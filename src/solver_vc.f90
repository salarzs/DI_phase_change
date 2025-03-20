module mod_solver_vc
  ! 
  use decomp_2d
  use mod_param     , only: dims,gacc
  use mod_common_mpi, only: coord,myid,ierr,comm_cart
  use mod_types     , only: rp
  !
  implicit none
  !
  integer, parameter :: HYPRESolverSMG      = 1, &
                        HYPRESolverPFMG     = 2, &
                        HYPRESolverGMRES    = 3, &
                        HYPRESolverBiCGSTAB = 4
  integer , parameter :: maxiter = 50
  real(rp), parameter :: tol = 1.0d-6, maxerror = tol
  !integer(8) ::  HYPRESolverType = HYPRESolverSMG  ! less efficient, but more robust
  integer(8) ::  HYPRESolverType = HYPRESolverPFMG ! more efficient, but less robust
  integer(8) ::  grid,stencil,solver,mat,vec,x,precond,precond_id
  !
  private
  public init_solver_mg,solver_mg,destroy_solver_mg
  !
  contains
  !
  subroutine init_solver_mg(n,dli,dzci,cbc,bc)
    !
    ! Helmholtz/Poisson solver for a Poisson equation with variable
    ! coefficients
    !
    implicit none
    !
    integer         , intent(in), dimension(3)     :: n
    real(rp)        , intent(in), dimension(3)     :: dli
    real(rp)        , intent(in), dimension(-2:)   :: dzci
    character(len=1), intent(in), dimension(0:1,3) :: cbc
    real(rp)        , intent(in), dimension(0:1,3) :: bc
    !
    ! hypre solver parameters/variables
    !
    real(rp) :: final_res_norm,tol
    integer :: extra,local_size
    integer :: mpi_comm
    !
    integer, parameter :: ndims = 3, nstencil = 7
    integer, dimension(ndims) :: ijlower, ijupper
    integer, dimension(ndims) :: periods
    integer, dimension(ndims,nstencil) :: offsets
    integer, dimension(nstencil) :: stencil_indices
    integer :: nvalues,num_iterations
    integer :: q,qq
    integer, dimension(3) :: ng
    !
    ng(1:2) = n(1:2)*dims(1:2)
    ng(3)   = n(3)
    periods(:)  = 0
    do q=1,3
      do qq=0,1
        select case(cbc(qq,q))
        case('N')
        case('D')
        case('P')
          periods(q) = ng(q)
        end select
      enddo
    enddo
    mpi_comm = comm_cart
    !
    ! initialize matrix
    !
    ! 1 - setup the grid
    !
    !   1.1 - create empty 2D grid object
    !
    call HYPRE_StructGridCreate(mpi_comm,ndims,grid,ierr)
    call HYPRE_StructGridSetPeriodic(grid,periods,ierr)
    !
    !   1.2 - set grid extents and assemble
    !
    ijlower(1:2) = (coord(:)  )*n(1:2)+1
    ijupper(1:2) = (coord(:)+1)*n(1:2)
    ijlower(3  ) = 1
    ijupper(3  ) = n(3)
    call HYPRE_StructGridSetExtents(grid,ijlower,ijupper,ierr)
    call HYPRE_StructGridAssemble(  grid,ierr)
    !
    ! 2 - setup the finite-difference stencil
    !
    call HYPRE_StructStencilCreate(ndims,nstencil,stencil,ierr)
    offsets = reshape((/ 0, 0, 0, &
                        -1, 0, 0, &
                         1, 0, 0, &
                         0,-1, 0, &
                         0, 1, 0, &
                         0, 0,-1, &
                         0, 0, 1 /),shape(offsets))
    do q = 1,nstencil
      call HYPRE_StructStencilSetElement(stencil,q-1,offsets(:,q),ierr)
    enddo
    !
    ! 3 - create coefficient matrix and rhs vector
    !
    call HYPRE_StructMatrixCreate(mpi_comm,grid,stencil,mat,ierr)
    call HYPRE_StructMatrixSetSymmetric(mat,1,ierr)                ! set the symmetric property of the matrix 
    call HYPRE_StructMatrixInitialize(mat,ierr)
    call HYPRE_StructVectorCreate(mpi_comm,grid,vec,ierr)
    call HYPRE_StructVectorInitialize(vec,ierr)
    call HYPRE_StructVectorCreate(mpi_comm,grid,x,ierr)
    call HYPRE_StructVectorInitialize(x,ierr)
    !
    ! setup solver, and solve
    ! note: this part was taken from the Paris Simulator code
    !
    if ( HYPRESolverType .eq. HYPRESolverSMG ) then 
      call HYPRE_StructSMGCreate(mpi_comm, solver, ierr)
      call HYPRE_StructSMGSetMaxIter(solver, maxiter, ierr)
      call HYPRE_StructSMGSetTol(solver, MaxError, ierr)
      call hypre_structSMGsetLogging(solver, 1, ierr)
      call HYPRE_StructSMGSetPrintLevel(solver,1,ierr) 
    elseif ( HYPRESolverType .eq. HYPRESolverPFMG ) then  
      call HYPRE_StructPFMGCreate(mpi_comm, solver, ierr)
      call HYPRE_StructPFMGSetMaxIter(solver, maxiter, ierr)
      call HYPRE_StructPFMGSetTol(solver, MaxError, ierr)
      call HYPRE_structPFMGsetLogging(solver, 1, ierr)
      call HYPRE_StructPFMGSetPrintLevel(solver,1,ierr) 
      call HYPRE_StructPFMGSetRelChange(solver, 1, ierr) 
      ! Relaxiation Method: 2 is the fastest if symm matrix 
      ! 0: Jacobi
      ! 1: Weighted Jacobi (default)
      ! 2: Red/Black Gauss-Seidel (symmetric: RB pre- and post-relaxation)
      ! 3: Red/Black Gauss-Seidel (nonsymmetric: RB pre- and post-relaxation)
      call HYPRE_StructPFMGSetRelaxType(solver,2, ierr) 
      call HYPRE_StructPFMGSetNumPreRelax(solver,2,ierr)
      call HYPRE_StructPFMGSetNumPostRelax(solver,2,ierr)
    elseif ( HYPRESolverType .eq. HYPRESolverGMRES .or. & 
             HYPRESolverType .eq. HYPRESolverBiCGSTAB   ) then
      if (HYPRESolverType .eq. HYPRESolverGMRES) then 
        call HYPRE_StructGMRESCreate(mpi_comm, solver,ierr)
        call HYPRE_StructGMRESSetMaxIter(solver, maxiter,ierr)
        call HYPRE_StructGMRESSetTol(solver, MaxError, ierr)
        !call HYPRE_StructGMRESSetLogging(solver, 1 ,ierr)
      elseif (HYPRESolverType .eq. HYPRESolverBiCGSTAB) then 
        !call HYPRE_StructBiCGSTABCreate(mpi_comm, solver,ierr)
        !call HYPRE_StructBiCGSTABSetMaxIter(solver, maxiter,ierr)
        !call HYPRE_StructBiCGSTABSetTol(solver, MaxError, ierr)
      endif
      ! Use PFMG as preconditioner
      call HYPRE_StructPFMGCreate(mpi_comm, precond, ierr)
      call HYPRE_StructPFMGSetMaxIter(precond,10, ierr)
      call HYPRE_StructPFMGSetTol(precond,0._rp,ierr)
      call HYPRE_StructPFMGSetZeroGuess(precond,ierr)
      call HYPRE_StructPFMGSetRelChange(precond,1,ierr) 
      call HYPRE_StructPFMGSetRelaxType(precond,2,ierr) 
      precond_id = 1   ! Set PFMG as preconditioner
      if (HYPRESolverType .eq. HYPRESolverGMRES) then 
        call HYPRE_StructGMRESSetPrecond(solver,precond_id,precond,ierr)
      elseif (HYPRESolverType .eq. HYPRESolverBiCGSTAB) then 
        call HYPRE_StructBiCGSTABSetPrecond(solver,precond_id,precond,ierr)
      endif
    endif
    !
    return 
  end subroutine init_solver_mg
  !
  subroutine solver_mg(n,dli,dzci,cbc,bc,alpha,csound2,dt,pg,p,is_poisson)
    !
    ! Helmholtz/Poisson solver for a Poisson equation with variable
    ! coeffcients
    !
    implicit none
    !
    integer         , intent(in   ), dimension(3)        :: n
    real(rp)        , intent(in   ), dimension(3)        :: dli
    real(rp)        , intent(in   ), dimension(-2:)      :: dzci
    character(len=1), intent(in   ), dimension(0:1,3)    :: cbc
    real(rp)        , intent(in   ), dimension(0:1,3)    :: bc
    real(rp)        , intent(in   ), dimension(-2:,-2:,-2:) :: alpha
    real(rp)        , intent(in   ), dimension(0:,0:,0:) :: csound2
    real(rp)        , intent(in   ), dimension(-2:,-2:,-2:) :: pg
    real(rp)        , intent(in   )                      :: dt
    real(rp)        , intent(inout), dimension(-2:,-2:,-2:) :: p
    logical         , intent(in   )                      :: is_poisson
    !
    integer , dimension(3) :: ng
    real(rp) :: cxp,cxm,cyp,cym,czp,czm,cc,rhs,cc_time
    real(rp) :: alphaxp,alphaxm,alphayp,alphaym,alphazp,alphazm
    real(rp), dimension(3) :: dl
    integer :: i,j,k,ii,jj,kk,q,qq
    real(rp), dimension(0:1,3) :: factor,sgn
    !
    ! hypre solver parameters/variables
    !
    real(rp) :: final_res_norm,tol
    integer :: extra,local_size
    integer :: mpi_comm
    !
    integer, parameter :: ndims = 3, nstencil = 7
    integer, dimension(ndims) :: ijlower, ijupper
    integer, dimension(ndims) :: periods
    integer, dimension(ndims,nstencil) :: offsets
    integer, dimension(nstencil) :: stencil_indices
    real(rp), allocatable, dimension(:) :: matvalues,vecvalues, &
                                          guessvalues
    integer :: nvalues,num_iterations
    !
    dl(:)=1._rp/dli(:)
    !
    ! decide if to solve a Poisson or an Helmholtz equation
    !
    if(is_poisson) then
      cc_time = 0._rp
    else
      cc_time = 1._rp
    endif 
    !
    ng(1:2) = n(1:2)*dims(1:2)
    ng(3)   = n(3)
    factor(:,:) = 0._rp
    sgn(   :,:) = 0._rp
    periods(:)  = 0
    do q=1,3
      do qq=0,1
        select case(cbc(qq,q))
        case('N')
          factor(qq,q) = 1._rp/dli(q)*bc(qq,q)
          sgn(   qq,q) = 1._rp
        case('D')
          factor(qq,q) = -2._rp*bc(qq,q)
          sgn(   qq,q) = -1._rp
        case('P')
          factor(qq,q) = 0._rp
          sgn(   qq,q) = 0._rp
          periods(q)   = ng(q)
        end select
      enddo
    enddo
    mpi_comm = comm_cart
    !
    ijlower(1:2) = (coord(:)  )*n(1:2)+1
    ijupper(1:2) = (coord(:)+1)*n(1:2)
    ijlower(3  ) = 1
    ijupper(3  ) = n(3)
    !
    ! preliminaries for setting up the coefficient matrix
    !
    stencil_indices = (/0,1,2,3,4,5,6/)
    nvalues = product(n(:))*nstencil ! number of grid points*number of stencil entries
    allocate(matvalues(nvalues))
    matvalues(:) = 0._rp
    nvalues = product(n(:))
    allocate(vecvalues(nvalues))
    vecvalues(:) = 0._rp
    allocate(guessvalues(nvalues))
    guessvalues(:) = 0._rp
    !
    ! compute stencil coefficients and rhs
    !
    q = 0
    do k=1,n(3)
      kk = k
      do j=1,n(2)
        jj = j+coord(2)*n(2)
        do i=1,n(1)
          ii = i+coord(1)*n(1)
          q = q + 1
          !
#ifdef TWOD
          alphaxp = 0._rp
          alphaxm = 0._rp
#else
          alphaxp = (0.5_rp*(alpha(i+1,j,k)+alpha(i,j,k)))**(-1)
          alphaxm = (0.5_rp*(alpha(i-1,j,k)+alpha(i,j,k)))**(-1)
#endif
          alphayp = (0.5_rp*(alpha(i,j+1,k)+alpha(i,j,k)))**(-1)
          alphaym = (0.5_rp*(alpha(i,j-1,k)+alpha(i,j,k)))**(-1)
          alphazp = (0.5_rp*(alpha(i,j,k+1)+alpha(i,j,k)))**(-1)
          alphazm = (0.5_rp*(alpha(i,j,k-1)+alpha(i,j,k)))**(-1)
          !
          cxm = -alphaxm*dli(1)**2
          cxp = -alphaxp*dli(1)**2
          cym = -alphaym*dli(2)**2
          cyp = -alphayp*dli(2)**2
          czm = -alphazm*dli(3)**2 
          czp = -alphazp*dli(3)**2 
          cc  = -(cxm+cxp+cym+cyp+czm+czp) + cc_time/(alpha(i,j,k)*csound2(i,j,k)*(dt**2))
          rhs = p(i,j,k)/(alpha(i,j,k)*csound2(i,j,k)*(dt**2))
          !
#ifndef TWOD
          if(periods(1).eq.0) then
            if(    ii.eq.  1) then
              rhs = rhs + cxm*factor(0,1)
              cc = cc + sgn(0,1)*cxm
              cxm = 0.d0
            elseif(ii.eq.ng(1)) then
              rhs = rhs + cxp*factor(1,1)
              cc = cc + sgn(1,1)*cxp
              cxp = 0.d0
            endif
          endif
#endif
          if(periods(2).eq.0) then
            if(    jj.eq.  1 ) then
              rhs = rhs + cym*factor(0,2)
              cc = cc + sgn(0,2)*cym
              cym = 0.d0
            elseif(jj.eq.ng(2)) then
              rhs = rhs + cyp*factor(1,2)
              cc = cc + sgn(1,2)*cyp
              cyp = 0.d0
            endif
          endif
          if(periods(3).eq.0) then
            if(    kk.eq.  1) then
              rhs = rhs + czm*factor(0,3)
              cc = cc + sgn(0,3)*czm
              czm = 0.d0
            elseif(kk.eq.ng(3)) then
              rhs = rhs + czp*factor(1,3)
              cc = cc + sgn(1,3)*czp
              czp = 0.d0
            endif
          endif
          !
          matvalues((q-1)*nstencil+1) = cc
          matvalues((q-1)*nstencil+2) = cxm
          matvalues((q-1)*nstencil+3) = cxp
          matvalues((q-1)*nstencil+4) = cym
          matvalues((q-1)*nstencil+5) = cyp
          matvalues((q-1)*nstencil+6) = czm
          matvalues((q-1)*nstencil+7) = czp
          vecvalues(q               ) = rhs
          guessvalues(q             ) = pg(i,j,k)
          !if (i.eq.32.and.j.eq.32) print*,k,cym,cyp,czm,czp,rhs,pg(i,j,k)          
          !
        enddo
      enddo
    enddo
    !
    call HYPRE_StructMatrixSetBoxValues(mat,ijlower,ijupper,nstencil, &
                                        stencil_indices,matvalues, &
                                        ierr)
    call HYPRE_StructMatrixAssemble(mat,ierr)
    call HYPRE_StructVectorSetBoxValues(vec,ijlower,ijupper, &
                                        vecvalues,ierr)
    call HYPRE_StructVectorAssemble(vec,ierr)
    !
    ! create soluction vector
    !
    call HYPRE_StructVectorSetBoxValues(x,ijlower,ijupper, &
                                        guessvalues,ierr)
    call HYPRE_StructVectorAssemble(x,ierr)
    deallocate(guessvalues)
    !
    ! setup solver, and solve
    !
    ! note: this part was based on the the Paris Simulator code
    !       freely available under a GPL license; see:
    !       http://www.ida.upmc.fr/~zaleski/paris/
    !
    if ( HYPRESolverType .eq. HYPRESolverSMG ) then 
      call HYPRE_StructSMGCreate(mpi_comm, solver,ierr)
      call HYPRE_StructSMGSetup(solver, mat, vec, x, ierr)
      call HYPRE_StructSMGSolve(solver, mat, vec, x, ierr)
      call HYPRE_StructSMGDestroy(solver,ierr)
      !call HYPRE_StructSMGGetNumIterations(solver, num_iterations,ierr)
    elseif ( HYPRESolverType .eq. HYPRESolverPFMG ) then  
      call HYPRE_StructPFMGCreate(mpi_comm, solver,ierr)
      call HYPRE_StructPFMGSetup(solver, mat, vec, x,ierr)
      call HYPRE_StructPFMGSolve(solver, mat, vec, x,ierr)
      call HYPRE_StructPFMGDestroy(solver,ierr)
      call HYPRE_StructPFMGGetNumIteration(solver,num_iterations,ierr)
      if(myid.eq.0) print*, "# iteration:", num_iterations
    elseif (HYPRESolverType .eq. HYPRESolverGMRES) then 
      call HYPRE_StructGMRESCreate(mpi_comm, solver,ierr)
      call HYPRE_StructGMRESSetup(solver, mat, vec, x, ierr)
      call HYPRE_StructGMRESSolve(solver, mat, vec, x, ierr)
      call HYPRE_StructGMRESDestroy(solver,ierr)
      !call HYPRE_StructGMRESGetNumIterations(solver, num_iterations,ierr)
    elseif (HYPRESolverType .eq. HYPRESolverBiCGSTAB) then 
      call HYPRE_StructBiCGSTABCreate(solver, mat, vec, x, ierr)
      call HYPRE_StructBiCGSTABSetup(solver, mat, vec, x, ierr)
      call HYPRE_StructBiCGSTABSolve(solver, mat, vec, x, ierr)
      call HYPRE_StructBiCGSTABDestroy(solver, mat, vec, x, ierr)
      !call HYPRE_StructBiCGSTABGetNumIterations(solver, num_iterations,ierr)
    endif ! HYPRESolverType
    !
    ! end of part based on the Paris Simulator code
    !
    ! fecth results
    !
    call HYPRE_StructVectorGetBoxValues(x,ijlower,ijupper,vecvalues,ierr)
    q = 0
    p(:,:,:) = 0.d0
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          q = q + 1
          p(i,j,k) = vecvalues(q)
        enddo
      enddo
    enddo
    deallocate(vecvalues,matvalues)
    !
    return 
  end subroutine solver_mg
  !
  subroutine destroy_solver_mg
    !
    implicit none
    !
    ! note: this part was based on the the Paris Simulator code
    !       freely available under a GPL license; see:
    !       http://www.ida.upmc.fr/~zaleski/paris/
    !
    if ( HYPRESolverType .eq. HYPRESolverSMG ) then 
      !call HYPRE_StructSMGDestroy(solver, ierr)
    elseif ( HYPRESolverType .eq. HYPRESolverPFMG ) then  
      !call HYPRE_StructPFMGDestroy(solver, ierr)
    elseif ( HYPRESolverType .eq. HYPRESolverGMRES ) then  
      !call HYPRE_StructGMRESDestroy(solver, ierr)
      call HYPRE_StructPFMGDestroy(precond, ierr)
    elseif ( HYPRESolverType .eq. HYPRESolverBiCGSTAB ) then  
      !call HYPRE_StructBiCGSTABDestroy(solver, ierr)
      call HYPRE_StructPFMGDestroy(precond, ierr)
    endif ! HYPRESolverType
    !
    ! end of part based on the Paris Simulator code
    !
    call HYPRE_StructGridDestroy(grid,ierr)
    call HYPRE_StructStencilDestroy(stencil,ierr)
    call HYPRE_StructMatrixDestroy(mat,ierr)
    call HYPRE_StructVectorDestroy(vec,ierr)
    call HYPRE_StructVectorDestroy(x,ierr)
    !
    return 
  end subroutine destroy_solver_mg
  !
end module mod_solver_vc
