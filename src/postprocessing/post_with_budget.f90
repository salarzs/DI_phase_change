module mod_post 
  !
  use mod_types
  use mod_bound     , only: boundp
  use mod_common_mpi, only: myid,ierr,ijk_start, &
                            MPI_PROC_NULL,       &
                            top,bottom,          &
                            front,back,          &
                            left,right,          &
                            comm_cart
  use mpi
  use mod_param,     only: lx,ly,lz,&
                           dx,dy,dz,&
                           n,dims_in,sigma,&
                           datadir_bal,datadir_cur1,datadir_cur2,datadir_cur3,datadir_cur4

#if defined(_OPENACC)
  use cudafor
  use mod_common_mpi, only: mydev
#endif
  !
  implicit none
  !
  private
  public  :: compute_vorticity,mixed_variables,budget,energy_balance_mf,&
             curv_pdf1,curv_pdf2,curv_pdf3,curv_pdf4
#if defined(_HEAT_TRANSFER)
  public  :: wall_avg
#endif
#if defined(_USE_VOF)
  public  :: time_tw_avg
#else
  public  :: time_sp_avg
#endif
  !
  contains
  !
#if defined(_HEAT_TRANSFER)
  subroutine wall_avg(idir,nx,ny,nz,ngx,ngy,ngz,dxi,dyi,dzi,nh_t,ka,tmp,time)
    !
    ! note: --> to be generalized for non-uniform grid along z (streched grid)
    !
    use mod_param, only: deltaT
#if defined(_USE_VOF)
    use mod_param, only: kappa2
#else
    use mod_param, only: kappa_sp
#endif
    !
    implicit none
    !
    integer , intent(in )                                      :: idir
    integer , intent(in )                                      :: nx,ny,nz
    integer , intent(in )                                      :: ngx,ngy,ngz
    real(rp), intent(in )                                      :: dxi,dyi,dzi
    integer , intent(in )                                      :: nh_t
    real(rp), intent(in ), dimension(     0:,     0:,     0:)  :: ka
    real(rp), intent(in ), dimension(1-nh_t:,1-nh_t:,1-nh_t:)  :: tmp
    real(rp), intent(in )                                      :: time
    !
    real(rp) :: dtdz,dtdy,dtdx,kap
    real(rp) :: lref,ka_ref
    real(rp) :: nusselt_up,nusselt_down
    integer  :: i,j,k,ip,jp,kp,im,jm,km
#if defined(_OPENACC)
    integer :: istat
    attributes(managed) :: ka, tmp
#endif
    !
#if defined(_USE_VOF)
    ka_ref = kappa2
#else
    ka_ref = kappa_sp
    kap    = kappa_sp
#endif
    !
    select case(idir)
    !
    ! Z-ORIENTATION
    !
    case(3)
      !
      lref = lz
      nusselt_up = 0._rp
      !
      if(top.eq.MPI_PROC_NULL) then ! z - top boundary 
        k  = nz+1
        km = nz
        !$acc kernels
        do j=1,ny
          do i=1,nx
            dtdz = (tmp(i,j,k)-tmp(i,j,km))*dzi
#if defined(_USE_VOF)
            kap  = 0.5_rp*(ka(i,j,k)+ka(i,j,km))
#endif
            nusselt_up = nusselt_up + kap*dtdz
          enddo
        enddo
        !$acc end kernels 
      endif
      call mpi_allreduce(MPI_IN_PLACE,nusselt_up,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      nusselt_up = (nusselt_up/(1._rp*ngx*ngy))*(lref/ka_ref/deltaT)
      !
      nusselt_down = 0._rp
      !
      if(bottom.eq.MPI_PROC_NULL) then ! z - bottom boundary 
        k  = 0
        kp = 1
        !$acc kernels
        do j=1,ny 
          do i=1,nx
            dtdz = (tmp(i,j,kp)-tmp(i,j,k))*dzi
#if defined(_USE_VOF)
            kap  = 0.5_rp*(ka(i,j,k)+ka(i,j,kp))
#endif
            nusselt_down = nusselt_down - kap*dtdz
          enddo
        enddo
        !$acc end kernels 
      endif
      call mpi_allreduce(MPI_IN_PLACE,nusselt_down,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      nusselt_down = (nusselt_down/(1._rp*ngx*ngy))*(lref/ka_ref/deltaT)
      ! 
      if(myid.eq.0) then
        open(94,file='data/post/wall/nusselt.out',position='append')
        write(94,'(3E15.7)') time,nusselt_down,nusselt_up
        close(94)
      endif
    !
    ! Y-ORIENTATION
    !
    case(2)
      !
      lref = ly
      nusselt_up = 0._rp
      !
      if(back.eq.MPI_PROC_NULL) then ! y - top boundary 
        j  = ny+1
        jm = ny
        !$acc kernels
        do k=1,nz
          do i=1,nx
            dtdy = (tmp(i,j,k)-tmp(i,jm,k))*dyi
#if defined(_USE_VOF)
            kap  = 0.5_rp*(ka(i,j,k)+ka(i,jm,k))
#endif
            nusselt_up = nusselt_up + kap*dtdy
          enddo
        enddo
        !$acc end kernels 
      endif
      call mpi_allreduce(MPI_IN_PLACE,nusselt_up,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      nusselt_up = (nusselt_up/(1._rp*ngx*ngz))*(lref/ka_ref/deltaT)
      !
      nusselt_down = 0._rp
      !
      if(front.eq.MPI_PROC_NULL) then ! y - bottom boundary 
        j  = 0
        jp = 1
        !$acc kernels
        do k=1,nz 
          do i=1,nx
            dtdy = (tmp(i,jp,k)-tmp(i,j,k))*dyi
            kap  = 0.5_rp*(ka(i,j,k)+ka(i,jp,k))
            nusselt_down = nusselt_down - kap*dtdy
          enddo
        enddo
        !$acc end kernels 
      endif
      call mpi_allreduce(MPI_IN_PLACE,nusselt_down,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      nusselt_down = (nusselt_down/(1._rp*ngx*ngz))*(lref/ka_ref/deltaT)
      ! 
      if(myid.eq.0) then
        open(94,file='data/post/wall/nusselt.out',position='append')
        write(94,'(3E15.7)') time,nusselt_down,nusselt_up
        close(94)
      endif
    !
    ! X-ORIENTATION
    !
    case(1)
      !
      lref = lx
      nusselt_up = 0._rp
      !
      if(right.eq.MPI_PROC_NULL) then ! x - top boundary 
        i  = nx+1
        im = nx
        !$acc kernels
        do k=1,nz
          do j=1,ny
            dtdx = (tmp(i,j,k)-tmp(im,j,k))*dxi
            kap  = 0.5_rp*(ka(i,j,k)+ka(im,j,k))
            nusselt_up = nusselt_up + kap*dtdx
          enddo
        enddo
        !$acc end kernels 
      endif
      call mpi_allreduce(MPI_IN_PLACE,nusselt_up,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      nusselt_up = (nusselt_up/(1._rp*ngy*ngz))*(lref/ka_ref/deltaT)
      !
      nusselt_down = 0._rp
      !
      if(left.eq.MPI_PROC_NULL) then ! x - bottom boundary 
        i  = 0
        ip = 1
        !$acc kernels
        do k=1,nz 
          do j=1,ny
            dtdx = (tmp(ip,j,k)-tmp(i,j,k))*dxi
            kap  = 0.5_rp*(ka(i,j,k)+ka(ip,j,k))
            nusselt_down = nusselt_down - kap*dtdx
          enddo
        enddo
        !$acc end kernels 
      endif
      call mpi_allreduce(MPI_IN_PLACE,nusselt_down,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      nusselt_down = (nusselt_down/(1._rp*ngy*ngz))*(lref/ka_ref/deltaT)
      ! 
      if(myid.eq.0) then
        open(94,file='data/post/wall/nusselt.out',position='append')
        write(94,'(3E15.7)') time,nusselt_down,nusselt_up
        close(94)
      endif
    end select
    !
    return
  end subroutine wall_avg
#endif
  !
#if defined(_USE_VOF)
  subroutine time_tw_avg(idir,do_avg,do_favre,is_stag,fname,n,ng,istep,istep_av,iout1d, &
                         nh_d,nh_v,nh_p,psi,rho_p,p,pout1,pout2,pvol1,pvol2)
    !
    ! writes the profile of a variable averaged
    ! over two domain directions (TWO-PHASE VERSION)
    !
    ! idir     -> the direction normal to the averaging plane
    ! do_avg   -> do or not averaging
    ! do_favre -> Favre or regular averaging
    ! is_stag  -> decide if "p" is a cell-centered or staggered variable
    ! fname    -> name of the file
    ! n        -> size of the input array
    ! ng       -> total size of computational domain
    ! istep    -> current time step
    ! istep_av -> size of statistical sample
    ! iout1d   -> print file every iout1d time steps
    ! idir     -> direction of the profile (all other directions are averaged)
    ! nh_      -> halos point
    ! psi      -> 3D vof field (0 --> liquid, 1 --> gas)
    ! rho_p    -> density of the phase (for compressible phase)
    ! p        -> 3D input scalar field to be averaged
    ! pout1    ->  first order time statistics of plane averaged field (mean)
    ! pout2    -> second order time statistics of plane averaged field (rms)
    ! pvol     ->  first order time statistics of volume averaged field (mean)
    ! pvol2    -> second order time statistics of volume averaged field (rms)
    !
    ! note: --> to be generalized for non-uniform grid along z (streched grid)
    !
    use mod_param, only: dl
    !
    implicit none
    !
    integer         , intent(in   )                                     :: idir
    logical         , intent(in   )                                     :: do_avg,do_favre
    integer         , intent(in   ), dimension(3)                       :: is_stag
    character(len=*), intent(in   )                                     :: fname
    integer         , intent(in   ), dimension(3)                       :: n,ng
    integer         , intent(in   )                                     :: istep,istep_av,iout1d
    integer         , intent(in   )                                     :: nh_d,nh_v,nh_p
    real(rp)        , intent(in   ), dimension(1-nh_v:,1-nh_v:,1-nh_v:) :: psi
    real(rp)        , intent(in   ), dimension(0     :,0     :,0     :) :: rho_p
    real(rp)        , intent(in   ), dimension(1-nh_p:,1-nh_p:,1-nh_p:) :: p
    real(rp)        , intent(inout), dimension(1:)                      :: pout1,pout2
    real(rp)        , intent(inout)                                     :: pvol1,pvol2
    !
    real(rp), allocatable, dimension(:) :: p1d1,p1d2,rhop
    real(rp) :: factor,factor2,p_dl_idir,p_avg
    integer  :: i,j,k,mm,qx,qy,qz
    integer  :: nx,ny,nz,ngx,ngy,ngz,ng_idir
    integer  :: start
    integer  :: iunit
#if defined(_OPENACC)
    integer :: istat
    attributes(managed) :: p,psi,rho_p,p1d1,p1d2,rhop,pout1,pout2
#endif
    !
    nx      = n(1)
    ny      = n(2)
    nz      = n(3)
    ngx     = ng(1)
    ngy     = ng(2)
    ngz     = ng(3)
    ng_idir = ng(idir)
    start   = ijk_start(idir)
    !
    qx = is_stag(1)
    qy = is_stag(2)
    qz = is_stag(3) 
    !
    allocate(p1d1(ng(idir)),p1d2(ng(idir)),rhop(ng(idir)))
    !
    iunit     = 10
    factor    = 1._rp*istep_av
    factor2   = 1._rp*ng_idir/(1._rp*ngx*ngy*ngz)
    p_dl_idir = dl(1)*dl(2)*dl(2)/dl(idir)
    !
    if(istep_av.eq.1) then  
      do mm=1,ng_idir
        pout1(mm) = 0._rp
        pout2(mm) = 0._rp
      enddo
    endif
    !
    !$acc kernels
    do mm=1,ng_idir
      p1d1(mm)  = 0._rp
      p1d2(mm)  = 0._rp
      rhop(mm)  = 0._rp
    enddo
    !$acc end kernels 
    !
    if(do_favre) then
      !
      ! Density-based averaging (Favre)
      !
      select case(idir)
      case(1)
        !$acc kernels
        do k=1,nz
          do j=1,ny
            do i=1,nx
              mm = start + i
              !
              p_avg    = 0.5_rp*(p(i,j,k)+p(i-qx,j-qy,k-qz))
              rhop(mm) = rhop(mm) + p_dl_idir*rho_p(i,j,k)
              p1d1(mm) = p1d1(mm) + p_dl_idir*rho_p(i,j,k)*(1._rp-psi(i,j,k))*p_avg 
              p1d2(mm) = p1d2(mm) + p_dl_idir*rho_p(i,j,k)*(1._rp-psi(i,j,k))*p_avg**2 
              !
            enddo
          enddo
        enddo
        !$acc end kernels 
      case(2)
        !$acc kernels
        do k=1,nz
          do j=1,ny
            do i=1,nx
              mm = start + j
              !
              p_avg    = 0.5_rp*(p(i,j,k)+p(i-qx,j-qy,k-qz))
              rhop(mm) = rhop(mm) + p_dl_idir*rho_p(i,j,k)
              p1d1(mm) = p1d1(mm) + p_dl_idir*rho_p(i,j,k)*(1._rp-psi(i,j,k))*p_avg 
              p1d2(mm) = p1d2(mm) + p_dl_idir*rho_p(i,j,k)*(1._rp-psi(i,j,k))*p_avg**2 
              !
            enddo
          enddo
        enddo
        !$acc end kernels 
      case(3)
        !$acc kernels
        do k=1,nz
          do j=1,ny
            do i=1,nx
              mm = start + k
              !
              p_avg    = 0.5_rp*(p(i,j,k)+p(i-qx,j-qy,k-qz))
              rhop(mm) = rhop(mm) + p_dl_idir*rho_p(i,j,k)
              p1d1(mm) = p1d1(mm) + p_dl_idir*rho_p(i,j,k)*(1._rp-psi(i,j,k))*p_avg 
              p1d2(mm) = p1d2(mm) + p_dl_idir*rho_p(i,j,k)*(1._rp-psi(i,j,k))*p_avg**2 
              !
            enddo
          enddo
        enddo
        !$acc end kernels 
      end select
      !
      call mpi_allreduce(MPI_IN_PLACE,rhop(1),ng_idir,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,p1d1(1),ng_idir,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,p1d2(1),ng_idir,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      !
      !$acc kernels
      do mm=1,ng_idir
        p1d1(mm) = p1d1(mm)/(rhop(mm))
        p1d2(mm) = p1d2(mm)/(rhop(mm))
      enddo
      !$acc end kernels 
      !
    else
      !
      ! Volumetric-based averaging (no Favre)
      !
      select case(idir)
      case(1)
        !$acc kernels
        do k=1,nz
          do j=1,ny
            do i=1,nx
              mm = start + i
              !
              p_avg    = 0.5_rp*(p(i,j,k)+p(i-qx,j-qy,k-qz))
              p1d1(mm) = p1d1(mm) + (1._rp-psi(i,j,k))*p_avg
              p1d2(mm) = p1d2(mm) + (1._rp-psi(i,j,k))*p_avg**2
              !
            enddo
          enddo
        enddo
        !$acc end kernels 
      case(2)
        !$acc kernels
        do k=1,nz
          do j=1,ny
            do i=1,nx
              mm = start + j
              !
              p_avg    = 0.5_rp*(p(i,j,k)+p(i-qx,j-qy,k-qz))
              p1d1(mm) = p1d1(mm) + (1._rp-psi(i,j,k))*p_avg
              p1d2(mm) = p1d2(mm) + (1._rp-psi(i,j,k))*p_avg**2
              !
            enddo
          enddo
        enddo
        !$acc end kernels 
      case(3)
        !$acc kernels
        do k=1,nz
          do j=1,ny
            do i=1,nx
              mm = start + k
              !
              p_avg    = 0.5_rp*(p(i,j,k)+p(i-qx,j-qy,k-qz))
              p1d1(mm) = p1d1(mm) + (1._rp-psi(i,j,k))*p_avg
              p1d2(mm) = p1d2(mm) + (1._rp-psi(i,j,k))*p_avg**2
              !
            enddo
          enddo
        enddo
        !$acc end kernels 
      end select
      !
      call mpi_allreduce(MPI_IN_PLACE,p1d1(1),ng_idir,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,p1d2(1),ng_idir,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      !
      !$acc kernels
      do mm=1,ng_idir
        p1d1(mm) = p1d1(mm)*factor2
        p1d2(mm) = p1d2(mm)*factor2
      enddo
      !$acc end kernels 
      !
    endif
    !
    ! decide or not to averaging 
    !
    if(.not.do_avg) then
      !$acc kernels
      do mm=1,ng_idir
        pout1(mm) = p1d1(mm)
        pout2(mm) = p1d2(mm)
      enddo
      !$acc end kernels 
    else
      !$acc kernels
      do mm=1,ng_idir
        pout1(mm) = ((factor-1._rp)*pout1(mm)+p1d1(mm))/factor
        pout2(mm) = ((factor-1._rp)*pout2(mm)+p1d2(mm))/factor
      enddo
      !$acc end kernels 
    endif
    pvol1 = sum(pout1)
    pvol2 = sum(pout2)
    !
    ! print 
    !  note: we put this condition on iout1d in order to ensure that the 
    !        averaging frequency (iout_av) is indenpendent 
    !        of the print frequency of the files (iout1d)
    !
    if(mod(istep,iout1d).eq.0) then
      if(myid.eq.0) then
        open(unit=iunit,file=fname)
        do mm=1,ng_idir
          write(iunit,'(3E15.7)') (mm-0.5_rp)*dl(idir),pout1(mm),pout2(mm)
        enddo
        close(iunit)
      endif
      pout1 = 0._rp
      pout2 = 0._rp
      pvol1 = 0._rp
      pvol2 = 0._rp
    endif
    deallocate(p1d1,p1d2,rhop)
    !
    return
  end subroutine time_tw_avg
#else
  subroutine time_sp_avg(idir,do_avg,do_favre,is_stag,fname,n,ng,istep,istep_av,iout1d, &
                         nh_d,nh_p,rho_p,p,pout1,pout2,pvol1,pvol2)
    !
    ! writes the profile of a variable averaged
    ! over two domain directions (SINGLE-PHASE VERSION)
    !
    ! idir     -> the direction normal to the averaging plane
    ! do_avg   -> do or not averaging
    ! do_favre -> Favre or regular averaging
    ! is_stag  -> decide if "p" is a cell-centered or staggered variable
    ! fname    -> name of the file
    ! n        -> size of the input array
    ! ng       -> total size of computational domain
    ! istep    -> current time step
    ! istep_av -> size of statistical sample
    ! iout1d   -> print file every iout1d time steps
    ! idir     -> direction of the profile (all other directions are averaged)
    ! nh_      -> halos point
    ! rho_p    -> density of the phase (for compressible phase)
    ! p        -> 3D input scalar field to be averaged
    ! pout1    ->  first order time statistics of plane averaged field (mean)
    ! pout2    -> second order time statistics of plane averaged field (rms)
    ! pvol     ->  first order time statistics of volume averaged field (mean)
    ! pvol2    -> second order time statistics of volume averaged field (rms)
    !
    ! note: --> to be generalized for non-uniform grid along z (streched grid)
    !
    use mod_param, only: dl
    !
    implicit none
    !
    integer         , intent(in   )                                     :: idir
    logical         , intent(in   )                                     :: do_avg,do_favre
    integer         , intent(in   ), dimension(3)                       :: is_stag
    character(len=*), intent(in   )                                     :: fname
    integer         , intent(in   ), dimension(3)                       :: n,ng
    integer         , intent(in   )                                     :: istep,istep_av,iout1d
    integer         , intent(in   )                                     :: nh_d,nh_p
    real(rp)        , intent(in   ), dimension(0     :,0     :,0     :) :: rho_p
    real(rp)        , intent(in   ), dimension(1-nh_p:,1-nh_p:,1-nh_p:) :: p
    real(rp)        , intent(inout), dimension(1:)                      :: pout1,pout2
    real(rp)        , intent(inout)                                     :: pvol1,pvol2
    !
    real(rp), allocatable, dimension(:) :: p1d1,p1d2,rhop
    real(rp) :: factor,factor2,p_dl_idir,p_avg
    integer  :: i,j,k,mm,qx,qy,qz
    integer  :: nx,ny,nz,ngx,ngy,ngz,ng_idir
    integer  :: start
    integer  :: iunit
#if defined(_OPENACC)
    integer :: istat
    attributes(managed) :: p,rho_p,p1d1,p1d2,rhop,pout1,pout2
#endif
    !
    nx      = n(1)
    ny      = n(2)
    nz      = n(3)
    ngx     = ng(1)
    ngy     = ng(2)
    ngz     = ng(3)
    ng_idir = ng(idir)
    start   = ijk_start(idir)
    !
    qx = is_stag(1)
    qy = is_stag(2)
    qz = is_stag(3) 
    !
    allocate(p1d1(ng(idir)),p1d2(ng(idir)),rhop(ng(idir)))
    !
    iunit     = 10
    factor    = 1._rp*istep_av
    factor2   = 1._rp*ng_idir/(1._rp*ngx*ngy*ngz)
    p_dl_idir = dl(1)*dl(2)*dl(2)/dl(idir)
    !
    if(istep_av.eq.1) then  
      do mm=1,ng_idir
        pout1(mm) = 0._rp
        pout2(mm) = 0._rp
      enddo
    endif
    !
    !$acc kernels
    do mm=1,ng_idir
      p1d1(mm)  = 0._rp
      p1d2(mm)  = 0._rp
      rhop(mm)  = 0._rp
    enddo
    !$acc end kernels 
    !
    if(do_favre) then
      !
      ! Density-based averaging (Favre)
      !
      select case(idir)
      case(1)
        !$acc kernels
        do k=1,nz
          do j=1,ny
            do i=1,nx
              mm = start + i
              !
              p_avg    = 0.5_rp*(p(i,j,k)+p(i-qx,j-qy,k-qz))
              rhop(mm) = rhop(mm) + p_dl_idir*rho_p(i,j,k)
              p1d1(mm) = p1d1(mm) + p_dl_idir*rho_p(i,j,k)*p_avg 
              p1d2(mm) = p1d2(mm) + p_dl_idir*rho_p(i,j,k)*p_avg**2 
              !
            enddo
          enddo
        enddo
        !$acc end kernels 
      case(2)
        !$acc kernels
        do k=1,nz
          do j=1,ny
            do i=1,nx
              mm = start + j
              !
              p_avg    = 0.5_rp*(p(i,j,k)+p(i-qx,j-qy,k-qz))
              rhop(mm) = rhop(mm) + p_dl_idir*rho_p(i,j,k)
              p1d1(mm) = p1d1(mm) + p_dl_idir*rho_p(i,j,k)*p_avg 
              p1d2(mm) = p1d2(mm) + p_dl_idir*rho_p(i,j,k)*p_avg**2 
              !
            enddo
          enddo
        enddo
        !$acc end kernels 
      case(3)
        !$acc kernels
        do k=1,nz
          do j=1,ny
            do i=1,nx
              mm = start + k
              !
              p_avg    = 0.5_rp*(p(i,j,k)+p(i-qx,j-qy,k-qz))
              rhop(mm) = rhop(mm) + p_dl_idir*rho_p(i,j,k)
              p1d1(mm) = p1d1(mm) + p_dl_idir*rho_p(i,j,k)*p_avg 
              p1d2(mm) = p1d2(mm) + p_dl_idir*rho_p(i,j,k)*p_avg**2 
              !
            enddo
          enddo
        enddo
        !$acc end kernels 
      end select
      !
      call mpi_allreduce(MPI_IN_PLACE,rhop(1),ng_idir,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,p1d1(1),ng_idir,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,p1d2(1),ng_idir,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      !
      !$acc kernels
      do mm=1,ng_idir
        p1d1(mm) = p1d1(mm)/(rhop(mm))
        p1d2(mm) = p1d2(mm)/(rhop(mm))
      enddo
      !$acc end kernels 
      !
    else
      !
      ! Volumetric-based averaging (no Favre)
      !
      select case(idir)
      case(1)
        !$acc kernels
        do k=1,nz
          do j=1,ny
            do i=1,nx
              mm = start + i
              !
              p_avg    = 0.5_rp*(p(i,j,k)+p(i-qx,j-qy,k-qz))
              p1d1(mm) = p1d1(mm) + p_avg
              p1d2(mm) = p1d2(mm) + p_avg**2
              !
            enddo
          enddo
        enddo
        !$acc end kernels 
      case(2)
        !$acc kernels
        do k=1,nz
          do j=1,ny
            do i=1,nx
              mm = start + j
              !
              p_avg    = 0.5_rp*(p(i,j,k)+p(i-qx,j-qy,k-qz))
              p1d1(mm) = p1d1(mm) + p_avg
              p1d2(mm) = p1d2(mm) + p_avg**2
              !
            enddo
          enddo
        enddo
        !$acc end kernels 
      case(3)
        !$acc kernels
        do k=1,nz
          do j=1,ny
            do i=1,nx
              mm = start + k
              !
              p_avg    = 0.5_rp*(p(i,j,k)+p(i-qx,j-qy,k-qz))
              p1d1(mm) = p1d1(mm) + p_avg
              p1d2(mm) = p1d2(mm) + p_avg**2
              !
            enddo
          enddo
        enddo
        !$acc end kernels 
      end select
      !
      call mpi_allreduce(MPI_IN_PLACE,p1d1(1),ng_idir,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,p1d2(1),ng_idir,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      !
      !$acc kernels
      do mm=1,ng_idir
        p1d1(mm) = p1d1(mm)*factor2
        p1d2(mm) = p1d2(mm)*factor2
      enddo
      !$acc end kernels 
      !
    endif
    !
    ! decide or not to averaging 
    !
    if(.not.do_avg) then
      !$acc kernels
      do mm=1,ng_idir
        pout1(mm) = p1d1(mm)
        pout2(mm) = p1d2(mm)
      enddo
      !$acc end kernels 
    else
      !$acc kernels
      do mm=1,ng_idir
        pout1(mm) = ((factor-1._rp)*pout1(mm)+p1d1(mm))/factor
        pout2(mm) = ((factor-1._rp)*pout2(mm)+p1d2(mm))/factor
      enddo
      !$acc end kernels 
    endif
    pvol1 = sum(pout1)
    pvol2 = sum(pout2)
    !
    ! print 
    !  note: we put this condition on iout1d in order to ensure that the 
    !        averaging frequency (iout_av) is indenpendent 
    !        of the print frequency of the files (iout1d)
    !
    if(mod(istep,iout1d).eq.0) then
      if(myid.eq.0) then
        open(unit=iunit,file=fname)
        do mm=1,ng_idir
          write(iunit,'(3E15.7)') (mm-0.5_rp)*dl(idir),pout1(mm),pout2(mm)
        enddo
        close(iunit)
      endif
      pout1 = 0._rp
      pout2 = 0._rp
      pvol1 = 0._rp
      pvol2 = 0._rp
    endif
    deallocate(p1d1,p1d2,rhop)
    !
    return
  end subroutine time_sp_avg
#endif
  !
  subroutine compute_vorticity(nx,ny,nz,dxi,dyi,dzi,nh_u,v,w,vor)
    !
    implicit none
    !
    integer , intent(in )                                     :: nx,ny,nz
    real(rp), intent(in )                                     :: dxi,dyi,dzi
    integer , intent(in )                                     :: nh_u
    real(rp), intent(in ), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: v,w
    real(rp), intent(out), dimension(     1:,     1:,     1:) :: vor
    !
    real(rp) :: vp,vm,wp,wm
    integer  :: i,j,k,im,jm,km,ip,jp,kp
#if defined(_OPENACC)
    integer :: istat
    attributes(managed) :: v, w
#endif
    !
    !$acc kernels
    do k=1,nz
      do j=1,ny
        do i=1,nx
          kp = k+1
          km = k-1
          jp = j+1
          jm = j-1
          ip = i+1
          im = i-1
          !
          vp = 0.25*(v(i,j,k)+v(i,jm,k)+v(i,j,kp)+v(i,jm,kp))
          vm = 0.25*(v(i,j,k)+v(i,jm,k)+v(i,j,km)+v(i,jm,km))
          wp = 0.25*(w(i,j,k)+w(i,j,km)+w(i,jp,k)+w(i,jp,km))
          wm = 0.25*(w(i,j,k)+w(i,j,km)+w(i,jm,k)+w(i,jm,km))
          !
          vor(i,j,k) = (wp-wm)*dyi-(vp-vm)*dzi
          !
        enddo
      enddo
    enddo
    !$acc end kernels 
    !
    return
  end subroutine compute_vorticity
  !
  subroutine mixed_variables(nx,ny,nz,dxi,dyi,dzi,nh_u,nh_s1, & 
                             u,v,w,s1, &
                             us1,vs1,ws1,uv,vw,wu)
    !
    implicit none
    !
    integer , intent(in )                                        :: nx,ny,nz
    real(rp), intent(in )                                        :: dxi,dyi,dzi
    integer , intent(in )                                        :: nh_u,nh_s1!,nh_s2
    real(rp), intent(in ), dimension(1-nh_u :,1-nh_u :,1-nh_u :) :: u,v,w
    real(rp), intent(in ), dimension(1-nh_s1:,1-nh_s1:,1-nh_s1:) :: s1 ! generic scalar
    real(rp), intent(out), dimension(      1:,      1:,      1:) :: us1,vs1,ws1
    real(rp), intent(out), dimension(      1:,      1:,      1:) :: uv ,vw ,wu
    !
    integer :: i,j,k,im,jm,km
#if defined(_OPENACC)
    integer :: istat
    attributes(managed) :: u,v,w,s1
#endif
    !
    !$acc kernels
    do k=1,nz
      do j=1,ny
        do i=1,nx
          km = k-1
          jm = j-1
          im = i-1
          !
          uv(i,j,k)  = 0.25_rp*(u(i,j,k)+u(im,j,k))*(v(i,j,k)+v(i,jm,k))
          vw(i,j,k)  = 0.25_rp*(v(i,j,k)+v(i,jm,k))*(w(i,j,k)+w(i,j,km))
          wu(i,j,k)  = 0.25_rp*(w(i,j,k)+w(i,j,km))*(u(i,j,k)+u(im,j,k))
          !
          us1(i,j,k) = 0.5_rp*(u(i,j,k)+u(im,j,k))*s1(i,j,k)
          vs1(i,j,k) = 0.5_rp*(v(i,j,k)+v(i,jm,k))*s1(i,j,k)
          ws1(i,j,k) = 0.5_rp*(w(i,j,k)+w(i,j,km))*s1(i,j,k)
          !
        enddo
      enddo
    enddo
    !$acc end kernels 
    !
    return
  end subroutine mixed_variables
  !
  subroutine budget(datadir,nx,ny,nz,ngx,ngy,ngz,rho1,rho2,dx,dy,dz, &
                    nh_u,u,v,w,rho,vof,time,istep)
    !
    implicit none
    !
    character(len=100), intent(in )                                     :: datadir
    integer           , intent(in )                                     :: nx ,ny ,nz
    integer           , intent(in )                                     :: ngx,ngy,ngz
    real(rp)          , intent(in )                                     :: rho1,rho2
    real(rp)          , intent(in )                                     :: dx,dy,dz
    integer           , intent(in )                                     :: nh_u
    real(rp)          , intent(in ), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w
    real(rp)          , intent(in ), dimension(     0:,     0:,     0:) :: rho
    real(rp)          , intent(in ), dimension(     0:,     0:,     0:) :: vof
    real(rp)          , intent(in )                                     :: time
    integer           , intent(in )                                     :: istep
    !
    real(rp) :: ke_t,ke_1,ke_2,dvol,volt,vol_t1,vol_t2
    integer  :: i,j,k,ip,jp,kp,im,jm,km
#if defined(_OPENACC)
    integer :: istat
    attributes(managed) :: u,v,w,rho,vof
#endif
    !
    real(rp), parameter :: eps = real(1e-12,rp)
    ! 
    volt = 1._rp*ngx*ngy*ngz
    dvol = dx*dy*dz
    !
    vol_t1 = 0._rp
    vol_t2 = 0._rp
    !
    ke_t = 0._rp
    ke_1 = 0._rp
    ke_2 = 0._rp
    !
    !$acc kernels
    do k=1,nz
      do j=1,ny
        do i=1,nx
          !
          ip = i + 1
          im = i - 1
          jp = j + 1
          jm = j - 1
          kp = k + 1
          km = k - 1
          !
          ! 1. volume
          ! 
          vol_t1 = vol_t1 + vof(i,j,k)
          vol_t2 = vol_t2 + (1._rp-vof(i,j,k))
          !
          ! 2. kinetic energy
          ! 
          ke_t = ke_t + rho(i,j,k)*&
          0.5_rp*(0.25_rp*(u(i,j,k)+u(im,j,k))**2 + 0.25_rp*(v(i,j,k)+v(i,jm,k))**2 + 0.25_rp*(w(i,j,k)+w(i,j,km))**2)
          ke_1 = ke_1 + rho1*vof(i,j,k)*&
          0.5_rp*(0.25_rp*(u(i,j,k)+u(im,j,k))**2 + 0.25_rp*(v(i,j,k)+v(i,jm,k))**2 + 0.25_rp*(w(i,j,k)+w(i,j,km))**2)
          ke_2 = ke_2 + rho2*(1._rp-vof(i,j,k))*&
          0.5_rp*(0.25_rp*(u(i,j,k)+u(im,j,k))**2 + 0.25_rp*(v(i,j,k)+v(i,jm,k))**2 + 0.25_rp*(w(i,j,k)+w(i,j,km))**2)
          !
        enddo
      enddo
    enddo
    !$acc end kernels 
    !
    call mpi_allreduce(MPI_IN_PLACE,vol_t1,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(MPI_IN_PLACE,vol_t2,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    !
    call mpi_allreduce(MPI_IN_PLACE,ke_t  ,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(MPI_IN_PLACE,ke_1  ,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(MPI_IN_PLACE,ke_2  ,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    !
    ke_t = ke_t/volt
    ke_1 = ke_1/(vol_t1+eps)
    ke_2 = ke_2/(vol_t2+eps)
    !
    if(myid.eq.0) then
      !
      ! a. post-processing one fluid
      !
      open(92,file=trim(datadir)//'ke_t.out',position='append')
      write(92,'(3E15.7)') 1._rp*istep,time,ke_t
      close(92)
      ! 
      ! b1. phase 1
      ! 
      open(93,file=trim(datadir)//'ke_1.out',position='append')
      write(93,'(4E15.7)') 1._rp*istep,time,vol_t1,ke_1
      close(93)
      ! 
      ! b2. phase 2
      !
      open(94,file=trim(datadir)//'ke_2.out',position='append')
      write(94,'(4E15.7)') 1._rp*istep,time,vol_t2,ke_2 
      close(94)
      !
    endif
    !
    return
  end subroutine budget
  ! 
  !Shahab
  subroutine energy_balance_mf(nh_d,nh_u,nh_v,dims,halo_v, &
                               dzc,dzf,time,istep,rho,mu,psi,kappa,p,u,v,w)
    !    
    ! energy balance in physical space
    !  note: --> for HIT (with minor modifications in the production term also HST);
    !        --> not yet on GPUs;
    !        --> it assumes that the production term is due to ABC forcing.
    !
    use mod_param, only: sigma,cbcpsi,bcpsi,rho1,mu1,rho2,mu2,lx,ly,lz,pi
#if defined(_TURB_FORCING)
    use mod_param, only: f0_t,abc_x, abc_y, abc_z ,k0_tx,k0_ty,k0_tz
#endif
    !
    implicit none
    !
    !integer         , intent(in  )                                     :: nx,ny,nz
    integer         , intent(in  )                                     :: nh_d,nh_u,nh_v
!    character(len=*), intent(in  )                                     :: datadir_bal
    integer         , intent(in  ), dimension(3)                       :: dims,halo_v
    real(rp)        , intent(in  ), dimension(1-nh_d:)                 :: dzc,dzf
    real(rp)        , intent(in  )                                     :: time
    integer         , intent(in  )                                     :: istep
    real(rp)        , intent(in  ), dimension(     0:,     0:,     0:) :: rho,mu,psi,kappa,p
    real(rp)        , intent(in  ), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w
    !
    real(rp), dimension(0:n(1)+1,0:n(2)+1,0:n(3)+1)   :: u_p,v_p,w_p
    real(rp), dimension(0:n(1)+1,0:n(2)+1,0:n(3)+1,6) :: S      
    real(rp),dimension(0:n(1)+1,0:n(2)+1,0:n(3)+1) :: surf_x, surf_y, surf_z
    real(rp) :: kappas
    integer  :: i,j,k,ip,im,jp,jm,kp,km,q             
    real(rp) :: dxi,dyi,dzi,vol
    real(rp) :: u_avg,v_avg,w_avg,S_S,E_E, &
                duTxx,duTxy,duTxz,dvTyx,dvTyy,dvTyz,dwTzx,dwTzy,dwTzz,t_prod
    real(rp) :: ur_t,vr_t,wr_t,ur1_v,vr1_v,wr1_v,ur2_v,vr2_v,wr2_v, & 
                um_t,vm_t,wm_t,um1_v,vm1_v,wm1_v,um2_v,vm2_v,wm2_v, & 
                eps_t ,eps1_v,eps2_v, & 
                prd_t ,prd1_v,prd2_v, & 
                surfi_t,surfp_t,surf1i_v,surf2i_v,surf1p_v,surf2p_v, &         
                t_surfp,t_surfi , &
                ke_t  ,ke1_v ,ke2_v , & 
                ens_t ,ens1_v,ens2_v, & 
                Psi_mf,Psi_nu,Psi_gf_v, &
                Tp1_v ,Tp2_v , &
                Tnu1_v,Tnu2_v, &
                vol1,vol2,lam,lambda,Relam,eta, area, ker, epsr
    real(rp) :: VGRux , VGRuy , VGRuz , &
                VGRvx , VGRvy , VGRvz , &
                VGRwx , VGRwy , VGRwz , &
                dupx  , dupy  , dupz  
    real(rp) :: fx_hit,fy_hit,fz_hit,xc,yc,zc
    !
    vol    = 1._rp*(n(1)*dims(1))*(n(2)*dims(2))*(n(3)*dims(3))
    dxi    = 1._rp/dx
    dyi    = 1._rp/dy
    dzi    = 1._rp/dz
    !
    E_E    = 0._rp 
    ens_t  = 0._rp
    ens1_v = 0._rp
    ens2_v = 0._rp
    !
    ur_t   = 0._rp 
    vr_t   = 0._rp
    wr_t   = 0._rp
    ur1_v  = 0._rp
    vr1_v  = 0._rp
    wr1_v  = 0._rp
    ur2_v  = 0._rp
    vr2_v  = 0._rp
    wr2_v  = 0._rp
    !
    um_t   = 0._rp 
    vm_t   = 0._rp
    wm_t   = 0._rp
    um1_v  = 0._rp
    vm1_v  = 0._rp
    wm1_v  = 0._rp
    um2_v  = 0._rp
    vm2_v  = 0._rp
    wm2_v  = 0._rp
    !
    epsr   = 0._rp
    eps_t  = 0._rp 
    eps1_v = 0._rp
    eps2_v = 0._rp
    !
    prd_t  = 0._rp
    prd1_v = 0._rp
    prd2_v = 0._rp
    !
    surfi_t= 0._rp
    surfp_t= 0._rp
    surf1i_v=0._rp
    surf2i_v=0._rp
    surf1p_v=0._rp
    surf2p_v=0._rp
    !
    ker    = 0._rp
    ke_t   = 0._rp
    ke1_v  = 0._rp
    ke2_v  = 0._rp
    !
    Psi_mf   = 0._rp 
    Psi_nu   = 0._rp
    Psi_gf_v = 0._rp
    !
    Tp1_v  = 0._rp 
    Tp2_v  = 0._rp
    !
    Tnu1_v = 0._rp 
    Tnu2_v = 0._rp
    !
    lam    = 0._rp 
    lambda = 0._rp
    !
    Relam  = 0._rp
    !
    eta    = 0._rp
    !
    area   = 0._rp
    !
    ! following Dodd and Ferrante, we compute (u_p,v_p,w_p)
    !
    u_avg = 0._rp
    v_avg = 0._rp
    w_avg = 0._rp
    !
    vol1 = 0._rp
    vol2 = 0._rp
    !
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          !
          ! average velocity
          !
          u_avg = u_avg + u(i,j,k)*dx*dy*dz
          v_avg = v_avg + v(i,j,k)*dx*dy*dz
          w_avg = w_avg + w(i,j,k)*dx*dy*dz
          !
          ! volume of the two phases
          !
          vol1 = vol1 + psi(i,j,k)
          vol2 = vol2 + (1._rp-psi(i,j,k))
          !
        enddo
      enddo
    enddo 
    call mpi_allreduce(MPI_IN_PLACE,u_avg,1,MPI_REAL_RP,mpi_sum,comm_cart,ierr)
    call mpi_allreduce(MPI_IN_PLACE,v_avg,1,MPI_REAL_RP,mpi_sum,comm_cart,ierr)
    call mpi_allreduce(MPI_IN_PLACE,w_avg,1,MPI_REAL_RP,mpi_sum,comm_cart,ierr)
    u_avg = u_avg/(lx*ly*lz)
    v_avg = v_avg/(lx*ly*lz)
    w_avg = w_avg/(lx*ly*lz)
    !
    call mpi_allreduce(MPI_IN_PLACE,vol1,1,MPI_REAL_RP,mpi_sum,comm_cart,ierr)
    call mpi_allreduce(MPI_IN_PLACE,vol2,1,MPI_REAL_RP,mpi_sum,comm_cart,ierr)
    vol1 = vol1 + 1.e-12 ! to avoid division by 0
    vol2 = vol2 + 1.e-12 ! to avoid division by 0
    !
    do k=0,n(3)+1
      do j=0,n(2)+1
        do i=0,n(1)+1
          u_p(i,j,k) = u(i,j,k)-u_avg
          v_p(i,j,k) = v(i,j,k)-v_avg
          w_p(i,j,k) = w(i,j,k)-w_avg
        enddo
      enddo
    enddo 
    !
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          !
          kp = k+1
          km = k-1
          jp = j+1
          jm = j-1
          ip = i+1
          im = i-1
          !
          ! a. MEAN KINETIC ENERGY BALANCE --> part 1
          !
          ! a1. calculation and storage (for the next loop) of the 
          !     strain rate tensor S_ij
          !
          ! along x
          !
          VGRux = ( u_p(i,j,k) - u_p(im,j,k) )*dxi
          VGRuy = 0.25_rp*( u_p(i,j,k) + u_p(i,jp,k) + u_p(im,j,k) + u_p(im,jp,k) )*dyi - &
                  0.25_rp*( u_p(i,j,k) + u_p(i,jm,k) + u_p(im,j,k) + u_p(im,jm,k) )*dyi  
          VGRuz = 0.25_rp*( u_p(i,j,k) + u_p(i,j,kp) + u_p(im,j,k) + u_p(im,j,kp) )*dzi - &
                  0.25_rp*( u_p(i,j,k) + u_p(i,j,km) + u_p(im,j,k) + u_p(im,j,km) )*dzi
          !
          ! along y
          !  
          VGRvx = 0.25_rp*( v_p(i,j,k) + v_p(ip,j,k) + v_p(i,jm,k) + v_p(ip,jm,k) )*dxi - &
                  0.25_rp*( v_p(i,j,k) + v_p(im,j,k) + v_p(i,jm,k) + v_p(im,jm,k) )*dxi
          VGRvy = ( v_p(i,j,k) - v_p(i,jm,k) )*dyi
          VGRvz = 0.25_rp*( v_p(i,j,k) + v_p(i,j,kp) + v_p(i,jm,k) + v_p(i,jm,kp) )*dzi - &
                  0.25_rp*( v_p(i,j,k) + v_p(i,j,km) + v_p(i,jm,k) + v_p(i,jm,km) )*dzi
          !
          ! along z
          ! 
          VGRwx = 0.25_rp*( w_p(i,j,k) + w_p(ip,j,k) + w_p(i,j,km) + w_p(ip,j,km) )*dxi - &
                  0.25_rp*( w_p(i,j,k) + w_p(im,j,k) + w_p(i,j,km) + w_p(im,j,km) )*dxi
          VGRwy = 0.25_rp*( w_p(i,j,k) + w_p(i,jp,k) + w_p(i,j,km) + w_p(i,jp,km) )*dyi - &
                  0.25_rp*( w_p(i,j,k) + w_p(i,jm,k) + w_p(i,j,km) + w_p(i,jm,km) )*dyi
          VGRwz = ( w_p(i,j,k) - w_p(i,j,km) )*dzi
          !
          ! incompressible part
          ! 
          !-----S11 = dU/dx
          S(i,j,k,1) = VGRux
          !-----S12 = 0.5*(dU/dy+dV/dx) = S21
          S(i,j,k,2) = 0.5_rp * (VGRuy + VGRvx)
          !-----S13 = 0.5*(dU/dz+dW/dx) = S31
          S(i,j,k,3) = 0.5_rp * (VGRuz + VGRwx)
          !-----S22 = dV/dy
          S(i,j,k,4) = VGRvy
          !-----S23 = 0.5*(dV/dz+dW/dy) = S32
          S(i,j,k,5) = 0.5_rp * (VGRvz + VGRwy)
          !-----S33 = dW/dz
          S(i,j,k,6) = VGRwz
          !
          ! a2. enstrophy
          !
          E_E    = (VGRwy-VGRvz)**2._rp + (VGRuz-VGRwx)**2._rp + (VGRvx-VGRuy)**2._rp 
          ens_t  = ens_t  + E_E
          ens1_v = ens1_v + E_E*psi(i,j,k)
          ens2_v = ens2_v + E_E*(1._rp-psi(i,j,k))
          !
        enddo
      enddo
    enddo
    ! 
    do q=1,6
      call boundp(cbcpsi,(/n(1),n(2),n(3)/),bcpsi,nh_d,nh_v,halo_v,(/dx,dy,dz/),dzc,dzf,S(:,:,:,q))
    enddo
    !
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          !
          kp = k+1
          km = k-1
          jp = j+1
          jm = j-1
          ip = i+1
          im = i-1
          !
          ! b1. RMS for one-fluid and two phases
          !
          ur_t = ur_t + (u_p(i,j,k))**2
          vr_t = vr_t + (v_p(i,j,k))**2
          wr_t = wr_t + (w_p(i,j,k))**2
          ! 
          ur1_v = ur1_v + ((u_p(i,j,k))**2)*psi(i,j,k)
          vr1_v = vr1_v + ((v_p(i,j,k))**2)*psi(i,j,k)
          wr1_v = wr1_v + ((w_p(i,j,k))**2)*psi(i,j,k)
          !
          ur2_v = ur2_v + ((u_p(i,j,k))**2)*(1._rp-psi(i,j,k))
          vr2_v = vr2_v + ((v_p(i,j,k))**2)*(1._rp-psi(i,j,k))
          wr2_v = wr2_v + ((w_p(i,j,k))**2)*(1._rp-psi(i,j,k))
          !
          ! b2. mean velocities for one-fluid and two phases
          !
          um_t = um_t + u_p(i,j,k)
          vm_t = vm_t + v_p(i,j,k)
          wm_t = wm_t + w_p(i,j,k)
          ! 
          um1_v = um1_v + u_p(i,j,k)*psi(i,j,k)
          vm1_v = vm1_v + v_p(i,j,k)*psi(i,j,k)
          wm1_v = wm1_v + w_p(i,j,k)*psi(i,j,k)
          !
          um2_v = um2_v + u_p(i,j,k)*(1._rp-psi(i,j,k))
          vm2_v = vm2_v + v_p(i,j,k)*(1._rp-psi(i,j,k))
          wm2_v = wm2_v + w_p(i,j,k)*(1._rp-psi(i,j,k))
          !
          ! c. MEAN KINETIC ENERGY BALANCE --> part 2
          !
          ! c1. dissipation
          !
          S_S =       S(i,j,k,1)*S(i,j,k,1) + & 
                2._rp*S(i,j,k,2)*S(i,j,k,2) + &
                2._rp*S(i,j,k,3)*S(i,j,k,3) + &
                      S(i,j,k,4)*S(i,j,k,4) + &
                2._rp*S(i,j,k,5)*S(i,j,k,5) + &
                      S(i,j,k,6)*S(i,j,k,6)
          !
          epsr   = 2._rp*S_S*(mu(i,j,k)/rho(i,j,k))
          eps_t  = eps_t  + epsr
          eps1_v = eps1_v + psi(i,j,k)*epsr
          eps2_v = eps2_v + (1._rp-psi(i,j,k))* epsr
          !
          ! c2. production
          !
#if defined(_TURB_FORCING)
          zc = (k+ijk_start(3)-0.5_rp)*dz/lz*2._rp*pi
          yc = (j+ijk_start(2)-0.5_rp)*dy/ly*2._rp*pi
          xc = (i+ijk_start(1)-0.5_rp)*dx/lx*2._rp*pi
          fx_hit = f0_t*(abc_x*sin(k0_tz*zc) + abc_z*cos(k0_ty*yc))
          fy_hit = f0_t*(abc_y*sin(k0_tx*xc) + abc_x*cos(k0_tz*zc))
          fz_hit = f0_t*(abc_z*sin(k0_ty*yc) + abc_y*cos(k0_tx*xc))
#else
          fx_hit = 0._rp
          fy_hit = 0._rp
          fz_hit = 0._rp
#endif
          t_prod = 0.5_rp*(u_p(i,j,k)+u_p(im,j,k))*fx_hit + &
                   0.5_rp*(v_p(i,j,k)+v_p(i,jm,k))*fy_hit + &
                   0.5_rp*(w_p(i,j,k)+w_p(i,j,km))*fz_hit
          prd_t  = prd_t  + t_prod
          prd1_v = prd1_v + t_prod*psi(i,j,k)
          prd2_v = prd2_v + t_prod*(1._rp-psi(i,j,k))
          !
          ! c3. surface tension work

          kappas = 0.5_rp*(kappa(ip,j,k)+kappa(i,j,k))
          surf_x(i,j,k) = dxi*sigma*kappas*(psi(ip,j,k)-psi(i,j,k))
          kappas = 0.5_rp*(kappa(i,jp,k)+kappa(i,j,k))
          surf_y(i,j,k) = dyi*sigma*kappas*(psi(i,jp,k)-psi(i,j,k)) 
          kappas = 0.5_rp*(kappa(i,j,kp)+kappa(i,j,k))
          surf_z(i,j,k) = dzi*sigma*kappas*(psi(i,j,kp)-psi(i,j,k)) 

          t_surfp = 0.5_rp*(u_p(i,j,k)+u_p(im,j,k))*surf_x(i,j,k) + &
                   0.5_rp*(v_p(i,j,k)+v_p(i,jm,k))*surf_y(i,j,k) + &
                   0.5_rp*(w_p(i,j,k)+w_p(i,j,km))*surf_z(i,j,k)
          surfp_t  = surfp_t  + t_surfp
          surf1p_v = surf1p_v + t_surfp*psi(i,j,k)
          surf2p_v = surf2p_v + t_surfp*(1._rp-psi(i,j,k))

          t_surfi = 0.5_rp*(u(i,j,k)+u(im,j,k))*surf_x(i,j,k) + &
                   0.5_rp*(v(i,j,k)+v(i,jm,k))*surf_y(i,j,k) + &
                   0.5_rp*(w(i,j,k)+w(i,j,km))*surf_z(i,j,k)
          surfi_t  = surfi_t  + t_surfi
          surf1i_v = surf1i_v + t_surfi*psi(i,j,k)
          surf2i_v = surf2i_v + t_surfi*(1._rp-psi(i,j,k))


          ! c3. kinetic energy
          !
          !ke_t  = ke_t  + rho(i,j,k)*&
          !0.5_rp*(0.25_rp*(u_p(i,j,k)+u_p(im,j,k))**2 + 0.25_rp*(v_p(i,j,k)+v_p(i,jm,k))**2 + 0.25_rp*(w_p(i,j,k)+w_p(i,j,km))**2)
          !ke1_v = ke1_v + rho1*psi(i,j,k)*&
          !0.5_rp*(0.25_rp*(u_p(i,j,k)+u_p(im,j,k))**2 + 0.25_rp*(v_p(i,j,k)+v_p(i,jm,k))**2 + 0.25_rp*(w_p(i,j,k)+w_p(i,j,km))**2)
          !ke2_v = ke2_v + rho2*(1._rp-psi(i,j,k))*&
          !0.5_rp*(0.25_rp*(u_p(i,j,k)+u_p(im,j,k))**2 + 0.25_rp*(v_p(i,j,k)+v_p(i,jm,k))**2 + 0.25_rp*(w_p(i,j,k)+w_p(i,j,km))**2)
          ker   =0.5_rp* &
          (0.25_rp*(u_p(i,j,k)+u_p(im,j,k))**2 + 0.25_rp*(v_p(i,j,k)+v_p(i,jm,k))**2 + 0.25_rp*(w_p(i,j,k)+w_p(i,j,km))**2)
          ke_t  = ke_t  + ker
          ke1_v = ke1_v + psi(i,j,k)* ker
          ke2_v = ke2_v + (1._rp-psi(i,j,k))* ker
 
!
          ! c4. external forces
          !
          Psi_nu = Psi_nu + sigma*kappa(i,j,k)*( &
                                  dxi*(0.5_rp*(psi(ip,j,k)-psi(im,j,k)))*0.5_rp*(u_p(i,j,k)+u_p(im,j,k)) + & 
                                  dyi*(0.5_rp*(psi(i,jp,k)-psi(i,jm,k)))*0.5_rp*(v_p(i,j,k)+v_p(i,jm,k)) + & 
                                  dzi*(0.5_rp*(psi(i,j,kp)-psi(i,j,km)))*0.5_rp*(w_p(i,j,k)+w_p(i,j,km))   &
                                               ) 
          !
          ! c7. pressure power (contribution per phase)
          !
          dupx = ( u_p(i,j,k)*0.5_rp*(p(ip,j,k)+p(i,j,k))-u_p(im,j,k)*0.5_rp*(p(i,j,k)+p(im,j,k)) )*dxi
          dupy = ( v_p(i,j,k)*0.5_rp*(p(i,jp,k)+p(i,j,k))-v_p(i,jm,k)*0.5_rp*(p(i,j,k)+p(i,jm,k)) )*dyi
          dupz = ( w_p(i,j,k)*0.5_rp*(p(i,j,kp)+p(i,j,k))-w_p(i,j,km)*0.5_rp*(p(i,j,k)+p(i,j,km)) )*dzi
          !
          Tp1_v = Tp1_v + (dupx+dupy+dupz)*psi(i,j,k) 
          Tp2_v = Tp2_v + (dupx+dupy+dupz)*(1._rp-psi(i,j,k))
          !
          ! c8. viscous power
          !
          duTxx = 0.500_rp*(u_p(i ,j,k)*(S(i,j,k,1)+S(ip,j,k,1))*(mu(i,j,k)+mu(ip,j,k))-&
                            u_p(im,j,k)*(S(i,j,k,1)+S(im,j,k,1))*(mu(i,j,k)+mu(im,j,k)))*dxi
          duTxy = 0.125_rp*(u_p(i,j,k)+u_p(i,jp,k)+u_p(im,j,k)+u_p(im,jp,k))*&
                           (S(i,jp,k,2)+S(i,j,k,2))*(mu(i,jp,k)+mu(i,j,k))*dyi-&
                  0.125_rp*(u_p(i,j,k)+u_p(i,jm,k)+u_p(im,j,k)+u_p(im,jm,k))*&
                           (S(i,jm,k,2)+S(i,j,k,2))*(mu(i,j,k)+mu(i,jm,k))*dyi  
          duTxz = 0.125_rp*(u_p(i,j,k)+u_p(i,j,kp)+u_p(im,j,k)+u_p(im,j,kp))*&
                           (S(i,j,kp,3)+S(i,j,k,3))*(mu(i,j,kp)+mu(i,j,k))*dzi-&
                  0.125_rp*(u_p(i,j,k)+u_p(i,j,km)+u_p(im,j,k)+u_p(im,j,km))*&
                           (S(i,j,k,3)+S(i,j,km,3))*(mu(i,j,k)+mu(i,j,km))*dzi
          !
          dvTyx = 0.125_rp*(v_p(i,j,k)+v_p(ip,j,k)+v_p(i,jm,k)+v_p(ip,jm,k))*&
                           (S(ip,j,k,2)+S(i,j,k,2))*(mu(ip,j,k)+mu(i,j,k))*dxi-&
                  0.125_rp*(v_p(i,j,k)+v_p(im,j,k)+v_p(i,jm,k)+v_p(im,jm,k))*&
                           ((S(i,j,k,2)+S(im,j,k,2))*(mu(i,j,k)+mu(im,j,k)))*dxi
          dvTyy = 0.500_rp*(v_p(i,j ,k)*(S(i,j,k,4)+S(i,jp,k,4))*(mu(i,j,k)+mu(i,jp,k))-&
                            v_p(i,jm,k)*(S(i,j,k,4)+S(i,jm,k,4))*(mu(i,j,k)+mu(i,jm,k)))*dyi
          dvTyz = 0.125_rp*(v_p(i,j,k)+v_p(i,j,kp)+v_p(i,jm,k)+v_p(i,jm,kp))*&
                           (S(i,j,k,5)+S(i,j,kp,5))*(mu(i,j,k)+mu(i,j,kp))*dzi-&
                  0.125_rp*(v_p(i,j,k)+v_p(i,j,km)+v_p(i,jm,k)+v_p(i,jm,km))*&
                           (S(i,j,k,5)+S(i,j,km,5))*(mu(i,j,k)+mu(i,j,km))*dzi
          !
          dwTzx = 0.125_rp*(w_p(i,j,k)+w_p(ip,j,k)+w_p(i,j,km)+w_p(ip,j,km))*&
                           (S(ip,j,k,3)+S(i,j,k,3))*(mu(ip,j,k)+mu(i,j,k))*dxi-&
                  0.125_rp*(w_p(i,j,k)+w_p(im,j,k)+w_p(i,j,km)+w_p(im,j,km))*& 
                           (S(i,j,k,3)+S(im,j,k,3))*(mu(i,j,k)+mu(im,j,k))*dxi
          dwTzy = 0.125_rp*(w_p(i,j,k)+w_p(i,jp,k)+w_p(i,j,km)+w_p(i,jp,km))*&
                           (S(i,j,k,5)+S(i,jp,k,5))*(mu(i,j,k)+mu(i,jp,k))*dyi-&
                  0.125_rp*(w_p(i,j,k)+w_p(i,jm,k)+w_p(i,j,km)+w_p(i,jm,km))*&    
                           (S(i,j,k,5)+S(i,jm,k,5))*(mu(i,j,k)+mu(i,jm,k))*dyi
          dwTzz = 0.500_rp*(w_p(i,j,k )*(S(i,j,k,6)+S(i,j,kp,6))*(mu(i,j,k)+mu(i,j,kp))-&
                            w_p(i,j,km)*(S(i,j,k,6)+S(i,j,km,6))*(mu(i,j,k)+mu(i,j,km)))*dzi
          !
          Tnu1_v = Tnu1_v + (duTxx+duTxy+duTxz+dvTyx+dvTyy+dvTyz+dwTzx+dwTzy+dwTzz)*psi(i,j,k)
          Tnu2_v = Tnu2_v + (duTxx+duTxy+duTxz+dvTyx+dvTyy+dvTyz+dwTzx+dwTzy+dwTzz)*(1._rp-psi(i,j,k))
          !
          !C9.lambda
          lam    =  (10.0*(mu(i,j,k)/rho(i,j,k))*ker/epsr)**0.5
 
          lambda = lambda + lam
         
          ! C10.Re_lambda
          
          Relam = Relam + (lam*rho(i,j,k)/mu(i,j,k))*(2.0/3.0*ker)**0.5
          
          ! C11. Kolmogorov length scale
          
          eta = eta + (mu(i,j,k)/rho(i,j,k))**(0.75) / (epsr**(0.25))   

          ! C12. area

          area = area + (((psi(i,j,k)-psi(im,j,k))/dx)**2 + & 
                 ((psi(i,j,k)-psi(i,jm,k))/dy)**2 + &
                 ((psi(i,j,k)-psi(i,j,km))/dz)**2)**0.5*dx*dy*dz

        enddo
      enddo
    enddo

    !
    call mpi_allreduce(MPI_IN_PLACE,E_E   ,1,MPI_REAL_RP,mpi_sum,comm_cart,ierr) ! a2
    call mpi_allreduce(MPI_IN_PLACE,ens1_v,1,MPI_REAL_RP,mpi_sum,comm_cart,ierr) ! a2
    call mpi_allreduce(MPI_IN_PLACE,ens2_v,1,MPI_REAL_RP,mpi_sum,comm_cart,ierr) ! a2
    E_E    = E_E   /vol
    ens1_v = ens1_v/vol1
    ens2_v = ens2_v/vol2
    !
    call mpi_allreduce(MPI_IN_PLACE,ur_t  ,1,MPI_REAL_RP,mpi_sum,comm_cart,ierr) ! b1
    call mpi_allreduce(MPI_IN_PLACE,vr_t  ,1,MPI_REAL_RP,mpi_sum,comm_cart,ierr) ! b1
    call mpi_allreduce(MPI_IN_PLACE,wr_t  ,1,MPI_REAL_RP,mpi_sum,comm_cart,ierr) ! b1
    call mpi_allreduce(MPI_IN_PLACE,ur1_v ,1,MPI_REAL_RP,mpi_sum,comm_cart,ierr) ! b1
    call mpi_allreduce(MPI_IN_PLACE,vr1_v ,1,MPI_REAL_RP,mpi_sum,comm_cart,ierr) ! b1
    call mpi_allreduce(MPI_IN_PLACE,wr1_v ,1,MPI_REAL_RP,mpi_sum,comm_cart,ierr) ! b1
    call mpi_allreduce(MPI_IN_PLACE,ur2_v ,1,MPI_REAL_RP,mpi_sum,comm_cart,ierr) ! b1
    call mpi_allreduce(MPI_IN_PLACE,vr2_v ,1,MPI_REAL_RP,mpi_sum,comm_cart,ierr) ! b1
    call mpi_allreduce(MPI_IN_PLACE,wr2_v ,1,MPI_REAL_RP,mpi_sum,comm_cart,ierr) ! b1
    ur_t  = sqrt(ur_t/vol)
    vr_t  = sqrt(vr_t/vol)
    wr_t  = sqrt(wr_t/vol)
    ur1_v = sqrt(ur1_v/vol1)
    vr1_v = sqrt(vr1_v/vol1)
    wr1_v = sqrt(wr1_v/vol1)
    ur2_v = sqrt(ur2_v/vol2)
    vr2_v = sqrt(vr2_v/vol2)
    wr2_v = sqrt(wr2_v/vol2)
    !
    call mpi_allreduce(MPI_IN_PLACE,um_t  ,1,MPI_REAL_RP,mpi_sum,comm_cart,ierr) ! b2
    call mpi_allreduce(MPI_IN_PLACE,vm_t  ,1,MPI_REAL_RP,mpi_sum,comm_cart,ierr) ! b2
    call mpi_allreduce(MPI_IN_PLACE,wm_t  ,1,MPI_REAL_RP,mpi_sum,comm_cart,ierr) ! b2
    call mpi_allreduce(MPI_IN_PLACE,um1_v ,1,MPI_REAL_RP,mpi_sum,comm_cart,ierr) ! b2
    call mpi_allreduce(MPI_IN_PLACE,vm1_v ,1,MPI_REAL_RP,mpi_sum,comm_cart,ierr) ! b2
    call mpi_allreduce(MPI_IN_PLACE,wm1_v ,1,MPI_REAL_RP,mpi_sum,comm_cart,ierr) ! b2
    call mpi_allreduce(MPI_IN_PLACE,um2_v ,1,MPI_REAL_RP,mpi_sum,comm_cart,ierr) ! b2
    call mpi_allreduce(MPI_IN_PLACE,vm2_v ,1,MPI_REAL_RP,mpi_sum,comm_cart,ierr) ! b2
    call mpi_allreduce(MPI_IN_PLACE,wm2_v ,1,MPI_REAL_RP,mpi_sum,comm_cart,ierr) ! b2
    um_t  = um_t/vol
    vm_t  = vm_t/vol
    wm_t  = wm_t/vol
    um1_v = um1_v/vol1
    vm1_v = vm1_v/vol1
    wm1_v = wm1_v/vol1
    um2_v = um2_v/vol2
    vm2_v = vm2_v/vol2
    wm2_v = wm2_v/vol2
    !
    call mpi_allreduce(MPI_IN_PLACE,eps_t ,1,MPI_REAL_RP,mpi_sum,comm_cart,ierr) ! c1
    call mpi_allreduce(MPI_IN_PLACE,eps1_v,1,MPI_REAL_RP,mpi_sum,comm_cart,ierr) ! c1
    call mpi_allreduce(MPI_IN_PLACE,eps2_v,1,MPI_REAL_RP,mpi_sum,comm_cart,ierr) ! c1
    eps_t  = eps_t /vol
    eps1_v = eps1_v/vol1
    eps2_v = eps2_v/vol2
    !
    call mpi_allreduce(MPI_IN_PLACE,prd_t ,1,MPI_REAL_RP,mpi_sum,comm_cart,ierr) ! c2
    call mpi_allreduce(MPI_IN_PLACE,prd1_v,1,MPI_REAL_RP,mpi_sum,comm_cart,ierr) ! c2
    call mpi_allreduce(MPI_IN_PLACE,prd2_v,1,MPI_REAL_RP,mpi_sum,comm_cart,ierr) ! c2
    prd_t  = prd_t /vol
    prd1_v = prd1_v/vol1
    prd2_v = prd2_v/vol2
    !
    call mpi_allreduce(MPI_IN_PLACE,surfp_t ,1,MPI_REAL_RP,mpi_sum,comm_cart,ierr) 
    call mpi_allreduce(MPI_IN_PLACE,surf1p_v,1,MPI_REAL_RP,mpi_sum,comm_cart,ierr) 
    call mpi_allreduce(MPI_IN_PLACE,surf2p_v,1,MPI_REAL_RP,mpi_sum,comm_cart,ierr) 
    call mpi_allreduce(MPI_IN_PLACE,surfi_t ,1,MPI_REAL_RP,mpi_sum,comm_cart,ierr) 
    call mpi_allreduce(MPI_IN_PLACE,surf1i_v,1,MPI_REAL_RP,mpi_sum,comm_cart,ierr) 
    call mpi_allreduce(MPI_IN_PLACE,surf2i_v,1,MPI_REAL_RP,mpi_sum,comm_cart,ierr) 
    surfp_t  = surfp_t /vol
    surf1p_v = surf1p_v/vol1
    surf2p_v = surf2p_v/vol2
    surfi_t  = surfi_t /vol
    surf1i_v = surf1i_v/vol1
    surf2i_v = surf2i_v/vol2
    !
    call mpi_allreduce(MPI_IN_PLACE,ke_t  ,1,MPI_REAL_RP,mpi_sum,comm_cart,ierr) ! c3
    call mpi_allreduce(MPI_IN_PLACE,ke1_v ,1,MPI_REAL_RP,mpi_sum,comm_cart,ierr) ! c3
    call mpi_allreduce(MPI_IN_PLACE,ke2_v ,1,MPI_REAL_RP,mpi_sum,comm_cart,ierr) ! c3
    ke_t  = ke_t /vol
    ke1_v = ke1_v/vol1
    ke2_v = ke2_v/vol2
    !
    call mpi_allreduce(MPI_IN_PLACE,Psi_mf  ,1,MPI_REAL_RP,mpi_sum,comm_cart,ierr) ! c4
    call mpi_allreduce(MPI_IN_PLACE,Psi_nu  ,1,MPI_REAL_RP,mpi_sum,comm_cart,ierr) ! c5
    call mpi_allreduce(MPI_IN_PLACE,Psi_gf_v,1,MPI_REAL_RP,mpi_sum,comm_cart,ierr) ! c6-a
    Psi_mf   = Psi_mf/vol
    Psi_nu   = Psi_nu/vol
    Psi_gf_v = Psi_gf_v/vol2
    !
    call mpi_allreduce(MPI_IN_PLACE,Tp1_v ,1,MPI_REAL_RP,mpi_sum,comm_cart,ierr) ! c7
    call mpi_allreduce(MPI_IN_PLACE,Tp2_v ,1,MPI_REAL_RP,mpi_sum,comm_cart,ierr) ! c7
    Tp1_v  = Tp1_v/vol1
    Tp2_v  = Tp2_v/vol2
    !
    call mpi_allreduce(MPI_IN_PLACE,Tnu1_v,1,MPI_REAL_RP,mpi_sum,comm_cart,ierr) ! c8
    call mpi_allreduce(MPI_IN_PLACE,Tnu2_v,1,MPI_REAL_RP,mpi_sum,comm_cart,ierr) ! c8
    Tnu1_v = Tnu1_v/vol1
    Tnu2_v = Tnu2_v/vol2
    !
    call mpi_allreduce(MPI_IN_PLACE,lambda,1,MPI_REAL_RP,mpi_sum,comm_cart,ierr) ! c8
    call mpi_allreduce(MPI_IN_PLACE,Relam,1,MPI_REAL_RP,mpi_sum,comm_cart,ierr) ! c8
    call mpi_allreduce(MPI_IN_PLACE,eta,1,MPI_REAL_RP,mpi_sum,comm_cart,ierr) ! c8
    call mpi_allreduce(MPI_IN_PLACE,area,1,MPI_REAL_RP,mpi_sum,comm_cart,ierr) ! c8
    lambda = lambda/vol
    Relam  = Relam/vol
    eta    = eta/vol

    if(myid.eq.0) then
      !
      ! overall budget
      !
      open(94,file=trim(datadir_bal)//'stats.out',position='append')
      write(94,'(19E19.7)') 1._rp*istep,time,lambda,Relam,area,eta,um_t,vm_t,wm_t,ur_t,vr_t,wr_t, &
                            ens_t,ke_t,prd_t,eps_t,psi_nu,psi_mf,psi_gf_v 
      close(94)
      !
      !
      !rms values    
!      open(94,file='data/post/balance/rms.out',position='append')
      open(94,file=trim(datadir_bal)//'rms.out',position='append')
      write(94,'(8E15.7)') 1._rp*istep,time,u_avg,v_avg,w_avg,ur_t,vr_t,wr_t
                      
      close(94)

      ! phase 1 - vol
      !
      open(94,file=trim(datadir_bal)//'budget_phase1_v.out',position='append')
      write(94,'(14E15.7)') 1._rp*istep,time,um1_v,vm1_v,wm1_v,ur1_v,vr1_v,wr1_v, &
                            ens1_v,ke1_v,prd1_v,eps1_v,Tp1_v,Tnu1_v 
      close(94)
      !
      ! phase 2 - vol
      !
      open(94,file=trim(datadir_bal)//'budget_phase2_v.out',position='append')
      write(94,'(14E15.7)') 1._rp*istep,time,um2_v,vm2_v,wm2_v,ur2_v,vr2_v,wr2_v, &
                            ens2_v,ke2_v,prd2_v,eps2_v,Tp2_v,Tnu2_v 
      close(94)
      !
      open(94,file=trim(datadir_bal)//'surf_work.out',position='append')
      write(94,'(19E19.7)') 1._rp*istep,time,surfi_t,surfp_t,surf1i_v,surf2i_v,surf1p_v,surf2p_v                    
      close(94)

    endif
    !
    return
  end subroutine energy_balance_mf
  !
  subroutine curv_pdf1(istep,phi,curv_field)
  implicit none
  integer,intent(in)                         :: istep
  real(rp),dimension(:,:,:),intent(in)       :: phi,curv_field
  real(rp),dimension(:),allocatable          :: curvature,curve_hist,curve_pdf
  integer                                    :: histogram,i,j,k,l
  real(rp)                                   :: hist_counter,coeff 
  real(rp)                                   :: curve_max,curve_min,delta_curve
  real(rp)                                   :: hist_begin,hist_end
  logical                                    :: interface_cond,cond_begin_hist,cond_end_hist
  character(len=8)                           :: fldnum
  !  
  histogram               = 500
  allocate(curvature(histogram))
  allocate(curve_hist(histogram))
  allocate(curve_pdf(histogram))
  !  
  curvature (1:histogram) = 0._rp
  curve_hist(1:histogram) = 0._rp
  curve_pdf(1:histogram)  = 0._rp
  !
  hist_counter            = 0._rp
  curve_max               =-1.e14
  curve_min               = 1.e14
  !
  do k=1,n(3)
     do j=1,n(2)
        do i=1,n(1)
           if (curv_field(i,j,k).gt.(-1e13)) then
              interface_cond =(phi(i,j,k).gt.(0.1)).and.(phi(i,j,k).lt.(0.9))
              if (interface_cond) then
                 curve_max = max(curve_max,curv_field(i,j,k))
                 curve_min = min(curve_min,curv_field(i,j,k))
              endif
           endif
        enddo
     enddo
  enddo
  !
  call mpi_allreduce(MPI_IN_PLACE,curve_max,1,mpi_real8,mpi_max,comm_cart,ierr)
  call mpi_allreduce(MPI_IN_PLACE,curve_min,1,mpi_real8,mpi_min,comm_cart,ierr)
  !
  delta_curve = (curve_max-curve_min)/(histogram)
  !
  do k=1,n(3)
     do j=1,n(2)
        do i=1,n(1)
           interface_cond =(phi(i,j,k).gt.(0.1)).and.(phi(i,j,k).lt.(0.9))
           if (interface_cond) then
              do l = 1,histogram
                 hist_begin                          = curve_min+ (l-1)*delta_curve
                 hist_end                            = curve_min+ l*delta_curve
                 cond_begin_hist                     = curv_field(i,j,k).ge.hist_begin
                 cond_end_hist                       = curv_field(i,j,k).lt.hist_end
                 if(l.eq.histogram) cond_end_hist    = curv_field(i,j,k).le.hist_end
                 if (cond_begin_hist.and.cond_end_hist) then
                    hist_counter  = hist_counter + 1._rp
                    curvature(l)  = curv_field(i,j,k)
                    curve_hist(l) = curve_hist(l)+1._rp
                    exit
                 endif
              enddo    
           endif
        enddo
     enddo 
  enddo
  !
  call mpi_allreduce(MPI_IN_PLACE,curvature,histogram,mpi_real8,mpi_sum,comm_cart,ierr)
  call mpi_allreduce(MPI_IN_PLACE,curve_hist,histogram,mpi_real8,mpi_sum,comm_cart,ierr)
  call mpi_allreduce(MPI_IN_PLACE,hist_counter,1,mpi_real8,mpi_sum,comm_cart,ierr)
  !
  do l=1,histogram
     curvature(l)=curvature(l)/(dims_in(1)*dims_in(2))
  enddo
  !
  coeff=1._rp/(delta_curve*sum(curve_hist))
  do l=1,histogram
     curve_pdf(l)=curve_hist(l)*coeff
  enddo
  !  
  if (myid.eq.0) then
     write(fldnum,'(i8.8)') istep
     open(1337,file=trim(datadir_cur1)//'cur_pdf'//fldnum//'.txt',position='append')
     do l= 1,histogram
        write(1337,'(6E16.8)' ) curve_min,curve_max,curvature(l),curve_hist(l),curve_pdf(l),hist_counter
     enddo
     close(1337)
  endif
  !
  deallocate(curvature)
  deallocate(curve_hist)
  deallocate(curve_pdf)
  return
  !
  end subroutine curv_pdf1
  !
  subroutine curv_pdf2(istep,phi,curv_field)
  implicit none
  integer,intent(in)                         :: istep
  real(rp),dimension(:,:,:),intent(in)       :: phi,curv_field
  real(rp),dimension(:),allocatable          :: curvature,curve_hist,curve_pdf
  integer                                    :: histogram,i,j,k,l
  real(rp)                                   :: hist_counter,coeff 
  real(rp)                                   :: curve_max,curve_min,delta_curve
  real(rp)                                   :: hist_begin,hist_end
  logical                                    :: interface_cond,cond_begin_hist,cond_end_hist
  character(len=8)                           :: fldnum
  !  
  histogram               = 500
  allocate(curvature(histogram))
  allocate(curve_hist(histogram))
  allocate(curve_pdf(histogram))
  !  
  curvature (1:histogram) = 0._rp
  curve_hist(1:histogram) = 0._rp
  curve_pdf(1:histogram)  = 0._rp
  !
  hist_counter            = 0._rp
  curve_max               =-1.e14
  curve_min               = 1.e14
  !
  do k=1,n(3)
     do j=1,n(2)
        do i=1,n(1)
           if (curv_field(i,j,k).gt.(-1e13)) then
              interface_cond =(phi(i,j,k).gt.(0.2)).and.(phi(i,j,k).lt.(0.8))
              if (interface_cond) then
                 curve_max = max(curve_max,curv_field(i,j,k))
                 curve_min = min(curve_min,curv_field(i,j,k))
              endif
           endif
        enddo
     enddo
  enddo
  !
  call mpi_allreduce(MPI_IN_PLACE,curve_max,1,mpi_real8,mpi_max,comm_cart,ierr)
  call mpi_allreduce(MPI_IN_PLACE,curve_min,1,mpi_real8,mpi_min,comm_cart,ierr)
  !
  delta_curve = (curve_max-curve_min)/(histogram)
  !
  do k=1,n(3)
     do j=1,n(2)
        do i=1,n(1)
           interface_cond =(phi(i,j,k).gt.(0.2)).and.(phi(i,j,k).lt.(0.8))
           if (interface_cond) then
              do l = 1,histogram
                 hist_begin                          = curve_min+ (l-1)*delta_curve
                 hist_end                            = curve_min+ l*delta_curve
                 cond_begin_hist                     = curv_field(i,j,k).ge.hist_begin
                 cond_end_hist                       = curv_field(i,j,k).lt.hist_end
                 if(l.eq.histogram) cond_end_hist    = curv_field(i,j,k).le.hist_end
                 if (cond_begin_hist.and.cond_end_hist) then
                    hist_counter  = hist_counter + 1._rp
                    curvature(l)  = curv_field(i,j,k)
                    curve_hist(l) = curve_hist(l)+1._rp
                    exit
                 endif
              enddo    
           endif
        enddo
     enddo 
  enddo
  !
  call mpi_allreduce(MPI_IN_PLACE,curvature,histogram,mpi_real8,mpi_sum,comm_cart,ierr)
  call mpi_allreduce(MPI_IN_PLACE,curve_hist,histogram,mpi_real8,mpi_sum,comm_cart,ierr)
  call mpi_allreduce(MPI_IN_PLACE,hist_counter,1,mpi_real8,mpi_sum,comm_cart,ierr)
  !
  do l=1,histogram
     curvature(l)=curvature(l)/(dims_in(1)*dims_in(2))
  enddo
  !
  coeff=1._rp/(delta_curve*sum(curve_hist))
  do l=1,histogram
     curve_pdf(l)=curve_hist(l)*coeff
  enddo
  !  
  if (myid.eq.0) then
     write(fldnum,'(i8.8)') istep
     open(1337,file=trim(datadir_cur2)//'cur_pdf'//fldnum//'.txt',position='append')
     do l= 1,histogram
        write(1337,'(6E16.8)' ) curve_min,curve_max,curvature(l),curve_hist(l),curve_pdf(l),hist_counter
     enddo
     close(1337)
  endif
  !
  deallocate(curvature)
  deallocate(curve_hist)
  deallocate(curve_pdf)
  return
  !
  end subroutine curv_pdf2
  !
  subroutine curv_pdf3(istep,phi,curv_field)
  implicit none
  integer,intent(in)                         :: istep
  real(rp),dimension(:,:,:),intent(in)       :: phi,curv_field
  real(rp),dimension(:),allocatable          :: curvature,curve_hist,curve_pdf
  integer                                    :: histogram,i,j,k,l
  real(rp)                                   :: hist_counter,coeff 
  real(rp)                                   :: curve_max,curve_min,delta_curve
  real(rp)                                   :: hist_begin,hist_end
  logical                                    :: interface_cond,cond_begin_hist,cond_end_hist
  character(len=8)                           :: fldnum
  !  
  histogram               = 500
  allocate(curvature(histogram))
  allocate(curve_hist(histogram))
  allocate(curve_pdf(histogram))
  !  
  curvature (1:histogram) = 0._rp
  curve_hist(1:histogram) = 0._rp
  curve_pdf(1:histogram)  = 0._rp
  !
  hist_counter            = 0._rp
  curve_max               =-1.e14
  curve_min               = 1.e14
  !
  do k=1,n(3)
     do j=1,n(2)
        do i=1,n(1)
           if (curv_field(i,j,k).gt.(-1e13)) then
              interface_cond =(phi(i,j,k).gt.(0.3)).and.(phi(i,j,k).lt.(0.7))
              if (interface_cond) then
                 curve_max = max(curve_max,curv_field(i,j,k))
                 curve_min = min(curve_min,curv_field(i,j,k))
              endif
           endif
        enddo
     enddo
  enddo
  !
  call mpi_allreduce(MPI_IN_PLACE,curve_max,1,mpi_real8,mpi_max,comm_cart,ierr)
  call mpi_allreduce(MPI_IN_PLACE,curve_min,1,mpi_real8,mpi_min,comm_cart,ierr)
  !
  delta_curve = (curve_max-curve_min)/(histogram)
  !
  do k=1,n(3)
     do j=1,n(2)
        do i=1,n(1)
           interface_cond =(phi(i,j,k).gt.(0.3)).and.(phi(i,j,k).lt.(0.7))
           if (interface_cond) then
              do l = 1,histogram
                 hist_begin                          = curve_min+ (l-1)*delta_curve
                 hist_end                            = curve_min+ l*delta_curve
                 cond_begin_hist                     = curv_field(i,j,k).ge.hist_begin
                 cond_end_hist                       = curv_field(i,j,k).lt.hist_end
                 if(l.eq.histogram) cond_end_hist    = curv_field(i,j,k).le.hist_end
                 if (cond_begin_hist.and.cond_end_hist) then
                    hist_counter  = hist_counter + 1._rp
                    curvature(l)  = curv_field(i,j,k)
                    curve_hist(l) = curve_hist(l)+1._rp
                    exit
                 endif
              enddo    
           endif
        enddo
     enddo 
  enddo
  !
  call mpi_allreduce(MPI_IN_PLACE,curvature,histogram,mpi_real8,mpi_sum,comm_cart,ierr)
  call mpi_allreduce(MPI_IN_PLACE,curve_hist,histogram,mpi_real8,mpi_sum,comm_cart,ierr)
  call mpi_allreduce(MPI_IN_PLACE,hist_counter,1,mpi_real8,mpi_sum,comm_cart,ierr)
  !
  do l=1,histogram
     curvature(l)=curvature(l)/(dims_in(1)*dims_in(2))
  enddo
  !
  coeff=1._rp/(delta_curve*sum(curve_hist))
  do l=1,histogram
     curve_pdf(l)=curve_hist(l)*coeff
  enddo
  !  
  if (myid.eq.0) then
     write(fldnum,'(i8.8)') istep
     open(1337,file=trim(datadir_cur3)//'cur_pdf'//fldnum//'.txt',position='append')
     do l= 1,histogram
        write(1337,'(6E16.8)' ) curve_min,curve_max,curvature(l),curve_hist(l),curve_pdf(l),hist_counter
     enddo
     close(1337)
  endif
  !
  deallocate(curvature)
  deallocate(curve_hist)
  deallocate(curve_pdf)
  return
  !
  end subroutine curv_pdf3
  !
  subroutine curv_pdf4(istep,phi,curv_field)
  implicit none
  integer,intent(in)                         :: istep
  real(rp),dimension(:,:,:),intent(in)       :: phi,curv_field
  real(rp),dimension(:),allocatable          :: curvature,curve_hist,curve_pdf
  integer                                    :: histogram,i,j,k,l
  real(rp)                                   :: hist_counter,coeff 
  real(rp)                                   :: curve_max,curve_min,delta_curve
  real(rp)                                   :: hist_begin,hist_end
  logical                                    :: interface_cond,cond_begin_hist,cond_end_hist
  character(len=8)                           :: fldnum
  !  
  histogram               = 500
  allocate(curvature(histogram))
  allocate(curve_hist(histogram))
  allocate(curve_pdf(histogram))
  !  
  curvature (1:histogram) = 0._rp
  curve_hist(1:histogram) = 0._rp
  curve_pdf(1:histogram)  = 0._rp
  !
  hist_counter            = 0._rp
  curve_max               =-1.e14
  curve_min               = 1.e14
  !
  do k=1,n(3)
     do j=1,n(2)
        do i=1,n(1)
           if (curv_field(i,j,k).gt.(-1e13)) then
              interface_cond =(phi(i,j,k).gt.(0.4)).and.(phi(i,j,k).lt.(0.6))
              if (interface_cond) then
                 curve_max = max(curve_max,curv_field(i,j,k))
                 curve_min = min(curve_min,curv_field(i,j,k))
              endif
           endif
        enddo
     enddo
  enddo
  !
  call mpi_allreduce(MPI_IN_PLACE,curve_max,1,mpi_real8,mpi_max,comm_cart,ierr)
  call mpi_allreduce(MPI_IN_PLACE,curve_min,1,mpi_real8,mpi_min,comm_cart,ierr)
  !
  delta_curve = (curve_max-curve_min)/(histogram)
  !
  do k=1,n(3)
     do j=1,n(2)
        do i=1,n(1)
           interface_cond =(phi(i,j,k).gt.(0.4)).and.(phi(i,j,k).lt.(0.6))
           if (interface_cond) then
              do l = 1,histogram
                 hist_begin                          = curve_min+ (l-1)*delta_curve
                 hist_end                            = curve_min+ l*delta_curve
                 cond_begin_hist                     = curv_field(i,j,k).ge.hist_begin
                 cond_end_hist                       = curv_field(i,j,k).lt.hist_end
                 if(l.eq.histogram) cond_end_hist    = curv_field(i,j,k).le.hist_end
                 if (cond_begin_hist.and.cond_end_hist) then
                    hist_counter  = hist_counter + 1._rp
                    curvature(l)  = curv_field(i,j,k)
                    curve_hist(l) = curve_hist(l)+1._rp
                    exit
                 endif
              enddo    
           endif
        enddo
     enddo 
  enddo
  !
  call mpi_allreduce(MPI_IN_PLACE,curvature,histogram,mpi_real8,mpi_sum,comm_cart,ierr)
  call mpi_allreduce(MPI_IN_PLACE,curve_hist,histogram,mpi_real8,mpi_sum,comm_cart,ierr)
  call mpi_allreduce(MPI_IN_PLACE,hist_counter,1,mpi_real8,mpi_sum,comm_cart,ierr)
  !
  do l=1,histogram
     curvature(l)=curvature(l)/(dims_in(1)*dims_in(2))
  enddo
  !
  coeff=1._rp/(delta_curve*sum(curve_hist))
  do l=1,histogram
     curve_pdf(l)=curve_hist(l)*coeff
  enddo
  !  
  if (myid.eq.0) then
     write(fldnum,'(i8.8)') istep
     open(1337,file=trim(datadir_cur4)//'cur_pdf'//fldnum//'.txt',position='append')
     do l= 1,histogram
        write(1337,'(6E16.8)' ) curve_min,curve_max,curvature(l),curve_hist(l),curve_pdf(l),hist_counter
     enddo
     close(1337)
  endif
  !
  deallocate(curvature)
  deallocate(curve_hist)
  deallocate(curve_pdf)
  return
  !
  end subroutine curv_pdf4
  !
end module mod_post
