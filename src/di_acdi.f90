! SPDX-License-Identifier: MIT
!
module mod_di_acdi
  !
  ! based on:
  ! --> Mirjalili, S. et. al. "A conservative diffuse interface method for two-phase flows
  !     with provable boundedness properties." 
  !     Journal of Computational Physics 401 (2020): 109006.
  ! --> Jain, S. S. "Accurate conservative phase-field method 
  !     for simulation of two-phase flows." 
  !     Journal of Computational Physics 469 (2022): 111529.
  !
  ! --> Attraction model developed by Salar Zamani & Armin Shahmardi 
  !
  use mod_types
  use mod_common_mpi, only: myid,ierr,ijk_start,comm_cart
  use mod_sanity    , only: flutas_error
  !
  implicit none
  !
  private
  public :: update_property,rk_psi,psi_to_ls,ls_to_psi,cmpt_norm,initls,initpsi,cmpt_umax
  public rk_c1,initc1, interfacial_cflux
  !
  public static_contact_angle_m
  !
  contains
    !
    !
  subroutine static_contact_angle_m(n,dl,nh_v,eps,theta,ap,ap_bound)
    !                             
    ! To do: we can treat the term multiplying grad_rho implicitly
    ! 
    implicit none
    !
    integer , intent(in ), dimension(3) :: n
    real(rp), intent(in ), dimension(3) :: dl
    integer , intent(in )               :: nh_v
    real(rp), intent(in )               :: eps,theta
    real(rp), intent(in ), dimension(1-nh_v:,1-nh_v:,1-nh_v:) :: ap
    real(rp), intent(out), dimension( 1:, 1:)     :: ap_bound
    real(rp) :: f,mv_wall,ap_wall
    integer  :: i,j
    !
    do j=1,n(2)
      do i=1,n(1)
        ap_wall  = 0.5*(ap(i,j,0)+ap(i,j,1))
        f        = ap_wall*(1.d0-ap_wall) ! non-linear term treated explicitly
        ap_bound(i,j) = ap(i,j,1)+(dl(3)/eps)*cos(theta)*f
        !
      enddo
    enddo
    return
  end subroutine static_contact_angle_m
    !
    subroutine rk_c1(f_t1,f_t2,n,dli,nh_d,nh_u,nh_v, &
                      dzc,dzf,diff1,diff2,u_reg,eps_int,u,v,w,c1,c2,phi,nor,dc1dtrko,dc2dtrko,id)     
                                                                          !phi is order paramter
                                                                          !c1 is scalar
    !
    ! low-storage 3rd-order Runge-Kutta scheme 
    ! for time integration of the c1 field.
    !
    implicit none
    !
    real(rp), intent(in   )                                        :: f_t1,f_t2
    integer , intent(in   ), dimension(3)                          :: n
    real(rp), intent(in   ), dimension(3)                          :: dli
    integer , intent(in   )                                        :: nh_d,nh_u,nh_v
    real(rp), intent(in   ), dimension(1-nh_d:)                    :: dzc,dzf
    real(rp), intent(in   )                                        :: diff1,diff2,eps_int,u_reg
    real(rp), intent(in   ), dimension(1-nh_u:,1-nh_u:,1-nh_u:)    :: u,v,w
    real(rp), intent(inout), dimension(1-nh_v:,1-nh_v:,1-nh_v:)    :: c1,c2  ! concentraion at current time level
    real(rp), intent(in   ), dimension(1-nh_v:,1-nh_v:,1-nh_v:)    :: phi ! level-set,sai
    real(rp), intent(in   ), dimension(1-nh_v:,1-nh_v:,1-nh_v:,1:) :: nor
    real(rp), intent(inout), dimension(     1:,     1:,     1:)    :: dc1dtrko,dc2dtrko
    integer , intent(in   )                                        :: id
    ! 
    real(rp), dimension(n(1),n(2),n(3)) :: dc1dtrk,dc2dtrk
    real(rp) :: factor1,factor2,factor3,tot_reg,int_den,q_t,tot_lag
    integer  :: i,j,k
    !
    factor1 = f_t1!rkpar(1)*dt
    factor2 = f_t2!rkpar(2)*dt
    factor3 = f_t1+f_t2!factor1+factor2
    !
    dc1dtrk = 0.d0
    dc2dtrk = 0.d0
    !
    ! 1. compute regularization
    !
    !
    if (id.eq.1)then
    call regularization_c1(n,dli,nh_d,nh_v,dzc,dzf,diff1,u_reg,eps_int,c1,Phi,nor,dc1dtrk)
    elseif(id.eq.-1)then 
    call regularization_c1(n,dli,nh_d,nh_v,dzc,dzf,diff2,u_reg,eps_int,c2,1.d0-Phi,-1.d0*nor,dc2dtrk)
    endif 
    !
    ! 2. compute interfacial scalar transfer
    !
    if (id.eq.1)then 
    call interfacial_cflux(n,dli,nh_d,nh_v,dzc,dzf,diff1,diff2,c1,c2,phi,dc1dtrk,id)
    elseif(id.eq.-1)then 
    call interfacial_cflux(n,dli,nh_d,nh_v,dzc,dzf,diff1,diff2,c1,c2,phi,dc2dtrk,id)
    endif
    !
    ! 3. Add advection 
    ! 
    if (id.eq.1)then 
    call moms_psi(n,dli,nh_d,nh_u,nh_v,u,v,w,c1,dc1dtrk)
    elseif(id.eq.-1)then 
    call moms_psi(n,dli,nh_d,nh_u,nh_v,u,v,w,c2,dc2dtrk)
    elseif(id.eq.-1)then 
    endif
    !
    !
    !
#if defined(_CONTACT_LINE)
   tot_reg = 0._rp
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          tot_reg = tot_reg + dc1dtrk(i,j,k)/(dli(1)*dli(2)*dli(3))
        enddo
      enddo
    enddo
    call mpi_allreduce(MPI_IN_PLACE,tot_reg,1,MPI_REAL_RP,MPI_SUM,comm_cart,ierr)
    if(myid.eq.0) print*, "Tot. reg.", tot_reg
    !
    int_den = 0._rp
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          int_den = int_den + 6._rp*c1(i,j,k)*(1._rp-c1(i,j,k))/(dli(1)*dli(2)*dli(3))
        enddo
      enddo
    enddo
    call mpi_allreduce(MPI_IN_PLACE,int_den,1,MPI_REAL_RP,MPI_SUM,comm_cart,ierr)
    q_t = tot_reg/int_den
#endif
    !
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP SHARED(n,f_t1,f_t2,psi,dpsidtrk,dpsidtrko)
    !
    !
    !
    !
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          !
#if defined(_CONTACT_LINE)
          tot_lag          = 6._rp*c1(i,j,k)*(1._rp-c1(i,j,k))*q_t
          c1(i,j,k)        = c1(i,j,k) + f_t1*( dc1dtrk(i,j,k)-tot_lag) + f_t2*dc1dtrko(i,j,k)
          dc1dtrko(i,j,k)  = dc1dtrk(i,j,k)-tot_lag
#else
   if(id.eq.1)then
          c1(i,j,k)        = c1(i,j,k) + f_t1*dc1dtrk(i,j,k) + f_t2*dc1dtrko(i,j,k)
          dc1dtrko(i,j,k)  = dc1dtrk(i,j,k)
   elseif(id.eq.-1) then 
          c2(i,j,k)        = c2(i,j,k) + f_t1*dc2dtrk(i,j,k) + f_t2*dc2dtrko(i,j,k)
          dc2dtrko(i,j,k)  = dc2dtrk(i,j,k)
   endif
#endif
          !
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
    return
  !
  end subroutine rk_c1
  !
  subroutine rk_psi(f_t1,f_t2,n,dli,nh_d,nh_u,nh_v, &
                    dzc,dzf,u_reg,eps_int,u,v,w,psi,ls,nor,kappa,dpsidtrko)
    !
    ! low-storage 3rd-order Runge-Kutta scheme 
    ! for time integration of the psi field.
    !
    implicit none
    !
    real(rp), intent(in   )                                        :: f_t1,f_t2
    integer , intent(in   ), dimension(3)                          :: n
    real(rp), intent(in   ), dimension(3)                          :: dli
    integer , intent(in   )                                        :: nh_d,nh_u,nh_v
    real(rp), intent(in   ), dimension(1-nh_d:)                    :: dzc,dzf
    real(rp), intent(in   )                                        :: u_reg,eps_int
    real(rp), intent(in   ), dimension(1-nh_u:,1-nh_u:,1-nh_u:)    :: u,v,w
    real(rp), intent(inout), dimension(1-nh_v:,1-nh_v:,1-nh_v:)    :: psi  ! volume-fraction at current time level
    real(rp), intent(in   ), dimension(1-nh_v:,1-nh_v:,1-nh_v:)    :: ls   ! level-set
    real(rp), intent(in   ), dimension(1-nh_v:,1-nh_v:,1-nh_v:,1:) :: nor
    real(rp), intent(in ), dimension( 0:, 0:, 0:)               :: kappa
    real(rp), intent(inout), dimension(     1:,     1:,     1:)    :: dpsidtrko
    ! 
    real(rp), dimension(n(1),n(2),n(3)) :: dpsidtrk
    real(rp) :: tot_reg,phase_change
#if defined(_CONTACT_LINE)
    real(rp) :: int_den,q_t,tot_lag
#endif
    integer  :: i,j,k
    !
    ! 1. compute regularization
    !
#if defined(_USE_ACDI)    
    call regularization_psi(n,dli,nh_d,nh_v,dzc,dzf,u_reg,eps_int,psi,ls,nor,dpsidtrk)
#else    
    call regularization_psi_cdi(n,dli,nh_d,nh_v,dzc,dzf,u_reg,eps_int,psi,nor,dpsidtrk)
#endif
    !
    ! 2. Compute lag mult
    !
#if defined(_CONTACT_LINE)
    tot_reg = 0._rp
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          tot_reg = tot_reg + dpsidtrk(i,j,k)/(dli(1)*dli(2)*dli(3))
        enddo
      enddo
    enddo
    call mpi_allreduce(MPI_IN_PLACE,tot_reg,1,MPI_REAL_RP,MPI_SUM,comm_cart,ierr)
    if(myid.eq.0) print*, "Tot. reg.", tot_reg
    !
    int_den = 0._rp
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          int_den = int_den + 6._rp*psi(i,j,k)*(1._rp-psi(i,j,k))/(dli(1)*dli(2)*dli(3))
        enddo
      enddo
    enddo
    call mpi_allreduce(MPI_IN_PLACE,int_den,1,MPI_REAL_RP,MPI_SUM,comm_cart,ierr)
    q_t = tot_reg/int_den
#endif
    !
    ! 3. Add advection 
    ! 
    call moms_psi(n,dli,nh_d,nh_u,nh_v,u,v,w,psi,dpsidtrk)
    !
    ! 4. Advance the color function
    !
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP SHARED(n,f_t1,f_t2,psi,dpsidtrk,dpsidtrko)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          !
#if defined(_CONTACT_LINE)
          tot_lag          = 6._rp*psi(i,j,k)*(1._rp-psi(i,j,k))*q_t
          psi(i,j,k)       = psi(i,j,k) + f_t1*( dpsidtrk(i,j,k)-tot_lag ) + f_t2*dpsidtrko(i,j,k)
          dpsidtrko(i,j,k) = dpsidtrk(i,j,k)-tot_lag
#else

          psi(i,j,k)       = psi(i,j,k) + f_t1*(dpsidtrk(i,j,k)) + f_t2*dpsidtrko(i,j,k)
          dpsidtrko(i,j,k) = dpsidtrk(i,j,k)
#endif
          !
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
    !
    return
  end subroutine rk_psi
  !
  subroutine moms_psi(n,dli,nh_d,nh_u,nh_v,u,v,w,s,dsdt)
    !
    implicit none
    !
    integer , intent(in   ), dimension(3)                       :: n
    real(rp), intent(in   ), dimension(3)                       :: dli
    integer , intent(in   )                                     :: nh_d,nh_u,nh_v
    real(rp), intent(in   ), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w
    real(rp), intent(in   ), dimension(1-nh_v:,1-nh_v:,1-nh_v:) :: s
    real(rp), intent(inout), dimension(     1:,     1:,     1:) :: dsdt
    !
    real(rp) :: s_flxp,s_flxm,s_flyp,s_flym,s_flzp,s_flzm
    integer  :: i,j,k
    !
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          !
          ! along x
          !
          s_flxp = 0.5_rp*( s(i+1,j,k)+s(i,j,k) )*u(i  ,j,k)
          s_flxm = 0.5_rp*( s(i-1,j,k)+s(i,j,k) )*u(i-1,j,k)
          !
          ! along y
          !
          s_flyp = 0.5_rp*( s(i,j+1,k)+s(i,j,k) )*v(i,j  ,k)
          s_flym = 0.5_rp*( s(i,j-1,k)+s(i,j,k) )*v(i,j-1,k)
          !
          ! along z
          !
          s_flzp = 0.5_rp*( s(i,j,k+1)+s(i,j,k) )*w(i,j,k  )
          s_flzm = 0.5_rp*( s(i,j,k-1)+s(i,j,k) )*w(i,j,k-1)
          !
          ! overall balance
          !
          dsdt(i,j,k) = - ( dli(1)*(s_flxp-s_flxm) + &
                            dli(2)*(s_flyp-s_flym) + &
                            dli(3)*(s_flzp-s_flzm) ) + dsdt(i,j,k)
          !
        enddo
      enddo
    enddo
    !
    return
  end subroutine moms_psi
  !
  !
  subroutine cmpt_umax(n,nh_u,u,v,w,gamma_v,u_reg)
    !
    implicit none
    !
    integer , intent(in ), dimension(3)                       :: n
    integer , intent(in )                                     :: nh_u
    real(rp), intent(in ), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w
    real(rp), intent(in )                                     :: gamma_v
    real(rp), intent(out)                                     :: u_reg
    !
    real(rp) :: abs_u_ijk 
    integer  :: i,j,k
    !
    u_reg = 0._rp 
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          !
          abs_u_ijk = sqrt(u(i,j,k)**2+v(i,j,k)**2+w(i,j,k)**2) 
          u_reg     = max(u_reg,abs_u_ijk)
          !
        enddo
      enddo
    enddo
    call mpi_allreduce(MPI_IN_PLACE,u_reg,1,MPI_REAL_RP,MPI_MAX,comm_cart,ierr)
    u_reg = gamma_v*u_reg
    !
    return
  end subroutine cmpt_umax
    !  
    !
    subroutine regularization_c1(n,dli,nh_d,nh_v, &
                               dzc,dzf,diff,gm,eps,s,phi,nor,dsdt) !s corresponds to c as concent.
                                                                   !phi corresponds to volume_fraction/phi
    !
    implicit none
    !
    integer , intent(in ), dimension(3)               :: n
    real(rp), intent(in ), dimension(3)               :: dli
    integer , intent(in )                             :: nh_d,nh_v
    real(rp), intent(in ), dimension(1-nh_d:)         :: dzc,dzf
    real(rp), intent(in )                             :: diff,gm,eps
    real(rp), intent(in ), dimension(1-nh_v:,1-nh_v:,1-nh_v:)    :: s
    real(rp), intent(in ), dimension(1-nh_v:,1-nh_v:,1-nh_v:)    :: Phi
    real(rp), intent(in ), dimension(1-nh_v:,1-nh_v:,1-nh_v:,1:) :: nor
    real(rp), intent(out), dimension(1:,1:,1:)                   :: dsdt
    real(rp),              dimension(0:n(1)+1,0:n(2)+1,0:n(3)+1) :: r
    !
    real(rp), dimension(3) :: dl
    real(rp) :: s_face,term_p,term_m,nor_face,s_min,s_max,rho_p,rho_m
    real(rp) :: Phi_face,A,D,urmax
    real(rp) :: gradrxp,gradrxm,gradryp,gradrym,gradrzp,gradrzm
    real(rp) :: ratioxp,ratioxm,ratioyp,ratioym,ratiozp,ratiozm
    real(rp) :: gradphixp,gradphixm,gradphiyp,gradphiym,gradphizp,gradphizm
    real(rp) :: phixp,phixm,phiyp,phiym,phizp,phizm
    real(rp) :: regxp,regxm,regyp,regym,regzp,regzm
    real(rp) :: diffxp,diffxm,diffyp,diffym,diffzp,diffzm
    real(rp) :: rhsxp,rhsxm,rhsyp,rhsym,rhszp,rhszm
    real(rp), parameter :: small = real(2.e-16,rp)
    integer  :: i,j,k
    integer  :: ip,im,jp,jm,kp,km
    !
    dl(:) = 1.d0/dli(:)
    dsdt(:,:,:) = 0.d0
    !
    !
    D = diff
    !
    !
 do k=0,n(3)+1
  do j=0,n(2)+1
   do i=0,n(1)+1
   r(i,j,k) = s(i,j,k)/(phi(i,j,k)+small)
   enddo
  enddo
 enddo

    !
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          !
          ip = i+1
          im = i-1
          jp = j+1
          jm = j-1
          kp = k+1
          km = k-1
          !
          ! along x-dir
          !
          !
#if defined(_TWOD)
          rhsxp = 0.d0
          rhsxm = 0.d0
#else
          gradrxp = ( r(ip,j,k) - r(i,j,k)   )*dli(1) !gradient of r=c/phi is computed at face i+1/2
          gradrxm = ( r(i,j,k)  - r(im,j,k)  )*dli(1) !gradient of r=c/phi is computed at face i-1/2

          ratioxp = ( r(i,j,k)   + r(ip,j,k) )*0.5    !phi/c is computed at face i+1/2
          ratioxm = ( r(i,j,k)   + r(im,j,k) )*0.5    !phi/c is computed at face i-1/2


          phixp   = ( phi(ip,j,k) + phi(i,j,k) )*0.5  !phi computed at i+1/2
          phixm   = ( phi(im,j,k) + phi(i,j,k) )*0.5  !phi computed at i-1/2

          gradphixp = ( phi(ip,j,k) - phi(i,j,k)   )*dli(1) !gradient phi computed at i+1/2
          gradphixm = ( phi(i,j,k)  - phi(im,j,k)  )*dli(1) !gradient phi computed at i-1/2

          regxp     = gm*(eps*gradphixp - ((phi(ip,j,k)*(1.d0-phi(ip,j,k))*nor(ip,j,k,1))  + &
                                           (phi(i,j,k)*(1.d0-phi(i,j,k))*nor(i,j,k,1)))*0.5)*ratioxp         !reg*ratio term computed at i+1/2

          regxm     = gm*(eps*gradphixm - ((phi(im,j,k)*(1.d0-phi(im,j,k))*nor(im,j,k,1))  + &
                                           (phi(i,j,k)*(1.d0-phi(i,j,k))*nor(i,j,k,1)))*0.5)*ratioxm         !reg*ratio term computed at i-1/2

          diffxp    = D*phixp*gradrxp                                                                        !diff term computed at i+1/2
          diffxm    = D*phixm*gradrxm                                                                        !diff term computed at i-1/2
          !
          rhsxp     = regxp + diffxp
          rhsxm     = regxm + diffxm
#endif
          !
          !
          !
          ! along y-dir
          !
          !
           gradryp = ( r(i,jp,k) - r(i,j,k)   )*dli(2) !gradient of r=c/phi is computed at face j+1/2
           gradrym = ( r(i,j,k)  - r(i,jm,k)  )*dli(2) !gradient of r=c/phi is computed at face j-1/2

           ratioyp = ( r(i,j,k)   + r(i,jp,k) )*0.5    !phi/c is computed at face j+1/2
           ratioym = ( r(i,j,k)   + r(i,jm,k) )*0.5    !phi/c is computed at face j-1/2

           phiyp   = ( phi(i,jp,k) + phi(i,j,k) )*0.5  !phi computed at j+1/2
           phiym   = ( phi(i,jm,k) + phi(i,j,k) )*0.5  !phi computed at j-1/2

           gradphiyp = ( phi(i,jp,k) - phi(i,j,k)   )*dli(2)  !gradient phi computed at j+1/2
           gradphiym = ( phi(i,j,k)  - phi(i,jm,k)  )*dli(2)  !gradient phi computed at j-1/2

           regyp     = gm*(eps*gradphiyp - ((phi(i,jp,k)*(1.d0-phi(i,jp,k))*nor(i,jp,k,2))  + &
                                            (phi(i,j,k)*(1.d0-phi(i,j,k))*nor(i,j,k,2)))*0.5)*ratioyp         !reg*ratio term computed at j+1/2

           regym     = gm*(eps*gradphiym - ((phi(i,jm,k)*(1.d0-phi(i,jm,k))*nor(i,jm,k,2))  + &
                                            (phi(i,j,k)*(1.d0-phi(i,j,k))*nor(i,j,k,2)))*0.5)*ratioym         !reg*ratio term computed at j-1/2

           diffyp    = D*phiyp*gradryp                                                                        !diff term computed at j+1/2
           diffym    = D*phiym*gradrym                                                                        !diff term computed at j-1/2

           rhsyp     = regyp + diffyp
           rhsym     = regym + diffym


          ! along z-dir
          !
          !
          gradrzp = ( r(i,j,kp) - r(i,j,k)   )*dli(3) !gradient of r=c/phi is computed at face k+1/2
          gradrzm = ( r(i,j,k)  - r(i,j,km)  )*dli(3) !gradient of r=c/phi is computed at face k-1/2

          ratiozp = ( r(i,j,k)   + r(i,j,kp) )*0.5    !phi/c is computed at face k+1/2
          ratiozm = ( r(i,j,k)   + r(i,j,km) )*0.5    !phi/c is computed at face k-1/2

          phizp   = ( phi(i,j,kp) + phi(i,j,k) )*0.5  !phi computed at k+1/2
          phizm   = ( phi(i,j,km) + phi(i,j,k) )*0.5  !phi computed at k-1/2

          gradphizp = ( phi(i,j,kp) - phi(i,j,k)   )*dli(3)  !gradient phi computed at k+1/2
          gradphizm = ( phi(i,j,k)  - phi(i,j,km)  )*dli(3)  !gradient phi computed at k-1/2

          regzp     = gm*(eps*gradphizp - ((phi(i,j,kp)*(1.d0-phi(i,j,kp))*nor(i,j,kp,3))  + &
                                           (phi(i,j,k)*(1.d0-phi(i,j,k))*nor(i,j,k,3)))*0.5)*ratiozp         !reg*ratio term computed at k+1/2

          regzm     = gm*(eps*gradphizm - ((phi(i,j,km)*(1.d0-phi(i,j,km))*nor(i,j,km,3))  + &
                                           (phi(i,j,k)*(1.d0-phi(i,j,k))*nor(i,j,k,3)))*0.5)*ratiozm         !reg*ratio term computed at k-1/2

          diffzp    = D*phizp*gradrzp                                                           !diff term computed at k+1/2
          diffzm    = D*phizm*gradrzm                                                           !diff term computed at k-1/2

          rhszp     = regzp + diffzp
          rhszm     = regzm + diffzm
          !
          dsdt(i,j,k)    =    (rhsxp - rhsxm) * dli(1) + &
                              (rhsyp - rhsym) * dli(2) + &
                              (rhszp - rhszm) * dli(3)
          !
        enddo
      enddo
    enddo
    !
    return
  end subroutine regularization_c1
  !
  !
    subroutine interfacial_cflux(n,dli,nh_d,nh_v, &
                                 dzc,dzf,diff1,diff2,c1,c2,phi,dsdt,id) 
    !
    use mod_param, only: A_m,k_eq 
    !
    implicit none
    !
    integer , intent(in ), dimension(3)               :: n
    real(rp), intent(in ), dimension(3)               :: dli 
    integer , intent(in )                             :: nh_d,nh_v
    real(rp), intent(in ), dimension(1-nh_d:)         :: dzc,dzf
    real(rp), intent(in )                             :: diff1,diff2
    real(rp), intent(in ), dimension(1-nh_v:,1-nh_v:,1-nh_v:)    :: c1,c2
    real(rp), intent(in ), dimension(1-nh_v:,1-nh_v:,1-nh_v:)    :: Phi
    real(rp), intent(inout), dimension(1:,1:,1:)                 :: dsdt
    integer , intent(in)                                         :: id
    !
    real(rp), dimension(3) :: dl
    real(rp) :: term1,term2_x,term2_y,term2_z,D_m
    integer  :: i,j,k
    !
    dl(:) = 1.d0/dli(:)
    term1 = 0.d0
    term2_x = 0.d0
    term2_y = 0.d0
    term2_z = 0.d0
    !
    !
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          !
          ! along x-dir
          !
          D_m = diff1*diff2/(k_eq*diff1*(1.d0-phi(i,j,k))+diff2*phi(i,j,k)+1e-30)
          !
          !First-term 
          term1 = A_m*D_m*(k_eq*c2(i,j,k)*phi(i,j,k)-c1(i,j,k)*(1.d0-phi(i,j,k)))
          !
          !Second-term along-x
          
          term2_x = -D_m*(phi(i+1,j,k)-phi(i-1,j,k))*0.5*dli(1)*(c1(i+1,j,k)-c1(i-1,j,k))*0.5*dli(1) - &
                     D_m*k_eq*(phi(i+1,j,k)-phi(i-1,j,k))*0.5*dli(1)*(c2(i+1,j,k)-c2(i-1,j,k))*0.5*dli(1) 
          !
          !Second-term along-y
          !
          term2_y = -D_m*(phi(i,j+1,k)-phi(i,j-1,k))*0.5*dli(2)*(c1(i,j+1,k)-c1(i,j-1,k))*0.5*dli(2) - &
                     D_m*k_eq*(phi(i,j+1,k)-phi(i,j-1,k))*0.5*dli(2)*(c2(i,j+1,k)-c2(i,j-1,k))*0.5*dli(2)
          !
          !      
          !Second-term along-z
          !
          term2_z = -D_m*(phi(i,j,k+1)-phi(i,j,k-1))*0.5*dli(3)*(c1(i,j,k+1)-c1(i,j,k-1))*0.5*dli(3) - &
                     D_m*k_eq*(phi(i,j,k+1)-phi(i,j,k-1))*0.5*dli(3)*(c2(i,j,k+1)-c2(i,j,k-1))*0.5*dli(3)
          !
          dsdt(i,j,k) = dsdt(i,j,k)+(term1+term2_x+term2_y+term2_z)*id
          !
        enddo
      enddo
    enddo
    !
    return
  end subroutine interfacial_cflux
  !
    subroutine regularization_psi_cdi(n,dli,nh_d,nh_v, &
                                     dzc,dzf,u_reg,eps_int,s,nor,dsdt)  !s is volume fraction/phi
    !
    implicit none
    !
    integer , intent(in ), dimension(3)               :: n
    real(rp), intent(in ), dimension(3)               :: dli 
    integer , intent(in )                             :: nh_d,nh_v
    real(rp), intent(in ), dimension(1-nh_d:)         :: dzc,dzf
    real(rp), intent(in )                             :: u_reg,eps_int
    real(rp), intent(in ), dimension(1-nh_v:,1-nh_v:,1-nh_v:)    :: s
    real(rp), intent(in ), dimension(1-nh_v:,1-nh_v:,1-nh_v:,1:) :: nor
    real(rp), intent(out), dimension(1:,1:,1:)                   :: dsdt
    !
    real(rp), dimension(3) :: dl
    real(rp) :: s_face,term_p,term_m,nor_face,s_min,s_max,rho_p,rho_m 
    real(rp) :: gradpxp,gradpxm,gradpyp,gradpym,gradpzp,gradpzm
    real(rp) :: psifxp,psifxm,psifyp,psifym,psifzp,psifzm
    real(rp) :: regxp,regxm,regyp,regym,regzp,regzm
    real(rp) :: abs_u_ijk
    real(rp) :: Psi_face
    integer  :: i,j,k,ip,im,jp,jm,kp,km
    !
    dsdt(:,:,:) = 0.d0
    !
    !
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
    !
        ip = i + 1
        im = i - 1
        jp = j + 1
        jm = j - 1
        kp = k + 1
        km = k - 1
#if defined(_TWOD)
      regxp   = 0.d0
      regxm   = 0.d0
#else
      gradpxp        =   ( s(ip,j,k) - s(i,j,k)  ) * dli(1)       !grad(phi) at i+0.5
      gradpxm        =   ( s(i,j,k)  - s(im,j,k) ) * dli(1)       !grad(phi) at i-0.5

      psifxp         =   ( s(i,j,k)*(1.d0-s(i,j,k))*nor(i,j,k,1) + s(ip,j,k)*(1.d0-s(ip,j,k))*nor(ip,j,k,1))*0.5 !phi*(1-phi)*nor at i+0.5 as p.15 Shahab.M
      psifxm         =   ( s(i,j,k)*(1.d0-s(i,j,k))*nor(i,j,k,1) + s(im,j,k)*(1.d0-s(im,j,k))*nor(im,j,k,1))*0.5 !phi*(1-phi)*nor at i-0.5 as p.15 Shahab.M

      regxp          =   u_reg*( eps_int*gradpxp - psifxp)        !regularization at i+0.5
      regxm          =   u_reg*( eps_int*gradpxm - psifxm)        !regularization at i-0.5
#endif

      gradpyp        =   ( s(i,jp,k) - s(i,j,k)  ) * dli(2)       !grad(phi) at j+0.5
      gradpym        =   ( s(i,j,k)  - s(i,jm,k) ) * dli(2)       !grad(phi) at j-0.5

      psifyp         =   ( s(i,j,k)*(1.d0-s(i,j,k))*nor(i,j,k,2) + s(i,jp,k)*(1.d0-s(i,jp,k))*nor(i,jp,k,2))*0.5 !phi*(1-phi)*nor at j+0.5 as p.15 Shahab.M
      psifym         =   ( s(i,j,k)*(1.d0-s(i,j,k))*nor(i,j,k,2) + s(i,jm,k)*(1.d0-s(i,jm,k))*nor(i,jm,k,2))*0.5 !phi*(1-phi)*nor at j-0.5 as p.15 Shahab.M

      regyp          =   u_reg*( eps_int*gradpyp - psifyp)        !regularization at j+0.5
      regym          =   u_reg*( eps_int*gradpym - psifym)        !regularization at j-0.5


      gradpzp        =   ( s(i,j,kp) - s(i,j,k)  ) * dli(3)       !grad(phi) at k+0.5
      gradpzm        =   ( s(i,j,k)  - s(i,j,km) ) * dli(3)       !grad(phi) at k-0.5

      psifzp         =   ( s(i,j,k)*(1.d0-s(i,j,k))*nor(i,j,k,3) + s(i,j,kp)*(1.d0-s(i,j,kp))*nor(i,j,kp,3))*0.5 !phi*(1-phi)*nor at k+0.5 as p.15 Shahab.M
      psifzm         =   ( s(i,j,k)*(1.d0-s(i,j,k))*nor(i,j,k,3) + s(i,j,km)*(1.d0-s(i,j,km))*nor(i,j,km,3))*0.5 !phi*(1-phi)*nor at k-0.5 as p.15 Shahab.M

      regzp          =   u_reg*( eps_int*gradpzp - psifzp)        !regularization at k+0.5
      regzm          =   u_reg*( eps_int*gradpzm - psifzm)        !regularization at k-0.5

 
      dsdt(i,j,k)    =    (regxp - regxm) * dli(1) + &
                          (regyp - regym) * dli(2) + &
                          (regzp - regzm) * dli(3) 


         !
        enddo
      enddo
    enddo
    !
    return
  end subroutine regularization_psi_cdi
  !
  !
  subroutine regularization_psi(n,dli,nh_d,nh_v, &
                                dzc,dzf,u_reg,eps_int,s,ls,nor,dsdt)
    !
    implicit none
    !
    integer , intent(in ), dimension(3)                          :: n
    real(rp), intent(in ), dimension(3)                          :: dli 
    integer , intent(in )                                        :: nh_d,nh_v
    real(rp), intent(in ), dimension(1-nh_d:)                    :: dzc,dzf
    real(rp), intent(in )                                        :: u_reg,eps_int
    real(rp), intent(in ), dimension(1-nh_v:,1-nh_v:,1-nh_v:)    :: s
    real(rp), intent(in ), dimension(1-nh_v:,1-nh_v:,1-nh_v:)    :: ls
    real(rp), intent(in ), dimension(1-nh_v:,1-nh_v:,1-nh_v:,1:) :: nor
    real(rp), intent(out), dimension(     1:,     1:,     1:)    :: dsdt
    !
    real(rp) :: ls_face,term_p,term_m,nor_face 
    integer  :: i,j,k
    !
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          !
          ! along x
          !
          ls_face  = 0.5_rp*(ls(i+1,j,k)+ls(i,j,k))
          nor_face = 0.5_rp*(nor(i+1,j,k,1)+nor(i,j,k,1))
          term_p   = u_reg*(eps_int*(s(i+1,j,k)-s(i,j,k))*dli(1) - &
                            0.25_rp*( 1._rp-tanh(ls_face/(2._rp*eps_int))**2 )*nor_face)
          ls_face  = 0.5_rp*(ls(i-1,j,k)+ls(i,j,k))
          nor_face = 0.5_rp*(nor(i-1,j,k,1)+nor(i,j,k,1))
          term_m   = u_reg*(eps_int*(s(i,j,k)-s(i-1,j,k))*dli(1) - &
                            0.25_rp*( 1._rp-tanh(ls_face/(2._rp*eps_int))**2 )*nor_face)
          !
          dsdt(i,j,k) = (term_p-term_m)*dli(1)
          !
          ! along y
          !
          ls_face  = 0.5_rp*(ls(i,j+1,k)+ls(i,j,k))
          nor_face = 0.5_rp*(nor(i,j+1,k,2)+nor(i,j,k,2))
          term_p   = u_reg*(eps_int*(s(i,j+1,k)-s(i,j,k))*dli(2) - &
                            0.25_rp*( 1._rp-tanh(ls_face/(2._rp*eps_int))**2 )*nor_face)
          ls_face  = 0.5*(ls(i,j-1,k)+ls(i,j,k))
          nor_face = 0.5*(nor(i,j-1,k,2)+nor(i,j,k,2))
          term_m   = u_reg*(eps_int*(s(i,j,k)-s(i,j-1,k))*dli(2) - &
                            0.25_rp*( 1._rp-tanh(ls_face/(2._rp*eps_int))**2 )*nor_face)
          !
          dsdt(i,j,k) = dsdt(i,j,k) + (term_p-term_m)*dli(2)
          !
          ! along z
          !
          ls_face  = 0.5_rp*(ls(i,j,k+1)+ls(i,j,k))
          nor_face = 0.5_rp*(nor(i,j,k+1,3)+nor(i,j,k,3))
          term_p   = u_reg*(eps_int*(s(i,j,k+1)-s(i,j,k))*dli(3) - &
                            0.25_rp*( 1._rp-tanh(ls_face/(2._rp*eps_int))**2 )*nor_face)
          ls_face  = 0.5_rp*(ls(i,j,k-1)+ ls(i,j,k))
          nor_face = 0.5_rp*(nor(i,j,k-1,3)+nor(i,j,k,3))
          term_m   = u_reg*(eps_int*(s(i,j,k)-s(i,j,k-1))*dli(3) - &
                            0.25_rp*( 1._rp-tanh(ls_face/(2._rp*eps_int))**2 )*nor_face)
          !
          dsdt(i,j,k) = dsdt(i,j,k) + (term_p-term_m)*dli(3)
          !
        enddo
      enddo
    enddo
    !
    return
  end subroutine regularization_psi
  !
  subroutine psi_to_ls(n,dl,nh_p,psi,eps_int,ls)
    !
    implicit none
    !
    integer , intent(in ), dimension(3) :: n
    real(rp), intent(in ), dimension(3) :: dl
    integer , intent(in )               :: nh_p
    real(rp), intent(in ), dimension(1-nh_p:,1-nh_p:,1-nh_p:) :: psi
    real(rp), intent(in )                                     :: eps_int
    real(rp), intent(out), dimension(1-nh_p:,1-nh_p:,1-nh_p:) :: ls
    !
    real(rp) :: clip_psi
    integer  :: i,j,k
    !
    real(rp), parameter :: small = real(2.e-16,rp)
    !
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          !
          clip_psi  = min(max(0.0_rp,psi(i,j,k)),1.0_rp) 
          ls(i,j,k) = eps_int*log( (clip_psi+small)/(1.d0-clip_psi+small) )
          !
        enddo
      enddo
    enddo
    !
    return
  end subroutine psi_to_ls
  !
  subroutine ls_to_psi(n,dl,nh_p,ls,eps_int,psi)
    !
    implicit none
    !
    integer , intent(in ), dimension(3) :: n
    real(rp), intent(in ), dimension(3) :: dl
    integer , intent(in )               :: nh_p
    real(rp), intent(in ), dimension(1-nh_p:,1-nh_p:,1-nh_p:) :: ls
    real(rp), intent(in )                                     :: eps_int
    real(rp), intent(out), dimension(1-nh_p:,1-nh_p:,1-nh_p:) :: psi
    !
    integer :: i,j,k
    !
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          psi(i,j,k) = 0.5_rp*( 1._rp+tanh(ls(i,j,k)/(2._rp*eps_int)) )
        enddo
      enddo
    enddo
    !
    return
  end subroutine ls_to_psi
  !
  subroutine cmpt_norm(n,dli,nh_p,phi,kappa,nor)
    !
    implicit none
    !
    integer , intent(in ), dimension(3) :: n
    real(rp), intent(in ), dimension(3) :: dli
    integer , intent(in )               :: nh_p
    real(rp), intent(in ), dimension(1-nh_p:,1-nh_p:,1-nh_p:)    :: phi
    real(rp), intent(out), dimension(     0:,     0:,     0:)    :: kappa
    real(rp), intent(out), dimension(1-nh_p:,1-nh_p:,1-nh_p:,1:) :: nor
    !
    real(rp), dimension(8) :: nx,ny,nz,mx,my,mz
    real(rp), dimension(3) :: dl
    real(rp) :: norm,cur_x,cur_y,cur_z
    integer  :: i,j,k,p
    !
    real(rp), parameter :: small = real(1.0e-12,rp)
    ! 
    dl(:) = dli(:)**(-1)
    !
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          !
          !i+1/2 j+1/2 k+1/2
          mx(1)=((phi(i+1,j  ,k  )+phi(i+1,j+1,k  )+phi(i+1,j  ,k+1)+phi(i+1,j+1,k+1))-&
                 (phi(i  ,j  ,k  )+phi(i  ,j+1,k  )+phi(i  ,j  ,k+1)+phi(i  ,j+1,k+1)))*dli(1)*0.25_rp
          !i+1/2 j-1/2 k+1/2
          mx(2)=((phi(i+1,j  ,k  )+phi(i+1,j-1,k  )+phi(i+1,j  ,k+1)+phi(i+1,j-1,k+1))-&
                 (phi(i  ,j  ,k  )+phi(i  ,j-1,k  )+phi(i  ,j  ,k+1)+phi(i  ,j-1,k+1)))*dli(1)*0.25_rp
          !i+1/2 j+1/2 k-1/2
          mx(3)=((phi(i+1,j  ,k  )+phi(i+1,j+1,k  )+phi(i+1,j  ,k-1)+phi(i+1,j+1,k-1))-&
                 (phi(i  ,j  ,k  )+phi(i  ,j+1,k  )+phi(i  ,j  ,k-1)+phi(i  ,j+1,k-1)))*dli(1)*0.25_rp
          !i+1/2 j-1/2 k-1/2
          mx(4)=((phi(i+1,j  ,k  )+phi(i+1,j-1,k  )+phi(i+1,j  ,k-1)+phi(i+1,j-1,k-1))-&
                 (phi(i  ,j  ,k  )+phi(i  ,j-1,k  )+phi(i  ,j  ,k-1)+phi(i  ,j-1,k-1)))*dli(1)*0.25_rp
          !i-1/2 j+1/2 k+1/2
          mx(5)=((phi(i  ,j  ,k  )+phi(i  ,j+1,k  )+phi(i  ,j  ,k+1)+phi(i  ,j+1,k+1))-&
                 (phi(i-1,j  ,k  )+phi(i-1,j+1,k  )+phi(i-1,j  ,k+1)+phi(i-1,j+1,k+1)))*dli(1)*0.25_rp
          !i-1/2 j-1/2 k+1/2
          mx(6)=((phi(i  ,j  ,k  )+phi(i  ,j-1,k  )+phi(i  ,j  ,k+1)+phi(i  ,j-1,k+1))-&
                 (phi(i-1,j  ,k  )+phi(i-1,j-1,k  )+phi(i-1,j  ,k+1)+phi(i-1,j-1,k+1)))*dli(1)*0.25_rp
          !i-1/2 j+1/2 k-1/2
          mx(7)=((phi(i  ,j  ,k  )+phi(i  ,j+1,k  )+phi(i  ,j  ,k-1)+phi(i  ,j+1,k-1))-&
                 (phi(i-1,j  ,k  )+phi(i-1,j+1,k  )+phi(i-1,j  ,k-1)+phi(i-1,j+1,k-1)))*dli(1)*0.25_rp
          !i-1/2 j-1/2 k-1/2
          mx(8)=((phi(i  ,j  ,k  )+phi(i  ,j-1,k  )+phi(i  ,j  ,k-1)+phi(i  ,j-1,k-1))-&
                 (phi(i-1,j  ,k  )+phi(i-1,j-1,k  )+phi(i-1,j  ,k-1)+phi(i-1,j-1,k-1)))*dli(1)*0.25_rp
          !
          !i+1/2 j+1/2 k+1/2
          my(1)=((phi(i  ,j+1,k  )+phi(i+1,j+1,k  )+phi(i  ,j+1,k+1)+phi(i+1,j+1,k+1))-&
                 (phi(i  ,j  ,k  )+phi(i+1,j  ,k  )+phi(i  ,j  ,k+1)+phi(i+1,j  ,k+1)))*dli(2)*0.25_rp
          !i+1/2 j-1/2 k+1/2
          my(2)=((phi(i  ,j  ,k  )+phi(i+1,j  ,k  )+phi(i  ,j  ,k+1)+phi(i+1,j  ,k+1))-&
                 (phi(i  ,j-1,k  )+phi(i+1,j-1,k  )+phi(i  ,j-1,k+1)+phi(i+1,j-1,k+1)))*dli(2)*0.25_rp
          !i+1/2 j+1/2 k-1/2
          my(3)=((phi(i  ,j+1,k  )+phi(i+1,j+1,k  )+phi(i  ,j+1,k-1)+phi(i+1,j+1,k-1))-&
                 (phi(i  ,j  ,k  )+phi(i+1,j  ,k  )+phi(i  ,j  ,k-1)+phi(i+1,j  ,k-1)))*dli(2)*0.25_rp
          !i+1/2 j-1/2 k-1/2
          my(4)=((phi(i  ,j  ,k  )+phi(i+1,j  ,k  )+phi(i  ,j  ,k-1)+phi(i+1,j  ,k-1))-&
                 (phi(i  ,j-1,k  )+phi(i+1,j-1,k  )+phi(i  ,j-1,k-1)+phi(i+1,j-1,k-1)))*dli(2)*0.25_rp
          !i-1/2 j+1/2 k+1/2
          my(5)=((phi(i  ,j+1,k  )+phi(i-1,j+1,k  )+phi(i  ,j+1,k+1)+phi(i-1,j+1,k+1))-&
                 (phi(i  ,j  ,k  )+phi(i-1,j  ,k  )+phi(i  ,j  ,k+1)+phi(i-1,j  ,k+1)))*dli(2)*0.25_rp
          !i-1/2 j-1/2 k+1/2
          my(6)=((phi(i  ,j  ,k  )+phi(i-1,j  ,k  )+phi(i  ,j  ,k+1)+phi(i-1,j  ,k+1))-&
                 (phi(i  ,j-1,k  )+phi(i-1,j-1,k  )+phi(i  ,j-1,k+1)+phi(i-1,j-1,k+1)))*dli(2)*0.25_rp
          !i-1/2 j+1/2 k-1/2
          my(7)=((phi(i  ,j+1,k  )+phi(i-1,j+1,k  )+phi(i  ,j+1,k-1)+phi(i-1,j+1,k-1))-&
                 (phi(i  ,j  ,k  )+phi(i-1,j  ,k  )+phi(i  ,j  ,k-1)+phi(i-1,j  ,k-1)))*dli(2)*0.25_rp
          !i-1/2 j-1/2 k-1/2
          my(8)=((phi(i  ,j  ,k  )+phi(i-1,j  ,k  )+phi(i  ,j  ,k-1)+phi(i-1,j  ,k-1))-&
                 (phi(i  ,j-1,k  )+phi(i-1,j-1,k  )+phi(i  ,j-1,k-1)+phi(i-1,j-1,k-1)))*dli(2)*0.25_rp
          !
          !i+1/2 j+1/2 k+1/2
          mz(1)=((phi(i  ,j  ,k+1)+phi(i+1,j  ,k+1)+phi(i  ,j+1,k+1)+phi(i+1,j+1,k+1))-&
                 (phi(i  ,j  ,k  )+phi(i+1,j  ,k  )+phi(i  ,j+1,k  )+phi(i+1,j+1,k  )))*dli(3)*0.25_rp
          !i+1/2 j-1/2 k+1/2
          mz(2)=((phi(i  ,j  ,k+1)+phi(i+1,j  ,k+1)+phi(i  ,j-1,k+1)+phi(i+1,j-1,k+1))-&
                 (phi(i  ,j  ,k  )+phi(i+1,j  ,k  )+phi(i  ,j-1,k  )+phi(i+1,j-1,k  )))*dli(3)*0.25_rp
          !i+1/2 j+1/2 k-1/2
          mz(3)=((phi(i  ,j  ,k  )+phi(i+1,j  ,k  )+phi(i  ,j+1,k  )+phi(i+1,j+1,k  ))-&
                 (phi(i  ,j  ,k-1)+phi(i+1,j  ,k-1)+phi(i  ,j+1,k-1)+phi(i+1,j+1,k-1)))*dli(3)*0.25_rp
          !i+1/2 j-1/2 k-1/2
          mz(4)=((phi(i  ,j  ,k  )+phi(i+1,j  ,k  )+phi(i  ,j-1,k  )+phi(i+1,j-1,k  ))-&
                 (phi(i  ,j  ,k-1)+phi(i+1,j  ,k-1)+phi(i  ,j-1,k-1)+phi(i+1,j-1,k-1)))*dli(3)*0.25_rp
          !i-1/2 j+1/2 k+1/2
          mz(5)=((phi(i  ,j  ,k+1)+phi(i-1,j  ,k+1)+phi(i  ,j+1,k+1)+phi(i-1,j+1,k+1))-&
                 (phi(i  ,j  ,k  )+phi(i-1,j  ,k  )+phi(i  ,j+1,k  )+phi(i-1,j+1,k  )))*dli(3)*0.25_rp
          !i-1/2 j-1/2 k+1/2
          mz(6)=((phi(i  ,j  ,k+1)+phi(i-1,j  ,k+1)+phi(i  ,j-1,k+1)+phi(i-1,j-1,k+1))-&
                 (phi(i  ,j  ,k  )+phi(i-1,j  ,k  )+phi(i  ,j-1,k  )+phi(i-1,j-1,k  )))*dli(3)*0.25_rp
          !i-1/2 j+1/2 k-1/2
          mz(7)=((phi(i  ,j  ,k  )+phi(i-1,j  ,k  )+phi(i  ,j+1,k  )+phi(i-1,j+1,k  ))-&
                 (phi(i  ,j  ,k-1)+phi(i-1,j  ,k-1)+phi(i  ,j+1,k-1)+phi(i-1,j+1,k-1)))*dli(3)*0.25_rp
          !i-1/2 j-1/2 k-1/2
          mz(8)=((phi(i  ,j  ,k  )+phi(i-1,j  ,k  )+phi(i  ,j-1,k  )+phi(i-1,j-1,k  ))-&
                 (phi(i  ,j  ,k-1)+phi(i-1,j  ,k-1)+phi(i  ,j-1,k-1)+phi(i-1,j-1,k-1)))*dli(3)*0.25_rp
          !
          do p=1,8
            norm  = sqrt(mx(p)**2+my(p)**2+mz(p)**2) + small
            nx(p) = mx(p)/norm
            ny(p) = my(p)/norm
            nz(p) = mz(p)/norm
          enddo
          !
          ! compute the normal vector
          !
          nor(i,j,k,1) = 0._rp
          nor(i,j,k,2) = 0._rp
          nor(i,j,k,3) = 0._rp
          do p=1,8
            nor(i,j,k,1) = nor(i,j,k,1) + 0.125_rp*mx(p)
            nor(i,j,k,2) = nor(i,j,k,2) + 0.125_rp*my(p)
            nor(i,j,k,3) = nor(i,j,k,3) + 0.125_rp*mz(p)
          enddo
          !
          norm         = sqrt(nor(i,j,k,1)**2+nor(i,j,k,2)**2+nor(i,j,k,3)**2) + small
          nor(i,j,k,1) = nor(i,j,k,1)/norm
          nor(i,j,k,2) = nor(i,j,k,2)/norm
          nor(i,j,k,3) = nor(i,j,k,3)/norm
          ! 
          cur_x = ((nx(1)+nx(2)+nx(3)+nx(4))-(nx(5)+nx(6)+nx(7)+nx(8)))*dl(1)*0.25_rp
          cur_y = ((ny(1)+ny(3)+ny(5)+ny(7))-(ny(2)+ny(4)+ny(6)+ny(8)))*dl(2)*0.25_rp
          cur_z = ((nz(1)+nz(2)+nz(5)+nz(6))-(nz(3)+nz(4)+nz(7)+nz(8)))*dl(3)*0.25_rp
          !     
          kappa(i,j,k) = -(cur_x*dli(1)**2+cur_y*dli(2)**2+cur_z*dli(3)**2)
          !
        enddo
      enddo
    enddo
    !
    return
  end subroutine cmpt_norm
  !
  subroutine update_property(n,prop12,psi,prop)
    !
    implicit none
    !
    integer , intent(in ), dimension(3)        :: n
    real(rp), intent(in ), dimension(2)        :: prop12
    real(rp), intent(in ), dimension(0:,0:,0:) :: psi
    real(rp), intent(out), dimension(0:,0:,0:) :: prop
    !
    real(rp) :: prop12_1, prop12_2
    integer  :: i,j,k
    integer  :: n1, n2, n3
    ! 
    n1 = n(1)
    n2 = n(2)
    n3 = n(3)
    !
    prop12_1 = prop12(1)
    prop12_2 = prop12(2)
    !
    !$acc parallel loop collapse(3)
    do k=1,n3
      do j=1,n2
        do i=1,n1
          prop(i,j,k) = psi(i,j,k)*prop12_1+(1.0_rp-psi(i,j,k))*prop12_2
        enddo
      enddo
    enddo
    !$acc end parallel loop 
    !
    return
  end subroutine update_property
  !
  subroutine initls(inils,n,dl,l,nh_p,ls)
    !
    use mod_param, only: xc,yc,zc,r,nbub
    !
    implicit none
    !
    character(len=3), intent(in )                                     :: inils
    integer         , intent(in ), dimension(3)                       :: n
    real(rp)        , intent(in ), dimension(3)                       :: dl
    real(rp)        , intent(in ), dimension(3)                       :: l
    integer         , intent(in )                                     :: nh_p
    real(rp)        , intent(out), dimension(1-nh_p:,1-nh_p:,1-nh_p:) :: ls
    !
    real(rp) :: phi
    real(rp) :: d,x,y,z,ry
    real(rp) :: rho1,rho2
    integer  :: i,j,k
    !
    !ry = 0.5d0*d0(1)
    ry = r(1)
    !
    select case(inils)
    case('tgn')
      do k=1,n(3)
        z = (k+ijk_start(3)-0.5)*dl(3)
        do j=1,n(2)
          y = (j+ijk_start(2)-0.5)*dl(2)
          do i=1,n(1)
            x = (i+ijk_start(1)-0.5)*dl(1)
            !
#if defined(_TWOD)
            d = sqrt(                (y-yc(1))**2.0+(z-zc(1))**2.0 )
#else
            d = sqrt( (x-xc(1))**2.0+(y-yc(1))**2.0+(z-zc(1))**2.0 )
#endif
            !
            ls(i,j,k) = ry-d
            !
          enddo
        enddo
      enddo
    case default
      call flutas_error('ERROR: invalid name for initial level-set field: check di_acdi.in')
    end select
    !
    return
  end subroutine initls
  !
  subroutine initpsi(inipsi,n,dl,l,nh_v,eps_int,psi)
    !
    use mod_param, only: xc,yc,zc,r,nbub
    !
    implicit none
    !
    character(len=3), intent(in )                                     :: inipsi
    integer         , intent(in ), dimension(3)                       :: n
    real(rp)        , intent(in ), dimension(3)                       :: dl
    real(rp)        , intent(in ), dimension(3)                       :: l
    integer         , intent(in )                                     :: nh_v
    real(rp)        , intent(in )                                     :: eps_int
    real(rp)        , intent(out), dimension(1-nh_v:,1-nh_v:,1-nh_v:) :: psi
    !
    real(rp) :: d,x,y,z
    integer  :: i,j,k,p
    !
    select case(inipsi)
    case('bub')
      do k=1,n(3)
        z = (k+ijk_start(3)-0.5_rp)*dl(3)
        do j=1,n(2)
          y = (j+ijk_start(2)-0.5_rp)*dl(2)
          do i=1,n(1)
            x = (i+ijk_start(1)-0.5_rp)*dl(1)
             !
             psi(i,j,k) = 0._rp
             do p=1,nbub
#if defined(_TWOD)
               d = sqrt(                  (y-yc(p))**2._rp+(z-zc(p))**2._rp )
#else
               d = sqrt( (x-xc(p))**2._rp+(y-yc(p))**2._rp+(z-zc(p))**2._rp )
#endif
               psi(i,j,k) = psi(i,j,k) + 0.5_rp*(1._rp+tanh( (r(p)-d)/(2._rp*eps_int) ))
             enddo
             !
          enddo
        enddo
      enddo
    case default
      call flutas_error('ERROR: invalid name for initial DI field: check di_acdi.in')
    end select
    !
    return
  end subroutine initpsi
    !
    !
    subroutine initc1(iniy1,n,dl,l,nh_v,eps_int,phi,c1,c2)
    !
    use mod_common_mpi, only: ijk_start
    !
    use mod_param, only: xc,yc,zc,r
    !
    implicit none
    !
    character(len=3), intent(in )                         :: iniy1
    integer         , intent(in ), dimension(3)           :: n
    real(rp)        , intent(in ), dimension(3)           :: dl
    real(rp)        , intent(in ), dimension(3)           :: l
    integer         , intent(in )                         :: nh_v
    real(rp)        , intent(in )                         :: eps_int
    real(rp)        , intent(in ), dimension(1-nh_v:,1-nh_v:,1-nh_v:)  :: phi
    real(rp)        , intent(out), dimension(1-nh_v:,1-nh_v:,1-nh_v:) :: c1,c2
    !
    real(rp) :: d,x,y,z,ry
    integer  :: i,j,k
    !
    ry = r(1)
    !
    select case(iniy1)
    case ('lin')
    do k=1,n(3)
      !
      z = (k+ijk_start(3)-0.5)*dl(3) 
      !
      do j=1,n(2)
      !
      y = (j+ijk_start(2)-0.5)*dl(2)
      !
        do i=1,n(1)
      !
      x = (i+ijk_start(1)-0.5)*dl(1)
      !
      c1(i,j,k) = phi(i,j,k)
          !
        enddo
      enddo
    enddo
    c2(:,:,:) = 0.d0
    end select
    !
    return
    !
  end subroutine initc1
    !
end module mod_di_acdi
