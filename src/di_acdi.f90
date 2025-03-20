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
  public :: update_property,rk_psi,rk_ent,psi_to_ls,ls_to_psi,cmpt_norm,initls,initpsi,cmpt_umax,inittmp
  public rk_c1,initc1,cmpt_rel_vel, interfacial_cflux,extended

  public static_contact_angle_m


  public cmpt_delta,marangoni_he
  !
  contains
    !
    !
  subroutine rk_ent(f_t1,f_t2,n,dli,nh_d,nh_u,nh_v,dzc,dzf,eps_int,u_reg, &
                      u,v,w,s,nor,tmp,cond,mfx,mfy,mfz,ent,dentdtrko)
    !
    ! low-storage 3rd-order Runge-Kutta scheme 
    ! for time integration of the Enthalpy field.
    !
    implicit none
    !
    real(rp), intent(in   )                                        :: f_t1,f_t2
    real(rp), intent(in   )                                        :: eps_int,u_reg
    integer , intent(in   ), dimension(3)                          :: n
    real(rp), intent(in   ), dimension(3)                          :: dli
    integer , intent(in   )                                        :: nh_d,nh_u,nh_v
    real(rp), intent(in   ), dimension(1-nh_d:)                    :: dzc,dzf
    real(rp), intent(in   ), dimension(1-nh_u:,1-nh_u:,1-nh_u:)    :: u,v,w,mfx,mfy,mfz
    real(rp), intent(in   ), dimension(1-nh_v:,1-nh_v:,1-nh_v:)    :: cond,s 
    real(rp), intent(inout   ), dimension(1-nh_v:,1-nh_v:,1-nh_v:)    :: tmp
    real(rp), intent(in ), dimension(0:,0:,0:,1:)                  :: nor
    real(rp), intent(inout), dimension(1-nh_v:,1-nh_v:,1-nh_v:)    :: ent  
    real(rp), intent(inout), dimension(     1:,     1:,     1:)    :: dentdtrko
    !
    real(rp) :: gradpxp,gradpxm,gradpyp,gradpym,gradpzp,gradpzm
    real(rp) :: psifxp,psifxm,psifyp,psifym,psifzp,psifzm
    real(rp) :: regxp,regxm,regyp,regym,regzp,regzm
    real(rp) :: tmpxp,tmpxm,tmpyp,tmpym,tmpzp,tmpzm
    !
    real(rp), dimension(n(1),n(2),n(3)) :: dentdtrk
    real(rp)  ::mxp,mxm,myp,mym,mzp,mzm
    integer  :: i,j,k,ip,im,jp,jm,kp,km,zp,zm
    !
    ! 1. compute diffusion
    !
    call diff_ent(n,dli,nh_d,nh_v,dzc,dzf,cond,tmp,dentdtrk)
    !
    !2 . Advection 
    call moms_reg(n,dli,nh_d,nh_u,nh_v,u,v,w,mfx,mfy,mfz,ent,dentdtrk)
    !
    ! 3. Advance enthalpy 
    
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
        !
        ent(i,j,k)       = ent(i,j,k) + f_t1*dentdtrk(i,j,k) + f_t2*dentdtrko(i,j,k)
        dentdtrko(i,j,k) = dentdtrk(i,j,k)
        !
        enddo
      enddo
    enddo
    !
    return
  end subroutine rk_ent
  !
    subroutine diff_ent(n,dli,nh_d,nh_v,dzc,dzf,cond,s,dsdt) 
    !
    implicit none
    !
    integer , intent(in ), dimension(3)               :: n
    real(rp), intent(in ), dimension(3)               :: dli
    integer , intent(in )                             :: nh_d,nh_v
    real(rp), intent(in ), dimension(1-nh_d:)         :: dzc,dzf
    real(rp), intent(in ), dimension(1-nh_v:,1-nh_v:,1-nh_v:)    :: s,cond
    real(rp), intent(out), dimension(1:,1:,1:)                   :: dsdt
    !
    real(rp), dimension(3) :: dl
    real(rp) :: sfxp,sfxm,sfyp,sfym,sfzp,sfzm
    integer  :: i,j,k,ip,im,jp,jm,kp,km
    !
    dl(:) = 1.d0/dli(:)
    dsdt(:,:,:) = 0.d0
    !
    !
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
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
          sfxp = (cond(ip,j,k)+cond(i,j,k))*(s(ip,j,k)-s(i,j,k))*dli(1)*0.5
          sfxm = (cond(im,j,k)+cond(i,j,k))*(s(i,j,k)-s(im,j,k))*dli(1)*0.5
          !
          dsdt(i,j,k) = (sfxp-sfxm)*dli(1)
          !
          !
          ! along y-dir
          !
          sfyp = (cond(i,jp,k)+cond(i,j,k))*(s(i,jp,k)-s(i,j,k))*dli(2)*0.5
          sfym = (cond(i,jm,k)+cond(i,j,k))*(s(i,j,k)-s(i,jm,k))*dli(2)*0.5
          !
          dsdt(i,j,k) = dsdt(i,j,k) + (sfyp-sfym)*dli(2) 
          !
          !
          ! along z-dir
          !
          sfzp = (cond(i,j,kp)+cond(i,j,k))*(s(i,j,kp)-s(i,j,k))*dli(3)*0.5
          sfzm = (cond(i,j,km)+cond(i,j,k))*(s(i,j,k)-s(i,j,km))*dli(3)*0.5
          !
          dsdt(i,j,k)  = dsdt(i,j,k)+ (sfzp-sfzm)*dli(3)
          !
          !
        enddo
      enddo
    enddo
    !
    return
  end subroutine diff_ent

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

  subroutine cmpt_delta(n,dli,phi,delta)
    !
    implicit none
    !
    integer , intent(in ), dimension(3) :: n
    real(rp), intent(in ), dimension(3) :: dli
    real(rp), intent(in ), dimension(0:,0:,0:) :: phi
    real(rp), intent(out), dimension(0:,0:,0:) :: delta
    !
    real(rp), dimension(8) :: nx,ny,nz,mx,my,mz
    real(rp), dimension(3) :: dl
    real(rp) :: dphidx,dphidy,dphidz
    real(rp) :: small = 1.e-12
    integer  :: i,j,k,p
    ! 
    dl(:) = dli(:)**(-1)
    !
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          !
#ifdef TWOD 
          mx(1:8) = 0.d0
#else
          !i+1/2 j+1/2 k+1/2
          mx(1)=((phi(i+1,j  ,k  )+phi(i+1,j+1,k  )+phi(i+1,j  ,k+1)+phi(i+1,j+1,k+1))-&
                 (phi(i  ,j  ,k  )+phi(i  ,j+1,k  )+phi(i  ,j  ,k+1)+phi(i  ,j+1,k+1)))*dli(1)*0.25
          !i+1/2 j-1/2 k+1/2
          mx(2)=((phi(i+1,j  ,k  )+phi(i+1,j-1,k  )+phi(i+1,j  ,k+1)+phi(i+1,j-1,k+1))-&
                 (phi(i  ,j  ,k  )+phi(i  ,j-1,k  )+phi(i  ,j  ,k+1)+phi(i  ,j-1,k+1)))*dli(1)*0.25
          !i+1/2 j+1/2 k-1/2
          mx(3)=((phi(i+1,j  ,k  )+phi(i+1,j+1,k  )+phi(i+1,j  ,k-1)+phi(i+1,j+1,k-1))-&
                 (phi(i  ,j  ,k  )+phi(i  ,j+1,k  )+phi(i  ,j  ,k-1)+phi(i  ,j+1,k-1)))*dli(1)*0.25
          !i+1/2 j-1/2 k-1/2
          mx(4)=((phi(i+1,j  ,k  )+phi(i+1,j-1,k  )+phi(i+1,j  ,k-1)+phi(i+1,j-1,k-1))-&
                 (phi(i  ,j  ,k  )+phi(i  ,j-1,k  )+phi(i  ,j  ,k-1)+phi(i  ,j-1,k-1)))*dli(1)*0.25
          !i-1/2 j+1/2 k+1/2
          mx(5)=((phi(i  ,j  ,k  )+phi(i  ,j+1,k  )+phi(i  ,j  ,k+1)+phi(i  ,j+1,k+1))-&
                 (phi(i-1,j  ,k  )+phi(i-1,j+1,k  )+phi(i-1,j  ,k+1)+phi(i-1,j+1,k+1)))*dli(1)*0.25
          !i-1/2 j-1/2 k+1/2
          mx(6)=((phi(i  ,j  ,k  )+phi(i  ,j-1,k  )+phi(i  ,j  ,k+1)+phi(i  ,j-1,k+1))-&
                 (phi(i-1,j  ,k  )+phi(i-1,j-1,k  )+phi(i-1,j  ,k+1)+phi(i-1,j-1,k+1)))*dli(1)*0.25
          !i-1/2 j+1/2 k-1/2
          mx(7)=((phi(i  ,j  ,k  )+phi(i  ,j+1,k  )+phi(i  ,j  ,k-1)+phi(i  ,j+1,k-1))-&
                 (phi(i-1,j  ,k  )+phi(i-1,j+1,k  )+phi(i-1,j  ,k-1)+phi(i-1,j+1,k-1)))*dli(1)*0.25
          !i-1/2 j-1/2 k-1/2
          mx(8)=((phi(i  ,j  ,k  )+phi(i  ,j-1,k  )+phi(i  ,j  ,k-1)+phi(i  ,j-1,k-1))-&
                 (phi(i-1,j  ,k  )+phi(i-1,j-1,k  )+phi(i-1,j  ,k-1)+phi(i-1,j-1,k-1)))*dli(1)*0.25
          !
#endif 
         my(1)=((phi(i  ,j+1,k  )+phi(i+1,j+1,k  )+phi(i  ,j+1,k+1)+phi(i+1,j+1,k+1))-&
                 (phi(i  ,j  ,k  )+phi(i+1,j  ,k  )+phi(i  ,j  ,k+1)+phi(i+1,j  ,k+1)))*dli(2)*0.25
          !i+1/2 j-1/2 k+1/2
          my(2)=((phi(i  ,j  ,k  )+phi(i+1,j  ,k  )+phi(i  ,j  ,k+1)+phi(i+1,j  ,k+1))-&
                 (phi(i  ,j-1,k  )+phi(i+1,j-1,k  )+phi(i  ,j-1,k+1)+phi(i+1,j-1,k+1)))*dli(2)*0.25
          !i+1/2 j+1/2 k-1/2
          my(3)=((phi(i  ,j+1,k  )+phi(i+1,j+1,k  )+phi(i  ,j+1,k-1)+phi(i+1,j+1,k-1))-&
                 (phi(i  ,j  ,k  )+phi(i+1,j  ,k  )+phi(i  ,j  ,k-1)+phi(i+1,j  ,k-1)))*dli(2)*0.25
          !i+1/2 j-1/2 k-1/2
          my(4)=((phi(i  ,j  ,k  )+phi(i+1,j  ,k  )+phi(i  ,j  ,k-1)+phi(i+1,j  ,k-1))-&
                 (phi(i  ,j-1,k  )+phi(i+1,j-1,k  )+phi(i  ,j-1,k-1)+phi(i+1,j-1,k-1)))*dli(2)*0.25
          !i-1/2 j+1/2 k+1/2
          my(5)=((phi(i  ,j+1,k  )+phi(i-1,j+1,k  )+phi(i  ,j+1,k+1)+phi(i-1,j+1,k+1))-&
                 (phi(i  ,j  ,k  )+phi(i-1,j  ,k  )+phi(i  ,j  ,k+1)+phi(i-1,j  ,k+1)))*dli(2)*0.25
          !i-1/2 j-1/2 k+1/2
          my(6)=((phi(i  ,j  ,k  )+phi(i-1,j  ,k  )+phi(i  ,j  ,k+1)+phi(i-1,j  ,k+1))-&
                 (phi(i  ,j-1,k  )+phi(i-1,j-1,k  )+phi(i  ,j-1,k+1)+phi(i-1,j-1,k+1)))*dli(2)*0.25
          !i-1/2 j+1/2 k-1/2
          my(7)=((phi(i  ,j+1,k  )+phi(i-1,j+1,k  )+phi(i  ,j+1,k-1)+phi(i-1,j+1,k-1))-&
                 (phi(i  ,j  ,k  )+phi(i-1,j  ,k  )+phi(i  ,j  ,k-1)+phi(i-1,j  ,k-1)))*dli(2)*0.25
          !i-1/2 j-1/2 k-1/2
          my(8)=((phi(i  ,j  ,k  )+phi(i-1,j  ,k  )+phi(i  ,j  ,k-1)+phi(i-1,j  ,k-1))-&
                 (phi(i  ,j-1,k  )+phi(i-1,j-1,k  )+phi(i  ,j-1,k-1)+phi(i-1,j-1,k-1)))*dli(2)*0.25
          !
          !i+1/2 j+1/2 k+1/2
          mz(1)=((phi(i  ,j  ,k+1)+phi(i+1,j  ,k+1)+phi(i  ,j+1,k+1)+phi(i+1,j+1,k+1))-&
                 (phi(i  ,j  ,k  )+phi(i+1,j  ,k  )+phi(i  ,j+1,k  )+phi(i+1,j+1,k  )))*dli(3)*0.25
          !i+1/2 j-1/2 k+1/2
          mz(2)=((phi(i  ,j  ,k+1)+phi(i+1,j  ,k+1)+phi(i  ,j-1,k+1)+phi(i+1,j-1,k+1))-&
                 (phi(i  ,j  ,k  )+phi(i+1,j  ,k  )+phi(i  ,j-1,k  )+phi(i+1,j-1,k  )))*dli(3)*0.25
          !i+1/2 j+1/2 k-1/2
          mz(3)=((phi(i  ,j  ,k  )+phi(i+1,j  ,k  )+phi(i  ,j+1,k  )+phi(i+1,j+1,k  ))-&
                 (phi(i  ,j  ,k-1)+phi(i+1,j  ,k-1)+phi(i  ,j+1,k-1)+phi(i+1,j+1,k-1)))*dli(3)*0.25
          !i+1/2 j-1/2 k-1/2
          mz(4)=((phi(i  ,j  ,k  )+phi(i+1,j  ,k  )+phi(i  ,j-1,k  )+phi(i+1,j-1,k  ))-&
                 (phi(i  ,j  ,k-1)+phi(i+1,j  ,k-1)+phi(i  ,j-1,k-1)+phi(i+1,j-1,k-1)))*dli(3)*0.25
          !i-1/2 j+1/2 k+1/2
          mz(5)=((phi(i  ,j  ,k+1)+phi(i-1,j  ,k+1)+phi(i  ,j+1,k+1)+phi(i-1,j+1,k+1))-&
                 (phi(i  ,j  ,k  )+phi(i-1,j  ,k  )+phi(i  ,j+1,k  )+phi(i-1,j+1,k  )))*dli(3)*0.25
          !i-1/2 j-1/2 k+1/2
          mz(6)=((phi(i  ,j  ,k+1)+phi(i-1,j  ,k+1)+phi(i  ,j-1,k+1)+phi(i-1,j-1,k+1))-&
                 (phi(i  ,j  ,k  )+phi(i-1,j  ,k  )+phi(i  ,j-1,k  )+phi(i-1,j-1,k  )))*dli(3)*0.25
          !i-1/2 j+1/2 k-1/2
          mz(7)=((phi(i  ,j  ,k  )+phi(i-1,j  ,k  )+phi(i  ,j+1,k  )+phi(i-1,j+1,k  ))-&
                 (phi(i  ,j  ,k-1)+phi(i-1,j  ,k-1)+phi(i  ,j+1,k-1)+phi(i-1,j+1,k-1)))*dli(3)*0.25
          !i-1/2 j-1/2 k-1/2
          mz(8)=((phi(i  ,j  ,k  )+phi(i-1,j  ,k  )+phi(i  ,j-1,k  )+phi(i-1,j-1,k  ))-&
                 (phi(i  ,j  ,k-1)+phi(i-1,j  ,k-1)+phi(i  ,j-1,k-1)+phi(i-1,j-1,k-1)))*dli(3)*0.25
          !
          ! compute the delta
          !
          dphidx = 0._rp
          dphidy = 0._rp
          dphidz = 0._rp
          do p=1,8
            dphidx = dphidx + 0.125_rp*mx(p)
            dphidy = dphidy + 0.125_rp*my(p)
            dphidz = dphidz + 0.125_rp*mz(p)
          enddo
          !
          delta(i,j,k) = sqrt(dphidx**2+dphidy**2+dphidz**2+small)
          !
        enddo
      enddo
    enddo
    !
    return
  end subroutine cmpt_delta
    !
    !
  subroutine marangoni_he(n,dli,dzc,dzf,tmp,nh_t,delta_arr,kappa,nor,a1,rho,mar_force)
    !
    use mod_param, only: sigma_t,sigma,tl0,tg0
    !
    implicit none
    !
    integer , intent(in ), dimension(3)                         :: n
    real(rp), intent(in ), dimension(3)                         :: dli
    real(rp), intent(in ), dimension(-2:)                       :: dzc,dzf
    integer , intent(in   )                                     :: nh_t
    real(rp), intent(in ), dimension(1-nh_t:,1-nh_t:,1-nh_t:)   :: tmp
    real(rp), intent(in ), dimension( 0:, 0:, 0:)               :: delta_arr,kappa
    real(rp), intent(in ), dimension(0:,0:,0:,1:)               :: nor
    real(rp), intent(in ), dimension(0:,0:,0:)                  :: rho
    real(rp), intent(in ), dimension(0:,0:,0:)                  :: a1
    real(rp), intent(out), dimension(-2:,-2:,-2:,1:)            :: mar_force
    !
    real(rp) :: tmpf,kappas
    real(rp) :: norx,nory,norz
    real(rp) :: dpsidx,dpsidy,dpsidz
    real(rp) :: dtmpdx,dtmpdy,dtmpdz,dt_nof
    real(rp) :: rhox, rhoy, rhoz, rhoxi, rhoyi, rhozi
    real(rp) :: tmp_min,tmp_max,sigma_min,sigma_max
    real(rp) :: am, ap, delta, deltas
    real(rp), dimension(3) :: dl
    real(rp), dimension(n(1),n(2),n(3)) :: sigma_arr
    real(rp):: tmp1
    integer  :: i,j,k,p
    !
    tmp1=tl0
    dl(:) = 1.d0/dli(:)
    !
    do k=1,n(3)
      do j=1,n(2)
         do i=1,n(1)
           !
           !
#ifdef TWOD
           mar_force(i,j,k,1) = 0.d0
#else
           !
           ! X-direction (three contributions)
           !
           tmpf   = 0.5*(tmp(  i+1,j,k  )+tmp(  i,j,k  ))
           kappas = 0.5*(kappa(i+1,j,k  )+kappa(i,j,k  ))
           norx   = 0.5*(nor(  i+1,j,k,1)+nor(  i,j,k,1))
           nory   = 0.5*(nor(  i+1,j,k,2)+nor(  i,j,k,2))
           norz   = 0.5*(nor(  i+1,j,k,3)+nor(  i,j,k,3))
           dpsidx = (a1( i+1,j,k)-a1( i,j,k))*dli(1)
           !
           rhox   = 0.5d0*(rho(i+1,j,k)+rho(i,j,k))
           rhoxi  = 1.d0 / rhox
           deltas = 0.5d0*(delta_arr(i+1,j,k)+delta_arr(i,j,k))
           !
           dtmpdx = (tmp(i+1,j,k)-tmp(i,j,k))*dli(1)
           dtmpdy = 0.25*( (tmp(i,j,k)+tmp(i+1,j,k)+tmp(i+1,j+1,k)+tmp(i,j+1,k)) - &
                           (tmp(i,j,k)+tmp(i+1,j,k)+tmp(i+1,j-1,k)+tmp(i,j-1,k)) )*dli(2)
           dtmpdz = 0.25*( (tmp(i,j,k)+tmp(i+1,j,k)+tmp(i+1,j,k+1)+tmp(i,j,k+1)) - &
                           (tmp(i,j,k)+tmp(i+1,j,k)+tmp(i+1,j,k-1)+tmp(i,j,k-1)) )*dli(3)
           dt_nof = dtmpdx*norx+dtmpdy*nory+dtmpdz*norz
           !
           ap     =  a1(i+1,j,k)
           am     =  a1(i  ,j,k)
           delta  =  ((ap-am)*dli(1))**2
           ap     =  0.25*(a1(i,j,k)+a1(i+1,j,k)+a1(i+1,j+1,k)+a1(i,j+1,k))
           am     =  0.25*(a1(i,j,k)+a1(i+1,j,k)+a1(i+1,j-1,k)+a1(i,j-1,k))
           delta  =  delta+((ap-am)*dli(2))**2
           ap     =  0.25*(a1(i,j,k)+a1(i+1,j,k)+a1(i+1,j,k+1)+a1(i,j,k+1))
           am     =  0.25*(a1(i,j,k)+a1(i+1,j,k)+a1(i+1,j,k-1)+a1(i,j,k-1))
           delta  =  sqrt(delta+((ap-am)*dli(3))**2)
           !
           !mar_force(i,j,k,1) = sigma_t*( (tmpf-tmp1)*kappas*dpsidx + dtmpdx*delta - dt_nof*dpsidx )*rhoxi
           mar_force(i,j,k,1) = sigma_t*( (tmpf-tmp1)*kappas*dpsidx + dtmpdx*deltas - dt_nof*dpsidx )*rhoxi
#endif
           !
           !
           ! Y-direction (three contribution)
           !
           tmpf   = 0.5*(tmp(  i,j+1,k  )+tmp(  i,j,k))
           kappas = 0.5*(kappa(i,j+1,k  )+kappa(i,j,k))
           norx   = 0.5*(nor(  i,j+1,k,1)+nor(  i,j,k,1))
           nory   = 0.5*(nor(  i,j+1,k,2)+nor(  i,j,k,2))
           norz   = 0.5*(nor(  i,j+1,k,3)+nor(  i,j,k,3))
           dpsidy = (a1( i,j+1,k)-a1( i,j,k))*dli(1)
           !
           rhoy   = 0.5d0*(rho(i,j+1,k)+rho(i,j,k))
           rhoyi  = 1.d0 / rhoy
           deltas = 0.5d0*(delta_arr(i,j+1,k)+delta_arr(i,j,k))
           !
           dtmpdx = 0.25*( (tmp(i,j,k)+tmp(i+1,j,k)+tmp(i+1,j+1,k)+tmp(i,j+1,k)) - &
                           (tmp(i,j,k)+tmp(i-1,j,k)+tmp(i-1,j+1,k)+tmp(i,j+1,k)) )*dli(1)
           dtmpdy = (tmp(i,j+1,k)-tmp(i,j,k))*dli(2)
           dtmpdz = 0.25*( (tmp(i,j,k)+tmp(i,j+1,k)+tmp(i,j+1,k+1)+tmp(i,j,k+1)) - &
                           (tmp(i,j,k)+tmp(i,j+1,k)+tmp(i,j+1,k+1)+tmp(i,j,k+1)) )*dli(3)
           dt_nof = dtmpdx*norx+dtmpdy*nory+dtmpdz*norz
           !
           ap     =  0.25*(a1(i,j,k)+a1(i+1,j,k)+a1(i+1,j+1,k)+a1(i,j+1,k))
           am     =  0.25*(a1(i,j,k)+a1(i-1,j,k)+a1(i-1,j+1,k)+a1(i,j+1,k))
           delta  =  ((ap-am)*dli(1))**2
           ap     =  a1(i,j+1,k)
           am     =  a1(i,j  ,k)
           delta  =  delta+((ap-am)*dli(2))**2
           ap     =  0.25*(a1(i,j,k)+a1(i,j+1,k)+a1(i,j+1,k+1)+a1(i,j,k+1))
           am     =  0.25*(a1(i,j,k)+a1(i,j+1,k)+a1(i,j+1,k-1)+a1(i,j,k-1))
           delta  =  sqrt(delta+((ap-am)*dli(3))**2)
           !
           !
           !mar_force(i,j,k,2) = sigma_t*( (tmpf-tmp1)*kappas*dpsidy + dtmpdy*delta - dt_nof*dpsidy )*rhoyi
           mar_force(i,j,k,2) = sigma_t*( (tmpf-tmp1)*kappas*dpsidy + dtmpdz*deltas - dt_nof*dpsidy )*rhoyi
           !
           ! Z-contribution (three contribution)
           !
           tmpf   = 0.5*(tmp(  i,j,k+1  )+tmp(  i,j,k))
           kappas = 0.5*(kappa(i,j,k+1  )+kappa(i,j,k))
           norx   = 0.5*(nor(  i,j,k+1,1)+nor(  i,j,k,1))
           nory   = 0.5*(nor(  i,j,k+1,2)+nor(  i,j,k,2))
           norz   = 0.5*(nor(  i,j,k+1,3)+nor(  i,j,k,3))
           dpsidz = (a1( i,j,k+1)-a1( i,j,k))*dli(1)
           !
           rhoz   = 0.5d0*(rho(i,j,k+1)+rho(i,j,k))
           rhozi  = 1.d0 / rhoz
           deltas = 0.5d0*(delta_arr(i,j,k+1)+delta_arr(i,j,k))
           !
           dtmpdx = 0.25*( (tmp(i,j,k)+tmp(i+1,j,k)+tmp(i+1,j,k+1)+tmp(i,j,k+1)) - &
                           (tmp(i,j,k)+tmp(i-1,j,k)+tmp(i-1,j,k+1)+tmp(i,j,k+1)) )*dli(1)
           dtmpdy = 0.25*( (tmp(i,j,k)+tmp(i,j+1,k)+tmp(i,j+1,k+1)+tmp(i,j,k+1)) - &
                           (tmp(i,j,k)+tmp(i,j-1,k)+tmp(i,j-1,k+1)+tmp(i,j,k+1)) )*dli(2)
           dtmpdz = (tmp(i,j,k+1)-tmp(i,j,k))*dli(3)
           dt_nof = dtmpdx*norx+dtmpdy*nory+dtmpdz*norz
           !
           ap     =  0.25*(a1(i,j,k)+a1(i+1,j,k)+a1(i+1,j,k+1)+a1(i,j,k+1))
           am     =  0.25*(a1(i,j,k)+a1(i-1,j,k)+a1(i-1,j,k+1)+a1(i,j,k+1))
           delta  =  ((ap-am)*dli(1))**2
           ap     =  0.25*(a1(i,j,k)+a1(i,j+1,k)+a1(i,j+1,k+1)+a1(i,j,k+1))
           am     =  0.25*(a1(i,j,k)+a1(i,j-1,k)+a1(i,j-1,k+1)+a1(i,j,k+1))
           delta  =  delta+((ap-am)*dli(2))**2
           ap     =  a1(i,j,k+1)
           am     =  a1(i,j,k)
           delta  =  sqrt(delta+((ap-am)*dli(3))**2)
           !
           !mar_force(i,j,k,3) = sigma_t*( (tmpf-tmp1)*kappas*dpsidz + dtmpdz*delta - dt_nof*dpsidz )*rhozi
           mar_force(i,j,k,3) = sigma_t*( (tmpf-tmp1)*kappas*dpsidz + dtmpdz*deltas - dt_nof*dpsidz )*rhozi
           !
         enddo
      enddo
    enddo
    !
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          sigma_arr(i,j,k) = sigma + sigma_t*(tmp(i,j,k)-tmp1)
        enddo
      enddo
    enddo
    !
    sigma_min = minval(sigma_arr(1:n(1),1:n(2),1:n(3)))
    call mpi_allreduce(MPI_IN_PLACE,sigma_min,1,mpi_real8,mpi_min,comm_cart,ierr)
    sigma_max = maxval(sigma_arr(1:n(1),1:n(2),1:n(3)))
    call mpi_allreduce(MPI_IN_PLACE,sigma_max,1,mpi_real8,mpi_max,comm_cart,ierr)
    !
    if(myid.eq.0) print*, "sigma_min", sigma_min
    if(myid.eq.0) print*, "sigma_max", sigma_max
    !
    tmp_min = minval(tmp(1:n(1),1:n(2),1:n(3)))
    call mpi_allreduce(MPI_IN_PLACE,tmp_min,1,mpi_real8,mpi_min,comm_cart,ierr)
    tmp_max = maxval(tmp(1:n(1),1:n(2),1:n(3)))
    call mpi_allreduce(MPI_IN_PLACE,tmp_max,1,mpi_real8,mpi_max,comm_cart,ierr)
    !
    if(myid.eq.0) print*, "tmp_min", tmp_min
    if(myid.eq.0) print*, "tmp_max", tmp_max
    !
    return
  end subroutine marangoni_he
    !
    !
    !
    subroutine rk_c1(f_t1,f_t2,n,dli,nh_d,nh_u,nh_v, &
                      dzc,dzf,diff1,diff2,u_reg,eps_int,u,v,w,ur,vr,wr,c1,c2,src,phi,nor,dc1dtrko,dc2dtrko,id)     
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
    real(rp), intent(in   ), dimension(1-nh_u:,1-nh_u:,1-nh_u:)    :: u,v,w,ur,vr,wr
    real(rp), intent(inout), dimension(1-nh_v:,1-nh_v:,1-nh_v:)    :: c1,c2  ! concentraion at current time level
    real(rp), intent(in   ), dimension(1-nh_v:,1-nh_v:,1-nh_v:)    :: phi,src ! level-set,sai
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
!   if (id.eq.1)then 
!   call interfacial_cflux(n,dli,nh_d,nh_v,dzc,dzf,diff1,diff2,c1,c2,phi,dc1dtrk,id)
!   elseif(id.eq.-1)then 
!   call interfacial_cflux(n,dli,nh_d,nh_v,dzc,dzf,diff1,diff2,c1,c2,phi,dc2dtrk,id)
!   endif
    !
    ! 3. Add advection 
    ! 
    if (id.eq.1)then 
    call moms_c1(n,dli,nh_d,nh_u,nh_v,u,v,w,ur,vr,wr,c1,Phi,nor,dc1dtrk)
    elseif(id.eq.-1)then 
    call moms_c1(n,dli,nh_d,nh_u,nh_v,u,v,w,ur,vr,wr,c2,1.d0-Phi,-1.d0*nor,dc2dtrk)
    !
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
        dc2dtrk(i,j,k) = dc2dtrk(i,j,k) + src(i,j,k)
        enddo
       enddo
      enddo
    endif   
    !
    !
    !
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
#if defined(_CONTACT_LINE)
    int_den = 0._rp
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
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          int_den = int_den + 6._rp*c1(i,j,k)*(1._rp-c1(i,j,k))/(dli(1)*dli(2)*dli(3))
        enddo
      enddo
    enddo
    call mpi_allreduce(MPI_IN_PLACE,int_den,1,MPI_REAL_RP,MPI_SUM,comm_cart,ierr)
    q_t = tot_reg/int_den
    !
#endif
    !
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP SHARED(n,f_t1,f_t2,psi,dpsidtrk,dpsidtrko)
    !
    !
    !
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
                    dzc,dzf,u_reg,eps_int,u,v,w,psi,src,ls,nor,kappa,dpsidtrko)
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
    real(rp), intent(in   ), dimension(1-nh_v:,1-nh_v:,1-nh_v:)    :: ls,src   ! level-set
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
    !
#if defined(_CONTACT_LINE)
    int_den = 0._rp
    tot_reg = 0._rp
    !
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
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          int_den = int_den + 6._rp*psi(i,j,k)*(1._rp-psi(i,j,k))/(dli(1)*dli(2)*dli(3))
        enddo
      enddo
    enddo
    call mpi_allreduce(MPI_IN_PLACE,int_den,1,MPI_REAL_RP,MPI_SUM,comm_cart,ierr)
    q_t = tot_reg/int_den
    !
#endif
    !
    ! 3. Add advection 
    ! 
    call moms_psi(n,dli,nh_d,nh_u,nh_v,u,v,w,psi,dpsidtrk)
    !
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
        dpsidtrk(i,j,k) = dpsidtrk(i,j,k)-src(i,j,k)
        enddo
      enddo
    enddo
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
          psi(i,j,k)       = psi(i,j,k) + f_t1*dpsidtrk(i,j,k) + f_t2*dpsidtrko(i,j,k)
          dpsidtrko(i,j,k) = dpsidtrk(i,j,k)-tot_lag
#else

          psi(i,j,k)       = psi(i,j,k) + f_t1*dpsidtrk(i,j,k) + f_t2*dpsidtrko(i,j,k)
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
  subroutine moms_reg(n,dli,nh_d,nh_u,nh_v,u,v,w,mfx,mfy,mfz,s,dsdt)
    !
    implicit none
    !
    integer , intent(in   ), dimension(3)                       :: n
    real(rp), intent(in   ), dimension(3)                       :: dli
    integer , intent(in   )                                     :: nh_d,nh_u,nh_v
    real(rp), intent(in   ), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w,mfx,mfy,mfz
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
          s_flxp = 0.5_rp*( s(i+1,j,k)+s(i,j,k) )*u(i  ,j,k)  -  mfx(i,j,k)
          s_flxm = 0.5_rp*( s(i-1,j,k)+s(i,j,k) )*u(i-1,j,k)  -  mfx(i-1,j,k)
          !
          ! along y
          !
          s_flyp = 0.5_rp*( s(i,j+1,k)+s(i,j,k) )*v(i,j  ,k)  -  mfy(i,j,k)
          s_flym = 0.5_rp*( s(i,j-1,k)+s(i,j,k) )*v(i,j-1,k)  -  mfy(i,j-1,k)
          !
          ! along z
          !
          s_flzp = 0.5_rp*( s(i,j,k+1)+s(i,j,k) )*w(i,j,k  )  -  mfz(i,j,k)
          s_flzm = 0.5_rp*( s(i,j,k-1)+s(i,j,k) )*w(i,j,k-1)  -  mfz(i,j,k-1)
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
  end subroutine moms_reg

  !
   subroutine moms_c1(n,dli,nh_d,nh_u,nh_v,u,v,w,ur,vr,wr,s,phi,nor,dsdt)
    !
    implicit none
    !
    integer , intent(in ), dimension(3)                       :: n
    real(rp), intent(in ), dimension(3)                       :: dli
    integer , intent(in )                                     :: nh_d,nh_u,nh_v
    real(rp), intent(in ), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w
    real(rp), intent(in ), dimension(1-nh_v:,1-nh_v:,1-nh_v:) :: s      !concentration
    real(rp), intent(in ), dimension(1-nh_v:,1-nh_v:,1-nh_v:)    :: phi !phase field 
    real(rp), intent(in ), dimension(1-nh_v:,1-nh_v:,1-nh_v:,1:) :: nor
    real(rp), intent(in ), dimension(0:,0:,0:)                   :: ur,vr,wr
    real(rp), intent(inout), dimension(      1:,     1:,     1:)    :: dsdt
    !
    integer  :: i,j,k,im,jm,km,ip,jp,kp,imm,ipp,jmm,jpp,kmm,kpp
    real(rp) :: s_flxp,s_flxm,s_flyp,s_flym,s_flzp,s_flzm      ,&
                Urx_p,Urx_m,Ury_p,Ury_m,Urz_p,Urz_m,A, &
                phi_face_xp,phi_face_xm,phi_face_yp,phi_face_ym, &
                phi_face_zp,phi_face_zm, nor_face_xp,nor_face_xm, &
                nor_face_yp,nor_face_ym,nor_face_zp,nor_face_zm, &
                c_face_xp,c_face_xm,c_face_yp,c_face_ym,c_face_zp, &
                c_face_zm,ur_face_xp,ur_face_xm,vr_face_yp,vr_face_ym, &
                wr_face_zp, wr_face_zm
    !    
    real(rp) :: abs_u_ijk,ureg_x,ureg_y,ureg_z,urmax
    !
    !
    do k=1,n(3)
      kp = k+1
      km = k-1
      do j=1,n(2)
        jp = j+1
        jm = j-1
        do i=1,n(1)
          ip  = i+1
          im  = i-1

          ! along x
          !
          ! interpolate relative velocity to face
          !
          ur_face_xp  = (ur(ip,j,k)+ur(i,j,k))*0.5
          ur_face_xm  = (ur(im,j,k)+ur(i,j,k))*0.5
          vr_face_yp  = (vr(i,jp,k)+vr(i,j,k))*0.5
          vr_face_ym  = (vr(i,jm,k)+vr(i,j,k))*0.5
          wr_face_zp  = (wr(i,j,kp)+wr(i,j,k))*0.5
          wr_face_zm  = (wr(i,j,km)+wr(i,j,k))*0.5
          !
          ! interpolate phi to face
          !
          phi_face_xp = (phi(ip,j,k)+phi(i,j,k))*0.5
          phi_face_xm = (phi(im,j,k)+phi(i,j,k))*0.5
          !
          ! interpolate normal to face
          !
          nor_face_xp = (nor(ip,j,k,1)+nor(i,j,k,1))*0.5
          nor_face_xm = (nor(im,j,k,1)+nor(i,j,k,1))*0.5
          ! 
          ! interpolate scalar to face 
          !
          c_face_xp   = (s(ip,j,k)+s(i,j,k))*0.5
          c_face_xm   = (s(im,j,k)+s(i,j,k))*0.5
          !
          ! compute u*c term at the face
          !
          s_flxp      = u(i,j,k)*c_face_xp
          s_flxm      = u(im,j,k)*c_face_xm
          !
          Urx_p       = ur_face_xp*c_face_xp
          Urx_m       = ur_face_xm*c_face_xm
          !
          !
          !
          ! along y
          !
          ! interpolate phi to face
          !
          phi_face_yp = (phi(i,jp,k)+phi(i,j,k))*0.5
          phi_face_ym = (phi(i,jm,k)+phi(i,j,k))*0.5
          !
          ! interpolate normal to face
          !
          nor_face_yp = (nor(i,jp,k,2)+nor(i,j,k,2))*0.5
          nor_face_ym = (nor(i,jm,k,2)+nor(i,j,k,2))*0.5
          ! 
          ! interpolate scalar to face 
          !
          c_face_yp   = (s(i,jp,k)+s(i,j,k))*0.5
          c_face_ym   = (s(i,jm,k)+s(i,j,k))*0.5
          !
          ! compute u*c term at the face
          !
          s_flyp      = v(i,j,k)*c_face_yp
          s_flym      = v(i,jm,k)*c_face_ym
          !
          !
          Ury_p       = vr_face_yp*c_face_yp
          Ury_m       = vr_face_ym*c_face_ym
          !        
          !
          ! along z
          !
          ! interpolate phi to face
          !
          phi_face_zp = (phi(i,j,kp)+phi(i,j,k))*0.5
          phi_face_zm = (phi(i,j,km)+phi(i,j,k))*0.5
          !
          ! interpolate normal to face
          !
          nor_face_zp = (nor(i,j,kp,3)+nor(i,j,k,3))*0.5
          nor_face_zm = (nor(i,j,km,3)+nor(i,j,k,3))*0.5
          ! 
          ! interpolate scalar to face 
          !
          c_face_zp   = (s(i,j,kp)+s(i,j,k))*0.5
          c_face_zm   = (s(i,j,km)+s(i,j,k))*0.5
          !
          ! compute u*c term at the face
          !
          s_flzp      = w(i,j,k)*c_face_zp
          s_flzm      = w(i,j,km)*c_face_zm
          !
          !        
          Urz_p        = wr_face_zp*c_face_zp
          Urz_m        = wr_face_zm*c_face_zm
          !
          ! compute the divergence between i+1/2 and i-1/2
          !
          ! overall balance
          !
#if defined(_TWOD)          
          s_flxp = 0.d0
          s_flxm = 0.d0
          Urx_p = 0.d0
          Urx_m = 0.d0
#endif
          dsdt(i,j,k) = - ( dli(1)*((s_flxp+Urx_p)-(s_flxm+Urx_m)) + &
                            dli(2)*((s_flyp+Ury_p)-(s_flym+Ury_m)) + &
                            dli(3)*((s_flzp+Urz_p)-(s_flzm+Urz_m)) )  + dsdt(i,j,k) 
        enddo
      enddo
    enddo

    return
  end subroutine moms_c1
  !
  !

   subroutine cmpt_rel_vel(n,dli,nh_d,nh_u,nh_v,u,v,w,s,phi,nor,ur,vr,wr)
    !
    !
    use mod_param, only: A_rel 
    !
    implicit none
    !
    integer , intent(in ), dimension(3)                       :: n
    real(rp), intent(in ), dimension(3)                       :: dli
    integer , intent(in )                                     :: nh_d,nh_u,nh_v
    real(rp), intent(in ), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w
    real(rp), intent(in ), dimension(1-nh_v:,1-nh_v:,1-nh_v:) :: s      !concentration
    real(rp), intent(in ), dimension(1-nh_v:,1-nh_v:,1-nh_v:)    :: phi !phase field 
    real(rp), intent(in ), dimension(1-nh_v:,1-nh_v:,1-nh_v:,1:) :: nor
    real(rp), intent(out), dimension(      0:,      0:,      0:) :: ur,vr,wr
    !
    integer  :: i,j,k,im,jm,km,ip,jp,kp,imm,ipp,jmm,jpp,kmm,kpp
    real(rp) :: s_flxp,s_flxm,s_flyp,s_flym,s_flzp,s_flzm      ,&
                Urx_p,Urx_m,Ury_p,Ury_m,Urz_p,Urz_m,c_sat, &
                phi_face_xp,phi_face_xm,phi_face_yp,phi_face_ym, &
                phi_face_zp,phi_face_zm, nor_face_xp,nor_face_xm, &
                nor_face_yp,nor_face_ym,nor_face_zp,nor_face_zm, &
                c_face_xp,c_face_xm,c_face_yp,c_face_ym,c_face_zp, &
                c_face_zm
    !    
    real(rp) :: abs_u_ijk,ureg_x,ureg_y,ureg_z,eps
    !
    eps   = 1e-30
    !
    c_sat = 1.0d0
    !
    !
    !Attraction Coefficient  
    !
    !compute rel_vel*phi at the cell centers
  do k=1,n(3)
   do j=1,n(2)
     do i=1,n(1)
          ur(i,j,k) =  A_rel*(0.5-phi(i,j,k))*abs(phi(i,j,k)-1.d0)/(abs(phi(i,j,k)-1.d0)+eps)* &
                          (c_sat - s(i,j,k))*s(i,j,k)**2*nor(i,j,k,1)*phi(i,j,k)
          !
          vr(i,j,k) =  A_rel*(0.5-phi(i,j,k))*abs(phi(i,j,k)-1.d0)/(abs(phi(i,j,k)-1.d0)+eps)* &
                          (c_sat - s(i,j,k))*s(i,j,k)**2*nor(i,j,k,2)*phi(i,j,k)
          !        
          wr(i,j,k) =  A_rel*(0.5-phi(i,j,k))*abs(phi(i,j,k)-1.d0)/(abs(phi(i,j,k)-1.d0)+eps)* &
                          (c_sat - s(i,j,k))*s(i,j,k)**2*nor(i,j,k,3)*phi(i,j,k)
         enddo
      enddo
    enddo

    return
  end subroutine cmpt_rel_vel
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
    !Inspired by S. Mirjalili, M. Khanwale AND A. Mani (2022)                                                                 
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
    subroutine interfacial_cflux(n,dli,nh_d,nh_v, &
                                 dzc,dzf,diff1,diff2,c1,c2,phi,dsdt,id) 
    !
    !Inspired by A computational model for interfacial heat and mass transfer
    !in two-phase flows using a phase field method (Shahab Mirjalili et.al)(2022)
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
          D_m = diff1*diff2/(k_eq*diff1*(1.d0-phi(i,j,k))+diff2*phi(i,j,k))
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
subroutine inittmp(initmp,n,dl,l,nh_v,eps_int,tmp)
    !
    use mod_common_mpi, only: ijk_start
    !
    use mod_param, only: xc,yc,zc,r,tl0,tg0,nbub
    !
    implicit none
    !
    character(len=3), intent(in )                         :: initmp
    integer         , intent(in ), dimension(3)           :: n
    real(rp)        , intent(in ), dimension(3)           :: dl
    real(rp)        , intent(in ), dimension(3)           :: l
    integer         , intent(in )                         :: nh_v
    real(rp)        , intent(in )                         :: eps_int
    real(rp)        , intent(out), dimension(1-nh_v:,1-nh_v:,1-nh_v:) :: tmp
    !
    real(rp) :: phi
    real(rp) :: d,x,y,z,ry
    real(rp) :: rho1,rho2
    integer  :: i,j,k,p
    !
    select case(initmp)
    case('ins')
      do k=1,n(3)
        z = (k+ijk_start(3)-0.5_rp)*dl(3)
        do j=1,n(2)
          y = (j+ijk_start(2)-0.5_rp)*dl(2)
          do i=1,n(1)
            x = (i+ijk_start(1)-0.5_rp)*dl(1)
             !
             tmp(i,j,k) = 0._rp
             do p=1,nbub
#if defined(_TWOD)
               d = sqrt(                  (y-yc(p))**2._rp+(z-zc(p))**2._rp )
#else
               d = sqrt( (x-xc(p))**2._rp+(y-yc(p))**2._rp+(z-zc(p))**2._rp )
#endif
               tmp(i,j,k) = tmp(i,j,k) + 0.5_rp*(1._rp+tanh( (r(p)-d)/(2._rp*eps_int) ))*tl0
             !  
             enddo
             !
          enddo
        enddo
      enddo

    case default
      call flutas_error('ERROR: invalid name for initial TMP field: check di_acdi.in')
    end select
    !
    return
    !
    end subroutine inittmp
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
      if (phi(i,j,k).lt.0.999d0.and.phi(i,j,k).gt.0.001d0)then
      c2(i,j,k) = 0.1        
      else
      c1(i,j,k) = 0.d0
      endif
      
        enddo
      enddo
    enddo
    c1(:,:,:) = 0.d0
    end select
    !
    return
    !
  end subroutine initc1
    !
    !
    !
    !
    subroutine extended(n,dli,nh_d,nh_v,halo_v,dzc,dzf,phi,nx,ny,nz,q,qext,rdir,delta) !phi is the level-set
    !
    use mod_bound, only: boundp
    use mod_param, only: bcpsi,cbcpsi
    !
    implicit none
    !
    integer, intent(in ), dimension(3)                    :: n
    real(8), intent(in ), dimension(3)                    :: dli
    integer , intent(in )                                 :: nh_d,nh_v
    real(8), intent(in ), dimension(1-nh_d:)                  :: dzc,dzf
    integer         , intent(in   ), dimension(3)         :: halo_v
    real(8), intent(in ), dimension(   1-nh_v:,   1-nh_v:,   1-nh_v:) :: phi
    real(8), intent(in ), dimension(    1-nh_v:,    1-nh_v:,    1-nh_v:) :: nx,ny,nz
    real(8), intent(in ), dimension(1-nh_v:,1-nh_v:,1-nh_v:) :: q
    real(8), intent(out), dimension(    0:,    0:,    0:) :: qext
    real(8), intent(in )                                  :: rdir,delta
    !
    real(8), dimension(0:n(1)+1,0:n(2)+1,0:n(3)+1) :: phi1
    real(8), dimension(3) :: dl
    integer :: i,j,k,im,jm,km,ip,jp,kp,ax,ay,az,s,p
    integer :: ps_step
    real(8) :: num,den,eps_mfx
    !
    real(8), parameter :: eps = 1e-08
    !
    if(delta.eq.0.d0) then
      ps_step = 5
      eps_mfx = 3.0d0
    else
      ps_step = 10
      eps_mfx = 6.0d0
    endif
    !
    dl(:) = dli(:)**(-1.d0)
    !
    ! 0. define the level-set curve:
    !    note: --> phi1, level 1;
    !          --> normal defined unless of an constant offset, so no need to be recomputed.
    !
    do k=0,n(3)+1
      do j=0,n(2)+1
        do i=0,n(1)+1
          !
          phi1(i,j,k) = phi(i,j,k)+rdir*delta*dl(2)
          !
        enddo
      enddo
    enddo
    !
    qext(:,:,:) = 0.d0 ! initialize to zero
    !
    !
    ! 1. first populate the cell cut at least in one direction by
    !    the interface (i.e., the interfacial cells)
    !
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          !
          ax = nint(sign(1.d0,nx  (i,j,k)))
          ay = nint(sign(1.d0,ny  (i,j,k)))
          az = nint(sign(1.d0,nz  (i,j,k)))
          p  = nint(sign(1.d0,phi1(i,j,k)))
          !
          if(phi1(i,j,k)*rdir.gt.0.d0) then
            !
            ! case 1 (cut along all the direction)
            !
            if(     phi1(i,j,k)*phi1(i-p*ax,j,k).lt.0.d0.and. &
                    phi1(i,j,k)*phi1(i,j-p*ay,k).lt.0.d0.and. &
                    phi1(i,j,k)*phi1(i,j,k-p*az).lt.0.d0 ) then
              num = ax*nx(i,j,k)*q(i-p*ax,j,k)*dli(1) + &
                    ay*ny(i,j,k)*q(i,j-p*ay,k)*dli(2) + &
                    az*nz(i,j,k)*q(i,j,k-p*az)*dli(3)
              den = ax*nx(i,j,k)*dli(1) + &
                    ay*ny(i,j,k)*dli(2) + &
                    az*nz(i,j,k)*dli(3)
              qext(i,j,k) = num/(den+eps)
            !
            ! case 2 (cut along only x,y, not z)
            !
            elseif( phi1(i,j,k)*phi1(i-p*ax,j,k).lt.0.d0.and. &
                    phi1(i,j,k)*phi1(i,j-p*ay,k).lt.0.d0.and. &
                    phi1(i,j,k)*phi1(i,j,k-p*az).gt.0.d0 ) then
              num = ax*nx(i,j,k)*q(i-p*ax,j,k)*dli(1) + &
                    ay*ny(i,j,k)*q(i,j-p*ay,k)*dli(2)
              den = ax*nx(i,j,k)*dli(1) + &
                    ay*ny(i,j,k)*dli(2)
              qext(i,j,k) = num/(den+eps)
            !
            ! case 3 (cut along only x,z, not y)
            ! 
            elseif( phi1(i,j,k)*phi1(i-p*ax,j,k).lt.0.d0.and. &
                    phi1(i,j,k)*phi1(i,j-p*ay,k).gt.0.d0.and. &
                    phi1(i,j,k)*phi1(i,j,k-p*az).lt.0.d0 ) then
              num = ax*nx(i,j,k)*q(i-p*ax,j,k)*dli(1) + &
                    az*nz(i,j,k)*q(i,j,k-p*az)*dli(3)
              den = ax*nx(i,j,k)*dli(1) + &
                    az*nz(i,j,k)*dli(3)
              qext(i,j,k) = num/(den+eps)
            !
            ! case 4 (cut along only y,z, not x)
            ! 
            elseif( phi1(i,j,k)*phi1(i-p*ax,j,k).gt.0.d0.and. &
                    phi1(i,j,k)*phi1(i,j-p*ay,k).lt.0.d0.and. &
                    phi1(i,j,k)*phi1(i,j,k-p*az).lt.0.d0 ) then
              num = ay*ny(i,j,k)*q(i,j-p*ay,k)*dli(2) + &
                    az*nz(i,j,k)*q(i,j,k-p*az)*dli(3)
              den = ay*ny(i,j,k)*dli(2) + &
                    az*nz(i,j,k)*dli(3)
              qext(i,j,k) = num/(den+eps)
            !
            ! case 5 (cut along only x, not y and z)
            !
            elseif( phi1(i,j,k)*phi1(i-p*ax,j,k).lt.0.d0.and. &
                    phi1(i,j,k)*phi1(i,j-p*ay,k).gt.0.d0.and. &
                    phi1(i,j,k)*phi1(i,j,k-p*az).gt.0.d0 ) then
              num = ax*nx(i,j,k)*q(i-p*ax,j,k)*dli(1)
              den = ax*nx(i,j,k)*dli(1)
              qext(i,j,k) = num/(den+eps)
            !
            ! case 6 (cut along only y, not x and z)
            !
            elseif( phi1(i,j,k)*phi1(i-p*ax,j,k).gt.0.d0.and. &
                    phi1(i,j,k)*phi1(i,j-p*ay,k).lt.0.d0.and. &
                    phi1(i,j,k)*phi1(i,j,k-p*az).gt.0.d0 ) then
              num = ay*ny(i,j,k)*q(i,j-p*ay,k)*dli(2)
              den = ay*ny(i,j,k)*dli(2)
              qext(i,j,k) = num/(den+eps)
            !
            ! case 7 (cut along only z, not y and z)
            !
            elseif( phi1(i,j,k)*phi1(i-p*ax,j,k).gt.0.d0.and. &
                    phi1(i,j,k)*phi1(i,j-p*ay,k).gt.0.d0.and. &
                    phi1(i,j,k)*phi1(i,j,k-p*az).lt.0.d0 ) then
              num = az*nz(i,j,k)*q(i,j,k-p*az)*dli(3)
              den = az*nz(i,j,k)*dli(3)
              qext(i,j,k) = num/(den+eps)
            !
            ! case 8 (no cut)
            !
            elseif( phi1(i,j,k)*phi1(i-p*ax,j,k).gt.0.d0.and. &
                    phi1(i,j,k)*phi1(i,j-p*ay,k).gt.0.d0.and. &
                    phi1(i,j,k)*phi1(i,j,k-p*az).gt.0.d0 ) then
              num = 0.d0
              den = 1.d0 ! arbitrary value
              qext(i,j,k) = num/(den+eps)
            !
            endif
            !
          else
            !
            ! case 8 (no cut)
            !
            if( phi1(i,j,k)*phi1(i-p*ax,j,k).gt.0.d0.and. &
                phi1(i,j,k)*phi1(i,j-p*ay,k).gt.0.d0.and. &
                phi1(i,j,k)*phi1(i,j,k-p*az).gt.0.d0 ) then
                qext(i,j,k) = 0.d0 !q(i,j,k) !num/den
            else
            !
            ! all the others (with at least one cut)
            !
                qext(i,j,k) = q(i,j,k)
            endif
            !
            !qext(i,j,k) = q(i,j,k)
            !qext(i,j,k) = 0.d0
            !
          endif
          !
        enddo
      enddo
    enddo
    call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,qext)
    !
    ! 2. iterate to fill all the cells at least within 3 grid-cells
    !    away from the interface
    !    Note: --> case 8 of previous list accounted as it is sufficient simply not to update qext;
    !          --> we extend in both directions (-/+ normal) using the interger p;
    !          --> ps_step defined above.
    !
    do s=1,ps_step
      !
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            !
            ax = nint(sign(1.d0,nx  (i,j,k)))
            ay = nint(sign(1.d0,ny  (i,j,k)))
            az = nint(sign(1.d0,nz  (i,j,k)))
            p  = nint(sign(1.d0,phi1(i,j,k))) ! so we account both direction (-/+ along the normal)
            !
            if(abs(phi1(i,j,k))*dli(2).lt.eps_mfx) then
              !
              ! case 1 (available neighbours: x,y,z)
              !
              if(     qext(i,j,k).eq.0.d0.and.qext(i-p*ax,j,k).ne.0.d0.and. &
                                              qext(i,j-p*ay,k).ne.0.d0.and. &
                                              qext(i,j,k-p*az).ne.0.d0      ) then
                num = ax*nx(i,j,k)*qext(i-p*ax,j,k)*dli(1) + &
                      ay*ny(i,j,k)*qext(i,j-p*ay,k)*dli(2) + &
                      az*nz(i,j,k)*qext(i,j,k-p*az)*dli(3)
                den = ax*nx(i,j,k)*dli(1) + &
                      ay*ny(i,j,k)*dli(2) + &
                      az*nz(i,j,k)*dli(3)
                qext(i,j,k) = num/(den+eps)
              !
              ! case 2 (available neighbours: x,y, not z)
              !
              elseif( qext(i,j,k).eq.0.d0.and.qext(i-p*ax,j,k).ne.0.d0.and. &
                                              qext(i,j-p*ay,k).ne.0.d0.and. &
                                              qext(i,j,k-p*az).eq.0.d0      ) then
                num = ax*nx(i,j,k)*qext(i-p*ax,j,k)*dli(1) + &
                      ay*ny(i,j,k)*qext(i,j-p*ay,k)*dli(2)
                den = ax*nx(i,j,k)*dli(1) + &
                      ay*ny(i,j,k)*dli(2)
                qext(i,j,k) = num/(den+eps)
              !
              ! case 3 (available neighbours: x,z, not y)
              !
              elseif( qext(i,j,k).eq.0.d0.and.qext(i-p*ax,j,k).ne.0.d0.and. &
                                              qext(i,j-p*ay,k).eq.0.d0.and. &
                                              qext(i,j,k-p*az).ne.0.d0      ) then
                num = ax*nx(i,j,k)*qext(i-p*ax,j,k)*dli(1) + &
                      az*nz(i,j,k)*qext(i,j,k-p*az)*dli(3)
                den = ax*nx(i,j,k)*dli(1) + &
                      az*nz(i,j,k)*dli(3)
                qext(i,j,k) = num/(den+eps)
              !
              ! case 4 (available neighbours: y,z, not x)
              !
              elseif( qext(i,j,k).eq.0.d0.and.qext(i-p*ax,j,k).eq.0.d0.and. &
                                              qext(i,j-p*ay,k).ne.0.d0.and. &
                                              qext(i,j,k-p*az).ne.0.d0      ) then
                num = ay*ny(i,j,k)*qext(i,j-p*ay,k)*dli(2) + &
                      az*nz(i,j,k)*qext(i,j,k-p*az)*dli(3)
                den = ay*ny(i,j,k)*dli(2) + &
                      az*nz(i,j,k)*dli(3)
                qext(i,j,k) = num/(den+eps)
              !
              ! case 5 (available neighbours: x, not y,z)
              !
              elseif( qext(i,j,k).eq.0.d0.and.qext(i-p*ax,j,k).ne.0.d0.and. &
                                              qext(i,j-p*ay,k).eq.0.d0.and. &
                                              qext(i,j,k-p*az).eq.0.d0      ) then
                num = ax*nx(i,j,k)*qext(i-p*ax,j,k)*dli(1)
                den = ax*nx(i,j,k)*dli(1)
                qext(i,j,k) = num/(den+eps)
              !
              ! case 6 (available neighbours: y, not x,z)
              !
              elseif( qext(i,j,k).eq.0.d0.and.qext(i-p*ax,j,k).eq.0.d0.and. &
                                              qext(i,j-p*ay,k).ne.0.d0.and. &
                                              qext(i,j,k-p*az).eq.0.d0      ) then
                num = ay*ny(i,j,k)*qext(i,j-p*ay,k)*dli(2)
                den = ay*ny(i,j,k)*dli(2)
                qext(i,j,k) = num/(den+eps)
              !
              ! case 7 (available neighbours: z, not x,y)
              !
              elseif( qext(i,j,k).eq.0.d0.and.qext(i-p*ax,j,k).eq.0.d0.and. &
                                              qext(i,j-p*ay,k).eq.0.d0.and. &
                                              qext(i,j,k-p*az).ne.0.d0      ) then
                num = az*nz(i,j,k)*qext(i,j,k-p*az)*dli(3)
                den = az*nz(i,j,k)*dli(3)
                qext(i,j,k) = num/(den+eps)
              !
              endif
              !
            else
              !
              !qext(i,j,k) = q(i,j,k)
              qext(i,j,k) = 0.d0
              !
            endif
            !
            !qext(i,j,k) = max(0.d0,qext(i,j,k))
            !
          enddo
        enddo
      enddo
      call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,qext)
      !
    enddo
    !
    return
  end subroutine extended

    !
end module mod_di_acdi
