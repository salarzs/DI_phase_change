module mod_phase_change
!
use mod_types
use mod_common_mpi, only: myid,ierr,ijk_start,comm_cart
use mod_bound, only:boundp,boundp_c,boundp_e
use mod_di_acdi, only:rk_c1
!
implicit none
!
private 
public :: cmpt_src,initvap,mass_fraction 
!
contains
!
!Antoine relation to compute 
!
real(8) function mass_fraction(pth,tmp)
!
use mod_param, only: m1,m2
!
real(8), intent(in) :: pth,tmp
!
real(8), parameter :: a = +8.07131, &  ! Water
                      b = +1730.63, &  ! Water    
                      c = +233.426, &  ! Water
                      d = +133.322     ! Water
!
real(8), parameter :: thr = 0.98d0
real(8), parameter :: eps = 1.0e-09
!
real(8)            :: pvap
!
pvap            = (10**(a-(b/(c+tmp-273.0))))*d
mass_fraction   = pvap/(pvap+(pth-pvap)*(m2/m1))
mass_fraction   = max(eps,min(mass_fraction,thr))
!
!
return
end function mass_fraction
!
!
subroutine initvap(inisca,n,dli,nh_d,nh_u,nh_v,halo_v,dzc,dzf,diff,rho,pth,tmp,phi,nor,c2,omega,src)
!
use mod_param, only: cbcpsi,bcpsi,cbcsca,bcsca,t_relax_i
!
! rho is the gas_density 
! computes initial conditions for the scalar field
!
implicit none
!
character(len=3), intent(in )                                  :: inisca
integer , intent(in   ), dimension(3)                          :: n
real(rp), intent(in   ), dimension(3)                          :: dli
integer , intent(in   )                                        :: nh_d,nh_u,nh_v
integer , intent(in   ), dimension(3)                          :: halo_v
real(rp), intent(in   ), dimension(1-nh_d:)                    :: dzc,dzf
real(rp), intent(in   )                                        :: diff
real(rp), intent(in   )                                        :: pth,rho
real(rp), intent(in   ), dimension(1-nh_v:,1-nh_v:,1-nh_v:)    :: phi,tmp
real(rp), intent(in   ), dimension(1-nh_v:,1-nh_v:,1-nh_v:,1:) :: nor
real(rp), intent(out)  , dimension(1-nh_v:,1-nh_v:,1-nh_v:)    :: c2,omega,src  
!
real(rp), dimension(1:n(1),1:n(2),1:n(3))                               :: dc1dtrko,dc2dtrko
real(rp), dimension(1-nh_v:n(1)+nh_v,1-nh_v:n(2)+nh_v,1-nh_v:n(3)+nh_v) :: u,v,w,ur,vr,wr 
real(rp), dimension(1-nh_v:n(1)+nh_v,1-nh_v:n(2)+nh_v,1-nh_v:n(3)+nh_v) :: omega_s,c1
real(rp), dimension(3) :: dl
real(rp)               :: dt
integer, parameter    :: nstep_sca = 10000
real(rp), parameter   :: small = real(2.e-16,rp)
integer               :: i,j,k,istep
!
dl(:) = dli(:)**(-1.d0)
!
omega(:,:,:)    = 0.d0
c2(:,:,:)       = 0.d0
!
!
select case(inisca)
!
case('std')
!
! b. calculation of the steady-state solution
!
if(myid.eq.0) print*, 'calculation of the steady state solution'
!
!
dt = 0.10d0/(min(dli(1)**2,dli(2)**2,dli(3)**2)*diff)
!
do i=1,n(1)
  do j=1,n(2)
    do k=1,n(3)

      omega_s(i,j,k) = mass_fraction(pth,tmp(i,j,k))

    enddo
  enddo
enddo
!
call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,omega_s)
!
do istep=1,nstep_sca
!
if(myid.eq.0) print*, 'Step = ', istep
!
!Necessary to be done this way
!
do i=1,n(1)
  do j=1,n(2)
    do k=1,n(3)

      src(i,j,k)=(6.0*(1.d0 - phi(i,j,k)) * phi(i,j,k)*((1.d0 - phi(i,j,k))*&
                  rho*(omega_s(i,j,k)-omega(i,j,k))))/(t_relax_i)

    enddo
  enddo
enddo
!
!call cmpt_src(n,dli,nh_d,nh_u,nh_v,dzc,dzf,phi,rho,omega,omega_s,src)
call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,src)
!
call rk_c1(dt,dt,n,dli,nh_d,nh_u,nh_v,dzc,dzf,0.d0,diff,0.d0,0.d0,u*0.d0,v*0.d0,w*0.d0, &
           ur*0.d0,vr*0.d0,wr*0.d0,c1,c2,src,phi,nor,dc1dtrko,dc2dtrko,-1)
!   
call boundp_c(n,nh_d,nh_v,halo_v,dl,dzc,dzf,phi,rho,omega,c2)
!
!
do i = 1,n(1)
  do j = 1,n(2)
    do k = 1,n(3)
        omega(i,j,k) = c2(i,j,k)/(rho*(1.0-phi(i,j,k))+small)
    enddo
  enddo
enddo
!
call boundp(cbcsca,n,bcsca,nh_d,nh_v,halo_v,dl,dzc,dzf,omega)
!
dc1dtrko(:,:,:) = 0.d0
dc2dtrko(:,:,:) = 0.d0
!
enddo
!
end select
!
return
end subroutine initvap
!
subroutine cmpt_src(n,dli,nh_d,nh_u,nh_v,dzc,dzf,phi,rho,t_relax,omega,omega_s,src)
!
!
!rho here is the gas density
!
! compute the interfacial mass-transfer rate 
!
implicit none
!
integer , intent(in   ), dimension(3)                          :: n
real(rp), intent(in   ), dimension(3)                          :: dli
integer , intent(in   )                                        :: nh_d,nh_u,nh_v
real(rp), intent(in   )                                        :: rho,t_relax
real(rp), intent(in   ), dimension(1-nh_d:)                    :: dzc,dzf
real(rp), intent(in   ), dimension(1-nh_v:,1-nh_v:,1-nh_v:)    :: phi
real(rp), intent(in   ), dimension(1-nh_v:,1-nh_v:,1-nh_v:)    :: omega,omega_s
real(rp), intent(out  ), dimension(1-nh_v:,1-nh_v:,1-nh_v:)    :: src
integer  :: i,j,k 
!
do i=1,n(1)
  do j=1,n(2)
    do k=1,n(3)

      src(i,j,k)=(6.0*(1.d0 - phi(i,j,k)) * phi(i,j,k)*((1.d0 - phi(i,j,k))*&
                  rho*(omega_s(i,j,k)-omega(i,j,k))))/(t_relax)

    enddo
  enddo
enddo
!
return
end subroutine cmpt_src
!
!
!
end module mod_phase_change
