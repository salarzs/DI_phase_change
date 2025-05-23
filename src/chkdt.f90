!
! SPDX-License-Identifier: MIT
!
module mod_chkdt
  !
  use mpi
  use mod_param     , only: gacc_x,gacc_y,gacc_z,small,cfl_c,cfl_d
  use mod_common_mpi, only: ierr
  use mod_types
  !@cuf use cudafor
  !
  implicit none
  !
  private
#if defined(_TWO_PHASE)
  public  :: chkdt_tw
#else
  public  :: chkdt_sp
#endif
  !
  contains
  !
#if defined(_TWO_PHASE)
  subroutine chkdt_tw(nx,ny,nz,dxi,dyi,dzi,nh_d,nh_u,dzci,dzfi,u,v,w,dtmax)
    !
    ! computes maximum allowed timestep for two-phase flow
    !  as in Kang, Fedkiw and Liu, 
    !  "A boundary condition Capturing Method for Multiphase Incompressible Flow",
    !  JSC 2000
    !
    use mod_param, only: rho1,rho2,mu1,mu2,sigma
#if defined(_HEAT_TRANSFER)
    use mod_param, only: cond1,cond2,cp1,cp2
#endif
    !
    implicit none
    !
    integer , intent(in )                                     :: nx,ny,nz
    real(rp), intent(in )                                     :: dxi,dyi,dzi
    integer , intent(in )                                     :: nh_d,nh_u
    real(rp), intent(in ), dimension(1-nh_d:)                 :: dzci,dzfi
    real(rp), intent(in ), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w
    real(rp), intent(out)                                     :: dtmax
    !
    real(rp) :: ux,uy,uz,vx,vy,vz,wx,wy,wz
    real(rp) :: dtix,dtiy,dtiz,dtic,dtiv,dtik,dtig,dti,dlmin,dlmini
#if defined(_HEAT_TRANSFER)
    real(rp) :: dtth
#endif
    !@cuf attributes(managed) :: u, v, w, dzci, dzfi
    integer :: i,j,k
    !
    dti = 0._rp
    !
#if defined(_OPENACC)
    !$acc parallel loop collapse(3) reduction(max:dti)
#else
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(n,u,v,w,dxi,dyi,dzi,dzci,dzfi) &
    !$OMP PRIVATE(i,j,k,ux,uy,uz,vx,vy,vz,wx,wy,wz,dtix,dtiy,dtiz) &
    !$OMP REDUCTION(max:dti)
#endif
    do k=1,nz
      do j=1,ny
        do i=1,nx
          !
          ux = abs(u(i,j,k))
          vx = 0.25_rp*abs( v(i,j,k)+v(i,j-1,k)+v(i+1,j,k)+v(i+1,j-1,k) )
          wx = 0.25_rp*abs( w(i,j,k)+w(i,j,k-1)+w(i+1,j,k)+w(i+1,j,k-1) )
          dtix = ux*dxi+vx*dyi+wx*dzfi(k)
          !
          uy = 0.25_rp*abs( u(i,j,k)+u(i,j+1,k)+u(i-1,j+1,k)+u(i-1,j,k) )
          vy = abs(v(i,j,k))
          wy = 0.25_rp*abs( w(i,j,k)+w(i,j+1,k)+w(i,j+1,k-1)+w(i,j,k-1) )
          dtiy = uy*dxi+vy*dyi+wy*dzfi(k)
          !
          uz = 0.25_rp*abs( u(i,j,k)+u(i-1,j,k)+u(i-1,j,k+1)+u(i,j,k+1) )
          vz = 0.25_rp*abs( v(i,j,k)+v(i,j-1,k)+v(i,j-1,k+1)+v(i,j,k+1) )
          wz = abs(w(i,j,k))
          dtiz = uz*dxi+vz*dyi+wz*dzci(k)
          !
          dti = max(dti,dtix,dtiy,dtiz)
          !
        enddo
      enddo
    enddo
#if defined(_OPENACC)
    !$acc end parallel loop 
#else
    !$OMP END PARALLEL DO
#endif
    call mpi_allreduce(MPI_IN_PLACE,dti,1,MPI_REAL_RP,MPI_MAX,MPI_COMM_WORLD,ierr)
    if(dti.eq.0._rp) dti = 1._rp
    dtic   = dti
    !
    dlmin  = min(1._rp/dxi,1._rp/dyi,1._rp/dzi)
    dlmin  = min(dlmin,minval(1._rp/dzfi(:))) ! minimum of dzf is an estimate on the safe side
    dlmini = dlmin**(-1)
    dtiv   = max(mu1/rho1,mu2/rho2)*dlmini**2/cfl_d
    if(sigma.eq.0._rp) then
      dtik = small 
    else
      dtik = sqrt(sigma/(min(rho1,rho2))*dlmini**3)
    endif
    dtig   = sqrt( max(abs(gacc_x),abs(gacc_y),abs(gacc_z))*dlmini )
    dtmax  = cfl_c*2._rp*(dtic+dtiv+sqrt((dtic+dtiv)**2+4._rp*(dtig**2+dtik**2)))**(-1)
    dtmax  = min(dtmax,1._rp/dtik) ! to ensure a satisfaction of the capillary time-step restriction if dtmax>dtk
    !
#if defined(_HEAT_TRANSFER)
    dtth   = cfl_d*(max(cond1/(rho1*cp1),cond2/(rho2*cp2))*dlmini**2)**(-1)
    dtmax  = min(dtmax,dtth)
#endif
    !
    return
  end subroutine chkdt_tw
#else
  subroutine chkdt_sp(nx,ny,nz,dxi,dyi,dzi,nh_d,nh_u,dzci,dzfi,u,v,w,dtmax)
    !
    ! computes maximum allowed timestep for single-phase flow 
    ! (optionally with heat transfer)
    ! 
    use mod_param, only: rho_sp,mu_sp
#if defined(_HEAT_TRANSFER)
    use mod_param, only: cond1,cond2,cp1,cp2!kappa_sp,cp_sp
#endif
    !
    implicit none
    !
    integer , intent(in )                                     :: nx,ny,nz
    real(rp), intent(in )                                     :: dxi,dyi,dzi
    integer , intent(in )                                     :: nh_d,nh_u
    real(rp), intent(in ), dimension(1-nh_d:)                 :: dzci,dzfi
    real(rp), intent(in ), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w
    real(rp), intent(out)                                     :: dtmax
    !
    real(rp) :: ux,uy,uz,vx,vy,vz,wx,wy,wz
    real(rp) :: dtix,dtiy,dtiz,dtic,dtiv,dtig,dti,dlmin,dlmini
#if defined(_HEAT_TRANSFER)
    real(rp) :: dtth
#endif
    !@cuf attributes(managed) :: u, v, w, dzci, dzfi
    integer :: i,j,k
    !
    dti = 0._rp
    !
#if defined(_OPENACC)
    !$acc parallel loop collapse(3) reduction(max:dti)
#else
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(n,u,v,w,dxi,dyi,dzi,dzci,dzfi) &
    !$OMP PRIVATE(i,j,k,ux,uy,uz,vx,vy,vz,wx,wy,wz,dtix,dtiy,dtiz) &
    !$OMP REDUCTION(max:dti)
#endif
    do k=1,nz
      do j=1,ny
        do i=1,nx
          !
          ux = abs(u(i,j,k))
          vx = 0.25_rp*abs( v(i,j,k)+v(i,j-1,k)+v(i+1,j,k)+v(i+1,j-1,k) )
          wx = 0.25_rp*abs( w(i,j,k)+w(i,j,k-1)+w(i+1,j,k)+w(i+1,j,k-1) )
          dtix = ux*dxi+vx*dyi+wx*dzfi(k)
          !
          uy = 0.25_rp*abs( u(i,j,k)+u(i,j+1,k)+u(i-1,j+1,k)+u(i-1,j,k) )
          vy = abs(v(i,j,k))
          wy = 0.25_rp*abs( w(i,j,k)+w(i,j+1,k)+w(i,j+1,k-1)+w(i,j,k-1) )
          dtiy = uy*dxi+vy*dyi+wy*dzfi(k)
          !
          uz = 0.25_rp*abs( u(i,j,k)+u(i-1,j,k)+u(i-1,j,k+1)+u(i,j,k+1) )
          vz = 0.25_rp*abs( v(i,j,k)+v(i,j-1,k)+v(i,j-1,k+1)+v(i,j,k+1) )
          wz = abs(w(i,j,k))
          dtiz = uz*dxi+vz*dyi+wz*dzci(k)
          !
          dti = max(dti,dtix,dtiy,dtiz)
          !
        enddo
      enddo
    enddo
#if defined(_OPENACC)
    !$acc end parallel loop 
#else
    !$OMP END PARALLEL DO
#endif
    call mpi_allreduce(MPI_IN_PLACE,dti,1,MPI_REAL_RP,MPI_MAX,MPI_COMM_WORLD,ierr)
    if(dti.eq.0._rp) dti = 1._rp
    dtic   = dti
    !
    dlmin  = min(1._rp/dxi,1._rp/dyi,1._rp/dzi)
    dlmin  = min(dlmin,minval(1._rp/dzfi(:))) ! minimum of dzf is an estimate on the safe side
    dlmini = dlmin**(-1)
    dtiv   = (mu_sp/rho_sp)*dlmini**2
    dtig   = sqrt( max(abs(gacc_x),abs(gacc_y),abs(gacc_z))*dlmini )
    dtmax  = min(cfl_c/dti,cfl_d/dtiv,1._rp/dtig)
#if defined(_HEAT_TRANSFER)
    dtth   = cfl_d*(max(cond1/(rho1*cp1),cond2/(rho2*cp2))*dlmini**2)**(-1)
    !dtth   = cfl_d*((max(cond1,cond2)/(rho_sp*min(cp1,cp2)))*dlmini**2)**(-1)
    dtmax  = min(dtmax,dtth)
#endif
    !
    return
  end subroutine chkdt_sp
#endif
  !
end module mod_chkdt
