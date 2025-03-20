! SPDX-License-Identifier: MIT
!
module mod_phase_change
  !

  !
  use mod_types
  use mod_common_mpi, only: myid,ierr,ijk_start,comm_cart
  use mod_sanity    , only: flutas_error
  use mod_bound     , only: boundp,updt_rhs_b, &
                            boundsb_array, & 
                            bounduvw, boundp_di
  use mod_di_acdi
  !
  implicit none
  !
  private
  public ::  src_evap, update_property_var,rk_m2, cmpt_rhs_incompressible,cmpt_rhs_low_mach,initvap

  contains
    !
    !
  !
  subroutine rk_m2(f_t1,f_t2,n,dli,nh_d,nh_u,nh_v, &
                    dzc,dzf,u_reg,eps_int,u,v,w,m2,psi2,src,ls,norm,kappa,dm2dtrko,rho_gas)
!    !
!    
!    !adam-bashforth scheme for time integration of m2 field
!    !
    implicit none
!    !
    real(rp), intent(in   )                                        :: f_t1,f_t2
    integer , intent(in   ), dimension(3)                          :: n
    real(rp), intent(in   ), dimension(3)                          :: dli
    integer , intent(in   )                                        :: nh_d,nh_u,nh_v
    real(rp), intent(in   ), dimension(1-nh_d:)                    :: dzc,dzf
    real(rp), intent(in   )                                        :: u_reg,eps_int
    real(rp), intent(in   ), dimension(1-nh_u:,1-nh_u:,1-nh_u:)    :: u,v,w,rho_gas, psi2
    real(rp), intent(inout), dimension(1-nh_v:,1-nh_v:,1-nh_v:)    :: m2  ! volume-fraction at current time level
    real(rp), intent(in   ), dimension(1-nh_v:,1-nh_v:,1-nh_v:)    :: ls,src   ! level-set
    real(rp), intent(in   ), dimension(1-nh_v:,1-nh_v:,1-nh_v:,1:) :: norm
    real(rp), intent(in ), dimension( 0:, 0:, 0:)               :: kappa
    real(rp), intent(inout), dimension(     1:,     1:,     1:)    :: dm2dtrko
!
!    ! 
    real(rp), dimension(n(1),n(2),n(3)) :: dm2dtrk
!    
    integer  :: i,j,k
!    !
!    ! 1. compute regularization  
!!   
!! 
#if defined(_REGU) 
    !print*, 'regu on m2'
    call regularization_m2_cdi(n,dli,nh_d,nh_v,dzc,dzf,u_reg,eps_int,psi2,norm,dm2dtrk,rho_gas)

#else
   ! print*, 'regu off m2'
    dm2dtrk(:,:,:) = 0.0

#endif
!!
!!    
!!    ! 2. Add advection 
!!    ! 
    call moms_m2(n,dli,nh_d,nh_u,nh_v,u,v,w,m2,dm2dtrk)
!!    !
!!    ! 3. Add source term
     do k=1,n(3)
      do j=1,n(2)
       do i=1,n(1)
        dm2dtrk(i,j,k) = dm2dtrk(i,j,k)-src(i,j,k)
        enddo
      enddo
    enddo
!!    !
!!    ! 4. Advance the color function
!!    !
!!    !$OMP PARALLEL DO DEFAULT(none) &
!!    !$OMP PRIVATE(i,j,k) &
!!    !$OMP SHARED(n,f_t1,f_t2,psi,dpsidtrk,dpsidtrko)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
!!   
          m2(i,j,k)       = m2(i,j,k) + f_t1*dm2dtrk(i,j,k) + f_t2*dm2dtrko(i,j,k)
          dm2dtrko(i,j,k) = dm2dtrk(i,j,k)
!!
!!          !
        enddo
      enddo
    enddo
!!    !$OMP END PARALLEL DO
!!    !
    return
  end subroutine rk_m2
!!  !
!!
!!  !!!!!!!!!!!!!!Advection of m2 field
!!
!!
  subroutine moms_m2(n,dli,nh_d,nh_u,nh_v,u,v,w,s,dsdt)
!!    !
    implicit none
!!    !
    integer , intent(in   ), dimension(3)                       :: n
    real(rp), intent(in   ), dimension(3)                       :: dli
    integer , intent(in   )                                     :: nh_d,nh_u,nh_v
    real(rp), intent(in   ), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w
    real(rp), intent(in   ), dimension(1-nh_v:,1-nh_v:,1-nh_v:) :: s
    real(rp), intent(inout), dimension(     1:,     1:,     1:) :: dsdt
!!    !
    real(rp) :: s_flxp,s_flxm,s_flyp,s_flym,s_flzp,s_flzm
    integer  :: i,j,k
!!    !
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
!!          !
!!          ! along x
!!          !
          s_flxp = 0.5_rp*( s(i+1,j,k)+s(i,j,k) )*u(i  ,j,k)
          s_flxm = 0.5_rp*( s(i-1,j,k)+s(i,j,k) )*u(i-1,j,k)
!!          !
!!          ! along y
!!          !
          s_flyp = 0.5_rp*( s(i,j+1,k)+s(i,j,k) )*v(i,j  ,k)
          s_flym = 0.5_rp*( s(i,j-1,k)+s(i,j,k) )*v(i,j-1,k)
!!          !
!!          ! along z
!!          !
          s_flzp = 0.5_rp*( s(i,j,k+1)+s(i,j,k) )*w(i,j,k  )
          s_flzm = 0.5_rp*( s(i,j,k-1)+s(i,j,k) )*w(i,j,k-1)

          ! overall balance
          !
          dsdt(i,j,k) = - ( dli(1)*(s_flxp-s_flxm) + &
                            dli(2)*(s_flyp-s_flym) + &
                            dli(3)*(s_flzp-s_flzm) ) + dsdt(i,j,k)

      enddo
    enddo
  enddo
   return
  end subroutine moms_m2
!!  !
!!    !  
!!  !!!!!!!!!!!!!!!!!!!!!!!1regularization part of the m2 field
!!
!!
    subroutine regularization_m2_cdi(n,dli,nh_d,nh_v, &
                                     dzc,dzf,u_reg,eps_int,s,nor,dsdt,rho_gas)  !s is volume fraction/phi
!!    !
    implicit none
!!    !
    integer , intent(in ), dimension(3)               :: n
    real(rp), intent(in ), dimension(3)               :: dli 
    integer , intent(in )                             :: nh_d,nh_v
    real(rp), intent(in ), dimension(1-nh_d:)         :: dzc,dzf
    real(rp), intent(in )                             :: u_reg,eps_int
    real(rp), intent(in ), dimension(1-nh_v:,1-nh_v:,1-nh_v:)    :: s,rho_gas
    real(rp), intent(in ), dimension(1-nh_v:,1-nh_v:,1-nh_v:,1:) :: nor
    real(rp), intent(out), dimension(1:,1:,1:)                   :: dsdt
!!    
!!    !
    real(rp), dimension(3) :: dl
    real(rp) :: s_face,term_p,term_m,nor_face,s_min,s_max,rho_p,rho_m 
    real(rp) :: gradpxp,gradpxm,gradpyp,gradpym,gradpzp,gradpzm
    real(rp) :: psifxp,psifxm,psifyp,psifym,psifzp,psifzm
    real(rp) :: regxp,regxm,regyp,regym,regzp,regzm
    real(rp) :: rho_xp,rho_xm,rho_yp,rho_ym,rho_zp,rho_zm
    real(rp) :: abs_u_ijk
    real(rp) :: Psi_face
    integer  :: i,j,k,ip,im,jp,jm,kp,km
!!    !
    dsdt(:,:,:) = 0.d0
!!    !
!!    !
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
!!    !
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
!
      psifxp         =   ( s(i,j,k)*(1.d0-s(i,j,k))*nor(i,j,k,1) + s(ip,j,k)*(1.d0-s(ip,j,k))*nor(ip,j,k,1))*0.5 !phi*(1-phi)*nor at i+0.5 as p.15 Shahab.M
      psifxm         =   ( s(i,j,k)*(1.d0-s(i,j,k))*nor(i,j,k,1) + s(im,j,k)*(1.d0-s(im,j,k))*nor(im,j,k,1))*0.5 !phi*(1-phi)*nor at i-0.5 as p.15 Shahab.M
!
!       
      rho_xp         =   (rho_gas(i+1,j,k)+rho_gas(i,j,k))*0.5
      rho_xm         =   (rho_gas(i,j,k)+rho_gas(i-1,j,k))*0.5
!!      
      regxp          =   rho_xp*u_reg*( eps_int*gradpxp - psifxp)        !regularization at i+0.5
      regxm          =   rho_xm*u_reg*( eps_int*gradpxm - psifxm)        !regularization at i-0.5
!!      
!!   
#endif
!!
      gradpyp        =   ( s(i,jp,k) - s(i,j,k)  ) * dli(2)       !grad(phi) at j+0.5
      gradpym        =   ( s(i,j,k)  - s(i,jm,k) ) * dli(2)       !grad(phi) at j-0.5
!!
      psifyp         =   ( s(i,j,k)*(1.d0-s(i,j,k))*nor(i,j,k,2) + s(i,jp,k)*(1.d0-s(i,jp,k))*nor(i,jp,k,2))*0.5 !phi*(1-phi)*nor at j+0.5 as p.15 Shahab.M
      psifym         =   ( s(i,j,k)*(1.d0-s(i,j,k))*nor(i,j,k,2) + s(i,jm,k)*(1.d0-s(i,jm,k))*nor(i,jm,k,2))*0.5 !phi*(1-phi)*nor at j-0.5 as p.15 Shahab.M
!!
      rho_yp         =   (rho_gas(i,j+1,k)+rho_gas(i,j,k))*0.5
      rho_ym         =   (rho_gas(i,j,k)+rho_gas(i,j-1,k))*0.5
!!      
!!    
      regyp          =   rho_yp*u_reg*( eps_int*gradpyp - psifyp)        !regularization at j+0.5
      regym          =   rho_ym*u_reg*( eps_int*gradpym - psifym)        !regularization at j-0.5
!!      
      gradpzp        =   ( s(i,j,kp) - s(i,j,k)  ) * dli(3)       !grad(phi) at k+0.5
      gradpzm        =   ( s(i,j,k)  - s(i,j,km) ) * dli(3)       !grad(phi) at k-0.5
!!
      psifzp         =   ( s(i,j,k)*(1.d0-s(i,j,k))*nor(i,j,k,3) + s(i,j,kp)*(1.d0-s(i,j,kp))*nor(i,j,kp,3))*0.5 !phi*(1-phi)*nor at k+0.5 as p.15 Shahab.M
      psifzm         =   ( s(i,j,k)*(1.d0-s(i,j,k))*nor(i,j,k,3) + s(i,j,km)*(1.d0-s(i,j,km))*nor(i,j,km,3))*0.5 !phi*(1-phi)*nor at k-0.5 as p.15 Shahab.M
!!      
      rho_zp         =   (rho_gas(i,j,k+1)+rho_gas(i,j,k))*0.5
      rho_zm         =   (rho_gas(i,j,k)+rho_gas(i,j,k-1))*0.5
!!      
!!      
      regzp          =   rho_zp*u_reg*( eps_int*gradpzp - psifzp)        !regularization at k+0.5
      regzm          =   rho_zm*u_reg*( eps_int*gradpzm - psifzm)        !regularization at k-0.5
!!      
      dsdt(i,j,k)    =    (regxp - regxm) * dli(1) + &
                          (regyp - regym) * dli(2) + &
                          (regzp - regzm) * dli(3) 
!!
!!     
!!         !
        enddo
      enddo
    enddo
!!    !
   return
  end subroutine regularization_m2_cdi
!  !
 
  subroutine src_evap(n,nh_v,tmp,src_final,src_c2,t_relax,psi,omega,rho_gas,m2,tmp_min,tmp_max,p0,c2,o_sat)  
    !
    implicit none
    
    integer , intent(in   ), dimension(3)                          :: n
    integer , intent(in   )                                        :: nh_v
    real(rp), intent(in   ), dimension(1-nh_v:,1-nh_v:,1-nh_v:)    ::psi,omega,rho_gas,tmp,m2,c2,o_sat 
    real(rp), intent(inout), dimension(1-nh_v:,1-nh_v:,1-nh_v:)    :: src_c2,src_final  !source term
    real(rp), intent(in)     :: t_relax,tmp_min,tmp_max,p0  !source term
    real(rp),  dimension(1:n(1),1:n(2),1:n(3))    :: src_mod,src
    real(rp) :: p_sat,omega_sat,grid_tot,src_tot,src_grid,src_mod_tot,temp
    integer  :: i,j,k
    !
    grid_tot = 0
    src_mod_tot = 0

    do i=1,n(1)
      do j=1,n(2)
        do k=1,n(3)
            !Anotine equation of state, coefficients are for water between 273 k -373 K
            !p_sat = (10**(8.07131-(1730.63/(233.426+tmp_min-273.0))))*133.322
            !p_sat = 5630.0
           ! omega_sat =  0.1!p_sat/(p_sat+(p0-p_sat)*(28.0/18.0))
           
             if(psi(i,j,k) .lt. 1.0 .and. psi(i,j,k) .gt. 0.0)then  
                      
                  !  temp = ((psi(i,j,k))*(1.0-psi(i,j,k))*10.0*rho_gas(i,j,k))/&
                  !          ((10.0*psi(i,j,k))+m2(i,j,k))
                     
                    temp = 6.0*psi(i,j,k)*(1.0-psi(i,j,k))

!                     temp = 30.0*psi(i,j,k)*psi(i,j,k)*(1.0-psi(i,j,k))*(1.0-psi(i,j,k))
                     
                   !  src(i,j,k)=t_relax*6.0*psi(i,j,k)*(1-psi(i,j,k))*&
                    !          rho_gas(i,j,k)*( c2(i,j,k) - (1.0-psi(i,j,k))*o_sat(i,j,k) )
                    src(i,j,k)=t_relax*temp*&
                              rho_gas(i,j,k)*( c2(i,j,k) - (1.0-psi(i,j,k))*o_sat(i,j,k) )
                    
                    src_c2(i,j,k) = t_relax*temp*&
                                  ( c2(i,j,k) - (1.0-psi(i,j,k))*o_sat(i,j,k) )

                 !   src(i,j,k)=t_relax*6.0*psi(i,j,k)*(1-psi(i,j,k))*&
                 !             m2(i,j,k)* ( omega(i,j,k) - o_sat(i,j,k) )

              else
                    src(i,j,k) = 0.0
                    src_c2(i,j,k) = 0.0
            endif

              if(src(i,j,k) .ne. 0) then 
                grid_tot = grid_tot + 1
              endif
    !          src(i,j,k) = 0.0

                
    
        enddo
      enddo
    enddo
    
    call mpi_allreduce(MPI_IN_PLACE,grid_tot ,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)

    src_tot  =  sum(src(1:n(1),1:n(2),1:n(3)))

    call mpi_allreduce(MPI_IN_PLACE, src_tot ,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)

    src_grid = src_tot/grid_tot

     do i=1,n(1)
      do j=1,n(2)
        do k=1,n(3)

         if(psi(i,j,k) .lt. 0.999 .and. psi(i,j,k) .gt. 0.001) then
           !  if(abs(src(i,j,k)) .gt. 1.0)then

             src_mod (i,j,k) =   6.0*psi(i,j,k)*(1.0-psi(i,j,k))*src_grid
          else
              
             src_mod(i,j,k) = 0.0
          endif
        
       enddo
      enddo
    enddo

    src_mod_tot  =  sum(src_mod(1:n(1),1:n(2),1:n(3)))
    call mpi_allreduce(MPI_IN_PLACE, src_mod_tot ,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)    

    do i=1,n(1)
      do j=1,n(2)
        do k=1,n(3)

          src_final(i,j,k) = min(src(i,j,k),0.0)!src_mod(i,j,k)*(src_tot/src_mod_tot)

         ! src_c2(i,j,k) = t_relax*&
          !               ( c2(i,j,k) - (1.0-psi(i,j,k))*o_sat(i,j,k) )

          src_c2(i,j,k) = min(src_c2(i,j,k),0.0)

        
       enddo
      enddo
    enddo
        


    return
  end subroutine src_evap


!  subroutine update_property_gas_mix(n,prop_vap,prop_air,prop_liq,omega,psi,psi2,prop_gas,prop)
!    !
!    implicit none
    !
!   integer , intent(in ), dimension(3)        :: n
!    real(rp), intent(in ), dimension(0:,0:,0:) :: omega,psi,psi2
!    real(rp), intent(out), dimension(0:,0:,0:) :: prop_gas,prop
!    real(rp), intent(in)  ::prop_vap,prop_air,prop_liq
!    integer  :: i,j,k
    
    ! 
   !$acc parallel loop collapse(3)
    
!    do k=1,n(3)
!      do j=1,n(2)
!         do i=1,n(1)
      
         ! prop_gas(i,j,k) = 1.0/((1.0-omega(i,j,k))/prop_air +&
                               ! omega(i,j,k)/prop_vap)
!          prop_gas(i,j,k) = omega(i,j,k)*prop_vap + &
!                           (1.0 - omega(i,j,k))*prop_air                     
!          prop(i,j,k) = psi(i,j,k)*prop_liq + &
!                     psi2(i,j,k)*prop_gas(i,j,k)
        
!        enddo
!      enddo
!    enddo
    !$acc end parallel loop 
    !
!    return
!  end subroutine update_property_gas_mix

 subroutine update_property_var(n,psi,prop_const,prop_var,prop)
    !
    implicit none
    !
    integer , intent(in ), dimension(3)        :: n
    real(rp), intent(in ), dimension(0:,0:,0:) :: prop_var,psi
    real(rp), intent(out), dimension(0:,0:,0:) :: prop
    real(rp), intent(in)  ::prop_const
    integer  :: i,j,k
    
    ! 
   !$acc parallel loop collapse(3)
    
    do k=1,n(3)
      do j=1,n(2)
         do i=1,n(1)
      
                           
          prop(i,j,k) = psi(i,j,k)*prop_const + &
                     (1.d0-psi(i,j,k))*prop_var(i,j,k)
        
        enddo
      enddo
    enddo
    !$acc end parallel loop 
    !
    return
  end subroutine update_property_var

  subroutine cmpt_rhs_incompressible(n,nh_v,dli,rho1,rho_air,rho_vap,omeg,m2,c2,diff2,src,psi,rhs_div_th)
    !
    implicit none
    
    integer , intent(in   ), dimension(3)                          :: n
    real(rp), intent(in ), dimension(3)               :: dli 
    integer , intent(in   )                                        :: nh_v
    real(rp), intent(in   ), dimension(1-nh_v:,1-nh_v:,1-nh_v:)    :: omeg,m2,src,psi,c2
    real(rp), intent(inout), dimension(1-nh_v:,1-nh_v:,1-nh_v:)    :: rhs_div_th  !source term
    real(rp), intent(in)     :: rho1,diff2,rho_air,rho_vap
  ! real(rp),  dimension(0:n(1)+1,0:n(2)+1,0:n(3)+1)    :: omeg

    real(rp) :: diffxp,diffxm,diffyp,diffym,diffzp,diffzm,vap_diff
    real(rp) :: gradxp,gradxm,gradyp,gradym,gradzp,gradzm
    integer  :: i,j,k,ip,im,jp,jm,kp,km
    !
    

  do i=1,n(1)
      do j=1,n(2)
        do k=1,n(3)

    
          ip = i+1
          im = i-1
          jp = j+1
          jm = j-1
          kp = k+1
          km = k-1


#if defined(_TWOD)
          diffxp = 0.d0
          diffxm = 0.d0
#else          
          gradxp = ( omeg(ip,j,k) - omeg(i,j,k)   )*dli(1) !gradient of omega computed at face i+1/2
          gradxm = ( omeg(i,j,k)  - omeg(im,j,k)  )*dli(1) !gradient of omega computed at face i-1/2

          diffxp = gradxp*(m2(i,j,k)+m2(i+1,j,k))*0.5
          diffxp = gradxp*(m2(i,j,k)+m2(i+1,j,k))*0.5
          diffxm = gradxm*(m2(i,j,k)+m2(i-1,j,k))*0.5                                                                                                      
#endif
        
          gradyp = ( omeg(i,jp,k) - omeg(i,j,k)   )*dli(2) !gradient of omega computed at face j+1/2
          gradym = ( omeg(i,j,k)  - omeg(i,jm,k)  )*dli(2) !gradient of omega computed at face j-1/2

          diffyp = gradyp*(m2(i,j,k)+m2(i,j+1,k))*0.5
          diffym = gradym*(m2(i,j,k)+m2(i,j-1,k))*0.5      


          gradzp = ( omeg(i,j,kp) - omeg(i,j,k)   )*dli(3) !gradient of omega computed at face k+1/2
          gradzm = ( omeg(i,j,k)  - omeg(i,j,km)  )*dli(3) !gradient of omega computed at face k-1/2

          diffzp = gradzp*(m2(i,j,k)+m2(i,j,k+1))*0.5
          diffzm = gradzm*(m2(i,j,k)+m2(i,j,k-1))*0.5 

          vap_diff   =  ((diffxp - diffxm) * dli(1) + &
                         (diffyp - diffym) * dli(2) + &
                         (diffzp - diffzm) * dli(3))*diff2  
          
          
          
          rhs_div_th(i,j,k) = src(i,j,k)*(1.d0/rho1 - 1.d0/rho_vap) + &
                            (1.0/rho_vap-1.0/rho_air)*vap_diff*0.0
       enddo
      enddo
    enddo
   
    return
  end subroutine cmpt_rhs_incompressible

  subroutine cmpt_rhs_low_mach(n,nh_v,dli,diff2,rho_gas,rho1,rho_air,rho_vap,cpp,cp1,cp_gas,cp_air,cp_vap,&
                              cv,cv1,cv_gas,cv_air,cv_vap,mol_mas,mol_mas_air,mol_mas_vap,& 
                              psi2,omega,m2,src,cond,tmp,u,v,w,dt,rhs_div_th,dp0dt,p0,h01,h02)
 
    implicit none
    
    integer , intent(in   ), dimension(3)                          :: n
    real(rp), intent(in ), dimension(3)               :: dli 
    integer , intent(in   )                                        :: nh_v
    real(rp), intent(in   ), dimension(1-nh_v:,1-nh_v:,1-nh_v:)    :: omega,m2,src,psi2,cond,tmp
    real(rp), intent(in   ), dimension(1-nh_v:,1-nh_v:,1-nh_v:)    :: cpp,cv,cp_gas,cv_gas,rho_gas, mol_mas
    real(rp), intent(in   ), dimension(1-nh_v:,1-nh_v:,1-nh_v:)    :: u, v, w
    real(rp), intent(inout), dimension(1-nh_v:n(1)+nh_v,1-nh_v:n(2)+nh_v,1-nh_v:n(3)+nh_v)    :: rhs_div_th  !source term
    real(rp), intent(in)     :: rho1,diff2,rho_air,rho_vap
    real(rp), intent(in)     :: cp1,cp_air,cp_vap,cv1,cv_air,cv_vap,mol_mas_air,mol_mas_vap, dt, h01, h02
    real(rp), intent(inout)     :: p0,dp0dt

    real(rp), dimension(0:n(1)+1,0:n(2)+1,0:n(3)+1)    ::  conduct_flux
    real(rp), dimension(1:n(1),1:n(2),1:n(3))    :: z, d_caps, vap_diff, divu, ks, b_hat
   
    real(rp) :: diffxp,diffxm,diffyp,diffym,diffzp,diffzm
    real(rp) :: gradxp,gradxm,gradyp,gradym,gradzp,gradzm
    real(rp) ::  gamma2, tot_dp0dt, tot_ks, h1, hv, ha
    real(rp) :: term1, term2, term3, term4, term5, term1_diff, term2_diff, term_new
    real(rp) :: sfxp, sfxm, sfyp, sfym, sfzp, sfzm 
    integer  :: i,j,k,ip,im,jp,jm,kp,km
    !


    tot_dp0dt = 0.0
    tot_ks = 0.0
   ! conduct_flux(:,:,:) = 0.0
    
    do i=1,n(1)
      do j=1,n(2)
        do k=1,n(3)
!         

          gamma2 = 1.4         

          ip = i+1
          im = i-1
          jp = j+1
          jm = j-1
          kp = k+1
          km = k-1

          
#if defined(_TWOD)
          diffxp = 0.d0
          diffxm = 0.d0
#else          
          gradxp = ( omega(ip,j,k) - omega(i,j,k)   )*dli(1) !gradient of omega computed at face i+1/2
          gradxm = ( omega(i,j,k)  - omega(im,j,k)  )*dli(1) !gradient of omega computed at face i-1/2

          diffxp =  gradxp*(m2(i,j,k)+m2(i+1,j,k)  )*0.5
          diffxm =  gradxm*(m2(i,j,k)+m2(i-1,j,k)  )*0.5

#endif
          gradyp = ( omega(i,jp,k) - omega(i,j,k)   )*dli(2) !gradient of omega computed at face j+1/2
          gradym = ( omega(i,j,k)  - omega(i,jm,k)  )*dli(2) !gradient of omega computed at face j-1/2


          diffyp = gradyp*(m2(i,j,k)+m2(i,j+1,k))*0.5
          diffym = gradym*(m2(i,j,k)+m2(i,j-1,k))*0.5      
!
!
          gradzp = ( omega(i,j,kp) - omega(i,j,k)   )*dli(3) !gradient of omega computed at face k+1/2
          gradzm = ( omega(i,j,k)  - omega(i,j,km)  )*dli(3) !gradient of omega computed at face k-1/2
!
          diffzp = gradzp*(m2(i,j,k)+m2(i,j,k+1))*0.5
          diffzm = gradzm*(m2(i,j,k)+m2(i,j,k-1))*0.5 
!
          vap_diff(i,j,k)   =  ((diffxp - diffxm) * dli(1) + &
                               (diffyp - diffym) * dli(2) + &
                               (diffzp - diffzm) * dli(3))*diff2  
                           
          
          term_new = tmp(i,j,k)*&
                     ((1.0-psi2(i,j,k))*rho1*cp1 + m2(i,j,k)*cp_gas(i,j,k))

          term1 = (1.d0/rho1 - 1.d0/rho_gas(i,j,k))

          term2 = (1.0-omega(i,j,k))*(mol_mas(i,j,k)/rho_gas(i,j,k))*(1.0/mol_mas_vap - 1.0/mol_mas_air)

          h1 = cp1*tmp(i,j,k)+ h01

          hv =  cp_vap*tmp(i,j,k)+h02

          ha =  cp_air*tmp(i,j,k)+h02
          
          term3 = (psi2(i,j,k)*(hv - h1))/term_new

          z(i,j,k) = term1 - term2 + term3

          term1_diff = (mol_mas(i,j,k)/rho_gas(i,j,k))*(1.0/mol_mas_vap - 1.0/mol_mas_air)

          term2_diff = (psi2(i,j,k)*(ha - hv))/term_new

         !if(psi2(i,j,k) .gt. 0.1)then

                d_caps(i,j,k) =  (term1_diff + term2_diff)
        !  else
        !        d_caps(i,j,k) = 0.0
        !  endif 
         
#if defined(_HEAT_TRANSFER)

           b_hat (i,j,k) = psi2(i,j,k)/term_new
         
          ! along x-dir
          !
          sfxp = (cond(i+1,j,k)+cond(i,j,k))*(tmp(i+1,j,k)-tmp(i,j,k))*dli(1)*0.5
          sfxm = (cond(i-1,j,k)+cond(i,j,k))*(tmp(i,j,k)-tmp(i-1,j,k))*dli(1)*0.5
          
          ! along y-dir
          !
          sfyp = (cond(i,j+1,k)+cond(i,j,k))*(tmp(i,j+1,k)-tmp(i,j,k))*dli(2)*0.5
          sfym = (cond(i,j-1,k)+cond(i,j,k))*(tmp(i,j,k)-tmp(i,j-1,k))*dli(2)*0.5
           
          ! along z-dir
          !
          sfzp = (cond(i,j,k+1)+cond(i,j,k))*(tmp(i,j,k+1)-tmp(i,j,k))*dli(3)*0.5
          sfzm = (cond(i,j,k-1)+cond(i,j,k))*(tmp(i,j,k)-tmp(i,j,k-1))*dli(3)*0.5
          !

          conduct_flux(i,j,k)  = (sfxp-sfxm)*dli(1)  + (sfyp-sfym)*dli(2) + (sfzp-sfzm)*dli(3)
#else
          b_hat(i,j,k) = 0.0
          conduct_flux(i,j,k) = 0.0
#endif

          ks(i,j,k) = psi2(i,j,k) * (1.0/p0 - 1.0/term_new)

       enddo
      enddo
    enddo


     

     do k=1,n(3)
      do j=1,n(2)
         do i=1,n(1)
      
           tot_dp0dt = tot_dp0dt + &
                            (z(i,j,k)*src(i,j,k) + d_caps(i,j,k)*vap_diff(i,j,k)+&
                             b_hat(i,j,k)*conduct_flux(i,j,k) )/&
                            (dli(1)*dli(2)*dli(3))
                            

           tot_ks = tot_ks + ks(i,j,k)/(dli(1)*dli(2)*dli(3))

        enddo
      enddo
    enddo


    call mpi_allreduce(MPI_IN_PLACE,tot_dp0dt,1,MPI_REAL_RP,MPI_SUM,comm_cart,ierr)
    call mpi_allreduce(MPI_IN_PLACE,tot_ks,1,MPI_REAL_RP,MPI_SUM,comm_cart,ierr)
    
    dp0dt = 0.0!tot_dp0dt/tot_ks
  

      do k=1,n(3)
       do j=1,n(2)
         do i=1,n(1)

         ! rhs_div_th(i,j,k) = conduct_flux(i,j,k)
         ! rhs_div_th(i,j,k) = cond(i,j,k)
         ! rhs_div_th(i,j,k) = tmp(i,j,k)

          rhs_div_th(i,j,k) = z(i,j,k)*src(i,j,k)*1.0 +  d_caps(i,j,k)*vap_diff(i,j,k) +&
                               b_hat(i,j,k)*conduct_flux(i,j,k) - &
                               1.0*ks(i,j,k)*dp0dt
 

        enddo
      enddo
    enddo


    return
  end subroutine cmpt_rhs_low_mach
!


subroutine initvap(n,dli,nh_d,nh_u,nh_v,halo_v,dzc,dzf,diff,rho,phi,nor,c2,omega,t_relax,o_sat)
!
use mod_param, only: cbcpsi,bcpsi,cbcsca,bcsca
!
! computes initial conditions for the scalar field
!
implicit none
!
integer , intent(in   ), dimension(3)                          :: n
real(rp), intent(in   ), dimension(3)                          :: dli
integer , intent(in   )                                        :: nh_d,nh_u,nh_v
integer , intent(in   ), dimension(3)                          :: halo_v
real(rp), intent(in   ), dimension(1-nh_d:)                    :: dzc,dzf
real(rp), intent(in   )                                        :: diff
real(rp), intent(in   )                                        :: rho,t_relax
real(rp), intent(in   ), dimension(1-nh_v:,1-nh_v:,1-nh_v:)    :: phi,o_sat
real(rp), intent(in   ), dimension(1-nh_v:,1-nh_v:,1-nh_v:,1:) :: nor
real(rp), intent(out)  , dimension(1-nh_v:,1-nh_v:,1-nh_v:)    :: c2,omega
!
real(rp), dimension(1:n(1),1:n(2),1:n(3))                               :: dc1dtrko,dc2dtrko
real(rp), dimension(1-nh_v:n(1)+nh_v,1-nh_v:n(2)+nh_v,1-nh_v:n(3)+nh_v) :: u,v,w,ur,vr,wr,src
real(rp), dimension(1-nh_v:n(1)+nh_v,1-nh_v:n(2)+nh_v,1-nh_v:n(3)+nh_v) :: omega_s,c1
real(rp), dimension(3) :: dl
real(rp)               :: dt
integer, parameter    :: nstep_sca = 10000
integer               :: i,j,k,istep
!
dl(:) = dli(:)**(-1.d0)
!
omega(:,:,:)    = 0.d0
c2(:,:,:)       = 0.d0
!
!
!select case(inisca)
!
!case('std')
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

      omega_s(i,j,k) = 0.1!mass_fraction(pth,tmp(i,j,k))

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
!call cmpt_src(n,dli,nh_d,nh_u,nh_v,dzc,dzf,phi,rho,omega,omega_s,src)

!call src_evap(n,nh_v,tmp,src,t_relax,psi,omega,rho_gas,m2,tmp_min,tmp_max,p0,c2)
!call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,src)
!

do i = 1,n(1)
  do j = 1,n(2)
    do k = 1,n(3)
        src(i,j,k) = t_relax*6.0*phi(i,j,k)*(1.0-phi(i,j,k))*&
                     (1.0-phi(i,j,k))*rho*(omega(i,j,k)-o_sat(i,j,k))
    enddo
  enddo
enddo
call boundp(cbcsca,n,bcsca,nh_d,nh_v,halo_v,dl,dzc,dzf,src)

! subroutine rk_c1(f_t1,f_t2,n,dli,nh_d,nh_u,nh_v, &
!                      dzc,dzf,diff1,diff2,diff,u_reg,eps_int,u,v,w,ur,vr,wr,c1,c2,src,phi,nor,dc1dtrko,dc2dtrko,id)

call rk_c1(dt,dt,n,dli,nh_d,nh_u,nh_v,dzc,dzf,0.d0,diff,0.d0,0.d0,u*0.d0,v*0.d0,w*0.d0, &
           ur*0.d0,vr*0.d0,wr*0.d0,c1,c2,src,phi,nor,dc1dtrko,dc2dtrko,-1)
!
call boundp(cbcsca,n,bcsca,nh_d,nh_v,halo_v,dl,dzc,dzf,c2)

!call boundp_c(n,nh_d,nh_v,halo_v,dl,dzc,dzf,phi,rho,omega,c2)
!
!
do i = 1,n(1)
  do j = 1,n(2)
    do k = 1,n(3)
        omega(i,j,k) = c2(i,j,k)/(rho*(1.0-phi(i,j,k))+1e-30)
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
!end select
!
return
end subroutine initvap

end module mod_phase_change
