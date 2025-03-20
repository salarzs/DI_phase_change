!
! SPDX-License-Identifier: MIT
!
!---------------------------------------------------------------------------------
! FluTAS -- Fluid Transport Accelerated Solver                            
!                                                                                 
!  a. Menu Title: FluTAS_two_phase_di;                                                 
!  b. Feautures of FluTAS_two_phase_di:                                                
!      --> two-fluid incompressible and adiabatic solver with Diffuse Interface;
!      --> allow for density and viscosity mismatch;                              
!      --> momentum equation advanced with Adams-Bashforth (explicit diffusion);  
!      --> pressure equation solved with a FFT-based direct solver.               
!---------------------------------------------------------------------------------
!
program flutas
  !
  ! module declaration 
  !  note: --> import what you really neeed 
  !
  use iso_c_binding , only: C_PTR
  use mpi
  use decomp_2d
  use mod_bound     , only: boundp,updt_rhs_b, &
                            boundsb_array, & 
                            bounduvw, boundp_di

  use mod_chkdiv    , only: chkdiv
#if defined(_TWO_PHASE)  
  use mod_chkdt     , only: chkdt_tw
#else
  use mod_chkdt     , only: chkdt_sp
#endif
  use mod_common_mpi, only: myid,ierr,comm_cart,n_z,ijk_start,ipencil
  use mod_correc    , only: correc
  use mod_debug     , only: cmpt_mean
  use mod_di_acdi   , only: rk_psi,initls,initpsi,psi_to_ls,cmpt_norm,ls_to_psi,cmpt_umax, &
                            static_contact_angle_m, inittmp,&
                            update_property 
  use mod_di_acdi   , only: rk_c1, rk_ent,initc1,cmpt_rel_vel 
  use mod_di_acdi   , only: cmpt_delta
  use mod_phase_change, only: src_evap,update_property_var,rk_m2,cmpt_rhs_incompressible, &
                              cmpt_rhs_low_mach,initvap


  use mod_fft       , only: fftini,fftend
  use mod_fillps    , only: fillps
  use mod_initflow  , only: initflow
  use mod_initgrid  , only: initgrid
  use mod_initmpi   , only: initmpi,halo_gen
#if defined(_OPENACC)
  use mod_initmpi   , only: alloc_buf
#endif
  use mod_initsolver, only: initsolver
  use mod_load      , only: load, load_scalar
  use mod_rk        , only: rk, cmpt_time_factors
  use mod_output    , only: out0d,out1d,out2d,out3d,write_visu_2d,write_visu_3d
  use mod_param     , only: lx,ly,lz,dx,dy,dz,dxi,dyi,dzi,small,is_wallturb, &
                            cbcvel,bcvel,cbcpre,bcpre, &
                            icheck,iout0d,iout1d,iout2d,iout3d,isave, &
                            nstep,time_max,tw_max,stop_type,restart, &
                            restart_sph,restart_mph,restart_mph_scalar,restart_mph_heat,cfl,    &
                            constant_dt,dt_input, &
                            inivel, &
                            itot,jtot,ktot,dims_in, &
                            nthreadsmax, &
                            gr, &
                            is_outflow,no_outflow,is_forced, &
                            rho1,rho2,rho0,mu1,mu2,cbcpsi,bcpsi,late_init,i_late_init, &
                            rho0,cbctmp,cbcent,bctmp,bcent,tl0,tg0,cp1,cp2,cond1,cond2, &
                            h01,h02,rho_air,rho_vap, t_relax,cv1,cv2,    &
                            cp_air,cp_vap,cv_air,cv_vap,cond_air,cond_vap, &
                            mol_mas_air, mol_mas_vap, &
                            sigma, sigma_t,theta, & 



                            cbcsca, bcsca,cbcrel,bcrel ,diff1,diff2,  &
                            n,ng,l,dl,dli, &
                            bulk_ftype,rkcoeff, &
                            time_scheme,space_scheme_mom,n_stage, &
                            gamma_v,eps_int,inipsi, initmp,&
                            cbcnor,bcnor, &
                            read_input
  !
#if defined(_DO_POSTPROC) && defined(_TWO_PHASE)
  use mod_param     , only: do_tagging,iout0d_ta
#endif
#if defined(_DO_POSTPROC) && defined(_TURB_FORCING)
  use mod_post      , only: budget
#endif
  use mod_sanity    , only: test_sanity
  use mod_source    , only: bulk_forcing_src
  use mod_source    , only: surft_src,grav_tw_src,pres_tw_src
#if defined(_TURB_FORCING)
  use mod_source    , only: forc_src 
#endif
  use mod_source    , only: mar_src
#if defined(_OPENACC)
  use mod_solver_gpu, only: solver_gpu
#else
  use mod_solver_cpu, only: solver_cpu
#endif
  use mod_types
  use profiler
#if defined(_DO_POSTPROC)
  use mod_tagging   , only: droplet_tagging
#endif
  !@cuf use mod_common_mpi, only: mydev
  !@cuf use cudafor
  !
  !$ use omp_lib
  !
  implicit none
  !
  ! Variables declaration
  !  note: --> first declare arrays, then the other variables;
  !        --> order of declaration: type, dimension, allocatable;
  !
  real(rp), dimension(:,:,:), allocatable :: u,v,w,div,p,mfx,mfy,mfz,uold,vold,wold
  real(rp), dimension(:,:,:), allocatable :: pold
  real(rp), dimension(:,:,:), allocatable :: mu,rho,muold,rhoold,rho_new
  real(rp), dimension(:,:,:), pointer :: dudtrko,dvdtrko,dwdtrko
  real(rp), dimension(:), allocatable :: dzc  ,dzf  ,zc  ,zf  ,dzci  ,dzfi
  real(rp), dimension(:), allocatable :: dzc_g,dzf_g,zc_g,zf_g,dzci_g,dzfi_g
  !
  real(rp), dimension(:,:,:)  , allocatable :: psi,kappa,psiold,kappaold,m2,psi2
  real(rp), dimension(:,:,:,:), allocatable :: nor,norold                   
  real(rp), dimension(:,:,:)  , allocatable :: ls,dpsidtrko,dm2dtrko,rhs_div_th,omeg
  real(rp), dimension(:,:,:)  , allocatable :: c1,c2,dc1dtrko,dc2dtrko,src_ph,omega_s,diff,omega,src_mod
  real(rp), dimension(:,:,:)  , allocatable :: ent,rhoh,tmp,cpp,cond,rhocp,dentdtrko,cv,conduction
  real(rp), dimension(:,:,:)  , allocatable :: rho_gas,cp_gas,cv_gas,rhocp_gas,cond_gas,mol_mas
  real(rp), dimension(:,:,:)  , allocatable :: cppold,condold,rhocpold,lh,o_sat
  real(rp), dimension(:,:,:)  , allocatable :: ur,vr,wr
  real(rp), dimension(:,:,:)  , allocatable :: regarx,regary,regarz,source
  !
  real(rp), dimension(:,:,:)  , allocatable :: delta

real(rp), dimension(:,:)  , allocatable :: psi_bound,c1_bound
  !
  type(C_PTR), dimension(2,2) :: arrplanp
  real(rp), allocatable, dimension(:,:) :: lambdaxyp
  real(rp), allocatable, dimension(:)   :: ap,bp,cp
  real(rp) :: normfftp
  ! 
  real(rp), allocatable, dimension(:,:,:) :: rhsbp_x, rhsbp_y, rhsbp_z
  !
  real(rp), dimension(3) :: f
  real(rp) :: dt,dto,dti,dtmax,time,dtrk,dtrki,divtot,divmax, &
              f_t1,f_t2,f_t12,f_t12_o,f_t12_i
#if defined(_TIMING)
  real(rp) :: dt12,dt12av,dt12min,dt12max
#endif
  !
  integer, dimension(3)   :: halo_u,halo_p,halo_d,halo_v
  integer, dimension(3)   :: dims
  integer, dimension(3,3) :: dims_xyz
  integer  :: nh_d,nh_u,nh_p,nh_v
  !
  integer  :: i,j,k,im,ip,jm,jp,km,kp
  integer  :: irk,istep
  character(len=100) :: datadir,datadir_ta,restart_dir
  character(len=1)   :: action_load
  logical  :: is_data
  !
  real(rp) :: meanvel,meanvelu,meanvelv,meanvelw
  real(rp), dimension(3) :: dpdl_c
  real(rp), dimension(10) :: var
  ! 
  character(len=9) :: fldnum
  real(rp) :: f1d,f2d,f3d
  real(rp) :: twi,tw,phi_wall
  integer :: kk
  logical :: is_done,kill
  real(rp) :: rho0i,u_reg,tmp_min,tmp_max, p0, dp0dt, vel_rise, v_sum
  real(rp) :: rhoxp,rhoyp,rhozp,rhox,rhoy,rhoz
  !
  !@cuf integer :: istat
  !@cuf integer(kind=cuda_count_kind) :: freeMem, totMem
  !@cuf attributes(managed) :: pold, kappa, mu, rho, psi, d_thinc, cur_t, nor
  !@cuf attributes(managed) :: u, v, w, p 
  !@cuf attributes(managed) :: dzc  , dzf  , zc  , zf  , dzci, dzfi
  !@cuf attributes(managed) :: zc_g, zf_g
  !@cuf attributes(managed) :: lambdaxyp, ap, bp, cp, rhsbp_x, rhsbp_y, rhsbp_z
  !@cuf attributes(managed) :: dudtrko, dvdtrko, dwdtrko
  !
  !if we don't use dropcheck.f90 we can comment the next line  
  real(rp) :: xd,yd,zd,ut,vt,wt,zcd,ycd,xcd,vol
  real(rp) :: psi_min,psi_max,vol_p1,vol_c1,vol_c2,vol_ct,ke_tot,omega_max,omega_min
  real(rp) ::  m_a, m_v,h1,h2, mass_gas, vol_gas, sum_rt, tmp_ini, p_sat
  real(rp) ::  sum_tmp, sum_grid, avg_tmp

  real(rp) :: gradpx, psifx, norfx, regx, &
              gradpy, psify, norfy, regy, &
              gradpz, psifz, norfz, regz, &
              tmpx, tmpy, tmpz

  real(rp) :: aux
  integer  :: ii,jj
  integer :: n1, n2, n3
  integer :: ng1, ng2, ng3
  !
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
  !
  call profiler_init()
  !
  ! read parameter file
  !
  call read_input(myid)
  !
  ! create data folder and subfolders for post-processing, if they do not exist.
  !
  inquire(file='data/',exist=is_data)
  if(.not.is_data.and.myid.eq.0) call execute_command_line('mkdir -p data')
  datadir = 'data/'
  !
  inquire(file='data/restart_dir/',exist=is_data)
  if(.not.is_data.and.myid.eq.0) call execute_command_line('mkdir -p data/restart_dir')
  restart_dir = 'data/restart_dir/'
  !
#if defined(_DO_POSTPROC)
  inquire(file='data/post/',exist=is_data)
  if(.not.is_data.and.myid.eq.0) call execute_command_line('mkdir -p data/post')
  inquire(file='data/post/tagging/',exist=is_data)
  if(.not.is_data.and.myid.eq.0.and.do_tagging) call execute_command_line('mkdir -p data/post/tagging')
  datadir_ta = 'data/post/tagging/'
#endif
  !
  ! initialize MPI/OpenMP
  !
  !$call omp_set_num_threads(nthreadsmax)
  call initmpi(.false.,ng,cbcpre,dims_in,dims_xyz,dims,n)
  !
  n1 = n(1)
  n2 = n(2)
  n3 = n(3)
  ng1 = ng(1)
  ng2 = ng(2)
  ng3 = ng(3)
  !
  twi = MPI_WTIME()
  !
  ! halo calculation
  !
  nh_u = 1
  nh_p = 1
  nh_v = 1
  !
  nh_d = max(nh_u,nh_p,nh_v) ! take the maximum of the previous ones
  !
  ! allocate memory
  !
  allocate(p(0:n(1)+1,0:n(2)+1,0:n(3)+1) , &
           u(1-nh_u:n(1)+nh_u,1-nh_u:n(2)+nh_u,1-nh_u:n(3)+nh_u) , &
           uold(1-nh_u:n(1)+nh_u,1-nh_u:n(2)+nh_u,1-nh_u:n(3)+nh_u) , &
           mfx(1-nh_u:n(1)+nh_u,1-nh_u:n(2)+nh_u,1-nh_u:n(3)+nh_u) , &
           regarx(1:n(1),1:n(2),1:n(3)) , &
           v(1-nh_u:n(1)+nh_u,1-nh_u:n(2)+nh_u,1-nh_u:n(3)+nh_u) , &
           vold(1-nh_u:n(1)+nh_u,1-nh_u:n(2)+nh_u,1-nh_u:n(3)+nh_u) , &
           mfy(1-nh_u:n(1)+nh_u,1-nh_u:n(2)+nh_u,1-nh_u:n(3)+nh_u) , &
           regary(1:n(1),1:n(2),1:n(3)) , &
           w(1-nh_u:n(1)+nh_u,1-nh_u:n(2)+nh_u,1-nh_u:n(3)+nh_u) , &
           wold(1-nh_u:n(1)+nh_u,1-nh_u:n(2)+nh_u,1-nh_u:n(3)+nh_u) , &
           mfz(1-nh_u:n(1)+nh_u,1-nh_u:n(2)+nh_u,1-nh_u:n(3)+nh_u) , &
           regarz(1:n(1),1:n(2),1:n(3)) , &
           source(1:n(1),1:n(2),1:n(3)) , &
           div(1:n(1),1:n(2),1:n(3)) , &
           pold(0:n(1)+1,0:n(2)+1,0:n(3)+1))
  allocate(dudtrko(n(1),n(2),n(3)), &
           dvdtrko(n(1),n(2),n(3)), &
           dwdtrko(n(1),n(2),n(3)))
  allocate(psi(      0:n(1)+1,0:n(2)+1,0:n(3)+1)  , &
           m2(      0:n(1)+1,0:n(2)+1,0:n(3)+1)  , &
           psi2(      0:n(1)+1,0:n(2)+1,0:n(3)+1)  , &
           psiold(       0:n(1)+1,0:n(2)+1,0:n(3)+1)  , &
           rhs_div_th(      0:n(1)+1,0:n(2)+1,0:n(3)+1)  , &
           kappa(    0:n(1)+1,0:n(2)+1,0:n(3)+1)  , &
           kappaold(    0:n(1)+1,0:n(2)+1,0:n(3)+1)  , &
           mu(       0:n(1)+1,0:n(2)+1,0:n(3)+1)  , &
           muold(       0:n(1)+1,0:n(2)+1,0:n(3)+1)  , &
           rhoold(      0:n(1)+1,0:n(2)+1,0:n(3)+1)  , &
           rho_new(      0:n(1)+1,0:n(2)+1,0:n(3)+1)  , &
           rho(      0:n(1)+1,0:n(2)+1,0:n(3)+1)  , &
           rho_gas(      0:n(1)+1,0:n(2)+1,0:n(3)+1)  , &
           cp_gas(      0:n(1)+1,0:n(2)+1,0:n(3)+1)  , &
           cv_gas(      0:n(1)+1,0:n(2)+1,0:n(3)+1)  , &
           rhocp_gas(      0:n(1)+1,0:n(2)+1,0:n(3)+1)  , &
           cond_gas(      0:n(1)+1,0:n(2)+1,0:n(3)+1)  , &
           mol_mas(      0:n(1)+1,0:n(2)+1,0:n(3)+1)  , &
           conduction(      0:n(1)+1,0:n(2)+1,0:n(3)+1)  , &
           nor(      0:n(1)+1,0:n(2)+1,0:n(3)+1,3), &
           norold(      0:n(1)+1,0:n(2)+1,0:n(3)+1,3), &
           ls(       0:n(1)+1,0:n(2)+1,0:n(3)+1)  , &
           dm2dtrko(0:n(1)+1,0:n(2)+1,0:n(3)+1), &
           dpsidtrko(0:n(1)+1,0:n(2)+1,0:n(3)+1))
  allocate(c1(0:n(1)+1,0:n(2)+1,0:n(3)+1)        , &
           c2(0:n(1)+1,0:n(2)+1,0:n(3)+1)        , &
           omega(0:n(1)+1,0:n(2)+1,0:n(3)+1)        , &
           o_sat(0:n(1)+1,0:n(2)+1,0:n(3)+1)        , &
           omeg(0:n(1)+1,0:n(2)+1,0:n(3)+1)        , &
           diff(0:n(1)+1,0:n(2)+1,0:n(3)+1)        , &
           src_ph(0:n(1)+1,0:n(2)+1,0:n(3)+1)        , &
           src_mod(0:n(1)+1,0:n(2)+1,0:n(3)+1)        , &
           ent(0:n(1)+1,0:n(2)+1,0:n(3)+1)        , &
           lh(0:n(1)+1,0:n(2)+1,0:n(3)+1)          ,&
           rhoh(0:n(1)+1,0:n(2)+1,0:n(3)+1)        , &
           tmp(0:n(1)+1,0:n(2)+1,0:n(3)+1)        , &
           cond(0:n(1)+1,0:n(2)+1,0:n(3)+1)        , &
           condold(0:n(1)+1,0:n(2)+1,0:n(3)+1)        , &
           cpp(0:n(1)+1,0:n(2)+1,0:n(3)+1)        , &
           cv(0:n(1)+1,0:n(2)+1,0:n(3)+1)        , &
           rhocp(0:n(1)+1,0:n(2)+1,0:n(3)+1)        , &
           rhocpold(0:n(1)+1,0:n(2)+1,0:n(3)+1)        , &
           cppold(0:n(1)+1,0:n(2)+1,0:n(3)+1)        , &
           dentdtrko(0:n(1)+1,0:n(2)+1,0:n(3)+1)  , &
           dc1dtrko(0:n(1)+1,0:n(2)+1,0:n(3)+1)  , &
           dc2dtrko(0:n(1)+1,0:n(2)+1,0:n(3)+1)  , &
           ur(0:n(1)+1,0:n(2)+1,0:n(3)+1), &
           vr(0:n(1)+1,0:n(2)+1,0:n(3)+1), &
           wr(0:n(1)+1,0:n(2)+1,0:n(3)+1) ) 
  allocate(delta(  0:n(1)+1,0:n(2)+1,0:n(3)+1))
  allocate(psi_bound(1:n(1),1:n(2))) 
  allocate(c1_bound(1:n(1),1:n(2))) 
  allocate(lambdaxyp(n_z(1),n_z(2)))
  allocate(ap(n_z(3)),bp(n_z(3)),cp(n_z(3)))
  allocate(dzc( 1-nh_d:n(3)+nh_d), &
           dzf( 1-nh_d:n(3)+nh_d), &
           zc(  1-nh_d:n(3)+nh_d), &
           zf(  1-nh_d:n(3)+nh_d), &
           dzci(1-nh_d:n(3)+nh_d), &
           dzfi(1-nh_d:n(3)+nh_d))
  allocate(dzc_g( 1-nh_d:ng(3)+nh_d), &
           dzf_g( 1-nh_d:ng(3)+nh_d), &
           zc_g(  1-nh_d:ng(3)+nh_d), &
           zf_g(  1-nh_d:ng(3)+nh_d), &
           dzci_g(1-nh_d:ng(3)+nh_d), &
           dzfi_g(1-nh_d:ng(3)+nh_d))
  allocate(rhsbp_x(n(2),n(3),0:1), &
           rhsbp_y(n(1),n(3),0:1), &
           rhsbp_z(n(1),n(2),0:1))
  !
  ! prefetching of the variables (TODO: remember to add the one of x-pencil!)
  !
  !@cuf istat = cudaMemAdvise(u, size(u), cudaMemAdviseSetPreferredLocation, mydev)
  !@cuf istat = cudaMemAdvise(v, size(v), cudaMemAdviseSetPreferredLocation, mydev)
  !@cuf istat = cudaMemAdvise(w, size(w), cudaMemAdviseSetPreferredLocation, mydev)
  !@cuf istat = cudaMemAdvise(p, size(p), cudaMemAdviseSetPreferredLocation, mydev)
  !@cuf istat = cudaMemAdvise(pold, size(pold), cudaMemAdviseSetPreferredLocation, mydev)
  !@cuf istat = cudaMemAdvise(dudtrko, size(dudtrko), cudaMemAdviseSetPreferredLocation, mydev)
  !@cuf istat = cudaMemAdvise(dvdtrko, size(dvdtrko), cudaMemAdviseSetPreferredLocation, mydev)
  !@cuf istat = cudaMemAdvise(dwdtrko, size(dwdtrko), cudaMemAdviseSetPreferredLocation, mydev)
  !@cuf istat = cudaMemAdvise(rhsbp_x, size(rhsbp_x), cudaMemAdviseSetPreferredLocation, mydev)
  !@cuf istat = cudaMemAdvise(rhsbp_y, size(rhsbp_y), cudaMemAdviseSetPreferredLocation, mydev)
  !@cuf istat = cudaMemAdvise(rhsbp_z, size(rhsbp_z), cudaMemAdviseSetPreferredLocation, mydev)
  !@cuf istat = cudaMemAdvise(lambdaxyp, size(lambdaxyp), cudaMemAdviseSetPreferredLocation, mydev)
  !@cuf istat = cudaMemAdvise(dzc, size(dzc), cudaMemAdviseSetReadMostly, 0)
  !@cuf istat = cudaMemAdvise(dzf, size(dzf), cudaMemAdviseSetReadMostly, 0)
  !@cuf istat = cudaMemAdvise(dzci, size(dzci), cudaMemAdviseSetReadMostly, 0)
  !@cuf istat = cudaMemAdvise(dzfi, size(dzfi), cudaMemAdviseSetReadMostly, 0)
  !
  !@cuf istat = cudaMemAdvise(zc_g , size(zc_g) , cudaMemAdviseSetPreferredLocation, mydev)
  !@cuf istat = cudaMemAdvise(zf_g , size(zf_g) , cudaMemAdviseSetPreferredLocation, mydev)
  !
#if defined(_USE_DI)
  !@cuf istat = cudaMemAdvise(mu, size(mu), cudaMemAdviseSetPreferredLocation, mydev)
  !@cuf istat = cudaMemAdvise(rho, size(rho), cudaMemAdviseSetPreferredLocation, mydev)
  !@cuf istat = cudaMemAdvise(kappa, size(kappa), cudaMemAdviseSetPreferredLocation, mydev)
  !@cuf istat = cudaMemAdvise(psi, size(psi), cudaMemAdviseSetPreferredLocation, mydev)
  !@cuf istat = cudaMemAdvise(nor, size(nor), cudaMemAdviseSetPreferredLocation, mydev)
  !@cuf istat = cudaMemAdvise(ls, size(ls), cudaMemAdviseSetPreferredLocation, mydev)
  !@cuf istat = cudaMemAdvise(dpsidtrko, size(dpsidtrko), cudaMemAdviseSetPreferredLocation, mydev)
#endif
  !
  if(myid.eq.0) print*, '************************************************'
  if(myid.eq.0) print*, '*** Beginning of simulation (TWO-PHASE mode) ***'
  if(myid.eq.0) print*, '************************************************'
  !
#if defined(_OPENACC)
  if(myid.eq.0) then
    print*, ' GPU accelerated version, grid size:', n(1)*dims(1), n(2)*dims(2), n(3)*dims(3)
  endif
#endif
  !
#if defined(_OPENACC)
  !
  ! Allocate buffers for halo communications (GPU-only)
  !
  call alloc_buf(n,nh_d)
  !
#else
  !
  ! halo generation using MPI derivate datatypes (CPU-only)
  !
  call halo_gen(n,nh_u ,halo_u )
  call halo_gen(n,nh_p ,halo_p )
  call halo_gen(n,nh_v ,halo_v )
  call halo_gen(n,nh_d ,halo_d )
  !
#endif
  !
  ! initialize the grid (using global variables along z)
  !
  call initgrid(inivel,ng(3),gr,lz,nh_d,dzc_g,dzf_g,zc_g,zf_g) 
  !
  if(myid.eq.0) then
    open(99,file=trim(datadir)//'grid.bin',access='direct',recl=4*ng(3)*sizeof(1._rp))
    write(99,rec=1) dzc_g(1:ng(3)),dzf_g(1:ng(3)),zc_g(1:ng(3)),zf_g(1:ng(3))
    close(99)
    open(99,file=trim(datadir)//'grid.out')
    do kk=1-nh_d,ng(3)+nh_d
      write(99,'(5E15.7)') 1._rp*kk,zf_g(kk),zc_g(kk),dzf_g(kk),dzc_g(kk)
    enddo
    close(99)
    open(99,file=trim(datadir)//'geometry.out')
      write(99,*) ng(1),ng(2),ng(3)
      write(99,*) l(1),l(2),l(3)
    close(99)
  endif
  !@cuf istat = cudaMemPrefetchAsync(zc_g , size( zc_g), mydev, 0) I did boiling with 2000 density ratio and evaporation with 1000
  !@cuf istat = cudaMemPrefetchAsync(zf_g , size( zf_g), mydev, 0)
  !
  do k=1-nh_d,ng3+nh_d
    dzfi_g(k) = 1._rp/dzf_g(k)
    dzci_g(k) = 1._rp/dzc_g(k)
  enddo
  !
  ! compute the spacing along z in local coordinates
  !
  do k=1-nh_d,n3+nh_d
    kk      = k + ijk_start(3)
    zc(k)   = zc_g(kk)
    zf(k)   = zf_g(kk) 
    dzf(k)  = dzf_g(kk)
    dzc(k)  = dzc_g(kk)
    dzfi(k) = 1._rp/dzf(k)
    dzci(k) = 1._rp/dzc(k)
  enddo
  !@cuf istat = cudaMemPrefetchAsync(dzci, size(dzci), mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(dzfi, size(dzfi), mydev, 0)
  !

  ! test input files before proceeding with the calculation
  !
  call test_sanity(ng,n,dims_xyz(:,3),ipencil,nh_d,nh_u,nh_p,halo_d,halo_u,halo_p,stop_type, &
                   cbcvel,cbcpre,bcvel,bcpre,is_outflow,is_forced,dli,dzci_g,dzfi_g)
  !
  if(.not.restart) then
    !
    istep = 0
    time  = 0._rp
      !
      ! Initialize DI 
      !
#if defined(_TWO_PHASE)
      call initpsi(inipsi,n,dl,l,nh_v,eps_int,psi)
#if defined(_CONTACT_LINE)
      call static_contact_angle_m(n,dl,nh_v,eps_int,theta,psi,psi_bound)
      call boundsb_array(n,psi_bound,nh_v,halo_v,psi)
#else
      call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,psi)
#endif
#if defined(_SCALAR_TRANS)      
      call initc1('lin',n,dl,l,nh_v,eps_int,psi,c1,c2)
      call boundp(cbcsca,n,bcsca,nh_d,nh_v,halo_v,dl,dzc,dzf,c1)
      call boundp(cbcsca,n,bcsca,nh_d,nh_v,halo_v,dl,dzc,dzf,c2)
#else
    c1(:,:,:) = 0.d0
    c2(:,:,:) = 0.d0
#endif      
#if defined(_USE_ACDI)
      call psi_to_ls(n,dl,nh_p,psi,eps_int,ls)
      call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,ls )
#else
    ls(:,:,:) = 0.d0
#endif
#else
    psi(:,:,:) = 0.d0
    ls(:,:,:)  = 0.d0
    c1(:,:,:)  = 0.d0
    c2(:,:,:)  = 0.d0
#endif
    !
    call initflow(inivel,n(1),n(2),n(3),dims,nh_d,nh_u,nh_p,rho2,mu2,zc/lz,dzc/lz,dzf/lz,u,v,w,p)
    !
    call bounduvw(cbcvel,n,bcvel,nh_d,nh_u,halo_u,is_outflow,dl,dzc,dzf,u,v,w)
    !
    ! set to zeros the rhs of momentum equation 
    ! (only for the first time-step, not for the restarting)
    ! 
    !$acc kernels
    do k=1,n3
      do j=1,n2
        do i=1,n1
          dudtrko(i,j,k) = 0._rp
          dvdtrko(i,j,k) = 0._rp
          dwdtrko(i,j,k) = 0._rp
          dpsidtrko(i,j,k) = 0.d0
          dm2dtrko(i,j,k) = 0.d0
          dc1dtrko(i,j,k)  = 0.d0
          dc2dtrko(i,j,k)  = 0.d0
          dentdtrko(i,j,k) = 0.d0
          rhs_div_th(i,j,k) = 0.d0
        enddo
      enddo
    enddo
    
    call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,rhs_div_th)

   p0 = 101325.0     !!background thermodynamic pressure initialization
   dp0dt = 0.0
   tmp_ini = tl0
  
#if defined(_EVAP)
    
    omega(:,:,:) = 0.0
    
    do i=1,n(1)
      do j=1,n(2)
        do k=1,n(3)

            !Anotine equation of state, coefficients are for water between 273 k -373 K
            p_sat = (10**(8.07131-(1730.63/(233.426+tmp_ini-273.0))))*133.322
            o_sat(i,j,k) =   p_sat/(p_sat+(p0-p_sat)*(mol_mas_air/mol_mas_vap))

       enddo
      enddo
    enddo

    call boundp(cbcsca,n,bcsca,nh_d,nh_v,halo_v,dl,dzc,dzf,o_sat)
    
   call initvap(n,dli,nh_d,nh_u,nh_v,halo_v,dzc,dzf,diff2,rho_air,psi,nor,c2,omega,t_relax,o_sat)

    do k=1,n3
      do j=1,n2
        do i=1,n1

           omega(i,j,k) =  min(omega(i,j,k),o_sat(i,j,k))
           omeg(i,j,k) =  omega(i,j,k)
        enddo
      enddo
    enddo
    
    call boundp(cbcsca,n,bcsca,nh_d,nh_v,halo_v,dl,dzc,dzf,omega)
    call boundp(cbcsca,n,bcsca,nh_d,nh_v,halo_v,dl,dzc,dzf,omeg)
#else

     omega(:,:,:) = 0.0
 
#endif  



#if defined(_LOW_MACH)

        !!!!!!!!!!!!!    LOW MACH FLAG START     

#if defined(_EVAP)

    print*, "Low Mach Evaporation"
   
       !!!!!!!!!!!!!    LOW MACH EVAPORATION FLAG START 
    
    
    !inside update property - 1 is arithmetic avg. , 2 is harmonic avg

    !molar mass initiation , based on composition

    call update_property(n,(/mol_mas_vap,mol_mas_air/),omega,mol_mas,2)
    call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,mol_mas)

    !rho_gas initiation, based on composition

    call update_property(n,(/rho_vap,rho_air/),omega,rho_gas,1)
    call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,rho_gas)

    !bulk density initiation, based on liquid and gas density(variable)

    call update_property_var(n,psi,rho1,rho_gas,rho)
    call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,rho)

    !!!!!!!!!!!!!    LOW MACH EVAPORATION FLAG END 

#else 
    print*, "Low Mach No Phase Change"
    !!!!!!!!!!!!!    LOW MACH NO EVAPORATION FLAG START 

    mol_mas(:,:,:) = mol_mas_air
    rho_gas(:,:,:) = min(rho1,rho2)

    !bulk density

    call update_property_var(n,psi,rho1,rho_gas,rho)
    call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,rho)

    !!!!!!!!!!!!!    LOW MACH NO EVAPORATION FLAG END

#endif

#else
    !!!!!!!!!!!!!    INCOMPRESSIBLE FLAG START

    mol_mas(:,:,:) = mol_mas_air
    rho_gas(:,:,:) = min(rho1,rho2)

    !!!!!!!!!!!!!!   INCOMPRESSIBLE EVAPORATION FLAG START

#if defined(_EVAP)


    print*, "Incompressible Evaporation"

    !bulk density initialization

    call update_property(n,(/rho1,rho2/),psi,rho,1)
    call boundp(cbctmp,n,bctmp,nh_d,nh_v,halo_v,dl,dzc,dzf,rho)

    !!!!!!!!!!!!!!   INCOMPRESSIBLE EVAPORATION FLAG END

#else
    
    print*, "Incompressible No Phase Change"

    !!!!!!!!!!!!!!!   INCOMPRESSIBLE NO EVAPORATION FLAG START

    !bulk density initialization

    call update_property(n,(/rho1,rho2/),psi,rho,1)
    call boundp(cbctmp,n,bctmp,nh_d,nh_v,halo_v,dl,dzc,dzf,rho)
    
    !!!!!!!!!!!!!!!    INCOMPRESSIBLE NO EVAPORATION FLAG START
#endif  


#endif



  do k=1,n3
    do j=1,n2
     do i=1,n1
           
        !c2(i,j,k) = rho_gas(i,j,k)*(1.d0-psi(i,j,k))*omega(i,j,k)
        c2(i,j,k) = (1.d0-psi(i,j,k))*omega(i,j,k)
        m2(i,j,k) = rho_gas(i,j,k)*(1.d0-psi(i,j,k))
        psi2(i,j,k) = (1.d0-psi(i,j,k))
       
      enddo
    enddo
  enddo
  !
  call boundp(cbcsca,n,bcsca,nh_d,nh_v,halo_v,dl,dzc,dzf,c2)      !This should be discussed
  call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,m2)
  call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,psi2)    
  
           

  call update_property(n,(/tl0,tg0/),psi,tmp,1)
  call boundp(cbctmp,n,bctmp,nh_d,nh_v,halo_v,dl,dzc,dzf,tmp)

#if defined(_HEAT_TRANSFER)     


#if defined (_EVAP)
     
     !thermal conductivity

     call update_property(n,(/cond_vap,cond_air/),omega,cond_gas,1)
     call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,cond_gas)

     call update_property_var(n,psi,cond1,cond_gas,cond)
     call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,cond)

     !cp

     call update_property(n,(/cp_vap,cp_air/),omega,cp_gas,1)
     call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,cp_gas)
     
     call update_property_var(n,psi,cp1,cp_gas,cpp)
     call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,cpp)

     !cv

     call update_property(n,(/cv_vap,cv_air/),omega,cv_gas,1)
     call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,cv_gas)

     call update_property_var(n,psi,cv1,cv_gas,cv)
     call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,cv)



     
#else
    
    !thermal conductivity
    
    call update_property(n,(/cond1,cond2/),psi,cond,1)
    call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,cond)
    
    !cp

    call update_property(n,(/cp1,cp2/),psi,cpp,1)
    call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,cpp)  
     
    cp_gas(:,:,:) = min(cp1,cp2)

    !cv
    call update_property(n,(/cv1,cv2/),psi,cv,1)
    call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,cv) 
    
    cv_gas(:,:,:) = min(cv1,cv2)

#endif

      !enthalpy of formation

     ! call update_property(n,(/h01*rho1,h02*rho2/),psi,lh,1)
     ! call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,lh)
 
      

      !rhocp calculation 

     ! call update_property(n,(/rho1*cp1,rho2*cp2/),psi,rhocp,1)
     ! call boundp(cbcent,n,bcent,nh_d,nh_v,halo_v,dl,dzc,dzf,rhocp)

#if defined(_LOW_MACH)

     do k=1,n(3)
       do j=1,n(2)
         do i=1,n(1)

           rhocp(i,j,k) = rho1*cp1*psi(i,j,k) + m2(i,j,k) * cp_gas(i,j,k)
           lh(i,j,k) = rho1*h01*psi(i,j,k) + m2(i,j,k)*h02
         end do
       end do
     end do
#else
     call update_property(n,(/rho1*cp1,rho2*cp2/),psi,rhocp,1)
     call update_property(n,(/h01*rho1,h02*rho2/),psi,lh,1)

#endif
     call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,rhocp)
     call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,lh)


      !temperature initialization

   !   call update_property(n,(/tl0,tg0/),psi,tmp,1)
   !   call boundp(cbctmp,n,bctmp,nh_d,nh_v,halo_v,dl,dzc,dzf,tmp)
      
     !for now just consider tl0 as t1 and tg0 as t2; will change the variable name later



      !enthalpy initialization 

    !  h1 = rho1*(cp1*tl0+h01)
    !  h2 = rho2*(cp2*tg0+h02)

   !   call update_property(n,(/h1,h2/),psi,ent,1)
   !   call boundp(cbctmp,n,bctmp,nh_d,nh_v,halo_v,dl,dzc,dzf,ent)
   !   call boundp_di(cbctmp,n,bctmp,nh_d,nh_p,halo_p,dl,dzc,dzf,ent,rhocp,lh,tmp) !impose b.c on "ent"
   

   do i=1,n(1)
      do j=1,n(2)
       do k=1,n(3)
       
        ent(i,j,k) = rhocp(i,j,k)*tmp(i,j,k) + lh(i,j,k)
   !    tmp(i,j,k) = (ent(i,j,k)-lh(i,j,k))/rhocp(i,j,k)

       enddo
      enddo
     enddo
    
    call boundp_di(cbctmp,n,bctmp,nh_d,nh_p,halo_p,dl,dzc,dzf,ent,rhocp,lh,tmp) ! b.c imposed on enthalpy


    do i=1,n(1)
      do j=1,n(2)
       do k=1,n(3)

       tmp(i,j,k) = (ent(i,j,k)-lh(i,j,k))/rhocp(i,j,k)

       enddo
      enddo
     enddo

     call boundp(cbctmp,n,bctmp,nh_d,nh_v,halo_v,dl,dzc,dzf,tmp)

#endif


      !!initial source term for evaporation

  tmp_min = minval(tmp(1:n1,1:n2,1:n3))
  call mpi_allreduce(MPI_IN_PLACE,tmp_min,1,MPI_REAL_RP,MPI_MIN,MPI_COMM_WORLD,ierr)
  tmp_max = maxval(tmp(1:n1,1:n2,1:n3))
  call mpi_allreduce(MPI_IN_PLACE,tmp_max,1,MPI_REAL_RP,MPI_MIN,MPI_COMM_WORLD,ierr)

#if defined(_EVAP)
  call src_evap(n,nh_v,tmp,src_ph,src_mod,t_relax,psi,omega,rho_gas,m2,tmp_min,tmp_max,p0,c2,o_sat)
#else
  src_ph(:,:,:) = 0.0
#endif
  call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,src_ph)
  call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,src_mod)

   vol_gas  =  sum(psi2(1:n1,1:n2,1:n3))*dl(1)*dl(2)*dl(3)
   call mpi_allreduce(MPI_IN_PLACE,vol_gas ,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)

  !  sum_rt  =  sum(psi2(1:n1,1:n2,1:n3)/(8.314*tmp(1:n1,1:n2,1:n3)))*dl(1)*dl(2)*dl(3)
  !  call mpi_allreduce(MPI_IN_PLACE,sum_rt ,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
   
!  mass_gas = (p0*vol_gas)/(287.0*tg) 
   mass_gas = vol_gas*rho_air
   p0 = (mass_gas*287.0*tg0)/vol_gas

!  print*,'the pressure is', p0
  !
    !$acc end kernels
    !
    if(myid.eq.0) print*, '*** Initial condition succesfully set ***'
    !
endif 
    !

if(restart_sph) then 
    !
    action_load = 'r'
    call load(action_load,trim(restart_dir)//'fldu.bin',n,      u(1:n(1),1:n(2),1:n(3)))
    call load(action_load,trim(restart_dir)//'fldv.bin',n,      v(1:n(1),1:n(2),1:n(3)))
    call load(action_load,trim(restart_dir)//'fldw.bin',n,      w(1:n(1),1:n(2),1:n(3)))
    call load(action_load,trim(restart_dir)//'flddu.bin',n,dudtrko(1:n(1),1:n(2),1:n(3)))
    call load(action_load,trim(restart_dir)//'flddv.bin',n,dvdtrko(1:n(1),1:n(2),1:n(3)))
    call load(action_load,trim(restart_dir)//'flddw.bin',n,dwdtrko(1:n(1),1:n(2),1:n(3)))
    call load(action_load,trim(restart_dir)//'fldp.bin',n,      p(1:n(1),1:n(2),1:n(3)))
    call load_scalar(action_load,trim(restart_dir)//'scalar.out',time,istep,dto)
    !
#if defined(_TWO_PHASE)
      call initpsi(inipsi,n,dl,l,nh_v,eps_int,psi)
#if defined(_CONTACT_LINE)
      call static_contact_angle_m(n,dl,nh_v,eps_int,theta,psi,psi_bound)
      call boundsb_array(n,psi_bound,nh_v,halo_v,psi)
#else
      call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,psi)
#endif
#if defined(_SCALAR_TRANS)      
      call initc1('lin',n,dl,l,nh_v,eps_int,psi,c1,c2)
      call boundp(cbcsca,n,bcsca,nh_d,nh_v,halo_v,dl,dzc,dzf,c1)
      call boundp(cbcsca,n,bcsca,nh_d,nh_v,halo_v,dl,dzc,dzf,c2)
#else
    c1(:,:,:) = 0.d0
    c2(:,:,:) = 0.d0
#endif      
#if defined(_USE_ACDI)
      call psi_to_ls(n,dl,nh_p,psi,eps_int,ls)
      call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,ls )
#else
    ls(:,:,:) = 0.d0
#endif
#else
    psi(:,:,:) = 0.d0
    ls(:,:,:)  = 0.d0
    c1(:,:,:)  = 0.d0
    c2(:,:,:)  = 0.d0
#endif
    !
    if(myid.eq.0) print*, '*** single phase Checkpoint loaded at time = ', time, 'time step = ', istep, '. ***'
    !
endif    
    !
if(restart_mph) then 

    action_load = 'r'
    call load(action_load,trim(restart_dir)//'fldu.bin',n,      u(1:n(1),1:n(2),1:n(3)))
    call load(action_load,trim(restart_dir)//'fldv.bin',n,      v(1:n(1),1:n(2),1:n(3)))
    call load(action_load,trim(restart_dir)//'fldw.bin',n,      w(1:n(1),1:n(2),1:n(3)))
    call load(action_load,trim(restart_dir)//'flddu.bin',n,dudtrko(1:n(1),1:n(2),1:n(3)))
    call load(action_load,trim(restart_dir)//'flddv.bin',n,dvdtrko(1:n(1),1:n(2),1:n(3)))
    call load(action_load,trim(restart_dir)//'flddw.bin',n,dwdtrko(1:n(1),1:n(2),1:n(3)))
    call load(action_load,trim(restart_dir)//'fldp.bin',n,      p(1:n(1),1:n(2),1:n(3)))
    call load(action_load,trim(restart_dir)//'fldpsi.bin',n,    psi(1:n(1),1:n(2),1:n(3)))
    call load(action_load,trim(restart_dir)//'flddpsi.bin',n,dpsidtrko(1:n(1),1:n(2),1:n(3)))
    call load_scalar(action_load,trim(restart_dir)//'scalar.out',time,istep,dto)

#if defined(_CONTACT_LINE)
      call static_contact_angle_m(n,dl,nh_v,eps_int,theta,psi,psi_bound)
      call boundsb_array(n,psi_bound,nh_v,halo_v,psi)
#else
      call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,psi)
#endif
#if defined(_SCALAR_TRANS)
      call initc1('lin',n,dl,l,nh_v,eps_int,psi,c1,c2)
      call boundp(cbcsca,n,bcsca,nh_d,nh_v,halo_v,dl,dzc,dzf,c1)
      call boundp(cbcsca,n,bcsca,nh_d,nh_v,halo_v,dl,dzc,dzf,c2)
#else
    c1(:,:,:) = 0.d0
    c2(:,:,:) = 0.d0
#endif
#if defined(_USE_ACDI)
      call psi_to_ls(n,dl,nh_p,psi,eps_int,ls)
      call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,ls )
#else
    ls(:,:,:) = 0.d0
#endif
    if(myid.eq.0) print*, '*** multiphase Checkpoint loaded at time = ', time, 'time step = ', istep, '. ***'

endif 
if(restart_mph_scalar) then 
    action_load = 'r'
    call load(action_load,trim(restart_dir)//'fldu.bin',n,      u(1:n(1),1:n(2),1:n(3)))
    call load(action_load,trim(restart_dir)//'fldv.bin',n,      v(1:n(1),1:n(2),1:n(3)))
    call load(action_load,trim(restart_dir)//'fldw.bin',n,      w(1:n(1),1:n(2),1:n(3)))
    call load(action_load,trim(restart_dir)//'flddu.bin',n,dudtrko(1:n(1),1:n(2),1:n(3)))
    call load(action_load,trim(restart_dir)//'flddv.bin',n,dvdtrko(1:n(1),1:n(2),1:n(3)))
    call load(action_load,trim(restart_dir)//'flddw.bin',n,dwdtrko(1:n(1),1:n(2),1:n(3)))
    call load(action_load,trim(restart_dir)//'fldp.bin',n,      p(1:n(1),1:n(2),1:n(3)))
    call load(action_load,trim(restart_dir)//'fldpsi.bin',n,    psi(1:n(1),1:n(2),1:n(3)))
    call load(action_load,trim(restart_dir)//'flddpsi.bin',n,dpsidtrko(1:n(1),1:n(2),1:n(3)))
    call load(action_load,trim(restart_dir)//'fldc1.bin',n,    c1(1:n(1),1:n(2),1:n(3)))
    call load(action_load,trim(restart_dir)//'fldc2.bin',n,    c2(1:n(1),1:n(2),1:n(3)))
    call load(action_load,trim(restart_dir)//'flddc1.bin',n,dc1dtrko(1:n(1),1:n(2),1:n(3)))
    call load(action_load,trim(restart_dir)//'flddc2.bin',n,dc2dtrko(1:n(1),1:n(2),1:n(3)))
    call load_scalar(action_load,trim(restart_dir)//'scalar.out',time,istep,dto)
#if defined(_CONTACT_LINE)
      call static_contact_angle_m(n,dl,nh_v,eps_int,theta,psi,psi_bound)
      call boundsb_array(n,psi_bound,nh_v,halo_v,psi)
#else
      call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,psi)
#endif
#if defined(_USE_ACDI)
      call psi_to_ls(n,dl,nh_p,psi,eps_int,ls)
      call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,ls )
#else
    ls(:,:,:) = 0.d0
#endif

    if(myid.eq.0) print*, '*** multiphase-scalar-transport Checkpoint loaded at time = ', time, 'time step = ', istep, '. ***'
endif
if(restart_mph_heat) then
    action_load = 'r'
    call load(action_load,trim(restart_dir)//'fldu.bin',n,      u(1:n(1),1:n(2),1:n(3)))
    call load(action_load,trim(restart_dir)//'fldv.bin',n,      v(1:n(1),1:n(2),1:n(3)))
    call load(action_load,trim(restart_dir)//'fldw.bin',n,      w(1:n(1),1:n(2),1:n(3)))
    call load(action_load,trim(restart_dir)//'flddu.bin',n,dudtrko(1:n(1),1:n(2),1:n(3)))
    call load(action_load,trim(restart_dir)//'flddv.bin',n,dvdtrko(1:n(1),1:n(2),1:n(3)))
    call load(action_load,trim(restart_dir)//'flddw.bin',n,dwdtrko(1:n(1),1:n(2),1:n(3)))
    call load(action_load,trim(restart_dir)//'fldp.bin',n,      p(1:n(1),1:n(2),1:n(3)))
    call load(action_load,trim(restart_dir)//'fldpsi.bin',n,    psi(1:n(1),1:n(2),1:n(3)))
    call load(action_load,trim(restart_dir)//'flddpsi.bin',n,dpsidtrko(1:n(1),1:n(2),1:n(3)))
    call load(action_load,trim(restart_dir)//'fldtmp.bin',n,    tmp(1:n(1),1:n(2),1:n(3)))
    call load(action_load,trim(restart_dir)//'flddtmp.bin',n,dentdtrko(1:n(1),1:n(2),1:n(3)))
    call load_scalar(action_load,trim(restart_dir)//'scalar.out',time,istep,dto)
#if defined(_CONTACT_LINE)
      call static_contact_angle_m(n,dl,nh_v,eps_int,theta,psi,psi_bound)
      call boundsb_array(n,psi_bound,nh_v,halo_v,psi)
#else
      call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,psi)
#endif
#if defined(_USE_ACDI)
      call psi_to_ls(n,dl,nh_p,psi,eps_int,ls)
      call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,ls )
#else
    ls(:,:,:) = 0.d0
#endif

    if(myid.eq.0) print*, '*** multiphase-scalar-transport Checkpoint loaded at time = ', time, 'time step = ', istep, '. ***'
endif
  !
  !@cuf istat = cudaMemPrefetchAsync(dudtrko, size(dudtrko), mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(dvdtrko, size(dvdtrko), mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(dwdtrko, size(dwdtrko), mydev, 0)
  !
  ! set boundary conditions on the initial/loaded fields
  !
  call bounduvw(cbcvel,n,bcvel,nh_d,nh_u,halo_u,is_outflow,dl,dzc,dzf,u,v,w)
  call boundp(cbcpre,n,bcpre,nh_d,nh_p,halo_p,dl,dzc,dzf,p)
  !
  ! for the first time-step and the restarting we use a 0th order extrapolation 
  ! in time-splitting of the pressure equation
  !
  !$acc kernels 
  do k=1,n3
    do j=1,n2
      do i=1,n1
        pold(i,j,k) = p(i,j,k)
      end do
    end do
  end do
  !$acc end kernels
  !
  call boundp(cbcpre,n,bcpre,nh_d,nh_p,halo_p,dl,dzc,dzf,pold)
  !
  !
  call update_property(n,(/mu1 ,mu2 /),psi,mu,1 ) 
  call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,mu)
  !
#if defined(_HEAT_TANSFER)
  call update_property(n,(/cond1,cond2/),psi,cond,1)
  call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,cond)
  call update_property(n,(/cp1,cp2/),psi,cpp,1)
  call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,cpp)
  call update_property(n,(/rho1*cp1,rho2*cp2/),psi,rhocp,1)
  call boundp(cbcent,n,bcent,nh_d,nh_v,halo_v,dl,dzc,dzf,rhocp)
  call update_property(n,(/h01*rho1,h02*rho2/),psi,lh,1)
  call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,lh)

#endif
  !
  !
  ! post-process and write initial condition
#if !defined(_BENCHMARK_NO_IO)
  !
  call profiler_start("OUT:initial", tag = .true., tag_color = COLOR_WHITE)
  !
  ! Prefetching back pre IO
  !@cuf istat = cudaMemPrefetchAsync(u, size(u), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(v, size(v), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(w, size(w), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(p, size(p), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(psi, size(psi), cudaCpuDeviceId, 0)
  !
#if defined(_TWO_PHASE)
#if defined(_USE_ACDI)  
  call cmpt_norm(n,dli,nh_p,ls,kappa,nor)
#else  
  call cmpt_norm(n,dli,nh_p,psi,kappa,nor)
#endif
  call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,kappa       )
  call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,nor(:,:,:,1))
  call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,nor(:,:,:,2))
  call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,nor(:,:,:,3))
#else
  nor(:,:,:,:) = 0.d0
  kappa(:,:,:) = 0.d0
#endif
  !
  write(fldnum,'(i9.9)') istep
  include 'out1d.h90'
  include 'out2d.h90'
  include 'out3d.h90'
  !
  ! Prefetching post IO
  !@cuf istat = cudaMemPrefetchAsync(u, size(u), mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(v, size(v), mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(w, size(w), mydev, 0)
  !
  call profiler_stop("OUT:initial")
#endif
  !
#if defined(_TWO_PHASE)
#if defined(_DO_POSTPROC)
  if(mod(istep,iout0d_ta).eq.0.and.do_tagging) then
    call droplet_tagging(n,dims,datadir_ta,dl,nh_d,nh_v,nh_u,halo_v,dzc,dzf,psi,u,v,w,istep,time)
  endif
#endif
#endif
  !
  ! compute an initial time-step
  !
  if(.not.constant_dt) then
#if defined(_TWO_PHASE)          
    call chkdt_tw(n(1),n(2),n(3),dxi,dyi,dzi,nh_d,nh_u,dzci,dzfi,u,v,w,dtmax)
#else    
    call chkdt_sp(n(1),n(2),n(3),dxi,dyi,dzi,nh_d,nh_u,dzci,dzfi,u,v,w,dtmax)
#endif    
    dt = cfl*dtmax
  else
    if(myid.eq.0) print*, 'the simulation is run at constant time-step'
    dtmax = dt_input
    dt    = dtmax
  endif
  if(istep.eq.0) dto = dt
  if(myid.eq.0) print*, 'dtmax = ', dtmax, 'dt = ', dt
  dti  = 1._rp/dt
  kill = .false.
  !
  ! Prefetching post IO
  !@cuf istat = cudaMemPrefetchAsync(psi, size(psi), mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(p, size(p), mydev, 0)
  !
  ! preliminary checks
  ! 
  psi_min = minval(psi(1:n1,1:n2,1:n3))
  call mpi_allreduce(MPI_IN_PLACE,psi_min,1,MPI_REAL_RP,MPI_MIN,MPI_COMM_WORLD,ierr)
  psi_max = maxval(psi(1:n1,1:n2,1:n3))
  call mpi_allreduce(MPI_IN_PLACE,psi_max,1,MPI_REAL_RP,MPI_MAX,MPI_COMM_WORLD,ierr)
  vol_p1  = sum(psi(1:n1,1:n2,1:n3))*dl(1)*dl(2)*dl(3)
  vol_c1  = sum(c1(1:n1,1:n2,1:n3))*dl(1)*dl(2)*dl(3)
  vol_c2  = sum(c2(1:n1,1:n2,1:n3))*dl(1)*dl(2)*dl(3)
  vol_ct  = sum(ent(1:n1,1:n2,1:n3))*dl(1)*dl(2)*dl(3)
  call mpi_allreduce(MPI_IN_PLACE,vol_p1 ,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
  call mpi_allreduce(MPI_IN_PLACE,vol_c1 ,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
  call mpi_allreduce(MPI_IN_PLACE,vol_c2 ,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
  call mpi_allreduce(MPI_IN_PLACE,vol_ct ,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
  !
  var(:) = 0._rp
  var(1) = 1._rp*istep
  var(2) = dt
  var(3) = time
  var(4) = psi_min
  var(5) = psi_max
  var(6) = vol_p1
  var(7) = vol_c1
  var(8) = vol_c2
  var(9) = p0
  var(10) = dp0dt
  call out0d(trim(datadir)//'psi_info.out',9,var)
  !
  ! initialize Poisson solver
  ! and deallocate global arrays (not needed anymore) 
  !
  call initsolver(n,dims,dims_xyz(:,3),dli,nh_d,dzci_g,dzfi_g,cbcpre,bcpre(:,:),(/'c','c','c'/),lambdaxyp, & 
                  ap,bp,cp,arrplanp,normfftp,rhsbp_x,rhsbp_y,rhsbp_z)
  !@cuf istat = cudaMemPrefetchAsync(rhsbp_x, size(rhsbp_x), mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(rhsbp_y, size(rhsbp_y), mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(rhsbp_z, size(rhsbp_z), mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(lambdaxyp, size(lambdaxyp), mydev, 0)
  deallocate(dzc_g,dzf_g,dzci_g,dzfi_g)
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!main loop!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  !                
  if(myid.eq.0) print*,'*** Calculation loop starts now ***'
  is_done = .false.
  !
  !@cuf istat = cudaMemGetInfo( freeMem, totMem )
  !@cuf if(myid.eq.0) print*, 'Used memory = ', totMem - freeMem
  !
  !call exit
  do while(.not.is_done)
    !
    !
#if defined(_TIMING)
    dt12  = MPI_WTIME()
#endif
    !
    istep = istep + 1
    !
    call profiler_start("STEP", tag = .true., tag_color = COLOR_WHITE)
    !
    time  = time + dt
    !
    if(myid.eq.0) print*, 'Timestep #', istep, 'Time = ', time
    !
    if(any(is_forced(:))) dpdl_c(1:3) = 0._rp
    !
    ! 0. compute the coefficients for the time advancement
    !
    call cmpt_time_factors(time_scheme,restart,istep,1,rkcoeff,dt,dto,f_t1,f_t2,f_t12)
    f_t12_i = 1._rp/f_t12
    f_t12_o = dto
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    1. Scalar advection and properties update --> c1,c2^(n+1) 
    !
    call profiler_start("DI", tag = .true., tag_color = COLOR_YELLOW)
    !
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  1a. Update properties at time-step (n)
   ! call update_property(n,(/rho1 ,rho2 /),psi,rho,1)
   ! call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,rho)

    call update_property(n,(/mu1 ,mu2 /),psi,mu,1)
    call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,mu)
    !
!#if defined(_HEAT_TRANSFER)
   ! call update_property(n,(/cond1,cond2/),psi,cond)
   ! call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,cond)
    !
   ! call update_property(n,(/cp1,cp2/),psi,cpp)
   ! call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,cpp)

!    call update_property(n,(/rho1*cp1,rho2*cp2/),psi,rhocp)
!    call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,rhocp)

!    call update_property_gas_mix(n,cp_vap,cp_air,cp1,omega,psi,psi2,cp_gas,cpp)
!    call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,cp_gas)
!    call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,cpp)

!    call update_property_gas_mix(n,cv_vap,cv_air,cv1,omega,psi,psi2,cv_gas,cv)
!    call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,cv_gas)
!    call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,cv)

   ! call update_property_gas_mix(n,rho_vap*cp_vap,rho_air*cp_air,rho1*cp1,omega,psi,rhocp_gas,rhocp)
   ! call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,rhocp_gas)
   ! call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,rhocp)

!    call update_property_gas_mix(n,cond_vap,cond_air,cond1,omega,psi,psi2,cond_gas,cond)
!    call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,cond_gas)
!    call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,cond)
    !
!    call update_property(n,(/h01*rho1,h02*rho2/),psi,lh)
!    call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,lh)
!#endif    

#if defined(_HEAT_TRANSFER)

#if defined (_EVAP)

     !thermal conductivity

     call update_property(n,(/cond_vap,cond_air/),omega,cond_gas,1)
     call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,cond_gas)

     call update_property_var(n,psi,cond1,cond_gas,cond)
     call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,cond)

     !cp

     call update_property(n,(/cp_vap,cp_air/),omega,cp_gas,1)
     call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,cp_gas)

     call update_property_var(n,psi,cp1,cp_gas,cpp)
     call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,cpp)

     !cv

     call update_property(n,(/cv_vap,cv_air/),omega,cv_gas,1)
     call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,cv_gas)

     call update_property_var(n,psi,cv1,cv_gas,cv)
     call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,cv)

#else

    !thermal conductivity

    call update_property(n,(/cond1,cond2/),psi,cond,1)
    call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,cond)

    !cp

    call update_property(n,(/cp1,cp2/),psi,cpp,1)
    call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,cpp)

    cp_gas(:,:,:) = min(cp1,cp2)
    call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,cp_gas)

    !cv

    call update_property(n,(/cv1,cv2/),psi,cv,1)
    call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,cv)
    
    cv_gas(:,:,:) = min(cv1,cv2)
    call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,cv_gas)

#endif    

    !enthalpy of formation

!   call update_property(n,(/h01*rho1,h02*rho2/),psi,lh,1)
!   call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,lh)

    !rhocp calculation

#if defined(_LOW_MACH)

     do k=1,n(3)
       do j=1,n(2)
         do i=1,n(1)

           rhocp(i,j,k) = rho1*cp1*psi(i,j,k) + m2(i,j,k) * cp_gas(i,j,k)
           lh(i,j,k) = rho1*h01*psi(i,j,k) + m2(i,j,k)*h02

         end do
       end do
     end do


#else
     call update_property(n,(/rho1*cp1,rho2*cp2/),psi,rhocp,1)
     call update_property(n,(/h01*rho1,h02*rho2/),psi,lh,1)

#endif
     call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,rhocp)
     call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,lh)

    
        !heat transfer flag end
#endif

    !
#if defined(_TWO_PHASE)    
    !
#if defined(_USE_ACDI)
    call cmpt_norm(n,dli,nh_p,ls,kappa,nor)
#else
    call cmpt_norm(n,dli,nh_p,psi,kappa,nor)
#endif
    call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,kappa       )
    call boundp(cbcnor,n,bcnor,nh_d,nh_v,halo_v,dl,dzc,dzf,nor(:,:,:,1))
    call boundp(cbcnor,n,bcnor,nh_d,nh_v,halo_v,dl,dzc,dzf,nor(:,:,:,2))
    call boundp(cbcnor,n,bcnor,nh_d,nh_v,halo_v,dl,dzc,dzf,nor(:,:,:,3))
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 1b.  time integration of c2^(n)--->c2^(n+1)

!#if defined(_EVAP)

!   sum_tmp = 0.0
!   sum_grid = 0.0

!   do i=1,n(1)
!      do j=1,n(2)
!        do k=1,n(3)
                    
!            if(psi(i,j,k)  .gt. 0.001 .and. psi(i,j,k) .lt. 0.009) then
!                    sum_tmp = sum_tmp + tmp(i,j,k)
!                    sum_grid = sum_grid + 1
!            endif
!       enddo
!      enddo
!    enddo
    
!     call mpi_allreduce(MPI_IN_PLACE,sum_tmp ,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
!     call mpi_allreduce(MPI_IN_PLACE,sum_grid ,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)

!     avg_tmp = sum_tmp/sum_grid


!     do i=1,n(1)
!      do j=1,n(2)
!        do k=1,n(3)

            !Anotine equation of state, coefficients are for water between 273 k -373 K
!            p_sat = (10**(8.07131-(1730.63/(233.426+avg_tmp-273.0))))*133.322
!            o_sat(i,j,k) =  p_sat/(p_sat+(p0-p_sat)*(mol_mas_air/mol_mas_vap))

!       enddo
!      enddo
!    enddo

!    call boundp(cbcsca,n,bcsca,nh_d,nh_v,halo_v,dl,dzc,dzf,o_sat)

!   call src_evap(n,nh_v,tmp,src_ph,src_mod,t_relax,psi,omega,rho_gas,m2,tmp_min,tmp_max,p0,c2,o_sat)
!#else
!    src_ph(:,:,:) = 0.0
!#endif
!    call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,src_ph)
!    call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,src_mod)


#if defined(_SCALAR_TRANS)
    call cmpt_rel_vel(n,dli,nh_d,nh_u,nh_v,u,v,w,c1,psi,nor,ur,vr,wr)     !change A inside the subroutine if you dont need attraction  
    call boundp(cbcrel,n,bcrel,nh_d,nh_v,halo_v,dl,dzc,dzf,ur)
    call boundp(cbcrel,n,bcrel,nh_d,nh_v,halo_v,dl,dzc,dzf,vr)
    call boundp(cbcrel,n,bcrel,nh_d,nh_v,halo_v,dl,dzc,dzf,wr)
    !call rk_c1(f_t1,f_t2,n,dli,nh_d,nh_u,nh_v, &
    !            dzc,dzf,diff1,diff2,u_reg,eps_int,u,v,w,ur,vr,wr,c1,c2,psi,nor,dc1dtrko,dc2dtrko,+1)!SCALAR TRANS. (Scalar inside the dispersed phase) 
    !call boundp(cbcsca,n,bcsca,nh_d,nh_v,halo_v,dl,dzc,dzf,c2)
    ! 
    !Update c1 at time-step "n+1"
    !
    call rk_c1(f_t1,f_t2,n,dli,nh_d,nh_u,nh_v, &
                dzc,dzf,diff1,diff2,u_reg,eps_int,u,v,w,ur,vr,wr,c1,c2,src_mod,psi,nor,dc1dtrko,dc2dtrko,-1)!SCALAR TRANS. (Scalar outside the dispersed phase) 
    ! 
    call boundp(cbcsca,n,bcsca,nh_d,nh_v,halo_v,dl,dzc,dzf,c2)    !NEEDS to be discussed 
    
    
!    call boundp(cbcsca,n,bcsca,nh_d,nh_v,halo_v,dl,dzc,dzf,c2)
    ! 
#endif
   !
   !
   psiold(:,:,:)    = psi(:,:,:)
   norold(:,:,:,:)  = nor(:,:,:,:)
   kappaold(:,:,:)  = kappa(:,:,:)
   

!for calculating dp0dt with the fields at nth time level  

#if defined(_LOW_MACH)
    call cmpt_rhs_low_mach(n,nh_v,dli,diff2,rho_gas,rho1,rho_air,rho_vap,cpp,cp1,cp_gas,cp_air,cp_vap,&
                          cv,cv1,cv_gas,cv_air,cv_vap,mol_mas,mol_mas_air,mol_mas_vap,&
                          psi2,omega,m2,src_ph,cond,tmp,u,v,w,dt,rhs_div_th,dp0dt,p0,h01,h02)
    call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,rhs_div_th)
    call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,conduction)

#endif


!   call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,rhs_div_th)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 2. Update the order-parameters at (psi, m2) to time-step (n+1)

   !
   call cmpt_umax(n,nh_u,u,v,w,gamma_v,u_reg)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 2a.incompressible field psi^(n) -------> psi^(n+1)

   call rk_psi(f_t1,f_t2,n,dli,nh_d,nh_u,nh_v, &
                dzc,dzf,u_reg,eps_int,u,v,w,psi,src_ph/rho1,ls,nor,kappa,dpsidtrko)

 
#if defined(_LOW_MACH)     
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 2b. compressible field  m2^(n)---------------->m2^(n+1)
 ! print*,'low mach on'
  
! do i = 1,n(1)
!     do j = 1,n(2)
!       do k = 1,n(3)

!         m2(i,j,k) = psi2(i,j,k)*rho_gas(i,j,k)


!        end do
!      end do
!    end do

!   call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,m2)


  call rk_m2(f_t1,f_t2,n,dli,nh_d,nh_u,nh_v, &
                dzc,dzf,u_reg,eps_int,u,v,w,m2,psi2,src_ph,ls,-nor,kappa,dm2dtrko,rho_gas)
 
  call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,m2)
  
   do i = 1,n(1)
     do j = 1,n(2)
       do k = 1,n(3)

         psi2(i,j,k) = 1.0-psi(i,j,k)
        

        end do
      end do
    end do
    
    call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,psi2)
#else
    do i = 1,n(1)
     do j = 1,n(2)
       do k = 1,n(3)
         
         psi2(i,j,k) = 1.0-psi(i,j,k)
         m2(i,j,k) = psi2(i,j,k)*rho_gas(i,j,k)

        end do
      end do
    end do
    ! print*,'low mach off'
    !call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,rho_gas)
#endif
   
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 2c. updating properties and calcultaing vapor mass fraction to n+1 with new fields
   ! 
   ! 
#if defined(_CONTACT_LINE) 
   call static_contact_angle_m(n,dl,nh_v,eps_int,theta,psi,psi_bound)
   call boundsb_array(n,psi_bound,nh_v,halo_v,psi)
#else  
   call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,psi)
   call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,m2)
   call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,psi2)
#endif
   !
#if defined(_USE_ACDI)
   call psi_to_ls(n,dl,nh_p,psi,eps_int,ls)
   call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,ls)
   call cmpt_norm(n,dli,nh_p,psi,kappa,nor)
#else
   call cmpt_norm(n,dli,nh_p,psi,kappa,nor)
   call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,kappa       )
   call boundp(cbcnor,n,bcnor,nh_d,nh_v,halo_v,dl,dzc,dzf,nor(:,:,:,1))
   call boundp(cbcnor,n,bcnor,nh_d,nh_v,halo_v,dl,dzc,dzf,nor(:,:,:,2))
   call boundp(cbcnor,n,bcnor,nh_d,nh_v,halo_v,dl,dzc,dzf,nor(:,:,:,3))
  !
#endif
#endif
   !
   rhoold(:,:,:)       = rho(:,:,:)
   muold(:,:,:)        = mu(:,:,:)
#if defined(_HEAT_TRANSFER)  
   cppold(:,:,:)       = cpp(:,:,:)
   condold(:,:,:)      = cond(:,:,:)
   rhocpold(:,:,:)     = rhocp(:,:,:)
#endif
   !
   !
   !Update properties at time-step (n+1)
   ! 


!#if defined(_EVAP)

!   sum_tmp = 0.0
!   sum_grid = 0.0

!    do i=1,n(1)
!      do j=1,n(2)
!        do k=1,n(3)

!            if(psi(i,j,k)  .gt. 0.001 .and. psi(i,j,k) .lt. 0.009) then
!                    sum_tmp = sum_tmp + tmp(i,j,k)
!                    sum_grid = sum_grid + 1
!            endif
!       enddo
!      enddo
!    enddo

!     call mpi_allreduce(MPI_IN_PLACE,sum_tmp ,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
!     call mpi_allreduce(MPI_IN_PLACE,sum_grid ,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)

!     avg_tmp = sum_tmp/sum_grid


#if defined(_EVAP)   
  
   do i = 1,n(1)
     do j = 1,n(2)
       do k = 1,n(3)
         
         omeg(i,j,k) = c2(i,j,k)/psi2(i,j,k)
         omega(i,j,k) =  c2(i,j,k)/psi2(i,j,k)
         
        ! if(omega(i,j,k) .gt. o_sat(i,j,k))then
         omega(i,j,k) =  min(omega(i,j,k),o_sat(i,j,k))
        ! endif
  
       enddo
     enddo
   enddo
   !
   call boundp(cbcsca,n,bcsca,nh_d,nh_v,halo_v,dl,dzc,dzf,omega)   !NEEDS to be discussed
   call boundp(cbcsca,n,bcsca,nh_d,nh_v,halo_v,dl,dzc,dzf,omeg)   
#else
    
    omega(:,:,:) = 0.0
    call boundp(cbcsca,n,bcsca,nh_d,nh_v,halo_v,dl,dzc,dzf,omega)

#endif 

   
   !
   call update_property(n,(/mu1 ,mu2 /),psi,mu,1)
   call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,mu)
   ! 

#if defined(_HEAT_TRANSFER)

#if defined (_EVAP)

     !thermal conductivity

     call update_property(n,(/cond_vap,cond_air/),omega,cond_gas,1)
     call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,cond_gas)

     call update_property_var(n,psi,cond1,cond_gas,cond)
     call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,cond)

     !cp

     call update_property(n,(/cp_vap,cp_air/),omega,cp_gas,1)
     call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,cp_gas)

     call update_property_var(n,psi,cp1,cp_gas,cpp)
     call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,cpp)

     !cv

     call update_property(n,(/cv_vap,cv_air/),omega,cv_gas,1)
     call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,cv_gas)

     call update_property_var(n,psi,cv1,cv_gas,cv)
     call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,cv)

#else

    !thermal conductivity

    call update_property(n,(/cond1,cond2/),psi,cond,1)
    call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,cond)
    !cp

    call update_property(n,(/cp1,cp2/),psi,cpp,1)
    call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,cpp)

    cp_gas(:,:,:) = min(cp1,cp2)
    call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,cp_gas)

    !cv

    call update_property(n,(/cv1,cv2/),psi,cv,1)
    call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,cv)
    
    cv_gas(:,:,:) = min(cv1,cv2)
    call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,cv_gas)

#endif

    !enthalpy of formation

!    call update_property(n,(/h01*rho1,h02*rho2/),psi,lh,1)
!    call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,lh)

    !rhocp calculation

#if defined(_LOW_MACH)

     do k=1,n(3)
       do j=1,n(2)
         do i=1,n(1)

           rhocp(i,j,k) = rho1*cp1*psi(i,j,k) + m2(i,j,k) * cp_gas(i,j,k)
           lh(i,j,k) = rho1*h01*psi(i,j,k) + m2(i,j,k)*h02
         end do
       end do
     end do
#else
     call update_property(n,(/rho1*cp1,rho2*cp2/),psi,rhocp,1)
     call update_property(n,(/h01*rho1,h02*rho2/),psi,lh,1)

#endif
     call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,rhocp)
     call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,lh)


        !heat transfer flag end
#endif

   !
   !
    call profiler_stop("DI") 
    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 3. calculation of regularization terms, time marching of enthalpy and getting back temperature

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 3a. calculation of regulaization terms 
#if defined(_HEAT_TRANSFER)
    gradpx = 0.d0
    psifx = 0.d0
    regx  = 0.d0

    gradpy = 0.d0
    psify = 0.d0
    regy  = 0.d0

    gradpz = 0.d0
    psifz = 0.d0
    regz  = 0.d0
   
    tmpx = 0.d0
    tmpy = 0.d0
    tmpz = 0.d0
  
    mfx(:,:,:) = 0.d0
    mfy(:,:,:) = 0.d0
    mfz(:,:,:) = 0.d0

    do k=1,n(3)
     do j=1,n(2)
      do i=1,n(1)

        ip = i + 1
        im = i - 1
        jp = j + 1
        jm = j - 1
        kp = k + 1
        km = k - 1
#if defined(_TWOD)
      regx   = 0.d0
#else
      gradpx =   ( psiold(ip,j,k) - psiold(i,j,k) ) * dli(1)
      psifx  =   ( psiold(i,j,k)*(1.d0-psiold(i,j,k))*norold(i,j,k,1)+psiold(ip,j,k)*(1.d0-psiold(ip,j,k))*norold(ip,j,k,1))*0.5
      regx   =   u_reg*( eps_int*gradpx - psifx)
      tmpx   =   (tmp(i,j,k)+tmp(ip,j,k))*0.5d0
#endif
      gradpy =   ( psiold(i,jp,k) - psiold(i,j,k) ) * dli(2)
      psify  =   ( psiold(i,j,k)*(1.d0-psiold(i,j,k))*norold(i,j,k,2)+psiold(i,jp,k)*(1.d0-psiold(i,jp,k))*norold(i,jp,k,2))*0.5
      regy   =   u_reg*( eps_int*gradpy - psify)
      tmpy   =   (tmp(i,jp,k)+tmp(i,j,k))*0.5d0

      gradpz =   ( psiold(i,j,kp) - psiold(i,j,k) ) * dli(3)
      psifz  =   ( psiold(i,j,k)*(1.d0-psiold(i,j,k))*norold(i,j,k,3)+psiold(i,j,kp)*(1.d0-psiold(i,j,kp))*norold(i,j,kp,3))*0.5
      regz   =   u_reg*( eps_int*gradpz - psifz)
      tmpz   =   (tmp(i,j,kp)+tmp(i,j,k))*0.5d0
      !
#if defined(_REGU)    
       
       !the formula is (rho1 - rho2), now rho_gas = min(rho1,rho2)
       !if rho1 is li1, rho1 is gt than rho2; so sho2 = rho_gas, we get back rho1-rho2


      if (rho1 .gt. rho2) then
      
       mfx(i,j,k)  = (rho1*(cp1*tmpx+h01)-rho_gas(i,j,k)*(cp2*tmpx+h02))*regx
       mfy(i,j,k)  = (rho1*(cp1*tmpy+h01)-rho_gas(i,j,k)*(cp2*tmpy+h02))*regy
       mfz(i,j,k)  = (rho1*(cp1*tmpz+h01)-rho_gas(i,j,k)*(cp2*tmpz+h02))*regz
      
      !Regularization is zero in case of density ratio 1

       else if(rho1 .eq. rho2) then
       
       mfx(i,j,k)  = 0.0
       mfy(i,j,k)  = 0.0
       mfz(i,j,k)  = 0.0
      
      !if rho1 is less than rho2; then rho1 = rho_gas, and we get back rho1-rho2

       else

       mfx(i,j,k)  = (rho_gas(i,j,k)*(cp1*tmpx+h01)-rho2*(cp2*tmpx+h02))*regx
       mfy(i,j,k)  = (rho_gas(i,j,k)*(cp1*tmpy+h01)-rho2*(cp2*tmpy+h02))*regy
       mfz(i,j,k)  = (rho_gas(i,j,k)*(cp1*tmpz+h01)-rho2*(cp2*tmpz+h02))*regz

       endif


#else
      ! print*, 'regu off enthalpy' 
       mfx(i,j,k)  = 0.0
       mfy(i,j,k)  = 0.0
       mfz(i,j,k)  = 0.0
#endif
      !
      enddo
     enddo
    enddo
    !
    !
    call bounduvw(cbcvel,n,bcvel,nh_d,nh_u,halo_u,no_outflow,dl,dzc,dzf,mfx,mfy,mfz)
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   3b. enthalpy advancement in time rho*h^(n)----------->rho*h^(n+1)
    
    !We pass rho*h and take back rho*h, so we need to apply boundary condition on (rho*h)
   

    call rk_ent(f_t1,f_t2,n,dli,nh_d,nh_u,nh_v,dzc,dzf,eps_int,u_reg, &
               u,v,w,psiold,norold,tmp,condold,mfx,mfy,mfz,psi2,dp0dt,ent,dentdtrko)
    !
     call boundp_di(cbctmp,n,bctmp,nh_d,nh_p,halo_p,dl,dzc,dzf,ent,rhocp,lh,tmp) ! b.c imposed on enthalpy
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   3c. recalculating back the temperature
 !   call boundp(cbctmp,n,bctmp,nh_d,nh_p,halo_p,dl,dzc,dzf,ent) 

    do i=1,n(1)
      do j=1,n(2)
       do k=1,n(3)

       tmp(i,j,k) = (ent(i,j,k)-lh(i,j,k))/rhocp(i,j,k)

       enddo
      enddo
     enddo


#else
    
    tmp(:,:,:) = tl0
    conduction(:,:,:) = 0.0

!   call boundp_di(cbctmp,n,bctmp,nh_d,nh_p,halo_p,dl,dzc,dzf,ent,rhocp,lh,tmp) ! b.c imposed on enthalpy
!
#endif

    call boundp(cbctmp,n,bctmp,nh_d,nh_p,halo_p,dl,dzc,dzf,conduction)
    call boundp(cbctmp,n,bctmp,nh_d,nh_p,halo_p,dl,dzc,dzf,tmp)

!updating some properties based on omega an m2

#if defined(_EVAP)
       
   do i = 1,n(1)
     do j = 1,n(2)
       do k = 1,n(3)
           
               !mol_mas (i,j,k) = 1.0 / (omega(i,j,k) / mol_mas_vap + &
               !                  (1.0-omega(i,j,k)) / mol_mas_air)
               mol_mas(i,j,k) = psi2(i,j,k)/&
                                (c2(i,j,k)*(1/mol_mas_vap - 1/mol_mas_air) +&
                                 psi2(i,j,k)/mol_mas_air + 1e-30)
               
         enddo
     enddo
   enddo
  
  call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,mol_mas)

#else
   
   mol_mas(:,:,:) = mol_mas_air
  
#endif

#if defined(_LOW_MACH)

    mass_gas  =  sum(m2(1:n1,1:n2,1:n3))*dl(1)*dl(2)*dl(3)
    call mpi_allreduce(MPI_IN_PLACE, mass_gas ,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)

    sum_rt = 0    

    do i = 1,n(1)
     do j = 1,n(2)
       do k = 1,n(3)
              
               
               sum_rt = sum_rt + &
                       ((1.0-psi(i,j,k))/(tmp(i,j,k)*(8314.0/mol_mas(i,j,k))))*dl(1)*dl(2)*dl(3)

      enddo
     enddo
   enddo

   call mpi_allreduce(MPI_IN_PLACE,sum_rt ,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)

 !  p0 = mass_gas/sum_rt

   do i = 1,n(1)
     do j = 1,n(2)
       do k = 1,n(3)

         rho(i,j,k) = psi(i,j,k)*rho1 + m2(i,j,k)
         rho_gas(i,j,k) = p0/(tmp(i,j,k)*(8314.0/mol_mas(i,j,k)))

         enddo
     enddo
   enddo
   
   call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,rho)
   call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,rho_gas)

#else
    
    do i = 1,n(1)
      do j = 1,n(2)
       do k = 1,n(3)

         rho_gas(i,j,k) = min(rho1,rho2)

         enddo
      enddo
    enddo

    call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,rho_gas)

    call update_property(n,(/rho1,rho2/),psi,rho,1)
    call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,rho)

#endif



      !!!!!We calculate p0 and dp0dt exactly after enthalpy advection - following Dalla-Barba and Scapin

#if defined(_LOW_MACH)
    call cmpt_rhs_low_mach(n,nh_v,dli,diff2,rho_gas,rho1,rho_air,rho_vap,cpp,cp1,cp_gas,cp_air,cp_vap,&
                          cv,cv1,cv_gas,cv_air,cv_vap,mol_mas,mol_mas_air,mol_mas_vap,&
                          psi2,omega,m2,src_ph,cond,tmp,u,v,w,dt,rhs_div_th,dp0dt,p0,h01,h02) 

    call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,rhs_div_th)               
#else
   call cmpt_rhs_incompressible(n,nh_v,dli,rho1,rho_air,rho_vap,omeg,m2,c2,diff2,src_ph,psi,rhs_div_th)
   call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,rhs_div_th)
#endif

!    call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,rhs_div_th)
     

#if defined(_EVAP)

    tmp_min = minval(tmp(1:n1,1:n2,1:n3))
    call mpi_allreduce(MPI_IN_PLACE,tmp_min,1,MPI_REAL_RP,MPI_MIN,MPI_COMM_WORLD,ierr)
    tmp_max = maxval(tmp(1:n1,1:n2,1:n3))
    call mpi_allreduce(MPI_IN_PLACE,tmp_max,1,MPI_REAL_RP,MPI_MIN,MPI_COMM_WORLD,ierr)




   sum_tmp = 0.0
   sum_grid = 0.0

    do i=1,n(1)
      do j=1,n(2)
        do k=1,n(3)

            if(psi(i,j,k)  .lt. 0.999 .and. psi(i,j,k) .gt. 0.001) then
                    sum_tmp = sum_tmp + tmp(i,j,k)
                    sum_grid = sum_grid + 1
            endif
       enddo
      enddo
    enddo

     call mpi_allreduce(MPI_IN_PLACE,sum_tmp ,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
     call mpi_allreduce(MPI_IN_PLACE,sum_grid ,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)

     avg_tmp = sum_tmp/sum_grid


     do i=1,n(1)
      do j=1,n(2)
        do k=1,n(3)

            !Anotine equation of state, coefficients are for water between 273 k -373 K
            p_sat = (10**(8.07131-(1730.63/(233.426+tmp(i,j,k)-273.0))))*133.322
            o_sat(i,j,k) =  p_sat/(p_sat+(p0-p_sat)*(mol_mas_air/mol_mas_vap))

       enddo
      enddo
    enddo

    call boundp(cbcsca,n,bcsca,nh_d,nh_v,halo_v,dl,dzc,dzf,o_sat)

     do i=1,n(1)
      do j=1,n(2)
        do k=1,n(3) 
           
      !      if(psi(i,j,k) .gt. 0.9) then

               omega(i,j,k) = min(omega(i,j,k),o_sat(i,j,k))
  !          endif  
        enddo
      enddo
    enddo

    call boundp(cbcsca,n,bcsca,nh_d,nh_v,halo_v,dl,dzc,dzf,omega)


    call src_evap(n,nh_v,tmp,src_ph,src_mod,t_relax,psi,omega,rho_gas,m2,tmp_min,tmp_max,p0,c2,o_sat)
#else
    src_ph(:,:,:) = 0.0
#endif
    call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,src_ph)
   call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,src_mod)


    !
#if defined(_NAVIER_STOKES)
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 4. two-fluid Navier-Stokes --> (u,v,w,p)^(n+1)
    !
    rho0i = 1._rp/rho0
    !
    !
    call profiler_start("RK", tag = .true., tag_color = COLOR_RED)
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 4a. Compute mass fluxes (regularization) at the faces which will be used in the momentum equation 
    !
    gradpx = 0.d0
    psifx = 0.d0
    regx  = 0.d0 

    gradpy = 0.d0
    psify = 0.d0
    regy  = 0.d0 

    gradpz = 0.d0
    psifz = 0.d0
    regz  = 0.d0
    mfx(:,:,:) = 0.d0
    mfy(:,:,:) = 0.d0
    mfz(:,:,:) = 0.d0

    do k=1,n(3)
     do j=1,n(2)
      do i=1,n(1)

        ip = i + 1
        im = i - 1
        jp = j + 1
        jm = j - 1
        kp = k + 1
        km = k - 1
#if defined(_TWOD)        
      regx   = 0.d0
#else
      gradpx =   ( psiold(ip,j,k) - psiold(i,j,k) ) * dli(1)
      psifx  =   ( psiold(i,j,k)*(1.d0-psiold(i,j,k))*norold(i,j,k,1)+psiold(ip,j,k)*(1.d0-psiold(ip,j,k))*norold(ip,j,k,1))*0.5
      regx   =   u_reg*( eps_int*gradpx - psifx) 
#endif
      gradpy =   ( psiold(i,jp,k) - psiold(i,j,k) ) * dli(2)
      psify  =   ( psiold(i,j,k)*(1.d0-psiold(i,j,k))*norold(i,j,k,2)+psiold(i,jp,k)*(1.d0-psiold(i,jp,k))*norold(i,jp,k,2))*0.5
      regy   =   u_reg*( eps_int*gradpy - psify) 

      gradpz =   ( psiold(i,j,kp) - psiold(i,j,k) ) * dli(3)
      psifz  =   ( psiold(i,j,k)*(1.d0-psiold(i,j,k))*norold(i,j,k,3)+psiold(i,j,kp)*(1.d0-psiold(i,j,kp))*norold(i,j,kp,3))*0.5
      regz   =   u_reg*( eps_int*gradpz - psifz) 
      !
#if defined(_REGU)
      
      if(rho1 .gt. rho2)then
    !  print*, 'regu on NSE' 
        mfx(i,j,k) = 0.5d0*(rhoold(ip,j,k)+rhoold(i,j,k))*u(i,j,k) - (rho1-rho_gas(i,j,k))*regx
        mfy(i,j,k) = 0.5d0*(rhoold(i,jp,k)+rhoold(i,j,k))*v(i,j,k) - (rho1-rho_gas(i,j,k))*regy
        mfz(i,j,k) = 0.5d0*(rhoold(i,j,kp)+rhoold(i,j,k))*w(i,j,k) - (rho1-rho_gas(i,j,k))*regz
      
      elseif(rho1 .eq. rho2) then
      
        mfx(i,j,k) = 0.0
        mfy(i,j,k) = 0.0
        mfz(i,j,k) = 0.0

      else
      
        mfx(i,j,k) = 0.5d0*(rhoold(ip,j,k)+rhoold(i,j,k))*u(i,j,k) - (rho_gas(i,j,k) - rho2)*regx
        mfy(i,j,k) = 0.5d0*(rhoold(i,jp,k)+rhoold(i,j,k))*v(i,j,k) - (rho_gas(i,j,k) - rho2)*regy
        mfz(i,j,k) = 0.5d0*(rhoold(i,j,kp)+rhoold(i,j,k))*w(i,j,k) - (rho_gas(i,j,k) - rho2)*regz

      endif
#else     
    !  print*, 'regu off NSE'   
       
        mfx(i,j,k) = 0.0
        mfy(i,j,k) = 0.0
        mfz(i,j,k) = 0.0
#endif      
      enddo
     enddo
    enddo
    !
    call bounduvw(cbcvel,n,bcvel,nh_d,nh_u,halo_u,no_outflow,dl,dzc,dzf,mfx,mfy,mfz)
    

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 4b. predictor equations - time integrations 

    !input of the prediction is velocities but the output is product of rho.u
    !
    !old density and viscoity are being used for diffusion and advective terms
    !
    call rk(space_scheme_mom,f_t1,f_t2,n(1),n(2),n(3),dxi,dyi,dzi,nh_d,nh_u,dzci,dzfi, &
            u,v,w,mfx,mfy,mfz,muold,rhoold,dudtrko,dvdtrko,dwdtrko)
    !
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 4c. add the source terms
    !
#if defined(_TWO_PHASE)    
    ! --> add gravity terms
    !
    !updated density will be used for surface forces and gravitational force 
    !
    call grav_tw_src(n(1),n(2),n(3),rho,dt,cbcpre,dxi,dyi,dzi,nh_d,dzfi,nh_u,u,v,w)  
    !
    ! --> add the surface tension forces
    !
    call surft_src(n(1),n(2),n(3),nh_d,nh_u,f_t12,dxi,dyi,dzi,dzci,kappa,psi,rho,u,v,w)
    !
#endif
!
!Marangoni is not yet validated 
!
!#if defined(_USE_MAR)
!    call cmpt_delta(n,dli,psi,delta)
!    call boundp(cbcpsi,n,bcpsi,nh_d,nh_p,halo_p,dl,dzc,dzf,delta)
!    !call mar_src(1./dzf,n(1),n(2),n(3),nh_d,nh_u,nh_t,f_t12,dxi,dyi,dzi,dzci,delta,nor,kappa,a1,tmp,rho,u,v,w)!(c1+1.d0)*300.d0,rho,u,v,w)
!    call mar_src(1./dzf,n(1),n(2),n(3),nh_d,nh_u,nh_v,f_t12,dxi,dyi,dzi,dzci,delta,nor,kappa,psi,(c1+1.d0)*300.d0,rho,u,v,w)
!#endif
!
!
!!Forcing is not yet validated 
!
!
#if defined(_TURB_FORCING)
    !
    ! --> add the forcing to sustain turbulence (triperiodic cases)
    !
    call forc_src(n(1),n(2),n(3),nh_d,nh_u,f_t12,dx,dy,dz,zc,rho,u,v,w)
    !
#endif
    !
    ! --> add the bulk velocity forcing 
    !     note: compute and add the bulk velocity forcing at the end so that the 
    !           f(1:3) accounts for all the terms in the prediction
    !
!    if(any(is_forced(:))) then
!      call bulk_forcing_src(bulk_ftype,is_forced,n(1),n(2),n(3),dx,dy,dz,f_t12, &
!                            nh_d,nh_u,dzc,dzf,rho,u,v,w,f)
!      dpdl_c(1:3) = dpdl_c(1:3) + f(1:3) 
!    endif
!
!
!Extracting velocity by dividing to updated density 
!
    do k=1,n(3)
     do j=1,n(2)
      do i=1,n(1)
         
       rhox = 0.5d0*(rho(i,j,k)+rho(i+1,j,k))
       rhoy = 0.5d0*(rho(i,j,k)+rho(i,j+1,k))
       rhoz = 0.5d0*(rho(i,j,k)+rho(i,j,k+1))

       u(i,j,k) = u(i,j,k)/rhox                !!!!!!!!Extracting velocity by dividing with the updated density 
       v(i,j,k) = v(i,j,k)/rhoy
       w(i,j,k) = w(i,j,k)/rhoz

    enddo
   enddo
  enddo
     !
     !impose the pressure gradient
     !
     !--> updated density will be used to compute add the pressure gradient term 
     !
#if defined(_TWO_PHASE)
    call pres_tw_src(n(1),n(2),n(3),dxi,dyi,dzi,nh_d,nh_u,dzci,rho0i,f_t12,f_t12_o,p,pold,rho,u,v,w)
#endif
    !
    !
    call bounduvw(cbcvel,n,bcvel,nh_d,nh_u,halo_u,no_outflow,dl,dzc,dzf,u,v,w) ! we impose bc at end (not valid for all cases)
    !
    call profiler_stop("RK")
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 4d. construct the Poisson equation - calculating the right hand side
    !
#if defined(_OPENACC)
    !$acc parallel loop collapse(3) present(p, pold)
#else
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP SHARED(p,pold)
#endif
    do k=1,n3
      do j=1,n2
        do i=1,n1
          pold(i,j,k) = p(i,j,k)
        enddo
      enddo
    enddo
#if defined(_OPENACC)
    !$acc end parallel loop 
#else
    !$OMP END PARALLEL DO
#endif
    call boundp(cbcpre,n,bcpre,nh_d,nh_p,halo_p,dl,dzc,dzf,pold)
    !
    !
    !Update the saturation pressure and saturation vapor mass at time-step "n+1"
    !
    tmp_min = minval(tmp(1:n1,1:n2,1:n3))
    call mpi_allreduce(MPI_IN_PLACE,tmp_min,1,MPI_REAL_RP,MPI_MIN,MPI_COMM_WORLD,ierr)
    tmp_max = maxval(tmp(1:n1,1:n2,1:n3))
    call mpi_allreduce(MPI_IN_PLACE,tmp_max,1,MPI_REAL_RP,MPI_MIN,MPI_COMM_WORLD,ierr)
    
    
    !
    !!!!calculating the source term for evaporation
  !  call src_evap(n,nh_v,tmp,src_ph,t_relax,psi,omega,rho_gas,m2,tmp_min,tmp_max,p0,c2)
   ! call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,src_ph)
    !
    !
    !The source-term at time-step (n+1) should be added to the continuity to avoid time-lag
    !
    !!for computation of z, RHS of div u in poisson solver, non zero value for phase change, different for incompressible and low mach
!conduction(:,:,:) = 0.0
#if defined(_LOW_MACH)
    call cmpt_rhs_low_mach(n,nh_v,dli,diff2,rho_gas,rho1,rho_air,rho_vap,cpp,cp1,cp_gas,cp_air,cp_vap,&
                          cv,cv1,cv_gas,cv_air,cv_vap,mol_mas,mol_mas_air,mol_mas_vap,&
                          psi2,omega,m2,src_ph,cond,tmp,u,v,w,dt,rhs_div_th,dp0dt,p0,h01,h02)
    call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,rhs_div_th)


#else
   call cmpt_rhs_incompressible(n,nh_v,dli,rho1,rho_air,rho_vap,omega,m2,c2,diff2,src_ph,psi,rhs_div_th)
   call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,rhs_div_th)
#endif


   call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,rhs_div_th)    

   
#if defined(_LOW_MACH)

   call fillps(n(1),n(2),n(3),nh_d,nh_u,dxi,dyi,dzi,dzfi,f_t12_i,rho0,u,v,w,p,rhs_div_th*1.0)

#else
   
   call fillps(n(1),n(2),n(3),nh_d,nh_u,dxi,dyi,dzi,dzfi,f_t12_i,rho0,u,v,w,p,rhs_div_th*1.0)

#endif

    call boundp(cbcpre,n,bcpre,nh_d,nh_p,halo_p,dl,dzc,dzf,p)
    !
    call updt_rhs_b(n(1),n(2),n(3),(/'c','c','c'/),cbcpre,nh_p,rhsbp_x,rhsbp_y,rhsbp_z,p)
    !
    call profiler_start("SOLVER", tag = .true., tag_color = COLOR_GREEN)
    !
#if defined(_OPENACC)
    call solver_gpu(n_z,dims_xyz(:,3),arrplanp,normfftp,lambdaxyp,ap,bp,cp,cbcpre(:,:),(/'c','c','c'/),p)
#else
    call solver_cpu(n,arrplanp,normfftp,lambdaxyp,ap,bp,cp,cbcpre(:,3),(/'c','c','c'/),p)
#endif
    call boundp(cbcpre,n,bcpre,nh_d,nh_p,halo_p,dl,dzc,dzf,p)
    !
    call profiler_stop("SOLVER")
    !
    call profiler_start("CORREC")
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 4e. correct the velocity and update the pressure
    !
    !u,v,w are the product of velocity and denisty as an input and output 
    !
    call correc(n(1),n(2),n(3),nh_d,nh_u,dxi,dyi,dzi,dzci,f_t12,rho0,p,u,v,w,rho)
    call bounduvw(cbcvel,n,bcvel,nh_d,nh_u,halo_u,is_outflow,dl,dzc,dzf,u,v,w)
    !
  
   

#if defined(_OPENACC)
    !$acc parallel loop collapse(3)
#else
    !$OMP WORKSHARE
#endif    
    do k=1,n3
      do j=1,n2
        do i=1,n1
          p(i,j,k) = pold(i,j,k) + p(i,j,k)
        enddo
      enddo
    enddo
#if defined(_OPENACC)
    !$acc end parallel loop
#else
    !$OMP END WORKSHARE
#endif
    !
    call boundp(cbcpre,n,bcpre,nh_d,nh_p,halo_p,dl,dzc,dzf,p)
    !
    call profiler_stop("CORREC")
    ! 
#endif
    ! 
    write(fldnum,'(i9.9)') istep
    !
    if(.not.late_init.and.istep.ge.i_late_init) then
      include 'dropcheck.h90'
    endif
    !
#if defined(_DO_POSTPROC) && defined(_TURB_FORCING)
    call budget(datadir,n(1),n(2),n(3),ng(1),ng(2),ng(3),rho1,rho2,dx,dy,dz,nh_u, &
                u,v,w,rho,psi,time,istep)
#endif
    !
    ! check simulation stopping criteria
    !
    if(stop_type(1)) then ! maximum number of time steps reached
      if(istep.ge.nstep   ) is_done = is_done.or..true.
    endif
    if(stop_type(2)) then ! maximum simulation time reached
      if(time .ge.time_max) is_done = is_done.or..true.
    endif
    if(stop_type(3)) then ! maximum wall-clock time reached
      tw = (MPI_WTIME()-twi)/3600._rp
      if(tw   .ge.tw_max  ) is_done = is_done.or..true.
    endif
    !
    ! check time-step size and velocity divergence
    !
    dto = dt  
    if(mod(istep,icheck).eq.0) then
      !
      if(myid.eq.0) print*, 'Checking stability and divergence...'
      !
      if(.not.constant_dt) then
        if(myid.eq.0) print*, 'updating the time-step size ...'
#if defined(_TWO_PHASE)        
        call chkdt_tw(n(1),n(2),n(3),dxi,dyi,dzi,nh_d,nh_u,dzci,dzfi,u,v,w,dtmax)
#else
        call chkdt_sp(n(1),n(2),n(3),dxi,dyi,dzi,nh_d,nh_u,dzci,dzfi,u,v,w,dtmax)
#endif        
        dt = cfl*dtmax
      else
        dtmax = dt_input
        dt    = dtmax
      endif
      !
      if(myid.eq.0) print*, 'dtmax = ', dtmax, 'dt = ',dt
      if(dtmax.lt.small) then
        if(myid.eq.0) print*, 'ERROR: timestep is too small.'
        if(myid.eq.0) print*, 'Aborting...'
        is_done = .true.
        kill = .true.
      endif
      dti = 1._rp/dt
      !
      if(myid.eq.0) print*, 'checking the velocity divergence ...'
      call chkdiv(n(1),n(2),n(3),dxi,dyi,dzi,nh_d,nh_u,dzfi,u,v,w,divtot,divmax)
    !  if(divmax.gt.small.or.divtot.ne.divtot) then
    !    if(myid.eq.0) print*, 'ERROR: maximum divergence is too large.'
    !    if(myid.eq.0) print*, 'Aborting...'
    !    is_done = .true.
    !    kill = .true.
    !  endif
    !  !
    endif
    !
    call profiler_stop("STEP")
    !
    ! output routines below
    !
    if(mod(istep,iout0d).eq.0) then
      !
      call profiler_start("OUT:iout0d", tag = .true., tag_color = COLOR_WHITE)
      !
      var(:) = 0._rp
      var(1) = 1._rp*istep
      var(2) = dt
      var(3) = time
      call out0d(trim(datadir)//'time.out',3,var)
      !
      
      
      psi_min = minval(psi(1:n(1),1:n(2),1:n(3)))
      call mpi_allreduce(MPI_IN_PLACE,psi_min,1,MPI_REAL_RP,MPI_MIN,MPI_COMM_WORLD,ierr)
      psi_max = maxval(psi(1:n(1),1:n(2),1:n(3)))
      call mpi_allreduce(MPI_IN_PLACE,psi_max,1,MPI_REAL_RP,MPI_MAX,MPI_COMM_WORLD,ierr)
      vol_p1  = sum(psi(1:n(1),1:n(2),1:n(3)))*dl(1)*dl(2)*dl(3)
      call mpi_allreduce(MPI_IN_PLACE,vol_p1 ,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
#if defined(_SCALAR_TRANS)      
      vol_c1  = sum(src_ph(1:n(1),1:n(2),1:n(3)))*dl(1)*dl(2)*dl(3)
      vol_c2  = sum(src_mod(1:n(1),1:n(2),1:n(3)))*dl(1)*dl(2)*dl(3)
      vol_ct  = sum(ent(1:n(1),1:n(2),1:n(3)))*dl(1)*dl(2)*dl(3)
      call mpi_allreduce(MPI_IN_PLACE,vol_c1 ,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,vol_c2 ,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,vol_ct ,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      
      omega_min = minval(omega(1:n(1),1:n(2),1:n(3)))
      call mpi_allreduce(MPI_IN_PLACE,omega_min,1,MPI_REAL_RP,MPI_MIN,MPI_COMM_WORLD,ierr)
      omega_max = maxval(omega(1:n(1),1:n(2),1:n(3)))
      call mpi_allreduce(MPI_IN_PLACE,omega_max,1,MPI_REAL_RP,MPI_MAX,MPI_COMM_WORLD,ierr)
      
      !
#endif
      !
      vol_gas  = sum(psi(1:n1,1:n2,1:n3)*rho_gas(1:n1,1:n2,1:n3))*dl(1)*dl(2)*dl(3)
      call mpi_allreduce(MPI_IN_PLACE,vol_gas ,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)


      v_sum =  sum(w(1:n1,1:n2,1:n3)*psi(1:n1,1:n2,1:n3)*rho_gas(1:n1,1:n2,1:n3))*dl(1)*dl(2)*dl(3)
      call mpi_allreduce(MPI_IN_PLACE,v_sum ,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)

      vel_rise = v_sum/vol_gas

      
      var(:) = 0._rp
      var(1) = 1._rp*istep
      var(2) = dt
      var(3) = time
      var(4) = vel_rise                  !mass of vapor
      var(5) = vol_p1*rho1+vol_c2       !mass of liquid + mass of vapor
      var(6) = vol_p1
      var(7) = vol_c1                  !total source term
      var(8) = vol_c2                  !total modified source term
      var(9) = p0              !max value of vapor mass fraction
      var(10) = avg_tmp
      call out0d(trim(datadir)//'psi_info.out',10,var)
      !
      if(any(is_forced(:))) then
        call cmpt_mean(n(1),n(2),n(3),nh_d,nh_u,dx,dy,dzf,lx,ly,lz,u,meanvelu)
        call cmpt_mean(n(1),n(2),n(3),nh_d,nh_u,dx,dy,dzf,lx,ly,lz,v,meanvelv)
        call cmpt_mean(n(1),n(2),n(3),nh_d,nh_u,dx,dy,dzf,lx,ly,lz,w,meanvelw)
        var(:)   = 0._rp
        var(1)   = time
        var(2)   = 1._rp*istep
        var(3:5) = -dpdl_c(1:3)*dti*rho2
        var(6)   = sqrt(maxval(abs(var(3:5)))*lz*0.5_rp)
        var(7)   = rho2*var(6)*0.5_rp*lz/mu2 ! Re_tau
        var(8)   = meanvelu
        var(9)   = meanvelv
        var(10)  = meanvelw
        call out0d(trim(datadir)//'forcing.out',10,var)
      endif 
      !
      call profiler_stop("OUT:iout0d")
      !
    endif
    !
#if defined(_DO_POSTPROC)
    if(mod(istep,iout0d_ta).eq.0.and.do_tagging) then
      call droplet_tagging(n,dims,datadir_ta,dl,nh_d,nh_v,nh_u,halo_v,dzc,dzf,psi,u,v,w,istep,time)
    endif
#endif
    !
    if(mod(istep,iout1d).eq.0) then
      call profiler_start("OUT:iout1d", tag = .true., tag_color = COLOR_WHITE)
      ! [TODO] Prefetching GPU->CPU (and viceversa)
      include 'out1d.h90'
      call profiler_stop("OUT:iout1d")
    endif
    if(mod(istep,iout2d).eq.0) then
      call profiler_start("OUT:iout2d", tag = .true., tag_color = COLOR_WHITE)
      ! [TODO] Prefetching GPU->CPU (and viceversa)
      include 'out2d.h90'
      call profiler_stop("OUT:iout2d")
    endif
    if(mod(istep,iout3d).eq.0) then
      call profiler_start("OUT:iout3d", tag = .true., tag_color = COLOR_WHITE)
      ! [TODO] Prefetching GPU->CPU (and viceversa)
      include 'out3d.h90'
      call profiler_stop("OUT:iout3d")
    endif
    !
#if !defined(_BENCHMARK_NO_IO)
    if ( (mod(istep,isave) .eq. 0) .or. is_done .and. (.not.kill) ) then
      !
      call profiler_start("OUT:isave", tag = .true., tag_color = COLOR_WHITE)
      !
      action_load = 'w'
      !
      !@cuf istat = cudaMemPrefetchAsync(u, size(u), cudaCpuDeviceId, 0)
      !@cuf istat = cudaMemPrefetchAsync(v, size(v), cudaCpuDeviceId, 0)
      !@cuf istat = cudaMemPrefetchAsync(w, size(w), cudaCpuDeviceId, 0)
      !@cuf istat = cudaMemPrefetchAsync(dudtrko, size(dudtrko), cudaCpuDeviceId, 0)
      !@cuf istat = cudaMemPrefetchAsync(dvdtrko, size(dvdtrko), cudaCpuDeviceId, 0)
      !@cuf istat = cudaMemPrefetchAsync(dwdtrko, size(dwdtrko), cudaCpuDeviceId, 0)
      !@cuf istat = cudaMemPrefetchAsync(p, size(p), cudaCpuDeviceId, 0)
      !@cuf istat = cudaMemPrefetchAsync(psi, size(psi), cudaCpuDeviceId, 0)
      !@cuf istat = cudaMemPrefetchAsync(dpsidtrko, size(dpsidtrko), cudaCpuDeviceId, 0)
      !
      inquire(file='data/fldu.bin', exist=is_data)
      if(myid.eq.0.and.is_data) call execute_command_line('mv data/fldu.bin   data/fldu_old.bin')
      call load(action_load,trim(restart_dir)//'fldu.bin',n,      u(1:n(1),1:n(2),1:n(3)))
      !@cuf if ( (mod(istep,isave)).eq. 0 .and. (.not. is_done) .and. (.not.kill) ) istat = cudaMemPrefetchAsync(u, size(u), mydev, 0)
      !
      inquire(file='data/fldv.bin', exist=is_data)
      if(myid.eq.0.and.is_data) call execute_command_line('mv data/fldv.bin   data/fldv_old.bin')
      call load(action_load,trim(restart_dir)//'fldv.bin',n,      v(1:n(1),1:n(2),1:n(3)))
      !@cuf if ( (mod(istep,isave)).eq. 0 .and. (.not. is_done) .and. (.not.kill) ) istat = cudaMemPrefetchAsync(v, size(v), mydev, 0)
      !
      inquire(file='data/fldw.bin', exist=is_data)
      if(myid.eq.0.and.is_data) call execute_command_line('mv data/fldw.bin   data/fldw_old.bin')
      call load(action_load,trim(restart_dir)//'fldw.bin',n,      w(1:n(1),1:n(2),1:n(3)))   
      !@cuf if ( (mod(istep,isave)).eq. 0 .and. (.not. is_done) .and. (.not.kill) ) istat = cudaMemPrefetchAsync(w, size(w), mydev, 0)
      !
      inquire(file='data/flddu.bin', exist=is_data)
      if(myid.eq.0.and.is_data) call execute_command_line('mv data/flddu.bin  data/flddu_old.bin')
      call load(action_load,trim(restart_dir)//'flddu.bin',n,dudtrko(1:n(1),1:n(2),1:n(3)))
      !@cuf if ( (mod(istep,isave)).eq. 0 .and. (.not. is_done) .and. (.not.kill) ) istat = cudaMemPrefetchAsync(dudtrko, size(dudtrko), mydev, 0)
      !
      inquire(file='data/flddv.bin', exist=is_data)
      if(myid.eq.0.and.is_data) call execute_command_line('mv data/flddv.bin  data/flddv_old.bin')
      call load(action_load,trim(restart_dir)//'flddv.bin',n,dvdtrko(1:n(1),1:n(2),1:n(3)))
      !@cuf if ( (mod(istep,isave)).eq. 0 .and. (.not. is_done) .and. (.not.kill) ) istat = cudaMemPrefetchAsync(dvdtrko, size(dvdtrko), mydev, 0)
      !
      inquire(file='data/flddw.bin', exist=is_data)
      if(myid.eq.0.and.is_data) call execute_command_line('mv data/flddw.bin  data/flddw_old.bin')
      call load(action_load,trim(restart_dir)//'flddw.bin',n,dwdtrko(1:n(1),1:n(2),1:n(3)))
      !@cuf if ( (mod(istep,isave)).eq. 0 .and. (.not. is_done) .and. (.not.kill) ) istat = cudaMemPrefetchAsync(dwdtrko, size(dwdtrko), mydev, 0)
      !
      inquire(file='data/flddp.bin', exist=is_data)
      if(myid.eq.0.and.is_data) call execute_command_line('mv data/fldp.bin   data/fldp_old.bin')
      call load(action_load,trim(restart_dir)//'fldp.bin',n,      p(1:n(1),1:n(2),1:n(3)))
      !@cuf if ( (mod(istep,isave)).eq. 0 .and. (.not. is_done) .and. (.not.kill) ) istat = cudaMemPrefetchAsync(p, size(p), mydev, 0)
      !
#if defined(_TWO_PHASE)      
      inquire(file='data/fldpsi.bin', exist=is_data)
      if(myid.eq.0.and.is_data) call execute_command_line('mv data/fldpsi.bin data/fldpsi_old.bin')
      call load(action_load,trim(restart_dir)//'fldpsi.bin',n,    psi(1:n(1),1:n(2),1:n(3)))
      !@cuf if ( (mod(istep,isave)).eq. 0 .and. (.not. is_done) .and. (.not.kill) ) istat = cudaMemPrefetchAsync(psi, size(psi), mydev, 0)
      !
      inquire(file='data/flddpsi.bin', exist=is_data)
      if(myid.eq.0.and.is_data) call execute_command_line('mv data/flddpsi.bin data/flddpsi_old.bin')
      call load(action_load,trim(restart_dir)//'flddpsi.bin',n, dpsidtrko(1:n(1),1:n(2),1:n(3)))
      !@cuf if ( (mod(istep,isave)).eq. 0 .and. (.not. is_done) .and. (.not.kill) ) istat = cudaMemPrefetchAsync(dpsidtrko, size(dpsidtrko), mydev, 0)
#endif
#if defined(_SCALAR_TRANS) 
      inquire(file='data/fldc1.bin', exist=is_data)
      if(myid.eq.0.and.is_data) call execute_command_line('mv data/fldc1.bin data/fldc1_old.bin')
      call load(action_load,trim(restart_dir)//'fldc1.bin',n, c1(1:n(1),1:n(2),1:n(3)))

      inquire(file='data/fldc2.bin', exist=is_data)
      if(myid.eq.0.and.is_data) call execute_command_line('mv data/fldc2.bin data/fldc2_old.bin')
      call load(action_load,trim(restart_dir)//'fldc2.bin',n, c2(1:n(1),1:n(2),1:n(3)))

      inquire(file='data/flddc1.bin', exist=is_data)
      if(myid.eq.0.and.is_data) call execute_command_line('mv data/flddc1.bin data/flddc1_old.bin')
      call load(action_load,trim(restart_dir)//'flddc1.bin',n, dc1dtrko(1:n(1),1:n(2),1:n(3)))

      inquire(file='data/flddc2.bin', exist=is_data)
      if(myid.eq.0.and.is_data) call execute_command_line('mv data/flddc2.bin data/flddc2_old.bin')
      call load(action_load,trim(restart_dir)//'flddc2.bin',n, dc2dtrko(1:n(1),1:n(2),1:n(3)))

#endif
#if defined(_HEAT_TRANSFER)
      inquire(file='data/fldtmp.bin', exist=is_data)
      if(myid.eq.0.and.is_data) call execute_command_line('mv data/fldtmp.bin data/fldtmp_old.bin')
      call load(action_load,trim(restart_dir)//'fldtmp.bin',n, tmp(1:n(1),1:n(2),1:n(3)))

      inquire(file='data/flddent.bin', exist=is_data)
      if(myid.eq.0.and.is_data) call execute_command_line('mv data/flddtmp.bin data/flddtmp_old.bin')
      call load(action_load,trim(restart_dir)//'flddtmp.bin',n, dentdtrko(1:n(1),1:n(2),1:n(3)))
#endif
      !
      inquire(file='data/scalar.out', exist=is_data)
      if(myid.eq.0.and.is_data) call execute_command_line('mv data/scalar.out data/scalar_old.out')
      call load_scalar(action_load,trim(restart_dir)//'scalar.out',time,istep,dto)
      !
      if(myid.eq.0) print*, '*** Checkpoint saved at time = ', time, 'time step = ', istep, '. ***'
      !
      call profiler_stop("OUT:isave")
      !
    endif
#endif
    !
#if defined(_TIMING)
    dt12 = MPI_WTIME()-dt12
    call MPI_ALLREDUCE(dt12,dt12av ,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(dt12,dt12min,1,MPI_REAL_RP,MPI_MIN,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(dt12,dt12max,1,MPI_REAL_RP,MPI_MAX,MPI_COMM_WORLD,ierr)
    !
    var(:) = 0._rp
    var(1) = 1._rp*istep
    var(2) = time
    var(3) = dt12av/(1._rp*product(dims))
    var(4) = dt12min
    var(5) = dt12max
    call out0d(trim(datadir)//'performance.out',5,var)
#endif
    !
  enddo
  !
  ! clear ffts
  !
  call fftend(arrplanp)
  !
  if(myid.eq.0.and.(.not.kill)) print*, '*** Fim ***'
  !
  call profiler_report()
  !
  call decomp_2d_finalize
  call MPI_FINALIZE(ierr)
  call exit
  !
end program flutas
