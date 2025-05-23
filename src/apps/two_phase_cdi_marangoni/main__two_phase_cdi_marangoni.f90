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
                            bounduvw

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
                            static_contact_angle_m,&
                            update_property 
  use mod_di_acdi   , only: rk_c1,initc1 

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
                            restart_sph,restart_mph,restart_mph_scalar,cfl,    &
                            constant_dt,dt_input, &
                            inivel, &
                            itot,jtot,ktot,dims_in, &
                            nthreadsmax, &
                            gr, &
                            is_outflow,no_outflow,is_forced, &
                            rho1,rho2,rho0,mu1,mu2,cbcpsi,bcpsi,late_init,i_late_init, &
                            rho0, &

                            sigma, sigma_t,theta, & 


                            cbcsca, bcsca,diff1,diff2,  &
                            n,ng,l,dl,dli, &
                            bulk_ftype,rkcoeff, &
                            time_scheme,space_scheme_mom,n_stage, &
                            gamma_v,eps_int,inipsi,&
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
  real(rp), dimension(:,:,:), allocatable :: u,v,w,p,mfx,mfy,mfz,uold,vold,wold
  real(rp), dimension(:,:,:), allocatable :: pold
  real(rp), dimension(:,:,:), allocatable :: mu,rho,muold,rhoold,rho_new
  real(rp), dimension(:,:,:), pointer :: dudtrko,dvdtrko,dwdtrko
  real(rp), dimension(:), allocatable :: dzc  ,dzf  ,zc  ,zf  ,dzci  ,dzfi
  real(rp), dimension(:), allocatable :: dzc_g,dzf_g,zc_g,zf_g,dzci_g,dzfi_g
  !
  real(rp), dimension(:,:,:)  , allocatable :: psi,kappa,psiold,lsold,kappaold 
  real(rp), dimension(:,:,:,:), allocatable :: nor,norold                   
  real(rp), dimension(:,:,:)  , allocatable :: ls,dpsidtrko
  real(rp), dimension(:,:,:)  , allocatable :: c1,c2,dc1dtrko,dc2dtrko
  real(rp), dimension(:,:,:)  , allocatable :: ent,entold,rhoh,tmp,cpp,cond,rhocp,dentdtrko
  real(rp), dimension(:,:,:)  , allocatable :: cppold,condold,rhocpold
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
  real(rp) :: rho0i,u_reg
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
  real(rp) :: psi_min,psi_max,vol_p1,vol_c1,vol_c2,vol_ct,ke_tot

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
           pold(0:n(1)+1,0:n(2)+1,0:n(3)+1))
  allocate(dudtrko(n(1),n(2),n(3)), &
           dvdtrko(n(1),n(2),n(3)), &
           dwdtrko(n(1),n(2),n(3)))
  allocate(psi(      0:n(1)+1,0:n(2)+1,0:n(3)+1)  , &
           psiold(       0:n(1)+1,0:n(2)+1,0:n(3)+1)  , &
           lsold(       0:n(1)+1,0:n(2)+1,0:n(3)+1)  , &
           kappa(    0:n(1)+1,0:n(2)+1,0:n(3)+1)  , &
           kappaold(    0:n(1)+1,0:n(2)+1,0:n(3)+1)  , &
           mu(       0:n(1)+1,0:n(2)+1,0:n(3)+1)  , &
           muold(       0:n(1)+1,0:n(2)+1,0:n(3)+1)  , &
           rhoold(      0:n(1)+1,0:n(2)+1,0:n(3)+1)  , &
           rho_new(      0:n(1)+1,0:n(2)+1,0:n(3)+1)  , &
           rho(      0:n(1)+1,0:n(2)+1,0:n(3)+1)  , &
           nor(      0:n(1)+1,0:n(2)+1,0:n(3)+1,3), &
           norold(      0:n(1)+1,0:n(2)+1,0:n(3)+1,3), &
           ls(       0:n(1)+1,0:n(2)+1,0:n(3)+1)  , &
           dpsidtrko(0:n(1)+1,0:n(2)+1,0:n(3)+1))
  allocate(c1(0:n(1)+1,0:n(2)+1,0:n(3)+1)        , &
           c2(0:n(1)+1,0:n(2)+1,0:n(3)+1)        , &
           ent(0:n(1)+1,0:n(2)+1,0:n(3)+1)        , &
           entold(0:n(1)+1,0:n(2)+1,0:n(3)+1)        , &
           rhoh(0:n(1)+1,0:n(2)+1,0:n(3)+1)        , &
           tmp(0:n(1)+1,0:n(2)+1,0:n(3)+1)        , &
           cond(0:n(1)+1,0:n(2)+1,0:n(3)+1)        , &
           condold(0:n(1)+1,0:n(2)+1,0:n(3)+1)        , &
           cpp(0:n(1)+1,0:n(2)+1,0:n(3)+1)        , &
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
  !@cuf istat = cudaMemPrefetchAsync(zc_g , size( zc_g), mydev, 0)
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
          dc1dtrko(i,j,k)  = 0.d0
          dc2dtrko(i,j,k)  = 0.d0
          dentdtrko(i,j,k) = 0.d0
        enddo
      enddo
    enddo
      !
      !
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
  !update the quantities derived from VoF using the lastest available VoF field
  !
  !call update_vof(n,dli,nh_d,dzc,dzf,nh_v,halo_v,psi,nor,cur_t,kappa,d_thinc)
  !
  call update_property(n,(/rho1,rho2/),psi,rho)
  call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,rho)
  call update_property(n,(/mu1 ,mu2 /),psi,mu ) 
  call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,mu)
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
  call boundp(cbcnor,n,bcnor,nh_d,nh_v,halo_v,dl,dzc,dzf,nor(:,:,:,1))
  call boundp(cbcnor,n,bcnor,nh_d,nh_v,halo_v,dl,dzc,dzf,nor(:,:,:,2))
  call boundp(cbcnor,n,bcnor,nh_d,nh_v,halo_v,dl,dzc,dzf,nor(:,:,:,3))
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
  vol_ct  = sum(c2(1:n1,1:n2,1:n3)+c1(1:n1,1:n2,1:n3))*dl(1)*dl(2)*dl(3)
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
  var(9) = vol_ct
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
  ! main loop
  !
  if(myid.eq.0) print*, '*** Calculation loop starts now ***'
  is_done = .false.
  !
  !@cuf istat = cudaMemGetInfo( freeMem, totMem )
  !@cuf if(myid.eq.0) print*, 'Used memory = ', totMem - freeMem
  !
  !call exit
  do while(.not.is_done)
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
    ! 1. DI advection and properties update --> psi^(n+1) 
    !
    call profiler_start("DI", tag = .true., tag_color = COLOR_YELLOW)
    !
    !
    call update_property(n,(/rho1,rho2/),psi,rho)
    call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,rho)
    !
    call update_property(n,(/mu1 ,mu2 /),psi,mu)
    call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,mu)
    !
    !
#if defined(_TWO_PHASE)    
    !
#if defined(_USE_ACDI)
    call psi_to_ls(n,dl,nh_p,psi,eps_int,ls)
    call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,ls)
    call cmpt_norm(n,dli,nh_p,ls,kappa,nor)
#else
    call cmpt_norm(n,dli,nh_p,psi,kappa,nor)
#endif
    call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,kappa       )
    call boundp(cbcnor,n,bcnor,nh_d,nh_v,halo_v,dl,dzc,dzf,nor(:,:,:,1))
    call boundp(cbcnor,n,bcnor,nh_d,nh_v,halo_v,dl,dzc,dzf,nor(:,:,:,2))
    call boundp(cbcnor,n,bcnor,nh_d,nh_v,halo_v,dl,dzc,dzf,nor(:,:,:,3))
    !
#if defined(_SCALAR_TRANS)
    call rk_c1(f_t1,f_t2,n,dli,nh_d,nh_u,nh_v, &
                dzc,dzf,diff1,diff2,u_reg,eps_int,u,v,w,c1,c2,psi,nor,dc1dtrko,dc2dtrko,+1)!SCALAR TRANS. (Scalar inside the dispersed phase) 
    call rk_c1(f_t1,f_t2,n,dli,nh_d,nh_u,nh_v, &
                dzc,dzf,diff1,diff2,u_reg,eps_int,u,v,w,c1,c2,psi,nor,dc1dtrko,dc2dtrko,-1)!SCALAR TRANS. (Scalar outside the dispersed phase) 
    call boundp(cbcsca,n,bcsca,nh_d,nh_v,halo_v,dl,dzc,dzf,c1)
    call boundp(cbcsca,n,bcsca,nh_d,nh_v,halo_v,dl,dzc,dzf,c2)
#endif
   !
   !
   psiold(:,:,:)    = psi(:,:,:)
   lsold(:,:,:)     = ls(:,:,:)
   norold(:,:,:,:)  = nor(:,:,:,:)
   kappaold(:,:,:)  = kappa(:,:,:)
   !
   !Advecting the order parameter
   !
   call cmpt_umax(n,nh_u,u,v,w,gamma_v,u_reg)
   call rk_psi(f_t1,f_t2,n,dli,nh_d,nh_u,nh_v, &
                dzc,dzf,u_reg,eps_int,u,v,w,psi,ls,nor,kappa,dpsidtrko)
   ! 
#if defined(_CONTACT_LINE)
   call static_contact_angle_m(n,dl,nh_v,eps_int,theta,psi,psi_bound)
   call boundsb_array(n,psi_bound,nh_v,halo_v,psi)
#else  
   call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,psi)
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
   !
   !
   call update_property(n,(/rho1,rho2/),psi,rho)
   call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,rho)
   !
   call update_property(n,(/mu1 ,mu2 /),psi,mu)
   call boundp(cbcpsi,n,bcpsi,nh_d,nh_v,halo_v,dl,dzc,dzf,mu)
   ! 
   !
   call profiler_stop("DI")
   !
   !
#if defined(_NAVIER_STOKES)
    !
    ! 2. two-fluid Navier-Stokes --> (u,v,w,p)^(n+1)
    !
    rho0i = 1._rp/rho0
    !
    !
    call profiler_start("RK", tag = .true., tag_color = COLOR_RED)
    !
    !Compute mass fluxes at the faces which will be used in the momentum equation 
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
#if defined(_USE_ACDI)
    do k=1,n(3)
     do j=1,n(2)
      do i=1,n(1)

        ip = i + 1
        im = i - 1
        jp = j + 1
        jm = j - 1
        kp = k + 1
        km = k - 1
      !
      ! along x-dir
      !  
#if defined(_TWOD)        
      regx   = 0.d0
#else
      gradpx =   ( psiold(ip,j,k) - psiold(i,j,k) ) * dli(1)

      psifx  =   ( 0.25d0*(1.d0 - tanh(lsold(i,j,k)/(2.d0*eps_int))**2)*norold(i,j,k,1) + &
                   0.25d0*(1.d0 - tanh(lsold(ip,j,k)/(2.d0*eps_int))**2)*norold(ip,j,k,1))*0.5d0 

      regx   =   u_reg*( eps_int*gradpx - psifx)
#endif
      !
      ! along y-dir
      !
      gradpy =   ( psiold(i,jp,k) - psiold(i,j,k) ) * dli(2)

      psify  =   ( 0.25d0*(1.d0 - tanh(lsold(i,j,k)/(2.d0*eps_int))**2)*norold(i,j,k,2) + &
                   0.25d0*(1.d0 - tanh(lsold(i,jp,k)/(2.d0*eps_int))**2)*norold(i,jp,k,2))*0.5d0

      regy   =   u_reg*( eps_int*gradpy - psify)
      !
      ! along z-dir
      !
      gradpz =   ( psiold(i,j,kp) - psiold(i,j,k) ) * dli(3)

      psifz  =   ( 0.25d0*(1.d0 - tanh(lsold(i,j,k)/(2.d0*eps_int))**2)*norold(i,j,k,3) + &
                   0.25d0*(1.d0 - tanh(lsold(i,j,kp)/(2.d0*eps_int))**2)*norold(i,j,kp,3))*0.5d0

      regz   =   u_reg*( eps_int*gradpz - psifz)
      !
      !
      mfx(i,j,k) = 0.5d0*(rhoold(ip,j,k)+rhoold(i,j,k))*u(i,j,k) - (rho1-rho2)*regx
      mfy(i,j,k) = 0.5d0*(rhoold(i,jp,k)+rhoold(i,j,k))*v(i,j,k) - (rho1-rho2)*regy
      mfz(i,j,k) = 0.5d0*(rhoold(i,j,kp)+rhoold(i,j,k))*w(i,j,k) - (rho1-rho2)*regz
      !
      !
      enddo
     enddo
    enddo
#else
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
      mfx(i,j,k) = 0.5d0*(rhoold(ip,j,k)+rhoold(i,j,k))*u(i,j,k) - (rho1-rho2)*regx
      mfy(i,j,k) = 0.5d0*(rhoold(i,jp,k)+rhoold(i,j,k))*v(i,j,k) - (rho1-rho2)*regy
      mfz(i,j,k) = 0.5d0*(rhoold(i,j,kp)+rhoold(i,j,k))*w(i,j,k) - (rho1-rho2)*regz
      !
      enddo
     enddo
    enddo
    !
#endif
    call bounduvw(cbcvel,n,bcvel,nh_d,nh_u,halo_u,no_outflow,dl,dzc,dzf,mfx,mfy,mfz)
    !
    !
    !input of the prediction is velocities but the output is product of rho.u
    !
    !old density and viscoity are being used for diffusion and advective terms
    call rk(space_scheme_mom,f_t1,f_t2,n(1),n(2),n(3),dxi,dyi,dzi,nh_d,nh_u,dzci,dzfi, &
            u,v,w,mfx,mfy,mfz,muold,rhoold,dudtrko,dvdtrko,dwdtrko)
    !
    !
    ! 2b. add the source terms
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
    ! 2c. construct the Poisson equation
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
    call fillps(n(1),n(2),n(3),nh_d,nh_u,dxi,dyi,dzi,dzfi,f_t12_i,rho0,u,v,w,p)
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
    ! 2d. correct the velocity and update the pressure
    !
    !u,v,w are the product of velocity and denisty as an input and output 
    !
    call correc(n(1),n(2),n(3),nh_d,nh_u,dxi,dyi,dzi,dzci,f_t12,rho0,p,u,v,w,rho)
    call bounduvw(cbcvel,n,bcvel,nh_d,nh_u,halo_u,is_outflow,dl,dzc,dzf,u,v,w)
    !
    !
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
      if(divmax.gt.small.or.divtot.ne.divtot) then
        if(myid.eq.0) print*, 'ERROR: maximum divergence is too large.'
        if(myid.eq.0) print*, 'Aborting...'
        is_done = .true.
        kill = .true.
      endif
      !
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
      vol_c1  = sum(c1(1:n(1),1:n(2),1:n(3)))*dl(1)*dl(2)*dl(3)
      vol_c2  = sum(c2(1:n(1),1:n(2),1:n(3)))*dl(1)*dl(2)*dl(3)
      vol_ct  = sum(c1(1:n(1),1:n(2),1:n(3))+c2(1:n(1),1:n(2),1:n(3)))*dl(1)*dl(2)*dl(3)
      call mpi_allreduce(MPI_IN_PLACE,vol_c1 ,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,vol_c2 ,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,vol_ct ,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      !
#endif
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
      var(9) = vol_ct
      var(10) = ke_tot
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
