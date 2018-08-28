!********************************* modules *********************************
! Keep all modules used in the program here
! To restore *.mod files use make filmods
!***************************************************************************

module dimM
       implicit none
       integer,parameter :: dimi = 10 ! should be equal to number of state index parameters
       integer :: dimj
       integer :: dimjp = 0 ! number of mass elements
       integer,allocatable,dimension(:,:) :: thr_domain
       integer :: lazy_thread
end module dimM

module stateIndex
       implicit none
       ! old labels x = 1, z = 2, vx = 3, vz = 4, d = 5, p = 6, t = 7, b = 8, mf = 9, pdg = 10
       ! labels for state vectors [see initialConditions subroutine]
       integer, parameter :: x = 1, z = 2, vx = 3, vz = 4, d = 5, p = 6, t = 7, b = 8, mf = 9, pdg = 10
       !x=1,y=2,z=3,  vx=4,vy=5,vz=6,  d=7,p=8,t=9,   b=10,mf=11,pdg=12, &
       !& x_int=13,y_int=14,z_int=15
end module

module stateMod
       use dimM
       implicit none
       real,allocatable,dimension(:,:) :: state,newState,oldState,stateTemp,statePlas,newStatePlas,oldStatePlas,statePlasTemp
       real :: time,tau ! time steps
       real :: tau_min_limit
       real :: tau_scale
end module stateMod

module rkVectors
       implicit none
       real,allocatable,dimension(:,:) :: s1,s2,s3,k1,k2,k3,k4,sa
end module rkVectors

module initCond
! used to produce different initial conditions
       implicit none
       real :: apex != 40.0 ! most tailward point (in Re)
       real :: scaleStep = 1.0 ! scale factor for initial mass element spacing
       real :: boundary ! ionospheric boundary (in Re)
       real :: step = 0.1408 ! element spacing (in Re) for dim = 499,apex = 40,boundary = 5,CW99, initial guesss for tracing
       real :: pressureScale ! ratio of filament pressure to background pressure
       real :: pressureScale_in
       real,parameter :: area0 = 1.0 ! cross sectional area of equitorial foot point [not used for now]
       real :: Phi ! flux constant [not used for now]
       logical :: plotPVG,useAnalyticModelGrid,use_dPdK,dontRunSim,useAM03ModelGrid,buildTsyData,useTsyModelGrid,use_bao_density
       logical :: rcmdata
       logical :: forceModel,dPdK_fixed,useBackgroundDensity,useKpk
       real :: Kpk,K_peak,fmsection
       real :: ampli
       logical :: buildKubyshData
       logical :: use_gallagherDensity,use_tm2003_density
       logical :: trace_scale_step
       real :: plotPVG_xmin,plotPVG_xmax
       integer :: plotPVG_grid
       logical :: save_external_file
       logical :: load_external_file
       logical :: save_after_sim
       real :: fac
end module initCond

module parametersBackground
       implicit none
       real :: A0
       real :: ht
       real :: alpha
       real :: alpha0
       real,parameter :: kp = 1.539 ! pressure scale parameter
       real,parameter :: ds = 10.0/4.0 ! dayside boundary (units of ht)
       integer :: nSum = 1000 !number of terms to use in potential sum (for now it's dynamically chosen)
       real :: tBack ! n * 1 keV/kb  background temp (K)
       real :: dBack
       logical :: useGrid
       integer :: interpType
       logical :: wallBound
       real :: txmax,txmin,tzmax,tzmin,tyslice,swpress,dst,byimf,bzimf,tsyg1,tsyg2,tilt,sw_density,sw_velocity
       integer :: tgridx,tgridz
       integer :: kpindex
       character(len=3) :: tsyVersion
       logical :: tecPlot3D
       real :: pressureCap
end module parametersBackground

module constants
       implicit none
       real,parameter :: nano  = 1.0e9                ! factor for unit conversion
       real,parameter :: atto  = 1.0e18                        ! factor for unit conversion
       real,parameter :: cm       = 100.0                               ! centimeter conversion
       real,parameter :: cmC       = 1.0e6                               ! centimeter cube conversion
       real,parameter :: Pi    = acos(-1.0)
       real,parameter :: Re    = 6.37104e+08          ! Earth radius (cm)
       real,parameter :: mu0   = nano*4.0*Pi*1.0e-7   ! vacuum permeability (nH/m)
       real,parameter :: tkev  = 1.6045e7                 ! 1 keV/kb  background temp (K)
       real,parameter :: q     = 1.602192e-19         ! unit charge (C)
       real,parameter :: mi    = 1.672622e-27         ! ion mass (kg)
       real,parameter :: kb    = cmC*nano*1.38065e-23 ! boltzman constant (1e6cm^3/m^3*nJ/K)
       real,parameter :: h     = 0.04                 ! step size for background gradients (in Re)
       real,parameter :: gamma = 5.0/3.0                        ! ratio of specific heats
       real,parameter :: condScale = 1.0/nano
       real,parameter :: jolt = 1.0e-1
       real,parameter :: jt   = 20.0
       real,parameter :: js   = jt
end module constants

module variables1 
! runtime variable for numerical methods and number of steps [see subroutine initialize]
       implicit none
       integer :: n,m,method,algIn,drvIn,cnfIn,run,runTemp
       real :: tMax,tInterv,nInterv
       logical :: adaptL,plotL,timePlotL,stepPlotL
       logical :: do_parallel,do_sweep
       integer :: nOutInterv
       integer :: threads_to_use
       character(len=100) :: sw_min,sw_max,sw_var
       integer :: sw_points
       logical :: rebuild_background
       logical :: write_binary_data
       logical :: write_full_filament
       logical :: write_midpoint
       logical :: write_boundary
       logical :: bbfrun
end module variables1

module dnams
       character(len=100) :: fileName,fileDir,datapath,tecPlotDir,tecPlotFileName,am03Model,am03ModelDir, &
       & tsyModelDir,tsyModel,fricDir,external_file
    real :: am03ModelTime
    integer :: upscalingFactor
end module dnams

module filConstants
       implicit none
       real :: totalM
       real :: Kfil
       real :: Ksum
       real :: wall
       real :: xEqual
       real :: ohm33
       real :: phase = 40 * 2.0*acos(-1.0)/360.0
       real :: Amp
       real,allocatable,dimension(:,:) :: s_0
end module filConstants

module errorMod
       implicit none
       logical :: errorFlag = .false.
       logical :: algLogic1 = .false.
       logical :: errorFlag2 = .false.
       logical :: firstbgSetup = .true.
end module errorMod

module boundMod
       implicit none
       integer :: bVelType = 1
       real :: Jthreshold = 0.0
       logical :: Jsticky = .false.
       logical :: Jswitch = .false.
       logical :: Jzjump = .false.
       logical :: radialDamping = .false.
       real :: sigmaP   ! Pederson conductivity, CW1999 has 4 mhos (Gmhos)
       real :: drive_tmin,drive_tmax,drag_tmin,drag_tmax
       logical:: use_drag_force,use_driver_force
       real :: drive_coeff,drag_coeff
       real :: drive_omega_0,drive_omega_f
       integer :: iopen
end module boundMod

module cont
! arrays to store the contour of the background grid, coni is x coordinates and conj is z coordinates
       implicit none
       real,allocatable,dimension(:) :: coni,conj
end module cont

module potential
! background potential and field grids
       implicit none
       integer :: gridi,gridj,gridk
       real,allocatable,dimension(:,:) :: Ag,Bxg,Byg,Bzg,Pg,Ng,Ptotg,testMesh
       real,allocatable,dimension(:,:,:) :: oAg,oBxg,oByg,oBzg,oPg
       real,allocatable,dimension(:) :: xg,yg,zg,oxg,oyg,ozg
       integer,allocatable,dimension(:,:) :: recallPos
       real,allocatable,dimension(:,:) :: dx_Bxg,dz_Bxg,dx_Bzg,dz_Bzg,dx_Ptotg,dz_Ptotg,dx_Pg,dz_Pg
       real,allocatable,dimension(:,:,:) :: Ag_d,Bxg_d,Byg_d,Bzg_d,Pg_d,Ng_d,Ptotg_d,testMesh_d
       real,allocatable,dimension(:) :: PVG_bg,dsob3_bg
       integer :: PVGgrid
       !*****************************!
       logical :: writeGrid = .true. ! calculates a new grid if true, otherwise use old grid from file
       logical :: saveGrid = .true.
       !*****************************!
       real :: midx = 0.01 ! stepoff from x = 0.0
       real :: xminl,xmaxl,yminl,ymaxl,zminl,zmaxl
       logical :: useAnalyticDipole
       integer :: plasGridi
       real,allocatable,dimension(:) :: plasGrid,plasGridx,plasGridz
       real :: dipole_moment = 30574.0
end module potential

module fricoptions
       implicit none
       integer :: initialization_option
       logical :: NSsymmetric_option
       ! (logical) flag to use explicit inner boundary conditions
       logical :: fric_use_inner_bc
       ! Viscosity for momentum equation (test only, not used, use 0)
       real :: fric_vis_mom           
       ! Limiter_choice (integer) for momentum (0--none, 1--limit3, 2--smoothf)
       integer :: fric_limiter             
       ! time interval to adjust pv^gamma (-ve is dont)
       integer :: fric_pvg_corr_inter       
       ! push_rho: -1 is no, 1 is yes
       integer :: fric_push_rho_opt            
       ! give sheath a high density? 1 = yes, 0 = no
       integer :: fric_sheathrho
       ! stopsign:       Set v=0 when K.E. peaks?
       logical :: fric_stopsign       
       ! tgamma:       Time to switch gamma to 5/3
       real :: fric_tgamma
       ! gamma:       Value of gamma to start with
       real :: fric_gammainp           
       ! rhomin:       (density floor)
       real :: fric_rhomin       
       ! pressmin0       (minimum pressure to allow)
       real :: fric_pressmin       
       ! cfl:              Courant number - default: 0.25
       real :: fric_cfl       
       ! force determination, 1 = volume integral, 0 = grid average
       integer :: fric_force_norm        
       ! max fric parameter - 10
       real :: fric_max_factor        
       ! min fric parameter 0.0001 
       real :: fric_min_factor               
       ! friction alog (1=orig,2=mod orig,3=velmax,4=mod velmax) 
       integer :: fric_algorithm        
       ! pressure increase factor
       real :: fric_pressure_factor               
       ! print diagnostic (screen/diag.files) every N iterations
       integer :: fric_print_diagnostic       
       ! minimum force imbalance
       real :: fric_min_avg_force        
       ! max number of iterations to run
       integer :: fric_max_steps            
       ! inwrite: datafile write interval
       integer :: fric_interval               
       ! starting friction value
       real :: fric_alpha       
end module
