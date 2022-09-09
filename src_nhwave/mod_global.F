     Module GLOBAL
# if defined (PARALLEL)
     use mpi
# endif
     implicit none

     ! define precision
# if defined (DOUBLE_PRECISION)
     integer, parameter :: SP=8
# if defined (PARALLEL)
     integer, parameter :: MPI_SP=MPI_DOUBLE_PRECISION
# endif
# else
     integer, parameter :: SP=SELECTED_REAL_KIND(6,30)
# if defined (PARALLEL)
     integer, parameter :: MPI_SP=MPI_REAL
# endif
# endif

     ! define parameters
     real(SP), parameter :: pi=3.141592653
     real(SP), parameter :: Small=0.000001
     real(SP), parameter :: Large=10000000.0
     real(SP), parameter :: Grav=9.81
     real(SP), parameter :: Zero=0.0
     real(SP), parameter :: One=1.0
     real(SP), parameter :: Rho0=1000.0
     real(SP), parameter :: RhoA=1.20

! fyshi change to integer 12/15/2011
     INTEGER, parameter :: MaxNumFreq=50
     INTEGER, parameter :: MaxNumDir=50

     ! ghost cells (>=1)
     integer, parameter :: Nghost=2

     ! define characters
     character(len=80) :: TITLE
     character(len=80) :: RESULT_FOLDER
     character(len=80) :: HIGH_ORDER
     character(len=80) :: TIME_ORDER
     character(len=80) :: WaveMaker
     character(len=80) :: DEPTH_TYPE
     character(len=80) :: dt_constraint

# if defined (PARALLEL)
     integer :: myid,ier
     integer :: comm2d
     integer :: n_west,n_east,n_suth,n_nrth
     integer :: npx,npy
     integer :: ndims=2
     integer :: NumP
     integer, dimension(2) :: dims,coords
     logical, dimension(2) :: periods
     logical :: reorder=.true.
# endif


# if defined (BUBBLE)
     integer, parameter :: Mg = 20
     real(SP), parameter :: Surface_Tension = 71.97e-3
     real(SP), parameter :: Con_b = 0.18
     real(SP), parameter :: TIM_B = 16.0
     real(SP), dimension(:), allocatable :: Rbg,Wbg,Entrain
     real(SP), dimension(:,:,:,:), allocatable :: Nbg,DNbg,DNbg0
     real(SP), dimension(:,:,:), allocatable :: Vbg
# endif

     ! define output logical parameters
     logical :: NON_HYDRO,VISCOUS_FLOW,SPONGE_ON,OUT_H,OUT_E,OUT_U,OUT_V,OUT_W,OUT_P, &
                OUT_K,OUT_D,OUT_S,OUT_C,OUT_B,PERIODIC,OUT_2DUV
! fyshi add tmp4preview

     LOGICAL :: OUT_preview

! DJN
     LOGICAL :: OUT_SLD     

     ! variables
     integer :: It_Order,Ibeg,Iend,Iend1,Jbeg,Jend,Jend1,Kbeg,Kend,Kend1,PX,PY,IVturb,IHturb,  &
                Mglob,Nglob,Kglob,Mloc,Nloc,Kloc,Mloc1,Nloc1,Kloc1,Icount,RUN_STEP,Ivgrd,SIM_STEPS,Ibot, &
                NumFreq,NumDir,NSTAT
     integer :: Bc_X0,Bc_Xn,Bc_Y0,Bc_Yn,Bc_Z0,Bc_Zn
     real(SP) :: dt,dt_old,dt_min,dt_max,dt_ini,dx,dy,Theta,CFL,VISCOUS_NUMBER,MinDep,TIME,TOTAL_TIME,Plot_Intv,  &
                 Screen_Intv,Screen_Count,Plot_Count,Visc,Cvs,Chs,Zob,Tke_min,Eps_min,Cd0,Plot_Start,Plot_Intv_Stat, &
                 Plot_Count_Stat
     real(SP) :: Amp_Wave,Per_Wave,Dep_Wave,Theta_Wave,Freq(MaxNumFreq),Dire(MaxNumDir),Wave_Spc2d(MaxNumDir,MaxNumFreq), &
                 Random_Phs(MaxNumDir,MaxNumFreq)
     real(SP) :: Sponge_West_Width,Sponge_East_Width,Sponge_South_Width,Sponge_North_Width,R_Sponge,A_Sponge, &
                 Xsource_West,Xsource_East,Ysource_Suth,Ysource_Nrth
     real(SP), dimension(3) :: ALPHA,BETA
     REAL(SP),DIMENSION(:),ALLOCATABLE :: XSTAT,YSTAT

! fyshi added time series boundary condition 12/17/2011
       CHARACTER(LEN=80) :: BoundaryFile
       INTEGER :: NumTimeData
       INTEGER :: icount_tide = 1
       REAL(SP),DIMENSION(:),ALLOCATABLE :: DataU_L,DataEta_L
       REAL(SP),DIMENSION(:),ALLOCATABLE :: DataU_R,DataEta_R
       REAL(SP),DIMENSION(:),ALLOCATABLE :: TimeData
       LOGICAL :: Nudging = .False.
! fyshi added bathymetry file 04/13/2012
       CHARACTER(LEN=80) :: Depth_File
! DJN added landslide thickness file 07/08/2013
       CHARACTER(LEN=80) :: Slide_File

     ! real arrays
     real(SP), dimension(:), allocatable :: x,xc,y,yc,sig,dsig,sigc,Ein_X0,Din_X0,Ein_Xn,Din_Xn
     real(SP), dimension(:,:), allocatable :: Ho,H,Hc,Hfx,Hfy,DeltH,DeltHo,Delt2H,DelxH,DelyH,D,D0,Eta,Eta0, &
                                              SourceX,SourceY,SourceC,DxL,DxR,DyL,DyR,EtaxL,EtaxR,EtayL,EtayR, &
                                              DelxEta,DelyEta,DelxD,DelyD,Uin_X0,Vin_X0,Win_X0,Uin_Xn,Vin_Xn, &
                                              Win_Xn,Bc_Prs,Sponge

!DJN	
     real(SP), dimension(:,:), allocatable :: SldM0, SldN0, SldZ0, SldDZ0
     real(SP), dimension(:,:), allocatable :: SldM1, SldN1, SldZ1, SldDZ1
     real(SP), dimension(:,:), allocatable :: SldMs, SldNs, SldZs, SldDZs
     real(SP), dimension(:,:), allocatable :: SldWZ
     real(SP), dimension(:,:), allocatable :: SldUs, SldVs

     real(SP) :: SLIDE_DENSITY
     real(SP) :: SLIDE_VISCOSITY
     real(SP) :: SLIDE_ROUGHNESS
     real(SP) :: SLIDE_MINTHICK 
     real(SP) :: SLIDE_DIFF 
     real(SP) :: FroudeCap

     real(SP) :: SLIDE_RKC1, SLIDE_RKC2, SLIDE_RKDT    
!DJN

     real(SP), dimension(:,:,:), allocatable :: U,V,W,Omega,P,DU,DV,DW,DU0,DV0,DW0, &
                                                UxL,UxR,VxL,VxR,WxL,WxR,DUxL,DUxR,DVxL,DVxR,DWxL, &
                                                DWxR,UyL,UyR,VyL,VyR,WyL,WyR,DUyL,DUyR,DVyL,DVyR,DWyL,DWyR, &
                                                UzL,UzR,VzL,VzR,WzL,WzR,OzL,OzR,SxL,SxR,SyL,SyR,ExL,ExR,FxL, &
                                                FxR,GxL,GxR,HxL,HxR,EyL,EyR,FyL,FyR,GyL,GyR,HyL,HyR,Ex,Ey,Fx, &
                                                Fy,Fz,Gx,Gy,Gz,Hx,Hy,Hz,DelxU,DelyU,DelzU,DelxV,DelyV,DelzV, &
                                                DelxW,DelyW,DelzW,DelxDU,DelyDU,DelxDV,DelyDV,DelxDW,DelyDW, &
                                                DelzO,Uf,Vf,Wf,Cmu,Cmuht,Cmuvt,Diffxx,Diffxy,Diffxz,Diffyx,  &
                                                Diffyy,Diffyz,Diffzx,Diffzy,Diffzz,DelxSc,DelySc,Rho,Tke,Eps, &
                                                DTke,DEps,DTke0,DEps0,Prod_s,Prod_b

     ! integer arrays
     integer, dimension(:,:), allocatable :: Mask,Mask_Struct,Mask9
     
     ! poisson solvers
     integer  :: itmax,isolver,neqns
     real(SP) :: tol
     real(SP), dimension(:),   allocatable :: Rhs
     integer,  dimension(:),   allocatable :: JCoef
     real(SP), dimension(:,:), allocatable :: Coef

! fyshi add tmp4preview
     REAL(SP),DIMENSION(:,:),ALLOCATABLE :: tmp4preview2D
     REAL(SP),DIMENSION(:,:,:),ALLOCATABLE :: tmp4preview3D

! add landslide parameter fyshi
!     real(SP) :: T_slide,b_slide,w_slide,e_slide,alpha_slide,kb_slide,kw_slide,ac_slide, &
!                 x0_slide,y0_slide,term_v,slope_slide,acceleration_lab

! add base bathymetry for landslide fyshi
     real(SP), dimension(:,:), allocatable :: Dep0,DepC0    

     REAL(SP),DIMENSION(:,:),ALLOCATABLE :: EleMax,EleMin,DepMaxFr,TimeMaxFr
     REAL(SP),DIMENSION(:,:,:),ALLOCATABLE :: FrMax,Froude,FrMax_total,Froude_total

     End Module GLOBAL
