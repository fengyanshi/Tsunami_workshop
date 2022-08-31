! --- For AMR modeling !ykchoi
MODULE LOCAL
       USE PARAM
       IMPLICIT NONE
       SAVE

! grid


       
         REAL(SP),DIMENSION(:,:,:),ALLOCATABLE :: DX_TILE,DY_TILE


! bathy
         REAL(SP),DIMENSION(:,:,:), ALLOCATABLE :: Depth_Tile, &
                 DepthX_Tile, DepthY_Tile
         REAL(SP),DIMENSION(:,:,:),ALLOCATABLE :: CD_Tile             !ykchoi 08/10/2018                 
         INTEGER,DIMENSION(:,:,:),ALLOCATABLE :: MASK_Tile,MASK9_Tile
         INTEGER,DIMENSION(:,:,:),ALLOCATABLE :: MASK_STRUC_Tile       !ykchoi 07/24/2018




         REAL(SP),DIMENSION(:,:,:),ALLOCATABLE :: Lat_Theta_Tile, Coriolis_Tile
         REAL(SP),DIMENSION(:,:,:),ALLOCATABLE :: SlopeX_Tile, SlopeY_Tile   !ykchoi 09/11/2018


! variable
         REAL(SP),DIMENSION(:,:,:),ALLOCATABLE :: U_Tile,V_Tile,Eta_Tile  
         REAL(SP),DIMENSION(:,:,:),ALLOCATABLE :: Ubar_Tile,Vbar_Tile       

       CONTAINS

SUBROUTINE ALLOCATE_VAR_TILE
       USE PARAM

       USE GLOBAL, ONLY: K,NumGrid,GridDimX,GridDimY,Nghost,ier,px,py,  &
                         MaxDimX,MaxDimX1,MaxDimY,MaxDimY1      




       IMPLICIT NONE


       MaxDimX = MAXVAL(GridDimX)/px+2*Nghost
       MaxDimY = MAXVAL(GridDimY)/py+2*Nghost




       MaxDimX1 = MaxDimX+1
       MaxDimY1 = MaxDimY+1
       
! static




       ALLOCATE ( DX_TILE(MaxDimX,MaxDimY,NumGrid) )
       ALLOCATE ( DY_TILE(MaxDimX,MaxDimY,NumGrid) )

       ALLOCATE ( Depth_Tile(MaxDimX,MaxDimY,NumGrid) )
       ALLOCATE ( DepthX_Tile(MaxDimX1,MaxDimY,NumGrid) )
       ALLOCATE ( DepthY_Tile(MaxDimX,MaxDimY1,NumGrid) )
       ALLOCATE ( MASK_Tile(MaxDimX,MaxDimY,NumGrid) )
       ALLOCATE ( MASK9_Tile(MaxDimX,MaxDimY,NumGrid) )
       ALLOCATE ( MASK_STRUC_Tile(MaxDimX,MaxDimY,NumGrid) )     !ykchoi 07/24/2018
       ALLOCATE ( CD_Tile(MaxDimX,MaxDimY,NumGrid) )             !ykchoi 08/10/2018



       ALLOCATE ( Lat_theta_Tile(MaxDimX,MaxDimY,NumGrid), &
	            Coriolis_Tile(MaxDimX,MaxDimY,NumGrid), &
				SlopeX_Tile(MaxDimX,MaxDimY,NumGrid),  &    !ykchoi 09/11/2018
				SlopeY_Tile(MaxDimX,MaxDimY,NumGrid) )      !ykchoi 09/11/2018

! dynamic
       ALLOCATE ( U_Tile(MaxDimX,MaxDimY,NumGrid))
       ALLOCATE ( V_Tile(MaxDimX,MaxDimY,NumGrid))
       ALLOCATE ( Ubar_Tile(MaxDimX,MaxDimY,NumGrid))
       ALLOCATE ( Vbar_Tile(MaxDimX,MaxDimY,NumGrid))

       ALLOCATE ( Eta_Tile(MaxDimX,MaxDimY,NumGrid)) 

! initialization
        MASK_TILE = 1
        MASK9_TILE = 1
        MASK_STRUC_TILE = 1   !ykchoi 07/24/2018
        U_TILE = ZERO
        V_TILE = ZERO
        Ubar_TILE = ZERO
        Vbar_TILE = ZERO
        Eta_TILE = ZERO
        CD_TILE = ZERO        !ykchoi 08/10/2018

END SUBROUTINE ALLOCATE_VAR_TILE

END MODULE LOCAL 
