!------------------------------------------------------------------------------------
!
!      FILE main.F
!
!      This file is part of the FUNWAVE-TVD program under the Simplified BSD license
!
!-------------------------------------------------------------------------------------
! 
!    Copyright (c) 2016, FUNWAVE Development Team
!
!    (See http://www.udel.edu/kirby/programs/funwave/funwave.html
!     for Development Team membership)
!
!    All rights reserved.
!
!    FUNWAVE_TVD is free software: you can redistribute it and/or modify
!    it under the terms of the Simplified BSD License as released by
!    the Berkeley Software Distribution (BSD).
!
!    Redistribution and use in source and binary forms, with or without
!    modification, are permitted provided that the following conditions are met:
!
!    1. Redistributions of source code must retain the above copyright notice, this
!       list of conditions and the following disclaimer.
!    2. Redistributions in binary form must reproduce the above copyright notice,
!    this list of conditions and the following disclaimer in the documentation
!    and/or other materials provided with the distribution.
!
!    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
!    ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
!    WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
!    DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
!    ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
!    (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
!    LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
!    ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
!    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
!    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!  
!    The views and conclusions contained in the software and documentation are those
!    of the authors and should not be interpreted as representing official policies,
!    either expressed or implied, of the FreeBSD Project.
!  
!-------------------------------------------------------------------------------------
!
!    FLUXES is subroutine to calculate fluxes at four sides
!
!    HISTORY: 
!        05/06/2010 Fengyan Shi
!
!-------------------------------------------------------------------------------------
SUBROUTINE FLUXES(ng)
     USE GLOBAL
     IMPLICIT NONE
     INTEGER,INTENT(IN) :: ng

      IF(HIGH_ORDER(1:3)=='FOU') THEN       !ykchoi (08/28/2016)

      ! for the fourth-order, Choi used 4th + combined minmod and van leer
      ! it is more stable than 4th + minmod

!  Choi s tests show this is unstable combination, switched off 08/31/2016 
!     IF(HIGH_ORDER(1:9)=='FOURTHMIN') THEN   !ykchoi (08/28/2016)
!       CALL CONSTRUCTION_HO_minmod
!       CALL WAVE_SPEED(Mloc,Nloc,Mloc1,Nloc1,UxL,UxR,VyL,VyR,HxL,HxR,HyL,HyR, &
!            SxL,SxR,SyL,SyR)

       CALL CONSTRUCTION_HO   ! note choi used the name construction_ho_vanleer
                              ! I changed back for consistency

       CALL WAVE_SPEED(Mloc,Nloc,Mloc1,Nloc1,  &
	      UxL(1:Mloc1,1:Nloc),UxR(1:Mloc1,1:Nloc),  &
		  VyL(1:Mloc,1:Nloc1),VyR(1:Mloc,1:Nloc1),  &
		  HxL(1:Mloc1,1:Nloc),HxR(1:Mloc1,1:Nloc),  &
		  HyL(1:Mloc,1:Nloc1),HyR(1:Mloc,1:Nloc1),  &
            SxL(1:Mloc1,1:Nloc),SxR(1:Mloc1,1:Nloc),  &
		  SyL(1:Mloc,1:Nloc1),SyR(1:Mloc,1:Nloc1))

      ELSE                                     

       CALL DelxyFun
       CALL CONSTRUCTION 
       CALL WAVE_SPEED(Mloc,Nloc,Mloc1,Nloc1,  &
	      UxL(1:Mloc1,1:Nloc),UxR(1:Mloc1,1:Nloc),  &
	 		  VyL(1:Mloc,1:Nloc1),VyR(1:Mloc,1:Nloc1),  &
	 		  HxL(1:Mloc1,1:Nloc),HxR(1:Mloc1,1:Nloc),  &
	 		  HyL(1:Mloc,1:Nloc1),HyR(1:Mloc,1:Nloc1),  &
	             SxL(1:Mloc1,1:Nloc),SxR(1:Mloc1,1:Nloc),  &
	 		  SyL(1:Mloc,1:Nloc1),SyR(1:Mloc,1:Nloc1))

     ENDIF

     IF(CONSTR(1:3)=='HLL')THEN
       CALL FLUX_AT_INTERFACE_HLL
     ELSE
       CALL FLUX_AT_INTERFACE
     ENDIF

     CALL BOUNDARY_CONDITION(ng) !For AMR routine 

END SUBROUTINE FLUXES

!-------------------------------------------------------------------------------------
!
!    FLUX_AT_INTERFACE is subroutine to calculate fluxes at four sides
!      using averaging approach (for predictor), dont use it for TVD
!
!    HISTORY: 
!      05/06/2010 Fengyan Shi
!
!-------------------------------------------------------------------------------------
SUBROUTINE FLUX_AT_INTERFACE
     USE GLOBAL

! for averaging approach for predictor
     P=0.5_SP*(PR+PL)
     Fx=0.5_SP*(FxR+FxL)
     Gx=0.5_SP*(GxR+GxL)
     Q=0.5_SP*(QR+QL)
     Fy=0.5_SP*(FyR+FyL)
     Gy=0.5_SP*(GyR+GyL)

END SUBROUTINE FLUX_AT_INTERFACE


!-------------------------------------------------------------------------------------
!
!    FLUX_AT_INTERFACE_HLL is subroutine to do HLL scheme
!
!    HISTORY: 
!    05/06/2010 Fengyan Shi
!
!-------------------------------------------------------------------------------------
SUBROUTINE FLUX_AT_INTERFACE_HLL
     USE GLOBAL
     REAL(SP)::SR,SL,FL,FR,UL,UR

     CALL HLL(Mloc1,Nloc,SxL(1:Mloc1,1:Nloc),SxR(1:Mloc1,1:Nloc),  &
              PL(1:Mloc1,1:Nloc),PR(1:Mloc1,1:Nloc),EtaRxL(1:Mloc1,1:Nloc),  &
			EtaRxR(1:Mloc1,1:Nloc),P(1:Mloc1,1:Nloc))
     CALL HLL(Mloc,Nloc1,SyL(1:Mloc,1:Nloc1),SyR(1:Mloc,1:Nloc1),  &
              QL(1:Mloc,1:Nloc1),QR(1:Mloc,1:Nloc1),EtaRyL(1:Mloc,1:Nloc1),  &
			EtaRyR(1:Mloc,1:Nloc1),Q(1:Mloc,1:Nloc1))
     CALL HLL(Mloc1,Nloc,SxL(1:Mloc1,1:Nloc),SxR(1:Mloc1,1:Nloc),  &
              FxL(1:Mloc1,1:Nloc),FxR(1:Mloc1,1:Nloc),HUxL(1:Mloc1,1:Nloc),  &
			HUxR(1:Mloc1,1:Nloc),Fx(1:Mloc1,1:Nloc))
     CALL HLL(Mloc,Nloc1,SyL(1:Mloc,1:Nloc1),SyR(1:Mloc,1:Nloc1),  &
              FyL(1:Mloc,1:Nloc1),FyR(1:Mloc,1:Nloc1),HUyL(1:Mloc,1:Nloc1),  &
			HUyR(1:Mloc,1:Nloc1),Fy(1:Mloc,1:Nloc1))
     CALL HLL(Mloc1,Nloc,SxL(1:Mloc1,1:Nloc),SxR(1:Mloc1,1:Nloc),  &
              GxL(1:Mloc1,1:Nloc),GxR(1:Mloc1,1:Nloc),HVxL(1:Mloc1,1:Nloc),  &
		    HVxR(1:Mloc1,1:Nloc),Gx(1:Mloc1,1:Nloc))
     CALL HLL(Mloc,Nloc1,SyL(1:Mloc,1:Nloc1),SyR(1:Mloc,1:Nloc1),  &
              GyL(1:Mloc,1:Nloc1),GyR(1:Mloc,1:Nloc1),HVyL(1:Mloc,1:Nloc1),  &
			HVyR(1:Mloc,1:Nloc1),Gy(1:Mloc,1:Nloc1))

END SUBROUTINE FLUX_AT_INTERFACE_HLL

!-------------------------------------------------------------------------------------
!
!    HLL is subroutine for HLL scheme
!
!    HISTORY: 
!      05/06/2010 Fengyan Shi
!
!-------------------------------------------------------------------------------------
SUBROUTINE HLL(M,N,SL,SR,FL,FR,UL,UR,FOUT)
     USE PARAM
     INTEGER,INTENT(IN)::M,N
     REAL(SP),INTENT(IN),DIMENSION(M,N)::SL,SR,FL,FR,UL,UR
     REAL(SP),INTENT(OUT),DIMENSION(M,N)::FOUT

      DO J=1,N
      DO I=1,M     
      IF(SL(I,J)>=ZERO) THEN
        FOUT(I,J)=FL(I,J)
      ELSEIF(SR(I,J)<=ZERO) THEN
        FOUT(I,J)=FR(I,J)
      ELSE
        FOUT(I,J)=SR(I,J)*FL(I,J)-SL(I,J)*FR(I,J)+SL(I,J)*SR(I,J)*(UR(I,J)-UL(I,J))
        IF((ABS(SR(I,J)-SL(I,J)))<SMALL)THEN
         FOUT(I,J)=FOUT(I,J)/SMALL
        ELSE
         FOUT(I,J)=FOUT(I,J)/(SR(I,J)-SL(I,J))
        ENDIF
      ENDIF
      ENDDO
      ENDDO

END SUBROUTINE HLL

!-------------------------------------------------------------------------------------
!
!    DelxyFun is subroutine to calculate derivative of x and y
!
!    HISTORY: 
!     05/06/2010 Fengyan Shi
!
!-------------------------------------------------------------------------------------
SUBROUTINE DelxyFun
     USE GLOBAL

! compute DelxFun


       

      CALL DelxFun(DX(1:Mloc,1:Nloc),Mloc,Nloc,Eta(1:Mloc,1:Nloc),DelxEtar(1:Mloc,1:Nloc))
      CALL DelxFun(DX(1:Mloc,1:Nloc),Mloc,Nloc,U(1:Mloc,1:Nloc),DelxU(1:Mloc,1:Nloc))
      CALL DelxFun(DX(1:Mloc,1:Nloc),Mloc,Nloc,V(1:Mloc,1:Nloc),DelxV(1:Mloc,1:Nloc))
      CALL DelxFun(DX(1:Mloc,1:Nloc),Mloc,Nloc,HU(1:Mloc,1:Nloc),DelxHU(1:Mloc,1:Nloc))
      CALL DelxFun(DX(1:Mloc,1:Nloc),Mloc,Nloc,HV(1:Mloc,1:Nloc),DelxHV(1:Mloc,1:Nloc))

! compute DelyFun
      CALL DelyFun(DY(1:Mloc,1:Nloc),Mloc,Nloc,Eta(1:Mloc,1:Nloc),DelyEtar(1:Mloc,1:Nloc))
      CALL DelyFun(DY(1:Mloc,1:Nloc),Mloc,Nloc,U(1:Mloc,1:Nloc),DelyU(1:Mloc,1:Nloc))
      CALL DelyFun(DY(1:Mloc,1:Nloc),Mloc,Nloc,V(1:Mloc,1:Nloc),DelyV(1:Mloc,1:Nloc))
      CALL DelyFun(DY(1:Mloc,1:Nloc),Mloc,Nloc,HV(1:Mloc,1:Nloc),DelyHV(1:Mloc,1:Nloc))
      CALL DelyFun(DY(1:Mloc,1:Nloc),Mloc,Nloc,HU(1:Mloc,1:Nloc),DelyHU(1:Mloc,1:Nloc))



END SUBROUTINE DelxyFun

!-------------------------------------------------------------------------------------
!
!    DelyFun is subroutine to calculate derivative of y
!
!    HISTORY: 
!      05/06/2010 Fengyan Shi
!
!-------------------------------------------------------------------------------------
SUBROUTINE DelyFun(DY,M,N,DIN,DOUT)
     USE PARAM
     IMPLICIT NONE



     REAL(SP),DIMENSION(M,N),INTENT(IN)::DY

     INTEGER,INTENT(IN)::M,N
     REAL(SP),INTENT(IN),DIMENSION(M,N)::DIN
     REAL(SP),INTENT(OUT),DIMENSION(M,N)::DOUT
!     REAL(SP) :: VANLEER_LIMITER

! van Leer Limiter     
     DO I=1,M
     DO J=2,N-1




      TMP1=(DIN(I,J+1)-DIN(I,J))/DY(I,J)
      TMP2=(DIN(I,J)-DIN(I,J-1))/DY(I,J-1)

      IF((ABS(TMP1)+ABS(TMP2))<SMALL)THEN
!       !WRITE(*,*)'WARNING É (a+b) in limiter x is too small','(I,J)',I,J
        DOUT(I,J)=ZERO
      ELSE
!        DOUT(I,J)=(TMP1*ABS(TMP2)+ABS(TMP1)*TMP2)/(ABS(TMP1)+ABS(TMP2))
!        DOUT(I,J)=VANLEER_LIMITER(TMP1,TMP2)
         DOUT(I,J)=(TMP1*ABS(TMP2)+ABS(TMP1)*TMP2)/(ABS(TMP1)+ABS(TMP2))
      ENDIF
     ENDDO
     ENDDO

     DO I=1,M




      DOUT(I,1)=(DIN(I,2)-DIN(I,1))/DY(I,1)
      DOUT(I,N)=(DIN(I,N)-DIN(I,N-1))/DY(I,N)

     ENDDO

END SUBROUTINE DelyFun

!-------------------------------------------------------------------------------------
!
!    DelxFun is subroutine to calculate derivative of x
!
!    HISTORY: 
!     05/06/2010 Fengyan Shi
!
!-------------------------------------------------------------------------------------
SUBROUTINE DelxFun(DX,M,N,DIN,DOUT)
     USE PARAM
     IMPLICIT NONE



     REAL(SP),DIMENSION(M,N),INTENT(IN)::DX

     INTEGER,INTENT(IN)::M,N
     REAL(SP),INTENT(IN),DIMENSION(M,N)::DIN
     REAL(SP),INTENT(OUT),DIMENSION(M,N)::DOUT
!     REAL(SP) :: VANLEER_LIMITER

! van Leer Limiter     
     DO I=2,M-1
     DO J=1,N




      TMP1=(DIN(I+1,J)-DIN(I,J))/DX(I,J)
      TMP2=(DIN(I,J)-DIN(I-1,J))/DX(I-1,J)

      IF((ABS(TMP1)+ABS(TMP2))<SMALL)THEN
!       !WRITE(*,*)'WARNING É (a+b) in limiter y is too small','(I,J)',I,J
        DOUT(I,J)=ZERO
      ELSE
!        DOUT(I,J)=VANLEER_LIMITER(TMP1,TMP2)
         DOUT(I,J)=(TMP1*ABS(TMP2)+ABS(TMP1)*TMP2)/(ABS(TMP1)+ABS(TMP2))
      ENDIF
     ENDDO
     ENDDO

     DO J=1,N




      DOUT(1,J)=(DIN(2,J)-DIN(1,J))/DX(1,J)
      DOUT(M,J)=(DIN(M,J)-DIN(M-1,J))/DX(M,J)

     ENDDO

END SUBROUTINE DelxFun

! ---------------------------------------------------
!
!    WAVE_SPEED is subroutine to calculate wave speed
!    no shear wave is calculated yet
!
!    HISTORY: 
!        01/21/2012 Fengyan Shi
!        10/29/2012 Fengyan Shi
!             Steve Brandt mentioned HxL and HxR not calculated 
!             inside ghost cells, it is not an issue though corrected
!             There's a bug found in Dmitry's case. In x direction, 
!             should use M1 and N1 for y direction.
!
! --------------------------------------------------
SUBROUTINE WAVE_SPEED(M,N,M1,N1,UL,UR,VL,VR,HxL,HxR,HyL,HyR,&
      SxL,SxR,SyL,SyR)
     USE PARAM
     USE GLOBAL, ONLY : Nghost
     IMPLICIT NONE
     INTEGER,INTENT(IN)::M,N,M1,N1
     REAL(SP),INTENT(IN),DIMENSION(M1,N)::UL,UR,HxL,HxR
     REAL(SP),INTENT(IN),DIMENSION(M,N1)::VL,VR,HyL,HyR
     REAL(SP),INTENT(OUT),DIMENSION(M1,N)::SxL,SxR
     REAL(SP),INTENT(OUT),DIMENSION(M,N1)::SyL,SyR         
     REAL(SP)::SQR_PHI_L,SQR_PHI_R,SQR_PHI_S,U_S


! Zhou et al., 2001 approach
! x interface
!     DO J=1,N
!     DO I=1,M
     DO J=1+Nghost,N-Nghost
     DO I=1+Nghost,M1-Nghost
       SQR_PHI_L=SQRT(GRAV*ABS(HxL(I,J)))
       SQR_PHI_R=SQRT(GRAV*ABS(HxR(I,J)))
       SQR_PHI_S=0.5*(SQR_PHI_L+SQR_PHI_R)+0.25*(UL(I,J)-UR(I,J))  
       U_S=0.5*(UL(I,J)+UR(I,J))+SQR_PHI_L-SQR_PHI_R
       SxL(I,J)=MIN(UL(I,J)-SQR_PHI_L,U_S-SQR_PHI_S)
       SxR(I,J)=MAX(UR(I,J)+SQR_PHI_R,U_S+SQR_PHI_S)
     ENDDO
     ENDDO

! ghost cells, this does not really matter
     DO J=1+Nghost,N-Nghost
      DO I=1,Nghost
       SxL(I,J)=SxL(Nghost+1,J)
       SxR(I,J)=SxR(Nghost+1,J)
      ENDDO
      DO I=M1-Nghost+1,M1
       SxL(I,J)=SxL(M1-Nghost,J)
       SxR(I,J)=SxR(M1-Nghost,J)       
      ENDDO
     ENDDO

     DO I=1,M1
       DO J=1,Nghost
         SxL(I,J)=SxL(I,Nghost+1)
         SxR(I,J)=SxR(I,Nghost+1)
       ENDDO
       DO J=N-Nghost+1,N
         SxL(I,J)=SxL(I,N-Nghost)
         SxR(I,J)=SxR(I,N-Nghost)
       ENDDO
     ENDDO

! y interface
     DO J=1+Nghost,N1-Nghost
     DO I=1+Nghost,M-Nghost
       SQR_PHI_L=SQRT(GRAV*ABS(HyL(I,J)))
       SQR_PHI_R=SQRT(GRAV*ABS(HyR(I,J)))
       SQR_PHI_S=0.5*(SQR_PHI_L+SQR_PHI_R)+0.25*(VL(I,J)-VR(I,J))  
       U_S=0.5*(VL(I,J)+VR(I,J))+SQR_PHI_L-SQR_PHI_R
       SyL(I,J)=MIN(VL(I,J)-SQR_PHI_L,U_S-SQR_PHI_S)
       SyR(I,J)=MAX(VR(I,J)+SQR_PHI_R,U_S+SQR_PHI_S)
     ENDDO
     ENDDO


! ghost cells
     DO I=1+Nghost,M-Nghost
      DO J=1,Nghost
       SyL(I,J)=SyL(I,Nghost+1)
       SyR(I,J)=SyR(I,Nghost+1)
      ENDDO
      DO J=N1-Nghost+1,N1
       SyL(I,J)=SyL(I,N1-Nghost)
       SyR(I,J)=SyR(I,N1-Nghost)       
      ENDDO
     ENDDO

     DO J=1,N1
       DO I=1,Nghost
         SyL(I,J)=SyL(Nghost+1,J)
         SyR(I,J)=SyR(Nghost+1,J)
       ENDDO
       DO I=M-Nghost+1,M
         SyL(I,J)=SyL(M-Nghost,J)
         SyR(I,J)=SyR(M-Nghost,J)
       ENDDO
     ENDDO

END SUBROUTINE WAVE_SPEED

!-------------------------------------------------------------------------------------
!
!    CONSTRUCTION is subroutine interface construction
!
!    HISTORY: 
!    05/06/2010 Fengyan Shi
!    08/06/2015 Young-Kwang Choi, added DelxU4 and DelxU4 terms
!    09/07/2015 Young-Kwang Choi, added MASK9u and MAST9v
!
!-------------------------------------------------------------------------------------
SUBROUTINE CONSTRUCTION
     USE GLOBAL
     IMPLICIT NONE




        

![ykchoi	(09.07.2015)
     REAL(SP),ALLOCATABLE :: MASK9u(:,:), MASK9v(:,:)
     ALLOCATE( MASK9u(1:Mloc1,1:Nloc), MASK9v(1:Mloc,1:Nloc1) )
	
     MASK9u(1:Mloc,1:Nloc)=MASK9(1:Mloc,1:Nloc)
     MASK9v(1:Mloc,1:Nloc)=MASK9(1:Mloc,1:Nloc)
     
     MASK9u(Mloc1,1:Nloc)=MASK9(Mloc,1:Nloc)
     MASK9v(1:Mloc,Nloc1)=MASK9(1:Mloc,Nloc)
!ykchoi]

! construct in x-direction

       
     CALL CONSTRUCT_X(Mloc,Mloc1,Nloc,DX(1:Mloc,1:Nloc),U(1:Mloc,1:Nloc),DelxU(1:Mloc,1:Nloc),  &
                      UxL(1:Mloc1,1:Nloc),UxR(1:Mloc1,1:Nloc),Kappa)
     CALL CONSTRUCT_X(Mloc,Mloc1,Nloc,DX(1:Mloc,1:Nloc),V(1:Mloc,1:Nloc),DelxV(1:Mloc,1:Nloc),  &
                      VxL(1:Mloc1,1:Nloc),VxR(1:Mloc1,1:Nloc),Kappa)
     CALL CONSTRUCT_X(Mloc,Mloc1,Nloc,DX(1:Mloc,1:Nloc),HU(1:Mloc,1:Nloc),DelxHU(1:Mloc,1:Nloc),  &
                      HUxL(1:Mloc1,1:Nloc),HUxR(1:Mloc1,1:Nloc),Kappa)
     CALL CONSTRUCT_X(Mloc,Mloc1,Nloc,DX(1:Mloc,1:Nloc),HV(1:Mloc,1:Nloc),DelxHV(1:Mloc,1:Nloc),  &
                      HVxL(1:Mloc1,1:Nloc),HVxR(1:Mloc1,1:Nloc),Kappa)
     CALL CONSTRUCT_X(Mloc,Mloc1,Nloc,DX(1:Mloc,1:Nloc),Eta(1:Mloc,1:Nloc),DelxEtar(1:Mloc,1:Nloc),  &
                      EtaRxL(1:Mloc1,1:Nloc),EtaRxR(1:Mloc1,1:Nloc),Kappa)

     HxL=EtaRxL+Depthx
     HxR=EtaRxR+Depthx

        

!ykchoi (09.07.2015)
!Following are modified for Ghost cell at Mloc1, Nloc1




     PL=HUxL 
     PR=HUxR






      FxL=Gamma3*PL*(UxL) &
        +0.5*GRAV*((EtaRxL)*(EtaRxL)*Gamma3+2.0_SP*(EtaRxL)*(Depthx))
      FxR=Gamma3*PR*(UxR) &
        +0.5*GRAV*((EtaRxR)*(EtaRxR)*Gamma3+2.0_SP*(EtaRxR)*(Depthx))





     GxL=HxL*UxL*VxL*Gamma3
     GxR=HxR*UxR*VxR*Gamma3




! construct in y-direction

       
     CALL CONSTRUCT_Y(Mloc,Nloc,Nloc1,DY(1:Mloc,1:Nloc),U(1:Mloc,1:Nloc),DelyU(1:Mloc,1:Nloc),  &
                      UyL(1:Mloc,1:Nloc1),UyR(1:Mloc,1:Nloc1),Kappa)
     CALL CONSTRUCT_Y(Mloc,Nloc,Nloc1,DY(1:Mloc,1:Nloc),V(1:Mloc,1:Nloc),DelyV(1:Mloc,1:Nloc),  &
                      VyL(1:Mloc,1:Nloc1),VyR(1:Mloc,1:Nloc1),Kappa)
     CALL CONSTRUCT_Y(Mloc,Nloc,Nloc1,DY(1:Mloc,1:Nloc),HV(1:Mloc,1:Nloc),DelyHV(1:Mloc,1:Nloc),  &
                      HVyL(1:Mloc,1:Nloc1),HVyR(1:Mloc,1:Nloc1),Kappa)
     CALL CONSTRUCT_Y(Mloc,Nloc,Nloc1,DY(1:Mloc,1:Nloc),HU(1:Mloc,1:Nloc),DelyHU(1:Mloc,1:Nloc),  &
                      HUyL(1:Mloc,1:Nloc1),HUyR(1:Mloc,1:Nloc1),Kappa)
     CALL CONSTRUCT_Y(Mloc,Nloc,Nloc1,DY(1:Mloc,1:Nloc),Eta(1:Mloc,1:Nloc),DelyEtar(1:Mloc,1:Nloc),  &
                      EtaRyL(1:Mloc,1:Nloc1),EtaRyR(1:Mloc,1:Nloc1),Kappa)

     HyL=EtaRyL+Depthy
     HyR=EtaRyR+Depthy






     QL=HVyL 
     QR=HVyR






     FyL=Gamma3*HyL*UyL*VyL*Gamma3
     FyR=Gamma3*HyR*UyR*VyR*Gamma3

! $$$ fyshi


     GyL=Gamma3*QL*(VyL) &
        +0.5*GRAV*((EtaRyL)*(EtaRyL)*Gamma3+2.0_SP*(EtaRyL)*(Depthy))
     GyR=Gamma3*QR*(VyR) &
        +0.5*GRAV*((EtaRyR)*(EtaRyR)*Gamma3+2.0_SP*(EtaRyR)*(Depthy))

        
	
	deallocate( MASK9u, MASK9v )
     
END SUBROUTINE CONSTRUCTION

!-------------------------------------------------------------------------------------
!
!    CONSTRUCT_Y is subroutine construct variable in Y direction
! 
!   HISTORY: 
!     06/20/2010 Fengyan Shi
!     06/01/2011 Jeff Harris, 
!                replaced limiter function, saving 1/3 time
!     08/31/2016 Choi removed kappa because 3rd order and 2nd order are identical
!
!-------------------------------------------------------------------------------------
SUBROUTINE CONSTRUCT_Y(M,N,N1,DY,Vin,Din,OutL,OutR,Kappa)
     USE PARAM
     IMPLICIT NONE
     INTEGER,INTENT(IN)::M,N,N1



     REAL(SP),DIMENSION(M,N),INTENT(IN)::DY

     REAL(SP),INTENT(IN),DIMENSION(M,N)::Vin,Din
     REAL(SP),INTENT(IN) :: Kappa
     REAL(SP),INTENT(OUT),DIMENSION(M,N1)::OutL,OutR


![ykchoi (08/28/2016)
! The "kappa" is meaningless in above routines.
! Thus, above routines can be simplified as follwoing.
! This is piecewise linear reconstruction using van Leer limiter (Zhou et al., 2001)
     DO I=1,M
        
        DO J=2,N



       
           OutL(I,J) = Vin(I,J-1) + 0.5_SP*DY(I,J-1)*Din(I,J-1)
           OutR(I,J) = Vin(I,J) - 0.5_SP*DY(I,J)*Din(I,J)

        ENDDO




       
        OutL(I,N1) = Vin(I,N) + 0.5_SP*DY(I,N)*Din(I,N)
        OutR(I,1) = Vin(I,1) - 0.5_SP*DY(I,1)*Din(I,1)

        OutL(I,1)=OutR(I,1)
        OutR(I,N1)=OutL(I,N1)

     ENDDO
!ykchoi]


END SUBROUTINE CONSTRUCT_Y

!-------------------------------------------------------------------------------------
!
!    CONSTRUCT_X is subroutine construct variable in X direction
!
!    HISTORY: 
!      06/20/2010 Fengyan Shi
!      06/01/2011 Jeff Harris, replaced limiter function, saving 1/3 time
! add third order in the code. 
!     08/31/2016 Choi removed kappa because 3rd order and 2nd order are identical
!
!-------------------------------------------------------------------------------------
SUBROUTINE CONSTRUCT_X(M,M1,N,DX,Vin,Din,OutL,OutR,Kappa)
     USE PARAM
     IMPLICIT NONE
     INTEGER,INTENT(IN)::M,M1,N



     REAL(SP),DIMENSION(M,N),INTENT(IN)::DX

     REAL(SP),INTENT(IN),DIMENSION(M,N)::Vin,Din
     REAL(SP),INTENT(IN)::Kappa
     REAL(SP),INTENT(OUT),DIMENSION(M1,N)::OutL,OutR

![ykchoi (08/28/2016)
! The "kappa" is meaningless in above routines.
! Thus, above routines can be simplified as follwoing.
! This is piecewise linear reconstruction using van Leer limiter (Zhou et al., 2001)
     DO J=1,N
        
	  DO I=2,M



       
          OutL(I,J) = Vin(I-1,J) + 0.5_SP*DX(I-1,J)*Din(I-1,J)
          OutR(I,J) = Vin(I,J) - 0.5_SP*DX(I,J)*Din(I,J)

        ENDDO



       
        OutL(M1,J) = Vin(M,J) + 0.5_SP*DX(M,J)*Din(M,J)
        OutR(1,J) = Vin(1,J) - 0.5_SP*DX(1,J)*Din(1,J)

        OutL(1,J) = OutR(1,J)
        OutR(M1,J) = OutL(M1,J)

     ENDDO
!ykchoi]

END SUBROUTINE CONSTRUCT_X

!-------------------------------------------------------------------------------------
!
!    CONSTRUCTION_HO is subroutine for 
!        high-order interface construction (Yamamoto et al.,1998)
!
!    HISTORY: 
!      09/07/2010 Fengyan Shi
!             use dummy variables
!      09/07/2015 Young-Kwang Choi, added MASK9u and MASK9v 
!      08/31/2016 Choi noticed 4th-order + minmod is not stable
!                 I changed construction_ho to _minmod to keep
!                 this subroutine for further tests
!
!-------------------------------------------------------------------------------------
SUBROUTINE CONSTRUCTION_HO_minmod
     USE GLOBAL
     IMPLICIT NONE



        
	
![ykchoi	(09.07.2015)
     REAL(SP),ALLOCATABLE :: MASK9u(:,:), MASK9v(:,:)
     ALLOCATE( MASK9u(1:Mloc1,1:Nloc), MASK9v(1:Mloc,1:Nloc1) )
	
     MASK9u(1:Mloc,1:Nloc)=MASK9(1:Mloc,1:Nloc)
     MASK9v(1:Mloc,1:Nloc)=MASK9(1:Mloc,1:Nloc)
     
     MASK9u(Mloc1,1:Nloc)=MASK9(Mloc,1:Nloc)
     MASK9v(1:Mloc,Nloc1)=MASK9(1:Mloc,Nloc)
!ykchoi]

! construct in x-direction

       
     CALL CONSTRUCT_HO_X(Mloc,Nloc,Mloc1,Ibeg,Iend,Jbeg,Jend,   &
           DX(1:Mloc,1:Nloc),MASK(1:Mloc,1:Nloc),U(1:Mloc,1:Nloc),UxL(1:Mloc1,1:Nloc),UxR(1:Mloc1,1:Nloc))
     CALL CONSTRUCT_HO_X(Mloc,Nloc,Mloc1,Ibeg,Iend,Jbeg,Jend,   &
           DX(1:Mloc,1:Nloc),MASK(1:Mloc,1:Nloc),V(1:Mloc,1:Nloc),VxL(1:Mloc1,1:Nloc),VxR(1:Mloc1,1:Nloc))
     CALL CONSTRUCT_HO_X(Mloc,Nloc,Mloc1,Ibeg,Iend,Jbeg,Jend,   &
           DX(1:Mloc,1:Nloc),MASK(1:Mloc,1:Nloc),HU(1:Mloc,1:Nloc),HUxL(1:Mloc1,1:Nloc),HUxR(1:Mloc1,1:Nloc))
     CALL CONSTRUCT_HO_X(Mloc,Nloc,Mloc1,Ibeg,Iend,Jbeg,Jend,   &
           DX(1:Mloc,1:Nloc),MASK(1:Mloc,1:Nloc),HV(1:Mloc,1:Nloc),HVxL(1:Mloc1,1:Nloc),HVxR(1:Mloc1,1:Nloc))
     CALL CONSTRUCT_HO_X(Mloc,Nloc,Mloc1,Ibeg,Iend,Jbeg,Jend,   &
           DX(1:Mloc,1:Nloc),MASK(1:Mloc,1:Nloc),Eta(1:Mloc,1:Nloc),EtaRxL(1:Mloc1,1:Nloc),EtaRxR(1:Mloc1,1:Nloc))



! dispersion
     HxL=EtaRxL+Depthx
     HxR=EtaRxR+Depthx

     PL=HUxL
     PR=HUxR
     FxL=Gamma3*PL*UxL &
        +0.5*GRAV*((EtaRxL)*(EtaRxL)*Gamma3+2.0_SP*(EtaRxL)*(Depthx))
     FxR=Gamma3*PR*UxR &
        +0.5*GRAV*((EtaRxR)*(EtaRxR)*Gamma3+2.0_SP*(EtaRxR)*(Depthx))






     GxL=Gamma3*HxL*UxL*VxL
     GxR=Gamma3*HxR*UxR*VxR


      
! construct in y-direction

       
     CALL CONSTRUCT_HO_Y(Mloc,Nloc,Nloc1,Ibeg,Iend,Jbeg,Jend,  &
       DY(1:Mloc,1:Nloc),MASK(1:Mloc,1:Nloc),U(1:Mloc,1:Nloc),UyL(1:Mloc,1:Nloc1),UyR(1:Mloc,1:Nloc1))
     CALL CONSTRUCT_HO_Y(Mloc,Nloc,Nloc1,Ibeg,Iend,Jbeg,Jend,  &
       DY(1:Mloc,1:Nloc),MASK(1:Mloc,1:Nloc),V(1:Mloc,1:Nloc),VyL(1:Mloc,1:Nloc1),VyR(1:Mloc,1:Nloc1))
     CALL CONSTRUCT_HO_Y(Mloc,Nloc,Nloc1,Ibeg,Iend,Jbeg,Jend,  &
       DY(1:Mloc,1:Nloc),MASK(1:Mloc,1:Nloc),HV(1:Mloc,1:Nloc),HVyL(1:Mloc,1:Nloc1),HVyR(1:Mloc,1:Nloc1))
     CALL CONSTRUCT_HO_Y(Mloc,Nloc,Nloc1,Ibeg,Iend,Jbeg,Jend,  &
       DY(1:Mloc,1:Nloc),MASK(1:Mloc,1:Nloc),HU(1:Mloc,1:Nloc),HUyL(1:Mloc,1:Nloc1),HUyR(1:Mloc,1:Nloc1))
     CALL CONSTRUCT_HO_Y(Mloc,Nloc,Nloc1,Ibeg,Iend,Jbeg,Jend,  &
       DY(1:Mloc,1:Nloc),MASK(1:Mloc,1:Nloc),Eta(1:Mloc,1:Nloc),EtaRyL(1:Mloc,1:Nloc1),EtaRyR(1:Mloc,1:Nloc1))



! dispersion
     HyL=EtaRyL+Depthy
     HyR=EtaRyR+Depthy



     QL=HVyL
     QR=HVyR
     GyL=Gamma3*QL*VyL &
        +0.5*GRAV*((EtaRyL)*(EtaRyL)*Gamma3+2.0_SP*(EtaRyL)*(Depthy))
     GyR=Gamma3*QR*VyR &
        +0.5*GRAV*((EtaRyR)*(EtaRyR)*Gamma3+2.0_SP*(EtaRyR)*(Depthy))

        




     FyL=Gamma3*HyL*UyL*VyL
     FyR=Gamma3*HyR*UyR*VyR


	
 deallocate( MASK9u, MASK9v )   !ykchoi 
     
END SUBROUTINE CONSTRUCTION_HO_minmod

!-------------------------------------------------------------------------------------
!
!    CONSTRUCT_HO_X is subroutine high-order construct variable in X direction
!
!    HISTORY: 
!    06/20/2010 Fengyan Shi, University of Delaware
!    02/14/2011 Jeff Harris, replaced limiter function, saving 1/3 time
!      08/31/2016 Choi noticed 4th-order + minmod is not stable
!                 I changed construction_ho to _minmod to keep
!                 this subroutine for further tests
!
!-------------------------------------------------------------------------------------
SUBROUTINE CONSTRUCT_HO_X_minmod(M,N,M1,Ibeg,Iend,Jbeg,Jend,DX,MASK,Vin,OutL,OutR)
     USE PARAM
     IMPLICIT NONE
     INTEGER, INTENT(IN)::M,N,M1,Ibeg,Iend,Jbeg,Jend



     REAL(SP),DIMENSION(M,N),INTENT(IN)::DX

     REAL(SP),INTENT(IN),DIMENSION(M,N)::Vin
     INTEGER, INTENT(IN),DIMENSION(M,N):: MASK
     REAL(SP),INTENT(OUT),DIMENSION(M1,N)::OutL,OutR

     REAL(SP),DIMENSION(M,N) :: Din
     REAL(SP) :: TXP1,TXP2,TXP3,TXP4,DVP1,DVP2,DVP3

     ! estimate Din first
     Din=0.0_SP
     DO J=Jbeg,Jend
     DO I=Ibeg-1,Iend+2
       TXP1=Vin(I-1,J)-Vin(I-2,J)
       TXP2=Vin(I,J)-Vin(I-1,J)
       TXP3=Vin(I+1,J)-Vin(I,J)

       if (TXP1.ge.0.0_SP) then
          DVP1=MAX(0.0_SP,MIN(TXP1,2.0*TXP2,2.0*TXP3))          
       else
          DVP1=MIN(0.0_SP,MAX(TXP1,2.0*TXP2,2.0*TXP3))          
       endif
       if (TXP2.ge.0.0_SP) then
          DVP2=MAX(0.0_SP,MIN(TXP2,2.0*TXP3,2.0*TXP1))          
       else
          DVP2=MIN(0.0_SP,MAX(TXP2,2.0*TXP3,2.0*TXP1))          
       endif
       if (TXP3.ge.0.0_SP) then
          DVP3=MAX(0.0_SP,MIN(TXP3,2.0*TXP1,2.0*TXP2))          
       else
          DVP3=MIN(0.0_SP,MAX(TXP3,2.0*TXP1,2.0*TXP2))          
       endif

! dry D-2 I-1 I D+1, lower-order
       IF(MASK(I-2,J)==0.OR.MASK(I+1,J)==0)THEN
       TXP2=Vin(I,J)-Vin(I-1,J)
       TXP1=TXP2
       TXP3=TXP2
       if (TXP1.ge.0.0_SP) then
          DVP1=MAX(0.0_SP,MIN(TXP1,2.0*TXP2,2.0*TXP3))          
       else
          DVP1=MIN(0.0_SP,MAX(TXP1,2.0*TXP2,2.0*TXP3))          
       endif
       if (TXP2.ge.0.0_SP) then
          DVP2=MAX(0.0_SP,MIN(TXP2,2.0*TXP3,2.0*TXP1))          
       else
          DVP2=MIN(0.0_SP,MAX(TXP2,2.0*TXP3,2.0*TXP1))          
       endif
       if (TXP3.ge.0.0_SP) then
          DVP3=MAX(0.0_SP,MIN(TXP3,2.0*TXP1,2.0*TXP2))          
       else
          DVP3=MIN(0.0_SP,MAX(TXP3,2.0*TXP1,2.0*TXP2))          
       endif
! here actually DVP1=DVP2=DVP3=TXP2
       ENDIF    
! dry I-2 D-1 D I+1, zero gradient
       IF(MASK(I-1,J)==0.OR.MASK(I,J)==0)THEN
       DVP1=ZERO
       DVP2=ZERO
       DVP3=ZERO
       ENDIF
       
       Din(I,J)=TXP2-1.0_SP/6.0_SP*(DVP3-2.0_SP*DVP2+DVP1)
     ENDDO
     ENDDO 

     DO J=Jbeg,Jend
     DO I=Ibeg,Iend+1
! Jeff modified the following statements 02/14/2011
       if (Din(I-1,J).ge.0.0_SP) then
          TXP1=MAX(0.0_SP,MIN(Din(I-1,J),4.0_SP*Din(I,J)))
       else
          TXP1=MIN(0.0_SP,MAX(Din(I-1,J),4.0_SP*Din(I,J)))
       endif
       if (Din(I,J).ge.0.0_SP) then
          TXP2=MAX(0.0_SP,MIN(Din(I,J),4.0_SP*Din(I-1,J)))
       else
          TXP2=MIN(0.0_SP,MAX(Din(I,J),4.0_SP*Din(I-1,J)))
       endif
! there was a HUGE bug here, 12 should versus 43, fixed. fyshi
       if (Din(I,J).ge.0.0_SP) then
          TXP4=MAX(0.0_SP,MIN(Din(I,J),4.0_SP*Din(I+1,J)))
       else
          TXP4=MIN(0.0_SP,MAX(Din(I,J),4.0_SP*Din(I+1,J)))
       endif
       if (Din(I+1,J).ge.0.0_SP) then
          TXP3=MAX(0.0_SP,MIN(Din(I+1,J),4.0_SP*Din(I,J)))
       else
          TXP3=MIN(0.0_SP,MAX(Din(I+1,J),4.0_SP*Din(I,J)))
       endif
       
       OutL(I,J)=Vin(I-1,J)+1.0_SP/6.0_SP*(TXP1+2.0_SP*TXP2)
       OutR(I,J)=Vin(I,J)-1.0_SP/6.0_SP*(TXP3+2.0_SP*TXP4)
     ENDDO
     ENDDO

END SUBROUTINE CONSTRUCT_HO_X_minmod


!-------------------------------------------------------------------------------------
!
!    CONSTRUCT_HO_Y is subroutine high-order construct variable in Y direction
!
!    HISTORY: 
!      06/20/2010 Fengyan Shi
!      02/14/2011 Jeff Harris, replaced limiter function, saving 1/3 time
!      08/31/2016 Choi noticed 4th-order + minmod is not stable
!                 I changed construction_ho to _minmod to keep
!                 this subroutine for further tests
!
!-------------------------------------------------------------------------------------
SUBROUTINE CONSTRUCT_HO_Y_minmod(M,N,N1,Ibeg,Iend,Jbeg,Jend,DY,MASK,Vin,OutL,OutR)
     USE PARAM
     IMPLICIT NONE
     INTEGER, INTENT(IN)::M,N,N1,Ibeg,Iend,Jbeg,Jend



     REAL(SP),DIMENSION(M,N),INTENT(IN)::DY

     REAL(SP),INTENT(IN),DIMENSION(M,N)::Vin
     INTEGER, INTENT(IN),DIMENSION(M,N):: MASK
     REAL(SP),INTENT(OUT),DIMENSION(M,N1)::OutL,OutR

     REAL(SP),DIMENSION(M,N) :: Din
     REAL(SP) :: TYP1,TYP2,TYP3,TYP4,DVP1,DVP2,DVP3

     ! estimate Din first
     Din=0.0_SP
     DO J=Jbeg-1,Jend+2
     DO I=Ibeg,Iend
       TYP1=Vin(I,J-1)-Vin(I,J-2)
       TYP2=Vin(I,J)-Vin(I,J-1)
       TYP3=Vin(I,J+1)-Vin(I,J)

       if (TYP1.ge.0.0_SP) then
          DVP1=MAX(0.0_SP,MIN(TYP1,2.0*TYP2,2.0*TYP3))          
       else
          DVP1=MIN(0.0_SP,MAX(TYP1,2.0*TYP2,2.0*TYP3))          
       endif
       if (TYP2.ge.0.0_SP) then
          DVP2=MAX(0.0_SP,MIN(TYP2,2.0*TYP3,2.0*TYP1))          
       else
          DVP2=MIN(0.0_SP,MAX(TYP2,2.0*TYP3,2.0*TYP1))          
       endif
       if (TYP3.ge.0.0_SP) then
          DVP3=MAX(0.0_SP,MIN(TYP3,2.0*TYP1,2.0*TYP2))          
       else
          DVP3=MIN(0.0_SP,MAX(TYP3,2.0*TYP1,2.0*TYP2))          
       endif

! dry D-2 J-1 J D+1, lower-order
       IF(MASK(I,J-2)==0.OR.MASK(I,J+1)==0)THEN
       TYP2=Vin(I,J)-Vin(I,J-1)
       TYP1=TYP2
       TYP3=TYP2
       if (TYP1.ge.0.0_SP) then
          DVP1=MAX(0.0_SP,MIN(TYP1,2.0*TYP2,2.0*TYP3))          
       else
          DVP1=MIN(0.0_SP,MAX(TYP1,2.0*TYP2,2.0*TYP3))          
       endif
       if (TYP2.ge.0.0_SP) then
          DVP2=MAX(0.0_SP,MIN(TYP2,2.0*TYP3,2.0*TYP1))          
       else
          DVP2=MIN(0.0_SP,MAX(TYP2,2.0*TYP3,2.0*TYP1))          
       endif
       if (TYP3.ge.0.0_SP) then
          DVP3=MAX(0.0_SP,MIN(TYP3,2.0*TYP1,2.0*TYP2))          
       else
          DVP3=MIN(0.0_SP,MAX(TYP3,2.0*TYP1,2.0*TYP2))          
       endif
! here actually DVP1=DVP2=DVP3=TYP2
       ENDIF    
! dry J-2 D-1 D J+1, zero gradient
       IF(MASK(I,J-1)==0.OR.MASK(I,J)==0)THEN
       DVP1=ZERO
       DVP2=ZERO
       DVP3=ZERO
       ENDIF

       Din(I,J)=TYP2-1.0_SP/6.0_SP*(DVP3-2.0_SP*DVP2+DVP1)
     ENDDO
     ENDDO     


     DO J=Jbeg,Jend+1
     DO I=Ibeg,Iend
! Jeff modified the following statements 02/14/2011
       if (Din(I,J-1).ge.0.0_SP) then
          TYP1=MAX(0.0_SP,MIN(Din(I,J-1),4.0_SP*Din(I,J)))
       else
          TYP1=MIN(0.0_SP,MAX(Din(I,J-1),4.0_SP*Din(I,J)))
       endif
       if (Din(I,J).ge.0.0_SP) then
          TYP2=MAX(0.0_SP,MIN(Din(I,J),4.0_SP*Din(I,J-1)))
       else
          TYP2=MIN(0.0_SP,MAX(Din(I,J),4.0_SP*Din(I,J-1)))
       endif
! there was a HUGE bug here, 12 should versus 43, fixed. fyshi
       if (Din(I,J).ge.0.0_SP) then
          TYP4=MAX(0.0_SP,MIN(Din(I,J),4.0_SP*Din(I,J+1)))
       else
          TYP4=MIN(0.0_SP,MAX(Din(I,J),4.0_SP*Din(I,J+1)))
       endif
       if (Din(I,J+1).ge.0.0_SP) then
          TYP3=MAX(0.0_SP,MIN(Din(I,J+1),4.0_SP*Din(I,J)))
       else
          TYP3=MIN(0.0_SP,MAX(Din(I,J+1),4.0_SP*Din(I,J)))
       endif

       OutL(I,J)=Vin(I,J-1)+1.0_SP/6.0_SP*(TYP1+2.0_SP*TYP2)
       OutR(I,J)=Vin(I,J)-1.0_SP/6.0_SP*(TYP3+2.0_SP*TYP4)
     ENDDO
     ENDDO     

END SUBROUTINE CONSTRUCT_HO_Y_minmod


!-------------------------------------------------------------------------------------
!    The following two functions are Van Leer and Minmod limiters.
!
!    HISTORY:   
!     05/27/2010 Gangfeng Ma
!     09/24/2012 Fengyan Shi
!    NOTE: some compiler complains about function explicit interface
!    put real in front of function
!-------------------------------------------------------------------------------------

REAL(SP) FUNCTION VANLEER_LIMITER(A,B)
    USE PARAM
    IMPLICIT NONE
    REAL(SP),INTENT(IN) :: A
    REAL(SP),OPTIONAL,INTENT(IN) :: B
!    REAL(SP) :: VANLEER_LIMITER

    IF(PRESENT(B)) THEN
      VANLEER_LIMITER=(A*ABS(B)+ABS(A)*B)/(ABS(A)+ABS(B))
    ELSE  
      VANLEER_LIMITER=(A+ABS(A))/(1.0+A)
    ENDIF

    RETURN
END FUNCTION VANLEER_LIMITER

FUNCTION MINMOD_LIMITER(A,B,C)
    USE PARAM
    IMPLICIT NONE
    REAL(SP),INTENT(IN) :: A,B
    REAL(SP),OPTIONAL,INTENT(IN) :: C
    REAL(SP) :: MINMOD_LIMITER

    IF(PRESENT(C)) THEN
      MINMOD_LIMITER=SIGN(1.0_SP,A)*MAX(0.0_SP,MIN(ABS(A),SIGN(1.0_SP,A)*B,SIGN(1.0_SP,A)*C))
    ELSE
      MINMOD_LIMITER=SIGN(1.0_SP,A)*MAX(0.0_SP,MIN(ABS(A),SIGN(1.0_SP,A)*B))
    ENDIF

    RETURN
END FUNCTION MINMOD_LIMITER


! --------------- ykchoi (08/28/2016)
! This is subroutine for the fourth order MUSCL-TVD scheme.
! Van-Leer limiter is used for the third-order part, 
! and the minmod limiter is used for the fourth order part in the fourth order scheme.
!
! Reference
! Hybrid finite-volume finite-difference scheme for the solution of Boussinesq equations
! Erduran et al. (2005)
!   08/31/21016  fyshi changed the subroutine name HO_vanleer to HO
!                to make a consistency
! ---------------------------------------------------------------------------------------
SUBROUTINE CONSTRUCTION_HO
     USE GLOBAL
     IMPLICIT NONE

     REAL(SP),ALLOCATABLE :: MASK9u(:,:), MASK9v(:,:)
     ALLOCATE( MASK9u(1:Mloc1,1:Nloc), MASK9v(1:Mloc,1:Nloc1) )
	
     MASK9u(1:Mloc,1:Nloc)=MASK9(1:Mloc,1:Nloc)
     MASK9v(1:Mloc,1:Nloc)=MASK9(1:Mloc,1:Nloc)
     
     MASK9u(Mloc1,1:Nloc)=MASK9(Mloc,1:Nloc)
     MASK9v(1:Mloc,Nloc1)=MASK9(1:Mloc,Nloc)

! construct in x-direction
     CALL CONSTRUCT_HO_X(Mloc,Nloc,Mloc1,Ibeg,Iend,Jbeg,Jend,  &
          MASK(1:Mloc,1:Nloc),U(1:Mloc,1:Nloc),UxL(1:Mloc1,1:Nloc),UxR(1:Mloc1,1:Nloc))
     CALL CONSTRUCT_HO_X(Mloc,Nloc,Mloc1,Ibeg,Iend,Jbeg,Jend,  &
          MASK(1:Mloc,1:Nloc),V(1:Mloc,1:Nloc),VxL(1:Mloc1,1:Nloc),VxR(1:Mloc1,1:Nloc))
     CALL CONSTRUCT_HO_X(Mloc,Nloc,Mloc1,Ibeg,Iend,Jbeg,Jend,  &
          MASK(1:Mloc,1:Nloc),HU(1:Mloc,1:Nloc),HUxL(1:Mloc1,1:Nloc),HUxR(1:Mloc1,1:Nloc))
     CALL CONSTRUCT_HO_X(Mloc,Nloc,Mloc1,Ibeg,Iend,Jbeg,Jend,  &
          MASK(1:Mloc,1:Nloc),HV(1:Mloc,1:Nloc),HVxL(1:Mloc1,1:Nloc),HVxR(1:Mloc1,1:Nloc))
     CALL CONSTRUCT_HO_X(Mloc,Nloc,Mloc1,Ibeg,Iend,Jbeg,Jend,  &
          MASK(1:Mloc,1:Nloc),Eta(1:Mloc,1:Nloc),EtaRxL(1:Mloc1,1:Nloc),EtaRxR(1:Mloc1,1:Nloc))

     HxL=EtaRxL+Depthx
     HxR=EtaRxR+Depthx


       

       
     PL=HUxL
     PR=HUxR
     FxL=Gamma3*PL*UxL &
        +0.5*GRAV*((EtaRxL)*(EtaRxL)*Gamma3+2.0_SP*(EtaRxL)*(Depthx))
     FxR=Gamma3*PR*UxR &
        +0.5*GRAV*((EtaRxR)*(EtaRxR)*Gamma3+2.0_SP*(EtaRxR)*(Depthx))





     GxL=Gamma3*HxL*UxL*VxL
     GxR=Gamma3*HxR*UxR*VxR


      
! construct in y-direction
     CALL CONSTRUCT_HO_Y(Mloc,Nloc,Nloc1,Ibeg,Iend,Jbeg,Jend,   &
          MASK(1:Mloc,1:Nloc),U(1:Mloc,1:Nloc),UyL(1:Mloc,1:Nloc1),UyR(1:Mloc,1:Nloc1))
     CALL CONSTRUCT_HO_Y(Mloc,Nloc,Nloc1,Ibeg,Iend,Jbeg,Jend,   &
          MASK(1:Mloc,1:Nloc),V(1:Mloc,1:Nloc),VyL(1:Mloc,1:Nloc1),VyR(1:Mloc,1:Nloc1))
     CALL CONSTRUCT_HO_Y(Mloc,Nloc,Nloc1,Ibeg,Iend,Jbeg,Jend,   &
          MASK(1:Mloc,1:Nloc),HV(1:Mloc,1:Nloc),HVyL(1:Mloc,1:Nloc1),HVyR(1:Mloc,1:Nloc1))
     CALL CONSTRUCT_HO_Y(Mloc,Nloc,Nloc1,Ibeg,Iend,Jbeg,Jend,   &
          MASK(1:Mloc,1:Nloc),HU(1:Mloc,1:Nloc),HUyL(1:Mloc,1:Nloc1),HUyR(1:Mloc,1:Nloc1))
     CALL CONSTRUCT_HO_Y(Mloc,Nloc,Nloc1,Ibeg,Iend,Jbeg,Jend,   &
          MASK(1:Mloc,1:Nloc),Eta(1:Mloc,1:Nloc),EtaRyL(1:Mloc,1:Nloc1),EtaRyR(1:Mloc,1:Nloc1))

     HyL=EtaRyL+Depthy
     HyR=EtaRyR+Depthy


       


       
     QL=HVyL
     QR=HVyR
     GyL=Gamma3*QL*VyL &
        +0.5*GRAV*((EtaRyL)*(EtaRyL)*Gamma3+2.0_SP*(EtaRyL)*(Depthy))
     GyR=Gamma3*QR*VyR &
        +0.5*GRAV*((EtaRyR)*(EtaRyR)*Gamma3+2.0_SP*(EtaRyR)*(Depthy))





     FyL=Gamma3*HyL*UyL*VyL
     FyR=Gamma3*HyR*UyR*VyR



	deallocate( MASK9u, MASK9v )
     
ENDSUBROUTINE CONSTRUCTION_HO


! --------------- ykchoi (08/28/2016)
! TVD-MUSCL scheme by Erduran et al. (2005)
! ---------------------------------------------------------------------------------------
SUBROUTINE CONSTRUCT_HO_X(M,N,M1,Ibeg,Iend,Jbeg,Jend,MASK,Vin,OutL,OutR)
     USE PARAM
     IMPLICIT NONE
     INTEGER,INTENT(IN) :: M,N,M1,Ibeg,Iend,Jbeg,Jend

     REAL(SP),INTENT(IN),DIMENSION(M,N) :: Vin
     INTEGER,INTENT(IN),DIMENSION(M,N) :: MASK
     REAL(SP),INTENT(OUT),DIMENSION(M1,N) :: OutL,OutR

     REAL(SP),DIMENSION(M,N) :: Din
     REAL(SP) :: TXP1,TXP2,TXP3,DVP1,DVP2,DVP3

     REAL(SP) :: VAN1, VAN2, RAT

     Din=0.0_SP
     DO J=Jbeg,Jend
        
	  DO I=Ibeg-1,Iend+2
           TXP1=Vin(I-1,J)-Vin(I-2,J)
           TXP2=Vin(I,J)-Vin(I-1,J)
           TXP3=Vin(I+1,J)-Vin(I,J)

           if (TXP1.ge.0.0_SP) then
              DVP1=MAX(0.0_SP,MIN(TXP1,2.0*TXP2,2.0*TXP3))          
           else
              DVP1=MIN(0.0_SP,MAX(TXP1,2.0*TXP2,2.0*TXP3))          
           endif
           if (TXP2.ge.0.0_SP) then
              DVP2=MAX(0.0_SP,MIN(TXP2,2.0*TXP3,2.0*TXP1))          
           else
              DVP2=MIN(0.0_SP,MAX(TXP2,2.0*TXP3,2.0*TXP1))          
           endif
           if (TXP3.ge.0.0_SP) then
              DVP3=MAX(0.0_SP,MIN(TXP3,2.0*TXP1,2.0*TXP2))          
           else
              DVP3=MIN(0.0_SP,MAX(TXP3,2.0*TXP1,2.0*TXP2))          
           endif

           ! dry D-2 I-1 I D+1, lower-order
           IF(MASK(I-2,J)==0.OR.MASK(I+1,J)==0)THEN
             TXP2=Vin(I,J)-Vin(I-1,J)
             TXP1=TXP2
             TXP3=TXP2
             if (TXP1.ge.0.0_SP) then
                DVP1=MAX(0.0_SP,MIN(TXP1,2.0*TXP2,2.0*TXP3))
             else
                DVP1=MIN(0.0_SP,MAX(TXP1,2.0*TXP2,2.0*TXP3))
             endif
             if (TXP2.ge.0.0_SP) then
                DVP2=MAX(0.0_SP,MIN(TXP2,2.0*TXP3,2.0*TXP1))
             else
                DVP2=MIN(0.0_SP,MAX(TXP2,2.0*TXP3,2.0*TXP1))
             endif
             if (TXP3.ge.0.0_SP) then
                DVP3=MAX(0.0_SP,MIN(TXP3,2.0*TXP1,2.0*TXP2))          
             else
                DVP3=MIN(0.0_SP,MAX(TXP3,2.0*TXP1,2.0*TXP2))          
             endif
           ! here actually DVP1=DVP2=DVP3=TXP2
           ENDIF
           ! dry I-2 D-1 D I+1, zero gradient
           IF(MASK(I-1,J)==0.OR.MASK(I,J)==0)THEN
             DVP1=ZERO
             DVP2=ZERO
             DVP3=ZERO
           ENDIF
       
           Din(I,J)=TXP2-1.0_SP/6.0_SP*(DVP3-2.0_SP*DVP2+DVP1)
        ENDDO

        DO I=Ibeg,Iend+1
           TMP1 = Din(I-1,J);   TMP2 = Din(I,J);
	  
           IF( ABS(TMP1).le.SMALL ) TMP1 = SMALL*SIGN( 1.0_SP, TMP1 )
           IF( ABS(TMP2).le.SMALL ) TMP2 = SMALL*SIGN( 1.0_SP, TMP2 )
           RAT = TMP2/TMP1
           VAN1 = 0.0_SP;
           IF( abs(1.0 + RAT) > SMALL ) VAN1 = ( RAT + ABS(RAT) )/( 1.0_SP + RAT)

           RAT = TMP1/TMP2
           VAN2 = 0.0_SP;
           IF( abs(1.0 + RAT) > SMALL ) VAN2 = ( RAT + ABS(RAT) )/( 1.0_SP + RAT)

           OutL(I,J) = Vin(I-1,J) + ( 1.0_SP/6.0_SP )*( VAN1*TMP1 + 2.0_SP*VAN2*TMP2 )

	     !!!!!!!!!!!!!!!

           TMP1 = Din(I,J);   TMP2 = Din(I+1,J);

           IF( ABS(TMP1).le.SMALL ) TMP1 = SMALL*SIGN( 1.0_SP, TMP1 )
           IF( ABS(TMP2).le.SMALL ) TMP2 = SMALL*SIGN( 1.0_SP, TMP2 )
           RAT = TMP2/TMP1
           VAN1 = 0.0_SP;
           IF( abs(1.0 + RAT) > SMALL ) VAN1 = ( RAT + ABS(RAT) )/( 1.0_SP + RAT)

           RAT = TMP1/TMP2
           VAN2 = 0.0_SP;
           IF( abs(1.0 + RAT) > SMALL ) VAN2 = ( RAT + ABS(RAT) )/( 1.0_SP + RAT)

           OutR(I,J) = Vin(I,J) - ( 1.0_SP/6.0_SP )*( 2.0_SP*VAN1*TMP1 + VAN2*TMP2 )
        ENDDO

     ENDDO 

END SUBROUTINE CONSTRUCT_HO_X

! --------------- ykchoi (08/28/2016)
! TVD-MUSCL scheme by Erduran et al. (2005)
! ---------------------------------------------------------------------------------------
SUBROUTINE CONSTRUCT_HO_Y(M,N,N1,Ibeg,Iend,Jbeg,Jend,MASK,Vin,OutL,OutR)
     USE PARAM
     IMPLICIT NONE
     INTEGER,INTENT(IN) :: M,N,N1,Ibeg,Iend,Jbeg,Jend

     REAL(SP),INTENT(IN),DIMENSION(M,N) :: Vin
     INTEGER,INTENT(IN),DIMENSION(M,N) :: MASK
     REAL(SP),INTENT(OUT),DIMENSION(M,N1) :: OutL,OutR

     REAL(SP),DIMENSION(M,N) :: Din
     REAL(SP) :: TYP1,TYP2,TYP3,DVP1,DVP2,DVP3

     REAL(SP) :: VAN1, VAN2, RAT

     Din=0.0_SP
     DO I=Ibeg,Iend

     	  DO J=Jbeg-1,Jend+2
           TYP1=Vin(I,J-1)-Vin(I,J-2)
           TYP2=Vin(I,J)-Vin(I,J-1)
           TYP3=Vin(I,J+1)-Vin(I,J)

           if (TYP1.ge.0.0_SP) then
              DVP1=MAX(0.0_SP,MIN(TYP1,2.0*TYP2,2.0*TYP3))          
           else
              DVP1=MIN(0.0_SP,MAX(TYP1,2.0*TYP2,2.0*TYP3))          
           endif
           if (TYP2.ge.0.0_SP) then
              DVP2=MAX(0.0_SP,MIN(TYP2,2.0*TYP3,2.0*TYP1))          
           else
              DVP2=MIN(0.0_SP,MAX(TYP2,2.0*TYP3,2.0*TYP1))          
           endif
           if (TYP3.ge.0.0_SP) then
              DVP3=MAX(0.0_SP,MIN(TYP3,2.0*TYP1,2.0*TYP2))          
           else
              DVP3=MIN(0.0_SP,MAX(TYP3,2.0*TYP1,2.0*TYP2))          
           endif

           ! dry D-2 J-1 J D+1, lower-order
           IF(MASK(I,J-2)==0.OR.MASK(I,J+1)==0)THEN
             TYP2=Vin(I,J)-Vin(I,J-1)
             TYP1=TYP2
             TYP3=TYP2
             if (TYP1.ge.0.0_SP) then
                DVP1=MAX(0.0_SP,MIN(TYP1,2.0*TYP2,2.0*TYP3))
             else
                DVP1=MIN(0.0_SP,MAX(TYP1,2.0*TYP2,2.0*TYP3))
             endif
             if (TYP2.ge.0.0_SP) then
                DVP2=MAX(0.0_SP,MIN(TYP2,2.0*TYP3,2.0*TYP1))
             else
                DVP2=MIN(0.0_SP,MAX(TYP2,2.0*TYP3,2.0*TYP1))
             endif
             if (TYP3.ge.0.0_SP) then
                DVP3=MAX(0.0_SP,MIN(TYP3,2.0*TYP1,2.0*TYP2))
             else
                DVP3=MIN(0.0_SP,MAX(TYP3,2.0*TYP1,2.0*TYP2))
             endif
           ! here actually DVP1=DVP2=DVP3=TYP2
           ENDIF    
           ! dry J-2 D-1 D J+1, zero gradient
           IF(MASK(I,J-1)==0.OR.MASK(I,J)==0)THEN
             DVP1=ZERO
             DVP2=ZERO
             DVP3=ZERO
           ENDIF

           Din(I,J)=TYP2-1.0_SP/6.0_SP*(DVP3-2.0_SP*DVP2+DVP1)
     	  ENDDO

	  DO J=Jbeg,Jend+1
           TMP1 = Din(I,J-1);   TMP2 = Din(I,J);

           IF( ABS(TMP1).le.SMALL ) TMP1 = SMALL*SIGN( 1.0_SP, TMP1 )
           IF( ABS(TMP2).le.SMALL ) TMP2 = SMALL*SIGN( 1.0_SP, TMP2 )
           RAT = TMP2/TMP1	   
           VAN1 = 0.0_SP;
           IF( abs(1.0 + RAT) > SMALL ) VAN1 = ( RAT + ABS(RAT) )/( 1.0_SP + RAT)

           RAT = TMP1/TMP2
           VAN2 = 0.0_SP;
           IF( abs(1.0 + RAT) > SMALL ) VAN2 = ( RAT + ABS(RAT) )/( 1.0_SP + RAT)

           OutL(I,J) = Vin(I,J-1) + ( 1.0_SP/6.0_SP )*( VAN1*TMP1 + 2.0_SP*VAN2*TMP2 )

	   !!!!!!!!!!!!!!!

           TMP1 = Din(I,J);   TMP2 = Din(I,J+1);

           IF( ABS(TMP1).le.SMALL ) TMP1 = SMALL*SIGN( 1.0_SP, TMP1 )
           IF( ABS(TMP2).le.SMALL ) TMP2 = SMALL*SIGN( 1.0_SP, TMP2 )
           RAT = TMP2/TMP1
           VAN1 = 0.0_SP;
           IF( abs(1.0 + RAT) > SMALL ) VAN1 = ( RAT + ABS(RAT) )/( 1.0_SP + RAT)

           RAT = TMP1/TMP2
           VAN2 = 0.0_SP;
           IF( abs(1.0 + RAT) > SMALL ) VAN2 = ( RAT + ABS(RAT) )/( 1.0_SP + RAT)

           OutR(I,J) = Vin(I,J) - ( 1.0_SP/6.0_SP )*( 2.0_SP*VAN1*TMP1 + VAN2*TMP2 )
        ENDDO

     ENDDO     

END SUBROUTINE CONSTRUCT_HO_Y