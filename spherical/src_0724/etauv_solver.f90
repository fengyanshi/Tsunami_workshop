!------------------------------------------------------------------------------------
!
!      FILE etauv_solver.F
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
!    ESTIMATE_HUV is subroutine to calculate eta, ubar and vbar
!      using 3rd-order LK scheme 
!
!    HISTORY: 
!      05/12/2011 Fengyan Shi
!
!-------------------------------------------------------------------------------------
SUBROUTINE ESTIMATE_HUV(ISTEP,ng)
     USE GLOBAL
     IMPLICIT NONE
     INTEGER,INTENT(IN)::ISTEP
     REAL(SP),PARAMETER::n_left=-1.0_SP,n_right=1.0_SP,n_bottom=-1.0_SP,n_top=1.0_SP
     REAL(SP)::F_left,F_right,F_bottom,F_top,WK_Source
     REAL(SP),DIMENSION(Ibeg:Iend,Jbeg:Jend)::R1,R2,R3
! now work for spherical # if defined (CARTESIAN)
     REAL(SP)::xmk,ymk
! now work for spherical # endif
     REAL(SP)::DXg,DYg

     INTEGER::kf,kd
     INTEGER,INTENT(IN) :: ng

! MUSCL-Hancock, Zhou et al., p. 7





! only for wavemaker
     DXg=DX(1,1)
     DYg=DY(1,1)


! solve eta
     DO J=Jbeg,Jend
     DO I=Ibeg,Iend
      F_left=P(I,J)
      F_right=P(I+1,J)
      F_bottom=Q(I,J)
      F_top=Q(I,J+1)
! now work for spherical # if defined (CARTESIAN)
      IF(WAVEMAKER(1:6)=='WK_IRR')THEN

            xmk=(I-Ibeg)*DXg+npx*(Mloc-2*Nghost)*DXg
            ymk=(J-Jbeg)*DYg+npy*(Nloc-2*Nghost)*DYg




         IF(ABS(xmk-Xc_WK)<Width_WK.AND. &
            ABS(ymk-Yc_WK)<Ywidth_WK/2.0_SP)THEN
          WK_Source=ZERO
          DO kf=1,Nfreq
           WK_Source=WK_Source+TANH(PI/(Time_ramp/FreqPeak)*TIME)*(Cm(I,J,kf) &
                       *COS(OMGN_IR(KF)*TIME) &
                       +Sm(I,J,kf)*SIN(OMGN_IR(KF)*TIME))
          ENDDO

          R1(I,J)=-1.0_SP/DXg*(F_right*n_right+F_left*n_left) &
                -1.0_SP/DYg*(F_top*n_top+F_bottom*n_bottom) &
        ! wavemaker
                +WK_Source      
         ELSE
         R1(I,J)=-1.0_SP/DXg*(F_right*n_right+F_left*n_left) &
                   -1.0_SP/DYg*(F_top*n_top+F_bottom*n_bottom)
         ENDIF
       ELSEIF(WAVEMAKER(1:6)=='WK_REG')THEN

            xmk=(I-Ibeg)*DXg+npx*(Mloc-2*Nghost)*DXg
            ymk=(J-Jbeg)*DYg+npy*(Nloc-2*Nghost)*DYg




         IF(ABS(xmk-Xc_WK)<Width_WK.AND. &
            ABS(ymk-Yc_WK)<Ywidth_WK/2.0_SP)THEN
          
          R1(I,J)=-1.0_SP/DXg*(F_right*n_right+F_left*n_left) &
                -1.0_SP/DYg*(F_top*n_top+F_bottom*n_bottom) &
        ! wavemaker 
                 +WaveMaker_Mass(I,J)    
         ELSE
         R1(I,J)=-1.0_SP/DXg*(F_right*n_right+F_left*n_left) &
                   -1.0_SP/DYg*(F_top*n_top+F_bottom*n_bottom)
         ENDIF
       ELSEIF(WAVEMAKER(1:7)=='WK_TIME')THEN

            xmk=(I-Ibeg)*DXg+npx*(Mloc-2*Nghost)*DXg
            ymk=(J-Jbeg)*DYg+npy*(Nloc-2*Nghost)*DYg




         IF(ABS(xmk-Xc_WK)<Width_WK.AND. &
            ABS(ymk-Yc_WK)<Ywidth_WK/2.0_SP)THEN
          
          R1(I,J)=-1.0_SP/DXg*(F_right*n_right+F_left*n_left) &
                -1.0_SP/DYg*(F_top*n_top+F_bottom*n_bottom) &
                +WaveMaker_Mass(I,J)      
         ELSE
         R1(I,J)=-1.0_SP/DXg*(F_right*n_right+F_left*n_left) &
                   -1.0_SP/DYg*(F_top*n_top+F_bottom*n_bottom)
         ENDIF   

       ELSEIF(WAVEMAKER(1:9)=='WK_DATA2D')THEN

            xmk=(I-Ibeg)*DXg+npx*(Mloc-2*Nghost)*DXg
            ymk=(J-Jbeg)*DYg+npy*(Nloc-2*Nghost)*DYg




         IF(ABS(xmk-Xc_WK)<Width_WK.AND. &
            ABS(ymk-Yc_WK)<Ywidth_WK/2.0_SP)THEN!
          
          R1(I,J)=-1.0_SP/DXg*(F_right*n_right+F_left*n_left) &
                -1.0_SP/DYg*(F_top*n_top+F_bottom*n_bottom) &
                +WaveMaker_Mass(I,J)      
         ELSE
         R1(I,J)=-1.0_SP/DXg*(F_right*n_right+F_left*n_left) &
                   -1.0_SP/DYg*(F_top*n_top+F_bottom*n_bottom)
         ENDIF 
   
      ELSE ! no wk_wavemaker, theres bug in version 1.1 Dxg,Dyg should be 
           ! replaced by Dxg() and Dy()




        R1(I,J)=-1.0_SP/DX(I,J)*(F_right*n_right+F_left*n_left) &
                   -1.0_SP/DY(I,J)*(F_top*n_top+F_bottom*n_bottom)

      ENDIF









! extra terms for depth-averaged u: 1/R tan theta HV
        R1(I,J)=R1(I,J)+1.0_SP/R_earth*TAN(Lat_theta(I,J))*HV(I,J)







      Eta(I,J)=ALPHA(ISTEP)*Eta0(I,J)+BETA(ISTEP)*(Eta(I,J)+DT*R1(I,J))

! eta_limiter is used for the case as wave touches the seabed within wavemaker
      IF(ETA_LIMITER)THEN
        IF(Eta(I,J)<TroughLimit)Eta(I,J)=TroughLimit
        IF(Eta(I,J)>CrestLimit)Eta(I,J)=CrestLimit
      ENDIF ! end eta_limiter


     ENDDO
     ENDDO

! solve ubar
     DO J=Jbeg,Jend
     DO I=Ibeg,Iend
      F_left=Fx(I,J)
      F_right=Fx(I+1,J)
      F_bottom=Fy(I,J)
      F_top=Fy(I,J+1)





      R2(I,J)=-1.0_SP/DX(I,J)*(F_right*n_right+F_left*n_left) &
                       -1.0_SP/DY(I,J)*(F_top*n_top+F_bottom*n_bottom) &
                        +SourceX(I,J)




      Ubar(I,J)=ALPHA(ISTEP)*Ubar0(I,J)+BETA(ISTEP)*(Ubar(I,J)+DT*R2(I,J))

     ENDDO
     ENDDO

! solve vbar
     DO J=Jbeg,Jend
     DO I=Ibeg,Iend
      F_left=Gx(I,J)
      F_right=Gx(I+1,J)
      F_bottom=Gy(I,J)
      F_top=Gy(I,J+1)





      R3(I,J)=-1.0_SP/DX(I,J)*(F_right*n_right+F_left*n_left) &
                       -1.0_SP/DY(I,J)*(F_top*n_top+F_bottom*n_bottom) &
                       +SourceY(I,J)




      Vbar(I,J)=ALPHA(ISTEP)*Vbar0(I,J)+BETA(ISTEP)*(Vbar(I,J)+DT*R3(I,J))

     ENDDO
     ENDDO

     CALL GET_Eta_U_V_HU_HV(ng)

END SUBROUTINE ESTIMATE_HUV

!-------------------------------------------------------------------------------------
!
!    GET_Eta_U_V_HU_HV is subroutine to obtain Eta, u,v,hu,hv
!
!  HISTORY: 
!       09/17/2010 Fengyan Shi
!       10/14/2012 Fengyan Shi, added nesting bc
!       08/06/2015 Young-Kwang Choi, modified U0 V0 shift 
!       01/27/2016 Fengyan Shi, made serial/parallel codes consistent
!                               added parallel code of periodic bc in y
!
!-------------------------------------------------------------------------------------
SUBROUTINE GET_Eta_U_V_HU_HV(ng)
     USE GLOBAL
     IMPLICIT NONE
     REAL(SP)::Fr,Utotal,Utheta,dep,depl,depr,reta,retal,retar
     REAL(SP),DIMENSION(Mloc,Nloc) :: myA,myC,myD,myF

     REAL(SP),DIMENSION(Nglob) :: AperG,BperG,CperG,DperG,VperG
     REAL(SP),DIMENSION(Mglob,Nglob) ::glbA,glbC,glbD
     REAL(SP),DIMENSION(Mglob+2*Nghost,Nglob+2*Nghost)::VG

      INTEGER :: IM
      INTEGER,INTENT(IN) :: ng

! calculate etar, u and vetar, HU, HV
     H=Eta*Gamma3+Depth


!   tridiagonal coefficient
! x direction

! shift U and V
!     U0=U    !ykchoi (15. 08. 06.) 
!     V0=V    !ykchoi
! ykchoi : U0, V0 are moved to the above part of "Do istage=1,3"

   IF(DISPERSION)THEN

     DO J=Jbeg,Jend
     DO I=Ibeg,Iend
       dep=Max(Depth(I,J),MinDepthFrc)
       depl=Max(Depth(I-1,J),MinDepthFrc)
       depr=Max(Depth(I+1,J),MinDepthFrc)








       tmp1=Gamma1*MASK9(I,J)*(b1/2.0_SP/DX(I,J)/DX(I,J)*dep*dep + b2/DX(I,J)/DX(I,J)*depl*dep)
       tmp2=1.0_SP+Gamma1*MASK9(I,J)*(-b1/DX(I,J)/DX(I,J)*dep*dep-2.0_SP*b2/DX(I,J)/DX(I,J)*dep*dep)
       tmp3=Gamma1*MASK9(I,J)*(b1/2.0_SP/DX(I,J)/DX(I,J)*dep*dep + b2/DX(I,J)/DX(I,J)*dep*depr)
       tmp4=Ubar(I,J)*MASK(I,J)/Max(H(I,J),MinDepthFrc)  &
            + Gamma1*MASK9(I,J)*( -b1/2.0_SP*dep*dep*Vxy(I,J)-b2*dep*DVxy(I,J)) 
! remember to document this part for spherical version 


! I added coupling condition 10/14/2012


        

!--------------------------------------------
! boundary condition for ng>1, AMR
     IF(ng.gt.1)THEN

	IF( n_west .eq. MPI_PROC_NULL ) THEN
        IF(I.eq.Ibeg)THEN
          tmp4=tmp4-tmp1*U(I-1,J)
        ENDIF
	ENDIF
	IF( n_east .eq. MPI_PROC_NULL ) THEN
        IF(I.eq.Iend)THEN
          tmp4=tmp4-tmp3*U(I+1,J)
        ENDIF
	ENDIF







        
     ENDIF
!--------------------------------------------

       IF(tmp2.NE.0.0_SP.OR.MASK(I,J).GT.0)THEN
          myA(I,J)=tmp1/tmp2
          myC(I,J)=tmp3/tmp2
          myD(I,J)=tmp4/tmp2
       ELSE
          myA(I,J)=ZERO
          myC(I,J)=ZERO
          myD(I,J)=ZERO
       ENDIF
     ENDDO
     ENDDO


     call TRIDx(myA,myC,myD,myF)
     U(Ibeg:Iend,Jbeg:Jend) = myF(Ibeg:Iend,Jbeg:Jend)





! y direction

     myA=ZERO
     myC=ZERO
     myD=ZERO

     DO I=Ibeg,Iend
     DO J=Jbeg,Jend
       dep=Max(Depth(I,J),MinDepthFrc)
       depl=Max(Depth(I,J-1),MinDepthFrc)
       depr=Max(Depth(I,J+1),MinDepthFrc)


       tmp1=Gamma1*MASK9(I,J)*(b1/2.0_SP/DY(I,J)/DY(I,J)*dep*dep + b2/DY(I,J)/DY(I,J)*depl*dep) 
       tmp2=1.0_SP+Gamma1*MASK9(I,J)*(-b1/DY(I,J)/DY(I,J)*dep*dep-2.0_SP*b2/DY(I,J)/DY(I,J)*dep*dep) 
       tmp3=Gamma1*MASK9(I,J)*(b1/2.0_SP/DY(I,J)/DY(I,J)*dep*dep + b2/DY(I,J)/DY(I,J)*dep*depr)
       tmp4=Vbar(I,J)*MASK(I,J)/Max(H(I,J),MinDepthFrc)  &
             + Gamma1*MASK9(I,J)*(-b1/2.0_SP*dep*dep*Uxy(I,J)-b2*dep*DUxy(I,J)) 
! remember to document this part for spherical version



        

!--------------------------------------------
! boundary condition for ng>1, AMR
     IF(ng.gt.1)THEN

	IF( n_suth .eq. MPI_PROC_NULL ) THEN
       IF(J.eq.Jbeg)THEN
         tmp4=tmp4-tmp1*V(I,J-1)
       ENDIF
	ENDIF
	IF( n_nrth .eq. MPI_PROC_NULL ) THEN      
	 IF(J.eq.Jend)THEN
         tmp4=tmp4-tmp3*V(I,J+1)
       ENDIF
	ENDIF







        
     ENDIF
!--------------------------------------------
  
       IF(tmp2.NE.0.0_SP.OR.MASK(I,J).GT.0)THEN
         myA(I,J)=tmp1/tmp2
         myC(I,J)=tmp3/tmp2
         myD(I,J)=tmp4/tmp2
       ELSE
         myA(I,J)=ZERO
         myC(I,J)=ZERO
         myD(I,J)=ZERO
       ENDIF
     ENDDO
     ENDDO ! end I

     IF(PERIODIC) THEN

!  Sherman-Morrison algorithm 01/27/2016 fyshi

      CALL TRIDy_periodic (myA,myC,myD,myF)  
       V(Ibeg:Iend,Jbeg:Jend) = myF(Ibeg:Iend,Jbeg:Jend)
     ELSE ! no periodic

       CALL TRIDy(myA,myC,myD,myF)
       V(Ibeg:Iend,Jbeg:Jend) = myF(Ibeg:Iend,Jbeg:Jend)





     ENDIF ! end if periodic

   ELSE  ! if no dispersion
     DO J=Jbeg,Jend
     DO I=Ibeg,Iend  
        U(I,J)=Ubar(I,J)/Max(H(I,J),MinDepthFrc)
        V(I,J)=Vbar(I,J)/Max(H(I,J),MinDepthFrc)
     ENDDO
     ENDDO   

   ENDIF  ! end dispersion

     DO J=Jbeg,Jend
     DO I=Ibeg,Iend   
       IF(MASK(I,J)<1)THEN
        Ubar(I,J)=ZERO
        Vbar(I,J)=ZERO
        U(I,J)=ZERO
        V(I,J)=ZERO
        HU(I,J)=ZERO
        HV(I,J)=ZERO
       ELSE
        HU(I,J)=Max(H(I,J),MinDepthFrc)*U(I,J)
        HV(I,J)=Max(H(I,J),MinDepthFrc)*V(I,J)
! apply Froude cap
        Utotal=SQRT(U(I,J)*U(I,J)+V(I,J)*V(I,J))
        Fr=SQRT(GRAV*Max(H(I,J),MinDepthFrc))
        IF(Utotal/Fr.gt.FroudeCap)THEN
          Utheta=ATAN2(V(I,J),U(I,J))
          U(I,J)=FroudeCap*Fr*COS(Utheta)
          V(I,J)=FroudeCap*Fr*SIN(Utheta)
          HU(I,J)=U(I,J)*Max(H(I,J),MinDepthFrc)
          HV(I,J)=V(I,J)*Max(H(I,J),MinDepthFrc)
        ENDIF
! end Froude cap
       ENDIF
     ENDDO
     ENDDO

!------------ykchoi 07/26/2016
!   --- is this for after correction of Froude cap?
     IF( .NOT. DISPERSION) THEN 			

	  Ubar = HU
	  Vbar = HV

     ENDIF

END SUBROUTINE GET_Eta_U_V_HU_HV

