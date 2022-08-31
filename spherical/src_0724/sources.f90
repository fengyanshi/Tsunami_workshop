!------------------------------------------------------------------------------------
!
!      FILE sources.F
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
!    SourceTerms is subroutine for all source terms
!
!    HISTORY: 
!       05/01/2010 Fengyan Shi
!       09/26/2013 Babak Tehranirad, added 2D Cd
!       08/18/2015 YoungKwang Choi, modified viscosity breaking
!       02/08/2016 Fengyan Shi, corrected wavemaker corresponding to 
!                               conservative form of momentum equations
!
! --------------------------------------------------
SUBROUTINE SourceTerms(ng)
     USE GLOBAL



     IMPLICIT NONE
     
     INTEGER,INTENT(IN) :: ng
     
     REAL,DIMENSION(Mloc,Nloc) :: nu_vis
     LOGICAL :: PQ_scheme = .TRUE.
     REAL(SP) :: xmk,ymk,Dxg,Dyg,WK_Source
     INTEGER :: kd,kf





! only for wavemaker
     DXg=DX(1,1)
     DYg=DY(1,1)



!ykchoi for nesting code     
     WaveMaker_Mass = ZERO
     IF(ng.Eq.1)THEN

! wavemaker fyshi 02/08/2016

          DO J=1,Nloc
          DO I=1,Mloc


            xmk=(I-Ibeg)*DXg+npx*(Mloc-2*Nghost)*DXg
            ymk=(J-Jbeg)*DYg+npy*(Nloc-2*Nghost)*DYg




         IF(ABS(xmk-Xc_WK)<Width_WK.AND. &
            ABS(ymk-Yc_WK)<Ywidth_WK/2.0_SP)THEN

       IF(WAVEMAKER(1:6)=='WK_REG')THEN          
          WaveMaker_Mass(I,J)=TANH(PI/(Time_ramp*Tperiod)*TIME)*D_gen &
                 *EXP(-Beta_gen*(xmk-Xc_WK)**2)&
                 *SIN(rlamda*(ymk-ZERO)-2.0_SP*PI/Tperiod*TIME)    
         ENDIF
       ENDIF

       IF(WAVEMAKER(1:6)=='WK_IRR')THEN
          WK_Source=ZERO
          DO kf=1,Nfreq
           WK_Source=WK_Source+TANH(PI/(Time_ramp/FreqPeak)*TIME)*(Cm(I,J,kf) &
                       *COS(OMGN_IR(KF)*TIME) &
                       +Sm(I,J,kf)*SIN(OMGN_IR(KF)*TIME))
          ENDDO
          WaveMaker_Mass(I,J)=WK_Source
       ENDIF

       IF(WAVEMAKER(1:7)=='WK_TIME')THEN
           WK_Source=ZERO
           DO kf=1,NumWaveComp
             WK_Source=WK_Source &
               +TANH(PI/(Time_ramp*PeakPeriod)*TIME)*D_genS(kf) &
                 *EXP(-Beta_genS(kf)*(xmk-Xc_WK)**2)&
                 *COS(2.0_SP*PI/WAVE_COMP(kf,1)*TIME-WAVE_COMP(kf,3)) 
           ENDDO
          WaveMaker_Mass(I,J)=WK_Source
       ENDIF

       IF(WAVEMAKER(1:9)=='WK_DATA2D')THEN

        ! make the efficient scheme using Cm and Sm 06/07/2016 fyshi
          WK_Source=ZERO

          DO kf=1,NumFreq
           WK_Source=WK_Source+TANH(PI/(Time_ramp/FreqPeak)*TIME)*(Cm(I,J,kf) &
                       *COS(OMGN2D(KF)*TIME) &
                       +Sm(I,J,kf)*SIN(OMGN2D(KF)*TIME))
          ENDDO
          WaveMaker_Mass(I,J)=WK_Source
       ENDIF

         ENDDO
         ENDDO

     ENDIF

! ----

     nu_vis=ZERO	
	! ykchoi(08.18.2015)
	! new variable :: WAVEMAKER_VIS
      !IF(VISCOSITY_BREAKING)THEN
      IF(VISCOSITY_BREAKING .OR. WAVEMAKER_VIS)THEN
       DO J=Jbeg,Jend
       DO I=Ibeg,Iend
         nu_vis(I,J)=nu_break(I,J)
       ENDDO
       ENDDO
      ENDIF

     IF(ng.Eq.1)THEN

     IF(DIFFUSION_SPONGE)THEN
       DO J=Jbeg,Jend
       DO I=Ibeg,Iend
         nu_vis(I,J)=nu_vis(I,J)+nu_sponge(I,J)
       ENDDO
       ENDDO          
     ENDIF

     ENDIF

107  format(500f12.6)

! depth gradient term
     DO J=Jbeg,Jend
     DO I=Ibeg,Iend



! second order, move the second term to left-hand side


       SourceX(I,J)=GRAV*(Eta(I,J))*SlopeX(I,J)*MASK(I,J) &
                       ! friction



                   -Cd(I,J)*U(I,J)*SQRT(U(I,J)*U(I,J)+V(I,J)*V(I,J)) &

                       ! dispersion
                        ! Ht(+V1p) = div(M)*(-U1p)
                    +Gamma1*MASK9(I,J)*((P(I+1,J)-P(I,J))/DX(I,J)+(Q(I,J+1)-Q(I,J))/DY(I,J)) &
                      *(-U1p(I,J)) &
                        ! Coriolis
                    +Coriolis(I,J)*0.5_SP*(Q(I,J)+Q(I,J+1))
          

       SourceY(I,J)=GRAV*(Eta(I,J))*SlopeY(I,J)*MASK(I,J) &
                          ! friction



                   -Cd(I,J)*V(I,J)*SQRT(U(I,J)*U(I,J)+V(I,J)*V(I,J)) &

                          ! dispersion
                          ! Ht(+V1p) = div(Q)*(-V1p)
                    +Gamma1*MASK9(I,J)*((P(I+1,J)-P(I,J))/DX(I,J)+(Q(I,J+1)-Q(I,J))/DY(I,J)) &
                      *(-V1p(I,J)) &
                        ! Coriolis
                    -Coriolis(I,J)*0.5_SP*(P(I,J)+P(I+1,J))




     ENDDO
     ENDDO







END SUBROUTINE SourceTerms


