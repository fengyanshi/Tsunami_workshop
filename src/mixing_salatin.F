!------------------------------------------------------------------------------------
!
!      FILE mixing.F
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
!-------------------------------------------------------------------------------------
!
!    MIXING_STUFF is subroutine to calculate mixing related, time-averaged properties
!    mean eta is also calculated.
!
!    HISTORY: 05/02/2011 Fengyan Shi
!
!-------------------------------------------------------------------------------------
SUBROUTINE MIXING_STUFF
     USE GLOBAL
     IMPLICIT NONE

! calculate mean for smagorinsky s mixing and wave height
	!ykchoi (be careful)
	!I think Umean, Vmean is not using in other routine.
      IF( time >= STEADY_TIME )THEN
      CALL CALCULATE_MEAN(T_INTV_mean,T_INTV_UVmean,T_sum,T_sumUVmean,&
           DT,Mloc,Nloc,U,V,ETA,ETA0,&
           Umean,Vmean,ETAmean,Usum,Vsum,ETAsum,WaveHeightRMS, &
           WaveHeightAve,Emax,Emin,Num_Zero_Up,Ibeg,Iend,Jbeg,Jend, &
           HrmsSum,HavgSum,ETA2,ETA2sum,SigWaveHeight,ETAvariance, &
           Vor_Instant, Vor_Mean, Vor_Intensity, NumDT1, NumDT2,Vor2,&
           Vor2sum,Vorsum)
	ENDIF    !ykchoi

END SUBROUTINE MIXING_STUFF


!-------------------------------------------------------------------------------------
!
!    CALCULATE_MEAN is subroutine to calculate mean u v required by
!      smagorinsky mixing and wave height
!      mean eta is also calculated.
!
!    HISTORY:
!      05/02/2011 Fengyan Shi
!                 Young-Kwang Choi added some time-averaging stuff
!
!-------------------------------------------------------------------------------------
SUBROUTINE CALCULATE_MEAN(T_INTV_mean,T_INTV_UVmean,T_sum,T_sumUVmean,&
                DT,M,N,U,V,ETA,ETA0,&
                Umean,Vmean,ETAmean,Usum,Vsum,ETAsum,&
                WaveHeightRMS, &
                WaveHeightAve,Emax,Emin,Num_Zero_Up,Ibeg,Iend,Jbeg,Jend, &
                HrmsSum,HavgSum,ETA2,ETA2sum,SigWaveHeight,ETAvariance, &
                Vor_Instant, Vor_Mean, Vor_Intensity, NumDT1, NumDT2,Vor2, &
                Vor2sum,Vorsum)

! calculate mean for smagorinsky mixing and wave height
      USE PARAM
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: M,N,Ibeg,Iend,Jbeg,Jend
      REAL(SP),DIMENSION(M,N),INTENT(IN)::U,V,ETA,ETA0
      REAL(SP),INTENT(IN) :: T_INTV_mean,DT,T_INTV_UVmean
      REAL(SP),DIMENSION(M,N),INTENT(OUT) :: Umean,Vmean
      REAL(SP),DIMENSION(M,N),INTENT(OUT) :: WaveHeightRMS,WaveHeightAve
      REAL(SP),DIMENSION(M,N),INTENT(INOUT) :: Usum,Vsum
      REAL(SP),DIMENSION(M,N),INTENT(INOUT) :: HrmsSum,HavgSum
      REAL(SP),INTENT(OUT) :: T_sum,T_sumUVmean
      REAL(SP)::Tmpe,Tmp_0
      REAL(SP),DIMENSION(M,N),INTENT(INOUT) :: Emax,Emin
      INTEGER,DIMENSION(M,N),INTENT(INOUT) :: Num_Zero_Up
      ! Chen, Salatin
      INTEGER,INTENT(INOUT) :: NumDT1, NumDT2
      REAL(SP),DIMENSION(M,N),INTENT(INOUT) :: Vor_Instant,Vor_Mean,Vor_Intensity
      REAL(SP),DIMENSION(M,N),INTENT(INOUT) :: Vor2,Vor2sum,Vorsum
      REAL(SP),DIMENSION(M,N) :: Vor_Variance
      REAL(SP),DIMENSION(M,N),INTENT(INOUT) :: ETAmean,ETAsum,ETA2,ETA2sum
      REAL(SP),DIMENSION(M,N) :: ETAvariance
      REAL(SP),DIMENSION(M,N),INTENT(INOUT) :: SigWaveHeight

      T_sum=T_sum+DT ! for wave height and MSL
      T_sumUVmean=T_sumUVmean+DT ! for velocity, speed, and vorticity

      IF(T_sum.GE.T_INTV_mean)THEN
        ! Significant Wave Height
        ! Chen, Salatin
        ! variance equation with initially unknown etamean
        ETAsum =ETA*DT+ETAsum
        ETAmean=ETAsum/T_sum
        DO J=1,N
        DO I=1,M
          ETAvariance(I,J) = 1/real(NumDT1-1)*ETA2sum(I,J)&
            -ETAmean(I,J)*ETAmean(I,J)*real(NumDT1)/real(NumDT1-1)
        ENDDO
        ENDDO
        SigWaveHeight = 4.0*SQRT(ETAvariance)
        ETAsum=ZERO
        ETA2sum=ZERO
        NumDT1=ZERO
        T_sum=T_sum-T_INTV_mean

! Average and RMS wave height
        DO J=1,N
        DO I=1,M
        IF(Num_Zero_Up(I,J)>=2)THEN
          WaveHeightAve(I,J)=HavgSum(I,J)/Num_Zero_Up(I,J)
          WaveHeightRMS(I,J)=SQRT(HrmsSum(I,J)/Num_Zero_Up(I,J))
        ENDIF
!        Num_Zero_Up(I,J)=ZERO
!        HavgSum(I,J)=ZERO
!        HrmsSum(I,J)=ZERO
        ENDDO
        ENDDO

        CALL PREVIEW_MEAN

      ELSE

       NumDT1 = NumDT1+1
       ETAsum=ETA*DT+ETAsum
       DO J=1,N
       DO I=1,M
         ETA2(I,J) = ETA(I,J)*ETA(I,J)
       ENDDO
       ENDDO
       ETA2sum=ETA2+ETA2sum

! wave height
       ! Chen, Salatin
       ETAmean = ETAsum/T_sum ! effective for long simulations,
       ! for short simulation just comment it
       DO J=1,N
       DO I=1,M
         if(Eta(i,j)>Emax(i,j)) Emax(i,j) = Eta(i,j)
         if(Eta(i,j)<Emin(i,j)) Emin(i,j) = Eta(i,j)
         Tmpe = Eta(i,j)-ETAmean(i,j)
         Tmp_0 = Eta0(i,j)-ETAmean(i,j)
         if(Tmpe>Tmp_0.and.Tmpe*Tmp_0<=Zero) then
           Num_Zero_Up(i,j) = Num_Zero_Up(i,j)+1
           if(Num_Zero_Up(i,j)>=2) then
               HavgSum(i,j) = HavgSum(i,j)+Emax(i,j)-Emin(i,j)
               HrmsSum(i,j) = HrmsSum(i,j)+(Emax(i,j)-Emin(i,j))**2
           endif
           ! reset Emax and Emin to find next wave
           Emax(i,j) = -1000.
           Emin(i,j) =  1000.
         endif
       ENDDO
       ENDDO

      ENDIF  ! end average time

!-------------------------------------------------------------------------------------
!
!    PREVIEW_UVMEAN is subroutine for print-out of mean U and V field data
!
!  HISTORY:
!    05/13/2020  Chen, Salatin UVmean
!    01/25/2020  Chen, Salatin Vorticity
!
! since U and V are preferably averaged over 2Tp,
! a separate function is assigned for that
!-------------------------------------------------------------------------------------

      IF(T_sumUVmean.GE.T_INTV_UVmean)THEN
        Usum=U*DT+Usum
        Vsum=V*DT+Vsum
        Umean=Usum/T_sumUVmean
        Vmean=Vsum/T_sumUVmean

        !Vorsum =Vor_Instant+Vorsum
        Vorsum =Vor_Instant*DT+Vorsum
        !Vor_Mean=Vorsum/NumDT2
        Vor_Mean=Vorsum/T_sumUVmean
        DO J=1,N
        DO I=1,M
          Vor_Variance(I,J) = 1/real(NumDT2-1)*Vor2sum(I,J)&
            -Vor_Mean(I,J)*Vor_Mean(I,J)*real(NumDT2)/real(NumDT2-1)
          Vor_Intensity(I,J) = SQRT(Vor_Variance(I,J))
        ENDDO
        ENDDO

        T_sumUVmean=T_sumUVmean-T_INTV_UVmean
        Usum=ZERO
        Vsum=ZERO
        Vorsum=ZERO
        Vor2sum=ZERO
        NumDT2=ZERO
        CALL PREVIEW_MEAN_UV
      ELSE
        Usum=U*DT+Usum
        Vsum=V*DT+Vsum
        !Vorticity Intensity
        NumDT2 = NumDT2+1
        !Vorsum =Vor_Instant+Vorsum
        Vorsum =Vor_Instant*DT+Vorsum
        DO J=1,N
        DO I=1,M
          Vor2(I,J) = Vor_Instant(I,J)*Vor_Instant(I,J)
        ENDDO
        ENDDO
        Vor2sum=Vor2+Vor2sum
      ENDIF

END SUBROUTINE CALCULATE_MEAN
