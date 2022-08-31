!------------------------------------------------------------------------------------
!
!      FILE misc.F
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
!    INDEX is subroutine to index for MPI
!
!    HISTORY: 
!    05/06/2010 Fengyan Shi
!    01/27/2016 Fengyan Shi, added ProcessorID for periodic bc
!
!-------------------------------------------------------------------------------------

!SUBROUTINE INDEX
!    USE GLOBAL
!    IMPLICIT NONE
!    INTEGER :: icp

! TEMP
! # if defined (1)

!    ALLOCATE(ProcessorID(PX,PY)) 

!    NumberProcessor = px*py
!    dims(1) = px
!    dims(2) = py
!    periods(1) = .false.
!    periods(2) = .false.

! dont use periodic because our TRID solver is not based on
!  the cyclic topology  
!    IF(PERIODIC) periods(2) = .true.  

!    coords(1) = 0
!    coords(2) = 0

!    call MPI_CART_CREATE( MPI_COMM_WORLD, ndims, dims, &
!         periods, reorder, comm2d, ier )
!    call MPI_CART_COORDS( comm2d, myid, 2, coords, ier)

!    npx = coords(1)
!    npy = coords(2)

!    call MPI_Cart_shift( comm2d, 0, 1, n_west, n_east, ier )
!    call MPI_Cart_shift( comm2d, 1, 1, n_suth, n_nrth, ier )

!    icp=0
!    DO I=1,PX
!    DO J=1,PY
!      ProcessorID(I,J) = icp
!      icp=icp+1
!    ENDDO
!    ENDDO

! check
! print*,myid, n_west,n_east,n_suth,n_nrth,ProcessorID(1,1),ProcessorID(1,2)
!      print*,myid,ProcessorID(1,1),ProcessorID(2,1),ProcessorID(1,2),ProcessorID(2,2)
!      call MPI_FINALIZE ( ier )

!# else
!  px=1
!  py=1
!# endif

! now for serial code
!    Mloc=Mglob/px+2*Nghost
!    Nloc=Nglob/py+2*Nghost
!    Mloc1=Mloc+1
!    Nloc1=Nloc+1

!    Ibeg=Nghost+1
!    Iend=Mloc-Nghost
!    Iend1=Mloc1-Nghost
!    Jbeg=Nghost+1
!    Jend=Nloc-Nghost
!    Jend1=Nloc1-Nghost

!END SUBROUTINE INDEX

!-------------------------------------------------------------------------------------
!
!    ESTIMATE_DT is subroutine evaluate dt based in CFL
!
!    HISTORY: 
!       05/06/2010 Fengyan Shi
!
!-------------------------------------------------------------------------------------

SUBROUTINE ESTIMATE_DT(M,N,DX,DY,U,V,H,MinDepthFrc,DT,CFL,TIME)
     USE PARAM

     USE GLOBAL, ONLY : ier

     IMPLICIT NONE
     INTEGER,INTENT(IN)::M,N

     REAL(SP) :: myvar





     REAL(SP),DIMENSION(M,N),INTENT(IN)::DX,DY

     REAL(SP),INTENT(IN),DIMENSION(M,N)::U,V,H
     REAL(SP),INTENT(IN)::CFL,MinDepthFrc
     REAL(SP),INTENT(OUT)::DT
     REAL(SP),INTENT(INOUT)::TIME

     TMP3=LARGE
     DO J=1,N
     DO I=1,M
! x direction
      TMP1=ABS(U(I,J))+SQRT(GRAV*MAX(H(I,J),MinDepthFrc))
      IF(TMP1<SMALL)THEN



       TMP2=DX(I,J)/SMALL

      ELSE



       TMP2=DX(I,J)/TMP1

      ENDIF
      IF(TMP2<TMP3)TMP3=TMP2
! y direction
      TMP1=ABS(V(I,J))+SQRT(GRAV*MAX(H(I,J),MinDepthFrc))
      IF(TMP1<SMALL)THEN



       TMP2=DY(I,J)/SMALL

      ELSE



       TMP2=DY(I,J)/TMP1

      ENDIF
      IF(TMP2<TMP3)TMP3=TMP2      
     ENDDO
     ENDDO

     call MPI_ALLREDUCE (TMP3,myvar,1,MPI_SP,MPI_MIN,&
          MPI_COMM_WORLD,ier)
     TMP3 = myvar

     DT=CFL*TMP3
! TEMP
     TIME=TIME+DT

END SUBROUTINE ESTIMATE_DT



!-------------------------------------------------------------------------------------
!
!    MAX_MIN_PROPERTY is subroutine to calculate Max and Min properties 
!        based on Chen et al., 2004
!
!    HISTORY: 
!        02/10/2016 Fengyan Shi
!
!-------------------------------------------------------------------------------------
SUBROUTINE MAX_MIN_PROPERTY
    USE GLOBAL
    IMPLICIT NONE
    REAL(SP) :: maxv, MaxAbsEta

      DO J=1,Nloc
      DO I=1,Mloc
       IF(OUT_Hmax)THEN
        IF(MASK(I,J).GT.0)THEN
        IF(Eta(I,J).GT.HeightMax(I,J)) HeightMax(I,J)=Eta(I,J)
        ENDIF
       ENDIF

       IF(OUT_Hmin)THEN
        IF(MASK(I,J).GT.0)THEN
        IF(Eta(I,J).LT.HeightMin(I,J)) HeightMin(I,J)=Eta(I,J)
        ENDIF
       ENDIF

       IF(OUT_Umax)THEN
        IF(MASK(I,J).GT.0)THEN
          maxv=SQRT(U(I,J)*U(I,J)+V(I,J)*V(I,J))
          IF(maxV.GT.VelocityMax(I,J)) VelocityMax(I,J)=maxV
        ENDIF
       ENDIF

       IF(OUT_MFmax)THEN
        IF(MASK(I,J).GT.0)THEN
          maxv=(U(I,J)*U(I,J)+V(I,J)*V(I,J))*(H(I,J))
          IF(maxv.GT.MomentumFluxMax(I,J)) MomentumFluxMax(I,J)=maxv
        ENDIF
       ENDIF


      ENDDO
      ENDDO


      IF(OUT_VORmax) THEN
       CALL phi_exch( VorticityMax(1:Mloc,1:Nloc) )
      ENDIF



END SUBROUTINE MAX_MIN_PROPERTY

!-------------------------------------------------------------------------------------
!
!    CHECK_BLOWUP is subroutine to check numerical stability 
!
!    HISTORY: 
!        01/23/2015 Young-Kwang Choi
!        02/15/2016 Fengyan Shi, added the threshold to 100*max_depth
!
!-------------------------------------------------------------------------------------
SUBROUTINE CHECK_BLOWUP
    USE GLOBAL
    IMPLICIT NONE
![ykchoi 15.01.23.
     REAL(SP) :: MaxAbsEta

     REAL(SP)::myvar_tmp


	MaxAbsEta=MAXVAL( abs(Eta(Ibeg:Iend,Jbeg:Jend)) )

      CALL MPI_ALLREDUCE(MaxAbsEta,myvar_tmp,1,MPI_SP,MPI_MAX,MPI_COMM_WORLD,ier)
      MaxAbsEta = myvar_tmp


	if (MaxAbsEta > EtaBlowVal) then

	   if (myid.eq.0) then

	      WRITE(*,*) "========================================="
		  WRITE(*,*) "BlowUp Time, MaxAbsEta=", Time, MaxAbsEta
	      WRITE(*,*) "========================================="

	   endif

	   ICOUNT=99998;
	   CALL PREVIEW_AMR(1)


     call MPI_FINALIZE ( ier )




	endif


END SUBROUTINE CHECK_BLOWUP

!-------------------------------------------------------------------------------------
!
!   wall_time_secs is used to calculate current wall time
!
!   HISTORY: 
!   Gangfeng Ma, 09/12/2011
!
!-------------------------------------------------------------------------------------
    SUBROUTINE wall_time_secs(tcurrent)
    use global, only: SP
    IMPLICIT NONE
    INTEGER, dimension(8) :: walltime
    real(SP), INTENT(OUT) :: tcurrent
    real(SP) :: msecs,secs,mins,hrs,days,months,mscale,years

    call date_and_time(VALUES=walltime)

    msecs = real(walltime(8))
    secs = real(walltime(7))
    mins = real(walltime(6))
    hrs = real(walltime(5))
    days = real(walltime(3))
    months = real(walltime(2))
    years = real(walltime(1))

    if((months.eq.1).or.(months.eq.3).or.(months.eq.5).or.  &
          (months.eq.7).or.(months.eq.8).or.(months.eq.10).or.  &                                                                                   
          (months.eq.12)) then
      mscale = 31.0
    elseif((months.eq.4).or.(months.eq.6).or.  &
          (months.eq.9).or.(months.eq.11)) then
      mscale = 30.0
    elseif(years.eq.4*int(years/4)) then
      mscale = 29.0
    else
      mscale = 28.0
    endif

    tcurrent = months*mscale*24.0*60.0*60.0+days*24.0*60.0*60.0+  &
         hrs*60.0*60.0+60.0*mins+secs+msecs/1000.0

    return
    end SUBROUTINE wall_time_secs




