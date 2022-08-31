!------------------------------------------------------------------------------------
!
!      FILE bc.F
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
!    BOUNDARY_CONDITION is subroutine to provide boundary conditions at edges of domain
!
!    HISTORY: 
!      05/06/2010 Fengyan Shi
!      04/05/2011  Jeff Harris, corrected bugs in do-loop
!
! -----------------------------------------------------------------------------------
SUBROUTINE BOUNDARY_CONDITION(ng)
     USE GLOBAL
     IMPLICIT NONE
     REAL(SP)::Xi,Deps
     INTEGER,INTENT(IN) :: ng

! four sides of computational domain


        if ( n_west .eq. MPI_PROC_NULL ) then

   IF(ng.EQ.1)THEN

       
     DO J=Jbeg,Jend
      P(Ibeg,J)=ZERO
      Xi=EtaRxR(Ibeg,J)
      Deps=Depthx(Ibeg,J)
      Fx(Ibeg,J)=0.5_SP*GRAV*(Xi*Xi*Gamma3+2.0_SP*Xi*Deps)
      Gx(Ibeg,J)=ZERO
      ENDDO

    ENDIF


      endif




        if ( n_east .eq. MPI_PROC_NULL ) then

   IF(ng.EQ.1)THEN

       
     DO J=Jbeg,Jend
      P(Iend1,J)=ZERO
      Xi=EtaRxL(Iend1,J)
      Deps=Depthx(Iend1,J)
      Fx(Iend1,J)=0.5_SP*GRAV*(Xi*Xi*Gamma3+2.0_SP*Xi*Deps)
      Gx(Iend1,J)=ZERO
     ENDDO

   ENDIF


      endif










      if ( n_suth .eq. MPI_PROC_NULL ) then

   IF(ng.EQ.1)THEN

       
     DO I=Ibeg,Iend
      Q(I,Jbeg)=ZERO
      Fy(I,Jbeg)=ZERO
      Xi=EtaRyR(I,Jbeg)
      Deps=Depthy(I,Jbeg)
      Gy(I,Jbeg)=0.5_SP*GRAV*(Xi*Xi*Gamma3+2.0_SP*Xi*Deps)
      ENDDO

   ENDIF


      endif


      if ( n_nrth .eq. MPI_PROC_NULL ) then

   IF(ng.EQ.1)THEN

       
     DO I=Ibeg,Iend
      Q(I,Jend1)=ZERO
      Fy(I,Jend1)=ZERO
      Xi=EtaRyL(I,Jend1)
      Deps=Depthy(I,Jend1)
      Gy(I,Jend1)=0.5_SP*GRAV*(Xi*Xi*Gamma3+2.0_SP*Xi*Deps)
     ENDDO

   ENDIF

     endif






! mask points
! Jeff pointed out the loop should be Jbeg-1, Jend+1
! The problem is that the fluxes on the inter-processor boundaries may be
!modified if the point next to the boundary (e.g., in the ghost cells,
!managed by a different processor) is land, but as is the routine doesnt
!check for this. 

     DO j=Jbeg-1,Jend+1
     DO i=Ibeg-1,Iend+1
      IF(MASK(I,J)<1)THEN
        P(I,J)=ZERO
! Jeff reported a bug here for parallel version

        IF((I/=Ibeg).or.(n_west.ne.MPI_PROC_NULL))THEN



!         Fx(I,J)=0.5_SP*GRAV*HxL(I,J)*HxL(I,J)*MASK(I-1,J)
!new splitting method
      Xi=EtaRxL(I,J)
      Deps=Depthx(I,J)
         Fx(I,J)=0.5_SP*GRAV*(Xi*Xi*Gamma3+2.0_SP*Xi*Deps)*MASK(I-1,J)
        ELSE
         Fx(I,J)=ZERO
        ENDIF
        Gx(I,J)=ZERO

        P(I+1,J)=ZERO
! Jeff also here

        IF((I/=Iend).or.(n_east.ne.MPI_PROC_NULL))THEN



!         Fx(I+1,J)=0.5_SP*GRAV*HxR(I+1,J)*HxR(I+1,J)*MASK(I+1,J)
! new splitting method
      Xi=EtaRxR(I+1,J)
      Deps=Depthx(I+1,J)
         Fx(I+1,J)=0.5_SP*GRAV*(Xi*Xi*Gamma3+2.0_SP*Xi*Deps)*MASK(I+1,J)
        ELSE
         Fx(I+1,J)=ZERO
        ENDIF
        Gx(I+1,J)=ZERO

        Q(I,J)=ZERO
        Fy(I,J)=ZERO
! Jeff also here

        IF((J/=Jbeg).or.(n_suth.ne.MPI_PROC_NULL))THEN



!         Gy(I,J)=0.5_SP*GRAV*HyL(I,J)*HyL(I,J)*MASK(I,J-1)
! new splitting method
      Xi=EtaRyL(I,J)
      Deps=Depthy(I,J)
         Gy(I,J)=0.5_SP*GRAV*(Xi*Xi*Gamma3+2.0_SP*Xi*Deps)*MASK(I,J-1)
        ELSE
         Gy(I,J)=ZERO
        ENDIF

        Q(I,J+1)=ZERO
        Fy(I,J+1)=ZERO
! Jeff also here

        IF((J/=Jend).or.(n_nrth.ne.MPI_PROC_NULL))THEN



!         Gy(I,J+1)=0.5_SP*GRAV*HyR(I,J+1)*HyR(I,J+1)*MASK(I,J+1)
! new splitting method
      Xi=EtaRyR(I,J+1)
      Deps=Depthy(I,J+1)
         Gy(I,J+1)=0.5_SP*GRAV*(Xi*Xi*Gamma3+2.0_SP*Xi*Deps)*MASK(I,J+1)
        ELSE
         Gy(I,J+1)=ZERO
        ENDIF
      ENDIF
     ENDDO
     ENDDO

END SUBROUTINE BOUNDARY_CONDITION


!-------------------------------------------------------------------
!
!   This subroutine is used to collect data into ghost cells                                                         
!
!   HISTORY:
!   07/09/2010 Fengyan Shi, use dummy variables 2) add vtype=3
!
!-------------------------------------------------------------------
SUBROUTINE EXCHANGE_DISPERSION(ng)
    USE GLOBAL
    IMPLICIT NONE
    INTEGER :: VTYPE
    INTEGER,INTENT(IN) :: ng


    VTYPE=2
    CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Uxx(1:Mloc,1:Nloc),VTYPE,ng)
    CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,DUxx(1:Mloc,1:Nloc),VTYPE,ng)
    VTYPE=3
    CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Vyy(1:Mloc,1:Nloc),VTYPE,ng)
    CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,DVyy(1:Mloc,1:Nloc),VTYPE,ng)

    VTYPE=1
    CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Uxy(1:Mloc,1:Nloc),VTYPE,ng)
    CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,DUxy(1:Mloc,1:Nloc),VTYPE,ng)
    CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Vxy(1:Mloc,1:Nloc),VTYPE,ng)
    CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,DVxy(1:Mloc,1:Nloc),VTYPE,ng)

      VTYPE=1  ! symetric in both direction
      CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Ux(1:Mloc,1:Nloc),VTYPE,ng)
      CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,DUx(1:Mloc,1:Nloc),VTYPE,ng)
      CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Vy(1:Mloc,1:Nloc),VTYPE,ng)
      CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,DVy(1:Mloc,1:Nloc),VTYPE,ng)
      VTYPE=3  !like v
      CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Uy(1:Mloc,1:Nloc),VTYPE,ng)
      CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,DUy(1:Mloc,1:Nloc),VTYPE,ng)
      Vtype=2  !like u
      CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Vx(1:Mloc,1:Nloc),VTYPE,ng)
      CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,DVx(1:Mloc,1:Nloc),VTYPE,ng)
      VTYPE=2  ! like u
      CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,ETAx(1:Mloc,1:Nloc),VTYPE,ng)
      VTYPE=3  ! like v
      CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,ETAy(1:Mloc,1:Nloc),VTYPE,ng)   

       
END SUBROUTINE EXCHANGE_DISPERSION

!---------------------------------------------------------------------------------------
!
!   EXCHANGE subroutine is used to collect data into ghost cells                                                         
!
!   HISTORY:
!   07/09/2010 Fengyan Shi 
!     1) use dummy variables 2) add vtype=3
!   08/19/2015 Choi, corrected segmentation fault, maybe memory leaking                                       
!
!---------------------------------------------------------------------------------------
SUBROUTINE EXCHANGE(ng)
    USE GLOBAL
    IMPLICIT NONE
    INTEGER :: VTYPE



    INTEGER,INTENT(IN) :: ng


    VTYPE=1
    CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Eta(1:Mloc,1:Nloc),VTYPE,ng)

    VTYPE=2
    CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,U(1:Mloc,1:Nloc),VTYPE,ng)
    CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,HU(1:Mloc,1:Nloc),VTYPE,ng)
    VTYPE=3
    CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,V(1:Mloc,1:Nloc),VTYPE,ng)
    CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,HV(1:Mloc,1:Nloc),VTYPE,ng)


! etaR x mask is a wrong idea
!    Eta=Eta*MASK

    U=U*MASK
    V=V*MASK
    HU=HU*MASK
    HV=HV*MASK


    
END SUBROUTINE EXCHANGE

!-----------------------------------------------------------------------------------
!
!   PHI_COLL_VARIABLE_LENGTH subroutine is used to collect data into ghost cells                                                         
!
!   HISTORY:
!   07/09/2010 Fengyan Shi                                      
!
!-----------------------------------------------------------------------------------
SUBROUTINE PHI_COLL_VARIABLE_LENGTH(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,PHI,VTYPE,ng)

    USE PARAM


    USE GLOBAL, ONLY : n_east,n_west,n_suth,n_nrth


    IMPLICIT NONE
    INTEGER,INTENT(IN) :: VTYPE
    INTEGER,INTENT(IN) :: Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost
    REAL(SP),INTENT(INOUT) :: PHI(Mloc,Nloc)
    INTEGER,INTENT(IN) :: ng

      ! x-direction

    if ( n_west .eq. MPI_PROC_NULL ) then

	IF(ng==1)THEN
      DO J=Jbeg,Jend  
      DO K=1,Nghost
        PHI(K,J)=PHI(Ibeg+Nghost-K,J)
      ENDDO
      ENDDO
	ENDIF


    endif



    if ( n_east .eq. MPI_PROC_NULL ) then

	IF(ng==1)THEN
      DO J=Jbeg,Jend  
      DO K=1,Nghost
        PHI(Iend+K,J)=PHI(Iend-K+1,J)
      ENDDO
      ENDDO
	ENDIF


    endif



      ! y-direction and corners

    if ( n_suth .eq. MPI_PROC_NULL ) then

	IF(ng==1)THEN
      DO I=1,Mloc
      DO K=1,Nghost
        PHI(I,K)=PHI(I,Jbeg+Nghost-K)
      ENDDO
      ENDDO
	ENDIF

    endif



    if ( n_nrth .eq. MPI_PROC_NULL ) then

	IF(ng==1)THEN
      DO I=1,Mloc
      DO K=1,Nghost
        PHI(I,Jend+K)=PHI(I,Jend-K+1)
      ENDDO
      ENDDO
	ENDIF

    endif



!    call phi_exch_variable_length (PHI,Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost)


END SUBROUTINE PHI_COLL_VARIABLE_LENGTH

!-------------------------------------------------------------------------------------
!
!   PHI_COLL subroutine is used to collect data into ghost cells
!
!   HISTORY:
!     05/01/2010  Fengyan Shi
!     09/07/2010 Fengyan Shi, fix:
!       1) u v symmetric problem, 2) remove use global 3) fix bug
!     05/27/2010 Gangfeng Ma, corrected some bugs
!
!-------------------------------------------------------------------------------------



SUBROUTINE PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,PHI,VTYPE,ng)

    USE PARAM








    USE GLOBAL, ONLY : n_east,n_west,n_suth,n_nrth, &
                     comm2d, ier,myid,PX,PY,&
                       NumberProcessor,ProcessorID


    IMPLICIT NONE
    INTEGER,INTENT(IN) :: VTYPE
    INTEGER,INTENT(IN) :: Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost
    REAL(SP),INTENT(INOUT) :: PHI(Mloc,Nloc)
    INTEGER,INTENT(IN) :: ng




    INTEGER,DIMENSION(1) :: req
    INTEGER :: nreq,len,ll,l,II,JJ
    integer,dimension(MPI_STATUS_SIZE,1) :: status
    REAL(SP),DIMENSION(Mloc,Nghost) :: xx,send2d
    REAL(sp) :: myvar,mybeta

	
 IF (ng==1) THEN
! periodic first because it is not related to VTYPE

! end cartesian
    ENDIF ! ng==1

! I added coupling condition 10/14/2012

! for Eta
    IF(VTYPE==1) THEN  ! for eta
      ! x-direction

    if ( n_west .eq. MPI_PROC_NULL ) then



        
	IF(ng==1)THEN
      DO J=Jbeg,Jend  
      DO K=1,Nghost
        PHI(K,J)=PHI(Ibeg+Nghost-K,J)
      ENDDO
      ENDDO
	ENDIF




    endif



    if ( n_east .eq. MPI_PROC_NULL ) then



        
	IF(ng==1)THEN
      DO J=Jbeg,Jend  
      DO K=1,Nghost
        PHI(Iend+K,J)=PHI(Iend-K+1,J)
      ENDDO
      ENDDO
	ENDIF




    endif



      ! y-direction and corners




    if ( n_suth .eq. MPI_PROC_NULL ) then



        
	IF(ng==1)THEN
      DO I=1,Mloc
      DO K=1,Nghost
        PHI(I,K)=PHI(I,Jbeg+Nghost-K)
      ENDDO
      ENDDO
	ENDIF




    endif



    if ( n_nrth .eq. MPI_PROC_NULL ) then



        
	IF(ng==1)THEN
      DO I=1,Mloc
      DO K=1,Nghost
        PHI(I,Jend+K)=PHI(I,Jend-K+1)
      ENDDO
      ENDDO
	ENDIF




    endif






     ENDIF ! end vtype=1

! for u
    IF(VTYPE==2) THEN  ! for u (x-mirror condition)
      ! x-direction

    if ( n_west .eq. MPI_PROC_NULL ) then



        
	IF(ng==1)THEN
      DO J=Jbeg,Jend
      DO K=1,Nghost
        PHI(K,J)=-PHI(Ibeg+Nghost-K,J)
      ENDDO
      ENDDO
	ENDIF




    endif



    if ( n_east .eq. MPI_PROC_NULL ) then



        
	IF(ng==1)THEN
      DO J=Jbeg,Jend
      DO K=1,Nghost
        PHI(Iend+K,J)=-PHI(Iend-K+1,J)
      ENDDO
      ENDDO
	ENDIF




    endif



      ! y-direction and corners





    if ( n_suth .eq. MPI_PROC_NULL ) then



        
	IF(ng==1)THEN
      DO I=1,Mloc
      DO K=1,Nghost
        PHI(I,K)=PHI(I,Jbeg+Nghost-K)
      ENDDO
      ENDDO
	ENDIF




    endif



    if ( n_nrth .eq. MPI_PROC_NULL ) then



        
	IF(ng==1)THEN
      DO I=1,Mloc
      DO K=1,Nghost
        PHI(I,Jend+K)=PHI(I,Jend-K+1)
      ENDDO
      ENDDO
	ENDIF




    endif






     ENDIF ! end vtype=2

    IF(VTYPE==3) THEN ! for v (y-mirror condition)
! for v
      ! x-direction

    if ( n_west .eq. MPI_PROC_NULL ) then



        
	IF(ng==1)THEN
      DO J=Jbeg,Jend
      DO K=1,Nghost
        PHI(K,J)=PHI(Ibeg+Nghost-K,J)
      ENDDO
      ENDDO
	ENDIF




    endif



    if ( n_east .eq. MPI_PROC_NULL ) then



        
	IF(ng==1)THEN
      DO J=Jbeg,Jend
      DO K=1,Nghost
        PHI(Iend+K,J)=PHI(Iend-K+1,J)
      ENDDO
      ENDDO
	ENDIF




    endif



      ! y-direction and corners



  ! end cartesian

    if ( n_suth .eq. MPI_PROC_NULL ) then



        
	IF(ng==1)THEN
      DO I=1,Mloc
      DO K=1,Nghost
        PHI(I,K)=-PHI(I,Jbeg+Nghost-K)
      ENDDO
      ENDDO
	ENDIF




    endif



    if ( n_nrth .eq. MPI_PROC_NULL ) then



        
	IF(ng==1)THEN
      DO I=1,Mloc
      DO K=1,Nghost
        PHI(I,Jend+K)=-PHI(I,Jend-K+1)
      ENDDO
      ENDDO
	ENDIF




    endif





     ENDIF ! end vtype=3

! for cross-derivatives
    IF(VTYPE==4) THEN ! VTYPE==4 for u and v cross-mirror
     ! x-direction

    if ( n_west .eq. MPI_PROC_NULL ) then



        
	IF(ng==1)THEN
      DO J=Jbeg,Jend
      DO K=1,Nghost
        PHI(K,J)=ZERO
      ENDDO
      ENDDO
	ENDIF




    endif



    if ( n_east .eq. MPI_PROC_NULL ) then



        
	IF(ng==1)THEN
      DO J=Jbeg,Jend
      DO K=1,Nghost
        PHI(Iend+K,J)=ZERO
      ENDDO
      ENDDO
	ENDIF




    endif


      ! y-direction and corners, this one is not an exact solution

    if ( n_suth .eq. MPI_PROC_NULL ) then



        
	IF(ng==1)THEN
      DO I=1,Mloc
      DO K=1,Nghost
        PHI(I,K)=ZERO
      ENDDO
      ENDDO
	ENDIF




    endif



    if ( n_nrth .eq. MPI_PROC_NULL ) then



        
	IF(ng==1)THEN
      DO I=1,Mloc
      DO K=1,Nghost
        PHI(I,Jend+K)=ZERO
      ENDDO
      ENDDO
	ENDIF




    endif



     ENDIF ! end vtype=4

! for symmetric
    IF(VTYPE==5)THEN
      ! x-direction

    if ( n_west .eq. MPI_PROC_NULL ) then



        
	IF(ng==1)THEN
      DO J=Jbeg,Jend
      DO K=1,Nghost
        PHI(K,J)=PHI(Ibeg+Nghost-K,J)
       ENDDO
      ENDDO
	ENDIF




    endif



    if ( n_east .eq. MPI_PROC_NULL ) then



        
	IF(ng==1)THEN
      DO J=Jbeg,Jend
      DO K=1,Nghost
        PHI(Iend+K,J)=PHI(Iend-K+1,J)
      ENDDO
      ENDDO
	ENDIF




    endif


      ! y-direction and corners


    if ( n_suth .eq. MPI_PROC_NULL ) then



        
	IF(ng==1)THEN
      DO I=1,Mloc
      DO K=1,Nghost
        PHI(I,K)=PHI(I,Jbeg+Nghost-K)
      ENDDO
      ENDDO
	ENDIF




    endif



    if ( n_nrth .eq. MPI_PROC_NULL ) then



        
	IF(ng==1)THEN
      DO I=1,Mloc
      DO K=1,Nghost
        PHI(I,Jend+K)=PHI(I,Jend-K+1)
      ENDDO
      ENDDO
	ENDIF




    endif


     ENDIF ! end vtype=5

! for anti-symmetric
     IF(VTYPE==6)THEN
      ! x-direction

    if ( n_west .eq. MPI_PROC_NULL ) then



        
	IF(ng==1)THEN
      DO J=Jbeg,Jend
      DO K=1,Nghost
        PHI(K,J)=-PHI(Ibeg+Nghost-K,J)
      ENDDO
      ENDDO
	ENDIF 




    endif



    if ( n_east .eq. MPI_PROC_NULL ) then



        
	IF(ng==1)THEN
      DO J=Jbeg,Jend
      DO K=1,Nghost
        PHI(Iend+K,J)=-PHI(Iend-K+1,J)
      ENDDO
      ENDDO
	ENDIF 




    endif


      ! y-direction and corners

    if ( n_suth .eq. MPI_PROC_NULL ) then



        
	IF(ng==1)THEN
      DO I=1,Mloc
      DO K=1,Nghost
        PHI(I,K)=-PHI(I,Jbeg+Nghost-K)
      ENDDO
      ENDDO
	ENDIF   




    endif



    if ( n_nrth .eq. MPI_PROC_NULL ) then



        
	IF(ng==1)THEN
      DO I=1,Mloc
      DO K=1,Nghost
        PHI(I,Jend+K)=-PHI(I,Jend-K+1)
      ENDDO
      ENDDO
	ENDIF     




    endif



    ENDIF ! end vtype=6


    call phi_exch (PHI(1:Mloc,1:Nloc))  !AMR 1:Mloc, 1:Nloc


END SUBROUTINE PHI_COLL