!------------------------------------------------------------------------------------
!
!      FILE dispersion.F
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
!    CAL_DISPERSION is subroutine to calculation dispersion terms
!    so far V^4 and V^1
!    called by
!       MAIN
!    call DERIVATIVE_XX
!         DERIVATIVE_XY
!    
!    HISTORY: 
!      05/01/2010 Fengyan Shi
!      10/14/2012 Fengyan Shi, added coupling bc,
!                              change derivative_xx_high to second order
!                              according to Harris suggestion
!      08/06/2015 - 08/18/2015 Young-Kwang Choi, modified t-derivatives
!                   corrected U1p,V1p, V2, V3 and omega_1 terms
!
!-------------------------------------------------------------------------------------
SUBROUTINE CAL_DISPERSION(ng)
     USE GLOBAL
     IMPLICIT NONE
     INTEGER, INTENT(IN) :: ng

     REAL(SP),Dimension(Mloc,Nloc) :: DU,DV,DUt,DVt
     REAL(SP) :: UxxVxy,UxyVyy,HUxxHVxy,HUxyHVyy, &
                 UxxVxy_x,UxxVxy_y,UxyVyy_x,UxyVyy_y, &
                 HUxxHVxy_x,HUxxHVxy_y,HUxyHVyy_x,HUxyHVyy_y, &
                 rh,rhx,rhy,reta,ken1,ken2,ken3,ken4,ken5




    

       
! uxx
    CALL DERIVATIVE_XX(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,   &
         MASK9(1:Mloc,1:Nloc),DX(1:Mloc,1:Nloc),U(1:Mloc,1:Nloc),Uxx(1:Mloc,1:Nloc))
! uxy
    CALL DERIVATIVE_XY(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,   &
         MASK9(1:Mloc,1:Nloc),DX(1:Mloc,1:Nloc),DY(1:Mloc,1:Nloc),  &
	   U(1:Mloc,1:Nloc),Uxy(1:Mloc,1:Nloc))
! vxy
    CALL DERIVATIVE_XY(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,   &
         MASK9(1:Mloc,1:Nloc),DX(1:Mloc,1:Nloc),DY(1:Mloc,1:Nloc),  &
	   V(1:Mloc,1:Nloc),Vxy(1:Mloc,1:Nloc))
! vyy
    CALL DERIVATIVE_YY(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,   &
         MASK9(1:Mloc,1:Nloc),DY(1:Mloc,1:Nloc),V(1:Mloc,1:Nloc),Vyy(1:Mloc,1:Nloc))



    IF(SHOW_BREAKING)THEN







       
     CALL DERIVATIVE_X(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,  &
          MASK9(1:Mloc,1:Nloc),DX(1:Mloc,1:Nloc),Eta(1:Mloc,1:Nloc),ETAx(1:Mloc,1:Nloc))
     CALL DERIVATIVE_Y(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,  &
          MASK9(1:Mloc,1:Nloc),DY(1:Mloc,1:Nloc),Eta(1:Mloc,1:Nloc),ETAy(1:Mloc,1:Nloc))

    ENDIF

! DU DV
      DO J=1,Nloc    !ykchoi for nesting 0526
      DO I=1,Mloc
       DU(I,J)=Max(Depth(I,J),MinDepthFrc)*U(I,J)
       DV(I,J)=Max(Depth(I,J),MinDepthFrc)*V(I,J)
      ENDDO
      ENDDO
! ykchoi (15. 08. 18)
! Computation of Etat for ( U1p )_t in conservative form of Shi et al. (2012)
! and for viscosity of wave maker

      DO J=1,Nloc-1   !ykchoi for nesting 0526
      DO I=1,Mloc-1


       
       ETAT(I,J)=-(P(I+1,J)-P(I,J))/DX(I,J)-(Q(I,J+1)-Q(I,J))/DY(I,J)

      ENDDO
      ENDDO 

! ETAT

    IF(SHOW_BREAKING .OR. WAVEMAKER_VIS)THEN


    ENDIF


       
! DUxx
    CALL DERIVATIVE_XX(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,  &
         MASK9(1:Mloc,1:Nloc),DX(1:Mloc,1:Nloc),DU(1:Mloc,1:Nloc),DUxx(1:Mloc,1:Nloc))
! DUxy
    CALL DERIVATIVE_XY(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,  &
         MASK9(1:Mloc,1:Nloc),DX(1:Mloc,1:Nloc),DY(1:Mloc,1:Nloc),  &
	   DU(1:Mloc,1:Nloc),DUxy(1:Mloc,1:Nloc))
! DVxy
    CALL DERIVATIVE_XY(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,  &
         MASK9(1:Mloc,1:Nloc),DX(1:Mloc,1:Nloc),DY(1:Mloc,1:Nloc),  &
	   DV(1:Mloc,1:Nloc),DVxy(1:Mloc,1:Nloc))
! DVyy
    CALL DERIVATIVE_YY(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,  &
         MASK9(1:Mloc,1:Nloc),DY(1:Mloc,1:Nloc),DV(1:Mloc,1:Nloc),DVyy(1:Mloc,1:Nloc))

      

      
! this may affect parallel version
! I added coupling boundary 10/14/2012

! BOUNDARY CONDITION FOR NG=1
    IF(ng.eq.1)THEN

!  boundary conditions

    if(n_west.eq.MPI_PROC_NULL) then




     DO J=1,Nloc
       Uxy(Ibeg,J)=ZERO
       DUxy(Ibeg,J)=ZERO
       Vxy(Ibeg,J)=ZERO
       DVxy(Ibeg,J)=ZERO
       Utxy(Ibeg,J)=ZERO
       DUtxy(Ibeg,J)=ZERO
       Vtxy(Ibeg,J)=ZERO
       DVtxy(Ibeg,J)=ZERO
     ENDDO



  

    endif  



    if(n_east.eq.MPI_PROC_NULL) then




     DO J=1,Nloc
       Uxy(Iend,J)=ZERO
       DUxy(Iend,J)=ZERO
       Vxy(Iend,J)=ZERO
       DVxy(Iend,J)=ZERO
       Utxy(Iend,J)=ZERO
       DUtxy(Iend,J)=ZERO
       Vtxy(Iend,J)=ZERO
       DVtxy(Iend,J)=ZERO
     ENDDO 




    endif  

  

    if(n_suth.eq.MPI_PROC_NULL) then




     DO I=1,Mloc
       Uxy(I,Jbeg)=ZERO
       DUxy(I,Jbeg)=ZERO
       Vxy(I,Jbeg)=ZERO
       DVxy(I,Jbeg)=ZERO
       Utxy(I,Jbeg)=ZERO
       DUtxy(I,Jbeg)=ZERO
       Vtxy(I,Jbeg)=ZERO
       DVtxy(I,Jbeg)=ZERO
     ENDDO   




    endif  



    if(n_nrth.eq.MPI_PROC_NULL) then




     DO I=1,Mloc
       Uxy(I,Jend)=ZERO
       DUxy(I,Jend)=ZERO
       Vxy(I,Jend)=ZERO
       DVxy(I,Jend)=ZERO
       Utxy(I,Jend)=ZERO
       DUtxy(I,Jend)=ZERO
       Vtxy(I,Jend)=ZERO
       DVtxy(I,Jend)=ZERO
     ENDDO 




    endif  

    
    ENDIF !end ng=1
    CALL EXCHANGE_DISPERSION(ng)
     
! calculate V1p  without nonlinear dispersion
     DO J=1,Nloc
     DO I=1,Mloc











! ykchoi( 15. 08. 06.)
! U1p, V1p terms are modified to 0.5_SP*(1.0_SP-Beta_1) --> 0.5_SP*(1.0_SP-Beta_1)*(1.0_SP-Beta_1)
               

       !U1p(I,J)=0.5_SP*(1.0_SP-Beta_1)  & !ykchoi
	 U1p(I,J)=0.5_SP*(1.0_SP-Beta_1)*(1.0_SP-Beta_1)  &
                *DEPTH(I,J)*DEPTH(I,J)*(Uxx(I,J)+Vxy(I,J)) &
               +(Beta_1-1.0_SP)*DEPTH(I,J)*(DUxx(I,J)+DVxy(I,J))

       !V1p(I,J)=0.5_SP*(1.0_SP-Beta_1)  & !ykchoi
       V1p(I,J)=0.5_SP*(1.0_SP-Beta_1)*(1.0_SP-Beta_1)  &
                *DEPTH(I,J)*DEPTH(I,J)  &
                *(Uxy(I,J)+Vyy(I,J)) &
               +(Beta_1-1.0_SP)*DEPTH(I,J)*(DUxy(I,J)+DVyy(I,J))

! end spherical extra dispersion terms which can be implemented in TIME_LEFT


     ENDDO
     ENDDO



END SUBROUTINE CAL_DISPERSION

