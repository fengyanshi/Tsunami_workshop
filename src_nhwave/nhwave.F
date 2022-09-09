     program NH_WAVE
!--------------------------------------------------------
!
!    Nonhydrostatic WAVE dynamics
!    
!    Code developer: Gangfeng Ma, University of Delaware
!    Last update: 14/04/2011
!    Last update: 14/08/2011, parallel implementation
!
!-------------------------------------------------------
     use global
     implicit none
     integer :: Istage,I,J,K
     real(SP) :: tbegin,tend

# if defined (PARALLEL)
     call MPI_INIT(ier)
     call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ier)
     call MPI_COMM_SIZE(MPI_COMM_WORLD,NumP,ier)
# endif

     ! record wall time
     call wall_time_secs(tbegin)

     ! read parameter input
     call read_input

# if defined (DEBUG)
print*,'read input OK'
# endif

     ! work index
     call index
 
# if defined (DEBUG)
print*,'index OK'
# endif

     ! allocate variables
     call allocate_variables

# if defined (DEBUG)
print*,'allocate variables OK'
# endif

     ! generate grids
     call generate_grid

# if defined (DEBUG)
print*,'generate grid OK'
# endif

     ! read bathymetry
     call read_bathymetry

# if defined (DEBUG)
print*,'read bathymetry OK'
# endif

     ! initialize model run
     call initial

# if defined (DEBUG)
print*,'initial OK'
# endif

     ! time integration
     do while (TIME<TOTAL_TIME.and.RUN_STEP<SIM_STEPS)

       ! time step     
       call estimate_dt

# if defined (DEBUG)
print*,'estimate dt OK'
# endif

# if defined (LANDSLIDE)
       ! run landslide-generated tsunami
       call update_bathymetry
# endif

# if defined (DEBUG)
print*,'update bathymetry OK'
# endif

       ! update boundary conditions       
       call update_wave_bc

# if defined (DEBUG)
print*,'update wave bc OK'
# endif

       ! update mask
       call update_mask

# if defined (DEBUG)
print*,'update mask OK'
# endif

       ! update vars
       call update_vars

# if defined (DEBUG)
print*,'update_vars OK'
# endif

       ! SSP Runge-Kutta time stepping
       do Istage = 1,It_Order

         ! well-balanced source terms
         call source_terms

# if defined (DEBUG)
print*,'source_terms OK'
# endif

         ! fluxes at cell faces
         call fluxes

# if defined (DEBUG)
print*,'fluxes OK'
# endif

         ! update all variables
         call eval_duvw(Istage)

# if defined (DEBUG)
print*,'eval_duvw OK'
# endif

         ! sponge layer
         call sponge_damping

# if defined (DEBUG)
print*,'sponge_damping OK'
# endif

         ! turbulence model
         if(VISCOUS_FLOW) call eval_turb(Istage)

# if defined (DEBUG)
print*,'eval_turb OK'
# endif

# if defined (BUBBLE)
         if(TIME>=TIM_B) then
           ! bubble rise velocity
           call bslip_velocity

           ! update bubble concentration
           call eval_bub(Istage)
         endif
# endif

       enddo


! find max and min elevation
      DO J=1,Nloc
      DO I=1,Mloc
      IF(MASK(I,J).GT.0)THEN
        IF(Eta(I,J).GT.EleMax(I,J)) EleMax(I,J)=Eta(I,J)
        IF(Eta(I,J).LT.EleMin(I,J)) EleMin(I,J)=Eta(I,J)
       DO K=1,Kloc
        IF(Froude(I,J,K).GT.FrMax(I,J,K)) THEN
           FrMax(I,J,K)=Froude(I,J,K)
        ENDIF
        IF(Froude_total(I,J,K).GT.FrMax_total(I,J,K)) THEN
           FrMax_total(I,J,K)=Froude_total(I,J,K)
           TimeMaxFr(I,J) = Time
           DepMaxFr(I,J) = Hc(I,J)
        ENDIF
       ENDDO
      ENDIF
      ENDDO
      ENDDO


       ! screen output
       Screen_Count = Screen_Count+dt
       if(Screen_Count>=Screen_Intv) then
         Screen_Count = Screen_Count-Screen_Intv
         call statistics
       endif

       ! probe output to files
       if(NSTAT>0) then
         Plot_Count_Stat = Plot_Count_Stat+dt
         if(Plot_Count_Stat>=Plot_Intv_Stat) then
           Plot_Count_Stat=Plot_Count_Stat-Plot_Intv_Stat
           call probes
         endif
       endif

       ! field output to files
       if(TIME>=Plot_Start) then
         Plot_Count = Plot_Count+dt
         if(Plot_Count>=Plot_Intv) then
           Plot_Count=Plot_Count-Plot_Intv
           call preview
         endif
       endif
      
     end do

# if defined (PARALLEL)
     if(myid.eq.0) write(*,*) 'Normal Termination!'
     if(myid.eq.0) write(3,*) 'Normal Termination!'
# else
     write(*,*) 'Normal Termination!'
     write(3,*) 'Normal Termination!'
# endif

     ! wall time at the end
     call wall_time_secs(tend)

# if defined (PARALLEL)
     if(myid.eq.0) write(*,*) 'Simulation takes',tend-tbegin,'seconds'
     if(myid.eq.0) write(3,*) 'Simulation takes',tend-tbegin,'seconds'
# else
     write(*,*) 'Simulation takes',tend-tbegin,'seconds'
     write(3,*) 'Simulation takes',tend-tbegin,'seconds'
# endif

# if defined (PARALLEL)
     call MPI_FINALIZE(ier)
# endif

     end

     
     subroutine projection_corrector
!-------------------------------------------
!    Correct the velocity field
!    Called by
!       eval_duvw
!    Last update: 25/03/2011, Gangfeng Ma
!-------------------------------------------
     use global, only: SP,Zero,D,DU,DV,DW,P,dx,dsig, &
                       sigc,DelxH,DelyH,dx,dy, &
                       dt,Rho0,Ibeg,Iend,Jbeg,Jend,Kbeg, &
                       Kend,Mloc,Nloc,Kloc,DelxSc,DelySc,Mask9
     implicit none
     integer :: i,j,k
     real(SP), dimension(:,:,:),allocatable :: DelxP,DelyP,DelzP
 
     allocate(DelxP(Mloc,Nloc,Kloc))
     allocate(DelyP(Mloc,Nloc,Kloc))
     allocate(DelzP(Mloc,Nloc,Kloc))

     DelxP = Zero
     DelyP = Zero
     DelzP = Zero
     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       DelxP(i,j,k) = 0.5*((P(i+1,j,k)-P(i-1,j,k))/(2.*dx)+(P(i+1,j,k+1)-P(i-1,j,k+1))/(2.*dx))
       DelyP(i,j,k) = 0.5*((P(i,j+1,k)-P(i,j-1,k))/(2.*dy)+(P(i,j+1,k+1)-P(i,j-1,k+1))/(2.*dy))
       DelzP(i,j,k) = (P(i,j,k+1)-P(i,j,k))/dsig(k)
     enddo
     enddo
     enddo

     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       if(Mask9(i,j)==0) cycle

       DU(i,j,k) = DU(i,j,k)-D(i,j)*dt/Rho0*(DelxP(i,j,k)+DelzP(i,j,k)*DelxSc(i,j,k))
       DV(i,j,k) = DV(i,j,k)-D(i,j)*dt/Rho0*(DelyP(i,j,k)+DelzP(i,j,k)*DelySc(i,j,k))
       DW(i,j,k) = DW(i,j,k)-dt/Rho0*DelzP(i,j,k)
     enddo
     enddo
     enddo

     deallocate(DelxP)
     deallocate(DelyP)
     deallocate(DelzP)

     return
     end subroutine projection_corrector
     

     subroutine poisson_solver
!--------------------------------------------
!    Solve poisson equation for dynamic pressure
!    Called by
!       eval_duvw
!    Last update: 24/03/2011, Gangfeng Ma
!----------------------------------------------
     use global
     implicit none
     integer :: i,j,k
# if !defined (PARALLEL)
     ! variables for serial computation
     real(SP), dimension(:), allocatable :: Wksp
     integer,  dimension(:), allocatable :: IWksp
     real(SP), dimension(neqns) :: Phi
     real(SP) :: RPARM(30),Pbar(1)
     integer :: IPARM(30),S(1),IS(1),nwksp,inwksp,Ndim,Mdim,N,Maxnz,ierr,neq
     external :: MIC3,IC3,SOR3,GMRES,CG,BCGS
# endif

     ! generate coefficient matrix and rhs
     call generate_coef_rhs

# if defined (PARALLEL)
     ! use HYPRE package for parallel computation
     call hypre_pres_solver
# else
     ! use NSPCG package for serial computation
     call dfault(IPARM,RPARM)

     ! reset default values
     IPARM(2) = itmax
! fyshi set IPARM(3)=0 to turn off write out
     IPARM(3) = 0
     IPARM(4) = 33
     RPARM(1) = tol

     Ndim = 5*neqns
     Mdim = 5*15
     N = neqns
     Maxnz = 15
     nwksp = 30*neqns
     inwksp = 10*neqns

     allocate(Wksp(nwksp))
     allocate(Iwksp(inwksp))
   
     ! initial guess
     neq = 0
     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       neq = neq+1
       Phi(neq) = P(i,j,k)
     enddo
     enddo
     enddo
     
     ! solve poisson equation using preconditioned CG or GMRES
     if(isolver==1) then
       call nspcg(MIC3,CG,Ndim,Mdim,N,Maxnz,Coef,JCoef,S,IS,  &
           Phi,Pbar,Rhs,Wksp,IWksp,nwksp,inwksp,IPARM,RPARM,ierr)
     elseif(isolver==2) then
       call nspcg(IC3,GMRES,Ndim,Mdim,N,Maxnz,Coef,JCoef,S,IS,  &
           Phi,Pbar,Rhs,Wksp,IWksp,nwksp,inwksp,IPARM,RPARM,ierr)  
     elseif(isolver==3) then
       call nspcg(SOR3,GMRES,Ndim,Mdim,N,Maxnz,Coef,JCoef,S,IS,  &
           Phi,Pbar,Rhs,Wksp,IWksp,nwksp,inwksp,IPARM,RPARM,ierr)
     endif   

     neq = 0
     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       neq = neq+1
       P(i,j,k) = Phi(neq)
     enddo
     enddo
     enddo

     deallocate(Wksp)
     deallocate(Iwksp)
# endif

     ! dry cells
     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
! gfma removed Mask9 12/19/2011
       if(Mask(i,j)==0) P(i,j,k) = Zero
     enddo
     enddo
     enddo

     ! collect into ghost cells
# if defined (PARALLEL)
     if(n_west.eq.MPI_PROC_NULL) then
# endif
     do k = Kbeg,Kend
     do j = Jbeg,Jend
       do i = 1,Nghost
         P(Ibeg-i,j,k) = P(Ibeg+i-1,j,k)
       enddo
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_east.eq.MPI_PROC_NULL) then
# endif
     do k = Kbeg,Kend
     do j = Jbeg,Jend 
       do i = 1,Nghost     
         P(Iend+i,j,k) = P(Iend-i+1,j,k)
       enddo
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_suth.eq.MPI_PROC_NULL) then
# endif
     do k = Kbeg,Kend
     do i = Ibeg,Iend
       do j = 1,Nghost
         P(i,Jbeg-j,k) = P(i,Jbeg+j-1,k)
       enddo
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_nrth.eq.MPI_PROC_NULL) then
# endif
     do k = Kbeg,Kend
     do i = Ibeg,Iend
       do j = 1,Nghost
         P(i,Jend+j,k) = P(i,Jend-j+1,k)
       enddo
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     call phi_3D_exch(P)
# endif

     end subroutine poisson_solver


     subroutine hypre_pres_solver
!---------------------------------------------
!    solve pressure using hypre package
!    called by
!       poisson_solver
!    Last update: 22/08/2011, Gangfeng Ma
!---------------------------------------------
# if defined (PARALLEL)
     use global
     implicit none
     integer, parameter :: ndim=3
     integer, parameter :: nentries=15
     integer :: i,j,k,n,ivalues,nvalues,neq,ientry,num_iterations,  &
                precond_id,n_pre,n_post,ierr
     integer*8 :: grid,stencil,matrix,vec_b,vec_x,solver,precond
     integer :: i_glob(Mloc),j_glob(Nloc),k_glob(Kloc)
     integer :: ilower(ndim),iupper(ndim),offsets(nentries,ndim),stencil_indices(nentries), &
                periodic_shift(ndim)
     real(SP) :: final_res_norm
     real(SP), dimension(:), allocatable :: values,Phi
     integer, dimension(:,:,:), allocatable :: indx 
     data ((offsets(i,j),j=1,ndim),i=1,nentries)/0,0,0,1,0,0,0,1,0,0,-1,1,-1,0,1,  &
             0,0,1,1,0,1,0,1,1,-1,0,0,0,-1,0,  &
             0,1,-1,1,0,-1,0,0,-1,-1,0,-1,0,-1,-1/

     ! set up a three dimensional grid
     call HYPRE_StructGridCreate(MPI_COMM_WORLD,ndim,grid,ierr)

     ! global indices
     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       i_glob(i) = npx*(Iend-Ibeg+1)+i-Nghost
       j_glob(j) = npy*(Jend-Jbeg+1)+j-Nghost
       k_glob(k) = k-Nghost
     enddo
     enddo
     enddo

     ilower(1) = i_glob(Ibeg)
     ilower(2) = j_glob(Jbeg)
     ilower(3) = k_glob(Kbeg)
     iupper(1) = i_glob(Iend)
     iupper(2) = j_glob(Jend)
     iupper(3) = k_glob(Kend)

     call HYPRE_StructGridSetExtents(grid,ilower,iupper,ierr)

     if(PERIODIC) then
       periodic_shift(1) = 0
       periodic_shift(2) = Nglob
       periodic_shift(3) = 0
       call HYPRE_StructGridSetPeriodic(grid,periodic_shift,ierr)
     endif

     call HYPRE_StructGridAssemble(grid,ierr)

     ! define the discretization stencil
     call HYPRE_StructStencilCreate(ndim,nentries,stencil,ierr)

     do ientry = 1,nentries
       call HYPRE_StructStencilSetElement(stencil,(ientry-1),offsets(ientry,:),ierr)
     enddo

     ! create matrix object
     call HYPRE_StructMatrixCreate(MPI_COMM_WORLD,grid,stencil,matrix,ierr)

     call HYPRE_StructMatrixInitialize(matrix,ierr)

     ! set the matrix coefficient
     do i = 1,nentries
       stencil_indices(i) = i-1
     enddo

     allocate(indx(Mloc,Nloc,Kloc))
 
     neq = 0
     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       neq = neq+1
       indx(i,j,k) = neq
     enddo
     enddo
     enddo
    
     nvalues = (Iend-Ibeg+1)*(Jend-Jbeg+1)*(Kend-Kbeg+1)*nentries
     allocate(values(nvalues))

     ivalues = 0
     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       do n = 1,nentries
         ivalues = ivalues+1
         values(ivalues) = Coef(indx(i,j,k),n)
       enddo
     enddo
     enddo
     enddo

     call HYPRE_StructMatrixSetBoxValues(matrix,ilower,iupper,nentries,  &
                                  stencil_indices,values,ierr) 
     call HYPRE_StructMatrixAssemble(matrix,ierr)
     !call HYPRE_StructMatrixPrint(matrix,zero,ierr)

     ! set up struct vectors for b and x
     call HYPRE_StructVectorCreate(MPI_COMM_WORLD,grid,vec_b,ierr)
     call HYPRE_StructVectorCreate(MPI_COMM_WORLD,grid,vec_x,ierr)

     call HYPRE_StructVectorInitialize(vec_b,ierr)
     call HYPRE_StructVectorInitialize(vec_x,ierr)

     ! set the vector coefficients
     call HYPRE_StructVectorSetBoxValues(vec_b,ilower,iupper,Rhs,ierr)   
     call HYPRE_StructVectorAssemble(vec_b,ierr)     
     !call HYPRE_StructVectorPrint(vec_b,zero,ierr)

     ! initial guess
     allocate(Phi(neqns))
     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       Phi(indx(i,j,k)) = P(i,j,k)
     enddo
     enddo
     enddo
     
     call HYPRE_StructVectorSetBoxValues(vec_x,ilower,iupper,Phi,ierr)
     call HYPRE_StructVectorAssemble(vec_x,ierr)
     !call HYPRE_StructVectorPrint(vec_x,zero,ierr)

     ! set up and use a solver
     call HYPRE_StructGMRESCreate(MPI_COMM_WORLD,solver,ierr)
     call HYPRE_StructGMRESSetMaxIter(solver,itmax,ierr)
     call HYPRE_StructGMRESSetTol(solver,tol,ierr)
     call HYPRE_StructGMRESSetPrintLevel(solver,0,ierr)
     call HYPRE_StructGMRESSetLogging(solver,0,ierr)

     ! use symmetric SMG as preconditioner
     n_pre = 1; n_post = 1
     call HYPRE_StructSMGCreate(MPI_COMM_WORLD,precond,ierr)
     call HYPRE_StructSMGSetMemoryUse(precond,0,ierr)
     call HYPRE_StructSMGSetMaxIter(precond,1,ierr)
     call HYPRE_StructSMGSetTol(precond,0.0,ierr)
     call HYPRE_StructSMGSetNumPreRelax(precond,n_pre,ierr)
     call HYPRE_StructSMGSetNumPostRelax(precond,n_post,ierr)
     call HYPRE_StructSMGSetLogging(precond,0,ierr)

     ! set up preconditioner
     precond_id = 0
     call HYPRE_StructGMRESSetPrecond(solver,precond_id,precond,ierr)
     
     ! do the setup
     call HYPRE_StructGMRESSetup(solver,matrix,vec_b,vec_x,ierr)
 
     ! do the solve
     call HYPRE_StructGMRESSolve(solver,matrix,vec_b,vec_x,ierr)

     ! get results
     call HYPRE_StructVectorGetBoxValues(vec_x,ilower,iupper,Phi,ierr)

     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       P(i,j,k) = Phi(indx(i,j,k))
     enddo
     enddo
     enddo

     ! get some info
     !call HYPRE_StructGMRESGetFinalRelati(solver,final_res_norm,ierr)
     !call HYPRE_StructGMRESGetNumIteratio(solver,num_iterations,ierr);
     !
     !if(myid.eq.0) then
     !  write(*,*)'Iterations = ',num_iterations
     !  write(*,*)'Final Relative Residual Norm = ',final_res_norm
     !endif

     ! free memory
     call HYPRE_StructGridDestroy(grid,ierr)
     call HYPRE_StructStencilDestroy(stencil,ierr)
     call HYPRE_StructMatrixDestroy(matrix,ierr)
     call HYPRE_StructVectorDestroy(vec_b,ierr)
     call HYPRE_StructVectorDestroy(vec_x,ierr)
     call HYPRE_StructGMRESDestroy(solver,ierr)
     call HYPRE_StructSMGDestroy(precond,ierr)

     deallocate(indx)
     deallocate(values)
     deallocate(Phi)

# endif
     return
     end subroutine hypre_pres_solver


     subroutine generate_coef_rhs
!---------------------------------------------
!    Generate coefficient matrix and rhs
!    Called by 
!       poisson_solver
!    Last update: 24/03/2011, Gangfeng Ma
!--------------------------------------------
     use global
     implicit none
     integer :: i,j,k,neq,n,ic
     real(SP), dimension(:,:,:), allocatable :: DelxS,DelyS,DelzS,A1
     integer,  dimension(:,:,:), allocatable :: indx

     allocate(DelxS(Mloc,Nloc,Kloc1))
     allocate(DelyS(Mloc,Nloc,Kloc1))
     allocate(DelzS(Mloc,Nloc,Kloc1))
     allocate(A1(Mloc,Nloc,Kloc1))
     allocate(indx(Mloc,Nloc,Kloc))

     DelxS = Zero
     DelyS = Zero
     DelzS = Zero
     A1 = Zero
     do k = Kbeg,Kend1
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       DelxS(i,j,k) = (1.-sig(k))/D(i,j)*DelxH(i,j)-sig(k)/D(i,j)*DelxEta(i,j)
       DelyS(i,j,k) = (1.-sig(k))/D(i,j)*DelyH(i,j)-sig(k)/D(i,j)*DelyEta(i,j)
       DelzS(i,j,k) = 1./D(i,j)

       A1(i,j,k) = DelxS(i,j,k)*DelxS(i,j,k)+DelyS(i,j,k)*DelyS(i,j,k)+  &
            DelzS(i,j,k)*DelzS(i,j,k)
     enddo
     enddo
     enddo
   
     ! generate coefficient matrix
     neq = 0
     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       neq = neq+1
       indx(i,j,k) = neq
     enddo
     enddo 
     enddo

     ! generate source term 
     Rhs = Zero
     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       Rhs(indx(i,j,k)) = -((Uf(i+1,j,k)-Uf(i-1,j,k))/(2.0*dx)+(U(i,j,k)-U(i,j,k-1))/(0.5*(dsig(k)+dsig(k-1)))*  &
              DelxS(i,j,k)+(Vf(i,j+1,k)-Vf(i,j-1,k))/(2.0*dy)+(V(i,j,k)-V(i,j,k-1))/(0.5*(dsig(k)+dsig(k-1)))*  &
              DelyS(i,j,k)+(W(i,j,k)-W(i,j,k-1))/(0.5*(dsig(k)+dsig(k-1)))*DelzS(i,j,k)-SourceC(i,j))*Rho0/dt
     enddo
     enddo
     enddo

     Coef = Zero
     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       Coef(indx(i,j,k),1) = (2./(dx*dx)+2./(dy*dy)+A1(i,j,k)/(0.5*(dsig(k)+dsig(k-1))*dsig(k))+  &
                A1(i,j,k)/(0.5*(dsig(k)+dsig(k-1))*dsig(k-1)))
       Coef(indx(i,j,k),2) = -1./(dx*dx)
       Coef(indx(i,j,k),3) = -1./(dy*dy)
       Coef(indx(i,j,k),4) = (DelyS(i,j-1,k)/(2.*dy*(dsig(k)+dsig(k-1)))+DelyS(i,j,k)/(2.*dy*(dsig(k)+dsig(k-1))))   
       Coef(indx(i,j,k),5) = (DelxS(i-1,j,k)/(2.*dx*(dsig(k)+dsig(k-1)))+DelxS(i,j,k)/(2.*dx*(dsig(k)+dsig(k-1))))
       Coef(indx(i,j,k),6) = -A1(i,j,k)/(0.5*(dsig(k)+dsig(k-1))*dsig(k))
       Coef(indx(i,j,k),7) = -(DelxS(i+1,j,k)/(2.*dx*(dsig(k)+dsig(k-1)))+DelxS(i,j,k)/(2.*dx*(dsig(k)+dsig(k-1))))
       Coef(indx(i,j,k),8) = -(DelyS(i,j+1,k)/(2.*dy*(dsig(k)+dsig(k-1)))+DelyS(i,j,k)/(2.*dy*(dsig(k)+dsig(k-1))))
       Coef(indx(i,j,k),9) = -1./(dx*dx)
       Coef(indx(i,j,k),10) = -1./(dy*dy)
       Coef(indx(i,j,k),11) = (DelyS(i,j+1,k)/(2.*dy*(dsig(k)+dsig(k-1)))+DelyS(i,j,k)/(2.*dy*(dsig(k)+dsig(k-1))))
       Coef(indx(i,j,k),12) = (DelxS(i+1,j,k)/(2.*dx*(dsig(k)+dsig(k-1)))+DelxS(i,j,k)/(2.*dx*(dsig(k)+dsig(k-1))))
       Coef(indx(i,j,k),13) = -A1(i,j,k)/(0.5*(dsig(k)+dsig(k-1))*dsig(k-1))
       Coef(indx(i,j,k),14) = -(DelxS(i-1,j,k)/(2.*dx*(dsig(k)+dsig(k-1)))+DelxS(i,j,k)/(2.*dx*(dsig(k)+dsig(k-1))))
       Coef(indx(i,j,k),15) = -(DelyS(i,j-1,k)/(2.*dy*(dsig(k)+dsig(k-1)))+DelyS(i,j,k)/(2.*dy*(dsig(k)+dsig(k-1))))
     enddo
     enddo
     enddo

     ! boundary conditions
     ! left side
# if defined (PARALLEL)
     if(n_west.eq.MPI_PROC_NULL) then
# endif
     i = Ibeg
     do k = Kbeg,Kend
     do j = Jbeg,Jend
       ic = indx(i,j,k)
       Coef(ic,1) = Coef(ic,1)+Coef(ic,9)
       Coef(ic,6) = Coef(ic,6)+Coef(ic,5)
       Coef(ic,13) = Coef(ic,13)+Coef(ic,14)
       Coef(ic,9) = Zero
       Coef(ic,5) = Zero
       Coef(ic,14) = Zero
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif

     ! right side
# if defined (PARALLEL)
     if(n_east.eq.MPI_PROC_NULL) then
# endif
     i = Iend
     do k = Kbeg,Kend
     do j = Jbeg,Jend
       ic = indx(i,j,k)
       Coef(ic,1) = Coef(ic,1)+Coef(ic,2)
       Coef(ic,6) = Coef(ic,6)+Coef(ic,7)
       Coef(ic,13) = Coef(ic,13)+Coef(ic,12)
       Coef(ic,2) = Zero
       Coef(ic,7) = Zero
       Coef(ic,12) = Zero
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif

     ! front side
# if defined (PARALLEL)
     if(n_suth.eq.MPI_PROC_NULL) then
# endif
     j = Jbeg
     do k = Kbeg,Kend
     do i = Ibeg,Iend
       ic = indx(i,j,k)         
       Coef(ic,1) = Coef(ic,1)+Coef(ic,10)
       Coef(ic,6) = Coef(ic,6)+Coef(ic,4)
       Coef(ic,13) = Coef(ic,13)+Coef(ic,15)
       Coef(ic,10) = Zero
       Coef(ic,4) = Zero
       Coef(ic,15) = Zero
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif

     ! back side
# if defined (PARALLEL)
     if(n_nrth.eq.MPI_PROC_NULL) then
# endif
     j = Jend
     do k = Kbeg,Kend
     do i = Ibeg,Iend
       ic = indx(i,j,k)
       Coef(ic,1) = Coef(ic,1)+Coef(ic,3)
       Coef(ic,6) = Coef(ic,6)+Coef(ic,8)
       Coef(ic,13) = Coef(ic,13)+Coef(ic,11)
       Coef(ic,3) = Zero
       Coef(ic,8) = Zero
       Coef(ic,11) = Zero
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif

     ! bottom side
     k = Kbeg
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       ic = indx(i,j,k)
# if defined (LANDSLIDE)
       Rhs(ic) = Rhs(ic)+Rho0*(dsig(Kbeg)+dsig(Kbeg-1))*(Coef(ic,13)*D(i,j)*Delt2H(i,j)+ &
            Coef(ic,12)*D(i+1,j)*Delt2H(i+1,j)+Coef(ic,11)*D(i,j+1)*Delt2H(i,j+1)+ &
            Coef(ic,14)*D(i-1,j)*Delt2H(i-1,j)+Coef(ic,15)*D(i,j-1)*Delt2H(i,j-1))
# endif
       Coef(ic,6) = Coef(ic,6)+Coef(ic,13)
       Coef(ic,7) = Coef(ic,7)+Coef(ic,12)
       Coef(ic,8) = Coef(ic,8)+Coef(ic,11)
       Coef(ic,5) = Coef(ic,5)+Coef(ic,14)
       Coef(ic,4) = Coef(ic,4)+Coef(ic,15)
       Coef(ic,13) = Zero
       Coef(ic,12) = Zero
       Coef(ic,11) = Zero
       Coef(ic,14) = Zero
       Coef(ic,15) = Zero
     enddo
     enddo

     ! top side (Dirichlet boundary)
     k = Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       ic = indx(i,j,k)
       Coef(ic,4) = Zero
       Coef(ic,5) = Zero
       Coef(ic,6) = Zero
       Coef(ic,7) = Zero
       Coef(ic,8) = Zero
     enddo
     enddo

     ! take (i=2,j=2,k=2) to obtain the diagonal information
     JCoef(1) = indx(Ibeg+1,Jbeg+1,Kbeg+1)-indx(Ibeg+1,Jbeg+1,Kbeg+1)  ! (i,j,k)
     JCoef(2) = indx(Ibeg+2,Jbeg+1,Kbeg+1)-indx(Ibeg+1,Jbeg+1,Kbeg+1)  ! (i+1,j,k) 
     JCoef(3) = indx(Ibeg+1,Jbeg+2,Kbeg+1)-indx(Ibeg+1,Jbeg+1,Kbeg+1)  ! (i,j+1,k)
     JCoef(4) = indx(Ibeg+1,Jbeg,Kbeg+2)-indx(Ibeg+1,Jbeg+1,Kbeg+1)    ! (i,j-1,k+1)
     JCoef(5) = indx(Ibeg,Jbeg+1,Kbeg+2)-indx(Ibeg+1,Jbeg+1,Kbeg+1)    ! (i-1,j,k+1)
     JCoef(6) = indx(Ibeg+1,Jbeg+1,Kbeg+2)-indx(Ibeg+1,Jbeg+1,Kbeg+1)  ! (i,j,k+1)
     JCoef(7) = indx(Ibeg+2,Jbeg+1,Kbeg+2)-indx(Ibeg+1,Jbeg+1,Kbeg+1)  ! (i+1,j,k+1)
     JCoef(8) = indx(Ibeg+1,Jbeg+2,Kbeg+2)-indx(Ibeg+1,Jbeg+1,Kbeg+1)  ! (i,j+1,k+1)
     JCoef(9) = indx(Ibeg,Jbeg+1,Kbeg+1)-indx(Ibeg+1,Jbeg+1,Kbeg+1)    ! (i-1,j,k)
     JCoef(10) = indx(Ibeg+1,Jbeg,Kbeg+1)-indx(Ibeg+1,Jbeg+1,Kbeg+1)   ! (i,j-1,k)
     JCoef(11) = indx(Ibeg+1,Jbeg+2,Kbeg)-indx(Ibeg+1,Jbeg+1,Kbeg+1)   ! (i,j+1,k-1)
     JCoef(12) = indx(Ibeg+2,Jbeg+1,Kbeg)-indx(Ibeg+1,Jbeg+1,Kbeg+1)   ! (i+1,j,k-1)
     JCoef(13) = indx(Ibeg+1,Jbeg+1,Kbeg)-indx(Ibeg+1,Jbeg+1,Kbeg+1)   ! (i,j,k-1)
     JCoef(14) = indx(Ibeg,Jbeg+1,Kbeg)-indx(Ibeg+1,Jbeg+1,Kbeg+1)     ! (i-1,j,k-1)
     JCoef(15) = indx(Ibeg+1,Jbeg,Kbeg)-indx(Ibeg+1,Jbeg+1,Kbeg+1)     ! (i,j-1,k-1)

     deallocate(DelxS)
     deallocate(DelyS)
     deallocate(DelzS)
     deallocate(A1) 
     deallocate(indx)

     return
     end subroutine generate_coef_rhs


     subroutine eval_duvw(ISTEP)
!-----------------------------------------------
!    Update all variables D,U,V,W,Omega
!    Called by
!       main
!    Last update: 25/12/2010, Gangfeng Ma
!----------------------------------------------
     use global
     implicit none
     integer,intent(in) :: ISTEP
     real(SP),dimension(:,:),allocatable :: R1
     real(SP),dimension(:,:,:),allocatable :: R2,R3,R4
     integer :: i,j,k,n,Ista

     !FS by DJN June 26, 2013
     REAL(SP) :: Fr,Utotal,Utheta,U_water_slide
!     REAL(SP) :: FroudeCap=0.5

     allocate(R1(Mloc,Nloc))
     allocate(R2(Mloc,Nloc,Kloc))
     allocate(R3(Mloc,Nloc,Kloc))
     allocate(R4(Mloc,Nloc,Kloc))

     ! tidal wave
# if defined (PARALLEL)
     if(n_west.eq.MPI_PROC_NULL) then
# endif
     if(WaveMaker(1:7)=='LEF_TID') then
       Ista = Ibeg+1
     else
       Ista = Ibeg
     endif
# if defined (PARALLEL)
     endif
# endif

     ! estimate diffusion terms
     if(VISCOUS_FLOW) call diffusion

     ! solve total water depth D
     R1 = Zero
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       do k = Kbeg,Kend
         R1(i,j) = R1(i,j)-1.0/dx*(Ex(i+1,j,k)-Ex(i,j,k))*dsig(k)  &
                        -1.0/dy*(Ey(i,j+1,k)-Ey(i,j,k))*dsig(k)
       enddo
       ! internal wavemaker
       R1(i,j) = R1(i,j)+D(i,j)*SourceC(i,j)
       D(i,j) = ALPHA(ISTEP)*D0(i,j)+BETA(ISTEP)*(D(i,j)+dt*R1(i,j))
     enddo
     enddo

     ! solve DU
     R2 = Zero
     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       R2(i,j,k) = -1.0/dx*(Fx(i+1,j,k)-Fx(i,j,k))-1.0/dy*(Fy(i,j+1,k)-Fy(i,j,k)) &
                         -1.0/dsig(k)*(Fz(i,j,k+1)-Fz(i,j,k))+SourceX(i,j) &
                         +Diffxx(i,j,k)+Diffxy(i,j,k)+Diffxz(i,j,k)
       DU(i,j,k) = DU(i,j,k)+dt*R2(i,j,k)
       if(Mask(i,j)==0) DU(i,j,k) = Zero
     enddo
     enddo
     enddo

     ! solve DV
     R3 = Zero
     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       R3(i,j,k) = -1.0/dx*(Gx(i+1,j,k)-Gx(i,j,k))-1.0/dy*(Gy(i,j+1,k)-Gy(i,j,k)) &
                        -1.0/dsig(k)*(Gz(i,j,k+1)-Gz(i,j,k))+SourceY(i,j)  &
                        +Diffyx(i,j,k)+Diffyy(i,j,k)+Diffyz(i,j,k) 
       DV(i,j,k) = DV(i,j,k)+dt*R3(i,j,k)
       if(Mask(i,j)==0) DV(i,j,k) = Zero
     enddo
     enddo
     enddo

     ! solve DW
     R4 = Zero
     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       R4(i,j,k) = -1.0/dx*(Hx(i+1,j,k)-Hx(i,j,k))-1.0/dy*(Hy(i,j+1,k)-Hy(i,j,k)) &
                        -1.0/dsig(k)*(Hz(i,j,k+1)-Hz(i,j,k))  &
                        +Diffzx(i,j,k)+Diffzy(i,j,k)+Diffzz(i,j,k) 
       DW(i,j,k) = DW(i,j,k)+dt*R4(i,j,k)
       if(Mask(i,j)==0) DW(i,j,k) = Zero
     enddo
     enddo
     enddo

     ! update D and Eta
     D = max(D,MinDep)
     call phi_2D_coll(D)
     Eta = D-Hc

     call delxFun_2D(Eta,DelxEta)
     call delyFun_2D(Eta,DelyEta)

     ! sigma transformation coefficient
     call sigma_transform

     ! run non-hydrostatic simulation  
     if(NON_HYDRO) then
       ! obtain hydrostatic velocity
       call get_UVW

       ! interpolate velocity into vertical faces
       call interpolate_velocity_to_faces

       ! solve dynamic pressure 
       call poisson_solver

       ! correct velocity field  
       call projection_corrector
     endif

     ! SSP Runge-Kutta time stepping
     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       DU(i,j,k) = ALPHA(ISTEP)*DU0(i,j,k)+BETA(ISTEP)*DU(i,j,k)
       DV(i,j,k) = ALPHA(ISTEP)*DV0(i,j,k)+BETA(ISTEP)*DV(i,j,k)
       DW(i,j,k) = ALPHA(ISTEP)*DW0(i,j,k)+BETA(ISTEP)*DW(i,j,k)

       if(Mask(i,j)==0) then
         DU(i,j,k) = Zero
         DV(i,j,k) = Zero
         DW(i,j,k) = Zero
       endif
     enddo
     enddo
     enddo  


     !FS by DJN June 26, 2013
     ! Froude number cap
     do j = Jbeg,Jend
     do i = Ibeg,Iend
        Fr=SQRT(GRAV*MAX(D(I,J),MinDep))*MAX(D(I,J),MinDep)
        
      IF(MASK(I,J)>ZERO)THEN
       do k = Kbeg,Kend
         Utotal=SQRT(DU(I,J,K)*DU(I,J,K)+DV(I,J,K)*DV(I,J,K) &
          +DW(I,J,K)*DW(I,J,K))
         U_water_slide= SQRT((DU(I,J,K)-SldUs(I,J))*(DU(I,J,K)-SldUs(I,J)) &
                           +(DV(I,J,K)-SldVs(I,J))*(DV(I,J,K)-SldVs(I,J)) &
                            )
         Froude(I,J,K) = Utotal/Fr
         Froude_total(I,J,K) = U_water_slide/Fr
         IF(Utotal/Fr.gt.FroudeCap)THEN
          Utheta=ATAN2(DV(I,J,K),DU(I,J,K))
          DU(I,J,K)=FroudeCap*Fr*COS(Utheta)
          DV(I,J,K)=FroudeCap*Fr*SIN(Utheta)
          DW(I,J,K)=ZERO
         ENDIF ! fr_cap
       enddo
      ENDIF ! mask
     enddo
     enddo



     ! obtain final velocity field
     call get_UVW

     ! update Omega
     do i = Ibeg,Iend
     do j = Jbeg,Jend
     do k = Kbeg,Kend
       Omega(i,j,k) = DeltH(i,j)-sigc(k)*R1(i,j)+U(i,j,k)*((1.0-sigc(k))*DelxH(i,j)-sigc(k)*DelxEta(i,j))+  &
              V(i,j,k)*((1.0-sigc(k))*DelyH(i,j)-sigc(k)*DelyEta(i,j))+W(i,j,k)
       if(Mask(i,j)==0) Omega(i,j,k) = Zero
     enddo
     do n = 1,Nghost
       Omega(i,j,Kbeg-n) = -Omega(i,j,Kbeg+n-1)
       Omega(i,j,Kend+n) = -Omega(i,j,Kend-n+1)
     enddo
     enddo
     enddo

     deallocate(R1)
     deallocate(R2)
     deallocate(R3)
     deallocate(R4)

     end subroutine eval_duvw


     subroutine interpolate_velocity_to_faces
!------------------------------------------------                                        
!    Interpolate U,V,W to vertical faces                                                 
!    Called by                                                                           
!       main                                                                             
!    Last Update: 19/03/2011, Gangfeng Ma                                                
!------------------------------------------------                                        
     use global, only: SP,U,V,W,Uf,Vf,Wf,dsig,sig,sigc,  &
                       Mloc,Nloc,Kloc,Kloc1,Nghost
     implicit none
     integer  :: i,j,k
     real(SP) :: Int_factor,Int_factor1,Int_factor2,Int_factor3
     logical  :: Linear_Interp=.True.

     if(Linear_Interp) then
       ! first-order linear interpolation
       do k = 2,Kloc
       do j = 1,Nloc
       do i = 1,Mloc
         Int_factor = dsig(k)/(dsig(k)+dsig(k-1))
         Uf(i,j,k) = (1.0-Int_factor)*U(i,j,k)+Int_factor*U(i,j,k-1)
         Vf(i,j,k) = (1.0-Int_factor)*V(i,j,k)+Int_factor*V(i,j,k-1)
         Wf(i,j,k) = (1.0-Int_factor)*W(i,j,k)+Int_factor*W(i,j,k-1)
       enddo
       enddo
       enddo

       do j = 1,Nloc
       do i = 1,Mloc
         Uf(i,j,1) = U(i,j,1)
         Vf(i,j,1) = V(i,j,1)
         Wf(i,j,1) = W(i,j,1)
         Uf(i,j,Kloc1) = U(i,j,Kloc)
         Vf(i,j,Kloc1) = V(i,j,Kloc)
         Wf(i,j,Kloc1) = W(i,j,Kloc)
       enddo
       enddo
     else
       ! second-order lagrange interpolation
       do k = 3,Kloc
       do j = 1,Nloc
       do i = 1,Mloc
         Int_factor1 = (sig(k)-sigc(k-1))*(sig(k)-sigc(k))/  &
             ((sigc(k-2)-sigc(k-1))*(sigc(k-2)-sigc(k)))
         Int_factor2 = (sig(k)-sigc(k-2))*(sig(k)-sigc(k))/  &
             ((sigc(k-1)-sigc(k-2))*(sigc(k-1)-sigc(k)))
         Int_factor3 = (sig(k)-sigc(k-2))*(sig(k)-sigc(k-1))/  &
             ((sigc(k)-sigc(k-2))*(sigc(k)-sigc(k-1)))
         Uf(i,j,k) = Int_factor1*U(i,j,k-2)+Int_factor2*U(i,j,k-1)+Int_factor3*U(i,j,k)
         Vf(i,j,k) = Int_factor1*V(i,j,k-2)+Int_factor2*V(i,j,k-1)+Int_factor3*V(i,j,k)
         Wf(i,j,k) = Int_factor1*W(i,j,k-2)+Int_factor2*W(i,j,k-1)+Int_factor3*W(i,j,k)
       enddo
       enddo
       enddo

       do j = 1,Nloc
       do i = 1,Mloc
         Int_factor1 = (sig(2)-sigc(2))*(sig(2)-sigc(3))/  &
             ((sigc(1)-sigc(2))*(sigc(1)-sigc(2)))
         Int_factor2 = (sig(2)-sigc(1))*(sig(2)-sigc(3))/  &
             ((sigc(2)-sigc(1))*(sigc(2)-sigc(3)))
         Int_factor3 = (sig(2)-sigc(1))*(sig(2)-sigc(2))/  &
             ((sigc(3)-sigc(1))*(sigc(3)-sigc(2)))
         Uf(i,j,2) = Int_factor1*U(i,j,1)+Int_factor2*U(i,j,2)+Int_factor3*U(i,j,3)
         Vf(i,j,2) = Int_factor1*V(i,j,1)+Int_factor2*V(i,j,2)+Int_factor3*V(i,j,3)
         Wf(i,j,2) = Int_factor1*W(i,j,1)+Int_factor2*W(i,j,2)+Int_factor3*W(i,j,3)
         Uf(i,j,1) = U(i,j,1)
         Vf(i,j,1) = V(i,j,1)
         Wf(i,j,1) = W(i,j,1)
         Uf(i,j,Kloc1) = U(i,j,Kloc)
         Vf(i,j,Kloc1) = V(i,j,Kloc)
         Wf(i,j,Kloc1) = W(i,j,Kloc)
       enddo
       enddo
     endif

     end subroutine interpolate_velocity_to_faces


     subroutine get_UVW
!------------------------------------------------
!    Obtain U,V,W
!    Called by
!       eval_duvw
!    Last update: 25/12/2010, Gangfeng Ma
!-----------------------------------------------
     use global
     implicit none
     integer :: i,j,k

     do j = Jbeg,Jend
     do i = Ibeg,Iend
     do k = Kbeg,Kend
       U(i,j,k) = DU(i,j,k)/D(i,j)
       V(i,j,k) = DV(i,j,k)/D(i,j)
       W(i,j,k) = DW(i,j,k)/D(i,j)
     enddo
     enddo
     enddo

     ! collect data into ghost cells
     call vel_bc
# if defined (PARALLEL)
     call phi_3D_exch(U)
     call phi_3D_exch(V)
     call phi_3D_exch(W)
     call phi_3D_exch(DU)
     call phi_3D_exch(DV)
     call phi_3D_exch(DW)
# endif

!     ! debug
!     do k = Kbeg,Kend
!     do j = Jbeg,Jend
!     do i = Ibeg,Iend
!       if(U(i,j,k)>3.0.or.V(i,j,k)>3.0.or.W(i,j,k)>3.0) then
!         write(11,'(I3,20E15.5)') k,x(i),y(j),D(i,j),U(i,j,k),V(i,j,k),W(i,j,k), &
!             CmuVt(i,j,k),Tke(i,j,k),Eps(i,j,k)
!       endif
!     enddo
!     enddo
!     enddo

     end subroutine get_UVW


     subroutine fluxes
!------------------------------------------------
!    This subroutine is used to calculate fluxes 
!    at cell faces
!    Called by
!       main
!    Last update: 23/12/2010, Gangfeng Ma
!------------------------------------------------
     use global
     implicit none

     ! second order construction
     call delxyzFun
     call construction  

     ! calculate wave speed
     call wave_speed

     ! calculate fluxes at faces
     call fluxes_at_faces_HLL

     ! impose boundary conditions
     call flux_bc

     end subroutine fluxes


     subroutine construction
!------------------------------------------
!    Second-order construction
!    Called by 
!       fluxes
!    Last update: 04/01/2011, Gangfeng Ma
!-----------------------------------------
     use global
     implicit none
     integer :: i,j,k

     call construct_2D_x(Eta,DelxEta,EtaxL,EtaxR)
     call construct_3D_x(U,DelxU,UxL,UxR)
     call construct_3D_x(V,DelxV,VxL,VxR)
     call construct_3D_x(W,DelxW,WxL,WxR)
     call construct_3D_x(DU,DelxDU,DUxL,DUxR)
     call construct_3D_x(DV,DelxDV,DVxL,DVxR)
     call construct_3D_x(DW,DelxDW,DWxL,DWxR)

     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend1
       DxL(i,j) = EtaxL(i,j)+hfx(i,j)
       DxR(i,j) = EtaxR(i,j)+hfx(i,j)
       ExL(i,j,k) = DUxL(i,j,k)
       ExR(i,j,k) = DUxR(i,j,k)
       FxL(i,j,k) = DUxL(i,j,k)*UxL(i,j,k)+0.5*Grav*(EtaxL(i,j)*EtaxL(i,j)+2.0*EtaxL(i,j)*hfx(i,j))
       FxR(i,j,k) = DUxR(i,j,k)*UxR(i,j,k)+0.5*Grav*(EtaxR(i,j)*EtaxR(i,j)+2.0*EtaxR(i,j)*hfx(i,j))
       GxL(i,j,k) = DxL(i,j)*UxL(i,j,k)*VxL(i,j,k)
       GxR(i,j,k) = DxR(i,j)*UxR(i,j,k)*VxR(i,j,k)
       HxL(i,j,k) = DxL(i,j)*UxL(i,j,k)*WxL(i,j,k)
       HxR(i,j,k) = DxR(i,j)*UxR(i,j,k)*WxR(i,j,k)
     enddo
     enddo
     enddo

     call construct_2D_y(Eta,DelyEta,EtayL,EtayR)
     call construct_3D_y(U,DelyU,UyL,UyR)
     call construct_3D_y(V,DelyV,VyL,VyR)
     call construct_3D_y(W,DelyW,WyL,WyR)
     call construct_3D_y(DU,DelyDU,DUyL,DUyR)
     call construct_3D_y(DV,DelyDV,DVyL,DVyR)
     call construct_3D_y(DW,DelyDW,DWyL,DWyR)

     do k = Kbeg,Kend
     do j = Jbeg,Jend1
     do i = Ibeg,Iend
       DyL(i,j) = EtayL(i,j)+hfy(i,j)
       DyR(i,j) = EtayR(i,j)+hfy(i,j)
       EyL(i,j,k) = DVyL(i,j,k)
       EyR(i,j,k) = DVyR(i,j,k)
       FyL(i,j,k) = DyL(i,j)*UyL(i,j,k)*VyL(i,j,k)
       FyR(i,j,k) = DyR(i,j)*UyR(i,j,k)*VyR(i,j,k)
       GyL(i,j,k) = DVyL(i,j,k)*VyL(i,j,k)+0.5*Grav*(EtayL(i,j)*EtayL(i,j)+2.0*EtayL(i,j)*hfy(i,j))
       GyR(i,j,k) = DVyR(i,j,k)*VyR(i,j,k)+0.5*Grav*(EtayR(i,j)*EtayR(i,j)+2.0*EtayR(i,j)*hfy(i,j))
       HyL(i,j,k) = DyL(i,j)*VyL(i,j,k)*WyL(i,j,k)
       HyR(i,j,k) = DyR(i,j)*VyR(i,j,k)*WyR(i,j,k) 
     enddo
     enddo
     enddo

     call construct_3D_z(U,DelzU,UzL,UzR)
     call construct_3D_z(V,DelzV,VzL,VzR)
     call construct_3D_z(W,DelzW,WzL,WzR)
     call construct_3D_z(Omega,DelzO,OzL,OzR)

     end subroutine construction


     subroutine construct_2D_x(Vin,Din,OutL,OutR)
!-------------------------------------------------
!    Construct 2D variables in x-direction
!    Called by
!       construction
!    Last update: 04/01/2011, Gangfeng Ma
!------------------------------------------------
     use global, only: SP,Zero,dx,Mloc,Nloc,Mloc1, &
                       Ibeg,Iend,Jbeg,Jend,Iend1
     implicit none
     real(SP),intent(in),dimension(Mloc,Nloc)   :: Vin,Din
     real(SP),intent(out),dimension(Mloc1,Nloc) :: OutL,OutR      
     integer :: i,j

     OutL = Zero
     OutR = Zero
     do i = Ibeg,Iend1
     do j = Jbeg,Jend
       OutL(i,j) = Vin(i-1,j)+0.5*dx*Din(i-1,j)
       OutR(i,j) = Vin(i,j)-0.5*dx*Din(i,j)
     enddo
     enddo

     end subroutine construct_2D_x


     subroutine construct_3D_x(Vin,Din,OutL,OutR)
!-------------------------------------------------
!    Construct 3D variables in x-direction
!    Called by 
!       construction 
!    Last update: 04/01/2011, Gangfeng Ma 
!------------------------------------------------
     use global, only: SP,Zero,dx,Mloc,Nloc,Kloc,Mloc1, &
                       Ibeg,Iend,Jbeg,Jend,Iend1,Kbeg,Kend
     implicit none
     real(SP),intent(in),dimension(Mloc,Nloc,Kloc)   :: Vin,Din
     real(SP),intent(out),dimension(Mloc1,Nloc,Kloc) :: OutL,OutR      
     integer :: i,j,k

     OutL = Zero
     OutR = Zero
     do i = Ibeg,Iend1
     do j = Jbeg,Jend
     do k = Kbeg,Kend
       OutL(i,j,k) = Vin(i-1,j,k)+0.5*dx*Din(i-1,j,k)
       OutR(i,j,k) = Vin(i,j,k)-0.5*dx*Din(i,j,k)
     enddo
     enddo
     enddo

     end subroutine construct_3D_x


     subroutine construct_2D_y(Vin,Din,OutL,OutR)
!-------------------------------------------------
!    Construct 2D variables in y-direction
!    Called by
!       construction 
!    Last update: 04/01/2011, Gangfeng Ma 
!------------------------------------------------ 
     use global, only: SP,Zero,dy,Mloc,Nloc,Nloc1, &
                       Ibeg,Iend,Jbeg,Jend,Jend1
     implicit none
     real(SP),intent(in),dimension(Mloc,Nloc)   :: Vin,Din
     real(SP),intent(out),dimension(Mloc,Nloc1) :: OutL,OutR      
     integer :: i,j

     OutL = Zero
     OutR = Zero
     do i = Ibeg,Iend
     do j = Jbeg,Jend1
       OutL(i,j) = Vin(i,j-1)+0.5*dy*Din(i,j-1)
       OutR(i,j) = Vin(i,j)-0.5*dy*Din(i,j)
     enddo
     enddo

     end subroutine construct_2D_y


     subroutine construct_3D_y(Vin,Din,OutL,OutR)
!-------------------------------------------------
!    Construct 3D variables in y-direction 
!    Called by
!       construction 
!    Last update: 04/01/2011, Gangfeng Ma 
!------------------------------------------------
     use global, only: SP,Zero,dy,Mloc,Nloc,Kloc,Nloc1, &
                       Ibeg,Iend,Jbeg,Jend,Jend1,Kbeg,Kend
     implicit none
     real(SP),intent(in),dimension(Mloc,Nloc,Kloc)   :: Vin,Din
     real(SP),intent(out),dimension(Mloc,Nloc1,Kloc) :: OutL,OutR      
     integer :: i,j,k

     OutL = Zero
     OutR = Zero
     do i = Ibeg,Iend
     do j = Jbeg,Jend1
     do k = Kbeg,Kend
       OutL(i,j,k) = Vin(i,j-1,k)+0.5*dy*Din(i,j-1,k)
       OutR(i,j,k) = Vin(i,j,k)-0.5*dy*Din(i,j,k)
     enddo
     enddo
     enddo

     end subroutine construct_3D_y


     subroutine construct_3D_z(Vin,Din,OutL,OutR)
!-------------------------------------------------
!    Construct 3D variables in z-direction
!    Called by
!       construction 
!    Last update: 04/01/2011, Gangfeng Ma
!------------------------------------------------
     use global, only: SP,Zero,dsig,Mloc,Nloc,Kloc,Kloc1, &
                       Ibeg,Iend,Jbeg,Jend,Kbeg,Kend,Kend1
     implicit none
     real(SP),intent(in),dimension(Mloc,Nloc,Kloc)   :: Vin,Din
     real(SP),intent(out),dimension(Mloc,Nloc,Kloc1) :: OutL,OutR      
     integer :: i,j,k

     OutL = Zero
     OutR = Zero
     do i = Ibeg,Iend
     do j = Jbeg,Jend
     do k = Kbeg,Kend1
       OutL(i,j,k) = Vin(i,j,k-1)+0.5*dsig(k-1)*Din(i,j,k-1)
       OutR(i,j,k) = Vin(i,j,k)-0.5*dsig(k)*Din(i,j,k)
     enddo
     enddo
     enddo

     end subroutine construct_3D_z


     subroutine delxyzFun 
!-------------------------------------------
!    Calculate variable derivatives 
!    Called by 
!       fluxes 
!    Last update: 04/01/2011, Gangfeng Ma
!------------------------------------------
     use global
     implicit none
     integer :: i
     
     call delxFun_2D(Eta,DelxEta)
     call delxFun_3D(U,DelxU)
     call delxFun_3D(V,DelxV)
     call delxFun_3D(W,DelxW)
     call delxFun_3D(DU,DelxDU)
     call delxFun_3D(DV,DelxDV)
     call delxFun_3D(DW,DelxDW)

     call delyFun_2D(Eta,DelyEta)
     call delyFun_3D(U,DelyU)
     call delyFun_3D(V,DelyV)
     call delyFun_3D(W,DelyW)
     call delyFun_3D(DU,DelyDU)
     call delyFun_3D(DV,DelyDV)
     call delyFun_3D(DW,DelyDW)

     call delzFun_3D(U,DelzU)
     call delzFun_3D(V,DelzV)
     call delzFun_3D(W,DelzW)
     call delzFun_3D(Omega,DelzO)

     end subroutine delxyzFun

     
     subroutine delxFun_2D(Din,Dout)
!-------------------------------------------
!    Second-order derivative in x
!    Called by
!       delxyzFun
!    Last update: 04/01/2011, Gangfeng Ma
!------------------------------------------
     use global, only: SP,Small,Zero,dx,Mloc,Nloc,Mask
     implicit none
     real(SP),intent(in),dimension(Mloc,Nloc)  :: Din
     real(SP),intent(out),dimension(Mloc,Nloc) :: Dout
     real(SP) :: TMP1,TMP2,VANLEER_LIMITER
     integer :: i,j
    
     do i = 2,Mloc-1
     do j = 1,Nloc
       if(Mask(i,j)==0) then
         Dout(i,j) = Zero
       else
         if(Mask(i+1,j)==0) then
           TMP1 = Zero
         else
           TMP1 = (Din(i+1,j)-Din(i,j))/dx
         endif
         if(Mask(i-1,j)==0) then
           TMP2 = Zero
         else
           TMP2 = (Din(i,j)-Din(i-1,j))/dx
         endif

         if((abs(TMP1)+abs(TMP2))<Small) then
           Dout(i,j) = Zero
         else
           Dout(i,j) = VANLEER_LIMITER(TMP1,TMP2)
         endif
       endif
     enddo
     enddo

     do j = 1,Nloc
       Dout(1,j) = (Din(2,j)-Din(1,j))/dx
       Dout(Mloc,j) = (Din(Mloc,j)-Din(Mloc-1,j))/dx
     enddo  

     return
     end subroutine delxFun_2D


     subroutine delxFun_3D(Din,Dout)
!------------------------------------------
!    Second-order derivative in x
!    Called by 
!       delxyzFun
!    Last update: 04/01/2011, Gangfeng Ma
!------------------------------------------
     use global, only: SP,Small,Zero,dx,Mloc,Nloc,Kloc,Mask
     implicit none
     real(SP),intent(in),dimension(Mloc,Nloc,Kloc)  :: Din
     real(SP),intent(out),dimension(Mloc,Nloc,Kloc) :: Dout
     real(SP) :: TMP1,TMP2,VANLEER_LIMITER
     integer :: i,j,k

     do i = 2,Mloc-1
     do j = 1,Nloc
     do k = 1,Kloc
       if(Mask(i,j)==0) then
         Dout(i,j,k) = Zero
       else
         if(Mask(i+1,j)==0) then
           TMP1 = Zero
         else
           TMP1 = (Din(i+1,j,k)-Din(i,j,k))/dx
         endif
         if(Mask(i-1,j)==0) then
           TMP2 = Zero
         else
           TMP2 = (Din(i,j,k)-Din(i-1,j,k))/dx
         endif

         if((abs(TMP1)+abs(TMP2))<Small) then
           Dout(i,j,k) = Zero
         else
           Dout(i,j,k) = VANLEER_LIMITER(TMP1,TMP2)
         endif
       endif
     enddo
     enddo
     enddo

     do j = 1,Nloc
     do k = 1,Kloc
       Dout(1,j,k) = (Din(2,j,k)-Din(1,j,k))/dx
       Dout(Mloc,j,k) = (Din(Mloc,j,k)-Din(Mloc-1,j,k))/dx
     enddo
     enddo

     return
     end subroutine delxFun_3D


     subroutine delyFun_2D(Din,Dout)
!-----------------------------------------
!    Second-order derivative in y
!    Called by 
!       delxyzFun  
!    Last update: 04/01/2011, Gangfeng Ma  
!------------------------------------------ 
     use global, only: SP,Small,Zero,dy,Mloc,Nloc,Mask
     implicit none
     real(SP),intent(in),dimension(Mloc,Nloc)  :: Din
     real(SP),intent(out),dimension(Mloc,Nloc) :: Dout
     real(SP) :: TMP1,TMP2,VANLEER_LIMITER
     integer :: i,j

     do i = 1,Mloc
     do j = 2,Nloc-1
       if(Mask(i,j)==0) then 
         Dout(i,j) = Zero
       else
         if(Mask(i,j+1)==0) then
           TMP1 = Zero
         else
           TMP1 = (Din(i,j+1)-Din(i,j))/dy
         endif
         if(Mask(i,j-1)==0) then
           TMP2 = Zero
         else
           TMP2 = (Din(i,j)-Din(i,j-1))/dy
         endif

         if((abs(TMP1)+abs(TMP2))<Small) then
           Dout(i,j) = Zero
         else
           Dout(i,j) = VANLEER_LIMITER(TMP1,TMP2)
         endif
       endif
     enddo
     enddo

     do i = 1,Mloc
       Dout(i,1) = (Din(i,2)-Din(i,1))/dy
       Dout(i,Nloc) = (Din(i,Nloc)-Din(i,Nloc-1))/dy
     enddo

     return
     end subroutine delyFun_2D


     subroutine delyFun_3D(Din,Dout)
!-------------------------------------------
!    Second-order derivative in y
!    Called by
!       delxyzFun 
!    Last update: 04/01/2011, Gangfeng Ma 
!-------------------------------------------
     use global, only: SP,Small,Zero,dy,Mloc,Nloc,Kloc,Mask
     implicit none
     real(SP),intent(in),dimension(Mloc,Nloc,Kloc)  :: Din
     real(SP),intent(out),dimension(Mloc,Nloc,Kloc) :: Dout
     real(SP) :: TMP1,TMP2,VANLEER_LIMITER
     integer :: i,j,k

     do i = 1,Mloc
     do j = 2,Nloc-1
     do k = 1,Kloc
       if(Mask(i,j)==0) then
         Dout(i,j,k) = Zero
       else
         if(Mask(i,j+1)==0) then
           TMP2 = Zero
         else
           TMP1 = (Din(i,j+1,k)-Din(i,j,k))/dy
         endif
         if(Mask(i,j-1)==0) then 
           TMP2 = Zero
         else
           TMP2 = (Din(i,j,k)-Din(i,j-1,k))/dy
         endif

         if((abs(TMP1)+abs(TMP2))<Small) then
           Dout(i,j,k) = Zero
         else
           Dout(i,j,k) = VANLEER_LIMITER(TMP1,TMP2)
         endif
       endif
     enddo
     enddo
     enddo

     do i = 1,Mloc
     do k = 1,Kloc
       Dout(i,1,k) = (Din(i,2,k)-Din(i,1,k))/dy
       Dout(i,Nloc,k) = (Din(i,Nloc,k)-Din(i,Nloc-1,k))/dy
     enddo
     enddo

     return 
     end subroutine delyFun_3D

     
     subroutine delzFun_3D(Din,Dout)
!-------------------------------------------
!    Second-order derivative in z
!    Called by
!       delxyzFun
!    Last update: 04/01/2011, Gangfeng Ma
!-------------------------------------------
     use global, only: SP,Small,Zero,dsig,Mloc,Nloc,Kloc
     implicit none
     real(SP),intent(in),dimension(Mloc,Nloc,Kloc)  :: Din
     real(SP),intent(out),dimension(Mloc,Nloc,Kloc) :: Dout
     real(SP) :: TMP1,TMP2,VANLEER_LIMITER
     integer :: i,j,k

     do i = 1,Mloc
     do j = 1,Nloc
     do k = 2,Kloc-1
       TMP1 = (Din(i,j,k+1)-Din(i,j,k))/dsig(k)
       TMP2 = (Din(i,j,k)-Din(i,j,k-1))/dsig(k)

       if((abs(TMP1)+abs(TMP2))<Small) then
         Dout(i,j,k) = Zero
       else
         Dout(i,j,k) = VANLEER_LIMITER(TMP1,TMP2)
       endif
     enddo
     enddo
     enddo

     do i = 1,Mloc
     do j = 1,Nloc
       Dout(i,j,1) = (Din(i,j,2)-Din(i,j,1))/(0.5*(dsig(1)+dsig(2)))
       Dout(i,j,Kloc) = (Din(i,j,Kloc)-Din(i,j,Kloc-1))/(0.5*(dsig(Kloc-1)+dsig(Kloc)))
     enddo
     enddo

     return
     end subroutine delzFun_3D


     subroutine flux_bc
!--------------------------------------------
!    This is subroutine to provide boundary conditions
!    Called by
!       fluxes
!    Last update: 25/12/2010, Gangfeng Ma
!--------------------------------------------
     use global
     implicit none
     integer :: i,j,k

     ! left and right side
# if defined (PARALLEL)
     if(n_west.eq.MPI_PROC_NULL) then
# endif
     do j = Jbeg,Jend
     do k = Kbeg,Kend
       if(Bc_X0==1.or.Bc_X0==2) then
         Ex(Ibeg,j,k) = Zero
         Fx(Ibeg,j,k) = 0.5*Grav*(EtaxR(Ibeg,j)*EtaxR(Ibeg,j)+2.0*EtaxR(Ibeg,j)*hfx(Ibeg,j))
         Gx(Ibeg,j,k) = Zero
         Hx(Ibeg,j,k) = Zero
       elseif(Bc_X0==3) then
         Ex(Ibeg,j,k) = Din_X0(j)*Uin_X0(j,k)
         Fx(Ibeg,j,k) = Din_X0(j)*Uin_X0(j,k)*Uin_X0(j,k)+0.5*Grav*(Ein_X0(j)*Ein_X0(j)+2.0*Ein_X0(j)*hfx(Ibeg,j))
         Gx(Ibeg,j,k) = Din_X0(j)*Uin_X0(j,k)*Vin_X0(j,k)
! gfma fixed the bug here 12/19/2011
         Hx(Ibeg,j,k) = Din_X0(j)*Uin_X0(j,k)*Win_X0(j,k)
       endif
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_east.eq.MPI_PROC_NULL) then
# endif
     do j = Jbeg,Jend
     do k = Kbeg,Kend
       if(Bc_Xn==1.or.Bc_Xn==2) then
         Ex(Iend1,j,k) = Zero
         Fx(Iend1,j,k) = 0.5*Grav*(EtaxL(Iend1,j)*EtaxL(Iend1,j)+2.0*EtaxL(Iend1,j)*hfx(Iend1,j))
         Gx(Iend1,j,k) = Zero
         Hx(Iend1,j,k) = Zero

! fyshi add influx condition

       elseif(Bc_Xn==3) then
         Ex(Iend1,j,k) = Din_Xn(j)*Uin_Xn(j,k)
         Fx(Iend1,j,k) = Din_Xn(j)*Uin_Xn(j,k)*Uin_Xn(j,k)+0.5*Grav*(Ein_Xn(j)*Ein_Xn(j)+2.0*Ein_Xn(j)*hfx(Iend1,j))
         Gx(Iend1,j,k) = Din_Xn(j)*Uin_Xn(j,k)*Vin_Xn(j,k)
! gfma also fixed here 12/19/2011
         Hx(Iend1,j,k) = Din_Xn(j)*Uin_Xn(j,k)*Win_Xn(j,k)
       endif
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif

     ! front and back side
# if defined (PARALLEL)
     if(n_suth.eq.MPI_PROC_NULL) then
# endif
     do i = Ibeg,Iend
     do k = Kbeg,Kend
       if(Bc_Y0==1.or.Bc_Y0==2) then
         Ey(i,Jbeg,k) = Zero
         Fy(i,Jbeg,k) = Zero
         Gy(i,Jbeg,k) = 0.5*Grav*(EtayR(i,Jbeg)*EtayR(i,Jbeg)+2.0*EtayR(i,Jbeg)*hfy(i,Jbeg))
         Hy(i,Jbeg,k) = Zero
       endif
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_nrth.eq.MPI_PROC_NULL) then
# endif
     do i = Ibeg,Iend
     do k = Kbeg,Kend
       if(Bc_Yn==1.or.Bc_Yn==2) then
         Ey(i,Jend1,k) = Zero
         Fy(i,Jend1,k) = Zero
         Gy(i,Jend1,k) = 0.5*Grav*(EtayL(i,Jend1)*EtayL(i,Jend1)+2.0*EtayL(i,Jend1)*hfy(i,Jend1))
         Hy(i,Jend1,k) = Zero
       endif
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif

     ! upper and bottom
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       Fz(i,j,Kbeg) = Zero
       Gz(i,j,Kbeg) = Zero
       Hz(i,j,Kbeg) = Zero
       Fz(i,j,Kend1) = Zero
       Gz(i,j,Kend1) = Zero
       Hz(i,j,Kend1) = Zero
     enddo
     enddo

     do k = 1,Kloc
     do j = 1,Nloc
     do i = 1,Mloc
       if(Mask(i,j)==0) then
         Ex(i,j,k) = Zero
         Fx(i,j,k) = 0.5*Grav*(EtaxL(i,j)*EtaxL(i,j)+2.0*EtaxL(i,j)*hfx(i,j))*Mask(i-1,j)
         Gx(i,j,k) = Zero
         Hx(i,j,k) = Zero

         Ex(i+1,j,k) = Zero
         Fx(i+1,j,k) = 0.5*Grav*(EtaxR(i+1,j)*EtaxR(i+1,j)+2.0*EtaxR(i+1,j)*hfx(i+1,j))*Mask(i+1,j)
         Gx(i+1,j,k) = Zero
         Hx(i+1,j,k) = Zero

         Ey(i,j,k) = Zero
         Fy(i,j,k) = Zero
         Gy(i,j,k) = 0.5*Grav*(EtayL(i,j)*EtayL(i,j)+2.0*EtayL(i,j)*hfy(i,j))*Mask(i,j-1)
         Hy(i,j,k) = Zero

         Ey(i,j+1,k) = Zero
         Fy(i,j+1,k) = Zero
         Gy(i,j+1,k) = 0.5*Grav*(EtayR(i,j+1)*EtayR(i,j+1)+2.0*EtayR(i,j+1)*hfy(i,j+1))*Mask(i,j+1)
         Hy(i,j+1,k) = Zero
       endif
     enddo
     enddo
     enddo

     end subroutine flux_bc


     subroutine fluxes_at_faces_HLL
!---------------------------------------------
!    Fluxes at cell faces estimated by HLL approximation
!    Called by 
!       fluxes
!    Last update: 24/12/2010, Gangfeng Ma
!---------------------------------------------
     use global
     implicit none
     integer  :: i,j,k
     real(SP) :: Omega_f
     real(SP), dimension(:,:,:), allocatable :: D3xL,D3xR,D3yL,D3yR

     ! temporary arrays
     allocate(D3xL(Mloc1,Nloc,Kloc))
     allocate(D3xR(Mloc1,Nloc,Kloc))
     do k = 1,Kloc
     do j = 1,Nloc
     do i = 1,Mloc1
       D3xL(i,j,k) = EtaxL(i,j)
       D3xR(i,j,k) = EtaxR(i,j)
     enddo
     enddo
     enddo

     allocate(D3yL(Mloc,Nloc1,Kloc))
     allocate(D3yR(Mloc,Nloc1,Kloc))
     do k = 1,Kloc
     do j = 1,Nloc1
     do i = 1,Mloc
       D3yL(i,j,k) = EtayL(i,j)
       D3yR(i,j,k) = EtayR(i,j)
     enddo
     enddo
     enddo

     ! horizontal fluxes
     call HLL(Mloc1,Nloc,Kloc,SxL,SxR,ExL,ExR,D3xL,D3xR,Ex)
     call HLL(Mloc,Nloc1,Kloc,SyL,SyR,EyL,EyR,D3yL,D3yR,Ey)
     call HLL(Mloc1,Nloc,Kloc,SxL,SxR,FxL,FxR,DUxL,DUxR,Fx)
     call HLL(Mloc,Nloc1,Kloc,SyL,SyR,FyL,FyR,DUyL,DUyR,Fy)
     call HLL(Mloc1,Nloc,Kloc,SxL,SxR,GxL,GxR,DVxL,DVxR,Gx)
     call HLL(Mloc,Nloc1,Kloc,SyL,SyR,GyL,GyR,DVyL,DVyR,Gy)
     call HLL(Mloc1,Nloc,Kloc,SxL,SxR,HxL,HxR,DWxL,DWxR,Hx)
     call HLL(Mloc,Nloc1,Kloc,SyL,SyR,HyL,HyR,DWyL,DWyR,Hy)     
 
     ! vertical fluxes
     do k = Kbeg+1,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       Omega_f = 0.5*(OzL(i,j,k)+OzR(i,j,k))
       Fz(i,j,k) = 0.5*(Omega_f*(UzL(i,j,k)+UzR(i,j,k))-abs(Omega_f)*(UzR(i,j,k)-UzL(i,j,k)))
       Gz(i,j,k) = 0.5*(Omega_f*(VzL(i,j,k)+VzR(i,j,k))-abs(Omega_f)*(VzR(i,j,k)-VzL(i,j,k)))
       Hz(i,j,k) = 0.5*(Omega_f*(WzL(i,j,k)+WzR(i,j,k))-abs(Omega_f)*(WzR(i,j,k)-WzL(i,j,k)))
!       Omega_f = 0.5*(Omega(i,j,k-1)+Omega(i,j,k))
!       Fz(i,j,k) = 0.5*Omega_f*(U(i,j,k-1)+U(i,j,k))
!       Gz(i,j,k) = 0.5*Omega_f*(V(i,j,k-1)+V(i,j,k))
!       Hz(i,j,k) = 0.5*Omega_f*(W(i,j,k-1)+W(i,j,k))
     enddo
     enddo
     enddo

     deallocate(D3xL)
     deallocate(D3xR)
     deallocate(D3yL)
     deallocate(D3yR)

     return
     end subroutine fluxes_at_faces_HLL


     subroutine HLL(M,N,L,SL,SR,FL,FR,UL,UR,FOUT)
!----------------------------------------------
!    HLLC reconstruction 
!    Called by
!       fluxes_at_faces_HLL
!    Last update: 24/12/2010, Gangfeng Ma
!---------------------------------------------
     use global, only: SP,ZERO,SMALL
     implicit none
     INTEGER,INTENT(IN)::M,N,L
     REAL(SP),INTENT(IN),DIMENSION(M,N,L)::SL,SR,FL,FR,UL,UR
     REAL(SP),INTENT(OUT),DIMENSION(M,N,L)::FOUT
     INTEGER :: I,J,K

     DO K = 1,L
     DO J = 1,N
     DO I = 1,M
       IF(SL(I,J,K)>=ZERO) THEN
         FOUT(I,J,K) = FL(I,J,K)
       ELSEIF(SR(I,J,K)<=ZERO) THEN
         FOUT(I,J,K) = FR(I,J,K)
       ELSE
         FOUT(I,J,K) = SR(I,J,K)*FL(I,J,K)-SL(I,J,K)*FR(I,J,K)+  &
               SL(I,J,K)*SR(I,J,K)*(UR(I,J,K)-UL(I,J,K))
         IF((ABS(SR(I,J,K)-SL(I,J,K)))<SMALL)THEN
           FOUT(I,J,K) = FOUT(I,J,K)/SMALL
         ELSE
           FOUT(I,J,K) = FOUT(I,J,K)/(SR(I,J,K)-SL(I,J,K))
         ENDIF
       ENDIF
     ENDDO
     ENDDO
     ENDDO

     return
     end subroutine HLL


     subroutine wave_speed
!----------------------------------------------
!    This subroutine is used to calculate wave speeds
!    Called by
!       fluxes
!    Last update: 24/12/2010, Gangfeng Ma
!    Last update: 12/04/2011, Gangfeng Ma, wetting-drying
!-----------------------------------------------
     use global, only: SP,Ibeg,Iend,Iend1,Jbeg,Jend,Jend1,Kbeg,Kend, &
                       DxL,DxR,DyL,DyR,UxL,UxR,VyL,VyR, &
                       SxL,SxR,SyL,SyR,Grav,Mask
     implicit none
     integer  :: i,j,k
     real(SP) :: SQR_PHI_L,SQR_PHI_R,SQR_PHI_S,U_S
     
     ! x-faces
     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend1
       if(Mask(i-1,j)==1.and.Mask(i,j)==1) then
         SQR_PHI_L = sqrt(Grav*abs(DxL(i,j)))
         SQR_PHI_R = sqrt(Grav*abs(DxR(i,j)))
         SQR_PHI_S = 0.5*(SQR_PHI_L+SQR_PHI_R)+0.25*(UxL(i,j,k)-UxR(i,j,k))
         U_S = 0.5*(UxL(i,j,k)+UxR(i,j,k))+SQR_PHI_L-SQR_PHI_R
         SxL(i,j,k) = min(UxL(i,j,k)-SQR_PHI_L,U_S-SQR_PHI_S)
         SxR(i,j,k) = max(UxR(i,j,k)+SQR_PHI_R,U_S+SQR_PHI_S)
       elseif(Mask(i-1,j)==0.and.Mask(i,j)==1) then
         ! left-side dry case
         SQR_PHI_R = sqrt(Grav*abs(DxR(i,j)))
         SxL(i,j,k) = UxR(i,j,k)-2.0*SQR_PHI_R
         SxR(i,j,k) = UxR(i,j,k)+SQR_PHI_R
       elseif(Mask(i-1,j)==1.and.Mask(i,j)==0) then
         ! right-side dry case
         SQR_PHI_L = sqrt(Grav*abs(DxL(i,j)))
         SxL(i,j,k) = UxL(i,j,k)-SQR_PHI_L
         SxR(i,j,k) = UxL(i,j,k)+2.0*SQR_PHI_L
       endif
     enddo
     enddo
     enddo

     ! y-faces
     do k = Kbeg,Kend
     do j = Jbeg,Jend1
     do i = Ibeg,Iend
       if(Mask(i,j-1)==1.and.Mask(i,j)==1) then
         SQR_PHI_L = sqrt(Grav*abs(DyL(i,j)))
         SQR_PHI_R = sqrt(Grav*abs(DyR(i,j)))
         SQR_PHI_S = 0.5*(SQR_PHI_L+SQR_PHI_R)+0.25*(VyL(i,j,k)-VyR(i,j,k))
         U_S = 0.5*(VyL(i,j,k)+VyR(i,j,k))+SQR_PHI_L-SQR_PHI_R
         SyL(i,j,k) = min(VyL(i,j,k)-SQR_PHI_L,U_S-SQR_PHI_S)
         SyR(i,j,k) = max(VyR(i,j,k)+SQR_PHI_R,U_S+SQR_PHI_S)
       elseif(Mask(i,j-1)==0.and.Mask(i,j)==1) then
         ! left-side dry case
         SQR_PHI_R = sqrt(Grav*abs(DyR(i,j)))
         SyL(i,j,k) = VyR(i,j,k)-2.0*SQR_PHI_R
         SyR(i,j,k) = VyR(i,j,k)+SQR_PHI_R
       elseif(Mask(i,j-1)==1.and.Mask(i,j)==0) then
         ! right-side dry case
         SQR_PHI_L = sqrt(Grav*abs(DyL(i,j)))
         SyL(i,j,k) = VyL(i,j,k)-SQR_PHI_L
         SyR(i,j,k) = VyL(i,j,k)+2.0*SQR_PHI_L
       endif
     enddo
     enddo
     enddo

     end subroutine wave_speed


     FUNCTION VANLEER_LIMITER(A,B)
     use global, only: SP
     IMPLICIT NONE
     REAL(SP),INTENT(IN) :: A,B
     REAL(SP) :: VANLEER_LIMITER

     VANLEER_LIMITER=(A*ABS(B)+ABS(A)*B)/(ABS(A)+ABS(B))

     RETURN
     END FUNCTION VANLEER_LIMITER


     FUNCTION MINMOD_LIMITER_3V(A,B,C)
     use global, only: SP,Zero,One
     IMPLICIT NONE
     REAL(SP),INTENT(IN) :: A,B,C
     REAL(SP) :: MINMOD_LIMITER_3V

     MINMOD_LIMITER_3V=SIGN(One,A)*MAX(Zero,MIN(ABS(A),SIGN(One,A)*B,SIGN(One,A)*C))

     RETURN
     END FUNCTION MINMOD_LIMITER_3V


     FUNCTION MINMOD_LIMITER_2V(A,B)
     use global, only: SP,Zero,One
     IMPLICIT NONE
     REAL(SP),INTENT(IN) :: A,B
     REAL(SP) :: MINMOD_LIMITER_2V

     MINMOD_LIMITER_2V=SIGN(One,A)*MAX(Zero,MIN(ABS(A),SIGN(One,A)*B))

     RETURN
     END FUNCTION MINMOD_LIMITER_2V


     subroutine source_terms
!------------------------------------------------
!    This subroutine is used to evaluate source
!    Called by
!       main
!    Last update: 23/12/2010, Gangfeng Ma
!------------------------------------------------
     use global
     implicit none
     integer :: i,j,Iter,nn,ndir,nfreq,nk
     real(SP) :: Segma,Celerity,Wave_Length,Wave_Number,Fk,Fkdif,Source_Area,myvar, &
                 WnumX,WnumY,Phs_lag,dfreq,ddir,Angle,tmp1,tmp2,tmp3,tmp4
     real(SP) :: Ytrough,Mod1,Zup,Zlow,Zmid,Xstart,Zero1,cnoidal_cn,cnoidal_ck,Atmp

     ! internal wavemaker for linear wave
     if(WaveMaker(1:7)=='INT_LIN') then
       ! Find wave number for linear wave (Newton-Ralphson Method)
       Segma = 2.0*pi/Per_Wave
       Celerity = sqrt(Grav*Dep_Wave)
       Wave_Length = Celerity*Per_Wave
       Wave_Number = 2.0*pi/Wave_Length
     
       Iter = 0
 55    Fk = Grav*Wave_Number*tanh(Wave_Number*Dep_Wave)-Segma**2
       if(abs(Fk)<=1.0e-8.or.Iter>1000) goto 65
       Fkdif = Grav*Wave_Number*Dep_Wave*(1.0-tanh(Wave_Number*Dep_Wave)**2)+  &
           Grav*tanh(Wave_Number*Dep_Wave) 
       Wave_Number = Wave_Number-Fk/Fkdif
       Iter = Iter+1
       goto 55
 65    continue
       Wave_Length = 2.0*pi/Wave_Number
       Celerity = Wave_Length/Per_Wave
       WnumX = Wave_Number*cos(Theta_Wave*pi/180.)
       WnumY = Wave_Number*sin(Theta_Wave*pi/180.)       

       Source_Area = 0.0    
       do i = Ibeg,Iend
         if(xc(i)>=Xsource_West.and.xc(i)<=Xsource_East) then
           Source_Area = Source_Area+dx*D(i,Jbeg+1)
         endif
       enddo

# if defined (PARALLEL)
       call MPI_ALLREDUCE(Source_Area,myvar,1,MPI_SP,MPI_SUM,MPI_COMM_WORLD,ier)
       Source_Area = myvar
# endif

       do j = Jbeg,Jend
       do i = Ibeg,Iend
         if(xc(i)>=Xsource_West.and.xc(i)<=Xsource_East.and. &
            yc(j)>=Ysource_Suth.and.yc(j)<=Ysource_Nrth) then
           Phs_lag = (j-Jbeg)*dy*WnumY
           SourceC(i,j) = Celerity*Amp_Wave/Source_Area*cos(pi/2-Segma*time+Phs_lag)
         endif
       enddo
       enddo
     endif

     ! internal wavemaker for random waves 
     if(WaveMaker(1:7)=='INT_SPC') then
       Source_Area = 0.0
       do i = Ibeg,Iend
         if(xc(i)>=Xsource_West.and.xc(i)<=Xsource_East) then
           Source_Area = Source_Area+dx*D(i,Jbeg+1)
         endif
       enddo

# if defined (PARALLEL)
       call MPI_ALLREDUCE(Source_Area,myvar,1,MPI_SP,MPI_SUM,MPI_COMM_WORLD,ier)  
       Source_Area = myvar
# endif

       do j = Jbeg,Jend
       do i = Ibeg,Iend
         SourceC(i,j) = Zero
         if(xc(i)>=Xsource_West.and.xc(i)<=Xsource_East.and. &
              yc(j)>=Ysource_Suth.and.yc(j)<=Ysource_Nrth) then
           do nfreq = 1,NumFreq
           do ndir = 1,NumDir
             Per_Wave = 1.0/Freq(nfreq)
             Segma = 2.0*pi/Per_Wave
             Celerity = sqrt(Grav*Dep_Wave)
             Wave_Length = Celerity*Per_Wave
             Wave_Number = 2.0*pi/Wave_Length
       
             Iter = 0
   75        Fk = Grav*Wave_Number*tanh(Wave_Number*Dep_Wave)-Segma**2
             if(abs(Fk)<=1.0e-8.or.Iter>1000) goto 85
             Fkdif = Grav*Wave_Number*Dep_Wave*(1.0-tanh(Wave_Number*Dep_Wave)**2)+  & 
               Grav*tanh(Wave_Number*Dep_Wave)
             Wave_Number = Wave_Number-Fk/Fkdif
             Iter = Iter+1
             goto 75
   85        continue
             Wave_Length = 2.0*pi/Wave_Number
             Celerity = Wave_Length/Per_Wave

             ! adjust wave direction for periodic bc
             Angle = Dire(ndir)*pi/180.
             if(Angle>zero) then
               tmp3 = zero
               tmp1 = Wave_Number
               nk = 0
               do while (tmp3<Angle)
                 nk = nk+1
                 tmp2 = nk*2.0*pi/(Nglob*dy)
                 if(tmp2>=tmp1) then
                   tmp3 = 0.5*pi-small
                 else
                   tmp3 = asin(tmp2/tmp1)
                 endif
               enddo

               ! judge between nk-1 and nk which is closer                                          
               tmp4 = asin((nk-1)*2.0*pi/(Nglob*dy)/tmp1)
               if(abs(tmp4-Angle)<abs(Angle-tmp3)) then
                 Angle = tmp4
               else
                 Angle = tmp3
               endif
             else
               tmp3 = zero
               tmp1 = Wave_Number
               nk = 0
               do while (tmp3>Angle)
                 nk = nk+1
                 tmp2 = nk*2.0*pi/(Nglob*dy)
                 if(tmp2>=tmp1) then
                   tmp3 = -0.5*pi+small
                 else
                   tmp3 = -asin(tmp2/tmp1)
                 endif
               enddo

               ! judge between nk-1 and nk which is closer                                          
               tmp4= asin((nk-1)*2.0*pi/(Nglob*dy)/tmp1)
               if(abs(tmp4-Angle)<abs(Angle-tmp3)) then
                 Angle = tmp4
               else
                 Angle = tmp3
               endif
             endif

             WnumX = Wave_Number*cos(Angle)
             WnumY = Wave_Number*sin(Angle)

             ! calculate root-mean-squre wave height for each component
             if(nfreq==1) then
               dfreq = Freq(2)-Freq(1)
             elseif(nfreq==NumFreq) then
               dfreq = Freq(NumFreq)-Freq(NumFreq-1)
             else
               dfreq = 0.5*(Freq(nfreq+1)-Freq(nfreq-1))
             endif
             dfreq = abs(dfreq)

             if(ndir==1) then
               ddir = Dire(2)-Dire(1)
             elseif(ndir==NumDir) then
               ddir = Dire(NumDir)-Dire(NumDir-1)
             else
               ddir = 0.5*(Dire(ndir+1)-Dire(ndir-1))
             endif
             ddir = abs(ddir)
         
             Amp_Wave = 2.0*sqrt(2.0*Wave_Spc2d(ndir,nfreq)*ddir*dfreq)

             Phs_lag = (j-Jbeg)*dy*WnumY
             SourceC(i,j) = SourceC(i,j)+Celerity*Amp_Wave/Source_Area*  &
                  cos(pi/2-Segma*time+Phs_lag+Random_Phs(ndir,nfreq)) 
           enddo
           enddo
         endif
       enddo
       enddo
     endif

     ! internal wavemaker for cnoidal wave
     if(WaveMaker(1:7)=='INT_CON') then
       call cnoidal(Amp_Wave,Dep_Wave,Per_Wave,Wave_Length,Celerity,Ytrough,Mod1)

       ! wave number
       Wave_Number = 2.0*pi/Wave_Length

!# if defined(PARALLEL)
!       if(myid.eq.0) write(*,*) 'Mod=',Mod1,'Ytrough=',Ytrough, &
!            'Wave_Number=',wave_number
!# endif      

       ! find zero start
       Zup = 1.0
       Zlow = 0.0
       Zmid= (Zup+Zlow)/2.0
       nn = 0
 200   nn = nn+1
       Zero1 = Ytrough+Amp_Wave*cnoidal_cn(Zmid*0.5*cnoidal_ck(Mod1),Mod1)**2                            

       if(abs(Zero1)<=1.0e-6) goto 210
       if(nn>1000) then
         write(*,*)'too many iterations; stop'
         stop
       endif
       if(Zero1<0.0) then
         Zup = Zmid
         Zmid = (Zup+Zlow)/2.0
         goto 200
       else
         Zlow = Zmid
         Zmid = (Zup+Zlow)/2.0
         goto 200
       endif
 210   continue
       Xstart = Zmid

       Source_Area = 0.0
       do i = Ibeg,Iend
         if(xc(i)>=Xsource_West.and.xc(i)<=Xsource_East) then
           Source_Area = Source_Area+dx*D(i,Jbeg+1)
         endif
       enddo

# if defined (PARALLEL)
       call MPI_ALLREDUCE(Source_Area,myvar,1,MPI_SP,MPI_SUM,MPI_COMM_WORLD,ier) 
       Source_Area = myvar
# endif

       do j = Jbeg,Jend
       do i = Ibeg,Iend
         if(xc(i)>=Xsource_West.and.xc(i)<=Xsource_East.and. &
            yc(j)>=Ysource_Suth.and.yc(j)<=Ysource_Nrth) then
           SourceC(i,j) = 2.0*Celerity/Source_Area*(Ytrough+Amp_Wave*cnoidal_cn(  &
               Xstart*0.5*cnoidal_ck(Mod1)+2.0*cnoidal_ck(Mod1)*(-TIME/Per_Wave),Mod1)**2)
         endif
       enddo
       enddo
     endif

     ! internal wavemaker for solitary wave
     if(WaveMaker(1:7)=='INT_SOL') then
       Celerity = sqrt(Grav*Dep_Wave*(1.0+Amp_Wave/Dep_Wave))
       Atmp = sqrt(0.75*Amp_Wave/Dep_Wave**3)
       Xstart = 4.0*Dep_Wave/sqrt(Amp_Wave/Dep_Wave)
       
       Source_Area = 0.0
       do i = Ibeg,Iend
         if(xc(i)>=Xsource_West.and.xc(i)<=Xsource_East) then
           Source_Area = Source_Area+dx*D(i,Jbeg+1)
         endif
       enddo

# if defined (PARALLEL)
       call MPI_ALLREDUCE(Source_Area,myvar,1,MPI_SP,MPI_SUM,MPI_COMM_WORLD,ier) 
       Source_Area = myvar
# endif

       do j = Jbeg,Jend
       do i = Ibeg,Iend
         if(xc(i)>=Xsource_West.and.xc(i)<=Xsource_East.and. &
            yc(j)>=Ysource_Suth.and.yc(j)<=Ysource_Nrth) then
           SourceC(i,j) = 2.0*Celerity/Source_Area*  &
               Amp_Wave/cosh(Atmp*(Xstart-Celerity*TIME))**2
         endif
       enddo
       enddo
     endif

     ! source terms for momentum eqs.
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       SourceX(i,j) = Grav*Eta(i,j)*DelxH(i,j)*Mask(i,j)
       SourceY(i,j) = Grav*Eta(i,j)*DelyH(i,j)*Mask(i,j)
     enddo
     enddo

     end subroutine source_terms


     subroutine statistics
!---------------------------------------------------
!    This subroutine is used to show statistics
!    Called by
!       main
!    Last update: 23/12/2010, Gangfeng Ma
!--------------------------------------------------
     use global
     implicit none
     real(SP) :: MassVolume,CellMass,Energy,MaxEta,MinEta,MaxU, &
                 MaxV,MaxW
     integer :: i,j,k
# if defined (PARALLEL)
     real(SP) :: myvar
# endif

     ! Vol = sum(D*dx*dy)
     ! Energy = sum(m*g*h+0.5*m*u^2), reference is at z = 0
     MassVolume = Zero
     Energy = Zero
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       MassVolume = MassVolume+D(i,j)*dx*dy
       do k = Kbeg,Kend
         CellMass = Rho0*dsig(k)*D(i,j)*dx*dy
         Energy = Energy+CellMass*Grav*(D(i,j)*sigc(k)-Hc(i,j))+  &
                    0.5*CellMass*(U(i,j,k)**2+V(i,j,k)**2+W(i,j,k)**2)
       enddo
     enddo
     enddo

     MaxEta = MAXVAL(Eta(Ibeg:Iend,Jbeg:Jend))
     MinEta = MINVAL(Eta(Ibeg:Iend,Jbeg:Jend))
     MaxU = MAXVAL(abs(U(Ibeg:Iend,Jbeg:Jend,Kbeg:Kend)))
     MaxV = MAXVAL(abs(V(Ibeg:Iend,Jbeg:Jend,Kbeg:Kend)))
     MaxW = MAXVAL(abs(W(Ibeg:Iend,Jbeg:Jend,Kbeg:Kend))) 

# if defined (PARALLEL)
     call MPI_ALLREDUCE(MassVolume,myvar,1,MPI_SP,MPI_SUM,MPI_COMM_WORLD,ier)        
     MassVolume = myvar
     call MPI_ALLREDUCE(Energy,myvar,1,MPI_SP,MPI_SUM,MPI_COMM_WORLD,ier)            
     Energy = myvar
     call MPI_ALLREDUCE(MaxEta,myvar,1,MPI_SP,MPI_MAX,MPI_COMM_WORLD,ier)            
     MaxEta = myvar
     call MPI_ALLREDUCE(MinEta,myvar,1,MPI_SP,MPI_MIN,MPI_COMM_WORLD,ier)            
     MinEta = myvar
     call MPI_ALLREDUCE(MaxU,myvar,1,MPI_SP,MPI_MAX,MPI_COMM_WORLD,ier)
     MaxU = myvar
     call MPI_ALLREDUCE(MaxV,myvar,1,MPI_SP,MPI_MAX,MPI_COMM_WORLD,ier)
     MaxV = myvar
     call MPI_ALLREDUCE(MaxW,myvar,1,MPI_SP,MPI_MAX,MPI_COMM_WORLD,ier)
     MaxW = myvar
# endif

# if defined (PARALLEL)
     if(myid.eq.0) then
# endif
     ! print screen
     WRITE(*,*),'----------------- STATISTICS ----------------'
     WRITE(*,*),' TIME        DT         DT_CONSTRAINT'
     WRITE(*,102) TIME,dt,TRIM(dt_constraint)
     WRITE(*,*),' MassVolume  Energy      MaxEta      MinEta      Max U       Max V       MaxW'
     WRITE(*,101) MassVolume,Energy,MaxEta,MinEta,MaxU,MaxV,MaxW

     ! print log file 
     WRITE(3,*),'----------------- STATISTICS ----------------'
     WRITE(3,*),' TIME        DT         DT_CONSTRAINT'
     WRITE(3,102) TIME,dt,TRIM(dt_constraint)
     WRITE(3,*),' MassVolume  Energy      MaxEta      MinEta      Max U       Max V       MaxW'
     WRITE(3,101), MassVolume,Energy,MaxEta,MinEta,MaxU,MaxV,MaxW
# if defined (PARALLEL)
     endif
# endif

101  FORMAT(10E12.4)
102  FORMAT(2E12.4,A8)

     end subroutine statistics

 
     subroutine probes_old
!--------------------------------------------------
!    This subroutine is used to output probes
!    Called by
!       main
!    Last update: 16/11/2011, Gangfeng Ma
!--------------------------------------------------
     use global
     implicit none
     integer :: n,iu,i,j,k
     character(len=80) :: STAT_FILE,FDIR,FILE_NUM

     FDIR = TRIM(RESULT_FOLDER)
     
     do n = 1,NSTAT
       iu = 100+n
       write(FILE_NUM(1:4),'(I4.4)') n
       STAT_FILE = TRIM(FDIR)//'probe_'//TRIM(FILE_NUM)
       open(iu,file=TRIM(STAT_FILE),access='APPEND')
       
       do j = Jbeg,Jend
       do i = Ibeg,Iend
         if(xstat(n)>=x(i).and.xstat(n)<=x(i+1).and.  &
            ystat(n)>=y(j).and.ystat(n)<=y(j+1)) then
           write(iu,'(100E12.4)') time,eta(i,j),(u(i,j,k),v(i,j,k),w(i,j,k),k=Kbeg,Kend)
         endif
       enddo
       enddo
       close(iu)
     enddo

     end subroutine probes_old

     subroutine probes
!--------------------------------------------------
!    This subroutine is used to output probes
!    Called by
!       main
!    Last update: 28/09/20113 Stephan Grilli
!    (xc,yc) center coordinate of cell of station n
!    (xs,ys) (etas,us,vs)  cooridnates and surface and velocity of station n
!    nfun(4) linear shape functon in intrinsic coord (xis,eis)
!    Works for odd or even number of sigma layers since kmid = int((Kglob+1)/2)
!           and sigc(kmid+Nghost) = 0.5 if odd or 0.5-dsgm/2 if even
!    
!--------------------------------------------------
     use global
     implicit none
     integer :: n,iu,i,j,k,kmid,iflag(1000)
     real(SP) :: xco,yco,xs,ys,xis,eis,dsgm,csts,cst1,nfun(4),etacf(4),ucf(4),vcf(4),etas,us,vs
     character(len=80) :: STAT_FILE,FDIR,FILE_NUM

     FDIR = TRIM(RESULT_FOLDER)

!    Find "middle" sigma layer and initialize the interpolation to z = 0.531h
!
     kmid = Nghost + (Kglob+1)/2
     dsgm = 1./kglob
     csts = (0.531 - sigc(kmid))/dsgm
     cst1 = 1 - csts

!    Initialize station print flag
!
     do n = 1,NSTAT
        iflag(n) = 0
     enddo

!    Loop on stations
!
     do n = 1,NSTAT
       xs = xstat(n)
       ys = ystat(n)
       if (xs >= x(Ibeg).and.xs <= x(Iend).and.ys >= y(Jbeg).and. &
          ys <= y(Jend).and.iflag(n) == 0) then
          i  = Ibeg + (xs - x(Ibeg))/dx
          j  = Jbeg + (ys - y(Jbeg))/dy
          xco = x(i) + 0.5*dx
          yco = y(j) + 0.5*dy
          xis = 2*(xs - xco)/dx
          eis = 2*(ys - yco)/dy

!         Linear interpolation in (x,y) and in sigma direciton
!
          call funf2(xis,eis,nfun)
          etacf(1) = eta(i,j)
          etacf(2) = eta(i+1,j)
          etacf(3) = eta(i+1,j+1)
          etacf(4) = eta(i,j+1)
          ucf(1) = cst1*u(i,j,kmid) + csts*u(i,j,kmid+1)
          ucf(2) = cst1*u(i+1,j,kmid) + csts*u(i+1,j,kmid+1)
          ucf(3) = cst1*u(i+1,j+1,kmid) + csts*u(i+1,j+1,kmid+1)
          ucf(4) = cst1*u(i,j+1,kmid) + csts*u(i,j+1,kmid+1)
          vcf(1) = cst1*v(i,j,kmid) + csts*v(i,j,kmid+1)
          vcf(2) = cst1*v(i+1,j,kmid) + csts*v(i+1,j,kmid+1)
          vcf(3) = cst1*v(i+1,j+1,kmid) + csts*v(i+1,j+1,kmid+1)
          vcf(4) = cst1*v(i,j+1,kmid) + csts*v(i,j+1,kmid+1)

          etas = 0.
          us   = 0.
          vs   = 0.
          do k = 1,4
            etas = etas + nfun(k)*etacf(k)
            us   = us   + nfun(k)*ucf(k)
            vs   = vs   + nfun(k)*vcf(k)
          enddo


!         Save interpolated station data
!
          iu = 100+n
          write(FILE_NUM(1:4),'(I4.4)') n
          STAT_FILE = TRIM(FDIR)//'probe_'//TRIM(FILE_NUM(1:4))
          open(iu,file=TRIM(STAT_FILE),access='APPEND')
          write(iu,'(5E12.4)') time,etas,us,vs,w(i,j,kmid)
          close(iu)
          iflag(n) = 1
       endif
     enddo

     end subroutine probes

     subroutine funf2(xis,eis,nfun)
!--------------------------------------------------
!    This subroutine is used to output probes
!    Called by
!       main
!    Last update: 16/11/2011, Gangfeng Ma
!--------------------------------------------------
     use global
     implicit none
     real(SP) :: nfun(4),xis,eis
     nfun(1) = 0.25*(1 - xis)*(1 - eis)
     nfun(2) = 0.25*(1 + xis)*(1 - eis)
     nfun(3) = 0.25*(1 + xis)*(1 + eis)
     nfun(4) = 0.25*(1 - xis)*(1 + eis)

     end subroutine funf2


     subroutine preview
!--------------------------------------------------- 
!    This subroutine is used to preview
!    Called by                         
!       main 
!    Last update: 23/12/2010, Gangfeng Ma 
! fyshi gave an initial string to avoid unexpected characters
!  12/16/2011
!--------------------------------------------------
     use global
     implicit none
     REAL(SP),DIMENSION(Mloc,Nloc) :: ETA_MASKOUT,UV2D
     REAL(SP):: dsgm,csts,cst1
     INTEGER :: kmid
     integer :: i,j,k,I1,I2,I3,I4
     character(len=80) :: FDIR=''
     character(len=80) :: FILE_NAME=''
     character(len=80) :: file=''

     ! file number
     Icount = Icount+1
   
     ! results directory
     FDIR = TRIM(RESULT_FOLDER)

# if defined (PARALLEL)
     if(myid.eq.0) write(*,102) 'Printing file No.',Icount,' TIME/TOTAL: ',TIME,'/',TOTAL_TIME
     if(myid.eq.0) write(3,102) 'Printing file No.',Icount,' TIME/TOTAL: ',TIME,'/',TOTAL_TIME     
# else
     write(*,102) 'Printing file No.',Icount,' TIME/TOTAL: ',TIME,'/',TOTAL_TIME
     write(3,102) 'Printing file No.',Icount,' TIME/TOTAL: ',TIME,'/',TOTAL_TIME
# endif
102  FORMAT(A20,I4,A14,F8.3,A2,F8.3)
100  FORMAT(5000E16.6)

     I1 = mod(Icount/1000,10)
     I2 = mod(Icount/100,10)
     I3 = mod(Icount/10,10)
     I4 = mod(Icount,10)

     write(FILE_NAME(1:1),'(I1)') I1
     write(FILE_NAME(2:2),'(I1)') I2
     write(FILE_NAME(3:3),'(I1)') I3
     write(FILE_NAME(4:4),'(I1)') I4

# if defined (PARALLEL)
     if(myid.eq.0) then
# endif
     open(5,file=TRIM(FDIR)//'time',position="append")
     write(5,*) TIME
     close(5)
# if defined (PARALLEL)
     endif
# endif

     if(Icount==1) then
       if(OUT_H) then
         file=TRIM(FDIR)//'depth'
         call putfile2D(file,Hc)
       endif
     endif

! fyshi add tmp4preview2d and 3d

     IF(OUT_preview)THEN
       file=TRIM(FDIR)//'tmp2d_'//TRIM(FILE_NAME)
       call putfile2D(file,tmp4preview2D)

       file=TRIM(FDIR)//'tmp3d_'//TRIM(FILE_NAME)
       call putfile3D(file,tmp4preview3D)      
     ENDIF

     if(OUT_E) then
       ETA_MASKOUT=ETA*MASK
       file=TRIM(FDIR)//'eta_'//TRIM(FILE_NAME)
       call putfile2D(file,ETA_MASKOUT)

       file=TRIM(FDIR)//'hmax_'//TRIM(FILE_NAME)
       call putfile2D(file,EleMax)

       file=TRIM(FDIR)//'hmin_'//TRIM(FILE_NAME)
       call putfile2D(file,EleMin)

       file=TRIM(FDIR)//'FrMax_'//TRIM(FILE_NAME)
       call putfile3D(file,FrMax)

       file=TRIM(FDIR)//'FrMaxT_'//TRIM(FILE_NAME)
       call putfile3D(file,FrMax_total)

       file=TRIM(FDIR)//'Depfm_'//TRIM(FILE_NAME)
       call putfile2D(file,DepMaxFr)

       file=TRIM(FDIR)//'Timefm_'//TRIM(FILE_NAME)
       call putfile2D(file,TimeMaxFr)

     endif

     if(OUT_SLD) then		                   !DJN

       file=TRIM(FDIR)//'slide_'//TRIM(FILE_NAME)  !DJN
       call putfile2D(file,SldDZ1)		   !DJN

       file=TRIM(FDIR)//'slideU_'//TRIM(FILE_NAME)  !fyshi
       call putfile2D(file,SldUs)

       file=TRIM(FDIR)//'slideV_'//TRIM(FILE_NAME)  !fyshi
       call putfile2D(file,SldVs)

     endif		                           !DJN

     if(OUT_U) then
       file=TRIM(FDIR)//'u_'//TRIM(FILE_NAME)
       call putfile3D(file,U)
     endif

     if(OUT_V) then
       file=TRIM(FDIR)//'v_'//TRIM(FILE_NAME)
       call putfile3D(file,V)
     endif

     if(OUT_W) then
       file=TRIM(FDIR)//'w_'//TRIM(FILE_NAME)
       call putfile3D(file,W)
     endif

     IF(OUT_2DUV) THEN

        kmid = Nghost + (Kglob+1)/2
        dsgm = 1./kglob
        csts = (0.531 - sigc(kmid))/dsgm
        cst1 = 1 - csts

        DO J=1,Nloc
        DO I=1,Mloc
          UV2D(I,J) = cst1*U(i,j,kmid) + csts*U(i,j,kmid+1)
        ENDDO
        ENDDO
       
       file=TRIM(FDIR)//'u2d_'//TRIM(FILE_NAME)
       call putfile2D(file,UV2D)


        DO J=1,Nloc
        DO I=1,Mloc
          UV2D(I,J) = cst1*V(i,j,kmid) + csts*V(i,j,kmid+1)
        ENDDO
        ENDDO

       file=TRIM(FDIR)//'v2d_'//TRIM(FILE_NAME)
       call putfile2D(file,UV2D)

     ENDIF  ! end uv2d

     if(OUT_P) then
       file=TRIM(FDIR)//'p_'//TRIM(FILE_NAME)
       call putfile3D(file,P)
     endif

     if(OUT_K) then
       file=TRIM(FDIR)//'k_'//TRIM(FILE_NAME)
       call putfile3D(file,Tke)
     endif

     if(OUT_D) then
       file=TRIM(FDIR)//'d_'//TRIM(FILE_NAME)
       call putfile3D(file,Eps)
     endif

     if(OUT_S) then
       file=TRIM(FDIR)//'s_'//TRIM(FILE_NAME)
       call putfile3D(file,Prod_s)
     endif

     if(OUT_C) then
       file=TRIM(FDIR)//'c_'//TRIM(FILE_NAME)
       call putfile3D(file,CmuVt)
     endif

# if defined (BUBBLE)
     if(OUT_B) then
       file=TRIM(FDIR)//'b_'//TRIM(FILE_NAME)
       call putfile3D(file,Vbg)
     endif
# endif

     end subroutine preview


# if defined (PARALLEL)
    subroutine putfile2D(file,phi)
    use global
    implicit none
    real(SP),dimension(Mloc,Nloc),intent(in) :: phi
    character(len=80) :: file
    integer,dimension(NumP) :: npxs,npys
    real(SP),dimension(NumP) :: xx
    real(SP),dimension(Mglob,Nglob) :: phiglob
    integer :: i,j,iglob,jglob,l

    call MPI_GATHER(npx,1,MPI_INTEGER,npxs,1,MPI_INTEGER,  &
           0,MPI_COMM_WORLD,ier)
    call MPI_GATHER(npy,1,MPI_INTEGER,npys,1,MPI_INTEGER,  &
           0,MPI_COMM_WORLD,ier)

    do j = Jbeg,Jend
    do i = Ibeg,Iend
      call MPI_GATHER(phi(i,j),1,MPI_SP, &
            xx,1,MPI_SP,0,MPI_COMM_WORLD,ier)
      if(i.eq.Ibeg) call MPI_Barrier(MPI_COMM_WORLD,ier)
      if(myid.eq.0) then
        do l = 1,NumP
          iglob = npxs(l)*(Iend-Ibeg+1)+i-Nghost
          jglob = npys(l)*(Jend-Jbeg+1)+j-Nghost
          phiglob(iglob,jglob) = xx(l)
        enddo
      endif
    enddo
    enddo

    if(myid.eq.0) then
      open(5,file=TRIM(file))
      do j = 1,Nglob
        write(5,100) (phiglob(i,j),i=1,Mglob)
      enddo
      close(5)
    endif
100 FORMAT(5000E16.6)
    end subroutine putfile2D

    subroutine putfile3D(file,phi)
    use global
    implicit none
    real(SP),dimension(Mloc,Nloc,Kloc),intent(in) :: phi
    character(len=80) :: file
    integer,dimension(NumP) :: npxs,npys
    real(SP),dimension(NumP) :: xx
    real(SP),dimension(Mglob,Nglob,Kglob) :: phiglob
    integer :: i,j,k,iglob,jglob,kk,l

    call MPI_GATHER(npx,1,MPI_INTEGER,npxs,1,MPI_INTEGER,  &
          0,MPI_COMM_WORLD,ier)
    call MPI_GATHER(npy,1,MPI_INTEGER,npys,1,MPI_INTEGER,  &
          0,MPI_COMM_WORLD,ier)

    do k = Kbeg,Kend
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      call MPI_GATHER(phi(i,j,k),1,MPI_SP, &
            xx,1,MPI_SP,0,MPI_COMM_WORLD,ier)
      if (i.eq.Ibeg) call MPI_Barrier(MPI_COMM_WORLD,ier)
      if(myid.eq.0) then
        do l = 1,NumP
          iglob = npxs(l)*(Iend-Ibeg+1)+i-Nghost
          jglob = npys(l)*(Jend-Jbeg+1)+j-Nghost
          kk = k-Nghost
          phiglob(iglob,jglob,kk) = xx(l)
        enddo
      endif
    enddo
    enddo
    enddo

    if(myid.eq.0) then
      open(5,file=TRIM(file))
      do k = 1,Kglob
      do j = 1,Nglob
        write(5,100) (phiglob(i,j,k),i=1,Mglob)
      enddo
      enddo
      close(5)
    endif
100 FORMAT(5000E16.6)
    end subroutine putfile3D
# else
     subroutine putfile2D(file,phi)
     use global
     implicit none
     real(SP),dimension(Mloc,Nloc),intent(in) :: phi
     character(len=80) :: file
     integer :: i,j

     open(5,file=trim(file))

     do j = Jbeg,Jend
       write(5,100) (phi(i,j),i=Ibeg,Iend)
     enddo
     close(5)
 100 FORMAT(5000E16.6)
     end subroutine putfile2D

     subroutine putfile3D(file,phi)
     use global
     implicit none
     real(SP),dimension(Mloc,Nloc,Kloc),intent(in) :: phi
     character(len=80) :: file
     integer :: i,j,k

     open(5,file=trim(file))
     do k = Kbeg,Kend
     do j = Jbeg,Jend
       write(5,100) (phi(i,j,k),i=Ibeg,Iend)
     enddo
     enddo
     close(5)
 100 FORMAT(5000E16.6)
     end subroutine putfile3D  
# endif

 
     subroutine estimate_dt
!----------------------------------------------------
!    This subroutine is used to estimate dt
!    Called by
!       main
!    Last update: 22/12/2010, Gangfeng Ma
!---------------------------------------------------
     use global
     implicit none
     integer :: i,j,k
     real(SP) :: tmp1,tmp2,dxonu,dyonv,dzonw,dt_growth,dt_courant,dt_viscous
     real(SP) :: dt_landslide  !DJN
# if defined (PARALLEL)
     real(SP) :: myvar
# endif

     ! save previous time step
     dt_old = dt
     dt_growth = 1.05*dt_old     

     tmp2 = Large
     do k = 1,Kloc
     do j = 1,Nloc
     do i = 1,Mloc
       tmp1 = abs(U(i,j,k))+sqrt(Grav*D(i,j))
       tmp1 = max(tmp1,Small)
       dxonu = dx/tmp1
       if(dxonu<tmp2) tmp2=dxonu

       tmp1 = abs(V(i,j,k))+sqrt(Grav*D(i,j))
       tmp1 = max(tmp1,Small)
       dyonv = dy/tmp1
       if(dyonv<tmp2) tmp2=dyonv

       tmp1 = max(abs(W(i,j,k)),Small)
       dzonw = dsig(k)*D(i,j)/tmp1
       if(dzonw<tmp2) tmp2=dzonw
     enddo
     enddo
     enddo
# if defined (PARALLEL)
     call MPI_ALLREDUCE(tmp2,myvar,1,MPI_SP,MPI_MIN,MPI_COMM_WORLD,ier)
     tmp2 = myvar
# endif
     dt_courant = CFL*tmp2

!!*DJN
     ! time step limit due to the landslide 
     tmp2 = Large
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       tmp1 = abs(SldUs(i,j))+sqrt(Grav*SldDZs(i,j))
       tmp1 = max(tmp1,Small)
       dxonu = dx/tmp1
       if(dxonu<tmp2) tmp2=dxonu

       tmp1 = abs(SldVs(i,j))+sqrt(Grav*SldDZs(i,j))
       tmp1 = max(tmp1,Small)
       dyonv = dy/tmp1
       if(dyonv<tmp2) tmp2=dyonv
     enddo
     enddo
# if defined (PARALLEL)
     call MPI_ALLREDUCE(tmp2,myvar,1,MPI_SP,MPI_MIN,MPI_COMM_WORLD,ier)
     tmp2 = myvar
# endif
     dt_landslide = 0.1d0*CFL*tmp2
!!*DJN

     ! time step limit due to explicit viscous stress terms
     dt_viscous = Large
     if(VISCOUS_FLOW) then
       tmp2 = Large
       do k = 1,Kloc
       do j = 1,Nloc
       do i = 1,Mloc
         tmp1 = dx**2/abs(CmuHt(i,j,k))
         if(tmp1<tmp2) tmp2 = tmp1

         tmp1 = dy**2/abs(CmuHt(i,j,k))
         if(tmp1<tmp2) tmp2 = tmp1

         tmp1 = (D(i,j)*dsig(k))**2/abs(CmuVt(i,j,k))
         if(tmp1<tmp2) tmp2 = tmp1
       enddo
       enddo
       enddo
# if defined (PARALLEL)
       call MPI_ALLREDUCE(tmp2,myvar,1,MPI_SP,MPI_MIN,MPI_COMM_WORLD,ier)
       tmp2 = myvar
# endif
       dt_viscous = VISCOUS_NUMBER*tmp2
     endif 

     ! get dt    
     dt = min(dt_growth,dt_courant,dt_viscous,dt_landslide,dt_max)
     if(dt<dt_min) then
# if defined (PARALLEL)
       if(myid.eq.0) then
         write(3,*) 'time step too small !!',dt,dt_courant,dt_viscous
         stop
       endif
# else
       write(3,*) 'time step too small !!',dt,dt_courant,dt_viscous
       stop
# endif
     endif
     TIME = TIME+dt
     RUN_STEP = RUN_STEP+1 
# if defined (PARALLEL)
     if(myid.eq.0) write(3,*) RUN_STEP,dt,TIME
# else 
     write(3,*) RUN_STEP,dt,TIME
# endif

     if(dt==dt_growth) then
       dt_constraint = 'GROWTH'
     elseif(dt==dt_courant) then
       dt_constraint = 'COURANT'
     elseif(dt==dt_viscous) then
       dt_constraint = 'VISCOUS'
     elseif(dt==dt_max) then
       dt_constraint = 'MAXIMUM'
     elseif(dt==dt_landslide) then  !DJN
       dt_constraint = 'LANDSLD'    !DJN	
     endif  

     end subroutine estimate_dt


     subroutine vel_bc
!----------------------------------------------------
!    Boundary conditions for velocity
!    Called by 
!       main and get_UVW
!    Last update: 01/02/2011, Gangfeng Ma
!---------------------------------------------------
     use global
     implicit none
     integer :: i,j,k
     real(SP) :: Wtop,Wbot,Cdrag,Phi,Kappa,Dz1

     ! left and right boundary
# if defined (PARALLEL)
     if(n_west.eq.MPI_PROC_NULL) then
# endif
     do k = Kbeg,Kend
     do j = Jbeg,Jend
       if(Bc_X0==1) then  ! free-slip wall
         do i = 1,Nghost
           U(Ibeg-i,j,k) = -U(Ibeg+i-1,j,k)
           V(Ibeg-i,j,k) = V(Ibeg+i-1,j,k)
           W(Ibeg-i,j,k) = W(Ibeg+i-1,j,k)
           DU(Ibeg-i,j,k) = -DU(Ibeg+i-1,j,k)
           DV(Ibeg-i,j,k) = DV(Ibeg+i-1,j,k)
           DW(Ibeg-i,j,k) = DW(Ibeg+i-1,j,k)
         enddo
       elseif(Bc_X0==2) then ! no-slip wall
         do i =1,Nghost
           U(Ibeg-i,j,k) = -U(Ibeg+i-1,j,k)
           V(Ibeg-i,j,k) = -V(Ibeg+i-1,j,k)
           W(Ibeg-i,j,k) = -W(Ibeg+i-1,j,k)
           DU(Ibeg-i,j,k) = -DU(Ibeg+i-1,j,k)
           DV(Ibeg-i,j,k) = -DV(Ibeg+i-1,j,k)
           DW(Ibeg-i,j,k) = -DW(Ibeg+i-1,j,k)
         enddo
       elseif(Bc_X0==3) then ! inflow
         do i = 1,Nghost
           U(Ibeg-i,j,k) = 2.0*Uin_X0(j,k)-U(Ibeg,j,k)
           V(Ibeg-i,j,k) = 2.0*Vin_X0(j,k)-V(Ibeg,j,k)
           W(Ibeg-i,j,k) = 2.0*Win_X0(j,k)-W(Ibeg,j,k)
           DU(Ibeg-i,j,k) = 2.0*Uin_X0(j,k)*Din_X0(j)-DU(Ibeg,j,k)
           DV(Ibeg-i,j,k) = 2.0*Vin_X0(j,k)*Din_X0(j)-DV(Ibeg,j,k)
           DW(Ibeg-i,j,k) = 2.0*Win_X0(j,k)*Din_X0(j)-DW(Ibeg,j,k)
         enddo
       elseif(Bc_X0==4) then ! outflow
         do i = 1,Nghost
           U(Ibeg-i,j,k) = U(Ibeg+i-1,j,k)
           V(Ibeg-i,j,k) = V(Ibeg+i-1,j,k)
           W(Ibeg-i,j,k) = W(Ibeg+i-1,j,k)
           DU(Ibeg-i,j,k) = DU(Ibeg+i-1,j,k)
           DV(Ibeg-i,j,k) = DV(Ibeg+i-1,j,k)
           DW(Ibeg-i,j,k) = DW(Ibeg+i-1,j,k)
         enddo
       endif
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_east.eq.MPI_PROC_NULL) then
# endif
     do k = Kbeg,Kend
     do j = Jbeg,Jend
       if(Bc_Xn==1) then  ! free-slip wall 
         do i = 1,Nghost
           U(Iend+i,j,k) = -U(Iend-i+1,j,k)
           V(Iend+i,j,k) = V(Iend-i+1,j,k)
           W(Iend+i,j,k) = W(Iend-i+1,j,k)
           DU(Iend+i,j,k) = -DU(Iend-i+1,j,k)
           DV(Iend+i,j,k) = DV(Iend-i+1,j,k)
           DW(Iend+i,j,k) = DW(Iend-i+1,j,k)
         enddo
       elseif(Bc_Xn==2) then ! no-slip wall
         do i = 1,Nghost
           U(Iend+i,j,k) = -U(Iend-i+1,j,k)
           V(Iend+i,j,k) = -V(Iend-i+1,j,k)
           W(Iend+i,j,k) = -W(Iend-i+1,j,k)
           DU(Iend+i,j,k) = -DU(Iend-i+1,j,k)
           DV(Iend+i,j,k) = -DV(Iend-i+1,j,k)
           DW(Iend+i,j,k) = -DW(Iend-i+1,j,k)
         enddo
       elseif(Bc_Xn==3) then ! inflow
         do i = 1,Nghost
           U(Iend+i,j,k) = 2.0*Uin_Xn(j,k)-U(Iend-i+1,j,k)
           V(Iend+i,j,k) = 2.0*Vin_Xn(j,k)-V(Iend-i+1,j,k)
           W(Iend+i,j,k) = 2.0*Win_Xn(j,k)-W(Iend-i+1,j,k)
           DU(Iend+i,j,k) = 2.0*Uin_Xn(j,k)*Din_Xn(j)-DU(Iend-i+1,j,k)
           DV(Iend+i,j,k) = 2.0*Vin_Xn(j,k)*Din_Xn(j)-DV(Iend-i+1,j,k)
           DW(Iend+i,j,k) = 2.0*Win_Xn(j,k)*Din_Xn(j)-DW(Iend-i+1,j,k)
         enddo
       elseif(Bc_Xn==4) then ! outflow
         do i = 1,Nghost
           U(Iend+i,j,k) = U(Iend-i+1,j,k)
           V(Iend+i,j,k) = V(Iend-i+1,j,k)
           W(Iend+i,j,k) = W(Iend-i+1,j,k)
           DU(Iend+i,j,k) = DU(Iend-i+1,j,k)
           DV(Iend+i,j,k) = DV(Iend-i+1,j,k)
           DW(Iend+i,j,k) = DW(Iend-i+1,j,k)
         enddo
       endif
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_suth.eq.MPI_PROC_NULL) then
# endif
     do k = Kbeg,Kend
     do i = Ibeg,Iend
       if(Bc_Y0==1) then  ! free-slip wall 
         do j = 1,Nghost
           U(i,Jbeg-j,k) = U(i,Jbeg+j-1,k)
           V(i,Jbeg-j,k) = -V(i,Jbeg+j-1,k)
           W(i,Jbeg-j,k) = W(i,Jbeg+j-1,k)
           DU(i,Jbeg-j,k) = DU(i,Jbeg+j-1,k)
           DV(i,Jbeg-j,k) = -DV(i,Jbeg+j-1,k)
           DW(i,Jbeg-j,k) = DW(i,Jbeg+j-1,k)
         enddo
       elseif(Bc_Y0==2) then ! no-slip wall 
         do j = 1,Nghost
           U(i,Jbeg-j,k) = -U(i,Jbeg+j-1,k)
           V(i,Jbeg-j,k) = -V(i,Jbeg+j-1,k)
           W(i,Jbeg-j,k) = -W(i,Jbeg+j-1,k)
           DU(i,Jbeg-j,k) = -DU(i,Jbeg+j-1,k)
           DV(i,Jbeg-j,k) = -DV(i,Jbeg+j-1,k)
           DW(i,Jbeg-j,k) = -DW(i,Jbeg+j-1,k)
         enddo
       endif
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_nrth.eq.MPI_PROC_NULL) then
# endif
     do k = Kbeg,Kend
     do i = Ibeg,Iend
       if(Bc_Yn==1) then  ! free-slip wall 
         do j = 1,Nghost
           U(i,Jend+j,k) = U(i,Jend-j+1,k)
           V(i,Jend+j,k) = -V(i,Jend-j+1,k)
           W(i,Jend+j,k) = W(i,Jend-j+1,k)
           DU(i,Jend+j,k) = DU(i,Jend-j+1,k)
           DV(i,Jend+j,k) = -DV(i,Jend-j+1,k)
           DW(i,Jend+j,k) = DW(i,Jend-j+1,k)
         enddo
       elseif(Bc_Yn==2) then ! no-slip wall 
         do j =1,Nghost
           U(i,Jend+j,k) = -U(i,Jend-j+1,k)
           V(i,Jend+j,k) = -V(i,Jend-j+1,k)
           W(i,Jend+j,k) = -W(i,Jend-j+1,k)
           DU(i,Jend+j,k) = -DU(i,Jend-j+1,k)
           DV(i,Jend+j,k) = -DV(i,Jend-j+1,k)
           DW(i,Jend+j,k) = -DW(i,Jend-j+1,k)
         enddo
       endif
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif

     ! top and bottom
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       Kappa = 0.41
       Dz1 = D(i,j)*dsig(Kbeg)
       if(ibot==1) then
         Cdrag = Cd0
       else
         Cdrag = 1./(1./Kappa*log(15.0*Dz1/Zob))**2
       endif
       Phi = Dz1*Cdrag*sqrt(U(i,j,Kbeg)**2+V(i,j,Kbeg)**2)/CmuVt(i,j,Kbeg)
       Phi = dmin1(Phi,2.0)

       if(Bc_Z0==1) then  ! free-slip
         Wbot = -DeltH(i,j)-U(i,j,Kbeg)*DelxH(i,j)-V(i,j,Kbeg)*DelyH(i,j)
         do k = 1,Nghost
           U(i,j,Kbeg-k) = U(i,j,Kbeg+k-1)
           V(i,j,Kbeg-k) = V(i,j,Kbeg+k-1)
           W(i,j,Kbeg-k) = 2.0*Wbot-W(i,j,Kbeg+k-1)
           DU(i,j,Kbeg-k) = D(i,j)*U(i,j,Kbeg-k)
           DV(i,j,Kbeg-k) = D(i,j)*V(i,j,Kbeg-k)
           DW(i,j,Kbeg-k) = D(i,j)*W(i,j,Kbeg-k)
         enddo
       elseif(Bc_Z0==2) then  ! no-slip
         Wbot = -DeltH(i,j)
         do k = 1,Nghost
           U(i,j,Kbeg-k) = -U(i,j,Kbeg+k-1)
           V(i,j,Kbeg-k) = -V(i,j,Kbeg+k-1)
           W(i,j,Kbeg-k) = 2.0*Wbot-W(i,j,Kbeg+k-1)
           DU(i,j,Kbeg-k) = D(i,j)*U(i,j,Kbeg-k)
           DV(i,j,Kbeg-k) = D(i,j)*V(i,j,Kbeg-k)
           DW(i,j,Kbeg-k) = D(i,j)*W(i,j,Kbeg-k)
         enddo
       elseif(Bc_Z0==5) then
         do k = 1,Nghost
           U(i,j,Kbeg-k) = (1.0-Phi)*U(i,j,Kbeg+k-1)
           V(i,j,Kbeg-k) = (1.0-Phi)*V(i,j,Kbeg+k-1)
           Wbot = -DeltH(i,j)-0.5*(U(i,j,Kbeg)+U(i,j,Kbeg-1))*DelxH(i,j)-  &
                    0.5*(V(i,j,Kbeg)+V(i,j,Kbeg-1))*DelyH(i,j)
           W(i,j,Kbeg-k) = 2.0*Wbot-W(i,j,Kbeg+k-1)
           DU(i,j,Kbeg-k) = D(i,j)*U(i,j,Kbeg-k)
           DV(i,j,Kbeg-k) = D(i,j)*V(i,j,Kbeg-k)
           DW(i,j,Kbeg-k) = D(i,j)*W(i,j,Kbeg-k)
         enddo
       endif

       ! at the surface (no stress)
       Wtop = (Eta(i,j)-Eta0(i,j))/dt+U(i,j,Kend)*DelxEta(i,j)+V(i,j,Kend)*DelyEta(i,j)
       do k = 1,Nghost
         U(i,j,Kend+k) = U(i,j,Kend-k+1)
         V(i,j,Kend+k) = V(i,j,Kend-k+1)
         W(i,j,Kend+k) = 2.0*Wtop-W(i,j,Kend-k+1)
         DU(i,j,Kend+k) = D(i,j)*U(i,j,Kend+k)
         DV(i,j,Kend+k) = D(i,j)*V(i,j,Kend+k)
         DW(i,j,Kend+k) = D(i,j)*W(i,j,Kend+k)
       enddo
     enddo
     enddo

     end subroutine vel_bc


     subroutine phi_2D_coll(phi)
!-----------------------------------------------------
!    This subroutine is used to collect data into ghost cells
!    Called by
!       eval_duvw
!    Last update: 22/12/2010, Gangfeng Ma
!-----------------------------------------------------
     use global
     implicit none
     real(SP), intent(inout) :: phi(Mloc,Nloc)
     integer :: i,j

     ! x-direction
# if defined (PARALLEL)
     if(n_west.eq.MPI_PROC_NULL) then
# endif
       do j = Jbeg,Jend
       do i = 1,Nghost
         phi(i,j) = phi(Ibeg+Nghost-i,j)
       enddo
       enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_east.eq.MPI_PROC_NULL) then
# endif
       do j = Jbeg,Jend
       do i = 1,Nghost
         phi(Iend+i,j) = phi(Iend-i+1,j)
       enddo
       enddo
# if defined (PARALLEL)
     endif
# endif
 
     ! y-direction and corners
# if defined (PARALLEL)
     if(n_suth.eq.MPI_PROC_NULL) then
# endif     
       do i = 1,Mloc
       do j = 1,Nghost
         phi(i,j) = phi(i,Jbeg+Nghost-j)
       enddo
       enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_nrth.eq.MPI_PROC_NULL) then
# endif
       do i = 1,Mloc
       do j = 1,Nghost
         phi(i,Jend+j) = phi(i,Jend-j+1)
       enddo
       enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     call phi_2D_exch(phi)
# endif    

     end subroutine phi_2D_coll


     subroutine update_mask
!------------------------------------------------------  
!    This subroutine is used to update mask for wetting-drying
!    Called by                                                
!       main
!    Last update: 22/12/2010, Gangfeng Ma 
!-----------------------------------------------------
     use global, only: Ibeg,Iend,Jbeg,Jend,Eta,Hc,D,MinDep,  &
                       Mask,Mask_Struct,Mask9
     implicit none
     integer :: i,j

     ! Mask at ghost cells keeps no change
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       if(Mask_Struct(i,j)==0) cycle
       
       ! flooding (dry->wet)
       if(Mask(i,j)==0) then
         if(Mask(i-1,j)==1.and.Eta(i-1,j)>Eta(i,j)) Mask(i,j)=1
         if(Mask(i+1,j)==1.and.Eta(i+1,j)>Eta(i,j)) Mask(i,j)=1
         if(Mask(i,j-1)==1.and.Eta(i,j-1)>Eta(i,j)) Mask(i,j)=1
         if(Mask(i,j+1)==1.and.Eta(i,j+1)>Eta(i,j)) Mask(i,j)=1
       else
         ! drying (wet->dry)
         if(abs(D(i,j)-MinDep)<=1.e-6) then
           Mask(i,j) = 0
           Eta(i,j) = MinDep-Hc(i,j)
           D(i,j) = Eta(i,j)+Hc(i,j)           
         endif
       endif
     enddo
     enddo

     Mask = Mask*Mask_Struct

# if defined (PARALLEL)
     ! collect mask into ghost cells
     call phi_int_exch(Mask)    
# endif

     do j = Jbeg,Jend
     do i = Ibeg,Iend
      Mask9(i,j) = Mask(i,j)*Mask(i-1,j)*Mask(i+1,j)  &
                *Mask(i+1,j+1)*Mask(i,j+1)*Mask(i-1,j+1) &
                *Mask(i+1,j-1)*Mask(i,j-1)*Mask(i-1,j-1)
     enddo
     enddo

! fyshi modify mask9 at domain boundaries


!# if defined (PARALLEL)
!      if(n_west.eq.MPI_PROC_NULL) then
!# endif
!      MASK9(Ibeg,:)=0
!# if defined (PARALLEL)
!      endif
!# endif

!# if defined (PARALLEL)
!      if(n_east.eq.MPI_PROC_NULL) then
!# endif
!      MASK9(Iend,:)=0
!# if defined (PARALLEL)
!      endif
!# endif

     end subroutine update_mask


     subroutine update_vars
!------------------------------------------------------ 
!    This subroutine is used to save variables at 
!    last time step
!    Called by   
!       main 
!    Last update: 22/12/2010, Gangfeng Ma 
!----------------------------------------------------- 
     use global        
     implicit none

     D0 = D
     Eta0 = Eta
     DU0 = DU
     DV0 = DV
     DW0 = DW
     DTke0 = DTke
     DEps0 = DEps

# if defined (BUBBLE)
     DNbg0 = DNbg
# endif

     end subroutine update_vars

  
     subroutine update_wave_bc
!------------------------------------------------------
!    This subroutine is used to update boundary conditions
!    Called by
!       main
!    Last update: 22/12/2010, Gangfeng Ma 
!-----------------------------------------------------
     use global
     implicit none

# if defined (PARALLEL)
     if(n_west.eq.MPI_PROC_NULL) then
# endif
     if(WaveMaker(1:7)=='LEF_SOL') then
       call solitary_wave_left_boundary
     elseif(WaveMaker(1:7)=='LEF_LIN') then
       call linear_wave_left_boundary
     elseif(WaveMaker(1:7)=='LEF_CON') then
       call cnoidal_wave_left_boundary
     elseif(WaveMaker(1:7)=='LEF_STK') then
       call stokes_wave_left_boundary
     elseif(WaveMaker(1:7)=='LEF_SPC') then
       ! not work yet
       call random_wave_left_boundary
     elseif(WaveMaker(1:7)=='LEF_TID') then
       call tidal_wave_left_boundary
     elseif(WaveMaker(1:10)=='TID_FLX_LR') then
       call tidal_flux_left_right
     elseif(WaveMaker(1:10)=='TID_ELE_LR') then
       call tidal_elevation_left_right
     endif
# if defined (PARALLEL)
     endif
# endif

     end subroutine update_wave_bc 

     subroutine stokes_wave_left_boundary
!-----------------------------------------------------
!    This subroutine is used to specify left boundary
!    Called by
!       update_wave_bc
!    Last update: 26/04/2011, Gangfeng Ma
!-----------------------------------------------------
     use global, only: SP,pi,Zero,Ibeg,Grav,TIME,Nloc,Kloc,Amp_Wave,Per_Wave,Dep_Wave, &
                       Ein_X0,Din_X0,Uin_X0,Vin_X0,Win_X0,Hfx,Jbeg,Jend,Kbeg,Kend,sigc
     implicit none
     integer  :: j,k,Iter
     real(SP) :: Segma,Celerity,Wave_Length,Wave_Number,Fk,Fkdif,Zlev

     ! Find wave number for linear wave (Newton-Ralphson Method)
     Segma = 2.0*pi/Per_Wave
     Celerity = sqrt(Grav*Dep_Wave)
     Wave_Length = Celerity*Per_Wave
     Wave_Number = 2.0*pi/Wave_Length
     
     Iter = 0
55   Fk = Grav*Wave_Number*tanh(Wave_Number*Dep_Wave)-Segma**2
     if(abs(Fk)<=1.0e-8.or.Iter>1000) goto 65
     Fkdif = Grav*Wave_Number*Dep_Wave*(1.0-tanh(Wave_Number*Dep_Wave)**2)+  &
        Grav*tanh(Wave_Number*Dep_Wave) 
     Wave_Number = Wave_Number-Fk/Fkdif
     Iter = Iter+1
     goto 55
65   continue
     
     Wave_Length = 2.0*pi/Wave_Number
     Celerity = Wave_Length/Per_Wave
     
     do j = 1,Nloc
       Ein_X0(j) = 0.5*Amp_Wave*cos(pi/2-Segma*TIME)+  &
                   Amp_Wave**2*Wave_Number/16.*cosh(Wave_Number*Dep_Wave)/  &
                   sinh(Wave_Number*Dep_Wave)**3*(2.0+cosh(2.*Wave_Number*Dep_Wave))* &
                   cos(2.*(pi/2-Segma*TIME))
       Din_X0(j) = Ein_X0(j)+Hfx(Ibeg,j)     
     enddo

     do k = Kbeg,Kend
     do j = Jbeg,Jend
       Zlev = sigc(k)*Din_X0(j)
       Uin_X0(j,k) = 0.5*Amp_Wave*Segma*cosh(Wave_Number*Zlev)/  &
           sinh(Wave_Number*Dep_Wave)*cos(pi/2-Segma*TIME)+  &
           3./16.*Amp_Wave**2*Segma*Wave_Number*cosh(2.*Wave_Number*Zlev)/  &
           sinh(Wave_Number*Dep_Wave)**4*cos(2.*(pi/2-Segma*TIME))
       Win_X0(j,k) = 0.5*Amp_Wave*Segma*sinh(Wave_Number*Zlev)/  &
           sinh(Wave_Number*Dep_Wave)*sin(pi/2-Segma*TIME)+  &
           3./16.*Amp_Wave**2*Segma*Wave_Number*sinh(2.*Wave_Number*Zlev)/  &
           sinh(Wave_Number*Dep_Wave)**4*sin(2.*(pi/2-Segma*TIME))
       Vin_X0(j,k) = 0.0
     enddo
     enddo

     end subroutine stokes_wave_left_boundary


     subroutine tidal_wave_left_boundary
!-----------------------------------------------------
!    This subroutine is used to specify left boundary
!    Called by
!       update_wave_bc
!    Last update: 06/02/2011, Gangfeng Ma
!-----------------------------------------------------
     use global, only: SP,pi,Zero,Grav,Ibeg,Nloc,Eta,Hc,D,TIME
     implicit none
     integer :: j

     do j = 1,Nloc
       Eta(Ibeg,j) = 1.0*sin(2.0*pi*(2.0*TIME/86400.0))
       D(Ibeg,j) = Eta(Ibeg,j)+Hc(Ibeg,j)
     enddo

     end subroutine tidal_wave_left_boundary

     subroutine tidal_flux_left_right
!-----------------------------------------------------
!    This subroutine is used to specify left boundary
!    Called by
!       update_wave_bc
!    Last update: 12/16/2011, Fengyan Shi
!-----------------------------------------------------
# if defined (PARALLEL)
     use global, only: SP,pi,Zero,Ibeg,TIME,Nloc,Kloc, &
                       Ein_X0,Din_X0,Uin_X0,Vin_X0,Win_X0,Hfx,Jbeg,Jend,Kbeg,Kend, &
                       Ein_Xn,Din_Xn,Uin_Xn,Vin_Xn,Win_Xn,Iend1,Eta,D,Iend,Hc,MASK,&
                       icount_tide, NumTimeData, &
                       DataU_L,DataEta_L,DataU_R,DataEta_R,TimeData, &
                       n_west,n_east,MPI_PROC_NULL, &
                       U,Iend,Nudging
# else
     use global, only: SP,pi,Zero,Ibeg,TIME,Nloc,Kloc, &
                       Ein_X0,Din_X0,Uin_X0,Vin_X0,Win_X0,Hfx,Jbeg,Jend,Kbeg,Kend, &
                       Ein_Xn,Din_Xn,Uin_Xn,Vin_Xn,Win_Xn,Iend1,Eta,D,Iend,Hc,MASK,&
                       icount_tide, NumTimeData, &
                       DataU_L,DataEta_L,DataU_R,DataEta_R,TimeData, &
                       U,Iend,Nudging
# endif
 
     implicit none
     integer  :: j,k,Iter
     real(SP) :: Segma,Celerity,Wave_Length,Wave_Number,Fk,Fkdif,Zlev,Stokes_Drift,Fact,Phs_Lag(Nloc)
     REAL(SP) :: TMP1,TMP2
     REAL(SP) :: Ztmp1,Ztmp2,Utmp1,Utmp2

! process time series
      IF(TIME>TimeData(icount_tide).AND.icount_tide<NumTimeData)THEN
       icount_tide=icount_tide+1
      ENDIF
     IF(icount_tide>1)THEN  ! tide start
       IF(TIME>TimeData(icount_tide))THEN
        tmp2=ZERO
       ELSE
        tmp2=(TimeData(icount_tide)-TIME) &
           /(TimeData(icount_tide)-TimeData(icount_tide-1))
       ENDIF

       Ztmp1= DataEta_L(icount_tide)*(1.0_SP-tmp2)   &
             +DataEta_L(icount_tide-1)*tmp2
       Ztmp2= DataEta_R(icount_tide)*(1.0_SP-tmp2)   &
             +DataEta_R(icount_tide-1)*tmp2

       Utmp1= DataU_L(icount_tide)*(1.0_SP-tmp2)   &
             +DataU_L(icount_tide-1)*tmp2
       Utmp2= DataU_R(icount_tide)*(1.0_SP-tmp2)   &
             +DataU_R(icount_tide-1)*tmp2
     ELSE

        Ztmp1=0.0
        Ztmp2=0.0
        Utmp1=0.0
        Utmp2=0.0

     ENDIF  ! end tide start

! NOTE always nudge surface elevation
     do j = 1,Nloc
       Ein_X0(j) = Ztmp1
       Din_X0(j) = Ein_X0(j)+Hfx(Ibeg,j)   
       Ein_Xn(j) = Ztmp2
       Din_Xn(j) = Ein_Xn(j)+Hfx(Iend1,j)  
     enddo

     do k = Kbeg,Kend
     do j = Jbeg,Jend
       IF(Nudging.OR.Utmp1>0.0)THEN
         Uin_X0(j,k)=Utmp1*MASK(Ibeg,J)
       ELSE
         Uin_X0(j,k)= U(Ibeg,j,k)
       ENDIF
       Win_X0(j,k) = 0.0
       Vin_X0(j,k) = 0.0
       IF(Nudging.OR.Utmp2<0.0)THEN
         Uin_Xn(j,k)= Utmp2*MASK(Iend,J)
       ELSE
         Uin_Xn(j,k) = U(Iend,j,k)
       ENDIF
       Win_Xn(j,k) = 0.0
       Vin_Xn(j,k) = 0.0
     enddo
     enddo  

! force elevation at left and right boundaries
! to avoid over-shooting 
# if defined (PARALLEL)
     if(n_west.eq.MPI_PROC_NULL) then
# endif
     do j = 1,Nloc
       IF(MASK(Ibeg,J)==1)THEN
         Eta(Ibeg,j) = Ein_X0(j)
         D(Ibeg,j) = Eta(Ibeg,j)+Hc(Ibeg,j)
       ENDIF
     enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_east.eq.MPI_PROC_NULL) then
# endif
     do j = 1,Nloc
       IF(MASK(Iend,J)==1)THEN
         Eta(Iend,j) = Ein_Xn(j)
         D(Iend,j) = Eta(Iend,j)+Hc(Iend,j)
       ENDIF
     enddo
# if defined (PARALLEL)
     endif
# endif

     end subroutine tidal_flux_left_right


     subroutine tidal_elevation_left_right
!-----------------------------------------------------
!    This subroutine is used to specify left and right boundary
!    Called by
!       update_wave_bc
!    Last update: 12/16/2011, Fengyan Shi
!    NOTE: only for large-scale domain because it may causing 
!          unreasonable upwelling and downwelling
!-----------------------------------------------------
# if defined (PARALLEL)
     use global, only: SP,pi,Zero,Ibeg,TIME,Nloc,Kloc, &
                       Ein_X0,Din_X0,Uin_X0,Vin_X0,Win_X0,Hfx,Jbeg,Jend,Kbeg,Kend, &
                       Ein_Xn,Din_Xn,Uin_Xn,Vin_Xn,Win_Xn,Iend1,Eta,D,Iend,Hc,MASK, &
                       icount_tide, NumTimeData, &
                       DataU_L,DataEta_L,DataU_R,DataEta_R,TimeData, &
                       n_west,n_east,MPI_PROC_NULL
# else
     use global, only: SP,pi,Zero,Ibeg,TIME,Nloc,Kloc, &
                       Ein_X0,Din_X0,Uin_X0,Vin_X0,Win_X0,Hfx,Jbeg,Jend,Kbeg,Kend, &
                       Ein_Xn,Din_Xn,Uin_Xn,Vin_Xn,Win_Xn,Iend1,Eta,D,Iend,Hc,MASK, &
                       icount_tide, NumTimeData, &
                       DataU_L,DataEta_L,DataU_R,DataEta_R,TimeData
# endif
     implicit none
     integer  :: j,k,Iter
     real(SP) :: Segma,Celerity,Wave_Length,Wave_Number,Fk,Fkdif,Zlev,Stokes_Drift,Fact,Phs_Lag(Nloc)

     REAL(SP) :: TMP1,TMP2
     REAL(SP) :: Ztmp1,Ztmp2

! process time series
      IF(TIME>TimeData(icount_tide).AND.icount_tide<NumTimeData)THEN
       icount_tide=icount_tide+1
      ENDIF
     IF(icount_tide>1)THEN  ! tide start
       IF(TIME>TimeData(icount_tide))THEN
        tmp2=ZERO
       ELSE
        tmp2=(TimeData(icount_tide)-TIME) &
           /(TimeData(icount_tide)-TimeData(icount_tide-1))
       ENDIF

       Ztmp1= DataEta_L(icount_tide)*(1.0_SP-tmp2)   &
             +DataEta_L(icount_tide-1)*tmp2
       Ztmp2= DataEta_R(icount_tide)*(1.0_SP-tmp2)   &
             +DataEta_R(icount_tide-1)*tmp2
     ELSE

        Ztmp1=0.0
        Ztmp2=0.0

     ENDIF  ! end tide start


# if defined (PARALLEL)
     if(n_west.eq.MPI_PROC_NULL) then
# endif
     do j = 1,Nloc
       IF(MASK(Ibeg,J)==1)THEN
         Eta(Ibeg,j) = Ztmp1
         D(Ibeg,j) = Eta(Ibeg,j)+Hc(Ibeg,j)
       ENDIF
     enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_east.eq.MPI_PROC_NULL) then
# endif
     do j = 1,Nloc
       IF(MASK(Iend,J)==1)THEN
         Eta(Iend,j) = Ztmp2
         D(Iend,j) = Eta(Iend,j)+Hc(Iend,j)
       ENDIF
     enddo
# if defined (PARALLEL)
     endif
# endif


     end subroutine tidal_elevation_left_right

     subroutine cnoidal_wave_left_boundary
!-----------------------------------------------------
!    This subroutine is used to specify left boundary 
!    Called by 
!       update_wave_bc 
!    Last update: 06/02/2011, Gangfeng Ma 
!----------------------------------------------------- 
     use global, only: SP,pi,Zero,Ibeg,Grav,TIME,Nloc,Kloc,Amp_Wave,Per_Wave,Dep_Wave, &
                       Ein_X0,Din_X0,Uin_X0,Vin_X0,Win_X0,Hfx,Jbeg,Jend,Kbeg,Kend,sigc
     implicit none
     integer :: j,k,nn
     real(SP) :: Wave_Length,Celerity,Ytrough,Mod1,Wave_Number,Zup,Zlow,Zmid,Xstart, &
                 Zero1,Zlev,Xcn,Xsn,Xdn,cnoidal_cn,cnoidal_ck,Stokes_Drift,Fact

     call cnoidal(Amp_Wave,Dep_Wave,Per_Wave,Wave_Length,Celerity,Ytrough,Mod1)
     
     ! wave number
     Wave_Number = 2.0*pi/Wave_Length

     ! find zero start
     Zup = 1.0
     Zlow = 0.0
     Zmid= (Zup+Zlow)/2.0
     nn = 0
200  nn = nn+1
     Zero1 = Ytrough+Amp_Wave*cnoidal_cn(Zmid*0.5*cnoidal_ck(Mod1),Mod1)**2                            

     if(abs(Zero1)<=1.0e-6) goto 210
     if(nn>1000) then
       write(*,*)'too many iterations; stop'
       stop
     endif
     if(Zero1<0.0) then
       Zup = Zmid
       Zmid = (Zup+Zlow)/2.0
       goto 200
     else
       Zlow = Zmid
       Zmid = (Zup+Zlow)/2.0
       goto 200
     endif
 210 continue
     Xstart = Zmid
!     write(*,*) "Mod=",Mod1,"Ytrough=",Ytrough,"Wave_Number=",wave_number,  &    
!                   "Xstart=",Xstart

     do j = 1,Nloc
       Ein_X0(j) = Ytrough+Amp_Wave*cnoidal_cn(Xstart*0.5*cnoidal_ck(Mod1)+  &
             2.0*cnoidal_ck(Mod1)*(-TIME/Per_Wave),Mod1)**2
       Din_X0(j) = Ein_X0(j)+Hfx(Ibeg,j)
     enddo

     ! mean mass transport
     Stokes_Drift = Grav*Amp_Wave**2/(Wave_Length/Per_Wave)/Dep_Wave/8.0
     Fact = 1.0

     do k = Kbeg,Kend
     do j = Jbeg,Jend
       Zlev = sigc(k)*Din_X0(j)

       Xcn = cnoidal_cn(Xstart*0.5*cnoidal_ck(Mod1)+2.0*cnoidal_ck(Mod1)*(-TIME/Per_Wave),Mod1)
       Xsn = sqrt(1.0-Xcn**2)
       Xdn = sqrt(1.0-Mod1*(1.0-Xcn**2)) 

       Uin_X0(j,k) = sqrt(Grav*Dep_Wave)*(-5./4.+3.*(Dep_Wave+Ytrough)/2.0/Dep_Wave-  &
            (Dep_Wave+Ytrough)**2/4./Dep_Wave**2+(3.*Amp_Wave/2./Dep_Wave-  &                       
            (Dep_Wave+Ytrough)*Amp_Wave/2./Dep_Wave**2)*Xcn**2-Amp_Wave**2/4./  &                      
            Dep_Wave**2*Xcn**4-8.0*Amp_Wave*cnoidal_ck(Mod1)**2/Wave_Length**2*(Dep_Wave/3-  &            
            Zlev**2/2/Dep_Wave)*(-Mod1**2*Xsn**2*Xcn**2+Xcn**2*Xdn**2-Xsn**2*Xdn**2))-  &
            ! substract mean mass transport
            Stokes_Drift*Fact
       Win_X0(j,k) = sqrt(Grav*Dep_Wave)*Zlev*2*Amp_Wave*cnoidal_ck(Mod1)/  &
            Wave_Length/Dep_Wave*(1+(Dep_Wave+Ytrough)/Dep_Wave+Amp_Wave/  &                        
            Dep_Wave*Xcn**2+32*cnoidal_ck(Mod1)**2/3/Wave_Length**2*(Dep_Wave**2-  &
            Zlev**2/2)*(Mod1**2*Xsn**2-Mod1**2*Xcn**2-Xdn**2))*Xsn*Xcn*Xdn
       Vin_X0(j,k) =0.0
     enddo
     enddo

     end subroutine cnoidal_wave_left_boundary


     subroutine cnoidal(height,depth,period,l,c,x2,mod)
!----------------------------------------------------------------
!    Cnoidal Function 
!    Called by
!       cnoidal_wave_left_boundary 
!    Last update: 06/02/2011, Gangfeng Ma 
!----------------------------------------------------------------
     use global, only: SP,Grav
     implicit none

     real(SP), intent(in) :: height, depth, period
     real(SP), intent(out) :: l, c, x2, mod
     real(SP) :: xa,xb,xtemp,cnoidal_cn,cnoidal_ck,cnoidal_ce
     integer :: n

     mod=0.99999999d0
     n=0
 40  n=n+1

     xa=mod*depth+2.0d0*height-mod*height-3.0d0*height*  &
        cnoidal_ce(mod)/cnoidal_ck(mod)-16.0d0*depth**3*mod**2*cnoidal_ck(mod)**2/  &
        3.0d0/Grav/height/period**2
     if(abs(xa).le.1.0e-8.or.n.gt.1000) goto 50
  
     xb=depth-height+3.0d0*height/2.0d0/mod/(1.0d0-mod)/cnoidal_ck(mod)**2*  &
        ((1.0d0-mod)*cnoidal_ck(mod)**2+cnoidal_ce(mod)**2-2.0d0*(1.0d0-mod)* &
        cnoidal_ck(mod)*cnoidal_ce(mod))-16.0d0*depth**3*mod*cnoidal_ck(mod)/3.0d0/Grav/ &
        (1.0d0-mod)/height/period**2*((1.0d0-mod)*cnoidal_ck(mod)+cnoidal_ce(mod))                          

     mod=mod-xa/xb
     goto 40
 50  continue

     ! sobet el at (1987, J. Waterway)                        
     l=4.0*cnoidal_ck(mod)*depth*sqrt(mod*depth/height/3.0)
     ! mei (1983) or simply c=L/T 
     xtemp=-mod+2.0-3.0*cnoidal_ce(mod)/cnoidal_ck(mod)
     c=sqrt(Grav*depth*(1.0+height/depth/mod*xtemp))
     x2=height/mod*(1.0-mod-cnoidal_ce(mod)/cnoidal_ck(mod))
 
     return
     end subroutine cnoidal


     function cnoidal_cn(u,mod)
!----------------------------------------------------------------
!    Cnoidal Function
!----------------------------------------------------------------
     use global, only: SP
     implicit none

     real(SP), intent(in) :: mod,u
     real(SP) ::  mod1,a0,a1,b0,b1,c0,c1,c(1000),a(1000),y(1000),cnoidal_cn
     integer :: n,i

     mod1 = 1.0-mod
     a0 = 1.0
     b0 = sqrt(mod1)
     c0 = sqrt(mod)
     n = 1
     a(n) = a0
     c(n) = c0

 15  if(abs(c0)<1.0e-15.or.n>1000) then
       goto 30
     else
       n = n+1
       a1 = (a0+b0)/2.0
       b1 = sqrt(a0*b0)
       c1 = (a0-b0)/2.0

       a0 = a1
       b0 = b1
       c0 = c1

       a(n) = a0
       c(n) = c0

       goto 15
     endif

 30  y(n) = 2.0**(n-1)*a(n)*u

     do i = n-1,1,-1
       y(i) = (y(i+1)+asin(c(i+1)/a(i+1)*sin(y(i+1))))/2.0
     enddo

     cnoidal_cn = cos(y(1))

     return 
     end function cnoidal_cn

     function cnoidal_ck(mod)
!----------------------------------------------------------------
!    Cnoidal Function
!----------------------------------------------------------------
     use global, only: SP
     implicit none

     real(SP), intent(in) :: mod
     real(SP) :: mod1,a0,a1,b0,b1,c0,c1,cnoidal_ck
     integer :: n

     mod1 = 1.0-mod
     cnoidal_ck = 0.0
     a0 = 1.0
     b0 = dsqrt(mod1)
     c0 = dsqrt(mod)
     n = 1

 15  if(abs(c0)<1.0e-15.or.n>1000) then
       goto 30
     else
       n = n+1
       a1 = (a0+b0)/2.0
       b1 = sqrt(a0*b0)
       c1 = (a0-b0)/2.0

       a0 = a1
       b0 = b1
       c0 = c1

       goto 15
     endif

 30  cnoidal_ck = 3.1415926535897932384626/2.0/a0

     return
     end function cnoidal_ck

     function cnoidal_ce(mod)
!----------------------------------------------------------------
!    Cnoidal Function
!----------------------------------------------------------------
     use global, only: SP
     implicit none

     real(SP), intent(in) :: mod
     real(SP) :: sum,mod1,a0,a1,b0,b1,c0,c1,ck,c(1000),cnoidal_ce
     integer :: n,k1

     cnoidal_ce = 0.0
     mod1 = 1.0-mod
     a0 = 1.0
     b0 = sqrt(mod1)
     c0 = sqrt(mod)
     n = 1
     c(n) = c0

 15  if(abs(c0)<1.0e-15.or.n>1000) then
       goto 30
     else
       n = n+1
       a1 = (a0+b0)/2.0
       b1 = sqrt(a0*b0)
       c1 = (a0-b0)/2.0

       a0 = a1
       b0 = b1
       c0 = c1
       c(n) = c0

       goto 15
     endif

 30  ck = 3.1415926535897932384626/2.0/a0

     sum = 0.0
     do k1 = 1,n
       sum = sum+2.0**(k1-2)*c(k1)**2
     enddo

     cnoidal_ce = ck*(1.0-sum)

     return
     end function cnoidal_ce


     subroutine random_wave_left_boundary
!-----------------------------------------------------    
!    This subroutine is used to specify left boundary 
!    Called by 
!       update_wave_bc  
!    Last update: 04/11/2011, Gangfeng Ma  
!    fyshi changed the option for parallel 12/16/2011
!-----------------------------------------------------  
# if defined (PARALLEL)
     use global, only: SP,pi,Zero,small,Ibeg,Grav,TIME,Nloc,Kloc,NumFreq,NumDir,Freq,Dire,Wave_Spc2d, & 
                       Ein_X0,Din_X0,Uin_X0,Vin_X0,Win_X0,Hfx,Jbeg,Jend,Kbeg,Kend,sigc,dy, &
                       MaxNumFreq,MaxNumDir,Nglob,myid,Dep_Wave,Random_Phs
# else
     use global, only: SP,pi,Zero,small,Ibeg,Grav,TIME,Nloc,Kloc,NumFreq,NumDir,Freq,Dire,Wave_Spc2d, & 
                       Ein_X0,Din_X0,Uin_X0,Vin_X0,Win_X0,Hfx,Jbeg,Jend,Kbeg,Kend,sigc,dy, &
                       MaxNumFreq,MaxNumDir,Nglob,Dep_Wave,Random_Phs
# endif
     implicit none
     integer :: nfreq,ndir,Iter,nk,i,j,k
     real(SP) :: Per_Wave,Segma,Celerity,Wave_Length,Wave_Number,Fk,Fkdif,Angle,tmp1,tmp2,tmp3,tmp4,Fact, &
                 Wtheta(MaxNumDir,MaxNumFreq),Wh(MaxNumDir,MaxNumFreq),Stokes_Drift(MaxNumDir,MaxNumFreq), &
                 Wnum(MaxNumDir,MaxNumFreq),Phs_Lag,dfreq,ddir,Zlev

     ! find the right angle for periodic bc
     do nfreq = 1,NumFreq
     do ndir = 1,NumDir
       Per_Wave = 1.0/Freq(nfreq)
       Segma = 2.0*pi/Per_Wave
       Celerity = sqrt(Grav*Dep_Wave)
       Wave_Length = Celerity*Per_Wave
       Wave_Number = 2.0*pi/Wave_Length

       Iter = 0
 75    Fk = Grav*Wave_Number*tanh(Wave_Number*Dep_Wave)-Segma**2
       if(abs(Fk)<=1.0e-8.or.Iter>1000) goto 85
         Fkdif = Grav*Wave_Number*Dep_Wave*(1.0-tanh(Wave_Number*Dep_Wave)**2)+  &                
            Grav*tanh(Wave_Number*Dep_Wave)
         Wave_Number = Wave_Number-Fk/Fkdif
         Iter = Iter+1
         goto 75
 85    continue
       Wave_Length = 2.0*pi/Wave_Number

       Angle = Dire(ndir)*pi/180.
       goto 100
       if(Angle>zero) then
         tmp3 = zero
         tmp1 = Wave_Number
         nk = 0
         do while (tmp3<Angle)
           nk = nk+1
           tmp2 = nk*2.0*pi/(Nglob*dy)
           if(tmp2>=tmp1) then
             tmp3 = 0.5*pi-small
           else
             tmp3 = asin(tmp2/tmp1)
           endif
         enddo

         ! judge between nk-1 and nk which is closer
         tmp4 = asin((nk-1)*2.0*pi/(Nglob*dy)/tmp1)
         if(abs(tmp4-Angle)<abs(Angle-tmp3)) then
           Angle = tmp4
         else
           Angle = tmp3
         endif
       else
         tmp3 = zero
         tmp1 = Wave_Number
         nk = 0
         do while (tmp3>Angle)
           nk = nk+1
           tmp2 = nk*2.0*pi/(Nglob*dy)
           if(tmp2>=tmp1) then
             tmp3 = -0.5*pi+small
           else
             tmp3 = -asin(tmp2/tmp1)
           endif
         enddo

         ! judge between nk-1 and nk which is closer 
         tmp4= asin((nk-1)*2.0*pi/(Nglob*dy)/tmp1)
         if(abs(tmp4-Angle)<abs(Angle-tmp3)) then
           Angle = tmp4
         else
           Angle = tmp3
         endif
       endif
 100   continue
       Wtheta(ndir,nfreq) = Angle
 
       if(nfreq==1) then
         dfreq = Freq(2)-Freq(1)
       elseif(nfreq==NumFreq) then
         dfreq = Freq(NumFreq)-Freq(NumFreq-1)
       else
         dfreq = 0.5*(Freq(nfreq+1)-Freq(nfreq-1))
       endif
       dfreq = abs(dfreq)

       if(ndir==1) then
         ddir = Dire(2)-Dire(1)
       elseif(ndir==NumDir) then
         ddir = Dire(NumDir)-Dire(NumDir-1)
       else
         ddir = 0.5*(Dire(ndir+1)-Dire(ndir-1))
       endif
       ddir = abs(ddir)
       
       ! save wave number and wave height
       Wnum(ndir,nfreq) = Wave_Number
       Wh(ndir,nfreq) = 2.0*sqrt(2.0*Wave_Spc2d(ndir,nfreq)*dfreq*ddir)

       ! Stokes Drift
       Stokes_Drift(ndir,nfreq) = Grav*Wh(ndir,nfreq)**2/(Wave_Length/Per_Wave)/Dep_Wave/8.0
!       if(myid.eq.0)write(12,'(8E20.10)') Freq(nfreq),Dire(ndir),dfreq,ddir,Wtheta(ndir,nfreq)*180./pi,  &
!               Wnum(ndir,nfreq),Wh(ndir,nfreq),Stokes_Drift(ndir,nfreq)
     enddo
     enddo

     Fact = 1.0
     do j = 1,Nloc
       Ein_X0(j) = 0.0
       do nfreq = 1,NumFreq
       do ndir = 1,NumDir
         Segma = 2.0*pi*Freq(nfreq)
         Phs_Lag = ((j-Jbeg)*dy)*sin(Wtheta(ndir,nfreq))*Wnum(ndir,nfreq)+Random_Phs(ndir,nfreq)
         Ein_X0(j) = Ein_X0(j)+0.5*Wh(ndir,nfreq)*cos(pi/2-Segma*TIME+Phs_Lag)
       enddo
       enddo
       Ein_X0(j)=Ein_X0(j)*tanh(TIME/20.0)
       Din_X0(j) = Ein_X0(j)+Hfx(Ibeg,j)
       !if(myid.eq.0) write(12,*) TIME,Ein_X0(Jbeg+Nglob/2)
     enddo

     do k = Kbeg,Kend
     do j = Jbeg,Jend
       Zlev = sigc(k)*Din_X0(j)
       Uin_X0(j,k) = 0.0
       Win_X0(j,k) = 0.0
       Vin_X0(j,k) = 0.0
       do nfreq = 1,NumFreq
       do ndir = 1,NumDir
         Segma = 2.0*pi*Freq(nfreq)
         Phs_Lag = ((j-Jbeg)*dy)*sin(Wtheta(ndir,nfreq))*Wnum(ndir,nfreq)+Random_Phs(ndir,nfreq)
         Uin_X0(j,k) = Uin_X0(j,k)+(0.5*Wh(ndir,nfreq)*Segma*cosh(Wnum(ndir,nfreq)*Zlev)/  &
             sinh(Wnum(ndir,nfreq)*Dep_Wave)*cos(pi/2-Segma*TIME+Phs_Lag)-Stokes_Drift(ndir,nfreq)*Fact)*  &
             cos(Wtheta(ndir,nfreq)) 
         Win_X0(j,k) = Win_X0(j,k)+0.5*Wh(ndir,nfreq)*Segma*sinh(Wnum(ndir,nfreq)*Zlev)/  &
             sinh(Wnum(ndir,nfreq)*Dep_Wave)*sin(pi/2-Segma*TIME+Phs_Lag)
         Vin_X0(j,k) = Vin_X0(j,k)+(0.5*Wh(ndir,nfreq)*Segma*cosh(Wnum(ndir,nfreq)*Zlev)/  & 
             sinh(Wnum(ndir,nfreq)*Dep_Wave)*cos(pi/2-Segma*TIME+Phs_Lag)-Stokes_Drift(ndir,nfreq)*Fact)*  & 
             sin(Wtheta(ndir,nfreq))
       enddo
       enddo
       Uin_X0(j,k)=Uin_X0(j,k)*tanh(TIME/20.0)
       Vin_X0(j,k)=Vin_X0(j,k)*tanh(TIME/20.0)
       Win_X0(j,k)=Win_X0(j,k)*tanh(TIME/20.0)
     enddo
     enddo

     end subroutine random_wave_left_boundary


     subroutine linear_wave_left_boundary
!-----------------------------------------------------
!    This subroutine is used to specify left boundary 
!    Called by
!       update_wave_bc 
!    Last update: 05/02/2011, Gangfeng Ma
!-----------------------------------------------------
     use global, only: SP,pi,Zero,Ibeg,Grav,TIME,Nloc,Kloc,Amp_Wave,Per_Wave,Dep_Wave,Theta_Wave, &
                       Ein_X0,Din_X0,Uin_X0,Vin_X0,Win_X0,Hfx,Jbeg,Jend,Kbeg,Kend,sigc,dy
     implicit none
     integer  :: j,k,Iter
     real(SP) :: Segma,Celerity,Wave_Length,Wave_Number,Fk,Fkdif,Zlev,Stokes_Drift,Fact,Phs_Lag(Nloc)

     ! Find wave number for linear wave (Newton-Ralphson Method)
     Segma = 2.0*pi/Per_Wave
     Celerity = sqrt(Grav*Dep_Wave)
     Wave_Length = Celerity*Per_Wave
     Wave_Number = 2.0*pi/Wave_Length
     
     Iter = 0
55   Fk = Grav*Wave_Number*tanh(Wave_Number*Dep_Wave)-Segma**2
     if(abs(Fk)<=1.0e-8.or.Iter>1000) goto 65
     Fkdif = Grav*Wave_Number*Dep_Wave*(1.0-tanh(Wave_Number*Dep_Wave)**2)+  &
        Grav*tanh(Wave_Number*Dep_Wave) 
     Wave_Number = Wave_Number-Fk/Fkdif
     Iter = Iter+1
     goto 55
65   continue
     Wave_Length = 2.0*pi/Wave_Number

     Stokes_Drift = Grav*Amp_Wave**2/(Wave_Length/Per_Wave)/Dep_Wave/8.0
     Fact = 1.0     

     do j = 1,Nloc
       Phs_Lag(j) = ((j-Jbeg)*dy)*sin(Theta_Wave*pi/180.)*Wave_Number
       Ein_X0(j) = 0.5*Amp_Wave*cos(pi/2-Segma*TIME+Phs_Lag(j))
       Din_X0(j) = Ein_X0(j)+Hfx(Ibeg,j)     
     enddo

     do k = Kbeg,Kend
     do j = Jbeg,Jend
       Zlev = sigc(k)*Din_X0(j)
       Uin_X0(j,k) = (0.5*Amp_Wave*Segma*cosh(Wave_Number*Zlev)/  &
           sinh(Wave_Number*Dep_Wave)*cos(pi/2-Segma*TIME+Phs_Lag(j))-Stokes_Drift*Fact)*cos(Theta_Wave*pi/180.)
       Win_X0(j,k) = 0.5*Amp_Wave*Segma*sinh(Wave_Number*Zlev)/  &
           sinh(Wave_Number*Dep_Wave)*sin(pi/2-Segma*TIME+Phs_Lag(j))
       Vin_X0(j,k) = (0.5*Amp_Wave*Segma*cosh(Wave_Number*Zlev)/  &
           sinh(Wave_Number*Dep_Wave)*cos(pi/2-Segma*TIME+Phs_Lag(j))-Stokes_Drift*Fact)*sin(Theta_Wave*pi/180.)
     enddo
     enddo         

     end subroutine linear_wave_left_boundary


     subroutine solitary_wave_left_boundary
!------------------------------------------------------
!    This subroutine is used to specify left boundary
!    Called by
!       update_bc
!    Last update: 22/12/2010, Gangfeng Ma 
!-----------------------------------------------------
     use global, only: SP,Zero,pi,Ibeg,Grav,TIME,Nloc,Kloc,Amp_Wave,Dep_Wave, &
                       Ein_X0,Din_X0,Uin_X0,Vin_X0,Win_X0,Hfx,Kbeg,Kend,sigc,BC_X0
     implicit none
     integer  :: j,k
     real(SP) :: Celerity,Atmp,Xstart,C2,D1,D2,D3,Zlev,Xc

     Uin_X0 = Zero
     Vin_X0 = Zero
     Win_X0 = Zero
     
     Celerity = sqrt(Grav*Dep_Wave*(1.0+Amp_Wave/Dep_Wave))
     do j = 1,Nloc
       Atmp = sqrt(0.75*Amp_Wave/Dep_Wave**3)
       Xstart = 4.0*Dep_Wave/sqrt(Amp_Wave/Dep_Wave)
       Xc = Xstart-Celerity*TIME
       Ein_X0(j) = Amp_Wave/cosh(Atmp*Xc)**2
       Din_X0(j) = Ein_X0(j)+Hfx(Ibeg,j)
 
       C2 = sqrt(Grav*Dep_Wave)
       D2 = 2.0*Amp_Wave*Atmp**2*(2.0*cosh(Atmp*Xc)**2-3.)/cosh(Atmp*Xc)**4
       D1 = -2.0*Amp_Wave*sinh(Atmp*Xc)*Atmp/cosh(Atmp*Xc)**3
       D3 = -8.0*Amp_Wave*sinh(Atmp*Xc)*Atmp**3*(cosh(Atmp*Xc)**2-3.)/cosh(Atmp*Xc)**5
       do k = Kbeg,Kend
         Zlev = sigc(k)*Din_X0(j)
         Uin_X0(j,k) = C2*Ein_X0(j)/Dep_Wave*(1.0-1.0/4.0*Ein_X0(j)/Dep_Wave+Dep_Wave/3.0*(Dep_Wave/Ein_X0(j))*  &
                           (1.0-3.0/2.0*Zlev**2/Dep_Wave**2)*D2)
         Win_X0(j,k) = -C2*Zlev/Dep_Wave*((1.0-1.0/2.0*Ein_X0(j)/Dep_Wave)*D1+1.0/3.0*Dep_Wave**2*  &
                           (1.0-1.0/2.0*Zlev**2/Dep_Wave**2)*D3)
         Vin_X0(j,k) = 0.0
       enddo
     enddo

     end subroutine solitary_wave_left_boundary


subroutine update_bathymetry
!------------------------------------------------------
!    This subroutine is used to update time-varying bathymetry 
!    Called by 
!       main 
!    update: 12/05/2011, Gangfeng Ma 
!    update: 04/13/2012, Fengyan Shi
!    update: 07/08/2013, Dmitry Nicolsky
!    DJN makes a viscous landslide application
!-----------------------------------------------------
     use global
     implicit none
     integer :: i,j,m,n,iter, RKIT

              
     ! save old bathymetry
     Ho = Hc

     ! update D and Eta
     !call phi_2D_coll(D)
     !call phi_2D_coll(Hc)
     !call phi_int_exch(Mask9)

     SldWZ = Eta !D-Hc


     SldZs=SldZ0;   SldDZs=SldDZ0
     SldNs=SldN0;   SldMs=SldM0
     DO RKIT=1,3
       IF (RKIT.EQ.1) THEN
	 SLIDE_RKC1=1.0d0; SLIDE_RKC2=0.0d0; SLIDE_RKDT=dt
       ELSEIF (RKIT.EQ.2) THEN
	 SLIDE_RKC1=3.0d0/4.0d0; SLIDE_RKC2=1.0d0/4.0d0; SLIDE_RKDT=1.0d0/4.0d0*dt
       ELSEIF (RKIT.EQ.3) THEN
	 SLIDE_RKC1=1.0d0/3.0d0; SLIDE_RKC2=2.0d0/3.0d0; SLIDE_RKDT=2.0d0/3.0d0*dt
       ENDIF

       call nlmassF90();      call nlmmntF90()

       SldZs=SldZ1;   SldDZs=SldDZ1
       SldNs=SldN1;   SldMs=SldM1
     ENDDO
     
     SldZ0=SldZ1;   SldDZ0=SldDZ1
     SldN0=SldN1;   SldM0=SldM1


     if(trim(adjustl(DEPTH_TYPE))=='CELL_GRID') then
       write(*,*) 'DEPTH_TYPE has to be cell_center'
     elseif(trim(adjustl(DEPTH_TYPE))=='CELL_CENTER') then

       ! interpolate into grid center
       do j=Jbeg,Jend
       do i=Ibeg,Iend
         Hc(i,j) = DepC0(i,j)-SldDZ1(i,j)
       enddo
       enddo
     endif ! end cell center

     ! ghost cells
     call phi_2D_coll(Hc)

     ! reconstruct depth at x-y faces
     do j = 1,Nloc
     do i = 2,Mloc
       Hfx(i,j) = 0.5*(Hc(i,j)+Hc(i+1,j))
     enddo
     Hfx(1,j) = Hfx(2,j)
     Hfx(Mloc1,j) = Hfx(Mloc,j)
     enddo

     do i = 1,Mloc
     do j = 2,Nloc
       Hfy(i,j) = 0.5*(Hc(i,j)+Hc(i,j+1))
     enddo
     Hfy(i,1) = Hfy(i,2)
     Hfy(i,Nloc1) = Hfy(i,Nloc)
     enddo

     ! derivatives of water depth at cell center
     do j = 1,Nloc
     do i = 1,Mloc
       DelxH(i,j) = (Hfx(i+1,j)-Hfx(i,j))/dx
       DelyH(i,j) = (Hfy(i,j+1)-Hfy(i,j))/dy
     enddo
     enddo

     ! time derivative of water depth
     DeltHo = DeltH

     DeltH = zero
     do j = 1,Nloc
     do i = 1,Mloc
       DeltH(i,j) = (Hc(i,j)-Ho(i,j))/dt
     enddo
     enddo

     ! second-order time derivative
     if(RUN_STEP>2) Delt2H = (DeltH-DeltHo)/dt

     end subroutine update_bathymetry




     subroutine read_bathymetry
!------------------------------------------------------ 
!    This subroutine is used to read bathymetry                                                
!    Called by  
!       main    
!    update: 21/12/2010, Gangfeng Ma    
!    update: 04/13/2012, Fengyan Shi
!    update: 07/08/2013, Dmitry Nicolsky
!    DJN makes a viscous landslide application
!-----------------------------------------------------
     use global
     implicit none
     integer :: i,j,m,n,iter,iglob,jglob
     integer :: Maskp(Mglob+1,Nglob+1)
     real(SP), dimension(Mglob+1,Nglob+1) :: HG
     real(SP), dimension(Mglob,Nglob) :: HCG

     real(SP), dimension(Mglob,Nglob) :: SldDZtmp


     if(trim(adjustl(DEPTH_TYPE))=='CELL_GRID') then
       write(*,*) 'DEPTH_TYPE has to be cell_center'
     endif
  
# if defined (DEBUG)
print*,'cell_grid OK'
# endif

     ! read bathymetry at cell center

     open(5,file=TRIM(Depth_File),status='old')
     do j = 1,Nglob
       read(5,*) (HCG(i,j),i=1,Mglob)
     enddo

# if defined (DEBUG)
print*,'open read bathymetry file OK'
# endif

# if defined (PARALLEL)
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       iglob = npx*(Mloc-2*Nghost)+i-Nghost
       jglob = npy*(Nloc-2*Nghost)+j-Nghost
       Hc(i,j) = HCG(iglob,jglob)
     enddo
     enddo
# else
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       iglob = i-Nghost
       jglob = j-Nghost
       Hc(i,j) = HCG(iglob,jglob)
     enddo
     enddo
# endif

# if defined (DEBUG)
print*,'assign HCG to Hc OK'
# endif

     call phi_2D_coll(Hc)
     DepC0 = Hc   ! save base bathymetry for field case

# if defined (DEBUG)
print*,'DepC0 OK'
# endif

     open(5001,file=TRIM(Slide_File),status='old')
     do j = 1,Nglob
       read(5001,*) (SldDZtmp(i,j),i=1,Mglob)
     enddo
     close(5001)

# if defined (DEBUG)
print*,'read slide OK'
# endif

     do i = Ibeg,Iend
     do j = Jbeg,Jend
# if defined (PARALLEL)
       iglob = npx*(Mloc-2*Nghost)+i-Nghost
       jglob = npy*(Nloc-2*Nghost)+j-Nghost
# else
       iglob = i-Nghost
       jglob = j-Nghost
# endif

       SldDZ0(i,j)= SldDZtmp(iglob,jglob)

       if(SldDZ0(i,j)<0.1) SldDZ0(i,j)=0.0

       SldZ0(i,j)=SldDZ0(i,j)-Hc(i,j)

       SldDZs(i,j)=SldDZ0(i,j)
       SldZs(i,j)=SldZ0(i,j)
       SldDZ1(i,j)=SldDZ0(i,j)
       SldZ1(i,j)=SldZ0(i,j)
     enddo
     enddo

# if defined (DEBUG)
print*,'SldZ OK'
# endif

     SldN0=Zero
     SldNs=Zero
     SldN1=Zero
     SldM0=Zero
     SldMs=Zero
     SldM1=Zero

! interpolate into grid center
     do j=Jbeg,Jend
     do i=Ibeg,Iend
       Hc(i,j) = Hc(i,j)-SldDZ0(i,j) 
     enddo
     enddo

# if defined (DEBUG)
print*,'interpolate into center OK'
# endif

     ! collect data into ghost cells
     call phi_2D_coll(Hc)
     call phi_2D_coll(SldDZ0) 
     call phi_2D_coll(SldDZs) 
     call phi_2D_coll(SldDZ1) 
     call phi_2D_coll(SldZ0)   
     call phi_2D_coll(SldZs)   
     call phi_2D_coll(SldZ1)  


# if defined (DEBUG)
print*,'fill ghost cells OK'
# endif

     ! find pernament dry cells
     Mask_Struct = 1
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       if(Hc(i,j)<-10.0) then
         Mask_Struct(i,j) = 0
       endif
     enddo
     enddo

     ! reconstruct depth at x-y faces

! gfma fixed the bug 12/20/2011
     do j = 1,Nloc
     do i = 2,Mloc
       Hfx(i,j) = 0.5*(Hc(i-1,j)+Hc(i,j))
     enddo
     Hfx(1,j) = Hc(1,j)
     Hfx(Mloc1,j) = Hc(Mloc,j)
     enddo

     do i = 1,Mloc
     do j = 2,Nloc
       Hfy(i,j) = 0.5*(Hc(i,j-1)+Hc(i,j))
     enddo
     Hfy(i,1) = Hc(i,1)
     Hfy(i,Nloc1) = Hc(i,Nloc)
     enddo


     ! derivatives of water depth at cell center
     do j = 1,Nloc
     do i = 1,Mloc
       DelxH(i,j) = (Hfx(i+1,j)-Hfx(i,j))/dx
       DelyH(i,j) = (Hfy(i,j+1)-Hfy(i,j))/dy
     enddo
     enddo

# if defined (DEBUG)
print*,'out of read_bathymetry OK'
# endif

     end subroutine read_bathymetry

    
     subroutine generate_grid
!------------------------------------------------------
!    This subroutine is used to generate grids
!    Called by
!       main
!    Last update: 20/12/2010, Gangfeng Ma
!-----------------------------------------------------
     use global
     implicit none
     integer :: i,j,k

     ! horizontal grid
# if defined (PARALLEL)
     x(Ibeg) = npx*(Mloc-2*Nghost)*dx
# else
     x(Ibeg) = Zero
# endif
     do i = Ibeg+1,Mloc1
       x(i) = x(i-1)+dx
       xc(i-1) = x(i-1)+0.5*dx
     enddo
     do i = Ibeg-1,Ibeg-Nghost,-1
       x(i) = x(i+1)-dx
       xc(i) = x(i+1)-0.5*dx
     enddo

# if defined (PARALLEL)
     y(Jbeg) = npy*(Nloc-2*Nghost)*dy
# else
     y(Jbeg) = Zero
# endif
     do j = Jbeg+1,Nloc1
       y(j) = y(j-1)+dy
       yc(j-1) = y(j-1)+0.5*dy
     enddo
     do j = Jbeg-1,Jbeg-Nghost,-1
       y(j) = y(j+1)-dy
       yc(j) = y(j+1)-0.5*dy
     enddo

     ! vertical grid
     if(Ivgrd==1) then
       do k = 1,Kloc
         dsig(k) = 1.0/float(Kglob)
       enddo
     elseif(Ivgrd==2) then
       write(*,*) 'under development'
     endif

     sig(Kbeg) = Zero
     do k = Kbeg+1,Kloc1
       sig(k) = sig(k-1)+dsig(k-1)
       sigc(k-1) = sig(k-1)+0.5*dsig(k-1)
     enddo
     do k = Kbeg-1,1,-1
       sig(k) = sig(k+1)-dsig(k)
       sigc(k) = sig(k+1)-0.5*dsig(k)
     enddo

     end subroutine generate_grid

     
     subroutine initial
!---------------------------------------------------  
!    This subroutine is used to initialize model run 
!    Called by                                       
!       main 
!    Last update: 21/12/2010, Gangfeng Ma  
!---------------------------------------------------
     use global
     implicit none
     integer  :: i,j,k,n,m,nmax
     real(SP) :: xsol(80),zsol(80),zmax,xmax,xterp,zterp,tmp,zc(Kloc), &
                 utmp1,wtmp1,utmp2,wtmp2,xk(321,16),zk(321,16),  &
                 uk(321,16),wk(321,16)

     ! simulation time
     TIME = Zero
     RUN_STEP = 0
     dt = dt_ini     
     Screen_Count = Zero
     Plot_Count = Zero
     Plot_Count_Stat = Zero
     Icount = 0
     
     ! working arrays
     D = Zero
     U = Zero
     V = Zero
     W = Zero
     P = Zero
     Omega = Zero
     DU = Zero
     DV = Zero
     DW = Zero
     D0 = Zero
     Eta0 = Zero
     DU0 = Zero
     DV0 = Zero
     DW0 = Zero
     Uf = Zero
     Vf = Zero
     Wf = Zero
     Rho = Rho0
     
     ! source terms
     SourceC = Zero
     SourceX = Zero
     SourceY = Zero

     ! fluxes
     DxL = Zero
     DxR = Zero
     DyL = Zero
     DyR = Zero
     UxL = Zero
     UxR = Zero
     UyL = Zero
     UyR = Zero
     UzL = Zero
     UzR = Zero
     VxL = Zero
     VxR = Zero
     VyL = Zero
     VyR = Zero
     VzL = Zero
     VzR = Zero
     WxL = Zero
     WxR = Zero
     WyL = Zero
     WyR = Zero
     WzL = Zero
     WzR = Zero
     DUxL = Zero
     DUxR = Zero
     DUyL = Zero
     DUyR = Zero
     DVxL = Zero
     DVxR = Zero
     DVyL = Zero
     DVyR = Zero
     DWxL = Zero
     DWxR = Zero
     DWyL = Zero
     DWyR = Zero
     OzL = Zero
     OzR = Zero
     SxL = Zero
     SxR = Zero
     SyL = Zero
     SyR = Zero
     Ex = Zero
     Ey = Zero
     Fx = Zero
     Fy = Zero
     Fz = Zero
     Gx = Zero
     Gy = Zero
     Gz = Zero
     Hx = Zero
     Hy = Zero
     Hz = Zero
     EtaxL = Zero
     EtaxR = Zero
     EtayL = Zero
     EtayR = Zero
     DelxEta = Zero
     DelyEta = Zero
     DeltH = Zero
     DeltHo = Zero
     Delt2H = Zero
     DelxD = Zero
     DelyD = Zero
     DelxU = Zero
     DelyU = Zero
     DelzU = Zero
     DelxV = Zero
     DelyV = Zero
     DelzV = Zero
     DelxW = Zero
     DelyW = Zero
     DelzW = Zero
     DelxDU = Zero
     DelyDU = Zero
     DelxDV = Zero
     DelyDV = Zero
     DelzO = Zero
     Sponge = One
     Cmu = Visc
     CmuHt = Visc
     CmuVt = Visc
     Diffxx = Zero
     Diffxy = Zero
     Diffxz = Zero
     Diffyx = Zero
     Diffyy = Zero
     Diffyz = Zero
     Diffzx = Zero
     Diffzy = Zero
     Diffzz = Zero

     Uin_X0 = Zero
     Vin_X0 = Zero
     Win_X0 = Zero
     Ein_X0 = Zero
     Din_X0 = Zero
     Uin_Xn = Zero
     Vin_Xn = Zero
     Win_Xn = Zero
     Ein_Xn = Zero
     Din_Xn = Zero   

     EleMax = Zero
     EleMin = Large
     FrMax = ZERO
     FrMax_total = ZERO
     TimeMaxFr = ZERO
     DepMaxFr = ZERO

# if defined (BUBBLE)
     Nbg = Zero
     DNbg = Zero
     DNbg0 = Zero
     Vbg = Zero
# endif

     Tke = Zero
     Eps = Zero
     DTke = Zero
     DEps = Zero
     DTke0 = Zero
     DEps0 = Zero

     ! pressure boundary
     Bc_Prs = Zero

     ! initial surface elevation (user-specified)
     Eta = Zero

     if(WaveMaker(1:3)=='INI') then

       goto 100

       ! solitary wave from Tanaka solution
       open(21,file='soliton.dat')
       do n = 1,80
         read(21,*) i,xsol(n),zsol(n),tmp,tmp
       enddo
       close(21)

       ! find the peak location                                                   
       zmax = -1.0e+10
       do n = 1,80
         if(zsol(n)>zmax) then
           zmax = zsol(n)
           xmax = xsol(n)
         endif
       enddo

       ! move the peak to x = 3.0m                                                
       do n = 1,80
         xsol(n) = xsol(n)+8.0-xmax
       enddo

       ! interpolate into computational grid
       do j = Jbeg,Jend
       do i = Ibeg,Iend
         if(xc(i)>xsol(80).and.xc(i)<xsol(1)) cycle
         do n = 2,80
           if(xc(i)>=xsol(n-1).and.xc(i)<xsol(n)) then
             xterp = (xc(i)-xsol(n-1))/(xsol(n)-xsol(n-1))
             Eta(i,j) = (1.0-xterp)*zsol(n-1)+xterp*zsol(n)
           endif
         enddo
       enddo
       enddo

       open(22,file='plotuv.dat')
       do n = 1,16
       do m = 1,321
         read(22,*) k,xk(m,n),zk(m,n),tmp,uk(m,n),wk(m,n)
         xk(m,n) = xk(m,n)+8.0
       enddo
       enddo
       close(22)
       
       do j = Jbeg,Jend
       do k = Kbeg,Kend
       do i = Ibeg,Iend
         if(xc(i)>xsol(80).and.xc(i)<xsol(1)) cycle
         zc(k) = (1.0+Eta(i,j))*sigc(k)-1.0

         do n = 2,16
         do m = 2,321
           if(xc(i)>=xk(m-1,n).and.xc(i)<xk(m,n).and.  &
                    zc(k)>=zk(m,n-1).and.zc(k)<zk(m,n)) then
             xterp = (xc(i)-xk(m-1,n))/(xk(m,n)-xk(m-1,n))
             zterp = (zc(k)-zk(m,n-1))/(zk(m,n)-zk(m,n-1)) 
             utmp1 = (1.0-xterp)*uk(m-1,n-1)+xterp*uk(m,n-1)
             wtmp1 = (1.0-xterp)*wk(m-1,n-1)+xterp*wk(m,n-1)
             utmp2 = (1.0-xterp)*uk(m-1,n)+xterp*uk(m,n)
             wtmp2 = (1.0-xterp)*wk(m-1,n)+xterp*wk(m,n)
             
             U(i,j,k) = (1.0-zterp)*utmp1+zterp*utmp2
             W(i,j,k) = (1.0-zterp)*wtmp1+zterp*wtmp2
           endif
         enddo
         enddo
       enddo
       enddo
       enddo
 100   continue
     endif    

     ! collect data into ghost cells
     call phi_2D_coll(Eta)
     Eta0 = Eta

     call vel_bc
# if defined (PARALLEL)
     call phi_3D_exch(U)
     call phi_3D_exch(V)
     call phi_3D_exch(W)
# endif

     ! wetting-drying mask
     ! Mask: 1 - wet; 0 - dry
     ! Mask_Struct: 0 - permanent dry point
     ! Mask9: mask for itself and 8 elements around
     Mask = 1
     do j = 1,Nloc
     do i = 1,Mloc
       if((Eta(i,j)+Hc(i,j))<=MinDep) then
         Mask(i,j) = 0
         Eta(i,j) = MinDep-Hc(i,j)
       else
         Mask(i,j) = 1
       endif
     enddo
     enddo
     Mask = Mask*Mask_Struct

# if defined (PARALLEL)
     ! collect mask into ghost cells
     call phi_int_exch(Mask)
# endif

     do j = Jbeg,Jend
     do i = Ibeg,Iend
      Mask9(i,j) = Mask(i,j)*Mask(i-1,j)*Mask(i+1,j)  &
                *Mask(i+1,j+1)*Mask(i,j+1)*Mask(i-1,j+1) &
                *Mask(i+1,j-1)*Mask(i,j-1)*Mask(i-1,j-1)
     enddo
     enddo


     ! total water depth and flux
     D = max(Hc+Eta, MinDep)

     do k = 1,Kloc
     do j = 1,Nloc
     do i = 1,Mloc
       DU(i,j,k) = D(i,j)*U(i,j,k)*Mask(i,j)
       DV(i,j,k) = D(i,j)*V(i,j,k)*Mask(i,j)
       DW(i,j,k) = D(i,j)*W(i,j,k)*Mask(i,j)
     enddo
     enddo
     enddo

     ! initial seeding values for turbulence
     Tke_min = 0.5*(1.4e-3)**2
     Eps_min = 0.09*Tke_min**2/(0.1*Visc)
     do k = 1,Kloc
     do j = 1,Nloc
     do i = 1,Mloc
       Tke(i,j,k) = Tke_min
       Eps(i,j,k) = Eps_min
       DTke(i,j,k) = D(i,j)*Tke(i,j,k)*Mask(i,j)
       DEps(i,j,k) = D(i,j)*Eps(i,j,k)*Mask(i,j)
     enddo
     enddo
     enddo

     ! SSP Runge-Kutta method parameters
     if(TIME_ORDER(1:3)=='THI') then
       It_Order = 3
       ALPHA(1) = 0.0
       ALPHA(2) = 3.0/4.0
       ALPHA(3) = 1.0/3.0
       BETA(1) = 1.0
       BETA(2) = 1.0/4.0
       BETA(3) = 2.0/3.0
     elseif(TIME_ORDER(1:3)=='SEC') then
       It_Order = 2
       ALPHA(1) = 0.0
       ALPHA(2) = 1.0/2.0
       BETA(1) = 1.0
       BETA(2) = 1.0/2.0
     else
       It_Order = 1
       ALPHA(1) = 0.0
       BETA(1) = 1.0
     endif

     ! sponge layer
     if(SPONGE_ON) then
       call calculate_sponge
     endif
    
     end subroutine initial

     subroutine allocate_variables 
!--------------------------------------------------- 
!    This subroutine is used to allocate variables
!    Called by                                    
!       main                                                                                  
!    Last update: 23/12/2010, Gangfeng Ma                                                   
!---------------------------------------------------
     use global
     implicit none

     ! one-dimensional vars
     ALLOCATE(x(Mloc1),xc(Mloc),y(Nloc1),yc(Nloc),sig(Kloc1),dsig(Kloc),sigc(Kloc),  &
              Ein_X0(Nloc),Din_X0(Nloc),Ein_Xn(Nloc),Din_Xn(Nloc))

     ! two-dimensional vars
     ALLOCATE(Ho(Mloc,Nloc),H(Mloc1,Nloc1),Hc(Mloc,Nloc),Hfx(Mloc1,Nloc),Hfy(Mloc,Nloc1),D(Mloc,Nloc),Eta0(Mloc,Nloc),  &
              D0(Mloc,Nloc),DeltH(Mloc,Nloc),DelxH(Mloc,Nloc),DelyH(Mloc,Nloc),Eta(Mloc,Nloc),Mask(Mloc,Nloc),  &
              Mask_Struct(Mloc,Nloc),Mask9(Mloc,Nloc),SourceC(Mloc,Nloc),SourceX(Mloc,Nloc),SourceY(Mloc,Nloc),DeltHo(Mloc,Nloc), &
              Uin_X0(Nloc,Kloc),Vin_X0(Nloc,Kloc),Win_X0(Nloc,Kloc),Uin_Xn(Nloc,Kloc),Delt2H(Mloc,Nloc),  &
              Vin_Xn(Nloc,Kloc),Win_Xn(Nloc,Kloc),Bc_Prs(Mloc,Nloc))
     ALLOCATE(DxL(Mloc1,Nloc),DxR(Mloc1,Nloc),DyL(Mloc,Nloc1),DyR(Mloc,Nloc1), &
              EtaxL(Mloc1,Nloc),EtaxR(Mloc1,Nloc),EtayL(Mloc,Nloc1),EtayR(Mloc,Nloc1), &
              DelxEta(Mloc,Nloc),DelyEta(Mloc,Nloc),DelxD(Mloc,Nloc),DelyD(Mloc,Nloc),Sponge(Mloc,Nloc))

! fyshi add tmp4preview
     ALLOCATE(tmp4preview2D(Mloc,Nloc))

!DJN add the viscous landslide
     ALLOCATE(SldM0(Mloc,Nloc), SldN0(Mloc,Nloc), SldZ0(Mloc,Nloc), SldDZ0(Mloc,Nloc), &
              SldM1(Mloc,Nloc), SldN1(Mloc,Nloc), SldZ1(Mloc,Nloc), SldDZ1(Mloc,Nloc), & 
              SldMs(Mloc,Nloc), SldNs(Mloc,Nloc), SldZs(Mloc,Nloc), SldDZs(Mloc,Nloc), & 
              SldWZ(Mloc,Nloc), SldUs(Mloc,Nloc), SldVs(Mloc,Nloc))
!DJN

     ALLOCATE(EleMax(Mloc,Nloc),EleMin(Mloc,Nloc),FrMax(Mloc,Nloc,Kloc),DepMaxFr(Mloc,Nloc), &
              Froude(Mloc,Nloc,Kloc),FrMax_total(Mloc,Nloc,Kloc),TimeMaxFr(Mloc,Nloc), &
              Froude_total(Mloc,Nloc,Kloc))

     ! three-dimensional vars

! fyshi add tmp4preview
     ALLOCATE(tmp4preview3D(Mloc,Nloc,Kloc))

! fyshi add base bathymetry for landslide
     ALLOCATE(Dep0(Mloc1,Nloc1),DepC0(Mloc,Nloc))

     ALLOCATE(U(Mloc,Nloc,Kloc),V(Mloc,Nloc,Kloc),W(Mloc,Nloc,Kloc),Omega(Mloc,Nloc,Kloc), &
              P(Mloc,Nloc,Kloc1),DU(Mloc,Nloc,Kloc),DV(Mloc,Nloc,Kloc),DW(Mloc,Nloc,Kloc),  &
              DU0(Mloc,Nloc,Kloc),DV0(Mloc,Nloc,Kloc),DW0(Mloc,Nloc,Kloc),Uf(Mloc,Nloc,Kloc1), &
              Vf(Mloc,Nloc,Kloc1),Wf(Mloc,Nloc,Kloc1),Cmu(Mloc,Nloc,Kloc), &
              Diffxx(Mloc,Nloc,Kloc),Diffxy(Mloc,Nloc,Kloc),Diffxz(Mloc,Nloc,Kloc), &
              Diffyx(Mloc,Nloc,Kloc),Diffyy(Mloc,Nloc,Kloc),Diffyz(Mloc,Nloc,Kloc),Diffzx(Mloc,Nloc,Kloc), &
              Diffzy(Mloc,Nloc,Kloc),Diffzz(Mloc,Nloc,Kloc),DelxSc(Mloc,Nloc,Kloc),DelySc(Mloc,Nloc,Kloc), &
              CmuHt(Mloc,Nloc,Kloc),CmuVt(Mloc,Nloc,Kloc),Rho(Mloc,Nloc,Kloc),Tke(Mloc,Nloc,Kloc), &
              Eps(Mloc,Nloc,Kloc),DTke(Mloc,Nloc,Kloc),DEps(Mloc,Nloc,Kloc),DTke0(Mloc,Nloc,Kloc), &
              DEps0(Mloc,Nloc,Kloc),Prod_s(Mloc,Nloc,Kloc),Prod_b(Mloc,Nloc,Kloc))

     ! fluxes for construction at cell faces    
     ALLOCATE(UxL(Mloc1,Nloc,Kloc),UxR(Mloc1,Nloc,Kloc),VxL(Mloc1,Nloc,Kloc),VxR(Mloc1,Nloc,Kloc), &
              WxL(Mloc1,Nloc,Kloc),WxR(Mloc1,Nloc,Kloc),DUxL(Mloc1,Nloc,Kloc),DUxR(Mloc1,Nloc,Kloc), &
              DVxL(Mloc1,Nloc,Kloc),DVxR(Mloc1,Nloc,Kloc),DWxL(Mloc1,Nloc,Kloc),DWxR(Mloc1,Nloc,Kloc), &
              UyL(Mloc,Nloc1,Kloc),UyR(Mloc,Nloc1,Kloc),VyL(Mloc,Nloc1,Kloc),VyR(Mloc,Nloc1,Kloc), &
              WyL(Mloc,Nloc1,Kloc),WyR(Mloc,Nloc1,Kloc),DUyL(Mloc,Nloc1,Kloc),DUyR(Mloc,Nloc1,Kloc), &
              DVyL(Mloc,Nloc1,Kloc),DVyR(Mloc,Nloc1,Kloc),DWyL(Mloc,Nloc1,Kloc),DWyR(Mloc,Nloc1,Kloc), &
              UzL(Mloc,Nloc,Kloc1),UzR(Mloc,Nloc,Kloc1),VzL(Mloc,Nloc,Kloc1),VzR(Mloc,Nloc,Kloc1), &
              WzL(Mloc,Nloc,Kloc1),WzR(Mloc,Nloc,Kloc1),OzL(Mloc,Nloc,Kloc1),OzR(Mloc,Nloc,Kloc1), &
              SxL(Mloc1,Nloc,Kloc),SxR(Mloc1,Nloc,Kloc),SyL(Mloc,Nloc1,Kloc),SyR(Mloc,Nloc1,Kloc), &
              ExL(Mloc1,Nloc,Kloc),ExR(Mloc1,Nloc,Kloc),FxL(Mloc1,Nloc,Kloc),FxR(Mloc1,Nloc,Kloc), &
              GxL(Mloc1,Nloc,Kloc),GxR(Mloc1,Nloc,Kloc),HxL(Mloc1,Nloc,Kloc),HxR(Mloc1,Nloc,Kloc), &
              EyL(Mloc,Nloc1,Kloc),EyR(Mloc,Nloc1,Kloc),FyL(Mloc,Nloc1,Kloc),FyR(Mloc,Nloc1,Kloc), &
              GyL(Mloc,Nloc1,Kloc),GyR(Mloc,Nloc1,Kloc),HyL(Mloc,Nloc1,Kloc),HyR(Mloc1,Nloc1,Kloc), &
              Ex(Mloc1,Nloc,Kloc),Ey(Mloc,Nloc1,Kloc),Fx(Mloc1,Nloc,Kloc),Fy(Mloc,Nloc1,Kloc), &
              Gx(Mloc1,Nloc,Kloc),Gy(Mloc,Nloc1,Kloc),Hx(Mloc1,Nloc,Kloc),Hy(Mloc,Nloc1,Kloc), &
              Fz(Mloc,Nloc,Kloc1),Gz(Mloc,Nloc,Kloc1),Hz(Mloc,Nloc,Kloc1),DelxU(Mloc,Nloc,Kloc), &
              DelyU(Mloc,Nloc,Kloc),DelzU(Mloc,Nloc,Kloc),DelxV(Mloc,Nloc,Kloc),DelyV(Mloc,Nloc,Kloc), &
              DelzV(Mloc,Nloc,Kloc),DelxW(Mloc,Nloc,Kloc),DelyW(Mloc,Nloc,Kloc),DelzW(Mloc,Nloc,Kloc), &
              DelxDU(Mloc,Nloc,Kloc),DelyDU(Mloc,Nloc,Kloc),DelxDV(Mloc,Nloc,Kloc),DelyDV(Mloc,Nloc,Kloc), &
              DelxDW(Mloc,Nloc,Kloc),DelyDW(Mloc,Nloc,Kloc),DelzO(Mloc,Nloc,Kloc)) 



# if defined (BUBBLE)
     ALLOCATE(Rbg(Mg),Wbg(Mg),Entrain(Mg))
     ALLOCATE(Nbg(Mloc,Nloc,Kloc,Mg),DNbg(Mloc,Nloc,Kloc,Mg),DNbg0(Mloc,Nloc,Kloc,Mg),Vbg(Mloc,Nloc,Kloc))
# endif

     ! poisson solver (for NSPCG use)
     neqns = (Iend-Ibeg+1)*(Jend-Jbeg+1)*(Kend-Kbeg+1)
     ALLOCATE(Coef(5*neqns,5*15),JCoef(5*15),Rhs(neqns))

     end subroutine allocate_variables
     

     subroutine index
!---------------------------------------------------
!    This subroutine is used to creat work index
!    Called by                   
!       main                                                                
!    Last update: 20/12/2010, Gangfeng Ma                                     
!---------------------------------------------------
     use global
     implicit none

# if defined (PARALLEL)
     dims(1)=PX
     dims(2)=PY
     periods(1)=.false.
     periods(2)=.false.
     if(PERIODIC) periods(2)=.true.
     coords(1)=0
     coords(2)=0

     call MPI_CART_CREATE(MPI_COMM_WORLD,ndims,dims, &
         periods,reorder,comm2d,ier)
     call MPI_CART_COORDS(comm2d,myid,2,coords,ier)

     npx=coords(1)
     npy=coords(2)

     call MPI_CART_SHIFT(comm2d,0,1,n_west,n_east,ier)
     call MPI_CART_SHIFT(comm2d,1,1,n_suth,n_nrth,ier)
# else
     PX=1
     PY=1
# endif

     ! local index
     Mloc = Mglob/PX+2*Nghost
     Nloc = Nglob/PY+2*Nghost
     Kloc = Kglob+2*Nghost
     Mloc1 = Mloc+1
     Nloc1 = Nloc+1
     Kloc1 = Kloc+1

     Ibeg = Nghost+1
     Iend = Mloc-Nghost
     Iend1 = Mloc1-Nghost
     Jbeg = Nghost+1
     Jend = Nloc-Nghost
     Jend1 = Nloc1-Nghost
     Kbeg = Nghost+1
     Kend = Kloc-Nghost
     Kend1 = Kloc1-Nghost

     end subroutine index


     subroutine read_input
!---------------------------------------------------
!    This subroutine is used to read input.txt
!    Called by 
!       main
!    Last update: 20/12/2010, Gangfeng Ma
!---------------------------------------------------
     use global
     use input_util
# if defined (INTEL)
     use IFPORT, only: rand
# endif
     implicit none
     character(len=80) :: FILE_NAME
     CHARACTER(LEN=80) :: WHAT
     real(SP) :: Segma,Celerity,Wave_Length,Wave_Number,  &
                 Fk,Fkdif,Theta_Calc,Wnumy,tmp,tmp1
     integer :: line,ierr,Iter,i,j
 
     ! log and error file
     open(3,file='log.txt')

     ! read from input.txt
     FILE_NAME='input.txt'

     ! title
     CALL GET_STRING_VAL(TITLE,FILE_NAME,'TITLE',line,ierr)
     IF(ierr==1)THEN
# if defined (PARALLEL)
       if(myid.eq.0) write(3,*) 'No TITLE in ', FILE_NAME, 'use default'
# else
       write(3,*) 'No TITLE in ', FILE_NAME, 'use default'
# endif
       TITLE='---TEST RUN---'
     ENDIF
# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,*)'---- LOG FILE ---'
     if(myid.eq.0) WRITE(3,*)TITLE
     if(myid.eq.0) WRITE(3,*)'--------------input start --------------'
# else
     WRITE(3,*)'---- LOG FILE ---'
     WRITE(3,*)TITLE
     WRITE(3,*)'--------------input start --------------'
# endif

     ! dimension                                             
     CALL GET_INTEGER_VAL(Mglob,FILE_NAME,'Mglob',line)
     CALL GET_INTEGER_VAL(Nglob,FILE_NAME,'Nglob',line)
     CALL GET_INTEGER_VAL(Kglob,FILE_NAME,'Kglob',line)
# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,'(A7,I5)')'Mglob= ',Mglob
     if(myid.eq.0) WRITE(3,'(A7,I5)')'Nglob= ',Nglob
     if(myid.eq.0) WRITE(3,'(A7,I5)')'Kglob= ',Kglob
# else
     WRITE(3,'(A7,I5)')'Mglob= ',Mglob
     WRITE(3,'(A7,I5)')'Nglob= ',Nglob
     WRITE(3,'(A7,I5)')'Kglob= ',Kglob
# endif

# if defined (PARALLEL)
     ! processor number
     CALL GET_INTEGER_VAL(PX,FILE_NAME,'PX',line)
     CALL GET_INTEGER_VAL(PY,FILE_NAME,'PY',line)
     if(myid.eq.0) WRITE(3,'(A4,I5)')'PX= ',PX
     if(myid.eq.0) WRITE(3,'(A4,I5)')'PY= ',PY
     if(PX*PY.ne.NumP) then
       if(myid.eq.0) WRITE(3,'(A6,I5)') 'NumP= ',NumP
       stop
     endif
# endif

     ! grid sizes
     CALL GET_Float_VAL(dx,FILE_NAME,'DX',line)
     CALL GET_Float_VAL(dy,FILE_NAME,'DY',line)
# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,'(A4,F8.4)')'DX= ',dx
     if(myid.eq.0) WRITE(3,'(A4,F8.4)')'DY= ',dy
# else
     WRITE(3,'(A4,F8.4)')'DX= ',dx
     WRITE(3,'(A4,F8.4)')'DY= ',dy
# endif

     ! vertical grid option
     call GET_INTEGER_VAL(Ivgrd,FILE_NAME,'IVGRD',line)
# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,'(A7,I3)')'Ivgrd= ',Ivgrd
# else
     WRITE(3,'(A7,I3)')'Ivgrd= ',Ivgrd
# endif

     ! time step
     CALL GET_Float_VAL(dt_ini,FILE_NAME,'DT_INI',line)
     CALL GET_Float_VAL(dt_min,FILE_NAME,'DT_MIN',line)
     CALL GET_Float_VAL(dt_max,FILE_NAME,'DT_MAX',line)
# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,'(A8,F6.4)')'DT_INI= ',dt_ini
     if(myid.eq.0) WRITE(3,'(A8,F6.4)')'DT_MIN= ',dt_min
     if(myid.eq.0) WRITE(3,'(A8,F6.4)')'DT_MAX= ',dt_max
# else
     WRITE(3,'(A8,F6.4)')'DT_INI= ',dt_ini
     WRITE(3,'(A8,F6.4)')'DT_MIN= ',dt_min
     WRITE(3,'(A8,F6.4)')'DT_MAX= ',dt_max
# endif

     ! result folder                                     
     CALL GET_STRING_VAL(RESULT_FOLDER,FILE_NAME,'RESULT_FOLDER',line,ierr)
# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,'(A15,A50)')'RESULT_FOLDER= ', RESULT_FOLDER
# else
     WRITE(3,'(A15,A50)')'RESULT_FOLDER= ', RESULT_FOLDER
# endif

     ! simulation steps and time
     call GET_INTEGER_VAL(SIM_STEPS,FILE_NAME,'SIM_STEPS',line)
     CALL GET_Float_VAL(TOTAL_TIME,FILE_NAME,'TOTAL_TIME',line)
     CALL GET_Float_VAL(Plot_Start,FILE_NAME,'PLOT_START',line)
     CALL GET_Float_VAL(Plot_Intv,FILE_NAME,'PLOT_INTV',line)
     CALL GET_Float_VAL(Screen_Intv,FILE_NAME,'SCREEN_INTV',line)
# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,'(A11,I12)')'SIM_STEPS= ', SIM_STEPS
     if(myid.eq.0) WRITE(3,'(A12,F8.2)')'TOTAL_TIME= ', TOTAL_TIME
     if(myid.eq.0) WRITE(3,'(A11,F8.2)')'PLOT_START= ', Plot_Start
     if(myid.eq.0) WRITE(3,'(A11,F8.2)')'PLOT_INTV= ', Plot_Intv
     if(myid.eq.0) WRITE(3,'(A13,F8.2)')'SCREEN_INTV= ', Screen_Intv
# else
     WRITE(3,'(A11,I12)')'SIM_STEPS= ', SIM_STEPS
     WRITE(3,'(A12,F8.2)')'TOTAL_TIME= ', TOTAL_TIME
     WRITE(3,'(A11,F8.2)')'PLOT_START= ', Plot_Start
     WRITE(3,'(A11,F8.2)')'PLOT_INTV= ', Plot_Intv
     WRITE(3,'(A13,F8.2)')'SCREEN_INTV= ', Screen_Intv
# endif

! fyshi added bathymetry file name 04/13/2012

     CALL GET_STRING_VAL(Depth_File,FILE_NAME,'DEPTH_FILE',line,ierr)
# if defined (PARALLEL)
       if(myid.eq.0) WRITE(3,'(A14,A20)')'Depth File= ', TRIM(Depth_File)
# else
       WRITE(3,'(A14,A20)')'Depth File= ', TRIM(Depth_File)
# endif


! DJN added landslide thickness file 07/08/2013

     CALL GET_STRING_VAL(Slide_File,FILE_NAME,'SLIDE_FILE',line,ierr)
     CALL GET_Float_VAL(SLIDE_DENSITY ,FILE_NAME,'SLIDE_DENSITY',line)
     CALL GET_Float_VAL(SLIDE_VISCOSITY ,FILE_NAME,'SLIDE_VISCOSITY',line)
     CALL GET_Float_VAL(SLIDE_DIFF ,FILE_NAME,'SLIDE_DIFF',line)
     CALL GET_Float_VAL(SLIDE_ROUGHNESS ,FILE_NAME,'SLIDE_ROUGHNESS',line)
     CALL GET_Float_VAL(SLIDE_MINTHICK ,FILE_NAME,'SLIDE_MINTHICK',line)
# if defined (PARALLEL)
       if(myid.eq.0) WRITE(3,'(A14,A20)')'Slide File= ', TRIM(Slide_File)
       if(myid.eq.0) WRITE(3,'(A12,F8.2)')'SLIDE_DENSITY= ', SLIDE_DENSITY
       if(myid.eq.0) WRITE(3,'(A12,F8.2)')'SLIDE_VISCOSITY= ', SLIDE_VISCOSITY
       if(myid.eq.0) WRITE(3,'(A12,F8.2)')'SLIDE_DIFF= ', SLIDE_DIFF
       if(myid.eq.0) WRITE(3,'(A12,F8.2)')'SLIDE_ROUGHNESS= ', SLIDE_ROUGHNESS
       if(myid.eq.0) WRITE(3,'(A12,F8.2)')'SLIDE_MINTHICK= ', SLIDE_MINTHICK
# else
       WRITE(3,'(A14,A20)')'Slide File= ', TRIM(Slide_File)
       WRITE(3,'(A12,F8.2)')'SLIDE_DENSITY= ', SLIDE_DENSITY
       WRITE(3,'(A12,F8.2)')'SLIDE_VISCOSITY= ', SLIDE_VISCOSITY
       WRITE(3,'(A12,F8.2)')'SLIDE_DIFF= ', SLIDE_DIFF
       WRITE(3,'(A12,F8.2)')'SLIDE_ROUGHNESS= ', SLIDE_ROUGHNESS
       WRITE(3,'(A12,F8.2)')'SLIDE_MINTHICK= ', SLIDE_MINTHICK
# endif

     ! courant number
     CALL GET_Float_VAL(CFL,FILE_NAME,'CFL',line)
# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,'(A5,F8.3)')'CFL= ',CFL
# else
     WRITE(3,'(A5,F8.3)')'CFL= ',CFL
# endif

     ! viscous number
     CALL GET_Float_VAL(VISCOUS_NUMBER,FILE_NAME,'VISCOUS_NUMBER',line)
# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,'(A16,F8.3)')'VISCOUS_NUMBER= ',VISCOUS_NUMBER
# else
     WRITE(3,'(A16,F8.3)')'VISCOUS_NUMBER= ',VISCOUS_NUMBER
# endif

     ! minimum depth
     CALL GET_Float_VAL(MinDep,FILE_NAME,'MinDep',line)
# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,'(A8,F8.3)')'MinDep= ',MinDep
# else
     WRITE(3,'(A8,F8.3)')'MinDep= ',MinDep
# endif

     ! laminar viscosity
     CALL GET_LOGICAL_VAL(VISCOUS_FLOW,FILE_NAME,'VISCOUS_FLOW',line)
     CALL GET_INTEGER_VAL(IVturb,FILE_NAME,'IVTURB',line)
     CALL GET_INTEGER_VAL(IHturb,FILE_NAME,'IHTURB',line)
     CALL GET_Float_VAL(Visc,FILE_NAME,'VISCOSITY',line)
     CALL GET_Float_VAL(Cvs,FILE_NAME,'Cvs',line)
     CALL GET_Float_VAL(Chs,FILE_NAME,'Chs',line)
# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,'(A14,L4)')'VISCOUS_FLOW= ',VISCOUS_FLOW
     if(myid.eq.0) WRITE(3,'(A8,I2)')'IVTURB= ',IVturb
     if(myid.eq.0) WRITE(3,'(A8,I2)')'IHTURB= ',IHturb
     if(myid.eq.0) WRITE(3,'(A11,F8.3)')'VISCOSITY= ',Visc
     if(myid.eq.0) WRITE(3,'(A6,F8.3)')'Cvs= ',Cvs
     if(myid.eq.0) WRITE(3,'(A6,F8.3)')'Chs= ',Chs
# else
     WRITE(3,'(A14,L4)')'VISCOUS_FLOW= ',VISCOUS_FLOW
     WRITE(3,'(A7,I2)')'IVTURB= ',IVturb
     WRITE(3,'(A7,I2)')'IHTURB= ',IHturb
     WRITE(3,'(A11,F8.3)')'VISCOSITY= ',Visc
     WRITE(3,'(A6,F8.3)')'Cvs= ',Cvs
     WRITE(3,'(A6,F8.3)')'Chs= ',Chs
# endif

     ! bathymetry     
     CALL GET_STRING_VAL(DEPTH_TYPE,FILE_NAME,'DEPTH_TYPE',line,ierr)
# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,'(A12,A50)')'DEPTH_TYPE= ',DEPTH_TYPE
# else
     WRITE(3,'(A12,A50)')'DEPTH_TYPE= ',DEPTH_TYPE
# endif

     ! bottom roughness
     CALL GET_INTEGER_VAL(Ibot,FILE_NAME,'Ibot',line)
     CALL GET_Float_VAL(Cd0,FILE_NAME,'Cd0',line)
     CALL GET_Float_VAL(Zob,FILE_NAME,'Zob',line)
# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,'(A6,I2)')'Ibot= ',Ibot
     if(myid.eq.0) WRITE(3,'(A5,F8.5)')'Cd0= ',Cd0
     if(myid.eq.0) WRITE(3,'(A5,F8.5)')'Zob= ',Zob
# else
     WRITE(3,'(A6,I2)')'Ibot= ',Ibot
     WRITE(3,'(A5,F8.5)')'Cd0= ',Cd0
     WRITE(3,'(A5,F8.5)')'Zob= ',Zob
# endif

     ! froudecap
     CALL GET_Float_VAL(FroudeCap,FILE_NAME,'FroudeCap',line)

# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,'(A12,F8.5)')'FroudeCap= ',FroudeCap
# else
     WRITE(3,'(A12,F8.5)')'FroudeCap= ',FroudeCap
# endif

     ! numerical scheme
     CALL GET_STRING_VAL(HIGH_ORDER,FILE_NAME,'HIGH_ORDER',line,ierr)
     CALL GET_STRING_VAL(TIME_ORDER,FILE_NAME,'TIME_ORDER',line,ierr)
     IF(ierr==1)THEN
# if defined (PARALLEL)
       if(myid.eq.0) WRITE(3,'(A12,A50)')'HIGH_ORDER', 'NOT DEFINED, USE DEFAULT'
# else
       WRITE(3,'(A12,A50)')'HIGH_ORDER', 'NOT DEFINED, USE DEFAULT'
# endif
       HIGH_ORDER='SECOND'
     ENDIF
# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,'(A12,A50)')'HIGH_ORDER= ', HIGH_ORDER
# else
     WRITE(3,'(A12,A50)')'HIGH_ORDER= ', HIGH_ORDER
# endif
     IF(ierr==1)THEN
# if defined (PARALLEL)
       if(myid.eq.0) WRITE(3,'(A12,A50)')'TIME_ORDER', 'NOT DEFINED, USE DEFAULT'
# else
       WRITE(3,'(A12,A50)')'TIME_ORDER', 'NOT DEFINED, USE DEFAULT'
# endif
       TIME_ORDER='THIRD'
     ENDIF
# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,'(A12,A50)')'TIME_ORDER= ', TIME_ORDER
# else
     WRITE(3,'(A12,A50)')'TIME_ORDER= ', TIME_ORDER
# endif

     ! if non-hydrostatic simulation
     CALL GET_LOGICAL_VAL(NON_HYDRO,FILE_NAME,'NON_HYDRO',line)
# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,'(A11,L4)')'NON_HYDRO= ',NON_HYDRO
# else
     WRITE(3,'(A11,L4)')'NON_HYDRO= ',NON_HYDRO
# endif

     ! poisson solver
     CALL GET_INTEGER_VAL(isolver,FILE_NAME,'ISOLVER',line)
     CALL GET_INTEGER_VAL(itmax,FILE_NAME,'ITMAX',line)
     CALL GET_Float_VAL(tol,FILE_NAME,'TOL',line)
# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,'(A9,I2)')'ISOLVER= ',isolver
     if(myid.eq.0) WRITE(3,'(A7,I5)')'ITMAX= ',itmax
     if(myid.eq.0) WRITE(3,'(A5,E12.3)')'TOL= ',tol
# else
     WRITE(3,'(A9,I2)')'ISOLVER= ',isolver
     WRITE(3,'(A7,I5)')'ITMAX= ',itmax
     WRITE(3,'(A5,E12.3)')'TOL= ',tol
# endif

     ! periodic bc
     CALL GET_LOGICAL_VAL(PERIODIC,FILE_NAME,'PERIODIC',line)
# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,'(A10,L4)')'PERIODIC= ',PERIODIC
# else
     WRITE(3,'(A10,L4)')'PERIODIC= ',PERIODIC
# endif    

     ! boundary type
     CALL GET_INTEGER_VAL(Bc_X0,FILE_NAME,'BC_X0',line)
     CALL GET_INTEGER_VAL(Bc_Xn,FILE_NAME,'BC_Xn',line)
     CALL GET_INTEGER_VAL(Bc_Y0,FILE_NAME,'BC_Y0',line)
     CALL GET_INTEGER_VAL(Bc_Yn,FILE_NAME,'BC_Yn',line)
     CALL GET_INTEGER_VAL(Bc_Z0,FILE_NAME,'BC_Z0',line)
     CALL GET_INTEGER_VAL(Bc_Zn,FILE_NAME,'BC_Zn',line)
# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,'(A7,I2)')'BC_X0= ',Bc_X0
     if(myid.eq.0) WRITE(3,'(A7,I2)')'BC_Xn= ',Bc_Xn
     if(myid.eq.0) WRITE(3,'(A7,I2)')'BC_Y0= ',Bc_Y0
     if(myid.eq.0) WRITE(3,'(A7,I2)')'BC_Yn= ',Bc_Yn
     if(myid.eq.0) WRITE(3,'(A7,I2)')'BC_Z0= ',Bc_Z0
     if(myid.eq.0) WRITE(3,'(A7,I2)')'BC_Zn= ',Bc_Zn
# else
     WRITE(3,'(A7,I2)')'BC_X0= ',Bc_X0
     WRITE(3,'(A7,I2)')'BC_Xn= ',Bc_Xn
     WRITE(3,'(A7,I2)')'BC_Y0= ',Bc_Y0
     WRITE(3,'(A7,I2)')'BC_Yn= ',Bc_Yn
     WRITE(3,'(A7,I2)')'BC_Z0= ',Bc_Z0
     WRITE(3,'(A7,I2)')'BC_Zn= ',Bc_Zn
# endif

     ! wavemaker 
     CALL GET_STRING_VAL(WaveMaker,FILE_NAME,'WAVEMAKER',line,ierr)
# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,'(A11,A50)')'WAVEMAKER= ', WAVEMAKER
# else
     WRITE(3,'(A11,A50)')'WAVEMAKER= ', WAVEMAKER
# endif
     IF(WaveMaker(1:10)=='TID_FLX_LR'.or.  &
        WaveMaker(1:10)=='TID_ELE_LR')THEN
       CALL GET_STRING_VAL(BoundaryFile,FILE_NAME,'BoundaryFile',line,ierr)
       CALL GET_LOGICAL_VAL(Nudging,FILE_NAME,'Nudging',line)

# if defined (PARALLEL)
       if(myid.eq.0) WRITE(3,'(A14,A20)')'BoundaryFile= ', TRIM(BoundaryFile)
# else
       WRITE(3,'(A14,A20)')'BoundaryFile= ', TRIM(BoundaryFile)
# endif
     ENDIF

     IF(WaveMaker(1:3)=='LEF'.or.WaveMaker(1:3)=='INT')THEN
       CALL GET_Float_VAL(Amp_Wave,FILE_NAME,'AMP',line)
       CALL GET_Float_VAL(Per_Wave,FILE_NAME,'PER',line)
       CALL GET_Float_VAL(Dep_Wave,FILE_NAME,'DEP',line)
       CALL GET_Float_VAL(Theta_Wave,FILE_NAME,'THETA',line)
# if defined (PARALLEL)
       if(myid.eq.0) WRITE(3,'(A9,F6.3)')'AMP_WAVE= ', Amp_Wave
       if(myid.eq.0) WRITE(3,'(A9,F6.3)')'PER_WAVE= ', Per_Wave
       if(myid.eq.0) WRITE(3,'(A9,F6.3)')'DEP_WAVE= ', Dep_Wave
       if(myid.eq.0) WRITE(3,'(A12,F6.3)')'THETA_WAVE= ', Theta_Wave
# else
       WRITE(3,'(A9,F6.3)')'AMP_WAVE= ', Amp_Wave
       WRITE(3,'(A9,F6.3)')'PER_WAVE= ', Per_Wave
       WRITE(3,'(A9,F6.3)')'DEP_WAVE= ', Dep_Wave
       WRITE(3,'(A12,F6.3)')'THETA_WAVE= ', Theta_Wave
# endif

       IF(WaveMaker(1:3)=='INT') then
         CALL GET_Float_VAL(Xsource_West,FILE_NAME,'Xsource_West',line)
         CALL GET_Float_VAL(Xsource_East,FILE_NAME,'Xsource_East',line)
         CALL GET_Float_VAL(Ysource_Suth,FILE_NAME,'Ysource_Suth',line)
         CALL GET_Float_VAL(Ysource_Nrth,FILE_NAME,'Ysource_Nrth',line)
# if defined (PARALLEL)
         if(myid.eq.0) WRITE(3,'(A14,F6.3)')'Xsource_West= ',Xsource_West
         if(myid.eq.0) WRITE(3,'(A14,F6.3)')'Xsource_East= ',Xsource_East
         if(myid.eq.0) WRITE(3,'(A14,F6.3)')'Ysource_Suth= ',Ysource_Suth
         if(myid.eq.0) WRITE(3,'(A14,F6.3)')'Ysource_Nrth= ',Ysource_Nrth
# else
         WRITE(3,'(A14,F6.3)')'Xsource_West= ',Xsource_West                                                                                    
         WRITE(3,'(A14,F6.3)')'Xsource_East= ',Xsource_East                                                                                    
         WRITE(3,'(A14,F6.3)')'Ysource_Suth= ',Ysource_Suth                                                                                    
         WRITE(3,'(A14,F6.3)')'Ysource_Nrth= ',Ysource_Nrth
# endif
       ENDIF

       ! test periodicity
       IF(PERIODIC.and.Theta_Wave.ne.Zero) then
         ! find wave number
         Segma = 2.0*pi/Per_Wave
         Celerity = sqrt(Grav*Dep_Wave)
         Wave_Length = Celerity*Per_Wave
         Wave_Number = 2.0*pi/Wave_Length

         Iter = 0
 75      Fk = Grav*Wave_Number*tanh(Wave_Number*Dep_Wave)-Segma**2
         if(abs(Fk)<=1.0e-8.or.Iter>1000) goto 85
         Fkdif = Grav*Wave_Number*Dep_Wave*(1.0-tanh(Wave_Number*Dep_Wave)**2)+  &                
            Grav*tanh(Wave_Number*Dep_Wave)
         Wave_Number = Wave_Number-Fk/Fkdif
         Iter = Iter+1
         goto 75
 85      continue

         if(Theta_Wave>Zero) then
           ! find right angle for periodic bc
           tmp = Large       
           do Iter = 1,10000
             Wnumy = Iter*2.0*pi/(Nglob*dy)
             if(WnumY<Wave_Number) then
               ! theta based on Ky = K*sin(theta)
               tmp1 = asin(Wnumy/Wave_Number)*180./pi
               if(abs(tmp1-Theta_Wave)<tmp) then
                 tmp = abs(tmp1-Theta_Wave)
                 Theta_Calc = tmp1
               endif
             endif
           enddo
         elseif(Theta_Wave<Zero) then
           ! find right angle for periodic bc 
           tmp = Large
           do Iter = 1,10000
             Wnumy = Iter*2.0*pi/(Nglob*dy)
             if(WnumY<Wave_Number) then
               ! theta based on Ky = K*sin(theta)
               tmp1 = -asin(Wnumy/Wave_Number)*180./pi
               if(abs(tmp1-Theta_Wave)<tmp) then
                 tmp = abs(tmp1-Theta_Wave)
                 Theta_Calc = tmp1
               endif
             endif
           enddo
         endif

# if defined (PARALLEL)
         if(myid.eq.0) then
           write(3,'(A20,F6.3)') 'Wave angle you set= ',Theta_Wave         
           write(3,'(A28,F6.3)') 'Wave angle for periodic bc= ',Theta_Calc
         endif
# else
         write(3,'(A20,F6.3)') 'Wave angle you set= ',Theta_Wave                                                      
         write(3,'(A28,F6.3)') 'Wave angle for periodic bc= ',Theta_Calc
# endif
         Theta_Wave = Theta_Calc
       ENDIF
     ENDIF

     ! random wave, read in 2d spectrum
     IF(WaveMaker(5:7)=='SPC') then
       open(14,file='spc2d.txt')
       read(14,*) NumFreq,NumDir
       if(NumFreq>MaxNumFreq) then
# if defined (PARALLEL)
         if(myid.eq.0) then
           write(3,'(A)') 'Please set a larger MaxNumFreq in mod_glob.F'
           stop
         endif
# else
         write(3,'(A)') 'Please set a larger MaxNumFreq in mod_glob.F'
         stop
# endif
       endif
       if(NumDir>MaxNumDir) then
# if defined (PARALLEL)
         if(myid.eq.0) then
           write(3,'(A)') 'Please set a larger MaxNumDir in mod_glob.F'
           stop
         endif
# else
         write(3,'(A)') 'Please set a larger MaxNumDir in mod_glob.F'
         stop
# endif
       endif
       do i = 1,NumFreq
         read(14,*) Freq(i)
       enddo
       do i = 1,NumDir
         read(14,*) Dire(i)
       enddo
       do j = 1,NumFreq
       do i = 1,NumDir
         read(14,*) Wave_Spc2d(i,j)
       enddo
       enddo
       close(14)

       ! random phase for each component
       do j = 1,NumFreq
       do i = 1,NumDir
# if defined (INTEL)
         Random_Phs(i,j) = rand()*2.0*pi
# else
         Random_Phs(i,j) = rand(0)*2.0*pi
# endif
       enddo
       enddo
     ENDIF

     ! sponge layer
     CALL GET_LOGICAL_VAL(SPONGE_ON,FILE_NAME,'SPONGE_ON',line)
# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,'(A11,L4)')'SPONGE_ON= ', SPONGE_ON
# else
     WRITE(3,'(A11,L4)')'SPONGE_ON= ', SPONGE_ON
# endif
     IF(SPONGE_ON)THEN
       CALL GET_Float_VAL(Sponge_West_Width,FILE_NAME,'Sponge_West_Width',line)
       CALL GET_Float_VAL(Sponge_East_Width,FILE_NAME,'Sponge_East_Width',line)
       CALL GET_Float_VAL(Sponge_South_Width,FILE_NAME,'Sponge_South_Width',line)
       CALL GET_Float_VAL(Sponge_North_Width,FILE_NAME,'Sponge_North_Width',line)
       CALL GET_Float_VAL(R_Sponge,FILE_NAME,'R_Sponge',line)
       CALL GET_Float_VAL(A_Sponge,FILE_NAME,'A_Sponge',line)
# if defined (PARALLEL)
       if(myid.eq.0) WRITE(3,'(A19,F6.3)')'Sponge_West_Width= ', Sponge_West_Width
       if(myid.eq.0) WRITE(3,'(A19,F6.3)')'Sponge_East_Width= ', Sponge_East_Width
       if(myid.eq.0) WRITE(3,'(A20,F6.3)')'Sponge_South_Width= ', Sponge_South_Width
       if(myid.eq.0) WRITE(3,'(A20,F6.3)')'Sponge_North_Width= ', Sponge_North_Width
       if(myid.eq.0) WRITE(3,'(A10,F6.3)')'R_Sponge= ', R_Sponge
       if(myid.eq.0) WRITE(3,'(A10,F6.3)')'A_Sponge= ', A_Sponge
# else
       WRITE(3,'(A19,F6.3)')'Sponge_West_Width= ', Sponge_West_Width
       WRITE(3,'(A19,F6.3)')'Sponge_East_Width= ', Sponge_East_Width
       WRITE(3,'(A20,F6.3)')'Sponge_South_Width= ', Sponge_South_Width
       WRITE(3,'(A20,F6.3)')'Sponge_North_Width= ', Sponge_North_Width
       WRITE(3,'(A10,F6.3)')'R_Sponge= ', R_Sponge
       WRITE(3,'(A10,F6.3)')'A_Sponge= ', A_Sponge
# endif
     ENDIF

!# if defined (LANDSLIDE)
!       CALL GET_Float_VAL(T_slide,FILE_NAME,'H_slide',line)
!       CALL GET_Float_VAL(b_slide,FILE_NAME,'L_slide',line)
!       CALL GET_Float_VAL(w_slide,FILE_NAME,'W_slide',line)
!       CALL GET_Float_VAL(e_slide,FILE_NAME,'e_slide',line)
!       CALL GET_Float_VAL(alpha_slide,FILE_NAME,'Angle_slide',line)
!       CALL GET_Float_VAL(x0_slide,FILE_NAME,'X0_slide',line)
!       CALL GET_Float_VAL(y0_slide,FILE_NAME,'Y0_slide',line)
!       CALL GET_Float_VAL(slope_slide,FILE_NAME,'Slope_slide',line)
!       CALL GET_Float_VAL(term_v,FILE_NAME,'TermV_slide',line)
!       IF(ANA_BATHY)THEN
!       CALL GET_Float_VAL(acceleration_lab,FILE_NAME,'ACC_lab',line)
!       ENDIF
!# if defined (PARALLEL)
!       if(myid.eq.0) WRITE(3,'(A15,F12.3)')'Slide height= ', T_slide
!       if(myid.eq.0) WRITE(3,'(A15,F12.3)')'Slide length= ', b_slide
!       if(myid.eq.0) WRITE(3,'(A15,F12.3)')'Slide width = ', w_slide
!       if(myid.eq.0) WRITE(3,'(A15,F12.3)')'Slide e para= ', e_slide
!       if(myid.eq.0) WRITE(3,'(A15,F12.3)')'Slide angle = ', alpha_slide
!       if(myid.eq.0) WRITE(3,'(A15,F12.3)')'Slide init x= ', x0_slide
!       if(myid.eq.0) WRITE(3,'(A15,F12.3)')'Slide init y= ', y0_slide
!       if(myid.eq.0) WRITE(3,'(A15,F12.3)')'Slide slope = ', slope_slide
!       if(myid.eq.0) WRITE(3,'(A15,F12.3)')'Slide term V= ', term_v
!# else
!       WRITE(3,'(A15,F12.3)')'Slide height= ', T_slide
!       WRITE(3,'(A15,F12.3)')'Slide length= ', b_slide
!       WRITE(3,'(A15,F12.3)')'Slide width = ', w_slide
!       WRITE(3,'(A15,F12.3)')'Slide e para= ', e_slide
!       WRITE(3,'(A15,F12.3)')'Slide angle = ', alpha_slide
!       WRITE(3,'(A15,F12.3)')'Slide init x= ', x0_slide
!       WRITE(3,'(A15,F12.3)')'Slide init y= ', y0_slide
!       WRITE(3,'(A15,F12.3)')'Slide slope = ', slope_slide
!       WRITE(3,'(A15,F12.3)')'Slide term V= ', term_v
!# endif
!# endif
     ! probe output
     CALL GET_INTEGER_VAL(NSTAT,FILE_NAME,'NSTAT',line)
     CALL GET_Float_VAL(Plot_Intv_Stat,FILE_NAME,'PLOT_INTV_STAT',line)
# if defined (PARALLEL)
       if(myid.eq.0) WRITE(3,'(A7,I3)')'NSTAT= ', NSTAT
       if(myid.eq.0) WRITE(3,'(A19,F6.3)')'Plot_Intv_Stat= ', Plot_Intv_Stat
# else
       WRITE(3,'(A7,I3)')'NSTAT= ', NSTAT
       WRITE(3,'(A19,F6.3)')'Plot_Intv_Stat= ', Plot_Intv_Stat
# endif

     if(NSTAT>0) then
       ALLOCATE(XSTAT(NSTAT),YSTAT(NSTAT))
       open(15,file='stat.txt',status='old')
       do i = 1,NSTAT
         read(15,*) xstat(i),ystat(i)
       enddo
       close(15)
     endif

     ! output parameters  
     CALL GET_LOGICAL_VAL(OUT_H,FILE_NAME,'OUT_H',line)
     CALL GET_LOGICAL_VAL(OUT_E,FILE_NAME,'OUT_E',line)
     CALL GET_LOGICAL_VAL(OUT_U,FILE_NAME,'OUT_U',line)
     CALL GET_LOGICAL_VAL(OUT_V,FILE_NAME,'OUT_V',line)
     CALL GET_LOGICAL_VAL(OUT_W,FILE_NAME,'OUT_W',line)
     CALL GET_LOGICAL_VAL(OUT_2DUV,FILE_NAME,'OUT_2DUV',line)
     CALL GET_LOGICAL_VAL(OUT_P,FILE_NAME,'OUT_P',line)
     CALL GET_LOGICAL_VAL(OUT_K,FILE_NAME,'OUT_K',line)
     CALL GET_LOGICAL_VAL(OUT_D,FILE_NAME,'OUT_D',line)
     CALL GET_LOGICAL_VAL(OUT_S,FILE_NAME,'OUT_S',line)
     CALL GET_LOGICAL_VAL(OUT_C,FILE_NAME,'OUT_C',line)
     CALL GET_LOGICAL_VAL(OUT_B,FILE_NAME,'OUT_B',line)
! fyshi add tmp4preview
     CALL GET_LOGICAL_VAL(OUT_preview,FILE_NAME,'OUT_preview',line)

     CALL GET_LOGICAL_VAL(OUT_SLD,FILE_NAME,'OUT_SLD',line) ! DJN 


# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,'(A7,L4)') 'OUT_H= ',OUT_H
     if(myid.eq.0) WRITE(3,'(A7,L4)') 'OUT_E= ',OUT_E
     if(myid.eq.0) WRITE(3,'(A7,L4)') 'OUT_U= ',OUT_U
     if(myid.eq.0) WRITE(3,'(A7,L4)') 'OUT_V= ',OUT_V
     if(myid.eq.0) WRITE(3,'(A7,L4)') 'OUT_W= ',OUT_W
     if(myid.eq.0) WRITE(3,'(A7,L4)') 'OUT_P= ',OUT_P
     if(myid.eq.0) WRITE(3,'(A7,L4)') 'OUT_K= ',OUT_K
     if(myid.eq.0) WRITE(3,'(A7,L4)') 'OUT_D= ',OUT_D
     if(myid.eq.0) WRITE(3,'(A7,L4)') 'OUT_S= ',OUT_S
     if(myid.eq.0) WRITE(3,'(A7,L4)') 'OUT_C= ',OUT_C
     if(myid.eq.0) WRITE(3,'(A7,L4)') 'OUT_B= ',OUT_B
     if(myid.eq.0) WRITE(3,'(A9,L4)') 'OUT_SLD= ',OUT_SLD   !DJN
# else
     WRITE(3,'(A7,L4)') 'OUT_H= ',OUT_H
     WRITE(3,'(A7,L4)') 'OUT_E= ',OUT_E
     WRITE(3,'(A7,L4)') 'OUT_U= ',OUT_U
     WRITE(3,'(A7,L4)') 'OUT_V= ',OUT_V
     WRITE(3,'(A7,L4)') 'OUT_W= ',OUT_W
     WRITE(3,'(A7,L4)') 'OUT_P= ',OUT_P
     WRITE(3,'(A7,L4)') 'OUT_K= ',OUT_K
     WRITE(3,'(A7,L4)') 'OUT_D= ',OUT_D
     WRITE(3,'(A7,L4)') 'OUT_S= ',OUT_S
     WRITE(3,'(A7,L4)') 'OUT_C= ',OUT_C
     WRITE(3,'(A7,L4)') 'OUT_B= ',OUT_B
     WRITE(3,'(A9,L4)') 'OUT_SLD= ',OUT_SLD   !DJN
# endif

# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,*)'--------------input end --------------'
# else
     WRITE(3,*)'--------------input end --------------'
# endif

! fyshi add some must-parameters 12/16/2011
     IF(WaveMaker(1:10)=='TID_FLX_LR')THEN
       BC_X0 = 3
       BC_Xn = 3
     ELSEIF(WaveMaker(1:10)=='TID_ELE_LR')THEN
       BC_X0 = 1
       BC_Xn = 1
     ENDIF
     IF(WaveMaker(1:10)=='TID_FLX_LR'.OR.  &
        WaveMaker(1:10)=='TID_ELE_LR') THEN
       OPEN(2,FILE=TRIM(BoundaryFile))
         READ(2,*)WHAT
         READ(2,*)NumTimeData
         ALLOCATE (TimeData(NumTimeData),  &
                   DataU_L(NumTimeData), DataU_R(NumTimeData), &   
                   DataEta_L(NumTimeData), DataEta_R(NumTimeData) )      
         DO I=1,NumTimeData
           READ(2,*,END=111) TimeData(I),DataEta_L(I),DataU_L(I), &
                                         DataEta_R(I),DataU_R(I) 
         ENDDO
111      CONTINUE
# if defined (PARALLEL)
        IF(myid==0) WRITE(*,*) 'Number of Data: ', I-1
# else
        WRITE(*,*) 'Number of Data: ', I-1
# endif
       CLOSE(2)
     ENDIF

     end subroutine read_input


     subroutine calculate_sponge
!-------------------------------------------------
!    Calculate sponge function
!    Called by
!      initial
!    Last update: 12/02/2011, Gangfeng Ma
!------------------------------------------------
     use global
     implicit none
     integer :: i,j,Iwidth
     real(SP) :: ri,lim

# if defined (PARALLEL)
     if(n_west.eq.MPI_PROC_NULL) then
# endif
     if(Sponge_West_Width>Zero) then
       Iwidth = int(Sponge_West_Width/dx)+Nghost
       do j = 1,Nloc
       do i = 1,Iwidth
         if(Sponge(i,j)>1.0) then
           lim = Sponge(i,j)
         else
           lim = 1.0
         endif
         ri = R_Sponge**(50*(i-1)/(Iwidth-1))
         Sponge(i,j) = max(A_Sponge**ri,lim)
       enddo
       enddo
     endif
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_east.eq.MPI_PROC_NULL) then
# endif
     if(Sponge_East_Width>Zero)then
       Iwidth = int(Sponge_East_Width/dx)+Nghost
       do j = 1,Nloc
       do i = Mloc-Iwidth+1,Mloc
         if(Sponge(i,j)>1.0) then
           lim = Sponge(i,j)
         else
           lim = 1.0
         endif
         ri = R_Sponge**(50*(Mloc-i)/(Iwidth-1))
         Sponge(i,j) = max(A_Sponge**ri,lim)
       enddo
       enddo
     endif
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_suth.eq.MPI_PROC_NULL) then
# endif
     if(Sponge_South_Width>Zero)then
       Iwidth = int(Sponge_South_Width/dy)+Nghost
       do j = 1,Iwidth
       do i = 1,Mloc
         if(Sponge(i,j)>1.0) then
           lim = Sponge(i,j)
         else
           lim = 1.0
         endif
         ri = R_Sponge**(50*(j-1)/(Iwidth-1))
         Sponge(i,j) = max(A_Sponge**ri,lim)
       enddo
       enddo
     endif
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_nrth.eq.MPI_PROC_NULL) then
# endif
     if(Sponge_North_Width>Zero)then
       Iwidth = int(Sponge_North_Width/dy)+Nghost
       do j = Nloc-Iwidth+1,Nloc
       do i = 1,Mloc
         if(Sponge(i,j)>1.0) then
           lim = Sponge(i,j)
         else
           lim = 1.0
         endif
         ri = R_Sponge**(50*(Nloc-j)/(Iwidth-1))
         Sponge(i,j) = max(A_Sponge**ri,lim)
       enddo
       enddo
     endif
# if defined (PARALLEL)
     endif
# endif

     end subroutine calculate_sponge

     
     subroutine sponge_damping
!---------------------------------------------------
!    This subroutine is used to damp waves using DHI type
!    sponge layer variables
!    Called by 
!      main
!    Last update: 12/02/2011, Gangfeng Ma
!--------------------------------------------------
     use global, only: Eta,Hc,D,U,V,W,Omega,Sponge,Mask, &
                       Mloc,Nloc,Kloc,DU,DV,DW
     implicit none
     integer :: i,j,k

     do j = 1,Nloc
     do i = 1,Mloc
       if(Mask(i,j)==1) then
         Eta(i,j) = Eta(i,j)/Sponge(i,j)
         D(i,j) = Eta(i,j)+Hc(i,j)

         do k = 1,Kloc
           U(i,j,k) = U(i,j,k)/Sponge(i,j)
           V(i,j,k) = V(i,j,k)/Sponge(i,j)
           W(i,j,k) = W(i,j,k)/Sponge(i,j)
           Omega(i,j,k) = Omega(i,j,k)/Sponge(i,j)
           DU(i,j,k) = D(i,j)*U(i,j,k)
           DV(i,j,k) = D(i,j)*V(i,j,k)
           DW(i,j,k) = D(i,j)*W(i,j,k)
         enddo
       endif
     enddo
     enddo
     
     end subroutine sponge_damping

 
     subroutine sigma_transform
!--------------------------------------------------- 
!    Calculate sigma transformation coefficient
!    Called by       
!      eval_duvw
!    Last update: 29/03/2011, Gangfeng Ma
!--------------------------------------------------
     use global, only: Zero,DelxSc,DelySc,D,DelxH,DelyH, &
                       DelxEta,DelyEta,sigc,Mloc,Nloc,Kloc
     implicit none
     integer :: i,j,k

     DelxSc = Zero
     DelySc = Zero
     do k = 1,Kloc
     do j = 1,Nloc
     do i = 1,Mloc
       DelxSc(i,j,k) = (1.0-sigc(k))/D(i,j)*DelxH(i,j)-sigc(k)/D(i,j)*DelxEta(i,j)
       DelySc(i,j,k) = (1.0-sigc(k))/D(i,j)*DelyH(i,j)-sigc(k)/D(i,j)*DelyEta(i,j)
     enddo
     enddo
     enddo

     end subroutine sigma_transform


     subroutine eval_turb(ISTEP)
!---------------------------------------------------
!    This subroutine is used to calculate viscosity
!    Called by                                                                                                             
!      main
!    Last update: 21/06/2011, Gangfeng Ma 
!--------------------------------------------------
     use global
     implicit none
     integer, intent(in) :: ISTEP
     integer :: i,j,k
     real(SP) :: DelsU,DelsV,Strxx,Stryy,Strxy

     ! laminar viscosity
     do k = 1,Kloc
     do j = 1,Nloc
     do i = 1,Mloc
       Cmu(i,j,k) = Visc
     enddo
     enddo
     enddo

     ! vertical turbulent viscosity
     if(IVturb==1) then
       ! constant vertical viscosity
       CmuVt = Cvs
     elseif(IVturb==2) then
       ! subgrid model
       do i = 1,Mloc
       do j = 1,Nloc
       do k = 2,Kloc-1
         DelsU = (U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))
         DelsV = (V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))
         CmuVt(i,j,k) = Cmu(i,j,k)+  &
            Cvs*D(i,j)*dsig(k)**2*sqrt(DelsU**2+DelsV**2)
       enddo
       CmuVt(i,j,1) = CmuVt(i,j,2)
       CmuVt(i,j,Kloc) = CmuVt(i,j,Kloc-1)
       enddo
       enddo
     elseif(IVturb==3) then
       ! k-epsilon turbulence model
       call kepsilon(ISTEP)
     elseif(IVturb==10) then
       ! 3D turbulence model
       call kepsilon_3D(ISTEP)
     elseif(IVturb==20) then
       call les_3D(ISTEP)
     endif
     
     ! horizontal turbulent viscosity
     if(IHturb==1) then
       ! constant viscosity
       CmuHt = Chs
     elseif(IHturb==2) then
       ! subgrid model
       do i = Ibeg,Iend
       do j = Jbeg,Jend
       do k = Kbeg,Kend
         Strxx = (U(i+1,j,k)-U(i-1,j,k))/(2.0*dx)+  &
               (U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k)
         Stryy = (V(i,j+1,k)-V(i,j-1,k))/(2.0*dy)+  &
               (V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k)
         Strxy = 0.5*((U(i,j+1,j)-U(i,j-1,k))/(2.0*dy)+  &
               (U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k)+  &
               (V(i+1,j,k)-V(i-1,j,k))/(2.0*dx)+  &
               (V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k))
         CmuHt(i,j,k) = Cmu(i,j,k)+  &
             Chs*dx*dy*sqrt(Strxx**2+2.0*Strxy**2+Stryy**2)
       enddo
       enddo
       enddo

       ! ghost cell
# if defined (PARALLEL)
       call phi_3D_exch(CmuHt)
# endif

# if defined (PARALLEL)
       if(n_west.eq.MPI_PROC_NULL) then
# endif
       do j = Jbeg,Jend
       do k = Kbeg,Kend
         do i = 1,Nghost
           CmuHt(Ibeg-i,j,k) = CmuHt(Ibeg+i-1,j,k)
         enddo
       enddo
       enddo
# if defined (PARALLEL)
       endif
# endif

# if defined (PARALLEL)
       if(n_east.eq.MPI_PROC_NULL) then
# endif
       do j = Jbeg,Jend
       do k = Kbeg,Kend
         do i = 1,Nghost
           CmuHt(Iend+i,j,k) = CmuHt(Iend-i+1,j,k)
         enddo
       enddo
       enddo
# if defined (PARALLEL)
       endif
# endif

# if defined (PARALLEL)
       if(n_suth.eq.MPI_PROC_NULL) then
# endif
       do i = Ibeg,Iend
       do k = Kbeg,Kend
         do j = 1,Nghost
           CmuHt(i,Jbeg-j,k) = CmuHt(i,Jbeg+j-1,k)
         enddo
       enddo
       enddo
# if defined (PARALLEL)
       endif
# endif

# if defined (PARALLEL)
       if(n_nrth.eq.MPI_PROC_NULL) then
# endif
       do i = Ibeg,Iend
       do k = Kbeg,Kend
         do j = 1,Nghost
           CmuHt(i,Jend+j,k) = CmuHt(i,Jend-j+1,k)
         enddo
       enddo
       enddo
# if defined (PARALLEL)
       endif
# endif

       do i = Ibeg,Iend
       do j = Jbeg,Jend
         do k = 1,Nghost
           CmuHt(i,j,Kbeg-k) = CmuHt(i,j,Kbeg+k-1)
         enddo
         do k = 1,Nghost
           CmuHt(i,j,Kend+k) = CmuHt(i,j,Kend-k+1)
         enddo
       enddo
       enddo  
     elseif(IHturb>=10) then
       ! use 3D turbulence model
       ! in this case, the length scales in all directions are
       ! in the same order
       CmuHt = CmuVt
     endif

     end subroutine eval_turb


     subroutine diffusion
!---------------------------------------------------
!    This subroutine is used to evaluate diffusion terms
!    Called by                                          
!      eval_duvw 
!    Last update: 29/03/2011, Gangfeng Ma 
!--------------------------------------------------
     use global
     implicit none
     integer :: i,j,k
     real(SP) :: Txx_L,Txx_R,Txx_U,Txx_D,Txy_L,Txy_R,Txy_F,Txy_B,  &
                 Txy_U,Txy_D,Txz_L,Txz_R,Txz_U,Txz_D,Tyy_F,Tyy_B,  &
                 Tyy_U,Tyy_D,Tyz_F,Tyz_B,Tyz_U,Tyz_D,Tzz_U,Tzz_D

     Diffxx = zero; Diffxy = zero; Diffxz = zero
     Diffyx = zero; Diffyy = zero; Diffyz = zero
     Diffzx = zero; Diffzy = zero; Diffzz = zero

     do k = Kbeg,Kend
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       if(Mask9(i,j)==0) cycle

       Txx_L = (CmuHt(i-1,j,k)+CmuHt(i,j,k))*((U(i,j,k)-U(i-1,j,k))/dx+  &
           0.25*(U(i,j,k+1)+U(i-1,j,k+1)-U(i,j,k-1)-U(i-1,j,k-1))/(sigc(k+1)-sigc(k-1))*  &
           (DelxSc(i,j,k)+DelxSc(i-1,j,k)))
       Txx_R = (CmuHt(i,j,k)+CmuHt(i+1,j,k))*((U(i+1,j,k)-U(i,j,k))/dx+  &  
           0.25*(U(i+1,j,k+1)+U(i,j,k+1)-U(i+1,j,k-1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))*  &
           (DelxSc(i+1,j,k)+DelxSc(i,j,k)))
       Txx_U = (CmuHt(i,j,k+1)+CmuHt(i,j,k))*(0.5*(U(i+1,j,k+1)+U(i+1,j,k)-U(i-1,j,k+1)-U(i-1,j,k))/(2.0*dx)+  &
           (U(i,j,k+1)-U(i,j,k))/(sigc(k+1)-sigc(k))*0.5*(DelxSc(i,j,k+1)+DelxSc(i,j,k)))
       Txx_D = (CmuHt(i,j,k)+CmuHt(i,j,k-1))*(0.5*(U(i+1,j,k)+U(i+1,j,k-1)-U(i-1,j,k)-U(i-1,j,k-1))/(2.0*dx)+  & 
           (U(i,j,k)-U(i,j,k-1))/(sigc(k)-sigc(k-1))*0.5*(DelxSc(i,j,k)+DelxSc(i,j,k-1)))
 
       Diffxx(i,j,k) = 0.5*((D(i+1,j)+D(i,j))*Txx_R-(D(i,j)+D(i-1,j))*Txx_L)/dx+  &
           D(i,j)*DelxSc(i,j,k)*(Txx_U-Txx_D)/dsig(k)

       Txy_L = 0.5*(CmuHt(i-1,j,k)+CmuHt(i,j,k))*(0.5*(U(i-1,j+1,k)+U(i,j+1,k)-U(i-1,j-1,k)-U(i,j-1,k))/(2.0*dy)+  &
           0.25*(U(i,j,k+1)+U(i-1,j,k+1)-U(i,j,k-1)-U(i-1,j,k-1))/(sigc(k+1)-sigc(k-1))*(DelySc(i,j,k)+DelySc(i-1,j,k))+  &
           (V(i,j,k)-V(i-1,j,k))/dx+0.25*(V(i,j,k+1)+V(i-1,j,k+1)-V(i,j,k-1)-V(i-1,j,k-1))/(sigc(k+1)-sigc(k-1))*  & 
           (DelxSc(i,j,k)+DelxSc(i-1,j,k)))
       Txy_R = 0.5*(CmuHt(i,j,k)+CmuHt(i+1,j,k))*(0.5*(U(i,j+1,k)+U(i+1,j+1,k)-U(i,j-1,k)-U(i+1,j-1,k))/(2.0*dy)+  &
           0.25*(U(i+1,j,k+1)+U(i,j,k+1)-U(i+1,j,k-1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))*(DelySc(i+1,j,k)+DelySc(i,j,k))+  &
           (V(i+1,j,k)-V(i,j,k))/dx+0.25*(V(i+1,j,k+1)+V(i,j,k+1)-V(i+1,j,k-1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))*  &
           (DelxSc(i+1,j,k)+DelxSc(i,j,k)))
       Txy_F = 0.5*(CmuHt(i,j-1,k)+CmuHt(i,j,k))*((U(i,j,k)-U(i,j-1,k))/dy+  &
           0.25*(U(i,j,k+1)+U(i,j-1,k+1)-U(i,j,k-1)-U(i,j-1,k-1))/(sigc(k+1)-sigc(k-1))*  &
           (DelySc(i,j,k)+DelySc(i,j-1,k))+0.5*(V(i+1,j,k)+V(i+1,j-1,k)-V(i-1,j,k)-V(i-1,j-1,k))/(2.0*dx)+  &
           0.25*(V(i,j,k+1)+V(i,j-1,k+1)-V(i,j,k-1)-V(i,j-1,k-1))/(sigc(k+1)-sigc(k-1))*  &
           (DelxSc(i,j,k)+DelxSc(i,j-1,k)))
       Txy_B = 0.5*(CmuHt(i,j,k)+CmuHt(i,j+1,k))*((U(i,j+1,k)-U(i,j,k))/dy+  & 
           0.25*(U(i,j+1,k+1)+U(i,j,k+1)-U(i,j+1,k-1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))*  & 
           (DelySc(i,j+1,k)+DelySc(i,j,k))+0.5*(V(i+1,j+1,k)+V(i+1,j,k)-V(i-1,j+1,k)-V(i-1,j,k))/(2.0*dx)+  & 
           0.25*(V(i,j+1,k+1)+V(i,j,k+1)-V(i,j+1,k-1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))*  & 
           (DelxSc(i,j+1,k)+DelxSc(i,j,k)))
       Txy_U = 0.5*(CmuHt(i,j,k)+CmuHt(i,j,k+1))*(0.5*(U(i,j+1,k+1)+U(i,j+1,k)-U(i,j-1,k+1)-U(i,j-1,k))/(2.0*dy)+ &
           0.5*(U(i,j,k+1)-U(i,j,k))/(sigc(k+1)-sigc(k))*(DelySc(i,j,k+1)+DelySc(i,j,k))+  &
           0.5*(V(i+1,j,k+1)+V(i+1,j,k)-V(i-1,j,k+1)-V(i-1,j,k))/(2.0*dx)+  &
           0.5*(V(i,j,k+1)-V(i,j,k))*(DelxSc(i,j,k+1)+DelxSc(i,j,k)))
       Txy_D = 0.5*(CmuHt(i,j,k-1)+CmuHt(i,j,k))*(0.5*(U(i,j+1,k)+U(i,j+1,k-1)-U(i,j-1,k)-U(i,j-1,k-1))/(2.0*dy)+ &
           0.5*(U(i,j,k)-U(i,j,k-1))/(sigc(k)-sigc(k-1))*(DelySc(i,j,k)+DelySc(i,j,k-1))+  & 
           0.5*(V(i+1,j,k)+V(i+1,j,k-1)-V(i-1,j,k)-V(i-1,j,k-1))/(2.0*dx)+  &  
           0.5*(V(i,j,k)-V(i,j,k-1))*(DelxSc(i,j,k)+DelxSc(i,j,k-1)))

       Diffxy(i,j,k) = 0.5*((D(i,j+1)+D(i,j))*Txy_B-(D(i,j)+D(i,j-1))*Txy_F)/dy+  &
           D(i,j)*DelySc(i,j,k)*(Txy_U-Txy_D)/dsig(k)
       Diffyx(i,j,k) = 0.5*((D(i+1,j)+D(i,j))*Txy_R-(D(i,j)+D(i-1,j))*Txy_L)/dx+  &
           D(i,j)*DelxSc(i,j,k)*(Txy_U-Txy_D)/dsig(k)

       Txz_L = 0.5*(CmuHt(i-1,j,k)+CmuHt(i,j,k))*((U(i-1,j,k+1)+U(i,j,k+1)-U(i-1,j,k-1)-U(i,j,k-1))/  &
           (sigc(k+1)-sigc(k-1))/(D(i-1,j)+D(i,j))+(W(i,j,k)-W(i-1,j,k))/dx+0.25*(W(i-1,j,k+1)+W(i,j,k+1)-  &
           W(i-1,j,k-1)-W(i,j,k-1))/(sigc(k+1)-sigc(k-1))*(DelxSc(i-1,j,k)+DelxSc(i,j,k)))
       Txz_R = 0.5*(CmuHt(i,j,k)+CmuHt(i+1,j,k))*((U(i,j,k+1)+U(i+1,j,k+1)-U(i,j,k-1)-U(i+1,j,k-1))/  & 
           (sigc(k+1)-sigc(k-1))/(D(i,j)+D(i+1,j))+(W(i+1,j,k)-W(i,j,k))/dx+0.25*(W(i,j,k+1)+W(i+1,j,k+1)-  &
           W(i,j,k-1)-W(i+1,j,k-1))/(sigc(k+1)-sigc(k-1))*(DelxSc(i,j,k)+DelxSc(i+1,j,k)))
       Txz_U = 0.5*(CmuVt(i,j,k+1)+CmuVt(i,j,k))*((U(i,j,k+1)-U(i,j,k))/(sigc(k+1)-sigc(k))/D(i,j)+  &
           0.5*(W(i+1,j,k+1)+W(i+1,j,k)-W(i-1,j,k+1)-W(i-1,j,k))/(2.0*dx)+  &
           0.5*(W(i,j,k+1)-W(i,j,k))/(sigc(k+1)-sigc(k))*(DelxSc(i,j,k)+DelxSc(i,j,k+1)))
       Txz_D = 0.5*(CmuVt(i,j,k)+CmuVt(i,j,k-1))*((U(i,j,k)-U(i,j,k-1))/(sigc(k)-sigc(k-1))/D(i,j)+  &
           0.5*(W(i+1,j,k)+W(i+1,j,k-1)-W(i-1,j,k)-W(i-1,j,k-1))/(2.0*dx)+  &  
           0.5*(W(i,j,k)-W(i,j,k-1))/(sigc(k)-sigc(k-1))*(DelxSc(i,j,k)+DelxSc(i,j,k-1)))

       Diffxz(i,j,k) = (Txz_U-Txz_D)/dsig(k)
       Diffzx(i,j,k) = 0.5*((D(i,j)+D(i+1,j))*Txz_R-(D(i-1,j)+D(i,j))*Txz_L)/dx+  &
           D(i,j)*DelxSc(i,j,k)*(Txz_U-Txz_D)/dsig(k)

       Tyy_F = (CmuHt(i,j-1,k)+CmuHt(i,j,k))*((V(i,j,k)-V(i,j-1,k))/dy+  &
           0.25*(V(i,j,k+1)+V(i,j-1,k+1)-V(i,j,k-1)-V(i,j-1,k-1))/(sigc(k+1)-sigc(k-1))*  &
           (DelySc(i,j,k)+DelySc(i,j-1,k)))
       Tyy_B = (CmuHt(i,j,k)+CmuHt(i,j+1,k))*((V(i,j+1,k)-V(i,j,k))/dy+  & 
           0.25*(V(i,j+1,k+1)+V(i,j,k+1)-V(i,j+1,k-1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))*  &  
           (DelySc(i,j+1,k)+DelySc(i,j,k)))
       Tyy_U = (CmuHt(i,j,k+1)+CmuHt(i,j,k))*(0.5*(V(i,j+1,k+1)+V(i,j+1,k)-V(i,j-1,k+1)-V(i,j-1,k))/(2.0*dy)+  &
           (V(i,j,k+1)-V(i,j,k))/(sigc(k+1)-sigc(k))*0.5*(DelySc(i,j,k+1)+DelySc(i,j,k)))  
       Tyy_D = (CmuHt(i,j,k)+CmuHt(i,j,k-1))*(0.5*(V(i,j+1,k)+V(i,j+1,k-1)-V(i,j-1,k)-V(i,j-1,k-1))/(2.0*dy)+  & 
           (V(i,j,k)-V(i,j,k-1))/(sigc(k)-sigc(k-1))*0.5*(DelySc(i,j,k)+DelySc(i,j,k-1)))

       Diffyy(i,j,k) = 0.5*((D(i,j+1)+D(i,j))*Tyy_B-(D(i,j)+D(i,j-1))*Tyy_F)/dy+  & 
           D(i,j)*DelySc(i,j,k)*(Tyy_U-Tyy_D)/dsig(k)
 
       Tyz_F = 0.5*(CmuHt(i,j-1,k)+CmuHt(i,j,k))*((V(i,j-1,k+1)+V(i,j,k+1)-V(i,j-1,k-1)-V(i,j,k-1))/  & 
           (sigc(k+1)-sigc(k-1))/(D(i,j-1)+D(i,j))+(W(i,j,k)-W(i,j-1,k))/dy+0.25*(W(i,j-1,k+1)+W(i,j,k+1)-  &
           W(i,j-1,k-1)-W(i,j,k-1))/(sigc(k+1)-sigc(k-1))*(DelySc(i,j-1,k)+DelySc(i,j,k)))
       Tyz_B = 0.5*(CmuHt(i,j,k)+CmuHt(i,j+1,k))*((V(i,j,k+1)+V(i,j+1,k+1)-V(i,j,k-1)-V(i,j+1,k-1))/  & 
           (sigc(k+1)-sigc(k-1))/(D(i,j)+D(i,j+1))+(W(i,j+1,k)-W(i,j,k))/dy+0.25*(W(i,j,k+1)+W(i,j+1,k+1)-  & 
           W(i,j,k-1)-W(i,j+1,k-1))/(sigc(k+1)-sigc(k-1))*(DelySc(i,j,k)+DelySc(i,j+1,k)))
       Tyz_U = 0.5*(CmuVt(i,j,k+1)+CmuVt(i,j,k))*((V(i,j,k+1)-V(i,j,k))/(sigc(k+1)-sigc(k))/D(i,j)+  & 
           0.5*(W(i,j+1,k+1)+W(i,j+1,k)-W(i,j-1,k+1)-W(i,j-1,k))/(2.0*dy)+  &
           0.5*(W(i,j,k+1)-W(i,j,k))/(sigc(k+1)-sigc(k))*(DelySc(i,j,k)+DelySc(i,j,k+1)))  
       Tyz_D = 0.5*(CmuVt(i,j,k)+CmuVt(i,j,k-1))*((V(i,j,k)-V(i,j,k-1))/(sigc(k)-sigc(k-1))/D(i,j)+  & 
           0.5*(W(i,j+1,k)+W(i,j+1,k-1)-W(i,j-1,k)-W(i,j-1,k-1))/(2.0*dy)+  & 
           0.5*(W(i,j,k)-W(i,j,k-1))/(sigc(k)-sigc(k-1))*(DelySc(i,j,k)+DelySc(i,j,k-1)))

       Diffyz(i,j,k) = (Tyz_U-Tyz_D)/dsig(k)
       Diffzy(i,j,k) = 0.5*((D(i,j)+D(i,j+1))*Tyz_B-(D(i,j-1)+D(i,j))*Tyz_F)/dy+  &
           D(i,j)*DelySc(i,j,k)*(Tyz_U-Tyz_D)/dsig(k)

       Tzz_U = (CmuVt(i,j,k+1)+CmuVt(i,j,k))*(W(i,j,k+1)-W(i,j,k))/(sigc(k+1)-sigc(k))/D(i,j)
       Tzz_D = (CmuVt(i,j,k)+CmuVt(i,j,k-1))*(W(i,j,k)-W(i,j,k-1))/(sigc(k)-sigc(k-1))/D(i,j)

       Diffzz(i,j,k) = (Tzz_U-Tzz_D)/dsig(k)

     enddo
     enddo
     enddo 


     end subroutine diffusion


# if defined (PARALLEL)
    subroutine phi_2D_exch(PHI)
    USE GLOBAL
    IMPLICIT NONE
    REAL(SP),INTENT(INOUT) :: PHI(Mloc,Nloc)

    INTEGER,DIMENSION(MPI_STATUS_SIZE,4) :: status
    INTEGER,DIMENSION(4) :: req
    INTEGER :: i,j,nreq,len
    REAL(SP),DIMENSION(Mloc,Nghost) :: rNmsg, sNmsg,rSmsg,sSmsg
    REAL(SP),DIMENSION(Nloc,Nghost) :: rWmsg, sWmsg,rEmsg,sEmsg

! for east-west

    len = Nloc * Nghost

    nreq = 0
    if ( n_west .ne. MPI_PROC_NULL ) then
       nreq = nreq + 1
       call MPI_IRECV( rWmsg, len, MPI_SP, &
            n_west, 0, comm2d, req(nreq), ier )
       do j = 1, Nloc
       do i = 1, Nghost
          sWmsg(j,i) = PHI(Ibeg+i-1,j)
       enddo
       enddo
       nreq = nreq +1
       call MPI_ISEND( sWmsg, len, MPI_SP, &
            n_west, 1, comm2d, req(nreq), ier )
    endif

    if ( n_east .ne. MPI_PROC_NULL ) then
       nreq = nreq + 1
       call MPI_IRECV( rEmsg, len, MPI_SP, &
            n_east, 1, comm2d, req(nreq), ier )
       do j = 1, Nloc
       do i = 1, Nghost
          sEmsg(j,i) = PHI(Iend-i+1,j)
       enddo
       enddo
       nreq = nreq +1
       call MPI_ISEND( sEmsg, len, MPI_SP, &
            n_east, 0, comm2d, req(nreq), ier )
    endif

    call MPI_WAITALL( nreq, req, status, ier )

    if ( n_west .ne. MPI_PROC_NULL ) then
       do j = 1, Nloc
       do i = 1, Nghost
          PHI(Ibeg-i,j) = rWmsg(j,i)
       enddo
       enddo
    endif

    if ( n_east .ne. MPI_PROC_NULL ) then
       do j = 1, Nloc
       do i = 1, Nghost
          PHI(Iend+i,j) = rEmsg(j,i)
       enddo
       enddo
    endif

! for nrth-suth

    len = Mloc * Nghost

    nreq = 0
    if ( n_suth .ne. MPI_PROC_NULL ) then
       nreq = nreq + 1
       call MPI_IRECV( rSmsg, len, MPI_SP, &
            n_suth, 0, comm2d, req(nreq), ier )
       do i = 1, Mloc
       do j = 1, Nghost
          sSmsg(i,j) = PHI(i,Jbeg+j-1)
       enddo
       enddo
       nreq = nreq +1
       call MPI_ISEND( sSmsg, len, MPI_SP, &
            n_suth, 1, comm2d, req(nreq), ier )
    endif

    if ( n_nrth .ne. MPI_PROC_NULL ) then
       nreq = nreq + 1
       call MPI_IRECV( rNmsg, len, MPI_SP, &
            n_nrth, 1, comm2d, req(nreq), ier )
       do i = 1, Mloc
       do j = 1, Nghost
          sNmsg(i,j) = PHI(i,Jend-j+1)
       enddo
       enddo
       nreq = nreq + 1
       call MPI_ISEND( sNmsg, len, MPI_SP, &
            n_nrth, 0, comm2d, req(nreq), ier )
    endif

    call MPI_WAITALL( nreq, req, status, ier )

    if ( n_suth .ne. MPI_PROC_NULL ) then
       do i = 1, Mloc
       do j = 1, Nghost
          PHI(i,Jbeg-j) = rSmsg(i,j)
       enddo
       enddo
    endif

    if ( n_nrth .ne. MPI_PROC_NULL ) then
       do i = 1, Mloc
       do j = 1, Nghost
          PHI(i,Jend+j) = rNmsg(i,j)
       enddo
       enddo
    endif

    return
    END SUBROUTINE phi_2D_exch
# endif


# if defined (PARALLEL)
    SUBROUTINE phi_3D_exch(PHI)
    USE GLOBAL
    IMPLICIT NONE
    REAL(SP),INTENT(INOUT) :: PHI(Mloc,Nloc,Kloc)

    INTEGER,DIMENSION(MPI_STATUS_SIZE,4) :: status
    INTEGER,DIMENSION(4) :: req
    INTEGER :: i,j,k,ik,jk,nreq,len
    REAL(SP),DIMENSION(Mloc*Kloc,Nghost) :: rNmsg, sNmsg,rSmsg,sSmsg
    REAL(SP),DIMENSION(Nloc*Kloc,Nghost) :: rWmsg, sWmsg,rEmsg,sEmsg

! for east-west

    len = Nloc * Kloc * Nghost

    nreq = 0
    if ( n_west .ne. MPI_PROC_NULL ) then
       nreq = nreq + 1
       call MPI_IRECV( rWmsg, len, MPI_SP, &
            n_west, 0, comm2d, req(nreq), ier )
       do k = 1, Kloc
       do j = 1, Nloc
       do i = 1, Nghost
          jk = (k-1)*Nloc+j
          sWmsg(jk,i) = PHI(Ibeg+i-1,j,k)
       enddo
       enddo
       enddo
       nreq = nreq + 1
       call MPI_ISEND( sWmsg, len, MPI_SP, &
            n_west, 1, comm2d, req(nreq), ier )
    endif

    if ( n_east .ne. MPI_PROC_NULL ) then
       nreq = nreq + 1
       call MPI_IRECV( rEmsg, len, MPI_SP, &
            n_east, 1, comm2d, req(nreq), ier )
       do k = 1, Kloc
       do j = 1, Nloc
       do i = 1, Nghost
          jk = (k-1)*Nloc+j
          sEmsg(jk,i) = PHI(Iend-i+1,j,k)
       enddo
       enddo
       enddo
       nreq = nreq +1
       call MPI_ISEND( sEmsg, len, MPI_SP, &
            n_east, 0, comm2d, req(nreq), ier )
    endif

    call MPI_WAITALL( nreq, req, status, ier )

    if ( n_west .ne. MPI_PROC_NULL ) then
       do k = 1, Kloc
       do j = 1, Nloc
       do i = 1, Nghost
          jk = (k-1)*Nloc+j
          PHI(Ibeg-i,j,k) = rWmsg(jk,i)
       enddo
       enddo
       enddo
    endif

    if ( n_east .ne. MPI_PROC_NULL ) then
       do k = 1, Kloc
       do j = 1, Nloc
       do i = 1, Nghost
          jk = (k-1)*Nloc+j
          PHI(Iend+i,j,k) = rEmsg(jk,i)
       enddo
       enddo
       enddo
    endif

! for nrth-suth

    len = Mloc * Kloc * Nghost

    nreq = 0
    if ( n_suth .ne. MPI_PROC_NULL ) then
       nreq = nreq + 1
       call MPI_IRECV( rSmsg, len, MPI_SP, &
            n_suth, 0, comm2d, req(nreq), ier )
       do k = 1, Kloc
       do i = 1, Mloc
       do j = 1, Nghost
          ik = (k-1)*Mloc+i
          sSmsg(ik,j) = PHI(i,Jbeg+j-1,k)
       enddo
       enddo
       enddo
       nreq = nreq +1
       call MPI_ISEND( sSmsg, len, MPI_SP, &
            n_suth, 1, comm2d, req(nreq), ier )
    endif

    if ( n_nrth .ne. MPI_PROC_NULL ) then
       nreq = nreq + 1
       call MPI_IRECV( rNmsg, len, MPI_SP, &
            n_nrth, 1, comm2d, req(nreq), ier )
       do k = 1, Kloc
       do i = 1, Mloc
       do j = 1, Nghost
          ik = (k-1)*Mloc+i
          sNmsg(ik,j) = PHI(i,Jend-j+1,k)
       enddo
       enddo
       enddo
       nreq = nreq + 1
       call MPI_ISEND( sNmsg, len, MPI_SP, &
            n_nrth, 0, comm2d, req(nreq), ier )
    endif

    call MPI_WAITALL( nreq, req, status, ier )

    if ( n_suth .ne. MPI_PROC_NULL ) then
       do k = 1, Kloc
       do i = 1, Mloc
       do j = 1, Nghost
          ik = (k-1)*Mloc+i
          PHI(i,Jbeg-j,k) = rSmsg(ik,j)
       enddo
       enddo
       enddo
    endif

    if ( n_nrth .ne. MPI_PROC_NULL ) then
       do k = 1, Kloc
       do i = 1, Mloc
       do j = 1, Nghost
          ik = (k-1)*Mloc+i
          PHI(i,Jend+j,k) = rNmsg(ik,j)
       enddo
       enddo
       enddo
    endif

    return
    END SUBROUTINE phi_3D_exch
# endif


# if defined(PARALLEL)
    ! Jeff added this subroutine to pass mask 02/14/2011
    SUBROUTINE phi_int_exch(PHI)
    USE GLOBAL
    IMPLICIT NONE
    INTEGER,INTENT(INOUT) :: PHI(Mloc,Nloc)

    INTEGER,DIMENSION(MPI_STATUS_SIZE,4) :: status
    INTEGER,DIMENSION(4) :: req
    INTEGER :: i,j,nreq,len
    INTEGER,DIMENSION(Mloc,Nghost) :: rNmsg, sNmsg,rSmsg,sSmsg
    INTEGER,DIMENSION(Nloc,Nghost) :: rWmsg, sWmsg,rEmsg,sEmsg

! for east-west

    len = Nloc * Nghost

    nreq = 0
    if ( n_west .ne. MPI_PROC_NULL ) then
       nreq = nreq + 1
       call MPI_IRECV( rWmsg, len, MPI_INTEGER, &
            n_west, 0, comm2d, req(nreq), ier )
       do j = 1, Nloc
       do i = 1, Nghost
          sWmsg(j,i) = PHI(Ibeg+i-1,j)
       enddo
       enddo
       nreq = nreq +1
       call MPI_ISEND( sWmsg, len, MPI_INTEGER, &
            n_west, 1, comm2d, req(nreq), ier )
    endif

    if ( n_east .ne. MPI_PROC_NULL ) then
       nreq = nreq + 1
       call MPI_IRECV( rEmsg, len, MPI_INTEGER, &
            n_east, 1, comm2d, req(nreq), ier )
       do j = 1, Nloc
       do i = 1, Nghost
          sEmsg(j,i) = PHI(Iend-i+1,j)
       enddo
       enddo
       nreq = nreq +1
       call MPI_ISEND( sEmsg, len, MPI_INTEGER, &
            n_east, 0, comm2d, req(nreq), ier )
    endif

    call MPI_WAITALL( nreq, req, status, ier )

    if ( n_west .ne. MPI_PROC_NULL ) then
       do j = 1, Nloc
       do i = 1, Nghost
          PHI(Ibeg-i,j) = rWmsg(j,i)
       enddo
       enddo
    endif

    if ( n_east .ne. MPI_PROC_NULL ) then
       do j = 1, Nloc
       do i = 1, Nghost
          PHI(Iend+i,j) = rEmsg(j,i)
       enddo
       enddo
    endif

! for nrth-suth

    len = Mloc * Nghost

    nreq = 0
    if ( n_suth .ne. MPI_PROC_NULL ) then
       nreq = nreq + 1
       call MPI_IRECV( rSmsg, len, MPI_INTEGER, &
            n_suth, 0, comm2d, req(nreq), ier )
       do i = 1, Mloc
       do j = 1, Nghost
          sSmsg(i,j) = PHI(i,Jbeg+j-1)
       enddo
       enddo
       nreq = nreq +1
       call MPI_ISEND( sSmsg, len, MPI_INTEGER, &
            n_suth, 1, comm2d, req(nreq), ier )
    endif

    if ( n_nrth .ne. MPI_PROC_NULL ) then
       nreq = nreq + 1
       call MPI_IRECV( rNmsg, len, MPI_INTEGER, &
            n_nrth, 1, comm2d, req(nreq), ier )
       do i = 1, Mloc
       do j = 1, Nghost
          sNmsg(i,j) = PHI(i,Jend-j+1)
       enddo
       enddo
       nreq = nreq + 1
       call MPI_ISEND( sNmsg, len, MPI_INTEGER, &
            n_nrth, 0, comm2d, req(nreq), ier )
    endif

    call MPI_WAITALL( nreq, req, status, ier )

    if ( n_suth .ne. MPI_PROC_NULL ) then
       do i = 1, Mloc
       do j = 1, Nghost
          PHI(i,Jbeg-j) = rSmsg(i,j)
       enddo
       enddo
    endif

    if ( n_nrth .ne. MPI_PROC_NULL ) then
       do i = 1, Mloc
       do j = 1, Nghost
          PHI(i,Jend+j) = rNmsg(i,j)
       enddo
       enddo
    endif
    END SUBROUTINE phi_int_exch
# endif


    subroutine adv_scalar_old(Phi,DPhi,R5,IVAR)
!--------------------------------------------------------
!   Subroutine for scalar convection and horizontal diffusion
!   IVAR: indication of different scalars
!    = 1: turbulent kinetic energy k
!    = 2: dissipation rate epsilon
!    = 3: salinity
!    = 4: temperature
!    = 5: bubble number density
!    = 6: sediment concentration
!   Last update: Gangfeng Ma, 31/08/2011
!-------------------------------------------------------
    use global
    implicit none
    integer, intent(in) :: IVAR
    real(SP), dimension(Mloc,Nloc,Kloc), intent(in) :: Phi,DPhi
    real(SP), dimension(Mloc,Nloc,Kloc), intent(inout) :: R5
    real(SP), dimension(:,:,:), allocatable :: DelxDPhi,DelyDPhi
    real(SP), dimension(:,:,:), allocatable :: DPhixL,DPhixR,DPhiyL,DPhiyR
    real(SP), dimension(:,:,:), allocatable :: ScalxL,ScalxR,Scalx,ScalyL,ScalyR,Scaly,  &
                                               Scalz,Sdiffx,Sdiffy,Sdiffz
    real(SP) :: omega_f,SchtH,SchtV,Sxx_L,Sxx_R,Sxx_U,Sxx_D,Syy_F,Syy_B,Syy_U,Syy_D
    integer :: i,j,k

    allocate(DelxDPhi(Mloc,Nloc,Kloc))
    allocate(DelyDPhi(Mloc,Nloc,Kloc))
    allocate(DPhixL(Mloc1,Nloc,Kloc))
    allocate(DPhixR(Mloc1,Nloc,Kloc))
    allocate(DPhiyL(Mloc,Nloc1,Kloc))
    allocate(DPhiyR(Mloc,Nloc1,Kloc))
    allocate(ScalxL(Mloc1,Nloc,Kloc))
    allocate(ScalxR(Mloc1,Nloc,Kloc))
    allocate(Scalx(Mloc1,Nloc,Kloc))
    allocate(ScalyL(Mloc,Nloc1,Kloc))
    allocate(ScalyR(Mloc,Nloc1,Kloc))
    allocate(Scaly(Mloc,Nloc1,Kloc))
    allocate(Scalz(Mloc,Nloc,Kloc1))
    allocate(Sdiffx(MLoc,Nloc,Kloc))
    allocate(Sdiffy(Mloc,Nloc,Kloc))
    allocate(Sdiffz(Mloc,Nloc,Kloc))

    call delxFun_3D(DPhi,DelxDPhi)
    call delyFun_3D(DPhi,DelyDPhi)

    call construct_3D_x(DPhi,DelxDPhi,DPhixL,DPhixR)
    call construct_3D_y(DPhi,DelyDPhi,DPhiyL,DPhiyR)

    ScalxL = DPhixL*UxL
    ScalxR = DPhixR*UxR
    ScalyL = DPhiyL*VyL
    ScalyR = DPhiyR*VyR

    call HLL(Mloc1,Nloc,Kloc,SxL,SxR,ScalxL,ScalxR,DPhixL,DPhixR,Scalx)
    call HLL(Mloc,Nloc1,Kloc,SyL,SyR,ScalyL,ScalyR,DPhiyL,DPhiyR,Scaly)

    Scalz = Zero
    do k = Kbeg+1,Kend
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      Omega_f = 0.5*(Omega(i,j,k-1)+Omega(i,j,k))
      Scalz(i,j,k) = 0.5*Omega_f*(Phi(i,j,k-1)+Phi(i,j,k))
    enddo
    enddo
    enddo

    ! at boundaries
    call flux_scalar_bc(Scalx,Scaly)

    ! Schmidt number
    if(IVAR==1) then  ! tke eq.
      SchtH = 1.0
      SchtV = 1.0
    elseif(IVAR==2) then  ! epsilon eq.
      SchtH = 1.0
      SchtV = 1.3
    else
      SchtH = 1.0
      SchtV = 1.0
    endif

    ! estimate diffusion
    Sdiffx = Zero; Sdiffy = Zero; Sdiffz = Zero
    do k = Kbeg,Kend
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      Sxx_L = 0.5*(CmuHt(i-1,j,k)+CmuHt(i,j,k))*((Phi(i,j,k)-Phi(i-1,j,k))/dx+  & 
           0.25*(Phi(i,j,k+1)+Phi(i-1,j,k+1)-Phi(i,j,k-1)-Phi(i-1,j,k-1))/(sigc(k+1)-sigc(k-1))*  &
           (DelxSc(i,j,k)+DelxSc(i-1,j,k)))
      Sxx_R = 0.5*(CmuHt(i,j,k)+CmuHt(i+1,j,k))*((Phi(i+1,j,k)-Phi(i,j,k))/dx+  & 
           0.25*(Phi(i+1,j,k+1)+Phi(i,j,k+1)-Phi(i+1,j,k-1)-Phi(i,j,k-1))/(sigc(k+1)-sigc(k-1))*  &
           (DelxSc(i+1,j,k)+DelxSc(i,j,k)))
      Sxx_U = 0.5*(CmuHt(i,j,k+1)+CmuHt(i,j,k))*(0.5*(Phi(i+1,j,k+1)+Phi(i+1,j,k)-Phi(i-1,j,k+1)-Phi(i-1,j,k))/(2.0*dx)+  & 
           (Phi(i,j,k+1)-Phi(i,j,k))/(sigc(k+1)-sigc(k))*0.5*(DelxSc(i,j,k+1)+DelxSc(i,j,k)))                                                           
      Sxx_D = 0.5*(CmuHt(i,j,k)+CmuHt(i,j,k-1))*(0.5*(Phi(i+1,j,k)+Phi(i+1,j,k-1)-Phi(i-1,j,k)-Phi(i-1,j,k-1))/(2.0*dx)+  &
           (Phi(i,j,k)-Phi(i,j,k-1))/(sigc(k)-sigc(k-1))*0.5*(DelxSc(i,j,k)+DelxSc(i,j,k-1)))
      Sdiffx(i,j,k) = 0.5*((D(i+1,j)+D(i,j))*Sxx_R-(D(i,j)+D(i-1,j))*Sxx_L)/dx+  &                                                                 
           D(i,j)*DelxSc(i,j,k)*(Sxx_U-Sxx_D)/dsig(k)

      Syy_F = 0.5*(CmuHt(i,j-1,k)+CmuHt(i,j,k))*((Phi(i,j,k)-Phi(i,j-1,k))/dy+  & 
           0.25*(Phi(i,j,k+1)+Phi(i,j-1,k+1)-Phi(i,j,k-1)-Phi(i,j-1,k-1))/(sigc(k+1)-sigc(k-1))*  & 
           (DelySc(i,j,k)+DelySc(i,j-1,k)))
      Syy_B = 0.5*(CmuHt(i,j,k)+CmuHt(i,j+1,k))*((Phi(i,j+1,k)-Phi(i,j,k))/dy+  & 
           0.25*(Phi(i,j+1,k+1)+Phi(i,j,k+1)-Phi(i,j+1,k-1)-Phi(i,j,k-1))/(sigc(k+1)-sigc(k-1))*  &
           (DelySc(i,j+1,k)+DelySc(i,j,k)))
      Syy_U = 0.5*(CmuHt(i,j,k+1)+CmuHt(i,j,k))*(0.5*(Phi(i,j+1,k+1)+Phi(i,j+1,k)-Phi(i,j-1,k+1)-Phi(i,j-1,k))/(2.0*dy)+  & 
           (Phi(i,j,k+1)-Phi(i,j,k))/(sigc(k+1)-sigc(k))*0.5*(DelySc(i,j,k+1)+DelySc(i,j,k)))                                                           
      Syy_D = 0.5*(CmuHt(i,j,k)+CmuHt(i,j,k-1))*(0.5*(Phi(i,j+1,k)+Phi(i,j+1,k-1)-Phi(i,j-1,k)-Phi(i,j-1,k-1))/(2.0*dy)+  & 
           (Phi(i,j,k)-Phi(i,j,k-1))/(sigc(k)-sigc(k-1))*0.5*(DelySc(i,j,k)+DelySc(i,j,k-1)))
      Sdiffy(i,j,k) = 0.5*((D(i,j+1)+D(i,j))*Syy_B-(D(i,j)+D(i,j-1))*Syy_F)/dy+  &                                                                 
           D(i,j)*DelySc(i,j,k)*(Syy_U-Syy_D)/dsig(k)

      Sdiffz(i,j,k) = 0.5*((CmuVt(i,j,k+1)+CmuVt(i,j,k))*(Phi(i,j,k+1)-Phi(i,j,k))/  &
              (sigc(k+1)-sigc(k))-(CmuVt(i,j,k)+CmuVt(i,j,k-1))*(Phi(i,j,k)-Phi(i,j,k-1))/ &
              (sigc(k)-sigc(k-1)))/SchtV/dsig(k)/D(i,j)
    enddo
    enddo
    enddo

    R5 = Zero
    do k = Kbeg,Kend
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      R5(i,j,k) = -1.0/dx*(Scalx(i+1,j,k)-Scalx(i,j,k))-1.0/dy*(Scaly(i,j+1,k)-Scaly(i,j,k))  &
                  -1.0/dsig(k)*(Scalz(i,j,k+1)-Scalz(i,j,k))+Sdiffx(i,j,k)+Sdiffy(i,j,k)+Sdiffz(i,j,k)
    enddo
    enddo
    enddo

    deallocate(DelxDPhi)
    deallocate(DelyDPhi)
    deallocate(DPhixL)
    deallocate(DPhixR)
    deallocate(DPhiyL)
    deallocate(DPhiyR)
    deallocate(ScalxL)
    deallocate(ScalxR)
    deallocate(Scalx)
    deallocate(ScalyL)
    deallocate(ScalyR)
    deallocate(Scaly)
    deallocate(Scalz)
    deallocate(Sdiffx)
    deallocate(Sdiffy)
    deallocate(Sdiffz)
     
    return
    end subroutine adv_scalar_old


    subroutine adv_scalar_tvd(Phi,R5,IVAR)
!--------------------------------------------------------
!   Subroutine for scalar convection and horizontal diffusion  
!   IVAR: indication of different scalars 
!    = 1: turbulent kinetic energy k
!    = 2: dissipation rate epsilon
!    = 3: salinity 
!    = 4: temperature
!    = 5: bubble number density 
!    = 6: sediment concentration 
!   Last update: Gangfeng Ma, 26/11/2011
!-------------------------------------------------------  
    use global
    implicit none
    integer, intent(in) :: IVAR
    integer, parameter :: TVD = 3
    real(SP), dimension(Mloc,Nloc,Kloc), intent(in) :: Phi
    real(SP), dimension(Mloc,Nloc,Kloc), intent(inout) :: R5
    real(SP), dimension(:,:,:), allocatable :: Scalx,Scaly,Scalz,Sdiffx,Sdiffy,Sdiffz 
    real(SP) :: Uface,Vface,Wface,Cp,Cm,Rp,Rm,Phifp,Phifm,DUf,DVf,SchtH,SchtV,Psi
    integer :: i,j,k

    allocate(Scalx(Mloc1,Nloc,Kloc))
    allocate(Scaly(Mloc,Nloc1,Kloc))
    allocate(Scalz(Mloc,Nloc,Kloc1))
    allocate(Sdiffx(MLoc,Nloc,Kloc))
    allocate(Sdiffy(Mloc,Nloc,Kloc))
    allocate(Sdiffz(Mloc,Nloc,Kloc))

    ! advection in x direction
    Scalx = Zero
    do k = Kbeg,Kend
    do j = Jbeg,Jend
    do i = Ibeg,Iend+1
      ! compute the courant number for TVD schemes
      Uface = 0.5*(U(i-1,j,k)+U(i,j,k))
      Cp = 0.5*(Uface+abs(Uface))*dt/dx
      Cm = 0.5*(Uface-abs(Uface))*dt/dx
     
      ! compute the upwind gradient ratio for TVD schemes
      Rp = (Phi(i-1,j,k)-Phi(i-2,j,k)+small)/(Phi(i,j,k)-Phi(i-1,j,k)+small)
      Rm = (Phi(i,j,k)-Phi(i+1,j,k)+small)/(Phi(i-1,j,k)-Phi(i,j,k)+small)

      Phifp = Phi(i-1,j,k)+0.5*Psi(Rp,TVD)*(1.0-Cp)*(Phi(i,j,k)-Phi(i-1,j,k))
      Phifm = Phi(i,j,k)-0.5*Psi(Rm,TVD)*(1.0+Cm)*(Phi(i,j,k)-Phi(i-1,j,k))

      DUf = 0.5*(DU(i-1,j,k)+DU(i,j,k))
      Scalx(i,j,k) = 0.5*(DUf+abs(DUf))*Phifp+0.5*(DUf-abs(DUf))*Phifm
    enddo
    enddo
    enddo

    ! advection in y direction
    Scaly = Zero
    do k = Kbeg,Kend
    do j = Jbeg,Jend+1
    do i = Ibeg,Iend
      ! compute the courant number for TVD schemes 
      Vface = 0.5*(V(i-1,j,k)+V(i,j,k))
      Cp = 0.5*(Vface+abs(Vface))*dt/dy
      Cm = 0.5*(Vface-abs(Vface))*dt/dy                                                                                                                                          
      ! compute the upwind gradient ratio for TVD schemes
      Rp = (Phi(i,j-1,k)-Phi(i,j-2,k)+small)/(Phi(i,j,k)-Phi(i,j-1,k)+small)
      Rm = (Phi(i,j,k)-Phi(i,j+1,k)+small)/(Phi(i,j-1,k)-Phi(i,j,k)+small)

      Phifp = Phi(i,j-1,k)+0.5*Psi(Rp,TVD)*(1.0-Cp)*(Phi(i,j,k)-Phi(i,j-1,k))
      Phifm = Phi(i,j,k)-0.5*Psi(Rm,TVD)*(1.0+Cm)*(Phi(i,j,k)-Phi(i,j-1,k)) 

      DVf = 0.5*(DV(i,j-1,k)+DV(i,j,k))
      Scaly(i,j,k) = 0.5*(DVf+abs(DVf))*Phifp+0.5*(DVf-abs(DVf))*Phifm
    enddo
    enddo
    enddo

    ! advection in z direction
    Scalz = Zero
    do k = Kbeg+1,Kend
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      ! compute the courant number for TVD schemes 
      Wface = 0.5*(Omega(i,j,k-1)+Omega(i,j,k))
      Cp = 0.5*(Wface+abs(Wface))*dt/(sigc(k)-sigc(k-1))
      Cm = 0.5*(Wface-abs(Wface))*dt/(sigc(k)-sigc(k-1))                                                                                                                                          
      ! compute the upwind gradient ratio for TVD schemes
      Rp = (Phi(i,j,k-1)-Phi(i,j,k-2)+small)/(Phi(i,j,k)-Phi(i,j,k-1)+small) 
      Rm = (Phi(i,j,k)-Phi(i,j,k+1)+small)/(Phi(i,j,k-1)-Phi(i,j,k)+small) 

      Phifp = Phi(i,j,k-1)+0.5*Psi(Rp,TVD)*(1.0-Cp)*(Phi(i,j,k)-Phi(i,j,k-1))
      Phifm = Phi(i,j,k)-0.5*Psi(Rm,TVD)*(1.0+Cm)*(Phi(i,j,k)-Phi(i,j,k-1))

      Scalz(i,j,k) = 0.5*(Wface+abs(Wface))*Phifp+0.5*(Wface-abs(Wface))*Phifm
    enddo
    enddo
    enddo

    ! at boundaries
    call flux_scalar_bc(Scalx,Scaly)

    ! Schmidt number
    if(IVAR==1) then  ! tke eq.
      SchtH = 1.0
      SchtV = 1.0
    elseif(IVAR==2) then  ! epsilon eq.
      SchtH = 1.0
      SchtV = 1.3
    elseif(IVAR==5) then ! bubble
      SchtH = 0.7
      SchtV = 0.7
    else
      SchtH = 1.0
      SchtV = 1.0
    endif

    ! estimate diffusion
    Sdiffx = Zero; Sdiffy = Zero; Sdiffz = Zero
    do k = Kbeg,Kend
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      Sdiffx(i,j,k) = 0.25*((CmuHt(i+1,j,k)+CmuHt(i,j,k))*(D(i+1,j)+D(i,j))*(Phi(i+1,j,k)-Phi(i,j,k))-  &
           (CmuHt(i,j,k)+CmuHt(i-1,j,k))*(D(i,j)+D(i-1,j))*(Phi(i,j,k)-Phi(i-1,j,k)))/SchtH/dx**2
      Sdiffy(i,j,k) = 0.25*((CmuHt(i,j+1,k)+CmuHt(i,j,k))*(D(i,j+1)+D(i,j))*(Phi(i,j+1,k)-Phi(i,j,k))-  & 
           (CmuHt(i,j,k)+CmuHt(i,j-1,k))*(D(i,j)+D(i,j-1))*(Phi(i,j,k)-Phi(i,j-1,k)))/SchtH/dy**2
      Sdiffz(i,j,k) = 0.5*((CmuVt(i,j,k+1)+CmuVt(i,j,k))*(Phi(i,j,k+1)-Phi(i,j,k))/  &   
           (sigc(k+1)-sigc(k))-(CmuVt(i,j,k)+CmuVt(i,j,k-1))*(Phi(i,j,k)-Phi(i,j,k-1))/ &  
           (sigc(k)-sigc(k-1)))/SchtV/dsig(k)/D(i,j)
    enddo
    enddo
    enddo

    R5 = Zero
    do k = Kbeg,Kend
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      R5(i,j,k) = -1.0/dx*(Scalx(i+1,j,k)-Scalx(i,j,k))-1.0/dy*(Scaly(i,j+1,k)-Scaly(i,j,k))  &
                  -1.0/dsig(k)*(Scalz(i,j,k+1)-Scalz(i,j,k))+Sdiffx(i,j,k)+Sdiffy(i,j,k)+Sdiffz(i,j,k)
    enddo
    enddo
    enddo

    deallocate(Scalx)
    deallocate(Scaly)
    deallocate(Scalz)
    deallocate(Sdiffx)
    deallocate(Sdiffy)
    deallocate(Sdiffz)

    end subroutine adv_scalar_tvd

    function Psi(r,TVD)
    use global, only: SP,Zero,One,Small
    real(SP), intent(in) :: r
    integer, intent(in) :: TVD
    real(SP) :: Psi

    if(TVD==1) then
      ! upwind
      Psi = Zero
    elseif(TVD==2) then
      ! Lax-Wendroff
      Psi = One
    elseif(TVD==3) then
      ! Superbee
      Psi = max(Zero,max(min(2*r,One),min(r,2.0)))
    elseif(TVD==4) then
      ! Van Leer
      Psi = (r+abs(r))/(1+r+small)
    else
      Psi = Zero
    endif

    return
    end function Psi

    subroutine adv_scalar_upwind(Phi,R5,IVAR)
!--------------------------------------------------------
!   Subroutine for scalar convection and horizontal diffusion
!   IVAR: indication of different scalars
!    = 1: turbulent kinetic energy k
!    = 2: dissipation rate epsilon
!    = 3: salinity
!    = 4: temperature
!    = 5: bubble number density
!    = 6: sediment concentration
!   Last update: Gangfeng Ma, 31/08/2011
!-------------------------------------------------------
    use global
    implicit none
    integer, intent(in) :: IVAR
    real(SP), dimension(Mloc,Nloc,Kloc), intent(in) :: Phi
    real(SP), dimension(Mloc,Nloc,Kloc), intent(inout) :: R5
    real(SP), dimension(:,:,:), allocatable :: DelxPhi,DelyPhi,DelzPhi
    real(SP), dimension(:,:,:), allocatable :: PhixL,PhixR,PhiyL,PhiyR,PhizL,PhizR
    real(SP), dimension(:,:,:), allocatable :: Scalx,Scaly,Scalz,Sdiffx,Sdiffy,Sdiffz
    real(SP) :: DUxf,DVyf,omega_f,SchtH,SchtV,Sxx_L,Sxx_R,Sxx_U,Sxx_D,Syy_F,Syy_B,Syy_U,Syy_D
    integer :: i,j,k

    allocate(DelxPhi(Mloc,Nloc,Kloc))
    allocate(DelyPhi(Mloc,Nloc,Kloc))
    allocate(DelzPhi(Mloc,Nloc,Kloc))
    allocate(PhixL(Mloc1,Nloc,Kloc))
    allocate(PhixR(Mloc1,Nloc,Kloc))
    allocate(PhiyL(Mloc,Nloc1,Kloc))
    allocate(PhiyR(Mloc,Nloc1,Kloc))
    allocate(PhizL(Mloc,Nloc,Kloc1))
    allocate(PhizR(Mloc,Nloc,Kloc1))
    allocate(Scalx(Mloc1,Nloc,Kloc))
    allocate(Scaly(Mloc,Nloc1,Kloc))
    allocate(Scalz(Mloc,Nloc,Kloc1))
    allocate(Sdiffx(MLoc,Nloc,Kloc))
    allocate(Sdiffy(Mloc,Nloc,Kloc))
    allocate(Sdiffz(Mloc,Nloc,Kloc))

    call delxFun_3D(Phi,DelxPhi)
    call delyFun_3D(Phi,DelyPhi)
    call delzFun_3D(Phi,DelzPhi)

    call construct_3D_x(Phi,DelxPhi,PhixL,PhixR)
    call construct_3D_y(Phi,DelyPhi,PhiyL,PhiyR)
    call construct_3D_z(Phi,DelzPhi,PhizL,PhizR)

    Scalx = Zero
    do k = Kbeg,Kend
    do j = Jbeg,Jend
    do i = Ibeg,Iend+1
      DUxf = 0.5*(DUxL(i,j,k)+DUxR(i,j,k))
      Scalx(i,j,k) = 0.5*(DUxf*(PhixL(i,j,k)+PhixR(i,j,k))-abs(DUxf)*(PhixR(i,j,k)-PhixL(i,j,k)))
    enddo
    enddo
    enddo

    Scaly = Zero
    do k = Kbeg,Kend
    do j = Jbeg,Jend+1
    do i = Ibeg,Iend
      DVyf = 0.5*(DVyL(i,j,k)+DVyR(i,j,k))
      Scaly(i,j,k) = 0.5*(DVyf*(PhiyL(i,j,k)+PhiyR(i,j,k))-abs(DVyf)*(PhiyR(i,j,k)-PhiyL(i,j,k))) 
    enddo
    enddo
    enddo

    Scalz = Zero
    do k = Kbeg+1,Kend
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      Omega_f = 0.5*(OzL(i,j,k)+OzR(i,j,k))
      Scalz(i,j,k) = 0.5*(Omega_f*(PhizL(i,j,k)+PhizR(i,j,k))-abs(Omega_f)*(PhizR(i,j,k)-PhizL(i,j,k)))
    enddo
    enddo
    enddo

    ! at boundaries
    call flux_scalar_bc(Scalx,Scaly)

    ! Schmidt number
    if(IVAR==1) then  ! tke eq.
      SchtH = 1.0
      SchtV = 1.0
    elseif(IVAR==2) then  ! epsilon eq.
      SchtH = 1.0
      SchtV = 1.3
    elseif(IVAR==5) then ! bubble
      SchtH = 0.7
      SchtV = 0.7
    else
      SchtH = 1.0
      SchtV = 1.0
    endif

    ! estimate diffusion
    Sdiffx = Zero; Sdiffy = Zero; Sdiffz = Zero
    do k = Kbeg,Kend
    do j = Jbeg,Jend
    do i = Ibeg,Iend
!      large errors are found when the free surface becomes steep using the following treatment
!      Sxx_L = 0.5*(CmuHt(i-1,j,k)+CmuHt(i,j,k))*((Phi(i,j,k)-Phi(i-1,j,k))/dx+  & 
!           0.25*(Phi(i,j,k+1)+Phi(i-1,j,k+1)-Phi(i,j,k-1)-Phi(i-1,j,k-1))/(sigc(k+1)-sigc(k-1))*  & 
!           (DelxSc(i,j,k)+DelxSc(i-1,j,k)))
!      Sxx_R = 0.5*(CmuHt(i,j,k)+CmuHt(i+1,j,k))*((Phi(i+1,j,k)-Phi(i,j,k))/dx+  &   
!           0.25*(Phi(i+1,j,k+1)+Phi(i,j,k+1)-Phi(i+1,j,k-1)-Phi(i,j,k-1))/(sigc(k+1)-sigc(k-1))*  &
!           (DelxSc(i+1,j,k)+DelxSc(i,j,k)))
!      Sxx_U = 0.5*(CmuHt(i,j,k+1)+CmuHt(i,j,k))*(0.5*(Phi(i+1,j,k+1)+Phi(i+1,j,k)-Phi(i-1,j,k+1)-Phi(i-1,j,k))/(2.0*dx)+  & 
!           (Phi(i,j,k+1)-Phi(i,j,k))/(sigc(k+1)-sigc(k))*0.5*(DelxSc(i,j,k+1)+DelxSc(i,j,k)))
!      Sxx_D = 0.5*(CmuHt(i,j,k)+CmuHt(i,j,k-1))*(0.5*(Phi(i+1,j,k)+Phi(i+1,j,k-1)-Phi(i-1,j,k)-Phi(i-1,j,k-1))/(2.0*dx)+  &  
!           (Phi(i,j,k)-Phi(i,j,k-1))/(sigc(k)-sigc(k-1))*0.5*(DelxSc(i,j,k)+DelxSc(i,j,k-1)))  
!      Sdiffx(i,j,k) = 0.5*((D(i+1,j)+D(i,j))*Sxx_R-(D(i,j)+D(i-1,j))*Sxx_L)/dx+  & 
!           D(i,j)*DelxSc(i,j,k)*(Sxx_U-Sxx_D)/dsig(k)
!
!      Syy_F = 0.5*(CmuHt(i,j-1,k)+CmuHt(i,j,k))*((Phi(i,j,k)-Phi(i,j-1,k))/dy+  &  
!           0.25*(Phi(i,j,k+1)+Phi(i,j-1,k+1)-Phi(i,j,k-1)-Phi(i,j-1,k-1))/(sigc(k+1)-sigc(k-1))*  & 
!           (DelySc(i,j,k)+DelySc(i,j-1,k)))
!      Syy_B = 0.5*(CmuHt(i,j,k)+CmuHt(i,j+1,k))*((Phi(i,j+1,k)-Phi(i,j,k))/dy+  & 
!           0.25*(Phi(i,j+1,k+1)+Phi(i,j,k+1)-Phi(i,j+1,k-1)-Phi(i,j,k-1))/(sigc(k+1)-sigc(k-1))*  & 
!           (DelySc(i,j+1,k)+DelySc(i,j,k)))
!      Syy_U = 0.5*(CmuHt(i,j,k+1)+CmuHt(i,j,k))*(0.5*(Phi(i,j+1,k+1)+Phi(i,j+1,k)-Phi(i,j-1,k+1)-Phi(i,j-1,k))/(2.0*dy)+  & 
!           (Phi(i,j,k+1)-Phi(i,j,k))/(sigc(k+1)-sigc(k))*0.5*(DelySc(i,j,k+1)+DelySc(i,j,k)))  
!      Syy_D = 0.5*(CmuHt(i,j,k)+CmuHt(i,j,k-1))*(0.5*(Phi(i,j+1,k)+Phi(i,j+1,k-1)-Phi(i,j-1,k)-Phi(i,j-1,k-1))/(2.0*dy)+  & 
!           (Phi(i,j,k)-Phi(i,j,k-1))/(sigc(k)-sigc(k-1))*0.5*(DelySc(i,j,k)+DelySc(i,j,k-1))) 
!      Sdiffy(i,j,k) = 0.5*((D(i,j+1)+D(i,j))*Syy_B-(D(i,j)+D(i,j-1))*Syy_F)/dy+  & 
!           D(i,j)*DelySc(i,j,k)*(Syy_U-Syy_D)/dsig(k)

      ! use the simplified form following FVCOM
      Sdiffx(i,j,k) = 0.25*((CmuHt(i+1,j,k)+CmuHt(i,j,k))*(D(i+1,j)+D(i,j))*(Phi(i+1,j,k)-Phi(i,j,k))-  &
           (CmuHt(i,j,k)+CmuHt(i-1,j,k))*(D(i,j)+D(i-1,j))*(Phi(i,j,k)-Phi(i-1,j,k)))/SchtH/dx**2
      Sdiffy(i,j,k) = 0.25*((CmuHt(i,j+1,k)+CmuHt(i,j,k))*(D(i,j+1)+D(i,j))*(Phi(i,j+1,k)-Phi(i,j,k))-  & 
           (CmuHt(i,j,k)+CmuHt(i,j-1,k))*(D(i,j)+D(i,j-1))*(Phi(i,j,k)-Phi(i,j-1,k)))/SchtH/dy**2
      Sdiffz(i,j,k) = 0.5*((CmuVt(i,j,k+1)+CmuVt(i,j,k))*(Phi(i,j,k+1)-Phi(i,j,k))/  &   
           (sigc(k+1)-sigc(k))-(CmuVt(i,j,k)+CmuVt(i,j,k-1))*(Phi(i,j,k)-Phi(i,j,k-1))/ &  
           (sigc(k)-sigc(k-1)))/SchtV/dsig(k)/D(i,j)
    enddo
    enddo
    enddo

    R5 = Zero
    do k = Kbeg,Kend
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      R5(i,j,k) = -1.0/dx*(Scalx(i+1,j,k)-Scalx(i,j,k))-1.0/dy*(Scaly(i,j+1,k)-Scaly(i,j,k))  &
                  -1.0/dsig(k)*(Scalz(i,j,k+1)-Scalz(i,j,k))+Sdiffx(i,j,k)+Sdiffy(i,j,k)+Sdiffz(i,j,k)
    enddo
    enddo
    enddo

    deallocate(DelxPhi)
    deallocate(DelyPhi)
    deallocate(DelzPhi)
    deallocate(PhixL)
    deallocate(PhixR)
    deallocate(PhiyL)
    deallocate(PhiyR)
    deallocate(PhizL)
    deallocate(PhizR)
    deallocate(Scalx)
    deallocate(Scaly)
    deallocate(Scalz)
    deallocate(Sdiffx)
    deallocate(Sdiffy)
    deallocate(Sdiffz)
     
    return
    end subroutine adv_scalar_upwind


    subroutine adv_scalar(Phi,R5,IVAR)
!--------------------------------------------------------
!   Subroutine for scalar convection and horizontal diffusion
!   IVAR: indication of different scalars
!    = 1: turbulent kinetic energy k
!    = 2: dissipation rate epsilon
!    = 3: salinity
!    = 4: temperature
!    = 5: bubble number density
!    = 6: sediment concentration
!   Last update: Gangfeng Ma, 31/08/2011
!-------------------------------------------------------
    use global
    implicit none
    integer, intent(in) :: IVAR
    real(SP), dimension(Mloc,Nloc,Kloc), intent(in) :: Phi
    real(SP), dimension(Mloc,Nloc,Kloc), intent(inout) :: R5
    real(SP), dimension(:,:,:), allocatable :: DelxPhi,DelyPhi,DelzPhi
    real(SP), dimension(:,:,:), allocatable :: PhixL,PhixR,PhiyL,PhiyR,PhizL,PhizR
    real(SP), dimension(:,:,:), allocatable :: Scalx,Scaly,Scalz,Sdiffx,Sdiffy,Sdiffz
    real(SP) :: DUxf,DVyf,omega_f,SchtH,SchtV,Sxx_L,Sxx_R,Sxx_U,Sxx_D,Syy_F,Syy_B,Syy_U,Syy_D
    integer :: i,j,k

    allocate(DelxPhi(Mloc,Nloc,Kloc))
    allocate(DelyPhi(Mloc,Nloc,Kloc))
    allocate(DelzPhi(Mloc,Nloc,Kloc))
    allocate(PhixL(Mloc1,Nloc,Kloc))
    allocate(PhixR(Mloc1,Nloc,Kloc))
    allocate(PhiyL(Mloc,Nloc1,Kloc))
    allocate(PhiyR(Mloc,Nloc1,Kloc))
    allocate(PhizL(Mloc,Nloc,Kloc1))
    allocate(PhizR(Mloc,Nloc,Kloc1))
    allocate(Scalx(Mloc1,Nloc,Kloc))
    allocate(Scaly(Mloc,Nloc1,Kloc))
    allocate(Scalz(Mloc,Nloc,Kloc1))
    allocate(Sdiffx(MLoc,Nloc,Kloc))
    allocate(Sdiffy(Mloc,Nloc,Kloc))
    allocate(Sdiffz(Mloc,Nloc,Kloc))

    call delxFun_3D(Phi,DelxPhi)
    call delyFun_3D(Phi,DelyPhi)
    call delzFun_3D(Phi,DelzPhi)

    call construct_3D_x(Phi,DelxPhi,PhixL,PhixR)
    call construct_3D_y(Phi,DelyPhi,PhiyL,PhiyR)
    call construct_3D_z(Phi,DelzPhi,PhizL,PhizR)

    Scalx = Zero
    do k = Kbeg,Kend
    do j = Jbeg,Jend
    do i = Ibeg,Iend+1
      DUxf = 0.5*(DUxL(i,j,k)+DUxR(i,j,k))
      Scalx(i,j,k) = 0.5*(DUxf*(PhixL(i,j,k)+PhixR(i,j,k))-abs(DUxf)*(PhixR(i,j,k)-PhixL(i,j,k)))
    enddo
    enddo
    enddo

    Scaly = Zero
    do k = Kbeg,Kend
    do j = Jbeg,Jend+1
    do i = Ibeg,Iend
      DVyf = 0.5*(DVyL(i,j,k)+DVyR(i,j,k))
      Scaly(i,j,k) = 0.5*(DVyf*(PhiyL(i,j,k)+PhiyR(i,j,k))-abs(DVyf)*(PhiyR(i,j,k)-PhiyL(i,j,k))) 
    enddo
    enddo
    enddo

    Scalz = Zero
    do k = Kbeg+1,Kend
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      Omega_f = 0.5*(OzL(i,j,k)+OzR(i,j,k))
      Scalz(i,j,k) = 0.5*(Omega_f*(PhizL(i,j,k)+PhizR(i,j,k))-abs(Omega_f)*(PhizR(i,j,k)-PhizL(i,j,k)))
    enddo
    enddo
    enddo

    ! at boundaries
    call flux_scalar_bc(Scalx,Scaly)

    ! Schmidt number
    if(IVAR==1) then  ! tke eq.
      SchtH = 1.0
      SchtV = 1.0
    elseif(IVAR==2) then  ! epsilon eq.
      SchtH = 1.0
      SchtV = 1.3
    else
      SchtH = 1.0
      SchtV = 1.0
    endif

    ! estimate diffusion
    Sdiffx = Zero; Sdiffy = Zero; Sdiffz = Zero
    do k = Kbeg,Kend
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      Sxx_L = 0.5*(CmuHt(i-1,j,k)+CmuHt(i,j,k))*((Phi(i,j,k)-Phi(i-1,j,k))/dx+  & 
           0.25*(Phi(i,j,k+1)+Phi(i-1,j,k+1)-Phi(i,j,k-1)-Phi(i-1,j,k-1))/(sigc(k+1)-sigc(k-1))*  & 
           (DelxSc(i,j,k)+DelxSc(i-1,j,k)))
      Sxx_R = 0.5*(CmuHt(i,j,k)+CmuHt(i+1,j,k))*((Phi(i+1,j,k)-Phi(i,j,k))/dx+  &   
           0.25*(Phi(i+1,j,k+1)+Phi(i,j,k+1)-Phi(i+1,j,k-1)-Phi(i,j,k-1))/(sigc(k+1)-sigc(k-1))*  &
           (DelxSc(i+1,j,k)+DelxSc(i,j,k)))
      Sxx_U = 0.5*(CmuHt(i,j,k+1)+CmuHt(i,j,k))*(0.5*(Phi(i+1,j,k+1)+Phi(i+1,j,k)-Phi(i-1,j,k+1)-Phi(i-1,j,k))/(2.0*dx)+  & 
           (Phi(i,j,k+1)-Phi(i,j,k))/(sigc(k+1)-sigc(k))*0.5*(DelxSc(i,j,k+1)+DelxSc(i,j,k)))
      Sxx_D = 0.5*(CmuHt(i,j,k)+CmuHt(i,j,k-1))*(0.5*(Phi(i+1,j,k)+Phi(i+1,j,k-1)-Phi(i-1,j,k)-Phi(i-1,j,k-1))/(2.0*dx)+  &  
           (Phi(i,j,k)-Phi(i,j,k-1))/(sigc(k)-sigc(k-1))*0.5*(DelxSc(i,j,k)+DelxSc(i,j,k-1)))  
      Sdiffx(i,j,k) = 0.5*((D(i+1,j)+D(i,j))*Sxx_R-(D(i,j)+D(i-1,j))*Sxx_L)/dx+  & 
           D(i,j)*DelxSc(i,j,k)*(Sxx_U-Sxx_D)/dsig(k)

      Syy_F = 0.5*(CmuHt(i,j-1,k)+CmuHt(i,j,k))*((Phi(i,j,k)-Phi(i,j-1,k))/dy+  &  
           0.25*(Phi(i,j,k+1)+Phi(i,j-1,k+1)-Phi(i,j,k-1)-Phi(i,j-1,k-1))/(sigc(k+1)-sigc(k-1))*  & 
           (DelySc(i,j,k)+DelySc(i,j-1,k)))
      Syy_B = 0.5*(CmuHt(i,j,k)+CmuHt(i,j+1,k))*((Phi(i,j+1,k)-Phi(i,j,k))/dy+  & 
           0.25*(Phi(i,j+1,k+1)+Phi(i,j,k+1)-Phi(i,j+1,k-1)-Phi(i,j,k-1))/(sigc(k+1)-sigc(k-1))*  & 
           (DelySc(i,j+1,k)+DelySc(i,j,k)))
      Syy_U = 0.5*(CmuHt(i,j,k+1)+CmuHt(i,j,k))*(0.5*(Phi(i,j+1,k+1)+Phi(i,j+1,k)-Phi(i,j-1,k+1)-Phi(i,j-1,k))/(2.0*dy)+  & 
           (Phi(i,j,k+1)-Phi(i,j,k))/(sigc(k+1)-sigc(k))*0.5*(DelySc(i,j,k+1)+DelySc(i,j,k)))  
      Syy_D = 0.5*(CmuHt(i,j,k)+CmuHt(i,j,k-1))*(0.5*(Phi(i,j+1,k)+Phi(i,j+1,k-1)-Phi(i,j-1,k)-Phi(i,j-1,k-1))/(2.0*dy)+  & 
           (Phi(i,j,k)-Phi(i,j,k-1))/(sigc(k)-sigc(k-1))*0.5*(DelySc(i,j,k)+DelySc(i,j,k-1))) 
      Sdiffy(i,j,k) = 0.5*((D(i,j+1)+D(i,j))*Syy_B-(D(i,j)+D(i,j-1))*Syy_F)/dy+  & 
           D(i,j)*DelySc(i,j,k)*(Syy_U-Syy_D)/dsig(k)

      Sdiffz(i,j,k) = 0.5*((CmuVt(i,j,k+1)+CmuVt(i,j,k))*(Phi(i,j,k+1)-Phi(i,j,k))/  &   
           (sigc(k+1)-sigc(k))-(CmuVt(i,j,k)+CmuVt(i,j,k-1))*(Phi(i,j,k)-Phi(i,j,k-1))/ &  
           (sigc(k)-sigc(k-1)))/SchtV/dsig(k)/D(i,j)
    enddo
    enddo
    enddo

    R5 = Zero
    do k = Kbeg,Kend
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      R5(i,j,k) = -1.0/dx*(Scalx(i+1,j,k)-Scalx(i,j,k))-1.0/dy*(Scaly(i,j+1,k)-Scaly(i,j,k))  &
                  -1.0/dsig(k)*(Scalz(i,j,k+1)-Scalz(i,j,k))+Sdiffx(i,j,k)+Sdiffy(i,j,k)+Sdiffz(i,j,k)
    enddo
    enddo
    enddo

    deallocate(DelxPhi)
    deallocate(DelyPhi)
    deallocate(DelzPhi)
    deallocate(PhixL)
    deallocate(PhixR)
    deallocate(PhiyL)
    deallocate(PhiyR)
    deallocate(PhizL)
    deallocate(PhizR)
    deallocate(Scalx)
    deallocate(Scaly)
    deallocate(Scalz)
    deallocate(Sdiffx)
    deallocate(Sdiffy)
    deallocate(Sdiffz)
     
    return
    end subroutine adv_scalar


    subroutine flux_scalar_bc(Scalx,Scaly)
!--------------------------------------------------------
!   Specify boundary conditions for scalar convection
!   Last update: Gangfeng Ma, 09/02/2011
!   fyshi add right boundary condition 12/16/2011
!-------------------------------------------------------
    use global
    implicit none
    real(SP), dimension(Mloc1,Nloc,Kloc), intent(inout) :: Scalx
    real(SP), dimension(Mloc,Nloc1,Kloc), intent(inout) :: Scaly
    real(SP), dimension(Nloc,Kloc) :: Scal_X0,Scal_Xn
    integer :: i,j,k

    ! temporarily set it here, why
    Scal_X0 = Zero
    Scal_Xn = Zero

    ! left and right side
# if defined (PARALLEL)
     if(n_west.eq.MPI_PROC_NULL) then
# endif
     do j = Jbeg,Jend
     do k = Kbeg,Kend
       if(Bc_X0==1.or.Bc_X0==2) then
         Scalx(Ibeg,j,k) = Zero
       elseif(Bc_X0==3) then
         Scalx(Ibeg,j,k) = Din_X0(j)*Uin_X0(j,k)*Scal_X0(j,k)
       endif
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif


# if defined (PARALLEL)
     if(n_east.eq.MPI_PROC_NULL) then
# endif
     do j = Jbeg,Jend
     do k = Kbeg,Kend
       if(Bc_Xn==1.or.Bc_Xn==2) then
         Scalx(Iend1,j,k) = Zero
!  fyshi add this inflow at right bounary
       elseif(Bc_Xn==3) then
         Scalx(Iend1,j,k) = Din_Xn(j)*Uin_Xn(j,k)*Scal_Xn(j,k)
       endif
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif      

     ! front and back side
# if defined (PARALLEL)
     if(n_suth.eq.MPI_PROC_NULL) then
# endif
     do i = Ibeg,Iend
     do k = Kbeg,Kend
       if(Bc_Y0==1.or.Bc_Y0==2) then
         Scaly(i,Jbeg,k) = Zero
       endif
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif


# if defined (PARALLEL)
     if(n_nrth.eq.MPI_PROC_NULL) then
# endif
     do i = Ibeg,Iend
     do k = Kbeg,Kend
       if(Bc_Yn==1.or.Bc_Yn==2) then
         Scaly(i,Jend1,k) = Zero
       endif
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif

    do k = Kbeg,Kend
    do j = Jbeg,Jend
    do i = Ibeg,Iend
      if(Mask(i,j)==0) then
        Scalx(i,j,k) = Zero
        Scalx(i+1,j,k) = Zero
        Scaly(i,j,k) = Zero
        Scaly(i,j+1,k) = Zero
      endif
    enddo
    enddo
    enddo

    return
    end subroutine flux_scalar_bc


    subroutine les_3D(ISTEP)
!------------------------------------------------------
!   large eddy simulation (LES)
!   Last update: Gangfeng Ma, 09/22/2011
!------------------------------------------------------
    use global
    implicit none
    integer, intent(in) :: ISTEP
    integer :: i,j,k,Iter
    real(SP) :: S11,S22,S33,S12,S13,S23,SijSij,Filter
    real(SP) :: Umag,Zdis,X0,Xa,Xn,FricU
 
    CmuVt = Cmu
    do i = Ibeg,Iend
    do j = Jbeg,Jend
    do k = Kbeg,Kend
      if(Mask9(i,j)==1) then
        S11 = (U(i+1,j,k)-U(i-1,j,k))/(2.0*dx)+(U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k) 
        S22 = (V(i,j+1,k)-V(i,j-1,k))/(2.0*dy)+(V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k)
        S33 = 1./D(i,j)*(W(i,j,k+1)-W(i,j,k-1))/(sigc(k+1)-sigc(k-1))
        S12 = 0.5*((U(i,j+1,k)-U(i,j-1,k))/(2.0*dy)+(U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k)+  &
              (V(i+1,j,k)-V(i-1,j,k))/(2.0*dx)+(V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k)) 
        S13 = 0.5*((U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))/D(i,j)+  &                            
              (W(i+1,j,k)-W(i-1,j,k))/(2.0*dx)+(W(i,j,k+1)-W(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k))
        S23 = 0.5*((V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))/D(i,j)+  &                
              (W(i,j+1,k)-W(i,j-1,k))/(2.0*dy)+(W(i,j,k+1)-W(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k))
        SijSij = S11**2+S22**2+S33**2+2.0*(S12**2+S13**2+S23**2)
        Filter = (dx*dy*dsig(k)*D(i,j))**(1./3.)
        CmuVt(i,j,k) = Cmu(i,j,k)+(Cvs*Filter)**2*sqrt(2.0*SijSij)
      endif
    enddo
    enddo
    enddo

    ! ghost cells
    ! at the bottom
    do i = Ibeg,Iend
    do j = Jbeg,Jend
      ! impose wall function
      Umag = sqrt(U(i,j,Kbeg)**2+V(i,j,Kbeg)**2)
      if(Umag<1.e-6) then
        CmuVt(i,j,Kbeg) = Cmu(i,j,Kbeg)
      else
        Zdis = 0.5*dsig(Kbeg)*D(i,j)
        if(Zob>=0.001) then
          ! rough wall
          FricU = Umag/(1/0.41*log(30.0*Zdis/Zob))
        else
          ! smooth wall
          X0 = 0.05
          Iter = 0
       
          Xa = dlog(9.0*Umag*Zdis/Visc)
 10       Xn = X0+(0.41-X0*(Xa+dlog(X0)))/(1.0+0.41/X0)
          if(Iter>=20) then
            write(*,*) 'Iteration exceeds 20 steps',i,j,Umag
          endif
          if(dabs((Xn-X0)/X0)>1.e-8.and.Xn>0.0) then
            X0 = Xn
            Iter = Iter+1
            goto 10
          else
            FricU = Xn*Umag
          endif
        endif

        CmuVt(i,j,Kbeg) = Cmu(i,j,Kbeg)+0.41*Zdis*FricU*  &
           (1.0-exp(-Zdis*FricU/Visc/19.0))**2
      endif
 100  continue

      do k = 1,Nghost
        CmuVt(i,j,Kbeg-k) = CmuVt(i,j,Kbeg+k-1)
      enddo
    enddo
    enddo

    ! at the free surface
    do i = Ibeg,Iend
    do j = Jbeg,Jend
      do k = 1,Nghost
        CmuVt(i,j,Kend+k) = CmuVt(i,j,Kend-k+1)
      enddo
    enddo
    enddo

# if defined (PARALLEL)
    call phi_3D_exch(CmuVt)
# endif

# if defined (PARALLEL)
     if(n_west.eq.MPI_PROC_NULL) then
# endif
     do j = Jbeg,Jend
     do k = Kbeg,Kend
       if(WaveMaker(1:3)=='LEF') then
         ! no turbulence at wave generation region
         CmuVt(Ibeg,j,k) = Cmu(Ibeg,j,k)
         do i = 1,Nghost
           CmuVt(Ibeg-i,j,k) = Cmu(Ibeg-i,j,k)
         enddo
       else       
         do i = 1,Nghost
           CmuVt(Ibeg-i,j,k) = CmuVt(Ibeg+i-1,j,k)
         enddo
       endif
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_east.eq.MPI_PROC_NULL) then
# endif
     do j = Jbeg,Jend
     do k = Kbeg,Kend
       do i = 1,Nghost
         CmuVt(Iend+i,j,k) = CmuVt(Iend-i+1,j,k)
       enddo
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif
    
# if defined (PARALLEL)
     if(n_suth.eq.MPI_PROC_NULL) then
# endif
     do i = Ibeg,Iend
     do k = Kbeg,Kend
       do j = 1,Nghost
         CmuVt(i,Jbeg-j,k) = CmuVt(i,Jbeg+j-1,k)
       enddo
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_nrth.eq.MPI_PROC_NULL) then
# endif
     do i = Ibeg,Iend
     do k = Kbeg,Kend
       do j = 1,Nghost
         CmuVt(i,Jend+j,k) = CmuVt(i,Jend-j+1,k)
       enddo
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif       

     ! no turbulence in the internal wavemaker region   
    if(WaveMaker(1:3)=='INT') then
      do k = 1,Kloc
      do j = 1,Nloc
      do i = 1,Mloc
        if(xc(i)>=Xsource_West.and.xc(i)<=Xsource_East.and. &
            yc(j)>=Ysource_Suth.and.yc(j)<=Ysource_Nrth) then
          CmuVt(i,j,k) = Cmu(i,j,k)
        endif
      enddo
      enddo
      enddo
    endif

    end subroutine les_3D


    subroutine kepsilon_3D(ISTEP)
!-------------------------------------------------------
!   k-epsilon turbulence model
!   Last update: Gangfeng Ma, 09/07/2011
!-------------------------------------------------------
    use global
    implicit none
    integer,  intent(in) :: ISTEP
    integer,  parameter :: ke_model = 2
    real(SP), parameter :: Dmin = 0.05
    real(SP), dimension(:,:,:), allocatable :: R5,DelzR,Tke_Old,Eps_Old,Uamp
    real(SP), dimension(:,:), allocatable :: VelGrad,ReynoldStress
    real(SP) :: c1e,c2e,c3e,cmiu,Umag,Zdis,X0,Xa,Xn,FricU,Tkeb,Epsb
    real(SP) :: smax,dmax,c_d,c_1,c_2,c_3,delta_nm
    real(SP) :: S11,S22,S33,S12,S13,S23,TKE_CAP
    integer :: i,j,k,n,m,l,g,IVAR,Iter
    
    allocate(R5(Mloc,Nloc,Kloc))
    allocate(DelzR(Mloc,Nloc,Kloc))
    allocate(Tke_Old(Mloc,Nloc,Kloc))
    allocate(Eps_Old(Mloc,Nloc,Kloc))
    allocate(VelGrad(3,3))
    allocate(ReynoldStress(3,3))
    allocate(Uamp(Mloc,Nloc,Kloc))

    ! some parameters
    c1e = 1.44
    c2e = 1.92
    c3e = -1.4
    cmiu = 0.09

    ! save old values
    Tke_Old = Tke
    Eps_Old = Eps

    Prod_s = Zero
    do i = Ibeg,Iend
    do j = Jbeg,Jend
    do k = Kbeg,Kend
      if(ke_model==1) then
        ! linear model
        if(D(i,j)>=Dmin) then
          S11 = (U(i+1,j,k)-U(i-1,j,k))/(2.0*dx)+(U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k)
          S22 = (V(i,j+1,k)-V(i,j-1,k))/(2.0*dy)+(V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k)
          S33 = 1./D(i,j)*(W(i,j,k+1)-W(i,j,k-1))/(sigc(k+1)-sigc(k-1)) 
          S12 = 0.5*((U(i,j+1,k)-U(i,j-1,k))/(2.0*dy)+(U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k)+  &
              (V(i+1,j,k)-V(i-1,j,k))/(2.0*dx)+(V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k))                   
          S13 = 0.5*((U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))/D(i,j)+  &                                            
              (W(i+1,j,k)-W(i-1,j,k))/(2.0*dx)+(W(i,j,k+1)-W(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k))                   
          S23 = 0.5*((V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))/D(i,j)+  &                                            
              (W(i,j+1,k)-W(i,j-1,k))/(2.0*dy)+(W(i,j,k+1)-W(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k))
          Prod_s(i,j,k) = 2.0*CmuVt(i,j,k)*(S11**2+S22**2+S33**2+2.0*(S12**2+S13**2+S23**2))
        endif
      elseif(ke_model==2) then
        ! nonlinear model (Lin and Liu, 1998)
        if(D(i,j)>=Dmin) then
          ! estimate gradient first
          VelGrad = Zero

          VelGrad(1,1) = (U(i+1,j,k)-U(i-1,j,k))/(2.0*dx)+(U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k)
          VelGrad(1,2) = (U(i,j+1,k)-U(i,j-1,k))/(2.0*dy)+(U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k)
          VelGrad(1,3) = 1./D(i,j)*(U(i,j,k+1)-U(i,j,k-1))/(sigc(k+1)-sigc(k-1))
          VelGrad(2,1) = (V(i+1,j,k)-V(i-1,j,k))/(2.0*dx)+(V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k)
          VelGrad(2,2) = (V(i,j+1,k)-V(i,j-1,k))/(2.0*dy)+(V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k) 
          VelGrad(2,3) = 1./D(i,j)*(V(i,j,k+1)-V(i,j,k-1))/(sigc(k+1)-sigc(k-1))
          VelGrad(3,1) = (W(i+1,j,k)-W(i-1,j,k))/(2.0*dx)+(W(i,j,k+1)-W(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelxSc(i,j,k)
          VelGrad(3,2) = (W(i,j+1,k)-W(i,j-1,k))/(2.0*dy)+(W(i,j,k+1)-W(i,j,k-1))/(sigc(k+1)-sigc(k-1))*DelySc(i,j,k)
          VelGrad(3,3) = 1./D(i,j)*(W(i,j,k+1)-W(i,j,k-1))/(sigc(k+1)-sigc(k-1))

          ! estimate Reynolds stress
          ReynoldStress = Zero
          
          smax = zero
          do n = 1,3
            if(abs(VelGrad(n,n))>smax) smax = abs(VelGrad(n,n))
          enddo
          smax = smax*Tke_Old(i,j,k)/Eps_Old(i,j,k)
          c_d = (1./3.)*(1./(3.7+smax))

          dmax = zero
          do n = 1,3
          do m = 1,3
            if(abs(VelGrad(n,m))>dmax) dmax = abs(VelGrad(n,m))
          enddo
          enddo
          dmax = dmax*Tke_Old(i,j,k)/Eps_Old(i,j,k)
          c_1 = 2./3./(123.5+2.0*dmax**2)
          c_2 = -2./3./(39.2+2.0*dmax**2)
          c_3 = 2./3./(246.9+2.0*dmax**2)

          do n = 1,3
          do m = 1,3
            if(n==m) then
              delta_nm = 1.
            else
              delta_nm = 0.
            endif

            ReynoldStress(n,m) = c_d*Tke_Old(i,j,k)**2/Eps_Old(i,j,k)*(VelGrad(n,m)+VelGrad(m,n))-  &
                   (2./3.)*Tke_Old(i,j,k)*delta_nm

            do l = 1,3
              ReynoldStress(n,m) = ReynoldStress(n,m)+  &
                   c_1*Tke_Old(i,j,k)**3/Eps_Old(i,j,k)**2*(VelGrad(n,l)*VelGrad(l,m)+  &
                   VelGrad(m,l)*VelGrad(l,n))
              do g = 1,3
                ReynoldStress(n,m) = ReynoldStress(n,m)-  &
                   (2./3.)*c_1*Tke_Old(i,j,k)**3/Eps_Old(i,j,k)**2*VelGrad(l,g)*VelGrad(g,l)*delta_nm
              enddo
            enddo

            do g = 1,3
              ReynoldStress(n,m) = ReynoldStress(n,m)+  &
                   c_2*Tke_Old(i,j,k)**3/Eps_Old(i,j,k)**2*VelGrad(n,g)*VelGrad(m,g)
              do l = 1,3
                ReynoldStress(n,m) = ReynoldStress(n,m)-  &
                   (1./3.)*c_2*Tke_Old(i,j,k)**3/Eps_Old(i,j,k)**2*VelGrad(l,g)*VelGrad(l,g)*delta_nm
              enddo
            enddo

            do g = 1,3
              ReynoldStress(n,m) = ReynoldStress(n,m)+  &
                   c_3*Tke_Old(i,j,k)**3/Eps_Old(i,j,k)**2*VelGrad(g,n)*VelGrad(g,m)
              do l = 1,3
                ReynoldStress(n,m) = ReynoldStress(n,m)-  &
                   (1./3.)*c_3*Tke_Old(i,j,k)**3/Eps_Old(i,j,k)**2*VelGrad(l,g)*VelGrad(l,g)*delta_nm 
              enddo
            enddo
          enddo
          enddo

          ! estimate shear production
          do n = 1,3
          do m = 1,3
            Prod_s(i,j,k) = Prod_s(i,j,k)+ReynoldStress(n,m)*VelGrad(n,m)*MASK9(I,J)
          enddo
          enddo

        endif
      endif
    enddo
    enddo
    enddo

tmp4preview3D=Prod_s

    ! buoyancy production
    Prod_b = Zero

    ! slove tke equation
    IVAR = 1
    call adv_scalar_tvd(Tke,R5,IVAR)

    do i = Ibeg,Iend
    do j = Jbeg,Jend
    do k = Kbeg,Kend
      if(D(i,j)>Dmin) then
        R5(i,j,k) = R5(i,j,k)+D(i,j)*(Prod_s(i,j,k)+Prod_b(i,j,k)-Eps_Old(i,j,k))
        DTke(i,j,k) = ALPHA(ISTEP)*DTke0(i,j,k)+BETA(ISTEP)*(DTke(i,j,k)+dt*R5(i,j,k))
        DTke(i,j,k) = max(DTke(i,j,k),D(i,j)*Tke_min)
      endif
    enddo
    enddo
    enddo

    ! solve epsilon equation
    IVAR = 2
    call adv_scalar_tvd(Eps,R5,IVAR)

    do i = Ibeg,Iend
    do j = Jbeg,Jend
    do k = Kbeg,Kend
      if(D(i,j)>=Dmin) then
        R5(i,j,k) = R5(i,j,k)+D(i,j)*(c1e*(Prod_s(i,j,k)+c3e*Prod_b(i,j,k))*Eps_Old(i,j,k)/Tke_Old(i,j,k)  &
                     -c2e*Eps_Old(i,j,k)**2/Tke_Old(i,j,k))                              
        DEps(i,j,k) = ALPHA(ISTEP)*DEps0(i,j,k)+BETA(ISTEP)*(DEps(i,j,k)+dt*R5(i,j,k))
        DEps(i,j,k) = max(DEps(i,j,k),D(i,j)*Eps_min)
      endif
    enddo
    enddo
    enddo

    ! ghost cells
    ! at the bottom
    do i = Ibeg,Iend
    do j = Jbeg,Jend
      ! impose wall function
      Umag = sqrt(U(i,j,Kbeg)**2+V(i,j,Kbeg)**2)
      if(Umag<1.e-6.or.D(i,j)<Dmin) cycle

      Zdis = 0.5*dsig(Kbeg)*D(i,j)
      X0 = 0.05
      Iter = 0
       
      Xa = dlog(9.0*Umag*Zdis/Visc)
 10   Xn = X0+(0.41-X0*(Xa+dlog(X0)))/(1.0+0.41/X0)
      if(Iter>=20) then
        write(*,*) 'Iteration exceeds 20 steps',i,j,Umag
      endif
      if(dabs((Xn-X0)/X0)>1.e-8.and.Xn>0.0) then
        X0 = Xn
        Iter = Iter+1
        goto 10
      else
        FricU = Xn*Umag
      endif

      Tkeb = FricU**2/sqrt(cmiu)
      Epsb = FricU**3/(0.41*Zdis)

      DTke(i,j,Kbeg) = D(i,j)*Tkeb
      DEps(i,j,Kbeg) = D(i,j)*Epsb

      do k = 1,Nghost
        DTke(i,j,Kbeg-k) = DTke(i,j,Kbeg+k-1)
        DEps(i,j,Kbeg-k) = DEps(i,j,Kbeg+k-1)
      enddo
    enddo
    enddo

    ! at the free surface
    do i = Ibeg,Iend
    do j = Jbeg,Jend
      do k = 1,Nghost
        DTke(i,j,Kend+k) = DTke(i,j,Kend-k+1)
        DEps(i,j,Kend+k) = DEps(i,j,Kend-k+1)
      enddo
    enddo
    enddo

# if defined (PARALLEL)
    call phi_3D_exch(DTke)
    call phi_3D_exch(DEps)
# endif


! fyshi apply tke cap because tke always exceed normal value, not enough epsilon?
     Uamp=sqrt(U**2+V**2+W**2)
     Umag=MAXVAL(Uamp)

     DO J=Jbeg,Jend
     DO I=Ibeg,Iend
     DO K=Kbeg,Kend
        TKE_CAP=0.5*Umag*Umag*D(I,J)
        IF(DTKE(I,J,K).GE.TKE_CAP)DTKE(I,J,K)=TKE_CAP
     ENDDO
     ENDDO
     ENDDO


# if defined (PARALLEL)
     if(n_west.eq.MPI_PROC_NULL) then
# endif
     do j = Jbeg,Jend
     do k = Kbeg,Kend
       if(WaveMaker(1:3)=='LEF') then
         DTke(Ibeg,j,k) = D(Ibeg,j)*Tke_min
         DEps(Ibeg,j,k) = D(Ibeg,j)*Eps_min
         do i = 1,Nghost
           DTke(Ibeg-i,j,k) = D(Ibeg-i,j)*Tke_min
           DEps(Ibeg-i,j,k) = D(Ibeg-i,j)*Eps_min
         enddo
       else
         do i = 1,Nghost
           DTke(Ibeg-i,j,k) = DTke(Ibeg+i-1,j,k)
           DEps(Ibeg-i,j,k) = DEps(Ibeg+i-1,j,k)
         enddo
       endif
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_east.eq.MPI_PROC_NULL) then
# endif
     do j = Jbeg,Jend
     do k = Kbeg,Kend
       do i = 1,Nghost
         DTke(Iend+i,j,k) = DTke(Iend-i+1,j,k)
         DEps(Iend+i,j,k) = DEps(Iend-i+1,j,k)
       enddo
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif
    
# if defined (PARALLEL)
     if(n_suth.eq.MPI_PROC_NULL) then
# endif
     do i = Ibeg,Iend
     do k = Kbeg,Kend
       do j = 1,Nghost
         DTke(i,Jbeg-j,k) = DTke(i,Jbeg+j-1,k)
         DEps(i,Jbeg-j,k) = DEps(i,Jbeg+j-1,k)
       enddo
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_nrth.eq.MPI_PROC_NULL) then
# endif
     do i = Ibeg,Iend
     do k = Kbeg,Kend
       do j = 1,Nghost
         DTke(i,Jend+j,k) = DTke(i,Jend-j+1,k)
         DEps(i,Jend+j,k) = DEps(i,Jend-j+1,k)
       enddo
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif

    do i = 1,Mloc
    do j = 1,Nloc
    do k = 1,Kloc
      if(D(i,j)>=Dmin) then
        Tke(i,j,k) = DTke(i,j,k)/D(i,j)
        Eps(i,j,k) = DEps(i,j,k)/D(i,j)
        CmuVt(i,j,k) = Cmu(i,j,k)+Cmiu*Tke(i,j,k)**2/Eps(i,j,k)
        CmuVt(i,j,k) = min(CmuVt(i,j,k),1.e-2)
      else
        Tke(i,j,k) = Tke_min
        Eps(i,j,k) = Eps_min
        DTke(i,j,k) = D(i,j)*Tke_min
        DEps(i,j,k) = D(i,j)*Eps_min
        CmuVt(i,j,k) = Cmu(i,j,k)
      endif
    enddo
    enddo
    enddo

    ! no turbulence in the internal wavemaker region
    if(WaveMaker(1:3)=='INT') then
      do k = 1,Kloc
      do j = 1,Nloc
      do i = 1,Mloc
        if(xc(i)>=Xsource_West.and.xc(i)<=Xsource_East.and. &
            yc(j)>=Ysource_Suth.and.yc(j)<=Ysource_Nrth) then
          Tke(i,j,k) = Tke_min
          Eps(i,j,k) = Eps_min
          DTke(i,j,k) = D(i,j)*Tke_min
          DEps(i,j,k) = D(i,j)*Eps_min
          CmuVt(i,j,k) = Cmu(i,j,k)
        endif  
      enddo
      enddo
      enddo
    endif

    deallocate(R5)
    deallocate(DelzR)
    deallocate(Tke_Old)
    deallocate(Eps_Old)
    deallocate(VelGrad)
    deallocate(ReynoldStress)

    end subroutine kepsilon_3D


    subroutine kepsilon(ISTEP)
!-------------------------------------------------------
!   k-epsilon turbulence model
!   Last update: Gangfeng Ma, 09/07/2011
!-------------------------------------------------------
    use global
    implicit none
    integer, intent(in) :: ISTEP
    real(SP), parameter :: Dmin = 0.05
    real(SP), dimension(:,:,:), allocatable :: R5,DelzR
    real(SP) :: c1e,c2e,c3e,cmiu,Umag,Zdis,X0,Xa,Xn,FricU,Tkeb,Epsb
    integer :: i,j,k,IVAR,Iter

    allocate(R5(Mloc,Nloc,Kloc))
    allocate(DelzR(Mloc,Nloc,Kloc))

    ! some parameters
    c1e = 1.44
    c2e = 1.92
    c3e = -1.4
    cmiu = 0.09

    ! shear production
    Prod_s = Zero
    do i = Ibeg,Iend
    do j = Jbeg,Jend
    do k = Kbeg,Kend
      if(Mask9(i,j)==1) then
        DelzU(i,j,k) = (Uf(i,j,k+1)-Uf(i,j,k))/dsig(k)
        DelzV(i,j,k) = (Vf(i,j,k+1)-Vf(i,j,k))/dsig(k)
        Prod_s(i,j,k) = CmuVt(i,j,k)*(DelzU(i,j,k)**2+DelzV(i,j,k)**2)/D(i,j)**2   
      endif
    enddo
    enddo
    enddo

    ! buoyancy production
    Prod_b = Zero
!    do i = Ibeg,Iend
!    do j = Jbeg,Jend
!    do k = Kbeg,Kend
!      Prod_b(i,j,k) = -Grav/Rho0*CmuVt(i,j,k)*DelzR(i,j,k)/D(i,j)
!    enddo
!    enddo
!    enddo

    ! slove tke equation
    IVAR = 1
    call adv_scalar(Tke,R5,IVAR)

    do i = Ibeg,Iend
    do j = Jbeg,Jend
    do k = Kbeg,Kend
      if(D(i,j)>Dmin) then
        R5(i,j,k) = R5(i,j,k)+D(i,j)*(Prod_s(i,j,k)+Prod_b(i,j,k)-Eps(i,j,k))
        DTke(i,j,k) = ALPHA(ISTEP)*DTke0(i,j,k)+BETA(ISTEP)*(DTke(i,j,k)+dt*R5(i,j,k))
        DTke(i,j,k) = max(DTke(i,j,k),D(i,j)*Tke_min)
      endif
    enddo
    enddo
    enddo

    ! solve epsilon equation
    IVAR = 2
    call adv_scalar(Eps,R5,IVAR)

    do i = Ibeg,Iend
    do j = Jbeg,Jend
    do k = Kbeg,Kend
      if(D(i,j)>Dmin) then
        R5(i,j,k) = R5(i,j,k)+D(i,j)*(c1e*(Prod_s(i,j,k)+c3e*Prod_b(i,j,k))*Eps(i,j,k)/(Tke(i,j,k)+1.e-16)  &
                     -c2e*Eps(i,j,k)**2/(Tke(i,j,k)+1.e-16))
        DEps(i,j,k) = ALPHA(ISTEP)*DEps0(i,j,k)+BETA(ISTEP)*(DEps(i,j,k)+dt*R5(i,j,k))
        DEps(i,j,k) = max(DEps(i,j,k),D(i,j)*Eps_min)
      endif
    enddo
    enddo
    enddo

    ! ghost cells
    ! at the bottom
     do i = Ibeg,Iend
     do j = Jbeg,Jend
       ! impose wall function (no-slip bc)
       Umag = sqrt(U(i,j,Kbeg)**2+V(i,j,Kbeg)**2)
       if(Umag<1.e-6) cycle

       Zdis = 0.5*dsig(Kbeg)*D(i,j)
       X0 = 0.05
       Iter = 0
       
       Xa = dlog(9.0*Umag*Zdis/Visc)
 10    Xn = X0+(0.41-X0*(Xa+dlog(X0)))/(1.0+0.41/X0)
       if(Iter>=20) then
         write(*,*) 'Iteration exceeds 20 steps',i,j,Umag
       endif
       if(dabs((Xn-X0)/X0)>1.e-8.and.Xn>0.0) then
         X0 = Xn
         Iter = Iter+1
         goto 10
       else
         FricU = Xn*Umag
       endif

       Tkeb = FricU**2/sqrt(cmiu)
       Epsb = FricU**2/(0.41*Zdis)

       DTke(i,j,Kbeg) = D(i,j)*Tkeb
       DEps(i,j,Kbeg) = D(i,j)*Epsb

       do k = 1,Nghost
         DTke(i,j,Kbeg-k) = DTke(i,j,Kbeg+k-1)
         DEps(i,j,Kbeg-k) = DEps(i,j,Kbeg+k-1)
       enddo
     enddo
     enddo

     ! at the free surface
     do i = Ibeg,Iend
     do j = Jbeg,Jend
       do k = 1,Nghost
         DTke(i,j,Kend+k) = DTke(i,j,Kend-k+1)
         DEps(i,j,Kend+k) = DEps(i,j,Kend-k+1)
       enddo
     enddo
     enddo

# if defined (PARALLEL)
    call phi_3D_exch(DTke)
    call phi_3D_exch(DEps)
# endif

# if defined (PARALLEL)
     if(n_west.eq.MPI_PROC_NULL) then
# endif
     do j = Jbeg,Jend
     do k = Kbeg,Kend
       do i = 1,Nghost
         DTke(Ibeg-i,j,k) = DTke(Ibeg+i-1,j,k)
         DEps(Ibeg-i,j,k) = DEps(Ibeg+i-1,j,k)
       enddo
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_east.eq.MPI_PROC_NULL) then
# endif
     do j = Jbeg,Jend
     do k = Kbeg,Kend
       do i = 1,Nghost
         DTke(Iend+i,j,k) = DTke(Iend-i+1,j,k)
         DEps(Iend+i,j,k) = DEps(Iend-i+1,j,k)
       enddo
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif
    
# if defined (PARALLEL)
     if(n_suth.eq.MPI_PROC_NULL) then
# endif
     do i = Ibeg,Iend
     do k = Kbeg,Kend
       do j = 1,Nghost
         DTke(i,Jbeg-j,k) = DTke(i,Jbeg+j-1,k)
         DEps(i,Jbeg-j,k) = DEps(i,Jbeg+j-1,k)
       enddo
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif

# if defined (PARALLEL)
     if(n_nrth.eq.MPI_PROC_NULL) then
# endif
     do i = Ibeg,Iend
     do k = Kbeg,Kend
       do j = 1,Nghost
         DTke(i,Jend+j,k) = DTke(i,Jend-j+1,k)
         DEps(i,Jend+j,k) = DEps(i,Jend-j+1,k)
       enddo
     enddo
     enddo
# if defined (PARALLEL)
     endif
# endif

    do i = 1,Mloc
    do j = 1,Nloc
    do k = 1,Kloc
      if(D(i,j)>Dmin) then
        Tke(i,j,k) = DTke(i,j,k)/D(i,j)
        Eps(i,j,k) = DEps(i,j,k)/D(i,j)
        CmuVt(i,j,k) = Cmu(i,j,k)+Cmiu*Tke(i,j,k)**2/Eps(i,j,k)
      else
        Tke(i,j,k) = Tke_min
        Eps(i,j,k) = Eps_min
        DTke(i,j,k) = D(i,j)*Tke_min
        DEps(i,j,k) = D(i,j)*Eps_min
        CmuVt(i,j,k) = Cmu(i,j,k)
      endif
    enddo
    enddo
    enddo

    ! no turbulence in the internal wavemaker region
    if(WaveMaker(1:3)=='INT') then
      do k = 1,Kloc
      do j = 1,Nloc
      do i = 1,Mloc
        if(xc(i)>=Xsource_West.and.xc(i)<=Xsource_East.and. &
            yc(j)>=Ysource_Suth.and.yc(j)<=Ysource_Nrth) then
          Tke(i,j,k) = Tke_min
          Eps(i,j,k) = Eps_min
          DTke(i,j,k) = D(i,j)*Tke_min
          DEps(i,j,k) = D(i,j)*Eps_min
          CmuVt(i,j,k) = Cmu(i,j,k)
        endif
      enddo
      enddo
      enddo
    endif

    deallocate(R5)
    deallocate(DelzR)

    end subroutine kepsilon


    subroutine wall_time_secs(tcurrent)
!--------------------------------------------------------
!   Calculate current wall time
!   Last update: Gangfeng Ma, 09/12/2011
!--------------------------------------------------------
    use global, only: SP
    implicit none
    integer, dimension(8) :: walltime
    real(SP), intent(out) :: tcurrent
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
    end subroutine wall_time_secs


# if defined (BUBBLE)
    subroutine bslip_velocity
!----------------------------------------------------------
!   Specify bubble radius and calculate rise velocity
!   Last update: Gangfeng Ma, 09/01/2011
!---------------------------------------------------------
    use global, only: SP,pi,Zero,Rho0,Mg,Rbg,Wbg,Entrain,Con_b,Surface_Tension
    implicit none
    integer :: g
    real(SP) :: rlogR0,rlogRn,rlogR1,rlogR2,rlogN0,alpha_b,beta_b,sum_e
    real(SP), dimension(Mg) :: DRbg,rlogN,specN,binN             

    ! specify bubble radius
    rlogR0 = -1.0
    rlogRn = 1.0
    do g = 1,Mg
      rlogR1 = rlogR0+(rlogRn-rlogR0)*(g-1)/float(Mg)
      rlogR2 = rlogR0+(rlogRn-rlogR0)*g/float(Mg)
      Rbg(g) = (0.5*(10**rlogR1+10**rlogR2))*0.001
      DRbg(g) = (10**rlogR2-10**rlogR1)*0.001
    enddo

    ! slip velocity
    do g = 1,Mg
      if(Rbg(g)<=7.0e-4) then
        Wbg(g) = 4474.*Rbg(g)**1.357
      elseif(Rbg(g)>5.1e-3) then
        Wbg(g) = 4.202*Rbg(g)**0.547
      else
        Wbg(g) = 0.23
      endif
    enddo

    ! read bubble size distribution in Deane and Stokes' paper
    rlogN0 = 4.3
    alpha_b = -3.0/2.0
    beta_b = -10.0/3.0

    rlogN(1) = rlogN0
    do g = 2,Mg/2
      rlogN(g) = rlogN(g-1)+alpha_b*(log10(Rbg(g))-log10(Rbg(g-1)))
    enddo
    do g = Mg/2+1,Mg
      rlogN(g) = rlogN(g-1)+beta_b*(log10(Rbg(g))-log10(Rbg(g-1)))
    enddo

    do g = 1,Mg
      specN(g) = 10**rlogN(g)*1.0e+6  ! convert to 1/(m^4)
      binN(g)  = specN(g)*DRbg(g)
    enddo

    ! entrainment coefficient
    sum_e = zero
    do g = 1,Mg
      sum_e = sum_e+Rbg(g)**2*binN(g)
    enddo

    do g = 1,Mg
      Entrain(g) = Con_b/(4.0*pi)*(Surface_Tension/Rho0)**(-1)*binN(g)/sum_e
    enddo

    end subroutine bslip_velocity

    
    subroutine eval_bub(ISTEP)
!---------------------------------------------------------------------
!   Update bubble concentration
!   Last update: Gangfeng Ma, 09/01/2011
!---------------------------------------------------------------------
    use global
    implicit none
    integer, intent(in) :: ISTEP
    real(SP), dimension(:,:,:), allocatable :: R5,Phi,DPhi
    integer :: i,j,k,g,IVAR 

    allocate(R5(Mloc,Nloc,Kloc))
    allocate(Phi(Mloc,Nloc,Kloc))
    allocate(DPhi(Mloc,Nloc,Kloc))

    do g = 1,Mg
      ! temporary arrays
      do k = 1,Kloc
      do j = 1,Nloc
      do i = 1,Mloc
        R5(i,j,k) = Zero 
        Phi(i,j,k) = Nbg(i,j,k,g)
        DPhi(i,j,k) = DNbg(i,j,k,g)
      enddo
      enddo
      enddo

      ! advection and diffusion 
      IVAR = 5     
      call adv_scalar(Phi,R5,IVAR)

      ! account for bubble rise
      do k = Kbeg,Kend
      do j = Jbeg,Jend
      do i = Ibeg,Iend
        R5(i,j,k) = R5(i,j,k)-  &
             Wbg(g)*(Phi(i,j,k+1)-Phi(i,j,k-1))/(sigc(k+1)-sigc(k-1))
      enddo
      enddo
      enddo

      ! bubble entrainment at the surface
      k = Kend
      do j = Jbeg,Jend
      do i = Ibeg,Iend
        if(Eps(i,j,k)>=0.01) then
          R5(i,j,k) = R5(i,j,k)+(1.0-Vbg(i,j,k))*Entrain(g)*D(i,j)*Eps(i,j,k)
        endif
      enddo
      enddo
 
      do k = Kbeg,Kend
      do j = Jbeg,Jend
      do i = Ibeg,Iend
        DNbg(i,j,k,g) = ALPHA(ISTEP)*DNbg0(i,j,k,g)+  &
                           BETA(ISTEP)*(DNbg(i,j,k,g)+dt*R5(i,j,k))
        if(Mask(i,j)==0) DNbg(i,j,k,g) = Zero
      enddo
      enddo
      enddo
    enddo

    ! collect data into ghost cells
# if defined (PARALLEL)
    do g = 1,Mg
      do k = 1,Kloc
      do j = 1,Nloc
      do i = 1,Mloc
        DPhi(i,j,k) = DNbg(i,j,k,g)
      enddo
      enddo
      enddo

      call phi_3D_exch(DPhi)

      do k = 1,Kloc
      do j = 1,Nloc
      do i = 1,Mloc
        DNbg(i,j,k,g) = DPhi(i,j,k)
      enddo
      enddo
      enddo
    enddo
# endif          

# if defined (PARALLEL)
    if(n_west.eq.MPI_PROC_NULL) then
# endif
    do g = 1,Mg
    do j = Jbeg,Jend
    do k = Kbeg,Kend
      do i = 1,Nghost
        DNbg(Ibeg-i,j,k,g) = DNbg(Ibeg+i-1,j,k,g)
      enddo
    enddo
    enddo
    enddo
# if defined (PARALLEL)
    endif
# endif

# if defined (PARALLEL)
    if(n_east.eq.MPI_PROC_NULL) then
# endif
    do g = 1,Mg
    do j = Jbeg,Jend
    do k = Kbeg,Kend
      do i = 1,Nghost
        DNbg(Iend+i,j,k,g) = DNbg(Iend-i+1,j,k,g)
      enddo
    enddo
    enddo
    enddo
# if defined (PARALLEL)
    endif
# endif

# if defined (PARALLEL)
    if(n_suth.eq.MPI_PROC_NULL) then
# endif
    do g = 1,Mg
    do i = Ibeg,Iend
    do k = Kbeg,Kend
      do j = 1,Nghost
        DNbg(i,Jbeg-j,k,g) = DNbg(i,Jbeg+j-1,k,g)
      enddo
    enddo
    enddo
    enddo
# if defined (PARALLEL)
    endif
# endif 

# if defined (PARALLEL)
    if(n_nrth.eq.MPI_PROC_NULL) then
# endif
    do g = 1,Mg
    do i = Ibeg,Iend
    do k = Kbeg,Kend
      do j = 1,Nghost
        DNbg(i,Jend+j,k,g) = DNbg(i,Jend-j+1,k,g)
      enddo
    enddo
    enddo
    enddo
# if defined (PARALLEL)
    endif
# endif

    do g = 1,Mg
    do i = Ibeg,Iend
    do j = Jbeg,Jend
      do k = 1,Nghost
        DNbg(i,j,Kbeg-k,g) = DNbg(i,j,Kbeg+k-1,g)
      enddo
      do k = 1,Nghost
        DNbg(i,j,Kend+k,g) = DNbg(i,j,Kend-k+1,g)
      enddo
    enddo
    enddo
    enddo

    ! update bubble number density
    Nbg = Zero
    do g = 1,Mg
    do k = 1,Kloc
    do j = 1,Nloc
    do i = 1,Mloc
      if(Mask(i,j)==1) then
        Nbg(i,j,k,g) = DNbg(i,j,k,g)/D(i,j)
      endif
    enddo
    enddo
    enddo
    enddo

    ! calculate void fraction
    Vbg = Zero
    do k = 1,Kloc
    do j = 1,Nloc
    do i = 1,Mloc
      do g = 1,Mg
        Vbg(i,j,k) = Vbg(i,j,k)+4./3.*pi*Rbg(g)**3*Nbg(i,j,k,g)
      enddo
    enddo
    enddo
    enddo

!    ! update cell density
!    do k = 1,Kloc
!    do j = 1,Nloc
!    do i = 1,Mloc
!      Rho(i,j,k) = (1.0-Vbg(i,j,k))*Rho0+Vbg(i,j,k)*RhoA
!    enddo
!    enddo
!    enddo

    deallocate(R5)
    deallocate(Phi)
    deallocate(DPhi)

    end subroutine eval_bub
# endif
