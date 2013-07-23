module RefineRemeshModule

use NumberKindsModule
use SphereGeomModule
use LoggerModule
use TracerAndVorticityDistributionModule
use ParticlesEdgesPanelsModule

! Must link to oldstripack.o and ssrfpack.o

implicit none
private
public InitRefine
public AdaptiveRemesh, DirectAdaptiveRemesh
public SINGLE_GAUSSIAN_VORTEX, RH4_WAVE, MULTIPLE_VORTICES, RH2_WAVE, &
	   STRATOSPHERE_MODEL, JET, ZONAL_MEAN, GAUSSIAN_HILLS, COSINE_BELLS, BLOCKM,&
	   SLOTTED_CYLINDERS, TWO_VORTS, TRIPOLE, SOLID_BODY
public SetMaxRefinementLimit, MAX_REFINEMENT

integer(kint), parameter :: logLevel=TRACE_LOGGING_LEVEL
type(Logger) :: log
character(len=28) :: logKey='RefineRemesh : '
character(len=256) :: logString
logical(klog), save :: logInit = .False.

integer(kint), parameter :: SINGLE_GAUSSIAN_VORTEX = 41, &
							RH4_WAVE = 42,&
							MULTIPLE_VORTICES = 43, &
							RH2_WAVE = 44, &
							STRATOSPHERE_MODEL = 45, &
							JET = 46, &
							ZONAL_MEAN = 47, &
							GAUSSIAN_HILLS = 48, &
							COSINE_BELLS = 49, &
							BLOCKM = 50,&
							SLOTTED_CYLINDERS = 51,&
							TWO_VORTS = 40, &
							TRIPOLE = 39, &
							SOLID_BODY = 52


integer(kint), save :: MAX_REFINEMENT = 2

contains

subroutine SetMaxRefinementLimit(newLimit)
	integer(kint), intent(in) :: newLimit
	MAX_REFINEMENT = newLimit
end subroutine


subroutine InitRefine(aParticles, anEdges, aPanels, circTol, varTol, problemID,&
	theta0,beta,perturbAmp,perturbWaveNum,&
	cent1,cent2,cent3, beta1,beta2,beta3, strength1, strength2, strength3, &
	procRank )
	type(Particles), intent(inout) :: aParticles
	type(Edges), intent(inout) :: anEdges
	type(Panels), intent(inout) :: aPanels
	real(kreal), intent(in) :: circTol, varTol
	integer(kint), intent(in) :: problemID
	integer(kint), intent(in), optional :: perturbWaveNum				! Options for Jet Test Case
	real(kreal), intent(in), optional :: beta, theta0, perturbAmp	! Options for Jet Test Case
	integer(kint), intent(in), optional :: procRank
	real(kreal), intent(in), optional :: cent1(3), cent2(3), cent3(3)
	real(kreal), intent(in), optional :: beta1, beta2, beta3, strength1, strength2, strength3
	integer(kint) :: panelKind, vertIndices(4), amrCount1, amrCount2, initNest
	integer(kint) :: amrLoopCounter, startIndex, maxNest, nOldParticles, nOldPanels, nOldPanels2
	real(kreal) :: maxx0(3), minx0(3), lagVar, dlambda0, dlambda1
	integer(kint) :: i, j, k, refineCount, spaceLeft
	logical(klog), allocatable :: refineFlag(:)
	logical(klog) :: maxRefinement, printMsg
	real(kreal), parameter :: zero = 0.0_kreal

	printMsg = .TRUE.
	if ( present(procRank) ) then
		write(logkey,'(A,I0.2,A)') 'RefineRemesh_',procRank,": "
		if ( procRank /= 0 ) printMsg = .FALSE.
	endif
	if ( .NOT. logInit) call InitLogger(log)

	call LogMessage(log,DEBUG_LOGGING_LEVEL,logKey,"Entering InitRefine.")

	maxRefinement = .FALSE.

	panelkind = GetPanelKind(aPanels)
	amrCount1 = 0
	amrCount2 = 0
	allocate(refineFlag(aPanels%N_Max))
	refineFlag = .False.
	initNest = maxVal(aPanels%nest)

	do j=1,aPanels%N
		if ( .NOT. aPanels%hasChildren(j) ) then
!!!!!
! AMR CRITERIA
!!!!!
			if ( (problemID >= GAUSSIAN_HILLS) .AND. ( problemID <= SLOTTED_CYLINDERS) ) then
! ADVECTION
				if ( abs(aPanels%tracer(j,1))*aPanels%area(j) > circTol ) then
					refineFlag(j) = .TRUE.
					amrCount1 = amrCount1 + 1
				endif
			else
! BVE		! AMR Criterion 1
				if ( abs(aPanels%relVort(j))*aPanels%area(j) > circTol ) then
					refineFlag(j) = .TRUE.
					amrCount1 = amrCount1 + 1
				else
				! AMR Criterion 2
					maxx0 = aPanels%x0(:,j)
					minx0 = aPanels%x0(:,j)
					vertIndices(1:panelKind) = aPanels%vertices(:,j)
					do k=1,panelKind
						if ( aParticles%x0(1,vertIndices(k)) < minX0(1) ) then
							minx0(1) = aParticles%x0(1,vertIndices(k))
						endif
						if ( aParticles%x0(1,vertIndices(k)) > maxx0(1) ) then
							maxx0(1) = aParticles%x0(1,vertIndices(k))
						endif
						if ( aParticles%x0(2,vertIndices(k)) < minx0(2) ) then
							minx0(2) = aParticles%x0(2,vertIndices(k))
						endif
						if ( aParticles%x0(2,vertIndices(k)) > maxx0(2) ) then
							maxx0(2) = aParticles%x0(2,vertIndices(k))
						endif
						if ( aParticles%x0(3,vertIndices(k)) < minx0(3) ) then
							minx0(3) = aParticles%x0(3,vertIndices(k))
						endif
						if ( aParticles%x0(3,vertIndices(k)) > maxx0(3) ) then
							maxx0(3) = aParticles%x0(3,vertIndices(k))
						endif
					enddo
					lagVar = sum(maxx0-minx0)
					if ( lagVar > varTol ) then
						refineFlag(j) = .TRUE.
						amrCount2 = amrCount2 + 1
					endif
				endif
			endif
		endif
	enddo

	refineCount = count(refineFlag)
	amrLoopCounter = 1
	startIndex = 1
	if ( refineCount == 0 ) then
		call LogMessage(log,TRACE_LOGGING_LEVEL,logKey,'AMR converged for these tolerances.')
		deallocate(refineFlag)
		return
	else
		do while (refineCount > 0 )
			spaceLeft = aPanels%N_Max - aPanels%N
			if ( refineCount < spaceLeft/4) then
				if (printMsg) then
					write(logString,'(A,I4,A,I8,A)') 'AMR loop ',amrLoopCounter,' : refining ',refineCount,' panels.'
					call LogMessage(log,TRACE_LOGGING_LEVEL,logKey,trim(logString))
				endif

				nOldPanels = aPanels%N
				do j=startIndex,nOldPanels
					if ( refineFlag(j) ) then
						nOldPanels2 = aPanels%N
						nOldParticles = aParticles%N
						! divide the panel
						call DividePanel(aParticles,anEdges,aPanels,j)
						refineFlag(j) = .FALSE.

!!!!!
! SET VORTICITY OR TRACER
!!!!!
!			Set absolute vorticity on new particles	and panels
						if ( problemID == SINGLE_GAUSSIAN_VORTEX) then
							do k=nOldParticles+1,aParticles%N
								aParticles%absVort(k) = GaussianVortexX(aParticles%x0(:,k)) - gaussConst + &
									2.0_kreal*Omega*aParticles%x0(3,k)
							enddo
							do k=nOldPanels2+1,aPanels%N
								aPanels%absVort(k) = GaussianVortexX(aPanels%x0(:,k)) - gaussConst + &
										2.0_kreal*Omega*aPanels%x0(3,k)
							enddo
						elseif (problemID == RH4_WAVE ) then
							do k=nOldParticles+1,aParticles%N
								aParticles%absVort(k) = HaurwitzStationary4RelVort(aParticles%x0(:,k)) + &
									 2.0_kreal*Omega*aParticles%x0(3,k)
							enddo
							do k=nOldPanels2+1,aPanels%N
								aPanels%absVort(k) = HaurwitzStationary4RelVort(aPanels%x0(:,k)) + &
									2.0_kreal*Omega*aPanels%x0(3,k)
							enddo
						elseif ( problemID == MULTIPLE_VORTICES ) then
							do k=nOldParticles+1,aParticles%N
								aParticles%absVort(k) = ManyVortsX(aParticles%x0(:,k)) - gaussConst + &
									2.0_kreal*Omega*aParticles%x0(3,k)
							enddo
							do k=nOldPanels2+1,aPanels%N
								aPanels%absVort(k) = ManyVortsX(aPanels%x0(:,k)) - gaussConst + &
									2.0_kreal*Omega*aPanels%x0(3,k)
							enddo
						elseif ( problemID == RH2_WAVE ) then
							do k=nOldParticles+1,aParticles%N
								aParticles%absVort(k) = RH2WaveX(aParticles%x0(:,k)) + &
									2.0_kreal*Omega*aParticles%x0(3,k)
							enddo
							do k=nOldPanels2+1,aPanels%N
								aPanels%absVort(k) = RH2WaveX(aPanels%x0(:,k)) + &
									2.0_kreal*Omega*aPanels%x0(3,k)
							enddo
						elseif ( problemID == STRATOSPHERE_MODEL) then
							do k=nOldParticles+1,aParticles%N
								aParticles%absVort(k) = StratosphereRelVortX(aParticles%x0(:,k)) - gaussConst &
									+ 2.0_kreal*Omega*aParticles%x0(3,k) + Juckes_Forcing(aparticles%x0(:,k),zero)
							enddo
							do k=nOldPanels2+1,aPanels%N
								aPanels%absVort(k) = StratosphereRelVortX(aPanels%x0(:,k)) - gaussConst &
									+ 2.0_kreal*Omega*aPanels%x0(3,k) + Juckes_Forcing(aPanels%x0(:,k),zero)
							enddo
						elseif ( problemID == JET) then
							if ( (( .NOT. present(perturbWaveNum) ) .OR. (.NOT. present(beta))) .OR. &
								 ( (.NOT. present(theta0)) .OR. (.NOT. present(perturbAmp)) ) ) then
								 call LogMessage(log,ERROR_LOGGING_LEVEL,logkey,'Remesh Jet ERROR : JetRelVortX arguments needed.')
								 return
							endif
							do k=nOldParticles+1,aParticles%N
								aParticles%absVort(k) = JetRelVortX(aParticles%x0(:,k),theta0,beta,perturbAmp,perturbWaveNum) + &
									2.0_kreal*Omega*aParticles%x0(3,k) - gaussConst
							enddo
							do k=nOldPanels2+1,aPanels%N
								aPanels%absVort(k) = JetRelVortX(aPanels%x0(:,k),theta0,beta,perturbAmp,perturbWaveNum) + &
									2.0_kreal*Omega*aPanels%x0(3,k) - gaussConst
							enddo
						elseif ( problemID == ZONAL_MEAN ) then
							do k=nOldParticles+1,aParticles%N
								aParticles%absVort(k) = ZonalMeanX(aParticles%x0(:,k)) + &
									2.0_kreal*Omega*aParticles%x0(3,k)
							enddo
							do k=nOldPanels2+1,aPanels%N
								aPanels%absVort(k) = ZonalMeanX(aPanels%x0(:,k)) + &
									2.0_kreal*Omega*aPanels%x0(3,k)
							enddo
						elseif ( problemID == GAUSSIAN_HILLS) then
							do k=nOldParticles+1,aParticles%N
								aParticles%tracer(k,1) = GaussianHillsX(aParticles%x0(:,k))
							enddo
							do k=nOldPanels2+1,aPanels%N
								aPanels%tracer(k,1) = GaussianHillsX(aPanels%x0(:,k))
							enddo
						elseif ( problemID == COSINE_BELLS) then
							do k=nOldParticles+1,aParticles%N
								aParticles%tracer(k,1) = CosineBellsX(aParticles%x0(:,k))
							enddo
							do k=nOldPanels2+1,aPanels%N
								aPanels%tracer(k,1) = CosineBellsX(aPanels%x0(:,k))
							enddo
						elseif ( problemID == BLOCKM) then
							do k=nOldParticles+1,aParticles%N
								aParticles%tracer(k,1) = BlockMX(aParticles%x0(:,k))
							enddo
							do k=nOldPanels2+1,aPanels%N
								aPanels%tracer(k,1) = BlockMX(aPanels%x0(:,k))
							enddo
						elseif ( problemID == SLOTTED_CYLINDERS ) then
							do k=nOldParticles+1,aParticles%N
								aParticles%tracer(k,1) = SlottedCylindersX(aParticles%x0(:,k))
							enddo
							do k=nOldPanels2+1,aPanels%N
								aPanels%tracer(k,1) = SlottedCylindersX(aPanels%x0(:,k))
							enddo
						elseif (problemID == TWO_VORTS) then
							do k=nOldParticles+1,aParticles%N
								aParticles%absVort(k) = TwoVortsX(aParticles%x0(:,k),cent1,cent2,beta1,beta2,strength1,strength2) &
									- gaussConst + 2.0_kreal*Omega*aParticles%x0(3,j)
							enddo
							do k=nOldPanels2+1,aPanels%N
								aPanels%absVort(k) = TwoVortsX(aPanels%x0(:,k),cent1,cent2,beta1,beta2,strength1,strength2) &
									- gaussConst + 2.0_kreal*Omega*aPanels%x0(3,j)
							enddo
						elseif (problemID == TRIPOLE) then
							do k=nOldParticles+1,aParticles%N
								aParticles%absVort(k) = TripoleX(aParticles%x0(:,k),cent1,cent2,cent3,beta1,beta2,beta3,strength1,strength2,strength3) &
									- gaussConst + 2.0_kreal*Omega*aParticles%x0(3,j)
							enddo
							do k=nOldPanels2+1,aPanels%N
								aPanels%absVort(k) = TripoleX(aPanels%x0(:,k),cent1,cent2,cent3,beta1,beta2,beta3,strength1,strength2,strength3) &
									- gaussConst + 2.0_kreal*Omega*aPanels%x0(3,j)
							enddo
						endif



!					Set relative vorticity on new particles
						do k=nOldParticles+1,aParticles%N
							aParticles%relVort(k) = aParticles%absVort(k) - 2.0_kreal*Omega*aParticles%x(3,k)
						enddo
!					Set relative vorticity on new panels
						do k=nOldPanels2+1,aPanels%N
							aPanels%relVort(k) = aPanels%absVort(k) - 2.0_kreal*Omega*aPanels%x(3,k)
						enddo
!					Reset divided panel vorticity to zero
						aPanels%relVort(j) = 0.0_kreal
						aPanels%absVort(j) = 0.0_kreal
						aPanels%area(j) = 0.0_kreal
! 					Check refinement criteria on new panels
!!!!!
! AMR CRITERIA
!!!!!
						if ( aPanels%nest(j) <= initNest + MAX_REFINEMENT) then
							do k=nOldPanels2+1,aPanels%N
								if ( (problemID >= GAUSSIAN_HILLS) .AND. ( problemID <= SLOTTED_CYLINDERS) ) then
! ADVECTION
									if ( abs(aPanels%tracer(k,1))*aPanels%area(k) > circTol ) then
										refineFlag(k) = .TRUE.
										amrCount1 = amrCount1 + 1
									endif
								else
! BVE						! AMR criterion 1
									if ( abs(aPanels%relVort(k))*aPanels%area(k) > circTol ) then
										refineFlag(k) = .True.
										amrCount1 = amrCount1 + 1
									else
								! AMR Criterion 2
										minx0 = aPanels%x0(:,k)
										maxx0 = aPanels%x0(:,k)
										vertIndices(1:panelKind) = aPanels%vertices(1:panelKind,k)
										do i=1,panelKind
											if ( aParticles%x0(1,vertIndices(i)) < minX0(1) ) then
												minX0(1) = aParticles%x0(1,vertIndices(i))
											endif
											if ( aParticles%x0(1,vertIndices(i)) > maxX0(1) ) then
												maxx0(1) = aParticles%x0(1,vertIndices(i))
											endif
											if ( aparticles%x0(2,vertIndices(i)) < minX0(2)) then
												minx0(2) = aParticles%x0(2,vertIndices(i))
											endif
											if ( aParticles%x0(2,vertIndices(i)) > maxX0(2) ) then
												maxx0(2) = aParticles%x0(2,vertIndices(i))
											endif
											if ( aparticles%x0(3,vertIndices(i)) < minx0(3) ) then
												minx0(3) = aParticles%x0(3,vertIndices(i))
											endif
											if (aParticles%x0(3,vertIndices(I)) > maxx0(3)) then
												maxX0(3) = aParticles%x0(3,vertIndices(i))
											endif
										enddo
										lagVar = sum(maxx0-minx0)
										if ( lagVar > varTol ) then
											refineFlag(k) = .TRUE.
											amrCount2 = amrCount2 + 1
										endif
									endif
								endif
							enddo ! new panels refinement criteria
						else
							if ( .NOT. maxRefinement ) then
								call LogMessage(log,WARNING_LOGGING_LEVEL,logKey,'Refinement WARNING : maximum refinement reached.')
								maxRefinement = .TRUE.
							endif
						endif
					endif ! refineFlag(j)
				enddo  ! one-time panels loop

				! set up for next loop through panels (if necessary)
				refineCount = count(refineFlag)
				if ( refineCount == 0 ) then
					if ( printMsg ) &
					call LogMessage(log,TRACE_LOGGING_LEVEL,logKey,'InitRefine converged for these tolerances.')
					maxNest = maxVal(aPanels%nest)
					if ( panelKind == TRI_PANEL) then
						dlambda0 = TriAvgMeshSize(initNest)
						dlambda1 = TriAvgMeshSize(maxNest)
					else
						dlambda0 = QuadAvgMeshSize(initNest)
						dlambda1 = QuadAvgMeshSize(maxNest)
					endif
					if ( printMsg ) then
						call LogMessage(log,TRACE_LOGGING_LEVEL,trim(logKey)//' N_Active = ',aPanels%N_Active)
						call LogMessage(log,TRACE_LOGGING_LEVEL,trim(logKey)//' initNest dLambda = ', dlambda0)
						call LogMessage(log,TRACE_LOGGING_LEVEL,trim(logKey)//' maxNest dLambda = ',dLambda1)
						call LogMessage(log,TRACE_LOGGING_LEVEL,trim(logKey)//' circ. AMR count = ',amrCount1)
						call LogMessage(log,TRACE_LOGGING_LEVEL,trim(logKey)//' var. AMR count = ',amrCount2)
					endif
				else
					amrLoopCounter = amrLoopCounter + 1
					startIndex = nOldPanels
				endif
			else ! not enough space left
				call LogMessage(log,WARNING_LOGGING_LEVEL,logKey,'Not enough memory to continue AMR.')
				deallocate(refineFlag)
				return
			endif! spaceleft
		enddo ! while refineCount > 0
	endif
	deallocate(refineFlag)
end subroutine



subroutine AdaptiveRemesh(aParticles,anEdges,aPanels,initNest,AMR,amrTol1, amrTol2,procRank, problemID &
	,theta0, beta, perturbAmp, perturbWaveNum, time,cent1,cent2,cent3,beta1,beta2,beta3,&
	strength1,strength2,strength3, rotationRate )
	type(Particles), pointer, intent(inout) :: aParticles
	type(Edges), pointer, intent(inout) :: anEdges
	type(Panels), pointer, intent(inout) :: aPanels
	integer(kint), intent(in) :: initNest, AMR, procRank, problemID
	real(kreal), intent(in) :: amrTol1, amrTol2
	integer(kint), intent(in), optional :: perturbWaveNum	! Options for Jet test case
	real(kreal), intent(in), optional :: beta, theta0, perturbAmp	! Options for Jet test case
	real(kreal), intent(in), optional :: time ! Options for stratosphere model
	real(kreal), intent(in), optional :: cent1(3), cent2(3), cent3(3), &
										 beta1, beta2,beta3, &
										 strength1, strength2, strength3, &
										 rotationRate

	! STRIPACK
	type(Panels), pointer :: activePanels, passivePanels
	integer(kint), allocatable :: activeMap(:), passiveMap(:)
	integer(kint) :: nActive, nPassive, panelKind, nTracer
	integer(kint) :: errCode
	integer(kint) :: n
	real(kreal), allocatable :: x(:), y(:), z(:), dist(:)
	integer(kint), allocatable :: list(:), lptr(:), lend(:), near(:), next(:)
	integer(kint) :: ltri(1,1), nb, ncol, lnew
	! SSRFPACK
	real(kreal), allocatable :: x0(:), y0(:), z0(:)
	real(kreal), allocatable :: gradx0(:,:), grady0(:,:), gradz0(:,:)
	real(kreal), allocatable :: sigmax0(:), sigmay0(:), sigmaz0(:)
	real(kreal) :: sigmaTol, dSigX, dSigY, dSigZ
	! Remesh variables
	type(Particles), pointer :: newParticles, tempParticles
	type(Edges), pointer :: newEdges, tempEdges
	type(Panels), pointer :: newPanels, tempPanels
	logical(klog), allocatable :: refineFlag(:)
	integer(kint) :: refineCount, spaceLeft, nOldPanels, nOldParticles, nOldPanels2
	integer(kint) :: amrLoopCounter, startIndex, sigmaFlag, gradFlag, startPanel
	real(kreal) :: newLat, newLong, maxX0(3), minX0(3), lagVar
	integer(kint) :: vertIndices(4), amrCount1, amrCount2
	integer(kint) :: j, k, i
	logical(klog) :: maxRefinement

	call LogMessage(log,DEBUG_LOGGING_LEVEL,logkey,'Entering AdaptiveRemesh.')
	maxRefinement = .False.
	!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
	!	Part 1 : Setup existing mesh as source for data interpolation  !
	!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
	! Sift low-level panels from old mesh structures
	allocate(activePanels)
	allocate(passivePanels)
	nActive = aPanels%N_Active
	nPassive = aPanels%N - nActive
	panelKind = GetPanelKind(aPanels)
	nTracer = GetNTracer(aPanels)
	call New(activePanels,nActive,panelKind,nTracer,BVE_SOLVER)
	call New(passivePanels,nPassive,panelKind,nTracer,BVE_SOLVER)
	allocate(activeMap(nActive))
	allocate(passiveMap(nPassive))
	activePanels%N = nActive
	activePanels%N_Active = nActive
	passivePanels%N = nPassive
	call GatherPanels(aPanels,activePanels,activeMap,passivePanels,passiveMap)
	! Allocate variables for STRIPACK's Delaunay triangulation program
	n = nActive + aparticles%N
	!n = nActive
	allocate(x(n))
	allocate(y(n))
	allocate(z(n))
	allocate(dist(n))
	allocate(near(n))
	allocate(lend(n))
	allocate(next(n))
	allocate(list(6*n-12))
	allocate(lptr(6*n-12))
	ltri = 0
	nb = 0
	ncol = 0
	! Allocate variables for SSRFPACK
	allocate(x0(n))
	allocate(y0(n))
	allocate(z0(n))
	allocate(gradx0(3,n))
	allocate(grady0(3,n))
	allocate(gradz0(3,n))
	allocate(sigmax0(6*n-12))
	allocate(sigmay0(6*n-12))
	allocate(sigmaz0(6*n-12))
	x0(1:nActive) = activePanels%x0(1,:)
	y0(1:nActive) = activePanels%x0(2,:)
	z0(1:nActive) = activePanels%x0(3,:)
	x0(nActive+1:n) = aParticles%x0(1,1:aParticles%N)
	y0(nActive+1:n) = aParticles%x0(2,1:aParticles%N)
	z0(nActive+1:n) = aParticles%x0(3,1:aParticles%N)
	do j=1,nActive
		activePanels%x(:,j) = activePanels%x(:,j)/&
					sqrt(sum(activePanels%x(:,j)*activePanels%x(:,j)))
		x(j) = activePanels%x(1,j)
		y(j) = activePanels%x(2,j)
		z(j) = activePanels%x(3,j)
	enddo
	do j=1,aParticles%N
		aParticles%x(:,j) = aParticles%x(:,j)/&
				sqrt(sum(aParticles%x(:,j)*aParticles%x(:,j)))
		x(nActive+j) = aParticles%x(1,j)
		y(nActive+j) = aParticles%x(2,j)
		z(nActive+j) = aparticles%x(3,j)
	enddo

	! Build the Delaunay triangulation
	call TRMESH(n,x,y,z,list,lptr,lend,lnew,near,next,dist,errCode)
	if ( errCode > 0 ) then
		call LogMessage(Log,ERROR_LOGGING_LEVEL,trim(logKey)//" TRMESH found a duplicate node at node ",errCode)
		return
	endif

	! Estimate the gradients for each component at each node
	do j=1,n
		call GRADL(n,j,x,y,z,x0,list,lptr,lend,gradx0(:,j),errCode)
		call GRADL(n,j,x,y,z,y0,list,lptr,lend,grady0(:,j),errCode)
		call GRADL(n,j,x,y,z,z0,list,lptr,lend,gradz0(:,j),errCode)
	enddo

	! Compute the shape-preserving tension factors
	sigmaTol = 0.01_kreal
	sigmax0 = 0.0_kreal
	sigmay0 = 0.0_kreal
	sigmaz0 = 0.0_kreal
	call GETSIG(n,x,y,z,x0,list,lptr,lend,gradx0,sigmaTol,sigmax0,dSigX,errCode)
	call GETSIG(n,x,y,z,y0,list,lptr,lend,grady0,sigmaTol,sigmay0,dSigY,errCode)
	call GETSIG(n,x,y,z,z0,list,lptr,lend,gradz0,sigmaTOl,sigmaz0,dSigZ,errCode)


	!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
	!	Part 2 : Build new mesh and use as destination for interpolation  !
	!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
	call New(newParticles,newEdges,newPanels,panelKind,initNest,AMR,nTracer,BVE_SOLVER)
	allocate(refineFlag(newPanels%N_Max))
	refineFlag = .False.
	! Interpolate x0 from old mesh to new mesh
	gradFlag = 1
	sigmaFlag = 1
	startPanel = 1
	amrCount1 = 0
	amrCount2 = 0
	do j=1,newPanels%N
		if ( .NOT. newPanels%hasChildren(j) ) then
			newLat = Latitude(newPanels%x(:,j))
			newLong = Longitude(newPanels%x(:,j))
			call INTRC1(n,newLat,newLong,x,y,z,x0,list,lptr,lend,&
				sigmaFlag,sigmax0,gradFlag,gradx0,startpanel,newPanels%x0(1,j),errCode)
			call INTRC1(n,newLat,newLong,x,y,z,y0,list,lptr,lend,&
				sigmaFlag,sigmay0,gradFlag,grady0,startpanel,newPanels%x0(2,j),errCode)
			call INTRC1(n,newLat,newLong,x,y,z,z0,list,lptr,lend,&
				sigmaFlag,sigmaz0,gradFlag,gradz0,startpanel,newPanels%x0(3,j),errCode)
		endif
	enddo
	startPanel = 1
	do j=1,newParticles%N
		newLat = Latitude(newParticles%x(:,j))
		newLong = Longitude(newParticles%x(:,j))
		call INTRC1(n,newLat,newLong,x,y,z,x0,list,lptr,lend,&
			sigmaFlag,sigmax0,gradFlag,gradx0,startPanel,newParticles%x0(1,j),errCode)
		call INTRC1(n,newLat,newLong,x,y,z,y0,list,lptr,lend,&
			sigmaFlag,sigmay0,gradFlag,grady0,startPanel,newParticles%x0(2,j),errCode)
		call INTRC1(n,newLat,newLong,x,y,z,z0,list,lptr,lend,&
			sigmaFlag,sigmaz0,gradFlag,gradz0,startPanel,newParticles%x0(3,j),errCode)
	enddo
	! Normalize interpolated coordinates to unit sphere
	do j=1,newPanels%N
		if (.NOT. newPanels%hasChildren(j) ) then
			newPanels%x0(:,j) = newPanels%x0(:,j)/sqrt(sum( newPanels%x0(:,j)*newPanels%x0(:,j)))
		endif
	enddo
	do j=1,newParticles%N
		newParticles%x0(:,j) = newParticles%x0(:,j)/sqrt(sum( newParticles%x0(:,j)*newParticles%x0(:,j)))
	enddo
!!!!!
! SET VORTICITY And TRACER (for advection problems)
!!!!!
! Set absolute vorticity on new grid
	if ( problemID == SINGLE_GAUSSIAN_VORTEX) then
		do k=1,newParticles%N
			newParticles%absVort(k) = GaussianVortexX(newParticles%x0(:,k)) - gaussConst + &
				2.0_kreal*Omega*newParticles%x0(3,k)
		enddo
		do k=1,newPanels%N
			if ( .NOT. newPanels%hasChildren(k) ) then
				newPanels%absVort(k) = GaussianVortexX(newPanels%x0(:,k)) - gaussConst + &
					2.0_kreal*Omega*newPanels%x0(3,k)
			endif
		enddo
	elseif (problemID == RH4_WAVE ) then
		do k=1,newParticles%N
			newParticles%absVort(k) = HaurwitzStationary4RelVort(newParticles%x0(:,k)) + &
				 2.0_kreal*Omega*newParticles%x0(3,k)
		enddo
		do k=1,newPanels%N
			if ( .NOT. newPanels%hasCHildren(k) ) then
				newPanels%absVort(k) = HaurwitzStationary4RelVort(newPanels%x0(:,k)) + &
					2.0_kreal*Omega*newPanels%x0(3,k)
			endif
		enddo
	elseif ( problemID == MULTIPLE_VORTICES ) then
		do k=1,newParticles%N
			newParticles%absVort(k) = ManyVortsX(newParticles%x0(:,k)) - gaussConst + &
				2.0_kreal*Omega*newParticles%x0(3,k)
		enddo
		do k=1,newPanels%N
			if ( .NOT. newPanels%hasChildren(k) ) then
				newPanels%absVort(k) = ManyVortsX(newPanels%x0(:,k)) - gaussConst + &
					2.0_kreal*Omega*newPanels%x0(3,k)
			endif
		enddo
	elseif ( problemID == RH2_WAVE) then
		do k=1,newParticles%N
			newParticles%absVort(k) = RH2WaveX(newParticles%x0(:,k)) + &
					2.0_kreal*Omega*newParticles%x0(3,k)
		enddo
		do k=1,newPanels%N
			if ( .NOT. newPanels%hasChildren(k) ) then
				newPanels%absVort(k) = RH2WaveX(newPanels%x0(:,k)) + &
					2.0_kreal*Omega*newPanels%x0(3,k)
			endif
		enddo
	elseif ( problemID == STRATOSPHERE_MODEL ) then
		if ( .NOT. present(time) ) then
			call LogMessage(log,ERROR_LOGGING_LEVEL,logkey,'Remesh Stratosphere ERROR : time argument missing.')
			return
		endif
		do k=1,newParticles%N
			newParticles%absVort(k) = StratosphereRelVortX(newParticles%x0(:,k)) - gaussConst &
				+ 2.0_kreal*Omega*newParticles%x0(3,k)
		enddo
		do k=1,newPanels%N
			if ( .NOT. newPanels%hasChildren(k)) then
				newPanels%absVort(k) = StratosphereRelVortX(newPanels%x0(:,k)) - gaussConst &
					+ 2.0_kreal*Omega*newPanels%x0(3,k)
			 endif
		enddo
	elseif ( problemID == JET ) then
		if ( (( .NOT. present(perturbWaveNum) ) .OR. (.NOT. present(beta))) .OR. &
			 ( (.NOT. present(theta0)) .OR. (.NOT. present(perturbAmp)) ) ) then
			 call LogMessage(log,ERROR_LOGGING_LEVEL,logkey,'Remesh Jet ERROR : JetRelVortX arguments needed.')
			 return
		endif
		do k=1,newParticles%N
			newParticles%absVort(k) = JetRelVortX(newParticles%x0(:,k),theta0,beta,perturbAmp,perturbWaveNum) + &
				2.0_kreal*Omega*newParticles%x0(3,k) - gaussConst
		enddo
		do k=1,newPanels%N
			if (.NOT. newPanels%hasChildren(k) ) then
				newPanels%absVort(k) = JetRelVortX(newPanels%x0(:,k),theta0,beta,perturbAmp,perturbWaveNum) + &
					2.0_kreal*Omega*newPanels%x0(3,k) - gaussConst
			endif
		enddo
	elseif ( problemID == ZONAL_MEAN) then
		do k=1,newParticles%N
			newParticles%absVort(k) = ZonalMeanX(newParticles%x0(:,k)) + &
				2.0_kreal*Omega*newParticles%x0(3,k)
		enddo
		do k=1,newPanels%N
			if ( .NOT. newPanels%hasChildren(k) ) then
				newPanels%absVort(k) = ZonalMeanX(newPanels%x0(:,k)) + &
					2.0_kreal*Omega*newPanels%x0(3,k)
			endif
		enddo
	elseif ( problemID == GAUSSIAN_HILLS) then
		do k=1,newParticles%N
			newParticles%tracer(k,1) = GaussianHillsX(newParticles%x0(:,k))
		enddo
		do k=1,newPanels%N
			if (.NOT. newPanels%hasChildren(k) ) then
				newPanels%tracer(k,1) = GaussianHillsX(newPanels%x0(:,k))
			endif
		enddo
	elseif ( problemID == COSINE_BELLS) then
		do k=1,newParticles%N
			newParticles%tracer(k,1) = CosineBellsX(newParticles%x0(:,k))
		enddo
		do k=1,newPanels%N
			if (.NOT. newPanels%hasChildren(k) ) then
				newPanels%tracer(k,1) = CosineBellsX(newPanels%x0(:,k))
			endif
		enddo
	elseif ( problemID == BLOCKM) then
		do k=1,newParticles%N
			newParticles%tracer(k,1) = BlockMX(newParticles%x0(:,k))
		enddo
		do k=1,newPanels%N
			if (.NOT. newPanels%hasChildren(k) ) then
				newPanels%tracer(k,1) = BlockMX(newPanels%x0(:,k))
			endif
		enddo
	elseif ( problemID == SLOTTED_CYLINDERS) then
		do k=1,newParticles%N
			newParticles%tracer(k,1) = SlottedCylindersX(newParticles%x0(:,k))
		enddo
		do k=1,newPanels%N
			if (.NOT. newPanels%hasChildren(k) ) then
				newPanels%tracer(k,1) = SlottedCylindersX(newPanels%x0(:,k))
			endif
		enddo
	elseif (problemID == TWO_VORTS) then
		do k=1,newParticles%N
			newParticles%absVort(k) = TwoVortsX(newParticles%x0(:,k),cent1,cent2,beta1,beta2,strength1,strength2) &
				- gaussConst + 2.0_kreal*Omega*newParticles%x0(3,k)
		enddo
		do k=1,newPanels%N
			if ( .NOT. newPanels%hasChildren(k)) then
				newPanels%absVort(k) = TwoVortsX(newPanels%x0(:,k),cent1,cent2,beta1,beta2,strength1,strength2)&
					-gaussConst + 2.0_kreal*Omega*newPanels%x0(3,k)
			endif
		enddo
	elseif (problemID == TRIPOLE) then
		do k=1,newParticles%N
			newParticles%absVort(k) = TripoleX(newParticles%x0(:,k),cent1,cent2,cent3,&
				beta1,beta2,beta3,strength1,strength2,strength3) - gaussConst &
				+ 2.0_kreal*Omega*newParticles%x0(3,k)
		enddo
		do k=1,newPanels%N
			if (.NOT. newPanels%hasChildren(k)) then
				newPanels%absVort(k) = TripoleX(newPanels%x0(:,k),cent1,cent2,cent3,&
					beta1,beta2,beta3,strength1,strength2,strength3) - gaussConst &
					+ 2.0_kreal*Omega*newPanels%x0(3,k)
			endif
		enddo
	elseif (problemID == SOLID_BODY) then
		do k=1,newParticles%N
			newParticles%absVort(k) = 0.0_kreal
			newParticles%relVort(k) = SolidBodyX(newParticles%x0(:,k),rotationRate)
		enddo
		do k=1,newPanels%N
			if ( .NOT. newPanels%hasChildren(k) ) then
				newPanels%absVort(k) = 0.0_kreal
				newPanels%relVort(k) = SolidBodyX(newPanels%x0(:,k),rotationRate)
			endif
		enddo
	endif

! Set relative vorticity on new grid
	if ( problemID == STRATOSPHERE_MODEL ) then
		do j=1,newParticles%N
			newParticles%relVort(j) = newParticles%absVort(j) - 2.0_kreal*Omega*newParticles%x(3,j) &
				- Juckes_Forcing(newParticles%x(:,j),time)
		enddo
		do j=1,newPanels%N
			if ( .NOT. newPanels%hasChildren(j) ) then
				newPanels%relVort(j) = newPanels%absVort(j) - 2.0_kreal*Omega*newPanels%x(3,j) &
					- Juckes_Forcing(newPanels%x(:,j),time)
			endif
		enddo
	elseif ( problemID /= SOLID_BODY) then
		do j=1,newParticles%N
			newParticles%relVort(j) = newParticles%absVort(j) - 2.0_kreal*Omega*newParticles%x(3,j)
		enddo
		do j=1,newPanels%N
			if ( .NOT. newPanels%hasChildren(j) ) then
				newPanels%relVort(j) = newPanels%absVort(j) -2.0_kreal*Omega*newPanels%x(3,j)
			endif
		enddo
	endif
!!!!!!!!!
!  AMR  !
!!!!!!!!!
	if ( AMR > 0 ) then
		amrCount1 = 0
		amrCount2 = 0
		! Find panels that exceed refinement tolerances
		do j=1,newPanels%N
			if ( .NOT. newPanels%hasChildren(j) ) then
!!!!!
! AMR CRITERIA
!!!!!
				if ( (problemID >= GAUSSIAN_HILLS) .AND. ( problemID <= SLOTTED_CYLINDERS) ) then
! ADVECTION
					if ( abs(newPanels%tracer(j,1))*newPanels%area(j) > amrTol1 ) then
						refineFlag(j) = .TRUE.
						amrCount1 = amrCount1+1
					endif
				else
! BVE			! AMR Criterion 1
					if ( abs(newPanels%relVort(j))*newPanels%area(j) > amrTol1 ) then
						refineFlag(j) = .TRUE.
						amrCount1 = amrCount1 + 1
					else
					! AMR Criterion 2
						maxX0 = newPanels%x0(:,j)
						minX0 = newPanels%x0(:,j)
						vertIndices(1:panelKind) = newPanels%vertices(:,j)
						do k=1,panelKind
							if ( newParticles%x0(1,vertIndices(k)) < minX0(1) ) then
								minx0(1) = newParticles%x0(1,vertIndices(k))
							endif
							if ( newParticles%x0(2,vertIndices(k)) < minX0(2) ) then
								minx0(2) = newParticles%x0(2,vertIndices(k))
							endif
							if ( newParticles%x0(3,vertIndices(k)) < minX0(3) ) then
								minx0(3) = newParticles%x0(3,vertIndices(k))
							endif
							if ( newParticles%x0(1,vertIndices(k)) > maxx0(1) ) then
								maxx0(1) = newParticles%x0(1,vertIndices(k))
							endif
							if ( newParticles%x0(2,vertIndices(k)) > maxx0(2) ) then
								maxx0(2) = newParticles%x0(2,vertIndices(k))
							endif
							if ( newParticles%x0(3,vertIndices(k)) > maxx0(3) ) then
								maxx0(3) = newParticles%x0(3,vertIndices(k))
							endif
						enddo
						lagVar = sum(maxx0-minx0)
						if ( lagVar > amrTol2) then
							refineFlag(j) = .TRUE.
							amrCount2 = amrCount2 + 1
						endif
					endif
				endif
			endif
		enddo !initial panels loop

		refineCount = count(refineFlag)
		startIndex = 1
		amrLoopCounter = 1
		startPanel = 1

		! Refine grid until all panels meet refinement criteria
		do while (refineCount > 0)
			spaceLeft = newPanels%N_Max - newPanels%N
			if ( spaceLeft/4 > refineCount ) then
				if ( procRank <= 0 ) then
					write(logString,'(A,I4,A,I8,A)') " AMR Loop ",amrLoopCounter,&
						" : refining ",refineCount," panels."
					call LogMessage(log,TRACE_LOGGING_LEVEL,trim(logKey),logString)
				endif

				! Loop over grid, refine flagged panels
				nOldPanels = newPanels%N
				do j=startIndex,nOldPanels
					if ( refineFlag(j) ) then
						! Divide panel
						nOldParticles = newParticles%N
						nOldPanels2 = newPanels%N
						call DividePanel(newParticles,newEdges,newPanels,j)
						refineFlag(j) = .False.
						newPanels%absVort(j) = 0.0_kreal
						newPanels%relVort(j) = 0.0_kreal
						newPanels%area(j) = 0.0_kreal
						!if ( newParticles%N > nOldParticles) then
							! new particles have been added to the mesh-- find their x0
							do k=nOldParticles+1,newParticles%N
								newLat = Latitude(newParticles%x(:,k))
								newLong = Longitude(newParticles%x(:,k))
								call INTRC1(n,newLat,newLong,x,y,z,x0,list,lptr,lend,&
									sigmaFlag,sigmax0,gradFlag,gradx0,startPanel,&
									newParticles%x0(1,k),errCode)
								call INTRC1(n,newLat,newLong,x,y,z,y0,list,lptr,lend,&
									sigmaFlag,sigmay0,gradFlag,grady0,startPanel,&
									newParticles%x0(2,k),errCode)
								call INTRC1(n,newLat,newLong,x,y,z,z0,list,lptr,lend,&
									sigmaFlag,sigmaz0,gradFlag,gradz0,startPanel,&
									newParticles%x0(3,k),errCode)
								newParticles%x0(:,k) = newParticles%x0(:,k)/&
									sqrt(sum(newParticles%x0(:,k)*newParticles%x0(:,k)))
!!!!!
! SET VORTICITY					Set vorticity on new particles
!!!!!
								if ( problemID == SINGLE_GAUSSIAN_VORTEX) then
									newParticles%absVort(k) = GaussianVortexX(newParticles%x0(:,k)) - &
										gaussConst + 2.0_kreal*Omega*newParticles%x0(3,k)
								elseif (problemID == RH4_WAVE ) then
									newParticles%absVort(k) =  HaurwitzStationary4RelVort(newParticles%x0(:,k)) +&
												 2.0_kreal*Omega*newParticles%x0(3,k)
								elseif ( problemID == MULTIPLE_VORTICES ) then
									newParticles%absVort(k) = ManyVortsX(newParticles%x0(:,k)) - &
										gaussConst + 2.0_kreal*Omega*newParticles%x0(3,k)
								elseif ( problemID == RH2_Wave ) then
									newParticles%absVort(k) = RH2WaveX(newParticles%x0(:,k)) + &
										2.0_kreal*Omega*newParticles%x0(3,k)
								elseif ( problemID == STRATOSPHERE_MODEL) then
									newParticles%absVort(k) = StratosphereRelVortX(newParticles%x0(:,k)) - gaussConst &
										+ 2.0_kreal*Omega*newParticles%x0(3,k)
								elseif ( problemID == JET) then
									newParticles%absVort(k) = JetRelVortX(newParticles%x0(:,k),theta0,beta,perturbAmp,perturbWaveNum) + &
										2.0_kreal*Omega*newParticles%x0(3,k) - gaussConst
								elseif ( problemID == ZONAL_MEAN ) then
									newParticles%absVort(k) = ZonalMeanX(newParticles%x0(:,k)) + &
										2.0_kreal*Omega*newParticles%x0(3,k)
								elseif ( problemID == GAUSSIAN_HILLS) then
									newParticles%tracer(k,1) = GaussianHillsX(newParticles%x0(:,k))
								elseif ( problemID == COSINE_BELLS) then
									newParticles%tracer(k,1) = CosineBellsX(newParticles%x0(:,k))
								elseif ( problemID == BLOCKM) then
									newParticles%tracer(k,1) = BlockMX(newParticles%x0(:,k))
								elseif ( problemID == SLOTTED_CYLINDERS) then
									newParticles%tracer(k,1) = SlottedCylindersX(newParticles%x0(:,k))
								elseif (problemID == TWO_VORTS) then
									newParticles%absVort(k) = TwoVortsX(newParticles%x0(:,k),cent1,cent2,beta1,beta2,strength1,strength2)&
										- gaussConst + 2.0_kreal*Omega*newParticles%x0(3,k)
								elseif (problemID == TRIPOLE) then
									newParticles%absVort(k) = TripoleX(newParticles%x0(:,k),cent1,cent2,cent3,beta1,beta2,beta3,strength1,strength2,strength3)&
										- gaussConst + 2.0_kreal*Omega*newParticles%x0(3,k)
								endif

								if ( problemID == STRATOSPHERE_MODEL) then
									newParticles%relVort(k) = newParticles%absVort(k) - &
										2.0_kreal*Omega*newParticles%x(3,k) &
										- Juckes_Forcing(newParticles%x(:,k),time)
								else
									newParticles%relVort(k) = newParticles%absVort(k) - &
											2.0_kreal*Omega*newParticles%x(3,k)
								endif
							enddo
						!endif ! particles added
						! Find x0 of new panels
						do k=nOldPanels2+1,newPanels%N
							newLat = Latitude(newPanels%x(:,k))
							newLong = Longitude(newPanels%x(:,k))
							call INTRC1(n,newLat,newLong,x,y,z,x0,list,lptr,lend,&
								sigmaFlag,sigmax0,gradFlag,gradx0,startPanel,&
								newPanels%x0(1,k),errCode)
							call INTRC1(n,newLat,newLong,x,y,z,y0,list,lptr,lend,&
								sigmaFlag,sigmay0,gradFlag,grady0,startPanel,&
								newPanels%x0(2,k),errCode)
							call INTRC1(n,newLat,newLong,x,y,z,z0,list,lptr,lend,&
								sigmaFlag,sigmaz0,gradFlag,gradz0,startPanel,&
								newPanels%x0(3,k),errCode)
							newPanels%x0(:,k) = newPanels%x0(:,k)/&
								sqrt(sum(newPanels%x0(:,k)*newPanels%x0(:,k)))
!!!!!
! SET VORTICITY				Set vorticity on new panels
!!!!!
							if ( problemID == SINGLE_GAUSSIAN_VORTEX) then
								newPanels%absVort(k) = GaussianVortexX(newPanels%x0(:,k)) - &
									gaussConst + 2.0_kreal*Omega*newPanels%x0(3,k)
							elseif (problemID == RH4_WAVE ) then
								newPanels%absVort(k) = HaurwitzStationary4RelVort(newPanels%x0(:,k)) + &
									2.0_kreal*Omega*newPanels%x0(3,k)
							elseif ( problemID == MULTIPLE_VORTICES ) then
								newPanels%absVort(k) = ManyVortsX(newPanels%x0(:,k)) - gaussConst + &
									2.0_kreal*Omega*newPanels%x0(3,k)
							elseif ( problemID == RH2_WAVE ) then
								newPanels%absVort(k) = RH2WaveX(newPanels%x0(:,k)) + &
									2.0_kreal*Omega*newPanels%x0(3,k)
							elseif ( problemID == STRATOSPHERE_MODEL ) then
								newPanels%absVort(k) = StratosphereRelVortX(newPanels%x0(:,k)) - gaussConst &
									+ 2.0_kreal*Omega*newPanels%x0(3,k)
							elseif ( problemID == JET ) then
								newPanels%absVort(k) = JetRelVortX(newPanels%x0(:,k),theta0,beta,perturbAmp,perturbWaveNum) + &
									2.0_kreal*Omega*newPanels%x0(3,k) - gaussConst
							elseif ( problemID == ZONAL_MEAN ) then
								newPanels%absVort(k) = ZonalMeanX(newPanels%x0(:,k)) + &
									2.0_kreal*Omega*newPanels%x0(3,k)
							elseif ( problemID == GAUSSIAN_HILLS) then
								newPanels%tracer(k,1) = GaussianHillsX(newPanels%x0(:,k))
							elseif ( problemID == COSINE_BELLS) then
								newPanels%tracer(k,1) = CosineBellsX(newPanels%x0(:,k))
							elseif ( problemID == BLOCKM) then
								newPanels%tracer(k,1) = BlockMX(newPanels%x0(:,k))
							elseif ( problemID == SLOTTED_CYLINDERS) then
								newPanels%tracer(k,1) = SlottedCylindersX(newPanels%x0(:,k))
							elseif (problemID == TWO_VORTS) then
								newPanels%absVort(k) = TwoVortsX(newPanels%x0(:,k),cent1,cent2,beta1,beta2,strength1,strength2)&
									- gaussConst + 2.0_kreal*Omega*newPanels%x0(3,k)
							elseif (problemID == TRIPOLE) then
								newPanels%absVort(k) = TripoleX(newPanels%x0(:,k),cent1,cent2,cent3,beta1,beta2,beta3,strength1,strength2,strength3) &
									- gaussConst + 2.0_kreal*Omega*newPanels%x0(3,k)
							endif

							if ( problemID == STRATOSPHERE_MODEL ) then
								newPanels%relVort(k) = newPanels%absVort(k) - &
										2.0_kreal*Omega*newPanels%x(3,k) &
										- Juckes_Forcing(newPanels%x(:,k),time)
							else
								newPanels%relVort(k) = newPanels%absVort(k) - &
										2.0_kreal*Omega*newPanels%x(3,k)
							endif

!!!!!
! AMR CRITERIA				Check refinement criteria on new panels
!!!!!
							if ( newPanels%nest(j) <= initNest + MAX_REFINEMENT) then
								if ( (problemID >= GAUSSIAN_HILLS) .AND. ( problemID <= SLOTTED_CYLINDERS) ) then
! ADVECTION
									if ( abs(newPanels%tracer(k,1))*newPanels%area(k) > amrTol1 ) then
										refineFlag(k) = .TRUE.
										amrCount1 = amrCount1 + 1
									endif
								else
! BVE
									if ( abs(newPanels%relVort(k))*newPanels%area(k) > amrTol1 ) then
										refineFlag(k) = .TRUE.
										amrCount1 = amrCount1 + 1
									else
										minx0 = newPanels%x0(:,k)
										maxx0 = newPanels%x0(:,k)
										vertIndices(1:panelKind) = newPanels%vertices(:,k)
										do i=1,panelKind
											if ( newParticles%x0(1,vertIndices(i)) < minx0(1) ) then
												minx0(1) = newParticles%x0(1,vertIndices(i))
											endif
											if ( newParticles%x0(2,vertIndices(i)) < minx0(2) ) then
												minx0(2) = newParticles%x0(2,vertIndices(i))
											endif
											if ( newParticles%x0(3,vertIndices(i)) < minx0(3) ) then
												minx0(3) = newParticles%x0(3,vertIndices(i))
											endif
											if ( newParticles%x0(1,vertIndices(i)) > maxx0(1) ) then
												maxx0(1) = newParticles%x0(1,vertIndices(i))
											endif
											if ( newParticles%x0(2,vertIndices(i)) > maxx0(2) ) then
												maxx0(2) = newParticles%x0(2,vertIndices(i))
											endif
											if ( newParticles%x0(3,vertIndices(i)) > maxx0(3) ) then
												maxx0(3) = newParticles%x0(3,vertIndices(i))
											endif
										enddo

										lagVar = sum(maxx0-minx0)
										if ( lagVar > amrTol2 ) then
											refineFlag(k) = .TRUE.
											amrCount2 = amrCount2 + 1
										endif
									endif
								endif
							else
								if ( .NOT. maxRefinement) then
									if ( procRank == 0 ) call LogMessage(log,WARNING_LOGGING_LEVEL,logKey,'Refinement WARNING : maximum refinement reached.')
									maxRefinement = .TRUE.
								endif
							endif
						enddo
					endif ! flagged panel
				enddo ! panels loop

				! Reset for next level of refinement
				refineCount = count(refineFlag)
				amrLoopCounter = amrLoopCounter + 1
				startIndex = nOldPanels

			else ! not enough memory for AMR
				 if ( procRank <= 0 ) call LogMessage(log,WARNING_LOGGING_LEVEL,trim(logKey),' WARNING: not enough memory for AMR.')
				 exit
			endif !spaceleft
		enddo ! while refineCount > 0

		! Log statistics
		if ( procRank <= 0) then
			write(logString,'(A,I8,A)') " AMR circulation criterion triggered ",amrCount1," times."
			call LogMessage(log,TRACE_LOGGING_LEVEL,trim(logKey),trim(logString))
			write(logString,'(A,I8,A)') " AMR Lag. coord. variation criterion triggered ",amrCount2," times."
			call LogMessage(log,TRACE_LOGGING_LEVEL,trim(logKey),trim(logString))
		endif
	endif ! AMR > 0

	!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
	!	Part 3 : Replace old grid with new grid  !
	!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

	tempParticles => aParticles
	tempEdges => anEdges
	tempPanels => aPanels

	aParticles => newParticles
	anEdges => newEdges
	aPanels => newPanels

	!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
	!	Delete old grid & Clean up	    !
	!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

	call Delete(tempParticles,tempEdges,tempPanels)
	nullify(newParticles)
	nullify(newEdges)
	nullify(newPanels)

	deallocate(refineFlag)
	deallocate(sigmaz0)
	deallocate(sigmay0)
	deallocate(sigmax0)
	deallocate(gradz0)
	deallocate(grady0)
	deallocate(gradx0)
	deallocate(x0)
	deallocate(y0)
	deallocate(z0)
	deallocate(x)
	deallocate(y)
	deallocate(z)
	deallocate(list)
	deallocate(lptr)
	deallocate(lend)
	deallocate(near)
	deallocate(next)
	deallocate(dist)
	deallocate(activeMap)
	deallocate(passiveMap)
	call Delete(activePanels)
	call Delete(passivePanels)
	deallocate(activePanels)
	deallocate(passivePanels)
	nullify(activePanels)
	nullify(passivePanels)
end subroutine


subroutine DirectAdaptiveRemesh(aParticles,anEdges,aPanels,initNest,AMR,amrTol1, amrTol2,procRank,problemID)
	type(Particles), pointer, intent(inout) :: aParticles
	type(Edges), pointer, intent(inout) :: anEdges
	type(Panels), pointer, intent(inout) :: aPanels
	integer(kint), intent(in) :: initNest, AMR
	real(kreal), intent(in) :: amrTol1, amrTol2
	integer(kint), intent(in) :: procRank, problemID
	! STRIPACK
	integer(kint) :: n, errCode
	real(kreal), allocatable :: x(:), y(:), z(:), dist(:)
	integer(kint), allocatable :: list(:), lptr(:), lend(:), near(:), next(:)
	integer(kint) :: ltri(1,1), nb, ncol, lnew
	! SSRFPACK
	real(kreal), allocatable :: absVort(:), tracer(:,:)
	real(kreal), allocatable :: gradAbsVort(:,:), gradTracer(:,:,:)
	real(kreal), allocatable :: sigmaAbsVort(:), sigmaTracer(:,:)
	real(kreal) :: sigmaTol, dSigAbsVort, dSigTracer
	integer(kint) :: sigmaFlag, gradFlag, startPanel
	! Remesh variables
	type(Panels), pointer :: activePanels, passivePanels
	integer(kint), allocatable :: activeMap(:), passiveMap(:)
	type(Particles), pointer :: newParticles, tempParticles
	type(Edges), pointer :: newEdges, tempEdges
	type(Panels), pointer :: newPanels, tempPanels
	logical(klog), allocatable :: refineFlag(:)
	integer(kint) :: refineCount, spaceLeft, nOldPanels, nOldPanels2, nOldParticles
	integer(kint) :: amrLoopCounter, amrCount1, amrCount2, startIndex
	integer(kint) :: nActive, nPassive, panelKind, nTracer, vertIndices(4)
	real(kreal) :: newLat, newLong, maxRelVort, minRelVort,  lagVar
	integer(kint) :: i,j,k
	logical(klog) :: maxRefinement

	if ( .NOT. logInit) call InitLogger(log)
	call LogMessage(log,DEBUG_LOGGING_LEVEL,logKey,"Entering DirectAdaptiveRemesh");
	maxRefinement = .FALSE.

	!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
	!	Part 1 : Setup existing mesh as source for data interpolation  !
	!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
	! Sift low-level panels from old mesh structures
	allocate(activePanels)
	allocate(passivePanels)
	nActive = aPanels%N_Active
	nPassive = aPanels%N - nActive
	panelKind = GetPanelKind(aPanels)
	nTracer = GetNTracer(aPanels)
	call New(activePanels,nActive,panelKind,nTracer,BVE_SOLVER)
	call New(passivePanels,nPassive,panelKind,nTracer,BVE_SOLVER)
	allocate(activeMap(nActive))
	allocate(passiveMap(nPassive))
	activePanels%N = nActive
	activePanels%N_Active = nActive
	passivePanels%N = nPassive
	call GatherPanels(aPanels,activePanels,activeMap,passivePanels,passiveMap)
	! Allocate variables for STRIPACK's Delaunay triangulation program
	n = nActive + aparticles%N
	allocate(x(n))
	allocate(y(n))
	allocate(z(n))
	allocate(dist(n))
	allocate(near(n))
	allocate(lend(n))
	allocate(next(n))
	allocate(list(6*n-12))
	allocate(lptr(6*n-12))
	ltri = 0
	nb = 0
	ncol = 0
	! Allocate variables for SSRFPACK
	allocate(absVort(n))
	allocate(tracer(n,nTracer))
	allocate(gradAbsVort(3,n))
	allocate(gradTracer(3,n,nTracer))
	allocate(sigmaAbsVort(6*n-12))
	allocate(sigmaTracer(6*n-12,nTracer))
	!!call LogMessage(log,DEBUG_LOGGING_LEVEL,trim(logKey)//" DirectRemesh :"," Allocation done.")

	! Reproject to sphere for STRIPACK
	do i=1,nActive
		activePanels%x(:,i) = activePanels%x(:,i)/sqrt(sum(activePanels%x(:,i)*activePanels%x(:,i)))
		x(i) = activePanels%x(1,i)
		y(i) = activePanels%x(2,i)
		z(i) = activePanels%x(3,i)
	enddo
	do i=1,aParticles%N
		aParticles%x(:,i) = aParticles%x(:,i)/sqrt(sum(aParticles%x(:,i)*aParticles%x(:,i)))
		x(nActive+i) = aParticles%x(1,i)
		y(nActive+i) = aParticles%x(2,i)
		z(nActive+i) = aParticles%x(3,i)
	enddo
	! Build the Delaunay triangulation
	call TRMESH(n,x,y,z,list,lptr,lend,lnew,near,next,dist,errCode)
	if ( errCode > 0 ) then
		call LogMessage(Log,ERROR_LOGGING_LEVEL,trim(logKey)//" TRMESH found a duplicate node at node ",errCode)
		return
	endif
	call LogMessage(log,DEBUG_LOGGING_LEVEL,trim(logKey)//" DirectRemesh :"," TRMESH Returned.")
	! Setup source data
	absVort(1:nActive) = activePanels%absVort
	absVort(nActive+1:n) = aParticles%absVort(1:aParticles%N)
	do i=1,nTracer
		tracer(1:nActive,i) = activePanels%tracer(:,i)
		tracer(nActive+1:n,i) = aParticles%tracer(1:aParticles%N,i)
	enddo
	! Estimate the gradients at each node
	do j=1,n
		call GRADL(n,j,x,y,z,absVort,list,lptr,lend,gradAbsVort(:,j),errCode)
	enddo
	do i=1,nTracer
		do j=1,n
			call GRADL(n,j,x,y,z,tracer(:,i),list,lptr,lend,gradTracer(:,j,i),errCode)
		enddo
	enddo
	call LogMessage(log,DEBUG_LOGGING_LEVEL,trim(logKey)//" DirectRemesh :"," Gradients done.")
	! Compute the shape-preserving tension factors
	sigmaTol = 0.01_kreal
	sigmaAbsVort = 0.0_kreal
	sigmaTracer = 0.0_kreal
	call GETSIG(n,x,y,z,absVort,list,lptr,lend,gradAbsVort,sigmaTol,sigmaAbsVort,dSigAbsVort,errCode)
	do i=1,nTracer
		call GETSIG(n,x,y,z,tracer(:,i),list,lptr,lend,gradTracer(:,:,i),sigmaTol,sigmaTracer(:,i),dSigTracer,errCode)
	enddo
	call LogMessage(log,DEBUG_LOGGING_LEVEL,trim(logKey)//" DirectRemesh :"," sigmas done.")
	!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
	!	Part 2 : Build new mesh and use as destination for interpolation  !
	!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
	call New(newParticles,newEdges,newPanels,panelKind,initNest,AMR,nTracer,BVE_SOLVER)
	allocate(refineFlag(newPanels%N_Max))
	!call LogMessage(log,DEBUG_LOGGING_LEVEL,trim(logKey)//" DirectRemesh :"," new uniform mesh returned.")
	refineFlag = .False.
	! Interpolate from old mesh to new mesh
	gradFlag = 1
	sigmaFlag = 1
	startPanel = 1
	amrCount1 = 0
	amrCount2 = 0
	do j=1,newPanels%N
		if ( .NOT. newPanels%hasChildren(j) ) then
			newLat = Latitude(newPanels%x(:,j))
			newLong = Longitude(newPanels%x(:,j))
			call INTRC1(n,newLat,newLong,x,y,z,absVort,list,lptr,lend,&
				sigmaFlag,sigmaAbsVort,gradFlag,gradAbsVort,startPanel,newPanels%absVort(j),errCode)
			do i=1,nTracer
				call INTRC1(n,newLat,newLong,x,y,z,tracer(:,i),list,lptr,lend,&
					sigmaFlag,sigmaTracer(:,i),gradFlag,gradTracer(:,:,i),startPanel,newPanels%tracer(j,i),errCode)
			enddo
		endif
	enddo
	call LogMessage(log,DEBUG_LOGGING_LEVEL,trim(logKey)//" DirectRemesh :"," newPanels initial interp done.")
	startPanel = 1
	do j=1,newParticles%N
		newLat = Latitude(newParticles%x(:,j))
		newLong = Longitude(newParticles%x(:,j))
		call INTRC1(n,newLat,newLong,x,y,z,absVort,list,lptr,lend,&
			sigmaFlag,sigmaAbsVort,gradFlag,gradAbsVort,startPanel,newParticles%absVort(j),errCode)
		do i=1,nTracer
			call INTRC1(n,newLat,newLong,x,y,z,tracer(:,i),list,lptr,lend,&
				sigmaFlag,sigmaTracer(:,i),gradFlag,gradTracer(:,:,i),startPanel,newParticles%tracer(j,i),errCode)
		enddo
	enddo
	call LogMessage(log,DEBUG_LOGGING_LEVEL,trim(logKey)//" DirectRemesh :"," newParticles initial interp done.")
	! Set relative vorticity on new grid
	do j=1,newPanels%N
		if ( .NOT. newPanels%hasChildren(j) ) then
			newPanels%relVort(j) = newPanels%absVort(j) - 2.0_kreal*Omega*newPanels%x(3,j)
		endif
	enddo
	do j=1,newParticles%N
		newParticles%relVort(j) = newParticles%absVort(j) - 2.0_kreal*Omega*newParticles%x(3,j)
	enddo
	!call LogMessage(log,DEBUG_LOGGING_LEVEL,trim(logKey)//" DirectRemesh :"," relative vorticity done.")
!!!!!!!!!
!  AMR  !
!!!!!!!!!
	if ( AMR > 0 ) then
		amrCount1 = 0
		amrCount2 = 0
		! Find panels that exceed refinement tolerances
		do j=1,newPanels%N
			if ( (problemID >= GAUSSIAN_HILLS) .AND. ( problemID <= SLOTTED_CYLINDERS) ) then
			! Advection AMR
				if ( abs(newPanels%tracer(j,1))*newPanels%area(j) > amrTol1 ) then
					refineFlag(j) = .TRUE.
					amrCount1 = amrCount1 + 1
				endif
			else
			! BVE AMR
				if ( abs(newPanels%relVort(j))*newPanels%area(j) > amrTol1 ) then
					refineFlag(j) = .TRUE.
					amrCount1 = amrCount1 + 1
				else
					maxRelVort = newPanels%relVort(j)
					minRelVort = newPanels%relVort(j)
					vertIndices(1:panelKind) = newPanels%vertices(:,j)
					do k=1,panelKind
						if ( newParticles%relVort(vertIndices(k)) < minRelVort ) &
							minRelVort = newParticles%relVort(vertIndices(k))
						if ( newParticles%relVort(vertIndices(k)) > maxRelVort ) &
							maxRelVort = newParticles%relVort(vertIndices(k))
					enddo
					lagVar = maxRelVort - minRelVort
					if ( lagVar > amrTol2) then
						refineFlag(j) = .TRUE.
						amrCount2 = amrCount2 + 1
					endif
				endif
			endif ! advection or bve
		enddo ! original panels loop

		refineCount = count(refineFlag)
		startIndex = 1
		amrLoopCounter = 1
		startPanel = 1

		! refine grid until all panels meet criteria or until refinement limits or memory limits are hit
		do while (refineCount > 0 )
			spaceLeft = newPanels%N_Max - newPanels%N
			if ( spaceleft/4 > refineCount ) then
				if ( procRank == 0 ) then
					write(logString,'(A,I4,A,I8,A)') " AMR Loop ",amrLoopCounter,&
						" : refining ",refineCount," panels."
					call LogMessage(log,TRACE_LOGGING_LEVEL,trim(logKey),logString)
				endif

				! Loop over grid, refine flagged panels
				nOldPanels = newPanels%N
				do j=startIndex,nOldPanels
					if ( refineFlag(j) ) then
						! Divide the panel
						nOldParticles = newParticles%N
						nOldPanels2 = newPanels%N
						call DividePanel(newParticles,newEdges,newPanels,j)
						refineFlag(j) = .FALSE.
						newPanels%absVort(j) = 0.0_kreal
						newPanels%relVort(j) = 0.0_kreal
						newPanels%area(j) = 0.0_kreal

						if ( GetLoggingLevel(log) == DEBUG_LOGGING_LEVEL) then
							write(logString,'(A,I5,A,I8,A,F12.8)') "j = ",j, " Na = ",newPanels%N_Active," surf. area = ",sum(newPanels%area(1:newPanels%N))
							call LogMessage(log,DEBUG_LOGGING_LEVEL,logkey,logString)
						endif

						do k=nOldParticles+1,newParticles%N
							newLat = Latitude(newParticles%x(:,k))
							newLong = Longitude(newParticles%x(:,k))
							call INTRC1(n,newLat,newLong,x,y,z,absVort,list,lptr,lend,&
								sigmaFlag,sigmaAbsVort,gradFlag,gradAbsVort,startPanel,&
								newParticles%absVort(k),errCode)
							newParticles%relVort(k) = newParticles%absVort(k) - 2.0_kreal*Omega*newParticles%x(3,k)
							do i=1,nTracer
								call INTRC1(n,newLat,newLong,x,y,z,tracer(:,i),list,lptr,lend,&
									sigmaFlag,sigmaTracer(:,i),gradFlag,gradTracer(:,:,i),startPanel,&
									newParticles%tracer(k,i),errCode)
							enddo
						enddo! new particles

						do k=nOldPanels2+1,newPanels%N
							newLat = Latitude(newPanels%x(:,k))
							newLong = Longitude(newPanels%x(:,k))
							call INTRC1(n,newLat,newLong,x,y,z,absVort,list,lptr,lend,&
								sigmaFlag,sigmaAbsVort,gradFlag,gradAbsVort,startPanel,&
								newPanels%absVort(k),errCode)
							newPanels%relVort(k) = newPanels%absVort(k) - 2.0_kreal*Omega*newPanels%x(3,k)
							do i=1,nTracer
								call INTRC1(n,newLat,newLong,x,y,z,tracer(:,i),list,lptr,lend,&
									sigmaFlag,sigmaTracer(:,i),gradFlag,gradTracer(:,:,i),startPanel,&
									newPanels%tracer(k,i),errCode)
							enddo

							! Check refinement criteria on new panels
							if ( newPanels%nest(j) <= initNest + MAX_REFINEMENT ) then
								if ( (problemID <= GAUSSIAN_HILLS ) .AND. (problemID <= SLOTTED_CYLINDERS) ) then
								! advection amr
									if ( abs(newPanels%tracer(k,1))*newPanels%area(k) > amrTol1) then
										refineFlag(k) = .TRUE.
										amrCount1 = amrCount1 + 1
									endif
								else
								! bve amr
									if ( abs(newPanels%relVort(k))*newPanels%area(k) > amrTol1 ) then
										refineFlag(k) = .TRUE.
										amrCount1 = amrCount1 + 1
									else
										minRelVort = newPanels%relVort(k)
										maxRelVort = newPanels%relVort(k)
										vertIndices(1:panelKind) = newPanels%vertices(:,k)
										do i=1,panelKind
											if ( newParticles%relVort(vertIndices(i)) < minRelVort ) &
												minRelVort = newParticles%relVort(vertIndices(i))
											if ( newParticles%relVort(vertIndices(i)) > maxRelVort ) &
												maxRelVort = newParticles%relVort(vertIndices(i))
										enddo
										lagVar = maxRelVort - minRelVort
										if ( lagVar > amrTol2 ) then
											refineFlag(k) = .TRUE.
											amrCount2 = amrCount2 + 1
										endif! amrTol2
									endif! bve criteria
								endif! advection or bve
							else
							! max refinement reached
								if ( .NOT. maxRefinement ) then
									if ( procRank == 0 ) then
										call LogMessage(log,WARNING_LOGGING_LEVEL,logKey,"Refinement AMR : maximum refinement reached.")
									endif
									maxRefinement = .TRUE.
								endif
							endif! <= max_refinement
						enddo! new Panels
					endif ! refineFlag(j)
				enddo! panels loop

				! Reset for next iteration of while loop
				refineCount = count(refineFlag)
				amrLoopCounter = amrLoopCounter + 1
				startIndex = nOldPanels
			else
			! not enough memory to continue AMR
				if ( procRank == 0 ) then
					call LogMessage(log,WARNING_LOGGING_LEVEL,trim(logKey)," WARNING: not enough memory for AMR.")
				endif
				exit
			endif! space left
		enddo !while refineCount > 0
		! Log statistics
		if ( procRank == 0) then
			write(logString,'(A,I8,A)') " AMR circulation criterion triggered ",amrCount1," times."
			call LogMessage(log,TRACE_LOGGING_LEVEL,trim(logKey),trim(logString))
			write(logString,'(A,I8,A)') " AMR Lag. coord. variation criterion triggered ",amrCount2," times."
			call LogMessage(log,TRACE_LOGGING_LEVEL,trim(logKey),trim(logString))
		endif
	endif ! AMR > 0

	!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
	!	Part 3 : Replace old grid with new grid  !
	!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

	tempParticles => aParticles
	tempEdges => anEdges
	tempPanels => aPanels

	aParticles => newParticles
	anEdges => newEdges
	aPanels => newPanels

	if (GetLoggingLevel(log) == DEBUG_LOGGING_LEVEL) then
		call LogMessage(log,DEBUG_LOGGING_LEVEL,logKey,'OLD MESH : ');
		call PrintStats(tempParticles)
		call PrintStats(tempEdges)
		call PrintStats(tempPanels)
		call LogMessage(log,DEBUG_LOGGING_LEVEL,logKey,'NEW MESH : ');
		call PrintStats(aParticles)
		call PrintStats(anEdges)
		call PrintStats(aPanels)
	endif

	!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
	!	Delete old grid & Clean up	    !
	!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
	call LogMessage(log,DEBUG_LOGGING_LEVEL,trim(logKey)//" DirectRemesh :"," Starting cleanup.")
	call Delete(tempParticles,tempEdges,tempPanels)
	nullify(newParticles)
	nullify(newEdges)
	nullify(newPanels)
	nullify(tempParticles)
	nullify(tempEdges)
	nullify(tempPanels)
	call LogMessage(log,DEBUG_LOGGING_LEVEL,trim(logKey)//" DirectRemesh :"," old grid deleted.")
	deallocate(refineFlag)
	deallocate(absVort)
	deallocate(tracer)
	deallocate(sigmaTracer)
	deallocate(sigmaAbsVort)
	deallocate(gradTracer)
	deallocate(gradAbsVort)
	call LogMessage(log,DEBUG_LOGGING_LEVEL,trim(logKey)//" DirectRemesh :"," break 1.")
	deallocate(lptr)
	deallocate(list)
	deallocate(next)
	deallocate(lend)
	deallocate(near)
	deallocate(dist)
	call LogMessage(log,DEBUG_LOGGING_LEVEL,trim(logKey)//" DirectRemesh :"," break 2.")
	deallocate(z)
	deallocate(y)
	deallocate(x)
	deallocate(passiveMap)
	deallocate(activeMap)
	call LogMessage(log,DEBUG_LOGGING_LEVEL,trim(logKey)//" DirectRemesh :"," break 3.")
	call Delete(activePanels)
	call Delete(passivePanels)
	!call LogMessage(log,DEBUG_LOGGING_LEVEL,trim(logKey)//" DirectRemesh :"," break 4.")
	deallocate(activePanels)
	deallocate(passivePanels)
	call LogMessage(log,DEBUG_LOGGING_LEVEL,trim(logKey)//" DirectRemesh :"," break 5.")
	nullify(activePanels)
	nullify(passivePanels)
	call LogMessage(log,DEBUG_LOGGING_LEVEL,trim(logKey)//" DirectRemesh :"," break 6.")
end subroutine


subroutine InitLogger(aLog)
	type(Logger), intent(inout) :: aLog
	integer(kint) :: logUnit
	logUnit = 6
	call New(aLog,logLevel,logUnit)
	logInit = .True.
end subroutine

end module
