module PanelVelocityParallelModule

use NumberKindsModule
use SphereGeomModule
use LoggerModule
use ParticlesEdgesPanelsModule
use TracerAndVorticityDistributionModule

implicit none

include 'mpif.h'

private
public totalKE
public kineticEnergy 	! KE at each panel
public BVERK4, AdvectionRK4, DivergentAdvectionRK4, StratosphereRK4 ! timestepping routines
public ResetRK4	 ! memory free routine
public ResetArea ! Set to true for divergent flow
public SetSmooth, GetSmooth
public InitializeMPIRK4, FinalizeMPIRK4, LoadBalance

!----------------
! Module variables
!----------------
real(kreal), allocatable, save :: kineticEnergy(:)
real(kreal), save :: totalKE
real(kreal), save :: smooth = 0.05_kreal
logical(klog), save :: resetAreas = .False.
logical(klog), save :: rk4IsReady = .False.

type(Panels), pointer, save :: activePanels=>null(), passivePanels=>null()
integer(kint), allocatable, save :: activeMap(:), passiveMap(:)

real(kreal), allocatable, save :: particlesInput(:,:),&
								  particlesStage1(:,:),&
								  particlesStage2(:,:),&
								  particlesStage3(:,:),&
								  particlesStage4(:,:),&
								  newParticlesX(:,:),&
!								  particlesVortInput(:), &
!								  particlesVortStage1(:), &
!								  particlesVortStage2(:), &
!								  particlesVortStage3(:), &
!								  particlesVortStage4(:), &
!								  newParticlesVort(:), &
								  activePanelsInput(:,:),&
								  activePanelsStage1(:,:),&
								  activePanelsStage2(:,:),&
								  activePanelsStage3(:,:),&
								  activePanelsStage4(:,:),&
								  newActivePanelsX(:,:),&
								  activeVortInput(:),&
								  activeVortStage1(:),&
								  activeVortStage2(:),&
								  activeVortStage3(:),&
								  activeVortStage4(:),&
								  newActiveVort(:),&
								  passivePanelsInput(:,:),&
								  passivePanelsStage1(:,:),&
								  passivePanelsStage2(:,:),&
								  passivePanelsStage3(:,:),&
								  passivePanelsStage4(:,:),&
								  newPassivePanelsX(:,:),&
!								  passiveVortTrash1(:),&
!								  passiveVortTrash2(:),&
								  area(:),&
								  areaStage1(:),&
								  areaStage2(:),&
								  areaStage3(:),&
								  areaStage4(:),&
								  newArea(:)

!----------------
! MPI Variables
!----------------
logical(klog), save :: mpiIsReady = .False.
integer(kint) :: mpiStatus(MPI_STATUS_SIZE), mpiErrorCode
integer(kint), parameter :: nullTag = 0
integer(kint), allocatable :: particlesIndexStart(:), particlesIndexEnd(:), particlesMessageSize(:), &
							  activePanelsIndexStart(:), activePanelsIndexEnd(:), activePanelsMessageSize(:),&
							  passivePanelsIndexStart(:), passivePanelsIndexEnd(:), passivePanelsMessageSize(:)

!----------------
! Logging
!----------------
logical(klog), save :: logInit = .False.
integer(kint), parameter :: logLevel = TRACE_LOGGING_LEVEL, logUnit = 6
character(len=28) :: logKey = 'PanelVelocityMPI : '
type(Logger) :: log

contains
!----------------
! User decisions
!----------------

subroutine SetSmooth(newSmooth)
	real(kreal), intent(in) :: newSmooth
	smooth = newSmooth
end subroutine


function GetSmooth()
	real(kreal) :: GetSmooth
	GetSmooth = smooth
end function


function GetOmega()
	real(kreal) :: GetOmega
	GetOmega = Omega
end function

subroutine ResetArea(onOff)
	logical(klog), intent(in), optional :: onOff
	if ( present(onOff)) then
		resetAreas = onOff
	else
		resetAreas = .True.
	endif
end subroutine

!----------------
! MPI Setup
!----------------

subroutine InitializeMPIRK4(aParticles,aPanels,procRank,numProcs)
	type(Particles), intent(in) :: aParticles
	type(Panels), intent(in) :: aPanels
	integer(kint), intent(in) :: procRank, numProcs

	if ( .NOT. logInit) then
		call InitLogger(log)
		write(logKey,'(A,I1,A)') 'PanelVelocityMPI ',procRank,' : '
	endif

	call LogMessage(log,DEBUG_LOGGING_LEVEL,trim(logKey)//" hello from process ",procRank)

	allocate(particlesIndexStart(0:numProcs-1))
	allocate(particlesIndexEnd(0:numProcs-1))
	allocate(particlesMessageSize(0:numProcs-1))
	particlesIndexStart = 0
	particlesIndexEnd = 0
	particlesMessageSize = 0

	allocate(activePanelsIndexStart(0:numProcs-1))
	allocate(activePanelsIndexEnd(0:numProcs-1))
	allocate(activePanelsMessageSize(0:numProcs-1))
	activePanelsIndexStart = 0
	activePanelsIndexEnd = 0
	activePanelsMessageSize = 0

	allocate(passivePanelsIndexStart(0:numProcs-1))
	allocate(passivePanelsIndexEnd(0:numProcs-1))
	allocate(passivePanelsMessageSize(0:numProcs-1))
	passivePanelsIndexStart = 0
	passivePanelsIndexEnd = 0
	passivePanelsMessageSize = 0

	mpiIsReady = .True.

	!call LoadBalance(aParticles,aPanels,numProcs)

end subroutine


subroutine LoadBalance(aParticles,aPanels,numProcs)
	type(Particles), intent(in) :: aParticles
	type(Panels), intent(in) :: aPanels
	integer(kint), intent(in) :: numProcs
	integer(kint) :: chunkSize, j
	integer(kint) :: nParticles, nActive, nPassive

	call LogMessage(log,DEBUG_LOGGING_LEVEL,logKey,'Entering LoadBalance.')

	nParticles = aParticles%N
	nActive = aPanels%N_Active
	nPassive = aPanels%N - aPanels%N_Active

	if (mpiIsReady ) then

		chunkSize = nParticles/numProcs
		do j=0,numProcs-1
			particlesIndexStart(j) = j*chunkSize + 1
			particlesIndexEnd(j) = (j+1)*chunkSize
		enddo
		particlesIndexEnd(numProcs-1) = nParticles
		particlesMessageSize = particlesIndexEnd - particlesIndexStart + 1

		chunkSize = nActive/numProcs
		do j=0,numProcs-1
			activePanelsIndexStart(j) = j*chunkSize + 1
			activePanelsIndexEnd(j) = (j+1)*chunkSize
		enddo
		activePanelsIndexEnd(numProcs-1) = nActive
		activePanelsMessageSize = activePanelsIndexEnd - activePanelsIndexStart + 1

		chunkSize = nPassive/numProcs
		do j=0,numProcs-1
			passivePanelsIndexStart(j) = j*chunkSize + 1
			passivePanelsIndexEnd(j) = (j+1)*chunkSize
		enddo
		passivePanelsIndexEnd(numProcs-1) = nPassive
		passivePanelsMessageSize = passivePanelsIndexEnd - passivePanelsIndexStart + 1

		!call LogMessage(log,DEBUG_LOGGING_LEVEL,logKey,'LoadBalance complete.')
	else

		call LogMessage(log,ERROR_LOGGING_LEVEL,logKey,'LoadBalance ERROR : MPI arrays not ready.')
		return

	endif
end subroutine


subroutine FinalizeMPIRK4()
	if (rk4isReady)	call ResetRK4()
	deallocate(particlesIndexStart)
	deallocate(particlesIndexEnd)
	deallocate(particlesMessageSize)
	deallocate(activepanelsIndexStart)
	deallocate(activePanelsIndexEnd)
	deallocate(activePanelsMessageSize)
	deallocate(passivePanelsIndexStart)
	deallocate(passivePanelsIndexend)
	deallocate(passivePanelsMessageSize)
	mpiIsReady = .False.
	if ( logInit) call Delete(log)
end subroutine



!----------------
! Memory handling
!----------------
subroutine InitRK4(aParticles,aPanels)
	type(Particles), intent(in) :: aParticles
	type(Panels), intent(in) :: aPanels
	integer(kint) :: nPassive, nActive, nParticles, nTracer, panelKind, problemKind

	call LogMessage(log,DEBUG_LOGGING_LEVEL,logKey,'Entering InitRK4.')

	nParticles = aParticles%N
	nActive = aPanels%N_Active
	nPassive = aPanels%N-aPanels%N_Active
	nTracer = GetNTracer(aPanels)
	panelKind = GetPanelKind(aPanels)
	problemKind = GetProblemKind(aPanels)

	allocate(activePanels)
	allocate(passivePanels)
	allocate(activeMap(nActive))
	allocate(passiveMap(nPassive))
	activeMap = 0
	passiveMap = 0
	call New(activePanels,nActive,panelKind,nTracer,problemKind)
	call New(passivePanels,nPassive,panelKind,nTracer,problemKind)
	call LogMessage(log,DEBUG_LOGGING_LEVEL,logKey,'InitRK4: active/passive panels ready.')
	activePanels%N = nActive
	activePanels%N_Active = nActive
	passivePanels%N = nPassive
	call GatherPanels(aPanels,activePanels,activeMap,passivePanels,passiveMap)
	call LogMessage(log,DEBUG_LOGGING_LEVEL,logKey,'InitRK4: Gather panels done.')

	allocate(particlesInput(3,nParticles))
	allocate(particlesStage1(3,nParticles))
	allocate(particlesStage2(3,nParticles))
	allocate(particlesStage3(3,nParticles))
	allocate(particlesStage4(3,nParticles))
	allocate(newParticlesX(3,nParticles))

	!allocate(particlesVortInput(nParticles))
!	allocate(particlesVortStage1(nParticles))
!	allocate(particlesVortStage2(nParticles))
!	allocate(particlesVortStage3(nParticles))
!	allocate(particlesVortStage4(nParticles))
!	allocate(newParticlesVort(nParticles))

	allocate(activePanelsInput(3,nActive))
	allocate(activePanelsStage1(3,nActive))
	allocate(activePanelsStage2(3,nActive))
	allocate(activePanelsStage3(3,nActive))
	allocate(activePanelsStage4(3,nActive))
	allocate(newActivePanelsX(3,nActive))

	allocate(passivePanelsInput(3,nPassive))
	allocate(passivePanelsStage1(3,nPassive))
	allocate(passivePanelsStage2(3,nPassive))
	allocate(passivePanelsStage3(3,nPassive))
	allocate(passivePanelsStage4(3,nPassive))
	allocate(newPassivePanelsX(3,nPassive))

!	allocate(passiveVortTrash1(nPassive))
!	allocate(passiveVortTrash2(nPassive))

	allocate(activeVortInput(nActive))
	allocate(activeVortStage1(nActive))
	allocate(activeVortStage2(nActive))
	allocate(activeVortStage3(nActive))
	allocate(activeVortStage4(nActive))
	allocate(newActiveVort(nActive))

	allocate(kineticEnergy(nActive))

	allocate(area(nActive))
	if ( resetAreas ) then
		allocate(areaStage1(nActive))
		allocate(areaStage2(nActive))
		allocate(areaStage3(nActive))
		allocate(areaStage4(nActive))
		allocate(newArea(nActive))
	endif

	rk4isReady = .TRUE.

	call LogMessage(log,DEBUG_LOGGING_LEVEL,trim(logKey),' RK4 memory available.')

end subroutine


subroutine ResetRK4()
	deallocate(particlesInput)
	deallocate(particlesStage1)
	deallocate(particlesStage2)
	deallocate(particlesStage3)
	deallocate(particlesStage4)
	deallocate(newParticlesX)

	!deallocate(particlesVortInput)
!	deallocate(particlesVortStage1)
!	deallocate(particlesVortStage2)
!	deallocate(particlesVortStage3)
!	deallocate(particlesVortStage4)
!	deallocate(newParticlesVort)

	deallocate(activePanelsInput)
	deallocate(activePanelsStage1)
	deallocate(activePanelsStage2)
	deallocate(activePanelsStage3)
	deallocate(activePanelsStage4)
	deallocate(newActivePanelsX)

	deallocate(passivePanelsInput)
	deallocate(passivePanelsStage1)
	deallocate(passivePanelsStage2)
	deallocate(passivePanelsStage3)
	deallocate(passivePanelsStage4)
	deallocate(newPassivePanelsX)

!	deallocate(passiveVortTrash1)
!	deallocate(passiveVortTrash2)

	deallocate(activeVortInput)
	deallocate(activeVortStage1)
	deallocate(activeVortStage2)
	deallocate(activeVortStage3)
	deallocate(activeVortStage4)
	deallocate(newActiveVort)

	deallocate(kineticEnergy)

	deallocate(area)
	if (allocated(areaStage1)) then
		deallocate(areaStage1)
		deallocate(areaStage2)
		deallocate(areaStage3)
		deallocate(areaStage4)
		deallocate(newArea)
	endif

	deallocate(activeMap)
	deallocate(passiveMap)
	call Delete(activePanels)
	call Delete(passivePanels)
	deallocate(activePanels)
	deallocate(passivePanels)
	rk4isReady = .False.

	call LogMessage(log,DEBUG_LOGGING_LEVEL,trim(logkey),' RK4 memory freed.')
end subroutine


subroutine ZeroRK4()
	particlesInput = 0.0_kreal
	particlesStage1 = 0.0_kreal
	particlesStage2 = 0.0_kreal
	particlesStage3 = 0.0_kreal
	particlesStage4 = 0.0_kreal
	newParticlesX = 0.0_kreal

!	particlesVortInput = 0.0_kreal
!	particlesVortStage1 = 0.0_kreal
!	particlesVortStage2 = 0.0_kreal
!	particlesVortStage3 = 0.0_kreal
!	particlesVortStage4 = 0.0_kreal
!	newParticlesVort = 0.0_kreal

	activePanelsInput = 0.0_kreal
	activePanelsStage1 = 0.0_kreal
	activePanelsStage2 = 0.0_kreal
	activePanelsStage3 = 0.0_kreal
	activePanelsStage4 = 0.0_kreal
	newActivePanelsX = 0.0_kreal

	passivePanelsInput = 0.0_kreal
	passivePanelsStage1 = 0.0_kreal
	passivePanelsStage2 = 0.0_kreal
	passivePanelsStage3 = 0.0_kreal
	passivePanelsStage4 = 0.0_kreal
	newPassivePanelsX = 0.0_kreal

	activeVortInput = 0.0_kreal
	activeVortStage1 = 0.0_kreal
	activeVortStage2 = 0.0_kreal
	activeVortStage3 = 0.0_kreal
	activeVortStage4 = 0.0_kreal
	newActiveVort = 0.0_kreal

	kineticEnergy = 0.0_kreal

	area = 0.0_kreal
	if (allocated(areaStage1) ) then
		areaStage1 = 0.0_kreal
		areaStage2 = 0.0_kreal
		areaStage3 = 0.0_kreal
		areaStage4 = 0.0_kreal
		newArea = 0.0_kreal
	endif
end subroutine

!----------------
! Timestepping
!----------------

subroutine BVERK4(aParticles,aPanels,dt,procRank,numProcs)
	type(Particles), intent(inout) :: aParticles
	type(Panels), intent(inout) :: aPanels
	real(kreal), intent(in) :: dt
	integer(kint), intent(in) :: procRank, numProcs
	integer(kint) :: j, nTracer

	call LogMessage(log,DEBUG_LOGGING_LEVEL,logKey,'Entering RK4')

	if ( rk4isReady ) then
		call ZeroRK4()
	else
		call InitRK4(aParticles,aPanels)
		call LogMessage(log,DEBUG_LOGGING_LEVEL,logKey,'Returned from InitRK4.')
		call LoadBalance(aParticles,aPanels,numProcs)
	endif

	nTracer = GetNTracer(apanels)
	call LogMessage(log,DEBUG_LOGGING_LEVEL,logKey,'Starting RK Stage 1')
	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! 		RK Stage 1       !
	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Set input arrays for stage 1
	particlesInput = aParticles%x(:,1:aParticles%N)
	activePanelsInput = activePanels%x
	passivePanelsInput = passivePanels%x
	area = activePanels%area
	activeVortInput = activePanels%relvort
	! Store velocity output in stage1 arrays
	! Portion work over processes
	call BVEActiveRHS(activePanelsStage1,activeVortStage1,& !output
					  activePanelsInput,activeVortInput, area, & !input
					  activePanelsIndexStart(procRank),activePanelsIndexEnd(procRank)) !indices
	call BVEPassiveRHS(particlesStage1,& !output
					   particlesInput, activePanelsInput,activeVortInput, area, & !input
					   particlesIndexStart(procRank), particlesIndexEnd(procRank)) ! indices
	call BVESmoothVelocity( passivePanelsStage1, & !output
							passivePanelsInput, activePanelsInput, activeVortInput, area, & ! input
							passivePanelsIndexStart(procRank), passivePanelsIndexEnd(procRank)) ! indices

	! broadcast results to all
	do j=0,numProcs-1
		call MPI_BCAST(activePanelsStage1(:,activePanelsIndexStart(j):activePanelsIndexEnd(j)),& ! send buffer
					   3*activePanelsMessageSize(j),MPI_DOUBLE_PRECISION,j, & ! message data and senderId
					   MPI_COMM_WORLD, mpiErrorCode)
		call MPI_BCAST(activeVortStage1(activePanelsIndexStart(j):activePanelsIndexEnd(j)),&
					   activePanelsMessageSize(j),MPI_DOUBLE_PRECISION,j, &
					   MPI_COMM_WORLD, mpiErrorCode)
		call MPI_BCAST(particlesStage1(:,particlesIndexStart(j):particlesIndexEnd(j)),&
					   3*particlesMessageSize(j),MPI_DOUBLE_PRECISION, j, &
					   MPI_COMM_WORLD, mpiErrorCode)
		call MPI_BCAST(passivePanelsStage1(:,passivePanelsIndexStart(j):passivePanelsIndexEnd(j)),&
					   3*passivePanelsMessageSize(j), MPI_DOUBLE_PRECISION, j, &
					   MPI_COMM_WORLD, mpiErrorCode)
	enddo


	!!!!
	!! Kinetic energy calculation (stage 1 only)
	!!!!
	do j=1,activePanels%N
		kineticEnergy(j) = sum(activePanelsStage1(:,j)*activePanelsStage1(:,j))
	enddo
	totalKE = 0.5_kreal*sum(kineticEnergy*area)
	if ( nTracer >= 2) then
		activePanels%tracer(:,2) = kineticEnergy
	endif

	activePanelsStage1 = dt*activePanelsStage1
	particlesStage1 = dt*particlesStage1
	passivePanelsStage1 = dt*passivePanelsStage1
	activeVortStage1 = dt*activeVortStage1

	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! 		RK Stage 2       !
	!!!!!!!!!!!!!!!!!!!!!!!!!!
	call LogMessage(log,DEBUG_LOGGING_LEVEL,logKey,'Starting RK Stage 2')
	! Set input arrays for stage 2
	particlesInput = aParticles%x(:,1:aParticles%N) + 0.5_kreal*particlesStage1
	activePanelsInput = activePanels%x + 0.5_kreal*activePanelsStage1
	passivePanelsINput = passivePanels%x + 0.5_kreal*passivePanelsStage1
	activeVortInput = activePanels%relVort + 0.5_kreal*activeVortStage1
	! Store velocity output in stage2 arrays
	! Portion work over processes
	call BVEActiveRHS(activePanelsStage2,activeVortStage2,& !output
					  activePanelsInput,activeVortInput, area, & !input
					  activePanelsIndexStart(procRank),activePanelsIndexEnd(procRank)) !indices
	call BVEPassiveRHS(particlesStage2,& !output
					   particlesInput, activePanelsInput,activeVortInput, area, & !input
					   particlesIndexStart(procRank), particlesIndexEnd(procRank)) ! indices
	call BVESmoothVelocity( passivePanelsStage2, & !output
							passivePanelsInput, activePanelsInput, activeVortInput, area, & ! input
							passivePanelsIndexStart(procRank), passivePanelsIndexEnd(procRank)) ! indices

	! broadcast results to all
	do j=0,numProcs-1
		call MPI_BCAST(activePanelsStage2(:,activePanelsIndexStart(j):activePanelsIndexEnd(j)),& ! send buffer
					   3*activePanelsMessageSize(j),MPI_DOUBLE_PRECISION,j, & ! message data and senderId
					   MPI_COMM_WORLD, mpiErrorCode)
		call MPI_BCAST(activeVortStage2(activePanelsIndexStart(j):activePanelsIndexEnd(j)),&
					   activePanelsMessageSize(j),MPI_DOUBLE_PRECISION,j, &
					   MPI_COMM_WORLD, mpiErrorCode)
		call MPI_BCAST(particlesStage2(:,particlesIndexStart(j):particlesIndexEnd(j)),&
					   3*particlesMessageSize(j),MPI_DOUBLE_PRECISION, j, &
					   MPI_COMM_WORLD, mpiErrorCode)
		call MPI_BCAST(passivePanelsStage2(:,passivePanelsIndexStart(j):passivePanelsIndexEnd(j)),&
					   3*passivePanelsMessageSize(j), MPI_DOUBLE_PRECISION, j, &
					   MPI_COMM_WORLD, mpiErrorCode)
	enddo

	activePanelsStage2 = dt*activePanelsStage2
	particlesStage2 = dt*particlesStage2
	passivePanelsStage2 = dt*passivePanelsStage2
	activeVortStage2 = dt*activeVortStage2

	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! 		RK Stage 3       !
	!!!!!!!!!!!!!!!!!!!!!!!!!!
	call LogMessage(log,DEBUG_LOGGING_LEVEL,logKey,'Starting RK Stage 3')
	! Set input arrays for stage 3
	particlesInput = aParticles%x(:,1:aParticles%N) + 0.5_kreal*particlesStage2
	activePanelsInput = activePanels%x + 0.5_kreal*activePanelsStage2
	passivePanelsINput = passivePanels%x + 0.5_kreal*passivePanelsStage2
	activeVortInput = activePanels%relVort + 0.5_kreal*activeVortStage2
	! Store velocity output in Stage3 arrays
	! Portion work over processes
	call BVEActiveRHS(activePanelsStage3,activeVortStage3,& !output
					  activePanelsInput,activeVortInput, area, & !input
					  activePanelsIndexStart(procRank),activePanelsIndexEnd(procRank)) !indices
	call BVEPassiveRHS(particlesStage3,& !output
					   particlesInput, activePanelsInput,activeVortInput, area, & !input
					   particlesIndexStart(procRank), particlesIndexEnd(procRank)) ! indices
	call BVESmoothVelocity( passivePanelsStage3, & !output
							passivePanelsInput, activePanelsInput, activeVortInput, area, & ! input
							passivePanelsIndexStart(procRank), passivePanelsIndexEnd(procRank)) ! indices

	! broadcast results to all
	do j=0,numProcs-1
		call MPI_BCAST(activePanelsStage3(:,activePanelsIndexStart(j):activePanelsIndexEnd(j)),& ! send buffer
					   3*activePanelsMessageSize(j),MPI_DOUBLE_PRECISION,j, & ! message data and senderId
					   MPI_COMM_WORLD, mpiErrorCode)
		call MPI_BCAST(activeVortStage3(activePanelsIndexStart(j):activePanelsIndexEnd(j)),&
					   activePanelsMessageSize(j),MPI_DOUBLE_PRECISION,j, &
					   MPI_COMM_WORLD, mpiErrorCode)
		call MPI_BCAST(particlesStage3(:,particlesIndexStart(j):particlesIndexEnd(j)),&
					   3*particlesMessageSize(j),MPI_DOUBLE_PRECISION, j, &
					   MPI_COMM_WORLD, mpiErrorCode)
		call MPI_BCAST(passivePanelsStage3(:,passivePanelsIndexStart(j):passivePanelsIndexEnd(j)),&
					   3*passivePanelsMessageSize(j), MPI_DOUBLE_PRECISION, j, &
					   MPI_COMM_WORLD, mpiErrorCode)
	enddo
	activePanelsStage3 = dt*activePanelsStage3
	particlesStage3 = dt*particlesStage3
	passivePanelsStage3 = dt*passivePanelsStage3
	activeVortStage3 = dt*activeVortStage3


	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! 		RK Stage 4       !
	!!!!!!!!!!!!!!!!!!!!!!!!!!
	call LogMessage(log,DEBUG_LOGGING_LEVEL,logKey,'Starting RK Stage 4')
	! Set input arrays for stage 4
	particlesInput = aParticles%x(:,1:aParticles%N) + particlesStage3
	activePanelsInput = activePanels%x + activePanelsStage3
	passivePanelsINput = passivePanels%x + passivePanelsStage3
	activeVortInput = activePanels%relVort + activeVortStage3
	! Store velocity output in Stage4 arrays
	! Portion work over processes
	call BVEActiveRHS(activePanelsStage4,activeVortStage4,& !output
					  activePanelsInput,activeVortInput, area, & !input
					  activePanelsIndexStart(procRank),activePanelsIndexEnd(procRank)) !indices
	call BVEPassiveRHS(particlesStage4,& !output
					   particlesInput, activePanelsInput,activeVortInput, area, & !input
					   particlesIndexStart(procRank), particlesIndexEnd(procRank)) ! indices
	call BVESmoothVelocity( passivePanelsStage4, & !output
							passivePanelsInput, activePanelsInput, activeVortInput, area, & ! input
							passivePanelsIndexStart(procRank), passivePanelsIndexEnd(procRank)) ! indices

	! broadcast results to all
	do j=0,numProcs-1
		call MPI_BCAST(activePanelsStage4(:,activePanelsIndexStart(j):activePanelsIndexEnd(j)),& ! send buffer
					   3*activePanelsMessageSize(j),MPI_DOUBLE_PRECISION,j, & ! message data and senderId
					   MPI_COMM_WORLD, mpiErrorCode)
		call MPI_BCAST(activeVortStage4(activePanelsIndexStart(j):activePanelsIndexEnd(j)),&
					   activePanelsMessageSize(j),MPI_DOUBLE_PRECISION,j, &
					   MPI_COMM_WORLD, mpiErrorCode)
		call MPI_BCAST(particlesStage4(:,particlesIndexStart(j):particlesIndexEnd(j)),&
					   3*particlesMessageSize(j),MPI_DOUBLE_PRECISION, j, &
					   MPI_COMM_WORLD, mpiErrorCode)
		call MPI_BCAST(passivePanelsStage4(:,passivePanelsIndexStart(j):passivePanelsIndexEnd(j)),&
					   3*passivePanelsMessageSize(j), MPI_DOUBLE_PRECISION, j, &
					   MPI_COMM_WORLD, mpiErrorCode)
	enddo
	activePanelsStage4 = dt*activePanelsStage4
	particlesStage4 = dt*particlesStage4
	passivePanelsStage4 = dt*passivePanelsStage4
	activeVortStage4 = dt*activeVortStage4

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! 		RK Position update   !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	call LogMessage(log,DEBUG_LOGGING_LEVEL,logKey,'Starting RK Update')
	newParticlesX = aParticles%x(:,1:aParticles%N) + particlesStage1/6.0_kreal +&
		particlesStage2/3.0_kreal + particlesStage3/3.0_kreal + particlesStage4/6.0_kreal
	newActivePanelsX = activePanels%x + activePanelsStage1/6.0_kreal + activePanelsStage2/3.0_kreal +&
		activePanelsStage3/3.0_kreal + activePanelsStage4/6.0_kreal
	newPassivePanelsX = passivePanels%x + passivePanelsStage1/6.0_kreal + passivePanelsStage2/3.0_kreal +&
		passivePanelsStage3/3.0_kreal + passivePanelsStage4/6.0_kreal
	newActiveVort = activePanels%relVort + activeVortStage1/6.0_kreal + activeVortStage2/3.0_kreal + &
		activeVortStage3/3.0_kreal + activeVortStage4/6.0_kreal

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! 	Copy to data structures   !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	call LogMessage(log,DEBUG_LOGGING_LEVEL,logKey,'Starting RK Copy')
	aParticles%x(:,1:aParticles%N) = newParticlesX
	aParticles%relVort(1:aParticles%N) = aParticles%absVort(1:aParticles%N) - 2.0_kreal*Omega*aParticles%x(3,1:aParticles%N)
	activePanels%x = newActivePanelsX
	activePanels%relVort = newActiveVort
	passivePanels%x = newPassivePanelsX
	call LogMessage(log,DEBUG_LOGGING_LEVEL,logKey,'Calling Scatter Panels')
	call ScatterPanels(aPanels,activePanels,activeMap,passivePanels,passiveMap)
	call LogMessage(log,DEBUG_LOGGING_LEVEL,logKey,'Returned from Scatter Panels')
end subroutine


subroutine StratosphereRK4(aParticles,aPanels,t,dt,procRank,numProcs)
	! Arguments
	type(Particles), intent(inout) :: aParticles
	type(Panels), intent(inout) :: aPanels
	real(kreal), intent(in) :: t, dt
	integer(kint), intent(in) :: procRank, numProcs
	! Locals
	integer(kint) :: nTracer, j

	call LogMessage(log,DEBUG_LOGGING_LEVEL,logKey,'Entering StratosphereRK4')

	if ( rk4isReady ) then
		call ZeroRK4()
	else
		call InitRK4(aParticles,aPanels)
		call LoadBalance(aParticles,aPanels,numProcs)
	endif

	nTracer = GetNTracer(apanels)

	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! 		RK Stage 1       !
	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Set input arrays for stage 1
	particlesInput = aParticles%x(:,1:aParticles%N)
	!particlesVortInput = aParticles%relVort(1:aParticles%N)
	activePanelsInput = activePanels%x
	passivePanelsInput = passivePanels%x
	area = activePanels%area
	!activeVortInput = activePanels%relvort
	do j=1,aPanels%N_Active
		activeVortInput(j) = activePanels%absVort(j) - 2.0_kreal*Omega*activePanelsInput(3,j) &
					- Juckes_Forcing(activePanelsInput(:,j),t)
	enddo

	! Store velocity output in stage1 arrays
	! Portion work over processes
	call BVEActiveRHS(activePanelsStage1,activeVortStage1,& !output
					  activePanelsInput,activeVortInput, area, & !input
					  activePanelsIndexStart(procRank),activePanelsIndexEnd(procRank)) !indices
	call BVEPassiveRHS(particlesStage1,& !output
					   particlesInput, activePanelsInput,activeVortInput, area, & !input
					   particlesIndexStart(procRank), particlesIndexEnd(procRank)) ! indices
	call BVESmoothVelocity( passivePanelsStage1, & !output
							passivePanelsInput, activePanelsInput, activeVortInput, area, & ! input
							passivePanelsIndexStart(procRank), passivePanelsIndexEnd(procRank)) ! indices

	! broadcast results to all
	do j=0,numProcs-1
		call MPI_BCAST(activePanelsStage1(:,activePanelsIndexStart(j):activePanelsIndexEnd(j)),& ! send buffer
					   3*activePanelsMessageSize(j),MPI_DOUBLE_PRECISION,j, & ! message data and senderId
					   MPI_COMM_WORLD, mpiErrorCode)
		call MPI_BCAST(activeVortStage1(activePanelsIndexStart(j):activePanelsIndexEnd(j)),&
					   activePanelsMessageSize(j),MPI_DOUBLE_PRECISION,j, &
					   MPI_COMM_WORLD, mpiErrorCode)
		call MPI_BCAST(particlesStage1(:,particlesIndexStart(j):particlesIndexEnd(j)),&
					   3*particlesMessageSize(j),MPI_DOUBLE_PRECISION, j, &
					   MPI_COMM_WORLD, mpiErrorCode)
		call MPI_BCAST(passivePanelsStage1(:,passivePanelsIndexStart(j):passivePanelsIndexEnd(j)),&
					   3*passivePanelsMessageSize(j), MPI_DOUBLE_PRECISION, j, &
					   MPI_COMM_WORLD, mpiErrorCode)
	enddo

	!!!!
	!! Kinetic energy calculation (stage 1 only)
	!!!!
	do j=1,activePanels%N
		kineticEnergy(j) = sum(activePanelsStage1(:,j)*activePanelsStage1(:,j))
	enddo
	totalKE = 0.5_kreal*sum(kineticEnergy*area)
	if ( nTracer >= 2) then
		activePanels%tracer(:,2) = kineticEnergy
	endif

	activePanelsStage1 = dt*activePanelsStage1
	activeVortStage1 = dt*activeVortStage1
	particlesStage1 = dt*particlesStage1
	!particlesVortStage1 = dt*particlesVortStage1
	passivePanelsStage1 = dt*passivePanelsStage1


	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! 		RK Stage 2       !
	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Set input arrays for stage 2
	activePanelsInput = activePanels%x + 0.5_kreal*activePanelsStage1
	!activeVortInput = activePanels%relVort + 0.5_kreal*activeVortStage1
	particlesInput = aParticles%x(:,1:aParticles%N) + 0.5_kreal*particlesStage1
	!particlesVortInput = aParticles%relVort(1:aParticles%N) + 0.5_kreal*particlesVortStage1
	passivePanelsINput = passivePanels%x + 0.5_kreal*passivePanelsStage1

	do j=1,aPanels%N_Active
		activeVortInput(j) = activePanels%absVort(j) - 2.0_kreal*Omega*activePanelsInput(3,j) &
					- Juckes_Forcing(activePanelsInput(:,j),t + 0.5_kreal*dt)
	enddo

	! Store velocity output in stage2 arrays
	! Portion work over processes
	call BVEActiveRHS(activePanelsStage2,activeVortStage2,& !output
					  activePanelsInput,activeVortInput, area, & !input
					  activePanelsIndexStart(procRank),activePanelsIndexEnd(procRank)) !indices
	call BVEPassiveRHS(particlesStage2,& !output
					   particlesInput, activePanelsInput,activeVortInput, area, & !input
					   particlesIndexStart(procRank), particlesIndexEnd(procRank)) ! indices
	call BVESmoothVelocity( passivePanelsStage2, & !output
							passivePanelsInput, activePanelsInput, activeVortInput, area, & ! input
							passivePanelsIndexStart(procRank), passivePanelsIndexEnd(procRank)) ! indices

	! broadcast results to all
	do j=0,numProcs-1
		call MPI_BCAST(activePanelsStage2(:,activePanelsIndexStart(j):activePanelsIndexEnd(j)),& ! send buffer
					   3*activePanelsMessageSize(j),MPI_DOUBLE_PRECISION,j, & ! message data and senderId
					   MPI_COMM_WORLD, mpiErrorCode)
		call MPI_BCAST(activeVortStage2(activePanelsIndexStart(j):activePanelsIndexEnd(j)),&
					   activePanelsMessageSize(j),MPI_DOUBLE_PRECISION,j, &
					   MPI_COMM_WORLD, mpiErrorCode)
		call MPI_BCAST(particlesStage2(:,particlesIndexStart(j):particlesIndexEnd(j)),&
					   3*particlesMessageSize(j),MPI_DOUBLE_PRECISION, j, &
					   MPI_COMM_WORLD, mpiErrorCode)
		call MPI_BCAST(passivePanelsStage2(:,passivePanelsIndexStart(j):passivePanelsIndexEnd(j)),&
					   3*passivePanelsMessageSize(j), MPI_DOUBLE_PRECISION, j, &
					   MPI_COMM_WORLD, mpiErrorCode)
	enddo

	activePanelsStage2 = dt*activePanelsStage2
	!activeVortStage2 = dt*activeVortStage2
	particlesStage2 = dt*particlesStage2
	!particlesVortStage2 = dt*particlesVortStage2
	passivePanelsStage2 = dt*passivePanelsStage2

	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! 		RK Stage 3       !
	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Set input arrays for stage 3
	activePanelsInput = activePanels%x + 0.5_kreal*activePanelsStage2
	!activeVortInput = activePanels%relVort + 0.5_kreal*activeVortStage2
	particlesInput = aParticles%x(:,1:aParticles%N) + 0.5_kreal*particlesStage2
	!particlesVortInput = aParticles%relVort(1:aParticles%N) + 0.5_kreal*particlesVortStage2
	passivePanelsINput = passivePanels%x + 0.5_kreal*passivePanelsStage2
	do j=1,aPanels%N_Active
		activeVortInput(j) = activePanels%absVort(j) - 2.0_kreal*Omega*activePanelsInput(3,j) &
					- Juckes_Forcing(activePanelsInput(:,j),t + 0.5_kreal*dt)
	enddo


	! Store velocity output in Stage3 arrays
	! Portion work over processes
	call BVEActiveRHS(activePanelsStage3,activeVortStage3,& !output
					  activePanelsInput,activeVortInput, area, & !input
					  activePanelsIndexStart(procRank),activePanelsIndexEnd(procRank)) !indices
	call BVEPassiveRHS(particlesStage3,& !output
					   particlesInput, activePanelsInput,activeVortInput, area, & !input
					   particlesIndexStart(procRank), particlesIndexEnd(procRank)) ! indices
	call BVESmoothVelocity( passivePanelsStage3, & !output
							passivePanelsInput, activePanelsInput, activeVortInput, area, & ! input
							passivePanelsIndexStart(procRank), passivePanelsIndexEnd(procRank)) ! indices

	! broadcast results to all
	do j=0,numProcs-1
		call MPI_BCAST(activePanelsStage3(:,activePanelsIndexStart(j):activePanelsIndexEnd(j)),& ! send buffer
					   3*activePanelsMessageSize(j),MPI_DOUBLE_PRECISION,j, & ! message data and senderId
					   MPI_COMM_WORLD, mpiErrorCode)
		call MPI_BCAST(activeVortStage3(activePanelsIndexStart(j):activePanelsIndexEnd(j)),&
					   activePanelsMessageSize(j),MPI_DOUBLE_PRECISION,j, &
					   MPI_COMM_WORLD, mpiErrorCode)
		call MPI_BCAST(particlesStage3(:,particlesIndexStart(j):particlesIndexEnd(j)),&
					   3*particlesMessageSize(j),MPI_DOUBLE_PRECISION, j, &
					   MPI_COMM_WORLD, mpiErrorCode)
		call MPI_BCAST(passivePanelsStage3(:,passivePanelsIndexStart(j):passivePanelsIndexEnd(j)),&
					   3*passivePanelsMessageSize(j), MPI_DOUBLE_PRECISION, j, &
					   MPI_COMM_WORLD, mpiErrorCode)
	enddo
	activePanelsStage3 = dt*activePanelsStage3
	particlesStage3 = dt*particlesStage3
	passivePanelsStage3 = dt*passivePanelsStage3
	activeVortStage3 = dt*activeVortStage3

	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! 		RK Stage 4       !
	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Set input arrays for stage 4
	activePanelsInput = activePanels%x + activePanelsStage3
	!activeVortInput = activePanels%relVort + activeVortStage3
	particlesInput = aParticles%x(:,1:aParticles%N) + particlesStage3
	!particlesVortInput = aParticles%relVort(1:aParticles%N) + particlesVortStage3
	passivePanelsINput = passivePanels%x + passivePanelsStage3
	do j=1,aPanels%N_Active
		activeVortInput(j) = activePanels%absVort(j) - 2.0_kreal*Omega*activePanelsInput(3,j) &
					- Juckes_Forcing(activePanelsInput(:,j),t + dt)
	enddo

	! Store velocity output in Stage4 arrays
	! Portion work over processes
	call BVEActiveRHS(activePanelsStage4,activeVortStage4,& !output
					  activePanelsInput,activeVortInput, area, & !input
					  activePanelsIndexStart(procRank),activePanelsIndexEnd(procRank)) !indices
	call BVEPassiveRHS(particlesStage4,& !output
					   particlesInput, activePanelsInput,activeVortInput, area, & !input
					   particlesIndexStart(procRank), particlesIndexEnd(procRank)) ! indices
	call BVESmoothVelocity( passivePanelsStage4, & !output
							passivePanelsInput, activePanelsInput, activeVortInput, area, & ! input
							passivePanelsIndexStart(procRank), passivePanelsIndexEnd(procRank)) ! indices

	! broadcast results to all
	do j=0,numProcs-1
		call MPI_BCAST(activePanelsStage4(:,activePanelsIndexStart(j):activePanelsIndexEnd(j)),& ! send buffer
					   3*activePanelsMessageSize(j),MPI_DOUBLE_PRECISION,j, & ! message data and senderId
					   MPI_COMM_WORLD, mpiErrorCode)
		call MPI_BCAST(activeVortStage4(activePanelsIndexStart(j):activePanelsIndexEnd(j)),&
					   activePanelsMessageSize(j),MPI_DOUBLE_PRECISION,j, &
					   MPI_COMM_WORLD, mpiErrorCode)
		call MPI_BCAST(particlesStage4(:,particlesIndexStart(j):particlesIndexEnd(j)),&
					   3*particlesMessageSize(j),MPI_DOUBLE_PRECISION, j, &
					   MPI_COMM_WORLD, mpiErrorCode)
		call MPI_BCAST(passivePanelsStage4(:,passivePanelsIndexStart(j):passivePanelsIndexEnd(j)),&
					   3*passivePanelsMessageSize(j), MPI_DOUBLE_PRECISION, j, &
					   MPI_COMM_WORLD, mpiErrorCode)
	enddo
	activePanelsStage4 = dt*activePanelsStage4
	particlesStage4 = dt*particlesStage4
	passivePanelsStage4 = dt*passivePanelsStage4
	activeVortStage4 = dt*activeVortStage4

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! 		RK Position update   !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	newParticlesX = aParticles%x(:,1:aParticles%N) + particlesStage1/6.0_kreal +&
		particlesStage2/3.0_kreal + particlesStage3/3.0_kreal + particlesStage4/6.0_kreal
	newActivePanelsX = activePanels%x + activePanelsStage1/6.0_kreal + activePanelsStage2/3.0_kreal +&
		activePanelsStage3/3.0_kreal + activePanelsStage4/6.0_kreal
	newPassivePanelsX = passivePanels%x + passivePanelsStage1/6.0_kreal + passivePanelsStage2/3.0_kreal +&
		passivePanelsStage3/3.0_kreal + passivePanelsStage4/6.0_kreal
	!newActiveVort = activePanels%relVort + activeVortStage1/6.0_kreal + activeVortStage2/3.0_kreal + &
	!	activeVortStage3/3.0_kreal + activeVortStage4/6.0_kreal
	!newParticlesVort = aParticles%relVort(1:aParticles%N) + particlesVortStage1/6.0_kreal + &
	!				   particlesVortStage2/3.0_kreal + particlesVortStage3/3.0_kreal + particlesVortStage4/6.0_kreal

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! 	Copy to data structures   !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	aParticles%x(:,1:aParticles%N) = newParticlesX
	!aParticles%relVort(1:aParticles%N) = newParticlesVort
	activePanels%x = newActivePanelsX
	!activePanels%relVort = newActiveVort
	passivePanels%x = newPassivePanelsX
	do j=1,activePanels%N
		activePanels%relVort(j) = activePanels%absVort(j) - 2.0_kreal*Omega*activePanels%x(3,j) &
			- Juckes_Forcing(activePanels%x(:,j),t+dt)
	enddo
	do j=1,aParticles%N
		aParticles%relVort(j) = aParticles%absVort(j) - 2.0_kreal*Omega*aParticles%x(3,j) &
			- Juckes_Forcing(aParticles%x(:,j),t+dt)
	enddo

	call ScatterPanels(aPanels,activePanels,activeMap,passivePanels,passiveMap)

end subroutine


subroutine AdvectionRK4(aParticles, aPanels, t, dt, procRank, numProcs)
	type(Particles), intent(inout) :: aParticles
	type(Panels), intent(inout) :: aPanels
	real(kreal), intent(in) :: t, dt
	integer(kint), intent(in) :: procRank, numProcs
	integer(kint) :: j

	call LogMessage(log,DEBUG_LOGGING_LEVEL,logKey,'Entering AdvectionRK4')

	if ( rk4isReady ) then
		call ZeroRK4()
	else
		call InitRK4(aParticles,aPanels)
		call LoadBalance(aParticles,aPanels,numProcs)
	endif

	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! 		RK Stage 1       !
	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Set input arrays for stage 1
	particlesInput = aParticles%x(:,1:aParticles%N)
	activePanelsInput = activePanels%x
	passivePanelsInput = passivePanels%x
	! Store velocity output in stage1 arrays
	! Portion work over processes
	do j=activePanelsIndexStart(procRank),activePanelsIndexEnd(procRank)
		activePanelsStage1(:,j) = NonDivergentWind(activePanelsInput(:,j),t)
	enddo
	do j=particlesIndexStart(procRank),particlesIndexEnd(procRank)
		particlesStage1(:,j) = NonDivergentWind(particlesInput(:,j),t)
	enddo
	do j=passivePanelsIndexStart(procRank),passivePanelsIndexEnd(procRank)
		passivePanelsStage1(:,j) = NonDivergentWind(passivePanelsInput(:,j),t)
	enddo

	! broadcast results to all
	do j=0,numProcs-1
		call MPI_BCAST(activePanelsStage1(:,activePanelsIndexStart(j):activePanelsIndexEnd(j)),& ! send buffer
					   3*activePanelsMessageSize(j),MPI_DOUBLE_PRECISION,j, & ! message data and senderId
					   MPI_COMM_WORLD, mpiErrorCode)
		call MPI_BCAST(particlesStage1(:,particlesIndexStart(j):particlesIndexEnd(j)),&
					   3*particlesMessageSize(j),MPI_DOUBLE_PRECISION, j, &
					   MPI_COMM_WORLD, mpiErrorCode)
		call MPI_BCAST(passivePanelsStage1(:,passivePanelsIndexStart(j):passivePanelsIndexEnd(j)),&
					   3*passivePanelsMessageSize(j), MPI_DOUBLE_PRECISION, j, &
					   MPI_COMM_WORLD, mpiErrorCode)
	enddo

	activePanelsStage1 = dt*activePanelsStage1
	particlesStage1 = dt*particlesStage1
	passivePanelsStage1 = dt*passivePanelsStage1

	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! 		RK Stage 2       !
	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Set input arrays for stage 2
	particlesInput = aParticles%x(:,1:aParticles%N) + 0.5_kreal*particlesStage1
	activePanelsInput = activePanels%x + 0.5_kreal*activePanelsStage1
	passivePanelsINput = passivePanels%x + 0.5_kreal*passivePanelsStage1
	! Store velocity output in stage2 arrays
	! Portion work over processes
	do j=activePanelsIndexStart(procRank),activePanelsIndexEnd(procRank)
		activePanelsStage2(:,j) = NonDivergentWind(activePanelsInput(:,j),t + 0.5_kreal*dt)
	enddo
	do j=particlesIndexStart(procRank),particlesIndexEnd(procRank)
		particlesStage2(:,j) = NonDivergentWind(particlesInput(:,j),t + 0.5_kreal*dt)
	enddo
	do j=passivePanelsIndexStart(procRank),passivePanelsIndexEnd(procRank)
		passivePanelsStage2(:,j) = NonDivergentWind(passivePanelsInput(:,j),t + 0.5_kreal*dt)
	enddo
	! broadcast results to all
	do j=0,numProcs-1
		call MPI_BCAST(activePanelsStage2(:,activePanelsIndexStart(j):activePanelsIndexEnd(j)),& ! send buffer
					   3*activePanelsMessageSize(j),MPI_DOUBLE_PRECISION,j, & ! message data and senderId
					   MPI_COMM_WORLD, mpiErrorCode)
		call MPI_BCAST(particlesStage2(:,particlesIndexStart(j):particlesIndexEnd(j)),&
					   3*particlesMessageSize(j),MPI_DOUBLE_PRECISION, j, &
					   MPI_COMM_WORLD, mpiErrorCode)
		call MPI_BCAST(passivePanelsStage2(:,passivePanelsIndexStart(j):passivePanelsIndexEnd(j)),&
					   3*passivePanelsMessageSize(j), MPI_DOUBLE_PRECISION, j, &
					   MPI_COMM_WORLD, mpiErrorCode)
	enddo

	activePanelsStage2 = dt*activePanelsStage2
	particlesStage2 = dt*particlesStage2
	passivePanelsStage2 = dt*passivePanelsStage2

	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! 		RK Stage 3       !
	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Set input arrays for stage 3
	particlesInput = aParticles%x(:,1:aParticles%N) + 0.5_kreal*particlesStage2
	activePanelsInput = activePanels%x + 0.5_kreal*activePanelsStage2
	passivePanelsINput = passivePanels%x + 0.5_kreal*passivePanelsStage2
	! Store velocity output in Stage3 arrays
	! Portion work over processes
	do j=activePanelsIndexStart(procRank),activePanelsIndexEnd(procRank)
		activePanelsStage3(:,j) = NonDivergentWind(activePanelsInput(:,j),t + 0.5_kreal*dt)
	enddo
	do j=particlesIndexStart(procRank),particlesIndexEnd(procRank)
		particlesStage3(:,j) = NonDivergentWind(particlesInput(:,j),t + 0.5_kreal*dt)
	enddo
	do j=passivePanelsIndexStart(procRank),passivePanelsIndexEnd(procRank)
		passivePanelsStage3(:,j) = NonDivergentWind(passivePanelsInput(:,j),t + 0.5_kreal*dt)
	enddo
	! broadcast results to all
	do j=0,numProcs-1
		call MPI_BCAST(activePanelsStage3(:,activePanelsIndexStart(j):activePanelsIndexEnd(j)),& ! send buffer
					   3*activePanelsMessageSize(j),MPI_DOUBLE_PRECISION,j, & ! message data and senderId
					   MPI_COMM_WORLD, mpiErrorCode)
		call MPI_BCAST(particlesStage3(:,particlesIndexStart(j):particlesIndexEnd(j)),&
					   3*particlesMessageSize(j),MPI_DOUBLE_PRECISION, j, &
					   MPI_COMM_WORLD, mpiErrorCode)
		call MPI_BCAST(passivePanelsStage3(:,passivePanelsIndexStart(j):passivePanelsIndexEnd(j)),&
					   3*passivePanelsMessageSize(j), MPI_DOUBLE_PRECISION, j, &
					   MPI_COMM_WORLD, mpiErrorCode)
	enddo

	activePanelsStage3 = dt*activePanelsStage3
	particlesStage3 = dt*particlesStage3
	passivePanelsStage3 = dt*passivePanelsStage3

	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! 		RK Stage 4       !
	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Set input arrays for stage 4
	particlesInput = aParticles%x(:,1:aParticles%N) + particlesStage3
	activePanelsInput = activePanels%x + activePanelsStage3
	passivePanelsINput = passivePanels%x + passivePanelsStage3
	! Store velocity output in Stage4 arrays
	! Portion work over processes
	do j=activePanelsIndexStart(procRank),activePanelsIndexEnd(procRank)
		activePanelsStage4(:,j) = NonDivergentWind(activePanelsInput(:,j),t + dt)
	enddo
	do j=particlesIndexStart(procRank),particlesIndexEnd(procRank)
		particlesStage4(:,j) = NonDivergentWind(particlesInput(:,j),t + dt)
	enddo
	do j=passivePanelsIndexStart(procRank),passivePanelsIndexEnd(procRank)
		passivePanelsStage4(:,j) = NonDivergentWind(passivePanelsInput(:,j),t + dt)
	enddo
	! broadcast results to all
	do j=0,numProcs-1
		call MPI_BCAST(activePanelsStage4(:,activePanelsIndexStart(j):activePanelsIndexEnd(j)),& ! send buffer
					   3*activePanelsMessageSize(j),MPI_DOUBLE_PRECISION,j, & ! message data and senderId
					   MPI_COMM_WORLD, mpiErrorCode)
		call MPI_BCAST(particlesStage4(:,particlesIndexStart(j):particlesIndexEnd(j)),&
					   3*particlesMessageSize(j),MPI_DOUBLE_PRECISION, j, &
					   MPI_COMM_WORLD, mpiErrorCode)
		call MPI_BCAST(passivePanelsStage4(:,passivePanelsIndexStart(j):passivePanelsIndexEnd(j)),&
					   3*passivePanelsMessageSize(j), MPI_DOUBLE_PRECISION, j, &
					   MPI_COMM_WORLD, mpiErrorCode)
	enddo
	activePanelsStage4 = dt*activePanelsStage4
	particlesStage4 = dt*particlesStage4
	passivePanelsStage4 = dt*passivePanelsStage4
	activeVortStage4 = dt*activeVortStage4

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! 		RK Position update   !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	newParticlesX = aParticles%x(:,1:aParticles%N) + particlesStage1/6.0_kreal +&
		particlesStage2/3.0_kreal + particlesStage3/3.0_kreal + particlesStage4/6.0_kreal
	newActivePanelsX = activePanels%x + activePanelsStage1/6.0_kreal + activePanelsStage2/3.0_kreal +&
		activePanelsStage3/3.0_kreal + activePanelsStage4/6.0_kreal
	newPassivePanelsX = passivePanels%x + passivePanelsStage1/6.0_kreal + passivePanelsStage2/3.0_kreal +&
		passivePanelsStage3/3.0_kreal + passivePanelsStage4/6.0_kreal

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! 	Copy to data structures   !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	aParticles%x(:,1:aParticles%N) = newParticlesX
	activePanels%x = newActivePanelsX
	passivePanels%x = newPassivePanelsX

	call ScatterPanels(aPanels,activePanels,activeMap,passivePanels,passiveMap)
end subroutine


subroutine DivergentAdvectionRK4(aParticles, aPanels, t, dt, procRank, numProcs)
	type(Particles), intent(inout) :: aParticles
	type(Panels), intent(inout) :: aPanels
	real(kreal), intent(in) :: t, dt
	integer(kint), intent(in) :: procRank, numProcs
	integer(kint) :: j

	call LogMessage(log,DEBUG_LOGGING_LEVEL,logKey,'Entering AdvectionRK4')

	if ( rk4isReady ) then
		call ZeroRK4()
	else
		call InitRK4(aParticles,aPanels)
		call LoadBalance(aParticles,aPanels,numProcs)
	endif

	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! 		RK Stage 1       !
	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Set input arrays for stage 1
	particlesInput = aParticles%x(:,1:aParticles%N)
	activePanelsInput = activePanels%x
	passivePanelsInput = passivePanels%x
	! Store velocity output in stage1 arrays
	! Portion work over processes
	do j=activePanelsIndexStart(procRank),activePanelsIndexEnd(procRank)
		activePanelsStage1(:,j) = DivergentWind(activePanelsInput(:,j),t)
	enddo
	do j=particlesIndexStart(procRank),particlesIndexEnd(procRank)
		particlesStage1(:,j) = DivergentWind(particlesInput(:,j),t)
	enddo
	do j=passivePanelsIndexStart(procRank),passivePanelsIndexEnd(procRank)
		passivePanelsStage1(:,j) = DivergentWind(passivePanelsInput(:,j),t)
	enddo

	! broadcast results to all
	do j=0,numProcs-1
		call MPI_BCAST(activePanelsStage1(:,activePanelsIndexStart(j):activePanelsIndexEnd(j)),& ! send buffer
					   3*activePanelsMessageSize(j),MPI_DOUBLE_PRECISION,j, & ! message data and senderId
					   MPI_COMM_WORLD, mpiErrorCode)
		call MPI_BCAST(particlesStage1(:,particlesIndexStart(j):particlesIndexEnd(j)),&
					   3*particlesMessageSize(j),MPI_DOUBLE_PRECISION, j, &
					   MPI_COMM_WORLD, mpiErrorCode)
		call MPI_BCAST(passivePanelsStage1(:,passivePanelsIndexStart(j):passivePanelsIndexEnd(j)),&
					   3*passivePanelsMessageSize(j), MPI_DOUBLE_PRECISION, j, &
					   MPI_COMM_WORLD, mpiErrorCode)
	enddo

	activePanelsStage1 = dt*activePanelsStage1
	particlesStage1 = dt*particlesStage1
	passivePanelsStage1 = dt*passivePanelsStage1

	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! 		RK Stage 2       !
	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Set input arrays for stage 2
	particlesInput = aParticles%x(:,1:aParticles%N) + 0.5_kreal*particlesStage1
	activePanelsInput = activePanels%x + 0.5_kreal*activePanelsStage1
	passivePanelsINput = passivePanels%x + 0.5_kreal*passivePanelsStage1
	! Store velocity output in stage2 arrays
	! Portion work over processes
	do j=activePanelsIndexStart(procRank),activePanelsIndexEnd(procRank)
		activePanelsStage2(:,j) = DivergentWind(activePanelsInput(:,j),t + 0.5_kreal*dt)
	enddo
	do j=particlesIndexStart(procRank),particlesIndexEnd(procRank)
		particlesStage2(:,j) = DivergentWind(particlesInput(:,j),t + 0.5_kreal*dt)
	enddo
	do j=passivePanelsIndexStart(procRank),passivePanelsIndexEnd(procRank)
		passivePanelsStage2(:,j) = DivergentWind(passivePanelsInput(:,j),t + 0.5_kreal*dt)
	enddo
	! broadcast results to all
	do j=0,numProcs-1
		call MPI_BCAST(activePanelsStage2(:,activePanelsIndexStart(j):activePanelsIndexEnd(j)),& ! send buffer
					   3*activePanelsMessageSize(j),MPI_DOUBLE_PRECISION,j, & ! message data and senderId
					   MPI_COMM_WORLD, mpiErrorCode)
		call MPI_BCAST(particlesStage2(:,particlesIndexStart(j):particlesIndexEnd(j)),&
					   3*particlesMessageSize(j),MPI_DOUBLE_PRECISION, j, &
					   MPI_COMM_WORLD, mpiErrorCode)
		call MPI_BCAST(passivePanelsStage2(:,passivePanelsIndexStart(j):passivePanelsIndexEnd(j)),&
					   3*passivePanelsMessageSize(j), MPI_DOUBLE_PRECISION, j, &
					   MPI_COMM_WORLD, mpiErrorCode)
	enddo

	activePanelsStage2 = dt*activePanelsStage2
	particlesStage2 = dt*particlesStage2
	passivePanelsStage2 = dt*passivePanelsStage2

	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! 		RK Stage 3       !
	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Set input arrays for stage 3
	particlesInput = aParticles%x(:,1:aParticles%N) + 0.5_kreal*particlesStage2
	activePanelsInput = activePanels%x + 0.5_kreal*activePanelsStage2
	passivePanelsINput = passivePanels%x + 0.5_kreal*passivePanelsStage2
	! Store velocity output in Stage3 arrays
	! Portion work over processes
	do j=activePanelsIndexStart(procRank),activePanelsIndexEnd(procRank)
		activePanelsStage3(:,j) = DivergentWind(activePanelsInput(:,j),t + 0.5_kreal*dt)
	enddo
	do j=particlesIndexStart(procRank),particlesIndexEnd(procRank)
		particlesStage3(:,j) = DivergentWind(particlesInput(:,j),t + 0.5_kreal*dt)
	enddo
	do j=passivePanelsIndexStart(procRank),passivePanelsIndexEnd(procRank)
		passivePanelsStage3(:,j) = DivergentWind(passivePanelsInput(:,j),t + 0.5_kreal*dt)
	enddo
	! broadcast results to all
	do j=0,numProcs-1
		call MPI_BCAST(activePanelsStage3(:,activePanelsIndexStart(j):activePanelsIndexEnd(j)),& ! send buffer
					   3*activePanelsMessageSize(j),MPI_DOUBLE_PRECISION,j, & ! message data and senderId
					   MPI_COMM_WORLD, mpiErrorCode)
		call MPI_BCAST(particlesStage3(:,particlesIndexStart(j):particlesIndexEnd(j)),&
					   3*particlesMessageSize(j),MPI_DOUBLE_PRECISION, j, &
					   MPI_COMM_WORLD, mpiErrorCode)
		call MPI_BCAST(passivePanelsStage3(:,passivePanelsIndexStart(j):passivePanelsIndexEnd(j)),&
					   3*passivePanelsMessageSize(j), MPI_DOUBLE_PRECISION, j, &
					   MPI_COMM_WORLD, mpiErrorCode)
	enddo

	activePanelsStage3 = dt*activePanelsStage3
	particlesStage3 = dt*particlesStage3
	passivePanelsStage3 = dt*passivePanelsStage3

	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! 		RK Stage 4       !
	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Set input arrays for stage 4
	particlesInput = aParticles%x(:,1:aParticles%N) + particlesStage3
	activePanelsInput = activePanels%x + activePanelsStage3
	passivePanelsINput = passivePanels%x + passivePanelsStage3
	! Store velocity output in Stage4 arrays
	! Portion work over processes
	do j=activePanelsIndexStart(procRank),activePanelsIndexEnd(procRank)
		activePanelsStage4(:,j) = DivergentWind(activePanelsInput(:,j),t + dt)
	enddo
	do j=particlesIndexStart(procRank),particlesIndexEnd(procRank)
		particlesStage4(:,j) = DivergentWind(particlesInput(:,j),t + dt)
	enddo
	do j=passivePanelsIndexStart(procRank),passivePanelsIndexEnd(procRank)
		passivePanelsStage4(:,j) = DivergentWind(passivePanelsInput(:,j),t + dt)
	enddo
	! broadcast results to all
	do j=0,numProcs-1
		call MPI_BCAST(activePanelsStage4(:,activePanelsIndexStart(j):activePanelsIndexEnd(j)),& ! send buffer
					   3*activePanelsMessageSize(j),MPI_DOUBLE_PRECISION,j, & ! message data and senderId
					   MPI_COMM_WORLD, mpiErrorCode)
		call MPI_BCAST(particlesStage4(:,particlesIndexStart(j):particlesIndexEnd(j)),&
					   3*particlesMessageSize(j),MPI_DOUBLE_PRECISION, j, &
					   MPI_COMM_WORLD, mpiErrorCode)
		call MPI_BCAST(passivePanelsStage4(:,passivePanelsIndexStart(j):passivePanelsIndexEnd(j)),&
					   3*passivePanelsMessageSize(j), MPI_DOUBLE_PRECISION, j, &
					   MPI_COMM_WORLD, mpiErrorCode)
	enddo
	activePanelsStage4 = dt*activePanelsStage4
	particlesStage4 = dt*particlesStage4
	passivePanelsStage4 = dt*passivePanelsStage4
	activeVortStage4 = dt*activeVortStage4

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! 		RK Position update   !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	newParticlesX = aParticles%x(:,1:aParticles%N) + particlesStage1/6.0_kreal +&
		particlesStage2/3.0_kreal + particlesStage3/3.0_kreal + particlesStage4/6.0_kreal
	newActivePanelsX = activePanels%x + activePanelsStage1/6.0_kreal + activePanelsStage2/3.0_kreal +&
		activePanelsStage3/3.0_kreal + activePanelsStage4/6.0_kreal
	newPassivePanelsX = passivePanels%x + passivePanelsStage1/6.0_kreal + passivePanelsStage2/3.0_kreal +&
		passivePanelsStage3/3.0_kreal + passivePanelsStage4/6.0_kreal

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! 	Copy to data structures   !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	aParticles%x(:,1:aParticles%N) = newParticlesX
	activePanels%x = newActivePanelsX
	passivePanels%x = newPassivePanelsX

	call ScatterPanels(aPanels,activePanels,activeMap,passivePanels,passiveMap)
end subroutine



!----------------
! Particle Velocity Subroutines
!----------------

subroutine BVEActiveRHS(dX, dVort, pointVorts, vort, area, indexStart, indexEnd)
!  Performs the Biot-Savart integral summation in the rotating frame using midpoint rule
!  quadrature for a set of active particles (panel centers).
	real(kreal), intent(out) :: dX(:,:)
	real(kreal), intent(out) :: dVort(:)
	real(kreal), intent(in) :: pointVorts(:,:)
	real(kreal), intent(in) :: vort(:)
	real(kreal), intent(in) :: area(:)
	integer(kint), intent(in) :: indexStart, indexEnd
		integer(kint) :: nn, j, k
	real(kreal) :: xyzJ(3), xyzK(3), denom
	! Check for size mismatch errors
	nn = size(pointVorts,2)
	if ( nn /= size(dX,2) ) then
		print *,"BVEActiveRHS ERROR : size mismatch 1."
		return
	endif
	if ( nn /= size(dVort) ) then
		print *,"BVEActiveRHS ERROR : size mismatch 2."
		return
	endif
	if ( nn /= size(vort) ) then
		print *,"BVEActiveRHS ERROR : size mismatch 3."
		return
	endif
	if ( nn /= size(area) ) then
		print *,"BVEActiveRHS ERROR : size mismatch 4."
		return
	endif
	dX = 0.0_kreal
	dVort = 0.0_kreal
	do j=indexStart,indexEnd
		do k=1,j-1
				denom = 1.0_kreal - sum(pointVorts(:,k)*pointVorts(:,j))
				dX(1,j) = dX(1,j) + (pointVorts(2,j)*pointVorts(3,k) - pointVorts(3,j)*pointVorts(2,k))*vort(k)*area(k)/denom
				dX(2,j) = dX(2,j) + (pointVorts(3,j)*pointVorts(1,k) - pointVorts(1,j)*pointVorts(3,k))*vort(k)*area(k)/denom
				dX(3,j) = dX(3,j) + (pointVorts(1,j)*pointVorts(2,k) - pointVorts(2,j)*pointVorts(1,k))*vort(k)*area(k)/denom
		enddo
		do k=j+1,nn
			denom = 1.0_kreal - sum(pointVorts(:,k)*pointVorts(:,j))
			dX(1,j) = dX(1,j) + (pointVorts(2,j)*pointVorts(3,k) - pointVorts(3,j)*pointVorts(2,k))*vort(k)*area(k)/denom
			dX(2,j) = dX(2,j) + (pointVorts(3,j)*pointVorts(1,k) - pointVorts(1,j)*pointVorts(3,k))*vort(k)*area(k)/denom
			dX(3,j) = dX(3,j) + (pointVorts(1,j)*pointVorts(2,k) - pointVorts(2,j)*pointVorts(1,k))*vort(k)*area(k)/denom
		enddo
	enddo
	dX = dX/(-4.0_kreal*PI)
	dVort = -2.0_kreal*Omega*dX(3,:)
end subroutine


subroutine BVEPassiveRHS(dX,X,pointVorts,vort,area,indexStart,indexEnd)
!  Performs the Biot-Savart integral summation in the rotating frame using midpoint rule
!  quadrature for a set of passive particles (panel vertices).
	real(kreal), intent(out) :: dX(:,:)
	real(kreal), intent(in) :: X(:,:)
	real(kreal), intent(in) :: pointVorts(:,:)
	real(kreal), intent(in) :: vort(:)
	real(kreal), intent(in) :: area(:)
	integer(kint), intent(in) :: indexStart, indexEnd
	integer(kint) :: mm, nn, j, k
	real(kreal) :: xyzJ(3), xyzK(3), denom
	! Check for size mismatch errors
	mm = size(x,2)
	nn = size(pointVorts,2)
	if ( mm /= size(dX,2) ) then
		print *,"BVEPassiveRHS ERROR : size mismatch 1."
		return
	endif
	if ( nn /= size(vort) ) then
		print *,"BVEPassiveRHS ERROR : size mismatch 2."
		return
	endif
	if ( nn /= size(area) ) then
		print *,"BVEPassiveRHS ERROR : size mismatch 3."
		return
	endif
	dX = 0.0_kreal
	do j=indexStart,indexEnd
		do k=1,nn
			denom = 1.0_kreal - sum(pointVorts(:,k)*X(:,j))
			dX(1,j) = dX(1,j) + (X(2,j)*pointVorts(3,k) - X(3,j)*pointVorts(2,k))*vort(k)*area(k)/denom
			dX(2,j) = dX(2,j) + (X(3,j)*pointVorts(1,k) - X(1,j)*pointVorts(3,k))*vort(k)*area(k)/denom
			dX(3,j) = dX(3,j) + (X(1,j)*pointVorts(2,k) - X(2,j)*pointVorts(1,k))*vort(k)*area(k)/denom
		enddo
	enddo
	dX = dX/(-4.0_kreal*PI)
end subroutine


subroutine BVESmoothVelocity(dX,X,pointVorts,vort,Area,indexStart,indexEnd)
!  Performs the regularized Biot-Savart integral summation in the rotating frame using midpoint rule
!  quadrature for a given smoothing parameter.
	real(kreal), intent(out) :: dX(:,:)
	real(kreal), intent(in) :: X(:,:)
	real(kreal), intent(in) :: pointVorts(:,:)
	real(kreal), intent(in) :: vort(:)
	real(kreal), intent(in) :: area(:)
	integer(kint), intent(in) :: indexStart, indexEnd
	integer(kint) :: j,k, nn, mm
	real(kreal) :: xyzJ(3), xyzK(3), denom
	! Check for size mismatch errors
	mm = size(x,2)
	if ( mm /= size(dx,2) ) then
		print *,"BVEInertialVelocity ERROR : size mismatch 1."
		return
	endif
	nn = size(pointVorts,2)
	if ( nn /= size(vort) ) then
		print *,"BVEInertialVelocity ERROR : size mismatch 2."
		return
	endif
	if ( nn /= size(area) ) then
		print *,"BVEInertialVelocity ERROR : size mismatch 3."
		return
	endif
	dX = 0.0_kreal
	do j=indexStart,indexEnd
		do k=1,nn
			denom = 1.0_kreal - sum(pointVorts(:,k)*x(:,j)) + smooth*smooth
			dX(1,j) = dX(1,j) + (X(2,j)*pointVorts(3,k) - X(3,j)*pointVorts(2,k))*vort(k)*area(k)/denom
			dX(2,j) = dX(2,j) + (X(3,j)*pointVorts(1,k) - X(1,j)*pointVorts(3,k))*vort(k)*area(k)/denom
			dX(3,j) = dX(3,j) + (X(1,j)*pointVorts(2,k) - X(2,j)*pointVorts(1,k))*vort(k)*area(k)/denom
		enddo
	enddo
	dX = dX/(-4.0_kreal*PI)
end subroutine


function NonDivergentWind(xyz,t)
	real(kreal), intent(in) :: xyz(3), t
	real(kreal) :: NonDivergentWind(3)
	real(kreal), parameter :: RR = 1.0_kreal, TT = 5.0_kreal
	real(kreal), parameter :: zeroTol = 1e-14
	real(kreal) :: u, v, lat, long, raxis
	lat = latitude(xyz)
	long = longitude(xyz)
	raxis = sqrt(xyz(1)*xyz(1) + xyz(2)*xyz(2))
	if ( raxis > zeroTol) then ! Check divide by zero
		u = 10.0_kreal*RR/TT*sin(long-2.0_kreal*PI*t/TT)*sin(long-2.0_kreal*PI*t/TT)*sin(2.0_kreal*lat)*cos(PI*t/TT) + &
			2.0_kreal*PI*RR/TT*cos(lat)
		v = 10.0_kreal*RR/TT*sin(2.0_kreal*(long-2.0_kreal*PI*t/TT))*cos(lat)*cos(PI*t/TT)

		NonDivergentWind(1) = -u*xyz(2)/raxis - v*xyz(1)*xyz(3)/raxis
		NonDivergentWind(2) =  u*xyz(1)/raxis - v*xyz(2)*xyz(3)/raxis
		NonDivergentWind(3) =  v*raxis
	else
		NonDivergentWind = 0.0_kreal
	endif
end function



function DivergentWind(xyz,t)
	real(kreal), intent(in) :: xyz(3), t
	real(kreal) :: DivergentWind(3)
	real(kreal), parameter :: RR = 1.0_kreal, TT = 5.0_kreal
	real(kreal), parameter :: zeroTol = 1e-14
	real(kreal) :: u, v, lat, long, raxis
	lat = latitude(xyz)
	long = longitude(xyz)
	raxis = sqrt(xyz(1)*xyz(1) + xyz(2)*xyz(2))
	if ( raxis > zeroTol) then ! Check divide by zero
		u = -5.0_kreal*RR/TT*sin((long-2.0_kreal*PI*t/TT)/2.0_kreal)*sin((long-2.0_kreal*PI*t/TT)/2.0_kreal)*sin(2.0_kreal*lat)*&
			cos(lat)*cos(lat)*cos(PI*t/TT) + 2.0_kreal*PI*RR*cos(lat)/TT
		v = 5.0_kreal*RR/(2.0_kreal*TT)*sin(long-2.0_kreal*PI*t/TT)*cos(lat)*cos(lat)*cos(lat)*cos(PI*t/TT)

		DivergentWind(1) = -u*xyz(2)/raxis - v*xyz(1)*xyz(3)/raxis
		DivergentWind(2) =  u*xyz(1)/raxis - v*xyz(2)*xyz(3)/raxis
		DivergentWind(3) =  v*raxis
	else
		DivergentWind = 0.0_kreal
	endif
end function


subroutine StratospherePassiveRHS(dX, dVort, X, pointVorts,vort,area,t,indexStart,indexEnd,smooth)
	real(kreal), intent(out) :: dX(:,:)
	real(kreal), intent(out) :: dVort(:)
	real(kreal), intent(in) :: X(:,:)   ! locations where velocity is needed
	real(kreal), intent(in) :: pointVorts(:,:) ! locations of active particles
	real(kreal), intent(in) :: vort(:)
	real(kreal), intent(in) :: area(:)
	real(kreal), intent(in) :: t
	real(kreal), intent(in), optional :: smooth
	integer(kint), intent(in) :: indexStart, indexEnd
	integer(kint) :: j, k, nn
	real(kreal) :: denom, sm
	real(kreal), parameter :: zeroTol = 1e-12
	real(kreal) :: raxis
	nn = size(pointVorts,2)
	dX = 0.0_kreal
	dVort = 0.0_kreal
	if ( present(smooth) ) then
		sm = smooth
	else
		sm = 0.0_kreal
	endif
	do j=indexStart,indexEnd
		do k=1,nn
			denom = 1.0_kreal - sum(pointVorts(:,k)*X(:,j)) + sm*sm
			dX(1,j) = dX(1,j) + (X(2,j)*pointVorts(3,k) - X(3,j)*pointVorts(2,k))*vort(k)*area(k)/denom
			dX(2,j) = dX(2,j) + (X(3,j)*pointVorts(1,k) - X(1,j)*pointVorts(3,k))*vort(k)*area(k)/denom
			dX(3,j) = dX(3,j) + (X(1,j)*pointVorts(2,k) - X(2,j)*pointVorts(1,k))*vort(k)*area(k)/denom
		enddo
	enddo
	dX = dX/(-4.0_kreal*PI)
	do j=indexStart,indexEnd
		dVort(j) = -2.0_kreal*Omega*dX(3,j) - Juckes_ForcingDerivative(X(:,j),dX(:,j),t)
	enddo
end subroutine


subroutine StratosphereActiveRHS(dX,dVort, pointVorts,vort,area,t,indexStart,indexEnd)
	real(kreal), intent(out) :: dX(:,:)
	real(kreal), intent(out) :: dVort(:)
	real(kreal), intent(in) :: pointVorts(:,:)
	real(kreal), intent(in) :: vort(:)
	real(kreal), intent(in) :: area(:)
	real(kreal), intent(in) :: t
	integer(kint), intent(in) :: indexStart, indexEnd
	integer(kint) :: j, k, nn
	real(kreal) :: denom
	real(kreal), parameter :: zeroTol = 1e-12
	real(kreal) :: raxis
	nn = size(pointVorts,2)
	dX = 0.0_kreal
	dVort = 0.0_kreal
	do j=indexStart,indexEnd
		do k=1,j-1
			denom = 1.0_kreal - sum(pointVorts(:,k)*pointVorts(:,j))
			dX(1,j) = dX(1,j) + (pointVorts(2,j)*pointVorts(3,k) - pointVorts(3,j)*pointVorts(2,k))*vort(k)*area(k)/denom
			dX(2,j) = dX(2,j) + (pointVorts(3,j)*pointVorts(1,k) - pointVorts(1,j)*pointVorts(3,k))*vort(k)*area(k)/denom
			dX(3,j) = dX(3,j) + (pointVorts(1,j)*pointVorts(2,k) - pointVorts(2,j)*pointVorts(1,k))*vort(k)*area(k)/denom
		enddo
		do k=j+1,nn
			denom = 1.0_kreal - sum(pointVorts(:,k)*pointVorts(:,j))
			dX(1,j) = dX(1,j) + (pointVorts(2,j)*pointVorts(3,k) - pointVorts(3,j)*pointVorts(2,k))*vort(k)*area(k)/denom
			dX(2,j) = dX(2,j) + (pointVorts(3,j)*pointVorts(1,k) - pointVorts(1,j)*pointVorts(3,k))*vort(k)*area(k)/denom
			dX(3,j) = dX(3,j) + (pointVorts(1,j)*pointVorts(2,k) - pointVorts(2,j)*pointVorts(1,k))*vort(k)*area(k)/denom
		enddo
	enddo
	dX = dX/(-4.0_kreal*PI)
	do j=indexStart,indexEnd
		dVort(j) = -2.0_kreal*Omega*dX(3,j) - Juckes_ForcingDerivative(pointVorts(:,j),dX(:,j),t)
	enddo
end subroutine


!----------------
! Logger
!----------------

subroutine InitLogger(aLog)
	type(Logger), intent(inout) :: aLog
	integer(kint) :: logUnit
	call New(aLog,logLevel,logUnit)
	logInit = .TRUE.
end subroutine

end module
