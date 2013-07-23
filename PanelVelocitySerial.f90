module PanelVelocitySerialModule

use NumberKindsModule
use SphereGeomModule
use LoggerModule
use ParticlesEdgesPanelsModule
use TracerAndVorticityDistributionModule

implicit none

private
public BVERK4, ResetRK4
public omega, smooth
public InitGaussianVortex, GaussianVortexX, gaussConst
public InitRH4, HaurwitzStationary4RelVort
public totalKE, kineticEnergy
public SetOmega, SetSmoothingParameter, ResetArea

!----------------
! Module variables
!----------------
real(kreal), allocatable, save :: kineticEnergy(:)
real(kreal), save :: totalKE
real(kreal), save :: omega = 2.0_kreal*PI
real(kreal), save :: smooth = 0.05_kreal
real(kreal), save :: gaussConst = 0.0_kreal

type(Panels), pointer, save :: activePanels=>null(), passivePanels=>null()
integer(kint), allocatable :: activeMap(:), passiveMap(:)

real(kreal), allocatable, save :: particlesInput(:,:),&
								  particlesStage1(:,:),&
								  particlesStage2(:,:),&
								  particlesStage3(:,:),&
								  particlesStage4(:,:),&
								  newParticlesX(:,:),&
								  activePanelsInput(:,:),&
								  activePanelsStage1(:,:),&
								  activePanelsStage2(:,:),&
								  activePanelsStage3(:,:),&
								  activePanelsStage4(:,:),&
								  newActivePanelsX(:,:),&
								  passivePanelsInput(:,:),&
								  passivePanelsStage1(:,:),&
								  passivePanelsStage2(:,:),&
								  passivePanelsStage3(:,:),&
								  passivePanelsStage4(:,:),&
								  newPassivePanelsX(:,:),&
								  activeVortInput(:),&
								  activeVortStage1(:),&
								  activeVortStage2(:),&
								  activeVortStage3(:),&
								  activeVortStage4(:),&
								  newActiveVort(:),&
								  area(:),&
								  areaStage1(:),&
								  areaStage2(:),&
								  areaStage3(:),&
								  areaStage4(:),&
								  newArea(:)



logical(klog), save :: rk4isReady = .False.
logical(klog), save :: resetAreas = .False.

logical(klog) :: logInit = .False.
integer(kint), parameter :: logLevel = DEBUG_LOGGING_LEVEL
character(len=28) :: logKey = "PanelVelocity : "
type(Logger) :: Log

contains
!----------------
! User decisions
!----------------
subroutine SetSmoothingParameter(newSmooth)
	real(kreal), intent(in) :: newSmooth
	smooth = newSmooth
end subroutine


subroutine SetOmega(newOmega)
	real(kreal), intent(in) :: newOmega
	Omega = newOmega
end subroutine


subroutine ResetArea(onOff)
	logical(klog), intent(in), optional :: onOff
	if ( present(onOff)) then
		resetAreas = onOff
	else
		resetAreas = .True.
	endif
end subroutine


!----------------
! Memory handling 
!----------------
subroutine InitRK4(aParticles,aPanels)
	type(Particles), intent(in) :: aParticles
	type(Panels), intent(in) :: aPanels
	integer(kint) :: nPassive, nActive, nParticles, nTracer, panelKind, problemKind
	
	if ( .NOT. logInit) call InitLogger(log)	
	
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
	activePanels%N = nActive
	activePanels%N_Active = nActive
	passivePanels%N = nPassive
	call GatherPanels(aPanels,activePanels,activeMap,passivePanels,passiveMap)	
	
	allocate(particlesInput(3,nParticles))
	allocate(particlesStage1(3,nParticles))
	allocate(particlesStage2(3,nParticles))
	allocate(particlesStage3(3,nParticles))
	allocate(particlesStage4(3,nParticles))
	allocate(newParticlesX(3,nParticles))
	
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

subroutine BVERK4(aParticles,aPanels,dt)
	type(Particles), intent(inout) :: aParticles
	type(Panels), intent(inout) :: aPanels
	real(kreal), intent(in) :: dt
	integer(kint) :: particlesIndexStart, particlesIndexEnd
	integer(kint) :: activeIndexStart, activeIndexEnd
	integer(kint) :: passiveIndexStart, passiveIndexEnd
	integer(kint) :: j, nTracer
	
	if ( rk4isReady ) then
		call ZeroRK4()
	else
		call InitRK4(aParticles,aPanels)
	endif
	
	nTracer = GetNTracer(apanels)
	
	particlesIndexStart = 1
	activeIndexStart = 1
	passiveIndexStart = 1
	particlesIndexEnd = aParticles%N
	activeIndexEnd = aPanels%N_Active
	passiveIndexEnd = aPanels%N - aPanels%N_Active
	
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
	call BVEActiveRHS(activePanelsStage1,activeVortStage1,& !output
					  activePanelsInput,activeVortInput, area, & !input
					  activeIndexStart,activeIndexEnd) !indices
	call BVEPassiveRHS(particlesStage1,& !output
					   particlesInput, activePanelsInput,activeVortInput, area, & !input
					   particlesIndexStart, particlesIndexEnd) ! indices 
	call BVESmoothVelocity( passivePanelsStage1, & !output
							passivePanelsInput, activePanelsInput, activeVortInput, area, & ! input
							passiveIndexStart, passiveINdexEnd) ! indices		
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
	! Set input arrays for stage 2
	particlesInput = aParticles%x(:,1:aParticles%N) + 0.5_kreal*particlesStage1
	activePanelsInput = activePanels%x + 0.5_kreal*activePanelsStage1
	passivePanelsINput = passivePanels%x + 0.5_kreal*passivePanelsStage1
	activeVortInput = activePanels%relVort + 0.5_kreal*activeVortStage1
	! Store velocity output in stage2 arrays
	call BVEActiveRHS(activePanelsStage2,activeVortStage2,& !output
					  activePanelsInput,activeVortInput, area, & !input
					  activeIndexStart,activeIndexEnd) !indices
	call BVEPassiveRHS(particlesStage2,& !output
					   particlesInput, activePanelsInput,activeVortInput, area, & !input
					   particlesIndexStart, particlesIndexEnd) ! indices 
	call BVESmoothVelocity( passivePanelsStage2, & !output
							passivePanelsInput, activePanelsInput, activeVortInput, area, & ! input
							passiveIndexStart, passiveINdexEnd) ! indices
	activePanelsStage2 = dt*activePanelsStage2
	particlesStage2 = dt*particlesStage2
	passivePanelsStage2 = dt*passivePanelsStage2
	activeVortStage2 = dt*activeVortStage2
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! 		RK Stage 3       !
	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Set input arrays for stage 3
	particlesInput = aParticles%x(:,1:aParticles%N) + 0.5_kreal*particlesStage2
	activePanelsInput = activePanels%x + 0.5_kreal*activePanelsStage2
	passivePanelsINput = passivePanels%x + 0.5_kreal*passivePanelsStage2
	activeVortInput = activePanels%relVort + 0.5_kreal*activeVortStage2
	! Store velocity output in Stage3 arrays
	call BVEActiveRHS(activePanelsStage3,activeVortStage3,& !output
					  activePanelsInput,activeVortInput, area, & !input
					  activeIndexStart,activeIndexEnd) !indices
	call BVEPassiveRHS(particlesStage3,& !output
					   particlesInput, activePanelsInput,activeVortInput, area, & !input
					   particlesIndexStart, particlesIndexEnd) ! indices 
	call BVESmoothVelocity( passivePanelsStage3, & !output
							passivePanelsInput, activePanelsInput, activeVortInput, area, & ! input
							passiveIndexStart, passiveINdexEnd) ! indices
	activePanelsStage3 = dt*activePanelsStage3
	particlesStage3 = dt*particlesStage3
	passivePanelsStage3 = dt*passivePanelsStage3
	activeVortStage3 = dt*activeVortStage3
	
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! 		RK Stage 4       !
	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Set input arrays for stage 4
	particlesInput = aParticles%x(:,1:aParticles%N) + particlesStage3
	activePanelsInput = activePanels%x + activePanelsStage3
	passivePanelsINput = passivePanels%x + passivePanelsStage3
	activeVortInput = activePanels%relVort + activeVortStage3
	! Store velocity output in Stage4 arrays
	call BVEActiveRHS(activePanelsStage4,activeVortStage4,& !output
					  activePanelsInput,activeVortInput, area, & !input
					  activeIndexStart,activeIndexEnd) !indices
	call BVEPassiveRHS(particlesStage4,& !output
					   particlesInput, activePanelsInput,activeVortInput, area, & !input
					   particlesIndexStart, particlesIndexEnd) ! indices 
	call BVESmoothVelocity( passivePanelsStage4, & !output
							passivePanelsInput, activePanelsInput, activeVortInput, area, & ! input
							passiveIndexStart, passiveINdexEnd) ! indices
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
	newActiveVort = activePanels%relVort + activeVortStage1/6.0_kreal + activeVortStage2/3.0_kreal + &
		activeVortStage3/3.0_kreal + activeVortStage4/6.0_kreal
		
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! 	Copy to data structures   !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		
	
	aParticles%x(:,1:aParticles%N) = newParticlesX
	activePanels%x = newActivePanelsX
	activePanels%relVort = newActiveVort
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


!----------------
! Vorticity functions
!----------------

subroutine InitRH4(aParticles,aPanels)
! Initializes the vorticity profiles (absolute and relative) of a stationary Rossby-Haurwitz wave with zonal 
! wavenumber 4 on a grid.
	type(Particles), intent(inout) :: aParticles
	type(Panels), intent(inout) :: aPanels
	integer(kint) :: j
	do j=1,aParticles%N
		aParticles%relVort(j) = HaurwitzStationary4RelVort(aParticles%x(:,j))
		aParticles%absVort(j) = aParticles%relVort(j) + 2.0_kreal*Omega*aParticles%x(3,j) 
	enddo
	do j=1,aPanels%N
		aPanels%relVort(j) = HaurwitzStationary4RelVort(aPanels%x(:,j))
		aPanels%absVort(j) = aPanels%relVort(j) + 2.0_kreal*Omega*aPanels%x(3,j)
	enddo
end subroutine


function Legendre54(z)
! Calculates the Legendre polynomial P_5^4
	real(kreal), intent(in):: z
	real(kreal) :: Legendre54
	Legendre54 = z*(-1.0_kreal + z*z)*(-1.0_kreal + z*z)
end function


function HaurwitzStationary4RelVort(xyz)
! Calculates the relative vorticity as a function of position for a Rossby-Haurwitz wave 
! with zonal wavenumber 4.
	real(kreal), intent(in) :: xyz(3)
	real(kreal) :: HaurwitzStationary4RelVort
	! prevailing zonal wind angular velocity
	real(kreal), parameter :: alphaS = PI/7.0_kreal
	real(kreal), parameter :: cc = -1.0_kreal
	HaurwitzStationary4RelVort = 2.0_kreal*alphaS*xyz(3) - 30.0_kreal*cc*cos(4.0_kreal*longitude(xyz))*Legendre54(xyz(3))
end function


subroutine InitGaussianVortex(aParticles, aPanels)
! Initializes the vorticity profiles (absolute and relative) of a Gaussian vortex on the grid.
! Calculates the constant required to keep total vorticity = 0 over the sphere ( a necessary
! condition for the Poisson problem to have a solution) and stores that value for use
! later in the calculation (i.e., at each remeshing).
	type(Particles), intent(inout) :: aParticles
	type(Panels), intent(inout) :: aPanels
	real(kreal), allocatable :: gvort(:)
	integer(kint) :: j
	allocate(gvort(aPanels%N))
	gvort = 0.0_kreal
	do j=1,aPanels%N
		if ( .NOT. aPanels%hasChildren(j) ) then
			gvort(j) = GaussianVortexX(aPanels%x(:,j))
		endif 
	enddo
	gaussConst = sum(gvort*aPanels%area(1:aPanels%N))/(4.0_kreal*PI)
	do j=1,aParticles%N
		aParticles%absVort(j) = GaussianVortexX(aParticles%x0(:,j)) - gaussConst + &
				 2.0_kreal*Omega*aParticles%x0(3,j)
		aParticles%relVort(j) = aParticles%absVort(j) - 2.0_kreal*Omega*aParticles%x(3,j)
	enddo
	do j=1,aPanels%N
		if (.NOT. aPanels%hasChildren(j) ) then
			aPanels%absVort(j) = GaussianVortexX(aPanels%x0(:,j)) - gaussConst  + &
				2.0_kreal*Omega*aPanels%x0(3,j)
			aPanels%relVort(j) = aPanels%absVort(j) - 2.0_kreal*Omega*aPanels%x(3,j)
		else
			aPanels%relVort(j) = 0.0_kreal	
			aPanels%absVort(j) = 0.0_kreal
		endif
	enddo
	deallocate(gvort)
end subroutine


function GaussianVortexX(xyz)
! Outputs the height of a Gaussian function used to calculate the relative vorticity 
! asssociated with a Gaussian vortex.  The appropriate constant must be subtracted from
! the output value of this function to give a valid vorticity profile.  
	real(kreal) :: GaussianVortexX
	real(kreal), intent(in) :: xyz(3)
	real(kreal), parameter :: 	lat0 = PI/20.0_kreal,&
								long0 = 0.0_kreal,&
								beta = 4.0_kreal
	real(kreal) :: xCenter(3)
	xCenter(1) = cos(lat0)*cos(long0)
	xCenter(2) = cos(lat0)*sin(long0)
	xCenter(3) = sin(lat0)
	GaussianVortexX = 4.0_kreal*PI*exp(-2.0_kreal*beta*beta*(1.0_kreal - sum(xyz*xCenter)))
end function



!----------------
! Logger
!----------------

subroutine InitLogger(aLog)
	type(Logger), intent(inout) :: aLog
	integer(kint) :: logUnit
	logUnit = 6
	call New(aLog,logLevel,logUnit)
	logInit = .TRUE.
end subroutine

end module
