module VoronoiVelocitySerialModule

use NumberKindsModule
use SphereGeomModule
use LoggerModule
use VoronoiPanelsModule

implicit none

private
public totalKE  		! kinetic energy integral
public kineticEnergy 	! KE at each panel
public BVERK4, ResetRK4	! Timestepping routines (and memory free routine)
public SetOmega, GetOmega	! background rotation rate
public ResetArea		! Set to true for divergent flow
public InitGaussianVortex, GaussianVortexX	! Vorticity profiles
public InitRH4Wave, HaurwitzStationary4Relvort ! Vorticity profiles
public UseSavedArea

!----------------
! Module variables
!----------------
real(kreal), pointer, save :: kineticEnergy(:)=>null()
real(kreal), pointer, save :: inputX(:)=>null(),&
							  stage1X(:)=>null(),&
							  stage2X(:)=>null(),&
							  stage3X(:)=>null(),&
							  stage4X(:)=>null(),&
							  newX(:)=>null(),&
							  inputY(:)=>null(),&
							  stage1Y(:)=>null(),&
							  stage2Y(:)=>null(),&
							  stage3Y(:)=>null(),&
							  stage4Y(:)=>null(),&
							  newY(:)=>null(),&
							  inputZ(:)=>null(),&
							  stage1Z(:)=>null(),&
							  stage2Z(:)=>null(),&
							  stage3Z(:)=>null(),&
							  stage4Z(:)=>null(),&
							  newZ(:)=>null()
real(kreal), pointer, save :: inputArea(:)=>null(),&
							  stage1Area(:)=>null(),&
							  stage2Area(:)=>null(),&
							  stage3Area(:)=>null(),&
							  stage4Area(:)=>null(),&
							  newArea(:)=>null(),&
							  savedArea(:)=>null()
real(kreal), pointer, save :: inputVort(:)=>null(),&
							  stage1Vort(:)=>null(),&
							  stage2Vort(:)=>null(),&
							  stage3Vort(:)=>null(),&
							  stage4Vort(:)=>null(),&
							  newVort(:)=>null()							  							  

logical(klog), save :: ResetAreas = .False.
logical(klog), save :: saveArea = .False.
logical(klog), save :: areasLoaded = .False.
integer(kint), save :: initCount = 0
logical(klog), save :: rk4IsReady = .False.

real(kreal), save :: gaussConst = 0.0_kreal
real(kreal), save :: totalKE = 0.0_kreal
real(kreal), save :: Omega = 2.0_kreal*PI

!----------------
! Logging
!----------------
logical(klog), save :: logInit = .False.
integer(kint), parameter :: logLevel = DEBUG_LOGGING_LEVEL
character(len=28) :: logKey = 'BVEVelocity : '
type(Logger) :: log

contains

!----------------
! User decisions
!----------------
subroutine SetOmega(newOmega)
	real(kreal), intent(in) :: newOmega
	Omega = newOmega
end subroutine


function GetOmega()
	real(kreal) :: GetOmega
	GetOmega = Omega
end function


subroutine UseSavedArea(onOff)
	logical(klog), intent(in), optional :: onOff
	if ( present(onOff) ) then
		saveArea = onOff
	else
		saveArea = .True.
	endif
end subroutine


subroutine ResetArea(onOff)
	logical(klog), intent(in), optional :: onOff
	if ( present(onOff) ) then
		resetAreas = onOff
	else
		resetAreas = .True.
	endif
end subroutine


!----------------
! Memory Management
!----------------

subroutine InitRK4(aPanels)
	type(VorPanels), intent(in) :: aPanels
	allocate(kineticEnergy(aPanels%N_Max))
	allocate(inputX(aPanels%N_Max))
	allocate(stage1X(aPanels%N_Max))
	allocate(stage2X(aPanels%N_Max))
	allocate(stage3X(aPanels%N_Max))
	allocate(stage4X(aPanels%N_Max))
	allocate(newX(aPanels%N_Max))
	allocate(inputY(aPanels%N_Max))
	allocate(stage1Y(aPanels%N_Max))
	allocate(stage2Y(aPanels%N_Max))
	allocate(stage3Y(aPanels%N_Max))
	allocate(stage4Y(aPanels%N_Max))
	allocate(newY(aPanels%N_Max))
	allocate(inputZ(aPanels%N_Max))
	allocate(stage1Z(aPanels%N_Max))
	allocate(stage2Z(aPanels%N_Max))
	allocate(stage3Z(aPanels%N_Max))
	allocate(stage4Z(aPanels%N_Max))
	allocate(newZ(aPanels%N_Max))
	allocate(inputArea(aPanels%N_Max))
	if ( (.NOT. areasLoaded) .AND. saveArea ) then
		allocate(savedArea(aPanels%N_Max))
		savedArea = aPanels%area
		areasLoaded = .True.
	endif
	if ( resetAreas ) then
		allocate(stage1Area(aPanels%N_Max))
		allocate(stage2Area(aPanels%N_Max))
		allocate(stage3Area(aPanels%N_Max))
		allocate(stage4Area(aPanels%N_Max))
		allocate(newArea(aPanels%N_Max))
	endif
	allocate(inputVort(aPanels%N_Max))
	allocate(stage1Vort(aPanels%N_Max))
	allocate(stage2Vort(aPanels%N_Max))
	allocate(stage3Vort(aPanels%N_Max))
	allocate(stage4Vort(aPanels%N_Max))
	allocate(newVort(aPanels%N_Max))
	
	call ZeroRK4()
	if ( .NOT. logInit) call InitLogger(log)
	rk4isReady = .True.
	call LogMessage(log,DEBUG_LOGGING_LEVEL,logkey,'RK4 initialized.')
end subroutine


subroutine ResetRK4()
	deallocate(kineticEnergy)
	deallocate(inputX)
	deallocate(stage1X)
	deallocate(stage2X)
	deallocate(stage3X)
	deallocate(stage4X)
	deallocate(newX)
	deallocate(inputY)
	deallocate(stage1Y)
	deallocate(stage2Y)
	deallocate(stage3Y)
	deallocate(stage4Y)
	deallocate(newY)
	deallocate(inputZ)
	deallocate(stage1Z)
	deallocate(stage2Z)
	deallocate(stage3Z)
	deallocate(stage4Z)
	deallocate(newZ)
	deallocate(inputArea)
	if ( associated(savedArea) ) deallocate(savedArea)
	if ( resetAreas ) then
		deallocate(stage1Area)
		deallocate(stage2Area)
		deallocate(stage3Area)
		deallocate(stage4Area)
		deallocate(newArea)
	endif
	deallocate(inputVort)
	deallocate(stage1Vort)
	deallocate(stage2Vort)
	deallocate(stage3Vort)
	deallocate(stage4Vort)
	deallocate(newVort)
end subroutine


subroutine ZeroRK4()
	kineticEnergy= 0.0_kreal
	inputX= 0.0_kreal
	stage1X= 0.0_kreal
	stage2X= 0.0_kreal
	stage3X= 0.0_kreal
	stage4X= 0.0_kreal
	newX= 0.0_kreal
	inputY= 0.0_kreal
	stage1Y= 0.0_kreal
	stage2Y= 0.0_kreal
	stage3Y= 0.0_kreal
	stage4Y= 0.0_kreal
	newY= 0.0_kreal
	inputZ= 0.0_kreal
	stage1Z= 0.0_kreal
	stage2Z= 0.0_kreal
	stage3Z= 0.0_kreal
	stage4Z= 0.0_kreal
	newZ= 0.0_kreal
	inputArea= 0.0_kreal
	if  (resetAreas ) then
		stage1Area= 0.0_kreal
		stage2Area= 0.0_kreal
		stage3Area= 0.0_kreal
		stage4Area= 0.0_kreal
		newArea= 0.0_kreal
	endif
	inputVort= 0.0_kreal
	stage1Vort= 0.0_kreal
	stage2Vort= 0.0_kreal
	stage3Vort= 0.0_kreal
	stage4Vort= 0.0_kreal
	newVort= 0.0_kreal
end subroutine


!----------------
! Velocities and time stepping
!----------------

subroutine BVERK4(aPanels, dt)
	type(VorPanels), intent(inout) :: aPanels
	real(kreal), intent(in) :: dt
	integer(kint) :: indexStart, indexEnd, j
	real(kreal) :: norm
	call LogMessage(log,DEBUG_LOGGING_LEVEL,logkey,'Entering BVERK4.')
	if ( .NOT. rk4isReady ) then
		call InitRK4(aPanels)
	else
		call ZeroRK4()	
	endif
	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! 		RK Stage 1       !
	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Set input arrays for stage 1
	indexStart = 1
	indexEnd = aPanels%N
	inputX = aPanels%x
	inputY = aPanels%y
	inputZ = aPanels%z
	if ( saveArea ) then
		inputArea = savedArea
	else
		inputArea = aPanels%area
	endif
	inputVort = aPanels%relVort
	! Store velocity output in stage1 arrays
	call BVEVelocity( stage1X,stage1Y,stage1Z,stage1Vort, & !output data 
					  inputX,inputY,inputZ, inputVort,inputArea,& !input data
					  indexStart, indexEnd) ! indices
	!!!!
	!! Kinetic energy calculation (stage 1 only)
	!!!!
	do j=1,aPanels%N
		kineticEnergy(j) = stage1X(j)*stage1X(j) + stage1Y(j)*stage1Y(j) + stage1Z(j)*stage1Z(j)
	enddo		
	totalKE = 0.5_kreal*sum(kineticEnergy(1:apanels%N)*aPanels%area(1:aPanels%N))
	
	stage1X = dt*stage1X
	stage1Y = dt*stage1Y
	stage1Z = dt*stage1Z
	stage1Vort = dt*stage1Vort
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! 		RK Stage 2       !
	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Set input arrays for stage 2
	inputX = aPanels%x + 0.5_kreal*stage1X
	inputY = aPanels%y + 0.5_kreal*stage1Y
	inputZ = aPanels%z + 0.5_kreal*stage1Z
	inputVort = aPanels%relVort + 0.5_kreal*stage1Vort
	! Store velocity output in stage2 arrays
	call BVEVelocity( 	stage2x, stage2y, stage2z, stage2vort, &
						inputX, inputY, inputZ, inputVort, inputArea, &
						indexStart, indexEnd)
	stage2x = dt*stage2x
	stage2y = dt*stage2y
	stage2z = dt*stage2z
	stage2vort = dt*stage2vort
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! 		RK Stage 3       !
	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Set input arrays for stage 3		
	inputX = aPanels%x + 0.5_kreal*stage2x
	inputY = aPanels%y + 0.5_kreal*stage2y
	inputZ = aPanels%z + 0.5_kreal*stage2z
	inputVort = aPanels%relVort + 0.5_kreal*stage2Vort
	! Store velocity output in stage3 arrays
	call BVEVelocity( 	stage3x, stage3y, stage3z, stage3vort, &
						inputX, inputY, inputZ, inputVort, inputArea, &
						indexStart, indexEnd)
	stage3x = dt*stage3x
	stage3y = dt*stage3y
	stage3z = dt*stage3z
	stage3vort = dt*stage3vort
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! 		RK Stage 4       !
	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Set input arrays for stage 4
	inputX = aPanels%x + stage3x
	inputY = aPanels%y + stage3y
	inputZ = aPanels%z + stage3z
	inputVort = aPanels%relVort + stage3vort
	! Store velocity output in stage4 arrays
	call BVEVelocity( 	stage4x, stage4y, stage4z, stage4vort, &
						inputX, inputY, inputZ, inputVort, inputArea, &
						indexStart, indexEnd)
	
	stage4x = dt*stage4x
	stage4y = dt*stage4y
	stage4z = dt*stage4z
	stage4vort = dt*stage4vort
	
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! 		RK Position update   !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		
	newX(1:aPanels%N) = aPanels%x(1:aPanels%N) + stage1X(1:aPanels%N)/6.0_kreal + &
						stage2X(1:aPanels%N)/3.0_kreal + stage3X(1:aPanels%N)/3.0_kreal + &
						stage4X(1:aPanels%N)/6.0_kreal
	newY(1:aPanels%N) = aPanels%y(1:aPanels%N) + stage1Y(1:aPanels%N)/6.0_kreal + &
						stage2Y(1:aPanels%N)/3.0_kreal + stage3Y(1:aPanels%N)/3.0_kreal +&
						stage4y(1:aPanels%N)/6.0_kreal
	newZ(1:aPanels%N) = aPanels%z(1:aPanels%N) + stage1z(1:aPanels%N)/6.0_kreal + &
						stage2Z(1:aPanels%N)/3.0_kreal + stage3Z(1:aPanels%N)/3.0_kreal +&
						stage4z(1:aPanels%N)/6.0_kreal
	newVort(1:aPanels%N) = aPanels%relVort(1:aPanels%N) + stage1Vort(1:aPanels%N)/6.0_kreal +&
						 stage2Vort(1:aPanels%N)/3.0_kreal + stage3Vort(1:aPanels%N)/3.0_kreal +&
						 stage4Vort(1:aPanels%N)/6.0_kreal

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! 	Copy to data structure   !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!							 
	aPanels%x(1:aPanels%N) = newX(1:aPanels%N)
	aPanels%y(1:aPanels%N) = newY(1:aPanels%N)
	aPanels%z(1:aPanels%N) = newZ(1:aPanels%N)
	aPanels%relVort(1:aPanels%N) = newVort(1:aPanels%N)
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! 	Update connectivity		 !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Renormalize for STRIPACK
	call NormalizeGrid(aPanels)
	call DelaunayTri(aPanels)
	call VoronoiGrid(aPanels)							 
end subroutine


subroutine BVEVelocity(dx,dy,dz,dVort,x,y,z,vort,area,indexStart,indexEnd)
	real(kreal), dimension(:), intent(out) :: dx, dy, dz, dVort
	real(kreal), dimension(:), intent(in) :: x,y,z,vort,area
	integer(kint), intent(in) :: indexStart, indexEnd
	integer(kint) :: nn, j, k
	real(kreal) :: denom
	nn = size(x)
	dX = 0.0_kreal

	do j=indexStart,indexEnd
		do k=1,j-1
			denom = 1.0_kreal - x(j)*x(k) - y(j)*y(k) - z(j)*z(k)
			dx(j) = dx(j) + ( y(j)*z(k) - z(j)*y(k))*vort(k)*area(k)/denom
			dy(j) = dy(j) + ( z(j)*x(k) - x(j)*z(k))*vort(k)*area(k)/denom
			dz(j) = dz(j) + ( x(j)*y(k) - y(j)*x(k))*vort(k)*area(k)/denom
		enddo
		do k=j+1,nn
			denom = 1.0_kreal - x(j)*x(k) - y(j)*y(k) - z(j)*z(k)
			dx(j) = dx(j) + ( y(j)*z(k) - z(j)*y(k))*vort(k)*area(k)/denom
			dy(j) = dy(j) + ( z(j)*x(k) - x(j)*z(k))*vort(k)*area(k)/denom
			dz(j) = dz(j) + ( x(j)*y(k) - y(j)*x(k))*vort(k)*area(k)/denom
		enddo
	enddo

	dx = dx/(-4.0_kreal*PI)
	dy = dy/(-4.0_kreal*PI)
	dz = dz/(-4.0_kreal*PI)
	dVort = -2.0_kreal*Omega*dz	
end subroutine


!----------------
! Vorticity functions
!----------------

function Legendre54(z)
! Calculates the Legendre polynomial P_5^4
	real(kreal), intent(in):: z
	real(kreal) :: Legendre54
	Legendre54 = z*(-1.0_kreal + z*z)*(-1.0_kreal + z*z)
end function


function HaurwitzStationary4Relvort(x,y,z)
	real(kreal), intent(in) :: x, y, z
	real(kreal) :: HaurwitzStationary4RelVort
	real(kreal), parameter :: alphaS = PI/7.0_kreal
	real(kreal), parameter :: cc = -1.0_kreal
	HaurwitzStationary4RelVort = 2.0_kreal*alphaS*z - 30.0_kreal*cc*&
		cos(4.0_kreal*Longitude(x,y,z))*Legendre54(z)
end function


subroutine InitRH4Wave(aPanels)
	type(VorPanels), intent(inout) :: aPanels
	integer(kint) :: j
	do j=1,aPanels%N
		aPanels%relVort(j) = HaurwitzStationary4Relvort(aPanels%x(j),aPanels%y(j),aPanels%z(j))
	enddo
	aPanels%absVort = 2.0_kreal*Omega*aPanels%z + aPanels%relVort
end subroutine


function GaussianVortexX(x,y,z)
	real(kreal), intent(in) :: x, y, z
	real(kreal) :: GaussianVortexX
	real(kreal), parameter :: lat0 = PI/20.0_kreal
	real(kreal), parameter :: long0 = 0.0_kreal
	real(kreal), parameter :: beta = 4.0_kreal
	real(kreal) :: xc, yc, zc
	xc = cos(lat0)*cos(long0)
	yc = cos(lat0)*sin(long0)
	zc = sin(lat0)
	GaussianVortexX = 4.0_kreal*PI*exp(-2.0_kreal*beta*beta*(1.0_kreal - (x*xc + y*yc + z*zc )))
end function


subroutine InitGaussianVortex(aPanels)
	type(VorPanels), intent(inout) :: aPanels
	integer(kint) :: j
	real(kreal), allocatable :: gVort(:)
	allocate(gVort(aPanels%N))
	gVort = 0.0_kreal
	do j=1,aPanels%N
		gVort(j) = GaussianVortexX(aPanels%x(j),aPanels%y(j),aPanels%z(j))
	enddo
	gaussConst = sum(gVort*aPanels%area(1:aPanels%N))/(4.0_kreal*PI)
	aPanels%relVort(1:aPanels%N) = gVort - gaussConst
	aPanels%absVort(1:aPanels%N) = aPanels%relVort(1:aPanels%N) + 2.0_kreal*Omega*aPanels%z(1:aPanels%N)
	deallocate(gVort)
end subroutine

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
