module VoronoiVelocitySerialModule

use NumberKindsModule
use SphereGeomModule
use LoggerModule
use VoronoiPanelsModule

implicit none

private
public totalKE, kineticEnergy
public BVERK4, ResetRK4
public SetOmega, GetOmega
public InitGaussianVortex, GaussianVortexX
public InitRH4Wave, HaurwitzStationary4Relvort
!----------------
! Types and Module variables
!----------------
type(VorPanels), pointer, save :: rkStage1=>null(),&
								  rkStage2=>null(),&
								  rkStage3=>null(),&
								  rkStage4=>null(),&
								  rkInput=>null(),&
								  rkOutput=>null()
real(kreal), pointer, save :: kineticEnergy(:)=>null()								  

real(kreal), save :: Omega = 2.0_kreal*PI
real(kreal), save :: totalKE = 0.0_kreal
logical(klog), save :: rk4isReady = .False.
real(kreal), save :: GaussConst = 0.0_kreal
type(Logger) :: log
integer(kint), parameter :: logLevel = DEBUG_LOGGING_LEVEL
logical(klog), save :: logInit=.False.
character(len=28) :: logKey='BVEVelocity : '

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

!----------------
! Memory Management
!----------------
subroutine InitRK4(aPanels)
	type(VorPanels), intent(in) :: aPanels
	call New(rkStage1,aPanels)
	call New(rkStage2,aPanels)
	call New(rkStage3,aPanels)
	call New(rkStage4,aPanels)
	call New(rkInput,aPanels)
	call New(rkOutput,aPanels)
	allocate(kineticEnergy(aPanels%N_Max))
	kineticEnergy = 0.0_kreal
	call ZeroRK4()
	call InitLogger(log)
	rk4isReady = .True.
	call LogMessage(log,DEBUG_LOGGING_LEVEL,logkey,'RK4 initialized.')
end subroutine


subroutine ResetRK4()
	call Delete(rkStage1)
	call Delete(rkStage2)
	call Delete(rkStage3)
	call Delete(rkStage4)
	call Delete(rkInput)
	call Delete(rkOutput)
	deallocate(kineticEnergy)
	call LogMessage(log,DEBUG_LOGGING_LEVEL,logkey,'RK4 memory released.')
	call Delete(log)
	rk4isReady = .False.
end subroutine


subroutine ZeroRK4()
	call LogMessage(log,DEBUG_LOGGING_LEVEL,logkey,'Entering ZeroRK4.')
	rkOutput%x = 0.0_kreal
	rkOutput%y = 0.0_kreal
	rkOutput%z = 0.0_kreal
	rkOutput%relVort = 0.0_kreal
	rkInput%x = 0.0_kreal
	rkInput%y = 0.0_kreal
	rkInput%z = 0.0_kreal
	rkInput%relVort = 0.0_kreal
	rkStage1%x = 0.0_kreal
	rkStage1%y = 0.0_kreal
	rkStage1%z = 0.0_kreal
	rkStage1%relVort = 0.0_kreal
	rkStage2%x = 0.0_kreal
	rkStage2%y = 0.0_kreal
	rkStage2%z = 0.0_kreal
	rkStage2%relVort = 0.0_kreal
	rkStage3%x = 0.0_kreal
	rkStage3%y = 0.0_kreal
	rkStage3%z = 0.0_kreal
	rkStage3%relVort = 0.0_kreal
	rkStage4%x = 0.0_kreal
	rkStage4%y = 0.0_kreal
	rkStage4%z = 0.0_kreal
	rkStage4%relVort = 0.0_kreal
end subroutine


!----------------
! Timestepping Subroutines
!----------------

subroutine BVERK4(aPanels,dt)
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
	rkInput = aPanels
	! Store output in rkStage1
	call BVEVelocity(rkStage1%x,rkStage1%y,rkStage1%z,rkStage1%relVort,&
			rkInput%x,rkInput%y,rkInput%z,rkInput%relVort,rkInput%area,indexStart,indexEnd)
	! STAGE 1 ONLY : Calculate kinetic energy
	do j=1,aPanels%N
		kineticEnergy(j) = rkStage1%x(j)*rkStage1%x(j) + rkStage1%y(j)*rkStage1%y(j) + &
						   rkStage1%z(j)*rkStage1%z(j)
	enddo
	totalKE = 0.5_kreal*sum(kineticEnergy(1:aPanels%N)*aPanels%area(1:aPanels%N))
			
	rkStage1%x = dt*rkStage1%x		
	rkStage1%y = dt*rkStage1%y		
	rkStage1%z = dt*rkStage1%z		
	rkStage1%relVort = dt*rkStage1%relVort
	!call LogMessage(log,DEBUG_LOGGING_LEVEL,logkey,'Stage 1 complete.')
	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! 		RK Stage 2       !
	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Set input for stage 2
	rkInput%x = aPanels%x + 0.5_kreal*rkStage1%x
	rkInput%y = aPanels%y + 0.5_kreal*rkStage1%y
	rkInput%z = aPanels%z + 0.5_kreal*rkStage1%z
	rkInput%relVort = rkInput%relVort + 0.5_kreal*rkStage1%relVort
	! Store output in rkStage2
	call BVEVelocity(rkStage2%x,rkStage2%y,rkStage2%z,rkStage2%relVort,&
			rkInput%x,rkInput%y,rkInput%z,rkInput%relVort,rkInput%area,indexStart,indexEnd)
	rkStage2%x = dt*rkStage2%x
	rkStage2%y = dt*rkStage2%y
	rkStage2%z = dt*rkStage2%z
	rkStage2%relVort = dt*rkStage2%relVort
	!call LogMessage(log,DEBUG_LOGGING_LEVEL,logkey,'Stage 2 complete.')
	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! 		RK Stage 3       !
	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Set input for stage 3
	rkInput%x = aPanels%x + 0.5_kreal*rkStage2%x
	rkInput%y = apanels%y + 0.5_kreal*rkStage2%y
	rkInput%z = aPanels%z + 0.5_kreal*rkStage2%z
	rkInput%relVort = apanels%relVort + 0.5_kreal*rkStage2%relVort
	! Store output in rkStage3
	call BVEVelocity(rkStage3%x,rkStage3%y,rkStage3%z,rkStage3%relVort,&
			rkInput%x,rkInput%y,rkInput%z,rkInput%relVort,rkInput%area,indexStart,indexEnd)
	rkStage3%x = dt*rkStage3%x
	rkStage3%y = dt*rkStage3%y
	rkStage3%z = dt*rkStage3%z
	rkStage3%relVort = dt*rkStage3%relVort
	!call LogMessage(log,DEBUG_LOGGING_LEVEL,logkey,'Stage 3 complete.')
	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! 		RK Stage 4       !
	!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Set input for stage 4
	rkInput%x = aPanels%x + rkStage3%x
	rkInput%y = aPanels%y + rkStage3%y
	rkInput%z = aPanels%z + rkStage3%z
	rkInput%relVort = aPanels%relVort + rkStage3%relVort
	! Store output in rkStage4
	call BVEVelocity(rkStage4%x,rkStage4%y,rkStage4%z,rkStage4%relVort,&
			rkInput%x,rkInput%y,rkInput%z,rkInput%relVort,rkInput%area,indexStart,indexEnd)
	rkStage4%x = dt*rkStage4%x
	rkStage4%y = dt*rkStage4%y
	rkStage4%z = dt*rkStage4%z
	rkStage4%relVort = dt*rkStage4%relVort
	!call LogMessage(log,DEBUG_LOGGING_LEVEL,logkey,'Stage 4 complete.')
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! RK Position update       !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	rkOutput%x = apanels%x + rkStage1%x/6.0_kreal + rkStage2%x/3.0_kreal +&
		 rkStage3%x/3.0_kreal + rkStage4%x/6.0_kreal
 	rkOutput%y = apanels%y + rkStage1%y/6.0_kreal + rkStage2%y/3.0_kreal +&
		 rkStage3%y/3.0_kreal + rkStage4%y/6.0_kreal
	rkOutput%z = apanels%z + rkStage1%z/6.0_kreal + rkStage2%z/3.0_kreal +&
		 rkStage3%z/3.0_kreal + rkStage4%z/6.0_kreal	 
	rkOutput%relVort = aPanels%relVort + rkStage1%relVort/6.0_kreal + rkStage2%relVort/3.0_kreal +&
		 rkStage3%relVort/3.0_kreal + rkStage4%relVort/6.0_kreal
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! 	Reset Delaunay / Voronoi connections !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!do j=1,rkOutput%N
!		norm = sqrt(rkOutput%x(j)*rkOutput%x(j) + rkOutput%y(j)*rkOutput%y(j) + rkOutput%z(j)*rkOutput%z(j))
!		rkOutput%x(j) = rkOutput%x(j)/norm
!		rkOutput%y(j) = rkOutput%y(j)/norm
!		rkOutput%z(j) = rkOutput%z(j)/norm
!	enddo
	call DelaunayTri(rkOutput)
	call VoronoiGrid(rkoutput)
	aPanels = rkOutput
	
end subroutine



subroutine BVEVelocity(dx,dy,dz,dVort,x,y,z,vort,area,indexStart,indexEnd)
	real(kreal), dimension(:), intent(out) :: dx, dy, dz, dVort
	real(kreal), dimension(:), intent(in) :: x,y,z,vort,area
	integer(kint), intent(in) :: indexStart, indexEnd
	integer(kint) :: nn, j, k
	real(kreal) :: denom
	!call LogMessage(log,DEBUG_LOGGING_LEVEL,logkey,'Entering BVEVelocity.')
	nn = size(x)
	dX = 0.0_kreal
	!dVort = 0.0_kreal
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


subroutine InitLogger(aLog)
	type(Logger), intent(inout) :: aLog
	integer(kint) :: logUnit
	logUnit = 6
	call New(aLog,logLevel,logUnit)
	logInit = .TRUE.
end subroutine

end module
