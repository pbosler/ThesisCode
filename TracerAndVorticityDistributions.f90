module TracerAndVorticityDistributionModule

use NumberKindsModule
use SphereGeomModule
use LoggerModule
use ParticlesEdgesPanelsModule
use VoronoiPanelsModule

implicit none

private
public gaussConst, omega
public SetOmega
public InitGaussianVortex, GaussianVortexX
public InitRH4Wave, HaurwitzStationary4RelVort
public InitManyVorts, ManyVortsX
public Juckes_A, Juckes_APrime
public Juckes_B, Juckes_BPrime
public Juckes_Forcing, Juckes_ForcingDerivative
public InitRH2Wave, RH2WaveX
public InitPolarVortex, PolarVortexX
public InitZonalMean, ZonalMeanX
public InitJet, JetRelVortX
public InitGaussianHills, GaussianHillsX
public InitSlottedCylinders, SlottedCylindersX
public InitCosineBells, CosineBellsX
public InitWilliamsonCosineBell, WilliamsonCosineBellX
public InitStratosphere, StratosphereRelVortX
public InitBlockM, BlockMX
public InitTwoVorts, TwoVortsX
public InitTripole, TripoleX
public InitSolidBody, SolidBodyX

real(kreal), save :: gaussConst = 0.0_kreal, &
					 omega = 2.0_kreal*PI

type(Logger) :: log
logical(klog), save :: logInit =.False.
character(len=28), save :: logKey = 'TracerVorticity : '
integer(kint), parameter :: logLevel = DEBUG_LOGGING_LEVEL, logUnit=6

interface InitGaussianVortex
	module procedure InitGaussianVortexVoronoi
	module procedure InitGaussianVortexPanels
end interface

interface GaussianVortexX
	module procedure GaussianVortexVector
	module procedure GaussianVortexComponents
end interface

interface InitRH4Wave
	module procedure InitRH4WaveVoronoi
	module procedure InitRH4WavePanels
end interface

interface InitRH2Wave
	module procedure InitRH2WaveVoronoi
	module procedure InitRH2WavePanels
end interface

interface HaurwitzStationary4RelVort
	module procedure HaurwitzStationary4RelVortVector
	module procedure HaurwitzStationary4RelVortComponents
end interface

interface RH2WaveX
	module procedure RHWave2RelVortVector
	module procedure RHWave2RelVortComponents
end interface

interface InitManyVorts
	module procedure InitManyVortsPanels
	module procedure InitManyVortsVoronoi
end interface

interface ManyVortsX
	module procedure ManyVortsVector
	module procedure ManyVortsComponents
end interface

interface InitPolarVortex
	module procedure InitPolarVortexPanels
	module procedure InitPolarVortexVoronoi
end interface

interface PolarVortexX
	module procedure PolarVortexVector
	module procedure PolarVortexComponents
end interface

interface InitZonalMean
	module procedure InitZonalMeanPanels
	module procedure InitZonalMeanVoronoi
end interface

interface ZonalMeanX
	module procedure ZonalMeanRelVortVector
	module procedure ZonalMeanRelVortComponents
end interface

interface InitJet
	module procedure InitJetPanels
	module procedure InitJetVoronoi
end interface

interface JetRelVortX
	module procedure JetRelVortVector
	module procedure JetRelVortComponents
end interface

interface Juckes_B
	module procedure Juckes_BVector
	module procedure Juckes_Bscalar
end interface

interface Juckes_BPrime
	module procedure Juckes_BPrimeVector
	module procedure Juckes_BPrimeScalar
end interface

interface BlockMX
	module procedure BlockMVector
end interface

interface TwoVortsX
	module procedure TwoGaussVortsVector
end interface

interface TripoleX
	module procedure TripoleVortVector
end interface

contains

subroutine SetOmega(newOmega)
	real(kreal), intent(in) :: newOmega
	Omega = newOmega
end subroutine


subroutine InitRH4WavePanels(aParticles,aPanels)
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


subroutine InitRH2WavePanels(aParticles,aPanels)
	type(Particles), intent(inout) :: aParticles
	type(Panels), intent(inout) :: aPanels
	integer(kint) :: j
	do j=1,aParticles%N
		aParticles%relVort(j) = RHWave2RelVortVector(aParticles%x0(:,j))
		aParticles%absVort(j) = aParticles%relVort(j) + 2.0_kreal*Omega*aParticles%x0(3,j)
	enddo
	do j=1,aPanels%N
		if ( .NOT. aPanels%hasChildren(j)) then
			aPanels%relvort(j) = RHWave2RelVortVector(aPanels%x0(:,j))
			aPanels%absVort(j) = aPanels%relVort(j) + 2.0_kreal*Omega*aPanels%x0(3,j)
		endif
	enddo
end subroutine


subroutine InitRH4WaveVoronoi(aPanels)
	type(VorPanels), intent(inout) :: aPanels
	integer(kint) :: j
	do j=1,aPanels%N
		aPanels%relVort(j) = HaurwitzStationary4Relvort(aPanels%x(j),aPanels%y(j),aPanels%z(j))
	enddo
	aPanels%absVort = 2.0_kreal*Omega*aPanels%z + aPanels%relVort
end subroutine


subroutine InitRH2WaveVoronoi(aPanels)
	type(VorPanels), intent(inout) :: aPanels
	integer(kint) :: J
	do j=1,aPanels%N
		aPanels%relVort(j) = RHWave2RelVortComponents(aPanels%x(j),aPanels%y(j),aPanels%z(j))
	enddo
	aPanels%absVort = 2.0_kreal*Omega*aPanels%z + aPanels%relVort
end subroutine


function HaurwitzStationary4RelvortComponents(x,y,z)
	real(kreal), intent(in) :: x, y, z
	real(kreal) :: HaurwitzStationary4RelVortComponents
	real(kreal), parameter :: alphaS = PI/7.0_kreal
	real(kreal), parameter :: cc = -1.0_kreal
	HaurwitzStationary4RelVortComponents = 2.0_kreal*alphaS*z - 30.0_kreal*cc*&
		cos(4.0_kreal*Longitude(x,y,z))*Legendre54(z)
end function



function HaurwitzStationary4RelVortVector(xyz)
! Calculates the relative vorticity as a function of position for a Rossby-Haurwitz wave
! with zonal wavenumber 4.
	real(kreal), intent(in) :: xyz(3)
	real(kreal) :: HaurwitzStationary4RelVortVector
	! prevailing zonal wind angular velocity
	real(kreal), parameter :: alphaS = PI/7.0_kreal
	real(kreal), parameter :: cc = -1.0_kreal
	HaurwitzStationary4RelVortVector = 2.0_kreal*alphaS*xyz(3) - 30.0_kreal*cc*&
		cos(4.0_kreal*longitude(xyz))*Legendre54(xyz(3))
end function


function RHWave2RelVortVector(xyz)
	real(kreal) :: RHWave2RelVortVector
	real(kreal), intent(in) :: xyz(3)
	real(kreal), parameter :: alphaS = PI/7.0_kreal, &
							  cc = -1.0_kreal

	RHWave2RelVortVector = 2.0_kreal*alphaS*xyz(3) - 30.0_kreal*cc*&
		cos(2.0_kreal*Longitude(xyz))*Legendre52(xyz(3))
end function


function RHWave2RelVortComponents(x,y,z)
	real(kreal) :: RHWave2RelVortComponents
	real(kreal), intent(in) :: x, y, z
	real(kreal), parameter :: alphaS = PI/7.0_kreal, &
							  cc = -1.0_kreal
	RHWave2RelVortComponents = 2.0_kreal*alphaS*z - 30.0_kreal*cc*&
		cos(2.0_kreal*Longitude(x,y,z))*Legendre52(z)
end function




subroutine InitGaussianVortexPanels(aParticles, aPanels)
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


subroutine InitPolarVortexPanels(aParticles,aPanels)
	type(Particles), intent(inout) :: aParticles
	type(Panels), intent(inout) :: aPanels
	real(kreal), allocatable :: pvort(:)
	integer(kint) :: j
	allocate(pvort(aPanels%N))
	pvort = 0.0_kreal
	do j=1,aPanels%N
		if ( .NOT. aPanels%hasChildren(j)) then
			pvort(j) = PolarVortexX(aPanels%x0(:,j))
		endif
	enddo
	gaussConst = sum(pvort*aPanels%area(1:aPanels%N))/(4.0_kreal*PI)
	do j=1,aParticles%N
		aParticles%relVort(j) = PolarVortexX(aParticles%x0(:,j)) - gaussConst
		aParticles%absVort(j) = aParticles%relVort(j) + 2.0_kreal*Omega*aParticles%x0(3,j)
	enddo
	do j=1,aPanels%N
		if ( .NOT. aPanels%hasChildren(j)) then
			aPanels%relVort(j) = PolarVortexX(aPanels%x0(:,j)) - gaussConst
			aPanels%absVort(j) = aPanels%relVort(j) + 2.0_kreal*Omega*aPanels%x0(3,j)
		endif
	enddo
	deallocate(pvort)
end subroutine


subroutine InitGaussianVortexVoronoi(aPanels)
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


subroutine InitTwoVorts(aParticles,aPanels,cent1,cent2,beta1,beta2,strength1,strength2)
	type(Particles), intent(inout) :: aParticles
	type(Panels), intent(inout) :: aPanels
	real(kreal), intent(in) :: cent1(3), cent2(3), beta1, beta2, strength1, strength2
	real(kreal), allocatable :: gVort(:)
	integer(kint) :: j
	allocate(gVort(aPanels%N))
	gVort = 0.0_kreal
	do j=1,aPanels%N
		if ( .NOT. aPanels%hasChildren(j) ) then
			gVort(j) = TwoVortsX(aPanels%x(:,j),cent1,cent2,beta1,beta2,strength1,strength2)
		endif
	enddo
	gaussConst = sum(gVort*aPanels%area(1:aPanels%N))/(4.0_kreal*PI)
	do j=1,aPanels%N
		if (.NOT. aPanels%hasChildren(j)) then
			aPanels%absVort(j) = gVort(j) + 2.0_kreal*Omega*aPanels%x(3,j) - gaussConst
			aPanels%relVort(j) = gVort(j) - gaussConst
		endif
	enddo
	do j=1,aParticles%N
		aParticles%absVort(j) = TwoVortsX(aParticles%x(:,j),cent1,cent2,beta1,beta2,strength1,strength2) &
			- gaussConst + 2.0_kreal*Omega*aParticles%x(3,j)
		aParticles%relVort(j) = TwoVortsX(aParticles%x(:,j),cent1,cent2,beta1,beta2,strength1,strength2) &
			- gaussConst
	enddo
	deallocate(gVort)
end subroutine


subroutine InitTripole(aParticles,aPanels,cent1,cent2,cent3,beta1,beta2,beta3,strength1,strength2,strength3)
	type(Particles), intent(inout) :: aParticles
	type(Panels), intent(inout) :: aPanels
	real(kreal), intent(in) :: cent1(3), cent2(3), cent3(3), beta1, beta2, beta3, strength1, strength2, strength3
	real(kreal), allocatable :: gVort(:)
	integer(kint) :: j
	allocate(gVort(aPanels%N))
	gVort = 0.0_kreal
	do j=1,aPanels%N
		if ( .NOT. aPanels%hasChildren(j) ) then
			gVort(j) = TripoleX(aPanels%x(:,j), cent1,cent2,cent3,beta1,beta2,beta3,strength1,strength2,strength3)
		endif
	enddo
	gaussConst = sum(gVort*aPanels%area(1:aPanels%N))/(4.0_kreal*PI)
	do j=1,apanels%N
		if ( .NOT. aPanels%hasChildren(j)) then
			aPanels%absVort(j) = gVort(j) + 2.0_kreal*Omega*aPanels%x(3,j) - gaussConst
			aPanels%relVort(j) = gVort(j) - gaussConst
		endif
	enddo
	do j=1,aParticles%N
		aParticles%absVort(j) = TripoleX(aParticles%x(:,j),cent1,cent2,cent3,&
			beta1,beta2,beta3,strength1,strength2,strength3) + 2.0_kreal*aParticles%x(3,j) - gaussConst
		aParticles%relVort(j) = TripoleX(aParticles%x(:,j),cent1,cent2,cent3,&
			beta1,beta2,beta3,strength1,strength2,strength3) - gaussConst
	enddo
	deallocate(gVort)
end subroutine


subroutine InitPolarVortexVoronoi(aPanels)
	type(VorPanels), intent(inout) :: aPanels
	integer(kint) :: j
	real(kreal), allocatable :: pVort(:)
	allocate(pvort(aPanels%N))
	pvort = 0.0_kreal
	do j=1,aPanels%N
		pvort(j) = PolarVortexX(aPanels%x0(j),aPanels%y0(j),aPanels%z0(j))
	enddo
	gaussConst = sum(pvort*aPanels%area(1:aPanels%N))/(4.0_kreal*PI)
	aPanels%relVort = pvort - gaussConst
	aPanels%absVort = aPanels%relVort(1:aPanels%N) + 2.0_kreal*Omega*aPanels%z0(1:aPanels%N)
	deallocate(pvort)
end subroutine


function GaussianVortexVector(xyz)
! Outputs the height of a Gaussian function used to calculate the relative vorticity
! asssociated with a Gaussian vortex.  The appropriate constant must be subtracted from
! the output value of this function to give a valid vorticity profile.
	real(kreal) :: GaussianVortexVector
	real(kreal), intent(in) :: xyz(3)
	real(kreal), parameter :: 	lat0 = PI/20.0_kreal,&
								long0 = 0.0_kreal,&
								beta = 4.0_kreal
	real(kreal) :: xCenter(3)
	xCenter(1) = cos(lat0)*cos(long0)
	xCenter(2) = cos(lat0)*sin(long0)
	xCenter(3) = sin(lat0)
	GaussianVortexVector = 4.0_kreal*PI*exp(-2.0_kreal*beta*beta*(1.0_kreal - sum(xyz*xCenter)))
end function


function GaussianVortexComponents(x,y,z)
	real(kreal), intent(in) :: x, y, z
	real(kreal) :: GaussianVortexComponents
	real(kreal), parameter :: lat0 = PI/20.0_kreal
	real(kreal), parameter :: long0 = 0.0_kreal
	real(kreal), parameter :: beta = 4.0_kreal
	real(kreal) :: xc, yc, zc
	xc = cos(lat0)*cos(long0)
	yc = cos(lat0)*sin(long0)
	zc = sin(lat0)
	GaussianVortexComponents= 4.0_kreal*PI*exp(-2.0_kreal*beta*beta*(1.0_kreal - (x*xc + y*yc + z*zc )))
end function



function PolarVortexVector(xyz)
	real(kreal) :: PolarVortexVector
	real(kreal), intent(in) :: xyz(3)
	real(kreal), parameter :: lat0 = PI/2.0_kreal,&
							  long0 = PI, &
							  beta = 2.0_kreal
	real(kreal) :: xCenter(3)
	xCenter(1) = cos(lat0)*cos(long0)
	xCenter(2) = cos(lat0)*sin(long0)
	xCenter(3) = sin(lat0)
	PolarVortexVector = 4.0_kreal*PI*exp(-2.0_kreal*beta*beta*(1.0_kreal - sum(xyz*xCenter)))
end function



function PolarVortexComponents(x,y,z)
	real(kreal) :: PolarVortexComponents
	real(kreal), intent(in) :: x,y,z
	real(kreal), parameter :: lat0 = PI/2.0_kreal,&
							  long0 = PI, &
							  beta = 2.0_kreal
	real(kreal) :: xc, yc, zc
	xc = cos(lat0)*cos(long0)
	yc = cos(lat0)*sin(long0)
	zc = sin(lat0)
	PolarVortexComponents = 4.0_kreal*PI*exp(-2.0_kreal*beta*beta*(1.0_kreal - (x*xc + y*yc + z*zc )))
end function


function Legendre54(z)
! Calculates the Legendre polynomial P_5^4
	real(kreal), intent(in):: z
	real(kreal) :: Legendre54
	Legendre54 = z*(-1.0_kreal + z*z)*(-1.0_kreal + z*z)
end function


function Legendre52(z)
	real(kreal) :: Legendre52
	real(kreal), intent(in) :: z
	Legendre52 = -(-1.0_kreal + z*z)*(-z+3.0_kreal*z*z*z)
end function


subroutine InitManyVortsPanels(aParticles, aPanels)
	type(Particles), intent(inout) :: aParticles
	type(Panels), intent(inout) :: aPanels
	real(kreal), allocatable :: gvort(:)
	integer(kint) :: j
	allocate(gvort(aPanels%N))
	gvort = 0.0_kreal
	do j=1,aPanels%N
		if ( .NOT. aPanels%hasChildren(j) ) then
			gvort(j) = ManyVortsVector(aPanels%x(:,j))
		endif
	enddo
	gaussConst = sum(gvort*aPanels%area(1:aPanels%N))/(4.0_kreal*PI)
	do j=1,aParticles%N
		aParticles%absVort(j) = ManyVortsVector(aParticles%x0(:,j)) - gaussConst + &
				 2.0_kreal*Omega*aParticles%x0(3,j)
		aParticles%relVort(j) = aParticles%absVort(j) - 2.0_kreal*Omega*aParticles%x(3,j)
	enddo
	do j=1,aPanels%N
		if (.NOT. aPanels%hasChildren(j) ) then
			aPanels%absVort(j) = ManyVortsVector(aPanels%x0(:,j)) - gaussConst  + &
				2.0_kreal*Omega*aPanels%x0(3,j)
			aPanels%relVort(j) = aPanels%absVort(j) - 2.0_kreal*Omega*aPanels%x(3,j)
		else
			aPanels%relVort(j) = 0.0_kreal
			aPanels%absVort(j) = 0.0_kreal
		endif
	enddo
	deallocate(gvort)
end subroutine


subroutine InitManyVortsVoronoi(aPanels)
	type(VorPanels), intent(inout) :: aPanels
	integer(kint) :: j
	real(kreal), allocatable :: gVort(:)
	allocate(gVort(aPanels%N))
	gVort = 0.0_kreal
	do j=1,aPanels%N
		gVort(j) = ManyVortsComponents(aPanels%x(j),aPanels%y(j),aPanels%z(j))
	enddo
	gaussConst = sum(gVort*aPanels%area(1:aPanels%N))/(4.0_kreal*PI)
	aPanels%relVort(1:aPanels%N) = gVort - gaussConst
	aPanels%absVort(1:aPanels%N) = aPanels%relVort(1:aPanels%N) + 2.0_kreal*Omega*aPanels%z(1:aPanels%N)
	deallocate(gVort)
end subroutine


function ManyVortsVector(xyz)
	real(kreal) :: ManyVortsVector
	real(kreal), intent(in) :: xyz(3)
	real(kreal) :: lat, lon
	lat = Latitude(xyz)
	lon = Longitude(xyz)
	ManyVortsVector = (4.0_kreal*PI*sin(4.0_kreal*lon)*sin(8.0_kreal*lat) + &
					   1.6_kreal*PI*cos(3.0_kreal*lon)*cos(6.0_kreal*lat) + &
					   1.2_kreal*PI*cos(5.0_kreal*lon)*cos(10.0_kreal*lat) + &
					   0.08_kreal*PI*sin(lon) + 0.08_kreal*PI*sin(2.0_kreal*lat))* &
					   exp(-0.5_kreal*lat**8)
end function


function ManyVortsComponents(x,y,z)
	real(kreal) :: ManyVortsComponents
	real(kreal), intent(in) :: x, y, z
	real(kreal) :: lat, lon
	lat = Latitude(x,y,z)
	lon = Longitude(x,y,z)
	ManyVortsComponents = (4.0_kreal*PI*sin(4.0_kreal*lon)*sin(8.0_kreal*lat) + &
					   1.6_kreal*PI*cos(3.0_kreal*lon)*cos(6.0_kreal*lat) + &
					   1.2_kreal*PI*cos(5.0_kreal*lon)*cos(10.0_kreal*lat) + &
					   0.08_kreal*PI*sin(lon) + 0.08_kreal*PI*sin(2.0_kreal*lat))*&
					   exp(-0.5_kreal*lat**8)
end function


function Juckes_A(t)
	real(kreal) :: Juckes_A
	real(kreal), intent(in) :: t
	real(kreal), parameter :: tfull = 4.0_kreal, &
							  tend = 15.0_kreal
	real(kreal) :: b

	b = PI/tfull

	if ( 0.0_kreal <= t .AND. t < tfull ) then
		Juckes_A = 0.5_kreal*(1.0_kreal - cos(b*t))
	elseif ( tfull <= t .AND. t < (tend-tfull) ) then
		Juckes_A = 1.0_kreal
	elseif ( (tend-tfull) <= t .AND. (t < tend) ) then
		Juckes_A = 0.5_kreal*(1.0_kreal - cos(b*(t-(tend-tfull)) + PI))
	else
		Juckes_A = 0.0_kreal
	endif
end function


function Juckes_APrime(t)
	real(kreal) :: Juckes_APrime
	real(kreal), intent(in) :: t
	real(kreal), parameter :: tfull = 4.0_kreal, &
							  tend = 15.0_kreal
	real(kreal) :: b

	b = PI/tfull

	if ( 0.0_kreal <= t .AND. t < tfull ) then
		Juckes_APrime = 0.5_kreal*b*sin(b*t)
	elseif ( tfull <= t .AND. t < (tend-tfull) ) then
		Juckes_APrime = 0.0_kreal
	elseif ( (tend-tfull) <= t .AND. t < tend ) then
		Juckes_APrime = 0.5_kreal*b*sin(b*(t-(tend-tfull)) + PI)
	else
		Juckes_APrime = 0.0_kreal
	endif
end function


function Juckes_RelVortVector(xyz)
	real(kreal) :: Juckes_RelVortVector
	real(kreal), intent(in) :: xyz(3)
	real(kreal) :: lat
	real(kreal), parameter :: 	beta = 3.0_kreal/2.0_kreal, &
								th0 = 15_kreal*PI/32.0_kreal
	lat = Latitude(xyz)
	Juckes_RelVortVector = -PI*(cos(lat)*(-2.0_kreal*beta*beta*( cos(th0)*sin(lat) - sin(th0)*cos(lat) )) &
		- sin(lat))*exp(-2.0_kreal*beta*beta*(1.0_kreal-cos(th0)*cos(lat)-sin(th0)*sin(lat)))
end function


function Juckes_BVector(xyz)
	real(kreal) :: Juckes_BVector
	real(kreal), intent(in) :: xyz(3)
	real(kreal), parameter :: th1 = PI/3.0_kreal
	real(kreal) :: lat, tan1sq
	lat = Latitude(xyz)
	if ( lat <= 0.0_kreal ) then
		Juckes_BVector = 0.0_kreal
	else
		tan1sq = tan(th1)*tan(th1)
		Juckes_BVector = tan1sq/(tan(lat)*tan(lat))*exp(1.0_kreal - tan1sq/(tan(lat)*tan(lat)))
	endif
end function

function Juckes_Bscalar(lat)
	real(kreal) :: Juckes_Bscalar
	real(kreal), intent(in) :: lat
	real(kreal), parameter :: th1 = PI/3.0_kreal
	real(kreal) :: tan1sq
	if ( lat <= 0.0_kreal ) then
		Juckes_BScalar = 0.0_kreal
	else
		tan1sq = tan(th1)*tan(th1)
		Juckes_BScalar = tan1sq/(tan(lat)*tan(lat))*exp(1.0_kreal - tan1sq/(tan(lat)*tan(lat)))
	endif
end function


function Juckes_BPrimeVector(xyz)
	real(kreal) :: Juckes_BPrimeVEctor
	real(kreal), intent(in) :: xyz(3)
	real(kreal) :: lat, tan1sq, fac1, fac2
	real(kreal), parameter :: th1 = PI/3.0_kreal

	lat = Latitude(xyz)
	if ( lat < 0.0_kreal ) then
		Juckes_BPrimeVector = 0.0_kreal
	else
		tan1sq = tan(th1)*tan(th1)
		fac1 = 1.0_kreal/(sin(lat)*sin(lat)*tan(lat))
		fac2 = 1.0_kreal/(tan(lat)*tan(lat)*tan(lat)*sin(lat)*sin(lat))
		Juckes_BPrimeVector = (-2.0_kreal*tan1sq*fac1 + 2.0_kreal*tan1sq*tan1sq*fac2)*&
			exp(1.0_kreal-tan1sq/(tan(lat)*tan(lat)))
	endif
end function


function Juckes_BPrimeScalar(lat)
	real(kreal) :: Juckes_BPrimeScalar
	real(kreal), intent (in) :: lat
	real(kreal) ::  tan1sq, fac1, fac2
	real(kreal), parameter :: th1 = PI/3.0_kreal
	if ( lat < 0.0_kreal ) then
		Juckes_BPrimeScalar = 0.0_kreal
	else
		tan1sq = tan(th1)*tan(th1)
		fac1 = 1.0_kreal/(sin(lat)*sin(lat)*tan(lat))
		fac2 = 1.0_kreal/(tan(lat)*tan(lat)*tan(lat)*sin(lat)*sin(lat))
		Juckes_BPrimeScalar = (-2.0_kreal*tan1sq*fac1 + 2.0_kreal*tan1sq*tan1sq*fac2)*&
			exp(1.0_kreal-tan1sq/(tan(lat)*tan(lat)))
	endif
end function


function ZonalMeanRelVortVector(xyz)
	real(kreal) :: ZonalMeanRelVortVector
	real(kreal), intent(in) :: xyz(3)
	real(kreal) :: lat
	lat = Latitude(xyz)
	ZonalMeanRelVortVector = 15.0_kreal*cos(lat)*(37.0_kreal*xyz(3) - 35.0_kreal*sin(3.0_kreal*lat))/24.0_kreal
end function


function ZonalMeanRelVortComponents(x,y,z)
	real(kreal) :: ZonalMeanRelVortComponents
	real(kreal), intent(in) :: x, y, z
	real(kreal) :: lat
	lat = Latitude(x,y,z)
	ZonalMeanRelVortComponents = 15.0_kreal*cos(lat)*(37.0_kreal*z - 35.0_kreal*sin(3.0_kreal*lat))/24.0_kreal
end function


subroutine InitZonalMeanPanels(aParticles,aPanels)
	type(Particles), intent(inout) :: aParticles
	type(Panels), intent(inout) :: aPanels
	integer(kint) :: j
	do j=1,aParticles%N
		aParticles%relVort(j) = ZonalMeanRelVortVector(aParticles%x0(:,j))
		aParticles%absVort(j) = aParticles%relVort(j) + 2.0_kreal*Omega*aparticles%x0(3,j)
	enddo
	do j=1,aPanels%N
		if ( .NOT. aPanels%hasChildren(j) ) then
			aPanels%relVort(j) = ZonalMeanRelVortVector(aPanels%x0(:,j))
			aPanels%absVort(j) = aPanels%relVort(j) + 2.0_kreal*Omega*aPanels%x0(3,j)
		endif
	enddo
end subroutine


subroutine InitZonalMeanVoronoi(aPanels)
	type(VorPanels), intent(inout) :: aPanels
	integer(kint) :: j
	do j=1,aPanels%N
		aPanels%relVort(j) = ZonalMeanRelVortComponents(aPanels%x0(j),aPanels%y0(j),aPanels%z0(j))
		aPanels%absVort(j) = aPanels%relVort(j) + 2.0_kreal*Omega*aPanels%z0(j)
	enddo
end subroutine


subroutine InitJetVoronoi(aPanels,theta0,beta,perturbAmp,perturbWaveNum)
	type(VorPanels), intent(inout) :: aPanels
	integer(kint), intent(in) :: perturbWaveNum
	real(kreal), intent(in) :: beta, theta0, perturbAmp
	real(kreal), allocatable :: relVortStar(:)
	integer(kint) :: j
	! Determine value of gaussConst
	allocate(relVortStar(aPanels%N))
	relVortStar = 0.0_kreal
	gaussConst = 0.0_kreal
	do j=1,aPanels%N
		relVortStar(j) = JetRelVortComponents(aPanels%x0(j),aPanels%y0(j),aPanels%z0(j),theta0,beta,perturbAmp,perturbWaveNum)
	enddo
	gaussConst = sum(relVortStar*aPanels%area(1:aPanels%N))/(4.0_kreal*PI)

	! Set vorticity on panels
	do j=1,aPanels%N
		aPanels%relVort(j) = JetRelVortComponents(aPanels%x0(j),aPanels%y0(j),aPanels%z0(j),theta0,beta,perturbAmp,perturbWaveNum) - gaussConst
		aPanels%absVort(j) = aPanels%relVort(j) + 2.0_kreal*Omega*aPanels%z0(j)
	enddo
	deallocate(relVortStar)
end subroutine


subroutine InitJetPanels(aParticles,aPanels,theta0,beta,perturbAmp,perturbWaveNum)
	type(Particles), intent(inout) :: aParticles
	type(Panels), intent(inout) :: aPanels
	integer(kint), intent(in) :: perturbWaveNum
	real(kreal), intent(in) :: beta, theta0, perturbAmp
	real(kreal), allocatable :: relVortStar(:)
	integer(kint) :: j, k

	! Determine the value of the Gauss Constant
	allocate(relVortStar(aPanels%N))
	relVortStar = 0.0_kreal
	gaussConst = 0.0_kreal
	do j=1,aPanels%N
		if (.NOT. aPanels%hasChildren(j) ) then
			relVortStar(j) = JetRelVortVector(aPanels%x0(:,j),theta0,beta,perturbAmp,perturbWaveNum)
		endif
	enddo
	gaussConst = sum(relVortStar*aPanels%area(1:aPanels%N))/(4.0_kreal*PI)

	! Set vorticity on panels
	do j=1,aPanels%N
		if (.NOT. aPanels%hasChildren(j) ) then
			aPanels%relVort(j) = relVortStar(j) - gaussConst
			aPanels%absVort(j) = relVortStar(j) - gaussConst + 2.0_kreal*Omega*aPanels%x0(3,j)
		endif
	enddo
	! Set vorticity on particles
	do j=1,aParticles%N
		aParticles%relVort(j) = JetRelVortVector(aParticles%x0(:,j),theta0,beta,perturbAmp,perturbWaveNum) - gaussConst
		aParticles%absVort(j) = aParticles%relVort(j) + 2.0_kreal*Omega*aParticles%x0(3,j)
	enddo

!	theta2 = theta1 + jetThickness
!	do j=1,aPanels%N
!		if ( .NOT. aPanels%hasChildren(j) ) then
!			relVortStar(j) = JetRelVortVector(aPanels%x0(:,j),nVorts,beta,theta1,jetThickness,latPerturb)
!		endif
!	enddo
!	gaussConst = sum(relVortStar*aPanels%area(1:aPanels%N))/(4.0_kreal*PI)
!
!	! Set vorticity on panels
!	do j=1,aPanels%N
!		if (.NOT. aPanels%hasChildren(j) ) then
!			aPanels%relVort(j) = relVortStar(j) - gaussConst
!			aPanels%absVort(j) = aPanels%relVort(j) + 2.0_kreal*Omega*aPanels%x0(3,j)
!		endif
!	enddo
!	! Set vorticity on particles
!	do j=1,aParticles%N
!		aParticles%relVort(j) = JetRelVortVector(aParticles%x0(:,j),nVorts,beta,theta1,jetThickness,latPerturb) - gaussConst
!	enddo
!	aParticles%absVort(1:aParticles%N) = aParticles%relVort(1:aParticles%N) + 2.0_kreal*Omega*aParticles%x0(3,1:aParticles%N)
	deallocate(relVortStar)
end subroutine


subroutine InitStratosphere(aParticles,aPanels)
	type(Particles), intent(inout) :: aParticles
	type(Panels), intent(inout) :: aPanels
	integer(kint) :: j
	real(kreal), allocatable :: relVortStar(:)

	! Determine the value of the Gauss constant
	allocate(relVortStar(aPanels%N))
	relVortStar = 0.0_kreal
	gaussConst = 0.0_kreal
	do j=1,aPanels%N
		if (.NOT. aPanels%hasChildren(j) ) then
			relVortStar(j) = StratosphereRelVortX(aPanels%x0(:,j))
		endif
	enddo
	gaussConst = sum(relVortStar*aPanels%area(1:aPanels%N))/(4.0_kreal*PI)

	! Set vorticity on panels
	do j=1,aPanels%N
		if ( .NOT. aPanels%hasChildren(j) ) then
			aPanels%relVort(j) = relVortStar(j) - gaussConst
			aPanels%absVort(j) = relVortStar(j) - gaussConst + 2.0_kreal*Omega*aPanels%x0(3,j)
		endif
	enddo
	! Set vorticity on particles
	do j=1,aParticles%N
		aParticles%relVort(j) = StratosphereRelVortX(aParticles%x0(:,j)) - gaussConst
		aParticles%absVort(j) = StratosphereRelVortX(aParticles%x0(:,j)) - gaussConst &
			+ 2.0_kreal*Omega*aParticles%x0(3,j)
	enddo

	deallocate(relVortStar)
end subroutine



function JetRelVortVector(xyz,theta0,beta,perturbAmp,perturbWaveNum,strengthPerturbAmp,strengthPerturbWaveNum)
	real(kreal) :: JetRelVortVector
	real(kreal), intent(in) :: xyz(3), theta0, beta, perturbAmp
	integer(kint), intent(in) :: perturbWaveNum
	real(kreal), intent(in), optional :: strengthPerturbAmp
	integer(kint), intent(in), optional :: strengthPerturbWaveNum
	real(kreal) :: lat, lon, thetaP, strAmp
	integer(kint) :: strWN
	if ( present(strengthPerturbAmp) ) then
		strAmp = strengthPerturbAmp
		strWN = strengthPerturbWaveNum
	else
		strAmp = 0.0_kreal
		strWN = 4
	endif
	lat = Latitude(xyz)
	lon = Longitude(xyz)
	thetaP = theta0 + perturbAmp*cos(real(perturbWaveNum,kreal)*lon)
!	JetRelVortVector = -PI/2.0_kreal*(1.0_kreal + strAmp*cos(real(strWN,kreal)*lon))*&
!		(cos(lat)*(-2.0_kreal*beta*beta*(cos(thetaP)*sin(lat) - &
!		sin(thetaP)*cos(lat))) - sin(lat))*exp(-2.0_kreal*beta*beta*&
!		(1.0_kreal - cos(lat)*cos(thetaP)-sin(lat)*sin(thetaP)))

	JetRelVortVector = 300.0_kreal*sin(lat - thetaP)*exp(-2.0_kreal*beta*beta*(1.0_kreal - cos(lat-thetaP)))
end function



!function JetRelVortVector(xyz,nVorts,beta,theta1,jetThickness,latPerturb)
!	real(kreal) :: JetRelVortVector
!	integer(kint), intent(in) :: nVorts
!	real(kreal), intent(in) :: xyz(3), beta,theta1, jetThickness, latPerturb
!	real(kreal) :: theta2
!	integer(kint) :: k
!	real(kreal) :: xCent1(3), xCent2(3), lat1, lat2
!
!	theta2 = theta1 + jetThickness
!	JetRelVortVector = 0.0_kreal
!	do k=0,nVorts-1
!		lat1 = theta1 + latPerturb*sin(12.0_kreal*k*2.0_kreal*PI/nVorts)
!		lat2 = theta2 + latPerturb*sin(12.0_kreal*k*2.0_kreal*PI/nVorts)
!
!		xCent1(1) = cos(lat1)*cos(k*2.0_kreal*PI/nVorts)
!		xCent1(2) = cos(lat1)*sin(k*2.0_kreal*PI/nVorts)
!		xCent1(3) = sin(lat1)
!
!		xCent2(1) = cos(lat2)*cos(k*2.0_kreal*PI/nVorts)
!		xCent2(2) = cos(lat2)*sin(k*2.0_kreal*PI/nVorts)
!		xCent2(3) = sin(lat2)
!
!		JetRelVortVector = JetRelVortVector - GeneralGaussian(xyz,xCent1,beta)  &
!			+ GeneralGaussian(xyz,xCent2,beta) - gaussConst
!	enddo
!end function

function JetRelVortComponents(x,y,z,theta0,beta,perturbAmp,perturbWaveNum)
	real(kreal) :: JetRelVortComponents
	integer(kint), intent(in) :: perturbWaveNum
	real(kreal), intent(in) :: x,y,z, beta, theta0, perturbAmp
	real(kreal) :: xyz(3)
	xyz = [x,y,z]
	JetRelVortComponents = JetRelVortVector(xyz,theta0,beta,perturbAmp,perturbWaveNum)
end function

function GeneralGaussian(xyz,xCenter,beta, strength)
	real(kreal) :: GeneralGaussian
	real(kreal), intent(in) :: xyz(3), xCenter(3), beta
	real(kreal), intent(in), optional :: strength
	real(kreal) :: amp
	if ( present(strength) ) then
		amp = strength
	else
		amp = 4.0_kreal*PI
	endif
	GeneralGaussian = amp*exp(-2.0_kreal*beta*beta*(1.0_kreal - sum(xyz*xCenter)))
end function




function TwoGaussVortsVector(xyz,cent1,cent2,beta1,beta2,strength1,strength2)
	real(kreal) :: TwoGaussVortsVector
	real(kreal), intent(in) :: xyz(3), cent1(3), cent2(3), beta1, beta2, strength1, strength2
	real(kreal) :: g1, g2
	g1 = GeneralGaussian(xyz,cent1,beta1,strength1)
	g2 = GeneralGaussian(xyz,cent2,beta2,strength2)
	TwoGaussVortsVector = g1 + g2
end function


function TripoleVortVector(xyz,cent1,cent2,cent3,beta1,beta2,beta3,strength1,strength2,strength3)
	real(kreal) :: TripoleVortVector
	real(kreal), intent(in) :: xyz(3), cent1(3), cent2(3), cent3(3)
	real(kreal), intent(in) :: beta1, beta2, beta3
	real(kreal), intent(in) :: strength1, strength2, strength3
	real(kreal) :: g1, g2, g3

	if ( strength3 >= 0.0_kreal) then
		print *, "Tripole error: strength3 must be negative."
		return
	endif

	g1 = GeneralGaussian(xyz,cent1,beta1,strength1)
	g2 = GeneralGaussian(xyz,cent2,beta2,strength2)
	g3 = GeneralGaussian(xyz,cent3,beta3,strength3)

	TripoleVortVector = g1 + g2 + g3

end function


subroutine InitSlottedCylinders(aParticles,aPanels)
	type(Particles), intent(inout) :: aParticles
	type(Panels), intent(inout) :: aPanels
	integer(kint) :: j
	if ( (associated(aParticles%tracer)) .AND. (associated(aPanels%tracer)) ) then
		do j=1,aParticles%N
			aParticles%tracer(j,1) = SlottedCylindersX(aParticles%x0(:,j))
		enddo
		do j=1,aPanels%N
			if (.NOT. aPanels%hasChildren(j)) then
				aPanels%tracer(j,1) = SlottedCylindersX(aPanels%x0(:,j))
			else
				aPanels%tracer(j,1) = 0.0_kreal
			endif
		enddo
	else
		print *, "SlottedCylinders ERROR : tracer memory not allocated."
		return
	endif
end subroutine


function SlottedCylindersX(xyz)
	! This function assigns a value to a coordinate on the unit sphere corresponding to the
	! slotted cylinders mass tracer distribution in Nair & Lauritzen (2010).
	real(kreal), intent(in) :: xyz(3)
	real(kreal) :: SlottedCylindersX
	real(kreal), parameter :: xx1 = -0.866025403784439_kreal, yy1 = 0.5_kreal, zz1 = 0.0_kreal
	real(kreal), parameter :: xx2 = -0.866025403784439_kreal, yy2 = -0.5_kreal, zz2 = 0.0_kreal
	real(kreal), parameter :: lat1 = 0.0_kreal, long1 = 5.0_kreal*PI/6.0_kreal
	real(kreal), parameter :: lat2 = 0.0_kreal, long2 = 7.0_kreal*PI/6.0_kreal
	real(kreal), parameter :: RR = 0.5_kreal, b = 0.1_kreal, c = 1.0_kreal
	real(kreal) :: r1, r2, lat, long
	lat = latitude(xyz)
	long = longitude(xyz)
	r1 = SphereDistance(xyz,[xx1,yy1,zz1])
	r2 = SphereDistance(xyz,[xx2,yy2,zz2])
	slottedCylindersX = b
	if ( r1 <= RR) then
		if ( abs(long-long1)>=RR/6.0_kreal) then
			slottedCylindersX = c
		else
			if ( lat - lat1 < -5.0_kreal*RR/12.0_kreal) slottedCylindersX = c
		endif
	endif
	if ( r2 <= RR ) then
		if ( abs(long-long2) >= RR/6.0_kreal) then
			slottedCylindersX = c
		else
			if ( lat-lat2 > 5.0_kreal*RR/12.0_kreal) slottedCylindersX = c
		endif
	endif
end function


subroutine InitBlockM(aParticles,aPanels)
	type(Particles), intent(inout) :: aParticles
	type(Panels), intent(inout) :: aPanels
	integer(kint) :: j
	if ( (associated(aParticles%tracer)) .AND. (associated(aPanels%tracer)) ) then
		do j=1,aParticles%N
			aParticles%tracer(j,1) = BlockMVector(aParticles%x0(:,j))
		enddo
		do j=1,aPanels%N
			if (.NOT. aPanels%hasChildren(j)) then
				aPanels%tracer(j,1) = BlockMVector(aPanels%x0(:,j))
			else
				aPanels%tracer(j,1) = 0.0_kreal
			endif
		enddo
	else
		print *, "InitBlockM ERROR : tracer memory not allocated."
		return
	endif
end subroutine


function BlockMVector(xyz)
	real(kreal) :: BlockMVector
	real(kreal), intent(in) :: xyz(3)
	real(kreal), parameter :: unit = PI/6.0_kreal
							  !rightLon = PI/6.0_kreal,&
							 !bottomLat = -7.0_kreal*PI/60.0_kreal ,&
							 ! topLat = 7.0_kreal*PI/60.0_kreal
	real(kreal) :: lat, lon, x, y

	lat = Latitude(xyz)
	lon = atan2(xyz(2),xyz(1))
	x = lon/unit
	y = lat/unit
	BlockMVector = ScaledBlockM(x,y)
end function


function ScaledBlockM(x,y)
	real(kreal) :: ScaledBlockM
	real(kreal), intent(in) :: x, y
	real(kreal), parameter :: m = 1.6625_kreal, &
							  vert1 = 0.38421_kreal, &
							  horiz1 = 0.315789_kreal, &
							  horiz2 = 0.8421_kreal,&
							  horiz3 = 0.4736842,&
							  top = 0.7_kreal, &
							  right = 1.0_kreal, &
							  bottom = -0.7_kreal, &
							  left = -1.0_kreal, &
							  one = 1.0_kreal, &
							  zero = 0.0_kreal

	ScaledBlockM = 0.1_kreal
	! Left Foot
	if ( (( y < -vert1) .AND. ( y > bottom )) .AND. ( (x < - horiz1) .AND. (x > left))) ScaledBlockM = one
	! Right Foot
	if ( ( (y < -vert1) .AND. ( y > bottom )) .AND. ( (x > horiz1) .AND. (x < right))) ScaledBlockM = one
	! Left Leg
	if ( ( (y > bottom) .AND. (y < top) ) .AND. ( (x > -horiz2) .AND. (x < -horiz3 ))) ScaledBlockM = one
	! Right Leg
	If ( ( (y > bottom) .AND. (y<top)) .AND. ( (x>horiz3) .AND. (x<horiz2))) ScaledBlockM = one
	! Left top
	If ( ( (y > vert1) .AND. (y<top)) .AND. ( (x<-horiz1) .AND. (x> left))) ScaledBlockM = one
	! Right top
	If ( ( (y > vert1) .AND. (y<top)) .AND. ( (x<right) .AND. (x>horiz1))) ScaledBlockM = one
	! Left triangle
	If ( (x <= zero) .AND. ( y < top) ) then
		if (( y > -m*x - top) .AND. ( y < -m*(x + horiz1) + top)) ScaledBlockM = one
	endif
	If ( (x >= zero) .AND. ( y<top)) then
		if ((y > m*x - top) .AND. (y<m*(x-horiz1)+top)) ScaledBlockM = one
	endif
end function


subroutine InitCosineBells(aParticles,aPanels)
	type(Particles), intent(inout) :: aParticles
	type(Panels), intent(inout) :: aPanels
	integer(kint) :: j
	if ( (associated(aParticles%tracer)) .AND. (associated(aPanels%tracer))) then
	do j=1,aParticles%N
		aParticles%tracer(j,1) = CosineBellsX(aParticles%x0(:,j))
	enddo
	do j=1,aPanels%N
		if (.NOT. aPanels%hasChildren(j)) then
			aPanels%tracer(j,1) = CosineBellsX(aPanels%x0(:,j))
		else
			aPanels%tracer(j,1) = 0.0_kreal
		endif
	enddo
	else
		print *,"CosineBells ERROR : tracer memory not allocated."
		return
	endif
end subroutine


function CosineBellsX(xyz)
	! This function assigns a value to a coordinate on the unit sphere corresponding to the
	! cosine bells mass tracer distribution in Lair & Lauritzen (2010).
	real(kreal), intent(in) :: xyz(3)
	real(kreal) :: CosineBellsX
	real(kreal), parameter :: xx1 = -0.866025403784439_kreal, yy1 = 0.5_kreal, zz1 = 0.0_kreal
	real(kreal), parameter :: xx2 = -0.866025403784439_kreal, yy2 = -0.5_kreal, zz2 = 0.0_kreal
	real(kreal), parameter :: hmax = 1.0_kreal, RR = 0.5_kreal, b = 0.1_kreal, c = 0.9_kreal
	real(kreal) :: r1, r2, h1, h2
	r1 = SphereDistance(xyz,[xx1,yy1,zz1])
	r2 = SphereDistance(xyz,[xx2,yy2,zz2])
	h1 = hmax*(1.0_kreal + cos(PI*r1/RR))/2.0_kreal
	h2 = hmax*(1.0_kreal + cos(PI*r2/RR))/2.0_kreal
	if ( r1 < RR) then
		CosineBellsX = b + c*h1
	elseif (r2 < RR) then
		CosineBellsX = b + c*h2
	else
		CosineBellsX = b
	endif
end function


subroutine InitGaussianHills(aParticles,aPanels)
	type(Particles), intent(inout) :: aParticles
	type(Panels), intent(inout) :: aPanels
	integer(kint) :: j
	if ( ( associated(aParticles%tracer)) .AND. ( associated(aPanels%tracer))) then
		do j=1,aParticles%N
			aParticles%tracer(j,1) = GaussianHillsX(aParticles%x0(:,j))
		enddo
		do j=1,aPanels%N
			if (.NOT. aPanels%hasChildren(j)) then
				aPanels%tracer(j,1) = GaussianHillsX(aPanels%x0(:,j))
			else
				aPanels%tracer(j,1) = 0.0_kreal
			endif
		enddo
	else
		print *,"GaussianHills ERROR : tracer memory not allocated."
		return
	endif
end subroutine


function GaussianHillsX(xyz)
	! This function assigns a value to a coordinate on the unit sphere corresponding to the
	! Gaussian hills mass tracer distribution in Lair & Lauritzen (2010).
	real(kreal), intent(in) :: xyz(3)
	real(kreal) :: GaussianHillsX
	real(kreal), parameter :: xx1 = -0.866025403784439_kreal, yy1 = 0.5_kreal, zz1 = 0.0_kreal
	real(kreal), parameter :: xx2 = -0.866025403784439_kreal, yy2 = -0.5_kreal, zz2 = 0.0_kreal
	real(kreal), parameter :: hmax = 0.95_kreal, b = 5.0_kreal
	real(kreal) :: h1, h2
	h1 = hmax*exp( - b * ( (xyz(1)-xx1)*(xyz(1)-xx1) + (xyz(2)-yy1)*(xyz(2)-yy1) + (xyz(3)-zz1)*(xyz(3)-zz1) ) )
	h2 = hmax*exp( - b * ( (xyz(1)-xx2)*(xyz(1)-xx2) + (xyz(2)-yy2)*(xyz(2)-yy2) + (xyz(3)-zz2)*(xyz(3)-zz2) ) )
	GaussianHillsX = h1 + h2
end function


subroutine InitWilliamsonCosineBell(aParticles,aPanels)
	type(Particles), intent(inout) :: aParticles
	type(Panels), intent(inout) :: aPanels
	integer(kint) :: j
	if ( ( associated(aParticles%tracer)) .AND. ( associated(aPanels%tracer))) then
		do j=1,aParticles%N
			aParticles%tracer(j,4) = WilliamsonCosineBellX(aParticles%x0(:,j))
		enddo
		do j=1,aPanels%N
			if (.NOT. aPanels%hasChildren(j)) then
				aPanels%tracer(j,4) = WilliamsonCosineBellX(aPanels%x0(:,j))
			else
				aPanels%tracer(j,4) = 0.0_kreal
			endif
		enddo
	else
		print *,"WilliamsonCosineBell ERROR : tracer memory not allocated."
		return
	endif
end subroutine


function WilliamsonCosineBellX(xyz)
	! This function outputs a value from the cosine bell distribution from Williamson et al. (1992)
	! test case 1, rescaled to the unit sphere.
	! Input : xyz = coordinate on the unit sphere
	! Output : WilliamsonCosineBell = tracer value at xyz
	! Parameters : 	lat0, long0 = latitude and longitude coordinates of cosine bell center
	!				RR = width parameter of cosine bell
	!				h0 = max heights of cosine bell
	real(kreal), intent(in) :: xyz(3)
	real(kreal) :: WilliamsonCosineBellX
	real(kreal), parameter :: lat0 = 0.0_KREAL, long0 = 0.0_KREAL, RR = 1.0_KREAL/3.0_KREAL, h0 = 0.0016_KREAL
	real(kreal) :: x0, y0, z0, r
	x0 = cos(lat0)*cos(long0)
	y0 = cos(lat0)*sin(long0)
	z0 = sin(lat0)
	r = SphereDistance(xyz,[x0,y0,z0]);
	if ( r < RR) then
		WilliamsonCosineBellX = (h0/2.0_KREAL)*(1.0_KREAL + cos(PI*r/RR))
	else
		WilliamsonCosineBellX = 0.0_kreal
	endif
end function





function StratosphereRelVortX(xyz)
	real(kreal) :: StratosphereRelVortX
	real(kreal), intent(in) :: xyz(3)
	real(kreal), parameter :: th0 = 15.0_kreal*PI/32.0_kreal, &
							  beta = 1.5_kreal
	real(kreal) :: lat
	lat = Latitude(xyz)
	StratosphereRelVortX = -PI*( cos(lat)*(-2.0_kreal*beta*beta*(cos(th0)*sin(lat)-sin(th0)*cos(lat))) &
		- sin(lat))*exp(-2.0_kreal*beta*beta*(1.0_kreal - cos(th0)*cos(lat) - sin(th0)*sin(lat)))
end function


function Juckes_Forcing(xyz,t)
	real(kreal) :: Juckes_Forcing
	real(kreal), intent(in) :: xyz(3), t
	real(kreal), parameter :: F0 = 0.3_kreal*4.0_kreal*PI
	real(kreal) :: lat, lon
	lat = Latitude(xyz)
	lon = Longitude(xyz)
	Juckes_Forcing = F0*Juckes_A(t)*Juckes_B(lat)*cos(lon)
end function


function Juckes_ForcingDerivative(xyz, dX, t)
	real(kreal) :: Juckes_ForcingDerivative
	real(kreal), intent(in) :: xyz(3), dX(3), t
	real(kreal), parameter :: F0 = 0.3_kreal*4.0_kreal*PI, &
							  zeroTol = 1e-12
	real(kreal) :: lat, lon, raxis2, u, v
	lat = Latitude(xyz)
	lon = Longitude(xyz)
	raxis2 = xyz(1)*xyz(1) + xyz(2)*xyz(2)
	if ( raxis2 <= zeroTol ) then
		Juckes_ForcingDerivative = 0.0_kreal
	else
		u = (-xyz(2)*dX(1) + xyz(1)*dX(2))/sqrt(raxis2)
		v = (-xyz(1)*xyz(3)*dX(1) - xyz(2)*xyz(3)*dX(2))/sqrt(raxis2) + sqrt(raxis2)*dX(3)
		Juckes_ForcingDerivative = F0*Juckes_APrime(t)*Juckes_B(lat)*cos(lon) - & ! dF/dt
			F0 * u * Juckes_A(t)*Juckes_B(lat)*sin(lon)/sqrt(raxis2) + & !
			F0 * v * Juckes_A(t)*Juckes_BPrime(lat)*cos(lon)
	endif
end function

subroutine InitSolidBody(aParticles,aPanels,rotRate)
	type(Particles), intent(inout) :: aParticles
	type(Panels), intent(inout) :: aPanels
	real(kreal), intent(in) :: rotRate
	integer(kint) :: j
	! exchange role of absVort and relVort for this test case only
	do j=1,aParticles%N
		aParticles%absVort(j) = 0.0_kreal
		aParticles%relVort(j) = 2.0_kreal*rotRate*aParticles%x0(3,j)
	enddo
	do j=1,aPanels%N
		if ( .NOT. aPanels%hasChildren(j)) then
			aPanels%absVort(j) = 0.0_kreal
			aPanels%relVort(j) = 2.0_kreal*rotRate*aPanels%x0(3,j)
		endif
	enddo
end subroutine

function SolidBodyX(xyz,rotationRate)
	real(kreal) :: solidBodyX
	real(kreal), intent(in) :: xyz(3), rotationRate
	SolidBodyX = 2.0_kreal*rotationRate*xyz(3)
end function

subroutine InitLogger(aLog)
	type(Logger), intent(inout) :: aLog
	call New(aLog,logLevel,logUnit)
	logInit = .TRUE.
end subroutine

end module
