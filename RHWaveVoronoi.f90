program VoronoiRHWave

use NumberKindsModule
use SphereGeomModule
use LoggerModule
use VoronoiPanelsModule
use VoronoiVelocitySerialModule

implicit none

! Grid variables
type(VorPanels), pointer :: panels=>null()
integer(kint) :: initNest, AMR, nTracer
integer(kint) :: j
real(kreal), allocatable :: totalKineticEnergy(:), totalEnstrophy(:),&
							 l1(:), l2(:), linf(:), relVortError(:)
integer(kint) :: remeshInterval
real(kreal) :: l1Denom, l2Denom, linfDenom

! Logger varialbes
type(Logger) :: exeLog
integer(kint) :: logOut
character(len=28) :: logKey
character(len=128) :: logString
real(kreal) :: programStart, programFinish

! Timestepping variables
real(kreal) :: t, dt, tfinal
integer(kint) :: timeJ, timesteps

! User input varialbes
character(len=128) :: jobPrefix
character(len=256) :: delTriVTKFile, voronoiVTKFile, dataFile
character(len=256) :: delTriVTKRoot, voronoiVTKroot
integer(kint), parameter :: readUnit = 12, writeUnit = 13
integer(kint) :: readStat, writeStat
namelist /gridInit/ initNest, AMR, remeshInterval
namelist /fileIO/ jobPrefix
namelist /time/ dt, tfinal
!
! Initialize computing environment
!
logOut = 6
call New(exeLog,DEBUG_LOGGING_LEVEL,logOut)
logKey = 'EXE_Log :'
call LogMessage(exeLog,TRACE_LOGGING_LEVEL,logKey,"*** Program START ***")

! Read user input from namelist file
open(unit=readUnit,file='RHWaveVoronoi.namelist',action='READ',status='OLD',iostat=readStat)
	if ( readStat /= 0) then
		call LogMessage(exeLog,ERROR_LOGGING_LEVEL,logKey,"ERROR : Unable to open namelist file.")
		stop
	endif
	read(readUnit,nml=gridInit)
	rewind(readUnit)
	read(readUnit,nml=time)
	rewind(readUnit)
	read(readUnit,nml=fileIO)
close(readUnit)	
! Init time stepping variables
t = 0.0_kreal
timeJ = 0
timesteps = floor(tfinal/dt)
if ( AMR == 0 ) then
	call UseSavedArea()
endif

! Prepare output files
write(dataFile,'(A,A,I1,A,F3.1,A,F5.3,A)') trim(jobPrefix),'_nest',initNest,'_rev',tfinal,'_dt',dt,'.dat'
write(delTriVTKRoot,'(A,A,A,I1,A,F3.1,A,F5.3,A)') 'vtkOut/',trim(jobPrefix),'delTri_nest',initNest,'_rev',tfinal,'_dt',dt,'_'
write(voronoiVTKRoot,'(A,A,A,I1,A,F3.1,A,F5.3,A)') 'vtkOut/',trim(jobPrefix),'voronoi_nest',initNest,'_rev',tfinal,'_dt',dt,'_'
write(delTriVTKFile,'(A,I0.4,A)') trim(delTriVTKRoot),timeJ,'.vtk'
write(voronoiVTKFile,'(A,I0.4,A)') trim(voronoiVTKRoot),timeJ,'.vtk'
call cpu_time(programStart)


allocate(totalKineticEnergy(0:timesteps))
totalKineticEnergy = 0.0_kreal
allocate(totalEnstrophy(0:timesteps))
totalEnstrophy = 0.0_kreal
allocate(l1(0:timesteps))
l1 = 0.0_kreal
allocate(l2(0:timesteps))
l2 = 0.0_kreal
allocate(linf(0:timesteps))
linf = 0.0_kreal


! Initialize grid
nTracer = 3
call New(panels,initNest,AMR,nTracer)
call InitRH4Wave(panels)
! Set initial latitude as tracer 1
do j=1,panels%N
	panels%tracer(j,1) = Latitude(panels%x(j),panels%y(j),panels%z(j))
enddo
totalEnstrophy(0) = 0.5_kreal*sum(panels%relVort(1:panels%N)*panels%relVort(1:panels%N)*panels%area(1:panels%N))
call VTKOutputDelaunayTriangulation(panels,delTriVTKFile)
call VTKOutputVoronoiGrid(panels,voronoiVTKFile)
linfDenom = maxval(abs(panels%relvort(1:panels%N)))
l2Denom = 2.0_kreal*totalEnstrophy(0)
l1Denom = sum(abs(panels%relVort(1:panels%N))*panels%area(1:panels%N))
!
!	Run the problem
!
do timeJ=0,timesteps-1
	if ( timeJ == 15 ) then
		print *,"breakpoint"
	endif
	! advance time
	call BVERK4(Panels,dt)
	t = real(timeJ+1,kreal)*dt
	! Store kinetic energy as tracer 2
	panels%tracer(1:panels%N,2) = kineticEnergy(1:panels%N)
	totalEnstrophy(timeJ+1) = 0.5_kreal*sum(panels%relVort(1:panels%N)*panels%relVort(1:panels%N)*panels%area(1:panels%N))
	totalKineticEnergy(timeJ+1) = totalKE
	
	! Calculate error with tracer 3
	do j=1,panels%N
		panels%tracer(j,3) = abs(panels%relVort(j) - &
				HaurwitzStationary4Relvort(panels%x(j),panels%y(j),panels%z(j)))
	enddo
	l2(timeJ+1) = sum(panels%tracer(1:panels%N,3)*panels%tracer(1:panels%N,3)*panels%area(1:panels%N))/l2denom
	l1(timeJ+1) = sum(panels%tracer(1:panels%N,3)*panels%area(1:panels%N))/l1denom
	linf(timeJ+1) = maxval(panels%tracer(1:panels%N,3))/linfDenom
	panels%tracer(1:panels%N,3) = panels%tracer(1:panels%N,3)/linfDenom
	! Write output
	call LogMessage(exeLog,DEBUG_LOGGING_LEVEL,trim(logKey)//" t = ",t)
	
	
	write(delTriVTKFile,'(A,I0.4,A)') trim(delTriVTKRoot),timeJ+1,'.vtk'
	write(voronoiVTKFile,'(A,I0.4,A)') trim(voronoiVTKRoot),timeJ+1,'.vtk'
	call VTKOutputDelaunayTriangulation(panels,delTriVTKFile)
	call VTKOutputVoronoiGrid(panels,voronoiVTKFile)
enddo

open(unit=writeUnit,file=dataFile,action='WRITE',status='REPLACE',iostat=writeStat)
if ( writeStat /= 0 ) call LogMessage(exeLog,ERROR_LOGGING_LEVEL,logkey,' ERROR opening dataFile.')
	write(writeUnit,'(5A24)') 'totalKE','totalEnstrophy','relVort l1','relVort l2','relVort linf'
	do j=0,timesteps
		write(writeUnit,'(5F24.15)') totalKineticEnergy(j), totalEnstrophy(j), &
			l1(j), l2(j), linf(j)
	enddo
close(writeUnit)

call cpu_time(programFinish)

call LogMessage(exeLog,TRACE_LOGGING_LEVEL,logKey,'Tracer1 = initial latitude.')
call LogMessage(exeLog,TRACE_LOGGING_LEVEL,logKey,'Tracer2 = kinetic energy.')
call LogMessage(exeLog,TRACE_LOGGING_LEVEL,logKey,'Tracer3 = vorticity error.')

write(logString,'(A,F9.2,A)') " elapsed time = ",(programFinish-programStart)/60.0_kreal," minutes."
call LogMessage(exeLog,TRACE_LOGGING_LEVEL,logKey,trim(logString))


call Delete(panels)
call ResetRK4()
deallocate(totalKineticEnergy)
deallocate(totalEnstrophy)
deallocate(l1)
deallocate(l2)
deallocate(linf)
call Delete(exeLog)

end program
