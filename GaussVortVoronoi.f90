program GaussVortVoronoiMPI

use NumberKindsModule
use SphereGeomModule
use OutputWriterModule
use LoggerModule
use VoronoiPanelsModule
use VoronoiVelocityParallelModule
use TracerAndVorticityDistributionModule

implicit none

include 'mpif.h'

! Grid variables
type(VorPanels), pointer :: panels=>null()
integer(kint) :: initNest, AMR, nTracer
integer(kint) :: j
real(kreal), allocatable :: totalKineticEnergy(:), totalEnstrophy(:)

! Machine & Logger variables
type(Logger) :: exeLog
integer(kint) :: logOut
character(len=28) :: logKey
character(len=128) :: logString
real(kreal) :: wTime, etime
integer(kint), parameter :: realBufferSize = 2, intBufferSize = 2
integer(kint) :: procRank, numProcs, mpiErrCode, intBuffer(intBufferSize)
real(kreal) :: realBuffer(realBufferSize)

! Timestepping variables
real(kreal) :: t, dt, tfinal
integer(kint) :: timeJ, timesteps

! User input varialbes
character(len=128) :: jobPrefix
character(len=256) :: delTriVTKFile, voronoiVTKFile, dataFile
character(len=256) :: delTriVTKRoot, voronoiVTKroot
integer(kint), parameter :: readUnit = 12, writeUnit = 13
integer(kint) :: readStat, writeStat
namelist /gridInit/ initNest, AMR
namelist /fileIO/ jobPrefix
namelist /time/ dt, tfinal

!
! Initialize computing environment
!

call MPI_INIT(mpiErrCode)
call MPI_COMM_SIZE(MPI_COMM_WORLD,numProcs,mpiErrCode)
call MPI_COMM_RANK(MPI_COMM_WORLD,procRank,mpiErrCode)

logOut = 6
call New(exeLog,DEBUG_LOGGING_LEVEL,logOut)
write(logKey,'(A,I0.2,A)') 'EXE_LOG_',procRank,' : '

if ( procRank == 0 ) then

	call LogMessage(exeLog,TRACE_LOGGING_LEVEL,logKey,"*** Program START ***")

	! Read user input from namelist file
	open(unit=readUnit,file='GaussVortVoronoi.namelist',action='READ',status='OLD',iostat=readStat)
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
	
	! Prepare output files
	write(dataFile,'(A,A,I1,A,F3.1,A,F5.3,A)') trim(jobPrefix),'_nest',initNest,'_rev',tfinal,'_dt',dt,'.dat'
	write(delTriVTKRoot,'(A,A,A,I1,A,F3.1,A,F5.3,A)') 'vtkOut/',trim(jobPrefix),'delTri_nest',initNest,'_rev',tfinal,'_dt',dt,'_'
	write(voronoiVTKRoot,'(A,A,A,I1,A,F3.1,A,F5.3,A)') 'vtkOut/',trim(jobPrefix),'voronoi_nest',initNest,'_rev',tfinal,'_dt',dt,'_'
	write(delTriVTKFile,'(A,I0.4,A)') trim(delTriVTKRoot),timeJ,'.vtk'
	write(voronoiVTKFile,'(A,I0.4,A)') trim(voronoiVTKRoot),timeJ,'.vtk'
	
	realBuffer(1) = dt
	realBuffer(2) = tfinal
	
	intBuffer(1) = initNest
	intBuffer(2) = AMR
endif
! Send user input data to all processes
call MPI_BCAST(realBuffer,realBufferSize,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiErrCode)
call MPI_BCAST(intBuffer,intBufferSize,MPI_INTEGER,0,MPI_COMM_WORLD,mpiErrCode)

dt = realBuffer(1)
tfinal = realBuffer(2)

initNest = intBuffer(1)
AMR = intBuffer(2)


! Init time stepping variables
t = 0.0_kreal
timeJ = 0
timesteps = floor(tfinal/dt)
if ( AMR == 0 ) then
	call UseSavedArea()
endif


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!	Part 2 : Initialize the grid		             !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

allocate(totalKineticEnergy(0:timesteps))
totalKineticEnergy = 0.0_kreal
allocate(totalEnstrophy(0:timesteps))
totalEnstrophy = 0.0_kreal

if (procRank == 0) wTime = MPI_WTIME()

nTracer = 2
call New(panels,initNest,AMR,nTracer)

call InitGaussianVortex(panels)
! Set initial latitude as tracer 1
do j=1,panels%N
	panels%tracer(j,1) = Latitude(panels%x(j),panels%y(j),panels%z(j))
enddo
totalEnstrophy(0) = 0.5_kreal*sum(panels%relVort(1:panels%N)*panels%relVort(1:panels%N)*panels%area(1:panels%N))
if ( procRank == 0) then
	call StartSection(exeLog%writer,logKey,' Grid Initialized :')
	call LogMessage(exeLog,TRACE_LOGGING_LEVEL,' n Faces = ',panels%N)
	call LogMessage(exeLog,TRACE_LOGGING_LEVEL,' surf. area = ',sum(panels%area))
	call LogMessage(exeLog,TRACE_LOGGING_LEVEL,' rel. vort. integral = ',sum(panels%relVort*panels%area))
	call VTKOutputDelaunayTriangulation(panels,delTriVTKFile)
	call VTKOutputVoronoiGrid(panels,voronoiVTKFile)
endif

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!	Part 3 : Run the problem								    !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
call MPI_BARRIER(MPI_COMM_WORLD,mpiErrCode)

call InitializeMPIRK4(panels%N,procRank,numProcs)

do timeJ=0,timesteps-1
	if ( timeJ == 0 ) etime = MPI_WTIME()
	! advance time
	call BVERK4(Panels,dt,procRank,numProcs)
	if ( timeJ == 1 .AND. procRank == 0 ) then
		etime = MPI_WTIME() - etime
		write(logString,'(A,F9.3,A)') 'Estimated time left = ',(timesteps-1)*etime/60.0,' minutes.'
		call LogMessage(exeLog,TRACE_LOGGING_LEVEL,logKey,logString)
	endif
	t = real(timeJ+1,kreal)*dt
	! Store kinetic energy as tracer 2
	panels%tracer(1:panels%N,2) = kineticEnergy(1:panels%N)
	totalEnstrophy(timeJ+1) = 0.5_kreal*sum(panels%relVort(1:panels%N)*panels%relVort(1:panels%N)*panels%area(1:panels%N))
	totalKineticEnergy(timeJ+1) = totalKE
	
	if ( procRank == 0) then
		! Write output
		call LogMessage(exeLog,DEBUG_LOGGING_LEVEL,trim(logKey)//" t = ",t)
		call LogMessage(exeLog,TRACE_LOGGING_LEVEL,' surf. area = ',sum(panels%area))
		call LogMessage(exeLog,TRACE_LOGGING_LEVEL,' rel. vort. integral = ',sum(panels%relVort*panels%area))
		write(delTriVTKFile,'(A,I0.4,A)') trim(delTriVTKRoot),timeJ+1,'.vtk'
		write(voronoiVTKFile,'(A,I0.4,A)') trim(voronoiVTKRoot),timeJ+1,'.vtk'
		call VTKOutputDelaunayTriangulation(panels,delTriVTKFile)
		call VTKOutputVoronoiGrid(panels,voronoiVTKFile)
	endif
enddo

if ( procRank == 0 ) then
	open(unit=writeUnit,file=dataFile,action='WRITE',status='REPLACE',iostat=writeStat)
	if ( writeStat /= 0 ) call LogMessage(exeLog,ERROR_LOGGING_LEVEL,logkey,' ERROR opening dataFile.')
		write(writeUnit,'(2A24)') 'totalKE','totalEnstrophy'
		do j=0,timesteps
			write(writeUnit,'(2F24.15)') totalKineticEnergy(j), totalEnstrophy(j)
		enddo
	close(writeUnit)
	wTime = MPI_WTIME() - wtime
	call LogMessage(exeLog,TRACE_LOGGING_LEVEL,logKey,"*** Program END ***")
	call LogMessage(exeLog,TRACE_LOGGING_LEVEL,trim(logKey)//' data file = ',trim(dataFile))
	write(delTriVTKFile,'(A,A)') trim(delTriVTKRoot),'XXXX.vtk'
	write(voronoiVTKFile,'(A,A)') trim(voronoiVTKRoot),'XXXX.vtk'
	call LogMessage(exeLog,TRACE_LOGGING_LEVEL,trim(logKey)//' Delaunay Tri. VTK Output = ',trim(delTriVTKFile))
	call LogMessage(exeLog,TRACE_LOGGING_LEVEL,trim(logKey)//' Voronoi grid VTK Output = ',trim(voronoiVTKFile))
	call LogMessage(exeLog,TRACE_LOGGING_LEVEL,logKey,'Tracer1 = initial latitude.')
	call LogMessage(exeLog,TRACE_LOGGING_LEVEL,logKey,'Tracer2 = kinetic energy.')
	call LogMessage(exeLog,TRACE_LOGGING_LEVEL,logKey,'Tracer3 = vorticity error.')
	call EndSection(exeLog%writer)
	write(logString,'(A,F9.2,A)') " elapsed time = ",wtime/60.0_kreal," minutes."
	call LogMessage(exeLog,TRACE_LOGGING_LEVEL,logKey,trim(logString))
endif

call FinalizeMPIRK4()
call Delete(panels)
deallocate(totalKineticEnergy)
deallocate(totalEnstrophy)
call Delete(exeLog)

call MPI_FINALIZE(mpiErrCode)	

end program