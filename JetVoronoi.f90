program JetVoronoi

use NumberKindsModule
use SphereGeomModule
use LoggerModule
use VoronoiPanelsModule
use VoronoiVelocityParallelModule
use VoronoiRemeshModule
use RefineRemeshModule
use TracerAndVorticityDistributionModule

implicit none

include 'mpif.h'

! Grid variables
type(VorPanels), pointer :: panels=>null()
integer(kint) :: initNest, AMR, nTracer
integer(kint) :: j
real(kreal), allocatable :: totalKineticEnergy(:), totalEnstrophy(:)
real(kreal) :: circTol, varTol, amrTol1, amrTol2
real(kreal) :: lat0, beta, perturbAmp
integer(kint) :: perturbWaveNum

! Machine & Logger variables
type(Logger) :: exeLog
integer(kint) :: logOut
character(len=28) :: logKey
character(len=128) :: logString
real(kreal) :: wTime, etime
integer(kint), parameter :: realBufferSize = 7, intBufferSize = 4
integer(kint) :: procRank, numProcs, mpiErrCode, intBuffer(intBufferSize)
real(kreal) :: realBuffer(realBufferSize)
real(kreal) :: longitudeBin, maxCirc, baseVar, minx0(3), maxx0(3)

! Timestepping variables
real(kreal) :: t, dt, tfinal
integer(kint) :: timeJ, timesteps

! Remeshing variables
integer(kint) :: remeshInterval
logical(klog) :: remeshFlag
integer(kint), parameter :: problemID = JET

! User input varialbes
character(len=128) :: jobPrefix
character(len=256) :: delTriVTKFile, voronoiVTKFile, dataFile
character(len=256) :: delTriVTKRoot, voronoiVTKroot
integer(kint), parameter :: readUnit = 12, writeUnit = 13
integer(kint) :: readStat, writeStat
namelist /gridInit/ initNest, AMR, remeshInterval, amrTol1, amrTol2
namelist /fileIO/ jobPrefix
namelist /time/ dt, tfinal
namelist /jetInit/ lat0, beta, perturbAmp, perturbWaveNum

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!	Part 1 : Initialize the computing environment / get user input  !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

call MPI_INIT(mpiErrCode)
call MPI_COMM_SIZE(MPI_COMM_WORLD,numProcs,mpiErrCode)
call MPI_COMM_RANK(MPI_COMM_WORLD,procRank,mpiErrCode)

logOut = 6
call New(exeLog,DEBUG_LOGGING_LEVEL,logOut)
write(logKey,'(A,I0.2,A)') 'EXE_LOG_',procRank,' : '
if ( procRank == 0 ) then

	call LogMessage(exeLog,TRACE_LOGGING_LEVEL,logKey,"*** Program START ***")
	
	! Read user input from namelist file
	open(unit=readUnit,file='JetVoronoi.namelist',action='READ',status='OLD',iostat=readStat)
		if ( readStat /= 0) then
			call LogMessage(exeLog,ERROR_LOGGING_LEVEL,logKey,"ERROR : Unable to open namelist file.")
			stop
		endif
		read(readUnit,nml=gridInit)
		rewind(readUnit)
		read(readUnit,nml=time)
		rewind(readUnit)
		read(readUnit,nml=fileIO)
		rewind(readUnit)
		read(readUnit,nml=jetInit)
	close(readUnit)	
	
	! Prepare output files
	write(dataFile,'(A,A,I1,A,F3.1,A,F5.3,A)') trim(jobPrefix),'_nest',initNest,'_rev',tfinal,'_dt',dt,'.dat'
	write(delTriVTKRoot,'(A,A,A,I1,A,F3.1,A,F5.3,A)') 'vtkOut/',trim(jobPrefix),'delTri_nest',initNest,'_rev',tfinal,'_dt',dt,'_'
	write(voronoiVTKRoot,'(A,A,A,I1,A,F3.1,A,F5.3,A)') 'vtkOut/',trim(jobPrefix),'voronoi_nest',initNest,'_rev',tfinal,'_dt',dt,'_'
	write(delTriVTKFile,'(A,I0.4,A)') trim(delTriVTKRoot),timeJ,'.vtk'
	write(voronoiVTKFile,'(A,I0.4,A)') trim(voronoiVTKRoot),timeJ,'.vtk'
	
	realBuffer(1) = dt
	realBuffer(2) = tfinal
	realBuffer(3) = circTol
	realBuffer(4) = varTol
	realBuffer(5) = lat0
	realBuffer(6) = beta
	realBuffer(7) = perturbAmp
	
	intBuffer(1) = initNest
	intBuffer(2) = AMR
	intBuffer(3) = remeshInterval
	intBuffer(4) = perturbWaveNum
	
endif
call MPI_BCAST(realBuffer,realBufferSize,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiErrCode)
call MPI_BCAST(intBuffer,intBufferSize,MPI_INTEGER,0,MPI_COMM_WORLD,mpiErrCode)

dt = realBuffer(1)
tfinal = realBuffer(2)
circTol = realBuffer(3)
varTol = realBuffer(4)
lat0 = realBuffer(5)
beta = realBuffer(6)
perturbAmp = realBuffer(7)

initNest = intBuffer(1)
AMR = intBuffer(2)
remeshInterval = intBuffer(3)
perturbWaveNum = intBuffer(4)
remeshFlag = .False.

if (procRank == 0) wTime = MPI_WTIME()

! Init time stepping variables
t = 0.0_kreal
timeJ = 0
timesteps = floor(tfinal/dt)


allocate(totalKineticEnergy(0:timesteps))
totalKineticEnergy = 0.0_kreal
allocate(totalEnstrophy(0:timesteps))
totalEnstrophy = 0.0_kreal

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!	Part 2 : Initialize the grid		             !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

nTracer = 3
call New(panels,initNest,AMR,nTracer)

call InitJet(panels,lat0,beta,perturbAmp,perturbWaveNum)

if ( AMR == 0 ) then
	call UseSavedArea()
else
	maxCirc = maxval(abs(panels%relVort(1:panels%N))*panels%area(1:panels%N))
	varTol = amrTol2 * sqrt(4.0_kreal*PI/panels%N)
	circTol = amrTol1 * maxCirc
	call Delete(panels)
	call New(panels,initNest,AMR,nTracer,circTol, varTol,problemID,lat0,beta,perturbAmp,perturbWaveNum)
endif

! Set initial latitude as tracer 1
do j=1,panels%N
	panels%tracer(j,1) = Latitude(panels%x(j),panels%y(j),panels%z(j))
enddo
! Use tracer 2 for kinetic energy
! Set initial Longitude as tracer 3
do j=1,panels%N
	longitudeBin = mod( Longitude(panels%x0(j),panels%y0(j),panels%z0(j))*180.0_kreal/PI, 90.0_kreal )
	if ( longitudeBin <= 30.0_kreal ) then
		panels%tracer(j,3) = -1.0_kreal
	elseif (longitudeBin > 60.0_kreal) then
		panels%tracer(j,3) = 1.0_kreal
	else
		panels%tracer(j,3) = 0.0_kreal
	endif
enddo

totalEnstrophy(0) = 0.5_kreal*sum(panels%relVort(1:panels%N)*panels%relVort(1:panels%N)*panels%area(1:panels%N))
if ( procRank == 0) then
	call StartSection(exeLog%writer,logKey,' Grid Initialized :')
	call LogMessage(exeLog,TRACE_LOGGING_LEVEL,' n Faces = ',panels%N)
	call LogMessage(exeLog,TRACE_LOGGING_LEVEL,' surf. area = ',sum(panels%area))
	call VTKOutputDelaunayTriangulation(panels,delTriVTKFile)
	call VTKOutputVoronoiGrid(panels,voronoiVTKFile)
endif

call InitializeMPIRK4(panels%N,procRank,numProcs)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!	Part 3 : Run the problem								    !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
do timeJ=0,timesteps-1

	! Remesh if necessary
	if ( mod(timeJ+1,remeshInterval) == 0 ) then
		remeshFlag = .TRUE.
		if ( procRank == 0 )&
			call LogMessage(exeLog,TRACE_LOGGING_LEVEL,logKey,"Remesh triggered by remeshInterval.")		
	endif
	if ( remeshFlag) then
		call RemeshVoronoi(panels,initNest,AMR,circTol,varTol,procRank,problemID,lat0,beta,perturbAmp,perturbWaveNum)
		remeshFlag = .False.
		do j=1,panels%N
			panels%tracer(j,1) = Latitude(panels%x0(j),panels%y0(j),panels%z0(j))
			longitudeBin = mod( Longitude(panels%x0(j),panels%y0(j),panels%z0(j))*180.0_kreal/PI, 90.0_kreal )
			if ( longitudeBin <= 30.0_kreal ) then
				panels%tracer(j,3) = -1.0_kreal
			elseif (longitudeBin > 60.0_kreal) then
				panels%tracer(j,3) = 1.0_kreal
			else
				panels%tracer(j,3) = 0.0_kreal
			endif
		enddo
	endif
	
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
		
		write(delTriVTKFile,'(A,I0.4,A)') trim(delTriVTKRoot),timeJ+1,'.vtk'
		write(voronoiVTKFile,'(A,I0.4,A)') trim(voronoiVTKRoot),timeJ+1,'.vtk'
		call VTKOutputDelaunayTriangulation(panels,delTriVTKFile)
		call VTKOutputVoronoiGrid(panels,voronoiVTKFile)
	endif
enddo

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!	Part 4 : Finish and clear memory							!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

call FinalizeMPIRK4()
call Delete(panels)
deallocate(totalKineticEnergy)
deallocate(totalEnstrophy)
call Delete(exeLog)

call MPI_FINALIZE(mpiErrCode)
end program