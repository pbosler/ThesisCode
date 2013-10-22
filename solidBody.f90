program SolidBody

use NumberKindsModule
use SphereGeomModule
use LoggerModule
use ParticlesEdgesPanelsModule
use PanelVelocityParallelModule
use TracerAndVorticityDistributionModule
use RefineRemeshModule

implicit none

include 'mpif.h'
!
! grid variables
!
type(Particles), pointer :: gridParticles=>null()
type(Edges), pointer :: gridEdges=>null()
type(Panels), pointer :: gridPanels=>null()

integer(kint) :: panelKind, initNest, AMR, nTracer, remeshInterval, reproject
integer(kint), parameter :: problemKind = BVE_SOLVER, problemID = SOLID_BODY
logical(klog) :: remeshFlag, renormalize

real(kreal) :: rotationRate

!
! Error calculation variables
!
real(kreal), allocatable :: particlesl2(:), particlesLinf(:), panelsL2(:), panelsLinf(:), maxVel(:), &
			totalKineticEnergy(:), totalEnstrophy(:), particlesExact(:,:), panelsExact(:,:),&
			particlesNorm(:), panelsNorm(:)
real(kreal) :: particleL2denom, particleLinfDenom, panelL2denom, panelLinfdenom

!
! Computation managment varialbes
!
type(logger) :: exeLog
integer(kint), parameter :: logOut = 6, REAL_INIT_SIZE = 3, INT_INIT_SIZE = 4
character(len=28) :: logKey
character(len=128) :: logString
real(kreal) :: wtime, etime, realBuffer(REAL_INIT_SIZE)
integer(kint) :: procRank, numProcs, mpiErrCode, intBuffer(INT_INIT_SIZE)
logical(klog) :: newEstimate

!
! timestepping variables
!
real(kreal) :: t, dt, tfinal
integer(kint) :: timeJ, timesteps

!
! user & I/O variables
!
character(len=128) :: jobPrefix
character(len=256) :: vtkRoot, vtkFile, dataFile
integer(kint), parameter :: readUnit=12, writeUnit=13
integer(kint) :: readStat, writeStat

namelist /gridInit/ panelKind, initNest, remeshInterval, rotationRate, reproject
namelist /time/ dt, tfinal
namelist /fileIO/ jobPrefix

!
! general variables
!
real(kreal), parameter :: zero = 0.0_kreal
integer(kint) :: j,k

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!	Part 1 : Initialize the computing environment / get user input  !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

call MPI_INIT(mpiErrCode)
call MPI_COMM_SIZE(MPI_COMM_WORLD,numProcs,mpiErrCode)
call MPI_COMM_RANK(MPI_COMM_WORLD,procRank,mpiErrCode)

call New(exeLog,DEBUG_LOGGING_LEVEL,logOut)
write(logKey,'(A,I0.2,A)') 'EXE_LOG_',procRank,' : '

if ( procRank == 0 ) then
	call LogMessage(exeLog,TRACE_LOGGING_LEVEL,logKey,"**** PROGRAM START ****")
	call LogMessage(exeLog,DEBUG_LOGGING_LEVEL,logKey//"numProcs = ",numProcs)
	call LogMessage(exeLog,DEBUG_LOGGING_LEVEL,logKey,"reading input file...")
	
	! read user input
	
	open(unit=readUnit,file='SolidBody.namelist',action='READ',status='OLD',iostat=readStat)
		if ( readStat /= 0 ) then
			call LogMessage(exeLog,ERROR_LOGGING_LEVEL,logKey,"ERROR opening namelist file.")
			stop
		endif
		read(readUnit,nml=gridInit)
		rewind(readUnit)
		read(readUnit,nml=time)
		rewind(readunit)
		read(readUnit,nml=fileIO)
		rewind(readUnit)
	close(readUnit)
	
	intBuffer(1) = panelKind
	intBuffer(2) = initNest
	intBuffer(3) = remeshInterval
	intBuffer(4) = reproject
	
	realBuffer(1) = dt
	realBuffer(2) = tfinal
	realBuffer(3) = rotationRate
	
	! prepare output files
	
	if ( panelKind == TRI_PANEL) then
		write(jobPrefix,'(A,A,I1,A,I0.2,A,F3.1,A,F7.5,A)') trim(jobPrefix),'_triNest',initNest,&
			'_LG',remeshInterval,'_rev',tfinal,'_dt',dt,'_'
	elseif (panelKind == QUAD_PANEL) then
		write(jobPrefix,'(A,A,I1,A,I0.2,A,F3.1,A,F7.5,A)') trim(jobPrefix),'_quadNest',initNest,&
			'_LG',remeshInterval,'_rev',tfinal,'_dt',dt,'_'
	endif
	write(dataFile,'(A,A)') trim(jobPrefix),'.dat'
	write(vtkRoot,'(A,A)') 'vtkOut/',trim(jobPrefix)
	write(vtkFile,'(A,I0.4,A)') trim(vtkRoot),0,'.vtk'
	
	call LogMessage(exeLog,DEBUG_LOGGING_LEVEL,logKey,'... done. Broadcasting input data...')
endif

call MPI_BCAST(intBuffer,INT_INIT_SIZE,MPI_INTEGER,0,MPI_COMM_WORLD,mpiErrCode)
call MPI_BCAST(realBuffer,REAL_INIT_SIZE,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiErrCode)
call LogMessage(exeLog,DEBUG_LOGGING_LEVEL,logKey,'... init broadcast complete.')

panelKind = intBuffer(1)
initNest = intBuffer(2)
remeshInterval = intBuffer(3)
reproject = intBuffer(4)

dt = realBuffer(1)
tfinal = realBuffer(2)
rotationRate = realBuffer(3)

remeshFlag = .False.
newEstimate = .FALSE.

AMR = 0
nTracer = 5
t = 0.0_kreal
timeJ = 0
timesteps = floor(tfinal/dt)

if ( reproject /= 0 ) then
	renormalize = .TRUE.
else
	renormalize = .FALSE.
endif

allocate(particlesL2(0:timesteps))
particlesL2 = 0.0_kreal
allocate(particlesLinf(0:timesteps))
particlesLinf = 0.0_kreal
allocate(panelsL2(0:timesteps))
panelsL2 = 0.0_kreal
allocate(panelsLinf(0:timesteps))
panelsLinf = 0.0_kreal
allocate(maxVel(0:timesteps))
maxVel = 0.0_kreal
allocate(totalKineticEnergy(0:timesteps))
totalKineticEnergy = 0.0_kreal
allocate(totalEnstrophy(0:timesteps))
totalEnstrophy = 0.0_kreal
allocate(particlesNorm(0:timesteps))
particlesNorm = 0.0_kreal
allocate(panelsNorm(0:timesteps))
panelsNorm = 0.0_kreal

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!	Part 2 : Initialize the grid		             !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

wtime = MPI_WTIME()

call New(gridParticles,gridEdges,gridPanels,panelKind,initNest,AMR,nTracer,problemKind)

call SetOmega(zero)

call InitWilliamsonCosineBell(gridParticles,gridPanels)

call InitSolidBody(gridParticles,gridPanels,rotationRate)

! store initial latitude in tracer 1
do j=1,gridParticles%N
	gridParticles%tracer(j,1) = Latitude(gridParticles%x0(:,j))
enddo
do j=1,gridPanels%N
	if ( .NOT. gridPanels%hasChildren(j) ) then
		gridPanels%tracer(j,1) = Latitude(gridPanels%x0(:,j))
	endif
enddo

totalEnstrophy(0) = 0.5_kreal*sum(gridPanels%area(1:gridPanels%N)* &
					gridPanels%relVort(1:gridPanels%N)*gridPanels%relVort(1:gridPanels%N))
particleL2denom = 0.0_kreal
do j=1,gridParticles%N
	particleL2denom = particleL2denom + sum(gridparticles%x(:,j)*gridparticles%x(:,j))*gridparticles%dualArea(j)
enddo
particleL2denom = sqrt(particleL2denom)
particleLinfDenom = 1.0_kreal

panelL2denom = 0.0_kreal
do j=1,gridPanels%N
	panelL2denom = panelL2denom + sum(gridpanels%x(:,j)*gridpanels%x(:,j))*gridpanels%area(j)
enddo
panelL2denom = sqrt(panelL2denom)
panelLinfdenom = 1.0_kreal

allocate(particlesExact(3,gridParticles%N))
do j=1,gridParticles%N
	particlesExact(:,j) = gridParticles%x(:,j)
enddo
allocate(panelsExact(3,gridPanels%N))
do j=1,gridPanels%N
	if ( gridPanels%hasChildren(j) ) then
		panelsExact(:,j) = 0.0_kreal
	else	
		panelsExact(:,j) = gridPanels%x(:,j)
	endif
enddo	

if ( procRank == 0 ) then
	call PrintStats(gridParticles)
	call PrintStats(gridEdges)
	call PrintStats(gridPanels)
	call vtkOutput(gridParticles,gridPanels,vtkFile)
endif

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!	Part 3 : Run the problem								    !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

call LogMessage(exeLog,DEBUG_LOGGING_LEVEL,logKey,'Setup complete. Starting time integration ...')

call InitializeMPIRK4(gridParticles,gridPanels,procRank,numProcs)

! MAIN Timestepping loop
do timeJ = 0,timesteps - 1
	! remesh if necessary
	if ( mod(timeJ+1,remeshInterval) == 0 ) then
		remeshFlag = .TRUE.
		newEstimate = .TRUE.
		if (procRank == 0 ) call LogMessage(exeLog,TRACE_LOGGING_LEVEL,logKey,"Remesh triggered by remeshInterval.")
	endif
	if ( remeshFlag) then
		call AdaptiveRemesh(gridParticles,gridEdges,gridPanels,&
							initNest,AMR, zero,zero, &
							procRank, problemID, &
							rotationRate=rotationRate)
		remeshFlag = .FALSE.
		
		call ResetRK4()
		
		do j=1,gridParticles%N
			gridParticles%tracer(j,1) = Latitude(gridParticles%x0(:,j))
		enddo
		do j=1,gridPanels%N
			if (.NOT. gridPanels%hasChildren(j) ) then
				gridPanels%tracer(j,1) = Latitude(gridPanels%x0(:,j))
			endif
		enddo
		if ( procRank == 0 ) then
			call PrintStats(gridParticles)
			call PrintStats(gridEdges)
			call PrintStats(gridPanels)
		endif
	endif
	
	! advance time
	call BVERK4(gridParticles,gridPanels,dt,procRank,numProcs)
	t = real(timeJ+1,kreal)*dt
	
	if ( renormalize) then 
		call RenormalizeGrid(gridParticles,gridPanels)
	else
		do j=1,gridParticles%N
			gridParticles%tracer(j,5) =sqrt(sum(gridParticles%x(:,j)*gridParticles%x(:,j)))
		enddo
		do j=1,gridPanels%N
			gridPanels%tracer(j,5) = sqrt(sum(gridPanels%x(:,j)*gridPanels%x(:,j)))
		enddo
		particlesNorm(timeJ+1) = maxval(abs(1.0_kreal - gridParticles%tracer(:,5)))
		panelsNorm(timeJ+1) = maxval(abs(1.0_kreal- gridPanels%tracer(:,5)))
	endif
	
	if ( timeJ==0) totalKineticEnergy(0) = totalKE
	totalKineticEnergy(timeJ+1) = totalKE
	maxVel(timeJ+1) = maxval(sqrt(gridPanels%tracer(1:gridPanels%N,2)))
	totalEnstrophy(timeJ+1) = 0.5_kreal*sum(gridPanels%relVort(1:gridPanels%N) * &
											gridPanels%relVort(1:gridPanels%N) * &
											gridPanels%area(1:gridPanels%N) )

	! Calculate 'exact' positions
	do j=1,gridParticles%N
		particlesExact(1,j) = gridParticles%x0(1,j)*cos(rotationRate*t)-gridParticles%x0(2,j)*sin(rotationRate*t)
		particlesExact(2,j) = gridParticles%x0(2,j)*cos(rotationRate*t)+gridParticles%x0(1,j)*sin(rotationRate*t)
		particlesExact(3,j) = gridParticles%x0(3,j)
	enddo
	do j=1,gridPanels%N
		if (.NOT. gridPanels%hasChildren(j) ) then
			panelsExact(1,j) = gridpanels%x0(1,j)*cos(rotationRate*t)-gridPanels%x0(2,j)*sin(rotationRate*t)
			panelsExact(2,j) = gridPanels%x0(2,j)*cos(rotationRate*t)+gridPanels%x0(1,j)*sin(rotationRate*t)
			panelsExact(3,j) = gridPanels%x0(3,j)
		endif
	enddo
	
	! Store position error in tracer 3
	do j=1,gridParticles%N
		gridParticles%tracer(j,3) = sqrt(sum((gridParticles%x(:,j)-particlesExact(:,j))*(gridParticles%x(:,j)-particlesExact(:,j))))
	enddo
	do j=1,gridPanels%N
		if (.NOT. gridPanels%hasChildren(j) ) then
			gridPanels%tracer(j,3) = sqrt(sum((gridPanels%x(:,j)-panelsExact(:,j))*(gridPanels%x(:,j)-panelsExact(:,j))))
		endif
	enddo
	
	! find error norms
	particlesL2(timeJ+1) = sum( gridParticles%tracer(:,3)*gridParticles%tracer(:,3)*gridParticles%dualArea)
	particlesL2(timeJ+1) = sqrt(particlesL2(timeJ+1))/particleL2denom
	panelsL2(timeJ+1) = sum(gridPanels%tracer(:,3)*gridPanels%tracer(:,3)*gridPanels%area)
	panelsL2(timeJ+1) = sqrt(panelsL2(timeJ+1))/panelL2denom
	particlesLinf(timeJ+1) = maxVal(gridParticles%tracer(:,3))
	panelsLinf(timeJ+1) = maxVal(gridPanels%tracer(:,3))
	
	
	if (newEstimate .AND. procRank == 0) then
		etime = MPI_WTIME() - wtime
		write(logString,'(A,F9.4,A)') 'Estimated time left = ', (timesteps-1-timeJ)/timeJ*etime/60.0_kreal, ' minutes.'
		call LogMessage(exeLog,TRACE_LOGGING_LEVEL,logKey,trim(logString))
		newEstimate = .FALSE.
	endif
	
	if (procRank== 0) then
		call LogMessage(exeLog,DEBUG_LOGGING_LEVEL,trim(logKey)//" t = ",t)
		write(vtkFile,'(A,I0.4,A)') trim(vtkRoot),timeJ+1,'.vtk'
		call vtkOutput(gridParticles,gridPanels,vtkFile)
	endif
	
enddo

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!	Part 4 : Finish and clear memory							!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

if (procRank == 0 ) then
	open(unit=writeUnit,file=dataFile,status='REPLACE',action='WRITE',iostat=writeStat)
	if (writeStat/=0) then
		call LogMessage(exeLog,ERROR_LOGGING_LEVEL,logKey,'ERROR opening datafile')
	else
		write(writeUnit,'(9A24)') 		'particlesL2',  'particlesLinf',  'panelsL2',  'panelsLinf',&
			'totalKE',             'totalEns',       'maxVel',      'particlesDist', 'panelsDist'
		do j=0,timesteps
			write(writeUnit,'(9F24.15)') particlesL2(j), particlesLinf(j), panelsL2(j), panelsLinf(j), &
			totalKineticEnergy(j), totalEnstrophy(j), maxVel(j), particlesNorm(j), panelsNorm(j)
		enddo
	endif
	close(writeUnit)
	
	wtime = MPI_WTIME()-wtime
	call LogMessage(exeLog,TRACE_LOGGING_LEVEL,logKey,'****** program end ****** ')
	call StartSection(exeLog,'-------- Calculation Summary -------------')
	call LogMessage(exeLog,TRACE_LOGGING_LEVEL,logKey//' dataFile = ',trim(dataFile))
	write(vtkFile,'(A,A)') trim(vtkRoot),'XXXX.vtk'
	call LogMessage(exeLog,TRACE_LOGGING_LEVEL,logKey//' vtkFiles = ',trim(vtkFile))
	call LogMessage(exeLog,TRACE_LOGGING_LEVEL,logKey,'tracer1 = initial latitude')
	call LogMessage(exeLog,TRACE_LOGGING_LEVEL,logKey,'tracer2 = kinetic energy')
	call LogMessage(exeLog,TRACE_LOGGING_LEVEL,logKey,'tracer3 = position error')
	call LogMessage(exeLog,TRACE_LOGGING_LEVEL,logKey//'panelKind = ',panelKind)
	call LogMessage(exeLog,TRACE_LOGGING_LEVEL,logKey//'remeshInterval = ',remeshInterval)
	call LogMessage(exeLog,TRACE_LOGGING_LEVEL,logKey//'dt = ',dt)
	call EndSection(exeLog)
	write(logString,'(A,F9.2,A)') " elapsed time = ",wtime/60.0_kreal," minutes."
	call LogMessage(exeLog,TRACE_LOGGING_LEVEL,logKey,trim(logString))
	
	
endif



deallocate(particlesL2)
deallocate(particlesLinf)
deallocate(panelsL2)
deallocate(panelsLinf)
deallocate(maxVel)
deallocate(particlesExact)
deallocate(panelsExact)
deallocate(totalKineticEnergy)
deallocate(totalEnstrophy)
call FinalizeMPIRK4()
call Delete(gridParticles,gridEdges,gridPanels)
call Delete(exeLog)

call MPI_FINALIZE(mpiErrCode)

end program
