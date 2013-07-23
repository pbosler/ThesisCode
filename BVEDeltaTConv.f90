program BVEDeltaTConv

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
type(Particles), pointer :: gridParticles=>null(), refParticles=>null()
type(Edges), pointer :: gridEdges=>null(), refEdges=>null()
type(Panels), pointer :: gridPanels=>null(), refPanels=>null()

integer(kint) :: panelKind, initNest, AMR, nTracer, remeshInterval, reproject
integer(kint), parameter :: problemKind = BVE_SOLVER, problemID = SOLID_BODY
logical(klog) :: remeshFlag, renormalize

!
! Error calculation variables
!
real(kreal), allocatable :: particlesl2(:), particlesLinf(:), panelsL2(:), panelsLinf(:)
real(kreal) :: particleL2denom, particleLinfDenom, panelL2denom, panelLinfdenom

!
! Computation managment varialbes
!
type(logger) :: exeLog
integer(kint), parameter :: logOut = 6, REAL_INIT_SIZE = 2, INT_INIT_SIZE = 3
character(len=28) :: logKey
character(len=128) :: logString
real(kreal) :: wtime, etime, realBuffer(REAL_INIT_SIZE)
integer(kint) :: procRank, numProcs, mpiErrCode, intBuffer(INT_INIT_SIZE)
logical(klog) :: newEstimate

!
! timestepping variables
!
real(kreal) :: t, dt, dtt, tfinal, refDt
integer(kint) :: timeJ, timesteps, dtj
integer(kint), parameter :: DT_NUMBER = 4
!
! user & I/O variables
!
character(len=128) :: jobPrefix
character(len=256) :: vtkRoot, vtkFile, dataFile
integer(kint), parameter :: readUnit=12, writeUnit=13
integer(kint) :: readStat, writeStat

namelist /gridInit/ panelKind, initNest, reproject
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
	
	open(unit=readUnit,file='BVEDTConv.namelist',action='READ',status='OLD',iostat=readStat)
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
	intBuffer(3) = reproject
	
	realBuffer(1) = dt
	realBuffer(2) = tfinal
	
	! prepare output files
	
	if ( panelKind == TRI_PANEL) then
		write(jobPrefix,'(A,A,I1,A,I0.2,A,F3.1,A,F7.5,A)') trim(jobPrefix),'_triNest',initNest,&
			'_LG',remeshInterval,'_rev',tfinal,'_dt',dt,'_'
	elseif (panelKind == QUAD_PANEL) then
		write(jobPrefix,'(A,A,I1,A,I0.2,A,F3.1,A,F7.5,A)') trim(jobPrefix),'_quadNest',initNest,&
			'_LG',remeshInterval,'_rev',tfinal,'_dt',dt,'_'
	endif
	write(dataFile,'(A,A)') trim(jobPrefix),'.dat'

	call LogMessage(exeLog,DEBUG_LOGGING_LEVEL,logKey,'... done. Broadcasting input data...')
endif

call MPI_BCAST(intBuffer,INT_INIT_SIZE,MPI_INTEGER,0,MPI_COMM_WORLD,mpiErrCode)
call MPI_BCAST(realBuffer,REAL_INIT_SIZE,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiErrCode)
call LogMessage(exeLog,DEBUG_LOGGING_LEVEL,logKey,'... init broadcast complete.')

panelKind = intBuffer(1)
initNest = intBuffer(2)
reproject = intBuffer(3)

dt = realBuffer(1)
tfinal = realBuffer(2)


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

allocate(particlesL2(DT_NUMBER))
particlesL2 = 0.0_kreal
allocate(particlesLinf(DT_NUMBER))
particlesLinf = 0.0_kreal
allocate(panelsL2(DT_NUMBER))
panelsL2 = 0.0_kreal
allocate(panelsLinf(DT_NUMBER))
panelsLinf = 0.0_kreal

wtime = MPI_WTIME()

! Compute reference solution
call New(refParticles,refEdges,refPanels,panelKind,initNest,AMR,nTracer,problemKind)

call SetOmega(zero)

call InitSolidBody(refParticles,refPanels,2.0_kreal*PI)

particleL2denom = 0.0_kreal
do j=1,refParticles%N
	particleL2denom = particleL2denom + sum(refParticles%x(:,j)*refParticles%x(:,j))*refParticles%dualArea(j)
enddo
panelL2denom = 0.0_kreal
do j=1,refPanels%N
	panelL2denom = panelL2denom + sum(refPanels%x(:,j)*refPanels%x(:,j))*refPanels%area(j)
enddo


if ( procRank == 0 ) then
	call PrintStats(refParticles)
	call PrintStats(refEdges)
	call PrintStats(refPanels)
endif


print *, " DEBUG 0"
!call MPI_BARRIER(MPI_COMM_WORLD)
print *, " DEBUG 1"
call LogMessage(exeLog,DEBUG_LOGGING_LEVEL,logKey,'Setup complete. Starting time integration ...')

call InitializeMPIRK4(refParticles,refPanels,procRank,numProcs)

refDT = dt/(2**(DT_NUMBER))

timesteps = floor(tfinal/refDt)

! MAIN Timestepping loop
do timeJ = 0,timesteps - 1
	! advance time
	call BVERK4(refParticles,refPanels,refdt,procRank,numProcs)
	t = real(timeJ+1,kreal)*refdt
	
	if ( renormalize) call RenormalizeGrid(refParticles,refPanels)

	if (procRank== 0) then
		call LogMessage(exeLog,DEBUG_LOGGING_LEVEL,trim(logKey)//" t = ",t)
	endif
enddo

call LogMessage(exeLog,DEBUG_LOGGING_LEVEL,logKey,'Reference solution done ...')

do dtJ = 1,DT_NUMBER
	call ResetRK4()
	call New(gridParticles,gridEdges,gridPanels,panelKind,initNest,AMR,nTracer,problemKind)

	call InitSolidBody(gridParticles,gridPanels,2.0_kreal*PI)
	
	if ( procRank == 0) then
		call PrintStats(gridParticles)
		call PrintStats(gridEdges)
		call PrintStats(gridPanels)
	endif
	
	dtt = dt/(2**(dtj-1))
	timesteps = floor(tfinal/dtt)
	
	do timeJ = 0,timesteps-1
		! advance time
		call BVERK4(gridParticles,gridPanels,dtt,procRank,numProcs)
		t = real(timeJ+1,kreal)*dtt
		if (procRank== 0) then
			call LogMessage(exeLog,DEBUG_LOGGING_LEVEL,trim(logKey)//" t = ",t)
		endif
		if ( renormalize) call RenormalizeGrid(gridParticles,gridPanels)		
	enddo
	
	do j=1,gridParticles%N
		gridParticles%tracer(j,3) = sqrt(sum( (gridParticles%x(:,j)-refParticles%x(:,j))*&
			(gridParticles%x(:,j)-refParticles%x(:,j))))
	enddo
	do j=1,gridPanels%N
		if ( .NOT. gridPanels%hasChildren(j) ) then
			gridPanels%tracer(j,3) = sqrt(sum( (gridPanels%x(:,j) - refPanels%x(:,j))*&
				(gridPanels%x(:,j) - refPanels%x(:,j))))
		else
			gridPanels%tracer(j,3) = 0.0_kreal
		endif
	enddo
	! find error norms
	particlesL2(dtJ) = sum(gridParticles%tracer(:,3)*gridParticles%tracer(:,3)*gridParticles%dualArea)
	particlesL2(dtJ) = sqrt(particlesL2(dtj))/particleL2denom
	panelsL2(dtj) = sum(gridPanels%tracer(:,3)*gridPanels%tracer(:,3)*gridPanels%area)
	panelsL2(dtj) = sqrt(panelsL2(dtj))/panelL2denom
	particlesLinf(dtj) = maxval(gridParticles%tracer(:,3))
	panelsLinf(dtj) = maxval(gridPanels%tracer(:,3))
	
	call Delete(gridParticles,gridEdges,gridPanels)
	write(logString,'(A,I6,A,I6,A)') '... finished trial ',dtJ,' of ',DT_NUMBER,'...'
	call LogMessage(exeLog,TRACE_LOGGING_LEVEL,logKey,logString)
enddo

if ( procRank == 0 ) then
	open(unit=writeUnit,file=dataFile,status='REPLACE',action='WRITE',iostat=writeStat)
	if (writeStat/=0) then
		call LogMessage(exeLog,ERROR_LOGGING_LEVEL,logKey,'ERROR opening dataFile.')
	else
		write(writeUnit,'(5A24)') ' dt','particleL2','particleLinf','panelL2','panelLinf'
		write(6,'(5A24)') ' dt','particleL2','particleLinf','panelL2','panelLinf'
		do j=1,DT_NUMBER
			write(writeUnit,'(5F24.15)') dt/(2**(j-1)), particlesL2(j), particlesLinf(j), panelsL2(j), panelsLinf(j)
			write(6,'(5F24.15)') dt/(2**(j-1)), particlesL2(j), particlesLinf(j), panelsL2(j), panelsLinf(j)
		enddo
	endif
	close(writeUnit)
endif

call LogMessage(exeLog,TRACE_LOGGING_LEVEL,logKey,'****** program end ****** ')
wtime = MPI_WTIME()-wtime
write(logString,'(A,F9.2,A)') " elapsed time = ",wtime/60.0_kreal," minutes."
call LogMessage(exeLog,TRACE_LOGGING_LEVEL,logKey,trim(logString))


deallocate(particlesL2)
deallocate(particlesLinf)
deallocate(panelsL2)
deallocate(panelsLinf)
call FinalizeMPIRK4()
call Delete(refParticles,refEdges,refPanels)
call Delete(exeLog)

call MPI_FINALIZE(mpiErrCode)


end program