program GaussianVortex

use NumberKindsModule
use SphereGeomModule
use LoggerModule
use ParticlesEdgesPanelsModule
use PanelVelocityParallelModule
use TracerAndVorticityDistributionModule
use RefineRemeshModule

implicit none

include 'mpif.h'

! Grid variables
type(Particles), pointer :: gridParticles=>null()
type(Edges), pointer :: gridEdges=>null()
type(Panels), pointer :: gridPanels=>null()
integer(kint) :: panelKind, initNest, AMR, nTracer, remeshInterval
integer(kint), parameter :: problemKind = BVE_SOLVER, problemID = SINGLE_GAUSSIAN_VORTEX
logical(klog) :: remeshFlag
real(kreal) :: amrTol1, circTol, amrTol2, baseVar, varTol, minX0(3), maxX0(3), maxCirc

! Data collection variables
real(kreal), allocatable :: totalKineticEnergy(:), totalEnstrophy(:)

! Logger and computation management variables
type(Logger) :: exeLog
integer(kint) :: logOut = 6
character(len=28) :: logKey
character(len=128) :: logString
real(kreal) :: wtime, etime
integer(kint), parameter :: REAL_INIT_BUFFER_SIZE = 4, INTEGER_INIT_BUFFER_SIZE = 4
integer(kint) :: procRank, numProcs, mpiErrCode, initIntegers(INTEGER_INIT_BUFFER_SIZE)
real(kreal) :: initReals(REAL_INIT_BUFFER_SIZE)
logical(klog) :: newEstimate

! Timestepping variables
real(kreal) :: t, dt, tfinal
integer(kint) :: timeJ, timesteps

! I/O and User variables
character(len=128) :: jobPrefix, outputDir
character(len=48) :: amrString
character(len=256) :: vtkRoot, vtkFile, dataFile, vtkRoot2, vtkFile2
integer(kint), parameter :: readUnit = 12, writeUnit = 13
integer(kint) :: readStat, writeStat, frameOut, frameCounter

namelist /gridInit/ panelKind, initNest, AMR, remeshInterval, amrTol1, amrTol2
namelist /time/ dt, tfinal
namelist /fileIO/ jobPrefix, outputDir, frameOut

! General variables
integer(kint) :: j, k

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!	Part 1 : Initialize the computing environment / get user input  !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

call MPI_INIT(mpiErrCode)
call MPI_COMM_SIZE(MPI_COMM_WORLD,numProcs,mpiErrCode)
call MPI_COMM_RANK(MPI_COMM_WORLD,procRank,mpiErrCode)

call New(exeLog,DEBUG_LOGGING_LEVEL,logOut)
write(logKey,'(A,I0.2,A)') 'EXE_LOG_',procRank,' : '

if (procRank == 0) then
	call LogMessage(exeLog,TRACE_LOGGING_LEVEL,logKey,"******  Program START ******")
	call LogMessage(exeLog,DEBUG_LOGGING_LEVEL,logKey//" numProcs = ",numProcs)

	! Read user input from namelist file
	open(unit=readUnit,file='GaussVortPanels.namelist',action='READ',status='OLD',iostat=readStat)
		if (readStat /= 0) then
			call LogMessage(exeLog,ERROR_LOGGING_LEVEL,logKey,' ERROR opening namelist file.')
			stop
		endif
		read(readUnit,nml=gridInit)
		rewind(readUnit)
		read(readUnit,nml=time)
		rewind(readUnit)
		read(readUnit,nml=fileIO)
		rewind(readUnit)
	close(readUnit)

	initIntegers(1) = panelKind
	initIntegers(2) = initNest
	initIntegers(3) = AMR
	initIntegers(4) = remeshInterval

	initReals(1) = amrTol1
	initReals(2) = amrTol2
	initReals(3) = dt
	initReals(4) = tfinal
endif

call MPI_BCAST(initIntegers,INTEGER_INIT_BUFFER_SIZE,MPI_INTEGER,0,MPI_COMM_WORLD,mpiErrCode)
call MPI_BCAST(initReals,REAL_INIT_BUFFER_SIZE,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiErrCode)

panelKind = initIntegers(1)
initNest = initIntegers(2)
AMR = initIntegers(3)
remeshInterval = initIntegers(4)

amrTol1 = initReals(1)
amrTol2 = initReals(2)
dt = initReals(3)
tfinal = initReals(4)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!	Part 2 : Initialize the grid		             !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

t = 0.0_kreal
timeJ = 0
timesteps = floor(tfinal/dt)

allocate(totalKineticEnergy(0:timesteps))
totalKineticEnergy = 0.0_kreal
allocate(totalEnstrophy(0:timesteps))
totalEnstrophy = 0.0_kreal

nTracer = 2

! build uniform grid to initNest
call New(gridParticles,gridEdges,gridPanels,panelKind,initNest,AMR,nTracer,problemKind)
call InitGaussianVortex(gridParticles,gridPanels)

if (AMR>0) then
	call SetMaxRefinementLimit(2)
	! Set AMR criteria relative to max and min rel vort and initial grid spacing
	maxCirc = maxval(abs(gridPanels%relVort(1:gridPanels%N))*gridPanels%area(1:gridPanels%N))
	baseVar = 0.0_kreal
	do j=1,gridPanels%N
		if ( .NOT. gridPanels%hasChildren(j) ) then
			minx0 = gridPanels%x0(:,j)
			maxx0 = gridPanels%x0(:,j)
			do k=1,panelKind
				if ( gridParticles%x0(1,gridPanels%vertices(k,j)) < minx0(1)) then
					minx0(1) = gridParticles%x0(1,gridPanels%vertices(k,j))
				endif
				if ( gridParticles%x0(1,gridPanels%vertices(k,j)) > maxx0(1)) then
					maxx0(1) = gridParticles%x0(1,gridPanels%vertices(k,j))
				endif
				if ( gridParticles%x0(2,gridPanels%vertices(k,j)) < minx0(2)) then
					minx0(2) = gridParticles%x0(2,gridPanels%vertices(k,j))
				endif
				if ( gridParticles%x0(2,gridPanels%vertices(k,j)) > maxx0(2)) then
					maxx0(2) = gridParticles%x0(2,gridPanels%vertices(k,j))
				endif
				if ( gridParticles%x0(3,gridPanels%vertices(k,j)) < minx0(3)) then
					minx0(3) = gridParticles%x0(3,gridPanels%vertices(k,j))
				endif
				if ( gridParticles%x0(3,gridPanels%vertices(k,j)) > maxx0(3)) then
					maxx0(3) = gridParticles%x0(3,gridPanels%vertices(k,j))
				endif
			enddo
			if ( sum(maxx0-minx0) > baseVar) baseVar = sum(maxx0-minx0)
		endif
	enddo

	circTol = amrTol1*maxCirc
	varTol = amrTol2*baseVar
	if ( procRank == 0 ) then
		call LogMessage(exeLog,TRACE_LOGGING_LEVEL,trim(logKey)//' circTol = ',circTol)
		call LogMessage(exeLog,TRACE_LOGGING_LEVEL,trim(logKey)//' varTol = ',varTol)
		call LogMessage(exeLog,TRACE_LOGGING_LEVEL,trim(logKey)//' maxNest will be ',initNest+MAX_REFINEMENT+1)

		write(amrString,'(A,I1,A,I0.2,A)') 'AMR_',initNest,'to',initNest+MAX_REFINEMENT+1,'_'

	endif

	call InitRefine(gridParticles,gridEdges,gridPanels,circTol,varTol,problemID,procRank=procRank)
	call InitGaussianVortex(gridParticles,gridPanels)
else
	write(amrString,'(A,I1)') 'nest',initNest
endif

! Store init latitude in tracer 1
do j=1,gridParticles%N
	gridParticles%tracer(j,1) = Latitude(gridParticles%x0(:,j))
enddo
do j=1,gridPanels%N
	if ( .NOT. gridPanels%hasChildren(j) ) then
		gridPanels%tracer(j,1) = Latitude(gridPanels%x0(:,j))
	endif
enddo

totalEnstrophy(0) = 0.5_kreal*sum(gridPanels%relVort(1:gridPanels%N) * &
								  gridPanels%relVort(1:gridPanels%N) * &
								  gridPanels%area(1:gridPanels%N) )

frameCounter = 0
if ( procRank == 0 ) then
	call PrintStats(gridParticles)
	call PrintStats(gridEdges)
	call PrintStats(gridPanels)

	! Prepare output files
	if ( panelKind == 3) then
		write(dataFile,'(A,A,A,A,A,F4.2,A,F6.4,A)') trim(outputDir),trim(jobPrefix),'_tri',trim(amrString),'_rev',tfinal,'_dt',dt,'_'
	elseif (panelKind == 4) then
		write(dataFile,'(A,A,A,A,F4.2,A,F6.4,A)') trim(outputDir),trim(jobPrefix),'_quad',trim(amrString),'_rev',tfinal,'_dt',dt,'_'
	endif
	write(vtkRoot,'(A,A,A)') trim(outputDir),'vtkOut/',trim(jobPrefix)
	write(vtkFile,'(A,I0.4,A)') trim(vtkRoot),frameCounter,'.vtk'
	write(dataFile,'(A,A)')trim(dataFile),'.dat'

	! output initial data
	call vtkOutput(gridParticles,gridPanels,vtkFile)


endif

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!	Part 3 : Run the problem						!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

call MPI_BARRIER(MPI_COMM_WORLD,mpiErrCode)
if ( procRank == 0 ) call LogMessage(exeLog,DEBUG_LOGGING_LEVEL,logKey,'Setup complete. Starting time integration.')

call InitializeMPIRK4(gridParticles,gridPanels,procRank,numProcs)

remeshFlag = .False.
newEstimate = .False.

if ( procRank == 0) wtime = MPI_WTIME()


do timeJ=0, timesteps-1
	if ( mod(timeJ+1,remeshInterval) == 0 ) then
		remeshFlag = .TRUE.
		newEstimate = .TRUE.
		if ( procRank == 0 ) call LogMessage(exeLog,TRACE_LOGGING_LEVEL,logKey," remesh triggered by remeshInterval.")
	endif
	if ( remeshFlag ) then
		call AdaptiveRemesh(gridParticles,gridEdges,gridPanels, &
							initNest, AMR, circTol, varTol, procRank, problemID)
		remeshFlag = .FALSE.

		call ResetRK4()

		do j=1,gridParticles%N
			gridParticles%tracer(j,1) = Latitude(gridParticles%x0(:,j))
		enddo
		do j=1,gridPanels%N
			if ( .NOT. gridPanels%hasChildren(j)) then
				gridPanels%tracer(j,1) = Latitude(gridPanels%x0(:,j))
			endif
		enddo
		if ( procRank == 0) then
			call PrintStats(gridParticles)
			call PrintStats(gridEdges)
			call PrintStats(gridPanels)
		endif
	endif

	! Advance time
	call BVERK4(gridParticles,gridPanels,dt,procRank,numProcs)
	t = real(timeJ+1,kreal)*dt

	if ( timeJ == 0) totalKineticEnergy(0) = totalKE
	totalKineticEnergy(timeJ+1) = totalKE
	totalEnstrophy(timeJ+1) = 0.5_kreal*sum(gridPanels%relVort(1:gridPanels%N) * &
										    gridPanels%relVort(1:gridPanels%N) * &
										    gridPanels%area(1:gridPanels%N) )

	if (newEstimate .AND. procRank == 0) then
		etime = MPI_WTIME() - wtime
		write(logString,'(A,F9.4,A)') 'Estimated time left = ', etime*real(timesteps,kreal)/(60.0_kreal*real(timeJ+1,kreal)), ' minutes.'
		call LogMessage(exeLog,TRACE_LOGGING_LEVEL,logKey,trim(logString))
		newEstimate = .FALSE.
	endif

	if (procRank== 0) then
		call LogMessage(exeLog,DEBUG_LOGGING_LEVEL,trim(logKey)//" t = ",t)
		if ( mod(timeJ+1,frameOut) == 0 ) then
			frameCounter = frameCounter + 1
			write(vtkFile,'(A,I0.4,A)') trim(vtkRoot),frameCounter,'.vtk'
			call vtkOutput(gridParticles,gridPanels,vtkFile)
		endif
	endif

enddo

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!	Part 4 : Finish and clear memory				!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

if ( procRank == 0) then
	open(unit=writeUnit,file=dataFile,status='REPLACE',action='WRITE',iostat=writeStat)
		if ( writeStat /= 0 ) then
			call LogMessage(exeLog,ERROR_LOGGING_LEVEL,logKey,"ERROR opening dataFile.")
		else
			write(writeUnit,'(2A24)') 'totalKE','totalEnstrophy'
			do j=0,timesteps
				write(writeUnit,'(2F24.15)') totalKineticEnergy(j), totalEnstrophy(j)
			enddo
		endif
	close(writeUnit)

	wtime = MPI_WTIME() - wtime
	call LogMessage(exeLog,TRACE_LOGGING_LEVEL,logKey,"****** program END ******")
	call LogMessage(exeLog,TRACE_LOGGING_LEVEL,trim(logKey)//" dataFile = ",trim(dataFile))
	write(vtkFile,'(A,A)') trim(vtkRoot),'XXXX.vtk'
	call LogMessage(exeLog,TRACE_LOGGING_LEVEL,trim(logKey)//' vtkFiles = ',trim(vtkFile))
	call LogMessage(exeLog,TRACE_LOGGING_LEVEL,logKey,' tracer 1 = initial latitude')
	call LogMEssage(exeLog,TRACE_LOGGING_LEVEL,logKey,' tracer 2 = kinetic energy ')
	call StartSection(exeLog,'-------- Calculation Summary -------------')
		if (panelKind == QUAD_PANEL .AND. AMR == 0) then
			write(logString,'(A,I8)') 'Cubed sphere nActive = ',gridPanels%N_Active
		elseif ( panelKind == TRI_PANEL .AND. AMR == 0) then
			write(logString,'(A,I8)') 'Icos. triangles nActive = ',gridPanels%N_Active
		elseif ( AMR > 0 ) then
			write(logString,'(A,F14.7,A,F14.7)') 'circTol = ',circTol,' varTol = ', varTol
		endif
		call LogMessage(exeLog,TRACE_LOGGING_LEVEL,logkey,trim(logString))
		write(logString,'(A,F12.8)') ' dt = ',dt
		call LogMessage(exeLog,TRACE_LOGGING_LEVEL,logkey,trim(logString))
		write(logString,'(A,I6)') 'remeshInterval = ',remeshInterval
		call LogMessage(exeLog,TRACE_LOGGING_LEVEL,logkey,trim(logString))
	call EndSection(exeLog)
	write(logString,'(A,F9.2,A)') " elapsed time = ",wtime/60.0_kreal," minutes."
	call LogMessage(exeLog,TRACE_LOGGING_LEVEL,logKey,trim(logString))
endif

deallocate(totalKineticEnergy)
deallocate(totalEnstrophy)
call FinalizeMPIRK4()
call Delete(gridParticles,gridEdges,gridPanels)
call Delete(exeLog)

call MPI_FINALIZE(mpiErrCode)

end program
