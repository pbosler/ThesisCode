program CosineBellsDirect

use NumberKindsModule
use SphereGeomModule
use LoggerModule
use ParticlesEdgesPanelsModule
use PanelVelocityParallelModule
use TracerAndVorticityDistributionModule
use RefineRemeshModule

implicit none

include 'mpif.h'

! Grid Variables
type(Particles), pointer :: gridParticles=>null()
type(Edges), pointer :: gridEdges=>null()
type(Panels), pointer :: gridPanels=>null()
integer(kint) :: panelKind, initNest, AMR, nTracer, remeshInterval
integer(kint), parameter :: problemKind = BVE_SOLVER
logical(klog) :: remeshFlag 
integer(kint) :: problemID=COSINE_BELLS
integer(kint) :: useRelativeTol, maxRefine
real(kreal) :: amrTol1, amrTol2, maxTracer

! Logger & Computation management variables
type(Logger) :: exeLog
integer(kint) :: logOut = 6
character(len=28) :: logKey 
character(len=128) :: logString
real(kreal) :: wtime, etime
logical(klog) :: newEst
integer(kint), parameter :: REAL_INIT_BUFFER_SIZE = 4, INT_INIT_BUFFER_SIZE = 6
integer(kint) :: procRank, numProcs, mpiErrCode, intBuffer(INT_INIT_BUFFER_SIZE)
real(kreal) :: realBuffer(REAL_INIT_BUFFER_SIZE)

! Global integral variables
real(kreal), allocatable :: totalMass(:)
real(kreal) :: finalL1, finalL2, finalLinf
real(kreal) :: filamentData(2,11), tau

! Timestepping variables
real(kreal) :: t, dt, tfinal
integer(kint) :: timeJ, timesteps

! I/O & User variables
character(len=128) :: jobPrefix
character(len=256) :: vtkRoot, vtkFile, dataFile, correlationFile, filamentFile
integer(kint), parameter :: readUnit = 12, writeUnit = 13
integer(kint) :: readStat, writeStat
namelist /gridInit/ panelKind, initnest, AMR, remeshInterval, amrTol1, amrTol2, &
					useRelativeTol, maxRefine
namelist /time/ dt, tfinal
namelist /fileIO/ jobPrefix

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
	
	if ( procRank == 0) then
		call LogMessage(exeLog,TRACE_LOGGING_LEVEL,logKey,"*** Program START ***")
		call LogMessage(exeLog,DEBUG_LOGGING_LEVEL,logKey//" numProcs = ",numProcs)
		! Read user input from namelist file
		open(unit=readUnit,file='AdvectCosineBells.namelist',action='READ',status='OLD',iostat=readStat)
			if (readstat /= 0 ) then
				call LogMessage(exeLog,ERROR_LOGGING_LEVEL,logKey," ERROR opening namelist file.")
				stop
			endif
			read(readunit,nml=gridinit)
			rewind(readunit)
			read(readunit,nml=time)
			rewind(readunit)
			read(readunit,nml=fileIO)
			rewind(readunit)
		close(readunit)	
	
		intBuffer(1) = panelKind
		intBuffer(2) = initNest
		intBuffer(3) = AMR
		intBuffer(4) = remeshInterval
		intBuffer(5) = useRelativeTol
		intBuffer(6) = maxRefine
	
		realBuffer(1) = amrTol1
		realBuffer(2) = amrTol2
		realBuffer(3) = dt
		realBuffer(4) = tfinal
		! Prepare output files
		if ( panelKind == TRI_PANEL ) then
			if ( AMR == 0) then
				write(dataFile,'(A,A,I1,A,F3.1,A,F5.3)') trim(jobPrefix),'_DIRECT_triNest',initNest,'_rev',tfinal,'_dt',dt
				write(vtkRoot,'(A,A,A,I1,A,F3.1,A,F5.3,A)') 'vtkOut/',trim(jobPrefix),'_DIRECT_triNest',initNest,'_rev',tfinal,'_dt',dt,'_'
			else
				write(dataFile,'(A,A,I1,A,F3.1,A,F5.3)') trim(jobPrefix),'_DIRECT_triAMR',initNest,'_rev',tfinal,'_dt',dt
				write(vtkRoot,'(A,A,A,I1,A,F3.1,A,F5.3,A)') 'vtkOut/',trim(jobPrefix),'_DIRECT_triAMR',initNest,'_rev',tfinal,'_dt',dt,'_'
			endif
		elseif (panelKind == QUAD_PANEL ) then
			if ( AMR == 0) then
				write(dataFile,'(A,A,I1,A,F3.1,A,F5.3)') trim(jobPrefix),'_DIRECT_quadNest',initNest,'_rev',tfinal,'_dt',dt
				write(vtkRoot,'(A,A,A,I1,A,F3.1,A,F5.3,A)') 'vtkOut/',trim(jobPrefix),'_DIRECT_quadNest',initNest,'_rev',tfinal,'_dt',dt,'_'
			else
				write(dataFile,'(A,A,I1,A,F3.1,A,F5.3)') trim(jobPrefix),'_DIRECT_quadAMR',initNest,'_rev',tfinal,'_dt',dt
				write(vtkRoot,'(A,A,A,I1,A,F3.1,A,F5.3,A)') 'vtkOut/',trim(jobPrefix),'_DIRECT_quadAMR',initNest,'_rev',tfinal,'_dt',dt,'_'
			endif
		else
			call LogMessage(exeLog,ERROR_LOGGING_LEVEL,logkey,' ERROR : Invalid panelKind.')
			stop
		endif
		write(vtkFile,'(A,I0.4,A)') trim(vtkRoot),timeJ,'.vtk'
		write(correlationFile,'(A,A)') trim(dataFile),'_corr.dat'
		write(filamentFile,'(A,A)') trim(dataFile),'_filament.dat'
		write(dataFile,'(A,A)') trim(dataFile),'.dat'
	endif
	call MPI_BCAST(intBuffer,INT_INIT_BUFFER_SIZE,MPI_INTEGER,0,MPI_COMM_WORLD,mpiErrCode)
	call MPI_BCAST(realBuffer,REAL_INIT_BUFFER_SIZE,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiErrCode)
	
	panelKind = intBuffer(1)
	initNest = intBuffer(2)
	AMR = intBuffer(3)
	remeshInterval = intBuffer(4)
	useRelativeTol = intBuffer(5)
	maxRefine = intBuffer(6)
	
	amrTol1 = realBuffer(1)
	amrTol2 = realBuffer(2)
	dt = realBuffer(3)
	tfinal = realBuffer(4)

	! Init time stepping variables
	t = 0.0_kreal
	timeJ = 0
	timesteps = floor(tfinal/dt)
	
	allocate(totalMass(0:timesteps))
	totalMass = 0.0_kreal
	
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!	Part 2 : Initialize the grid		             !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
	
	if (procRank == 0) wtime = MPI_WTIME()
	
	nTracer = 3
	call New(gridParticles,gridEdges,gridPanels,panelKind,initNest,AMR,nTracer,problemKind)
	call InitCosineBells(gridParticles,gridPanels)
	
	if ( AMR> 0 ) then
		call SetMaxRefinementLimit(maxRefine)
		if ( useRelativeTol > 0 ) then
			maxTracer = maxVal(abs(gridPanels%tracer(1:gridPanels%N,1))*gridPanels%area(1:gridPanels%N))
			amrTol1 = maxTracer*amrTol1	
			if ( procRank == 0 ) then
				call LogMessage(exeLog,TRACE_LOGGING_LEVEL,trim(logKey)//' amrTol1 = ',amrTol1)
			endif
		endif
		call InitRefine(gridParticles,gridEdges,gridPanels,amrTol1,amrTol2, problemID, procRank=procRank)
		call InitCosineBells(gridParticles,gridPanels)	
	endif ! AMR
	
	! Initialize correlated tracer
	do j=1,gridParticles%N
		gridParticles%tracer(j,3) = -0.8_kreal*gridParticles%tracer(j,1)*gridParticles%tracer(j,1) + 0.9_kreal
	enddo
	do j=1,gridPanels%N
		if ( .NOT. gridPanels%hasChildren(j)) then
		gridPanels%tracer(j,3) = -0.8_kreal*gridPanels%tracer(j,1)*gridPanels%tracer(j,1) + 0.9_Kreal
		endif
	enddo
	
	totalMass(0) = sum(gridPanels%tracer(1:gridPanels%N,1)*gridPanels%area(1:gridPanels%N))
	
	! Calculate initial filament diagnostics
	filamentData = 0.0_kreal
	do k=1,11
		tau = (k-1)*0.1_kreal
		do j=1,gridPanels%N
			if ( .NOT. gridPanels%hasChildren(j) ) then
				if ( gridPanels%tracer(j,1) > tau) then
					filamentData(1,k) = filamentData(1,k) + gridPanels%area(j)
				endif
			endif
		enddo
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

	call MPI_BARRIER(MPI_COMM_WORLD,mpiErrCode)
	call LogMessage(exeLog,DEBUG_LOGGING_LEVEL,logKey,'Setup complete. Starting time integration.')
	call InitializeMPIRK4(gridParticles,gridPanels,procRank,numProcs)

	remeshFlag = .False.
	newEst = .False.
	do timeJ = 0,timesteps-1
		! Remesh if necessary
		if ( mod(timeJ+1,remeshInterval) == 0 ) then
			remeshFlag = .True.
			if ( procRank == 0 ) &
			call LogMessage(exeLog,TRACE_LOGGING_LEVEL,logKey," Remesh triggered by remeshInterval.")
		endif
		if ( remeshFlag ) then
		
			call DirectAdaptiveRemesh( gridParticles,gridEdges,gridPanels,&
									initNest, AMR, amrTol1, amrTol2, &
									procRank, problemID)
!			do j=1,gridParticles%N
!				gridParticles%tracer(j,3) = -0.8_kreal*gridParticles%tracer(j,1)*gridParticles%tracer(j,1) + 0.9_kreal
!			enddo
!			do j=1,gridPanels%N
!				if ( .NOT. gridPanels%hasChildren(j)) then
!				gridPanels%tracer(j,3) = -0.8_kreal*gridPanels%tracer(j,1)*gridPanels%tracer(j,1) + 0.9_Kreal
!				endif
!			enddo									
		
			remeshFlag = .False.

			call ResetRK4()
			
			
			if ( procRank == 0 ) then
				call PrintStats(gridParticles)
				call PrintStats(gridEdges)
				call PrintStats(gridPanels)
			endif
		
			newEst = .TRUE.
		endif
	
		if ( ( procRank == 0 .AND. timeJ == 1) .OR. (procRank == 0 .AND. newEst ) ) etime = MPI_WTIME()
		! Advance time
		call AdvectionRK4(gridParticles,gridPanels,t,dt,procRank,numProcs)
		t = real(timeJ+1,kreal)*dt
		
		totalMass(timeJ+1) = sum(gridPanels%tracer(1:gridPanels%N,1)*gridPanels%area(1:gridPanels%N))
		if (timeJ+1 == timesteps/2) then
			do k=1,11
				tau = (k-1)*0.1_kreal
				do j=1,gridPanels%N
					if ( .NOT. gridPanels%hasChildren(j) ) then
						if ( gridPanels%tracer(j,1) > tau) then
							filamentData(2,k) = filamentData(2,k) + gridPanels%area(j)
						endif
					endif
				enddo
			enddo
		endif

		if ( procRank == 0 ) then
			if ( ( timeJ == 1) .OR. ( newEst) ) then
				etime = MPI_WTIME() - etime
				write(logString,'(A,F9.3,A)') 'Estimated time left = ',(timesteps-timeJ)*etime/60.0,' minutes.'
				call LogMessage(exeLog,TRACE_LOGGING_LEVEL,logKey,logString)
				newEst = .False.
			endif
	
			call LogMessage(exeLog,DEBUG_LOGGING_LEVEL,trim(logKey)//" t = ",t)
			write(vtkFile,'(A,I0.4,A)') trim(vtkRoot),timeJ+1,'.vtk'
			call vtkOutput(gridParticles,gridPanels,vtkFile)
		endif
	enddo	
	
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!	Part 4 : Finish and clear memory				!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

	! Calculate error
	do j=1,gridPanels%N
		if ( .NOT. gridPanels%hasChildren(j) ) then
			gridPanels%tracer(j,2) = abs( CosineBellsX(gridPanels%x(:,j)) - gridPanels%tracer(j,1))
		endif
	enddo
	
	finalL1 = sum( gridPanels%tracer(1:gridPanels%N,2)*gridPanels%area(1:gridPanels%N))/&
			sum(abs(gridPanels%tracer(1:gridPanels%N,1))*gridPanels%area(1:gridPanels%N))
	finalL2 = sqrt(sum( gridPanels%tracer(1:gridPanels%N,2)*gridPanels%tracer(1:gridPanels%N,2)*gridPanels%area(1:gridPanels%N)))&
			/ sqrt(sum( gridPanels%tracer(1:gridPanels%N,1)*gridPanels%tracer(1:gridPanels%N,1)*gridPanels%area(1:gridPanels%N)))
	finalLinf = maxval(gridPanels%tracer(1:gridPanels%N,2))/maxVal(abs(gridPanels%tracer(1:gridPanels%N,1)))			
	
	
	
	if ( procRank == 0 ) then
		write(6,'(A,F24.18)') ' L1 Error at t = 5 : ', finalL1
		write(6,'(A,F24.18)') ' L2 Error at t = 5 : ', finalL2
		write(6,'(A,F24.18)') 'Linf Error at t = 5: ', finalLinf
	
		open(unit=writeUnit,file=dataFile,status='REPLACE',action='WRITE',iostat=writeStat)
			if ( writeStat /= 0 ) call LogMessage(exeLog,ERROR_LOGGING_LEVEL,logkey,'ERROR opening dataFile.')
		write(writeUnit,'(A,F24.18)') ' L1 Error at t = 5 : ', finalL1
		write(writeUnit,'(A,F24.18)') ' L2 Error at t = 5 : ', finalL2
		write(writeUnit,'(A,F24.18)') 'Linf Error at t = 5: ', finalLinf
		write(writeUnit,'(A)') ' '
		write(writeUnit,'(2A24)') 'time', 'totalMass'
		do j=0,timesteps
			write(writeUnit,'(2F24.18)') j*dt, totalMass(j)
		enddo
		close(writeUnit)
		
		open(unit=writeUnit,file=filamentFile,status='REPLACE',action='WRITE',iostat=writeStat)
			if ( writeStat /=0 ) call LogMessage(exeLog,ERROR_LOGGING_LEVEL,logKey,'ERROR opening filamentFile')
			do j=1,11
				write(writeUnit,'(2F24.16)') filamentData(:,j)
			enddo
		close(writeUnit)
		
		open(unit=writeUnit,file=correlationFile,status='REPLACE',action='WRITE',iostat=writeStat)
			if (writeStat /= 0) call LogMessage(exeLog,ERROR_LOGGING_LEVEL,logKey,'ERROR opening correlationFile')
			do j=1,gridPanels%N
				if ( .NOT. gridPanels%hasChildren(j) ) then
					write(writeUnit,'(2F24.16)') gridPanels%tracer(j,1), gridPanels%tracer(j,3)
				endif
			enddo
		close(writeUnit)
		
		wTime = MPI_WTIME() - wtime
		call LogMessage(exeLog,TRACE_LOGGING_LEVEL,logKey,"*** Program END ***")
		call LogMessage(exeLog,TRACE_LOGGING_LEVEL,trim(logKey)//' data file = ',trim(dataFile))
		write(vtkFile,'(A,A)') trim(vtkRoot),'XXXX.vtk'
		call LogMessage(exeLog,TRACE_LOGGING_LEVEL,trim(logKey)//' vtkFiles = ',trim(vtkFile))
		if ( AMR > 0 ) then
			if ( useRelativeTOl > 0 ) call LogMessage(exeLog,TRACE_LOGGING_LEVEL,logKey,"AMR used relative tolerances.")
			call LogMessage(exeLog,TRACE_LOGGING_LEVEL,logKey//'panelKind = ',panelKind)
			call LogMessage(exeLog,TRACE_LOGGING_LEVEL,logkey//'initNest = ',initNest)
		endif
		write(logString,'(A,F9.2,A)') " elapsed time = ",wtime/60.0_kreal," minutes."
		call LogMessage(exeLog,TRACE_LOGGING_LEVEL,logKey,trim(logString))
	endif
	
	deallocate(totalMass)
	call FinalizeMPIRK4()
	call Delete(gridParticles,gridEdges,gridPanels)
	call Delete(exeLog)
	
	call MPI_FINALIZE(mpiErrCode)	
end program
