program ZonalMean

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
integer(kint), parameter :: problemID = ZONAL_MEAN
integer(kint) :: useRelativeTol
real(kreal) :: amrTol1, amrTol2, newOmega
real(kreal) :: circTol, varTol, maxCirc, baseVar, maxx0(3), minx0(3)

! Error calculation variables
real(kreal), allocatable :: totalKineticEnergy(:), totalEnstrophy(:)

! Logger & Computation management variables
type(Logger) :: exeLog
integer(kint) :: logOut = 6
character(len=28) :: logKey 
character(len=128) :: logString
real(kreal) :: wtime, etime
logical(klog) :: newEst
integer(kint), parameter :: REAL_INIT_BUFFER_SIZE = 4, INT_INIT_BUFFER_SIZE = 5
integer(kint) :: procRank, numProcs, mpiErrCode, intBuffer(INT_INIT_BUFFER_SIZE)
real(kreal) :: realBuffer(REAL_INIT_BUFFER_SIZE)

! Timestepping variables
real(kreal) :: t, dt, tfinal
integer(kint) :: timeJ, timesteps

! I/O & User variables
character(len=128) :: jobPrefix
character(len=256) :: vtkRoot, vtkFile, dataFile
integer(kint), parameter :: readUnit = 12, writeUnit = 13
integer(kint) :: readStat, writeStat
namelist /gridInit/ panelKind, initnest, AMR, remeshInterval, amrTol1, amrTol2, useRelativeTol, newOmega
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
		open(unit=readUnit,file='ZonalMeanPanels.namelist',action='READ',status='OLD',iostat=readStat)
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
	
		realBuffer(1) = amrTol1
		realBuffer(2) = amrTol2
		realBuffer(3) = dt
		realBuffer(4) = tfinal
	
	
		! Prepare output files
		if ( panelKind == TRI_PANEL ) then
			if ( AMR == 0) then
				write(dataFile,'(A,A,I1,A,F3.1,A,F5.3,A)') trim(jobPrefix),'triNest',initNest,'_rev',tfinal,'_dt',dt,'.dat'
				write(vtkRoot,'(A,A,A,I1,A,F3.1,A,F5.3,A)') 'vtkOut/',trim(jobPrefix),'triNest',initNest,'_rev',tfinal,'_dt',dt,'_'
			else
				write(dataFile,'(A,A,I1,A,F3.1,A,F5.3,A)') trim(jobPrefix),'triAMR',initNest,'_rev',tfinal,'_dt',dt,'.dat'
				write(vtkRoot,'(A,A,A,I1,A,F3.1,A,F5.3,A)') 'vtkOut/',trim(jobPrefix),'triAMR',initNest,'_rev',tfinal,'_dt',dt,'_'
			endif
		elseif (panelKind == QUAD_PANEL ) then
			if ( AMR == 0) then
				write(dataFile,'(A,A,I1,A,F3.1,A,F5.3,A)') trim(jobPrefix),'quadNest',initNest,'_rev',tfinal,'_dt',dt,'.dat'
				write(vtkRoot,'(A,A,A,I1,A,F3.1,A,F5.3,A)') 'vtkOut/',trim(jobPrefix),'quadNest',initNest,'_rev',tfinal,'_dt',dt,'_'
			else
				write(dataFile,'(A,A,I1,A,F3.1,A,F5.3,A)') trim(jobPrefix),'quadAMR',initNest,'_rev',tfinal,'_dt',dt,'.dat'
				write(vtkRoot,'(A,A,A,I1,A,F3.1,A,F5.3,A)') 'vtkOut/',trim(jobPrefix),'quadAMR',initNest,'_rev',tfinal,'_dt',dt,'_'
			endif
		else
			call LogMessage(exeLog,ERROR_LOGGING_LEVEL,logkey,' ERROR : Invalid panelKind.')
			stop
		endif
		write(vtkFile,'(A,I0.4,A)') trim(vtkRoot),timeJ,'.vtk'
	endif
	call MPI_BCAST(intBuffer,INT_INIT_BUFFER_SIZE,MPI_INTEGER,0,MPI_COMM_WORLD,mpiErrCode)
	call MPI_BCAST(realBuffer,REAL_INIT_BUFFER_SIZE,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiErrCode)

	panelKind = intBuffer(1)
	initNest = intBuffer(2)
	AMR = intBuffer(3)
	remeshInterval = intBuffer(4)
	useRelativeTol = intBuffer(5)

	amrTol1 = realBuffer(1)
	amrTol2 = realBuffer(2)
	dt = realBuffer(3)
	tfinal = realBuffer(4)

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
	if (procRank == 0) wtime = MPI_WTIME()

	nTracer = 2
	call New(gridParticles,gridEdges,gridPanels,panelKind,initNest,AMR,nTracer,problemKind)

	call InitZonalMean(gridParticles,gridPanels)
	if ( AMR > 0 ) then
		if ( useRelativeTol > 0 ) then
			maxCirc = maxval(abs(gridPanels%relVort(1:gridPanels%N))*gridPanels%area(1:gridPanels%N))
			baseVar = 0.0_kreal	
			do j=1,gridPanels%N
				if ( .NOT. gridPanels%hasChildren(j) ) then
					maxx0 = gridPanels%x0(:,j)
					minx0 = gridPanels%x0(:,j)
					do k=1,panelKind
						if ( gridParticles%x0(1,gridPanels%vertices(k,j)) > maxx0(1) ) then
							maxx0(1) = gridParticles%x0(1,gridPanels%vertices(k,j))
						endif
						if ( gridParticles%x0(1,gridPanels%vertices(k,j)) < minx0(1) ) then
							minx0(1) = gridParticles%x0(1,gridPanels%vertices(k,j))
						endif
						if ( gridParticles%x0(2,gridPanels%vertices(k,j)) > maxx0(2) ) then
							maxx0(2) = gridParticles%x0(2,gridPanels%vertices(k,j))
						endif
						if ( gridParticles%x0(2,gridPanels%vertices(k,j)) < minx0(2) ) then
							minx0(2) = gridParticles%x0(2,gridPanels%vertices(k,j))
						endif
						if ( gridParticles%x0(3,gridPanels%vertices(k,j)) > maxx0(3) ) then
							maxx0(3) = gridParticles%x0(3,gridPanels%vertices(k,j))
						endif
						if ( gridParticles%x0(3,gridPanels%vertices(k,j)) < minx0(3) ) then
							minx0(3) = gridParticles%x0(3,gridPanels%vertices(k,j))
						endif
					enddo
					if ( sum(maxx0 - minx0) > baseVar ) baseVar = sum(maxx0 - minx0)
				endif
			enddo

			circTol = amrTol1 * maxCirc
			varTol = amrTol2 * baseVar
			if ( procRank == 0 ) then
				call LogMessage(exeLog,TRACE_LOGGING_LEVEL,trim(logKey)//' circTol = ',circTol)
				call LogMessage(exeLog,TRACE_LOGGING_LEVEL,trim(logKey)//' varTol = ',varTol)
			endif
		else
			circTol = amrTol1
			varTol = amrTol2
		endif
		call InitRefine(gridParticles,gridEdges,gridPanels,circTol,varTol, problemID, procRank=procRank)
		call InitZonalMean(gridParticles,gridPanels)	
	endif
	
	!
	! Store initial latitude in tracer 1
	!
	do j=1,gridParticles%N
		gridParticles%tracer(j,1) = Latitude(gridParticles%x0(:,j))
	enddo
	do j=1,gridPanels%N
		if ( .NOT. gridPanels%hasChildren(j) )  then
			gridpanels%tracer(j,1) = Latitude(gridPanels%x0(:,j))
		endif
	enddo
	!
	! Tracer 2 used for kinetic energy
	!
	totalEnstrophy(0) = 0.5_kreal*sum(gridPanels%area(1:gridpanels%N)*&
					gridPanels%relVort(1:gridPanels%N)*gridPanels%relVort(1:gridPanels%N))

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
		
			call AdaptiveRemesh( gridParticles,gridEdges,gridPanels,&
									initNest, AMR, circTol, varTol, &
									procRank, problemID)
		
			remeshFlag = .False.

			call ResetRK4()

			do j=1,gridParticles%N
				gridParticles%tracer(j,1) = Latitude(gridParticles%x0(:,j))
			enddo
			do j=1,gridPanels%N
				if ( .NOT. gridPanels%hasChildren(j) )  then
					gridpanels%tracer(j,1) = Latitude(gridPanels%x0(:,j))
				endif
			enddo
			
			if ( procRank == 0 ) then
				call PrintStats(gridParticles)
				call PrintStats(gridEdges)
				call PrintStats(gridPanels)
			endif
		
			newEst = .TRUE.
		endif
	
	
	
		if ( ( procRank == 0 .AND. timeJ == 1) .OR. (procRank == 0 .AND. newEst ) ) etime = MPI_WTIME()
		! Advance time
		call BVERK4(gridParticles,gridPanels,dt,procRank,numProcs)
		t = real(timeJ+1,kreal)*dt
	
		totalKineticEnergy(timeJ+1) = totalKE
		totalEnstrophy(timeJ+1) = 0.5_kreal*sum(gridPanels%area(1:gridPanels%N)*&
				gridPanels%relVort(1:gridPanels%N)*gridPanels%relVort(1:gridPanels%N))
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
	if ( procRank == 0 ) then
		open(unit=writeUnit,file=dataFile,status='REPLACE',action='WRITE',iostat=writeStat)
		if ( writeStat /= 0 ) call LogMessage(exeLog,ERROR_LOGGING_LEVEL,logkey,'ERROR opening dataFile.')
			write(writeUnit,'(2A24)') 'totalKE','totalEnstrophy'
			do j=0,timesteps
				write(writeUnit,'(2F24.15)') totalKineticEnergy(j), totalEnstrophy(j)
			enddo
		close(writeUnit)
	
		wTime = MPI_WTIME() - wtime
		call LogMessage(exeLog,TRACE_LOGGING_LEVEL,logKey,"*** Program END ***")
		call LogMessage(exeLog,TRACE_LOGGING_LEVEL,trim(logKey)//' data file = ',trim(dataFile))
		write(vtkFile,'(A,A)') trim(vtkRoot),'XXXX.vtk'
		call LogMessage(exeLog,TRACE_LOGGING_LEVEL,trim(logKey)//' vtkFiles = ',trim(vtkFile))
		call LogMessage(exeLog,TRACE_LOGGING_LEVEL,logKey,'Tracer1 = initial latitude.')
		call LogMessage(exeLog,TRACE_LOGGING_LEVEL,logKey,'Tracer2 = kinetic energy.')
		if ( AMR > 0 ) then
			if ( useRelativeTOl > 0 ) call LogMessage(exeLog,TRACE_LOGGING_LEVEL,logKey,"AMR used relative tolerances.")
			call LogMessage(exeLog,TRACE_LOGGING_LEVEL,logKey//'panelKind = ',panelKind)
			call LogMessage(exeLog,TRACE_LOGGING_LEVEL,logkey//'initNest = ',initNest)
			call LogMessage(exeLog,TRACE_LOGGING_LEVEL,logKey//'circTol = ',circTol)
			call LogMessage(exeLog,TRACE_LOGGING_LEVEL,logkey//'varTol = ',varTol)
		endif
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
