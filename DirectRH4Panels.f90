program DirectRHWave

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
integer(kint), parameter :: problemID = RH4_WAVE
real(kreal) :: amrTol1, amrTol2

! Error calculation variables
real(kreal), allocatable :: totalKineticEnergy(:), totalEnstrophy(:),&
							l1(:), l2(:), linf(:)
real(kreal) :: l1Denom, l2Denom, linfDenom

! Logger & Computation management variables
type(Logger) :: exeLog
integer(kint) :: logOut = 6
character(len=28) :: logKey 
character(len=128) :: logString
real(kreal) :: wtime, etime
integer(kint), parameter :: REAL_INIT_BUFFER_SIZE = 4, INT_INIT_BUFFER_SIZE = 4
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
namelist /gridInit/ panelKind, initnest, AMR, remeshInterval, amrTol1, amrTol2
namelist /time/ dt, tfinal
namelist /fileIO/ jobPrefix

! General variables
integer(kint) :: j


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
	open(unit=readUnit,file='RHWavePanels.namelist',action='READ',status='OLD',iostat=readStat)
		if (readstat /= 0 ) then
			call LogMessage(exeLog,ERROR_LOGGING_LEVEL,logKey," ERROR opening namelist file.")
			stop
		endif
		read(readunit,nml=gridinit)
		rewind(readunit)
		read(readunit,nml=time)
		rewind(readunit)
		read(readunit,nml=fileIO)
	close(readunit)	
	
	intBuffer(1) = panelKind
	intBuffer(2) = initNest
	intBuffer(3) = AMR
	intBuffer(4) = remeshInterval
	
	realBuffer(1) = amrTol1
	realBuffer(2) = amrTol2
	realBuffer(3) = dt
	realBuffer(4) = tfinal
	
	! Prepare output files
	if ( panelKind == TRI_PANEL ) then
		if ( AMR == 0) then
			write(dataFile,'(A,A,I1,A,F3.1,A,F5.3,A)') trim(jobPrefix),'_DIRECT_triNest',initNest,'_rev',tfinal,'_dt',dt,'.dat'
			write(vtkRoot,'(A,A,A,I1,A,F3.1,A,F5.3,A)') 'vtkOut/',trim(jobPrefix),'_DIRECT_triNest',initNest,'_rev',tfinal,'_dt',dt,'_'
		else
			write(dataFile,'(A,A,I1,A,F3.1,A,F5.3,A)') trim(jobPrefix),'_DIRECT_triAMR',initNest,'_rev',tfinal,'_dt',dt,'.dat'
			write(vtkRoot,'(A,A,A,I1,A,F3.1,A,F5.3,A)') 'vtkOut/',trim(jobPrefix),'_DIRECT_triAMR',initNest,'_rev',tfinal,'_dt',dt,'_'
		endif
	elseif (panelKind == QUAD_PANEL ) then
		if ( AMR == 0) then
			write(dataFile,'(A,A,I1,A,F3.1,A,F5.3,A)') trim(jobPrefix),'DIRECT_quadNest',initNest,'_rev',tfinal,'_dt',dt,'.dat'
			write(vtkRoot,'(A,A,A,I1,A,F3.1,A,F5.3,A)') 'vtkOut/',trim(jobPrefix),'DIRECT_quadNest',initNest,'_rev',tfinal,'_dt',dt,'_'
		else
			write(dataFile,'(A,A,I1,A,F3.1,A,F5.3,A)') trim(jobPrefix),'DIRECT_quadAMR',initNest,'_rev',tfinal,'_dt',dt,'.dat'
			write(vtkRoot,'(A,A,A,I1,A,F3.1,A,F5.3,A)') 'vtkOut/',trim(jobPrefix),'DIRECT_quadAMR',initNest,'_rev',tfinal,'_dt',dt,'_'
		endif
	else
		call LogMessage(exeLog,ERROR_LOGGING_LEVEL,logkey,' ERROR : Invalid panelKind.')
		stop
	endif
	write(vtkFile,'(A,I0.4,A)') trim(vtkRoot),timeJ,'.vtk'
endif
call MPI_BCAST(intBuffer,INT_INIT_BUFFER_SIZE,MPI_INTEGER,0,MPI_COMM_WORLD,mpiErrCode)
call MPI_BCAST(realBuffer,REAL_INIT_BUFFER_SIZE,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiErrCode)

!call LogMessage(exeLog,DEBUG_LOGGING_LEVEL,logKey,'initial broadcast complete.')

panelKind = intBuffer(1)
initNest = intBuffer(2)
AMR = intBuffer(3)
remeshInterval = intBuffer(4)

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
allocate(l1(0:timesteps))
l1 = 0.0_kreal
allocate(l2(0:timesteps))
l2 = 0.0_kreal
allocate(linf(0:timesteps))
linf = 0.0_kreal


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!	Part 2 : Initialize the grid		             !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!


if (procRank == 0) wtime = MPI_WTIME()

nTracer = 3
call New(gridParticles,gridEdges,gridPanels,panelKind,initNest,AMR,nTracer,problemKind)
call InitRH4Wave(gridParticles,gridPanels)
if ( AMR > 0 ) then
	call InitRefine(gridParticles,gridEdges,gridPanels,amrTol1,amrTol2, problemID,procRank=procRank)
	call InitRH4Wave(gridParticles,gridPanels)
endif

!call LogMessage(exeLog,DEBUG_LOGGING_LEVEL,logKey,'initial grid returned.')

! Store initial latitude in tracer 1
do j=1,gridParticles%N
	gridParticles%tracer(j,1) = Latitude(gridParticles%x0(:,j))
enddo
do j=1,gridPanels%N
	if ( .NOT. gridPanels%hasChildren(j) )  then
		gridpanels%tracer(j,1) = Latitude(gridPanels%x0(:,j))
	endif
enddo

totalEnstrophy(0) = 0.5_kreal*sum(gridPanels%area(1:gridpanels%N)*&
					gridPanels%relVort(1:gridPanels%N)*gridPanels%relVort(1:gridPanels%N))
l1Denom = sum(abs(gridPanels%relVort(1:gridPanels%N))*gridPanels%area(1:gridPanels%N))
l2Denom = 2.0_kreal*totalEnstrophy(0)					
linfDenom = maxval(abs(gridPanels%relVort(1:gridPanels%N)))

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
		remeshFlag = .False.
		call ResetRK4()
!		do j=1,gridParticles%N
!			gridParticles%tracer(j,1) = Latitude(gridParticles%x0(:,j))
!		enddo
!		do j=1,gridPanels%N
!			if ( .NOT. gridPanels%hasChildren(j) )  then
!				gridpanels%tracer(j,1) = Latitude(gridPanels%x0(:,j))
!			endif
!		enddo
		if ( procRank == 0 ) then
			print *, "PROC0 End remesh"
			call PrintStats(gridParticles)
			call PrintStats(gridEdges)
			call PrintStats(gridPanels)
		endif
	endif
!	call LogMessage(exeLog,DEBUG_LOGGING_LEVEL,logKey,'REMESH DONE.')
!	call MPI_BARRIER(MPI_COMM_WORLD,mpiErrCode)
	
	if ( procRank == 0 .AND. timeJ == 1) etime = MPI_WTIME()
	
	! Advance time
	call BVERK4(gridParticles,gridPanels,dt,procRank,numProcs)
	t = real(timeJ+1,kreal)*dt
	!print *, "DEBUG t = ",t
	totalKineticEnergy(timeJ+1) = totalKE
	totalEnstrophy(timeJ+1) = 0.5_kreal*sum(gridPanels%area(1:gridPanels%N)*&
			gridPanels%relVort(1:gridPanels%N)*gridPanels%relVort(1:gridPanels%N))
	!print *, "DEBUG Break 2"
	if ( (procRank == 0) .AND. (timeJ == 1)) then
		etime = MPI_WTIME() - etime
		write(logString,'(A,F9.3,A)') 'Estimated time left = ',(timesteps-1)*etime/60.0,' minutes.'
		call LogMessage(exeLog,TRACE_LOGGING_LEVEL,logKey,logString)
	endif
	
	! Calculate error
	gridPanels%tracer(:,3) = 0.0_kreal
	gridParticles%tracer(:,3) = 0.0_kreal
	!print *,"DEBUG Break 3"
	do j=1,gridParticles%N
		gridParticles%tracer(j,3) = abs(gridParticles%relVort(j) - HaurwitzStationary4RelVort(gridParticles%x(:,j)))
	enddo
	gridParticles%tracer(:,3) = gridParticles%tracer(:,3)/linfDenom
	!print *,"DEBUG Break 4"
	do j=1,gridPanels%N
		if ( .NOT. gridPanels%hasChildren(j) ) then
			gridPanels%tracer(j,3) = abs(gridPanels%relVort(j) - &
						HaurwitzStationary4RelVort(gridPanels%x(:,j)))
		endif
	enddo
	
	l1(timeJ+1) = sum(gridPanels%tracer(1:gridPanels%N,3)*gridPanels%area(1:gridPanels%N))/l1Denom
	l2(timeJ+1) = sum(gridPanels%tracer(1:gridPanels%N,3)*gridPanels%tracer(1:gridPanels%N,3)*&
				gridPanels%area(1:gridPanels%N))/l2Denom
	linf(timeJ+1) = maxVal(gridPanels%tracer(1:gridPanels%N,3))/linfDenom
	gridPanels%tracer(1:gridPanels%N,3) = gridPanels%tracer(1:gridPanels%N,3)/linfDenom
	
	if ( procRank == 0 ) then
		call LogMessage(exeLog,DEBUG_LOGGING_LEVEL,trim(logKey)//" t = ",t)
		write(vtkFile,'(A,I0.4,A)') trim(vtkRoot),timeJ+1,'.vtk'
		call vtkOutput(gridParticles,gridPanels,vtkFile)
	endif
enddo	

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!	Part 4 : Finish and clear memory
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

if ( procRank == 0 ) then
	open(unit=writeUnit,file=dataFile,status='REPLACE',action='WRITE',iostat=writeStat)
	if ( writeStat /= 0 ) call LogMessage(exeLog,ERROR_LOGGING_LEVEL,logkey,'ERROR opening dataFile.')
		write(writeUnit,'(5A24)') 'totalKE','totalEnstrophy','relVort l1','relVort l2','relVort linf'
		do j=0,timesteps
			write(writeUnit,'(5F24.15)') totalKineticEnergy(j), totalEnstrophy(j), &
				l1(j), l2(j), linf(j)
		enddo
	close(writeUnit)
	
	wTime = MPI_WTIME() - wtime
	call LogMessage(exeLog,TRACE_LOGGING_LEVEL,logKey,"*** Program END ***")
	call LogMessage(exeLog,TRACE_LOGGING_LEVEL,trim(logKey)//' data file = ',trim(dataFile))
	write(vtkFile,'(A,A)') trim(vtkRoot),'XXXX.vtk'
	call LogMessage(exeLog,TRACE_LOGGING_LEVEL,trim(logKey)//' vtkFiles = ',trim(vtkFile))
	call LogMessage(exeLog,TRACE_LOGGING_LEVEL,logKey,'Tracer1 = initial latitude.')
	call LogMessage(exeLog,TRACE_LOGGING_LEVEL,logKey,'Tracer2 = kinetic energy.')
	call LogMessage(exeLog,TRACE_LOGGING_LEVEL,logKey,'Tracer3 = vorticity error.')
	!call EndSection(exeLog%writer)
	call StartSection(exeLog, ' ------ Calculation Summary ----------- ')
		if (panelKind == QUAD_PANEL .AND. AMR == 0) then
			write(logString,'(A,I8)') 'Cubed sphere nActive = ',gridPanels%N_Active
		elseif ( panelKind == TRI_PANEL .AND. AMR == 0) then
			write(logString,'(A,I8)') 'Icos. triangles nActive = ',gridPanels%N_Active
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
deallocate(l1)
deallocate(l2)
deallocate(linf)
call FinalizeMPIRK4()
call Delete(gridParticles,gridEdges,gridPanels)
call Delete(exeLog)

call MPI_FINALIZE(mpiErrCode)

end program
