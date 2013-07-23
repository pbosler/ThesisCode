program RHWavePanels

use NumberKindsModule
use SphereGeomModule
use LoggerModule
use ParticlesEdgesPanelsModule
use PanelVelocitySerialModule
use TracerAndVorticityDistributionModule
use RefineRemeshModule

implicit none

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
character(len=28) :: logKey = 'EXE_LOG : '
character(len=128) :: logString
real(kreal) :: programStart, programFinish
integer(kint) :: procRank = -1

! Timestepping variables
real(kreal) :: t, dt, tfinal
integer(kint) :: timeJ, timesteps

! I/O & User variables
character(len=128) :: jobPrefix, outputDir
character(len=256) :: vtkRoot, vtkFile, dataFile
integer(kint), parameter :: readUnit = 12, writeUnit = 13
integer(kint) :: readStat, writeStat
namelist /gridInit/ panelKind, initnest, AMR, remeshInterval, amrTol1, amrTol2
namelist /time/ dt, tfinal
namelist /fileIO/ jobPrefix, outputDir

! General variables
integer(kint) :: j

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!	Part 1 : Initialize the computing environment / get user input  !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

call New(exeLog,DEBUG_LOGGING_LEVEL,logOut)
call LogMessage(exeLog,TRACE_LOGGING_LEVEL,logKey,"*** Program START ***")

! Read user input from namelist file
open(unit=readUnit,file='RHWavePanels.namelist',action='READ',status='OLD',iostat=readStat)
	if ( readStat /= 0 ) then
		call LogMessage(exeLog,ERROR_LOGGING_LEVEL,logKey,"ERROR : Unable to open namelist file.")
		stop
	endif
	read(readUnit,nml=gridInit)
	rewind(readUnit)
	read(readUnit,nml=time)
	rewind(readUnit)
	read(readUnit,nml=fileIO)
	rewind(readUnit)
close(readUnit)

! Init timestepping variables
t = 0.0_kreal
timeJ = 0
timesteps = floor(tfinal/dt)
! Allocate error vs. time arrays
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

! Prepare output files
if ( panelKind == TRI_PANEL ) then
	if ( AMR == 0) then
		write(dataFile,'(A,A,A,I1,A,F3.1,A,F5.3,A)') trim(outputDir),trim(jobPrefix),'triNest',initNest,'_rev',tfinal,'_dt',dt,'.dat'
		write(vtkRoot,'(A,A,A,A,I1,A,F3.1,A,F5.3,A)') trim(outputDir),'vtkOut/',trim(jobPrefix),'triNest',initNest,'_rev',tfinal,'_dt',dt,'_'
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


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!	Part 2 : Initialize the grid		             !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

call cpu_time(programStart)

nTracer = 3
call New(gridParticles,gridEdges,gridPanels,panelKind,initNest,AMR,nTracer,problemKind)
call InitRH4(gridParticles,gridPanels)
if ( AMR > 0 ) then
	call InitRefine(gridParticles,gridEdges,gridPanels,amrTol1,amrTol2, problemID)
	call InitRH4(gridParticles,gridPanels)
endif

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

call PrintStats(gridParticles)
call PrintStats(gridEdges)
call PrintStats(gridPanels)
call vtkOutput(gridParticles,gridPanels,vtkFile)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!	Part 3 : Run the problem								    !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
remeshFlag = .False.
do timeJ = 0,timesteps-1

	! Remesh if necessary

	if ( mod(timeJ+1,remeshInterval) == 0 ) then
		remeshFlag = .True.
		call LogMessage(exeLog,TRACE_LOGGING_LEVEL,logKey," Remesh triggered by remeshInterval.")
	endif
	if ( remeshFlag ) then
		call AdaptiveRemesh( gridParticles,gridEdges,gridPanels,&
								initNest, AMR, amrTol1, amrTol2, &
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
		call PrintStats(gridParticles)
		call PrintStats(gridEdges)
		call PrintStats(gridPanels)
	endif

	! Advance time
	call BVERK4(gridParticles,gridPanels,dt)
	t = real(timeJ+1,kreal)*dt
	! Calculate error
	gridPanels%tracer(:,3) = 0.0_kreal
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

	! Write output
	call LogMessage(exeLog,DEBUG_LOGGING_LEVEL,trim(logKey)//" t = ",t)
	write(vtkFile,'(A,I0.4,A)') trim(vtkRoot),timeJ+1,'.vtk'
	call vtkOutput(gridParticles,gridPanels,vtkFile)

enddo
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!	Part 4 : Finish and clear memory
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

! Output Remaining Data
open(unit=writeUnit,file=dataFile,status='REPLACE',action='WRITE',iostat=writeStat)
if ( writeStat /= 0 ) call LogMessage(exeLog,ERROR_LOGGING_LEVEL,logkey,'ERROR opening dataFile.')
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

!
! Clean up & Free memory
!
deallocate(totalKineticEnergy)
deallocate(totalEnstrophy)
deallocate(l1)
deallocate(l2)
deallocate(linf)
call ResetRK4()
call Delete(gridParticles,gridEdges,gridPanels)
call Delete(exeLog)


end program
