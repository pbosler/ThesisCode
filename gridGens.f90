program gridGens

use NumberKindsModule
use VoronoiPanelsModule
use LoggerModule


implicit none

! Grid variables
type(VorPanels), pointer :: panels
integer(kint) :: initNest, AMR, nTracer
character(len=128) :: delTriFilename, vorGridFilename
integer(kint) :: i,j

! Logger variables
type(Logger) :: exeLog
integer(kint) :: logOut
character(len=28) :: logKey

! User input variables
integer(kint), parameter :: readUnit = 12
integer(kint) :: readStat, maxNest
namelist /gridInit/ initNest, AMR, nTracer
namelist /fileIO/ delTriFilename, vorGridFilename

! output variables
character(len=128) :: delTriOut, vorGridOut

! INITIALIZE Logger
logOut = 6
call New(exeLog,DEBUG_LOGGING_LEVEL,logOut)
logKey = "EXE log :"
call LogMessage(exeLog,DEBUG_LOGGING_LEVEL,logKey,"initialized.")

! READ User Input
open(unit=readUnit,file='gridGens.namelist',action='READ',status='OLD',iostat=readStat)
	if ( readStat /= 0 ) then
		call LogMessage(exeLog,ERROR_LOGGING_LEVEL,logKey,'ERROR : cannot open namelist file.')
	endif
	read(readUnit,nml=gridInit)
	rewind(readUnit)
	read(readUnit,nml=fileIO)
	call LogMessage(exeLog,TRACE_LOGGING_LEVEL,trim(logkey)//' Building Icos grid to nest ',initNest)
	call LogMessage(exeLog,TRACE_LOGGING_LEVEL,trim(logkey)//' Delaunay triangulation output file : ',trim(delTriFilename))
	call LogMessage(exeLog,TRACE_LOGGING_LEVEL,trim(logkey)//' Voronoi grid output file : ',trim(vorGridFilename))
close(readUnit)
maxNest = initNest

do i=0,maxNest
	allocate(panels)
	write(delTriOut,'(A,A,I1,A)') trim(delTriFilename),'_nest',i,'.vtk'
	write(vorGridOut,'(A,A,I1,A)') trim(vorGridFilename),'_nest',i,'.vtk'
	! BUILD New Grid
	call New(panels,initNest,AMR,nTracer)
	call LogMessage(exeLog,TRACE_LOGGING_LEVEL,logkey,' New Vor Grid returned.')
	do j=1,panels%N
		panels%relVort(j) = (-1.0_kreal)**j*real(mod(j,39),kreal)
	enddo
	! OUTPUT Grid Data
	call VTKOutputDelaunayTriangulation(panels,delTriOut)
	call LogMessage(exeLog,TRACE_LOGGING_LEVEL,trim(logKey)//' VTK Output : Del Tri ','Returned.')
	call VTKOutputVoronoiGrid(panels,vorGridOut)
	call LogMessage(exeLog,TRACE_LOGGING_LEVEL,trim(logKey)//' VTK Output : Vor. grid ','Returned.')
	! DELETE / FREE MEMORY
	call Delete(panels)
	deallocate(panels)
	call LogMessage(exeLog,TRACE_LOGGING_LEVEL,trim(logKey)//'Delete VorGrid',' returned.')
enddo


call Delete(exeLog)

end program
