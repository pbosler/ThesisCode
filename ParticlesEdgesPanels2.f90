module ParticlesEdgesPanelsModule
!	Pete Bosler
!   University of Michigan, Department of Mathematics
!   July 18, 2012.
use NumberKindsModule
use SphereGeomModule
use IntegerListModule
use LoggerModule

implicit none

!----------------
! USAGE :  This module provides the basic data structures for particle-panel methods on the 
! unit sphere.  Memory is preallocated statically for speed, assuming the column-priority indexing
! typical to Fortran.  This means that information associated with panel j will be contained in
! column j of the panels' structure's various arrays.  
!		For example : the physical position of panel j is stored at panels%x(:,j), and the area
!				associated with panel j is stored at panels%area(j).  Vertices that make up the corners
!				of panel j are stored at panels%vertices(:,j), etc...  The same is true for particles and
!				edges, so that the indices to the particles that make up edge k are stored at
!				edges%verts(:,k).  
!		Memory is handled with the "New" and "Delete" interfaces and procedures; grids are initialized
! 				with subroutines labeled "Init*".  
!----------------
private
public Particles, Edges, Panels
public New, Delete, Copy
public QUAD_PANEL, TRI_PANEL, VORONOI_PANEL, ADVECTION_SOLVER, BVE_SOLVER
public GetNTracer, GetPanelKind, GetProblemKind
public vtkOutput, PrintStats
public GatherPanels, ScatterPanels, ResetGridArea, RenormalizeGrid
public FindRootPanel, AdjacentPanelSearch, CCWDualVertList
public DividePanel
public TriAvgMeshSize, QuadAvgMeshSize

integer(kint), parameter :: QUAD_PANEL = 4, &
							TRI_PANEL = 3, &
							VORONOI_PANEL = 10, &
							ADVECTION_SOLVER = 60, &
							BVE_SOLVER = 61
!----------------
! Types
!----------------

type Particles
	! Grid variables (used always)
	real(kreal), dimension(:,:), pointer :: x		! Cartesian coordinate vectors
	real(kreal), dimension(:,:), pointer :: x0      ! Cartesian coordinate vectors for Lagrange coordinate
	integer(kint), dimension(:,:), pointer :: edges ! Edges that begin or end at this particle
	integer(kint) :: N								! Current number of particles
	integer(kint) :: N_Max							! Max number of particles allowed in memory
	real(kreal), dimension(:), pointer :: dualArea  ! Area of dual mesh panel centered at particle
	! Data variables (context dependent)
	real(kreal), dimension(:,:), pointer :: tracer  ! Passive tracer field
	real(kreal), dimension(:), pointer :: absVort   ! absolute vorticity
	real(kreal), dimension(:), pointer :: relVort   ! relative vorticity
end type

type Edges
	! Grid variables (used always)
	integer(kint), dimension(:,:), pointer :: verts 	! 2-dim arrays. Edge points from verts(1) to verts(2)
	integer(kint), dimension(:), pointer :: leftPanel 	! panel to left of directed edge
	integer(kint), dimension(:), pointer :: rightPanel  ! panel to right of directed edge
	logical(klog), dimension(:), pointer :: hasChildren ! true if edge has been divided
	integer(kint), dimension(:,:), pointer :: children  ! indices to child edges
	integer(kint) :: N									! current number of edges
	integer(kint) :: N_Max								! max number of edges allowed in memory
	real(kreal), dimension(:), pointer :: length	    ! current length of edge
	real(kreal), dimension(:), pointer :: length0       ! initial edge length
end type

type Panels
	! Grid variables (used always)
	real(kreal), dimension(:,:), pointer :: x			! Cartesian coordinate vectors
	real(kreal), dimension(:,:), pointer :: x0			! Cartesian coordinate vectors for Lagrange coordinate
	real(kreal), dimension(:), pointer ::  area			! area of low-level panels (divided panels have area = 0.0)
	integer(kint) :: N									! current number of panels in memory
	integer(kint) :: N_Active							! current number of low-level panels
	integer(kint) :: N_Max								! max number of panels allowed in memory
	integer(kint), dimension(:,:), pointer :: vertices  ! indices to particles at panel corners
	integer(kint), dimension(:,:), pointer :: edges     ! indices to edges of panel
	logical(klog), dimension(:), pointer :: hasChildren ! true if panel has been divided
	integer(kint), dimension(:,:), pointer :: children  ! indices to child panels
	integer(kint), dimension(:), pointer :: nest		! nest level of panel
	! Data variables (context dependent)
	real(kreal), dimension(:,:), pointer :: tracer		! Passive tracer field
	real(kreal), dimension(:), pointer :: absVort		! absolute vorticity
	real(kreal), dimension(:), pointer :: relVort		! relative vorticity
end type

logical(klog), save :: logInit = .FALSE.
type(Logger) :: log
character(len=28), save :: logkey = 'PPE2 :'
integer(kint), parameter :: logLevel = TRACE_LOGGING_LEVEL

!----------------
! Interfaces
!----------------

interface New
	module procedure NewGrid
	module procedure NewPrivatePanels
	module procedure NewPrivateParticles
	module procedure NewPrivateEdges
end interface

interface Delete
	module procedure DeleteGrid
	module procedure DeletePrivatePanels
	module procedure DeletePrivateParticles
end interface

interface GetNTracer
	module procedure NTracerParticles
	module procedure NTracerPanels
end interface

interface GetProblemKind
	module procedure ParticlesProblemKind
	module procedure PanelsProblemKind
end interface

interface PrintStats
    module procedure PrintParticlesStats
    module procedure PrintPanelsStats
    module procedure PrintEdgesStats
end interface

interface Copy
	module procedure CopyPanelByIndex
end interface

interface vtkOutput
	module procedure vtkOutputUnif
	module procedure vtkOutputAMR
	module procedure vtkOutputActivePassive
end interface

contains

subroutine InitLogger(aLog)
	type(Logger), intent(inout) :: aLog
	integer(kint) :: logUnit
	logUnit = 6
	call New(aLog,logLevel,logUnit)
	logInit = .TRUE.
end subroutine


!----------------
! Standard Methods : Construction, destruction, copy
!----------------
subroutine NewPrivateParticles(aParticles, nMax, panelKind, nTracer, problemKind)
! Allocates memory for the Particles data structure, and sets initial
! state to zero/null.
	type(Particles), intent(out) :: aParticles
	integer(kint), intent(in) :: nMax, panelKind, nTracer, problemKind
	if (nMax <= 0) stop "Particles: Invalid nMax."
	if (  (problemKind /= ADVECTION_SOLVER) .AND. (problemKind /= BVE_SOLVER) ) then
		stop "NewParticles ERROR : Invalid problemKind."
	else
		if ( ((panelKind /= TRI_PANEL) .AND. (panelKind /= QUAD_PANEL )) .AND. (panelKind /= VORONOI_PANEL) ) then
			stop "NewParticles ERROR : Invalid panelKind."
		else
			allocate(aParticles%x(3,nMax))
				aParticles%x = 0.0_kreal
			allocate(aParticles%x0(3,nMax))
				aParticles%x0 = 0.0_kreal
			if ( panelKind == TRI_PANEL) then
				allocate(aParticles%edges(6,nMax))
			elseif (panelKind == QUAD_PANEL ) then
				allocate(aParticles%edges(4,nMax))
			elseif (panelKind == VORONOI_PANEL) then
				allocate(aParticles%edges(10,nMax))
			endif
			aParticles%edges = 0
			allocate(aParticles%dualArea(nMax))
			aParticles%dualArea = 0.0_kreal
			if (nTracer > 0) then
				allocate(aParticles%tracer(nMax,nTracer))
				aParticles%tracer = 0.0_kreal
			else
				nullify(aParticles%tracer)
			endif
			if (problemKind == ADVECTION_SOLVER) then
				nullify(aParticles%absVort)
				nullify(aParticles%relVort)
			elseif (problemKind == BVE_SOLVER) then
				allocate(aParticles%absVort(nMax))
				aParticles%absVort = 0.0_kreal
				allocate(aParticles%relVort(nMax))
				aParticles%relVort = 0.0_kreal
			endif
			aParticles%N = 0
			aParticles%N_Max = nMax
		endif
	endif
end subroutine


subroutine DeletePrivateParticles(aParticles)
! Frees memory associated with (deletes) a Particles data type.
	type(Particles), intent(inout) :: aParticles
	deallocate(aParticles%x)
	deallocate(aParticles%x0)
	deallocate(aParticles%edges)
	deallocate(aParticles%dualArea)
	if ( associated(aParticles%tracer)) deallocate(aParticles%tracer)
	aParticles%N = 0
	aParticles%N_Max = 0
	if (associated(aParticles%absVort)) deallocate(aParticles%absVort)
	if (associated(aParticles%relVort)) deallocate(aParticles%relVort)
end subroutine


subroutine NewPrivateEdges(anEdges, nMax )
! Allocates memory for the Edges data structure, and sets initial
! state to zero/null.
	type(Edges), intent(out) :: anEdges
	integer(kint), intent(in) :: nMax
	if (nMax <= 0) stop "NewEdges : Invalid nMax."
	allocate(anEdges%verts(2,nMax))
	anEdges%verts = 0
	allocate(anEdges%leftPanel(nMax))
	anEdges%leftPanel = 0
	allocate(anEdges%rightPanel(nMax))
	anEdges%rightPanel = 0
	allocate(anEdges%hasChildren(nMax))
	anEdges%hasChildren = .False.
	allocate(anEdges%children(2,nMax))
	anEdges%children = 0
	anEdges%N = 0
	anEdges%N_Max = nMax
	allocate(anEdges%length(nMax))
	anEdges%length = 0.0_kreal
	allocate(anEdges%length0(nMax))
	anEdges%length0 = 0.0_kreal
end subroutine


subroutine DeletePrivateEdges(anEdges)
! Frees memory asssociated with an Edges data structure
	type(Edges), intent(inout) :: anEdges
	deallocate(anEdges%verts)
	deallocate(anEdges%leftPanel)
	deallocate(anEdges%rightPanel)
	deallocate(anEdges%length)
	deallocate(anEdges%length0)
	deallocate(anEdges%hasChildren)
	deallocate(anEdges%children)
	anEdges%N = 0
	anEdges%N_Max = 0
end subroutine


subroutine NewPrivatePanels(aPanels, nMax, panelKind, nTracer, problemKind)
! Allocates memory for the Panels data structure, and sets initial
! state to zero/null.
	type(Panels), intent(out) :: aPanels
	integer(kint), intent(in) :: nMax, panelKind, nTracer, problemKind
	
	if (nMax <= 0 ) stop "Panels : Invalid nMax."
	if ( ((panelKind == TRI_PANEL) .OR. (panelKind == QUAD_PANEL)) ) then
		if ( (problemKind == ADVECTION_SOLVER) .OR. (problemKind == BVE_SOLVER)) then
			allocate(aPanels%x(3,nMax))
			aPanels%x = 0.0_kreal
			allocate(aPanels%x0(3,nMax))
			aPanels%x0 = 0.0_kreal
			allocate(aPanels%area(nMax))
			aPanels%area = 0.0_kreal
			if (panelKind == TRI_PANEL) then
				allocate(aPanels%vertices(3,nMax))
				allocate(aPanels%edges(3,nMax))
			elseif (panelKind == QUAD_PANEL ) then
				allocate(aPanels%vertices(4,nMax))
				allocate(aPanels%edges(4,nMax))
			elseif (panelKind == VORONOI_PANEL) then
				allocate(aPanels%vertices(10,nMax))
				allocate(aPanels%edges(10,nMax))
			endif
			aPanels%vertices = 0
			aPanels%edges = 0
			allocate(aPanels%hasChildren(nMax))
			aPanels%hasChildren = .False.
			allocate(aPanels%children(4,nMax))
			aPanels%children = 0
			allocate(aPanels%nest(nMax))
			aPanels%nest = -1
			aPanels%N = 0
			aPanels%N_Active = 0
			aPanels%N_Max = nMax
			if (nTracer > 0) then
				allocate(aPanels%tracer(nMax,nTracer))
				aPanels%tracer = 0.0_kreal
			else
				nullify(aPanels%tracer)
			endif
			if (problemKind == ADVECTION_SOLVER) then
				nullify(aPanels%relVort)
				nullify(aPanels%absVort)
			elseif (problemKind == BVE_SOLVER) then
				allocate(aPanels%relVort(nMax))
				allocate(aPanels%absVort(nMax))
				aPanels%relVort = 0.0_kreal
				aPanels%absVort = 0.0_kreal
			endif
		else
			stop "NewPanels ERROR : Invalid problemKind."
		endif
	else
		stop "NewPanels ERROR : Invalid panelKind."
	endif
end subroutine


subroutine DeletePrivatePanels(aPanels)
! Frees memory (deletes) a Panels data structure
	type(Panels), intent(inout) :: aPanels
	deallocate(aPanels%x)
	deallocate(aPanels%x0)
	deallocate(aPanels%vertices)
	deallocate(aPanels%area)
	deallocate(aPanels%edges)
	deallocate(aPanels%hasChildren)
	deallocate(aPanels%children)
	deallocate(aPanels%nest)
	if (associated(aPanels%tracer)) deallocate(aPanels%tracer)
	if (associated(aPanels%relVort)) deallocate(aPanels%relVort)
	if (associated(aPanels%absVort)) deallocate(aPanels%absVort)
	aPanels%N = 0
	aPanels%N_Max = 0
	aPanels%N_Active = 0
end subroutine


subroutine CopyPanelByIndex(newPanels, newIndex, oldPanels, oldIndex)
! Copies the information for the panel and oldIndex of the oldPanels
! data structure to the newIndex location of the newPanels data structure.
	type(Panels), intent(inout) :: newPanels
	type(Panels), intent(in) :: oldPanels
	integer(kint), intent(in) :: newIndex, oldIndex
	newPanels%x(:,newIndex) = oldPanels%x(:,oldIndex)
	newPanels%x0(:,newIndex) = oldPanels%x0(:,oldIndex)
	newPanels%area(newIndex) = oldPanels%area(oldIndex)
	newPanels%vertices(:,newIndex) = oldPanels%vertices(:,oldIndex)
	newPanels%edges(:,newIndex) = oldPanels%edges(:,oldIndex)
	newPanels%hasChildren(newIndex) = oldPanels%hasChildren(oldIndex)
	newPanels%children(:,newIndex) = oldPanels%children(:,oldIndex)
	newPanels%nest(newIndex) = oldPanels%nest(oldIndex)
	if (associated(oldPanels%tracer)) then	
		newPanels%tracer(newIndex,:) = oldPanels%tracer(oldIndex,:)
	endif
	if (associated(oldPanels%absVort)) then
		newPanels%absVort(newIndex) = oldPanels%absVort(oldIndex)
	endif
	if (associated(oldPanels%relVort)) then
		newPanels%relVort(newIndex) = oldPanels%relVort(oldIndex)
	endif
end subroutine


subroutine NewGrid(aParticles, anEdges, aPanels, panelKind, initNest, AMR, nTracer, problemKind)
! This is the primary user interface function for grid generation.  A grid is made of of 
! particles, edges, and panels.  This routine allocates memory for a grid and calls the appropriate
! initialization subroutine based on the input parameters.
	type(Particles), pointer, intent(out) :: aParticles
	type(Edges), pointer, intent(out) :: anEdges
	type(Panels), pointer, intent(out) :: aPanels
	integer(kint), intent(in) :: panelKind, initNest, AMR, nTracer, problemKind
	integer(kint) :: nPanels, nParticles, nEdges
	
	if ( .NOT. logInit ) then
		call InitLogger(Log)
	endif
	call LogMessage(log,DEBUG_LOGGING_LEVEL,logKey," entering NewGrid.")
	
	
	! allocate the pointers
	allocate(aParticles)
	allocate(anEdges)
	allocate(aPanels)
	! compute the array sizes
	nPanels = PanelMax(panelKind,initNest + AMR)
	nParticles = ParticleMax(panelKind,initNest + AMR)
	nEdges = EdgeMax(panelKind,initNest+AMR)
	! allocate the data structures
	call NewPrivateParticles(aParticles,nParticles,panelKind,nTracer,problemKind)
	call NewPrivateEdges(anEdges,nEdges)
	call NewPrivatePanels(aPanels,nPanels,panelKind,nTracer,problemKind)
	! Initialize the grid
	if (panelKind == TRI_PANEL ) then
		call InitIcosTriGrid(aParticles,anEdges,aPanels,initNest)
	elseif (panelKind == QUAD_PANEL ) then
		call InitCubedSphereGrid(aParticles,anEdges,aPanels,initNest)
	endif
end subroutine


subroutine DeleteGrid(aParticles, anEdges, aPanels)
! This is the primary Delete subroutine.  
! Deletes a grid and frees memory with each of its sub-components.
	type(Particles), pointer, intent(inout) :: aParticles
	type(Edges), pointer, intent(inout) :: anEdges
	type(Panels), pointer, intent(inout) :: aPanels
	! Delete the data structure arrays
	call DeletePrivateParticles(aParticles)
	call DeletePrivateEdges(anEdges)
	call DeletePrivatePanels(aPanels)
	! Deallocate the pointers
	deallocate(aParticles)
	deallocate(anEdges)
	deallocate(aPanels)
end subroutine

!----------------
! Grid Methods : 
!----------------

subroutine InitIcosTriGrid(aParticles, anEdges, aPanels, initNest)
! Initializes an icosahedral triangle grid.  It assumes memory has been
! allocated properly for the particles, edges, and panels data structures previously.
	type(Particles), intent(inout) :: aParticles
	type(Edges), intent(inout) :: anEdges
	type(Panels), intent(inout) :: aPanels
	integer(kint), intent(in) :: initNest
	! local variables
	integer(kint) :: j,k, startIndex, nOldPanels, panelVerts(3)
	!
	! set up the icosahedron
	!
	call LogMessage(log,DEBUG_LOGGING_LEVEL,logKey," initializing icosTri grid.")
	! Particles
	aParticles%N = 12
	aParticles%x(:,1) = [0.0_kreal, 0.0_kreal, 1.0_kreal]
	aParticles%x(:,2) = [0.723606797749978969640917366873_kreal,&
						 0.525731112119133606025669084848_kreal,&
						 0.447213595499957939281834733746_kreal]
	aParticles%x(:,3) =[-0.276393202250021030359082633126_kreal,&
						 0.850650808352039932181540497063_kreal,&
						 0.447213595499957939281834733746_kreal]
	aParticles%x(:,4) =[-0.894427190999915878563669467492_kreal,&
						 0.0_kreal,&
						 0.447213595499957939281834733746_kreal]
	aParticles%x(:,5) =[-0.276393202250021030359082633127_kreal,&
						-0.850650808352039932181540497063_kreal,&
						 0.447213595499957939281834733746_kreal]
	aParticles%x(:,6) = [0.723606797749978969640917366873_kreal,&
						-0.525731112119133606025669084848_kreal,&
						 0.447213595499957939281834733746_kreal]
	aParticles%x(:,7) = [0.894427190999915878563669467492_kreal,&
						 0.0_kreal,&
						-0.447213595499957939281834733746_kreal]
	aParticles%x(:,8) = [0.276393202250021030359082633127_kreal,&
						 0.850650808352039932181540497063_kreal,&
						-0.447213595499957939281834733746_kreal]
	aParticles%x(:,9) =[-0.723606797749978969640917366873_kreal,&
						 0.525731112119133606025669084848_kreal,&
						-0.447213595499957939281834733746_kreal]
	aParticles%x(:,10)=[-0.723606797749978969640917366873_kreal,&
						-0.525731112119133606025669084848_kreal,&
						-0.447213595499957939281834733746_kreal]
	aParticles%x(:,11)= [0.276393202250021030359082633127_kreal,&
						-0.850650808352039932181540497063_kreal,&
						-0.447213595499957939281834733746_kreal]
	do j=1,11
		aParticles%x(:,j) = aParticles%x(:,j)/sqrt(sum(aParticles%x(:,j)*aParticles%x(:,j)))
	enddo					
	aParticles%x(:,12)= [0.0_kreal,0.0_kreal,-1.0_kreal]
	aParticles%edges(:,1) = [1,3,5,7,9,0]
	aParticles%edges(:,2) = [1,10,24,11,2,0]
	aParticles%edges(:,3) = [3,2,12,14,4,0]
	aParticles%edges(:,4) = [5,4,15,17,6,0]
	aParticles%edges(:,5) = [7,6,18,20,8,0]
	aParticles%edges(:,6) = [9,8,21,23,10,0]
	aParticles%edges(:,7) = [23,22,30,25,24,0]
	aParticles%edges(:,8) = [11,25,26,13,12,0]
	aParticles%edges(:,9) = [14,13,27,16,15,0]
	aParticles%edges(:,10)= [17,16,28,19,18,0]
	aParticles%edges(:,11)= [20,19,29,22,21,0]
	aParticles%edges(:,12)= [26,30,29,28,27,0]
	aParticles%x0(:,1:12) = aParticles%x(:,1:12)
	
	
	! EDGES
	
	anEdges%N = 30
	anEdges%verts(:,1) = [1,2]
	anEdges%verts(:,2) = [2,3]
	anEdges%verts(:,3) = [3,1]
	anEdges%verts(:,4) = [3,4]
	anEdges%verts(:,5) = [4,1]
	anEdges%verts(:,6) = [4,5]
	anEdges%verts(:,7) = [5,1]
	anEdges%verts(:,8) = [5,6]
	anEdges%verts(:,9) = [6,1]
	anEdges%verts(:,10)= [6,2]
	anEdges%verts(:,11)= [2,8]
	anEdges%verts(:,12)= [8,3]
	anEdges%verts(:,13)= [8,9]
	anEdges%verts(:,14)= [9,3]
	anEdges%verts(:,15)= [9,4]
	anEdges%verts(:,16)= [9,10]
	anEdges%verts(:,17)= [10,4]
	anEdges%verts(:,18)= [10,5]
	anEdges%verts(:,19)= [10,11]
	anEdges%verts(:,20)= [11,5]
	anEdges%verts(:,21)= [11,6]
	anEdges%verts(:,22)= [11,7]
	anEdges%verts(:,23)= [7,6]
	anEdges%verts(:,24)= [7,2]
	anEdges%verts(:,25)= [7,8]
	anEdges%verts(:,26)= [8,12]
	anEdges%verts(:,27)= [12,9]
	anEdges%verts(:,28)= [12,10]
	anEdges%verts(:,29)= [12,11]
	anEdges%verts(:,30)= [12,7]
	anEdges%leftPanel(1) = 1
	anEdges%leftPanel(2) = 1
	anEdges%leftPanel(3) = 1
	anEdges%leftPanel(4) = 2
	anEdges%leftPanel(5) = 2
	anEdges%leftPanel(6) = 3
	anEdges%leftPanel(7) = 3
	anEdges%leftPanel(8) = 4
	anEdges%leftPanel(9) = 4
	anEdges%leftPanel(10)= 5
	anEdges%leftPanel(11)= 6
	anEdges%leftPanel(12)= 6
	anEdges%leftPanel(13)= 7
	anEdges%leftPanel(14)= 7
	anEdges%leftPanel(15)= 8
	anEdges%leftPanel(16)= 9
	anEdges%leftPanel(17)= 9
	anEdges%leftPanel(18)= 10
	anEdges%leftPanel(19)= 11
	anEdges%leftPanel(20)= 11
	anEdges%leftPanel(21)= 12
	anEdges%leftPanel(22)= 13
	anEdges%leftPanel(23)= 13
	anEdges%leftPanel(24)= 14
	anEdges%leftPanel(25)= 15
	anEdges%leftPanel(26)= 16
	anEdges%leftPanel(27)= 16
	anEdges%leftPanel(28)= 17
	anEdges%leftPanel(29)= 18
	anEdges%leftPanel(30)= 19
	anEdges%rightPanel(1) = 5
	anEdges%rightPanel(2) = 6
	anEdges%rightPanel(3) = 2
	anEdges%rightPanel(4) = 8
	anEdges%rightPanel(5) = 3
	anEdges%rightPanel(6) = 10
	anEdges%rightPanel(7) = 4
	anEdges%rightPanel(8) = 12
	anEdges%rightPanel(9) = 5
	anEdges%rightPanel(10)= 14
	anEdges%rightPanel(11)= 15
	anEdges%rightPanel(12)= 7
	anEdges%rightPanel(13)= 16
	anEdges%rightPanel(14)= 8
	anEdges%rightPanel(15)= 9
	anEdges%rightPanel(16)= 17
	anEdges%rightPanel(17)= 10
	anEdges%rightPanel(18)= 11
	anEdges%rightPanel(19)= 18
	anEdges%rightPanel(20)= 12
	anEdges%rightPanel(21)= 13
	anEdges%rightPanel(22)= 19
	anEdges%rightPanel(23)= 14
	anEdges%rightPanel(24)= 15
	anEdges%rightPanel(25)= 20
	anEdges%rightPanel(26)= 20
	anEdges%rightPanel(27)= 17
	anEdges%rightPanel(28)= 18
	anEdges%rightPanel(29)= 19
	anEdges%rightPanel(30)= 20
	
	! PANELS
	
	aPanels%N = 20
	aPanels%N_Active = 20
	aPanels%edges(:,1) = [1,2,3]
	aPanels%edges(:,2) = [3,4,5]
	aPanels%edges(:,3) = [5,6,7]
	aPanels%edges(:,4) = [7,8,9]
	aPanels%edges(:,5) = [9,10,1]
	aPanels%edges(:,6) = [12,2,11]
	aPanels%edges(:,7) = [12,13,14]
	aPanels%edges(:,8) = [15,4,14]
	aPanels%edges(:,9) = [15,16,17]
	aPanels%edges(:,10)= [18,6,17]
	aPanels%edges(:,11)= [18,19,20]
	aPanels%edges(:,12)= [21,8,20]
	aPanels%edges(:,13)= [21,22,23]
	aPanels%edges(:,14)= [24,10,23]
	aPanels%edges(:,15)= [24,25,11]
	aPanels%edges(:,16)= [27,13,26]
	aPanels%edges(:,17)= [28,16,27]
	aPanels%edges(:,18)= [29,19,28]
	aPanels%edges(:,19)= [30,22,29]
	aPanels%edges(:,20)= [26,25,30]
	aPanels%vertices(:,1) = [1,2,3]
	aPanels%vertices(:,2) = [1,3,4]
	aPanels%vertices(:,3) = [1,4,5]
	aPanels%vertices(:,4) = [1,5,6]
	aPanels%vertices(:,5) = [1,6,2]
	aPanels%vertices(:,6) = [8,3,2]
	aPanels%vertices(:,7) = [3,8,9]
	aPanels%vertices(:,8) = [9,4,3]
	aPanels%vertices(:,9) = [4,9,10]
	aPanels%vertices(:,10)= [10,5,4]
	aPanels%vertices(:,11)= [5,10,11]
	aPanels%vertices(:,12)= [11,6,5]
	aPanels%vertices(:,13)= [6,11,7]
	aPanels%vertices(:,14)= [7,2,6]
	aPanels%vertices(:,15)= [2,7,8]
	aPanels%vertices(:,16)= [12,9,8]
	aPanels%vertices(:,17)= [12,10,9]
	aPanels%vertices(:,18)= [12,11,10]
	aPanels%vertices(:,19)= [12,7,11]
	aPanels%vertices(:,20)= [12,8,7]
	! panel centers and areas
	do j=1,20
		panelVerts = aPanels%vertices(:,j)
		aPanels%x(:,j) = SphereTriCenter(aParticles%x(:,panelVerts(1)),aParticles%x(:,panelVerts(2)),aParticles%x(:,panelVerts(3)))
		aPanels%area(j) = TriPanelArea(aParticles,anEdges,aPanels,j)
	enddo
	aPanels%x0(:,1:20) = aPanels%x(:,1:20)
	aPanels%nest(1:20) = 0
	
	! ICOSAHEDRON is complete.  Recursively divide it to the desired resolution
	startIndex = 1
	do k=1,initNest
		nOldPanels = aPanels%N
		do j=startIndex,nOldPanels
			if ( .NOT. aPanels%hasChildren(j)) then
				call DivideTriPanel(aParticles,anEdges,aPanels,j)
			endif
		enddo
		startIndex = nOldPanels
	enddo
	do j=1,aParticles%N
		aParticles%dualArea(j) = DualPanelArea(aParticles,anEdges,aPanels,j)
	enddo	
end subroutine


subroutine InitCubedSphereGrid(aParticles,anEdges,aPanels,initNest)
! Initializes a cubed sphere (quadrilateral) grid for the unit sphere.  
! Memory must have been allocated for the grid data structures prior to calling this routine.
	type(Particles), intent(inout) :: aParticles
	type(Edges), intent(inout) :: anEdges
	type(Panels), intent(inout) :: aPanels
	integer(kint), intent(in) :: initNest
	! local variables
	real(kreal), parameter :: a = 1.0_kreal/dsqrt(3.0_kreal)
	integer(kint) :: j, k, startIndex, nOldPanels
	
	! setup the root cube
	aParticles%N = 8
	aParticles%x(:,1) = [a,-a,a]
	aParticles%x(:,2) = [ a,-a,-a]
	aParticles%x(:,3) = [ a, a,-a]
	aParticles%x(:,4) = [ a, a, a]
	aParticles%x(:,5) = [-a, a,-a]
	aParticles%x(:,6) = [-a, a, a]
	aParticles%x(:,7) = [-a,-a,-a]
	aParticles%x(:,8) = [-a,-a, a]
	aParticles%edges(:,1) = [0,12,1,4]
	aParticles%edges(:,2) = [1,11,0,2]
	aParticles%edges(:,3) = [3,2,0,5]
	aParticles%edges(:,4) = [0,4,3,7]
	aParticles%edges(:,5) = [6,5,0,8]
	aParticles%edges(:,6) = [0,7,6,10]
	aParticles%edges(:,7) = [9,8,0,11]
	aParticles%edges(:,8) = [0,10,9,12]
	aParticles%x0(:,1:8) = aParticles%x(:,1:8)
	
	anEdges%N = 12
	anEdges%verts(:,1) = [1,2]
	anEdges%verts(:,2) = [2,3]
	anEdges%verts(:,3) = [3,4]
	anEdges%verts(:,4) = [4,1]
	anEdges%verts(:,5) = [3,5]
	anEdges%verts(:,6) = [5,6]
	anEdges%verts(:,7) = [6,4]
	anEdges%verts(:,8) = [5,7]
	anEdges%verts(:,9) = [7,8]
	anEdges%verts(:,10)= [8,6]
	anEdges%verts(:,11)= [7,2]
	anEdges%verts(:,12)= [1,8]
	anEdges%leftPanel(1) = 1
	anEdges%leftPanel(2) = 1
	anEdges%leftPanel(3) = 1
	anEdges%leftPanel(4) = 1
	anEdges%leftPanel(5) = 2
	anEdges%leftPanel(6) = 2
	anEdges%leftPanel(7) = 2
	anEdges%leftPanel(8) = 3
	anEdges%leftPanel(9) = 3
	anEdges%leftPanel(10)= 3
	anEdges%leftPanel(11)= 4
	anEdges%leftPanel(12)= 4	
	anEdges%rightPanel(1) = 4
	anEdges%rightPanel(2) = 6
	anEdges%rightPanel(3) = 2
	anEdges%rightPanel(4) = 5
	anEdges%rightPanel(5) = 6
	anEdges%rightPanel(6) = 3
	anEdges%rightPanel(7) = 5
	anEdges%rightPanel(8) = 6
	anEdges%rightPanel(9) = 4
	anEdges%rightPanel(10)= 5
	anEdges%rightPanel(11)= 6
	anEdges%rightPanel(12)= 5

	aPanels%nest(1:6) = 0
	aPanels%N = 6
	aPanels%N_Active = 6
	aPanels%x(:,1) = [1.0_kreal,0.0_kreal,0.0_kreal]
	aPanels%x(:,2) = [0.0_kreal,1.0_kreal,0.0_kreal]
	aPanels%x(:,3) = [-1.0_kreal,0.0_kreal,0.0_kreal]
	aPanels%x(:,4) = [0.0_kreal,-1.0_kreal,0.0_kreal]
	aPanels%x(:,5) = [0.0_kreal,0.0_kreal,1.0_kreal]
	aPanels%x(:,6) = [0.0_kreal,0.0_kreal,-1.0_kreal]
	aPanels%edges(:,1) = [1,2,3,4]
	aPanels%edges(:,2) = [3,5,6,7]
	aPanels%edges(:,3) = [6,8,9,10]
	aPanels%edges(:,4) = [9,11,1,12]
	aPanels%edges(:,5) = [12,4,7,10]
	aPanels%edges(:,6) = [11,8,5,2]
	aPanels%vertices(:,1) = [1,2,3,4]
	aPanels%vertices(:,2) = [4,3,5,6]
	aPanels%vertices(:,3) = [6,5,7,8]
	aPanels%vertices(:,4) = [8,7,2,1]
	aPanels%vertices(:,5) = [8,1,4,6]
	aPanels%vertices(:,6) = [2,7,5,3]
	apanels%x0(:,1:6) = apanels%x(:,1:6)
	do j=1,6
		apanels%area(j) = QuadPanelArea(aParticles,anEdges,aPanels,j)
	enddo
	
	startIndex = 1
	do k=1,initNest
		nOldPanels = aPanels%N
		do j=startIndex,nOldPanels
			if ( .NOT. aPanels%hasChildren(j) ) then
				call DivideQuadPanel(aParticles,anEdges,aPanels,j)
			endif
		enddo
		startIndex = nOldPanels
	enddo
	
	do j=1,aParticles%N
		aParticles%dualArea(j) = DualPanelArea(aParticles,anEdges,aPanels,j)
	enddo	
end subroutine


subroutine DivideTriPanel(aParticles,anEdges,aPanels,panelIndex)
! Divides a triangular panel into 4 subpanels.
! Subpanel areas are assigned in this routine.
! Lagrangian coordinates of new panels are found by averaging Lagrangian coordinates
! of the divided panel.  
! No physical data (i.e., vorticity/tracer) are assigned here.  These values must be
! set separately.
	type(Particles), intent(inout) :: aParticles
	type(Edges), intent(inout) :: anEdges
	type(Panels), intent(inout) :: aPanels
	integer(kint), intent(in) :: panelIndex
	! local variables
	integer(kint) :: vertexIndices(3), edgeIndices(3), nestLevel, panelVerts(3)
	integer(kint) :: nParticles, nEdges, nPanels
	integer(kint) :: j, child1, child2
	
	if ( aPanels%N_Max - aPanels%N  < 4 ) then
		print *,"DivideTriPanel ERROR : not enough memory."
		return
	endif
	! Get info from state of current grid
	vertexIndices = aPanels%vertices(:,panelIndex)
	edgeIndices = aPanels%edges(:,panelIndex)
	nestLevel = aPanels%nest(panelIndex)
	nParticles = aParticles%N
	nEdges = anEdges%N
	nPanels = aPanels%N
	! connect parent vertices to subpanels
	aPanels%vertices(1,nPanels+1) = vertexIndices(1)
	aPanels%vertices(2,nPanels+2) = vertexIndices(2)
	aPanels%vertices(3,nPanels+3) = vertexIndices(3)
	! Parent edge 1
	if ( .NOT. anEdges%hasChildren(edgeIndices(1)) ) then ! divide edge 1
		! create new particle
		aParticles%x(:,nParticles+1) = SphereMidpoint(aParticles%x(:,vertexIndices(1)),aParticles%x(:,vertexIndices(2)))
		aParticles%x0(:,nParticles+1) = SphereMidpoint(aParticles%x0(:,vertexIndices(1)),aParticles%x0(:,vertexIndices(2)))
		! connect new particle to subpanels
		aPanels%vertices(2,nPanels+1) = nParticles+1
		aPanels%vertices(1,nPanels+2) = nParticles+1
		aPanels%vertices(3,nPanels+4) = nParticles+1
		! create new edges
		anEdges%hasChildren(edgeIndices(1)) = .True.
		anEdges%children(:,edgeIndices(1)) = [nEdges+1,nEdges+2]
		if ( panelIndex == anEdges%leftPanel(edgeIndices(1)) ) then ! edge 1 has positive orientation
			anEdges%verts(:,nEdges+1) = [vertexIndices(1),nParticles+1]
			anEdges%verts(:,nEdges+2) = [nParticles+1,vertexIndices(2)]
			aParticles%edges(:,nParticles+1) = [0,0,nEdges+1,0,0,nEdges+2]
			anEdges%leftPanel(nEdges+1) = nPanels+1
			anEdges%leftPanel(nEdges+2) = nPanels+2
			anEdges%rightPanel(nEdges+1) = anEdges%rightPanel(edgeIndices(1))
			anEdges%rightPanel(nEdges+2) = anEdges%rightPanel(edgeIndices(1))
			aPanels%edges(1,nPanels+1) = nEdges + 1
			aPanels%edges(1,nPanels+2) = nEdges + 2
			do j=1,6
				if (aParticles%edges(j,vertexIndices(1)) == edgeIndices(1)) then
					aParticles%edges(j,vertexIndices(1)) = nEdges+1
				endif
				if (aParticles%edges(j,vertexIndices(2)) == edgeIndices(1)) then
					aParticles%edges(j,vertexIndices(2)) = nEdges + 2
				endif
			enddo
		else ! edge 1 has negative orientation
			anEdges%verts(:,nEdges+1) = [vertexIndices(2),nParticles+1]
			anEdges%verts(:,nEdges+2) = [nParticles+1,vertexIndices(1)]
			aParticles%edges(:,nParticles+1) = [0,0,nEdges+1,0,0,nEdges+2]
			anEdges%leftPanel(nEdges+1) = anEdges%leftPanel(edgeIndices(1))
			anEdges%leftPanel(nEdges+2) = anEdges%leftPanel(edgeIndices(1))
			anEdges%rightPanel(nEdges+1) = nPanels+2
			anEdges%rightPanel(nEdges+2) = nPanels+1
			aPanels%edges(1,nPanels+1) = nEdges + 2
			aPanels%edges(1,nPanels+2) = nEdges+1
			do j=1,6
				if (aParticles%edges(j,vertexIndices(1)) == edgeIndices(1)) then
					aParticles%edges(j,vertexIndices(1)) = nEdges + 2
				endif
				if (aParticles%edges(j,vertexIndices(2)) == edgeIndices(1)) then
					aParticles%edges(j,vertexIndices(2)) = nEdges + 1
				endif
			enddo
		endif
		nParticles = nParticles+1
		nEdges = nEdges + 2
	else ! edge has already been divided
		child1 = anEdges%children(1,edgeIndices(1))
		child2 = anEdges%children(2,edgeIndices(1))
		! connect subpanels to parent edge midpoint
		aPanels%vertices(2,nPanels+1) = anEdges%verts(2,child1)
		aPanels%vertices(1,nPanels+2) = anEdges%verts(2,child1)
		aPanels%vertices(3,nPanels+4) = anEdges%verts(2,child1)
		if ( panelIndex == anEdges%leftPanel(edgeIndices(1))) then ! edge 1 has positive orientation
			anEdges%leftPanel(child1) = nPanels+1
			anEdges%leftPanel(child2) = nPanels+2
			aPanels%edges(1,nPanels+1) = child1
			aPanels%edges(1,nPanels+2) = child2
		else ! edge 1 has negative orientation
			anEdges%rightPanel(child1) = nPanels+2
			anEdges%rightPanel(child2) = nPanels+1
			aPanels%edges(1,nPanels+1) = child2
			aPanels%edges(1,nPanels+2) = child1
		endif
	endif
	
	! Parent Edge 2
	if ( .NOT. anEdges%hasChildren(edgeIndices(2)) ) then ! divide edge 2
		! create new particle
		aParticles%x(:,nParticles+1) = SphereMidpoint(aParticles%x(:,vertexIndices(2)),aParticles%x(:,vertexIndices(3)))
		aParticles%x0(:,nParticles+1) = SphereMidpoint(aParticles%x0(:,vertexIndices(2)),aParticles%x0(:,vertexIndices(3)))
		! connect subpanels to new particle
		aPanels%vertices(3,nPanels+2) = nParticles+1
		aPanels%vertices(2,nPanels+3) = nParticles+1
		aPanels%vertices(1,nPanels+4) = nParticles+1
		! create new edges
		anEdges%hasChildren(edgeIndices(2)) = .True.
		anEdges%children(:,edgeIndices(2)) = [nEdges + 1,nEdges+2]
		if ( panelIndex == anEdges%leftPanel(edgeIndices(2)) ) then ! edge 2 has positive orientation
			anEdges%verts(:,nEdges+1) = [vertexIndices(2),nParticles+1]
			anEdges%verts(:,nEdges+2) = [nParticles+1,vertexIndices(3)]
			aParticles%edges(:,nParticles+1) = [0,nEdges+1,0,0,nEdges+2,0]
			anEdges%leftPanel(nEdges+1) = nPanels+2
			anEdges%leftPanel(nEdges+2) = nPanels+3
			anEdges%rightPanel(nEdges+1) = anEdges%rightPanel(edgeIndices(2))
			anEdges%rightPanel(nEdges+2) = anEdges%rightPanel(edgeIndices(2))
			aPanels%edges(2,nPanels+2) = nEdges + 1
			aPanels%edges(2,nPanels+3) = nEdges + 2
			do j=1,6
				if (aParticles%edges(j,vertexIndices(2)) == edgeIndices(2)) then
					aParticles%edges(j,vertexIndices(2)) = nEdges+1
				endif
				if ( aParticles%edges(j,vertexIndices(3)) == edgeIndices(2)) then
					aParticles%edges(j,vertexIndices(3)) = nEdges + 2
				endif
			enddo
		else ! edge 2 has negative orientation
			anEdges%verts(:,nEdges + 1) = [vertexIndices(3),nParticles+1]
			anEdges%verts(:,nEdges + 2) = [nParticles+1,vertexIndices(2)]
			aParticles%edges(:,nParticles+1) = [0,nEdges+2,0,0,nEdges+1,0]
			anEdges%leftPanel(nEdges+1) = anEdges%leftPanel(edgeIndices(2))
			anEdges%leftPanel(nEdges+2) = anEdges%leftPanel(edgeIndices(2))
			anEdges%rightPanel(nEdges+1) = nPanels+3
			anEdges%rightPanel(nEdges+2) = nPanels+2
			aPanels%edges(2,nPanels+2) = nEdges + 2
			aPanels%edges(2,nPanels+3) = nEdges + 1
			do j=1,6
				if (aParticles%edges(j,vertexIndices(2)) == edgeIndices(2)) then
					aParticles%edges(j,vertexIndices(2)) = nEdges+2
				endif
				if (aParticles%edges(j,vertexIndices(3)) == edgeIndices(2)) then
					aParticles%edges(j,vertexIndices(3)) = nEdges + 1
				endif
			enddo
		endif
		nParticles = nParticles + 1
		nEdges = nEdges + 2
	else ! edge 2 has already been divided
		child1 = anEdges%children(1,edgeIndices(2))
		child2 = anEdges%children(2,edgeIndices(2))
		! connect subpanels to parent edge midpoint
		aPanels%vertices(3,nPanels+2) = anEdges%verts(2,child1)
		aPanels%vertices(2,nPanels+3) = anEdges%verts(2,child1)
		aPanels%vertices(1,nPanels+4) = anEdges%verts(2,child1)
		! connect subpanels to child edges
		if (panelIndex == anEdges%leftPanel(edgeIndices(2))) then! edge 2 has positive orientation
			aPanels%edges(2,nPanels+2) = child1
			aPanels%edges(2,nPanels+3) = child2
			anEdges%leftPanel(child1) = nPanels+2
			anEdges%leftPanel(child2) = nPanels+3
		else ! edge 2 has negative orientation
			aPanels%edges(2,nPanels+2) = child2
			aPanels%edges(2,nPanels+3) = child1
			anEdges%rightPanel(child1) = nPanels+3
			anEdges%rightPanel(child2) = nPanels+2
		endif
	endif
	
	! Parent Edge 3
	if ( .NOT. anEdges%hasChildren(edgeIndices(3)) ) then ! divide edge 3
		! create new particle
		aParticles%x(:,nParticles+1) = SphereMidpoint(aParticles%x(:,vertexIndices(3)),aParticles%x(:,vertexIndices(1)))
		aParticles%x0(:,nParticles+1) = SphereMidpoint(aParticles%x0(:,vertexIndices(3)),aParticles%x0(:,vertexIndices(1)))
		! connect subpanels to new particle
		aPanels%vertices(1,nPanels+3) = nParticles+1
		aPanels%vertices(2,nPanels+4) = nParticles+1
		aPanels%vertices(3,nPanels+1) = nParticles+1
		! create new edges
		anEdges%hasChildren(edgeIndices(3)) = .True.
		anEdges%children(:,edgeIndices(3)) = [nEdges + 1,nEdges + 2]
		if ( panelIndex == anEdges%leftPanel(edgeIndices(3)) ) then ! edge 3 has positive orientation
			anEdges%verts(:,nEdges+1) = [vertexIndices(3),nParticles+1]
			anEdges%verts(:,nEdges+2) = [nParticles+1,vertexIndices(1)]
			aParticles%edges(:,nParticles+1) = [nEdges+2,0,0,nEdges+1,0,0]
			anEdges%leftPanel(nEdges+1) = nPanels+3
			anEdges%leftPanel(nEdges+2) = nPanels+1
			anEdges%rightPanel(nEdges+1) = anEdges%rightPanel(edgeIndices(3))
			anEdges%rightPanel(nEdges+2) = anEdges%rightPanel(edgeIndices(3))
			aPanels%edges(3,nPanels+3) = nEdges + 1
			aPanels%edges(3,nPanels+1) = nEdges + 2
			do j=1,6
				if (aParticles%edges(j,vertexIndices(3)) == edgeIndices(3)) then
					aParticles%edges(j,vertexIndices(3)) = nEdges+1
				endif
				if (aParticles%edges(j,vertexIndices(1)) == edgeIndices(3)) then
					aParticles%edges(j,vertexIndices(1)) = nEdges+2
				endif
			enddo
		else ! edge 3 has negative orientation
			anEdges%verts(:,nEdges+1) = [vertexIndices(1),nParticles+1]
			anEdges%verts(:,nEdges+2) = [nParticles+1,vertexIndices(3)]
			aParticles%edges(:,nParticles+1) = [nEdges+1,0,0,nEdges+2,0,0]
			anEdges%leftPanel(nEdges+1) = anEdges%leftPanel(edgeIndices(3))
			anEdges%leftPanel(nEdges+2) = anEdges%leftPanel(edgeIndices(3))
			anEdges%rightPanel(nEdges+1) = nPanels+1
			anEdges%rightPanel(nEdges+2) = nPanels+3
			aPanels%edges(3,nPanels+3) = nEdges + 2
			aPanels%edges(3,nPanels+1) = nEdges + 1
			do j=1,6
				if (aParticles%edges(j,vertexIndices(3)) == edgeIndices(3)) then
					aParticles%edges(j,vertexIndices(3)) = nEdges+2
				endif
				if (aParticles%edges(j,vertexIndices(1)) == edgeIndices(3)) then
					aParticles%edges(j,vertexIndices(1)) = nEdges+1
				endif
			enddo
		endif
		nParticles = nParticles + 1
		nEdges = nEdges + 2
	else ! edge 3 has already been divided
		child1 = anEdges%children(1,edgeIndices(3))
		child2 = anEdges%children(2,edgeIndices(3))
		! connect subpanels to parent edge midpoint
		aPanels%vertices(1,nPanels+3) = anEdges%verts(2,child1)
		aPanels%vertices(2,nPanels+4) = anEdges%verts(2,child1)
		aPanels%vertices(3,nPanels+1) = anEdges%verts(2,child1)
		if (panelIndex == anEdges%leftPanel(edgeIndices(3)) ) then ! edge 3 has positive orientation
			aPanels%edges(3,nPanels+3) = child1
			aPanels%edges(3,nPanels+1) = child2
			anEdges%leftPanel(child1) = nPanels+3
			anEdges%leftPanel(child2) = nPanels+1
		else ! edge 3 has negative orientaiton
			aPanels%edges(3,nPanels+3) = child2
			aPanels%edges(3,nPanels+1) = child1
			anEdges%rightPanel(child1) = nPanels+1
			anEdges%rightPanel(child2) = nPanels+3
		endif
	endif
	
	! Interior edges
	anEdges%verts(:,nEdges+1) = [aPanels%vertices(1,nPanels+4),aPanels%vertices(2,nPanels+4)]
	anEdges%leftPanel(nEdges+1) = nPanels+4
	anEdges%rightPanel(nEdges+1) = nPanels+3
	aPanels%edges(1,nPanels+4) = nEdges+1
	aPanels%edges(1,nPanels+3) = nEdges+1
	do j=1,6
	   if (aParticles%edges(j,aPanels%vertices(1,nPanels+4))==0) then
		  aParticles%edges(j,aPanels%vertices(1,nPanels+4)) = nEdges + 1
		  exit
	   endif
	enddo
	do j=1,6
	   if (aParticles%edges(j,aPanels%vertices(2,nPanels+4))==0) then
		  aParticles%edges(j,aPanels%vertices(2,nPanels+4)) = nEdges + 1
		  exit
	   endif
	enddo

	anEdges%verts(:,nEdges+2) = [aPanels%vertices(2,nPanels+4),aPanels%vertices(3,nPanels+4)]
	anEdges%leftPanel(nEdges+2) = nPanels+4
	anEdges%rightPanel(nEdges+2) = nPanels+1
	aPanels%edges(2,nPanels+4) = nEdges + 2
	aPanels%edges(2,nPanels+1) = nEdges + 2
	do j=1,6
	   if (aParticles%edges(j,aPanels%vertices(2,nPanels+4))==0) then
		  aParticles%edges(j,aPanels%vertices(2,nPanels+4)) = nEdges + 2
		  exit
	   endif
	enddo
	do j=1,6
	   if (aParticles%edges(j,aPanels%vertices(3,nPanels+4))==0) then
		  aParticles%edges(j,aPanels%vertices(3,nPanels+4)) = nEdges + 2 
		  exit
	   endif
	enddo
	
	anEdges%verts(:,nEdges+3) = [aPanels%vertices(3,nPanels+4),aPanels%vertices(1,nPanels+4)]
	anEdges%leftPanel(nEdges+3) = nPanels+4
	anEdges%rightPanel(nEdges+3) = nPanels+2
	aPanels%edges(3,nPanels+4) = nEdges+3
	aPanels%edges(3,nPanels+2) = nEdges+3
	do j=1,6
	   if (aParticles%edges(j,aPanels%vertices(3,nPanels+4)) == 0) then
		  aParticles%edges(j,aPanels%vertices(3,nPanels+4)) = nEdges + 3
		  exit
	   endif
	enddo
	do j=1,6
	   if (aParticles%edges(j,aPanels%vertices(1,nPanels+4)) == 0) then
		  aParticles%edges(j,aPanels%vertices(1,nPanels+4)) = nEdges + 3
		  exit
	   endif
	enddo
	
	nEdges = nEdges + 3
	
	! Subpanel centers and areas
	do j=1,4
		panelVerts = aPanels%vertices(:,nPanels+j)
		aPanels%x(:,nPanels+j) = SphereTriCenter(aParticles%x(:,panelVerts(1)),aParticles%x(:,panelVerts(2)),&
			aParticles%x(:,panelVerts(3)))
		aPanels%x0(:,nPanels+j) = SphereTriCenter(aParticles%x0(:,panelVerts(1)),aParticles%x0(:,panelVerts(2)),&
			aParticles%x0(:,panelVerts(3)))
		aPanels%area(nPanels+j) = TriPanelArea(aparticles,anedges,apanels,nPanels+j)
	enddo
	
	! Center particle
	aPanels%x(:,panelIndex) = aPanels%x(:,nPanels+4)
	aPanels%x0(:,panelIndex) = aPanels%x0(:,nPanels+4)
	! TO DO : Don't duplicate center panel's particle
	
	
	! Update data structures and return
	aPanels%area(panelIndex) = 0.0_kreal ! exclude divided panel from dynamics
	aPanels%hasChildren(panelIndex) = .True.
	aPanels%nest(nPanels+1:nPanels+4) = nestLevel+1
	aPanels%children(:,panelIndex) = [1,2,3,4] + nPanels
	aPanels%N = nPanels+4
	aPanels%N_Active = aPanels%N_Active + 3
	
	anEdges%N = nEdges
	
	aParticles%N = nParticles
	
end subroutine


subroutine DivideQuadPanel(aParticles,anEdges,aPanels,panelIndex)
! Divides a quadrilateral panel into 4 subpanels.
! Subpanel areas are assigned in this routine.
! Lagrangian coordinates of new panels are found by averaging Lagrangian coordinates
! of the divided panel.  
! No physical data (i.e., vorticity/tracer) are assigned here.  These values must be
! set separately.
	type(Particles), intent(inout) :: aParticles
	type(Edges), intent(inout) :: anEdges
	type(Panels), intent(inout) :: aPanels
	integer(kint), intent(in) :: panelIndex
	! local variables
	integer(kint) :: vertexIndices(4), edgeIndices(4), nestLevel, panelVerts(4)
	integer(kint) :: nParticles, nPanels, nEdges
	integer(kint) :: j, child1, child2, noEdgeSlot
	
	if ( aPanels%N_Max - aPanels%N < 4) then
		print *,"DivideQuadPanel ERROR : not enough memory."
		return
	endif
	! Get info on the grid's current state
	vertexIndices = aPanels%vertices(:,panelIndex)
	edgeIndices = aPanels%edges(:,panelIndex)
	nestLevel = aPanels%nest(panelIndex)
	nParticles = aParticles%N
	nPanels = aPanels%N
	nEdges = anEdges%N
	! connect parent vertices to subpanels
	apanels%vertices(1,nPanels+1) = vertexIndices(1)
	apanels%vertices(2,nPanels+2) = vertexIndices(2)
	apanels%vertices(3,nPanels+3) = vertexIndices(3)
	apanels%vertices(4,nPanels+4) = vertexIndices(4)
	! Parent Edge 1
	if ( .NOT. anEdges%hasChildren(edgeIndices(1)) ) then ! divide edge 1
		! create new particle
		aParticles%x(:,nParticles+1) = SphereMidpoint(aParticles%x(:,vertexIndices(1)),aParticles%x(:,vertexIndices(2)))
		aParticles%x0(:,nParticles+1) = SphereMidpoint(aParticles%x0(:,vertexIndices(1)),aParticles%x0(:,vertexIndices(2)))
		! connect subpanels to new particle
		aPanels%vertices(2,nPanels+1) = nParticles+1
		aPanels%vertices(1,nPanels+2) = nParticles+1
		! create new edges
		anEdges%hasChildren(edgeIndices(1)) = .True.
		anEdges%children(:,edgeIndices(1)) = [nEdges+1,nEdges+2]
		if ( panelIndex == anEdges%leftPanel(edgeIndices(1))) then ! edge 1 has positive orientation
			! connect children edges to particles
			anEdges%verts(:,nEdges+1) = [vertexIndices(1),nParticles+1]
			anEdges%verts(:,nEdges+2) = [nParticles+1,vertexIndices(2)]
			! connect new particle to new edges
			aParticles%edges(:,nParticles+1) = [nEdges+1,0,nEdges+2,0]
			! connect new edges to panels
			anEdges%leftPanel(nEdges+1) = nPanels+1
			anEdges%leftPanel(nEdges+2) = nPanels+2
			anEdges%rightPanel(nEdges+1) = anEdges%rightPanel(edgeIndices(1))
			anEdges%rightPanel(nEdges+2) = anEdges%rightPanel(edgeIndices(1))
			aPanels%edges(1,nPanels+1) = nEdges+1
			aPanels%edges(1,nPanels+2) = nEdges+2
			! replace parent edge with children at parent vertices
			do j=1,4
			   if (aParticles%edges(j,vertexIndices(1)) == edgeIndices(1)) then
				  aParticles%edges(j,vertexIndices(1)) = nEdges+1 
			   endif
			   if (aParticles%edges(j,vertexIndices(2)) == edgeIndices(1)) then
				   aParticles%edges(j,vertexIndices(2)) = nEdges+2
			   endif
			enddo	
		else	! edge 1 has negative orientation 
			! record edge vertices
		 	anEdges%verts(:,nEdges+1) = [vertexIndices(2),nParticles+1]
			anEdges%verts(:,nEdges+2) = [nParticles+1,vertexIndices(1)]
			! connect edges to panels
			aParticles%edges(:,nParticles+1) = [nEdges+2,0,nEdges+1,0]
			anEdges%leftPanel(nEdges+1) = anEdges%leftPanel(edgeIndices(1))
			anEdges%leftPanel(nEdges+2) = anEdges%leftPanel(edgeIndices(1))
			anEdges%rightPanel(nEdges+1) = nPanels+2
			anEdges%rightPanel(nEdges+2) = nPanels+1
			aPanels%edges(1,nPanels+1) = nEdges+2
			aPanels%edges(1,nPanels+2) = nEdges+1
			! replace parent edge with child edges at parent vertices
			do j=1,4
			   if (aParticles%edges(j,vertexIndices(1)) == edgeIndices(1)) then
				  aParticles%edges(j,vertexIndices(1)) = nEdges+2 
			   endif
			   if (aParticles%edges(j,vertexIndices(2)) == edgeIndices(1)) then
				   aParticles%edges(j,vertexIndices(2)) = nEdges+1
			   endif
			enddo
		endif
		nParticles = nParticles + 1
		nEdges = nEdges + 2
	else ! edge 1 has been divided previously
		child1 = anEdges%children(1,edgeIndices(1))
		child2 = anEdges%children(2,edgeIndices(1))
		! connect subpanels to divided edge midpoint
		aPanels%vertices(2,nPanels+1) = anEdges%verts(2,child1)
		aPanels%vertices(1,nPanels+2) = anEdges%verts(2,child1)
		! connect panels and child edges
		if ( panelIndex == anEdges%leftPanel(edgeIndices(1)) ) then ! edge 1 has positive orientation
			aPanels%edges(1,nPanels+1) = child1
			aPanels%edges(1,nPanels+2) = child2
			anEdges%leftPanel(child1) = nPanels+1
			anEdges%leftPanel(child2) = nPanels+2
		else ! edge 1 has negative orientation
			aPanels%edges(1,nPanels+1) = child2
			aPanels%edges(1,nPanels+2) = child1
			anEdges%rightPanel(child1) = nPanels+2
			anEdges%rightPanel(child2) = nPanels+1
		endif
	endif
	! Parent Edge 2
	if ( .NOT. anEdges%hasChildren(edgeIndices(2)) ) then ! divide edge 2
		! create new particle
		aParticles%x(:,nParticles+1) = SphereMidpoint(aParticles%x(:,vertexIndices(2)),aParticles%x(:,vertexIndices(3)))
		aParticles%x0(:,nParticles+1) = SphereMidpoint(aParticles%x0(:,vertexIndices(2)),aParticles%x0(:,vertexIndices(3)))
		! connect subpanels to new particle
		aPanels%vertices(3,nPanels+2) = nParticles+1
		aPanels%vertices(2,nPanels+3) = nParticles+1
		! create new edges
		anEdges%hasChildren(edgeIndices(2)) = .True.
		anEdges%children(:,edgeIndices(2)) = [nEdges+1,nEdges+2]
		if ( panelIndex == anEdges%leftPanel( edgeIndices(2)) ) then ! edge 2 has positive orientation
			anEdges%verts(:,nEdges+1) = [vertexIndices(2),nParticles+1]
			anEdges%verts(:,nEdges+2) = [nParticles+1,vertexIndices(3)]
			aParticles%edges(:,nParticles+1) = [0,nEdges+1,0,nEdges+2]
			anEdges%leftPanel(nEdges+1) = nPanels+2
			anEdges%leftPanel(nEdges+2) = nPanels+3
			anEdges%rightPanel(nEdges+1) = anEdges%rightPanel(edgeIndices(2))
			anEdges%rightPanel(nEdges+2) = anEdges%rightPanel(edgeIndices(2))
			aPanels%edges(2,nPanels+2) = nEdges+1
			aPanels%edges(2,nPanels+3) = nEdges+2
			do j=1,4
			   if (aParticles%edges(j,vertexIndices(2)) == edgeIndices(2)) then
				   aParticles%edges(j,vertexIndices(2)) = nEdges + 1
			   endif
			   if (aParticles%edges(j,vertexIndices(3)) == edgeIndices(2)) then
				   aParticles%edges(j,vertexIndices(3)) = nEdges + 2
			   endif
			enddo
		else ! edge 2 has negative orientation
			anEdges%verts(:,nEdges+1) = [vertexIndices(3),nParticles+1]
			anEdges%verts(:,nEdges+2) = [nParticles+1,vertexIndices(2)]
			aParticles%edges(:,nParticles+1) = [0,nEdges+2,0,nEdges+1]
			anEdges%leftPanel(nEdges+1) = anEdges%leftPanel(edgeIndices(2))
			anEdges%leftPanel(nEdges+2) = anEdges%leftPanel(edgeIndices(2))
			anEdges%rightPanel(nEdges+1) = nPanels+3
			anEdges%rightPanel(nEdges+2) = nPanels+2
			aPanels%edges(2,nPanels+2) = nEdges+2
			aPanels%edges(2,nPanels+3) = nEdges+1
			do j=1,4
			   if (aParticles%edges(j,vertexIndices(2)) == edgeIndices(2)) then
				   aParticles%edges(j,vertexIndices(2)) = nEdges + 2
			   endif
			   if (aParticles%edges(j,vertexIndices(3)) == edgeIndices(2)) then
				   aParticles%edges(j,vertexIndices(3)) = nEdges + 1
			   endif
			enddo
		endif
		nParticles = nParticles+1
		nEdges = nEdges + 2
	else ! edge 2 has already been divided
		child1 = anEdges%children(1,edgeIndices(2))
		child2 = anEdges%children(2,edgeIndices(2))
		! connect subpaneles to divided edge midpoint
		aPanels%vertices(3,nPanels+2) = anEdges%verts(2,child1)
		aPanels%vertices(2,nPanels+3) = anEdges%verts(2,child1)
		if ( panelIndex == anEdges%leftPanel(edgeIndices(2)) ) then !edge 2 has positive orientation
			aPanels%edges(2,nPanels+2) = child1
			aPanels%edges(2,nPanels+3) = child2
			anEdges%leftPanel(child1) = nPanels+2
			anEdges%leftPanel(child2) = nPanels+3
		else ! edge 2 has negative orientation
			aPanels%edges(2,nPanels+2) = child2
			aPanels%edges(2,nPanels+3) = child1
			anEdges%rightPanel(child1) = nPanels+3
			anEdges%rightPanel(child2) = nPanels+2
		endif
	endif
	
	! Parent Edge 3
	if ( .NOT. anEdges%hasChildren(edgeIndices(3)) ) then !divide edge 3
		! create new particle
		aParticles%x(:,nParticles+1) = SphereMidpoint(aParticles%x(:,vertexIndices(3)),aParticles%x(:,vertexIndices(4)))
		aParticles%x0(:,nParticles+1) = SphereMidpoint(aParticles%x0(:,vertexIndices(3)),aParticles%x0(:,vertexIndices(4)))
		! connect subpanels to new particle
		aPanels%vertices(4,nPanels+3) = nParticles+1
		aPanels%vertices(3,nPanels+4) = nParticles+1
		! create new edges
		anEdges%hasChildren(edgeIndices(3)) = .TRUE.
		anEdges%children(:,edgeIndices(3)) = [nEdges+1,nEdges+2]
		if ( panelIndex == anEdges%leftPanel(edgeIndices(3)) ) then ! edge 3 has positive orientation
			anEdges%verts(:,nEdges+1) = [vertexIndices(3),nParticles+1]
			anEdges%verts(:,nEdges+2) = [nParticles+1,vertexIndices(4)]
			aParticles%edges(:,nParticles+1) = [nEdges+2,0,nEdges+1,0]
			anEdges%leftPanel(nEdges+1) = nPanels+3
			anEdges%leftPanel(nEdges+2) = nPanels+4
			anEdges%rightPanel(nEdges+1) = anEdges%rightPanel(edgeIndices(3))
			anEdges%rightPanel(nEdges+2) = anEdges%rightPanel(edgeIndices(3))
			aPanels%edges(3,nPanels+3) = nEdges+1
			aPanels%edges(3,nPanels+4) = nEdges+2
			do j=1,4
			   if (aParticles%edges(j,vertexIndices(3)) == edgeIndices(3)) then
				   aParticles%edges(j,vertexIndices(3)) = nEdges + 1
			   endif
			   if (aParticles%edges(j,vertexIndices(4)) == edgeIndices(3)) then
				   aParticles%edges(j,vertexIndices(4)) = nEdges + 2
			   endif
			enddo
		else ! edge 3 has negative orientation
			anEdges%verts(:,nEdges+1) = [vertexIndices(4),nParticles+1]
			anEdges%verts(:,nEdges+2) = [nParticles+1,vertexIndices(3)]
			aParticles%edges(:,nParticles+1) = [nEdges+1,0,nEdges+2,0]
			anEdges%leftPanel(nEdges+1) = anEdges%leftPanel(edgeIndices(3))
			anEdges%leftPanel(nEdges+2) = anEdges%leftPanel(edgeIndices(3))
			anEdges%rightPanel(nEdges+1) = nPanels+4
			anEdges%rightPanel(nEdges+2) = nPanels+3
			aPanels%edges(3,nPanels+3) = nEdges+2
			aPanels%edges(3,nPanels+4) = nEdges+1
			do j=1,4
			   if (aParticles%edges(j,vertexIndices(3)) == edgeIndices(3)) then
				   aParticles%edges(j,vertexIndices(3)) = nEdges + 2
			   endif
			   if (aParticles%edges(j,vertexIndices(4)) == edgeIndices(3)) then
				   aParticles%edges(j,vertexIndices(4)) = nEdges + 1
			   endif
			enddo
		endif
		nParticles = nParticles+1
		nEdges = nEdges + 2
	else ! edge 3 has already been divided
		child1 = anEdges%children(1,edgeIndices(3))
		child2 = anEdges%children(2,edgeIndices(3))
		! connect subpanels to divided edge midpoint
		aPanels%vertices(4,nPanels+3) = anEdges%verts(2,child1)
		aPanels%vertices(3,nPanels+4) = anEdges%verts(2,child1)
		if (panelIndex == anEdges%leftPanel(edgeIndices(3)) ) then! edges 3 has positive orientation
			aPanels%edges(3,nPanels+3) = child1
			aPanels%edges(3,nPanels+4) = child2
			anEdges%leftPanel(child1) = nPanels+3
			anEdges%leftPanel(child2) = nPanels+4
		else
			aPanels%edges(3,nPanels+3) = child2
			aPanels%edges(3,nPanels+4) = child1
			anEdges%rightPanel(child1) = nPanels+4
			anEdges%rightPanel(child2) = nPanels+3
		endif
	endif
	! Parent Edge 4
	if ( .NOT. anEdges%hasChildren(edgeIndices(4)) ) then ! divide edge 4
		! create new particle
		aParticles%x(:,nParticles+1) = SphereMidpoint(aParticles%x(:,vertexIndices(4)),aParticles%x(:,vertexIndices(1)))
		aParticles%x0(:,nParticles+1) = SphereMidpoint(aParticles%x0(:,vertexIndices(4)),aParticles%x0(:,vertexIndices(1)))
		! connect new particle to subpanels
		aPanels%vertices(1,nPanels+4) = nParticles+1
		aPanels%vertices(4,nPanels+1) = nParticles+1
		! create new edges
		anEdges%hasChildren(edgeIndices(4)) = .TRUE.
		anEdges%children(:,edgeIndices(4)) = [nEdges+1,nEdges+2]
		if ( panelIndex == anEdges%leftPanel(edgeIndices(4)) ) then ! edge 4 has positive orientation
			anEdges%verts(:,nEdges+1) = [vertexIndices(4),nParticles+1]
			anEdges%verts(:,nEdges+2) = [nParticles+1,vertexIndices(1)]
			aParticles%edges(:,nParticles+1) = [0,nEdges+2,0,nEdges+1]
			anEdges%leftPanel(nEdges+1) = nPanels+4
			anEdges%leftPanel(nEdges+2) = nPanels+1
			anEdges%rightPanel(nEdges+1) = anEdges%rightPanel(edgeIndices(4))
			anEdges%rightPanel(nEdges+2) = anEdges%rightPanel(edgeIndices(4))
			aPanels%edges(4,nPanels+4) = nEdges+1
			aPanels%edges(4,nPanels+1) = nEdges+2
			do j=1,4
			   if (aParticles%edges(j,vertexIndices(4)) == edgeIndices(4)) then
				  aParticles%edges(j,vertexIndices(4)) = nEdges + 1 
			   endif
			   if (aParticles%edges(j,vertexIndices(1)) == edgeIndices(4)) then
				   aParticles%edges(j,vertexIndices(1)) = nEdges + 2
			   endif
			enddo
		else ! edge 4 has negative orientation
			anEdges%verts(:,nEdges+1) = [vertexIndices(1),nParticles+1]
			anEdges%verts(:,nEdges+2) = [nParticles+1,vertexIndices(4)]
			aParticles%edges(:,nParticles+1) = [0,nEdges+1,0,nEdges+2]
			anEdges%leftPanel(nEdges+1) = anEdges%leftPanel(edgeIndices(4))
			anEdges%leftPanel(nEdges+2) = anEdges%leftPanel(edgeIndices(4))
			anEdges%rightPanel(nEdges+1) = nPanels+1
			anEdges%rightPanel(nEdges+2) = nPanels+4
			aPanels%edges(4,nPanels+4) = nEdges+2
			aPanels%edges(4,nPanels+1) = nEdges+1
			do j=1,4
			   if (aParticles%edges(j,vertexIndices(4)) == edgeIndices(4)) then
				  aParticles%edges(j,vertexIndices(4)) = nEdges + 2 
			   endif
			   if (aParticles%edges(j,vertexIndices(1)) == edgeIndices(4)) then
				   aParticles%edges(j,vertexIndices(1)) = nEdges + 1
			   endif
			enddo
		endif
		nParticles = nParticles+1
		nEdges = nEdges + 2
	else ! edge 4 has aleady been divided
		child1 = anEdges%children(1,edgeIndices(4))
		child2 = anEdges%children(2,edgeIndices(4))
		! connect subpanels to midpoint of parent edge
		aPanels%vertices(1,nPanels+4) = anEdges%verts(2,child1)
		aPanels%vertices(4,nPanels+1) = anEdges%verts(2,child1)
		if ( panelIndex == anEdges%leftPanel(edgeIndices(4)) ) then ! edge 4 has positive orientation
		   aPanels%edges(4,nPanels+4) = child1
		   aPanels%edges(4,nPanels+1) = child2
		   anEdges%leftPanel(child1) = nPanels+4
		   anEdges%leftPanel(child2) = nPanels+1
		else ! edge 4 has negative orientation
			aPanels%edges(4,nPanels+4) = child2
			aPanels%edges(4,nPanels+1) = child1
			anEdges%rightPanel(child1) = nPanels+1
			anEdges%rightPanel(child2) = nPanels+4
		endif
	endif
	
	! Center particle
	aParticles%x(:,nParticles+1) = aPanels%x(:,panelIndex)
	aParticles%x0(:,nParticles+1) = aPanels%x0(:,panelIndex)
	aParticles%edges(:,nParticles+1) = [4,1,3,2]+nEdges
	aPanels%vertices(3,nPanels+1) = nParticles+1
	aPanels%vertices(4,nPanels+2) = nParticles+1
	aPanels%vertices(1,nPanels+3) = nParticles+1
	aPanels%vertices(2,nPanels+4) = nParticles+1
	
	nParticles = nParticles + 1
	
	! Interior edges
	anEdges%verts(:,nEdges+1) = [aPanels%vertices(2,nPanels+1),aPanels%vertices(3,nPanels+1)]
	anEdges%leftPanel(nEdges+1) = nPanels+1
	anEdges%rightPanel(nEdges+1)= nPanels+2
	aPanels%edges(2,nPanels+1) = nEdges+1
	aPanels%edges(4,nPanels+2) = nEdges+1
	if ( aParticles%edges(4,aPanels%vertices(2,nPanels+1)) == 0 ) then
		aParticles%edges(4,aPanels%vertices(2,nPanels+1)) = nEdges+1
	else
		noEdgeSlot = 0
		do j=1,4
			if (aParticles%edges(j,aPanels%vertices(2,nPanels+1)) == 0) then
				aParticles%edges(j,aPanels%vertices(2,nPanels+1)) = nEdges+1
				exit
			else
				noEdgeSlot = noEdgeSlot + 1
			endif
		enddo
		if (noEdgeSlot >= 4) then
			print *,"DivideQuadPanel ERROR: no space for interior edge 1 at parent edge 1 midpoint."
			return
		endif
	endif
	
	anEdges%verts(:,nEdges+2) = [aPanels%vertices(2,nPanels+4),aPanels%vertices(3,nPanels+4)]
	anEdges%leftPanel(nEdges+2) = nPanels+4
	anEdges%rightPanel(nEdges+2)= nPanels+3
	aPanels%edges(2,nPanels+4) = nEdges+2
	aPanels%edges(4,nPanels+3) = nEdges+2
	if (aParticles%edges(2,aPanels%vertices(3,nPanels+4)) == 0) then
      	aParticles%edges(2,aPanels%vertices(3,nPanels+4)) = nEdges + 2
	else
		noEdgeSlot = 0
		do j=1,4
			if (aParticles%edges(j,aPanels%vertices(3,nPanels+4)) == 0) then
				aParticles%edges(j,aPanels%vertices(3,nPanels+4)) = nEdges + 2
				exit
			else
				noEdgeSlot = noEdgeSlot + 1
			endif
		enddo
		if (noEdgeSlot >= 4) then
			print *, 'DivideQuadPanel ERROR : no space for interior edge 2 at parent edge 3 midpoint.'
			return
		endif
	endif
	
	anEdges%verts(:,nEdges+3) = [aPanels%vertices(2,nPanels+3),aPanels%vertices(1,nPanels+3)]
	anEdges%leftPanel(nEdges+3) = nPanels+2
	anEdges%rightPanel(nEdges+3) = nPanels+3
	aPanels%edges(3,nPanels+2) = nEdges+3
	aPanels%edges(1,nPanels+3) = nEdges+3
	if (aParticles%edges(1,aPanels%vertices(2,nPanels+3)) == 0) then
		aParticles%edges(1,aPanels%vertices(2,nPanels+3)) = nEdges + 3
	else
		noEdgeSlot = 0
		do j=1,4
			if (aParticles%edges(j,aPanels%vertices(2,nPanels+3)) == 0) then
				aParticles%edges(j,aPanels%vertices(2,nPanels+3)) = nEdges+3
				exit
			else
				noEdgeSlot = noEdgeSlot + 1
			endif
		enddo
		if (noEdgeSlot >= 4) then
			print *, 'DivideQuadPanel ERROR : no space for interior edge 3 at parent edge 2 midpoint.'
			return
		endif
	endif
	
	anEdges%verts(:,nEdges+4) = [aPanels%vertices(2,nPanels+4),aPanels%vertices(1,nPanels+4)]
	anEdges%leftPanel(nEdges+4) = nPanels+1
	anEdges%rightPanel(nEdges+4) = nPanels+4
	aPanels%edges(1,nPanels+4) = nEdges+4
	aPanels%edges(3,nPanels+1) = nEdges+4
	if (aParticles%edges(3,aPanels%vertices(4,nPanels+1)) == 0) then
		aParticles%edges(3,aPanels%vertices(4,nPanels+1)) = nEdges + 4
	else
		noEdgeSlot = 0
		do j=1,4
			if (aParticles%edges(j,aPanels%vertices(4,nPanels+1)) == 0) then
				aParticles%edges(j,aPanels%vertices(4,nPanels+1)) = nEdges + 4
				exit
			else
				noEdgeSlot = noEdgeSlot + 1
			endif
		enddo
		if (noEdgeSlot >= 4) then
			print *,'DivideQuadPanel ERROR: no space for interior edge 4 at parent edge 4 midpoint.'
			return
		endif
	endif
	
	nEdges = nEdges + 4
	
	! Subpanel centers and areas
	do j=1,4
		panelVerts = aPanels%vertices(:,nPanels+j)
		aPanels%x(:,nPanels+j) = SphereQuadCenter(aParticles%x(:,panelVerts(1)),&
			aParticles%x(:,panelVerts(2)),aParticles%x(:,panelVerts(3)),aParticles%x(:,panelVerts(4)))
 	    aPanels%x0(:,nPanels+j) = SphereQuadCenter(aParticles%x0(:,panelVerts(1)),&
 	    	aParticles%x0(:,panelVerts(2)),aParticles%x0(:,panelVerts(3)),aParticles%x0(:,panelVerts(4)))
	    aPanels%area(nPanels+j) = QuadPanelArea(aparticles,anedges,apanels,nPanels+j)
	enddo
	
	! Update data structures and return
	aPanels%area(panelIndex) = 0.0_kreal ! exclude divided panel from future dynamics
	aPanels%hasChildren(panelIndex) = .True.
	aPanels%nest(nPanels+1:nPanels+4) = nestLevel+1
	aPanels%children(:,panelIndex) = [1,2,3,4] + nPanels
	aPanels%N = nPanels+4
	aPanels%N_Active = aPanels%N_Active + 3
	
	anEdges%N = nEdges
	
	aParticles%N = nParticles	
	
end subroutine


subroutine DividePanel(aParticles,anEdges,aPanels,panelindex)
! Simple interface for calling appropriate DividePanel routine, based on input panels
! kind (TRI_PANEL or QUAD_PANEL).
	type(Particles), intent(inout) :: aParticles
	type(Edges), intent(inout) :: anEdges
	type(Panels), intent(inout) :: aPanels
	integer(kint), intent(in) :: panelIndex
	if ( GetPanelKind(aPanels) == TRI_PANEL ) then
		call DivideTriPanel(aParticles,anEdges,aPanels,panelIndex)
	elseif (GetPanelKind(aPanels) == QUAD_PANEL ) then
		call DivideQuadPanel(aParticles,anEdges,aPanels,panelIndex)
	else
		print *,"DividePanel ERROR : invalid panelKind."
		return
	endif
end subroutine


subroutine GatherPanels(aPanels, activePanels, activeMap, passivePanels, passiveMap)
! Separates divided panels (passive panels) from low-level panels (active panels).
! It keeps track of the memory location of each panel and stores them in integer arrays
! activeMap and passiveMap, for later use in the ScatterPanels subroutine.
!
! Must be called any time the low-level panels are needed without their tree structure.
!
	type(Panels), intent(in) :: aPanels
	type(Panels), intent(out) :: activePanels, passivePanels
	integer(kint), intent(out) :: activeMap(:), passiveMap(:)
	integer(kint) :: j, activeK, passiveK, nActive, nPassive
	nActive = aPanels%N_Active
	nPassive = aPanels%N - nActive
	if ( size(activeMap) /= nActive) then
		call LogMessage(log,ERROR_LOGGING_LEVEL,logKey,"GatherPanels ERROR : active size 1")
!		print *,"GatherPanels ERROR : active size mismatch 1."
		return
	endif
	if ( nActive /= activePanels%N) then
		call LogMessage(log,ERROR_LOGGING_LEVEL,logKey,"GatherPanels ERROR : active size 2")
!		print *,"GatherPanels ERROR : active size mismatch 2."
		return
	endif
	if ( nActive /= activePanels%N_Active) then
		call LogMessage(log,ERROR_LOGGING_LEVEL,logKey,"GatherPanels ERROR : active size 3")
!		print *,"GatherPanels ERROR : active size mismatch 3."
		return
	endif
	if ( size(passiveMap) /= nPassive) then
		call LogMessage(log,ERROR_LOGGING_LEVEL,logKey,"GatherPanels ERROR : passive size 1")
!		print *,"GatherPanels ERROR : passive size mismatch 1."
		return
	endif
	if ( nPassive /= passivePanels%N) then
		call LogMessage(log,ERROR_LOGGING_LEVEL,logKey,"GatherPanels ERROR : passive size 1")
!		print *, "GatherPanels ERROR : passive size mismatch 2."
		return
	endif	
	call LogMessage(log,DEBUG_LOGGING_LEVEL,logKey,' entering GatherPanels .')
	
	activeK = 1
	passiveK = 1
	do j=1,aPanels%N
		if ( .NOT. aPanels%hasChildren(j) ) then ! panel j is active
			call Copy(activePanels,activeK,aPanels,j)
			activeMap(activeK) = j
			activeK = activeK + 1
		else ! panel j is passive
			call Copy(passivePanels,passiveK,aPanels,j)
			passiveMap(passiveK) = j
			passiveK = passiveK + 1
		endif
	enddo
	
	if ( GetLoggingLevel(log) == DEBUG_LOGGING_LEVEL ) then
		call LogMessage(log,DEBUG_LOGGING_LEVEL,logKey,' full panels data :')
		call PrintStats(aPanels)
		call LogMessage(log,DEBUG_LOGGING_LEVEL,logKey,' activePanels data :')
		call PrintStats(activePanels)
		call LogMessage(log,DEBUG_LOGGING_LEVEL,logKey,' passivePanels data:')
		call PrintStats(passivePanels)
	endif
	call LogMessage(log,DEBUG_LOGGING_LEVEL,logKey,' GatherPanels : work done.')
end subroutine


subroutine ScatterPanels(aPanels,activePanels,activeMap,passivePanels,passiveMap)
! Recombines the activePanels and passivePanels separated by the 
! GatherPanels subroutine back into the original data structure, through the use of 
! the activeMap and passiveMap arrays.
	type(Panels), intent(inout) :: aPanels
	type(Panels), intent(in) :: activePanels, passivePanels
	integer(kint), intent(in) :: activeMap(:), passiveMap(:)
	integer(kint) :: j
	call LogMessage(log,DEBUG_LOGGING_LEVEL,logKey,'Entering ScatterPanels.')
	do j=1,activePanels%N
		call Copy(aPanels,activeMap(j),activePanels,j)
	enddo
	do j=1,passivePanels%N
		call Copy(aPanels,passiveMap(j),passivePanels,j)
	enddo
end subroutine


function TriPanelArea0(aParticles,anEdges,aPanels,panelIndex)
! Returns panel area assuming all neighbor panels are at the same resolution.  If used
! globally, grids using this function will not conserve mass.  
! Primary use for this function is as a first guess at panel area for AMR convergence criteria.
! Once grid has returned from AMR functions, users should relcalculate grid area using the 
! normal area functions.
	real(kreal) :: TriPanelArea0
	type(Particles), intent(in) :: aParticles
	type(Edges), intent(in) :: anEdges
	type(Panels), intent(in) :: aPanels
	integer(kint), intent(in) :: panelIndex
	integer(kint) :: vertIndices(3), j
	real(kreal) :: centerX(3)
	vertIndices = aPanels%vertices(:,panelIndex)
	centerX = aPanels%x(:,panelIndex)
	TriPanelArea0 = 0.0_kreal
	do j=1,2
		TriPanelArea0 = TriPanelArea0 + SphereTriArea(aParticles%x(:,vertIndices(j)),centerX,aParticles%x(:,vertIndices(j+1)))
	enddo
	TriPanelArea0 = TriPanelArea0 + SphereTriArea(aParticles%x(:,vertIndices(3)),centerX,aParticles%x(:,vertIndices(1)))
end function


subroutine ResetGridArea(aParticles,anEdges,aPanels)
! Resets area on the whole grid.  Must be called at each time step for divergent flow
! problems.
	type(Particles), intent(in) :: aParticles
	type(Edges), intent(in) :: anEdges
	type(Panels), intent(inout) :: aPanels
	integer(kint) :: panelKind, j
	panelKind = GetPanelKind(aPanels)
	if (panelKind == TRI_PANEL) then
		do j=1,aPanels%N
			if (.NOT. aPanels%hasChildren(j) ) then
				aPanels%area(j) = TriPanelArea(aParticles,anEdges,aPanels,j)
			else
				aPanels%area(j) = 0.0_kreal
			endif
		enddo
	elseif (panelKind == QUAD_PANEL ) then
		do j=1,aPanels%N
			if ( .NOT. aPanels%hasChildren(j) ) then
				aPanels%area(j) = QuadPanelArea(aParticles,anEdges,aPanels,j)
			else
				aPanels%area(j) = 0.0_kreal
			endif
		enddo
	endif
end subroutine 

subroutine RenormalizeGrid(aParticles,aPanels)
	type(Particles), intent(inout) :: aParticles
	type(Panels), intent(inout) :: aPanels
	integer(kint) :: j

	do j=1,aParticles%N
		aParticles%x(:,j) = aParticles%x(:,j)/sqrt(sum(aParticles%x(:,j)*aParticles%x(:,j)))	
	enddo
	
	do j=1,aPanels%N
		aPanels%x(:,j) = aPanels%x(:,j)/sqrt(sum(aPanels%x(:,j)*aPanels%x(:,j)))
	enddo
end subroutine

function QuadPanelArea0(aParticles,anEdges,aPanels,panelIndex)
! Returns panel area assuming all neighbor panels are at the same resolution.  If used
! globally, grids using this function will not conserve mass.  
! Primary use for this function is as a first guess at panel area for AMR convergence criteria.
! Once grid has returned from AMR functions, users should relcalculate grid area using the 
! normal area functions.
	real(kreal) :: QuadPanelArea0
	type(Particles), intent(in) :: aParticles
	type(Edges), intent(in) :: anEdges
	type(Panels), intent(in) :: aPanels
	integer(kint), intent(in) :: panelIndex
	integer(kint) :: vertIndices(4), j
	real(kreal) :: centerX(3)
	vertIndices = aPanels%vertices(:,panelIndex)
	centerX = aPanels%x(:,panelIndex)
	QuadPanelArea0 = 0.0_kreal
	do j=1,3
		QuadPanelArea0 = QuadPanelArea0 + SphereTriArea(aParticles%x(:,vertIndices(j)),centerX,aParticles%x(:,vertIndices(j+1)))
	enddo
	QuadPanelArea0 = QuadPanelArea0 + SphereTriArea(aParticles%x(:,vertIndices(4)),centerX,aParticles%x(:,vertIndices(1)))
end function


function TriPanelArea(aParticles,anEdges,aPanels,panelIndex)
! This subroutine computes the area of a triangular panel, whose neighbors may have been
! divided further than this one.
! 
!  This is the primary area function for triangular grids that use AMR.  Use of this function
! guarantees mass conservation between remeshings.
	real(kreal) :: TriPanelArea
	type(Particles), intent(in) :: aParticles
	type(Edges), intent(in) :: anEdges
	type(Panels), intent(in) :: aPanels
	integer(kint), intent(in) :: panelIndex
	! local variables
	type(IntegerListNode), pointer :: edge1list, edge2list, edge3list, current
	integer(kint), allocatable :: vertList1(:), vertList2(:), vertList3(:)
	integer(kint) :: j, parentEdges(3), verts(2)
	real(kreal) :: centerX(3)
	
	centerX = aPanels%x(:,panelIndex)
	parentEdges = aPanels%edges(:,panelIndex)
	TriPanelArea = 0.0_kreal
	! initialize empty edge lists
	nullify(edge1list)
	nullify(edge2list)
	nullify(edge3list)
	call New(edge1list)
	call New(edge2list)
	call New(edge3list)
	! Find leaf edges for all parent edges
	call FindLeafEdges(edge1list,anEdges,parentEdges(1))
	call FindLeafEdges(edge2list,anEdges,parentEdges(2))
	call FindLeafEdges(edge3list,anEdges,parentEdges(3))
	! Allocate memory for vertex lists
	allocate(vertList1(edge1list%listSize+1))
	allocate(vertList2(edge2list%listSize+1))
	allocate(vertList3(edge3list%listSize+1))	
	! populate vertex lists
	current => edge1list
	do j=1,edge1list%listSize
		verts = anEdges%verts(:,current%int)
		vertList1(j) = verts(1)
		current => current%next
	enddo
	vertList1(edge1list%listSize+1) = verts(2)
	current => edge2list
	do j=1,edge2list%listSize
		verts = anEdges%verts(:,current%int)
		vertList2(j) = verts(1)
		current=>current%next
	enddo
	vertList2(edge2list%listSize+1) = verts(2)
	current => edge3list
	do j=1,edge3list%listSize
		verts = anEdges%verts(:,current%int)
		vertList3(j) = verts(1)
		current=>current%next
	enddo
	vertList3(edge3list%listSize+1) = verts(2)
	! Calculate area
	do j=1,(size(vertList1)-1)
		TriPanelArea = TriPanelArea + SphereTriArea(aparticles%x(:,vertList1(j)),centerX,aparticles%x(:,vertList1(j+1)))
	enddo
	do j=1,(size(vertList2)-1)
		TriPanelArea = TriPanelArea + SphereTriArea(aparticles%x(:,vertList2(j)),centerX,aparticles%x(:,vertList2(j+1)))
	enddo
	do j=1,(size(vertList3)-1)
		TriPanelArea = TriPanelArea + SphereTriArea(aparticles%x(:,vertList3(j)),centerX,aparticles%x(:,vertList3(j+1)))
	enddo
	! clean up
	deallocate(vertList1)
	deallocate(vertList2)
	deallocate(vertList3)
	call Delete(edge1list)
	call Delete(edge2list)
	call Delete(edge3list)
end function


function QuadPanelArea(aParticles,anEdges,aPanels,panelIndex)
! Calculates the area of a quadlilateral panel whose neighbors may be at a different level of
! refinement.
! 
!  This is the primary area function for quadrilateral grids that use AMR.  Use of this function
! guarantees mass conservation between remeshings.
	! Calling parameters
	real(kreal) :: QuadPanelArea
	type(Particles), intent(in) :: aParticles
	type(Edges), intent(in) :: anEdges
	type(Panels), intent(in) :: aPanels
	integer(kint), intent(in) :: panelIndex
	! Local variables
	type(IntegerListNode), pointer :: edge1list, edge2list, edge3list, edge4list, current
	integer(kint), allocatable :: vertlist1(:), vertList2(:), vertList3(:), vertList4(:)
	integer(kint) :: j, parentEdges(4), verts(2)
	real(kreal) :: centerX(3)
	centerX = aPanels%x(:,panelIndex)
	QuadPanelArea = 0.0_kreal
	! Initialize empty edges lists
	nullify(edge1list)
	nullify(edge2list)
	nullify(edge3list)
	nullify(edge4list)
	call New(edge1list)
	call New(edge2list)
	call New(edge3list)
	call New(edge4list)
	! Find leaf edges for all four parent edges
	parentEdges = aPanels%edges(:,panelIndex)
	call FindLeafEdges(edge1list,anEdges,parentEdges(1))
	call FindLeafEdges(edge2list,anEdges,parentEdges(2))
	call FindLeafEdges(edge3list,anEdges,parentEdges(3))
	call FindLeafEdges(edge4list,anEdges,parentEdges(4))
	! allocate memory for vertex lists associated with leaf edges
	allocate(vertList1(edge1list%listSize+1))
	allocate(vertList2(edge2list%listSize+1))
	allocate(vertList3(edge3list%listSize+1))
	allocate(vertList4(edge4list%listSize+1))
	! populate vertex lists
	current => edge1list
	do j=1,edge1list%listSize
		verts = anEdges%verts(:,current%int)
		vertList1(j) = verts(1)
		current=>current%next
	enddo
	vertList1(edge1list%listSize+1) = verts(2)
	current => edge2list
	do j=1,edge2list%listSize
		verts = anEdges%verts(:,current%int)
		vertList2(j) = verts(1)
		current=>current%next
	enddo
	vertList2(edge2list%listSize+1) = verts(2)
	current => edge3list
	do j=1,edge3list%listSize
		verts = anEdges%verts(:,current%int)
		vertList3(j) = verts(1)
		current=>current%next
	enddo
	vertList3(edge3list%listSize+1) = verts(2)
	current => edge4list
	do j=1,edge4list%listSize
		verts = anEdges%verts(:,current%int)
		vertList4(j) = verts(1)
		current=>current%next
	enddo
	vertList4(edge4list%listSize+1) = verts(2)

	do j=1,(size(vertList1)-1)
		QuadPanelArea = QuadPanelArea + SphereTriArea(aparticles%x(:,vertList1(j)),centerX,aparticles%x(:,vertList1(j+1)))
	enddo
	do j=1,(size(vertList2)-1)
		QuadPanelArea = QuadPanelArea + SphereTriArea(aparticles%x(:,vertList2(j)),centerX,aparticles%x(:,vertList2(j+1)))
	enddo
	do j=1,(size(vertList3)-1)
		QuadPanelArea = QuadPanelArea + SphereTriArea(aparticles%x(:,vertList3(j)),centerX,aparticles%x(:,vertList3(j+1)))
	enddo
	do j=1,(size(vertList4)-1)
		QuadPanelArea = QuadPanelArea + SphereTriArea(aparticles%x(:,vertList4(j)),centerX,aparticles%x(:,vertList4(j+1)))
	enddo
	! Clean up
	deallocate(vertList1)
	deallocate(vertList2)
	deallocate(vertList3)
	deallocate(vertList4)
	call Delete(edge1list)
	call Delete(edge2list)
	call Delete(edge3list)
	call Delete(edge4list)
end function


function DualPanelArea(aParticles,anEdges,aPanels,particleIndex)
! Calculates the area of a dual mesh panel centered on the particle at index particle index.
! Note that this routine allows for polygons of varying kinds (triangles, quadrilaterals, 
! pentagons, etc.)
	real(kreal) :: DualPanelArea
	type(Particles), intent(in) :: aParticles
	type(Edges), intent(in) :: anEdges
	type(Panels), intent(in) :: aPanels
	integer(kint), intent(in) :: particleIndex
	integer(kint) :: j, triDualList(6), quadDualList(4), panelCount, panelKind
	triDualList = 0
	quadDualList = 0
	panelCount = 0
	panelKind = GetPanelKind(aPanels)
	DualPanelArea = 0.0_kreal
	if ( panelKind == TRI_PANEL) then
		call CCWDualVertList(triDualList,aparticles,anEdges,aPanels,particleIndex)
		do j=1,6
			if (triDualList(j) > 0) panelCount = panelCount + 1
		enddo
		do j=1,panelCount -1
			DualPanelArea = DualPanelArea + SphereTriArea(aPanels%x(:,triDualList(j)),&
							aParticles%x(:,particleIndex),aPanels%x(:,triDualList(j+1)))
		enddo
		DualPanelArea = DualPanelArea + SphereTriArea(aPanels%x(:,triDualList(panelCount)),&
								aParticles%x(:,particleIndex),aPanels%x(:,triDualList(1)))
	elseif (panelKind == QUAD_PANEL) then 
		call CCWDualVertList(quadDualList,aParticles,anEdges,aPanels,particleIndex)
		do j=1,4
			if (quadDualList(j) > 0) panelCount = panelCount + 1
		enddo
		do j=1,panelCount-1
			DualPanelArea = DualPanelArea + SphereTriArea(apanels%x(:,quadDualList(j)),&
				aParticles%x(:,particleIndex),aPanels%x(:,quadDualList(j+1)))
		enddo
		DualPanelArea = DualPanelArea + SphereTriArea(aPanels%x(:,quadDualList(panelCount)),&
					aParticles%x(:,particleIndex),aPanels%x(:,quadDualList(1)))
	endif
end function


subroutine CCWDualVertList(vertList, aParticles, anEdges, aPanels, particleIndex)
! Returns the list of panel centers adjacent to a particle in CCW order.
! Necessary for calculating dual panel areas and handling polygons of differing types.
	integer(kint), intent(out) :: vertList(:)
	type(Particles), intent(in) :: aParticles
	type(Edges), intent(in) :: anEdges
	type(Panels), intent(in) :: aPanels
	integer(kint), intent(in) :: particleIndex
	integer(kint) :: panelKind, j, k, quadEdgeList(4), triEdgeList(6)
	integer(kint) :: startEdge, currentEdge, currentPanel, startPanel, nextEdge, nextPanel
	integer(kint) :: vertNum, edgeNum, edgeCounter, panelK
	panelKind = GetPanelKind(aPanels)
	if (panelKind == TRI_PANEL) then
		if ( size(vertList) /= 6 ) then
			print *, " CCWDualVertList ERROR : vertList must have size = 6."
			return
		endif
		vertList = 0
		triEdgeList = aParticles%edges(:,particleIndex)
		startEdge = triEdgeList(1)
		! choose first edge for starting point
		if ( anEdges%verts(1,triEdgeList(1)) == particleIndex) then
			vertList(1) = anEdges%leftPanel(triEdgeList(1))
		else
			vertList(1) = anEdges%rightPanel(triEdgeList(1))
		endif
		startPanel = vertList(1)
		currentEdge = startEdge
		currentPanel = startPanel
		do j=2,6
			do k=1,3
				if (aPanels%vertices(k,currentPanel) == particleIndex) then
					vertNum = k
				endif
				if (aPanels%edges(k,currentPanel) == currentEdge) then
					edgeNum = k
				endif
			enddo
			nextEdge = NextTriEdgeCCW(vertNum,edgeNum)
			nextEdge = aPanels%edges(nextEdge,currentPanel)
			edgeCounter = 0
			do k=1,6
				if ( triEdgeList(k) == nextEdge) then
					edgeCounter = edgeCounter + 1
				endif
			enddo
			if (edgeCounter /= 1) then
				print *,"CCWDualVertList ERROR : next edge not found at current vertex."
				print *,"edgeCounter = ",edgeCounter
				print *,"Edge index = ",nextEdge
				print *,"Particle index = ",particleIndex
				print *,"Particle edges = ",aParticles%edges(:,particleIndex)
				print *,"TriEdgeList = ",triEdgeList
				return
			else
				if (anEdges%verts(1,nextEdge) == particleIndex) then
					nextPanel = anEdges%leftPanel(nextEdge)
				else
					nextPanel = anEdges%rightPanel(nextEdge)
				endif
			endif
			if ( nextPanel == startPanel) then
				exit
			else
				vertList(j) = nextPanel
				currentEdge = nextEdge
				currentPanel = nextPanel
			endif
		enddo
	elseif (panelKind == QUAD_PANEL) then
		if ( size(vertList) /= 4) then
			print *, "CCWDualvertList ERROR : vertList must have size = 4."
			return
		endif
		vertList = 0
		quadEdgeList = aParticles%edges(:,particleIndex)
		panelK = 1
		do k=1,4
			if (quadEdgeList(k) /= 0) then
				if (anEdges%verts(1,quadEdgeList(k)) == particleIndex) then
					vertList(panelK) = anEdges%leftPanel(quadEdgeList(k))
				else
					vertList(panelK) = anEdges%rightPanel(quadEdgeList(k))
				endif
				panelK = panelK + 1
			endif
		enddo
	endif
end subroutine


function NextTriEdgeCCW(currVert,currEdge)
! Used by CCWDualVertList.  Returns a pointer to the next CCW panel in a triangular mesh.
! Outputs local triangle information, valid relative to each panel independently.  
	integer(kint) :: NextTriEdgeCCW
	integer(kint), intent(in) :: currVert, currEdge
	if (currVert == 1) then
		if ( currEdge == 1) then
			NextTriEdgeCCW = 3
		else
			NextTriEdgeCCW = 1
		endif
	elseif (currVert == 2) then
		if ( currEdge == 1) then
			NextTriEdgeCCW = 2
		else
			NextTriEdgeCCW = 1
		endif
	elseif (currVert == 3) then
		if ( currEdge == 2 ) then
			NextTriEdgeCCW = 3
		else
			NextTriEdgeCCW = 2
		endif
	else
		print *,"NextTriEdgeCCW ERROR : bad vertex input."
		NextTriEdgeCCW = -1
	endif	
end function

!----------------
! Output Methods : Plotting and console
!----------------
subroutine Printparticlesstats(self)
	type(Particles), intent(in) :: self
	print *,"PARTICLES Stats:"
	print *, "  N_Max = ",self%N_Max
	print *, "  N     = ",self%N
	print *, "  Dual mesh surf. area = ",sum(self%dualArea)
	if ( associated(self%tracer)) then
		print *,"  Max tracer1 = ",maxval(self%tracer(1:self%N,1))
		print *,"  Min tracer1 = ",minval(self%tracer(1:self%N,1))
	endif
	if ( associated(self%relVort)) then
		print *,"  Max relVort = ",maxval(self%relVort(1:self%N))
		print *,"  Min relVort = ",minval(self%relVort(1:self%N))
	endif
	if ( associated(self%absVort)) then
		print *,"  Max absVort = ",maxVal(self%absVort(1:self%N))
		print *,"  Min absVort = ",minval(self%absVort(1:self%N))
	endif
end subroutine


subroutine PrintPanelsStats(self)
	type(Panels), intent(in) :: self
	print *,"PANELS Stats : "
	print *,"  N_Max    = ",self%N_Max
	print *,"  N        = ",self%N
	print *,"  N_Active = ",self%N_Active
	print *,"  MaxNest  = ",maxVal(self%nest(1:self%N))
	print *,"  Surface Area = ", sum(self%area(1:self%N))
	if ( associated(self%tracer)) then
		print *,"  Max tracer1 = ",maxval(self%tracer(1:self%N,1))
		print *,"  Min tracer1 = ",minval(self%tracer(1:self%N,1))
		print *,"  tracer integral = ",sum(self%tracer(1:self%N,1)*self%area(1:self%N))
	endif
	if ( associated(self%relVort)) then
		print *,"  Max relVort = ",maxval(self%relVort(1:self%N))
		print *,"  Min relVort = ",minval(self%relVort(1:self%N))
		print *,"  relvort integral = ",sum(self%relVort*self%area)
	endif
	if ( associated(self%absVort)) then
		print *,"  Max absVort = ",maxVal(self%absVort(1:self%N))
		print *,"  Min absVort = ",minval(self%absVort(1:self%N))
		print *,"  absVort integral = ",sum(self%absVort(1:self%N)*self%area(1:self%N))
	endif
end subroutine


subroutine PrintEdgesStats(self)
	type(Edges), intent(in) :: self
	print *,"EDGES Stats : "
	print *, "   N     = ",self%N
	print *, "   N_Max = ",self%N_Max
end subroutine


subroutine vtkOutputUnif(aParticles,aPanels,filename)
	! Outputs a vtkPolyData file suitable for input into VTK C++ programs or 
	! ParaView.  
	type(Particles), intent(in) :: aParticles
	type(Panels), intent(in) :: aPanels
	character(len=*), intent(in) :: filename
	! Local variables
	integer(KINT) :: cellListSize,errorCode, nScalars, nPoints, j,k, nTracers, nActive, nParticles
	integer(kint) :: nVerts, panelKind
	integer(KINT), parameter :: unitNumber = 11
	character(len=28) :: dataString
	! Format statements
	103 format (A)
	104 format (A,I8,A)
	105 format (F24.18, F24.18, F24.18)
	106 format (A, I8, I8)
	107 format (I8,3I8)
	108 format (F24.18)
	109 format (A,I8)
	110 format (I8,4I8)
	! VTK data size integers
	nTracers = 0
	nActive = aPanels%N_Active
	nParticles = aParticles%N
	if (associated(aParticles%tracer)) then
		nTracers = size(aParticles%tracer,2)
	endif
	nScalars = nTracers
	if ( associated(aParticles%relVort)) then
		nScalars = nScalars + 1
	endif
	if ( associated(aparticles%absVort)) then
		nScalars = nScalars + 1
	endif
	panelKind = GetPanelKind(aPanels)	
	if ( panelKind /= VORONOI_PANEL) then
		cellListSize = (panelKind+1)*nActive
	else
		cellListSize = 0
		do j=1,aPanels%N
			do k=1,10
				if (aPanels%vertices(k,j) > 0 ) cellListSize = cellListSize + 1
			enddo
		enddo
		cellListSize = cellListSize + 2*nactive
	endif
	nPoints = aParticles%N
	
	open(unit=unitNumber,file=filename,status='REPLACE',action='WRITE',iostat=errorCode)
		if (errorCode /= 0 ) stop 'vtkOutput ERROR: cannot open file.'
	! VTK HEADER
	write(unitNumber,103) '# vtk DataFile Version 2.0'
	write(unitNumber,104) 'Sphere Grid : ',nActive,' panels.'
	write(unitNumber,103) 'ASCII'
	write(unitNumber,103) 'DATASET POLYDATA'
	
	! VTK POINTS
	write(unitNumber,104) 'POINTS ',nParticles,'   double'
	do j=1,nParticles
		write(unitNumber,105) aParticles%x(1,j), aParticles%x(2,j), aParticles%x(3,j)
	enddo
	
	! VTK CELLS
	write(unitNumber,106) 'POLYGONS', nActive, cellListSize
	if (GetPanelKind(aPanels) == TRI_PANEL) then
		do j=1,aPanels%N
			if (.NOT. aPanels%hasChildren(j)) then
				write(unitNumber,107) 3, aPanels%vertices(:,j)-1
			endif
		enddo
	elseif (GetPanelKind(aPanels) == QUAD_PANEL) then
		do j=1,aPanels%N
			if (.NOT. aPanels%hasChildren(j)) then
				write(unitNumber,110) 4, aPanels%vertices(:,j)-1
			endif
		enddo
	elseif ( GetPanelKind(aPanels) == VORONOI_PANEL) then
		do j=1,aPanels%n
			nVerts = 0
			do k=1,VORONOI_PANEL
				if ( aPanels%vertices(k,j) >0 ) nVerts = nVerts + 1
			enddo
			write(unitNumber,*) nVerts,aPanels%vertices(1:nVerts,j) - 1
		enddo	
	endif
	! POINTS SCALARS	
	write(unitNumber,109) 'POINT_DATA   ',nPoints
	if (nTracers > 0) then
		do k=1,nTracers
			write(dataString,'(A6,I1)') trim('Tracer'),k
			write(unitNumber,*) 'SCALARS  ',trim(dataString),'  double 1'
			write(unitNumber,103) 'LOOKUP_TABLE default'
			do j=1,nParticles
				write(unitNumber,108) aParticles%tracer(j,k)
			enddo
		enddo
	endif	
	if ( associated(aParticles%relVort)) then
		write(unitNumber,103) 'SCALARS   relVort  double 1'
		write(unitNumber,103) 'LOOKUP_TABLE default'
		do j=1,nParticles
			write(unitNumber,108) aparticles%relVort(j)
		enddo	
	endif
	if ( associated(aParticles%absVort)) then
		write(unitNumber,103) 'SCALARS   absVort  double 1'
		write(unitNumber,103) 'LOOKUP_TABLE default'
		do j=1,nParticles
			write(unitNumber,108) aparticles%absVort(j)
		enddo	
	endif
	! CELL SCALARS
	write(unitNumber,109) 'CELL_DATA   ',nActive
	if ( associated(aPanels%tracer)) then
		do k=1,nTracers
			write(dataString,'(A,I1)') trim('TracerPanel'),k
			write(unitNumber,'(A,A,A)') 'SCALARS ',trim(dataString),' double 1'
			write(unitNumber,103) 'LOOKUP_TABLE default'
			do j=1,aPanels%N
				if ( .NOT. aPanels%hasChildren(j)) then
					write(unitNumber,108) aPanels%tracer(j,k)		
				endif
			enddo
		enddo
	endif
	if ( associated(aPanels%relVort) ) then
		write(unitNumber,103) 'SCALARS  relVortPanel  double 1'
		write(unitNumber,103) 'LOOKUP_TABLE default'
		do j=1,aPanels%N
			if ( .NOT. aPanels%hasChildren(j)) then
				write(unitNumber,108) aPanels%relVort(j)
			endif
		enddo
	endif
	if ( associated(aPanels%absVort) ) then
		write(unitNumber,103) 'SCALARS  absVortPanel  double 1'
		write(unitNumber,103) 'LOOKUP_TABLE default'
		do j=1,aPanels%N
			if ( .NOT. aPanels%hasChildren(j)) then
				write(unitNumber,108) aPanels%absVort(j)
			endif
		enddo
	endif
	close(unitNumber)
end subroutine


subroutine vtkOutputActivePassive(aParticles,anEdges,aPanels,passiveFilename,activeFilename)
	type(Particles), intent(in) :: aParticles
	type(Edges), intent(in) :: anEdges
	type(Panels), intent(in) :: aPanels
	character(len=*), intent(in) :: activeFilename, passiveFilename
	integer(kint) :: nTracer, i,j,k, nPoints, errCode
	integer(kint), parameter :: unitNumber=15
	character(len=28) :: dataString
	
	call vtkOutputUnif(aParticles,aPanels,passiveFilename)
	
	nTracer = GetNTracer(aPanels)
	nPoints = aPanels%N_Active
	
	open(unit=unitNumber,file=activeFilename,action='WRITE',status='REPLACE',iostat=errCode)
		if (errCode /= 0 ) stop 'vtkOutput ERROR: cannot open file.'
	! VTK HEADER
	write(unitNumber,'(A)') '# vtk DataFile Version 2.0'
	write(unitNumber,'(A,I8,A)') 'Sphere Grid ActivePanels : ',nPoints,' panels.'
	write(unitNumber,'(A)') 'ASCII'
	write(unitNumber,'(A)') 'DATASET POLYDATA'
	
	! VTK POINTS
	write(unitNumber,'(A,I8,A)') 'POINTS ',nPoints,' double'
	do j=1,aPanels%N
		if (.NOT. aPanels%hasChildren(j) ) then
			write(unitNumber,'(3F24.15)') aPanels%x(:,j)
		endif
	enddo
	
	! VTK TOPOLOGY
	write(unitNumber,'(A,I8,I8)') 'VERTICES ', nPoints, 2*nPoints
	do j=1,aPanels%N
		if ( .NOT. aPanels%hasChildren(j) ) then
			write(unitNumber,'(I4,I8)') 1, j-1
		endif	
	enddo
	
	! POINTS SCALARS
	write(unitNumber,'(A,I8)') 'POINT_DATA ',nPoints
	if ( associated(aParticles%relVort) ) then
		write(unitNumber,'(A)') 'SCALARS relVort double 1'
		write(unitNumber,'(A)') 'LOOKUP_TABLE default'
		do j=1,aPanels%N
			if ( .NOT. aPanels%hasChildren(j) ) then
				write(unitNumber,'(F24.15)') aPanels%relVort(j)
			endif
		enddo
	endif
	if ( associated(aParticles%absVort) ) then
		write(unitNumber,'(A)') 'SCALARS absVort double 1'
		write(unitNumber,'(A)') 'LOOKUP_TABLE default'
		do j=1,aPanels%N
			if ( .NOT. aPanels%hasChildren(j) ) then
				write(unitNumber,'(F24.15)') aPanels%absVort(j)
			endif
		enddo
	endif
	if ( nTracer > 0 ) then
		do i=1,nTracer
			write(dataString,'(A,I1)') trim('Tracer'),i
			write(unitNumber,'(A,A,A)') 'SCALARS ',trim(dataString),' double 1'
			write(unitNumber,'(A)') 'LOOKUP_TABLE default'
			do j=1,aPanels%N
				if (.NOT. aPanels%hasChildren(j) ) then
					write(unitNumber,'(F24.15)') aPanels%tracer(j,i)
				endif
			enddo
		enddo
	endif	
	close(unitNumber)
	
end subroutine


subroutine vtkOutputAMR(aParticles,anEdges, aPanels,filename)
	! Outputs a vtkPolyData file suitable for input into VTK C++ programs or 
	! ParaView.  Use this routine for data sets with adaptive meshes.
	! Arguments
	type(Particles), intent(in) :: aParticles
	type(Edges), intent(in) :: anEdges
	type(Panels), intent(in) :: aPanels
	character(len=*), intent(in) :: filename
	! Locals
	integer(kint) :: nTracer, errCode, j, k, nPoints, cellListSize
	integer(kint) :: panelKind, nActive, nParticles
	integer(kint), parameter :: unitNumber = 11, MAX_EDGES = 20
	character(len=28) :: dataString
	type(IntegerListNode), pointer :: leafEdgeList
	integer(kint) :: numberOfEdges, edgeList(MAX_EDGES)
	
	nParticles = aParticles%N
	nActive = aPanels%N_Active
	panelKind = GetPanelKind(aPanels)
	nTracer = GetNTracer(aPanels)
	
	! Determine cellListSize
	cellListSize = 0
	do j=1,aPanels%N
		if (.NOT. aPanels%hasChildren(j) ) then
			call New(leafEdgeList)
			do k=1,panelKind
				call FindLeafEdges(leafEdgeList,anEdges,aPanels%edges(k,j))
			enddo
			cellListSize = cellListSize + leafEdgeList%listSize + 1
			call Delete(leafEdgeList)	
		endif
	enddo
	
	open(unit=unitNumber,file=filename,status='REPLACE',action='WRITE',iostat=errCode)
		if ( errCode /= 0 ) then
			call LogMessage(log,ERROR_LOGGING_LEVEL,logKey,'vtkOutput ERROR: cannot open .vtk file')
			return
		endif
		! .vtk File HEADER
		write(unitNumber,'(A)') '# vtk DataFile Version 2.0'
		write(unitNumber,'(A,I8,A)') 'Sphere Grid : ',nActive,' panels.'
		write(unitNumber,'(A)') 'ASCII'
		write(unitNumber,'(A)') 'DATASET POLYDATA'
		
		! VTK POINTS
		write(unitNumber,'(A,I8,A)') 'POINTS ',nParticles,' double'
		do j=1,nParticles
			write(unitNumber,'(F24.18,F24.18,F24.18)') aParticles%x(1,j), aParticles%x(2,j), aParticles%x(3,j)
		enddo
		
		! VTK CELLS
		write(unitNumber,'(A,I8,I8)') 'POLYGONS ', nActive, cellListSize
		do j=1,aPanels%N
			if (.NOT. aPanels%hasChildren(j) ) then
				edgeList = -1
				call New(leafEdgeList)
				do k=1,panelKind
					call FindLeafEdges(leafEdgeList,anEdges,aPanels%edges(k,j))
				enddo
				numberOfEdges = leafEdgeList%listSize
				if ( numberOfEdges > MAX_EDGES) then
					call LogMessage(log,ERROR_LOGGING_LEVEL,logKey,'vtkOutput ERROR : found panel with too many edges')
					return
				endif
				call ReturnIntegerArrayFromList(edgeList(1:numberOfEdges),leafEdgeList)
				
				write(unitNumber,'(I6)',advance='NO') numberOfEdges
				
				! Generate CCW vertex list from edge list
				do k=1,numberOfEdges
					if ( anEdges%leftPanel(edgeList(k)) == j) then
						write(unitNumber,'(I8)',advance='NO') anEdges%verts(1,edgeList(k)) - 1
					else
						write(unitNumber,'(I8)',advance='NO') anEdges%verts(2,edgeList(k)) - 1
					endif
				enddo
				write(unitNumber,'(A)',advance='YES') ' '
				call Delete(leafEdgeList)
			endif
		enddo
		
		write(unitNumber,'(A,I8)') 'POINT_DATA ',nParticles
		if ( associated(aParticles%relVort) ) then
			write(unitNumber,'(A)') 'SCALARS relVort double 1'
			write(unitNumber,'(A)') 'LOOKUP_TABLE default'
			do j=1,nParticles
				write(unitNumber,'(F24.18)') aParticles%relVort(j)
			enddo
		endif
		if ( associated(aParticles%absVort) ) then
			write(unitNumber,'(A)') 'SCALARS absVort double 1'
			write(unitNumber,'(A)') 'LOOKUP_TABLE default'
			do j=1,nParticles
				write(unitNumber,'(F24.18)') aParticles%absVort(j)
			enddo
		endif
		if (nTracer > 0 ) then
			do k=1,nTracer
				write(dataString,'(A,I1)') 'Tracer',k
				write(unitNumber,'(A,A,A)') 'SCALARS ',trim(dataString),' double 1'
				write(unitNumber,'(A)') 'LOOKUP_TABLE default'
				do j=1,nParticles
					write(unitNumber,'(F24.18)') aParticles%tracer(j,k)
				enddo
			enddo
		endif
		
		write(unitNumber,'(A,I8)') 'CELL_DATA ',nActive
		if ( associated(aPanels%relVort) ) then
			write(unitNumber,'(A)') 'SCALARS relVortPanel double 1'
			write(unitNumber,'(A)') 'LOOKUP_TABLE default'
			do j=1,aPanels%N
				if ( .NOT. aPanels%hasChildren(j) ) then
					write(unitNumber,'(F24.18)') aPanels%relVort(j)
				endif
			enddo
		endif
		if ( associated(aPanels%absVort) ) then
			write(unitNumber,'(A)') 'SCALARS absVortPanel double 1'
			write(unitNumber,'(A)') 'LOOKUP_TABLE default'
			do j=1,aPanels%N
				if ( .NOT. aPanels%hasChildren(j) ) then
					write(unitNumber,'(F24.18)') aPanels%absVort(j)
				endif
			enddo
		endif
		if ( nTracer > 0 ) then
			do k=1,nTracer
				write(dataString,'(A,I1)') 'TracerPanel',k
				write(unitNumber,'(A,A,A)') 'SCALARS ',trim(dataString),' double 1'
				write(unitNumber,'(A)') 'LOOKUP_TABLE default'
				do j=1,aPanels%N
					if ( .NOT. aPanels%hasChildren(j) ) then
						write(unitNumber,'(F24.18)') aPanels%tracer(j,k)
					endif
				enddo
			enddo
		endif
		
		close(unitNumber)
end subroutine


!----------------
! Other Methods :
!----------------
function GetPanelKind(aPanels)
! Returns the panelKind variable from a panels data structure.
	type(Panels), intent(in) :: aPanels
	integer(kint) :: GetPanelKind
	GetPanelKind = size(apanels%vertices,1)
end function


function NTracerParticles(aParticles)
! Returns the number of tracers carried by a grid.
	type(Particles), intent(in) :: aParticles
	integer(kint) :: NTracerParticles
	if (associated(aParticles%tracer)) then
		NTracerParticles = size(aParticles%tracer,2)
	else
		NTracerParticles = 0
	endif
end function


function NTracerPanels(aPanels)
! Returns the number of tracers carried by a grid.
	type(Panels), intent(in) :: aPanels
	integer(kint) :: NTracerPanels
	if (associated(aPanels%tracer)) then
		NTracerPanels = size(aPanels%tracer,2)
	else
		NTracerPanels = 0
	endif
end function


function ParticlesProblemKind(aParticles)
! Returns the problemKind associated with a grid.
	integer(kint) :: ParticlesProblemKind
	type(Particles), intent(in) :: aParticles
	if ( associated(aParticles%absVort) ) then
		ParticlesProblemKind = BVE_SOLVER
	else
		ParticlesProblemKind = ADVECTION_SOLVER
	endif
end function


function PanelsProblemKind(aPanels)
! Returns the problemKind associated with a grid.
	integer(kint) :: PanelsProblemKind
	type(Panels), intent(in) :: aPanels
	if ( associated(aPanels%absVort) ) then
		PanelsProblemKind = BVE_SOLVER
	else
		PanelsProblemKind = ADVECTION_SOLVER
	endif
end function


function CubedSphereMaxMeshSize( nestLevel)
! Returns the maximum edge length of a cubed sphere grid recursively refined down to an initial 
! nest level nestLevel.
!
! Mesh size output is given in degrees
!
	! Calling parameters
	integer(KINT):: nestLevel
	real(KREAL) :: cubedSphereMaxMeshSize
	! Local variables
	integer :: p
	
	p = 4 - nestLevel
	cubedSphereMaxMeshSize = 5.6464d0*2**p	
end function


function IcosTriMaxMeshSize( nestLevel )
! Returns the maximum edge length of an icosahedral triangle grid recursively refined down
! to nestLevel.
!
! Mesh size output is given in degrees
!
	! Calling parameters
	integer(KINT) :: nestLevel
	real(KREAL) :: icosTriMaxMeshSize
	! Local variables
	integer :: p
	
	p = 4 - nestLevel
	icosTriMaxMeshSize = 4.7342d0*2**p	
end function


function QuadAvgMeshSize(initNest)
! Returns the average edge length of a quadrilateral mesh refined to nest level initNest.
	real(kreal) :: QuadAvgMeshSize
	integer(kint), intent(in) :: initNest
	if ( initNest >= 3) then
		QuadAvgMeshSize = 5.20_kreal*2.0_kreal**(4-initNest)
	elseif (initNest == 1) then
		QuadAvgMeshSize = 40.13_kreal
	elseif (initNest == 2) then
		QuadAvgMeshSize = 20.62_kreal
	elseif ( initNest == 0) then
		QuadAvgMeshSize = 70.13_kreal
	else
		QuadAvgMeshSize = 0.0_kreal
		stop "ERROR : initNest must be nonnegative."
	endif
end function


function TriAvgMeshSize(initNest)
! Returns the average edge length of an icosahedral triangular mesh refined to nestLevel 
! initNest.
	real(kreal) :: TriAvgMeshSize
	integer(kint), intent(in) :: initNest
	if ( initNest >= 3 ) then
		TriAvgMeshSize = 4.33_kreal*2.0_kreal**(4-initNest)
	elseif ( initNest == 2) then
		TriAvgMeshSize = 17.22_kreal
	elseif (initNest == 1) then
		TriAvgMeshSize = 33.89_kreal
	elseif (initNest == 0) then
		TriAvgMeshSize = 63.44_kreal
	else
		TriAvgMeshSize = 0.0_kreal
		stop "ERROR : initNest must be nonnegative."
	endif
end function


function PanelMax(panelKind,maxNest)
! Returns the maximum number of panels needed at a given nest for memory allocation.
	integer(kint) :: panelMax
	integer(kint), intent(in) :: panelKind, maxNest
	integer(kint) :: j
	panelMax = 0
	if ( panelKind == TRI_PANEL) then
		do j=0,maxNest
			panelMax = panelMax + 20*4**j
		enddo
	elseif (panelKind == QUAD_PANEL) then
		do j=0,maxNest
			panelMax = panelMax + 6*4**j
		enddo
	endif
	
end function


function ParticleMax(panelKind,maxNest)
! Returns the maximum number of particles needed for memory allocation.
	integer(kint) :: particleMax
	integer(kint), intent(in) :: panelKind, maxNest
	if ( panelKind == TRI_PANEL ) then
		particleMax = 2 + 10*4**maxNest
	elseif (panelKind == QUAD_PANEL ) then
		particleMax = 2 + 6*4**maxNest	
	endif
end function


function EdgeMax(panelKind,maxNest)
! Returns the maximum number of edges needed for memory allocation.
	integer(kint) :: edgeMax
	integer(kint), intent(in) :: panelKind, maxNest
	integer(kint) :: j
	edgeMax = 0
	if ( panelKind == TRI_PANEL ) then
		do j=0,maxNest
			edgeMax = edgeMax + 30*4**j
		enddo	
	elseif (panelKind == QUAD_PANEL ) then
		do j=0,maxNest
			edgeMax = edgeMax + 12*4**j
		enddo
	elseif (panelKind == VORONOI_PANEL ) then
		EdgeMax = 6*PanelMax(panelKind,maxNest) - 12
	endif
end function


recursive subroutine FindLeafEdges(leafList,anEdges,edgeIndex)
! Returns an integer linked list of leaf edges descended from a parent edge.
	type(IntegerListNode), pointer, intent(inout) :: leafList
	type(Edges), intent(in) :: anEdges
	integer(kint), intent(in) :: edgeIndex
	integer(kint) :: child1, child2
	if (.NOT. anEdges%hasChildren(edgeIndex)) then	
		call AddIntegerToList(leafList,edgeIndex)
	else
		child1 = anEdges%children(1,edgeIndex)
		child2 = anEdges%children(2,edgeIndex)
		call FindLeafEdges(leafList,anEdges,child1)
		call FindLeafEdges(leafList,anEdges,child2) 
	endif
end subroutine


recursive subroutine PanelTreeSearch(nearestPanel,xyz,rootPanel,aPanels)
! Searches the tree-structure of panels for the nearest panel to input point xyz.  Note that
! this routine only sees children of its starting panel, and therefore may not converge
! to the actual closest panel on the grid.  
!
! Use this routine to provide a good starting point for the AdjacentPanelSearch subroutine.
!
	integer(kint), intent(out) :: nearestPanel
	real(kreal), intent(in) :: xyz(3)
	integer(kint), intent(in) :: rootPanel
	type(Panels), intent(in) :: aPanels
	real(kreal) :: currentMin, testDist
	integer(kint) :: j, nextRoot
	if ( aPanels%hasChildren(rootPanel) ) then
		currentMin = PI
		do j=1,4
			testDist = ChordDistance(xyz,aPanels%x(:,aPanels%children(j,rootPanel)))
			if ( testDist < currentMin ) then
				currentMin = testDist
				nextRoot = aPanels%children(j,rootPanel)
			endif
		enddo
		call PanelTreeSearch(nearestPanel,xyz,nextRoot,aPanels)
	else
		nearestPanel = rootPanel
	endif
end subroutine


recursive subroutine AdjacentPanelSearch(inPanel,xyz,startPanel,aParticles,anEdges,aPanels)
! Searches a grid for the panel that contains the point xyz, using edges to move through the 
! mesh, by searching adjacent panels beginning at startPanel.
!
! TO DO : This subroutine sometimes creates an infinite loop due to finite-precision roundoff
! errors.  Change the convergence criteria or restart search from a random starting point
! to get around this.  This is particularly a problem in AMR grids.  
!	Idea1 : Change convergence criteria to distance from centroid instead of leftTest
!				What if polygons are non-convex?
!
	integer(kint), intent(out) :: inPanel
	real(kreal), intent(in) :: xyz(3)
	integer(kint), intent(in) :: startPanel
	type(Particles), intent(in) :: aParticles
	type(Edges), intent(in) :: anEdges
	type(Panels), intent(in) :: aPanels
	real(kreal) :: currentMin, testDist, leftCheck
	type(IntegerListNode), pointer :: leafEdgeList, current, next
	integer(kint) :: currentPanel, panelKind, j, edgeIndex, edgeStartIndex, edgeEndIndex, nextPanel
	integer(kint) :: edgeCounter, numberOfEdges
	nullify(leafEdgeList)
	nullify(current)
	nullify(next)
	! check input for errors
	if ( (startPanel <= 0) .OR. (startPanel >= aPanels%N)  ) then
		print *,"AdjacentPanelSearch ERROR : invalid initial guess."
		return
	endif
	panelKind = GetPanelKind(aPanels)
	! Case 1 : panel has no children
	!print *,"DEBUG : entering AdjacentPanelSearch : startPanel = ",startPanel," newX = ",xyz
	if ( .NOT. aPanels%hasChildren(startPanel) ) then
		! Find low-level edges of this panel
		call New(leafEdgeList)
		do j=1,panelKind
			call FindLeafEdges(leafEdgeList,anEdges,aPanels%edges(j,startPanel))
		enddo
		numberOfEdges = leafEdgeList%listSize
		!print *,"DEBUG : leafEdgesListSize = ",leafEdgeList%listSize
		! use left test to determine if xyz lies in this panel
		edgeCounter = 0
		current=>leafEdgeList
		do j=1,leafEdgeList%listSize
			edgeIndex = current%int
			! put edge vertices in ccw order
			if ( startPanel == anEdges%leftPanel(edgeIndex) ) then
				edgeStartIndex = anEdges%verts(1,edgeIndex)
				edgeEndIndex = anEdges%verts(2,edgeIndex)
			else
				edgeStartIndex = anEdges%verts(2,edgeIndex)
				edgeEndIndex = anEdges%verts(1,edgeIndex)
			endif
			leftCheck = Determinant(xyz,aParticles%x(:,edgeStartIndex),aParticles%x(:,edgeEndIndex))
			!print *, "DEBUG : leftCheck = ",leftCheck
			if ( leftCheck < -0.001_kreal) then
				! flip the edge, and restart search
				if ( startPanel == anEdges%leftPanel(edgeIndex) ) then
					nextPanel = anEdges%rightPanel(edgeIndex)
				else
					nextPanel = anEdges%leftPanel(edgeIndex)
				endif
				call Delete(leafEdgeList)
				nullify(leafEdgeList)
				!print *, "DEBUG returned from Delete(EdgeList)."
				call AdjacentPanelSearch(inPanel,xyz,nextPanel,aParticles,anEdges,aPanels)
				exit
			else
				edgeCounter = edgeCounter + 1
			endif
			current=>current%next
		enddo
		!print *,"DEBUG : edgeCounter = ",edgeCounter
		!print *,"DEBUG : panel center = ",aPanels%x(:,startPanel)
		if ( edgeCounter == numberOfEdges) then
			call Delete(leafEdgeList)
			inPanel = startPanel
			return			
		endif
	else ! Case 2 : panel has children; use tree search
		! Find the nearest child
		currentMin = PI
		currentPanel = 0
		do j=1,4
			testDist = ChordDistance(xyz,apanels%x(:,aPanels%children(j,startPanel)))
			if ( testDist < currentMin ) then
				currentMin = testDist
				currentPanel = aPanels%children(j,startPanel)
			endif
		enddo
		! Restart search with nearest child
		call AdjacentPanelSearch(inPanel,xyz,currentPanel,aParticles,anEdges,aPanels)
	endif
end subroutine


function FindRootPanel(xyz,apanels)
! Finds the closest panel center of the root polyhedron (cubed sphere or spherical icosahedron)
! to the point xyz.  Used to initialize the tree search subroutine.
	integer(kint) :: FindRootPanel
	real(kreal), intent(in) :: xyz(3)
	type(Panels), intent(in) :: aPanels
	integer(kint) :: panelKind, nRootFaces, j
	real(kreal) :: currentMin, testDist
	panelKind = GetPanelKind(aPanels)
	if ( panelKind == TRI_PANEL) then
		nRootFaces = 20
	elseif (panelKind == QUAD_PANEL) then
		nRootFaces = 6
	else
		print *,"FindRootPanel ERROR : invalid panelkind."
	endif
	currentMin = PI
	FindRootPanel = 1
	do j=1,nRootFaces
		testDist = ChordDistance(xyz,aPanels%x(:,j))
		if ( testDist < currentMin) then
			currentMin = testDist
			FindRootPanel = j
		endif
	enddo
end function


end module
