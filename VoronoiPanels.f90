module VoronoiPanelsModule

use NumberKindsModule
use SphereGeomModule
use ParticlesEdgesPanelsModule
use LoggerModule
use STRIPACKModule
!use RefineRemeshModule

implicit none

private
public VorPanels
public New, Delete
public DelaunayTri, VoronoiGrid
public VTKOutputDelaunayTriangulation, VTKOutputVoronoiGrid
public GetNTracer
public NormalizeGrid
public numThreads

!----------------
! Types and Module variables
!----------------

type VorPanels
	! Voronoi generators (active points)
	real(kreal), dimension(:), pointer :: x, y, z 		! Physical coordinate
	real(kreal), dimension(:), pointer :: x0, y0, z0	! Lagrangian coordinate
	real(kreal), dimension(:), pointer :: area			! Area of Voronoi polygons
	real(kreal), dimension(:), pointer :: absVort		! absolute vorticity 
	real(kreal), dimension(:), pointer :: relVort		! Relative vorticity
	real(kreal), dimension(:,:), pointer :: tracer		! passive tracers
	! Topological data
	integer(kint), dimension (:), pointer :: numVerts	! number of vertices in each voronoi polygon
	integer(kint), dimension(:,:), pointer :: vertIndices ! Indices to xc, yc, zc for vornoi corners
	integer(kint), dimension(:,:), pointer :: triangles ! Indices to other generators of Delaunay triangulation
	integer(kint) :: N									! number of generators in calculation
	integer(kint) :: N_Max								! Max number of generators allowed in memory
	integer(kint) :: N_Tri								! Number of triangles = number of voronoi corners
	
	! Voronoi vertices (passive points)
	real(kreal), dimension(:), pointer :: xc, yc, zc	! Coordinates of voronoi corners
	
	! STRIPACK Variables
	real(kreal), dimension(:), pointer :: dist, rc		! workspaces
	integer(kint), dimension(:), pointer :: list		! pointers for STRIPACK triangulation 
	integer(kint), dimension(:), pointer :: lptr		! pointers for STRIPACK triangulation 
	integer(kint), dimension(:), pointer :: lend		! pointers for STRIPACK triangulation 
	integer(kint), dimension(:), pointer :: near		! workspace
	integer(kint), dimension(:), pointer :: next		! workspace
	integer(kint), dimension(:), pointer :: listc		! pointer to xc,yc,zc for STRIPACK
	integer(kint), dimension(:,:), pointer :: ltri		! pointer to boundary nodes (not used here)
	integer(kint) :: lnew								! index of next free space in STRIPACK arrays
	integer(kint) :: stripackErrCode					! Error code
	integer(kint) :: nB									! workspace
	integer(kint) :: nCol								! workspace
end type

integer(kint), parameter :: MAX_POLY_SIDES = 20

integer(kint), save :: numThreads
integer(kint) :: logLevel = TRACE_LOGGING_LEVEL
type(Logger) :: Log
logical(klog), save :: logInit = .False.
integer(kint), save :: refCount = 0
character(len=28), save :: logKey = "VorPanels : "
!----------------
! Interfaces
!----------------
interface New
	module procedure NewPrivate
	module procedure NewCopyPrivate
	module procedure NewFromRemesh
end interface

interface Delete
	module procedure DeletePrivatePointer
end interface

interface Assignment(=)
	module procedure AssignPrivate
end interface

interface GetNTracer
	module procedure GetNTracerPrivate
end interface

contains
!----------------
! Standard Methods : Construction, destruction, copy, assignment
!----------------

subroutine NewCopyPrivate(self,other)
	type(VorPanels), intent(in) :: other
	type(VorPanels), pointer, intent(out) :: self
	integer(kint) :: nTracer, nMax
	
	nTracer = GetNTracer(other)
	nMax = other%N_Max
	
	allocate(self)
	call AllocMemory(self, nMax, nTracer)
	call CopyPrivate(self,other)
	refCount = refCount + 1
end subroutine


subroutine AssignPrivate(self,other)
	type(VorPanels), intent(in) :: other
	type(VorPanels), intent(out) :: self
	call CopyPrivate(self,other)
end subroutine


subroutine DeletePrivatePointer(self)
	type(VorPanels), pointer, intent(inout) :: self
	call DeletePrivate(self)
	deallocate(self)
	nullify(self)
end subroutine


subroutine CopyPrivate(self,other)
	type(VorPanels), intent(out) :: self
	type(VorPanels), intent(in) :: other
	self%x = other%x
	self%y = other%y
	self%z = other%z
	self%x0 = other%x0
	self%y0 = other%y0
	self%z0 = other%z0
	self%area = other%area
	self%relVort = other%relVort
	self%absVort = other%absVort
	if ( associated(other%tracer) ) then
		self%tracer = other%tracer
	else
		nullify(self%tracer)
	endif
	self%numVerts = other%numVerts
	self%vertIndices = other%vertIndices
	self%triangles = other%triangles
	self%N = other%N
	self%N_Max = other%N_Max
	self%N_Tri = other%N_Tri
	self%xc = other%xc
	self%yc = other%yc
	self%zc = other%zc
	self%rc = other%rc
	self%dist = other%dist
	self%list = other%list
	self%lptr = other%lptr
	self%lend = other%lend
	self%lnew = other%lnew
	self%near = other%near
	self%next = other%next
	self%listc = other%listc
	self%ltri = other%ltri
	self%stripackErrCode = other%stripackErrCode
	self%nB = other%nB
	self%nCol = other%nCol
end subroutine


subroutine AllocMemory(self, nGens, nTracer)
	type(VorPanels), intent(out) :: self
	integer(kint), intent(in) :: nGens, nTracer
	integer(kint) :: nTri, nArcs
	
	if ( nGens <= 0 ) then
		call LogMessage(log,ERROR_LOGGING_LEVEL,logKey,"AllocMemory ERROR : nGens <= 0.")
		return
	endif
	
	nTri = 2*nGens - 4
	nArcs = 3*nGens - 6
	
	! Reserve memory for n=nGens Voronoi generators (active particles)
	self%N = 0
	self%N_Max = nGens
	allocate(self%x(nGens))
	allocate(self%y(nGens))
	allocate(self%z(nGens))
	allocate(self%x0(nGens))
	allocate(self%y0(nGens))
	allocate(self%z0(nGens))
	allocate(self%area(nGens))
	allocate(self%absVort(nGens))
	allocate(self%relVort(nGens))
	if ( nTracer > 0 ) then
		allocate(self%tracer(nGens,nTracer))
		self%tracer = 0.0_kreal
	else
		nullify(self%tracer)
	endif
	self%x = 0.0_kreal
	self%y = 0.0_kreal
	self%z = 0.0_kreal
	self%x0 = 0.0_kreal
	self%y0 = 0.0_kreal
	self%z0 = 0.0_kreal
	self%area = 0.0_kreal
	self%absVort = 0.0_kreal
	self%relVort = 0.0_kreal
	allocate(self%numVerts(nGens))
	allocate(self%vertIndices(MAX_POLY_SIDES,nGens))
	self%numVerts = 0
	self%vertIndices = 0
	allocate(self%dist(nGens))
	allocate(self%near(nGens))
	self%dist = 0.0_kreal
	self%near = 0
	allocate(self%lend(nGens))
	allocate(self%next(nGens))
	allocate(self%ltri(1,1))
	self%lend = 0
	self%next = 0
	self%ltri = 0
	
	! Reserve memory for nTri = 2*n-4 Voronoi corners (passive particles)
	! This is equivalent to the number of triangles in the Delaunay triangulation.
	allocate(self%triangles(3,nTri))
	self%triangles = 0
	allocate(self%xc(nTri))
	allocate(self%yc(nTri))
	allocate(self%zc(nTri))
	allocate(self%rc(nTri))
	self%xc = 0.0_kreal
	self%yc = 0.0_kreal
	self%zc = 0.0_kreal
	self%rc = 0.0_kreal
	allocate(self%list(2*nArcs))
	allocate(self%lptr(2*nArcs))
	allocate(self%listc(2*nArcs))
	self%list = 0
	self%lptr = 0
	self%listc = 0
	self%stripackErrCode = 0
	
	self%nCol = 0
	self%nB = 0
end subroutine


subroutine NewPrivate(self, initNest,AMR,nTracer,amrTol1, amrTol2, problemID,&
	 lat0, beta, perturbAmp, perturbWaveNum)
	type(VorPanels), pointer, intent(out) :: self
	integer(kint), intent(in) :: initNest,AMR, nTracer
	real(kreal), intent(in), optional :: amrTol1, amrTol2, lat0, beta, perturbAmp
	integer(kint), intent(in), optional :: problemID, perturbWaveNum
	type(Particles), pointer :: tempParticles
	type(Edges), pointer :: tempEdges
	type(Panels), pointer :: tempPanels
	integer(kint), parameter :: panelKind = TRI_PANEL, problemKind=BVE_SOLVER
	integer(kint) :: i, j, k, MemLimit, nGens, nActive, nPassive
	integer(kint) :: nCol, nB
	logical(klog) :: useAMR
	
	
	allocate(self)
	
	! Start logger for this module
	if ( .NOT. logInit) then
		call InitLogger(Log)
		refCount = 1
	else
		refCount = refCount + 1
	endif
	call LogMessage(Log,DEBUG_LOGGING_LEVEL,logKey, "Started Logger.")
	
	!numThreads = OMP_GET_NUM_THREADS()
	
	useAMR = .False.
	! Use icosahedron to generate initial grid
	call New(tempParticles,tempEdges,tempPanels,panelKind,initNest,AMR,nTracer,problemKind)
!!!!!!!!!!
!   AMR		: Use old grids to do amr
!!!!!!!!!!	
!	if ( AMR > 0 ) then
!		if ( (present(amrTol1) .AND. present(amrTol2)) .And. present(problemID) ) useAMR = .TRUE.
!	endif
!	if ( useAMR) then
!		if ( problemID /= JET ) then
!			!call InitRefine(tempParticles,tempEdges,tempPanels,amrTol1,amrTol2,problemID)
!		else
!			!call InitRefine(tempParticles,tempEdges,tempPanels,amrTol1,amrTol2,problemID,&
!			!	lat0,beta,perturbAmp,perturbWaveNum)
!		endif
!	endif
!	
	
	
	call LogMessage(log,DEBUG_LOGGING_LEVEL,logKey," old grid returned.")
	nGens = tempParticles%N
	if ( AMR > 0 ) then
		memLimit = tempParticles%N_Max
	else
		memlimit = nGens
	endif
	call AllocMemory(self,memLimit,nTracer)
	do i = 1,tempParticles%N
		self%x(i) = tempParticles%x(1,i)
		self%y(i) = tempParticles%x(2,i)
		self%z(i) = tempParticles%x(3,i)
	enddo
	self%x0(1:tempParticles%N) = self%x(1:tempParticles%N)
	self%y0(1:tempParticles%N) = self%y(1:tempParticles%N)
	self%z0(1:tempParticles%N) = self%z(1:tempParticles%N)
	self%N = tempParticles%N
	
	! Get Delaunay triangulation from STRIPACK
	call DelaunayTri(self)
	! Make Voronoi polygons with STRIPACK
	call VoronoiGrid(self)
end subroutine


subroutine NewFromRemesh(self,nPanels,nTracer,triParticles)
	type(VorPanels), intent(out) :: self
	type(Particles), intent(in) :: triParticles
	integer(kint), intent(in) :: nPanels, nTracer
	integer(kint) :: j
	call AllocMemory(self,nPanels,nTracer)
	do j=1,triParticles%N
		self%x(j) = triParticles%x(1,j)
		self%x0(j) = triParticles%x0(1,j)
		self%y(j) = triParticles%x(2,j)
		self%y0(j) = triParticles%x0(2,j)
		self%z(j) = triParticles%x(3,j)
		self%z0(j) = triParticles%x0(3,j)
		self%absVort(j) = triParticles%absVort(j)
		self%relVort(j) = triParticles%relVort(j)
	enddo
	self%N = triParticles%N
	call DelaunayTri(self)
	call VoronoiGrid(self)
end subroutine


subroutine DeletePrivate(self)
	type(VorPanels), intent(inout) :: self
	deallocate(self%x)
	deallocate(self%y)
	deallocate(self%z)
	deallocate(self%x0)
	deallocate(self%y0)
	deallocate(self%z0)
	deallocate(self%absVort)
	deallocate(self%relVort)
	deallocate(self%area)
	if ( associated(self%tracer) ) deallocate(self%tracer)
	deallocate(self%numVerts)
	deallocate(self%vertIndices)
	deallocate(self%triangles)
	deallocate(self%xc)
	deallocate(self%yc)
	deallocate(self%zc)
	deallocate(self%dist)
	deallocate(self%rc)
	deallocate(self%list)
	deallocate(self%lptr)
	deallocate(self%lend)
	deallocate(self%near)
	deallocate(self%next)
	deallocate(self%listc)
	deallocate(self%ltri)

	call LogMessage(Log,TRACE_LOGGING_LEVEL,logKey,"Deleting : VorPanels and Logger.")
	if ( refCount == 1 ) then
		refCount = 0
		call Delete(Log)
		logInit = .False.
	else 
		refCount = refCount -1
	endif
end subroutine


function GetNTracerPrivate(self)
	type(VorPanels), intent(in) :: self
	integer(kint) :: GetNTracerPrivate
	if ( associated(self%tracer) ) then
		GetNTracerPrivate = size(self%tracer,2)
	else
		GetNTracerPrivate = 0
	endif
end function


!----------------
! STRIPACK Interfaces : Transform Delaunay triangulation and voronoi grids into VorPanels
!						data type.
!----------------

subroutine NormalizeGrid(self)
	type(VorPanels), intent(inout) :: self
	integer(kint) :: j
	real(kreal) :: norm, newx, newy, newz
	do j=1,self%N
		norm = sqrt( self%x(j)*self%x(j) + self%y(j)*self%y(j) + self%z(j)*self%z(j) )
		self%x(j) = self%x(j)/norm
		self%y(j) = self%y(j)/norm
		self%z(j) = self%z(j)/norm
	enddo
end subroutine


subroutine DelaunayTri(self)
	type(VorPanels), intent(inout) :: self
	self%stripackErrCode = 0
	call TRMESH(self%N,self%x,self%y,self%z,self%list,self%lptr,self%lend,&
		self%lnew,self%near,self%next,self%dist,self%stripackErrCode)
	if ( self%stripackErrCode > 0 ) then
		call LogMessage(Log,ERROR_LOGGING_LEVEL,trim(logKey)//" TRMESH a duplicate node at node ",self%stripackErrCode)
		return
	endif
	self%stripackErrCode = 0
	call TRLIST2(self%N,self%list,self%lptr,self%lend,self%N_Tri,self%triangles,self%stripackErrCode)
	if ( self%stripackErrCode > 0 ) then
		call LogMessage(log,ERROR_LOGGING_LEVEL,trim(logKey)//" TRLIST2 errorCode = ",self%stripackErrCode)
		return
	endif
end subroutine


subroutine VoronoiGrid(self)
	type(VorPanels), intent(inout) :: self
	integer(kint) :: faceCount, edgeCount, vertCount
	integer(kint) :: i, j, k
	integer(kint) :: lp, lpl, kt, nv
	logical(klog) :: keepGoing
	self%stripackErrCode = 0
	call CRList(self%N,self%nCol,self%x,self%y,self%z,self%list,self%lend,self%lptr,&
		self%lnew,self%ltri,self%listc,self%nb,self%xc,self%yc,self%zc,self%rc,self%stripackErrCode)
	if (self%stripackErrCode == 1) then
		call LogMessage(Log,ERROR_LOGGING_LEVEL,trim(logKey)//" CRLIST ERROR : ","Less than 3 input nodes.")
	elseif (self%stripackErrCode == 2 ) then
		call LogMessage(Log,ERROR_LOGGING_LEVEL,trim(logKey)//" CRLIST ERROR : ","nCol < nB -2")
	elseif (self%stripackErrCode == 3 ) then
		call LogMessage(Log,ERROR_LOGGING_LEVEL,trim(logKey)//" CRLIST ERROR : ","found degenerate triangle.")
	endif
	faceCount = self%N
	vertCount = 0
	edgeCount = 0
	do i=1,self%N	
		lpl = self%lend(i) 	! initialize lpl to last index of triangle i
		lp = lpl		  	! initialize lp to first index of triangle i
		nv = 0				! count the number of vertices
		keepGoing = .True.	! loop termination condition
		do while (keepGoing)
			lp = self%lptr(lp)	! move to next node
			kt = self%listc(lp) ! kt points to indices in xc,yc,zc stripack arrays of voronoi corners
			nv = nv + 1			! add another vertex to polygon i
			self%vertIndices(nv,i) = kt ! keep track of kt in VorPanels data structure
			vertCount = max(vertCount,kt) ! This prevents vertices from being counted twice
			edgeCount = edgeCount + 1 ! edges are counted twice
			if ( nV > MAX_POLY_SIDES ) then ! check for invalid state
				call LogMessage(log,ERROR_LOGGING_LEVEL,trim(logKey)//"VorGrid ERROR: exceeded max poly size at node ",i)
				return
			endif
			! Stop when polygon has closed.
			if ( lp == lpl ) keepGoing = .False.
		enddo
		self%numVerts(i) = nv
	enddo
	edgeCount = edgeCount/2
	
	if ( logLevel == DEBUG_LOGGING_LEVEL) then
		! Check Euler's formula
		call LogMessage(Log,DEBUG_LOGGING_LEVEL,trim(logKey)//" n Faces = ",faceCount)
		call LogMessage(Log,DEBUG_LOGGING_LEVEL,trim(logKey)//" n Vertices =",vertCount)
		call LogMessage(Log,DEBUG_LOGGING_LEVEL,trim(logKey)//" n Edges = ",edgeCount)
		call LogMessage(Log,DEBUG_LOGGING_LEVEL,trim(logKey)//" F + V - E -2 = ",faceCount + vertCount - edgeCount - 2)
	endif
	
	! Set areas on polygons
	do i=1,self%N
		self%area(i) = 0.0_kreal
		do j=1,self%numVerts(i) ! loop over subtriangles
			self%area(i) = self%area(i) + SphereTriArea( & 
			self%xc(self%vertIndices(j,i)),self%yc(self%vertIndices(j,i)),self%zc(self%vertIndices(j,i)), & ! first vertex of polygon subtriangle
			self%x(i), self%y(i), self%z(i), & ! second vertex is polygon's generator (the "center" of the polygon)
			self%xc(self%vertIndices(mod(j,self%numVerts(i))+1,i)),& ! third vertex is next one in vertIndices list
			self%yc(self%vertIndices(mod(j,self%numVerts(i))+1,i)),&
			self%zc(self%vertIndices(mod(j,self%numVerts(i))+1,i)) )
		enddo
	enddo
	call LogMessage(log,DEBUG_LOGGING_LEVEL,trim(logKey)//" Voronoi surf. area = ",sum(self%area(1:self%N)))
end subroutine


!----------------
! Output Methods : Console and VTK
!----------------


subroutine VTKOutputDelaunayTriangulation(self,filename)
	type(VorPanels), intent(in) :: self
	character(len=*), intent(in) :: filename
	
	integer(kint) :: nPoints, nScalars, cellListSize, errCode
	integer(kint), parameter :: unitNumber = 11
	integer(kint) :: i,j,k
	character(len=28) :: dataString
	
	open(unit=unitNumber,file=filename,status='REPLACE',action='WRITE',iostat=errCode)
		if ( errCode /= 0 ) then
			call LogMessage(Log,ERROR_LOGGING_LEVEL,"VTKOutputDelaunay ERROR : ","cannot open file.")
			return
		endif
	! WRITE VTK HEADER
	write(unitNumber,'(A)') '# vtk DataFile Version 2.0'
	write(unitNumber,'(A,I8,A)') 'Voronoi Diagram : ',self%N,' generators'
	write(unitNumber,'(A)') 'ASCII'
	write(unitNumber,'(A)') 'DATASET POLYDATA'	
	
	! WRITE VTK POINTS DATA (voronoi generators)
	write(unitNumber,'(A,I8,A)') 'POINTS ',self%N,'  double'
	do j=1,self%N
		write(unitNumber,'(3F24.18)') self%x(j), self%y(j), self%z(j)
	enddo
	
	! WRITE VTK POLYGON DATA
	! N generators give 2*N-4 triangles. Each triangle requires 4 integers to record.
	cellListSize = 4*(2*self%N-4)
	write(unitNumber,'(A,I8,I8)') 'POLYGONS ',2*self%N-4,cellListSize
	do j=1,2*self%N-4
		write(unitNumber,'(4I8)') 3, self%triangles(1,j)-1, self%triangles(2,j)-1, self%triangles(3,j)-1
	enddo
	
	! WRITE SCALARS DATA
	if ( associated(self%tracer) ) then
		nScalars = size(self%tracer,2)
	else
		nScalars = 0
	endif
	write(unitNumber,'(A,I8)') 'POINT_DATA ',self%N
	write(unitNumber,'(A)') 'SCALARS relVort double 1'
	write(unitNumber,'(A)') 'LOOKUP_TABLE default'
	do j=1,self%N
		write(unitNumber,'(F24.15)') self%relVort(j)
	enddo
	write(unitNumber,'(A)') 'SCALARS absVort double 1'
	write(unitNumber,'(A)') 'LOOKUP_TABLE default'
	do j=1,self%N
		write(unitNumber,'(F24.15)') self%absVort(j)
	enddo
	do k=1,nScalars
		write(dataString,'(A,I1)') trim('tracer'),k
		write(unitNumber,'(A,A,A)') 'SCALARS ',trim(dataString),' double 1'
		write(unitNumber,'(A)') 'LOOKUP_TABLE default'
		do j=1,self%N
			write(unitNumber,'(F24.15)') self%tracer(j,k)
		enddo
	enddo
	
	close(unitNumber)	
end subroutine





subroutine VTKOutputVoronoiGrid(self,filename)
	type(VorPanels), intent(in) :: self
	character(len=*), intent(in) :: filename
	!local variables
	integer(kint) :: nPoints, nScalars, cellListSize, errCode
	integer(kint), parameter :: unitNumber = 11
	integer(kint) :: j, k
	character(len=28) :: dataString
	
	open(unit=unitNumber,file=filename,status='REPLACE',action='WRITE',iostat=errCode)
		if ( errCode /= 0 ) then
			call LogMessage(Log,ERROR_LOGGING_LEVEL,"VTKOutputVoronoi ERROR : ","cannot open file.")
			return
		endif
	! WRITE VTK HEADER
	write(unitNumber,'(A)') '# vtk DataFile Version 2.0'
	write(unitNumber,'(A,I8,A)') 'Voronoi Diagram : ',self%N,' generators'
	write(unitNumber,'(A)') 'ASCII'
	write(unitNumber,'(A)') 'DATASET POLYDATA'
	
	
	! WRITE VTK POINTS DATA (vertices of voronoi polygons)
	write(unitNumber,'(A,I8,A)') 'POINTS ',2*self%N-4,'  double'
	do j=1,2*self%N-4
		write(unitNumber,'(3F24.18)') self%xc(j), self%yc(j), self%zc(j)
	enddo
	
	! WRITE VTK POLYGON DATA
	cellListSize = self%N + sum(self%numVerts(1:self%N))
	write(unitNumber,'(A,I8,I8)') 'POLYGONS ',self%N,cellListSize
	do j=1,self%N
		write(unitNumber,'(I8)',advance='NO') self%numVerts(j)
		do k=1,self%numVerts(j)
			write(unitNumber,'(I8)',advance='NO') self%vertIndices(k,j)-1
		enddo
		write(unitNumber,'(A)',advance='YES') ' '
	enddo
	
	! WRITE SCALARS DATA
	if ( associated(self%tracer) ) then
		nScalars = size(self%tracer,2)
	else
		nScalars = 0
	endif
	write(unitNumber,'(A,I8)') 'CELL_DATA ',self%N
	write(unitNumber,'(A)') 'SCALARS relVortPanel double 1'
	write(unitNumber,'(A)') 'LOOKUP_TABLE default'
	do j=1,self%N
		write(unitNumber,'(F24.15)') self%relVort(j)
	enddo
	write(unitNumber,'(A)') 'SCALARS absVortPanel double 1'
	write(unitNumber,'(A)') 'LOOKUP_TABLE default'
	do j=1,self%N
		write(unitNumber,'(F24.15)') self%absVort(j)
	enddo
	do k=1,nScalars
		write(dataString,'(A,I1)') trim('tracerPanel'),k
		write(unitNumber,'(A,A,A)') 'SCALARS ',trim(dataString),' double 1'
		write(unitNumber,'(A)') 'LOOKUP_TABLE default'
		do j=1,self%N
			write(unitNumber,'(F24.15)') self%tracer(j,k)
		enddo
	enddo
	
	close(unitNumber)
end subroutine

!----------------
! Logging
!----------------
subroutine InitLogger(aLog)
	type(Logger), intent(inout) :: aLog
	integer(kint) :: logUnit
	logUnit = 6
	call New(aLog,logLevel,logUnit)
	logInit = .TRUE.
end subroutine

!subroutine ExtractDelaunayTriangulation(self)
!	type(VorPanels), intent(inout) :: self
!	integer(kint) :: i, nt
!	integer(kint) :: nrow
!	nrow = size(self%ltri,1)
!	
!	call LogMessage(Log,TRACE_LOGGING_LEVEL,logKey," Entering : ExtractDelaunayTriangulation.")
!		
!	call TRList2(self%N,self%list,self%lptr,self%lend,nt,self%triangles,self%stripackErrCode)
!	call LogMessage(Log,DEBUG_LOGGING_LEVEL,trim(logKey)//" TRLIST2 returned errCode = ",self%stripackErrCode)
!	
!end subroutine

!subroutine GenerateVoronoiGrid(self)
!	type(VorPanels), intent(inout) :: self
!	integer(kint) :: i, j, k
!	integer(kint) :: lp, lpl, kt, nv
!	logical(klog) :: keepGoing
!	integer(kint) :: vertCount, edgeCount, faceCount
!	
!	call LogMessage(Log,TRACE_LOGGING_LEVEL," Entering : "," GenerateVoronoiGrid")
!	
!	faceCount = self%N
!	vertCount = 0
!	edgeCount = 0
!
!	do i = 1,self%N
!		lpl = self%lend(i)
!		lp = lpl
!		nv = 0
!		keepGoing = .True.
!		do while (keepGoing)
!			lp = self%lptr(lp)
!			kt = self%listc(lp)
!			nv = nv + 1
!			self%vertIndices(nv,i) = kt
!			vertCount = max(vertCount,kt)
!			edgeCount = edgeCount + 1
!			if ( nv > MAX_POLY_SIDES ) then
!				call LogMessage(Log,ERROR_LOGGING_LEVEL,'Max Sides exceeded at node ',i)
!				return
!			endif
!			if ( lp == lpl ) keepGoing = .False.
!		enddo
!		self%numVerts(i) = nv
!		
!	enddo
!	edgeCount = edgeCount/2
!	
!	if ( logLevel == DEBUG_LOGGING_LEVEL) then
!		! Check Euler's formula
!		call LogMessage(Log,DEBUG_LOGGING_LEVEL,trim(logKey)//" n Faces = ",faceCount)
!		call LogMessage(Log,DEBUG_LOGGING_LEVEL,trim(logKey)//" n Vertices =",vertCount)
!		call LogMessage(Log,DEBUG_LOGGING_LEVEL,trim(logKey)//" n Edges = ",edgeCount)
!		call LogMessage(Log,DEBUG_LOGGING_LEVEL,trim(logKey)//" F + V - E -2 = ",faceCount + vertCount - edgeCount - 2)
!	endif
!	
!	! Set areas from Voronoi polygons
!	do i=1,self%N
!		self%area(i) = 0.0_kreal
!		do j=1,self%numVerts(i)
!			self%area(i) = self%area(i) + SphereTriArea( &
!				self%xc(self%vertIndices(j,i)), self%yc(self%vertIndices(j,i)), self%zc(self%vertIndices(j,i)), &
!				self%x(i),self%y(i),self%z(i), &
!				self%xc(self%vertIndices(mod(j,self%numVerts(i))+1,i)),&
!				self%yc(self%vertIndices(mod(j,self%numVerts(i))+1,i)),&
!				self%zc(self%vertIndices(mod(j,self%numVerts(i))+1,i)))
!		enddo
!	enddo	
!	call LogMessage(Log,DEBUG_LOGGING_LEVEL,trim(logKey)//"surf area = ",sum(self%area(1:self%N)))
!	
!end subroutine

!subroutine NewFromPanelsPrivate( self, initNest, AMR, nTracer )
!	type(VorPanels), intent(out) :: self
!	integer(kint), intent(in) :: initNest, AMR, nTracer
!	! Local variables
!	type(Particles), pointer :: tempParticles
!	type(Edges), pointer :: tempEdges
!	type(Panels), pointer :: tempPanels, activePanels, passivePanels
!	integer(kint), parameter :: panelKind = TRI_PANEL, problemKind=BVE_SOLVER
!	integer(kint), allocatable :: activeMap(:), passiveMap(:)
!	integer(kint) :: i, j, k, MemLimit, nGens, nActive, nPassive
!	integer(kint) :: nCol, nB
!	
!	if ( .NOT. logInit) call InitLogger(Log)
!	
!!	if ( associated(self)) then
!!		call LogMessage(Log,ERROR_LOGGING_LEVEL,logkey,&
!!			"New VorPanels ERROR : pointer already associated.")
!!		return
!!	endif
!	
!	call LogMessage(Log,TRACE_LOGGING_LEVEL,logKey, "Started Logger.")
!
!	! Use icosahedron to generate initial grid
!
!	call New(tempParticles,tempEdges,tempPanels,panelKind,initNest,AMR,nTracer,problemKind)
!	nActive = tempPanels%N_Active
!	nPassive = tempPanels%N-tempPanels%N_Active
!	call LogMessage(Log,TRACE_LOGGING_LEVEL,logKey,"old Grid returned.")
!	if ( nPassive > 0 ) then
!		allocate(activePanels)
!		call New(activePanels,nActive,panelKind,nTracer,problemKind)
!		activePanels%N = tempPanels%N_Active
!		activePanels%N_Active = activePanels%N
!		allocate(activeMap(nActive))
!
!		allocate(passivePanels)
!		call New(passivePanels,nPassive,panelKind,nTracer,problemKind)	
!		allocate(passiveMap(nPassive))
!		passivePanels%N = tempPanels%N-tempPanels%N_Active
!		call GatherPanels(tempPanels,activePanels,activeMap,passivePanels,passiveMap)
!	else
!		activePanels=>tempPanels
!	endif
!	
!
!	call LogMessage(log,TRACE_LOGGING_LEVEL,logKey,"old grid gathered.")
!
!	
!	nGens = tempPanels%N_Active
!	if ( AMR > 0 ) then
!		memLimit = tempPanels%N_Max
!	else
!		memLimit = tempPanels%N_Active
!	endif
!	
!	! Generator data	
!!	allocate(self)
!	allocate(self%x(memLimit))
!	allocate(self%y(memLimit))
!	allocate(self%z(memLimit))
!	self%x = 0.0_kreal
!	self%y = 0.0_kreal
!	self%z = 0.0_kreal
!	allocate(self%x0(memLimit))
!	allocate(self%y0(memLimit))
!	allocate(self%z0(memLimit))
!	self%x0 = 0.0_kreal
!	self%y0 = 0.0_kreal
!	self%z0 = 0.0_kreal
!	allocate(self%numVerts(memLimit))
!	self%numVerts = 0
!	allocate(self%absVort(memLimit))
!	self%absVort = 0.0_kreal
!	allocate(self%relVort(memLimit))
!	self%relVort = 0.0_Kreal
!	allocate(self%area(memLimit))
!	self%area = 0.0_kreal
!	allocate(self%triangles(3,2*memLimit-4))
!	self%triangles = 0
!	if (nTracer > 0 ) then
!		allocate(self%tracer(memLimit,nTracer))
!		self%tracer = 0.0_kreal
!	else
!		nullify(self%tracer)
!	endif
!	
!	! Voronoi corner data
!	allocate(self%vertIndices(MAX_POLY_SIDES,2*memLimit-4))
!	self%vertIndices = 0
!	allocate(self%xc(2*memLimit-4))
!	allocate(self%yc(2*memLimit-4))
!	allocate(self%zc(2*memLimit-4))
!	self%xc = 0.0_kreal
!	self%yc = 0.0_kreal
!	self%zc = 0.0_kreal
!	
!	
!	! Set up STRIPACK workspaces and arrays
!	allocate(self%list(6*memLimit-12))
!	allocate(self%lptr(6*memLimit-12))
!	allocate(self%lend(memLimit))
!	allocate(self%near(memLimit))
!	allocate(self%next(memLimit))
!	allocate(self%dist(memLimit))
!	allocate(self%listc(6*memLimit-12))
!	allocate(self%rc(2*memLimit-4))
!	allocate(self%ltri(1,1))
!	self%list = 0
!	self%lptr = 0
!	self%lend = 0
!	self%near = 0
!	self%next = 0
!	self%dist = 0.0_kreal
!	self%listc = 0
!	self%rc = 0.0_kreal
!	self%ltri = 0
!	self%stripackErrCode = 0
!	
!	do i=1,tempPanels%N_Active
!		self%x(i) = activePanels%x(1,i)
!		self%y(i) = activePanels%x(2,i)
!		self%z(i) = activePanels%x(3,i)
!		self%x0(i) = activePanels%x0(1,i)
!		self%y0(i) = activePanels%x0(2,i)
!		self%z0(i) = activePanels%x0(3,i)
!	enddo
!	self%N = tempPanels%N_Active
!	self%N_Max = memLimit
!	
!	! Generate the Delaunay triangulation with STRIPACK
!	call TRMesh(self%N,self%x,self%y,self%z,self%list,self%lptr,self%lend, &
!				self%lnew,self%near,self%next,self%dist,self%stripackErrCode)
!	call LogMessage(Log,DEBUG_LOGGING_LEVEL,trim(logKey)//" TRMESH returned errCode = ",self%stripackErrCode)
!	if ( self%stripackErrCode > 0 ) then
!		call LogMessage(Log,ERROR_LOGGING_LEVEL,trim(logKey)//" TRMESH a duplicate node at node ",self%stripackErrCode)
!		return
!	endif
!	! Generate the Voronoi vertices with STRIPACK
!	nCol = 0
!	nB = 0	
!	call CRList(self%N,nCol,self%x,self%y,self%z,self%list,self%lend,self%lptr,&
!			self%lnew,self%ltri,self%listc,nb,self%xc,self%yc,self%zc,self%rc,self%stripackErrCode)
!	call LogMessage(Log,DEBUG_LOGGING_LEVEL,trim(logKey)//"CRList returned errCode = ",self%stripackErrCode)
!	if (self%stripackErrCode == 1) then
!		call LogMessage(Log,ERROR_LOGGING_LEVEL,trim(logKey)//" CRLIST ERROR : ","Less than 3 input nodes.")
!	elseif (self%stripackErrCode == 2 ) then
!		call LogMessage(Log,ERROR_LOGGING_LEVEL,trim(logKey)//" CRLIST ERROR : ","nCol < nB -2")
!	elseif (self%stripackErrCode == 3 ) then
!		call LogMessage(Log,ERROR_LOGGING_LEVEL,trim(logKey)//" CRLIST ERROR : ","found degenerate triangle.")
!	endif
!	
!	! Generate the area weights for each generator from the voronoi diagram
!	call GenerateVoronoiGrid(self)
!	call ExtractDelaunayTriangulation(self)
!	
!	! Clean up
!	call Delete(tempParticles,tempEdges,tempPanels)
!	if ( associated(passivePanels) ) then
!		call Delete(activePanels)
!		call Delete(passivePanels)
!		deallocate(activeMap)
!		deallocate(passiveMap)
!		deallocate(activePanels)
!		deallocate(passivePanels)
!	else
!		nullify(activePanels)	
!	endif
!	
!end subroutine

end module
