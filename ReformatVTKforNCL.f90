program ReformatVTKtoNCL

use NumberKindsModule
use SphereGeomModule

implicit none

! File Processing
character(len=128) :: namelistFile
character(len=256) :: vtkDataFile,  dataFileHeader, outputPrefix
character(len=256) :: scalarOutputFile
character(len=128) :: lineIn, scratch1, scratch2
integer(kint), parameter :: readUnit = 11
integer(kint), parameter :: llWriteUnit = 12
real(kreal), allocatable :: vertexXYZ(:,:), pointsScalars(:,:), cellScalars(:,:), cellXYZ(:,:)
integer(kint), allocatable :: cellData(:,:)
integer(kint) :: nPoints, lineNumber,  scalarCount
integer(kint) :: readStat, llWriteStat
integer(kint) :: nVTKPoints, nVTKCells, polygonHeaderLine, pointDataLine, cellDataLine
character(len=56) :: scalarName, foundName, formatString
integer(kint) :: i, j, k, i1, i2, i3
logical(klog) :: keepGoing
real(kreal) :: xA(3), xB(3), xc(3), cellX(3)

!
real(kreal), allocatable :: lat(:), lon(:), llgrid(:,:)
real(kreal) :: degreeSpacing, radSpacing, leftLonDeg, leftLonRad
integer(kint) :: nLat, nLon, totalLLpoints, llCounter

!
!	stripack variables
!
real(kreal), allocatable :: x(:), y(:), z(:), dist(:)
integer(kint), allocatable :: list(:), lptr(:), lend(:), near(:), next(:)
integer(kint) :: n, lnew, errCode, startTriangle
!
!	ssrfpack variables
!
real(kreal), allocatable :: sourceData(:), sourceGrad(:,:), sigma(:)
real(kreal) :: dsig, sigTol
integer(kint) :: SIGMA_FLAG = 1, GRAD_FLAG = 1


namelist /filenames/ vtkDataFile, outputPrefix
namelist /interp/ degreeSpacing, leftLonDeg

print *, "This program takes particles and panels data from a .vtk file and "
print *, "uses nearest SSRFPACK to assign values to a uniform "
print *, "lat-long grid to use for contour plot purposes. "
print *, "User input via 'conversion.namelist' file."


namelistFile = 'smoothConversion.namelist'
! Get user input for filenames
open(unit=readUnit,file=namelistFile,status='OLD',delim='APOSTROPHE')
	read(readUnit,nml=filenames)
	read(readUnit,nml=interp)
close(readUnit)

!
! Read VTK Metadata
!
open(unit=readUnit,file=vtkDataFile,status='OLD',action='READ',iostat=readStat)
if ( readStat /= 0 ) stop "ERROR opening vtk input file."
! Preprocessing : Find number of points, number of cells, and number of scalars in data file
keepGoing = .TRUE.
lineNumber = 0
scalarCount = 0
polygonHeaderLine = -1
pointDataLine = -1
print '("Getting VTK File MetaData...")'
print '("VTK FILE INFO :")'
do while (keepGoing)
	! Read entire line from input file into string 'lineIn'
	read(readUnit,'(A)',iostat=readStat) lineIn
	lineNumber = lineNumber + 1
	if ( readStat > 0 ) then
		print '("ERROR : problem encountered at line ",I8," .")' , lineNumber
		keepGoing = .FALSE.
		stop
	elseif ( readStat < 0 ) then
		print '("End of file found at line ",I8,".")', lineNumber
		keepGoing = .FALSE.
	else
		if ( lineNumber == 5 ) then
			! Get number of vertices
			lineIn = adjustl(lineIn)
			lineIn = trim(lineIn)
			i1 = scan(lineIn,' ')
			i2 = scan(lineIn,' ',back=.True.)
			scratch1 = lineIn(i1+1:i2)
			read(scratch1,'(I8)') nVTKPoints
			print *,' nPoints = ',nVTKPoints
			polygonHeaderLine = 5 + nVTKPoints + 1
		elseif ( lineNumber == polygonHeaderLine ) then
			! Get number of panels
			lineIn = adjustl(lineIn)
			lineIn = trim(lineIn)
			i1 = scan(lineIn,' ')
			i2 = scan(trim(lineIn),' ',back=.TRUE.)
			scratch1 = lineIn(i1:i2)
			scratch1 = adjustl(scratch1)
			scratch1 = trim(scratch1)
			read(scratch1,'(I8)') nVTKCells
			print *,' nVTKCells = ',nVTKCells
			pointDataLine = polygonHeaderLine + nVTKCells + 1
		elseif ( lineNumber > pointDataLine .AND. lineNumber > 5) then
			!Count the number of scalars
			i1 = 0
			i2 = 0
			lineIn = adjustl(lineIn)
			i1 = scan(lineIn,'C')
			i2 = scan(lineIn,'S')
			if (i1 == 1) then ! Stop at 'CELL_DATA' line
				keepGoing = .FALSE.
				cellDataLine = lineNumber
			else
				if ( i2 == 1) scalarCount = scalarCount + 1
			endif
		endif
	endif
enddo

print *, " File contains ", scalarCount, " scalars."
print *, " polygonDataLine = ", polygonHeaderLine
print *, " pointDataLine = ",pointDataLine
print *, " cellDataLine = ", cellDataLine

!
!	allocate memory
!
allocate(vertexXYZ(3,nVTKPoints))
vertexXYZ = 0.0_kreal
allocate(pointsScalars(nVTKPoints,scalarCount))
pointsScalars = 0.0_kreal
allocate(cellData(4,nVTKCells))
cellData = -1
allocate(cellScalars(nVTKCells,scalarCount))
cellScalars= 0.0_kreal
allocate(cellXYZ(3,nVTKCElls))
cellXYZ = 0.0_kreal

!
!	read vtk data
!
print '("Reading VTK File data ...")'

write(datafileheader,'(2A24)') 'Longitude','Latitude'

rewind(readUnit)
keepGoing = .TRUE.
lineNumber = 0
j = 1
k = 1
do while (keepGoing)
	read(readUnit,'(A)',iostat=readStat) lineIn
	lineNumber = lineNumber + 1
	if ( readStat > 0 ) then
		print '("ERROR : problem encountered at line ",I8," .")' , lineNumber
		keepGoing = .FALSE.
		stop
	elseif ( readStat < 0 ) then
		print '("End of file found at line ",I8,".")', lineNumber
		keepGoing = .FALSE.
	else
		if ( lineNumber > 5 .AND. lineNumber < polygonHeaderLine) then
		! READ VERTEX COORDINATE DATA
			read(lineIn,'(3D)') vertexXYZ(:,j)
			j = j+1
		endif
		if ( lineNumber == polygonHeaderLine) j = 1 ! reset counter
		if ( lineNumber > polygonHeaderLine .AND. lineNumber < pointDataLine) then
		! READ CONNECTIVITY DATA
			lineIn = adjustl(lineIn)
			i1 = scan(lineIn,' ')
			scratch1 = lineIn(i1+1:len_trim(lineIn))
			read(scratch1,'(4I8)') cellData(:,k)
			k = k+1
		endif
		if ( lineNumber == pointDataLine ) k = 0 ! reset counter
		if ( lineNumber > pointDataLine .AND. lineNumber < cellDataLine) then
		! READ VERTEX SCALAR DATA
			lineIn = adjustl(lineIn)
			i2 = scan(lineIn,'S')
			if ( i2 == 1) then
				k = k+1
				j = 1
				i1 = scan(lineIn,' ')
				scratch1 = lineIn(i1+1:len_trim(lineIN))
				scratch1 = adjustl(scratch1)
				i3 = scan(scratch1,' ')
				foundName = scratch1(1:i3+1)
				write(datafileHeader,'(A,A24)') trim(dataFileHeader),trim(foundName)
			endif
			i1 = scan(lineIn,'L')
			if ( i1 == 0 ) then
				read(lineIn,'(D)') pointsScalars(j,k)
				j = j+1
			endif
		endif
		if ( lineNumber == cellDataLine ) k = 0
		if ( lineNumber > cellDataLine ) then
			lineIn = adjustl(lineIN)
			i2 = scan(lineIn,'S')
			if ( i2 == 1) then
				j = 1
				k = k+1
			endif
			i1 = scan(lineIn,'L')
			if ( i1 == 0) then
				read(lineIn,'(D)') cellScalars(j,k)
				j = j+1
			endif
		endif
	endif
enddo

! normalize
do j=1,nVTKPoints
	vertexXYZ(:,j) = vertexXYZ(:,j)/sqrt(sum(vertexXYZ(:,j)*vertexXYZ(:,j)))
enddo

print '("... VTK data loaded. Preparing to interpolate.")'
close(readUnit)

nlon = floor(360.0/degreeSpacing)
nlat = nlon/2 + 1
allocate(lat(nLat))
allocate(lon(nLon))
allocate(llGrid(nLat,nLon))

do j=1,nLon
	lon(j) = (j-1)*(2.0_kreal*PI/nLon)
enddo
do j=1,nLat
	lat(j) = -PI/2.0_kreal + (j-1)*(2.0_kreal*PI/nLon)
enddo
print '("... lat long vectors ready.")'
allocate(list(6*nVTKPoints - 12))
allocate(lptr(6*nVTKPoints - 12))
allocate(lend(nVTKPoints))
allocate(near(nVTKPoints))
allocate(next(nVTKPoints))
allocate(dist(nVTKPoints))

call TRMESH(nVTKPoints, vertexXYZ(1,:), vertexXYZ(2,:), vertexXYZ(3,:), &
			list, lptr, lend, lnew, near, next, dist, errCode)
print '("Delaunay triangulation complete.")'
allocate(sourceGrad(3,nVTKPoints))
allocate(sigma(6*nVTKPoints-12))
allocate(sourceData(nVTKPoints))
print '("ssrfpack memory ready.")'
print '(" size(sourcedata) = ", I8, " nVTKPoints = ", I8 )', size(sourceData), nVTKPoints

sigTol = 0.01_kreal
startTriangle = 1


print '("scalarCount = ",I8)', scalarCount

do k=1,scalarCount
	print '("size source data = ", I8 ," nVTKPoints = ", I8 )', size(sourceData), nVTKPoints
	sourceData = pointsScalars(:,k)
	print '("source data ready")'
	do j=1,nVTKPoints
		call GRADL(nVTKPoints,j,vertexXYZ(1,:),vertexXYZ(2,:),vertexXYZ(3,:), sourceData, &
				   list, lptr, lend, sourceGrad(:,j),errCode)
	enddo
	print '(" gradients done.")'
	call GETSIG(nVTKPoints,vertexXYZ(1,:),vertexXYZ(2,:), vertexXYZ(3,:), sourceData, &
				list, lptr, lend, sourceGrad, sigTol, sigma, dsig, errCode)


	do j=1,nLon
		do i=1,nLat
			call INTRC1(nVTKPoints, lat(i), lon(j), vertexXYZ(1,:), vertexXYZ(2,:), vertexxYZ(3,:), &
						sourceData, list, lptr, lend, SIGMA_FLAG, sigma, GRAD_FLAG, sourceGrad, startTriangle, llGrid(i,j), errCode)
		enddo
	enddo
	print '("... scalar interp complete ...")'
	! Get scalar name from dataFileHeader
	scalarName = adjustl(dataFileHeader(49+(k-1)*24:48+k*24))
	write(scalarOutputFile,'(A,A,A,A)') trim(outputPrefix),'_',trim(scalarName),'.dat'
	! Write interpolated data to output file
	open(unit=llWriteUnit,file=scalarOutputFile,status='REPLACE',action='WRITE',iostat=llWriteStat)
	if (llWriteStat /= 0 ) then
		print '("WRITE ERROR opening outputScalarFile.")'
		stop
	endif

	do i=1,nlat
		do j=1,nLon-1
			write(llWriteUnit,'(F24.18)',advance='NO') llgrid(i,j)
		enddo
		write(llWriteUnit,'(F24.18)',advance='YES') llgrid(i,nlon)
	enddo

	close(llWriteUnit)
enddo

print '("Interpolation complete.")'
deallocate(sourceData)
deallocate(sourcegrad)
deallocate(sigma)
deallocate(list)
deallocate(lptr)
deallocate(lend)
deallocate(near)
deallocate(next)
deallocate(dist)
deallocate(lat)
deallocate(lon)
deallocate(llgrid)
deallocate(vertexXYZ)
deallocate(pointsScalars)
deallocate(cellScalars)
deallocate(cellData)


end program
