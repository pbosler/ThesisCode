program FindPotVort

use NumberKindsModule
use SphereGeomModule
use TracerAndVorticityDistributionModule
use LoggerModule

implicit none

! File Processing
character(len=128) :: namelistFile
character(len=256) :: vtkDataFile,  dataFileHeader, outputPrefix
character(len=256) :: outputFile
character(len=128) :: lineIn, scratch1, scratch2
integer(kint), parameter :: readUnit = 11
integer(kint), parameter :: writeUnit = 12
real(kreal), allocatable :: vertexXYZ(:,:), absVort(:), potVort(:)
integer(kint) :: nPoints, lineNumber,  scalarCount
integer(kint) :: readStat, writeStat
integer(kint) :: nVTKPoints, nVTKCells, polygonHeaderLine, pointDataLine, cellDataLine
integer(kint) :: absVortLine
character(len=56) :: scalarName, foundName, formatString
integer(kint) :: i, j, k, i1, i2, i3, i4, i5, i6, i7
logical(klog) :: keepGoing
real(kreal) :: programStart, programEnd, time, lat, lon, norm

type(Logger) :: log

namelist /inputdata/ vtkDataFile, outputPrefix, time

call cpu_time(programStart)

print *, "This program takes particles and panels data from a stratModel .vtk file and "
print *, "converts absolute vorticity to potential vorticity."
print *, "User input via 'findPotVort.namelist' file."

i = 6
call New(log,DEBUG_LOGGING_LEVEL,i)

namelistFile = 'findPotVort.namelist'
! Get user input for filenames
open(unit=readUnit,file=namelistFile,status='OLD',delim='APOSTROPHE')
	read(readUnit,nml=inputdata)
close(readUnit)

write(outputFile,'(A,A,F7.4,A)') trim(outputPrefix),'potVort_t',time,'.dat'

! Read VTK Metadata
open(unit=readUnit,file=vtkDataFile,status='OLD',action='READ',iostat=readStat)
if ( readStat /= 0 ) stop "ERROR opening vtk input file."
!
! Preprocessing : Find number of points, number of cells, and number of scalars in data file
!
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
				if ( i2 == 1) then 
					scalarCount = scalarCount + 1
					i1 = scan(lineIn,' ')
					scratch1 = lineIn(i1+1:len_trim(lineIn))
					scratch1 = adjustl(scratch1)
					i3 = scan(scratch1,' ')
					foundName = scratch1(1:i3+1)
					i4 = scan(foundName,'a')
					i5 = scan(foundName,'b')
					i6 = scan(foundName,'s')
					i7 = scan(foundName,'P')
					! FIND ABS VORT SCALAR DATA
					if ( i4 + i5 + i6 + i7 == 6 ) then ! absVort data found
						absVortLine = lineNumber
					endif
				endif
			endif
		endif
	endif
enddo

print *, " File contains ", scalarCount, " scalars."
print *, " polygonDataLine = ", polygonHeaderLine
print *, " pointDataLine = ",pointDataLine
print *, " cellDataLine = ", cellDataLine
print *, " absVortLine = ", absVortLine

allocate(vertexXYZ(3,nVTKPoints))
vertexXYZ = 0.0_kreal
allocate(absVort(nVTKPoints))
absVort = 0.0_kreal
allocate(potVort(nVTKPoints))
potVort = 0.0_kreal

call LogMessage(log,TRACE_LOGGING_LEVEL,'status: ', "Reading VTK File Data ...")
!
!	Read Data
!
rewind(readUnit)
keepGoing = .TRUE.
lineNumber = 0
j=1
k=1
do while (keepGoing)
	read(readUnit,'(A)',iostat=readStat) lineIn	! Read entire line of file
	lineNumber = lineNumber + 1
	if ( readStat > 0 ) then
		call LogMessage(log,ERROR_LOGGING_LEVEL,'Read ERROR at line ',lineNumber)
		keepGoing = .FALSE.
		stop
	elseif (readStat < 0 ) then
		call LogMessage(log,TRACE_LOGGING_LEVEL,'End of file found at line',lineNumber)
		keepGoing = .FALSE.
	else
		if ( lineNumber > 5 .AND. lineNumber < polygonHeaderLine ) then	
		! READ VERTEX COORDINATES
			read(lineIn,'(3D)') vertexXYZ(:,j)
			! DEBUG
			norm = sqrt(sum(vertexXYZ(:,j)*vertexXYZ(:,j)))
!			if ( abs(norm - 1.0_kreal) > 1.0E-12) then
!				call LogMessage(log,WARNING_LOGGING_LEVEL,"WARNING : norm = ",norm)
!				call LogMessage(log,WARNING_LOGGING_LEVEL," at lineNumber = ",lineNumber)
!				call LogMessage(log,WARNING_LOGGING_LEVEL," and j = ",j)
!			endif
			! END DEBUG
			j=j+1
		endif
		if ( lineNumber == polygonHeaderLine) j=1
		if ( lineNumber > absVortLine + 1 .AND. lineNumber < absVortLine + nVTKPoints + 1) then
			read(lineIn,'(D)') absVort(j)
			j=j+1
		endif
		if (lineNumber == absVortLine + nVTKPoints + 1) keepGoing = .FALSE.
	endif
enddo
call LogMessage(log,TRACE_LOGGING_LEVEL,'status: ', "... done.")

call LogMessage(log,TRACE_LOGGING_LEVEL,'status: ', "Converting absVort to potVort...")
	do j=1,nVTKPoints
		potVort(j) = absVort(j) + Juckes_Forcing(vertexXYZ(:,j),time)
!		if ( isnan(potVort(j)) ) then
!			call LogMessage(log,ERROR_LOGGING_LEVEL,"ERROR: Nan found at lineNumber = ",lineNumber)
!			print '("particle coordinates = ",3F24.18)', vertexXYZ(:,j)
!			print '("Lat, Lon = ",2F24.18)', Latitude(vertexXYZ(:,j)), Longitude(vertexXYZ(:,j))
!			print '("absVort = ",F24.18)', absVort(j)
!			print '("F(x,t)=",F24.18)', Juckes_Forcing(vertexXYZ(:,j),time)
!			print '("lineNumber = ",I8)', lineNumber
!			print '("j = ", I8)', j
!		endif
	enddo
call LogMessage(log,TRACE_LOGGING_LEVEL,'status: ','...done.')

open(unit=writeUnit,file=outputFile,status='REPLACE',action='WRITE',iostat=writeStat)
if (writeStat /= 0 ) then
	call LogMessage(log,ERROR_LOGGING_LEVEL,'ERROR :',' cannot open output file.')
	stop
endif
write(writeUnit,'(A,I8)') 'POINT_DATA  ',nVTKPoints
write(writeUnit,'(A)') 'SCALARS potVort double 1'
write(writeUnit,'(A)') 'LOOKUP_TABLE default'
do j=1,nVTKPoints
	write(writeUnit,'(F24.16)') potVort(j)
enddo
close(writeUnit)

deallocate(vertexXYZ)
deallocate(potVort)
deallocate(absVort)
call Delete(log)
end program