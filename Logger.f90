module LoggerModule
!	Based on : McCormack, D. "Scientific Software Development with Fortran," 2009.
use NumberKindsModule
use OutputWriterModule

implicit none
private

public Logger
public New, Delete
public DEBUG_LOGGING_LEVEL, TRACE_LOGGING_LEVEL, WARNING_LOGGING_LEVEL, ERROR_LOGGING_LEVEL
public LogMessage
public GetLoggingLevel
public StartSection, EndSection

integer(kint), parameter :: DEBUG_LOGGING_LEVEL = 1, &
							TRACE_LOGGING_LEVEL = 2, &
							WARNING_LOGGING_LEVEL = 3, &
							ERROR_LOGGING_LEVEL = 4

type Logger
	integer(kint) :: level
	type(OutputWriter) :: writer
	character(len=256) :: logFileName
	logical(klog) :: outputToFile
end type	

interface New
	module procedure NewPrivate
end interface

interface Delete
	module procedure DeletePrivate
end interface

interface LogMessage
	module procedure LogMessageString
	module procedure LogMessageInteger
	module procedure LogMessageReal
end interface

interface StartSection
	module procedure StartSectionLogger
end interface

interface EndSection
	module procedure EndSectionLogger
end interface

contains

function GetLoggingLevel(self)
	type(Logger), intent(in) :: self
	integer(kint) :: GetLoggingLevel
	GetLoggingLevel = self%level
end function


subroutine NewPrivate(self,level,fileUnit,filename)
	type(Logger), intent(out) :: self
	integer(kint), intent(in) :: level, fileUnit
	character(len=*), intent(in), optional :: filename
	integer(kint) :: logStat, stdOut
	
	stdOut = 6
	self%level = level
	
	self%outputToFile = .False.
	
	if ( present(filename) ) then
		self%logfileName = filename
		if ( level <= 6 ) then
			print *, "Logger WARNING : File output requires a fileUnit > 6."
		else	
			self%outputToFile = .True.
			open(unit=fileUnit,file=filename,status='REPLACE',action='WRITE',iostat=logStat)
			if ( logStat /= 0 ) then
				print *,"ERROR opening log file : ",filename
				self%outputToFile = .False.
			endif
		endif
	else
		self%logfileName = ' '
	endif
	
	if ( self%outputToFile ) then
		call New(self%writer,fileUnit)
	else
		call New(self%writer,stdOut)	
	endif
	
end subroutine


subroutine DeletePrivate(self)
	type(Logger), intent(inout) :: self
	call Delete(self%writer)
end subroutine


subroutine LogMessageString(self,msgLevel,key,str)
	type(Logger), intent(in) :: self
	integer(kint), intent(in) :: msgLevel
	character(len=*), intent(in) :: key,str
	if ( msgLevel >= self%level ) then
		call Write(self%writer,key,str)
	endif
end subroutine


subroutine LogMessageInteger(self,msgLevel,key,int)
	type(Logger), intent(in) :: self
	integer(kint), intent(in) :: msgLevel, int
	character(len=*), intent(in) :: key
	if ( msgLevel >= self%level) then
		call Write(self%writer,key,int)
	endif
end subroutine


subroutine LogMessageReal(self,msgLevel,key,float)
	type(Logger), intent(in) :: self
	integer(kint), intent(in) :: msgLevel
	character(len=*), intent(in) :: key
	real(kreal), intent(in) :: float
	if ( msgLevel >= self%level ) then
		call Write(self%writer,key,float)
	endif
end subroutine


subroutine StartSectionLogger(self,sectionName,description)
	type(Logger), intent(inout) :: self
	character(len=*), intent(in) :: sectionName
	character(len=*), intent(in), optional ::  description
	call StartSection(self%writer,sectionName,description)
end subroutine


subroutine EndSectionLogger(self)
	type(Logger), intent(inout) :: self
	call EndSection(self%writer)
end subroutine

end module
