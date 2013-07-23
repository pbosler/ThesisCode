module OutputWriterModule
	!	Based on : McCormack, D. "Scientific Software Development with Fortran," 2009.

	use NumberKindsModule
	
	implicit none
	private
	
	public OutputWriter
	public New, Delete
	public StartSection, EndSection
	public Write
	
	type OutputWriter
		private
		integer(kint) :: fileUnit
		integer(kint) :: indentLevel
	end type
	
	interface New
		module procedure NewPrivate
	end interface
		
	interface Delete
		module procedure DeletePrivate
	end interface	
	
	interface Write
		module procedure WriteString
		module procedure WriteInteger
		module procedure WriteReal
	end interface
	
	interface StartSection
		module procedure StartSectionWriter
	end interface

	interface EndSection
		module procedure EndSectionWriter
	end interface
contains

	subroutine NewPrivate(self, fileUnit)
		type(OutputWriter), intent(out) :: self
		integer(kint), intent(in), optional :: fileUnit
		if ( present(fileUnit) ) then
			self%fileUnit = fileUnit
		else
			self%fileUnit = 6
		endif
		self%indentLevel = 0
	end subroutine

	
	subroutine DeletePrivate(self)
		type(OutputWriter), intent(inout) :: self
	end subroutine


function FormatWithIndent(self,formatString)
	type(OutputWriter), intent(in) :: self
	character(len=*), intent(in) :: formatString
	character(len=len(formatString)+10) :: FormatWithIndent
	if ( self%indentLevel > 0 ) then
		write(FormatWithIndent,'(A,I2,2A)') '(',4*self%indentLevel,'X,',formatString(2:)
	else
		FormatWithIndent = formatString
	endif
end function

	
subroutine WriteString(self,key,str)
	type(OutputWriter), intent(in) :: self
	character(len=*), intent(in) :: key, str
	character(len=32) :: form
	form = FormatWithIndent(self,'(A,2X,A)')
	write(self%fileUnit,form) trim(key), str
end subroutine	
	
subroutine WriteInteger(self,key,val)
	type(OutputWriter), intent(in) :: self
	character(len=*), intent(in) :: key
	integer(kint), intent(in) :: val
	character(len=32) :: form
	form = FormatWithIndent(self,'(A,2X,I8)')
	write(self%fileUnit,form) trim(key), val
end subroutine


subroutine WriteReal(self,key,val)
	type(OutputWriter), intent(in) :: self
	character(len=*), intent(in) :: key
	real(kreal), intent(in) :: val
	character(len=32) :: form
	form = FormatWithIndent(self,'(A,2X,F24.15)')
	write(self%fileUnit,form) trim(key), val
end subroutine


subroutine StartSectionWriter(self,sectionName,description)
		type(OutputWriter), intent(inout) :: self
		character(len=*), intent(in) :: sectionName
		character(len=*), intent(in), optional ::  description
		character(len=32) :: form
		form = FormatWithIndent(self,'(A)')
		write(self%fileUnit,form) sectionName
		self%indentLevel = self%indentLevel + 1
		if ( present(description) ) then
			form = FormatWithIndent(self,'(A)')
			write(self%fileUnit,form) description
		endif
end subroutine


subroutine EndSectionWriter(self)
		type(OutputWriter), intent(inout) :: self
		if ( self%indentLevel == 0 ) then
			print *, "EndSection WARNING : Indentation Level is already zero."
		else
			self%indentLevel = self%indentLevel - 1
		endif
end subroutine




	
end module
