module IntegerVectorModule

use NumberKindsModule

implicit none
private
public IntegerVector, New, Delete
public AddIntDataToEnd, AddIntDataToFront

type IntegerVector
	integer(kint), pointer :: data
	integer(kint) :: n
end type

interface New
	module procedure NewPrivate
end interface

interface Delete
	module procedure DeletePrivate
end interface

interface AddIntDataToEnd
	module procedure AddIntToEnd
	module procedure AddIntArrayToEnd
end interface

interface AddIntDataToFront
	module procedure AddIntToFront
	module procedure AddIntArrayToFront
end interface

contains

subroutine NewPrivate(self,reserve)
	type(IntegerVector), intent(out) :: self
	integer(kint), intent(in) :: reserve
	self%n = 0
	allocate(self%data(reserve))
	self%data = 0
end subroutine


subroutine DeletePrivate(self)
	type(IntegerVector), intent(inout) :: self
	deallocate(self%data)
	self%n = 0
end subroutine


recursive subroutine AddIntToEnd(self,newInt)
	type(IntegerVector), intent(inout) :: self
	integer(kint), intent(in) :: newInt
	if ( size(self%data) >= self%n + 1) then
		self%data(self%n+1) = newInt
		self%n = self%n+1
	else
		call ResizeIntVector(self)
		call AddIntToEnd(self,newInt)
	endif
end subroutine


recursive subroutine AddIntToFront(self,newInt)
	type(IntegerVector), intent(inout) :: self
	integer(kint), intent(in) :: newInt
	integer(kint) :: i, n
	integer(kint), allocatable :: tempData(:)
	if ( size(self%data) >= self%n + 1) then
		n = self%n
		allocate(tempData(n))
		tempData = self%data(1:n)
		self%data(1) = newInt
		do i=1,self%n
			self%data(i+1) = tempData(i)
		enddo
		self%n = self%n+1
		deallocate(tempData)
	else
		call ResizeIntVector(self)
		call AddIntToEnd(self,newInt)
	endif
end subroutine


recursive subroutine AddIntArrayToEnd(self, newInts)
	type(IntegerVector), intent(inout) :: self
	integer(kint), intent(in) :: newInts(:)
	integer(kint) :: i,m
	m = size(newInts)
	if ( size(self%data) >= self%n+m) then
		do i = 1,m
			self%data(self%n+i) = newInts(i)
		enddo
		self%n = self%n + m
	else
		call ResizeIntVector(self)
		call AddIntArrayToEnd(self,newInts)
	endif
end subroutine


recursive subroutine AddIntArrayToFront(self, newInts)
	type(IntegerVector), intent(inout) :: self
	integer(kint), intent(in) :: newInts(:)
	integer(kint) :: i,m
	integer(kint), allocatable :: tempData(:)
	m = size(newInts)
	if ( size(self%data) >= self%n+m) then
		allocate(tempData(self%n))
		tempData = self%data(1:self%n)
		do i = 1,m
			self%data(i) = newInts(i)
		enddo
		do i=1,self%n
			self%data(i+m) = tempData(i)
		enddo
		self%n = self%n + m
		deallocate(tempData)
	else
		call ResizeIntVector(self)
		call AddIntArrayToEnd(self,newInts)
	endif
end subroutine


subroutine ResizeIntVector(self)
	type(IntegerVector), intent(inout) :: self
	integer(kint), allocatable :: tempData(:)
	integer(kint) :: n
	n = self%n
	allocate(tempData(n))
	tempData = self%data(1:n)
	deallocate(self%data)
	allocate(self%data(2*n))
	self%data(1:n) = tempData
	deallocate(tempData)
end subroutine

	
end module
