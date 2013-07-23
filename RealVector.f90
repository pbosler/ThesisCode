module DoubleVectorModule

use NumberKindsModule

implicit none
private

public DoubleVector, New, Delete
public AddRealDataToEnd, AddRealDataToFront

type RealVector
	real(kreal), pointer :: data
	integer(kint) :: n
end type

interface New
	module procedure NewPrivate
end interface

interface Delete
	module procedure DeletePrivate
end interface

interface AddRealDataToEnd
	module procedure AddRealToEnd
	module procedure AddRealArrayToEnd
end interface

interface AddRealDataToFront
	module procedure AddRealToFront
	module procedure AddRealArrayToFront
end interface

contains

subroutine NewPrivate(self,reserve)
	type(RealVector), intent(out) :: self
	integer(kint), intent(in) :: reserve
	self%n = 0
	allocate(self%data(reserve))
	self%data = 0.0_kreal
end subroutine


subroutine DeletePrivate(self)
	type(RealVector), intent(inout) :: self
	deallocate(self%data)
	self%n
end subroutine


recursive subroutine AddRealToEnd(self,newReal)
	type(RealVector), intent(inout) :: self
	real(kreal), intent(in) :: newReal
	if ( size(self%data) >= self%n + 1) then
		self%data(self%n+1) = newReal
		self%n = self%n+1
	else
		call ResizeRealVector(self)
		call AddRealToEnd(self,newReal)
	endif
end subroutine


recursive subroutine AddRealToFront(self,newReal)
	type(RealVector), intent(inout) :: self
	real(kreal), intent(in) :: newReal
	integer(kint) :: i, n
	real(kreal), allocatable :: tempData(:)
	if ( size(self%data) >= self%n+1) then
		n = self%n
		allocate(tempData(n))
		tempData = self%data(1:n)
		self%data(1) = newReal
		do i=1,n
			self%data(i+1) = tempData(i)
		enddo
		self%n = n+1
		deallocate(tempData)
	else
		call ResizeRealVector(self)
		call AddRealToFront(self,newReal)
	endif
end subroutine


recursive subroutine AddRealArrayToEnd(self,newReals)
	type(RealVector), intent(inout) :: self
	real(kreal), intent(in) :: newReals(:)
	integer(kint) :: i, m
	m = size(newReals)
	if ( size(self%data) >= self%n + m ) then
		do i=1,m
			self%data(self%n+i) = newReals(i)
		enddo
		self%n = self%n+m
	else
		call ResizeRealVector(self)
		call AddRealArrayToEnd(self,newReals)
	endif
end subroutine


subroutine AddRealArrayToFront(self,newReals)
	type(RealVector), intent(inout) :: self
	real(kreal), intent(in) :: newReals(:)
	integer(kint) :: i,n,m
	real(kreal), allocatable :: tempData(:)
	m = size(newReals)
	if ( size(self%data) >= self%n + m) then
		allocate(tempData(self%n))
		tempData = self%data(1:self%n)
		self%data(1:m) = newReals
		self%data(m+1:m+self%n) = tempData
		deallocate(tempData)	
	else
		call ResizeRealVector(self)
		call AddRealArrayToFront(self,newReals)
	endif
end subroutine


subroutine ResizeRealVector(self)
	type(RealVector), intent(inout) :: self
	real(kreal), allocatable :: tempData(:)
	allocate(tempData(self%n))
	tempData = self%data(1:self%n)
	deallocate(self%data)
	allocate(self%data(2*self%n))
	self%data(1:self%n) = tempData
	deallocate(tempData)
end subroutine


end module