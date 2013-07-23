module IntegerListModule
! This module implements a simple linked list data structure for lists of integers

use NumberKindsModule

implicit none
private
public IntegerListNode
public New, Delete
public AddIntegerToList, ReturnIntegerArrayFromList

type IntegerListNode
	integer(kint) :: int
	integer(kint) :: listSize
	type(IntegerListNode), pointer :: next
end type

interface New
	module procedure NewPrivate
end interface

interface Delete
	module procedure DeletePrivate
end interface

contains

subroutine NewPrivate(listRoot,number)
	type(IntegerListNode), pointer, intent(out) :: listRoot
	integer(kint), intent(in), optional :: number
	allocate(listRoot)
	if ( present(number) ) then
		listRoot%int = number
		listRoot%listSize = 1
	else
		listRoot%listSize = 0
	endif
	nullify(listRoot%next)
end subroutine


recursive subroutine DeletePrivate(listRoot)
	type(IntegerListNode), pointer, intent(inout) :: listRoot
	type(IntegerListNode), pointer :: next
	next => listRoot%next
	deallocate(listRoot)
	if ( associated(next)) then
		call DeletePrivate(next)
	endif
end subroutine


subroutine AddIntegerToList(listRoot,number)
	type(IntegerListNode), pointer, intent(inout) :: listRoot
	integer(kint), intent(in) :: number
	type(IntegerListNode), pointer :: current, next
	! check special case of empty list
	if ( listRoot%listSize == 0 ) then
		listRoot%int = number
		listRoot%listSize = 1
		nullify(listRoot%next)
	else
		listRoot%listSize = listRoot%listSize + 1
		current => listRoot
		next => listRoot%next
		! find end of list
		do while (associated(next)) 
			current=>current%next
			next=>current%next
		enddo
		! allocate a new node
		allocate(next)
		! connect new node to end of list
		current%next => next
		next%int = number
		nullify(next%next)
	endif
end subroutine


subroutine ReturnIntegerArrayFromList(intArray,listRoot)
	integer(kint), intent(out) :: intArray(:)
	type(IntegerListNode), pointer, intent(in) :: listRoot
	type(IntegerListNode), pointer :: current
	integer(kint) :: j
	if ( size(intArray) /= listRoot%listSize) then
		print *, "ReturnIntegerArrayFromList ERROR : size mismatch."
		return
	endif
	current => listRoot
	j = 1
	do while (associated(current))
		intArray(j) = current%int
		current=>current%next
		j = j+1
	enddo
end subroutine

end module
