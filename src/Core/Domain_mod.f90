module Domain_mod
    
use, intrinsic :: iso_c_binding
    
use Node_mod
use Element_mod

use stlMap_mod

implicit none

integer(kind=4), private :: numNodes, numElements
! For associative access
type(c_ptr), private :: nodeMap$, elementMap$ 
! For sequential access
type(NodePointer_t), allocatable :: nodearray(:)
type(ElementPointer_t), allocatable :: elemarray(:)

contains

!--------------------------------------------------------
!                      Build Routines
!--------------------------------------------------------

subroutine make_domain

numNodes = 0 
numElements = 0

nodeMap$ = create_map()
elementMap$ = create_map()

end subroutine make_domain

!--------------------------------------------------------

subroutine destroy_domain

! There were put in here simply to test this functionality
! and kept here for memory, has to be deleted
!type(c_ptr) :: iter$
!type(ElementPointer_t), pointer :: elemptr
!iter$ = get_begin_iterator(elementMap$)
!call delete_iterator(iter$)
!call c_f_pointer(get_val_for_key(elementMap$, 1001), elemptr)

integer(kind=4) :: i

do i=1,numNodes
    call destroy_Node(nodearray(i)%ptr)
    deallocate(nodearray(i)%ptr)
enddo
if (allocated(nodearray)) deallocate(nodearray)
call delete_map(nodeMap$)

call delete_map(elementMap$)

end subroutine destroy_domain

!--------------------------------------------------------

subroutine startDomainBuild
! Do nothing
end subroutine startDomainBuild

!--------------------------------------------------------

subroutine endDomainBuild

type(c_ptr) :: iterator$, node$

integer(kind=4) :: i

! Allocate node array for iterative access
if (numNodes .gt. 0) then
    allocate( nodearray(numNodes) )
    ! Add nodes to array container
    iterator$ = get_begin_iterator(nodeMap$)
    do i=1,numNodes
        call c_f_pointer(get_value_for_iterator(iterator$), nodearray(i)%ptr)
        call iterate_next(iterator$)
    enddo
    call delete_iterator(iterator$)
endif

!if (numElements .gt. 0) allocate( elemarray(numElements) )


end subroutine endDomainBuild

!--------------------------------------------------------

subroutine add_node_to_domain(node_)

type(Node_t), target :: node_

numNodes = numNodes + 1
call add_to_map(nodeMap$, node_%ID, c_loc(node_))

end subroutine add_node_to_domain

!--------------------------------------------------------

subroutine add_element_to_domain(elem_)

class(Element_t), target :: elem_

type(ElementPointer_t), pointer :: elemptr
allocate(elemptr)
elemptr%ptr => elem_

numElements = numElements + 1
call add_to_map(elementMap$, elem_%ID, c_loc(elemptr))

end subroutine add_element_to_domain

!--------------------------------------------------------
    
end module Domain_mod