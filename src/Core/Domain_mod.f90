module Domain_mod
    
use, intrinsic :: iso_c_binding
    
use Node_mod
use Element_mod
use RigidOffset_mod

use stlMap_mod

implicit none

integer(kind=4), private :: numNodes, numElements, numRigidOffsets
! For associative access
type(c_ptr), private :: nodeMap$, elementMap$, rigidoffsetMap$
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
numRigidOffsets = 0

nodeMap$ = create_map()
elementMap$ = create_map()
rigidoffsetMap$ = create_map()

end subroutine make_domain

!--------------------------------------------------------

subroutine destroy_domain

integer(kind=4) :: i

type(c_ptr) :: iterator$
type(RigidOffset_t), pointer :: offset

! Nodes
do i=1,numNodes
    call destroy_Node(nodearray(i)%ptr)
    deallocate(nodearray(i)%ptr)
enddo
if (allocated(nodearray)) deallocate(nodearray)
call delete_map(nodeMap$)

! Elements
call delete_map(elementMap$)

! Rigid offsets
if (numRigidOffsets .gt. 0) then
    iterator$ = get_begin_iterator(rigidoffsetMap$)
    do i=1,numRigidOffsets
        call c_f_pointer(get_value_for_iterator(iterator$), offset)
        call destroy_rigidoffset(offset)
        deallocate(offset)
        call iterate_next(iterator$)
    enddo
    call delete_iterator(iterator$)
endif
call delete_map(rigidoffsetMap$)

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

subroutine add_rigidoffset_to_domain(offset_)

type(RigidOffset_t), target :: offset_

numRigidOffsets = numRigidOffsets + 1
call add_to_map(rigidoffsetMap$, offset_%ID, c_loc(offset_))

end subroutine add_rigidoffset_to_domain

!--------------------------------------------------------!--------------------------------------------------------
    
end module Domain_mod