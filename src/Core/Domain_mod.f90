module Domain_mod
    
use, intrinsic :: iso_c_binding
    
use Node_mod
use Element_mod
use RigidOffset_mod

use stlMap_mod

implicit none

real(kind=8), dimension(3) :: gravity_u = (/0.0, 0.0, 0.0/)

! Nodes and elements
integer(kind=4), private :: numNodes, numElements
! For associative access
type(c_ptr), private :: nodeMap$, elementMap$
! For sequential access
type(NodePointer_t), allocatable, private :: nodearray(:)
type(ElementPointer_t), allocatable, private :: elemarray(:)

! Rigid Offsets
integer(kind=4), private :: numRigidOffsets
type(c_ptr), private :: rigidoffsetMap$

! Boundary conditions, DOF constraints and Element sets
integer(kind=4), private :: numBCs, numElemSets, numDofConstraints
type(c_ptr), private :: bcMap$, elementSetMap$, dofConstraintMap$

contains

!--------------------------------------------------------
!                      Build Routines
!--------------------------------------------------------

subroutine make_domain

numNodes = 0
nodeMap$ = create_map()

numElements = 0
elementMap$ = create_map()

numRigidOffsets = 0
rigidoffsetMap$ = create_map()

numBCs = 0
bcMap$ = create_map()

numElemSets = 0
elementSetMap$ = create_map()

numDofConstraints = 0
dofConstraintMap$ = create_map()

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
do i=1,numElements
    deallocate(elemarray(i)%ptr)
enddo
if (allocated(elemarray)) deallocate(elemarray)
call delete_map(elementMap$)

! Rigid offsets
call destroy_RigidOffsetMap(rigidoffsetMap$, numRigidOffsets)

! Boundary conditions

! Element sets

! DOF constraints

end subroutine destroy_domain

!--------------------------------------------------------
! Generic construction of routines to destroy maps

#define SPECIFIC_DESTROY_MAP destroy_RigidOffsetMap
#define SPECIFIC_COMPONENT_TYPE RigidOffset_t
#define SPECIFIC_DESTROY_ROUTINE destroy_RigidOffset
#include "destroy_map_generic.inc"
#undef SPECIFIC_DESTROY_MAP
#undef SPECIFIC_COMPONENT_TYPE
#undef SPECIFIC_DESTROY_ROUTINE

!--------------------------------------------------------

subroutine startDomainBuild
! Do nothing
end subroutine startDomainBuild

!--------------------------------------------------------

subroutine endDomainBuild

type(c_ptr) :: iterator$, node$
integer(kind=4) :: i, nodeind
type(ElementPointer_t), pointer :: tempElemptr
type(RigidOffset_t), pointer :: offsetptr

! Add nodes to array container
if (numNodes .gt. 0) then
    allocate( nodearray(numNodes) )
    iterator$ = get_begin_iterator(nodeMap$)
    do i=1,numNodes
        call c_f_pointer(get_value_for_iterator(iterator$), nodearray(i)%ptr)
        call iterate_next(iterator$)
    enddo
    call delete_iterator(iterator$)
endif

! Add elements to array container, and associate with nodes
! and rigid offsets
if (numElements .gt. 0) then
    allocate( elemarray(numElements) )
    iterator$ = get_begin_iterator(elementMap$)
    do i=1,numElements
        call c_f_pointer(get_value_for_iterator(iterator$), tempElemptr)
        elemarray(i)%ptr => tempElemptr%ptr
        do nodeind = 1, elemarray(i)%ptr%numNodes
            call c_f_pointer(get_val_for_key(nodeMap$,elemarray(i)%ptr%nodeIDlist(nodeind)), &
                             elemarray(i)%ptr%nodes(nodeind)%ptr)
            call c_f_pointer(get_val_for_key(rigidoffsetMap$,&
                      elemarray(i)%ptr%rigidoffsetIDlist(nodeind)), offsetptr)
            elemarray(i)%ptr%offsets(nodeind)%vector = offsetptr%offset
        enddo
        call iterate_next(iterator$)
    enddo
    call delete_iterator(iterator$)
endif

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