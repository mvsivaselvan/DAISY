module Domain_mod
    
use, intrinsic :: iso_c_binding
    
use Node_mod
use Element_mod

use stlMap_mod

implicit none

integer(kind=4), private :: numNodes, numElements
type(c_ptr), private :: nodeMap$, elementMap$

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

type(c_ptr) :: iter$
type(ElementPointer_t), pointer :: elemptr

iter$ = get_begin_iterator(elementMap$)
call delete_iterator(iter$)

call c_f_pointer(get_val_for_key(elementMap$, 1001), elemptr)

call delete_map(nodeMap$)
call delete_map(elementMap$)

end subroutine destroy_domain

!--------------------------------------------------------

subroutine startDomainBuild
! Do nothing
end subroutine startDomainBuild

!--------------------------------------------------------

subroutine endDomainBuild
! Do nothing
end subroutine endDomainBuild

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