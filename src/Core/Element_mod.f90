module Element_mod

use Node_mod, only : NodePointer_t

implicit none
    
type, abstract :: Element_t
    integer(kind=4) :: ID
    integer(kind=4) :: numNodes ! number of nodes in the element
    integer(kind=4) :: nodeIDlist
    type(NodePointer_t), allocatable, dimension(:) :: nodes ! pointers to nodes 
    integer(kind=4) :: numEqs ! number of element-level DOF
    integer(kind=4), allocatable, dimension(:) :: equationNumber
    logical :: active
contains
    procedure :: make_element
    procedure :: destroy_element
    procedure(setState_element), deferred :: setState
end type Element_t

abstract interface
    subroutine setState_element(this, dyn, u, x, xd, xdd)
        import :: Element_t
        class(Element_t), intent(out) :: this
        integer, intent(in) :: dyn
        real(kind=8), dimension(3), intent(in) :: u
        real(kind=8), dimension(:), intent(in) :: x
        real(kind=8), dimension(:), intent(in), optional :: xd, xdd
    end subroutine setState_element
end interface


contains

!--------------------------------------------------------

subroutine make_element(this, ID, numNodes, nodeIDlist, numEqs)

class(Element_t), intent(out) :: this
integer(kind=4), intent(in) :: ID
integer(kind=4), intent(in) :: numNodes
integer(kind=4), dimension(:), intent(in) :: nodeIDlist
integer(kind=4), intent(in) :: numEqs

this%ID = ID
this%numNodes = numNodes
allocate(this%nodes(numNodes))
this%numEqs = numEqs
allocate(this%equationNumber(numEqs))

end subroutine make_element
    
!--------------------------------------------------------

subroutine destroy_element(this)

class(Element_t), intent(out) :: this

if (allocated(this%nodes)) deallocate(this%nodes)
if (allocated(this%equationNumber)) deallocate(this%equationNumber)

end subroutine destroy_element
    
!--------------------------------------------------------
    
end module Element_mod