module Element_mod

use Node_mod, only : NodePointer_t

implicit none

type Offset_t
    real(kind=8), dimension(3) :: vector = (/0.0,0.0,0.0/)
end type Offset_t
    
type, abstract :: Element_t
    integer(kind=4) :: ID
    integer(kind=4) :: numNodes ! number of nodes in the element
    integer(kind=4), dimension(:), allocatable :: nodeIDlist
    integer(kind=4), dimension(:), allocatable :: rigidoffsetIDlist
    type(NodePointer_t), dimension(:), allocatable :: nodes ! pointers to nodes 
    type(Offset_t), dimension(:), allocatable :: offsets
    integer(kind=4) :: numEqs ! number of element-level DOF
    integer(kind=4), dimension(:), allocatable :: equationNumber
    real(kind=8), dimension(:), allocatable :: x, xd, xdd ! State
    logical :: active
contains
    procedure :: make_element
    procedure :: destroy_element
    procedure(setState_element), deferred :: setState
    procedure :: setDOF
end type Element_t

abstract interface
    subroutine setState_element(this, dyn, u, x, xd, xdd)
        import :: Element_t
        class(Element_t), intent(inout) :: this
        integer, intent(in) :: dyn
        real(kind=8), dimension(3), intent(in) :: u
        real(kind=8), dimension(:), intent(in) :: x
        real(kind=8), dimension(:), intent(in), optional :: xd, xdd
    end subroutine setState_element
end interface

type ElementPointer_t
    class(Element_t), pointer :: ptr => null()
end type ElementPointer_t

contains

!--------------------------------------------------------

subroutine make_element(this, ID, numNodes, nodeIDlist, rigidoffsetIDList, &
                        numEqs, active)
class(Element_t), intent(inout) :: this
integer(kind=4), intent(in) :: ID
integer(kind=4), intent(in) :: numNodes
integer(kind=4), dimension(numNodes), intent(in) :: nodeIDlist
integer(kind=4), dimension(numNodes), intent(in) :: rigidOffsetIDlist
integer(kind=4), intent(in) :: numEqs
logical, intent(in) :: active

this%ID = ID
this%numNodes = numNodes

allocate(this%nodeIDlist(numNodes))
this%nodeIDlist = nodeIDlist
allocate(this%nodes(numNodes))

allocate(this%rigidoffsetIDlist(numNodes))
this%rigidoffsetIDlist = rigidoffsetIDlist
allocate(this%offsets(numNodes))                        

this%numEqs = numEqs
if (numEqs .gt. 0) then
    allocate(this%equationNumber(numEqs))
    allocate(this%x(numEqs))
    allocate(this%xd(numEqs))
    allocate(this%xdd(numEqs))
endif

this%active = active

end subroutine make_element
    
!--------------------------------------------------------

subroutine destroy_element(this)

class(Element_t), intent(inout) :: this

if (allocated(this%nodeIDlist)) deallocate(this%nodeIDlist)
if (allocated(this%nodes)) deallocate(this%nodes)
if (allocated(this%rigidoffsetIDlist)) deallocate(this%rigidoffsetIDlist)
if (allocated(this%offsets)) deallocate(this%offsets)
if (allocated(this%equationNumber)) deallocate(this%equationNumber)
if (allocated(this%x)) deallocate(this%x)
if (allocated(this%xd)) deallocate(this%xd)
if (allocated(this%xdd)) deallocate(this%xdd)

end subroutine destroy_element
    
!--------------------------------------------------------

subroutine setDOF(this, x, xd, xdd)

class(Element_t), intent(inout) :: this
real(kind=8), dimension(:), intent(in) :: x, xd, xdd

if (this%numEqs .eq. 0) return

this%x = x(this%equationNumber)
this%xd = xd(this%equationNumber)
this%xdd = xdd(this%equationNumber)

end subroutine setDOF

!--------------------------------------------------------
    
end module Element_mod