module Node_mod

implicit none

integer, parameter :: NUM_NODE_DOF = 6 ! number of DOF per node

!--------------------------------------------------------

type Node_t
    integer(kind=4) :: ID
    real(kind=8), dimension(3) :: x0 ! reference position
    real(kind=8), dimension(9) :: RJ ! joint coordinate frame orientation
    integer(kind=4) :: constraints(NUM_NODE_DOF)
    integer(kind=4) :: equationNumber(NUM_NODE_DOF)
    integer(kind=4) :: reactNumber(NUM_NODE_DOF)
    double precision, dimension(NUM_NODE_DOF) :: x, xd, xdd ! displacement, velocity, accel
    logical :: active
end type Node_t

type NodePointer_t
    type(Node_t), pointer :: ptr
end type NodePointer_t

!--------------------------------------------------------

contains

!--------------------------------------------------------

subroutine make_Node(this, ID, x, y, z, RJ, constraints, active)

type(Node_t), intent(out) :: this
integer(kind=4), intent(in) :: ID
real(kind=8), intent(in) :: x, y, z
real(kind=8), dimension(9), intent(in) :: RJ
integer(kind=4), intent(in) :: constraints(NUM_NODE_DOF)
logical, intent(in) :: active

this%ID = ID
this%x0(1) = x
this%x0(2) = y
this%x0(3) = z
this%RJ = RJ
this%constraints = constraints
this%equationNumber = 0
this%reactNumber = 0
this%x = 0.d0
this%xd = 0.d0
this%xdd = 0.d0
this%active = active

end subroutine make_Node

!--------------------------------------------------------

subroutine destroy_Node(this)
type(Node_t), intent(out) :: this

this%active = .false.

end subroutine destroy_Node

!--------------------------------------------------------
!    Utility functions
!--------------------------------------------------------

function distance_between(this, that, deformed) result(dist)

type(Node_t), intent(in) :: this, that
logical, intent(in) :: deformed ! deformed or reference distance
real(kind=8) :: dist

real(kind=8), dimension(3) :: x1, x2

x1 = this%x0
x2 = that%x0
if (deformed) then ! add displacement
    x1 = x1 + this%x(1:3)
    x2 = x2 + that%x(1:3)
endif

dist = sqrt((x2(1)-x1(1))**2.d0+(x2(2)-x1(2))**2.d0+(x2(3)-x1(3))**2.d0)

end function distance_between
	
!--------------------------------------------------------

end module Node_mod