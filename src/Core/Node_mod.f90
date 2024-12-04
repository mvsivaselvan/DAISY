module Node

implicit none

integer, parameter :: NUM_NODE_DOF = 6 ! number of DOF per node

!===============================================

type Node_t
    integer(kind=4) :: ID
    real(kind=8), dimension(3) :: pos ! reference position
    integer(kind=4) :: constraints(NUM_NODE_DOF)
    integer(kind=4) :: equationNumber(NUM_NODE_DOF)
    integer(kind=4) :: reactNumber(NUM_NODE_DOF)
    double precision, dimension(NUM_NODE_DOF) :: x, xd, xdd ! displacement, velocity, accel
    logical :: active
end type Node_t

!===============================================

contains

!===============================================

subroutine make_Node(this, ID, x, y, z, constraints, active)

type(Node_t), intent(out) :: this
integer(kind=4), intent(in) :: ID
real(kind=8), intent(in) :: x, y, z
integer(kind=4), intent(in) :: constraints(NUM_NODE_DOF)
logical, intent(in) :: active

this%ID = ID
this%pos(1) = x
this%pos(2) = y
this%pos(3) = z
this%constraints = constraints
this%equationNumber = 0
this%reactNumber = 0
this%x = 0.d0
this%xd = 0.d0
this%xdd = 0.d0
this%active = active

end subroutine make_Node

!===============================================

subroutine destroy_Node(this)
type(Node_t), intent(out) :: this

this%active = .false.

end subroutine destroy_Node

!============================================================

!============================================================
!    Utility functions
!============================================================

function distance_between(this, that, deformed) result(dist)

type(Node_t), intent(in) :: this, that
logical, intent(in) :: deformed ! deformed or reference distance
real(kind=8) :: dist

real(kind=8), dimension(3) :: x1, x2

x1 = this%pos
x2 = that%pos
if (deformed) then ! add displacement
    x1 = x1 + this%x(1:3)
    x2 = x2 + that%x(1:3)
endif

dist = sqrt((x2(1)-x1(1))**2.d0+(x2(2)-x1(2))**2.d0+(x2(3)-x1(3))**2.d0)

end function distance_between
	
!===============================================

end module Node