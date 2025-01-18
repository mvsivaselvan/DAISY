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
    real(kind=8), dimension(NUM_NODE_DOF) :: x, xd, xdd ! displacement, velocity, accel
end type Node_t

type NodePointer_t
    type(Node_t), pointer :: ptr => null()
end type NodePointer_t

!--------------------------------------------------------

contains

!--------------------------------------------------------

subroutine make_Node(this, ID, x, y, z, RJ, constraints)

type(Node_t), intent(out) :: this
integer(kind=4), intent(in) :: ID
real(kind=8), intent(in) :: x, y, z
real(kind=8), dimension(9), intent(in) :: RJ
integer(kind=4), intent(in) :: constraints(NUM_NODE_DOF)

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

end subroutine make_Node

!--------------------------------------------------------

subroutine destroy_Node(this)
type(Node_t), intent(inout) :: this

! Nothing to deallocate

end subroutine destroy_Node

!--------------------------------------------------------

subroutine setDOF(this, x, xd, xdd)

type(Node_t), intent(inout) :: this
real(kind=8), dimension(:), intent(in) :: x, xd, xdd

integer(kind=4) :: n

do n = 1,NUM_NODE_DOF
    if (this%constraints(n) .eq. 0) then
        this%x(n) = x(this%equationNumber(n))
        this%xd(n) = xd(this%equationNumber(n))
        this%xdd(n) = xdd(this%equationNumber(n))
    endif
enddo

end subroutine setDOF

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