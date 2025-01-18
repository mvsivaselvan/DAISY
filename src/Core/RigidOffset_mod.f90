module RigidOffset_mod

implicit none
    
type RigidOffset_t
    integer(kind=4) :: ID
    real(kind=8), dimension(3) :: offset = (/0.0,0.0,0.0/)
end type RigidOffset_t

contains
    
!-------------------------------------------

subroutine make_RigidOffset(this, ID, x, y, z)

type(RigidOffset_t), intent(out) :: this
integer(kind=4), intent(in) :: ID
real(kind=8), intent(in) :: x, y, z

this%ID = ID
this%offset(1) = x
this%offset(2) = y
this%offset(3) = z

end subroutine make_RigidOffset

!-------------------------------------------

subroutine destroy_RigidOffset(this)

type(RigidOffset_t), intent(inout) :: this

! Nothing to do

end subroutine destroy_RigidOffset

!-------------------------------------------
    
end module RigidOffset_mod