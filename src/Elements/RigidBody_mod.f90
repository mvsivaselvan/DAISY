module RigidBody_mod

use emxArray

implicit none

type :: RigidBody_t
    !--------------------------------------------------------------
    ! GEOMETRY
    !--------------------------------------------------------------
    ! x0 = reference position of joint 1; (3x1) vector
    ! RJ = rotation of joint 1 coordinate frame with respect to global
    ! rr = end offset, 3x1 vector
    real(kind=8), dimension(3) :: x0, rr
    real(kind=8), dimension(9) :: RJ
    ! MASS PROPERTIES
    !--------------------------------------------------------------
    real(kind=8) :: mass ! mass
    real(kind=8), dimension(9) :: II ! 3x3 mass moment of inertia
    !--------------------------------------------------------------
    ! STIFFNESS AND DAMPING PROPERTIES
    !--------------------------------------------------------------
    real(kind=8), dimension(9) :: KR, KT, CR, CT
    !--------------------------------------------------------------
    ! STATE
    !--------------------------------------------------------------
    real(kind=8), dimension(3) :: d, phi ! displacement and rotation
    real(kind=8), dimension(3) :: ddot, phidot
    real(kind=8), dimension(3) :: dddot, phiddot
    real(kind=8), dimension(6) :: F ! Force
    real(kind=8), dimension(36) :: M, C, K
    real(kind=8), dimension(18) :: B
end type RigidBody_t

contains

!--------------------------------------------------------

subroutine rigidbody_setup(this, x0, RJ, rr,&
                           mass, II, KT, KR, &
                           alph_m, bet_KT, alph_II, bet_KR)

type (RigidBody_t), intent(out) :: this
real(kind=8), dimension(3), intent(in) :: x0, rr
real(kind=8), dimension(9), intent(in) :: RJ
real(kind=8), intent(in) :: mass
real(kind=8) :: alph_m ! mass-proportional damping
real(kind=8) :: bet_KT ! translational stiffness-proportional damping
real(kind=8) :: alph_II ! moment of inertia-proportional damping
real(kind=8) :: bet_KR ! rotational stiffness-proportional damping
real(kind=8), dimension(9), intent(in) :: II, KT, KR

this%x0 = x0
this%RJ = RJ
this%rr = rr

this%mass = mass
this%II = II
this%KT = KT
this%KR = KR

this%CT = bet_KT*KT
this%CT(1) = alph_m*mass
this%CT(5) = alph_m*mass
this%CT(9) = alph_m*mass
this%CR = alph_II*II + bet_KR*KR

end subroutine rigidbody_setup

!--------------------------------------------------------

subroutine rigidbody_destroy(this)

type(RigidBody_t), intent(out) :: this

end subroutine rigidbody_destroy

!--------------------------------------------------------

subroutine rigidbody_setState(this, dyn, u, x, xd, xdd)

use MATLABelements, only : RigidBodyForce

type(RigidBody_t), intent(out) :: this
integer, intent(in) :: dyn
real(kind=8), dimension(3), intent(in) :: u
real(kind=8), dimension(:), intent(in) :: x
real(kind=8), dimension(:), intent(in), optional :: xd, xdd

this%d = x(1:3)
this%phi = x(4:6)

if (present(xd)) then
  this%ddot = xd(1:3)
  this%phidot = xd(4:6)
endif

if (present(xdd)) then
  this%dddot = xdd(1:3)
  this%phiddot = xdd(4:6)
endif

call RigidBodyForce(this%d, this%phi, &
                    this%ddot, this%phidot, this%dddot, this%phiddot,  &
                    this%mass, this%II, this%KT, this%CT, this%KR, this%CR, &
                    this%rr, this%RJ, u, &
                    this%F, this%K, this%C, this%M, this%B)

end subroutine rigidbody_setState

!--------------------------------------------------------

end module RigidBody_mod
    