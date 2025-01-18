module RotationUtils_mod

implicit none
    
contains
    
!------------------------------------------------
    
subroutine invRodrigues(R, phi)

use lapack95

real(kind=8), dimension(3,3), intent(in) :: R
real(kind=8), dimension(3), intent(out) :: phi

real(kind=8), dimension(2,2) :: A = 0.d0
real(kind=8), dimension(3) :: x = 0.d0
integer(kind=4), dimension(2) :: ipiv
real(kind=8) :: traceR, thet, tracexhatR

! based on matlab function invRodrigues

traceR = R(1,1)+R(2,2)+R(3,3)

if (dabs(traceR-3.d0) < 1.d-10) then ! R is the identity
    phi = 0.d0
    return
endif

! find eigenvector of R corresponding to eigenvalue 1
A(1,1) = 1.d0 - R(1,1)
A(1,2) = -R(1,2)
A(2,1) = -R(2,1)
A(2,2) = 1.d0 - R(2,2)
x(1:2) = R(1:2,3)
call gesv(A, x, ipiv)
x(3) = 1.d0

x = x/dsqrt(x(1)**2.d0 + x(2)**2.d0 + x(3)**2.d0)

tracexhatR = x(1)*(R(2,3)-R(3,2)) + x(2)*(R(3,1)-R(1,3)) &
           + x(3)*(R(1,2)-R(2,1))

thet = datan2(-tracexhatR,traceR-1.d0)
phi = thet*x;

end subroutine invRodrigues

!------------------------------------------------

end module RotationUtils_mod