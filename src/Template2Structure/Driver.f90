program main
    
use GaussQuad
    
implicit none

integer :: N = 10 ! number of control points
integer :: d = 3  ! Degree -- cubic spline
integer :: dbrev = 2 ! Degree for twist
integer :: dbar = 1 ! Degree for strain projection basis
integer :: Ng = 3 ! Number of Gauss points per element

real(kind=8), dimension(3) :: xg, wg

call GaussQuadrature(Ng, xg, wg)

print*,xg
print*,wg

print*,'Hello world Template 2'
    
end program main