program main
    
use CableElement

    
implicit none

integer :: N = 10 ! number of control points
integer :: d = 3  ! Degree -- cubic spline
integer :: dbrev = 2 ! Degree for twist
integer :: dbar = 1 ! Degree for strain projection basis
integer :: Ng = 3 ! Number of Gauss points per element

type(Cable) :: cable1

call cable_setup(cable1, N, d, dbrev, dbar, Ng)

print*,'Cable created.'

call cable_destroy(cable1)

print*,'Cable destroyed.'
    
end program main