program main
    
use CableElement

    
implicit none

integer :: N = 10 ! number of control points
integer :: d = 3  ! Degree -- cubic spline
integer :: dbrev = 2 ! Degree for twist
integer :: dbar = 1 ! Degree for strain projection basis
integer :: Ng = 3 ! Number of Gauss points per element

real(kind=8), dimension(101,3) :: refGeom
integer :: i

type(Cable) :: cable1

! Read reference geometry from file
open(file='refCircle.txt', unit=100, status='old')
do i = 1,101
    read(100,*)refGeom(i,1:3)
enddo
close(unit=100)

! Set up cable
call cable_setup(cable1, N, d, dbrev, dbar, Ng, refGeom)
print*,'Cable created.'

! Write P0 to file
open(file='P0.dat', unit=100, status='unknown')
do i = 1,N
    write(100,*)cable1%P0%data_(:,i)
enddo
close(unit=100)

! Destroy cable
call cable_destroy(cable1)
print*,'Cable destroyed.'
    
end program main