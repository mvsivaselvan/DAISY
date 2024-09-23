program main
    
use CableElement
    
implicit none

real(kind=8), parameter :: pi = 3.14159265358979d0
real(kind=8), parameter :: g = 386.4d0 ! in/s^2

integer :: N = 10 ! number of control points
integer :: d = 3  ! Degree -- cubic spline
integer :: dbrev = 2 ! Degree for twist
integer :: dbar = 1 ! Degree for strain projection basis
integer :: Ng = 3 ! Number of Gauss points per element

real(kind=8), dimension(101,3) :: refGeom ! Coordinates pf reference geometry
integer :: i

! Material and section properties
real(kind=8) :: Acable = 1.57 ! in^2
real(kind=8) :: density = 2.54d-4 ! (lb-s^2/in)/in^3
real(kind=8) :: Ecable = 1e7 ! psi
real(kind=8) :: Icable = 2.18e-3 ! in^4
real(kind=8) :: rho ! mass per unit length
real(kind=8) :: EI 
real(kind=8) :: EA 
real(kind=8) :: GJ 
real(kind=8) :: betBEND = 0.01 ! damping coefficient for bending
real(kind=8) :: betAX = 0.01 ! damping coefficient for axial
real(kind=8) :: betTOR = 0.01 ! damping coefficient for torsion

real(kind=8), dimension(3,3) :: II 

type(Cable) :: cable1

! Read reference geometry from file
open(file='refCircle.txt', unit=100, status='old')
do i = 1,101
    read(100,*)refGeom(i,1:3)
enddo
close(unit=100)

! Section properties
rho = density*Acable
EI = Ecable*Icable
EA = Ecable*Acable
GJ = (2.d0*EI)/(2.d0*(1.d0+0.3d0)) ! I is Imin, i.e., wire-wise, so J = 2*I;
                                   ! and G = E/2/(1+nu), taking n = 0.3
II = reshape((/1.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 2.d0/),shape(II))
II = II*(density*(Acable/pi)**2.d0/64.d0)

! Set up cable
call cable_setup(cable1, N, d, dbrev, dbar, Ng, refGeom, &
                 rho, EI, EA, GJ, betBEND, betAX, betTOR)
print*,'Cable created.'

! Write R0 to file
open(file='Kbar.dat', unit=100, status='unknown')
do i = 1,N-d+dbar
    write(100,'(8e13.5)')cable1%Kbar11%data_(:,i)
enddo
close(unit=100)

! Destroy cable
call cable_destroy(cable1)
print*,'Cable destroyed.'
    
end program main