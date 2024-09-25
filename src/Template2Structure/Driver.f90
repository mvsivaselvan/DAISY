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
real(kind=8) :: alph0 = 0.d0 ! mass proportional damping factor

real(kind=8), dimension(3,3) :: II 

real(kind=c_double), dimension(3) :: x01, r1, x02, r2
real(kind=c_double), dimension(9) :: RJ1, RE1, RJ2, RE2

type(Cable) :: cable1

real(kind=8), dimension(:), allocatable :: x ! displacement
real(kind=8), dimension(3) :: u ! input (gravity+ground acceleration

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

x01 = (/0.d0, 0.d0, 0.d0/)
RJ1 = (/1.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 1.d0/)
RE1 = (/1.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 1.d0/)
r1 = (/0.d0, 0.d0, 0.d0/)
x02 = (/121.276066636024d0, 0.d0, 0.d0/)
RJ2 = (/1.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 1.d0/)
RE2 = (/1.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 1.d0/)
r2 = (/0.d0, 0.d0, 0.d0/)

! Set up cable
call cable_setup(cable1, x01, RJ1, RE1, r1, x02, RJ2, RE2, r2, &
                 N, d, dbrev, dbar, Ng, refGeom, &
                 rho, EI, EA, GJ, betBEND, betAX, betTOR, alph0, II)
print*,'Cable created.'

! A trial displacement
allocate(x(3*cable1%N+cable1%Nbrev))
x(1:3) = 0.d0
x(4:6) = 0.d0
x(7) = sqrt((cable1%P0%data_(1,2)-cable1%P0%data_(1,1))**2.d0 &
          + (cable1%P0%data_(2,2)-cable1%P0%data_(2,1))**2.d0 &
          + (cable1%P0%data_(3,2)-cable1%P0%data_(3,1))**2.d0)
x(1:3) = 0.d0
x(4:6) = 0.d0
x(7) = sqrt((cable1%P0%data_(1,2)-cable1%P0%data_(1,1))**2.d0 &
          + (cable1%P0%data_(2,2)-cable1%P0%data_(2,1))**2.d0 &
          + (cable1%P0%data_(3,2)-cable1%P0%data_(3,1))**2.d0)
x(8) = 125.d0 - cable1%P0%data_(1,cable1%N)
x(9:10) = 0.d0
x(11:13) = 0.d0
X(14) = -sqrt((cable1%P0%data_(1,N)-cable1%P0%data_(1,N-1))**2.d0 &
           + (cable1%P0%data_(2,N)-cable1%P0%data_(2,N-1))**2.d0 &
           + (cable1%P0%data_(3,N)-cable1%P0%data_(3,N-1))**2.d0)
x(15:(14+3*(N-4))) = reshape(cable1%P0%data_(1:3,3:N-2),[3*(N-4)])
x((15+3*(N-4)):(14+3*(N-4)+cable1%Nbrev-2)) = 0.d0

u = (/0.d0, 0.d0, -386.4d0/)

call cable_setState(cable1, 0, u, x)
print*,'Done setting state.'

! Write Fb and Kb to files
open(file='Fb.dat', unit=100, status='unknown')
do i = 1,3*N+cable1%Nbrev
    write(100,'(1e13.5)')cable1%Fb%data_(i)
enddo
close(unit=100)

open(file='Kb.dat', unit=100, status='unknown')
do i = 1,3*N+cable1%Nbrev
    write(100,'(39e13.5)')cable1%Kb%data_(i,:)
enddo
close(unit=100)

! Destroy cable
call cable_destroy(cable1)
deallocate(x)
print*,'Cable destroyed.'
    
end program main