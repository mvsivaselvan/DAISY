program main
    
use blas95
use lapack95
use flib_dom
use Cable_mod
    
implicit none

real(kind=8), parameter :: pi = 3.14159265358979d0
real(kind=8), parameter :: g = 386.4d0 ! in/s^2

character(80) :: inpFilename
integer :: cmdLineStatus

integer :: N ! number of control points
integer :: d ! Degree -- cubic spline
integer :: dbrev ! Degree for twist
integer :: dbar ! Degree for strain projection basis
integer :: Ng ! Number of Gauss points per element

real(kind=8), dimension(101,3) :: refGeom ! Coordinates pf reference geometry
integer :: i

! Material and section properties
real(kind=8) :: rho ! mass per unit length
real(kind=8) :: EI 
real(kind=8) :: EA 
real(kind=8) :: GJ 
real(kind=8) :: betBEND ! damping coefficient for bending
real(kind=8) :: betAX ! damping coefficient for axial
real(kind=8) :: betTOR ! damping coefficient for torsion
real(kind=8) :: alph0 ! mass proportional damping factor

real(kind=8), dimension(3,3) :: II 

real(kind=c_double), dimension(3) :: x01, r1, x02, r2
real(kind=c_double), dimension(9) :: RJ1, RE1, RJ2, RE2

class(Cable_t), pointer :: cable1

real(kind=8), dimension(:), allocatable :: x ! displacement
real(kind=8), dimension(:), allocatable :: x_ ! displacement during line search

! input (gravity+ground acceleration
real(kind=8), dimension(3) :: u = (/0.d0, 0.d0, -386.4d0/) 

! line search parameters (see page 33 of Nocedal and Wright)
real(kind=8), parameter :: armijo_gamm = 0.9d0
real(kind=8), parameter :: armijo_c1 = 1.d-2
real(kind=8) :: armijo_alpha
integer :: k ! Newton iteration index
integer :: l ! Line search step index
integer :: nDof
real(kind=8), dimension(:), allocatable :: gk ! residual force
real(kind=8) :: normgk
real(kind=8), dimension(:), allocatable :: Dx ! Newton direction
real(kind=8), dimension(:), allocatable :: gk_ ! residual force during line search
real(kind=8), dimension(:,:), allocatable :: Dg ! stiffness matrix
integer, dimension(:), allocatable :: ipiv ! pivot for linear solve

integer, parameter :: ndefpts = 101
real(kind=8), dimension(ndefpts,3) :: shpe

type(fNode), pointer :: doc
type(fNodeList), pointer :: domNodes
type(fNode), pointer :: cableDOMnode
type(fNode), pointer :: refgeomDOMnode
type(fNode), pointer :: propertiesDOMnode
type(fNode), pointer :: sectionInertiaDOMnode
type(fNode), pointer :: splineparamDOMnode
type(fNode), pointer :: stepDOMnode
type(fNode), pointer :: endpositionDOMnode
character(20) :: refgeomFile
character(20) :: attribString
character(20) :: defgeomFile

real(kind=8) :: endpositionX

call getarg(1, inpFilename, cmdLineStatus)

! doc => parsefile('Circle2Catenary.xml')
doc => parsefile(inpFilename)
domNodes => getElementsByTagName(doc,'Cable')
cableDOMnode => item(domNodes,0)
domNodes => getElementsByTagName(cableDOMnode,'ReferenceGeometry')
refgeomDOMnode => item(domNodes,0)
refgeomFile = getAttribute(refgeomDOMnode,'file')
domNodes => getElementsByTagName(cableDOMnode,'Properties')
propertiesDOMnode => item(domNodes,0)
attribString = getAttribute(propertiesDOMnode,'rho')
rho = dnum(attribString)
attribString = getAttribute(propertiesDOMnode,'EA')
EA = dnum(attribString)
attribString = getAttribute(propertiesDOMnode,'EI')
EI = dnum(attribString)
attribString = getAttribute(propertiesDOMnode,'GJ')
GJ = dnum(attribString)
domNodes => getElementsByTagName(propertiesDOMnode,'SectionMassMomentOfInertia')
sectionInertiaDOMnode => item(domNodes,0)
attribString = getAttribute(sectionInertiaDOMnode,'I11')
II(1,1) = dnum(attribString)
attribString = getAttribute(sectionInertiaDOMnode,'I12')
II(1,2) = dnum(attribString)
attribString = getAttribute(sectionInertiaDOMnode,'I13')
II(1,3) = dnum(attribString)
attribString = getAttribute(sectionInertiaDOMnode,'I22')
II(2,2) = dnum(attribString)
attribString = getAttribute(sectionInertiaDOMnode,'I23')
II(2,3) = dnum(attribString)
attribString = getAttribute(sectionInertiaDOMnode,'I33')
II(3,3) = dnum(attribString)
II(2,1) = II(1,2)
II(3,1) = II(1,3)
II(3,2) = II(2,3)
domNodes => getElementsByTagName(cableDOMnode,'SplineParameters')
splineparamDOMnode => item(domNodes,0)
attribString = getAttribute(splineparamDOMnode,'N')
N = jnum(attribString)
attribString = getAttribute(splineparamDOMnode,'d')
d = jnum(attribString)
attribString = getAttribute(splineparamDOMnode,'dtwist')
dbrev = jnum(attribString)
attribString = getAttribute(splineparamDOMnode,'dstrain')
dbar = jnum(attribString)
attribString = getAttribute(splineparamDOMnode,'NGauss')
Ng = jnum(attribString)

domNodes => getElementsByTagName(doc,'Step')
stepDOMnode => item(domNodes,0)
defgeomFile = getAttribute(stepDOMnode, 'OutputShapeFile')
domNodes => getElementsByTagName(stepDOMnode,'EndPosition')
endpositionDOMnode => item(domNodes,0)
attribString = getAttribute(endpositionDOMnode,'X')
endpositionX = dnum(attribString)

call destroyNode(doc)

! Read reference geometry from file
open(file=trim(refgeomFile), unit=100, status='old')
do i = 1,101
    read(100,*)refGeom(i,1:3)
enddo
close(unit=100)

x01 = (/0.d0, 0.d0, 0.d0/)
RJ1 = (/1.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 1.d0/)
RE1 = (/1.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 1.d0/)
r1 = (/0.d0, 0.d0, 0.d0/)
x02 = (/121.276066636024d0, 0.d0, 0.d0/)
RJ2 = (/1.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 1.d0/)
RE2 = (/1.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 1.d0/)
r2 = (/0.d0, 0.d0, 0.d0/)

! Set up cable
allocate(cable1)
call make_cable(cable1, x01, RJ1, RE1, r1, x02, RJ2, RE2, r2, &
                 N, d, dbrev, dbar, Ng, refGeom, &
                 rho, EI, EA, GJ, betBEND, betAX, betTOR, alph0, II)
print*,'Cable created.'

! A trial displacement
nDof = 3*cable1%N + cable1%Nbrev - 6 ! since the ends are fixed
allocate(x(nDof+6))
allocate(x_(nDof+6))
x(1:3) = 0.d0
x(4:6) = 0.d0
x(7) = sqrt((cable1%P0%data_(1,2)-cable1%P0%data_(1,1))**2.d0 &
          + (cable1%P0%data_(2,2)-cable1%P0%data_(2,1))**2.d0 &
          + (cable1%P0%data_(3,2)-cable1%P0%data_(3,1))**2.d0)
x(8) = endpositionX - cable1%P0%data_(1,cable1%N)
x(9:10) = 0.d0
x(11:13) = 0.d0
x(14) = -sqrt((cable1%P0%data_(1,N)-cable1%P0%data_(1,N-1))**2.d0 &
           + (cable1%P0%data_(2,N)-cable1%P0%data_(2,N-1))**2.d0 &
           + (cable1%P0%data_(3,N)-cable1%P0%data_(3,N-1))**2.d0)
x(15:(14+3*(N-4))) = reshape(cable1%P0%data_(1:3,3:N-2),[3*(N-4)])
x((15+3*(N-4)):(14+3*(N-4)+cable1%Nbrev-2)) = 0.d0

!----------------------------------------------------------------
! Newton's method
!----------------------------------------------------------------
allocate(gk(nDof))
allocate(Dx(nDof))
allocate(gk_(nDof))
allocate(Dg(nDof,nDof))
allocate(ipiv(nDof))
do k = 1, 200
    call cable1%setState(0, u, x)
    
    gk = [cable1%Fb%data_(7), cable1%Fb%data_(14:nDof+6)]
    normgk = nrm2(gk)
    if (normgk <= 1e-8) exit
    
    Dg = cable1%Kb%data_([7,14:nDof+6],[7,14:nDof+6])
    Dx = -gk
    call gesv(Dg, Dx, ipiv) ! this is the Newton direction
    
    armijo_alpha = 1.d0
    do l = 1,50
        x_(1:6) = x(1:6)
        x_(8:13) = x(8:13)
        x_([7,14:nDof+6]) = x([7,14:nDof+6]) + armijo_alpha*Dx
        call cable1%setState(0, u, x_)
        gk_ = [cable1%Fb%data_(7), cable1%Fb%data_(14:nDof+6)]
        if (0.5d0*nrm2(gk_)**2.d0 <= (0.5d0 - armijo_alpha*armijo_c1)*normgk**2.d0) exit
        armijo_alpha = armijo_alpha*armijo_gamm;
    enddo

    x = x_
    write(*,'(a4,i3,a9,e13.5,a10,e13.5,a16,e13.5)')"k = ",k, &
                 ",||g|| = ",normgk, &
                 ",||Dx|| = ",nrm2(Dx), &
                 ",armijo_alpha = ",armijo_alpha     
enddo
deallocate(gk)
deallocate(gk_)
deallocate(Dx)
deallocate(Dg)
deallocate(ipiv)

! Cable deformed shape
call cable_getShape(cable1, ndefpts, shpe)
open(file=defgeomFile, unit=100, status='unknown')
do i = 1,ndefpts
    write(100,'(3e13.5)')shpe(i,:)
enddo
close(unit=100)

! Destroy cable
call destroy_cable(cable1)
deallocate(cable1)
deallocate(x)
deallocate(x_)
print*,'Cable destroyed.'

!! Write Fb and Kb to files
!open(file='Fb.dat', unit=100, status='unknown')
!do i = 1,3*N+cable1%Nbrev
!    write(100,'(1e13.5)')cable1%Fb%data_(i)
!enddo
!close(unit=100)
!
!open(file='Kb.dat', unit=100, status='unknown')
!do i = 1,3*N+cable1%Nbrev
!    write(100,'(39e13.5)')cable1%Kb%data_(i,:)
!enddo
!close(unit=100)

end program main