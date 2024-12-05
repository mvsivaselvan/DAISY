module Cable_mod

use emxArray
use Element_mod, only : Element_t

implicit none

type, extends(Element_t) :: Cable_t
    !--------------------------------------------------------------
    ! GEOMETRY
    !--------------------------------------------------------------
    ! x01 = reference position of joint 1; (3x1) vector
    ! RJ1 = rotation of joint 1 coordinate frame with respect to global
    ! RE1 = rotation of cable end 1 with respect to joint 1 coordinate frame
    ! r1 = position of cable end 1 relative to joint 1 in joint 1 coordinate 
    !       system (end offset, 3x1 vector)
    ! x02, RJ2, RE2, r2 - same at end 2 (s=1) 
    real(kind=8), dimension(3) :: x01, r1, x02, r2
    real(kind=8), dimension(9) :: RJ1, RE1, RJ2, RE2
    
    !--------------------------------------------------------------
    ! SECTION PROPERTIES
    !--------------------------------------------------------------
    real(kind=8) :: rho, EA, EI, GJ, betAX, betBEND, betTOR, alph0
    real(kind=8), dimension(9) :: II

    !--------------------------------------------------------------
    ! SPLINE DISCRETIZATION
    !--------------------------------------------------------------
    integer :: N, Nbrev, Nbar ! Number of control points for displacement, twist and strain proj
    integer :: d ! degree of spline (ex, 3 for cubic) for displacement
    integer :: dbrev ! degree of spline for twist
    integer :: dbar ! degree of spline for strain projection
    integer :: Ng ! number of Gauss pts per element (knot interval)
    
    type(emxArray_1d_wrapper) :: knots, xg, wg
    type(emxArray_2d_wrapper) :: colmat, colmat_brev, colmat_bar
    ! nel(n) = index of element to which quadrature pt xg(n) belongs to
    type(emxArray_1d_wrapper) :: nel
    
    type(emxArray_2d_wrapper) :: P0, R0
    
    type(emxArray_2d_wrapper) :: Mbar, Kbar11, Dbar11
    
    !--------------------------------------------------------------
    ! STATE
    !--------------------------------------------------------------
    real(kind=8), dimension(3) :: d1, d2, phi1, phi2
    real(kind=8) :: gamma1, gamma2
    type(emxArray_2d_wrapper) :: Pmid
    type(emxArray_1d_wrapper) :: varThetamid
    real(kind=8), dimension(3) :: d1dot, d2dot, phi1dot, phi2dot
    real(kind=8) :: gamma1dot, gamma2dot
    type(emxArray_2d_wrapper) :: Pmiddot
    type(emxArray_1d_wrapper) :: varThetamiddot
    real(kind=8), dimension(3) :: d1ddot, d2ddot, phi1ddot, phi2ddot
    real(kind=8) :: gamma1ddot, gamma2ddot
    type(emxArray_2d_wrapper) :: Pmidddot
    type(emxArray_1d_wrapper) :: varThetamidddot
    
    type(emxArray_1d_wrapper) :: Fb
    type(emxArray_2d_wrapper) :: Kb, Cb, Mb, Bb
contains
    procedure :: Cable_t => make_cable
    procedure :: destroy_cable
    procedure :: setState => setState_cable
end type Cable_t

contains

!--------------------------------------------------------

subroutine make_cable(this, &
                      x01, RJ1, RE1, r1, x02, RJ2, RE2, r2, &
                      N, d, dbrev, dbar, Ng, refGeom, &
                      rho, EI, EA, GJ, betBEND, betAX, betTOR, alph0, II)

use GaussQuad
use BSpline
use MATLABelements, only: SplineApproximation, getBishopFrame, CableMbar

class(Cable_t), intent(out) :: this

real(kind=8), dimension(3), intent(in) :: x01, r1, x02, r2
real(kind=8), dimension(9), intent(in) :: RJ1, RE1, RJ2, RE2

integer, intent(in) :: N, d, dbrev, dbar, Ng
real(kind=8), dimension(:,:), intent(in) :: refGeom ! (x,y,z) coord of reference geometry
real(kind=8), intent(in) :: rho, EI, EA, GJ, betBEND, betAX, betTOR, alph0 ! section properties
real(kind=8), dimension(9), intent(in) :: II ! section mass moment of inertia

integer :: Nbrev, Nbar, Nelem
real(kind=8) :: Lelem
real(kind=8), dimension(:), allocatable :: xg, wg
real(kind=8), dimension(:), allocatable :: knots_brev, knots_bar
type(emxArray_2d_wrapper) :: gamm_
type(emxArray_1d_wrapper) :: J_
real(kind=8), dimension(3) :: P0col1
integer :: i
real(kind=8) :: err
integer :: nDof

this%x01 = x01
this%RJ1 = RJ1
this%RE1 = RE1
this%r1 = r1
this%x02 = x02
this%RJ2 = RJ2
this%RE2 = RE2
this%r2 = r2

this%rho = rho
this%EA = EA
this%EI = EI
this%GJ = GJ
this%betAx = betAX
this%betBEND = betBEND
this%betTOR = betTOR
this%alph0 = alph0
this%II = II

this%N = N
this%d = d
this%dbrev = dbrev
this%dbar = dbar
this%Ng = Ng

! Create knot vector
call emxArray_1d_create(this%knots, N+d+1)
call getKnotVector(N, d+1, this%knots%data_)

! Create quadrature points
Nelem = N - d
Lelem = 1.d0/Nelem
allocate(xg(Ng))
allocate(wg(Ng))
call GaussQuadrature(Ng,xg,wg) 
xg = (1.d0+xg)/2.d0*Lelem ! Gauss points \in (0,Lelem)
wg = wg*(Lelem/2.d0) ! weights multiplied by appropriate jacobian 
call emxArray_1d_create(this%xg, Nelem*Ng)
call emxArray_1d_create(this%wg, Nelem*Ng)
call emxArray_1d_create(this%nel, Nelem*Ng)
do i = 1,Nelem
    this%xg%data_((i-1)*Ng+1:i*Ng) = xg + (i-1)*Lelem
    this%wg%data_((i-1)*Ng+1:i*Ng) = wg
    this%nel%data_((i-1)*Ng+1:i*Ng) = i
enddo
deallocate(xg)
deallocate(wg)

! Create collocation matrix for displacement, 
! containing B-Splines and derivatives
call emxArray_2d_create(this%colmat, Ng*Nelem*3, N)
call spcol(this%knots%data_, d+1, this%xg%data_, 3, this%colmat%data_)

! Create collocation matrix for twist
allocate(knots_brev(N-d+1+2*dbrev))
Nbrev = N - d + dbrev
this%Nbrev = Nbrev
call getKnotVector(Nbrev, dbrev+1, knots_brev)
call emxArray_2d_create(this%colmat_brev, Ng*Nelem*2, Nbrev)
call spcol(knots_brev, dbrev+1, this%xg%data_, 2, this%colmat_brev%data_)
deallocate(knots_brev)

! Create collocation matrix for strain projection
allocate(knots_bar(N-d+1+2*dbar))
Nbar = N - d + dbar
this%Nbar = Nbar
call getKnotVector(Nbar, dbar+1, knots_bar)
call emxArray_2d_create(this%colmat_bar, Ng*Nelem, Nbar)
call spcol(knots_bar, dbar+1, this%xg%data_, 1, this%colmat_bar%data_)
deallocate(knots_bar)

! Construct Spline approximation to given data
call emxArray_2d_create(this%P0, 3, N)
call emxArray_2d_create(gamm_, 3, Ng*Nelem+2)
call emxArray_1d_create(J_, Ng*Nelem+2)
call PiecewiseLinearCurve(refGeom, this%xg%data_, gamm_%data_, J_%data_)
call SplineApproximation(gamm_%emx, J_%emx, real(N,8), this%xg%emx, this%wg%emx, &
                         this%colmat%emx, this%P0%emx, err)
call emxArray_2d_destroy(gamm_)
call emxArray_1d_destroy(J_)

! Move start point of curve to origin
P0col1 = this%P0%data_(:,1)
do i = 1, N
    this%P0%data_(:,i) = this%P0%data_(:,i) - P0col1
enddo

! Construct Bishop frame
call emxArray_2d_create(this%R0, 3, 3*Ng*Nelem+6) ! a rotation matrix per Gauss pt + ends
call getBishopFrame(this%P0%emx, this%knots%emx, real(d,8), this%xg%emx, this%R0%emx)

! Compute strain projection related matrices
call emxArray_2d_create(this%Mbar, Nbar, Nbar)
call emxArray_2d_create(this%Kbar11, Nbar, Nbar)
call emxArray_2d_create(this%Dbar11, Nbar, Nbar)
call CableMbar(this%P0%emx, EA, betAX, this%colmat%emx, this%colmat_bar%emx, this%wg%emx, &
               this%Mbar%emx, this%Kbar11%emx, this%Dbar11%emx)

! allocate and initialize state variables
call emxArray_2d_create(this%Pmid, 3, N-4)
call emxArray_1d_create(this%varThetamid, Nbrev-2)
call emxArray_2d_create(this%Pmiddot, 3, N-4)
call emxArray_1d_create(this%varThetamiddot, Nbrev-2)
call emxArray_2d_create(this%Pmidddot, 3, N-4)
call emxArray_1d_create(this%varThetamidddot, Nbrev-2)
this%d1 = 0.d0
this%phi1 = 0.d0 ! NOT REALLY
this%gamma1 = 0.d0
this%d2 = 0.d0
this%phi2 = 0.d0 ! NOT REALLY
this%gamma2 = 0.d0
this%Pmid%data_ = 0.d0
this%varThetamid%data_ = 0.d0
this%d1dot = 0.d0
this%phi1dot = 0.d0
this%gamma1dot = 0.d0
this%d2dot = 0.d0
this%phi2dot = 0.d0
this%gamma2dot = 0.d0
this%Pmiddot%data_ = 0.d0
this%varThetamiddot%data_ = 0.d0
this%d1ddot = 0.d0
this%phi1ddot = 0.d0
this%gamma1ddot = 0.d0
this%d2ddot = 0.d0
this%phi2ddot = 0.d0
this%gamma2ddot = 0.d0
this%Pmidddot%data_ = 0.d0
this%varThetamidddot%data_ = 0.d0

nDof = 3*N + Nbrev
call emxArray_1d_create(this%Fb, nDof)
call emxArray_2d_create(this%Kb, nDof, nDof)
call emxArray_2d_create(this%Cb, nDof, nDof)
call emxArray_2d_create(this%Mb, nDof, nDof)
call emxArray_2d_create(this%Bb, nDof, 3)

end subroutine make_cable
    
!--------------------------------------------------------

subroutine destroy_cable(this)

class(Cable_t), intent(out) :: this

call emxArray_1d_destroy(this%knots)
call emxArray_1d_destroy(this%xg)
call emxArray_1d_destroy(this%wg)
call emxArray_1d_destroy(this%nel)
call emxArray_2d_destroy(this%colmat)
call emxArray_2d_destroy(this%colmat_brev)
call emxArray_2d_destroy(this%colmat_bar)
call emxArray_2d_destroy(this%P0)
call emxArray_2d_destroy(this%R0)
call emxArray_2d_destroy(this%Mbar)
call emxArray_2d_destroy(this%Kbar11)
call emxArray_2d_destroy(this%Dbar11)
call emxArray_2d_destroy(this%Pmid)
call emxArray_1d_destroy(this%varThetamid)
call emxArray_2d_destroy(this%Pmiddot)
call emxArray_1d_destroy(this%varThetamiddot)
call emxArray_2d_destroy(this%Pmidddot)
call emxArray_1d_destroy(this%varThetamidddot)
call emxArray_1d_destroy(this%Fb)
call emxArray_2d_destroy(this%Kb)
call emxArray_2d_destroy(this%Cb)
call emxArray_2d_destroy(this%Mb)
call emxArray_2d_destroy(this%Bb)

end subroutine destroy_cable
                       
!--------------------------------------------------------

subroutine setState_cable(this, dyn, u, x, xd, xdd)

use MATLABelements, only : CableForceRotBCinCoord

class(Cable_t), intent(inout) :: this
integer, intent(in) :: dyn
real(kind=8), dimension(3), intent(in) :: u
real(kind=8), dimension(:), intent(in) :: x
real(kind=8), dimension(:), intent(in), optional :: xd, xdd

this%d1 = x(1:3)
this%phi1 = x(4:6)
this%gamma1 = x(7)
this%d2 = x(8:10)
this%phi2 = x(11:13)
this%gamma2 = x(14)
this%Pmid%data_ = reshape(x(15:(14+(3*(this%N-4)))),[3,(3*(this%N-4))])
this%varThetamid%data_ = x((15+(3*(this%N-4))):(14+(3*(this%N-4))+this%Nbrev-2))

if (present(xd)) then
  this%d1dot = xd(1:3)
  this%phi1dot = xd(4:6)
  this%gamma1dot = xd(7)
  this%d2dot = xd(8:10)
  this%phi2dot = xd(11:13)
  this%gamma2dot = xd(14)
  this%Pmiddot%data_ = reshape(xd(15:(14+(3*(this%N-4)))),[3,(3*(this%N-4))])
  this%varThetamiddot%data_ = xd((15+(3*(this%N-4))):(14+(3*(this%N-4))+this%Nbrev-2))
endif

if (present(xdd)) then
  this%d1ddot = xdd(1:3)
  this%phi1ddot = xdd(4:6)
  this%gamma1ddot = xdd(7)
  this%d2ddot = xdd(8:10)
  this%phi2ddot = xdd(11:13)
  this%gamma2ddot = xdd(14)
  this%Pmidddot%data_ = reshape(xdd(15:(14+(3*(this%N-4)))),[3,(3*(this%N-4))])
  this%varThetamidddot%data_ = xdd((15+(3*(this%N-4))):(14+(3*(this%N-4))+this%Nbrev-2))
endif

call CableForceRotBCinCoord &
        (this%d1, this%phi1, this%gamma1, this%d2, this%phi2, this%gamma2, &
         this%Pmid%emx, this%varThetamid%emx, this%P0%emx, &
         this%d1dot, this%phi1dot, this%gamma1dot, &
         this%d2dot, this%phi2dot, this%gamma2dot, &
         this%Pmiddot%emx, this%varThetamiddot%emx, &
         this%d1ddot, this%phi1ddot, this%gamma1ddot,&
         this%d2ddot, this%phi2ddot, this%gamma2ddot, &
         this%Pmidddot%emx, this%varThetamidddot%emx, &
         this%x01, this%RJ1, this%RE1, this%r1, this%x02, this%RJ2, this%RE2, this%r2, &
         this%R0%emx, this%II, &
         this%rho, this%EA, this%EI, this%GJ, this%betAX, this%betBEND, this%betTOR, &
         this%xg%emx, this%wg%emx, this%nel%emx, &
         this%colmat%emx, this%colmat_brev%emx, this%colmat_bar%emx, &
         real(this%d,8), real(this%dbrev,8), real(this%dbar,8), &
         this%Mbar%emx, u, this%Kbar11%emx, this%Dbar11%emx, real(dyn,8), this%alph0, &
         this%Fb%emx, this%Kb%emx, this%Cb%emx, this%Mb%emx, this%Bb%emx)

end subroutine setState_cable

!--------------------------------------------------------

subroutine cable_getShape(this, npts, shpe)

use BSpline
use MATLABelements, only : getCableDeformedShape

type(Cable_t), intent(in) :: this
integer, intent(in) :: npts ! number of points in shaps
real(kind=8), dimension(:,:), intent(out) :: shpe

real(kind=8), dimension(:), allocatable :: s
real(kind=8), dimension(:,:), allocatable :: colmat

real(kind=8), dimension(9) :: Rb10, Rb20
type(emxArray_2d_wrapper) :: P

integer :: Nelem, Rlastind
integer :: i, j

Nelem = this%N - this%d
Rlastind = 3*(1+Nelem*this%Ng)
do i = 1,3
    do j = 1,3
        Rb10((j-1)*3+i) = this%R0%data_(i,j)
        Rb20((j-1)*3+i) = this%R0%data_(i,Rlastind+j)
    enddo
enddo

call emxArray_2d_create(P, 3, this%N)
call getCableDeformedShape(this%d1, this%phi1, this%gamma1, &
                           this%d2, this%phi2, this%gamma2, &
                           this%Pmid%emx, &
                           this%x01, this%RJ1, this%RE1, this%r1, Rb10, &
                           this%x02, this%RJ2, this%RE2, this%r2, Rb20, &
                           P%emx)

allocate(s(npts))
allocate(colmat(npts,this%N))

s(1) = 0.d0
do i = 2,npts
    s(i) = s(i-1) + 1.d0/(npts-1)
enddo

call spcol(this%knots%data_, this%d+1, s, 1, colmat)

shpe = matmul(colmat,transpose(P%data_)) ! USE BLAS FOR THIS

deallocate(s)
deallocate(colmat)
call emxArray_2d_destroy(P)

end subroutine cable_getShape

!--------------------------------------------------------

subroutine PiecewiseLinearCurve(coord, s, x, jac)

! Piecewise linear curve connecting given points coord
! The parametrization s \in [0,1] \mapsto x is such that 
! x(0) = the first point in coord
! x(s) is such that the length of the curve from the first point 
!      to x(s) is equal to s*(total length)
! This implies x(1) = last point in coord
! J(s) is the jacobina at x(s), which is simply L!

! Length of piecewise linear curve
! MATLAB CODE
! dL = sqrt(diff(coord(:,1)).^2 + diff(coord(:,2)).^2 + diff(coord(:,3)).^2);
! L = [0; cumsum(dL)]; % so that L(n) = distance from the first point 
!                      % to the nth point
! x = interp1(L, coord, s*L(end))'; % transpose, so that it is columnwise
! J = ones(1,size(x,2))*L(end);
! END MATLAB CODE

real(kind=8), dimension(:,:), intent(in) :: coord
real(kind=8), dimension(:), intent(in) :: s
real(kind=8), dimension(:,:), intent(out) :: x
real(kind=8), dimension(:), intent(out) :: jac

integer :: n, ns
real(kind=8), dimension(:), allocatable :: L
real(kind=8) :: dL, Ls
integer :: i, j, k, ind

n = size(coord,1) ! number of points in the reference geometry
ns = size(s) ! number of parameter values (inside (0,1))

! cummulative length of line segments
allocate(L(n))
L = 0.d0 
do i = 1,n-1
    dL = (coord(i+1,1)-coord(i,1))**2.d0 &
        +(coord(i+1,2)-coord(i,2))**2.d0 &
        +(coord(i+1,3)-coord(i,3))**2.d0
    dL = sqrt(dL)
    L(i+1) = L(i) + dL
enddo

x(:,1) = coord(1,:)
do j = 1,ns
    Ls = s(j)*L(n)
    do i = 1,n
        if (L(i) > Ls) exit
    enddo
    ind = i-1 ! left point for interpolation
    do k = 1,3
        x(k,j+1) = coord(ind,k) + &
            (coord(ind+1,k)-coord(ind,k))*((Ls-L(ind))/(L(ind+1)-L(ind)))
    enddo
enddo
x(:,ns+2) = coord(n,:)

jac = L(n)

deallocate(L)

end subroutine PiecewiseLinearCurve

!--------------------------------------------------------
    
end module Cable_mod