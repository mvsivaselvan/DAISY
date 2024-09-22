module CableElement

use emxArray

implicit none

type :: Cable
    integer :: N ! Number of control points
    integer :: d ! degree of spline (ex, 3 for cubic) for displacement
    integer :: dbrev ! degree of spline for twist
    integer :: dbar ! degree of spline for strain projection
    integer :: Ng ! number of Gauss pts per element (knot interval)
    
    type(emxArray_1d_wrapper) :: knots, xg, wg
    type(emxArray_2d_wrapper) :: colmat, colmat_brev, colmat_bar
    ! nel(n) = index of element to which quadrature pt xg(n) belongs to
    type(emxArray_1d_wrapper) :: nel
    
    type(emxArray_2d_wrapper) :: P0
end type Cable

contains

!--------------------------------------------------------

subroutine cable_setup(this, N, d, dbrev, dbar, Ng, &
                       refGeom)

use GaussQuad
use BSpline
use MATLABelements, only: SplineApproximation

type(Cable), intent(out) :: this
integer, intent(in) :: N, d, dbrev, dbar, Ng
real(kind=8), dimension(:,:), intent(in) :: refGeom ! (x,y,z) coord of reference geometry

integer :: Nelem
real(kind=8) :: Lelem
real(kind=8), dimension(:), allocatable :: xg, wg
real(kind=8), dimension(:), allocatable :: knots_brev, knots_bar
type(emxArray_2d_wrapper) :: gamm_
type(emxArray_1d_wrapper) :: J_
integer :: i
real(kind=8) :: err

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
allocate(knots_brev(N+dbrev+1))
call getKnotVector(N, dbrev+1, knots_brev)
call emxArray_2d_create(this%colmat_brev, Ng*Nelem*2, N)
call spcol(knots_brev, dbrev+1, this%xg%data_, 2, this%colmat_brev%data_)
deallocate(knots_brev)

! Create collocation matrix for strain projection
allocate(knots_bar(N+dbar+1))
call getKnotVector(N, dbar+1, knots_bar)
call emxArray_2d_create(this%colmat_bar, Ng*Nelem, N)
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

end subroutine cable_setup
    
!--------------------------------------------------------

subroutine cable_destroy(this)

type(Cable), intent(out) :: this

call emxArray_1d_destroy(this%knots)
call emxArray_1d_destroy(this%xg)
call emxArray_1d_destroy(this%wg)
call emxArray_1d_destroy(this%nel)
call emxArray_2d_destroy(this%colmat)
call emxArray_2d_destroy(this%colmat_brev)
call emxArray_2d_destroy(this%colmat_bar)
call emxArray_2d_destroy(this%P0)

end subroutine cable_destroy
                       
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
    
end module CableElement