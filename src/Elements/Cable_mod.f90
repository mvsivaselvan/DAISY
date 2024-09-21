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
end type Cable

contains

!--------------------------------------------------------

subroutine cable_setup(this, N, d, dbrev, dbar, Ng)

use GaussQuad
use BSpline

type(Cable), intent(out) :: this
integer, intent(in) :: N, d, dbrev, dbar, Ng
!real(kind=8), dimension(:,:), intent(in) :: refGeom ! (x,y,z) coord of reference geometry

integer :: Nelem
real(kind=8) :: Lelem
real(kind=8), dimension(:), allocatable :: xg, wg
real(kind=8), dimension(:), allocatable :: knots_brev, knots_bar
integer :: i

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

end subroutine cable_destroy
                       
!--------------------------------------------------------
    
end module CableElement