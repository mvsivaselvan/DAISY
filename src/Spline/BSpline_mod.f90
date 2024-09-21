module BSpline
    
use, intrinsic :: iso_c_binding

implicit none
    
contains

!-------------------------------------------------------------    
    
subroutine getKnotVector(N, k, knots)

! Get knot vector for clamped B-Splines
! INPUTS
! N = number basis functions
! k = order (degree+1)
! knots = knot vector
integer, intent(in) :: N
integer, intent(in) :: k
real(kind=8), dimension(:), intent(out) :: knots

integer :: i

knots(1:k-1) = 0.d0
do i = k,N+1
    knots(i) = (real(i,8)-k)/(N-k+1)
enddo
knots(N+2:N+k) = 1.d0

end subroutine getKnotVector
    
!-------------------------------------------------------------   
    
subroutine spcol(knots, k, colpts, nderiv, colmat)

! Computes values and derivatives of B-Spline basis functions
! at given locations
! INPUTS
! knots = knot vector
! k = order of the B-Spline basis (degree+1, for example 4 for cubic)
! colpts = positions at which to compute basis and derivatives
! nderiv = number of derivatives to compute 
!          (1 for just value, 2 for value and 1st derivative etc.)
! OUTPUTS
! colmat = (nderiv*length(colpts)) x N matrix, 
!          where N = number of basis functions (or control points)

real(kind=8), dimension(:), intent(in) :: knots
integer, intent(in) :: k
real(kind=8), dimension(:), intent(in) :: colpts
integer, intent(in) :: nderiv
real(kind=8), dimension(:,:), intent(out) :: colmat

real(kind=8), dimension(:,:), allocatable :: awork ! work area for bsplvd
                                                   ! see doc for bsplvd 
real(kind=8), dimension(:,:), allocatable :: dbiatx ! temp arrray to get
                                                    ! output from bsplvd

integer :: N ! number of basis functions

integer :: nknots, ncolpts
integer :: left

integer :: i, j, p, q

nknots = size(knots)
N = nknots - k;

ncolpts = size(colpts)

allocate( awork(k,k) )
allocate( dbiatx(k, nderiv) )

colmat = 0.d0
do i = 1,ncolpts
    left = 0 ! index of knot to the left of this colpt
    do j = 1,nknots
        if (colpts(i) >= knots(j)) then ! larget index with knot value < colpts(i)
            left = j
        else
            exit
        endif
    enddo

    ! Compute the value and derivatives of non-zero basis functions
    call bsplvd(knots, k, colpts(i), left, awork, dbiatx, nderiv)
    
    ! Populate these in colmat
    do p = 1, nderiv
        do q = 1, k
            colmat((i-1)*nderiv+p, left-k+q) = dbiatx(q, p)
        enddo
    enddo
    
enddo

deallocate(awork)
deallocate(dbiatx)

end subroutine spcol

!-------------------------------------------------------------

subroutine spcolC(knots, nknots, k, colpts, ncolpts, nderiv, colmat, ncolmat) bind(C, name="spcolC")

integer(c_int), intent(in), value :: nknots, ncolpts, ncolmat
real(kind=c_double), dimension(nknots), intent(in) :: knots
integer(c_int), intent(in), value :: k
real(kind=c_double), dimension(ncolpts), intent(in) :: colpts
integer(c_int), intent(in), value :: nderiv
real(kind=c_double), dimension(ncolmat), intent(out), target :: colmat
real(kind=c_double), dimension(:,:), pointer :: colmat_

integer :: nbasis

nbasis = nknots - k;

call c_f_pointer(c_loc(colmat),colmat_,[ncolpts*nderiv,nbasis])

call spcol(knots, k, colpts, nderiv, colmat_)

end subroutine spcolC

!-------------------------------------------------------------

end module BSpline