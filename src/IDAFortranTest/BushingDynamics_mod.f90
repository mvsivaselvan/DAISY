module BushingDynamics

use, intrinsic :: iso_c_binding
use fsundials_core_mod
use fida_mod
use fnvector_serial_mod
    
implicit none

real(kind=8), parameter :: g = 386.4d0

type :: Bushing
    real(kind=8) :: m  ! bushing mass
    real(kind=8) :: II ! bushing mass moment of inertia about pivot point
    real(kind=8) :: kv ! vertical stiffness of cover plate
    real(kind=8) :: kr ! rotational stiffness of cover plate
    real(kind=8) :: cv ! damping coefficient for vertial motion
    real(kind=8) :: cr ! damping coefficient for rotational motion
    real(kind=8) :: h ! height of bushing center of mass from pivot point
    real(kind=8) :: theta0 ! misalignment from vertical (radians)
    real(kind=8), dimension(:), allocatable :: ux, uz ! ground motion components
    real(kind=8) :: dtu ! sampling time for ux and uz
end type Bushing

contains

!----------------------------------------------------------

subroutine bushing_create(this, m , II, kv, kr, cv, cr, h, theta0, ufile, dtu, &
                          uxscale, uzscale)

type(Bushing), intent(out) :: this
real(kind=8), intent(in) :: m, II, kv, kr, cv, cr, h, theta0, dtu
character(*), intent(in) :: ufile
real(kind=8), intent(in) :: uxscale, uzscale ! scale factor for input

integer, parameter :: unit_num = 100
integer :: numpoints ! number of points in the acceleration file
real(kind=8), dimension(3) :: udata
integer :: i

this%m = m
this%II = II
this%kv = kv
this%kr = kr
this%cv = cv
this%cr = cr
this%h = h
this%theta0 = theta0
this%dtu = dtu

open (unit=unit_num, file=ufile, status='old',action='read')
do while (.not. eof(unit_num) )
	numpoints = numpoints + 1
	read(unit_num,*)udata
enddo
allocate(this%ux(numpoints))
allocate(this%uz(numpoints))
rewind(unit_num)
do i=1,numpoints
	read(unit_num,*)udata
	this%ux(i) = udata(1) * uxscale
    this%uz(i) = udata(3) * uzscale
enddo
close(unit=unit_num)

end subroutine bushing_create

!----------------------------------------------------------

subroutine bushing_destroy(this)

type(Bushing), intent(out) :: this

if (allocated(this%ux)) deallocate(this%ux)
if (allocated(this%uz)) deallocate(this%uz)

end subroutine bushing_destroy

!----------------------------------------------------------

subroutine bushing_getInitialDisplacement(this, y)

type(Bushing), intent(in) :: this
real(kind=8), dimension(2), intent(out) :: y

real(kind=8) :: f, fp
integer :: niter

y(1) = -this%m*g/this%kv

! to get y(2), solve
! f(thet) = kr*thet_-m*g*h*sin(thet_+thet0) = 0
! Newton's method
y(2) = 0.d0 ! initial guess
f = this%kr*y(2)-this%m*g*this%h*sin(y(2)+this%theta0)
niter = 0
do while (abs(f) > 1e-10)
    fp = this%kr-this%m*g*this%h*cos(y(2)+this%theta0)
    y(2) = y(2) - f/fp
    f = this%kr*y(2)-this%m*g*this%h*sin(y(2)+this%theta0)
    niter = niter + 1
    if (niter > 20) then
        print*,'Newton for initial displcement did not converge'
        exit
    endif
enddo

end subroutine bushing_getInitialDisplacement

!----------------------------------------------------------

integer(c_int) function bushing_res(t, sunvec_y, sunvec_yp, sunvec_r, user_data) &
    result(ierr) bind(C, name='bushing_res')

real(c_double), value, intent(in) :: t
! although y and yp are really intent in, they cannot be declared as such because 
! they is argument to the FN_VGetArrayPointer call
type(N_Vector), intent(inout) :: sunvec_y 
type(N_Vector), intent(inout) :: sunvec_yp
type(N_Vector), intent(out) :: sunvec_r
type(c_ptr), value, intent(in) :: user_data

type(Bushing), pointer :: this

real(kind=8), dimension(:), pointer :: y, yp, r

real(kind=8), dimension(2,2) :: Mmat, Cmat, Kmat
real(kind=8), dimension(2) :: F

y => FN_VGetArrayPointer(sunvec_y)
yp => FN_VGetArrayPointer(sunvec_yp)
r => FN_VGetArrayPointer(sunvec_r)
call c_f_pointer(user_data, this)

call bushing_ForceAndTangent(this, t, y, yp, F, Mmat, Cmat, Kmat)

! Residual
r(1:2) = yp(1:2)- y(3:4)
r(3:4) = F

ierr = 0

end function bushing_res

!----------------------------------------------------------

subroutine bushing_MatrixEmbeddedLS(ida_mem, ctx, LS) &
    bind(C, name='bushing_MatrixEmbeddedLS')

type(c_ptr), intent(in) :: ida_mem
type(c_ptr), intent(in) :: ctx
type(SUNLinearSolver), pointer, intent(out) :: LS

type(SUNLinearSolver_Ops), pointer :: ops

LS => FSUNLinSolNewEmpty(ctx)

call c_f_pointer(LS%ops, ops)
ops%gettype = c_funloc(bushing_MatrixEmbeddedLSType)
ops%solve = c_funloc(bushing_MatrixEmbeddedLSSolve)
ops%free = c_funloc(bushing_MatrixEmbeddedLSFree)

LS%content = ida_mem

end subroutine bushing_MatrixEmbeddedLS
    
!----------------------------------------------------------

integer(SUNLinearSolver_Type) function bushing_MatrixEmbeddedLSType(LS) &
    result(lstype) bind(C, name='bushing_MatrixEmbeddedLSType')

type(SUNLinearSolver), intent(in) :: LS

lstype = SUNLINEARSOLVER_MATRIX_EMBEDDED

end function bushing_MatrixEmbeddedLSType

!----------------------------------------------------------

integer(c_int) function bushing_MatrixEmbeddedLSSolve(LS, sunmat_A, sunvec_x, sunvec_b, tol) &
    result(ierr) bind(C, name='bushing_MatrixEmbeddedLSSolve')

use lapack95

type(SUNLinearSolver), intent(in) :: LS
type(SUNMatrix), intent(in) :: sunmat_A ! This is a SUNMatrix object that is not used
type(N_Vector), intent(inout) :: sunvec_x ! This is actually intent out
type(N_Vector), intent(inout) :: sunvec_b ! This is actually intent in
real(c_double), value, intent(in) :: tol

real(kind=8), dimension(:), pointer :: x, b

integer, parameter :: neq = 4
real(kind=8), dimension(neq), target :: ypred, yppred, yn, ypn, res
real(kind=8), dimension(1) :: t, cj
type(c_ptr) :: user_data
type(Bushing), pointer :: this
integer :: retval

real(kind=8), dimension(2) :: F
real(kind=8), dimension(2,2) :: Mmat, Cmat, Kmat, Kbar
real(kind=8), dimension(2) :: rhs
integer, dimension(2) :: ipiv

x => FN_VGetArrayPointer_Serial(sunvec_x)
b => FN_VGetArrayPointer_Serial(sunvec_b)

retval = FIDAGetNonlinearSystemData(LS%content, t, &
                                    c_loc(ypred), c_loc(yppred), &
                                    c_loc(yn), c_loc(ypn), &
                                    c_loc(res), cj, user_data)
call c_f_pointer(user_data, this)

call bushing_ForceAndTangent(this, t(1), ypred, yppred, F, Mmat, Cmat, Kmat)

rhs(1) = b(3) + (Cmat(1,1)+cj(1)*Mmat(1,1))*b(1) + (Cmat(1,2)+cj(1)*Mmat(1,2))*b(2) 
rhs(2) = b(4) + (Cmat(2,1)+cj(1)*Mmat(2,1))*b(1) + (Cmat(2,2)+cj(1)*Mmat(2,2))*b(2) 

Kbar = Kmat + cj(1)*Cmat + Cj(1)*cj(1)*Mmat

call gesv(Kbar, rhs, ipiv)

x(1:2) = rhs
x(3:4) = cj(1)*x(1:2) - b(1:2)

ierr = SUN_SUCCESS

end function bushing_MatrixEmbeddedLSSolve
    
!----------------------------------------------------------

integer(c_int) function bushing_MatrixEmbeddedLSFree(LS) &
    result(ierr) bind(C, name='bushing_MatrixEmbeddedLSFree')

type(SUNLinearSolver), target :: LS

LS%content = c_null_ptr
call FSUNLinSolFreeEmpty(LS)

ierr = SUN_SUCCESS

end function bushing_MatrixEmbeddedLSFree
    
!----------------------------------------------------------

subroutine bushing_ForceAndTangent(this, t, y, yp, F, Mmat, Cmat, Kmat)

type(Bushing), intent(in) :: this
real(kind=8), intent(in) :: t
real(kind=8), dimension(4), intent(in) :: y, yp
real(kind=8), dimension(2), intent(out) :: F
real(kind=8), dimension(2,2), intent(out) :: Mmat, Cmat, Kmat

real(kind=8) :: Delt, thet, dDelt, dthet, ddDelt, ddThet
real(kind=8) :: c, s
real(kind=8) :: ux_, uz_
integer :: n

Delt = y(1)
thet = y(2)
dDelt = y(3)
dthet = y(4)
ddDelt = yp(3)
ddthet = yp(4)

c = cos(thet+this%theta0)
s = sin(thet+this%theta0)

! Mass matrix
Mmat(1,1) = this%m
Mmat(2,1) = -this%m*this%h*s
Mmat(1,2) = Mmat(2,1)
Mmat(2,2) = this%II

! Input ground acceleration
n = floor(t/this%dtu) + 1
if (n > size(this%ux)-1) then
    ux_ = 0.d0
    uz_ = 0.d0
else
    ux_ = this%ux(n) + (this%ux(n+1)-this%ux(n))/this%dtu*(t - (n-1)*this%dtu)
    uz_ = this%uz(n) + (this%uz(n+1)-this%uz(n))/this%dtu*(t - (n-1)*this%dtu)
endif

F(1) = Mmat(1,1)*ddDelt + Mmat(1,2)*ddthet &
     - this%m*this%h*c*dthet**2.d0 &
     + this%kv*Delt &
     + this%m*g*(1.d0+uz_) &
     + this%cv*dDelt
F(2) = Mmat(2,1)*ddDelt + Mmat(2,2)*ddthet &
     + this%kr*thet &
     + this%m*this%h*g*(-(1.d0+uz_)*s+ux_*c) &
     + this%cr*dthet

Cmat(1,1) = this%cv
Cmat(1,2) = -2.d0*this%m*this%h*c*dthet
Cmat(2,1) = 0.d0
Cmat(2,2) = this%cr

Kmat(1,1) = this%kv
Kmat(1,2) = this%m*this%h*(-c*ddthet+s*dthet**2.d0)
Kmat(2,1) = 0.d0
Kmat(2,2) = this%kr-this%m*this%h*(c*(g*(1.d0+uz_)+ddDelt)+s*g*ux_)

end subroutine bushing_ForceAndTangent 

!----------------------------------------------------------
    
end module BushingDynamics