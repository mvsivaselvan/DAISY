program main

use, intrinsic :: iso_c_binding
use fsundials_core_mod
use fida_mod
use fnvector_serial_mod
use fsunnonlinsol_newton_mod

use BushingDynamics

implicit none

real(kind=8), parameter :: pi = 3.14159265358979d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Bushing properties
type(Bushing), target :: bush
real(kind=8), parameter :: m = 2.41d0 ! lb-s^2/in
real(kind=8), parameter :: II = 4500.d0 ! lb-s^2/in in^2
real(kind=8), parameter :: kv = 6950.d0 ! lb/in
real(kind=8), parameter :: kr = 3245000.d0 ! lb-in/rad
real(kind=8), parameter :: zetav = 0.006d0 ! damping ratio for vertical motion
real(kind=8), parameter :: zetar = 0.008d0 ! damping ratio for rocking motion
real(kind=8) :: cv
real(kind=8) :: cr
real(kind=8), parameter :: h = 20.d0 ! in
real(kind=8), parameter :: theta0 = 20.d0/180.d0*pi
character(20), parameter :: ufile = 'eqdata.dat'
real(kind=8), parameter :: dtu = 1.d0/256.d0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! State variables and tolerances
integer, parameter :: neq = 4
real(kind=8), dimension(neq) :: y, yp
real(kind=8), parameter :: rtol = 1.d-6
real(kind=8), dimension(neq) :: atol = (/1.d-8, 1.d-8, 1.d-8, 1.d-8/)
real(kind=8) :: t0 = 0.d0
real(kind=8) :: tout, tret(1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SUNDIALS stuff
integer(c_int) :: retval
type(c_ptr) :: sunctx ! SUNDIALS simulation context
type(N_Vector), pointer :: sunvec_y, sunvec_yp
type(N_Vector), pointer :: sunvec_atol
type(c_ptr) :: ida_mem 
type(SUNLinearSolver), pointer :: linsolver
type(SUNMatrix), pointer :: amat ! not used since using custom linear solver
type(SUNNonlinearSolver), pointer :: NLS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialize bushing object
cv = zetav*2*sqrt(kv*m);
cr = zetar*2*sqrt(kr*II);
call bushing_create(bush, m , II, kv, kr, cv, cr, h, theta0, ufile, dtu, &
                    1.d0, 1.d0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Setup IDA
! STEP 1
retval = FSUNContext_Create(SUN_COMM_NULL, sunctx)

! STEP 2 - initial conditions + create SUNvectors
call bushing_getInitialDisplacement(bush, y(1:2))
print*,y(2)
y(3:4) = 0.d0 ! initial velocity
yp = 0.d0 ! initial velocity and acceleration
sunvec_y => FN_VMake_serial(neq, y, sunctx)
sunvec_yp => FN_VMake_serial(neq, yp, sunctx)
sunvec_atol => FN_VMake_serial(neq, atol, sunctx)

! STEP 3 - IDACreate, IDAInit, set tolerances, set user data
ida_mem = FIDACreate(sunctx)
retval = FIDAInit(ida_mem, c_funloc(bushing_res), t0, sunvec_y, sunvec_yp)
retval = FIDASVtolerances(ida_mem, rtol, sunvec_atol)
retval = FIDASetUserData(ida_mem, c_loc(bush))

! Step 4 - Construct and attach custom linear solver
call bushing_MatrixEmbeddedLS(ida_mem, sunctx, linsolver)
nullify(amat)
retval = FIDASetLinearSolver(ida_mem, linsolver, amat)

! Step 5 - Create and set nonlinear solver (Newton is used by default, so
!          not necessary to do this
NLS => FSUNNonlinSol_Newton(sunvec_y, sunctx)
retval = FIDASetNonlinearSolver(ida_mem, NLS)

! Step 6 - Take one step
open(unit=100, file='idaout.dat', status='unknown')
tout = 0.d0
write(100,'(9e13.5)')tout,y(1:4),yp(1:4)
do while (tout<=50.d0)
    tout = tout + 0.01d0
    retval = FIDASolve(ida_mem, tout, tret, sunvec_y, sunvec_yp, IDA_NORMAL)
    write(100,'(9e13.5)')tout,y(1:4),yp(1:4)
    print*,tout
enddo
close(unit=100)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Start cleanup
call bushing_destroy(bush)

call FIDAFree(ida_mem)
retval = FSUNNonlinSolFree(NLS)
retval = FSUNLinSolFree(linsolver)
call FN_VDestroy(sunvec_y)
call FN_VDestroy(sunvec_yp)
call FN_VDestroy(sunvec_atol)
retval = FSUNContext_Free(sunctx)
! End cleanup
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

print*,"Hello world"
    
end program main