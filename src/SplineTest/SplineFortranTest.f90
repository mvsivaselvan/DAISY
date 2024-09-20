!  Spline.f90 
!
!  FUNCTIONS:
!  Spline - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: Spline
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program Spline

    use BSpline
    
    implicit none

    integer :: N, k
    real(kind=8), dimension(14) :: t
    real(kind=8) :: x
    integer :: nderiv
    integer :: left
    real(kind=8), dimension(4,4) :: a
    real(kind=8), dimension(4,3) :: dbiatx
    
    real(kind=8), dimension(21) :: xg
    real(kind=8), dimension(63,10) :: colmat
    
    integer :: i, j
    
    N = 10 ! number of control points
    k = 4 ! order for cubic B-spline
    
    ! Construct knot vector
    t(1:3) = 0.d0
    do i = 4,11
        t(i) = (i-4.d0)/7.d0
    enddo
    t(12:14) = 0.d0
    
    x = 0.5d0
    nderiv = 3 ! Function, derivative, second derivative
    left = 7
    
    ! Test bsplvd
    print*,'Testing bsplvd ...'
    call bsplvd(t, k, x, left, a, dbiatx, nderiv)
    do i = 1,4
        write(*,*)(dbiatx(i,j),j=1,3)
    enddo

    ! Test spcol
    print*,'Testing spcol ...'
    xg = (/0.0161, 0.0714, 0.1268, 0.1590, 0.2143, 0.2696, 0.3018, 0.3571, 0.4125, &
           0.4447, 0.5000, 0.5553, 0.5875, 0.6429, 0.6982, 0.7304, 0.7857, 0.8410, &
           0.8732, 0.9286, 0.9839/)
    call spcol(t, k, xg, nderiv, colmat)
    open(file='colmat1.dat',unit=100,status='unknown')
    do i = 1,size(xg)*nderiv
        write(100,'(10E13.5)')(colmat(i,j),j=1,size(t)-k)
    enddo
    close(unit=100)
    
    ! Body of Spline
    print *, 'Hello World'

    end program Spline

