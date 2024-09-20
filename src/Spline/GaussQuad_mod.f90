module GaussQuad
    
implicit none

    contains

!-----------------------------------------------------    
    
subroutine GaussQuadrature(n, xg, wg)

! Gauss Legendre quadrature points \in (-1,1) and weights

integer, intent(in) :: n
real(kind=8), dimension(n), intent(out) :: xg, wg

integer :: lwork ! unused dummy variable
real(kind=8) :: w ! unused dummy variable 

integer :: ierror

call gaqd(n, xg, wg, w, lwork, ierror) 

xg = cos(xg)
xg(1:n) = xg(n:1:-1) ! flipud

end subroutine GaussQuadrature
    
!-----------------------------------------------------
    
end module GaussQuad    