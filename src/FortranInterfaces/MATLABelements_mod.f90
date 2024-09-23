module MATLABelements
    
use, intrinsic :: iso_c_binding
use emxArray

implicit none

!-------------------------------------------------------------

interface
    subroutine CableForceRotBCinCoord &
        (d1, phi1, gamm1, d2, phi2, gamm2, Pmid, varThetamid, P0, &
         d1dot, phi1dot, gamm1dot, d2dot, phi2dot, gamm2dot, &
         Pmiddot, varThetamiddot, &
         d1ddot, phi1ddot, gamm1ddot,d2ddot, phi2ddot, gamm2ddot, &
         Pmidddot, varThetamidddot, &
         x01, RJ1, RE1, r1, x02, RJ2, RE2, r2, R0, &
         II, rho, EA, EI, GJ, betAX, betBEND, betTOR, &
         sg, wg, nel, colmat, colmat_brev, colmat_bar, d, dbrev, dbar, &
         Mbar, u, Kbar11, Dbar11, dyn, alph0, &
         Fb, Kb, Cb, Mb, Bb) bind(C, name="CableForceRotBCinCoord")
        use, intrinsic :: iso_c_binding
        use emxArray
        implicit none
        real(kind=c_double), dimension(3), intent(in) :: d1, phi1, d2, phi2
        real(kind=c_double), value, intent(in) :: gamm1, gamm2
        type(emxArray_real_T), intent(in) :: Pmid, varThetamid, P0
        real(kind=c_double), dimension(3), intent(in) :: d1dot, phi1dot, d2dot, phi2dot
        real(kind=c_double), value, intent(in) :: gamm1dot, gamm2dot
        type(emxArray_real_T), intent(in) :: Pmiddot, varThetamiddot
        real(kind=c_double), dimension(3), intent(in) :: d1ddot, phi1ddot, d2ddot, phi2ddot
        real(kind=c_double), value, intent(in) :: gamm1ddot, gamm2ddot
        type(emxArray_real_T), intent(in) :: Pmidddot, varThetamidddot
        real(kind=c_double), dimension(3), intent(in) :: x01, r1, x02, r2
        real(kind=c_double), dimension(9), intent(in) :: RJ1, RE1, RJ2, RE2
        type(emxArray_real_T), intent(in) :: R0
        real(kind=c_double), dimension(9), intent(in) :: II
        real(kind=c_double), value, intent(in) :: rho, EA, EI, GJ, betAX, betBEND, betTOR
        type(emxArray_real_T), intent(in) :: sg, wg, nel, colmat, colmat_brev, colmat_bar
        real(kind=c_double), value, intent(in) :: d, dbrev, dbar
        type(emxArray_real_T), intent(in) :: Mbar, Kbar11, Dbar11
        real(kind=c_double), dimension(3), intent(in) :: u
        real(kind=c_double), value, intent(in) :: dyn, alph0
        type(emxArray_real_T), intent(out) :: Fb, Kb, Cb, Mb, Bb
    end subroutine CableForceRotBCinCoord
end interface 

!-------------------------------------------------------------

interface
    subroutine getBishopFrame (P, knots, d, xg, R) &
            bind(C, name="getBishopFrame")
        use, intrinsic :: iso_c_binding
        use emxArray
        implicit none
        type(emxArray_real_T), intent(in) :: P, knots, xg
        real(kind=c_double), value, intent(in) :: d
        type(emxArray_real_T), intent(out) :: R
    end subroutine getBishopFrame
end interface

!-------------------------------------------------------------

interface
    subroutine CableMbar(P0, EA, betAX, colmat, colmat_bar, wg, &
        Mbar, Kbar11, Dbar11) bind(C, name="CableMbar")
        use, intrinsic :: iso_c_binding
        use emxArray
        implicit none
        type(emxArray_real_T), intent(in) :: P0
        real(kind=c_double), value, intent(in) :: EA, betAX
        type(emxArray_real_T), intent(in) :: colmat, colmat_bar, wg
        type(emxArray_real_T), intent(out) :: Mbar, Kbar11, Dbar11
    end subroutine CableMbar
end interface

!-------------------------------------------------------------!-------------------------------------------------------------

interface
    subroutine SplineApproximation(gamm_, J_, N, xg, wg, colmat, p,err) &
        bind(C, name="SplineApproximation")
        use, intrinsic :: iso_c_binding
        use emxArray
        implicit none
        type(emxArray_real_T), intent(in) :: gamm_, J_
        real(kind=c_double), value, intent(in) :: N
        type(emxArray_real_T), intent(in) :: xg, wg, colmat
        type(emxArray_real_T), intent(out) :: p
        real(kind=c_double), intent(out) :: err
    end subroutine SplineApproximation
end interface

!-------------------------------------------------------------!-------------------------------------------------------------

interface
    subroutine RigidBodyForce(d, phi, dd, phid, ddd, phidd,  &
        m, II, KT, CT, KR, CR, rr, RJ, u, &
        F, K, C, MM, B) bind(C, name="RigidBodyForce")
        use, intrinsic :: iso_c_binding
        use emxArray
        implicit none
        real(kind=c_double), dimension(3), intent(in) :: d, phi, dd, phid, ddd, phidd
        real(kind=c_double), value, intent(in) :: m
        real(kind=c_double), dimension(9), intent(in) :: II, KT, CT, KR, CR
        real(kind=c_double), dimension(3), intent(in) :: rr
        real(kind=c_double), dimension(9), intent(in) :: RJ
        real(kind=c_double), dimension(3), intent(in) :: u
        real(kind=c_double), dimension(6), intent(out) :: F
        real(kind=c_double), dimension(36), intent(out) :: K, C, MM
        real(kind=c_double), dimension(18), intent(out) :: B
    end subroutine RigidBodyForce
end interface

!-------------------------------------------------------------

end module MATLABelements