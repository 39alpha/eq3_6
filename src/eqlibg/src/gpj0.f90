subroutine gpj0(dpj0,d2pj0,pj0,x)
    !! This subroutine evaluates the function J0(x) (pj0) and its
    !! first and second derivatives with respect to x (dpj0 and
    !! d2pj0). The evaluation here uses an approximation given by
    !! Pitzer (1975, eq. 47). Most modern usage of Pitzer's equations
    !! is consistent with the later approximation of Harvie (1981).
    !! These functions are used to compute E-lambda and E-theta
    !! coefficients and their derivatives with respect to ionic
    !! strength. They describe the higher-order electrostatic effects
    !! in Pitzer's equations.
    !! It is important to note that the tabulation of J(x) and J'(x)
    !! in Table II of Pitzer (1975) is not for the approximation
    !! represented by eq. 47. Pitzer (1975) does not provide such a
    !! table for the eq. 47 approximation.
    !!                        References
    !! Harvie, C. E. 1981. Theoretical Investigations in Geochemistry
    !!   and Atom Surface Scattering. Ph.D. dissertation, University
    !!   of California, San Diego (Available as #8203026 from
    !!   University Microfilms International, Ann Arbor, Michigan).
    !! Pitzer, K.S. 1975. Thermodynamics of electrolytes. V. Effects
    !!   of higher-order electrostatic terms. Journal of Solution
    !!   Chemistry, v. 4, p. 249-265.
    !! This subroutine is called by:
    !!   EQLIBG/elmdd.f
    !!   EQLIBG/cwrpjt.f
    !! Input:
    !!   x      = independent variable
    !! Output:
    !!   dpj0   = the derivative dJ0(x)/dx
    !!   d2pj0  = the derivative d2 J0(x)/dx2
    !!   pj0    = the function J0(x)
    implicit none

    ! Calling sequence variable declarations.
    real(kind=8) :: dpj0
    real(kind=8) :: d2pj0
    real(kind=8) :: pj0
    real(kind=8) :: x

    ! Local variable declarations.
    real(kind=8) :: alx
    real(kind=8) :: bk2
    real(kind=8) :: bk3
    real(kind=8) :: c1
    real(kind=8) :: c2
    real(kind=8) :: c3
    real(kind=8) :: c4
    real(kind=8) :: xc2
    real(kind=8) :: c3xc4
    real(kind=8) :: ec3xc4
    real(kind=8) :: t
    real(kind=8) :: xc2p

    ! These parameters are those for eq. 47 given in Table III of
    ! Pitzer (1975).
    data c1,c2,c3,c4 /4.581,0.7237,0.0120,0.528/

    pj0 = 0.
    dpj0 = 0.
    d2pj0 = 0.

    if (x .le. 0.) then
        go to 999
    end if

    alx = log(x)

    xc2 = exp(-c2 * alx)
    c3xc4 = c3*exp(c4 * alx)
    ec3xc4 = exp(-c3xc4)

    pj0 = x / (4. + c1*xc2*ec3xc4)

    t = c4*c3xc4
    xc2p = exp((-c2 - 1.)*alx)
    bk2 = c1*xc2p *(c2 + t)*ec3xc4

    dpj0 = (pj0 / x)*(1. + pj0*bk2)

    bk3 = (t/x)*((c4 / (c2 + t)) - (c2 + 1.) / t)
    d2pj0 = ((dpj0 / pj0) - 1. / x)*(2.*dpj0 + pj0*bk3)

999 continue
end subroutine gpj0