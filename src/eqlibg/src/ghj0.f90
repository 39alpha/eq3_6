subroutine ghj0(dhj0,d2hj0,hj0,hj1,hj2,x)
    !! This subroutine evaluates the function J0(x) (hj0), its first
    !! two derivatives with respect to x (dhj0 and d2hj0), and
    !! the related functions J1(x) (hj1) and J2(x) (hj2). The
    !! algorithm used here is based on Chebyshev polynomials and is
    !! taken from the Appendix B of Harvie (1981). The treatment of
    !! the second derivative and of J2(x) is an extension of Harvie's
    !! work. The Harvie (1981) approximation is considered to be
    !! more accurate than the earlier approximation given by Pitzer
    !! (1975) and is used in almost all modern work involving
    !! Pitzer's equations.
    !! These functions are used to compute E-lambda and E-theta
    !! coefficients and their derivatives with respect to ionic
    !! strength. They describe the higher-order electrostatic effects
    !! in Pitzer's equations.
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
    !!   dhj0   = the derivative dJ0(x)/dx
    !!   d2hj0  = the derivative d2 J0(x)/dx2
    !!   hj0    = the function J0(x)
    !!   hj1    = the function J1(x)
    !!   hj2    = the function J2(x)
    implicit none

    ! Calling sequence variable declarations.
    real(kind=8) :: dhj0
    real(kind=8) :: d2hj0
    real(kind=8) :: hj0
    real(kind=8) :: hj1
    real(kind=8) :: hj2
    real(kind=8) :: x

    ! Local variable declarations.
    integer :: k
    integer :: ki

    real(kind=8) :: ai(0:20)
    real(kind=8) :: aii(0:20)
    real(kind=8) :: b(0:22)
    real(kind=8) :: d(0:22)
    real(kind=8) :: e(0:22)
    real(kind=8) :: bdif
    real(kind=8) :: ddif
    real(kind=8) :: edif
    real(kind=8) :: dz
    real(kind=8) :: d2z
    real(kind=8) :: z

    data b(21) /0./,b(22) /0./
    data d(21) /0./,d(22) /0./
    data e(21) /0./,e(22) /0./

    data (ai(k),k = 0,9) /1.925154014814667,-.060076477753119,-.029779077456514,-.007299499690937,.000388260636404,.000636874599598,.000036583601823,-.000045036975204,-.000004537895710,.000002937706971/

    data (ai(k), k = 10,20) /.000000396566462,-.000000202099617,-.000000025267769,.000000013522610,.000000001229405,-.000000000821969,-.000000000050847,.000000000046333,.000000000001943,-.000000000002563,-.000000000010991/

    data (aii(k), k = 0,9) /.628023320520852,.462762985338493,.150044637187895,-.028796057604906,-.036552745910311,-.001668087945272,.006519840398744,.001130378079086,-.000887171310131,-.000242107641309/

    data (aii(k), k = 10,20) /.000087294451594,.000034682122751,-.000004583768938,-.000003548684306,-.000000250453880,.000000216991779,.000000080779570,.000000004558555,-.000000006944757,-.000000002849257,.000000000237816/

    hj0 = 0.
    dhj0 = 0.
    d2hj0 = 0.

    if (x .le. 0.) then
        go to 999
    end if

    if (x .le. 1.) then
        ! Case 1.
        z = 4.*(x**0.2) - 2.
        dz = 0.8*(x**(-0.8))
        d2z = 0.64*(x**(-1.8))

        do ki = 0,20
            k = 20 - ki
            b(k) = z*b(k + 1) - b(k + 2) + ai(k)
            d(k) = b(k + 1) + z*d(k + 1) - d(k + 2)
            e(k) = 2.*d(k + 1) + z*e(k + 1) - e(k + 2)
        end do
    else
        ! Case 2.
        z = 4.44444444444444*(x**(-0.1)) - 2.44444444444444
        dz = -0.444444444444444*(x**(-1.1))
        d2z = 0.488888888888889*(x**(-2.1))

        do ki = 0,20
            k = 20 - ki
            b(k) = z*b(k + 1) - b(k + 2) + aii(k)
            d(k) = b(k + 1) + z*d(k + 1) - d(k + 2)
            e(k) = 2.*d(k + 1) + z*e(k + 1) - e(k + 2)
        end do
    end if

    ! Compute (b0 - b2), (d0 - d2), and (e0 - e2).
    bdif = b(0) - b(2)
    ddif = d(0) - d(2)
    edif = e(0) - e(2)

    hj0 = 0.25*x - 1.0 + bdif/2.
    dhj0 = 0.25 + 0.5*dz*ddif
    d2hj0 = 0.5*(edif*(dz**2) + ddif*d2z)

    hj1 = x*dhj0
    hj2 = x*d2hj0

999 continue
end subroutine ghj0
