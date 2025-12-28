subroutine elmdd(aphi,el,elp,elpp,fxi,ijz,qpit75)
    !! This subroutine calculates the E-lambda function (el) and its
    !! first two derivatives (elp and elpp) with respect to ionic
    !! strength (fxi) for the charge pair product ijz. The E-lambda
    !! function and its derivatives are used to evaluate Pitzer's
    !! equations. Here aphi is the Debye-Huckel A(phi) parameter.
    !! Note: el = elp = elpp = 0 if ijz is less than or equal to 0.
    !! This subroutine is called by:
    !!   EQLIBG/gelam.f.
    !! Principal input:
    !!   aphi   = Debye-Huckel A(phi) parameter
    !!   fxi    = the ionic strength (the 2nd-order electrostatic
    !!              moment function I)
    !!   ijz    = input charge product
    !! Principal output:
    !!   el     = array of values of E-lambda functions
    !!   elp    = array of first ionic strength derivatives of
    !!              E-lambda functions
    !!   elpp   = array of second ionic strength derivatives of
    !!              E-lambda functions
    implicit none

    ! Calling sequence variable declarations.
    integer :: ijz

    logical :: qpit75

    real(kind=8) :: aphi
    real(kind=8) :: el
    real(kind=8) :: elp
    real(kind=8) :: elpp
    real(kind=8) :: fxi

    ! Local variable declarations.
    real(kind=8) :: ck
    real(kind=8) :: cx
    real(kind=8) :: dhj0
    real(kind=8) :: d2hj0
    real(kind=8) :: hj0
    real(kind=8) :: hj1
    real(kind=8) :: hj2
    real(kind=8) :: ri
    real(kind=8) :: sjp
    real(kind=8) :: sjpp
    real(kind=8) :: x
    real(kind=8) :: fxisqt
    real(kind=8) :: xp
    real(kind=8) :: xpp

    el = 0.
    elp = 0.
    elpp = 0.

    if (ijz .gt. 0.) then
        fxisqt = sqrt(fxi)
        cx = 3.*ijz*aphi
        x = 2.*cx*fxisqt

        ! Compute x' and x''.
        xp = cx/fxisqt
        xpp = -(0.5*xp)/fxi

        ! Get J0(x) and related functions.
        if (.not.qpit75) then
            ! Use the Chebyshev polynomial approximation
            ! of Harvie (1981). This is what should normally
            ! be used. Note that hj1 and hj2 are outputs
            ! that are not used.
            call ghj0(dhj0,d2hj0,hj0,hj1,hj2,x)
        else
            ! Use the less accurate approximation of
            ! Pitzer (1975).
            ! Calling sequence substitutions:
            !   dhj0 for dpj0
            !   d2hj0 for d2pj0
            !   hj0 for pj0
            call gpj0(dhj0,d2hj0,hj0,x)
        end if

        ! Convert these to derivatives with respect to I.
        sjp = dhj0*xp
        sjpp = dhj0*xpp + d2hj0*xp*xp

        ck = 0.25*ijz
        ri = 1./fxi

        el = ri*ck*hj0
        elp = ri*(ck*sjp - el)
        elpp = ri*(ck*sjpp - 2.*elp)
    end if
end subroutine elmdd