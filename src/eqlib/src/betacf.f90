subroutine betacf(acflg,acflgo,bacfmx,nst,nstmax,ubacmx,uspec)
    !! This subroutine finds the activity coefficient residual with
    !! the largest magnitude (bacfmx). The range of activity coefficient
    !! residuals covers all species.
    !! This subroutine is called by:
    !!   EQLIB/ngcadv.f
    !! Principal input:
    !!   uspec  = name array for aqueous species
    !!   acflg  = array of current values of the log activity
    !!              coefficients
    !!   acflgo = array of the previous values of the log activity
    !!              coefficients
    !! Principal output:
    !!   bacfmx = the largest activity coefficient residual
    !!   ubacmx = the name of the species with the largest activity
    !!              coefficient residual is bacfmx
    implicit none

    ! Calling sequence variable declarations.
    integer :: nstmax

    integer :: nst

    character(len=48) :: uspec(nstmax)
    character(len=48) :: ubacmx

    real(kind=8) :: acflg(nstmax)
    real(kind=8) :: acflgo(nstmax)
    real(kind=8) :: bacfmx

    ! Local variable declarations.
    integer :: ileft
    integer :: ns

    character(len=48) :: ux48

    real(kind=8) :: adacfl
    real(kind=8) :: bx

    ! Note the use of local variables (bx and ux48) within the loop.
    ! Note also that the loop is unrolled.
    bacfmx = 0.
    ubacmx = 'None'
    bx = 0.
    ux48 = 'None'
    ileft = (nst/8)*8

    do ns = 1,ileft,8
        adacfl = abs(acflg(ns) - acflgo(ns))

        if (adacfl .gt. bx) then
            bx = adacfl
            ux48 = uspec(ns)
        end if

        adacfl = abs(acflg(ns + 1) - acflgo(ns + 1))

        if (adacfl .gt. bx) then
            bx = adacfl
            ux48 = uspec(ns + 1)
        end if

        adacfl = abs(acflg(ns + 2) - acflgo(ns + 2))

        if (adacfl .gt. bx) then
            bx = adacfl
            ux48 = uspec(ns + 2)
        end if

        adacfl = abs(acflg(ns + 3) - acflgo(ns + 3))

        if (adacfl .gt. bx) then
            bx = adacfl
            ux48 = uspec(ns + 3)
        end if

        adacfl = abs(acflg(ns + 4) - acflgo(ns + 4))

        if (adacfl .gt. bx) then
            bx = adacfl
            ux48 = uspec(ns + 4)
        end if

        adacfl = abs(acflg(ns + 5) - acflgo(ns + 5))

        if (adacfl .gt. bx) then
            bx = adacfl
            ux48 = uspec(ns + 5)
        end if

        adacfl = abs(acflg(ns + 6) - acflgo(ns + 6))

        if (adacfl .gt. bx) then
            bx = adacfl
            ux48 = uspec(ns + 6)
        end if

        adacfl = abs(acflg(ns + 7) - acflgo(ns + 7))

        if (adacfl .gt. bx) then
            bx = adacfl
            ux48 = uspec(ns + 7)
        end if
    end do

    do ns = ileft + 1,nst
        adacfl = abs(acflg(ns) - acflgo(ns))

        if (adacfl .gt. bx) then
            bx = adacfl
            ux48 = uspec(ns)
        end if
    end do

    bacfmx = bx
    ubacmx = ux48
end subroutine betacf