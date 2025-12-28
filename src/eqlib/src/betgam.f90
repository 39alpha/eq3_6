subroutine betgam(acflg,acflgo,bgamx,narn1,narn2,nstmax,ubgamx,uspec)
    !! This subroutine finds the aqueous species activity coefficient
    !! residual with the largest magnitude (bgamx).
    !! This subroutine is called by:
    !!   EQLIB/ngcadv.f
    !!   EQ6/optmzr.f
    !! Principal input:
    !!   uspec  = name array for aqueous species
    !!   acflg  = array of current values of the log activity
    !!              coefficients
    !!   acflgo = array of the previous values of the log activity
    !!              coefficients
    !! Principal output:
    !!   bgamx  = the largest activity coefficient residual
    !!   ubgamx = the name of the species with the largest activity
    !!              coefficient residual is bgamx
    implicit none

    ! Calling sequence variable declarations.
    integer :: nstmax

    integer :: narn1
    integer :: narn2

    character(len=48) :: uspec(nstmax)
    character(len=48) :: ubgamx

    real(kind=8) :: acflg(nstmax)
    real(kind=8) :: acflgo(nstmax)
    real(kind=8) :: bgamx

    ! Local variable declarations.
    integer :: ileft
    integer :: ns
    integer :: nval

    character(len=48) :: ux48

    real(kind=8) :: adacfl
    real(kind=8) :: bx

    ! Note the use of local variables (bx and ux48) within the loop.
    ! Note also that the loop is unrolled.
    bgamx = 0.
    ubgamx = 'None'
    bx = 0.
    ux48 = 'None'
    nval = narn2 - narn1 + 1
    ileft = (nval/8)*8 + narn1 -1

    do ns = narn1,ileft,8
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

    do ns = ileft + 1,narn2
        adacfl = abs(acflg(ns) - acflgo(ns))

        if (adacfl .gt. bx) then
            bx = adacfl
            ux48 = uspec(ns)
        end if
    end do

    bgamx = bx
    ubgamx = ux48
end subroutine betgam