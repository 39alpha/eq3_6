subroutine suprdx(jflag,jsflag,narn1a,narn2a,ndrsd,ndrsmx,ndrsrd,nrdxsp,nsta,nstmax,uspeca)
    !! This subroutine executes the option to suppress all redox
    !! reactions. This is not the same as suppressing all redox
    !! species. An auxiliary basis species (say Oxalate-) that is
    !! in the active basis set (and hence has its own mass balance
    !! relation) is detached from any other species in the active
    !! basis set, including the redox species. It should not be
    !! suppressed as part of a general suppression of redox. For
    !! example, if Oxalate- is given its own mass balance relation,
    !! then equilibirum of its associated reaction, which is of redox
    !! type linking it to HCO3-, would be overridden. Tn effect,
    !! Oxalate- would be treated as being composed of a pseudo-element.
    !! This subroutine is called by:
    !!   EQ6/eq6.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: ndrsmx
    integer :: nstmax

    integer :: jflag(nstmax)
    integer :: jsflag(nstmax)
    integer :: ndrsd(ndrsmx)
    integer :: ndrsrd(2,nstmax)

    integer :: narn1a
    integer :: narn2a
    integer :: nrdxsp
    integer :: nsta

    character(len=48) :: uspeca(nstmax)

    ! Local variable declarations.
    integer :: kcount
    integer :: n
    integer :: ncount
    integer :: nr1
    integer :: nr2
    integer :: ns
    integer :: nse

    ncount = 0

    ! First, suppress any species whose associated reaction involves
    ! the redox species, unless equilibrium for that reaction has been
    ! over-ridden. Currently, such an over-ride is possible only for
    ! certain aqueous species. For them, this condition is marked by
    ! jflag is not 30.
    do ns = 1,narn1a - 1
        ! Species preceding the block of aqueous species. There
        ! should be none.
        if (jsflag(ns) .lt. 2) then
            nr1 = ndrsrd(1,ns)
            nr2 = ndrsrd(2,ns)

            do n = nr1 + 1,nr2
                nse = ndrsd(n)

                if (nse .eq. nrdxsp) then
                    jsflag(ns) = 2
                    ncount = ncount + 1
                    go to 100
                end if
            end do

100 continue
        end if
    end do

    do ns = narn1a,narn2a
        ! Aqueous species.
        if (jflag(ns) .eq. 30) then
            if (jsflag(ns) .lt. 2) then
                nr1 = ndrsrd(1,ns)
                nr2 = ndrsrd(2,ns)

                do n = nr1 + 1,nr2
                    nse = ndrsd(n)

                    if (nse .eq. nrdxsp) then
                        jsflag(ns) = 2
                        ncount = ncount + 1
                        go to 110
                    end if
                end do

110 continue
            end if
        end if
    end do

    do ns = narn2a + 1,nsta
        ! Species following the block of aqueous species.
        if (jsflag(ns) .lt. 2) then
            nr1 = ndrsrd(1,ns)
            nr2 = ndrsrd(2,ns)

            do n = nr1 + 1,nr2
                nse = ndrsd(n)

                if (nse .eq. nrdxsp) then
                    jsflag(ns) = 2
                    ncount = ncount + 1
                    go to 120
                end if
            end do

120 continue
        end if
    end do

    if (ncount .le. 0) then
        go to 999
    end if

    ! Now check to insure that any species that are actively linked
    ! to species that are now suppressed are themselves suppressed.
210 continue
    kcount = 0

    do ns = 1,narn1a - 1
        ! Species preceding the block of aqueous species. There
        ! should be none.
        if (jsflag(ns) .lt. 2) then
            nr1 = ndrsrd(1,ns)
            nr2 = ndrsrd(2,ns)

            do n = nr1 + 1,nr2
                nse = ndrsd(n)

                if (jsflag(nse) .ge. 2) then
                    jsflag(ns) = 2
                    kcount = kcount + 1
                    ncount = ncount + 1
                    go to 220
                end if
            end do

220 continue
        end if
    end do

    do ns = narn1a,narn2a
        ! Aqueous species.
        if (jflag(ns) .eq. 30) then
            if (jsflag(ns) .lt. 2) then
                nr1 = ndrsrd(1,ns)
                nr2 = ndrsrd(2,ns)

                do n = nr1 + 1,nr2
                    nse = ndrsd(n)

                    if (jsflag(nse) .ge. 2) then
                        jsflag(ns) = 2
                        kcount = kcount + 1
                        ncount = ncount + 1
                        go to 230
                    end if
                end do

230 continue
            end if
        end if
    end do

    do ns = narn2a + 1,nsta
        ! Species following the block of aqueous species.
        if (jsflag(ns) .lt. 2) then
            nr1 = ndrsrd(1,ns)
            nr2 = ndrsrd(2,ns)

            do n = nr1 + 1,nr2
                nse = ndrsd(n)

                if (jsflag(nse) .ge. 2) then
                    jsflag(ns) = 2
                    kcount = kcount + 1
                    ncount = ncount + 1
                    go to 240
                end if
            end do

240 continue
        end if
    end do

    if (kcount .gt. 0) then
        go to 210
    end if

999 continue
end subroutine suprdx