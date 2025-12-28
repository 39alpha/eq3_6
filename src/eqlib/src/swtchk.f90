subroutine swtchk(cdrs,jflag,jsflag,nbaspx,nbt,nbtmax,ndrs,ndrsmx,ndrsr,noutpt,ns1,ns2,nstmax,nttyo,uspec)
    !! This subroutine checks a proposed basis switch for EQLIB/switch.f
    !! to see if doing the switch is ok or not.
    !! This subroutine is called by:
    !!   EQLIB/switch.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: nbtmax
    integer :: ndrsmx
    integer :: nstmax

    integer :: noutpt
    integer :: nttyo

    integer :: jflag(nstmax)
    integer :: jsflag(nstmax)
    integer :: nbaspx(nbtmax)
    integer :: ndrs(ndrsmx)
    integer :: ndrsr(2,nstmax)
    integer :: nbt
    integer :: ns1
    integer :: ns2

    character(len=48) :: uspec(nstmax)

    real(kind=8) :: cdrs(ndrsmx)

    ! Local variable declarations.
    integer :: jlen1
    integer :: jlen2
    integer :: nb2
    integer :: ncmpt2
    integer :: nerr

    integer :: nbasis

    character(len=56) :: usp156
    character(len=56) :: usp256

    real(kind=8) :: cx12

    real(kind=8) :: coefdr

    nerr = 0

    ! Calling sequence substitutions:
    !   jlen1 for jlen
    !   uspec(ns1) for unam48
    !   usp156 for uspn56
    call fmspnx(jlen1,uspec(ns1),usp156)

    ! Calling sequence substitutions:
    !   jlen2 for jlen
    !   uspec(ns2) for unam48
    !   usp256 for uspn56
    call fmspnx(jlen2,uspec(ns2),usp256)

    if (ns1 .eq. ns2) then
        write (noutpt,1000) usp156(1:jlen1)
        write (nttyo,1000) usp156(1:jlen1)
1000 format(/" * Error - (EQLIB/swtchk) Can't replace the species ",a,/7x,'with itself in the basis set.')

        nerr = nerr + 1
    end if

    ! Make sure that the first species is not suppressed
    if (jsflag(ns1) .gt. 0) then
        write (noutpt,1010) usp156(1:jlen1),usp256(1:jlen2)
        write (nttyo,1010) usp156(1:jlen1),usp256(1:jlen2)
1010 format(/" * Error - (EQLIB/swtchk) Can't replace the species ",a,/7x,'in the basis set with ',a,', because the former is',' suppressed.')

        nerr = nerr + 1
    end if

    ! Make sure that the second species is not suppressed
    if (jsflag(ns2) .gt. 0) then
        write (noutpt,1020) usp156(1:jlen1),usp256(1:jlen2)
        write (nttyo,1020) usp156(1:jlen1),usp256(1:jlen2)
1020 format(/" * Error - (EQLIB/swtchk) Can't replace the species ",a,/7x,'in the basis set with ',a,', because the latter is',' suppressed.')

        nerr = nerr + 1
    end if

    ! Check the linking reaction.
    ncmpt2 = ndrsr(2,ns2) - ndrsr(1,ns2) + 1

    if (ncmpt2 .lt. 2) then
        ! The second species must not be in the strict basis, because
        ! there is then no possibility of a linking reaction.
        write (noutpt,1030) usp156(1:jlen1),usp256(1:jlen2)
        write (nttyo,1030) usp156(1:jlen1),usp256(1:jlen2)
1030 format(/" * Error - (EQLIB/swtchk) Can't replace the species ",a,/7x,'in the basis set with ',a,', because the latter is',' already in',/7x,'the strict basis set, so there is no',' possibility of a valid',/7x,'linking reaction.')
    else
        ! Make sure that the first species appears as a product in the
        ! reaction belonging to the second species.
        ! Calling sequence substitutions:
        !   ns1 for nse
        !   ns2 for ns
        cx12 = coefdr(cdrs,ndrs,ndrsmx,ndrsr,ns1,ns2,nstmax)

        if (cx12 .eq. 0.) then
            write (noutpt,1040) usp156(1:jlen1),usp256(1:jlen2)
            write (nttyo,1040) usp156(1:jlen1),usp256(1:jlen2)
1040 format(/" * Error - (EQLIB/swtchk) Can't replace the species ",a,/7x,'in the basis set with ',a,', because the former does',' not appear',/7x,'in the linking reaction.')

            ! Calling sequence substitutions:
            !   noutpt for nf
            !   ns2 for ns
            call prreac(cdrs,ndrs,ndrsmx,ndrsr,noutpt,ns2,nstmax,uspec)
            nerr = nerr + 1
        end if
    end if

    if (jflag(ns2) .ne. 30) then
        ! Calling sequence substitutions:
        !   nbaspx for nbasp
        !   ns2 for ns
        nb2 = nbasis(nbaspx,nbt,nbtmax,ns2)

        if (nb2 .gt. 0) then
            write (noutpt,1050) usp156(1:jlen1),usp256(1:jlen2)
            write (nttyo,1050) usp156(1:jlen1),usp256(1:jlen2)
1050 format(/" * Error - (EQLIB/swtchk) Can't replace the",' species ',a,/7x,'in the basis set with ',a,', because',' the latter is already in',/7x,'the active basis set.')

            nerr = nerr + 1
        end if
    end if

    if (nerr .gt. 0) then
        stop
    end if
end subroutine swtchk