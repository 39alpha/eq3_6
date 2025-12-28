subroutine intbsw(nbasp,nbaspx,nbt,nbtmax,nobswt,noutpt,nst,nstmax,nttyo,uobsw,uspec)
    !! This subroutine interprets ordinary basis switching directives
    !! specified on the input file.
    !! This subroutine is called by:
    !!   EQ3NR/eq3nr.f
    !!   EQ6/eq6.f
    !! Principal input:
    !!   nbasp  = the indices of the species in the basis set
    !!              (modified by this subroutine)
    !!   nbt    = the number of species in the basis set
    !!   nobswt = the number of basis switches to make
    !!   uobsw  = the name pairs for the species to be switched
    !!   uspec  = array of species names
    !! Principal output:
    !!   nbasp  = the indices of the species in the basis set
    !!              (modified by this subroutine)
    implicit none

    ! Calling sequence variable declarations.
    integer :: nbtmax
    integer :: nstmax

    integer :: noutpt
    integer :: nttyo

    integer :: nbasp(nbtmax)
    integer :: nbaspx(nbtmax)
    integer :: nobswt
    integer :: nbt
    integer :: nst

    character(len=48) :: uobsw(2,nbtmax)
    character(len=48) :: uspec(nstmax)

    ! Local variable declarations.
    integer :: jlen
    integer :: n
    integer :: nb
    integer :: nerr
    integer :: ns1
    integer :: ns2

    character(len=56) :: uspn56
    character(len=48) :: unam48

    nerr = 0

    do n = 1,nobswt
        ! Find the first species.
        do nb = 1,nbt
            ns1 = nbaspx(nb)

            if (uspec(ns1)(1:48) .eq. uobsw(1,n)(1:48)) then
                go to 110
            end if
        end do

        nerr = nerr + 1

        unam48 = uobsw(1,n)
        call fmspnx(jlen,unam48,uspn56)
        write (noutpt,1000) uspn56(1:jlen)
        write (nttyo,1000) uspn56(1:jlen)
1000 format(/' * Error - (EQLIB/intbsw) The species ',a,/7x,'is specified to be replaced in an ordinary basis switch,',/7x,'but it is not in the active basis set.')

        go to 150

110 continue

        ! Find the second species.
        do ns2 = 1,nst
            if (uspec(ns2)(1:48) .eq. uobsw(2,n)(1:48)) then
                go to 130
            end if
        end do

        nerr = nerr + 1

        unam48 = uobsw(2,n)
        call fmspnx(jlen,unam48,uspn56)
        write (noutpt,1010) uspn56(1:jlen)
        write (nttyo,1010) uspn56(1:jlen)
1010 format(/' * Error - (EQLIB/intbsw) The species ',a,/7x,'is specified to be put into the basis set by an ordinary',/7x,'basis switch but it is not in the set of active',' species.')

        go to 150

130 continue

        ! Mark to switch.
        nbasp(nb) = ns2

150 continue
    end do

    ! Stop if name match errors have been encountered.
    if (nerr .gt. 0) then
        write (noutpt,1020) nerr
        write (nttyo,1020) nerr
1020 format(/' * Error - (EQLIB/intbsw) ',i4,' errors were',/7x,'encountered interpreting ordinary basis switches that',/7x,'were directed on the input file.')

        stop
    end if
end subroutine intbsw