subroutine intsbs(nb1,nb2,nbaspd,nbtd,nbtmax,noutpt,ns1,ns2,nsbsw,nstmax,nttyo,usbsw,uspeca)
    !! This subroutine interprets the nsbsw-th directive read from the
    !! input file for effecting special basis switching. The roles of
    !! a strict basis species and an auxiliary basis species are
    !! switched by this process. The name of the former species is
    !! usbsw(1,nsbsw), that of the latter is usbsw(2,nsbsw). Here nb1
    !! is the basis index of the strict basis species and nb2 is that
    !! of the auxiliary basis species. Here also ns1 and ns2 are the
    !! corresponding species indices.
    !! This subroutine is called by:
    !!   EQ3NR/eq3nr.f
    !!   EQ6/eq6.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: nbtmax
    integer :: nstmax

    integer :: noutpt
    integer :: nttyo

    integer :: nbaspd(nbtmax)
    integer :: nb1
    integer :: nb2
    integer :: nbtd
    integer :: ns1
    integer :: ns2
    integer :: nsbsw

    character(len=48) :: usbsw(2,nbtmax)
    character(len=48) :: uspeca(nstmax)

    ! Local variable declarations.
    integer :: jlen
    integer :: j2
    integer :: nerr

    integer :: ilnobl

    character(len=48) :: ublk48
    character(len=48) :: unam48
    character(len=56) :: uspn56
    character(len=8) :: ux8

    data ublk48/'                                                '/

    nerr = 0

    do nb1 = 1,nbtd
        ns1 = nbaspd(nb1)

        if (usbsw(1,nsbsw)(1:48) .eq. uspeca(ns1)(1:48)) then
            go to 120
        end if
    end do

    nerr = nerr + 1

    unam48 = usbsw(1,nsbsw)
    call fmspnx(jlen,unam48,uspn56)
    write (noutpt,1000) uspn56(1:jlen)
    write (nttyo,1000) uspn56(1:jlen)
1000 format(/' * Error - (EQLIB/intsbs) Invalid special basis',' directive on the',/7x,'input file: ',a,' is not in the basis',' set and therefore',/7x,'can not be specified in a special',' basis switch.')

120 continue

    if (usbsw(2,nsbsw)(1:48) .eq. ublk48(1:48)) then
        nb2 = nb1
        ns2 = ns1
    else
        do nb2 = 1,nbtd
            ns2 = nbaspd(nb2)

            if (usbsw(2,nsbsw)(1:48) .eq. uspeca(ns2)(1:48)) then
                go to 140
            end if
        end do

        nerr = nerr + 1

        unam48 = usbsw(2,nsbsw)
        call fmspnx(jlen,unam48,uspn56)
        write (noutpt,1000) uspn56(1:jlen)
        write (nttyo,1000) uspn56(1:jlen)

140 continue
    end if

    ! Stop if name match errors have been encountered.
    if (nerr .gt. 0) then
        write (ux8,'(i5)') nerr
        call lejust(ux8)
        j2 = ilnobl(ux8)
        write (noutpt,1020) ux8(1:j2)
        write (nttyo,1020) ux8(1:j2)
1020 format(/' * Error - (EQLIB/intsbs) ',a,' errors were',/7x,'encountered in interpreting special basis switches',/7x,'specified on the input file.')

        stop
    end if
end subroutine intsbs