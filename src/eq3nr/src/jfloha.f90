subroutine jfloha(jflgi,nbtd,nbti,nbtmax,noutpt,nttyo,nxmdmx,nxmod,ubasp,uspeci,uxmod)
    !! This subroutine modifies the jflgi array to insure that jflag
    !! defaults to 27, not 30, for the species H2(aq) and O2(aq).
    !! These defaults are desirable because they provide some
    !! poising (buffering of the fO2 or Eh) in EQ6 mode calculations.
    !! However, the jflag is defaulted to -1 if a species is on the
    !! suppression list (nxmod option).
    !! This subroutine is called by:
    !!   EQ3NR/eq3nr.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: nbtmax
    integer :: nxmdmx

    integer :: noutpt
    integer :: nttyo

    integer :: jflgi(nbtmax)
    integer :: nbtd
    integer :: nbti
    integer :: nxmod

    character(len=48) :: ubasp(nbtmax)
    character(len=48) :: uspeci(nbtmax)
    character(len=48) :: uxmod(nxmdmx)

    ! Local variable declarations.
    integer :: n
    integer :: nb
    integer :: nbi

    character(len=24) :: ux24

    ! Process H2(aq).
    ux24 = 'H2(aq)'

    ! Is it on the input file?
    do nbi = 1,nbti
        if (uspeci(nbi)(1:24) .eq. ux24(1:24)) then
            go to 110
        end if
    end do

    ! Is it on the data file?
    do nb = 1,nbtd
        if (ubasp(nb)(1:24) .eq. ux24(1:24)) then
            go to 100
        end if
    end do

    go to 110

    ! See if the input can be extended.
100 continue
    if (nbti .ge. nbtmax) then
        write (noutpt,1000) nbtmax,ux24(1:6)
        write (nttyo,1000) nbtmax,ux24(1:6)
1000 format(/' * Error - (EQ3NR/jfloha) Have exceeded the maximum',/7x,i4,' basis species trying to add ',a,/7x,'to the species read from the input file.',/7x,'Increase the dimensioning parameter nbtpar.')

        stop
    end if

    ! Extend the input with the desired default.
    nbti = nbti + 1
    uspeci(nbti) = ux24

    ! Is it on the suppression list?
    do n = 1,nxmod
        if (uxmod(n)(1:24) .eq. ux24(1:24)) then
            jflgi(nbti) = -1
            go to 110
        end if
    end do

    jflgi(nbti) = 27

110 continue

    ! Process O2(aq).
    ux24 = 'O2(aq)'

    ! Is it on the input file?
    do nbi = 1,nbti
        if (uspeci(nbi)(1:24) .eq. ux24(1:24)) then
            go to 210
        end if
    end do

    ! Is it on the data file?
    do nb = 1,nbtd
        if (ubasp(nb)(1:24) .eq. ux24(1:24)) then
            go to 200
        end if
    end do

    go to 210

    ! See if the input can be extended.
200 continue
    if (nbti .ge. nbtmax) then
        write (noutpt,1000) nbtmax,ux24(1:6)
        write (nttyo,1000) nbtmax,ux24(1:6)
        stop
    end if

    ! Extend the input with the desired default.
    nbti = nbti + 1
    uspeci(nbti) = ux24

    ! Is it on the suppression list?
    do n = 1,nxmod
        if (uxmod(n)(1:24) .eq. ux24(1:24)) then
            jflgi(nbti) = -1
            go to 210
        end if
    end do

    jflgi(nbti) = 27

210 continue
end subroutine jfloha