subroutine intieb(iebal,iebal3,ier,jsflag,nbasp,nbt,nbtmax,noutpt,nstmax,nttyo,uebal,uspec,zchar)
    !! This subroutine finds the basis index (iebal) of the species to be
    !! adjusted for electrical balance, if any.
    !! This subroutine is called by:
    !!   EQ3NR/eq3nr.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: nbtmax
    integer :: nstmax

    integer :: jsflag(nstmax)
    integer :: nbasp(nbtmax)
    integer :: iebal
    integer :: iebal3
    integer :: ier
    integer :: nbt
    integer :: noutpt
    integer :: nttyo

    character(len=48) :: uspec(nstmax)
    character(len=24) :: uebal

    real(kind=8) :: zchar(nstmax)

    ! Local variable declarations.
    integer :: j2
    integer :: nb
    integer :: nebal
    integer :: ns

    integer :: ilnobl

    ! if iebal3 is 0 or iebal = 1 but uebal is 'none' or blank, do no
    ! electrical balancing.
    iebal = 0

    if (iebal3 .le. 0) then
        uebal(1:24) = 'None'
    end if

    if (uebal(1:8) .eq. '        ') then
        uebal(1:24) = 'None'
    end if

    if (uebal(1:8) .eq. 'none    ') then
        uebal(1:24) = 'None'
    end if

    if (uebal(1:8) .eq. 'NONE    ') then
        uebal(1:24) = 'None'
    end if

    if (uebal(1:5) .eq. 'None ') then
        go to 999
    end if

    j2 = ilnobl(uebal)

    ! Find the basis index of the species specified on the input file.
    do nb = 1,nbt
        ns = nbasp(nb)

        if (uebal(1:24).eq.uspec(ns)(1:24)) then
            iebal = nb
            nebal = ns
            go to 100
        end if
    end do

    write (noutpt,1000) uebal(1:j2)
    write (nttyo,1000) uebal(1:j2)
1000 format(/' * Error - (EQ3NR/intieb) The input file specifies',/7x,'that the concentration of ',a,' is to be adjusted',/7x,'to achieve electrical balance. However, this species',/7x,"isn't in the basis set.")

    ier = 1
100 continue

    ! Test the specified ion to see if it is okay.
    ! Does the ion selected for electrical balancing have charge?
    if (zchar(nebal) .eq. 0.) then
        write (noutpt,1010) uebal(1:j2)
        write (nttyo,1010) uebal(1:j2)
1010 format(/' * Warning - (EQ3NR/intieb) The input file specifies',/7x,'that the concentration of ',a,' is to be adjusted',/7x,'to achieve electrical balance. However, this species',/7x,"doesn't have an electrical charge. The success of this",/7x,'calculation will depend on the concentration adjustment',/7x,'influencing the concentration of one or more charged',/7x,'dependent species.')
    end if

    ! Is the ion selected for electrical balancing in the model?
    if (jsflag(nebal) .gt. 0) then
        write (noutpt,1020) uebal(1:j2)
        write (nttyo,1020) uebal(1:j2)
1020 format(/' * Error - (EQ3NR/intieb) The input file specifies',/7x,'that the concentration of ',a,' is to be adjusted',/7x,'to achieve electrical balance. However, this species',/7x,"is effectively suppressed and thus can't be present",/7x,'in the model.')

        ier = 1
    end if

999 continue
end subroutine intieb