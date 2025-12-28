subroutine grdxsp(nbasp,nbt,nbtmax,nct,ndrsr,noutpt,nrdxsp,nstmax,nttyo,uspec)
    !! This subroutine finds the redox basis species (nrdxsp).
    !! This subroutine is called by:
    !!   EQ3NR/eq3nr.f
    !!   EQ6/eq6.f
    !! Principal input:
    !!   nbasp  = array of indices of basis species
    !! Principal output:
    !!   nrdxsp = species index of the redox basis species (= 0 if
    !!            there is none)
    implicit none

    ! Calling sequence variable declarations.
    integer :: nbtmax
    integer :: nstmax

    integer :: noutpt
    integer :: nttyo

    integer :: nbasp(nbtmax)
    integer :: ndrsr(2,nstmax)
    integer :: nbt
    integer :: nct
    integer :: nrdxsp

    character(len=48) :: uspec(nstmax)

    ! Local variable declarations.
    integer :: jlen
    integer :: nb
    integer :: nr1
    integer :: nr2
    integer :: ns
    integer :: nt1

    character(len=56) :: uspn56

    ! The redox basis species, if it exists, folllows the other strict
    ! basis species, each of which must associate one-to-one with a
    ! chemical element. The strategy is to test the basis species after
    ! those to see if it is a strict basis species (has no associated
    ! real reaction). If that is the case, then it is the redox species
    ! (species index nrdxsp). If the redox species doesn't exist, nrdxsp
    ! is returned as zero.
    nrdxsp = 0

    ! Get the candidate basis index (nb).
    nb = nct + 1

    if (nb .le. nbt) then
        ! This index is in range.
        ns = nbasp(nb)
        nr1 = ndrsr(1,ns)
        nr2 = ndrsr(2,ns)
        nt1 = nr2 - nr1 + 1

        if (nt1 .lt. 2) then
            ! The candidate species is a strict basis species.
            ! It must be the redox species.
            nrdxsp = ns
        end if
    end if

    if (nrdxsp .gt. 0) then
        ! Calling sequence substitutions:
        !   uspec(nrdxsp) for unam48
        call fmspnx(jlen,uspec(nrdxsp),uspn56)

        write (noutpt,1000) uspn56(1:jlen)
        write (nttyo,1000) uspn56(1:jlen)
1000 format(/' The redox basis species is ',a,'.')
    else
        write (noutpt,1010)
        write (nttyo,1010)
1010 format(' No redox basis species was found.')
    end if
end subroutine grdxsp