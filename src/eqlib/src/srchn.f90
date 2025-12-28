subroutine srchn(nrn1a,nrn2a,ns,nstmax,unam,uspeca)
    !! This subroutine matches the species name unam with the
    !! corresponding entry in the nrn1a-th through nrn2a-th range of
    !! the species name array uspeca. Only the first 24 characters are
    !! compared (uspeca has 48 characters, unam only 24). This
    !! subroutine returns the species index ns. If there is no match,
    !! ns is returned with a value of 0.
    !! This subroutine is called by:
    !!   EQLIB/inbdot.f
    !!   EQLIB/inupt.f
    !!   EQLIB/srchne.f
    !! Principal input:
    !!   nrn1a  = start of the range of species to search
    !!   nrn2a  = end of the range of species to search
    !!   unam   = name of the species whose index is to be found
    !!   uspeca = array of species names
    !! Principal output:
    !!   ns     = index of the species whose name is unam
    implicit none

    ! Calling sequence variable declarations.
    integer :: nstmax

    integer :: nrn1a
    integer :: nrn2a
    integer :: ns

    character(len=48) :: uspeca(nstmax)
    character(len=24) :: unam

    ! Local variable declarations.
    do ns = nrn1a,nrn2a
        if (uspeca(ns)(1:24) .eq. unam(1:24)) then
            go to 999
        end if
    end do

    ! No match was found.
    ns = 0

999 continue
end subroutine srchn