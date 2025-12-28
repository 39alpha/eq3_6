subroutine srchne(nrn1a,nrn2a,ns,nsta_asv,unam,uspeca)
    !! This subroutine matches the species name unam with the
    !! corresponding entry in the nrn1a-th through nrn2a-th range of
    !! the species name array uspeca. Only the first 24 characters are
    !! compared (uspeca has 48 characters, unam only 24). This
    !! subroutine returns the species index ns. If there is no match,
    !! the species name is converted to all lower case and the code again
    !! searches for a match. If none is found, the species name is
    !! coverted to all upper case and the code again searches for a
    !! match. If no match is found, ns is returned with a value of 0.
    !! This subroutine is much like EQLIB/srchn.f, except that if no
    !! match is found, the above-noted case conversions are tried in
    !! an effort to find a match.
    !! This subroutine is called by:
    !!   None
    !! Principal input:
    !!   nrn1a  = start of the range of species to search
    !!   nrn2a  = end of the range of species to search
    !!   unam   = name of the species whose index is to be found
    !!   uspeca = array of species names
    !! Principal output:
    !!   ns     = index of the species whose name is unam
    implicit none

    ! Calling sequence variable declarations.
    integer :: nsta_asv

    integer :: nrn1a
    integer :: nrn2a
    integer :: ns

    character(len=48) :: uspeca(nsta_asv)
    character(len=24) :: unam

    ! Local variable declarations.
    integer :: k
    integer :: nfa
    integer :: nfacap
    integer :: nfdif
    integer :: nfz
    integer :: nfzcap
    integer :: nlettr

    character(len=24) :: unamlc
    character(len=24) :: unamuc
    character(len=1) :: ulettr

    ! Determine the ranges of upper and lower case letters.
    nfa = ichar('a')
    nfz = ichar('z')
    nfacap = ichar('A')
    nfzcap = ichar('Z')
    nfdif = ichar('A') - ichar('a')

    ! Try the species name as is.
    call srchn(nrn1a,nrn2a,ns,nsta_asv,unam,uspeca)

    if (ns .le. 0) then
        ! Convert to lower case.
        unamlc = unam

        do k = 1,24
            ulettr = unam(k:k)
            nlettr = ichar(ulettr)

            if (nlettr.ge.nfacap .and. nlettr.le.nfzcap) then
                unamlc(k:k) = char(nlettr - nfdif)
            end if
        end do

        ! Test to see if this is different. If so, try it.
        if (unamlc .ne. unam) then
            ! Calling sequence substitutions:
            !   unamlc for unam
            call srchn(nrn1a,nrn2a,ns,nsta_asv,unamlc,uspeca)
        end if
    end if

    if (ns .le. 0) then
        ! Convert to upper case.
        unamuc = unam

        do k = 1,24
            ulettr = unam(k:k)
            nlettr = ichar(ulettr)

            if (nlettr.ge.nfa .and. nlettr.le.nfz) then
                unamuc(k:k) = char(nlettr + nfdif)
            end if
        end do

        ! Test to see if this is different. If so, try it.
        if (unamuc .ne. unam) then
            ! Calling sequence substitutions:
            !   unamuc for unam
            call srchn(nrn1a,nrn2a,ns,nsta_asv,unamuc,uspeca)
        end if
    end if
end subroutine srchne