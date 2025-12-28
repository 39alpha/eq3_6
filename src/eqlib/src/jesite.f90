integer function jesite(jern1,jern2,jetmax,jgext,ne,netmax,ns)
    !! This subroutine returns the index of the site which contains
    !! the ns-th species, a species of the ne-th generic ion exchanger
    !! phase. If a match is not found, a zero value is returned.
    !! This subroutine is called by:
    !!   None
    !! Input:
    !!   jern1  = array giving the starts of the species ranges of sites
    !!            of generic ion exchanger phases
    !!   jern2  = array giving the ends of the species ranges of sites
    !!            of generic ion exchanger phases
    !!   jgext  = array of the number of sites in generic ion exchanger
    !!              phases
    !!   ne     = the index of the desired generic ion exchanger phase
    !!   ns     = the index of the desired species
    !! Output:
    !!   jesite = the index of the site containing the species whose
    !!              index is ns, and which is a species in the ne-th
    !!              generic ion exchanger phase
    implicit none

    ! Calling sequence variable declarations.
    integer :: jetmax
    integer :: netmax

    integer :: jern1(jetmax,netmax)
    integer :: jern2(jetmax,netmax)
    integer :: jgext(netmax)
    integer :: ne
    integer :: ns

    ! Local variable declarations.
    integer :: je
    integer :: nrr1
    integer :: nrr2

    jesite = 0

    do je = 1,jgext(ne)
        nrr1 = jern1(je,ne)
        nrr2 = jern2(je,ne)

        if (ns .le. nrr2) then
            if (ns .ge. nrr1) then
                jesite = je
                go to 999
            end if
        end if
    end do

999 continue
end function jesite