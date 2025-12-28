subroutine ptzchk(narn1,narn2,natmax,nmxi,noutpt,nstmax,nsxi,nttyo,uspec)
    !! This subroutine checks each aqueous species in the current model
    !! to see if it has any S-lambda or mu Pitzer coefficients. If
    !! not, a warning is printed. If a species has at least one
    !! coefficient of either type but no S-lambda coefficients (i.e.,
    !! has only one or more mu coefficients), a warning is printed.
    !! No warning is printed if a species has one or more S-lambda
    !! coefficients but no mu coefficients.
    !! This subroutine is called by:
    !!   EQ3NR/eq3nr.f
    !!   EQ6/eq6.f
    !! Principal input:
    !!   narn1  = start of aqueous species range
    !!   narn2  = end of aqueous species range
    !!   uspec  = names of all species
    !!   nmxi   = range pointer array into the nmxx array:
    !!              nmxi(1,na) and nmxi(2,na) are the first and last
    !!              values of the second subscript (nmx) of the nmxx
    !!              array for entries pertaining to the na-th
    !!              aqueous species (ns-th species)
    !!   nsxi   = range pointer array into the nsxx array:
    !!              nsxi(1,na) and nsxi(2,na) are the first and last
    !!              values of the second subscript (nsx) of the nsxx
    !!              array for entries pertaining to the na-th
    !!              aqueous species (ns-th species)
    !! Principal output:
    !!    None
    !! Not used here, but referenced above:
    !!   nmxx   = pointer array:
    !!              nmxx(1,nmx) = the species index of the second
    !!              species in the nmu-th triplet, nmxx(2,nmx) is the
    !!              species index of the third species, and
    !!              nmxx(3,nmx) = nmu
    !!   nsxx   = pointer array:
    !!              nsxx(1,nsx) = the species index of the second
    !!              species in the nsl-th pair, where
    !!              nsxx(2,nsx) = nsl
    implicit none

    ! Calling sequence variable declarations.
    integer :: natmax
    integer :: nstmax

    integer :: nmxi(2,natmax)
    integer :: nsxi(2,natmax)
    integer :: narn1
    integer :: narn2
    integer :: noutpt
    integer :: nttyo

    character(len=48) :: uspec(nstmax)

    ! Local variable declarations.
    integer :: j2
    integer :: islt
    integer :: imut
    integer :: na
    integer :: ns
    integer :: ilnobl

    ! The following loop assumes that water is the narn1-th species.
    do ns = narn1 + 1,narn2
        na = ns - narn1 + 1

        islt = nsxi(2,na) - nsxi(1,na) + 1
        imut = nmxi(2,na) - nmxi(1,na) + 1

        if (islt.eq.0 .and. imut.eq.0) then
            if (uspec(ns)(1:6).ne.'O2(g) ' .and.      uspec(ns)(1:3).ne.'e- ') then
                j2 = ilnobl(uspec(ns)(1:24))
                write (noutpt,2000) uspec(ns)(1:j2)
                write (nttyo ,2000) uspec(ns)(1:j2)
2000 format(/' * Warning - (EQLIBG/ptzchk) The species ',a,/7x,"has no Pitzer interaction coefficients.")
            end if
        else if (islt .eq. 0) then
            j2 = ilnobl(uspec(ns)(1:24))
            write (noutpt,2010) uspec(ns)(1:j2)
            write (nttyo ,2010) uspec(ns)(1:j2)
2010 format(/' * Warning - (EQLIBG/ptzchk) The species ',a,/7x,"has no Pitzer S-lambda interaction coefficients.")
        end if
    end do
end subroutine ptzchk