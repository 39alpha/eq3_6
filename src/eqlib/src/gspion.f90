subroutine gspion(narn1,narn2,nchlor,nelect,nhydr,nhydx,noutpt,no2gaq,nstmax,nttyo,uspec)
    !! This subroutine finds the species indices of H+, OH-, Cl-,
    !! aqueous O2(g), and aqueous e-.
    !! This subroutine is called by:
    !!   EQ3NR/eq3nr.f
    !!   EQ6/eq6.f
    !!  Principal input:
    !!    narn1  = start of the range of aqueous species; also the
    !!               species index of solvent water
    !!    narn2  = end of the range of aqueous species
    !!    uspec  = array of species names
    !!  Principal output:
    !!    nchlor = species index of the Cl- ion
    !!    nhydr  = species index of the H+ ion
    !!    nhydx  = species index of the OH- ion
    !!    no2gaq = species index of the fictive species, aqueous O2(g)
    !!    nelect = species index of the fictive species, aqueous e-
    implicit none

    ! Calling sequence variable declarations.
    integer :: nstmax

    integer :: noutpt
    integer :: nttyo

    integer :: narn1
    integer :: narn2
    integer :: nchlor
    integer :: nelect
    integer :: nhydr
    integer :: nhydx
    integer :: no2gaq

    character(len=48) :: uspec(nstmax)

    ! Local variable declarations.
    integer :: j2
    integer :: nerr
    integer :: ns

    integer :: ilnobl

    character(len=8) :: uhydr
    character(len=8) :: uhydx
    character(len=8) :: uchlor
    character(len=8) :: uo2gaq
    character(len=8) :: uelect

    data uhydr  /'H+      '/,uhydx  /'OH-     '/,uchlor /'Cl-     '/,uo2gaq /'O2(g)   '/,uelect /'e-      '/

    nhydr = 0
    nchlor = 0
    nhydx = 0
    no2gaq = 0
    nelect = 0
    nerr = 0

    ! Find H+.
    do ns = narn1,narn2
        if (uspec(ns)(1:8) .eq. uhydr(1:8)) then
            nhydr = ns
            go to 100
        end if
    end do

    j2 = ilnobl(uhydr)
    write (noutpt,1000) uhydr(1:j2)
    write (nttyo,1000) uhydr(1:j2)
1000 format(/' * Error - (EQLIB/gspion) The species ',a," isn't",' present in the model',/7x,'system, as is required. Check to',' see that this species is not missing',/7x,'from the',' supporting data file. If it is there, and the present input',/7x,'file is an EQ6 input file, check to see that this file',' was created',/7x,'by a previous EQ3NR or EQ6 run using the',' same data file.',/7x,'Otherwise,the input file may have been',' corrupted.')

    nerr = nerr + 1

100 continue

    ! Find OH-.
    do ns = narn1,narn2
        if (uspec(ns)(1:8) .eq. uhydx(1:8)) then
            nhydx = ns
            go to 110
        end if
    end do

    j2 = ilnobl(uhydx)
    write (noutpt,1000) uhydx(1:j2)
    write (nttyo,1000) uhydx(1:j2)
    nerr = nerr + 1

110 continue

    ! Find Cl-.
    do ns = narn1,narn2
        if (uspec(ns)(1:8) .eq. uchlor(1:8)) then
            nchlor = ns
            go to 120
        end if
    end do

    j2 = ilnobl(uchlor)
    write (noutpt,1000) uchlor(1:j2)
    write (nttyo,1000) uchlor(1:j2)
    nerr = nerr + 1

120 continue

    ! Find the fictive, aqueous O2(g).
    do ns = narn1,narn2
        if (uspec(ns)(1:5) .eq. uo2gaq(1:5)) then
            no2gaq = ns
            go to 130
        end if
    end do

    j2 = ilnobl(uo2gaq)
    write (noutpt,1000) uo2gaq(1:j2)
    write (nttyo,1000) uo2gaq(1:j2)
    write (noutpt,1010)
    write (nttyo,1010)
1010 format(/7x,'Note- The species O2(g) (the fictive aqueous',' species, not',/7x,'the real gas species) is required to be',' present even if a problem',/7x,'has no redox aspect. The',' species is then inactive.')

    nerr = nerr + 1

130 continue

    ! Find the fictive, aqueous e-.
    do ns = narn1,narn2
        if (uspec(ns)(1:2) .eq. uelect(1:2)) then
            nelect = ns
            go to 140
        end if
    end do

    ! Note: it is currently not an error for this species to not be
    ! present. In future development, it may serve as an alternative
    ! to the fictive, aqueous O2(g).
140 continue

    if (nerr .gt. 0) then
        stop
    end if
end subroutine gspion