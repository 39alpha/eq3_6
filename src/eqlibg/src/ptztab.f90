subroutine ptztab(iopr,narn1,narn2,natmax,nmutmx,nmux,nmxi,nmxmax,nmxx,noprmx,noutpt,nsltmx,nslx,nstmax,nsxi,nsxmax,nsxx,uspec)
    !! This subroutine tabulates the species combinations corresponding
    !! to coefficients for Pitzer's equations. The tabulation is
    !! controlled by iopr(10):
    !!    0 = Don't print
    !!    1 = Print a summary of the names of the species present and
    !!          the number of Pitzer interaction coefficients
    !!    2 = Print a summary of the names of the species present and
    !!          the number of Pitzer interaction coefficients
    !! This subroutine is called by:
    !!   EQ3NR/eq3nr.f
    !!   EQ6/eq6.f
    !! Principal input:
    !!   iopr   = array of print option switches
    !!   narn1  = start of aqueous species range
    !!   narn2  = end of aqueous species range
    !!   nmux   = array of aqueous species indices defining triplets
    !!              of species for which mu data are present
    !!   nmutmx = the maximum number of triplets of aqueous species
    !!              for which Pitzer interaction parameters are
    !!              defined; the second dimension of the nmux array
    !!   nmxi   = range pointer array into the nmxx array:
    !!              nmxi(1,na) and nmxi(2,na) are the first and last
    !!              values of the second subscript (nmx) of the nmxx
    !!              array for entries pertaining to the na-th
    !!              aqueous species (ns-th species)
    !!   nmxx   = pointer array:
    !!              nmxx(1,nmx) = the species index of the second
    !!              species in the nmu-th triplet, nmxx(2,nmx) is the
    !!              species index of the third species, and
    !!              nmxx(3,nmx) = nmu
    !!   nslx   = array of aqueous species indices defining pairs
    !!              of species for which S-lambda data are present
    !!   nsltmx = the maximum number of pairs of aqueous species
    !!              for which Pitzer interaction parameters are
    !!              defined; the second dimension of the nslx array
    !!   nsxi   = range pointer array into the nsxx array:
    !!              nsxi(1,na) and nsxi(2,na) are the first and last
    !!              values of the second subscript (nsx) of the nsxx
    !!              array for entries pertaining to the na-th
    !!              aqueous species (ns-th species)
    !!   nsxx   = pointer array:
    !!              nsxx(1,nsx) = the species index of the second
    !!              species in the nsl-th pair, and nsxx(2,nsx) = nsl
    !!   uspec  = names of species
    !! Principal output:
    !!    None
    implicit none

    ! Calling sequence variable declarations.
    integer :: natmax
    integer :: nmutmx
    integer :: nmxmax
    integer :: noprmx
    integer :: nsltmx
    integer :: nstmax
    integer :: nsxmax

    integer :: noutpt

    integer :: iopr(noprmx)
    integer :: nmux(3,nmutmx)
    integer :: nmxi(2,natmax)
    integer :: nmxx(3,nmxmax)
    integer :: nslx(2,nsltmx)
    integer :: nsxi(2,natmax)
    integer :: nsxx(2,nsxmax)
    integer :: narn1
    integer :: narn2

    character(len=48) :: uspec(nstmax)

    ! Local variable declarations.
    integer :: j2
    integer :: j3
    integer :: j4
    integer :: islt
    integer :: imut
    integer :: ixf
    integer :: ixl
    integer :: na
    integer :: nmu
    integer :: nmx
    integer :: ns
    integer :: nsl
    integer :: nsx
    integer :: ns1
    integer :: ns2
    integer :: ns3

    integer :: ilnobl

    if (iopr(10) .le. 0) then
        go to 999
    end if

    ! The following loop assumes that water is the narn1-th species.
    if (iopr(10) .eq. 1) then
        ! Write a table giving the S-lambda and mu tallies for each
        ! aqueous species in the current model.
        write (noutpt,1000)
1000 format(/11x,'--- Pitzer Interaction Coefficient Summary ---',/4x,'Species',15x,'S-lambda Sets',3x,'Mu Sets',/)

        ! The following loop assumes that water is the narn1-th species.
        do ns = narn1 + 1, narn2
            if (uspec(ns)(1:6).ne.'O2(g) ' .and.      uspec(ns)(1:3).ne.'e- ') then
                na = ns - narn1 + 1
                islt = nsxi(2,na) - nsxi(1,na) + 1
                imut = nmxi(2,na) - nmxi(1,na) + 1
                write (noutpt,1010) uspec(ns),islt,imut
1010 format(2x,a24,2x,i4,7x,i4)
            end if
        end do

        write (noutpt,1020)
1020 format(/1x)
    end if

    if (iopr(10) .ge. 2) then
        ! Write a table giving the S-lambda and mu tallies for eachi
        ! aqueous species in the current model. Also list the other
        ! species involved in the S-lambda and mu combinations.
        write (noutpt,1100)
1100 format(/15x,'--- Pitzer Interaction Coefficient Sets ---')

        ! The following loop assumes that water is the narn1-th species.
        do ns = narn1 + 1, narn2
            if (uspec(ns)(1:6).ne.'O2(g) ' .and.      uspec(ns)(1:3).ne.'e- ') then
                na = ns - narn1 + 1
                j2 = ilnobl(uspec(ns)(1:24))
                write (noutpt,1110) uspec(ns)(1:j2)
1110 format(//' Coefficients for ',a,':',/)

                ixf = nsxi(1,na)
                ixl = nsxi(2,na)
                islt = ixl - ixf + 1
                write (noutpt, 1120) islt
1120 format(3x,'No. of S-lambda sets= ',i4,':',/)

                if (islt .gt. 0) then
                    do nsx = ixf,ixl
                        nsl = nsxx(2,nsx)
                        ns1 = nslx(1,nsl)
                        ns2 = nslx(2,nsl)
                        j2 = ilnobl(uspec(ns1)(1:24))
                        j3 = ilnobl(uspec(ns2)(1:24))
                        write (noutpt,1130) uspec(ns1)(1:j2),uspec(ns2)(1:j3)
1130 format(5x,a,', ',a)
                    end do
                end if

                ixf = nmxi(1,na)
                ixl = nmxi(2,na)
                imut = ixl - ixf + 1
                write (noutpt, 2530) imut
2530 format(/3x,'No. of mu sets= ',i4,':',/)

                if (imut .gt. 0) then
                    do nmx = ixf,ixl
                        nmu = nmxx(3,nmx)
                        ns1 = nmux(1,nmu)
                        ns2 = nmux(2,nmu)
                        ns3 = nmux(3,nmu)
                        j2 = ilnobl(uspec(ns1)(1:24))
                        j3 = ilnobl(uspec(ns2)(1:24))
                        j4 = ilnobl(uspec(ns3)(1:24))
                        write (noutpt,2540) uspec(ns1)(1:j2),uspec(ns2)(1:j3),uspec(ns3)(1:j4)
2540 format(5x,a,', ',a,', ',a)
                    end do
                end if
            end if
        end do

        write (noutpt,1020)
    end if

999 continue
end subroutine ptztab