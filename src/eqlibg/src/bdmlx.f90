subroutine bdmlx(narn1,narn2,natmax,nmut,nmutmx,nmux,nmxi,nmxmax,nmxx,noutpt,nttyo)
    !! This subroutine builds the nmxi and nmxx arrays. These are
    !! pointer arrays used in connection with the mu parts of
    !! Pitzer's equations.
    !! This subroutine is called by:
    !!   EQ3NR/eq3nr.f
    !!   EQ6/eq6.f
    !! Principal input:
    !!   narn1  = start of range of aqueous species
    !!   narn2  = end of range of aqueous species
    !!   natmax = the maximum number of aqeuous species
    !!   nmux   = array of aqueous species indices defining triplets
    !!              of species for which mu data are present
    !!   nmutmx = the maximum number of triplets of aqueous species
    !!              for which Pitzer interaction parameters are
    !!              defined; the second dimension of the nmux array
    !!   nmxmax = the second dimension of the nmxx pointer array
    !! Principal output:
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
    !! Note: the usage of the nmxi and nmxx arrays is illustrated by
    !! the following pseudo-code to evaluate SUM(jk) mu(ijk)m(j)m(k),
    !! where i, j, and k are species indices and na is the aqueous
    !! species index of the i-th species:
    !!   sum = 0.
    !!   i = na + narn1 - 1
    !!   do nmx = nmxi(1,na),nmxi(2,na)
    !!     j = nmxx(1,nmx)
    !!     k = nmxx(2,nmx)
    !!     nmu = nmxx(3,nmx)
    !!     sum = sum + mu(nmu)*m(j)*m(k)
    !!   enddo
    !!   sum = 2.*sum
    !! This sum is actually evaluated in EQLIBG/gmdsm.f.
    implicit none

    ! Calling sequence variable declarations.
    integer :: natmax
    integer :: nmutmx
    integer :: nmxmax

    integer :: nmux(3,nmutmx)
    integer :: nmxi(2,natmax)
    integer :: nmxx(3,nmxmax)

    integer :: narn1
    integer :: narn2
    integer :: nmut
    integer :: noutpt
    integer :: nttyo

    ! Local variable declarations.
    integer :: n
    integer :: na
    integer :: nmu
    integer :: nmx
    integer :: ns

    nmx = 1

    ! The following loop assumes that water is the narn1-th species.
    do ns = narn1 + 1,narn2
        na = ns - narn1 + 1

        ! Set the beginning of the range for the current species.
        nmxi(1,na) = nmx

        ! Search column 1 of the nmux array for the na-th aqueous species.
        do nmu = 1,nmut
            if (nmux(1,nmu) .eq. ns) then
                ! Found one. Get the other two aqueous species indices and the
                ! nmx index of the triplet.
                nmxx(1,nmx) = nmux(2,nmu)
                nmxx(2,nmx) = nmux(3,nmu)
                nmxx(3,nmx) = nmu
                nmx = nmx + 1

                if (nmx .gt. nmxmax) then
                    n = 3*nmxmax
                    write (noutpt,1010) nmxmax,n
                    write (nttyo,1010) nmxmax,n
1010 format(/' * Error - (EQLIBG/bdmlx) Have overflow of the',/7x,'mu pointer array nmxx. Increase the dimensioning',/7x,'parameter nmxpar from ',i5,' to no more than ',i5,'.')

                    stop
                end if
            end if
        end do

        ! Search column 2.
        do nmu = 1,nmut
            if (nmux(2,nmu) .eq. ns) then
                ! Found one. Skip if the same species is also in the first
                ! column.
                if (ns .ne. nmux(1,nmu)) then
                    nmxx(1,nmx) = nmux(1,nmu)
                    nmxx(2,nmx) = nmux(3,nmu)
                    nmxx(3,nmx) = nmu
                    nmx = nmx + 1

                    if (nmx .gt. nmxmax) then
                        n = 3*nmxmax
                        write (noutpt,1010) nmxmax,n
                        write (nttyo,1010) nmxmax,n
                        stop
                    end if
                end if
            end if
        end do

        ! Search column 3.
        do nmu = 1,nmut
            if (nmux(3,nmu) .eq. ns) then
                ! Found one. Skip if the species is also in column 1 or
                ! column 2.
                if (ns.ne.nmux(1,nmu) .and. ns.ne.nmux(2,nmu)) then
                    nmxx(1,nmx) = nmux(1,nmu)
                    nmxx(2,nmx) = nmux(2,nmu)
                    nmxx(3,nmx) = nmu
                    nmx = nmx + 1

                    if (nmx .gt. nmxmax) then
                        n = 3*nmxmax
                        write (noutpt,1010) nmxmax,n
                        write (nttyo,1010) nmxmax,n
                        stop
                    end if
                end if
            end if
        end do

        ! Set the end of the range for the current species.
        nmxi(2,na) = nmx - 1
    end do
end subroutine bdmlx