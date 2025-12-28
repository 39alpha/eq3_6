subroutine bdslx(narn1,narn2,natmax,noutpt,nslt,nsltmx,nslx,nsxi,nsxx,nsxmax,nttyo)
    !! This subroutine builds the nsxi and nsxx arrays. These are
    !! pointer arrays used in connection with the S-lambda
    !! parts of Pitzer's equations.
    !! This subroutine is called by:
    !!   EQ3NR/eq3nr.f
    !!   EQ6/eq6.f
    !! Principal input:
    !!   narn1  = start of aqueous species range
    !!   narn2  = end of aqueous species range
    !!   natmax = the maximum number of aqeuous species
    !!   nslx   = array of aqueous species indices defining pairs
    !!              of species for which S-lambda data are present
    !!   nsltmx = the maximum number of pairs of aqueous species
    !!              for which Pitzer interaction parameters are
    !!              defined; the second dimension of the nslx array
    !!   nsxmax = the second dimension of the nsxx pointer array
    !! Principal output:
    !!   nsxi   = range pointer array into the nsxx array:
    !!              nsxi(1,na) and nsxi(2,na) are the first and last
    !!              values of the second subscript (nsx) of the nsxx
    !!              array for entries pertaining to the na-th
    !!              aqueous species (ns-th species)
    !!   nsxt   = the total number of entries in the nsxx array
    !!   nsxx   = pointer array:
    !!              nsxx(1,nsx) = the species index of the second
    !!              species in the nsl-th pair, and nsxx(2,nsx) = nsl
    !! Note: the usage of the nsxi and nsxx arrays is illustrated by
    !! the following pseudo-code to evaluate the sum:
    !! SUM(j) S-lambda(ij)m(j)
    !!   sum = 0.
    !!   do nsx = nsxi(1,i),nsxi(2,i)
    !!     j = nsxx(1,nsx)
    !!     nsl = nsxx(2,nsx)
    !!     sum = sum + S-lambda(nsl)*m(j)
    !!   enddo
    !! This sum is actually evaluated in EQLIBG/gsgsm.f.
    implicit none

    ! Calling sequence variable declarations.
    integer :: natmax
    integer :: nsltmx
    integer :: nsxmax

    integer :: nslx(2,nsltmx)
    integer :: nsxi(2,natmax)
    integer :: nsxx(2,nsxmax)

    integer :: narn1
    integer :: narn2
    integer :: noutpt
    integer :: nslt
    integer :: nttyo

    ! Local variable declarations.
    integer :: n
    integer :: na
    integer :: nsx
    integer :: ns
    integer :: nsl

    nsx = 1

    ! The following loop assumes that water is the narn1-th species.
    do ns = narn1 + 1,narn2
        na = ns - narn1 + 1

        ! Set the beginning of the range for the current species.
        nsxi(1,na) = nsx

        ! Search column 1 of the nslx array for the ns-th species.
        do nsl = 1,nslt
            if (nslx(1,nsl) .eq. ns) then
                ! Found one. Get the other aqueous species index and the
                ! nsl index of the pair.
                nsxx(1,nsx) = nslx(2,nsl)
                nsxx(2,nsx) = nsl
                nsx = nsx + 1

                if (nsx .gt. nsxmax) then
                    n = 2*nsltmx
                    write (noutpt,1010) nsxmax,n
                    write (nttyo,1010) nsxmax,n
1010 format(/' * Error - (EQLIBG/bdslx) Have overflow of the',/7x,'S-lambda pointer array nsxx. Increase the ','dimensioning',/7x,'parameter nsxpar from ',i5,' to no more than ',i5,'.')

                    stop
                end if
            end if
        end do

        ! Search column 2.
        do nsl = 1,nslt
            if (nslx(2,nsl) .eq. ns) then
                ! Found one. Skip if the species is also in the first column.
                if (ns .ne. nslx(1,nsl)) then
                    nsxx(1,nsx) = nslx(1,nsl)
                    nsxx(2,nsx) = nsl
                    nsx = nsx + 1

                    if (nsx .gt. nsxmax) then
                        n = 2*nsltmx
                        write (noutpt,1010) nsxmax,n
                        write (nttyo,1010) nsxmax,n
                        stop
                    end if
                end if
            end if
        end do

        ! Set the end of the range for the current species.
        nsxi(2,na) = nsx - 1
    end do
end subroutine bdslx