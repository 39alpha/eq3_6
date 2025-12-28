subroutine inttfx(narn1a,narn2a,noutpt,nsta_asv,ntfxa,ntfxmx,ntfxta,nttyo,tfxa,uspeca,utfxxd)
    !! This subroutine sets up arrays for handling alkalinity
    !! coefficients Species names used to tag alkalinity coefficients
    !! are matched against the aqueous species names read from the data
    !! file. If a match is not found, the corresponding alkalinity
    !! coefficient is ignored.
    !! This subroutine is called by:
    !!   EQ3NR/eq3nr.f
    !!   EQ6/eq6.f
    !! Principal input:
    !!   uspeca = array of names of species read from the data file
    !!   utfxxd = array of string entries containing the names of
    !!              of species and corresponding alkalinity factors
    !! Principal output:
    !!   ntfxa  = pointer array containing the indices of species
    !!              corresponding to the titration factors in the
    !!              tfxa array
    !!   ntfxta = number of species with alkalinity factors
    !!   tfxa   = array of titration factors
    implicit none

    ! Calling sequence variable declarations.
    integer :: noutpt
    integer :: nttyo

    integer :: nsta_asv
    integer :: ntfxmx

    integer :: ntfxa(ntfxmx)
    integer :: narn1a
    integer :: narn2a
    integer :: ntfxta

    character(len=48) :: uspeca(nsta_asv)
    character(len=32) :: utfxxd(ntfxmx)

    real(kind=8) :: tfxa(ntfxmx)

    ! Local variable declarations.
    integer :: i
    integer :: j
    integer :: j2
    integer :: n
    integer :: nerr
    integer :: nn
    integer :: ns

    integer :: ilnobl

    real(kind=8) :: tfx

    character(len=32) :: uentry
    character(len=32) :: uscr1
    character(len=32) :: uscr2
    character(len=24) :: unam
    character(len=8) :: ux8

    nerr = 0

    ntfxta = 0

    do nn = 1,ntfxmx
        ! Decode the current entry in the pseudo-data file. The
        ! content of an entry is a string of the form: '"CaCO3(aq)", 2.0'.
        ! Double-quote marks surround the name. A comma separates the
        ! name from the number. All blanks are optional.
        uentry = utfxxd(nn)
        i = index(uentry,'"')

        if (i .le. 0) then
            j2 = ilnobl(uentry)
            write (noutpt,1000) uentry(1:j2)
            write (nttyo,1000) uentry(1:j2)
1000 format(/' * Error - (EQ3NR/inttfx) Encountered an error',' while reading a',/7x,'pseudo-data file of alkalinity',' factors. The following entry does not',/7x,'contain an',' initial double-quote mark:',/9x,"'",a,"'")

            nerr = nerr + 1
        end if

        if (nerr .gt. 4) then
            stop
        end if

        ! Okay, found the first double-quote mark.
        uscr1 = uentry(i + 1:32)
        call lejust(uscr1)
        j = index(uscr1,'"')

        if (j .le. 0) then
            j2 = ilnobl(uentry)
            write (noutpt,1010) uentry(1:j2)
            write (nttyo,1010) uentry(1:j2)
1010 format(/' * Error - (EQ3NR/inttfx) Encountered an error',' while reading a',/7x,'pseudo-data file of alkalinity',' factors. The following entry does not',/7x,'contain a',' second double-quote mark:',/9x,"'",a,"'")

            nerr = nerr + 1
        end if

        if (nerr .gt. 4) then
            stop
        end if

        ! Okay, found the second double-quote mark.
        if (j .eq. 1) then
            j2 = ilnobl(uentry)
            write (noutpt,1020) uentry(1:j2)
            write (nttyo,1020) uentry(1:j2)
1020 format(/' * Error - (EQ3NR/inttfx) Encountered an error',' while reading a',/7x,'pseudo-data file of alkalinity',' factors. The following entry does not',/7x,'contain a',' species name between the double-quote marks:',/9x,"'",a,"'")

            nerr = nerr + 1
        end if

        if (nerr .gt. 4) then
            stop
        end if

        unam = uscr1(1:j - 1)

        ! Test for the end of data.
        if (unam(1:6).eq.'endit.') then
            go to 150
        end if

        uscr2 = uscr1(j:32)
        i = index(uscr2,',')

        if (i .le. 0) then
            j2 = ilnobl(uentry)
            write (noutpt,1030) uentry(1:j2)
            write (nttyo,1030) uentry(1:j2)
1030 format(/' * Error - (EQ3NR/inttfx) Encountered an error',' while reading a',/7x,'pseudo-data file of alkalinity',' factors. The following entry does not',/7x,'contain a',' comma after the species name:',/9x,"'",a,"'")

            nerr = nerr + 1
        end if

        if (nerr .gt. 4) then
            stop
        end if

        uscr1 = uscr2(i + 1:32)
        call lejust(uscr1)
        ux8 = uscr1(1:8)

        ! Read the titration factor.
        read (ux8,'(f5.1)',err=100,end=100) tfx
        go to 110

100 continue
        j2 = ilnobl(uentry)
        write (noutpt,1040) uentry(1:j2)
        write (nttyo,1040) uentry(1:j2)
1040 format(/' * Error - (EQ3NR/inttfx) Encountered a read error',' while reading a',/7x,'pseudo-data file of alkalinity',' factors. The following entry does not',/7x,'contain a',' readable titration number:',/9x,"'",a,"'")

        nerr = nerr + 1

        if (nerr .gt. 4) then
            stop
        end if

110 continue

        ! Search for a matching species name.
        ! Calling sequence substitutions:
        !   narn1a for nrn1a
        !   narn2a for nrn2a
        call srchn(narn1a,narn2a,ns,nsta_asv,unam,uspeca)

        if (ns .gt. 0) then
            ! If a previous match was found, skip.
            do n = 1,ntfxta
                if (ns .eq. ntfxa(n)) then
                    go to 140
                end if
            end do

            ntfxta = ntfxta + 1
            ntfxa(ntfxta) = ns
            tfxa(ntfxta) = tfx
140 continue
        end if
    end do

150 continue
    if (nerr .gt. 0) then
        stop
    end if
end subroutine inttfx