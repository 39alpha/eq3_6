subroutine rdpni(abeta,ipbtmx,jpfcmx,nat,natmax,ndat0s,nerr,noutpt,npxni,npx2mx,npx2t,nttyo,nwarn,uaqsp,upair,zaqsp)
    !! This subroutine reads from the DATA1 file the coefficients
    !! required to compute those Pitzer interaction parameters
    !! associated with neutral-cation (nc) and neutral-anion (na)
    !! pairs. This particular set of interaction coefficients
    !! includes the following:
    !!   1. lambda(nc), lambda(na)
    !! Coefficients required to compute Pitzer interaction coefficients
    !! associated with other pair types or with triplets are not read
    !! by this subroutine. This subroutine is one of several that
    !! replace rdpz2.f and rdpz3.f.
    !! This subroutine is called by:
    !!   EQPT/eqpt.f
    !! Principal input:
    !!   ndat0s = unit number of the stripped DATA0 file
    !!   ndat1f = unit number of the DATA1F file
    !! Principal output:
    !!   abeta  = array of coefficients for computing beta parameters
    !!   nerr   = cumulative error counter
    !!   npxni  = the number of nc and na pairs for which Pitzer
    !!              parameters were read from the data file
    !!   npx2t  = the number of pairs of all species types for which
    !!              Pitzer parameters were read from the data file.
    !!   nwarn  = cumulative warning counter
    !!   uaqsp  = array of names of aqueous species
    !!   upair  = array of species names for species pairs
    !!   uaqsp  = array of electrical charge numbers of aqueous species
    implicit none

    ! Calling sequence variable declarations.
    integer :: ipbtmx
    integer :: jpfcmx
    integer :: natmax
    integer :: npx2mx

    integer :: ndat0s
    integer :: noutpt
    integer :: nttyo

    integer :: nat
    integer :: nerr
    integer :: npxni
    integer :: npx2t
    integer :: nwarn

    character(len=24) :: uaqsp(natmax)
    character(len=24) :: upair(2,npx2mx)

    real(kind=8) :: abeta(jpfcmx,0:ipbtmx,npx2mx)
    real(kind=8) :: zaqsp(natmax)

    ! Local variable declarations.
    integer :: ier
    integer :: iz1
    integer :: iz2
    integer :: j
    integer :: j2
    integer :: j3
    integer :: j4
    integer :: j5
    integer :: j6
    integer :: npx2
    integer :: npx2te
    integer :: n1
    integer :: n2

    integer :: ilnobl

    character(len=80) :: udastr
    character(len=80) :: uline
    character(len=80) :: ux80
    character(len=24) :: unam1
    character(len=24) :: unam2
    character(len=24) :: ux24a
    character(len=24) :: ux24b
    character(len=24) :: uhdni
    character(len=24) :: uhdnxt
    character(len=16) :: ustr16
    character(len=8) :: uterm
    character(len=8) :: ux8

    real(kind=8) :: var
    real(kind=8) :: z1
    real(kind=8) :: z2

    data uterm  /'+-------'/
    data uhdni  /"nc and na combinations: "/
    data uhdnxt /"nn combinations: lambda("/

    npx2te = npx2t
    npx2 = npx2te
    j2 = 0
    j3 = 0

    ! Advance to the start of the superblock for these parameters.
100 continue
    read (ndat0s,1000,end=990,err=995) uline
1000 format(a80)

    if (uline(1:24) .ne. uhdni(1:24)) then
        go to 100
    end if

    ! Skip the terminator.
    read (ndat0s,1000,end=990,err=995) uline

    ! Get the names of the species comprising a pair.
110 continue
    read (ndat0s,1010,end=990,err=995) unam1,unam2
1010 format(a24,2x,a24)

    if (unam1(1:24) .eq. uhdnxt(1:24)) then
        ! Found the header for the next superblock of Pitzer parameters.
        ! Back up one line and terminate reading the current superblock.
        npxni = npx2 - npx2t
        npx2t = npx2
        backspace(ndat0s)
        go to 999
    end if

    npx2 = npx2 + 1

    ! Setup up some string data for error messages.
    j2 = ilnobl(unam1)

    if (j2 .le. 0) then
        unam1 = '<blank>'
        j2 = ilnobl(unam1)
    end if

    j3 = ilnobl(unam2)

    if (j3 .le. 0) then
        unam2 = '<blank>'
        j3 = ilnobl(unam2)
    end if

    write (ux8,'(i5)') npx2
    call lejust(ux8)
    j4 = ilnobl(ux8)

    ! Check for blank names.
    if (j2.le.0 .or. j3.le.0) then
        write (noutpt,1020) ux8(1:j4),unam1(1:j2),unam2(1:j3)
        write (nttyo,1020) ux8(1:j4),unam1(1:j2),unam2(1:j3)
1020 format(/' * Error - (EQPT/rdpni) Have encountered a blank',' species name',/7x,'in the species pair for block ',a,' of the superblock for',/7x,'Pitzer nc and na parameter data.',' The offending pair shows as ',a,', ',a,'.')

        if (npx2 .gt. (npx2te + 1)) then
            ux24a = upair(1,npx2 - 1)
            ux24b = upair(2,npx2 - 1)

            if (ux24a(1:7).ne.'<blank>' .or.      ux24b(1:7).ne.'<blank>') then
                j5 = ilnobl(ux24a)
                j6 = ilnobl(ux24b)
                write (noutpt,1030) ux24a(1:j5),ux24b(1:j6)
                write (nttyo,1030) ux24a(1:j5),ux24b(1:j6)
1030 format(7x,'This block follows the one for ',a,', ',a,'.')
            end if
        end if

        nerr = nerr + 1
    end if

    ! Get the indices of the two species (n1, n2). These
    ! indices point to the data read from the previously
    ! read species blocks.
    ! Calling sequence substitutions:
    !   n1 for na
    !   unam1 for unams
    call gspidx(ier,n1,nat,natmax,uaqsp,unam1)

    if (ier .gt. 0) then
        if (unam1(1:7).ne.'<blank>' .and.    unam2(1:7).ne.'<blank>') then
            write (noutpt,1100) unam1(1:j2),unam1(1:j2),unam2(1:j3)
            write (nttyo,1100) unam1(1:j2),unam1(1:j2),unam2(1:j3)
1100 format(/' * Error - (EQPT/rdpni) The aqueous species ',a,' appearing in',/7x,'the Pitzer nc and na parameter data',' block for the species pair',/7x,a,', ',a,' does not',' appear in an aqueous species',/7x,'data block.')

            nerr = nerr + 1
        end if
    end if

    ! Calling sequence substitutions:
    !   n2 for na
    !   unam2 for unams
    call gspidx(ier,n2,nat,natmax,uaqsp,unam2)

    if (ier .gt. 0) then
        if (unam1(1:7).ne.'<blank>' .and.    unam2(1:7).ne.'<blank>') then
            write (noutpt,1100) unam2(1:j3),unam1(1:j2),unam2(1:j3)
            write (nttyo,1100) unam2(1:j3),unam1(1:j2),unam2(1:j3)
            nerr = nerr + 1
        end if
    end if

    ! Check for illegal appearance of solvent water.
    if (unam1(1:4).eq.'H2O ' .or. unam2(1:4).eq.'H2O ' .or.  n1.eq.1 .or. n2.eq.1) then
        write (noutpt,1200) unam1(1:j2),unam2(1:j3)
        write (nttyo,1200) unam1(1:j2),unam2(1:j3)
1200 format(/' * Error - (EQPT/rdpni) Have found an illegal',' data block for the',/7x,'species pair ',a,', ',a,' in the superblock of Pitzer data',/7x,"for nc and na",' pairs. Solvent water may not appear in such a pair.')

        nerr = nerr + 1
    end if

    ! Get the electrical charges of the species.
    z1 = 0.
    z2 = 0.
    iz1 = 0
    iz2 = 0

    if (n1 .gt. 0) then
        z1 = zaqsp(n1)
        iz1 = nint(z1)
    end if

    if (n2 .gt. 0) then
        z2 = zaqsp(n2)
        iz2 = nint(z2)
    end if

    if (n1.gt.0 .and. n2.gt.0) then
        ! Make sure that the species pair is appropriate for the
        ! current superblock. It must be of type nc or na. Other
        ! pair types such as ca, cc', aa', nn, and nn' are not
        ! permitted here.
        if ((iz1.eq.0 .and. iz2.eq.0) .or.    (iz1.ne.0 .and. iz2.ne.0)) then
            ! The current data block is in the wrong superblock.
            write (noutpt,1230) unam1(1:j2),unam2(1:j3)
            write (nttyo,1230) unam1(1:j2),unam2(1:j3)
1230 format(/' * Error - (EQPT/rdpni) Have found an',' illegal data block for the',/7x,'species pair ',a,', ',a,' in the',/7x,'superblock of Pitzer data for'," nc and na pairs.",/7x,'The present pair',' type may not appear in the current superblock.')

            nerr = nerr + 1
        end if
    end if

    ! Species pair storage rules:
    !   1. neutral < cation < anion
    !   2. If both species are of the same charge type,
    !      store alphabetically.
    ! Here: the neutral comes first.
    if (iz1 .eq. 0) then
        ! The first species is the neutral.
        upair(1,npx2) = unam1
        upair(2,npx2) = unam2
    else
        ! The second species is the neutral.
        upair(1,npx2) = unam2
        upair(2,npx2) = unam1
    end if

    ! Read the coefficients for the lambda(nc) and lambda(na)
    ! parameters. Store the lambda data as beta(0) data.
    ! Check the required parameter header.
    read (ndat0s,1000,end=990,err=995) uline
    j4 = ilnobl(uline)
    j4 = min(j4,70)
    ux80 = uline
    call lejust(ux80)
    ustr16 = 'lambda'
    j5 = 6

    if (ux80(1:j5) .ne. ustr16(1:j5)) then
        write (noutpt,1330) uline(1:j4),ustr16(1:j5),unam1(1:j2),unam2(1:j3)
        write (nttyo,1330) uline(1:j4),ustr16(1:j5),unam1(1:j2),unam2(1:j3)
1330 format(/' * Error - (EQPT/rdpni) Have found a line',' starting with:',/7x,'"',a,'"',/7x,'where "',a,'" was',' expected in the block for the species',/7x,'pair ',a,', ',a,' in the superblock of Pitzer data',/7x,"for nc and na pairs.")

        nerr = nerr + 1
        go to 999
    end if

    ! Read the coefficients for the associated temperature function.
    do j = 1,jpfcmx
        read (ndat0s,1000,end=990,err=995) uline
        j4 = ilnobl(uline)
        j4 = min(j4,70)
        ux80 = uline
        call lejust(ux80)
        ustr16 = 'a  ='
        j5 = 4
        write (ustr16(2:2),'(i1)') j

        if (ux80(1:j5) .eq. ustr16(1:j5)) then
            udastr = ux80
            udastr(1:j5) = ' '
            call g1dat(ier,noutpt,nttyo,udastr,var)

            if (ier .gt. 0) then
                write (noutpt,1320) uline(1:j4),unam1(1:j2),unam2(1:j3)
                write (nttyo,1320) uline(1:j4),unam1(1:j2),unam2(1:j3)
1320 format(/' * Error - (EQPT/rdpni) Have found a line',' starting with:',/7x,'"',a,'"',/7x,'containing an',' expected numerical field that could not be read.',/7x,'This occurred in the block for the species pair',/7x,a,', ',a,' in the superblock of Pitzer data',/7x,"for nc and na pairs.")

                nerr = nerr + 1
                go to 999
            end if

            abeta(j,0,npx2) = var
        else
            write (noutpt,1330) uline(1:j4),ustr16(1:j5),unam1(1:j2),unam2(1:j3)
            write (nttyo,1330) uline(1:j4),ustr16(1:j5),unam1(1:j2),unam2(1:j3)
            nerr = nerr + 1
            go to 999
        end if
    end do

    ! Skip past the block delimiter line.
    read (ndat0s,1000,end=990,err=995) uline

    if (uline(1:8) .ne. uterm(1:8)) then
        j4 = ilnobl(uline)
        j4 = min(j4,70)
        write (noutpt,1470) uline(1:j4),unam1(1:j2),unam2(1:j3)
        write (nttyo,1470) uline(1:j4),unam1(1:j2),unam2(1:j3)
1470 format(/' * Error - (EQPT/rdpni) Have found a line starting',' with:',/7x,'"',a,'"',/7x,' in the data block for the species',' pair ',a,', ',a,' in the',/7x,'superblock of Pitzer data for'," nc and na pairs. Should have",/7x,'found the',' delimiter line marking the end of the block.')

        nerr = nerr + 1
        go to 999
    end if

    ! Process the next data block.
    go to 110

990 continue
    write (noutpt,2010)
    write (nttyo,2010)
2010 format(/' * Error - (EQPT/rdpni) Unexpectedly encountered',/7x,"an end-of-file error while reading the nc and na",/7x,'(lambda) superblock of the DATA0 file.')

    stop

995 continue
    write (noutpt,2020)
    write (nttyo,2020)
2020 format(/' * Error - (EQPT/rdpni) Encountered a read format',/7x,"error while reading the nc and na (lambda) superblock",/7x,'of the DATA0 file.')

    stop

999 continue
end subroutine rdpni