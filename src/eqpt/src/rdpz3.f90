subroutine rdpz3(apsi,atheta,jpfcmx,nat,natmax,ndat0s,nerr,noutpt,npx3mx,npx3t,nthdt,nttyo,nwarn,uaqsp,uethfl,uthdtr,utripl,zaqsp)
    !! This suboutine reads the following Pitzer interaction parameter
    !! data from the DATA0 file:
    !!   S-theta(MM'X) (and S-theta (MXX')) and their temperature
    !!     derivatives
    !!   psi(MM'X) (and psi (MXX')) and their temperature derivatives
    !! (Here M = a cation, M' = a second cation, X = an anion, and
    !! X' = a second anion).
    !! This suboutine is called by:
    !!   EQPT/eqpt.f
    !! Principal input:
    !!   ndat0s = unit number of the stripped DATA0 file
    !! Principal output:
    !!   apsi   = array of coefficients for computing psi parameters
    !!   atheta = array of coefficients for computing S-theta parameters
    !!   nerr   = cumulative error counter
    !!   nwarn  = cumulative warning counter
    !!   nthdt  = the number of distinct S-theta parameter entries
    !!              (i.e., not counting duplicates)
    !!   npx3t  = the number of triplets of species for which Pitzer
    !!              interaction parameters are defined
    !!   utripl = array of names in a species triplet
    !!   uthdtr = array of names of cc'a and aa'c triplets which
    !!              provide the values used for coefficients for
    !!              computing S-theta parameters for the cc' and aa'
    !!              pairs
    implicit none

    ! Calling sequence variable declarations.
    integer :: jpfcmx
    integer :: natmax
    integer :: npx3mx

    integer :: ndat0s
    integer :: noutpt
    integer :: nttyo

    integer :: nat
    integer :: nerr
    integer :: npx3t
    integer :: nthdt
    integer :: nwarn

    character(len=24) :: uaqsp(natmax)
    character(len=24) :: uthdtr(3,npx3mx)
    character(len=24) :: utripl(3,npx3mx)
    character(len=8) :: uethfl

    real(kind=8) :: atheta(jpfcmx,npx3mx)
    real(kind=8) :: apsi(jpfcmx,npx3mx)
    real(kind=8) :: zaqsp(natmax)

    ! Local variable declarations.
    integer :: i
    integer :: ier
    integer :: ifound
    integer :: iz1
    integer :: iz2
    integer :: iz3
    integer :: j
    integer :: j2
    integer :: j3
    integer :: j4
    integer :: j5
    integer :: j6
    integer :: j7
    integer :: j8
    integer :: na
    integer :: nc
    integer :: nn
    integer :: npx3
    integer :: nthd
    integer :: n1
    integer :: n2
    integer :: n3

    integer :: ilnobl

    logical :: qdup12
    logical :: qdup13
    logical :: qdup23

    character(len=80) :: uline
    character(len=24) :: unam1
    character(len=24) :: unam2
    character(len=24) :: unam3
    character(len=24) :: uthd1
    character(len=24) :: uthd2
    character(len=24) :: uthd3
    character(len=24) :: ux24a
    character(len=24) :: ux24b
    character(len=24) :: ux24c
    character(len=24) :: u1
    character(len=24) :: u2
    character(len=24) :: u3
    character(len=8) :: uelemt
    character(len=8) :: uterm
    character(len=8) :: ux8

    real(kind=8), dimension(:), allocatable :: athetx

    real(kind=8) :: z1
    real(kind=8) :: z2
    real(kind=8) :: z3

    data uelemt / 'elements'/
    data uterm  / '+-------'/

    ! Allocate work array.
    ALLOCATE(athetx(jpfcmx))

    nthdt = 0
    npx3 = 0

    ! Read the E-theta flag.
    uethfl = ' '
    read (ndat0s,1000,end=990,err=995) uline
1000 format(a80)

    i = index(uline,'E-theta flag')

    if (i .eq. 0) then
        write (noutpt,1010)
        write (nttyo,1010)
1010 format(/' * Error - (EQPT/rdpz3) The E-theta flag line is',' missing from',/7x,'the DATA0 file. This line is specific to DATA0 files',/7x,"tied to Pitzer's equations. It should appear after",/7x,'the line containing "mixture parameters". The E-theta',/7x,'flag line should contain "E-theta flag = on" or',/7x,'"E-theta flag = off". A value of "on" causes inclusion',/7x,'of the higher-order electrostatic (E-theta) terms in',/7x,"Pitzer's",' equations. A value of "off" causes the',/7x,"exclusion of these terms. Most models based on Pitzer's",/7x,'equations are consistent with the inclusion of these',/7x,'terms, but some are not. The values of the mixture',/7x,'coefficients theta and psi are consistent with either',/7x,'the inclusion or the exclusion of these terms.')

        stop
    end if

    i = index(uline,'= on')

    if (i .gt. 0) then
        uethfl = uline(i + 2:i + 3)
    else
        i = index(uline,'= off')

        if (i .gt. 0) then
            uethfl = uline(i + 2:i + 4)
        else
            write (noutpt,1020) uline
            write (nttyo,1020) uline
1020 format(/' * Error - (EQPT/rdpz3) The value of the E-theta',' flag on the line:',//7x,'"',a32,'"',//7x,"isn't ",'"on" or "off", as is required.')

            stop
        end if
    end if

    ! LOOP HERE for each species block
100 continue

    ! Skip past the block terminator line.
    read (ndat0s,1000,end=990,err=995) uline

    if (uline(1:8) .ne. uterm(1:8)) then
        go to 100
    end if

    ! Check for the end of the present superblock.
    read (ndat0s,1030,end=990,err=995) unam1,unam2,unam3
1030 format(a24,2x,a24,2x,a24)

    if (unam1(1:8) .eq. uelemt(1:8)) then
        ! Found the header for the elements block.
        npx3t = npx3
        go to 999
    end if

    ! Have another species triplet.
    npx3 = npx3 + 1

    j2 = ilnobl(unam1)
    j3 = ilnobl(unam2)
    j4 = ilnobl(unam3)

    ! Check for blank names.
    if (j2.le.0 .or. j3.le.0 .or. j4.le.0) then
        write (ux8,'(i5)') npx3
        call lejust(ux8)
        j5 = ilnobl(ux8)

        if (j2 .le. 0) then
            unam1 = '<blank>'
            j2 = ilnobl(unam1)
        end if

        if (j3 .le. 0) then
            unam2 = '<blank>'
            j3 = ilnobl(unam2)
        end if

        if (j4 .le. 0) then
            unam3 = '<blank>'
            j4 = ilnobl(unam3)
        end if

        write (noutpt,1040) ux8(1:j5),unam1(1:j2),unam2(1:j3),unam3(1:j4)
        write (nttyo,1040) ux8(1:j5),unam1(1:j2),unam2(1:j3),unam3(1:j4)
1040 format(/' * Error - (EQPT/rdpz3) Have encountered a blank',' species name',/7x,'in the species triplet for block ',a,' of the superblock for',/7x,'Pitzer parameter data',' corresponding to species triplets.',/7x,'The offending',' triplet shows as ',a,', ',a,', ',a,'.')

        if (npx3 .gt. 1) then
            ux24a = utripl(1,npx3 - 1)
            ux24b = utripl(2,npx3 - 1)
            ux24c = utripl(3,npx3 - 1)

            if (ux24a(1:7).ne.'<blank>' .or.      ux24b(1:7).ne.'<blank>' .or.      ux24c(1:7).ne.'<blank>') then
                j6 = ilnobl(ux24a)
                j7 = ilnobl(ux24b)
                j8 = ilnobl(ux24c)
                write (noutpt,1050) ux24a(1:j6),ux24b(1:j7),ux24c(1:j8)
                write (nttyo,1050) ux24a(1:j6),ux24b(1:j7),ux24c(1:j8)
1050 format(7x,'This block follows the one for ',a,', ',a,', ',a,'.')
            end if
        end if

        nerr = nerr + 1
    end if

    ! Get the indices of the three species (n1, n2, n3). These
    ! indices point to the data read from the previously read
    ! species blocks.
    ! Get the index for the first species name.
    ! Calling sequence substitutions:
    !   n1 for na
    !   unam1 for unams
    call gspidx(ier,n1,nat,natmax,uaqsp,unam1)

    if (ier .gt. 0) then
        if (unam1(1:7).ne.'<blank>' .and.    unam2(1:7).ne.'<blank>' .and. unam3(1:7).ne.'<blank>') then
            write (noutpt,1100) unam1(1:j2),unam1(1:j2),unam2(1:j3),unam3(1:j4)
            write (nttyo,1100) unam1(1:j2),unam1(1:j2),unam2(1:j3),unam3(1:j4)
1100 format(/' * Error - (EQPT/rdpz3) The aqueous species ',a,' appearing in',/7x,'the Pitzer parameter data block',' for the species triplet',/7x,a,', ',a,', ',a,' does',' not appear in an aqueous',/7x,'species data block.')

            nerr = nerr + 1
        end if
    end if

    ! Get the index for the second species name.
    qdup12 = .false.

    if (unam2(1:24) .eq. unam1(1:24)) then
        qdup12 = .true.
        n2 = n1
    else
        ! Calling sequence substitutions:
        !   n2 for na
        !   unam2 for unams
        call gspidx(ier,n2,nat,natmax,uaqsp,unam2)

        if (ier .gt. 0) then
            if (unam1(1:7).ne.'<blank>' .and.      unam2(1:7).ne.'<blank>' .and. unam3(1:7).ne.'<blank>') then
                write (noutpt,1100) unam2(1:j3),unam1(1:j2),unam2(1:j3),unam3(1:j4)
                write (nttyo,1100) unam2(1:j3),unam1(1:j2),unam2(1:j3),unam3(1:j4)
                nerr = nerr + 1
            end if
        end if
    end if

    ! Get the index for the third species name.
    qdup13 = .false.
    qdup23 = .false.

    if (unam3(1:24) .eq. unam1(1:24)) then
        qdup13 = .true.

        if (qdup12) then
            qdup23 = .true.
        end if

        n3 = n1
    else if (unam3(1:24) .eq. unam2(1:24)) then
        qdup23 = .true.

        if (qdup12) then
            qdup13 = .true.
        end if

        n3 = n2
    else
        ! Calling sequence substitutions:
        !   n3 for na
        !   unam3 for unams
        call gspidx(ier,n3,nat,natmax,uaqsp,unam3)

        if (ier .gt. 0) then
            if (unam1(1:7).ne.'<blank>' .and.      unam2(1:7).ne.'<blank>' .and. unam3(1:7).ne.'<blank>') then
                write (noutpt,1100) unam3(1:j4),unam1(1:j2),unam2(1:j3),unam3(1:j4)
                write (nttyo,1100) unam3(1:j4),unam1(1:j2),unam2(1:j3),unam3(1:j4)
                nerr = nerr + 1
            end if
        end if
    end if

    if (n1.ne.0 .and. n2.ne.0 .and. n3.ne.0) then
        ! Get the charges.
        z1 = zaqsp(n1)
        z2 = zaqsp(n2)
        z3 = zaqsp(n3)

        iz1 = nint(z1)
        iz2 = nint(z2)
        iz3 = nint(z3)

        ! Count the neutral species.
        nn = 0

        if (iz1 .eq. 0) then
            nn = 1
        end if

        if (iz2 .eq. 0) then
            nn = nn + 1
        end if

        if (iz3 .eq. 0) then
            nn = nn + 1
        end if

        ! Count the cations.
        nc = 0

        if (iz1 .gt. 0) then
            nc = 1
        end if

        if (iz2 .gt. 0) then
            nc = nc + 1
        end if

        if (iz3 .gt. 0) then
            nc = nc + 1
        end if

        ! Count the anions.
        na = 0

        if (iz1 .lt. 0) then
            na = 1
        end if

        if (iz2 .lt. 0) then
            na = na + 1
        end if

        if (iz3 .lt. 0) then
            na = na + 1
        end if

        ! Check for illegal species combinations. These checks depend
        ! on the charge combinations; hence, they should not be made
        ! unless the all the charges have been determined.
        call tripck(na,nc,nerr,nn,noutpt,nttyo,n1,n2,n3,qdup12,qdup13,qdup23,unam1,unam2,unam3)
    end if

    ! Copy the names before rearranging for storage.
    u1 = unam1
    u2 = unam2
    u3 = unam3

    ! Rearrange for storage. The rules are:
    !   If exactly one neutral is present:
    !     neutral, cation, anion
    !   If exactly two neutrals are present:
    !     neutral 1, neutral 1, neutral2
    !   If two cations are present:
    !     cation 1, cation 2, anion
    !   If two anions are present:
    !     anion 1, anion 2, cation
    !   If two species are of the same charge type,
    !   store alphabetically.
    call artrip(iz1,iz2,iz3,na,nc,nn,n1,n2,n3,u1,u2,u3,z1,z2,z3)

    ! Store the triplet names in the order required by the
    ! storage rules. The triplet in the original order is
    ! preserved in unam1, unam2, unam3 for use in error messages
    ! and such.
    utripl(1,npx3) = u1
    utripl(2,npx3) = u2
    utripl(3,npx3) = u3

    ! Read the S-theta and psi parameters.
    read (ndat0s,1160,end=990,err=995) athetx(1),apsi(1,npx3)
1160 format(13x,f8.5,13x,f8.5)

    ! Read the first source.
    read (ndat0s,1170,end=990,err=995) uline
1170 format(13x,a24)

    ! Read the first and second temperature derivatives of the
    ! S-theta and psi parameters.
    read (ndat0s,1180,end=990,err=995) athetx(2),athetx(3)
    read (ndat0s,1180,end=990,err=995) apsi(2,npx3),apsi(3,npx3)
1180 format(13x,e10.3,13x,e10.3)

    ! Read the second source.
    read (ndat0s,1170,end=990,err=995) uline

    ! Check for illegal inputs in the S-theta data fields.
    call thetck(athetx,jpfcmx,na,nc,nerr,nn,noutpt,nttyo,n1,n2,n3,unam1,unam2,unam3)

    if ((nc.eq.2 .and. na.eq.1) .or.  (nc.eq.1 .and. na.eq.2)) then
        ! Find the theta index for the two cations or two anions.
        ifound = 0

        do nthd = 1,nthdt
            uthd1 = uthdtr(1,nthd)
            uthd2 = uthdtr(2,nthd)
            uthd3 = uthdtr(3,nthd)

            if (u1.eq.uthd1 .and. u2.eq.uthd2) then
                ifound = 1
                j2 = ilnobl(u1)
                j3 = ilnobl(u2)
                j4 = ilnobl(u3)
                j5 = ilnobl(uthd1)
                j6 = ilnobl(uthd2)
                j7 = ilnobl(uthd3)

                do j = 1,jpfcmx
                    write (ux8,'(i5)') j
                    call lejust(ux8)
                    j8 = ilnobl(ux8)

                    if (abs(atheta(j,nthd) - athetx(j)) .gt. 1.e-6) then
                        write (noutpt,1200) ux8(1:j8),u1(1:j2),u2(1:j3),u3(1:j4),athetx(j),uthd1(1:j5),uthd2(1:j6),uthd3(1:j7),atheta(j,nthd)
                        write (nttyo,1200) ux8(1:j8),u1(1:j2),u2(1:j3),u3(1:j4),athetx(j),uthd1(1:j5),uthd2(1:j6),uthd3(1:j7),atheta(j,nthd)
1200 format(/' * Warning - (EQPT/rdpz3) The S-theta',' a(',a,') coefficient for the',/7x,'triplet ',a,', ',a,', ',a,' is ',e12.5,',',/7x,'whereas the',' value previously read for the triplet',/7x,a,', ',a,', ',a,' is ',e12.5,'. The former',' value will be ignored.')

                        nwarn = nwarn + 1
                    end if
                end do

                go to 110
            end if
        end do

110 continue

        if (ifound .eq. 0) then
            ! Found a new pair defining a theta parameter.
            nthdt = nthdt + 1
            nthd = nthdt
            uthdtr(1,nthd) = u1
            uthdtr(2,nthd) = u2
            uthdtr(3,nthd) = u3

            do j = 1,jpfcmx
                atheta(j,nthd) = athetx(j)
            end do
        end if
    end if

    ! Go back to read the next species block.
    go to 100

990 continue
    write (noutpt,2000)
    write (nttyo,2000)
2000 format(/' * Error - (EQPT/rdpz3) Unexpectedly encountered',/7x,'end-of-file while reading the DATA0 file.')

    stop

995 continue
    write (noutpt,2010)
    write (nttyo,2010)
2010 format(/' * Error - (EQPT/rdpz3) Encountered a read format',/7x,'error while reading the DATA0 file.')

    stop

999 continue

    ! Deallocate work array.
    DEALLOCATE(athetx)
end subroutine rdpz3