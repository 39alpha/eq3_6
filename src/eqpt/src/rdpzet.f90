subroutine rdpzet(apsi,jpfcmx,nat,natmax,ndat0s,nerr,noutpt,npxzet,npx3mx,npx3t,nttyo,nwarn,uaqsp,utripl,zaqsp)
    !! This subroutine reads from the DATA1 file the coefficients
    !! required to compute those Pitzer interaction parameters
    !! associated with neutral-cation-anion (nca) triplets. This
    !! particular set of interaction coefficients includes the
    !! following:
    !!   1. zeta(nca)
    !! Coefficients required to compute Pitzer interaction coefficients
    !! associated with other triplet types or with pairs are not read
    !! by this subroutine. This subroutine is one of several that
    !! replace rdpz2.f and rdpz3.f.
    !! This suboutine is called by:
    !!   EQPT/eqpt.f
    !! Principal input:
    !!   ndat0s = unit number of the stripped DATA0 file
    !!   ndat1f = unit number of the DATA1F file
    !! Principal output:
    !!   apsi   = array of coefficients for computing psi parameters
    !!   nerr   = cumulative error counter
    !!   npxzet = the number of nca triplets for which Pitzer parameters
    !!             were read from the data file
    !!   npx3t  = the number of triplets of all species types for which
    !!              Pitzer parameters were read from the data file.
    !!   nwarn  = cumulative warning counter
    !!   uaqsp  = array of names of aqueous species
    !!   utripl = array of species names for species triplets
    !!   uaqsp  = array of electrical charge numbers of aqueous species
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
    integer :: npxzet
    integer :: npx3t
    integer :: nwarn

    character(len=24) :: uaqsp(natmax)
    character(len=24) :: utripl(3,npx3mx)

    real(kind=8) :: apsi(jpfcmx,npx3mx)
    real(kind=8) :: zaqsp(natmax)

    ! Local variable declarations.
    integer :: ier
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
    integer :: j9
    integer :: na
    integer :: nc
    integer :: nn
    integer :: npx3
    integer :: npx3te
    integer :: n1
    integer :: n2
    integer :: n3

    integer :: ilnobl

    logical :: qdupan
    logical :: qdup12
    logical :: qdup13
    logical :: qdup23

    character(len=80) :: udastr
    character(len=80) :: uline
    character(len=80) :: ux80
    character(len=24) :: unam1
    character(len=24) :: unam2
    character(len=24) :: unam3
    character(len=24) :: ux24a
    character(len=24) :: ux24b
    character(len=24) :: ux24c
    character(len=24) :: u1
    character(len=24) :: u2
    character(len=24) :: u3
    character(len=24) :: uhdzet
    character(len=24) :: uhdnxt
    character(len=16) :: ustr16
    character(len=8) :: uterm
    character(len=8) :: ux8

    real(kind=8) :: var
    real(kind=8) :: z1
    real(kind=8) :: z2
    real(kind=8) :: z3

    data uterm  /'+-------'/
    data uhdzet /"nca combinations: zeta(n"/
    data uhdnxt /"nnn' combinations: mu(nn"/

    npx3te = npx3t
    npx3 = npx3te
    j2 = 0
    j3 = 0
    j4 = 0

    ! Advance to the start of the superblock for these parameters.
100 continue
    read (ndat0s,1000,end=990,err=995) uline
1000 format(a80)

    if (uline(1:24) .ne. uhdzet(1:24)) then
        go to 100
    end if

    ! Skip the terminator.
    read (ndat0s,1000,end=990,err=995) uline

    ! Get the names of the species comprising a triplet.
110 continue
    read (ndat0s,1010,end=990,err=995) unam1,unam2,unam3
1010 format(a24,2x,a24,2x,a24)

    if (unam1(1:24) .eq. uhdnxt(1:24)) then
        ! Found the header for the next superblock of Pitzer parameters.
        ! Back up one line and terminate reading the current superblock.
        npxzet = npx3 - npx3te
        npx3t = npx3
        backspace(ndat0s)
        go to 999
    end if

    npx3 = npx3 + 1

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

    j4 = ilnobl(unam3)

    if (j4 .le. 0) then
        unam3 = '<blank>'
        j4 = ilnobl(unam3)
    end if

    write (ux8,'(i5)') npx3
    call lejust(ux8)
    j5 = ilnobl(ux8)

    ! Check for blank names.
    if (j2.le.0 .or. j3.le.0 .or. j4.le.0) then
        write (noutpt,1020) ux8(1:j5),unam1(1:j2),unam2(1:j3),unam3(1:j3)
        write (nttyo,1020) ux8(1:j5),unam1(1:j2),unam2(1:j3),unam3(1:j3)
1020 format(/' * Error - (EQPT/rdpzet) Have encountered a blank',' species name',/7x,'in the species triplet for block ',a,' of the superblock for',/7x,"Pitzer nca parameter data.",' The offending triplet shows as',/7x,a,', ',a,', ',a,'.')

        if (npx3 .gt. (npx3te + 1)) then
            ux24a = utripl(1,npx3 - 1)
            ux24b = utripl(2,npx3 - 1)
            ux24c = utripl(3,npx3 - 1)

            if (ux24a(1:7).ne.'<blank>' .or.      ux24b(1:7).ne.'<blank>' .or. ux24c(1:7).ne.'<blank>') then
                j7 = ilnobl(ux24a)
                j8 = ilnobl(ux24b)
                j9 = ilnobl(ux24c)
                write (noutpt,1030) ux24a(1:j7),ux24b(1:j8),ux24c(1:j9)
                write (nttyo,1030) ux24a(1:j7),ux24b(1:j8),ux24c(1:j9)
1030 format(7x,'This block follows the one for ',a,', ',a,', ',a,'.')
            end if
        end if

        nerr = nerr + 1
    end if

    ! XXXXXXXXX
    !      Get the indices of the three species (n1, n2, n3). These
    !      indices point to the data read from the previously read
    !      species blocks.
    !      Calling sequence substitutions:
    !        n1 for na
    !        unam1 for unams
    call gspidx(ier,n1,nat,natmax,uaqsp,unam1)

    if (ier .gt. 0) then
        if (unam1(1:7).ne.'<blank>' .and.    unam2(1:7).ne.'<blank>' .and. unam3(1:7).ne.'<blank>') then
            write (noutpt,1100) unam1(1:j2),unam1(1:j2),unam2(1:j3),unam3(1:j4)
            write (nttyo,1100) unam1(1:j2),unam1(1:j2),unam2(1:j3),unam3(1:j4)
1100 format(/' * Error - (EQPT/rdpzet) The aqueous species ',a,' appearing in',/7x,"the Pitzer nca parameter data",' block for the species triplet',/7x,a,', ',a,', ',a,' does not appear in an',/7x,'aqueous species data block.')

            nerr = nerr + 1
        end if
    end if

    ! Calling sequence substitutions:
    !   n2 for na
    !   unam2 for unams
    call gspidx(ier,n2,nat,natmax,uaqsp,unam2)

    if (ier .gt. 0) then
        if (unam1(1:7).ne.'<blank>' .and.    unam2(1:7).ne.'<blank>' .and. unam3(1:7).ne.'<blank>') then
            write (noutpt,1100) unam2(1:j3),unam1(1:j2),unam2(1:j3),unam3(1:j4)
            write (nttyo,1100) unam2(1:j3),unam1(1:j2),unam2(1:j3),unam3(1:j4)
            nerr = nerr + 1
        end if
    end if

    ! Calling sequence substitutions:
    !   n3 for na
    !   unam3 for unams
    call gspidx(ier,n3,nat,natmax,uaqsp,unam3)

    if (ier .gt. 0) then
        if (unam1(1:7).ne.'<blank>' .and.    unam2(1:7).ne.'<blank>' .and. unam3(1:7).ne.'<blank>') then
            write (noutpt,1100) unam3(1:j3),unam1(1:j2),unam2(1:j3),unam3(1:j4)
            write (nttyo,1100) unam3(1:j3),unam1(1:j2),unam2(1:j3),unam3(1:j4)
            nerr = nerr + 1
        end if
    end if

    ! Check for illegal appearance of solvent water.
    if (unam1(1:4) .eq.'H2O ' .or. unam2(1:4) .eq.'H2O ' .or.  unam3(1:4).eq.'H2O ' .or. n1.eq.1 .or. n2.eq.1 .or.  n3.eq.1) then
        write (noutpt,1200) unam1(1:j2),unam2(1:j3),unam3(1:j4)
        write (nttyo,1200) unam1(1:j2),unam2(1:j3),unam3(1:j4)
1200 format(/' * Error - (EQPT/rdpzet) Have found an illegal data',' block for the',/7x,'species triplet ',a,', ',a,', ',a,' in the superblock',/7x,"of Pitzer data for nca triplets.",' Solvent water may not',/7x,'appear in such a triplet.')

        nerr = nerr + 1
    end if

    ! Get the electrical charges of the species.
    z1 = 0.
    z2 = 0.
    z3 = 0.
    iz1 = 0
    iz2 = 0
    iz3 = 0

    if (n1 .gt. 0) then
        z1 = zaqsp(n1)
        iz1 = nint(z1)
    end if

    if (n2 .gt. 0) then
        z2 = zaqsp(n2)
        iz2 = nint(z2)
    end if

    if (n3 .gt. 0) then
        z3 = zaqsp(n3)
        iz3 = nint(z3)
    end if

    if (n1.ne.0 .and. n2.ne.0 .and. n3.ne.0) then
        ! Make sure that the species triplet is appropriate for the
        ! current superblock. It must be of type nca. Other triplet
        ! types such as cc'a, aa'c, nnn, nnn', and nn'n' are not
        ! permitted here.
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

        ! Check for duplications.
        qdup12 = unam1(1:24) .eq. unam2(1:24)
        qdup13 = unam1(1:24) .eq. unam3(1:24)
        qdup23 = unam2(1:24) .eq. unam3(1:24)
        qdupan = qdup12 .or. qdup23 .or. qdup13

        if (nc.eq.3 .or. na.eq.3    .or. (nc.eq.2 .and. nn.eq.1)    .or. (na.eq.2 .and. nn.eq.1)    .or. (nn.eq.2 .and. nc.eq.1)    .or. (nn.eq.2 .and. na.eq.1)    .or. (nc.eq.2 .and. na.eq.1 .and. qdupan)    .or. (na.eq.2 .and. nc.eq.1 .and. qdupan)) then
            ! The current data block refers to a triplet combination
            ! that is not valid in the Pitzer framework.
            write (noutpt,1220) unam1(1:j2),unam2(1:j3),unam3(1:j4)
            write (nttyo,1220) unam1(1:j2),unam2(1:j3),unam3(1:j4)
1220 format(/' * Error - (EQPT/rdpzet Have found an',' illegal data block for the',/7x,'species triplet ',a,', ',a,', ',a,' in the',/7x,'superblock of Pitzer'," data for nca triplets.",/7x,'The present triplet',' type is not valid in the Pitzer framework.')

            nerr = nerr + 1
        else if (.not.(nn.eq.1 .and. nc.eq.1 .and. na.eq.1)) then
            ! The current data block is in the wrong superblock.
            write (noutpt,1230) unam1(1:j2),unam2(1:j3),unam3(1:j4)
            write (nttyo,1230) unam1(1:j2),unam2(1:j3),unam3(1:j4)
1230 format(/' * Error - (EQPT/rdpzet Have found an',' illegal data block for the',/7x,'species triplet ',a,', ',a,', ',a,' in the',/7x,'superblock of Pitzer'," data nca triplets. The present",/7x,'triplet type',' may not appear in the current superblock.')

            nerr = nerr + 1
        end if
    end if

    ! Species triplet storage rules:
    !   nnn : neutral, neutral, neutral
    !   nnn': neutral 1, neutral1, neutral 2
    !   nca : neutral, cation, anion
    !   cc'a: cation 1, cation 2, anion
    !   aa'c: anion 1, anion 2, cation
    !   A duplicated species appears first.
    !   Otherwise, if two species are of the same charge type,
    !   store alphabetically.
    ! Here: neutral, cation, anion.
    ! Copy the names before rearranging for storage.
    u1 = unam1
    u2 = unam2
    u3 = unam3

    ! Rearrange according to the storage rules. Note that
    ! n1, n2, n3, z1, z2, z3, iz1, iz2, iz3 are all changed
    ! in addtion to u1, u2, u3. Note that unam1, unam2, and
    ! unam3 are not changed.
    call artrip(iz1,iz2,iz3,na,nc,nn,n1,n2,n3,u1,u2,u3,z1,z2,z3)

    ! Store the triplet names in the order required by the
    ! storage rules. The triplet in the original order is
    ! preserved in unam1, unam2, unam3 for use in error messages
    ! and such.
    utripl(1,npx3) = u1
    utripl(2,npx3) = u2
    utripl(3,npx3) = u3

    ! Read the coefficients for the zeta(nca) parameters. Store the
    ! zeta data as psi data.
    ! Check the required parameter header.
    read (ndat0s,1000,end=990,err=995) uline
    j5 = ilnobl(uline)
    j5 = min(j5,70)
    ux80 = uline
    call lejust(ux80)
    ustr16 = 'zeta'
    j6 = 4

    if (ux80(1:j6) .ne. ustr16(1:j6)) then
        write (noutpt,1330) uline(1:j5),ustr16(1:j6),unam1(1:j2),unam2(1:j3),unam3(1:j4)
        write (nttyo,1330) uline(1:j5),ustr16(1:j6),unam1(1:j2),unam2(1:j3),unam3(1:j4)
1330 format(/' * Error - (EQPT/rdpzet) Have found a line',' starting with:',/7x,'"',a,'"',/7x,'where "',a,'" was',' expected in the block for the species',/7x,'triplet ',a,', ',a,', ',a,' in the superblock of Pitzer data',/7x,"for nca triplets.")

        nerr = nerr + 1
        go to 999
    end if

    ! Read the coefficients for the associated temperature function.
    do j = 1,jpfcmx
        read (ndat0s,1000,end=990,err=995) uline
        j5 = ilnobl(uline)
        j5 = min(j5,70)
        ux80 = uline
        call lejust(ux80)
        ustr16 = 'a  ='
        j6 = 4
        write (ustr16(2:2),'(i1)') j

        if (ux80(1:j6) .eq. ustr16(1:j6)) then
            udastr = ux80
            udastr(1:j6) = ' '
            call g1dat(ier,noutpt,nttyo,udastr,var)

            if (ier .gt. 0) then
                write (noutpt,1320) uline(1:j5),unam1(1:j2),unam2(1:j3),unam3(1:j4)
                write (nttyo,1320) uline(1:j5),unam1(1:j2),unam2(1:j3),unam3(1:j4)
1320 format(/' * Error - (EQPT/rdpzet) Have found a line',' starting with:',/7x,'"',a,'"',/7x,'containing an',' expected numerical field that could not be read.',/7x,'This occurred in the block for the species triplet',/7x,a,', ',a,', ',a,' in the superblock of Pitzer data',/7x,"for nca triplets.")

                nerr = nerr + 1
                go to 999
            end if

            apsi(j,npx3) = var
        else
            write (noutpt,1330) uline(1:j5),ustr16(1:j6),unam1(1:j2),unam2(1:j3),unam3(1:j4)
            write (nttyo,1330) uline(1:j5),ustr16(1:j6),unam1(1:j2),unam2(1:j3),unam3(1:j4)
            nerr = nerr + 1
            go to 999
        end if
    end do

    ! Skip past the block delimiter line.
    read (ndat0s,1000,end=990,err=995) uline

    if (uline(1:8) .ne. uterm(1:8)) then
        j5 = ilnobl(uline)
        j5 = min(j5,70)
        write (noutpt,1470) uline(1:j5),unam1(1:j2),unam2(1:j3),unam3(1:j4)
        write (nttyo,1470) uline(1:j5),unam1(1:j2),unam2(1:j3),unam3(1:j4)
1470 format(/' * Error - (EQPT/rdpzet) Have found a line starting',' with:',/7x,'"',a,'"',/7x,' in the data block for the species',' triplet ',a,', ',a,', ',a,/7x,'in the superblock of'," Pitzer data for nca triplets. Should have",/7x,'found the delimiter line marking the end of the block.')

        nerr = nerr + 1
        go to 999
    end if

    ! Process the next data block.
    go to 110

990 continue
    write (noutpt,2010)
    write (nttyo,2010)
2010 format(/' * Error - (EQPT/rdpzet) Unexpectedly encountered',/7x,"an end-of-file error while reading the nca (zeta)",/7x,'superblock of the DATA0 file.')

    stop

995 continue
    write (noutpt,2020)
    write (nttyo,2020)
2020 format(/' * Error - (EQPT/rdpzet) Encountered a read format',/7x,"error while reading the nca (zeta) superblock of",/7x,'the DATA0 file.')

    stop

999 continue
end subroutine rdpzet