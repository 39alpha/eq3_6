subroutine pcrss(apx,bpx,iapxmx,ibpxmx,iktmax,issot,nbtmx1,ndata1,ndat0s,ndat1f,nerr,nmodwr,nmt,nmtmax,noutpt,nslist,nttyo,nxt,nxtmax,uminsp,ussosp,ussoph)
    !! This subroutine reads data on solid solutions from the stripped
    !! DATA0 file, processes this data, and writes the results
    !! on the DATA1 and DATA1F files. The counter "nerr" is incremented
    !! by one for each error encountered.
    !! This subroutine is called by:
    !!   EQPT/eqpt.f
    !! Principal input:
    !!   ndata1 = unit number of the DATA1 file
    !!   ndat0s = unit number of the stripped DATA0 file
    !!   ndat1f = unit number of the DATA1F file
    !!   nerr   = cumulative error counter
    !! Principal output:
    !!   apx    = array of interaction coefficients for computing
    !!              activity coefficients in solid solutions
    !!   bpx    = array of site-mixing parameters for computing
    !!              activity coefficients in solid solutions
    !!   nerr   = cumulative error counter
    implicit none

    ! Calling sequence variable declarations.
    integer :: iapxmx
    integer :: ibpxmx
    integer :: iktmax
    integer :: nbtmx1
    integer :: nmtmax
    integer :: nxtmax

    integer :: ndata1
    integer :: ndat0s
    integer :: ndat1f
    integer :: noutpt
    integer :: nslist
    integer :: nttyo
    integer :: nxt

    integer :: nerr
    integer :: nmodwr
    integer :: nmt

    integer :: issot(nxtmax)

    character(len=24) :: uminsp(nmtmax)
    character(len=24) :: ussosp(iktmax,nxtmax)
    character(len=24) :: ussoph(nxtmax)

    real(kind=8) :: apx(iapxmx)
    real(kind=8) :: bpx(ibpxmx)

    ! Local variable declarations.
    integer :: i
    integer :: iapxt
    integer :: ibpxt
    integer :: ii
    integer :: ikt
    integer :: j
    integer :: jsol
    integer :: j2
    integer :: j3
    integer :: j4
    integer :: j5
    integer :: k
    integer :: n
    integer :: nmodx

    integer :: ilnobl

    character(len=80) :: uline
    character(len=72) :: uterm
    character(len=72) :: utermc
    character(len=24) :: usolid
    character(len=24) :: ux24
    character(len=24) :: ux24a
    character(len=24) :: uend24
    character(len=24) :: ublk24
    character(len=24) :: unone
    character(len=8) :: uendit
    character(len=8) :: uref
    character(len=8) :: ux8
    character(len=8) :: ux8a

    data ublk24 /'                        '/
    data usolid /'solid solutions         '/
    data unone  /'none                    '/
    data uendit /'endit.  '/
    data uref  /'referenc'/

    uterm(1:48) = '+-----------------------------------------------'
    uterm(49:72) = '------------------------'
    utermc = uterm
    utermc(1:1) = '*'

    uend24(1:8) = uendit(1:8)
    uend24(9:24) = ublk24(9:24)
    nxt = 0
    ikt = 0
    nmodx = 0

    j2 =ilnobl(usolid)
    write (ndata1) usolid
    write (ndat1f,1000) usolid(1:j2)
1000 format(a)

    write (ndat1f,1000) utermc(1:72)
    write (noutpt,1010) usolid(1:j2)
    write (nttyo,1010) usolid(1:j2)
    write (nslist,1010) usolid(1:j2)
1010 format(//1x,a,/)

    ! Note: the main loop returns here.
100 continue
    read (ndat0s,1020,end=990,err=995) uline
1020 format(a)

    ux24 = uline(1:24)
    ux8 = uline(1:8)
    call locase(ux8)

    if (ux8(1:8) .eq. uref(1:8)) then
        ! Found the 'references' superblock header. Finish and exit here.
        if (nxt .le. 0) then
            ! No solid solutions were found on the data file.
            write (nttyo,1030) nxt,unone(1:4)
            write (nslist,1030) nxt,unone(1:4)
            write (noutpt,1030) nxt,unone(1:4)
1030 format(1x,i5,2x,a)
        end if

        if (nxt.gt.0 .and. nmodx.ne.1) then
            ! The last solid solution was not previously echoed
            ! to the screen. Do this now. The echo consists of the
            ! name of the solid solution and the names of its
            ! end-member components.
            j2 = ilnobl(ussoph(nxt))
            write (nttyo,1030) nxt,ussoph(nxt)(1:j2)

            i = 1
120 continue

            if (i .eq. ikt) then
                j2 = ilnobl(ussosp(i,nxt))
                write (nttyo,1040) ussosp(i,nxt)
1040 format(12x,a)

                i = i + 1
            else
                j2 = ilnobl(ussosp(i + 1,nxt))
                write (nttyo,1050) ussosp(i,nxt),ussosp(i + 1,nxt)(1:j2)
1050 format(12x,a24,2x,a)

                i = i + 2
            end if

            if (i .le. ikt) then
                go to 120
            end if
        end if

        write (ndata1) uend24,ublk24,ublk24
        j2 = ilnobl(uend24)
        write (ndat1f,1000) uend24(1:j2)
        write (ndat1f,1000) utermc(1:72)

        ! Skip the terminator line.
        read (ndat0s,1020,end=990,err=995) uline
        go to 999
    end if

    nxt = nxt + 1

    ! Check for a blank solid solution name.
    j2 = ilnobl(ux24)

    if (j2 .le. 0) then
        write (ux8a,'(i5)') nxt
        call lejust(ux8a)
        j3 = ilnobl(ux8a)
        j4 = ilnobl(usolid)
        write (noutpt,1080) ux8a(1:j3),usolid(1:j4)
        write (nttyo,1080) ux8a(1:j3),usolid(1:j4)
1080 format(/' * Error - (EQPT/pcrss) Have encountered a blank',' phase name',/7x,'for phase block ',a,' of the ',a,' superblock.')

        if (nxt .gt. 1) then
            ux24a = ussoph(nxt - 1)

            if (ux24a(1:7) .ne. '<blank>') then
                j5 = ilnobl(ux24a)
                write (noutpt,1090) ux24a(1:j5)
                write (nttyo,1090) ux24a(1:j5)
1090 format(7x,'This block follows the one for ',a,'.')
            end if
        end if

        ux24 = '<blank>'
        nerr = nerr + 1
    end if

    ussoph(nxt) = uline(1:24)

    nmodx = mod(nxt,nmodwr)

    if (nmodwr .eq. 1) then
        nmodx = 1
    end if

    j2 = ilnobl(ussoph(nxt))

    if (nmodx .eq. 1) then
        write (nttyo,1030) nxt,ussoph(nxt)(1:j2)
    end if

    write (nslist,1030) nxt,ussoph(nxt)(1:j2)

    ! Skip the 'entered by' and 'keys' lines.
    read (ndat0s,1020,end=990,err=995) uline
    read (ndat0s,1020,end=990,err=995) uline

    ! Read number of end members.
    read (ndat0s,1020,end=990,err=995) uline
    read (uline,1150,err=995) ikt
1150 format(i3)

    if (ikt .le. 0) then
        j2 = ilnobl(ussoph(nxt))
        write (noutpt,1160) ussoph(nxt)(1:j2)
        write (nttyo,1160) ussoph(nxt)(1:j2)
1160 format(/' * Error - (EQPT/pcrss) Solid solution "',a,'" is',' not',/7x,'comprised of any end-members.')

        nerr = nerr + 1
    end if

    if (ikt .eq. 1) then
        j2 = ilnobl(ussoph(nxt))
        write (noutpt,1170) ussoph(nxt)(1:j2)
        write (nttyo,1170) ussoph(nxt)(1:j2)
1170 format(/' * Error - (EQPT/pcrss) Solid solution "',a,'" is',' comprised',/7x,'of only one end-member. Two or more',' end-members are required.')

        nerr = nerr + 1
    end if

    if (ikt .gt. iktmax) then
        j2 = ilnobl(ussoph(nxt))
        write (noutpt,1180) ussoph(nxt)(1:j2),ikt,iktmax
        write (nttyo,1180) ussoph(nxt)(1:j2),ikt,iktmax
1180 format(/' * Error - (EQPT/pcrss) Solid solution ',a,/7x,'is composed of ',i3,' end-members, which exceeds the',/7x,'maximum dimension (iktpar) of ',i3,'.')

        nerr = nerr + 1
    end if

    issot(nxt) = ikt

    if (ikt .gt. 0) then
        n = 0

        do ii = 1,ikt/2
            read (ndat0s,1020,end=990,err=995) uline
            read (uline,1190,err=995) (ussosp(n + k,nxt), k = 1,2)
1190 format((13x,a24,13x,a24))

            n = n + 2
        end do

        j = mod(ikt,2)

        if (j .gt. 0) then
            read (ndat0s,1020,end=990,err=995) uline
            read (uline,1190,err=995) (ussosp(n + k,nxt), k = 1,j)
            n = n + j
        end if

        i = 1
130 continue

        if (i .eq. ikt) then
            j2 = ilnobl(ussosp(i,nxt))

            if (nmodx .eq. 1) then
                write (nttyo,1040) ussosp(i,nxt)(1:j2)
            end if

            write (nslist,1040) ussosp(i,nxt)(1:j2)
            i = i + 1
        else
            j2 = ilnobl(ussosp(i + 1,nxt))

            if (nmodx .eq. 1) then
                write (nttyo,1050) ussosp(i,nxt),ussosp(i + 1,nxt)(1:j2)
            end if

            write (nslist,1050) ussosp(i,nxt),ussosp(i + 1,nxt)(1:j2)
            i = i + 2
        end if

        if (i .le. ikt) then
            go to 130
        end if
    end if

    ! Check the end-member names for blanks.
    do i = 1,ikt
        ux24 = ussosp(i,nxt)
        j3 = ilnobl(ux24)

        if (j3 .le. 0) then
            j2 = ilnobl(ussoph(nxt))
            write (ux8,'(i5)') i
            call lejust(ux8)
            j4 = ilnobl(ux8)
            write (noutpt,1210) ussoph(nxt)(1:j2),ux8(1:j4)
            write (nttyo,1210) ussoph(nxt)(1:j2),ux8(1:j4)
1210 format(/' * Error - (EQPT/pcrss) Solid solution "',a,'"',/7x,'has a blank name for end-member number ',a,'.')

            if (i .gt. 1) then
                ux24a = ussosp(i - 1,nxt)

                if (ux24a(1:7) .ne. '<blank>') then
                    j5 = ilnobl(ux24a)
                    write (noutpt,1220) ux24a(1:j5)
                    write (nttyo,1220) ux24a(1:j5)
1220 format(7x,'This follows the end-member ',a,'.')
                end if
            end if

            ussosp(i,nxt) = '<blank>'
            nerr = nerr + 1
        end if
    end do

    ! Note: the end-members names are checked for duplicates by
    ! nxspck.f, and validated against pure mineral names by vxspck.f.
    ! These subroutines are called by eqpt.f.
    ! Read the activity coefficient code.
    read (ndat0s,1020,end=990,err=995) uline
    read (uline,1250,err=995) jsol
1250 format(10x,i1)

    call initaz(apx,iapxmx)
    call initaz(bpx,ibpxmx)

    ! Read the number of interaction coefficients for computing
    ! activity coefficients in this solid solution.
    read (ndat0s,1020,end=990,err=995) uline
    read (uline,1270,err=995) iapxt
1270 format(i3)

    if (iapxt .gt. 0) then
        n = 0

        do ii = 1,iapxt/6
            read (ndat0s,1020,end=990,err=995) uline
            read (uline,1280,err=995) (apx(n + k), k = 1,6)
1280 format(6f6.3)

            n = n + 6
        end do

        j = mod(ikt,6)

        if (j .gt. 0) then
            read (ndat0s,1020,end=990,err=995) uline
            read (uline,1280,err=995) (apx(n + k), k = 1,j)
            n = n + j
        end if
    end if

    ! Read the number of site-mixing parameters for computing
    ! activity coefficients in this solid solution.
    read (ndat0s,1270,end=990,err=995) ibpxt

    if (ibpxt .gt. 0) then
        n = 0

        do ii = 1,ibpxt/6
            read (ndat0s,1020,end=990,err=995) uline
            read (uline,1280,err=995) (bpx(n + k), k = 1,6)
            n = n + 6
        end do

        j = mod(ikt,6)

        if (j .gt. 0) then
            read (ndat0s,1020,end=990,err=995) uline
            read (uline,1280,err=995) (bpx(n + k), k = 1,j)
            n = n + j
        end if
    end if

    ! Write block on the DATA1 and DATA1F files.
    write (ndata1) ussoph(nxt),ikt,jsol
    write (ndata1) (ussosp(i,nxt), i = 1,ikt)
    write (ndata1) iapxt

    if (iapxt .gt. 0) then
        write (ndata1) (apx(i), i = 1,iapxt)
    end if

    write (ndata1) ibpxt

    if (ibpxt .gt. 0) then
        write (ndata1) (bpx(i), i = 1,ibpxt)
    end if

    write (ndat1f,1330) ussoph(nxt),ikt,jsol
1330 format(a24,2x,2i5)

    i = 1
250 continue

    if (i .eq. ikt) then
        j2 = ilnobl(ussosp(i,nxt))
        write (ndat1f,1340) ussosp(i,nxt)(1:j2)
1340 format(5x,a)

        i = i + 1
    else
        j2 = ilnobl(ussosp(i + 1,nxt))
        write (ndat1f,1350) ussosp(i,nxt),ussosp(i + 1,nxt)(1:j2)
1350 format(5x,a24,5x,a)

        i = i + 2
    end if

    if (i .le. ikt) then
        go to 250
    end if

    write (ndat1f,1360) iapxt
1360 format(i3)

    if (iapxt .gt. 0) then
        write (ndat1f,1370) (apx(i), i = 1,iapxt)
    end if

1370 format(6f6.3)

    write (ndat1f,1360) ibpxt

    if (ibpxt .gt. 0) then
        write (ndat1f,1370) (bpx(i), i = 1,ibpxt)
    end if

    j2 = ilnobl(ussoph(nxt))
    write (noutpt,1380) ussoph(nxt)(1:j2),ikt,jsol
1380 format(1x,a,/5x,'Number of components= ',i3,2x,'Activity coefficient code= ',i3)

    i = 1
260 continue

    if (i .eq. ikt) then
        j2 = ilnobl(ussosp(i,nxt))
        write (noutpt,1040) ussosp(i,nxt)(1:j2)
        i = i + 1
    else
        j2 = ilnobl(ussosp(i + 1,nxt))
        write (noutpt,1050) ussosp(i,nxt),ussosp(i + 1,nxt)(1:j2)
        i = i + 2
    end if

    if (i .le. ikt) then
        go to 260
    end if

    write (ndat1f,1000) utermc(1:72)

    ! Skip to the terminator.
170 continue

    read (ndat0s,1020,end=990,err=995) uline
    ux8 = uline(1:8)

    if(ux8(1:8) .ne. uterm(1:8)) then
        go to 170
    end if

    go to 100

990 continue
    write (noutpt,2000)
    write (nttyo,2000)
2000 format(/' * Error - (EQPT/pcrss) Unexpectedly encountered',/7x,'end-of-file while reading the DATA0 file.')

    write (noutpt,2010) nxt
    write (nttyo,2010) nxt
2010 format(7x,'The value of the solid solution block counter is ',i3,'.')

    if (nxt .gt. 0) then
        j2 = ilnobl(ussoph(nxt))
        write (noutpt,2020) ussoph(nxt)(1:j2)
        write (nttyo,2020) ussoph(nxt)(1:j2)
2020 format(7x,'The last solid solution name read was "',a,'".')
    else
        write (noutpt,2030) unone(1:j2)
        write (nttyo,2030) unone(1:j2)
2030 format(7x,'No solid solutions were read.')
    end if

    j2 = ilnobl(uline)

    if (j2 .gt. 0) then
        j2 = min(j2,70)
        write (noutpt,2040) uline(1:j2)
        write (nttyo,2040) uline(1:j2)
2040 format(7x,'The last line read was the following:',/7x,'"',a,'"')
    end if

    stop

995 continue
    write (noutpt,2050)
    write (nttyo,2050)
2050 format(/' * Error - (EQPT/pcrss) Encountered a read format',/7x,'error while reading the DATA0 file.')

    write (noutpt,2010) nxt
    write (nttyo,2010) nxt

    if (nxt .gt. 0) then
        j2 = ilnobl(ussoph(nxt))
        write (noutpt,2020) ussoph(nxt)(1:j2)
        write (nttyo,2020) ussoph(nxt)(1:j2)
    else
        write (noutpt,2030)
        write (nttyo,2030)
    end if

    j2 = ilnobl(uline)

    if (j2 .gt. 0) then
        j2 = min(j2,70)
        write (noutpt,2040) uline(1:j2)
        write (nttyo,2040) uline(1:j2)
    end if

    stop

999 continue
end subroutine pcrss