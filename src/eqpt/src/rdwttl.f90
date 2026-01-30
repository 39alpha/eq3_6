subroutine rdwttl(ipch,ipcv,jpdblo,jpfcmx,jptffl,narxt,ndata1,ndat0s,ndat1f,noutpt,nslist,ntitld,ntidmx,ntprmx,ntprt,nttyo,uakey,utitld)
    !! This suboutine reads the title on the DATA0 file and writes it
    !! on the DATA1 and DATA1F files. It checks the DATA0 file
    !! header for validity. It also searches the title for embedded
    !! flags and data indicating special treatment.
    !! Possible embedded flags and data deal with the following:
    !!   1. Defining the temperature grid used on the data file.
    !!      By default, the classic eight-temperature grid is used.
    !!   2. Indicating the presence of data grids for enthalpy
    !!      and volume functions. By default, these data grids
    !!      are assumed to be absent.
    !!   3. Defining the temperature function used to describe
    !!      Pitzer interaction parameters. By default, this function
    !!      is the classic second-order Taylor's series centered
    !!      at 25C.
    !! This suboutine is called by:
    !!   EQPT/eqpt.f
    !! Principal input:
    !!   ndat0s = unit number of the stripped DATA0 file
    !! Principal output:
    !!   ipch   = enthalpy functions data grid flag:
    !!              -1 = no enthalpy grids
    !!               0 = enthalpy grids present
    !!               1 = grids present for the enthalpy and its first
    !!                     pressure derivative
    !!               2 = grids present for the enthalpy and its first
    !!                     and second pressure derivatives
    !!   ipcv   = volume functions data flag:
    !!              -1 = no volume grids
    !!               0 = volume grids present
    !!               1 = grids present for the volume and its first
    !!                     pressure derivative
    !!               2 = grids present for the volume and its first
    !!                     and second pressure derivatives
    !!   ntitld = the number of lines in the title
    !!   utitld = array of title lines
    implicit none

    ! Calling sequence variable declarations.
    integer :: jpfcmx
    integer :: ntidmx
    integer :: ntprmx

    integer :: narxt(ntprmx)

    integer :: ndata1
    integer :: ndat0s
    integer :: ndat1f
    integer :: noutpt
    integer :: nslist
    integer :: ntitld
    integer :: nttyo

    integer :: ipch
    integer :: ipcv
    integer :: jpdblo
    integer :: jptffl
    integer :: ntprt

    character(len=80) :: utitld(ntidmx)
    character(len=8) :: uakey

    ! Local variable declarations.
    integer :: i
    integer :: j
    integer :: j2
    integer :: j3
    integer :: j4
    integer :: k
    integer :: n
    integer :: narxti
    integer :: nerr
    integer :: ntpr
    integer :: ntpri

    integer :: ilnobl

    logical :: qrderr

    character(len=56) :: ux56

    character(len=80) :: ulbufa
    character(len=80) :: ulbufb
    character(len=72) :: uterm
    character(len=72) :: utermc
    character(len=32) :: upchst
    character(len=32) :: upcvst
    character(len=8) :: ustr
    character(len=8) :: ux8a
    character(len=8) :: ux8b

    uterm(1:48) = '+-----------------------------------------------'
    uterm(49:72) = '------------------------'
    utermc = uterm
    utermc(1:1) = '*'

    ! Read the data file title.
    ntitld = ntidmx

    do n = 1,ntitld
        read (ndat0s,1000,end=990,err=995) utitld(n)
1000 format(a)

        if (utitld(n)(1:8) .eq. uterm(1:8)) then
            go to 110
        end if
    end do

110 continue

    ! Write title information.
    do n = 1,ntitld
        j2 = ilnobl(utitld(n))
        j2 = min(j2,79)
        write (ndata1) utitld(n)
        write (ndat1f,1000) utitld(n)(1:j2)
        write (nttyo,1140) utitld(n)(1:j2)
        write (noutpt,1140) utitld(n)(1:j2)
        write (nslist,1140) utitld(n)(1:j2)
1140 format(' ',a)
    end do

    ! Get the parameters for the temperature grid used for log K values
    ! and standard state enthalpy and volume functions.
    !   ntprt       = the number of temperature ranges
    !   narxt(ntpr) = the number of points or coefficients in the
    !                   ntpr-th temperature range
    ! The last point of the ntpr-th temperature range is also the first
    ! point of the following range, if any. This insures continuity.
    ! The default values correspond to the original EQ3/6 data grid:
    !   ntprt = 2
    !   narxt(1) = 4  (0, 25, 60 100 C)
    !   narxt(2) = 5  (100, 150, 200, 250, 300 C)
    ! Note: the number of temperature ranges (ntprt) is equal to
    ! the dimensioned limit (ntpr_asv), which has been determined
    ! by a prior scan. The prior scan has also provided narx_asv,
    ! the dimenionsed limit on the number of values in a temperature
    ! range. This was found as the greatest value of any range.
    ! Now, it is necessary to get the actual value for each range.
    ntprt = ntprmx
    call initiz(narxt,ntprmx)

    if (ntprt .eq. 2) then
        narxt(1) = 4
        narxt(2) = 5
    end if

    do ntpr = 1,ntprt
        do n = 1,ntitld
            ulbufa = ' '
            ulbufb = ' '
            j = 0

            ! Check for a keystring.
            i = index(utitld(n),'NO. OF POINTS IN RANGE')

            if (i .gt. 0) then
                j = i + 22
            else
                i = index(utitld(n),'NUMBER OF POINTS IN RANGE')

                if (i .gt. 0) then
                    j = i + 25
                else
                    i = index(utitld(n),'NARXT')

                    if (i .gt. 0) then
                        j = i + 5
                    end if
                end if
            end if

            if (j .gt. 0) then
                ! Extract a number acting as a subscript to the keystring.
                ulbufb = utitld(n)(j:80)
                call lejust(ulbufb)
                j = index(ulbufb,' ')
                j = min(j,9)
                k = j - 1

                if (k .gt. 0) then
                    ustr = ulbufb(1:k)
                    call chrint(ntpri,nttyo,qrderr,ustr)

                    if (qrderr) then
                        go to 997
                    end if

                    if (ntpri .ne. ntpr) then
                        j = 0
                    end if
                else
                    j = 0
                end if
            end if

            if (j .gt. 0) then
                ! Check for an equal sign following the keystring.
                ulbufa = ulbufb(j:80)
                call lejust(ulbufa)
                i = index(ulbufa,'=')
                j = 2

                if (i .ne. 1) then
                    j = 0
                end if
            end if

            if (j .gt. 0) then
                ! Extract a number matching the keystring.
                ulbufb = ulbufa(j:80)
                call lejust(ulbufb)
                ustr = ulbufb(1:8)
                call chrint(narxti,nttyo,qrderr,ustr)

                if (qrderr) then
                    go to 997
                end if

                narxt(ntpr) = narxti
                go to 220
            end if
        end do

220 continue
    end do

    nerr = 0
    ntpr = 1

    if (narxt(ntpr) .lt. 2) then
        write (ux8a,'(i5)') ntpr
        call lejust(ux8a)
        j3 = ilnobl(ux8a)
        write (ux8b,'(i5)') narxt(ntpr)
        call lejust(ux8b)
        j4 = ilnobl(ux8b)
        write (noutpt,1320) ux8a(1:j3),ux8b(1:j4)
        write (nttyo,1320) ux8a(1:j3),ux8b(1:j4)
1320 format(/' * Error - (EQPT/rdwttl) The number of points in',/7x,'temperature range ',a,' is ',a,'. The number of',/7x,'points in the first or last temperature range must be',/7x,'at least 2.')

        nerr = nerr + 1
    end if

    do ntpr = 2,ntprt - 1
        if (narxt(ntpr) .lt. 3) then
            write (ux8a,'(i5)') ntpr
            call lejust(ux8a)
            j3 = ilnobl(ux8a)
            write (ux8b,'(i5)') narxt(ntpr)
            call lejust(ux8b)
            j4 = ilnobl(ux8b)
            write (noutpt,1330) ux8a(1:j3),ux8b(1:j4)
            write (nttyo,1330) ux8a(1:j3),ux8b(1:j4)
1330 format(/' * Error - (EQPT/rdwttl) The number of points in',/7x,'temperature range ',a,' is ',a,'. The number of',/7x,'points in an interior temperature range must be at',/7x,'least 3.')

            nerr = nerr + 1
        end if
    end do

    ntpr = ntprt

    if (ntpr.gt.1 .and. narxt(ntpr).lt.2) then
        write (ux8a,'(i5)') ntpr
        call lejust(ux8a)
        j3 = ilnobl(ux8a)
        write (ux8b,'(i5)') narxt(ntpr)
        call lejust(ux8b)
        j4 = ilnobl(ux8b)
        write (noutpt,1320) ux8a(1:j3),ux8b(1:j4)
        write (nttyo,1320) ux8a(1:j3),ux8b(1:j4)
        nerr = nerr + 1
    end if

    if (nerr .gt. 0) then
        stop
    end if

    ! Write the grid parameters.
    ux56 = 'Number of ranges in the logK temperature grid'
    write (ndata1) ux56
    j2 = ilnobl(ux56)
    write (ndat1f,1350) ux56(1:j2)
1350 format(a)

    write (ndata1) ntprt
    write (ndat1f,1360) ntprt
1360 format(i5)

    do ntpr = 1,ntprt
        ux56 = 'Number of points in range '
        j2 = ilnobl(ux56)
        write (ux8a,'(i5)') ntpr
        call lejust(ux8a)
        j3 = ilnobl(ux8a)
        ux56(j2 + 2: j2 + 1 + j3) = ux8a(1:j3)

        write (ndata1) ux56
        j2 = ilnobl(ux56)
        write (ndat1f,1350) ux56(1:j2)

        n = narxt(ntpr)
        write (ndata1) n
        write (ndat1f,1360) n
    end do

    write (ux8a,'(i5)') ntprt
    call lejust(ux8a)
    j3 = ilnobl(ux8a)
    write (noutpt,1380) ux8a(1:j3)
    write (nttyo,1380) ux8a(1:j3)
1380 format(/' Number of logK temperature grid ranges= ',a)

    do ntpr = 1,ntprt
        write (ux8a,'(i5)') ntpr
        call lejust(ux8a)
        j2 = ilnobl(ux8a)
        n = narxt(ntpr)
        write (ux8b,'(i5)') n
        call lejust(ux8a)
        j3 = ilnobl(ux8a)
    end do

    ! Set flags for enthalpy and volume functions data grids.
    ipch = -1

    do n = 1,ntitld
        ulbufa = ' '
        j = 0

        ! Check for a keystring.
        i = index(utitld(n),'ENTHALPY')

        if (i .gt. 0) then
            j = i + 8
        end if

        if (j .gt. 0) then
            ! Check for an equal sign following the keystring.
            ulbufa = utitld(n)(j:80)
            call lejust(ulbufa)
            i = index(ulbufa,'=')
            j = 2

            if (i .ne. 1) then
                j = 0
            end if
        end if

        if (j .gt. 0) then
            ! Extract a string input matching the keystring.
            ustr = ulbufa(2:9)
            call lejust(ustr)
            k = index(ustr,'ON')

            if (k .le. 0) then
                k = index(ustr,'PRESENT')
            end if

            if (k .le. 0) then
                k = index(ustr,'ACTIVE')
            end if

            if (k .gt. 0) then
                ipch = 0
            end if

            go to 300
        end if
    end do

300 continue

    do n = 1,ntitld
        ulbufa = ' '
        j = 0

        ! Check for a keystring.
        i = index(utitld(n),'dH/dP')

        if (i .gt. 0) then
            j = i + 5
        end if

        if (j .gt. 0) then
            ! Check for an equal sign following the keystring.
            ulbufa = utitld(n)(j:80)
            call lejust(ulbufa)
            i = index(ulbufa,'=')
            j = 2

            if (i .ne. 1) then
                j = 0
            end if
        end if

        if (j .gt. 0) then
            ! Extract a string input matching the keystring.
            ustr = ulbufa(2:9)
            call lejust(ustr)
            k = index(ustr,'ON')

            if (k .le. 0) then
                k = index(ustr,'PRESENT')
            end if

            if (k .le. 0) then
                k = index(ustr,'ACTIVE')
            end if

            if (k .gt. 0) then
                if (ipch .eq. 0) then
                    ipch = 1
                else
                    write (noutpt,1400)
                    write (nttyo,1400)
1400 format(/" * Error - (EQPT/rdwttl) Can't use dH/dP data",' without',/7x,'enthalpy (H) data.')

                    stop
                end if
            end if

            go to 310
        end if
    end do

310 continue

    do n = 1,ntitld
        ulbufa = ' '
        j = 0

        ! Check for a keystring.
        i = index(utitld(n),'d2H/dP2')

        if (i .gt. 0) then
            j = i + 7
        end if

        if (j .gt. 0) then
            ! Check for an equal sign following the keystring.
            ulbufa = utitld(n)(j:80)
            call lejust(ulbufa)
            i = index(ulbufa,'=')
            j = 2

            if (i .ne. 1) then
                j = 0
            end if
        end if

        if (j .gt. 0) then
            ! Extract a string input matching the keystring.
            ustr = ulbufa(2:9)
            call lejust(ustr)
            k = index(ustr,'ON')

            if (k .le. 0) then
                k = index(ustr,'PRESENT')
            end if

            if (k .le. 0) then
                k = index(ustr,'ACTIVE')
            end if

            if (k .gt. 0) then
                if (ipch .eq. 1) then
                    ipch = 2
                else
                    write (noutpt,1410)
                    write (nttyo,1410)
1410 format(/" * Error - (EQPT/rdwttl) Can't use d2H/dP2",' data without',/7x,'enthalpy (H) and dH/dP data.')

                    stop
                end if
            end if

            go to 320
        end if
    end do

320 continue

    ipcv = -1

    do n = 1,ntitld
        ulbufa = ' '
        j = 0

        ! Check for a keystring.
        i = index(utitld(n),'VOLUME')

        if (i .gt. 0) then
            j = i + 6
        end if

        if (j .gt. 0) then
            ! Check for an equal sign following the keystring.
            ulbufa = utitld(n)(j:80)
            call lejust(ulbufa)
            i = index(ulbufa,'=')
            j = 2

            if (i .ne. 1) then
                j = 0
            end if
        end if

        if (j .gt. 0) then
            ! Extract a string input matching the keystring.
            ustr = ulbufa(2:9)
            call lejust(ustr)
            k = index(ustr,'ON')

            if (k .le. 0) then
                k = index(ustr,'PRESENT')
            end if

            if (k .le. 0) then
                k = index(ustr,'ACTIVE')
            end if

            if (k .gt. 0) then
                ipcv = 0
            end if

            go to 330
        end if
    end do

330 continue

    do n = 1,ntitld
        ulbufa = ' '
        j = 0

        ! Check for a keystring.
        i = index(utitld(n),'dV/dP')

        if (i .gt. 0) then
            j = i + 5
        end if

        if (j .gt. 0) then
            ! Check for an equal sign following the keystring.
            ulbufa = utitld(n)(j:80)
            call lejust(ulbufa)
            i = index(ulbufa,'=')
            j = 2

            if (i .ne. 1) then
                j = 0
            end if
        end if

        if (j .gt. 0) then
            ! Extract a string input matching the keystring.
            ustr = ulbufa(2:9)
            call lejust(ustr)
            k = index(ustr,'ON')

            if (k .le. 0) then
                k = index(ustr,'PRESENT')
            end if

            if (k .le. 0) then
                k = index(ustr,'ACTIVE')
            end if

            if (k .gt. 0) then
                if (ipcv .eq. 0) then
                    ipcv = 1
                else
                    write (noutpt,1420)
                    write (nttyo,1420)
1420 format(/" * Error - (EQPT/rdwttl) Can't use dV/dP data",' without',/7x,'volume (V) data.')

                    stop
                end if
            end if

            go to 340
        end if
    end do

340 continue

    do n = 1,ntitld
        ulbufa = ' '
        j = 0

        ! Check for a keystring.
        i = index(utitld(n),'d2V/dP2')

        if (i .gt. 0) then
            j = i + 7
        end if

        if (j .gt. 0) then
            ! Check for an equal sign following the keystring.
            ulbufa = utitld(n)(j:80)
            call lejust(ulbufa)
            i = index(ulbufa,'=')
            j = 2

            if (i .ne. 1) then
                j = 0
            end if
        end if

        if (j .gt. 0) then
            ! Extract a string input matching the keystring.
            ustr = ulbufa(2:9)
            call lejust(ustr)
            k = index(ustr,'ON')

            if (k .le. 0) then
                k = index(ustr,'PRESENT')
            end if

            if (k .le. 0) then
                k = index(ustr,'ACTIVE')
            end if

            if (k .gt. 0) then
                if (ipcv .eq. 1) then
                    ipcv = 2
                else
                    write (noutpt,1430)
                    write (nttyo,1430)
1430 format(/" * Error - (EQPT/rdwttl) Can't use d2V/dP2",' data without',/7x,'volume (V) and dV/dP data.')

                    nerr = nerr + 1
                end if
            end if

            go to 350
        end if
    end do

350 continue

    if (ipch.gt.0 .and. ipcv.lt.0) then
        write (noutpt,1440)
        write (nttyo,1440)
1440 format(/" * Error - (EQPT/rdwttl) Can't have pressure",' corrections',/7x,'for enthalpy functions without having such corrections',/7x,'for log K functions. At a minimum, volume functions',/7x,'must be present.')

        nerr = nerr + 1
    end if

    ! Write flags for enthalpy and volume functions grids.
    ux56 = 'Enthalpy functions flag'
    write (ndata1) ux56
    j2 = ilnobl(ux56)
    write (ndat1f,1450) ux56(1:j2)
1450 format(a)

    write (ndata1) ipch
    write (ndat1f,1460) ipch
1460 format(i5)

    upchst = 'Not present'

    if (ipch .eq. 0) then
        upchst = 'Present, no dH/dP data'
    else if (ipch .eq. 1) then
        upchst = 'Present, with dH/dP data'
    else if (ipch .ge. 2) then
        upchst = 'Present, with dH/dP-dnH/dPn data'
        write (upchst,'(21x,i1,4x,i1)') ipch
    end if

    j2 = ilnobl(upchst)

    write (noutpt,1470) ipch,upchst(1:j2)
    write (nttyo,1470) ipch,upchst(1:j2)
1470 format(/' Enthalpy functions flag= ',i2,' (',a,')')

    ux56 = 'Volume functions flag'
    write (ndata1) ux56
    j2 = ilnobl(ux56)
    write (ndat1f,1450) ux56(1:j2)

    write (ndata1) ipcv
    write (ndat1f,1460) ipcv

    upcvst = 'Not present'

    if (ipcv .eq. 0) then
        upcvst = 'Present, no dV/dP data'
    else if (ipcv .eq. 1) then
        upcvst = 'Present, with dV/dP data'
    else if (ipcv .ge. 2) then
        upcvst = 'Present, with dV/dP-dnV/dPn data'
        write (upcvst,'(21x,i1,4x,i1)') ipcv
    end if

    j2 = ilnobl(upcvst)

    write (noutpt,1480) ipcv,upcvst(1:j2)
    write (nttyo,1480) ipcv,upcvst(1:j2)
1480 format(' Volume functions flag  = ',i2,' (',a,')')

    if (uakey(1:8) .eq. 'Pitzer  ') then
        ! Write the flag for Pitzer data block organization.
        ux56 = 'Pitzer data block organization flag'
        write (ndata1) ux56
        j2 = ilnobl(ux56)
        write (ndat1f,1450) ux56(1:j2)

        write (ndata1) jpdblo
        write (ndat1f,1460) jpdblo
    end if

    ! Check the temperature function used to represent Pitzer
    ! interaction parameters. This is determined by the value
    ! of jptffl, which was read previously by ggridp.f.
    if (uakey(1:8) .eq. 'Pitzer  ') then
        if (jpdblo .eq. -1) then
            if (jptffl .ne. -1) then
                write (noutpt,1510)
                write (nttyo,1510)
1510 format(/' * Note - (EQPT/rdwttl) Resetting the jptffl',' Pitzer data',/7x,'temperature function flag to -1'," (classic 25C-centric Taylor's series",/7x,'truncated',' at second order), because the Pitzer data block',/7x,'organization flag (jpdblo) is set to -1',' ("Classical").',/7x,'Other Pitzer data temperature',' functions are not permitted.')

                jptffl = -1
            end if
        end if

        if (jptffl .eq. -1) then
            if (jpfcmx .ne. 3) then
                write (noutpt,1520)
                write (nttyo,1520)
1520 format(/' * Error - (EQPT/rdwttl) The number of terms',' (jpfcmx) used in',/7x,'the Pitzer data temperature',' function is not 3, as is required',/7x,'in the case'," of the classical 25C-centric Taylor's series",/7x,'truncated at second order (jptffl = -1).')

                nerr = nerr + 1
            end if
        end if

        if (jptffl .eq. 1) then
            if (jpfcmx .ne. 8) then
                write (noutpt,1530)
                write (nttyo,1530)
1530 format(/' * Error - (EQPT/rdwttl) The number of terms',' (jpfcmx) used in',/7x,'the Pitzer data temperature',' function is not 8, as is required',/7x,'in the case',' of the case of the Greenberg-Moller combination',/7x,'temperature function (jptffl = 1).')

                nerr = nerr + 1
            end if
        end if

        if (jptffl .eq. 0) then
            if (jpfcmx.lt.1 .or. jpfcmx.gt.5) then
                write (noutpt,1540)
                write (nttyo,1540)
1540 format(/' * Error - (EQPT/rdwttl) The number of terms',' terms (jpfcmx)',/7x,'used in the LLNL maximal five-term',' temperature equation (jptffl = 0)',/7x,'is not in the',' allowed range of of 1-5. ')

                nerr = nerr + 1
            end if
        end if
    end if

    if (uakey(1:8) .eq. 'Pitzer  ') then
        ! Write the flag for the temperature function used to represent
        ! Pitzer interaction coefficients.
        ux56 = 'Pitzer parameter temperature function'
        write (ndata1) ux56
        j2 = ilnobl(ux56)
        write (ndat1f,1450) ux56(1:j2)

        write (ndata1) jptffl
        write (ndat1f,1460) jptffl
    end if

    write (ndat1f,1220) utermc(1:72)
1220 format(a)

    if (nerr .gt. 0) then
        stop
    end if

    go to 999

990 continue
    write (noutpt,2000)
    write (nttyo,2000)
2000 format(/' * Error - (EQPT/rdwttl) Unexpectedly encountered',/7x,'end-of-file while reading the title of the DATA0 file.')

    if (n .le. 1) then
        write (noutpt,2010)
        write (nttyo,2010)
2010 format(/7x,'This occurred while trying to read the first line',/7x,'of the title.')
    else
        ulbufa = utitld(n - 1)
        call lejust(ulbufa)
        j2 = ilnobl(ulbufa)
        j2 = min(j2,70)
        write (noutpt,2020) ulbufa(1:j2)
        write (nttyo,2020) ulbufa(1:j2)
2020 format(/7x,'This occurred while trying to read the line',' following',/7x,'the line:',/7x,'"',a,'"')
    end if

    stop

995 continue
    write (noutpt,2030)
    write (nttyo,2030)
2030 format(/' * Error - (EQPT/rdwttl) Encountered a read format',/7x,'error while reading the DATA0 file.')

    if (n .le. 1) then
        write (noutpt,2010)
        write (nttyo,2010)
    else
        ulbufa = utitld(n - 1)
        call lejust(ulbufa)
        j2 = ilnobl(ulbufa)
        j2 = min(j2,70)
        write (noutpt,2020) ulbufa(1:j2)
    end if

    stop

997 continue
    ulbufa = utitld(n)
    call lejust(ulbufa)
    j2 = ilnobl(ulbufa)
    j2 = min(j2,70)
    write (noutpt,2050) ulbufa(1:j2)
    write (nttyo,2050) ulbufa(1:j2)
2050 format(/' * Error - (EQPT/rdwttl) Encountered a read format',/7x,'error while reading data embedded in the title of the',/7x,'DATA0 file. This occurred while attempting to process',' the line:',/7x,'"',a,'"')

    stop

999 continue
end subroutine rdwttl
