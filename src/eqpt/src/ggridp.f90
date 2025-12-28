subroutine ggridp(ipch_asv,ipcv_asv,itgenf,jpdblo,jpfc_asv,jptffl,narx_asv,ndat0s,ndb_asv,noutpt,ntid_asv,ntpr_asv,nttyo,q500fl,uakey)
    !! This suboutine gets the necessary dimensioning parameters for the
    !! temperature grid on which thermodynamic data are represented.
    !! These parameters are presently embedded in the data file title,
    !! in order to maintain maximum forward and backward compatibility
    !! of data files between Versions 7 and 8 of EQ3/6. There is
    !! compatibility if the temperature grid matches the old "standard"
    !! form and if there are no data grids for volume and enthalpy data.
    !! This suboutine is called by:
    !!   EQPT/eqpt.f
    !! Principal input:
    !!   ndat0s = unit number of the stripped DATA0 file
    !!   uakey  = string identifying the type of aqueous species
    !!            activity coefficient model ("SEDH" or "Pitzer")
    !! Principal output:
    !!   ipch_asv = the maximum order of pressure corrections to the
    !!                enthalpy functions; ipch_asv = min(1,ipch)
    !!   ipcv_asv = the maximum order of pressure corrections to the
    !!   itgenf   = temperature grid enforcement parameter:
    !!                -1 = do nothing
    !!                 0 = warn if any sparsely filled ranges in
    !!                     "log K" temperature grid data
    !!                 1 = declare error condition if any sparsely
    !!                       filled ranges in "log K" temperature
    !!                       grid data
    !!   jpdblo   = Pitzer data block organization flag; -1 =
    !!                "classical", 0 = new
    !!                volume functions; ipcv_asv = min(1,ipcv)
    !!   jpfc_asv = the number of terms in the temperature function
    !!                used to represent Pitzer interaction parameters;
    !!                must be at least 1
    !!   jptffl   = integer flag denoting the temperature function
    !!                used to represent Pitzer interaction parameters;
    !!                -1 = classical truncated Taylor's series,
    !!                0 = LLNL5TERM, 1 = Greenberg and Moller (1989)
    !!                eight-term equation
    !!   narx_asv = the maximum number of points in any range of the
    !!                temperature grid
    !!   ntid_asv = the number of lines in the data file title
    !!   ntpr_asv = the number of temperature ranges on the temperature
    !!                grid (ntpr_asv = ntprt)
    !!   ndb_asv  = the maximum number of distinct points on the
    !!                temperature grid;
    !!                ndb_asv = ntpr_asv*(narx_asv - 1) + 1
    !!   q500fl   = logical flag; if .true., treat instances of "500."
    !!                in "log K" temperature grids as a "no data"
    !!                condition"
    implicit none

    ! Calling sequence variable declarations.
    integer :: ndat0s
    integer :: noutpt
    integer :: nttyo

    integer :: ipch_asv
    integer :: ipcv_asv
    integer :: jpfc_asv
    integer :: narx_asv
    integer :: ndb_asv
    integer :: ntid_asv
    integer :: ntpr_asv

    integer :: itgenf
    integer :: jpdblo
    integer :: jptffl

    character(len=8) :: uakey

    ! Local variable declarations.
    integer :: i
    integer :: j
    integer :: j2
    integer :: k
    integer :: n
    integer :: narxti
    integer :: narx_x
    integer :: ntpr
    integer :: ntpri

    integer :: ilnobl

    logical :: q500fl
    logical :: qnarxr
    logical :: qrderr

    character(len=80), dimension(:), allocatable :: ux80ar

    character(len=80) :: ulbufa
    character(len=80) :: ulbufb
    character(len=80) :: uline
    character(len=16) :: ustr16
    character(len=16) :: uterm
    character(len=8) :: ustr

    data uterm  /'+---------------'/

    ! Determine the number of lines in the data file title.
    ! Count that first terminator line as the last line of the title.
    n = 0
100 continue
    read (ndat0s,1000,end=990,err=995) uline
1000 format(a)

    n = n + 1

    if (uline(1:16) .ne. uterm(1:16)) then
        go to 100
    end if

    ntid_asv = n
    rewind(ndat0s)

    ! Allocate an array to hold the lines in the title.
    ALLOCATE(ux80ar(ntid_asv))

    ! Read the title into the holding array.
    do n = 1,ntid_asv
        read (ndat0s,1000,end=990,err=995) ux80ar(n)
    end do

    rewind(ndat0s)

    ! Set parameters for the temperature grid used to represent log K
    ! values and standard state enthalpy and volume functions.
    ! Values corresponding to the "classical" EQ3/6 temperature grid
    ! will be used as defaults. These are intended to support the use
    ! of version 7 format data files by the version 8 and higher
    ! software.
    !   ntprt       = the number of temperature ranges
    !   narxt(ntpr) = the number of points or coefficients in the
    !                   ntpr-th temperature range
    ! The old standard grid corresponds to:
    !   ntprt = 2
    !   narxt(1) = 4  (0, 25, 60 100 C)
    !   narxt(2) = 5  (100, 150, 200, 250, 300 C)
    !   ipch = -1   (no enthalpy function grids are present)
    !   ipcv = -1   (no volume function grids are present)
    ! Note that the last point of the ntpr-th temperature range
    ! is also the first point of the following range, if any.
    ! This insures continuity.
    ! Hence the matching default dimensioning is as follows. Note
    ! that if ipch < 1, ipch_asv must be set to one, as an array
    ! can't have a null or negative dimension. The same applies to
    ! ipcv and ipcv_asv.
    ntpr_asv = 2
    narx_asv = 5
    ipch_asv = 1
    ipcv_asv = 1

    ! Get the number of temperature ranges (ntpr_asv).
    do n = 1,ntid_asv
        ulbufa = ' '
        ulbufb = ' '
        j = 0

        ! Check for a keystring.
        i = index(ux80ar(n),'NO. OF TEMPERATURE RANGES')

        if (i .gt. 0) then
            j = i + 25
        else
            i = index(ux80ar(n),'NUMBER OF TEMPERATURE RANGES')

            if (i .gt. 0) then
                j = i + 28
            else
                i = index(ux80ar(n),'NTPRT')

                if (i .gt. 0) then
                    j = i + 5
                end if
            end if
        end if

        if (j .gt. 0) then
            ! Check for an equal sign following the keystring.
            ulbufa = ux80ar(n)(j:80)
            call lejust(ulbufa)
            i = index(ulbufa,'=')
            j = 2

            if (i .ne. 1) then
                j = 0
            end if
        end if

        if (j .gt. 0) then
            ! Extract a number matching the keystring.
            ulbufb = ulbufa(2:80)
            call lejust(ulbufb)
            ustr = ulbufb(1:8)
            call chrint(ntpr_asv,nttyo,qrderr,ustr)

            if (qrderr) then
                go to 997
            end if

            go to 200
        end if
    end do

200 continue

    ! Get the maximum number of points in any temperature range
    ! (narx_asv). Note that this is done by finding the number of
    ! points in each range and determining the maximum. However, an
    ! array to hold the actual values has not yet been allocated.
    ! The actual values will be redetermined and stored by rdwttl.f.
    qnarxr = .false.
    narx_x = 0
    narxti = 0

    do ntpr = 1,ntpr_asv
        do n = 1,ntid_asv
            ulbufa = ' '
            ulbufb = ' '
            j = 0

            ! Check for a keystring.
            i = index(ux80ar(n),'NO. OF POINTS IN RANGE')

            if (i .gt. 0) then
                j = i + 22
            else
                i = index(ux80ar(n),'NUMBER OF POINTS IN RANGE')

                if (i .gt. 0) then
                    j = i + 25
                else
                    i = index(ux80ar(n),'NARXT')

                    if (i .gt. 0) then
                        j = i + 5
                    end if
                end if
            end if

            if (j .gt. 0) then
                ! Extract a number acting as a subscript to the keystring.
                ulbufb = ux80ar(n)(j:80)
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

                qnarxr = .true.
                narx_x = max(narx_x,narxti)
                go to 220
            end if
        end do

220 continue
    end do

    ! If a value was read for any temperature range, use it. But
    ! do not use a number less than 1, as an array can not have a
    ! null or negative dimension.
    if (qnarxr) then
        narx_asv = narx_x
    end if

    narx_asv = max(1,narx_asv)

    ! Calculate the maximum possible number of points on the temperature
    ! grid for the current grid parameters.
    ndb_asv = ntpr_asv*(narx_asv - 1) + 1

    ! Set flag for dealing with sparsely filled ranges in the
    ! "log K" temperature grids.
    !   itgenf:
    !     -1 = do nothing
    !      0 = generate warnings for sparsely filled ranges
    !      1 = generate errors for sparsely filled ranges
    itgenf = -1

    do n = 1,ntid_asv
        ulbufa = ' '
        j = 0

        ! Check for a keystring.
        i = index(ux80ar(n),'SPARSE GRID RANGE CONDITION')

        if (i .gt. 0) then
            j = i + 27
        end if

        if (j .gt. 0) then
            ! Check for an equal sign following the keystring.
            ulbufa = ux80ar(n)(j:80)
            call lejust(ulbufa)
            i = index(ulbufa,'=')
            j = 2

            if (i .ne. 1) then
                j = 0
            end if
        end if

        if (j .gt. 0) then
            ! Extract a string input matching the keystring.
            ustr16 = ulbufa(2:17)
            call lejust(ustr16)
            k = index(ustr16,'IGNORE')

            if (k .gt. 0) then
                itgenf = -1
            else
                k = index(ustr16,'WARN')

                if (k .gt. 0) then
                    itgenf = 0
                else
                    k = index(ustr16,'ERROR')

                    if (k .gt. 0) then
                        itgenf = 1
                    end if
                end if
            end if

            go to 300
        end if
    end do

300 continue

    ! Set flag for dealing with values of "500." on the "log K"
    ! temperature grids.
    !   q500fl:
    !     .false. = do nothing
    !      .true. = treat values of 500. as "no data"
    q500fl = .true.

    do n = 1,ntid_asv
        ulbufa = ' '
        j = 0

        ! Check for a keystring.
        i = index(ux80ar(n),'INTERPRET 500 AS NO DATA')

        if (i .gt. 0) then
            j = i + 24
        end if

        if (j .gt. 0) then
            ! Check for an equal sign following the keystring.
            ulbufa = ux80ar(n)(j:80)
            call lejust(ulbufa)
            i = index(ulbufa,'=')
            j = 2

            if (i .ne. 1) then
                j = 0
            end if
        end if

        if (j .gt. 0) then
            ! Extract a string input matching the keystring.
            ustr16 = ulbufa(2:17)
            call lejust(ustr16)
            k = index(ustr16,'NO')

            if (k .le. 0) then
                k = index(ustr16,'FALSE')
            end if

            if (k .gt. 0) then
                q500fl = .false.
            else
                k = index(ustr16,'YES')

                if (k .le. 0) then
                    k = index(ustr16,'TRUE')
                end if

                if (k .gt. 0) then
                    q500fl = .true.
                end if
            end if

            go to 310
        end if
    end do

310 continue

    ! Get the dimension corresponding to the order of the pressure
    ! dependence of the enthalpy functions. An order of -1 means that
    ! No enthalpy functions data are present on the data file. An order
    ! of 0 means that enthalpy functions data are present, but no data
    ! for the pressure dependence of the same are present. However,
    ! the dimensioning variable ipch_asv must be set to a value of
    ! at least one, as an array can not have a null or negative
    ! dimension. Thus, it is necessary here to check only for data
    ! corresponding to order 2 or higher.
    do n = 1,ntid_asv
        ulbufa = ' '
        j = 0

        ! Check for a keystring.
        i = index(ux80ar(n),'d2H/dP2')

        if (i .gt. 0) then
            j = i + 7
        end if

        if (j .gt. 0) then
            ! Check for an equal sign following the keystring.
            ulbufa = ux80ar(n)(j:80)
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
                ipch_asv = 2
            end if

            go to 320
        end if
    end do

320 continue

    ! Get the dimension corresponding to the order of the pressure
    ! dependence of the volume functions. This process and the dimension
    ! itself are analogous to those described above for the pressure
    ! dependence of the enthalpy functions.
    do n = 1,ntid_asv
        ulbufa = ' '
        j = 0

        ! Check for a keystring.
        i = index(ux80ar(n),'d2V/dP2')

        if (i .gt. 0) then
            j = i + 7
        end if

        if (j .gt. 0) then
            ! Check for an equal sign following the keystring.
            ulbufa = ux80ar(n)(j:80)
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
                ipcv_asv = 2
            end if

            go to 350
        end if
    end do

350 continue

    ! Set flag for Pitzer data block organization.
    !   jpdblo:
    !     = -1 Classical: two superblocks, one for pure aqueous
    !          electrolytes, one for mixtures of two such electrolytes
    !          having a common ion. A given theta parameter can appear
    !          multiple times. Data for parameters involving electrically
    !          neutral species are folded into this structure. Checks are
    !          in place to sort things out.
    !     =  0 New: multiple superblocks, each corresponding to a
    !          charge-based type of species pair or triplet. This
    !          structure is more rational. Theta parameters for example
    !          have their own superblock.
    jpdblo = -1

    if (uakey(1:8) .eq. 'Pitzer  ') then
        do n = 1,ntid_asv
            ulbufa = ' '
            j = 0

            ! Check for a keystring.
            i = index(ux80ar(n),'PITZER DATA BLOCK ORG.')

            if (i .gt. 0) then
                j = i + 22
            end if

            if (j .gt. 0) then
                ! Check for an equal sign following the keystring.
                ulbufa = ux80ar(n)(j:80)
                call lejust(ulbufa)
                i = index(ulbufa,'=')
                j = 2

                if (i .ne. 1) then
                    j = 0
                end if
            end if

            if (j .gt. 0) then
                ! Extract a string input matching the keystring.
                ustr16 = ulbufa(2:17)
                call lejust(ustr16)
                k = index(ustr16,'NEW')

                if (k .le. 0) then
                    k = index(ustr16,'BYSPECIES')
                end if

                if (k .gt. 0) then
                    jpdblo = 0
                end if

                go to 400
            end if
        end do

400 continue
    end if

    ! Get the identification of the temperature function used to
    ! represent Pitzer interaction parameters. Note that the number
    ! of terms or coefficients (jpfcmx) is a dimensioning parameter
    ! that may also be set by an option embedded in the data file
    ! title.
    ! Temperature function flag.
    !   jptffl:
    !     = -1 Classic: a Taylor's series centered at 25C and
    !          truncated at second order.
    !            a1 = value at 25C
    !            a2 = first derivative at 25C
    !            a3 = second derivative at 25C
    !          This is currently required if the Pitzer data block
    !          organization flag is set to "Classical" (jpdblo = -1).
    !          The number of terms used (jpfcmx) is fixed at 3.
    !     =  0 The function:
    !            x(T) = a1 + a2/(T - Tr) + a3*ln(T/Tr)
    !                        a4(T - Tr) + a5(T**2 - Tr**2)
    !          This is also 25C centric (Tr = 25C).
    !          The number of terms used (jpfcmx) may vary from
    !          1 to 5. This always refers to the first jpfcmx terms.
    !     =  1 The Greenberg-Moller (1988) combination temperature
    !          function. This is not 25C centric. The number of terms
    !          used (mx) is fixed at 8.
    jptffl = -1

    if (uakey(1:8) .eq. 'Pitzer  ') then
        do n = 1,ntid_asv
            ulbufa = ' '
            j = 0

            ! Check for a keystring.
            i = index(ux80ar(n),'PITZER TEMP FUNCTION')

            if (i .gt. 0) then
                j = i + 20
            end if

            if (j .gt. 0) then
                ! Check for an equal sign following the keystring.
                ulbufa = ux80ar(n)(j:80)
                call lejust(ulbufa)
                i = index(ulbufa,'=')
                j = 2

                if (i .ne. 1) then
                    j = 0
                end if
            end if

            if (j .gt. 0) then
                ! Extract a string input matching the keystring.
                ustr16 = ulbufa(2:17)
                call lejust(ustr16)
                k = index(ustr16,'NEW')

                if (k .le. 0) then
                    k = index(ustr16,'Livermore')
                end if

                if (k .le. 0) then
                    k = index(ustr16,'LIVERMORE')
                end if

                if (k .le. 0) then
                    k = index(ustr16,'LLNL')
                end if

                if (k .le. 0) then
                    k = index(ustr16,'L5TERM')
                end if

                if (k .le. 0) then
                    k = index(ustr16,'5TERM')
                end if

                if (k .gt. 0) then
                    jptffl = 0
                    go to 410
                end if

                k = index(ustr16,'TEQUIL')

                if (k .le. 0) then
                    k = index(ustr16,'T8TERM')
                end if

                if (k .le. 0) then
                    k = index(ustr16,'GM8TERM')
                end if

                if (k .le. 0) then
                    k = index(ustr16,'8TERM')
                end if

                if (k .gt. 0) then
                    jptffl = 1
                    go to 410
                end if
            end if
        end do

410 continue
    end if

    ! Get the number of terms or coefficients (jpfc_asv) in the
    ! temperature function used to represent Pitzer parameters.
    jpfc_asv = 1

    if (uakey(1:8) .eq. 'Pitzer  ') then
        ! Three is the default value for a Pitzer model.
        jpfc_asv = 3

        do n = 1,ntid_asv
            ulbufa = ' '
            j = 0

            ! Check for a keystring.
            i = index(ux80ar(n),'NO. OF PITZER TEMP FUNC TERMS')

            if (i .gt. 0) then
                j = i + 29
            end if

            if (j .gt. 0) then
                ! Check for an equal sign following the keystring.
                ulbufa = ux80ar(n)(j:80)
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
                read (ustr,'(i8)') jpfc_asv
                go to 420
            end if
        end do

420 continue
    end if

    ! Deallocate the array holding the data file title.
    DEALLOCATE(ux80ar)
    go to 999

    ! Write messages for errors associated with reading the DATA0
    ! file.
990 continue
    write (noutpt,2000)
    write (nttyo,2000)
2000 format(/' * Error - (EQPT/ggridp) Unexpectedly encountered',/7x,'end-of-file while reading the title of the DATA0 file.')

    if (n .le. 1) then
        write (noutpt,2010)
        write (nttyo,2010)
2010 format(/7x,'This occurred while trying to read the first line',/7x,'of the title.')
    else
        ulbufa = ux80ar(n - 1)
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
2030 format(/' * Error - (EQPT/ggridp) Encountered a read format',/7x,'error while reading the DATA0 file.')

    if (n .le. 1) then
        write (noutpt,2010)
        write (nttyo,2010)
    else
        ulbufa = ux80ar(n - 1)
        call lejust(ulbufa)
        j2 = ilnobl(ulbufa)
        j2 = min(j2,70)
        write (noutpt,2020) ulbufa(1:j2)
    end if

    stop

997 continue
    ulbufa = ux80ar(n)
    call lejust(ulbufa)
    j2 = ilnobl(ulbufa)
    j2 = min(j2,70)
    write (noutpt,2050) ulbufa(1:j2)
    write (nttyo,2050) ulbufa(1:j2)
2050 format(/' * Error - (EQPT/ggridp) Encountered a read format',/7x,'error while reading data embedded in the title of the',/7x,'DATA0 file. This occurred while attempting to process',' the line:',/7x,'"',a,'"')

    stop

999 continue
end subroutine ggridp