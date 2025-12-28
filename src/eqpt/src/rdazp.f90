subroutine rdazp(azero,insgf,nazt,naztmx,ndat0s,nerr,noutpt,nttyo,uazp)
    !! This subroutine reads hard core diamaters and related parameters
    !! used in the B-dot equation from the DATA0 file. This data
    !! consists of aqueous species names, hard core diameter
    !! values, and an integer flag to indicate whether a species that
    !! is electrically neutral should be considered polar (set log
    !! gamma = 0) or non-polar (log gamma = log gamma for CO2(aq)).
    !! The data are returned in the uazp, azero, and insgf arrays. They
    !! are later written on the DATA1 and DATA1F files by EQPT/wrazp.f.
    !! This subroutine is called by:
    !!   EQPT/eqpt.f
    !! Principal input:
    !!   ndat0s = unit number of the stripped DATA0 file
    !! Principal output:
    !!   nazt   = the number of specified hard core diameters
    !!   uazp   = array of aqueous species names used to specify
    !!              hard core diamters on the data file
    !!   azero  = array of corresponding hard core diameters
    !!   insgf  = array of corresponding neutral species
    !!              activity coefficient flags
    implicit none

    ! Calling sequence variable declarations.
    integer :: naztmx

    integer :: ndat0s
    integer :: noutpt
    integer :: nttyo

    integer :: nazt
    integer :: nerr

    integer :: insgf(naztmx)

    character(len=24) :: uazp(naztmx)

    real(kind=8) :: azero(naztmx)

    ! Local variable declarations.
    integer :: j2
    integer :: j3
    integer :: j4

    integer :: ilnobl

    character(len=80) :: uline
    character(len=24) :: ux24
    character(len=8) :: ubdot
    character(len=8) :: uterm
    character(len=8) :: ux8

    data ubdot / 'bdot par' /
    data uterm / '+-------' /

    nazt = 0

    ! Advance to bdot header.
100 continue

    read (ndat0s,1000,end=990,err=995) uline
1000 format(a)

    ux8 = uline(1:8)

    if (ux8(1:8) .ne. ubdot(1:8)) then
        go to 100
    end if

    ! Skip the terminator line.
    read (ndat0s,1000,end=990,err=995) uline

    ! Read in the azero ('bdot') data.
110 continue

    read (ndat0s,1000,end=990,err=995) uline
    ux8 = uline(1:8)

    if (ux8(1:8) .eq. uterm(1:8)) then
        ! Have found the block terminator. Skip past the element
        ! header.
        read (ndat0s,1000,end=990,err=995) uline

        go to 999
    end if

    ! Found another entry. Read the data from the line.
    nazt = nazt + 1
    read(uline,1010,err=998) uazp(nazt),azero(nazt),insgf(nazt)
1010 format(a24,7x,f7.1,4x,i2)

    j2 = ilnobl(uazp(nazt))

    if (j2 .le. 0) then
        write (ux8,'(i5)') nazt
        call lejust(ux8)
        j3 = ilnobl(ux8)
        write (noutpt,1030) ux8(1:j3)
        write (nttyo,1030) ux8(1:j3)
1030 format(/' * Error - (EQPT/rdazp) Have encountered a blank',' species name on line ',a,/7x,'of the block of hard core',' diameter values for aqueous species.')

        if (nazt .gt. 1) then
            ux24 = uazp(nazt - 1)

            if (ux24(1:7) .ne. '<blank>') then
                j4 = ilnobl(ux24)
                write (noutpt,1040) ux24(1:j4)
                write (nttyo,1040) ux24(1:j4)
1040 format(7x,'This line follows the one for ',a,'.')
            end if
        end if

        uazp(nazt) = '<blank>'
        nerr = nerr + 1
    end if

    go to 110

990 continue
    write (noutpt,2000)
    write (nttyo,2000)
2000 format(/' * Error - (EQPT/rdazp) Unexpectedly encountered',/7x,'end-of-file while reading the DATA0 file.')

    write (noutpt,2010) nazt
    write (nttyo,2010) nazt
2010 format(7x,'The value of the local block line counter is ',i4,'.')

    j2 = ilnobl(uline)

    if (j2 .gt. 0) then
        j2 = min(j2,70)
        write (noutpt,2030) uline(1:j2)
        write (nttyo,2030) uline(1:j2)
2030 format(7x,'The last line read was the following:',/7x,'"',a,'"')
    end if

    stop

995 continue
    write (noutpt,2040)
    write (nttyo,2040)
2040 format(/' * Error - (EQPT/rdazp) Encountered a read format',/7x,'error while reading the DATA0 file.')

    write (noutpt,2010) nazt
    write (nttyo,2010) nazt
    j2 = ilnobl(uline)

    if (j2 .gt. 0) then
        j2 = min(j2,70)
        write (noutpt,2030) uline(1:j2)
        write (nttyo,2030) uline(1:j2)
    end if

    stop

998 continue
    write (noutpt,2050)
    write (nttyo,2050)
2050 format(/' * Error - (EQPT/rdazp) Encountered a read format error',/7x,'while reading data from lines read from the DATA0 file.')

    j2 = ilnobl(uline)

    if (j2 .gt. 0) then
        j2 = min(j2,70)
        write (noutpt,2060) uline(1:j2)
        write (nttyo,2060) uline(1:j2)
2060 format(7x,'The line with the problem is:',/7x,'"',a,'"')
    end if

    stop

999 continue
end subroutine rdazp