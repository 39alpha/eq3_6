subroutine initim(iexec0,jexec0,texec0,noutpt,nttyo,udate0,utime0)
    !! This subroutine gets the time and date of the start of execution.
    !! It also calculates the data in the form needed to compute the
    !! actual run time at the end of execution (see EQLIBU/runtim.f).
    !! This routine should be called at the start of execution.
    !! This subroutine is called by:
    !!   EQPT/eqpt.f
    !!   EQ3NR/eq3nr.f
    !!   EQ6/eq6.f
    !! Input:
    !!   noutpt = unit number of the output file
    !!   nttyo  = unit number of the screen file
    !! Output:
    !!   texec0 = time at the start of execution, since midnight,
    !!              in seconds
    !!   iexec0 = day of the year at the start of execution
    !!   jexec0 = the last two digits of the year at the start of
    !!              execution
    !!   udate0 = the date at the start of execution
    !!   utime0 = the time at the start of execution
    implicit none

    ! Calling sequence variable declarations.
    integer :: noutpt
    integer :: nttyo

    integer :: iexec0
    integer :: jexec0

    character(len=11) :: utime0
    character(len=9) :: udate0

    real(kind=8) :: texec0

    ! Local variable declarations.
    integer :: ndays(12)

    integer :: mo
    integer :: mm
    integer :: dd
    integer :: yy
    integer :: hr
    integer :: hs
    integer :: mi
    integer :: se
    integer :: j2
    integer :: j3

    integer :: ilnobl

    logical :: qleapy

    character(len=8) :: umonth(12)
    character(len=8) :: ux
    character(len=3) :: um

    data (umonth(mo), mo = 1,12) /'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'/

    data ndays(1) /31/
    data (ndays(mo), mo = 3,12) /31,30,31,30,31,31,30,31,30,31/

    ! Get the time and date.
    ! Calling sequence substitutions:
    !   udate0 for udate
    !   utime0 for utime
    call timdat(udate0,utime0)

    ! Compute the data needed to calculate the actual run time (this
    ! will be used by EQLIBU/runtim.f).
    ux = udate0(1:2)
    read (ux,1000,err=900) dd
1000 format(i2)

    um = udate0(3:5)

    do mo = 1,12
        if (umonth(mo)(1:3) .eq. um(1:3)) then
            mm = mo
            go to 110
        end if
    end do

    j2 = ilnobl(um)
    j3 = ilnobl(udate0)
    write (noutpt,1010) um(1:j2),udate0(1:j3)
    write (nttyo,1010) um(1:j2),udate0(1:j3)
1010 format(/' * Error - (EQLIBU/initim) The string "',a,'"'," doesn't match the standard",/7x,'abbreviation for any',' month. This was extracted from the string "',a,'",',/7x,'which should contain the date at the start of',' execution. Check',/7x,'porting modifications in',' EQLIBU/timdat.f.')

    stop

110 continue
    if (udate0(8:9) .ne. '  ') then
        ux = udate0(6:9)
        read (ux,1020,err=900) yy
1020 format(i4)
    else
        ux = udate0(6:7)
        read (ux,1000,err=900) yy
    end if

    ux = utime0(1:2)
    read (ux,1000,err=910) hr

    ux = utime0(4:5)
    read (ux,1000,err=910) mi

    ux = utime0(7:8)
    read (ux,1000,err=910) se

    ux = utime0(10:11)
    read (ux,1000,err=910) hs

    texec0 = hr*3600. + mi*60. + se*1. + hs*0.01

    ! Get the number of days in February of the current year.
    ndays(2) = 28
    call tleapy(yy,qleapy)

    if (qleapy) then
        ndays(2) = 29
    end if

    iexec0 = dd

    if (mm .gt. 1) then
        do mo = 1, mm - 1
            iexec0 = iexec0 + ndays(mo)
        end do
    end if

    jexec0 = yy

    go to 999

900 continue
    j3 = ilnobl(udate0)
    write (noutpt,1030) udate0(1:j3)
    write (nttyo,1030) udate0(1:j3)
1030 format(/' * Error - (EQLIBU/initim) The date at the start of',' execution (udate0)'/7x,'must be a string in the format "','ddMmmccyy" (e.g., 23Sep1997) or',/7x,'"ddmmyy", (e.g.,',' 23Sep97). The actual string "',a,'"'," isn't",/7x,'in',' compliance with this requirement. Check porting',' modifications.',/7x,'in EQLIBU/timdat.f.')

    stop

910 continue
    j3 = ilnobl(utime0)
    write (noutpt,1040) utime0(1:j3)
    write (nttyo,1040) utime0(1:j3)
1040 format(/' * Error - (EQLIBU/initim) The time at the start of',' execution (utime0)'/7x,'must be a string in the format "','hh:mm:ss:hs" (e.g., 15:32:19:28).',7x,'The actual string "',a,'"'," isn't in compliance with this requirement.",/7x, 'Check porting modifications in EQLIBU/timdat.f.')

    stop

999 continue
end subroutine initim