subroutine runtim(iexec0,jexec0,texec0,noutpt,nttyo,trun,tuser,tcpu,udate1,utime1)
    !! This subroutine gets the time and date at the end of execution,
    !! and also the run, user, and cpu times. The user and cpu times
    !! are not defined on all machines. In such cases, zero values are
    !! returned. In order to use this subroutine, EQLIBU/initim.f must
    !! be called at the start of execution to get values for iexec0,
    !! jexec0, and texec0. This routine should be called at the end of
    !! execution.
    !! This subroutine is called by:
    !!   EQPT/eqpt.f
    !!   EQ3NR/eq3nr.f
    !!   EQ6/eq6.f
    !! Input:
    !!   iexec0 = day of the year at the start of execution
    !!   jexec0 = the last two digits of the year at the start of
    !!              execution
    !!   texec0 = time at the start of execution, since midnight,
    !!            in seconds
    !!   noutpt = unit number of the output file
    !!   nttyo  = unit number of the screen file
    !! Output:
    !!   tcpu   = cpu time, in seconds
    !!   trun   = run time, in seconds
    !!   tuser  = user time, in seconds
    !!   udate1 = the date at the end of execution
    !!   utime1 = the time at the end of execution
    implicit none

    ! Calling sequence variable declarations.
    integer :: noutpt
    integer :: nttyo

    integer :: iexec0
    integer :: jexec0

    character(len=11) :: utime1
    character(len=9) :: udate1

    real(kind=8) :: texec0
    real(kind=8) :: trun
    real(kind=8) :: tuser
    real(kind=8) :: tcpu

    ! Local variable declarations.
    integer :: ndays(12)

    integer :: idays
    integer :: iexec1
    integer :: iy
    integer :: jexec1
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

    real(kind=8) :: texec1

    ! om   BEGIN_UNIX_DEPENDENT_CODE (Optional)
    ! om
    !        real*4 timear(2),tx
    ! om
    !        real*8 etime
    ! om
    !        external etime
    ! om
    ! om   END_UNIX_DEPENDENT_CODE
    data (umonth(mo), mo = 1,12) /'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'/

    data ndays(1) /31/
    data (ndays(mo), mo = 3,12) /31,30,31,30,31,31,30,31,30,31/

    ! Get the time and date.
    ! Calling sequence substitutions:
    !   udate1 for udate
    !   utime1 for utime
    call timdat(udate1,utime1)

    ! Compute the data needed to calculate the actual run time.
    ux = udate1(1:2)
    read (ux,1000,err=900) dd
1000 format(i2)

    um = udate1(3:5)

    do mo = 1,12
        if (umonth(mo)(1:3) .eq. um(1:3)) then
            mm = mo
            go to 110
        end if
    end do

    j2 = ilnobl(um)
    j3 = ilnobl(udate1)
    write (noutpt,1010) um(1:j2),udate1(1:j3)
    write (nttyo,1010) um(1:j2),udate1(1:j3)
1010 format(/' * Error - (EQLIBU/runtim) The string "',a,'"'," doesn't match the standard",/7x,'abbreviation for any',' month. This was extracted from the string "',a,'",',/7x,'which should contain the date at the end of',' execution. Check',/7x,'porting modifications in',' EQLIBU/timdat.f.')

    stop

110 continue
    if (udate1(8:9) .ne. '  ') then
        ux = udate1(6:9)
        read (ux,1020,err=900) yy
1020 format(i4)
    else
        ux = udate1(6:7)
        read (ux,1000,err=900) yy
    end if

    ux = utime1(1:2)
    read (ux,1000,err=910) hr

    ux = utime1(4:5)
    read (ux,1000,err=910) mi

    ux = utime1(7:8)
    read (ux,1000,err=910) se

    ux = utime1(10:11)
    read (ux,1000,err=910) hs

    texec1 = hr*3600. + mi*60. + se*1. + hs*0.01

    ! Get the number of days in February of the current year.
    ndays(2) = 28
    call tleapy(yy,qleapy)

    if (qleapy) then
        ndays(2) = 29
    end if

    iexec1 = dd

    if (mm .gt. 1) then
        do mo = 1, mm - 1
            iexec1 = iexec1 + ndays(mo)
        end do
    end if

    jexec1 = yy

    idays = iexec1 - iexec0

    if (jexec1 .gt. jexec0) then
        do iy = jexec0, jexec1 - 1
            idays = idays + 365

            if (mod(iy,4) .eq. 4) then
                if (mod(iy,100) .ne. 0) then
                    idays = idays + 1
                end if
            end if
        end do
    end if

    trun = idays*86400. + (texec1 - texec0)

    ! Normally, actual values are not returned for the user and cpu
    ! times. These times are not defined on all systems. The user time
    ! is usually not useful. The cpu time gives an idea of how much
    ! time a code is actually carrying out calculations, as opposed to
    ! executing other functions such as I/O.
    tuser = 0.
    tcpu = 0.

    ! om   BEGIN_UNIX_DEPENDENT_CODE (Optional)
    ! om
    ! om     Activate this block only if you have a UNIX system and want
    ! om     to obtain the user and cpu times on the output files.
    ! om
    ! om     Note: the variable declarations in the corresponding block
    ! om     above must be activated in order to successfully activate
    ! om     the present block.
    ! om
    ! om     The subroutine ETIME called below is a UNIX subroutine that
    ! om     gets the elapsed user and cpu times, returning the former in
    ! om     TIMEAR(1), the latter in TIMEAR(2). The array TIMEAR is real*4.
    ! om     On some systems, the name of this subroutine may be "ETIME_"
    ! om     instead of "ETIME".
    ! om
    !        tx = etime(timear)
    !        tuser = timear(1)
    !        tcpu = timear(2)
    ! om
    ! om     Note: the following statement doesn't really do anything except
    ! om     cause the compiler not to complain that tx is not used.
    ! om
    !        timear(1) = tx
    ! om
    ! om   END_UNIX_DEPENDENT_CODE
    go to 999

900 continue
    j3 = ilnobl(udate1)
    write (noutpt,1030) udate1(1:j3)
    write (nttyo,1030) udate1(1:j3)
1030 format(/' * Error - (EQLIBU/runtim) The date at the end of',' execution (udate1)'/7x,'must be a string in the format "','ddMmmccyy" (e.g., 23Sep1997) or',/7x,'"ddmmyy", (e.g.,',' 23Sep97). The actual string "',a,'"'," isn't",/7x,'in',' compliance with this requirement. Check porting',' modifications.',/7x,'in EQLIBU/timdat.f.')

    stop

910 continue
    j3 = ilnobl(utime1)
    write (noutpt,1040) utime1(1:j3)
    write (nttyo,1040) utime1(1:j3)
1040 format(/' * Error - (EQLIBU/runtim) The time at the end of',' execution (utime1)'/7x,'must be a string in the format "','hh:mm:ss:hs" (e.g., 15:32:19:28).',7x,'The actual string "',a,'"'," isn't in compliance with this requirement.",/7x, 'Check porting modifications in EQLIBU/timdat.f.')

    stop

999 continue
end subroutine runtim