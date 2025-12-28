subroutine timdat(udate,utime)
    !! This subroutine gets the date and clock time in ASCII for time and
    !! date stamping on output files. The variable "udate" is returned in
    !! the form "ddMmmccyy" (e.g., 23Sep1997) under Fortran 90, in the
    !! form 'ddMmmyy" (e.g., 23Sep97) under Fortran 77. The variable
    !! "utime" is returned in the form "hh:mm:ss:hs", though in some
    !! cases the hundred-th second part ("hs") may be blank.
    !! See also EQLIBU/initim.f and EQLIBU/runtim.f.
    !! This subroutine is called by:
    !!   EQLIBU/initim.f
    !!   EQLIBU/runtim.f
    !! Input:
    !!   None
    !! Output:
    !!   udate  = string containing the date
    !!   utime  = string containing the clock time
    implicit none

    ! Calling sequence variable declarations.
    character(len=11) :: utime
    character(len=9) :: udate

    integer :: values(8)

    character(len=8) :: date
    character(len=10) :: time
    character(len=6) :: zone
    character(len=8) :: ux
    character(len=8) :: umonth(12)

    integer :: mo
    integer :: mm

    data (umonth(mo), mo = 1,12) /'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'/

    ! Note: date_and_time is a Fortran 90 intrinsic subroutine.
    call date_and_time(date,time,zone,values)

    ! Note: date has the form ccyymmdd; udate has the form ddMmmccyy.
    udate(1:2) = date(7:8)
    ux = date(5:6)
    read (ux,'(i5)') mm
    udate(3:5) = umonth(mm)(1:3)
    udate(6:9) = date(1:4)

    ! Note: time has the form hhmmss.sss (sss = milliseconds);
    ! utime has the form hh:mm:ss:hs.
    utime(1:2) = time(1:2)
    utime(3:3) = ':'
    utime(4:5) = time(3:4)
    utime(6:6) = ':'
    utime(7:8) = time(5:6)
    utime(9:9) = ':'
    utime(10:11) = time(8:9)

    ! The following is nonsense intended to keep the compiler from
    ! complaining that zone and values are not used.
    ux = zone
    mm = values(1)

999 continue
end subroutine timdat