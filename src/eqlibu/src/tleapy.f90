subroutine tleapy(yy,qleapy)
    !! This subroutine determines if the year yy is a leap year. The
    !! logical flag qleapy is returned as true if this is the case.
    !! Modern era rules apply:
    !!   1. Any year not evenly divisible by 4 is not a leap year.
    !!   2. If a year is evenly divisible by 100 (i.e., is a "century
    !!      year"), it is not a leap year, unless it is also evenly
    !!      divisible by 400.
    !! Examples:
    !!    1. The year 1998 is not a leap year because it is not evenly
    !!       divisible by 4.
    !!    2. Century years like 1700, 1800, and 1900 are not leap years
    !!       because they are not evenly divisible by 400.
    !!    3. Century years like 1600 and 2000 are leap years because
    !!       they are evenly divisible by 400.
    !! The modern era rules (part of the Gregorian calendar) first
    !! applied in 1582. The first leap year of the Gregorian calendar
    !! was therefore 1584. This subroutine will not return any prior
    !! year as a leap year. Thus, 1580 is not a leap year.
    !! This subroutine is called by:
    !!   EQLIB/initim.f
    !!   EQLIB/runtim.f
    !! Input:
    !!   yy     = year, assumed to be an integer up to four digits
    !! Output:
    !!   qleapy = logical flag, true if yy is a leap year.
    implicit none

    ! Calling sequence variable declarations.
    integer :: yy

    logical :: qleapy

    ! Local variable declarations.
    ! None
    qleapy = .false.

    if (yy .ge. 1584) then
        ! Have a year sufficiently recent to possibly be a leap year.
        if (mod(yy,4) .eq. 0) then
            ! Have a year evenly divisible by 4. It might be a leap year.
            if (mod(yy,100) .eq. 0) then
                ! Have a century year. It might be a leap year.
                if (mod(yy,400) .eq. 0) then
                    ! Have a century year that is evenly divisible by 400.
                    ! This is a leap year.
                    qleapy = .true.
                end if
            else
                ! Have a year evenly divisible by 4, but it is not a
                ! century year. This is a leap year.
                qleapy = .true.
            end if
        end if
    end if
end subroutine tleapy