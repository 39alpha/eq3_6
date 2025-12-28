subroutine gntpr(ntpr,ntprmx,ntprt,tempc,tempcu)
    !! This subroutine finds the value of ntpr, the temperature range
    !! flag.
    !! This subroutine is called by:
    !!   EQ3NR/eq3nr.f
    !!   EQ6/eq6.f
    !!   EQ6/tpadv.f
    !! Principal input:
    !!   ntprt  = the number of temperature ranges
    !!   tempc  = temperature, C
    !!   tempcu = array of temperatures, C, defining the upper boundaries
    !!              of the temperature ranges
    !! Principal output:
    !!   ntpr   = temperature range flag
    implicit none

    ! Calling sequence variable declarations.
    integer :: ntprmx

    integer :: ntpr
    integer :: ntprt

    real(kind=8) :: tempcu(ntprmx)

    real(kind=8) :: tempc

    ! Local variable declarations.
    ! None
    ! Determine the temperature range flag (ntpr).
    do ntpr = 1,ntprt
        if (tempc .le. tempcu(ntpr)) then
            go to 100
        end if
    end do

    ntpr = ntprt
100 continue
end subroutine gntpr