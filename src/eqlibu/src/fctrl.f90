real(kind=8) function fctrl(i)
    !! This subroutine returns the factorial function. The factorials
    !! from 0! to 10! are programmed into a data statement.
    !! This subroutine is called by:
    !!   Any
    !! Input:
    !!   i     = argument
    !! Output:
    !!   fctrl = factorial
    implicit none

    ! Calling sequence variable declarations.
    integer :: i

    ! Local variable declarations.
    integer :: ifac(0:10)
    integer :: ifactr
    integer :: j

    ! Set low factorials.
    data (ifac(j), j = 0,10) /1,1,2,6,24,120,720,5040,40320,362880,3628800/

    if ( i.ge.0 .and. i.le.10) then
        fctrl = ifac(i)
    else if (i.gt.10) then
        ifactr = ifac(10)

        do j = 10 + 1,i
            ifactr = ifactr*j
        end do

        fctrl = ifactr
    end if
end function fctrl