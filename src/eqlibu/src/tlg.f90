real(kind=8) function tlg(x)
    !! This subroutine is the function log10, except that the log10(0.)
    !! is returned with a value of -99999.
    !! This subroutine is called by:
    !!   Any
    !! Input:
    !!   x      = argument
    !! Output:
    !!   tlg    = log10 x (except that tlg(0.) = -99999.)
    implicit none

    ! Calling sequence variable declarations.
    real(kind=8) :: x

    ! Local variable declarations.
    !   None
    if (x .eq. 0.) then
        tlg = -99999.
    else
        tlg = log10(x)
    end if
end function tlg