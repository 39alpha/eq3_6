real(kind=8) function texp(x)
    !! This subroutine is the function 10**x. The argument is tested in
    !! order to avoid overflow.
    !! This subroutine is called by:
    !!   Any
    !! Input:
    !!   x      = argument
    !! Output:
    !!   texp   = 10**x (maximum returned value is 1.e+125)
    implicit none

    ! Calling sequence variable declarations.
    real(kind=8) :: x

    ! Local variable declarations.
    !   None
    ! Test for overflow. The recommended truncation limit on the
    ! argument x is 125. On most systems, the real*8 exponent
    ! range is at least 308, which provides sufficient cushion.
    ! In considering how much cushion is enough, remember that
    ! exponentiated numbers may be subsequently be squared (in which
    ! case 1.e+125 gives rise to 1.e+250) or otherwise raised to some
    ! positive power. If your machine has a smaller but still acceptable
    ! exponent range (no less than 100), you should adjust the the limit
    ! used below accordingly. For example, if your machine's exponent
    ! range is 100, try a truncation limit of about 1.e+40. Your
    ! machine's exponent range is actually computed by EQLIBU/flpars.f.
    if (x .ge. 125.) then
        !        Avoid overflow.
        ! om     BEGIN_MACHINE_DEPENDENT_CODE
        ! om
        ! om       Note: Some compilers may not allow use of the construction
        ! om       "1.e+125" in the following setting. They may infer from the
        ! om       "e" part that this is a REAL*4 number, and conclude that the
        ! om       exponent is out of range. The "d" construction is appropriate
        ! om       for any 32-bit machine. It shouldn't cause a problem on a
        ! om       64-bit machine, such as a CRAY. If it should, however, change
        ! om       to the "e" construction.
        ! om
        texp = 1.d+125

        ! xx       texp = 1.e+125
        ! om
        ! om     END_MACHINE_DEPENDENT_CODE
    else
        texp = 10.**x
    end if
end function texp