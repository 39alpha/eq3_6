real(kind=8) function arrmxn(array,nmax)
    !! This subroutine finds the max norm of the real*8 array "array".
    !! Note that array could actually have more than one dimension. For
    !! example, if the array is really "a", which has dimensions of i1
    !! and i2, this subroutine could be used by calling it in the
    !! following manner: call arrmxn(a,(i1*i2)). Use this subroutine
    !! with caution if nmax is not the product of the true (declared)
    !! dimensions.
    !! To find the index of the first element corresponding to the max
    !! norm, call EQLIBU/iarmxn.f instead. Then just use the index to
    !! get the max norm if this is also desired.
    !! This subroutine is called by:
    !!   Any
    !! Input:
    !!   array  = the array
    !!   nmax   = the dimension or pseudo-dimension of the array
    !! Output:
    !!   arrmxn = the max norm of the array
    implicit none

    ! Calling sequence variable declarations.
    integer :: nmax

    real(kind=8) :: array(nmax)

    ! Local variable declarations.
    integer :: i
    integer :: ileft

    real(kind=8) :: aax
    real(kind=8) :: ax

    arrmxn = 0.

    if (nmax .le. 0) then
        go to 999
    end if

    ! Note the use of a local variable (ax) within the loop.
    ! Note also that the loop is unrolled.
    ax = 0.
    ileft = (nmax/8)*8

    do i = 1,ileft,8
        aax = abs(array(i))

        if (aax .gt. ax) then
            ax = aax
        end if

        aax = abs(array(i + 1))

        if (aax .gt. ax) then
            ax = aax
        end if

        aax = abs(array(i + 2))

        if (aax .gt. ax) then
            ax = aax
        end if

        aax = abs(array(i + 3))

        if (aax .gt. ax) then
            ax = aax
        end if

        aax = abs(array(i + 4))

        if (aax .gt. ax) then
            ax = aax
        end if

        aax = abs(array(i + 5))

        if (aax .gt. ax) then
            ax = aax
        end if

        aax = abs(array(i + 6))

        if (aax .gt. ax) then
            ax = aax
        end if

        aax = abs(array(i + 7))

        if (aax .gt. ax) then
            ax = aax
        end if
    end do

    do i = ileft + 1,nmax
        aax = abs(array(i))

        if (aax .gt. ax) then
            ax = aax
        end if
    end do

    arrmxn = ax

999 continue
end function arrmxn