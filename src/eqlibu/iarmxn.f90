integer function iarmxn(array,nmax)
    !! This subroutine finds the index of the element corresponding to
    !! the max norm of the real*8 array "array". If there is more than
    !! one such element, the index returned corresponds to the first one.
    !! Note that array could actually have more than one dimension. For
    !! example, if the array is really "a", which has dimensions of i1
    !! and i2, this subroutine could be used by calling it in the
    !! following manner: call iarmxn(a,(i1*i2)). Use this subroutine
    !! with caution if nmax is not the product of the true (declared)
    !! dimensions.
    !! To just find the max norm itself, use EQLIBU/arrmxn.f instead.
    !! This subroutine is nearly identical to EQLIBU/idamax.f, for which
    !! the calling sequence arguments are reversed.
    !! This subroutine is called by:
    !!   Any
    !! Input:
    !!   array  = the array
    !!   nmax   = dimension or pseudo-dimension of the array
    !! Output:
    !!   iarmxn = index of the element corresponding to the max norm
    !!              (zero if nmax .le. 0)
    implicit none

    ! Calling sequence variable declarations.
    integer :: nmax

    real(kind=8) :: array(nmax)

    ! Local variable declarations.
    integer :: i
    integer :: ileft
    integer :: ix
    integer :: j

    real(kind=8) :: aax
    real(kind=8) :: arrmxn

    iarmxn = 0
    arrmxn = 0.

    if (nmax .le. 0) then
        go to 999
    end if

    ! Note the use of a local variable (ix) within the loop.
    ! Note also that the loop is unrolled.
    ix = 0
    ileft = (nmax/8)*8

    do i = 1,ileft,8
        aax = abs(array(i))

        if (aax .gt. arrmxn) then
            ix = i
            arrmxn = aax
        end if

        j = i + 1
        aax = abs(array(j))

        if (aax .gt. arrmxn) then
            ix = j
            arrmxn = aax
        end if

        j = i + 2
        aax = abs(array(j))

        if (aax .gt. arrmxn) then
            ix = j
            arrmxn = aax
        end if

        j = i + 3
        aax = abs(array(j))

        if (aax .gt. arrmxn) then
            ix = j
            arrmxn = aax
        end if

        j = i + 4
        aax = abs(array(j))

        if (aax .gt. arrmxn) then
            ix = j
            arrmxn = aax
        end if

        j = i + 5
        aax = abs(array(j))

        if (aax .gt. arrmxn) then
            ix = j
            arrmxn = aax
        end if

        j = i + 6
        aax = abs(array(j))

        if (aax .gt. arrmxn) then
            ix = j
            arrmxn = aax
        end if

        j = i + 7
        aax = abs(array(j))

        if (aax .gt. arrmxn) then
            ix = j
            arrmxn = aax
        end if
    end do

    do i = ileft + 1,nmax
        aax = abs(array(i))

        if (aax .gt. arrmxn) then
            ix = i
            arrmxn = aax
        end if
    end do

    iarmxn = ix

999 continue
end function iarmxn