integer function idamax(nmax,array,incx)
    !! This subroutine finds the index of the element corresponding to
    !! the max norm of the real*8 array "array". If there is more than
    !! one such element, the index returned corresponds to the first one.
    !! This subroutine is an adaptation of the 1979 Linpack subroutine of
    !! the same name. The size and order in the original calling sequence
    !! has been preserved: (n,dx,incx) == (nmax,array,incx).
    !! The increment incx is not used here (it is taken as having
    !! a value of 1; input of any other value constitutes an error).
    !! This is a pseudo-Linpack BLAS (Basic Linear Algebra Subsystem)
    !! subroutine.
    !! This subroutine is nearly identical to EQLIBU/iarmnx.f, for which
    !! the calling sequence arguments are reversed.
    !! This subroutine is called by:
    !!   EQLIBU/dgefa.f
    !! Input:
    !!   nmax   = dimension or pseudo-dimension of array
    !!   array  = the array
    !!   incx   = the increment (must have a value of 1)
    !! Output:
    !!   idamax = index of the element corresponding to the max norm
    !!              (zero if nmax .le. 0)
    implicit none

    ! Calling sequence variable declarations.
    integer :: nmax
    integer :: incx

    real(kind=8) :: array(nmax)

    ! Local variable declarations.
    integer :: i
    integer :: ileft
    integer :: ix
    integer :: j

    real(kind=8) :: aax
    real(kind=8) :: arrmxn

    ! Trap illegal value of incx.
    if (incx .ne. 1) then
        write (6,1000) incx
1000 format(/' * Error - (EQLIBU/idamax) The argument incx has a',/7x,'value of ',i5,'. This version of idamax only allows this',/7x,'argument to have a value of 1. This is a programming',/7x,'error or a linking error.')

        stop
    end if

    idamax = 0
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

    idamax = ix

999 continue
end function idamax