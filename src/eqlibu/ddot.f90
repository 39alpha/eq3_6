real(kind=8) function ddot(nmax,array1,incx,array2,incy)
    !! This subroutine computes the dot product (ddot) of the two real*8
    !! arrays "array1" and "array2". This subroutine is an adaptation of
    !! the 1979 Linpack subroutine of the same name. The size and order
    !! in the original calling sequence has been preserved:
    !! (n,dx,incx,dy,incy) = (nmax,array1,incx,array2,incy).
    !! The increments incx and incy are not used here (each is taken
    !! as having a value of 1; input of any other value constitutes
    !! an error). This is a pseudo-Linpack BLAS (Basic Linear Algebra
    !! Subsystem) subroutine.
    !! This subroutine is nearly identical to EQLIBU/dotpra.f, which has
    !! the same function, but a different calling sequence.
    !! This subroutine is called by:
    !!   Any
    !! Input:
    !!   nmax   = dimension or pseudo-dimension of array
    !!   array1 = the first array
    !!   incx   = the x increment (must have a value of 1)
    !!   array2 = the second array
    !!   incy   = the y increment (must have a value of 1)
    !! Output:
    !!   ddot   = the dot product
    implicit none

    ! Calling sequence variable declarations.
    integer :: nmax
    integer :: incx
    integer :: incy

    real(kind=8) :: array1(nmax)
    real(kind=8) :: array2(nmax)

    ! Local variable declarations.
    integer :: i
    integer :: ileft

    real(kind=8) :: dx

    ! Trap illegal values of incx and inc6.
    if (incx .ne. 1) then
        write (6,1000) incx
1000 format(/' * Error - (EQLIBU/ddot) The argument incx has a',/7x,'value of ',i5,'. This version of ddot only allows this',/7x,'argument to have a value of 1. This is a programming',/7x,'error or a linking error.')

        stop
    end if

    if (incy .ne. 1) then
        write (6,1010) incy
1010 format(/' * Error - (EQLIBU/ddot) The argument incy has a',/7x,'value of ',i5,'. This version of ddot only allows this',/7x,'argument to have a value of 1. This is a programming',/7x,'error or a linking error.')

        stop
    end if

    ! Note the use of a local variable (dx) within the loop.
    ! Note also that the loop is unrolled.
    dx = 0.
    ileft = (nmax/8)*8

    do i = 1,ileft,8
        dx = dx + array1(i)*array2(i) + array1(i + 1)*array2(i + 1)  + array1(i + 2)*array2(i + 2) + array1(i + 3)*array2(i + 3)  + array1(i + 4)*array2(i + 4) + array1(i + 5)*array2(i + 5)  + array1(i + 6)*array2(i + 6) + array1(i + 7)*array2(i + 7)
    end do

    do i = ileft + 1,nmax
        dx = dx + array1(i)*array2(i)
    end do

    ddot = dx
end function ddot