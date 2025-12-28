subroutine daxpy(nmax,scfact,array1,incx,array2,incy)
    !! This subroutine multiplies the real*8 array "array1" by the
    !! real*8 constant scfact and adds that to the real*8 array "array2".
    !! The result is returned in "array2". The relevant equation is
    !! (y) = a*(x) + (y). This subroutine is an adaptation of the 1979
    !! Linpack subroutine of the same name. The size and order in the
    !! original calling sequence has been preserved:
    !! (n,da,dx,incx,dy,incy) == (nmax,scfact,array1,incx,array2,incy).
    !! The increments incx and incy are not used here (each is taken
    !! as having a value of 1; input of any other value constitutes
    !! an error). This is a pseudo-Linpack BLAS (Basic Linear Algebra
    !! Subsystem) subroutine.
    !! This subroutine is called by:
    !!   EQLIBU/dgefa.f
    !!   EQLIBU/dgesl.f
    !! Input:
    !!   nmax   = dimension or pseudo-dimension of array
    !!   scfact = the scale factor
    !!   array1 = the first array (x)
    !!   array2 = the second array (y)
    !!   incx   = the x increment (must have a value of 1)
    !!   incy   = the y increment (must have a value of 1)
    !! Output:
    !!   array2 = the original array multiplied by the scale factor.
    implicit none

    ! Calling sequence variable declarations.
    integer :: nmax
    integer :: incx
    integer :: incy

    real(kind=8) :: array1(nmax)
    real(kind=8) :: array2(nmax)

    real(kind=8) :: scfact

    ! Local variable declarations.
    integer :: i
    integer :: ileft

    real(kind=8) :: scx

    ! Trap illegal values of incx and incy.
    if (incx .ne. 1) then
        write (6,1000) incx
1000 format(/' * Error - (EQLIBU/daxpy) The argument incx has a',/7x,'value of ',i5,'. This version of daxpy only allows this',/7x,'argument to have a value of 1. This is a programming',/7x,'error or a linking error.')

        stop
    end if

    if (incy .ne. 1) then
        write (6,1010) incy
1010 format(/' * Error - (EQLIBU/daxpy) The argument incy has a',/7x,'value of ',i5,'. This version of daxpy only allows this',/7x,'argument to have a value of 1. This is a programming',/7x,'error or a linking error.')

        stop
    end if

    ! Note that the loop is unrolled.
    scx = scfact
    ileft = (nmax/8)*8

    do i = 1,ileft,8
        array2(i) = array2(i) + scx*array1(i)
        array2(i + 1) = array2(i + 1) + scx*array1(i + 1)
        array2(i + 2) = array2(i + 2) + scx*array1(i + 2)
        array2(i + 3) = array2(i + 3) + scx*array1(i + 3)
        array2(i + 4) = array2(i + 4) + scx*array1(i + 4)
        array2(i + 5) = array2(i + 5) + scx*array1(i + 5)
        array2(i + 6) = array2(i + 6) + scx*array1(i + 6)
        array2(i + 7) = array2(i + 7) + scx*array1(i + 7)
    end do

    do i = ileft + 1,nmax
        array2(i) = array2(i) + scx*array1(i)
    end do
end subroutine daxpy