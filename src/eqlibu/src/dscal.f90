subroutine dscal(nmax,scfact,array,incx)
    !! This subroutine multiplies the real*8 array "array" by the real*8
    !! constant scfact. The result is returned in "array". This
    !! subroutine is an adaptation of the 1979 Linpack subroutine of the
    !! same name. The size and order in the original calling sequence has
    !! been preserved: (n,da,dx,incx) == (nmax,scfact,array,incx).
    !! The increment incx is not used here (it is taken as having
    !! a value of 1; input of any other value constitutes an error).
    !! This is a pseudo-Linpack BLAS (Basic Linear Algebra Subsystem)
    !! subroutine.
    !! This subroutine is called by:
    !!   EQLIBU/dgefa.f
    !! Input:
    !!   nmax   = dimension or pseudo-dimension of array
    !!   scfact = the scale factor
    !!   array  = the array
    !!   incx   = the increment (must have a value of 1)
    !! Output:
    !!   array  = the original array multiplied by the scale factor.
    implicit none

    ! Calling sequence variable declarations.
    integer :: nmax
    integer :: incx

    real(kind=8) :: array(nmax)
    real(kind=8) :: scfact

    ! Local variable declarations.
    integer :: i
    integer :: ileft

    real(kind=8) :: scx

    ! Trap illegal value of incx.
    if (incx .ne. 1) then
        write (6,1000) incx
1000 format(/' * Error - (EQLIBU/dscal) The argument incx has a',/7x,'value of ',i5,'. This version of dscal only allows this',/7x,'argument to have a value of 1. This is a programming',/7x,'error or a linking error.')

        stop
    end if

    ! Note that the loop is unrolled.
    ileft = (nmax/8)*8
    scx = scfact

    do i = 1,ileft,8
        array(i) = scx*array(i)
        array(i + 1) = scx*array(i + 1)
        array(i + 2) = scx*array(i + 2)
        array(i + 3) = scx*array(i + 3)
        array(i + 4) = scx*array(i + 4)
        array(i + 5) = scx*array(i + 5)
        array(i + 6) = scx*array(i + 6)
        array(i + 7) = scx*array(i + 7)
    end do

    do i = ileft + 1,nmax
        array(i) = scx*array(i)
    end do
end subroutine dscal