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

    ! Trap illegal value of incx.
    if (incx .ne. 1) then
        write (6,1000) incx
1000 format(/' * Error - (EQLIBU/idamax) The argument incx has a',/7x,'value of ',i5,'. This version of idamax only allows this',/7x,'argument to have a value of 1. This is a programming',/7x,'error or a linking error.')

        stop
    end if

    idamax = maxloc(abs(array), 1)
end function idamax
