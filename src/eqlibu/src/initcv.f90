subroutine initcv(uarray,nmax,uvalue)
    !! This subroutine initializes the character array uarray to the
    !! string contained in uvalue over the first nmax positions.
    !! Normally, nmax would be the dimension of the 1D uarray. However,
    !! nmax could be less than the true dimension. Also, uarray could
    !! actually have more than one dimension. For example, if the array
    !! is really ua, which has dimensions of i1 and i2, this subroutine
    !! could be used by calling it in the following manner:
    !! call initcv((ua,(i1*i2),uv). Use this subroutine with caution if
    !! nmax is not the product of the true (declared) dimensions.
    !! To initialize a character array to blanks, use EQLIBU/initcb.f
    !! instead.
    !! This subroutine is called by:
    !!   Any
    !! Input:
    !!   uarray = the array
    !!   nmax   = dimension or pseudo-dimension of uarray
    !!   uvalue = the assigned value
    !! Output:
    !!   uarray = uarray, with the first nmax positions set to the
    !!              string contained in uv
    implicit none

    ! Calling sequence variable declarations.
    integer :: nmax

    character(len=*) :: uarray(nmax)
    character(len=*) :: uvalue

    ! Local variable declarations.
    integer :: i
    integer :: ileft

    character(len=80) :: uv

    ! Note the assignment of the value to a local variable.
    ! Note that the loop is unrolled.
    uv = uvalue
    ileft = (nmax/8)*8

    do i = 1,ileft,8
        uarray(i) = uv
        uarray(i + 1) = uv
        uarray(i + 2) = uv
        uarray(i + 3) = uv
        uarray(i + 4) = uv
        uarray(i + 5) = uv
        uarray(i + 6) = uv
        uarray(i + 7) = uv
    end do

    do i = ileft + 1,nmax
        uarray(i) = uv
    end do
end subroutine initcv