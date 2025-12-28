subroutine initii(iarray,nmax)
    !! This subroutine initializes each element of the integer array
    !! iarray its index value over the first nmax positions. Normally,
    !! nmax would be the dimension of the 1D iarray. However, nmax could
    !! be less than the true dimension. Also, iarray could actually have
    !! more than one dimension. For example, if the array is really ib,
    !! which has dimensions of i1 and i2, this subroutine could be used
    !! by calling it in the following manner: call initii(ib,(i1*i2)).
    !! Use this subroutine with caution if nmax is not the product of the
    !! true (declared) dimensions.
    !! This subroutine is called by:
    !!   Any
    !! Input:
    !!   iarray = the array
    !!   nmax   = dimension or pseudo-dimension of iarray
    !! Output:
    !!   iarray = iarray, with the first nmax positions set to the
    !!              corresponding index values
    implicit none

    ! Calling sequence variable declarations.
    integer :: nmax

    integer :: iarray(nmax)

    ! Local variable declarations.
    integer :: i
    integer :: ileft

    ! Note that the loop is unrolled.
    ileft = (nmax/8)*8

    do i = 1,ileft,8
        iarray(i) = i
        iarray(i + 1) = i + 1
        iarray(i + 2) = i + 2
        iarray(i + 3) = i + 3
        iarray(i + 4) = i + 4
        iarray(i + 5) = i + 5
        iarray(i + 6) = i + 6
        iarray(i + 7) = i + 7
    end do

    do i = ileft + 1,nmax
        iarray(i) = i
    end do
end subroutine initii