subroutine multav(array1,nmax,avalue,array2)
    !! This subroutine multiplies the first nmax elements of the real*8
    !! array array1 by the real*8 value avalue. The result is placed
    !! in the real*8 array array2. Normally, nmax would be the dimension
    !! of both arrays. However, nmax could be less than the true
    !! dimension. Also, the arrays could actually have more than one
    !! dimension. For example, if array1 and array2 are really a and b,
    !! respectively, having common dimensions of i1 and i2, this
    !! subroutine could be used by calling it in the following manner:
    !! call multav(a,(i1*i2),avalue,b). Use this subroutine with caution
    !! if the two arrays do not have identical dimensioning or nmax is
    !! not the product of the true (declared) dimensions.
    !! This subroutine is called by:
    !!   Any
    !! Input:
    !!   array1 = the first array
    !!   nmax   = dimension or pseudo-dimension of array1 and array2
    !!   avalue = the value by which to multipy array1
    !! Output:
    !!   array2 = the second array, the product of avalue and array1
    implicit none

    ! Calling sequence variable declarations.
    integer :: nmax

    real(kind=8) :: array1(nmax)
    real(kind=8) :: array2(nmax)

    real(kind=8) :: avalue

    ! Local variable declarations.
    integer :: i
    integer :: ileft

    real(kind=8) :: av

    ! Note the assignment of the value to a local variable.
    ! Note also that the loop is unrolled.
    av = avalue
    ileft = (nmax/8)*8

    do i = 1,ileft,8
        array2(i) = av*array1(i)
        array2(i + 1) = av*array1(i + 1)
        array2(i + 2) = av*array1(i + 2)
        array2(i + 3) = av*array1(i + 3)
        array2(i + 4) = av*array1(i + 4)
        array2(i + 5) = av*array1(i + 5)
        array2(i + 6) = av*array1(i + 6)
        array2(i + 7) = av*array1(i + 7)
    end do

    do i = ileft + 1,nmax
        array2(i) = av*array1(i)
    end do
end subroutine multav