subroutine copyaa(array1,array2,nmax)
    !! This subroutine copies the first nmax elements of the real*8
    !! array array1 into the real*8 array array2. Normally, nmax would
    !! be the dimension of both arrays. However, nmax could be less
    !! than the true dimension. Also, the arrays could actually have
    !! more than one dimension. For example, if array1 and array2
    !! are really a and b, respectively, having common dimensions of
    !! i1 and i2, this subroutine could be used by calling it in the
    !! following manner: call copyaa(a,b,(i1*i2)). Use this subroutine
    !! with caution if the two arrays do not have identical dimensioning
    !! or nmax is not the product of the true (declared) dimensions.
    !! This subroutine is called by:
    !!   Any
    !! Input:
    !!   array1 = the first array
    !!   nmax   = dimension or pseudo-dimension of array1 and array2
    !! Output:
    !!   array2 = the second array
    implicit none

    ! Calling sequence variable declarations.
    integer :: nmax

    real(kind=8) :: array1(nmax)
    real(kind=8) :: array2(nmax)

    ! Local variable declarations.
    integer :: i
    integer :: ileft

    ! Note that the loop is unrolled.
    ileft = (nmax/8)*8

    do i = 1,ileft,8
        array2(i) = array1(i)
        array2(i + 1) = array1(i + 1)
        array2(i + 2) = array1(i + 2)
        array2(i + 3) = array1(i + 3)
        array2(i + 4) = array1(i + 4)
        array2(i + 5) = array1(i + 5)
        array2(i + 6) = array1(i + 6)
        array2(i + 7) = array1(i + 7)
    end do

    do i = ileft + 1,nmax
        array2(i) = array1(i)
    end do
end subroutine copyaa