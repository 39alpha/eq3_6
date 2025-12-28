subroutine copyia(iarr1,iarr2,nmax)
    !! This subroutine copies the first nmax elements of the integer
    !! array iarr1 into the integer array iarr2. Normally, nmax would
    !! be the dimension of both arrays. However, nmax could be less
    !! than the true dimension. Also, the arrays could actually have
    !! more than one dimension. For example, if iarr1 and iarr2
    !! are really ia and ib, respectively, having common dimensions of
    !! i1 and i2, this subroutine could be used by calling it in the
    !! following manner: call copyia(ia,ib,(i1*i2)). Use this subroutine
    !! with caution if the two arrays do not have identical dimensioning
    !! or nmax is not the product of the actual (declared) dimensions.
    !! This subroutine is called by:
    !!   Any
    !! Input:
    !!   iarr1  = the first array
    !!   nmax   = dimension or pseudo-dimension of iarr1 and iarr2
    !! Output:
    !!   iarr2  = the second array
    implicit none

    ! Calling sequence variable declarations.
    integer :: nmax

    integer :: iarr1(nmax)
    integer :: iarr2(nmax)

    ! Local variable declarations.
    integer :: i
    integer :: ileft

    ! Note that the loop is unrolled.
    ileft = (nmax/8)*8

    do i = 1,ileft,8
        iarr2(i) = iarr1(i)
        iarr2(i + 1) = iarr1(i + 1)
        iarr2(i + 2) = iarr1(i + 2)
        iarr2(i + 3) = iarr1(i + 3)
        iarr2(i + 4) = iarr1(i + 4)
        iarr2(i + 5) = iarr1(i + 5)
        iarr2(i + 6) = iarr1(i + 6)
        iarr2(i + 7) = iarr1(i + 7)
    end do

    do i = ileft + 1,nmax
        iarr2(i) = iarr1(i)
    end do
end subroutine copyia