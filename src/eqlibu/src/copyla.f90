subroutine copyla(qarr1,qarr2,nmax)
    !! This subroutine copies the first nmax elements of the logical
    !! array qarr1 into the logical array qarr2. Normally, nmax would
    !! be the dimension of both arrays. However, nmax could be less
    !! than the true dimension. Also, the arrays could actually have
    !! more than one dimension. For example, if qarr1 and qarr2
    !! are really qa and qb, respectively, having common dimensions of
    !! i1 and i2, this subroutine could be used by calling it in the
    !! following manner: call copyla(qa,qb,(i1*i2)). Use this subroutine
    !! with caution if the two arrays do not have identical dimensioning
    !! or nmax is not the product of the actual (declared) dimensions.
    !! This subroutine is called by:
    !!   Any
    !! Input:
    !!   qarr1  = the first array
    !!   nmax   = dimension or pseudo-dimension of qarr1 and qarr2
    !! Output:
    !!   qarr2  = the second array
    implicit none

    ! Calling sequence variable declarations.
    integer :: nmax

    logical :: qarr1(nmax)
    logical :: qarr2(nmax)

    ! Local variable declarations.
    integer :: i
    integer :: ileft

    ! Note that the loop is unrolled.
    ileft = (nmax/8)*8

    do i = 1,ileft,8
        qarr2(i) = qarr1(i)
        qarr2(i + 1) = qarr1(i + 1)
        qarr2(i + 2) = qarr1(i + 2)
        qarr2(i + 3) = qarr1(i + 3)
        qarr2(i + 4) = qarr1(i + 4)
        qarr2(i + 5) = qarr1(i + 5)
        qarr2(i + 6) = qarr1(i + 6)
        qarr2(i + 7) = qarr1(i + 7)
    end do

    do i = ileft + 1,nmax
        qarr2(i) = qarr1(i)
    end do
end subroutine copyla