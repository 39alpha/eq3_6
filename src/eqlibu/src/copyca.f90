subroutine copyca(uarr1,uarr2,nmax)
    !! This subroutine copies the first nmax elements of the character
    !! array uarr1 into the character array uarr2. Normally, nmax would
    !! be the dimension of both arrays. However, nmax could be less
    !! than the true dimension. Also, the arrays could actually have
    !! more than one dimension. For example, if uarr1 and uarr2
    !! are really a and b, respectively, having common dimensions of
    !! i1 and i2, this subroutine could be used by calling it in the
    !! following manner: call copyca(a,b,(i1*i2)). Use this subroutine
    !! with caution if the two arrays do not have identical dimensioning,
    !! nmax is not the product of the actual (declared) dimensions, or
    !! the character lengths of the two arrays are not identical.
    !! This subroutine is called by:
    !!   Any
    !! Input:
    !!   uarr1  = the first array
    !!   nmax   = dimension or pseudo-dimension of uarr1 and uarr2
    !! Output:
    !!   uarr2  = the second array
    implicit none

    ! Calling sequence variable declarations.
    integer :: nmax

    character(len=*) :: uarr1(nmax)
    character(len=*) :: uarr2(nmax)

    ! Local variable declarations.
    integer :: i
    integer :: ileft

    ! Note that the loop is unrolled.
    ileft = (nmax/8)*8

    do i = 1,ileft,8
        uarr2(i) = uarr1(i)
        uarr2(i + 1) = uarr1(i + 1)
        uarr2(i + 2) = uarr1(i + 2)
        uarr2(i + 3) = uarr1(i + 3)
        uarr2(i + 4) = uarr1(i + 4)
        uarr2(i + 5) = uarr1(i + 5)
        uarr2(i + 6) = uarr1(i + 6)
        uarr2(i + 7) = uarr1(i + 7)
    end do

    do i = ileft + 1,nmax
        uarr2(i) = uarr1(i)
    end do
end subroutine copyca