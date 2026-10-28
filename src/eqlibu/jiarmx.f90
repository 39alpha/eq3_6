integer function jiarmx(iarray,nmax)
    !! This subroutine finds the maximum element of the integer array
    !! iarray. Note that iarray could actually have more than one
    !! dimension. For example, if the array is really "ia", which has
    !! dimensions of i1 and i2, this subroutine could be used by calling
    !! it in the following manner: call jiarmx(ia,(i1*i2)). Use this
    !! subroutine with caution if nmax is not the product of the true
    !! (declared) dimensions.
    !! To find the corresponding minimum value, use EQLIBU/jiarmn.f
    !! instead.
    !! This subroutine is called by:
    !!   Any
    !! Input:
    !!   iarray = the array
    !!   nmax   = dimension or pseudo-dimension of the array
    !! Output:
    !!   jiarmx = the maximum element of the array (zero if nmax .le. 0)
    implicit none

    ! Calling sequence variable declarations.
    integer :: nmax

    integer :: iarray(nmax)

    ! Local variable declarations.
    integer :: i
    integer :: ileft
    integer :: ix

    jiarmx = 0

    if (nmax .le. 0) then
        go to 999
    end if

    ! Note the use of a local variable (ix) within the loop.
    ! Note also that the loop is unrolled.
    ix = iarray(1)
    ileft = (nmax/8)*8

    do i = 1,ileft,8
        ix = max(ix,iarray(i),iarray(i + 1),iarray(i + 2),iarray(i + 3),iarray(i + 4),iarray(i + 5),iarray(i + 6),iarray(i + 7))
    end do

    do i = ileft + 1,nmax
        ix = max(ix,iarray(i))
    end do

    jiarmx = ix

999 continue
end function jiarmx