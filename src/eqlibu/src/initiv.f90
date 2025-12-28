subroutine initiv(iarray,nmax,ivalue)
    !! This subroutine initializes the integer array iarray to ivalue
    !! over the first nmax positions. Normally, nmax would be the
    !! dimension of the 1D iarray. However, nmax could be less than the
    !! true dimension. Also, iarray could actually have more than one
    !! dimension. For example, if the array is really ia, which has
    !! dimensions of i1 and i2, this subroutine could be used by calling
    !! it in the following manner: call initiv(ia,(i1*i2),ivalue).
    !! Use this subroutine with caution if nmax is not the product of the
    !! true (declared) dimensions.
    !! To initialize an integer array to a zero value, use
    !! EQLIBU/initiz.f instead.
    !! This subroutine is called by:
    !!   Any
    !! Input:
    !!   iarray = the array
    !!   nmax   = dimension or pseudo-dimension of array
    !!   ivalue = the initialization value
    !! Output:
    !!   iarray = iarray, with the first nmax positions set to the
    !!            initialization value
    implicit none

    ! Calling sequence variable declarations.
    integer :: nmax

    integer :: iarray(nmax)
    integer :: ivalue

    ! Local variable declarations.
    integer :: i
    integer :: ileft
    integer :: iv

    ! Note the assignment of the value to a local variable.
    ! Note that the loop is unrolled.
    iv = ivalue
    ileft = (nmax/8)*8

    do i = 1,ileft,8
        iarray(i) = iv
        iarray(i + 1) = iv
        iarray(i + 2) = iv
        iarray(i + 3) = iv
        iarray(i + 4) = iv
        iarray(i + 5) = iv
        iarray(i + 6) = iv
        iarray(i + 7) = iv
    end do

    do i = ileft + 1,nmax
        iarray(i) = iv
    end do
end subroutine initiv