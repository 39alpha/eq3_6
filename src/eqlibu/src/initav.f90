subroutine initav(array,nmax,value)
    !! This subroutine initializes the real*8 array "array" to the value
    !! of "value" over the first nmax positions. Normally, nmax would
    !! be the dimension of the 1D "array". However, nmax could be less
    !! than the true dimension. Also, array could actually have more
    !! than one dimension. For example, if the array is really "a",
    !! having dimensions of i1 and i2, this subroutine could be used by
    !! calling it in the following manner: call initar(a,(i1*i2),value).
    !! Use this subroutine with caution if the two arrays do not have
    !! identical dimensioning or nmax is not the product of the true
    !! (declared) dimensions.
    !! To initialize a real*8 array to a zero value, use EQLIBU/initaz.f
    !! instead.
    !! This subroutine is called by:
    !!   Any
    !! Input:
    !!   array  = the array
    !!   nmax   = dimension or pseudo-dimension of array
    !!   value  = the initialization value
    !! Output:
    !!   array  = array, with the first nmax positions set to the
    !!            initialization value
    implicit none

    ! Calling sequence variable declarations.
    integer :: nmax

    real(kind=8) :: array(nmax)
    real(kind=8) :: value

    ! Local variable declarations.
    integer :: i
    integer :: ileft

    real(kind=8) :: vx

    ! Note the assignment of the value to a local variable.
    ! Note that the loop is unrolled.
    vx = value
    ileft = (nmax/8)*8

    do i = 1,ileft,8
        array(i) = vx
        array(i + 1) = vx
        array(i + 2) = vx
        array(i + 3) = vx
        array(i + 4) = vx
        array(i + 5) = vx
        array(i + 6) = vx
        array(i + 7) = vx
    end do

    do i = ileft + 1,nmax
        array(i) = vx
    end do
end subroutine initav