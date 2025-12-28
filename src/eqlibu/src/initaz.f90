subroutine initaz(array,nmax)
    !! This subroutine initializes the real*8 array "array" to zero over
    !! the first nmax positions. Normally, nmax would be the dimension
    !! of the 1D "array". However, nmax could be less than the true
    !! dimension. Also, array could actually have more than one
    !! dimension. For example, if the array is really b, which has
    !! dimensions of i1 and i2, this subroutine could be used by calling
    !! it in the following manner: call initaz(b,(i1*i2)). Use this
    !! subroutine with caution if the two arrays do not have identical
    !! dimensioning or nmax is not the product of the true (declared)
    !! dimensions.
    !! To initialize a real*8 array to a non-zero value, use
    !! EQLIBU/initav.f instead.
    !! NOTE: It may be more efficient to just initialize the array
    !! to zero in a DO loop in the calling subroutine.
    !! This subroutine is called by:
    !!   Any
    !! Input:
    !!   array  = the array
    !!   nmax   = dimension or pseudo-dimension of array
    !! Output:
    !!   array  = array, with the first nmax positions set to zero
    implicit none

    ! Calling sequence variable declarations.
    integer :: nmax

    real(kind=8) :: array(nmax)

    ! Local variable declarations.
    integer :: i

    ! Caution: efficiency may be best served by not unrolling the
    ! following loop.
    do i = 1,nmax
        array(i) = 0.0
    end do
end subroutine initaz