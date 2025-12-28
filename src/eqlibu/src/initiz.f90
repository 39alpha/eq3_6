subroutine initiz(iarray,nmax)
    !! This subroutine initializes the integer array iarray to zero over
    !! the first nmax positions. Normally, nmax would be the dimension
    !! of the 1D iarray. However, nmax could be less than the true
    !! dimension. Also, iarray could actually have more than one
    !! dimension. For example, if the array is really ib, which has
    !! dimensions of i1 and i2, this subroutine could be used by calling
    !! it in the following manner: call initiz(ib,(i1*i2)). Use this
    !! subroutine with caution if nmax is not the product of the true
    !! (declared) dimensions.
    !! To initialize an integer array to a non-zero value, use
    !! EQLIBU/initav.f instead.
    !! NOTE: It may be more efficient to just initialize the array
    !! to zero in a DO loop in the calling subroutine.
    !! This subroutine is called by:
    !!   Any
    !! Input:
    !!   iarray = the array
    !!   nmax   = dimension or pseudo-dimension of iarray
    !! Output:
    !!   iarray = iarray, with the first nmax positions set to zero
    implicit none

    ! Calling sequence variable declarations.
    integer :: nmax

    integer :: iarray(nmax)

    ! Local variable declarations.
    integer :: i

    ! Caution: efficiency may be best served by not unrolling the
    ! following loop.
    do i = 1,nmax
        iarray(i) = 0
    end do
end subroutine initiz