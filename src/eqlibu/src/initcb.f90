subroutine initcb(uarray,nmax)
    !! This subroutine initializes the character array uarray to blanks
    !! over the first nmax positions. Normally, nmax would be the
    !! dimension of the 1D uarray. However, nmax could be less than the
    !! true dimension. Also, uarray could actually have more than one
    !! dimension. For example, if the array is really ub, having
    !! dimensions of i1 and i2, this subroutine could be used by calling
    !! it in the following manner: call initcb(ub,i1*i2)). Use this
    !! subroutine with caution if nmax is not the product of the true
    !! (declared) dimensions.
    !! To initialize a character array to a non-blank string, use
    !! EQLIBU/initcv.f instead.
    !! This subroutine is called by:
    !!   Any
    !! Input:
    !!   uarray = the array
    !!   nmax   = dimension or pseudo-dimension of uarray
    !! Output:
    !!   uarray = uarray, with the first nmax positions set to zero
    implicit none

    ! Calling sequence variable declarations.
    integer :: nmax

    character(len=*) :: uarray(nmax)

    ! Local variable declarations.
    integer :: i

    ! Caution: efficiency may be best served by not unrolling the
    ! following loop.
    do i = 1,nmax
        uarray(i) = ' '
    end do
end subroutine initcb