subroutine scalx1(avx,avxmax,avxs,ier,nmax)
    !! This subroutine scales the elements in the avx array to the
    !! interval (-1,1). The results are placed in the avxs array. This
    !! subroutine is normally used in conjunction with EQLIBU/rscaly.f
    !! and EQLIBU/polfit.f to fitting interpolating polynomials to
    !! data in which avx contains values of the independent variable.
    !! This subroutine is called by:
    !!   EQPT/eqpt.f
    !!   Any
    !! Input:
    !!   avx    = array whose contents are to be scaled
    !!   avxmax = the max norm of avx
    !!   n      = the number of elements in avx
    !! Output:
    !!   avxs   = array whose contents are to be scaled
    !!   ier    = error flag
    !!              = 0  okay
    !!              = 1  all the elements of avx are zero
    implicit none

    ! Calling sequence variable declarations.
    integer :: ier
    integer :: nmax

    real(kind=8) :: avx(nmax)
    real(kind=8) :: avxs(nmax)
    real(kind=8) :: avxmax

    ! Local variable declarations.
    integer :: i

    real(kind=8) :: ax

    ! Initialize error flag.
    ier = 0

    ! Note- nmax is expected to be relatively small. Therefore, loop
    ! unrolling or calling of subroutines which employ loop unrolling
    ! is not used here.
    ! Find the max norm of avx.
    avxmax = 0.

    do i = 1,nmax
        ax = abs(avx(i))

        if (ax .gt. avxmax) then
            avxmax = ax
        end if
    end do

    if (avxmax .le. 0.) then
        ier = 1
        go to 999
    end if

    ! Scale the array.
    do i = 1,nmax
        avxs(i) = avx(i)/avxmax
    end do

999 continue
end subroutine scalx1