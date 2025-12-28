subroutine rscaly(avxmax,avy,avys,eps100,nmax)
    !! This subroutine rescales the elements in the avys array. The
    !! results are placed in the avy array. This subroutine is normally
    !! used in conjunction with EQLIBU/scalx1.f and EQLIBU/polfit.f
    !! to fitting interpolating polynomials to data in which avx
    !! contains values of the independent variable, avxs contains
    !! values of the independent variable scaled to the interval
    !! (-1,1), avys contains values of the dependent variable as
    !! fit with the independent variable scaled, and avy contains
    !! values of the dependent variable corrected (rescaled) for
    !! scaling of the in dependent variable.
    !! This subroutine is called by:
    !!   EQPT/intrp.f
    !!   Any
    !! Input:
    !!   avys   = array whose contents are to be rescaled
    !!   avxmax = the max norm of avx
    !!   n      = the number of elements in avys
    !! Output:
    !!   avy    = array whose contents are rescaled
    implicit none

    ! Calling sequence variable declarations.
    integer :: nmax

    real(kind=8) :: avy(nmax)
    real(kind=8) :: avys(nmax)
    real(kind=8) :: avxmax
    real(kind=8) :: eps100

    ! Local variable declarations.
    integer :: i

    real(kind=8) :: ax
    real(kind=8) :: ay

    ! Note- nmax is expected to be relatively small. Therefore, loop
    ! unrolling or calling of subroutines which employ loop unrolling
    ! is not used here.
    ! Test avxmax.
    ax = avxmax - 1.0

    if (abs(ax) .le. eps100) then
        do i = 1,nmax
            avy(i) = avys(i)
        end do

        go to 999
    end if

    ! Rescale the array.
    ay = 1.
    avy(1) = avys(1)

    do i = 2,nmax
        ay = ay/avxmax
        avy(i) = avys(i)*ay
    end do

999 continue
end subroutine rscaly