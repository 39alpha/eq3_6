subroutine gslam(dgpit,dpslm,gpit,ipbtmx,nalpha,napmax,nslt,nsltmx,pslamn,pslm)
    !! This subroutine computes the S-lambda coefficients and their
    !! ionic strength derivatives. These coefficients are used
    !! in Pitzer's equations.
    !! This subroutine is called by:
    !!   EQLIBG/gcoeff.f
    !! Principal input:
    !!   gpit   = array of values for the g(x) function
    !!   dgpit  = array of values of the ionic strength derivatives
    !!              of g(x)
    !!   nslt   = number of S-lambda coeffcient parameter sets
    !!   pslamn = array of S-lambda(n) coefficient parameters
    !!   nalpha = array that gives the index in the palpha array
    !!              corresponding to a given set of S-lambda coefficient
    !!              parameters.
    !! Principal output:
    !!   pslm   = array of the S-lambda functions
    !!   dpslm  = array of ionic strength derivatives of the
    !!              S-lambda functions
    implicit none

    ! Calling sequence variable declarations.
    integer :: ipbtmx
    integer :: napmax
    integer :: nsltmx

    integer :: nalpha(nsltmx)
    integer :: nslt

    real(kind=8) :: dgpit(2,ipbtmx,napmax)
    real(kind=8) :: dpslm(2,nsltmx)
    real(kind=8) :: gpit(ipbtmx,napmax)
    real(kind=8) :: pslamn(0:ipbtmx,nsltmx)
    real(kind=8) :: pslm(nsltmx)

    ! Local variable declarations.
    integer :: i
    integer :: k
    integer :: nap
    integer :: nsl

    real(kind=8) :: dpx
    real(kind=8) :: px

    do nsl = 1,nslt
        nap = nalpha(nsl)

        px = pslamn(0,nsl)

        do i = 1,ipbtmx
            px = px + gpit(i,nap)*pslamn(i,nsl)
        end do

        pslm(nsl) = px

        do k = 1,2
            dpx = 0.

            do i = 1,ipbtmx
                dpx = dpx + dgpit(k,i,nap)*pslamn(i,nsl)
            end do

            dpslm(k,nsl) = dpx
        end do
    end do
end subroutine gslam