subroutine gdd(dgpit,fxi,gpit,ipbtmx,napmax,napt,palpha)
    !! This subroutine computes the functions g(xi) and their ionic
    !! strength derivatives (xi = alpha(i) * I**1/2). These functions
    !! are used in Pitzer's equations. Note: alpha(i) = 0 implies
    !! that g(xi) and its derivatives are zero.
    !! This subroutine is called by:
    !!   EQLIBG/gcoeff.f
    !! Principal input:
    !!   fxi    = the ionic strength (the 2nd-order electrostatic
    !!              moment function I)
    !!   napt   = number of alpha coefficient pairs in the palpha array
    !! Principal output:
    !!   gpit   = array of values for the g(x) function
    !!   dgpit  = array of values for the ionic strength derivatives
    !!              of the g(x) function
    implicit none

    ! Calling sequence variable declarations.
    integer :: ipbtmx
    integer :: napmax

    integer :: napt

    real(kind=8) :: dgpit(2,ipbtmx,napmax)
    real(kind=8) :: gpit(ipbtmx,napmax)
    real(kind=8) :: palpha(ipbtmx,napmax)
    real(kind=8) :: fxi

    ! Local variable declarations.
    integer :: i
    integer :: nap

    real(kind=8) :: alphai
    real(kind=8) :: dxdi
    real(kind=8) :: emx
    real(kind=8) :: gpx
    real(kind=8) :: gppx
    real(kind=8) :: x
    real(kind=8) :: xsq

    do nap = 1,napt
        do i = 1,ipbtmx
            gpit(i,nap) = 0.
            dgpit(1,i,nap) = 0.
            dgpit(2,i,nap) = 0.
            alphai = palpha(i,nap)
            x = alphai * sqrt(fxi)

            if (x .ne. 0.) then
                xsq = x*x
                emx = exp(-x)

                ! Calculate g(x).
                gpit(i,nap) = (2./xsq) * (1. - emx*(1. + x))

                ! Calculate the derivatives with respect to x.
                gpx = (-2./(xsq*x)) * (2. - emx*(2. + x*(2. + x)))
                gppx = (-2./(xsq*xsq))      * (-6. + emx*(6. + x*(6. + x*(3. + x))))

                ! Calculate the derivatives with respect to ionic strength.
                dxdi = alphai/(2.*sqrt(fxi))
                dgpit(1,i,nap) = dxdi * gpx
                dgpit(2,i,nap) = dxdi * (gppx*dxdi - 0.5*gpx/fxi)
            end if
        end do
    end do
end subroutine gdd