subroutine gwhcfa(akmat1,delxi,dxsm11,hhcvec,nord,nordmx,nrd1mx,whcfac,xhcvec)
    !! This subroutine computes the w factor (whcfac), which is
    !! required for higher-order (stiff) ODE corrections. Logically,
    !! w = g dot FCh where g and h are vectors, and (FC) is the matrix
    !! required to convert finite differences utilizing the new point
    !! of reaction progress to derivatives at the base point. This
    !! factor depends only on the recent step size history, including
    !! the current step size (delxi).
    !! This subroutine is called by:
    !!   EQ6/path.f
    !! Principal input:
    !!   akmat1 = the (FC) matrix
    !!   delxi  = the current step size
    !!   dxsm11 = the array of cumulative step sizes going back from
    !!              new point, centered on the new point
    !!   nord   = the order of the finite difference method
    !!   nordmx = the maximum order of the finite difference method
    !!   nrd1mx = nordmx + 1
    !!   hhcvec = the h vector, here effectively work space
    !!   xhcvec = the FCh vector, here effectively work space
    !! Principal output:
    !!   whcfac = the w factor
    implicit none

    ! Calling sequence variable declarations.
    integer :: nordmx
    integer :: nrd1mx

    integer :: nord

    real(kind=8) :: akmat1(nrd1mx,nrd1mx)
    real(kind=8) :: dxsm11(nrd1mx)
    real(kind=8) :: hhcvec(nordmx)
    real(kind=8) :: xhcvec(nordmx)

    real(kind=8) :: delxi
    real(kind=8) :: whcfac

    ! Local variable declarations.
    integer :: i
    integer :: j

    real(kind=8) :: gx

    ! First calculate h (hhcvec).
    hhcvec(1) = 1./delxi

    do j = 2,nord
        hhcvec(j) = hhcvec(j - 1)/dxsm11(j)
    end do

    ! Now calculate the (FCh) vector (xhcvec).
    do i = 1,nord
        xhcvec(i) = 0.

        do j = 1,nord
            xhcvec(i) = xhcvec(i) + akmat1(i,j)*hhcvec(j)
        end do
    end do

    ! Finally, emulate w = g dot (FCh).
    whcfac = 0.
    gx = delxi

    do j = 1,nord
        gx = gx*delxi/(j + 1)
        whcfac = whcfac + gx*xhcvec(j)
    end do
end subroutine gwhcfa