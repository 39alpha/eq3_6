subroutine xderiv(akmat0,dxx0,fdxx0,nord,nordmx,nrd1mx)
    !! This subroutine computes estimates of the derivatives (dxx0)
    !! for some quantity (such as pH) from the corresponding finite
    !! differences (fdxx0). Note that (dxx0) = (akmat0)(fdxx0).
    !! Compare with:
    !!   EQ6/bderiv.f
    !!   EQ6/pderiv.f
    !!   EQ6/rderiv.f
    !!   EQ6/zderiv.f
    !! This subroutine is called by:
    !!   EQ6/path.f
    !! Principal input:
    !! Principal output:
    implicit none

    integer :: nordmx
    integer :: nrd1mx

    integer :: nord

    real(kind=8) :: akmat0(nrd1mx,nrd1mx)
    real(kind=8) :: dxx0(nordmx)
    real(kind=8) :: fdxx0(nordmx)

    ! Local variable declarations.
    integer :: k
    integer :: n

    real(kind=8) :: dx

    do n = 1,nord
        dx = 0.

        do k = n,nord
            dx = dx + fdxx0(k)*akmat0(n,k)
        end do

        dxx0(n) = dx
    end do
end subroutine xderiv