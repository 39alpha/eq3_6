subroutine aderiv(akmat0,daffp0,fdaf0,nord,nordmx,npt,nptmax,nrd1mx)
    !! This subroutine computes estimates of the derivatives of the phase
    !! affinities (affp) from the corresponding finite differences.
    !! Note that (daffp0) = (akmat0)(fdaf0).
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
    integer :: nptmax
    integer :: nrd1mx

    integer :: nord
    integer :: npt

    real(kind=8) :: akmat0(nrd1mx,nrd1mx)
    real(kind=8) :: daffp0(nordmx,nptmax)
    real(kind=8) :: fdaf0(nordmx,nptmax)

    ! Local variable declarations.
    integer :: k
    integer :: n
    integer :: np

    real(kind=8) :: dx

    do np = 1,npt
        do n = 1,nord
            dx = 0.

            do k = n,nord
                dx = dx + fdaf0(k,np)*akmat0(n,k)
            end do

            daffp0(n,np) = dx
        end do
    end do
end subroutine aderiv