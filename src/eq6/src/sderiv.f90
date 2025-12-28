subroutine sderiv(akmat0,demos0,fdse0,nord,nordmx,nrd1mx,nset,nsetmx)
    !! This subroutine computes estimates of the derivatives of the
    !! numbers of moles of selected species in the equilibrium system.
    !! These species include H2O(l) for the aqueous solution and all
    !! species of the non-aqueous phases. These derivatives (demos0)
    !! are computed from the corresponding finite differences (fdse0).
    !! The relation is: (demos0) = (akmat0)(fdse0).
    !! Compare with:
    !!   EQ6/aderiv.f
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
    integer :: nsetmx

    integer :: nord
    integer :: nset

    real(kind=8) :: akmat0(nrd1mx,nrd1mx)
    real(kind=8) :: demos0(nordmx,nsetmx)
    real(kind=8) :: fdse0(nordmx,nsetmx)

    ! Local variable declarations.
    integer :: k
    integer :: n
    integer :: nse

    real(kind=8) :: dx

    do nse = 1,nset
        do n = 1,nord
            dx = 0.

            do k = n,nord
                dx = dx + fdse0(k,nse)*akmat0(n,k)
            end do

            demos0(n,nse) = dx
        end do
    end do
end subroutine sderiv