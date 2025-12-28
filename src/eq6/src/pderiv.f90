subroutine pderiv(akmat0,demop0,fdpe0,nord,nordmx,npet,npetmx,nrd1mx)
    !! This subroutine computes estimates of the derivatives of the
    !! numbers of moles of the phases present in the equilibrium
    !! system. These derivatives (demop0) are computed from the
    !! corresponding finite differences (fdpe0). The relation is:
    !! (demop0) = (akmat0)(fdpe0).
    !! Compare with:
    !!   EQ6/aderiv.f
    !!   EQ6/bderiv.f
    !!   EQ6/rderiv.f
    !!   EQ6/sderiv.f
    !!   EQ6/zderiv.f
    !! This subroutine is called by:
    !!   EQ6/path.f
    !! Principal input:
    !! Principal output:
    implicit none

    integer :: nordmx
    integer :: npetmx
    integer :: nrd1mx

    integer :: nord
    integer :: npet

    real(kind=8) :: akmat0(nrd1mx,nrd1mx)
    real(kind=8) :: demop0(nordmx,npetmx)
    real(kind=8) :: fdpe0(nordmx,npetmx)

    ! Local variable declarations.
    integer :: k
    integer :: n
    integer :: npe

    real(kind=8) :: dx

    do npe = 1,npet
        do n = 1,nord
            dx = 0.

            do k = n,nord
                dx = dx + fdpe0(k,npe)*akmat0(n,k)
            end do

            demop0(n,npe) = dx
        end do
    end do
end subroutine pderiv