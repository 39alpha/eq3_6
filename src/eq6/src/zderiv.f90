subroutine zderiv(akmat0,dzvc0,fdzv0,kdim,kmax,nord,nrd1mx)
    !! This subroutine computes estimates of the derivatives of the
    !! master algebraic variables (z vector) from the corresponding
    !! finite differences. Note that (dzvc0) = (akmat0)(fdzv0).
    !! Compare with:
    !!   EQ6/aderiv.f
    !!   EQ6/bderiv.f
    !!   EQ6/pderiv.f
    !!   EQ6/rderiv.f
    !! This subroutine is called by:
    !!   EQ6/path.f
    !! Principal input:
    !! Principal output:
    implicit none

    integer :: kmax
    integer :: nrd1mx

    integer :: kdim
    integer :: nord

    real(kind=8) :: akmat0(nrd1mx,nrd1mx)
    real(kind=8) :: dzvc0(nrd1mx,kmax)
    real(kind=8) :: fdzv0(nrd1mx,kmax)

    ! Local variable declarations.
    integer :: k
    integer :: kcol
    integer :: n

    real(kind=8) :: dx

    do kcol = 1,kdim
        do n = 1,nord
            dx = 0.

            do k = n,nord
                dx = dx + fdzv0(k,kcol)*akmat0(n,k)
            end do

            dzvc0(n,kcol) = dx
        end do
    end do
end subroutine zderiv