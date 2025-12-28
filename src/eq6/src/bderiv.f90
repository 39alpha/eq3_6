subroutine bderiv(akmat0,dafrc0,fdar0,jreac,nord,nordmx,nrct,nrctmx,nrd1mx)
    !! This subroutine computes estimates of the derivatives of the
    !! affinities of reactants (afrc1) from the corresponding finite
    !! differences. Note that (dafrc0) = (akmat0)(fdar0).
    !! Compare with:
    !!   EQ6/aderiv.f
    !!   EQ6/pderiv.f
    !!   EQ6/rderiv.f
    !!   EQ6/zderiv.f
    !! This subroutine is called by:
    !!   EQ6/path.f
    !! Principal input:
    !! Principal output:
    implicit none

    integer :: nordmx
    integer :: nrctmx
    integer :: nrd1mx

    integer :: nord
    integer :: nrct

    integer :: jreac(nrctmx)

    real(kind=8) :: akmat0(nrd1mx,nrd1mx)
    real(kind=8) :: dafrc0(nordmx,nrctmx)
    real(kind=8) :: fdar0(nordmx,nrctmx)

    ! Local variable declarations.
    integer :: k
    integer :: n
    integer :: nrc

    real(kind=8) :: dx

    ! Recall the following jreac reactant status flag
    ! conventions:
    !   jreac =  0: set to react
    !   jreac =  1: exhausted
    !   jreac = -1: saturated, but the remaining reactant mass
    !                 continues to react irreversibly
    !   jreac =  2: saturated; the status of any remaining reactant
    !                 mass is changed to that of a product phase
    ! For jreac(nrc) = -1 or 2, the affinity is fixed at zero.
    ! However, the actual calculated affinity values may be non-zero
    ! owing to convergence tolerances. In EQ6/stepfd.f, the
    ! corresponding finite differences are set to zero. In the present
    ! subroutine, the corresponding derivatives are set to zero.
    do nrc = 1,nrct
        if (jreac(nrc).eq.-1 .or. jreac(nrc).eq.2) then
            do n = 1,nord
                dafrc0(n,nrc) = 0.
            end do
        else
            do n = 1,nord
                dx = 0.

                do k = n,nord
                    dx = dx + fdar0(k,nrc)*akmat0(n,k)
                end do

                dafrc0(n,nrc) = dx
            end do
        end if
    end do
end subroutine bderiv