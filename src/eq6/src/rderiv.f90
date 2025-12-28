subroutine rderiv(akmat0,drer0,drir0,fdri0,fdrr0,jreac,nord,nrct,nrctmx,nrd1mx)
    !! This subroutine transforms finite differences of the rates of
    !! irreversible reactions into the corresponding derivatives.
    !! Note that (drer0) = (akmat0)(fdrr0). The inverse rate is treated
    !! similarly; i.e., (drir0) = (akmat0)(fdri0).
    !! Compare with:
    !!   EQ6/aderiv.f
    !!   EQ6/bderiv.f
    !!   EQ6/pderiv.f
    !!   EQ6/zderiv.f
    !! This subroutine is called by:
    !!   EQ6/chksti.f
    !!   EQ6/path.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: nrctmx
    integer :: nrd1mx

    integer :: jreac(nrctmx)

    integer :: nord
    integer :: nrct

    real(kind=8) :: akmat0(nrd1mx,nrd1mx)
    real(kind=8) :: drer0(nrd1mx,nrctmx)
    real(kind=8) :: drir0(nrd1mx)
    real(kind=8) :: fdri0(nrd1mx)
    real(kind=8) :: fdrr0(nrd1mx,nrctmx)

    ! Local variable declarations.
    integer :: k
    integer :: n
    integer :: nrc

    real(kind=8) :: dxx

    ! Calculate the derivatives of the inverse rate.
    do n = 1,nord
        dxx = 0.

        do k = n,nord
            dxx = dxx + fdri0(k)*akmat0(n,k)
        end do

        drir0(n) = dxx
    end do

    ! Calculate the derivatives of the rates of irreversible reactions.
    do nrc = 1,nrct
        if (jreac(nrc).eq.0 .or. jreac(nrc).eq.-1)  then
            do n = 1,nord
                dxx = 0.

                do k = n,nord
                    dxx = dxx + fdrr0(k,nrc)*akmat0(n,k)
                end do

                drer0(n,nrc) = dxx
            end do
        end if
    end do
end subroutine rderiv