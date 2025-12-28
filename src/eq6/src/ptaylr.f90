subroutine ptaylr(delxi,demop0,emop0,emop,nord,nordmx,npet,npetmx)
    !! This subroutine evaluates Taylor's series expansions for the
    !! number of moles of phases in the ES. These expansions are used
    !! to find phase boundaries at which phases disappear from the ES.
    !! Compare with:
    !!   EQ6/ataylr.f
    !!   EQ6/rtaylr.f
    !!   EQ6/ztaylr.f
    !! See also:
    !!   EQ6/d1ptay.f
    !!   EQ6/d2ptay.f
    !! This subroutine is called by:
    !!   EQ6/fpbdpp.f
    !!   EQ6/fpbflo.f
    !! Principal input:
    !!   delxi  = the step size (in reaction progress)
    !!   demop0 = the number of moles derivative vector for phases in
    !!              the ES at the base point
    !!   nord   = the order of the truncated Taylor's series
    !!   nordmx = the maximum order of the truncated Taylor's series
    !!   npet   = the number of phases in the ES
    !!   npetmx = the maximum number of phases in the ES
    !!   emop0  = the number of moles vector for phases in the ES at
    !!              at the base point
    !! Principal output:
    !!   emop   = the number of moles vector for phases in the ES at
    !!              at the new point
    implicit none

    ! Calling sequence variable declarations.
    integer :: nordmx
    integer :: npetmx

    integer :: nord
    integer :: npet

    real(kind=8) :: emop(npetmx)
    real(kind=8) :: emop0(npetmx)
    real(kind=8) :: demop0(nordmx,npetmx)

    real(kind=8) :: delxi

    ! Local variable declarations.
    integer :: n
    integer :: npe

    real(kind=8) :: dxp
    real(kind=8) :: ex

    real(kind=8) :: fctrl

    do npe = 1,npet
        ex = emop0(npe)
        dxp = 1.

        do n = 1,nord
            dxp = dxp*delxi
            ex = ex + ( demop0(n,npe)/fctrl(n) )*dxp
        end do

        emop(npe) = ex
    end do
end subroutine ptaylr