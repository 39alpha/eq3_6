subroutine rtaylr(delxi,drer0,drir0,jreac,nord,nrct,nrctmx,nrd1mx,rirec0,rirecp,rrelr0,rrelrp)
    !! This subroutine evaluates Taylor's series expansions for the
    !! inverse rate and the relative rates of all irreversible reactions.
    !! Compare with:
    !!   EQ6/ataylr.f
    !!   EQ6/ptaylr.f
    !!   EQ6/ztaylr.f
    !! This subroutine is called by:
    !!   EQ6/path.f
    !! Principal input:
    !!   delxi  = the step size (in reaction progress)
    !!   drer0  = the relative rate derivative vector at the base point
    !!   drir0  = the inverse rate derivative vector at the base point
    !!   nord   = the order of the truncated Taylor's series
    !!   nrd1mx = the maximum order of the truncated Taylor's series + 1
    !!   nrct   = the number of reactants
    !!   nrctmx = the maximum number of reactants
    !!   rirec0 = the inverse rate at the base point
    !!   rrelr0 = the relative rate vector at the base point
    !! Principal output:
    !!   rirecp = the inverse rate at the new point
    !!   rrelrp = the relative rate vector at the new point
    implicit none

    ! Calling sequence variable declarations.
    integer :: nrctmx
    integer :: nrd1mx

    integer :: jreac(nrctmx)

    integer :: nord
    integer :: nrct

    real(kind=8) :: drer0(nrd1mx,nrctmx)
    real(kind=8) :: drir0(nrd1mx)
    real(kind=8) :: rrelr0(nrctmx)
    real(kind=8) :: rrelrp(nrctmx)

    real(kind=8) :: delxi
    real(kind=8) :: rirecp
    real(kind=8) :: rirec0

    ! Local variable declarations.
    integer :: n
    integer :: nrc

    real(kind=8) :: dxp
    real(kind=8) :: rx

    real(kind=8) :: fctrl

    rirecp = rirec0
    dxp = 1.

    do n = 1,nord
        dxp = dxp*delxi
        rirecp = rirecp + ( drir0(n)/fctrl(n) )*dxp
    end do

    do nrc = 1,nrct
        rx = 0.

        if (jreac(nrc) .eq. 0) then
            rx = rrelr0(nrc)
            dxp = 1.

            do n = 1,nord
                dxp = dxp*delxi
                rx = rx + ( drer0(n,nrc)/fctrl(n) )*dxp
            end do
        end if

        rrelrp(nrc) = rx
    end do
end subroutine rtaylr