subroutine integr(delxi,dlxrct,drer0,nord,nrc,nrctmx,nrd1mx,rrelr0)
    !! This subroutine integrates a Taylor's series for the relative
    !! rate of a reaction to calculate the advancement in the
    !! corresponding reaction progress variable (dlxrct).
    !! This subroutine is called by:
    !!   EQ6/reacts.f
    !! Principal input:
    !!   delxi  = overall reaction progress variable
    !!   drer0  = array of derivatives of relative rates for individual
    !!              irreversible reactions
    !!   nord   = order of the finite difference approximation
    !! Principal output:
    !!   dlxrct = array of increments of reaction progress for
    implicit none

    ! Calling sequence variable declarations.
    integer :: nrctmx
    integer :: nrd1mx

    integer :: nord
    integer :: nrc

    real(kind=8) :: drer0(nrd1mx,nrctmx)
    real(kind=8) :: rrelr0(nrctmx)
    real(kind=8) :: delxi
    real(kind=8) :: dlxrct

    ! Local variable declarations.
    integer :: n

    real(kind=8) :: dxp

    real(kind=8) :: fctrl

    dlxrct = rrelr0(nrc)*delxi

    if (nord .gt. 0) then
        dxp = delxi

        do n = 1,nord
            dxp = dxp*delxi
            dlxrct = dlxrct + ( drer0(n,nrc)/fctrl(n + 1) )*dxp
        end do
    end if
end subroutine integr