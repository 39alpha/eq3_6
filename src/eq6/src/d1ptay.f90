subroutine d1ptay(delxi,demop0,d1emp1,nord,nordmx,npet,npetmx)
    !! This subroutine evaluates the first derivative of Taylor's series
    !! expansions for the number of moles of phases in the ES. These
    !! expansions are used to find points of reaction progress at which
    !! such variables maximize. At such points, for example in the
    !! fluid-centered flow-through open system model,it may be necessary
    !! to transfer mass of a corresponding phase from the ES to the PRS.
    !! See also:
    !!   EQ6/d2ptay.f
    !!   EQ6/ztaylr.f
    !! This subroutine is called by:
    !!   EQ6/fpbflo.f
    !! Principal input:
    !!   delxi  = the step size (in reaction progress)
    !!   demop0 = the number of moles derivative vector for phases in
    !!              the ES at the base point
    !!   nord   = the order of the truncated Taylor's series
    !!   nordmx = the maximum order of the truncated Taylor's series
    !!   npet   = the number of phases in the ES
    !!   npetmx = the maximum number of phases in the ES
    !! Principal output:
    !!   d1emp1 = the first derivative of the number of moles vector
    !!              for phases in the ES at at the new point
    implicit none

    ! Calling sequence variable declarations.
    integer :: nordmx
    integer :: npetmx

    integer :: nord
    integer :: npet

    real(kind=8) :: demop0(nordmx,npetmx)
    real(kind=8) :: d1emp1(npetmx)

    real(kind=8) :: delxi

    ! Local variable declarations.
    integer :: j
    integer :: jk
    integer :: npe

    real(kind=8) :: dxp
    real(kind=8) :: d1p

    real(kind=8) :: fctrl

    if (nord .gt. 0) then
        do npe = 1,npet
            d1p = demop0(1,npe)

            if (nord .gt. 1) then
                jk = nord - 1
                dxp = 1.

                do j = 1,jk
                    dxp = dxp*delxi
                    d1p = d1p + ( demop0(j + 1,npe)/fctrl(j) )*dxp
                end do
            end if

            d1emp1(npe) = d1p
        end do
    end if
end subroutine d1ptay