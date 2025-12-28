subroutine d2ptay(delxi,demop0,d2emp1,nord,nordmx,npet,npetmx)
    !! This subroutine computes the Taylor's series expansion for the
    !! second derivative of the number of moles of the npe-th phase
    !! in the ES. This second derivative is used to test whether or not
    !! a critical point corresponds to a maximum.
    !! See also:
    !!   EQ6/d1ptay.f
    !!   EQ6/ptaylr.f
    !! This subroutine is called by:
    !!   EQ6/fpbflo.f
    !! Principal input:
    !!   delxi  = the step size (in reaction progress)
    !!   demop0 = the number of moles derivative vector for phases in
    !!              the ES at the base point
    !!   nord   = the order of the truncated Taylor's series
    !!   nordmx = the maximum order of the truncated Taylor's series
    !!   npe    = the ES index of the phase for which the second
    !!              derivative of the number of moles is to be computed
    !!   npetmx = the maximum number of phases in the ES
    !! Principal output:
    !!   d2emp1 = the second derivative of the number of moles vector
    !!              for phases in the ES at the new point
    implicit none

    ! Calling sequence variable declarations.
    integer :: nordmx
    integer :: npetmx

    integer :: nord
    integer :: npet

    real(kind=8) :: demop0(nordmx,npetmx)
    real(kind=8) :: d2emp1(npetmx)

    real(kind=8) :: delxi

    ! Local variable declarations.
    integer :: j
    integer :: jk
    integer :: npe

    real(kind=8) :: dxp
    real(kind=8) :: d2p

    real(kind=8) :: fctrl

    do npe = 1,npet
        d2emp1(npe) = 0.

        if (nord .gt. 1) then
            d2p = demop0(2,npe)

            if (nord .gt. 2) then
                jk = nord - 2
                dxp = 1.

                do j = 1,jk
                    dxp = dxp*delxi
                    d2p = d2p + ( demop0(j + 2,npe )/fctrl(j) )*dxp
                end do
            end if

            d2emp1(npe) = d2p
        end if
    end do
end subroutine d2ptay