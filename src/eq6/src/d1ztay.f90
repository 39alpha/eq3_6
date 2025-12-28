subroutine d1ztay(delxi,dzvc0,d1zvc1,kdim,kmax,nord,nrd1mx)
    !! This subroutine computes the Taylor's series expansions for the
    !! first derivatives of the master variables. These expansions are
    !! used to find points of reaction progress at which such variables
    !! maximize. At such points, for example in the fluid-centered
    !! flow-through open system model,it may be necessary to transfer
    !! mass of a corresponding phase from the ES to the PRS.
    !! See also:
    !!   EQ6/d2ztay.f
    !!   EQ6/ztaylr.f
    !! This subroutine is called by:
    !!   EQ6/eqcalc.f
    !! Principal input:
    !!   delxi  = the step size (in reaction progress)
    !!   dzvc0  = the dz/d(xi) vector at the base point
    !!   kdim   = the number of elements in the z vector
    !!   kmax   = the maximum number of elements in the z vector
    !!   nord   = the order of the truncated Taylor's series
    !!   nrd1mx = the maximum order of the truncated Taylor's series + 1
    !! Principal output:
    !!   d1zvc1 = the dz/d(xi) vector at the new point
    implicit none

    ! Calling sequence variable declarations.
    integer :: kmax
    integer :: nrd1mx

    integer :: kdim
    integer :: nord

    real(kind=8) :: dzvc0(nrd1mx,kmax)
    real(kind=8) :: d1zvc1(kmax)

    real(kind=8) :: delxi

    ! Local variable declarations.
    integer :: j
    integer :: jk
    integer :: kcol

    real(kind=8) :: dxp
    real(kind=8) :: d1z

    real(kind=8) :: fctrl

    if (nord .gt. 0) then
        do kcol = 1,kdim
            d1z = dzvc0(1,kcol)

            if (nord .gt. 1) then
                jk = nord - 1
                dxp = 1.

                do j = 1,jk
                    dxp = dxp*delxi
                    d1z = d1z + ( dzvc0(j + 1,kcol)/fctrl(j) )*dxp
                end do
            end if

            d1zvc1(kcol) = d1z
        end do
    end if
end subroutine d1ztay