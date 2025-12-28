subroutine d2ztay(delxi,dzvc0,d2zvc1,kcol,kmax,nord,nrd1mx)
    !! This subroutine computes the Taylor's series expansion for the
    !! second derivative of the kcol-th master variables. This second
    !! derivative is used to test whether or not a critical point
    !! corresponds to a maximum.
    !! See also:
    !!   EQ6/d1ztay.f
    !!   EQ6/ztaylr.f
    !! This subroutine is called by:
    !!   None
    !! Principal input:
    !!   delxi  = the step size (in reaction progress)
    !!   dzvc0  = the dz/d(xi) vector at the base point
    !!   kcol   = the index of the master variable whose second
    !!              derivative is to be calculated
    !!   kmax   = the maximum number of elements in the z vector
    !!   nord   = the order of the truncated Taylor's series
    !!   nrd1mx = the maximum order of the truncated Taylor's series + 1
    !! Principal output:
    !!   d2zvc1 = the d2z/d(xi)2 element for the kcol-th master
    !!              variable
    implicit none

    ! Calling sequence variable declarations.
    integer :: kmax
    integer :: nrd1mx

    integer :: kcol
    integer :: nord

    real(kind=8) :: dzvc0(nrd1mx,kmax)

    real(kind=8) :: delxi
    real(kind=8) :: d2zvc1

    ! Local variable declarations.
    integer :: j
    integer :: jk

    real(kind=8) :: dxp
    real(kind=8) :: d2z

    real(kind=8) :: fctrl

    d2zvc1 = 0.

    if (nord .gt. 1) then
        d2z = dzvc0(2,kcol)

        if (nord .gt. 2) then
            jk = nord - 2
            dxp = 1.

            do j = 1,jk
                dxp = dxp*delxi
                d2z = d2z + ( dzvc0(j + 2,kcol )/fctrl(j) )*dxp
            end do
        end if

        d2zvc1 = d2z
    end if
end subroutine d2ztay