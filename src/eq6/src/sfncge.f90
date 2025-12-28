subroutine sfncge(delxi,xval0,xtargv,dxval0,nord,nrd1mx,resx)
    !! This subroutine computes a general search function for
    !! EQ6/search.f. This search function is a residual defined as the
    !! difference between the calculated value of a function (xvalc)
    !! described by a truncated Taylor's series and a target value
    !! (xtargv).
    !! This subroutine is called by:
    !!   EQ6/search.f
    !! Principal input:
    !!   delxi  = step size
    !!   xtargv = target value
    !!   xval0  = the function value at delxi = 0
    !!   dxval0 = the derivatives of this function at delxi = 0
    !!   nord   = the order of the truncated Taylor's series
    !! Principal output:
    !!   resx   = value of residual function
    implicit none

    ! Calling sequence variable declarations.
    integer :: nrd1mx

    integer :: nord

    real(kind=8) :: dxval0(nrd1mx)

    real(kind=8) :: delxi
    real(kind=8) :: resx
    real(kind=8) :: xtargv
    real(kind=8) :: xval0

    real(kind=8) :: fctrl

    ! Local variable declarations.
    integer :: n

    real(kind=8) :: dxp
    real(kind=8) :: xvalc

    xvalc = xval0

    if (nord .gt. 0) then
        dxp = 1.

        do n = 1,nord
            dxp = dxp*delxi
            xvalc = xvalc + ( dxval0(n)/fctrl(n) )*dxp
        end do
    end if

    resx = xvalc - xtargv
end subroutine sfncge