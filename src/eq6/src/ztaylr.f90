subroutine ztaylr(delxi,dzvc0,kdim,kmax,km1,kxt,nord,nrd1mx,qztayl,zklogu,zvclg0,zvclg1,zvec0,zvec1)
    !! This subroutine calculates new values for algebraic master
    !! variables (the z vector elements) from the finite-difference-
    !! based truncated Taylor's series. If qztayl = .true., change
    !! limits are applied to the results. The purpose of these limits
    !! is to assist the hybrid Newton-Raphson iteration either to
    !! converge or to generate useful divergence diagnostics.
    !! Compare with:
    !!   EQ6/ataylr.f
    !!   EQ6/ptaylr.f
    !!   EQ6/rtaylr.f
    !!  See also:
    !!   EQ6/d1ztay.f
    !!   EQ6/d2ztay.f
    !! This subroutine is called by:
    !!   EQ6/eqshel.f
    !!   EQ6/ldlxrc.f
    !!   EQ6/path.f
    !! Principal input:
    !!   delxi  = the step size (in reaction progress)
    !!   dzvc0  = the dz/d(xi) vector at the base point
    !!   kdim   = the number of elements in the z vector
    !!   kmax   = the maximum number of elements in the z vector
    !!   nord   = the order of the truncated Taylor's series
    !!   nrd1mx = the maximum order of the truncated Taylor's series + 1
    !!   qztayl = flag to apply change limits
    !!   zvclg0 = the log z vector at the base point
    !!   zvec0  = the z vector at the base point
    !! Principal output:
    !!   zvec1  = the z vector at the new point
    !!   zvclg1 = the log z vector at the new point
    implicit none

    ! Calling sequence variable declarations.
    integer :: kmax
    integer :: nrd1mx

    integer :: kdim
    integer :: km1
    integer :: kxt
    integer :: nord

    logical :: qztayl

    real(kind=8) :: dzvc0(nrd1mx,kmax)
    real(kind=8) :: zvclg0(kmax)
    real(kind=8) :: zvclg1(kmax)
    real(kind=8) :: zvec0(kmax)
    real(kind=8) :: zvec1(kmax)

    real(kind=8) :: delxi
    real(kind=8) :: zklogu

    ! Local variable declarations.
    integer :: kcol
    integer :: n

    real(kind=8) :: dxp
    real(kind=8) :: zx
    real(kind=8) :: zx0
    real(kind=8) :: zxl
    real(kind=8) :: zxu

    real(kind=8) :: fctrl
    real(kind=8) :: tlg

    ! Compute the expansions from the Taylor's series.
    do kcol = 1,kdim
        zx = zvec0(kcol)
        dxp = 1.

        do n = 1,nord
            dxp = dxp*delxi
            zx = zx + ( dzvc0(n,kcol)/fctrl(n) )*dxp
        end do

        zvec1(kcol) = zx

        ! Compute the corresponding logarithmic variable. Provide
        ! protection if zx is negative. If this is the case, the
        ! logarithmic variable won't be used.
        if (zx .lt. 0) then
            zx = 0.
        end if

        zvclg1(kcol) = tlg(zx)
    end do

    if (qztayl) then
        ! Apply change limits.
        do kcol = 1,kdim
            zx0 = zvec0(kcol)
            zxl = 1.e-20*zx0
            zxu = 1.e+20*zx0
            zx = zvec1(kcol)
            zx = max(zx,zxl)
            zx = min(zx,zxu)

            if (zx .lt. 0.) then
                zx = 0.
            end if

            zvec1(kcol) = zx
            zvclg1(kcol) = tlg(zx)
        end do

        do kcol = km1,kxt
            if (zvclg0(kcol) .lt. zklogu) then
                zvclg1(kcol) = zvclg0(kcol)
                zvec1(kcol) = zvec0(kcol)
            end if
        end do
    end if
end subroutine ztaylr