subroutine ataylr(delxi,daffp0,nord,nordmx,npt,nptmax,affp0,affp)
    !! This subroutine evaluates Taylor's series expansions for phase
    !! affinities. These expansions are used to find phase boundaries
    !! at which new phases appear in the ES.
    !! Compare with:
    !!   EQ6/ptaylr.f
    !!   EQ6/rtaylr.f
    !!   EQ6/ztaylr.f
    !! This subroutine is called by:
    !!   EQ6/fpbnpp.f
    !! Principal input:
    !!   delxi  = the step size (in reaction progress)
    !!   daffp0 = the phase affinity derivative vector at the base point
    !!   nord   = the order of the truncated Taylor's series
    !!   nordmx = the maximum order of the truncated Taylor's series
    !!   npt    = the number of phases
    !!   nptmax = the maximum number of phases
    !!   affp0  = the phase affinity vector at the base point
    !! Principal output:
    !!   affp   = the phase affinity vector at the new point
    implicit none

    ! Calling sequence variable declarations.
    integer :: nordmx
    integer :: nptmax

    integer :: nord
    integer :: npt

    real(kind=8) :: affp(nptmax)
    real(kind=8) :: affp0(nptmax)
    real(kind=8) :: daffp0(nordmx,nptmax)

    real(kind=8) :: delxi

    ! Local variable declarations.
    integer :: n
    integer :: np

    real(kind=8) :: ax
    real(kind=8) :: dxp

    real(kind=8) :: fctrl

    do np = 1,npt
        ax = affp0(np)
        dxp = 1.

        do n = 1,nord
            dxp = dxp*delxi
            ax = ax + ( daffp0(n,np)/fctrl(n) )*dxp
        end do

        affp(np) = ax
    end do
end subroutine ataylr