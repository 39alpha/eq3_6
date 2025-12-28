subroutine corrfd(delxi,dxsm11,fdlim,fdre0,fdre1,fdri0,fdri1,fdrr0,fdrr1,iodb,iopt,jreac,nodbmx,noptmx,nord,nordmx,noutpt,npts,nrct,nrctmx,nrd1mx,rirec0,rirec1,rreac0,rreac1,rrelr0,rrelr1)
    !! This subroutine computes finite differences for use in ODE
    !! corrector iteration. These finite differences are based at the
    !! new point stepped to, as opposed to the point stepped from.
    !! Finite differences based on the point stepped from are used to
    !! generate predictor functions. Those finite differences are
    !! computed by EQ6/stepfd.f and comprise a larger set. Here finite
    !! differences are computed only for rate functions.
    !! This subroutine is called by:
    !!   EQ6/path.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: nodbmx
    integer :: noptmx
    integer :: nordmx
    integer :: nrctmx
    integer :: nrd1mx

    integer :: noutpt

    integer :: iodb(nodbmx)
    integer :: iopt(noptmx)
    integer :: jreac(nrctmx)

    integer :: nord
    integer :: npts
    integer :: nrct

    real(kind=8) :: dxsm11(nrd1mx)
    real(kind=8) :: fdre0(nordmx,nrctmx)
    real(kind=8) :: fdre1(nordmx,nrctmx)
    real(kind=8) :: fdri0(nrd1mx)
    real(kind=8) :: fdri1(nrd1mx)
    real(kind=8) :: fdrr0(nrd1mx,nrctmx)
    real(kind=8) :: fdrr1(nrd1mx,nrctmx)

    real(kind=8) :: rreac0(nrctmx)
    real(kind=8) :: rreac1(nrctmx)
    real(kind=8) :: rrelr0(nrctmx)
    real(kind=8) :: rrelr1(nrctmx)

    real(kind=8) :: delxi
    real(kind=8) :: fdlim
    real(kind=8) :: rirec0
    real(kind=8) :: rirec1

    ! Local variable declarations.
    integer :: itrunc
    integer :: j
    integer :: k
    integer :: nmax
    integer :: nrc

    real(kind=8) :: dfx
    real(kind=8) :: dfxl
    real(kind=8) :: dfxu

    if (npts .eq. 1) then
        ! Zero all finite differences.
        nmax = nrd1mx*nrctmx
        call initaz(fdrr1,nmax)

        if (iopt(2) .gt. 0) then
            call initaz(fdri1,nrd1mx)

            nmax = nrd1mx*nrctmx
            call initaz(fdre1,nmax)
        end if

        go to 999
    end if

    ! Overflow protection is automatically activated in the code blocks
    ! below in order to run on VAX machines with a small exponent range
    ! (+/- 38) for real*8.
    itrunc = nord

    ! Compute finite differences for relative rates.
    do nrc = 1,nrct
        if (jreac(nrc).eq.1 .or. jreac(nrc).eq.2) then
            ! The relative rate must be a constant zero, because no
            ! reactant mass remains:
            !   jreac =  1: exhausted
            !   jreac =  2: saturated; any remaining reactant mass is
            !                 converted to the corresponding product
            !                 phase, so the "reactant" is effectively
            !                 exhausted.
            ! Hence the corresponding finite differences must also be
            ! zero. Avoid calculating finite differences from relative
            ! rate values that may by non-zero owing to convergence
            ! tolerances and the like.
            do j = 1,nord + 1
                fdrr1(j,nrc) = 0.
            end do
        else
            ! The relative rate is calculated for the reactant, which
            ! is actively reacting according to a rate law:
            !   jreac =  0: set to react
            !   jreac = -1: saturated, but the remaining reactant mass
            !                 continues to react irreversibly
            ! Hence calculate the finite differences.
            fdrr1(1,nrc) = (rrelr1(nrc) - rrelr0(nrc))/delxi

            do j = 2,nord + 1
                dfxu = fdlim*dxsm11(j)
                dfxl = -dfxu
                k = j - 1
                dfx = fdrr1(k,nrc) - fdrr0(k,nrc)
                dfx = min(dfxu,dfx)
                dfx = max(dfxl,dfx)

                if (abs(dfx) .ge. dfxu) then
                    itrunc = min(itrunc,k)
                end if

                fdrr1(j,nrc) = dfx/dxsm11(j)
            end do
        end if
    end do

    if (iopt(2) .gt. 0) then
        ! Compute finite differences for the inverse rate.
        fdri1(1) = (rirec1 - rirec0)/delxi

        do j = 2,nord + 1
            dfxu = fdlim*dxsm11(j)
            dfxl = -dfxu
            k = j - 1
            dfx = fdri1(k) - fdri0(k)
            dfx = min(dfxu,dfx)
            dfx = max(dfxl,dfx)

            if (abs(dfx) .ge. dfxu) then
                itrunc = min(itrunc,k)
            end if

            fdri1(j) = dfx/dxsm11(j)
        end do

        ! Compute finite differences for reaction rates.
        do nrc = 1,nrct
            if (jreac(nrc).eq.1 .or. jreac(nrc).eq.2) then
                ! The reaction rate must be a constant zero, because no
                ! reactant mass remains:
                !   jreac =  1: exhausted
                !   jreac =  2: saturated; any remaining reactant mass is
                !                 converted to the corresponding product
                !                 phase, so the "reactant" is effectively
                !                 exhausted.
                ! Hence the corresponding finite differences must also be
                ! zero. Avoid calculating finite differences from reaction
                ! rate values that may by non-zero owing to convergence
                ! tolerances and the like.
                do j = 1,nord + 1
                    fdre1(j,nrc) = 0.
                end do
            else
                ! The reaction rate is calculated for the reactant, which
                ! is actively reacting according to a rate law:
                !   jreac =  0: set to react
                !   jreac = -1: saturated, but the remaining reactant mass
                !                 continues to react irreversibly
                ! Hence calculate the finite differences.
                fdre1(1,nrc) = (rreac1(nrc) - rreac0(nrc))/delxi

                do j = 2,nord + 1
                    dfxu = fdlim*dxsm11(j)
                    dfxl = -dfxu
                    k = j - 1
                    dfx = fdre1(k,nrc) - fdre0(k,nrc)
                    dfx = min(dfxu,dfx)
                    dfx = max(dfxl,dfx)

                    if (abs(dfx) .ge. dfxu) then
                        itrunc = min(itrunc,k)
                    end if

                    if (j .le. nord) then
                        fdre1(j,nrc) = dfx/dxsm11(j)
                    end if
                end do
            end if
        end do
    end if

    if (itrunc .lt. nord) then
        if (iodb(1) .ge. 1) then
            write (noutpt,1000) itrunc
        end if

1000 format(/' * Note - (EQ6/corrfd) Cutting the order to ',i2,/7x,'to stay within the finite difference limit.')

        nord = itrunc
    end if

999 continue
end subroutine corrfd