subroutine gabswx(beta,ibswx,iindx1,kbt,kmax,nbt,nbtmax)
    !! This subroutine supports automatic basis switching as a
    !! pre-Newton-Raphson optimization technique. It resolves conflicts
    !! in the ibswx array, which contains candidate switches. This array
    !! could call for the same non-basis species to be switched with
    !! more than one basis species. The approach here is to use the size
    !! of the associated relative mass balance residual to any resolve
    !! conflicts. This subroutine is somewhat similar to EQLIB/gbfac.f,
    !! which resolves conflicts affecting continued fraction
    !! calculations.
    !! This subroutine is called by:
    !!   EQLIB/absswa.f
    !! Principal input:
    !!   ibswx = array of indices of species not in the active
    !!             master set that are candidates for switching into
    !!             the basis set
    !!   beta  = array of normalized Newton-Raphson residual functions
    !!   kbt   = number of species in the active basis set
    !! Principal output:
    !!   ibswx = input array modified to remove conflicts.
    implicit none

    ! Calling sequence variable declarations.
    integer :: kmax
    integer :: nbtmax

    integer :: ibswx(nbtmax)
    integer :: iindx1(kmax)

    integer :: kbt
    integer :: nbt

    real(kind=8) :: beta(kmax)

    ! Local variable declarations.
    integer :: kcol
    integer :: kk
    integer :: krow
    integer :: nb
    integer :: nbig
    integer :: nb1
    integer :: nb2
    integer :: ns1
    integer :: ns2

    real(kind=8) :: bbig
    real(kind=8) :: frac

    ! Find the largest residual among potential basis switching cases.
    nbig = 0
    bbig = 0.

    do kcol = 1,kbt
        nb = iindx1(kcol)

        if (ibswx(nb) .gt. 0) then
            if (beta(kcol) .gt. bbig) then
                nbig = nb
                bbig = beta(kcol)
            end if
        end if
    end do

    ! Eliminate candidates not involving largest residuals.
    if (bbig .le. 0.5) then
        ! If the largest residual is not very large, wipe out all other
        ! proposed switches.
        do nb = 1,nbt
            if (nb .ne. nbig) then
                ibswx(nb) = 0
            end if
        end do

        go to 999
    else
        ! If the largest residual is rather large, wipe out all other
        ! proposed switches with residuals that are not within two
        ! orders of magnitude of the largest residual.
        do kcol = 1,kbt
            nb = iindx1(kcol)

            if (ibswx(nb) .gt. 0) then
                frac = beta(kcol)/bbig

                if (frac .le. 1.e-2) then
                    ibswx(nb) = 0
                end if
            end if
        end do
    end if

    ! Eliminate any remaining conflicts.
    do krow = 1,kbt - 1
        nb1 = iindx1(krow)
        ns1 = ibswx(nb1)

        if (ns1 .gt. 0) then
            kk = krow + 1

            do kcol = kk,kbt
                nb2 = iindx1(kcol)
                ns2 = ibswx(nb2)

                if (ns2 .eq. ns1) then
                    if (beta(krow) .gt. beta(kcol)) then
                        ibswx(nb2) = 0
                    else
                        ibswx(nb1) = 0
                        go to 25
                    end if
                end if
            end do
        end if

25 continue
    end do

999 continue
end subroutine gabswx