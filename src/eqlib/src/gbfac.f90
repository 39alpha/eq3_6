subroutine gbfac(beta,bfac,efac,iindx1,kbt,kmax,nbt,nbtmax,nfac)
    !! This subroutine calculates the bfac array, which is used in making
    !! continued fraction corrections. It resolves conflicts when the
    !! same aqueous species dominates more than one mass balance (the
    !! species dominating a given mass balance is determined by
    !! EQLIB/fdomsp). The continued fraction algorithm can be applied
    !! to the master species associated with only one mass balance in
    !! such a set, otherwise, oscillatory behavior will occur. In each
    !! set of mass balances with a common dominant species, this
    !! subroutine finds the mass balance with the greatest residual and
    !! completes the calculation of its bfac factor by doing the
    !! appropriate exponentiation. It sets bfac to unity for the others
    !! in the set.
    !! This subroutine is called by:
    !!   EQ3NR/arrset.f
    !!   EQ6/optmzr.f
    !! Principal input:
    !!   nfac   = array of indices of dominant aqueous species
    !!   beta   = array of normalized Newton-Raphson residual functions
    !!   efac   = array of corresponding reciprocal stoichiometric
    !!              weights
    !!   kbt    = the number of active basis species
    !!   nbt    = the number of active basis species
    !! Principal output:
    !!   bfac   = array (in terms of matrix indexing) of the quantity
    !!              used in making a continued fraction correction
    !!              (e.g., conc (new) = conc (old) / bfac ); in general,
    !!              bfac = ( beta + 1. )**efac, but bfac may be set to
    !!              unity for a species in resolving a conflict
    implicit none

    ! Calling sequence variable declarations.
    integer :: kmax
    integer :: nbtmax

    integer :: iindx1(kmax)
    integer :: nfac(nbtmax)
    integer :: kbt
    integer :: nbt

    real(kind=8) :: beta(kmax)
    real(kind=8) :: bfac(nbtmax)
    real(kind=8) :: efac(nbtmax)

    ! Local variable declarations.
    integer :: krow
    integer :: nb
    integer :: nb1
    integer :: nb2

    real(kind=8) :: bpx1

    ! Calculate the bfac correction factors: bfac = (beta + 1)**efac.
    do krow = 1,kbt
        nb = iindx1(krow)
        bpx1 = beta(krow) + 1.

        if (bpx1 .le. 0.) then
            ! Protect against a singularity in case (beta + 1) is not
            ! a positive number. The value set in such a case allows a
            ! generous enough one-step increase in the value of the
            ! corrected concentration, number of moles, etc. (20 orders
            ! of magnitude). A more restrictive limit may be applied
            ! later when the actual correction is made.
            bfac(nb) = 1.e-20
        else
            ! Calculate bfac from the usual formula.
            bfac(nb) = bpx1**efac(nb)
        end if
    end do

    ! Eliminate conflicts.
    do nb1 = 1,nbt - 1
        if (nfac(nb1) .gt. 0) then
            do nb2 = nb1 + 1,nbt
                if (nfac(nb2) .eq. nfac(nb1)) then
                    if (bfac(nb1) .gt. bfac(nb2)) then
                        nfac(nb2) = 0
                        bfac(nb2) = 1.
                    else
                        nfac(nb1) = 0
                        bfac(nb1) = 1.
                        go to 15
                    end if
                end if
            end do
        end if

15 continue
    end do
end subroutine gbfac