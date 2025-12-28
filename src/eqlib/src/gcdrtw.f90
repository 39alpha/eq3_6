subroutine gcdrtw(cdrs,cdrtw,narn1,narn2,ndrs,ndrsmx,ndrsr,nelect,no2gaq,nst,nstmax)
    !! This subroutine computes the cdrtw array. Each element of this
    !! array contains the sum of the reaction coefficients of the
    !! aqueous solute species in the corresponding reaction. The sum
    !! represents the number of times that the number of moles of
    !! solvent water is implicitly represented in the reaction via
    !! solute molalities. Note that generic ion exchanger species do
    !! not count in these calculations, even though molalities can
    !! be calculated for such species. That is because their
    !! thermodynamic activities are defined by mole fractions.
    !! Do not confuse this array with the cdrw array computed by
    !! EQLIB/gcdrw.f. The present array is really only used by
    !! EQ6. It presently exists in EQ3NR only for the sake of
    !! consistency.
    !! This subroutine is called by:
    !!   EQLIB/absswa.f
    !!   EQ3NR/eq3nr.f
    !!   EQ6/path.f
    !!   EQ6/absswb.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: ndrsmx
    integer :: nstmax

    integer :: ndrs(ndrsmx)
    integer :: ndrsr(2,nstmax)

    integer :: narn1
    integer :: narn2
    integer :: nelect
    integer :: no2gaq
    integer :: nst

    real(kind=8) :: cdrs(ndrsmx)
    real(kind=8) :: cdrtw(nstmax)

    ! Local variable declarations.
    integer :: n
    integer :: nn
    integer :: nr1
    integer :: nr2
    integer :: ns

    real(kind=8) :: cx

    do ns = 1,nst
        nr1 = ndrsr(1,ns)
        nr2 = ndrsr(2,ns)
        cx = 0.

        do n = nr1,nr2
            nn = ndrs(n)

            if (nn.gt.narn1 .and. nn.le.narn2) then
                if (nn.ne.no2gaq .and. nn.ne.nelect) then
                    cx = cx + cdrs(n)
                end if
            end if
        end do

        cdrtw(ns) = cx
    end do
end subroutine gcdrtw