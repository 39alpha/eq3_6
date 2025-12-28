real(kind=8) function coefdr(cdrs,ndrs,ndrsmx,ndrsr,nse,ns,nstmax)
    !! This subroutine finds the coefficient of the nse-th species
    !! in the reaction for the destruction of the ns-th species.
    !! This subroutine is called by:
    !!   EQLIB/elim.f
    !!   EQLIB/ncmpex.f
    !!   EQLIB/switch.f
    !!   EQLIB/swtchb.f
    !!   EQLIB/swtchk.f
    !!   EQ3NR/arrset.f
    !!   EQ3NR/arrsim.f
    !!   EQ3NR/dawfix.f
    !!   EQ3NR/balcon.f
    !!   EQ3NR/betas.f
    !!   EQ3NR/eq3nr.f
    !!   EQ3NR/matrix.f
    !!   EQ3NR/scripx.f
    !!   EQ6/balcmz.f
    !!   EQ6/jgibbs.f
    !!   EQ6/matrxz.f
    !!   EQ6/mincsp.f
    !!   EQ6/scanlm.f
    !!   EQ6/tstrdx.f
    !! Principal input:
    !!   cdrs   = array of reaction coefficients
    !!   ndrs   = array of indices of the species corresponding to the
    !!              coefficients in the cdrs array
    !!   ndrsr  = array giving the range in the cdrs and ndrs arrays
    !!              containing the reaction for a given species
    !! Principal output:
    !!   coefdr = the reaction coefficient of the nse-th species in
    !!              the reaction for the ns-th species
    implicit none

    ! Calling sequence variable declarations.
    integer :: ndrsmx
    integer :: nstmax

    integer :: ndrs(ndrsmx)
    integer :: ndrsr(2,nstmax)
    integer :: nse
    integer :: ns

    real(kind=8) :: cdrs(ndrsmx)

    ! Local variable declarations.
    integer :: n
    integer :: n1
    integer :: n2

    coefdr = 0.

    n1 = ndrsr(1,ns)
    n2 = ndrsr(2,ns)

    do n = n1,n2
        if (nse .eq. ndrs(n)) then
            coefdr = cdrs(n)
            go to 999
        end if
    end do

999 continue
end function coefdr