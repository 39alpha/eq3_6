subroutine afcalc(actlg,af,afcnst,cdrs,jflag,jsflag,ndrs,ndrsmx,ndrsr,ns,nstmax,si,xlks)
    !! This subroutine computes the affinity function:
    !!   A = 2.303 RT log Q/K
    !! and the saturation index:
    !!   SI = RT log Q/K
    !! for the ns-th reaction. Note that this affinity is the affinity
    !! for the forward direction (A(+)), hence it has the same sign as
    !! the saturation index (SI).
    !! EQLIB/gafsir.f makes the same calculations for all reactions.
    !! This subroutine is called by:
    !!   EQLIB/betas.f
    !!   EQLIB/hpsat.f
    !!   EQ6/raff.f
    !!   EQ6/satchk.f
    !! Principal input:
    !!   actlg  = array of log activities of species
    !!   afcnst = the factor 2.303 RT
    !!   cdrs   = array of reaction coefficients
    !!   jflag  = flag array denoting whether species are in the active
    !!              basis set (jflag is not 30) or out of it (jflag = 30)
    !!   ndrs   = array of indices of species whose reaction coefficients
    !!              are in the cdrs array
    !!   ndrsr  = pointer array for the cdrs/ndrs arrays, denoting the
    !!              range of entries corresponding to a given reaction
    !!   ns     = index of the species for whose associated reaction
    !!              the affinity and saturation index are to be
    !!              calculated
    !!   xlks   = array of equilibrium constants
    !! Principal output:
    !!   af     = affinity
    !!   si     = saturation index (SI)
    implicit none

    ! Calling sequence variable declarations.
    integer :: ndrsmx
    integer :: nstmax

    integer :: jflag(nstmax)
    integer :: jsflag(nstmax)
    integer :: ndrs(ndrsmx)
    integer :: ndrsr(2,nstmax)
    integer :: ns

    real(kind=8) :: actlg(nstmax)
    real(kind=8) :: cdrs(ndrsmx)
    real(kind=8) :: xlks(nstmax)
    real(kind=8) :: af
    real(kind=8) :: afcnst
    real(kind=8) :: si

    ! Local variable declarations.
    integer :: n
    integer :: nr1
    integer :: nr2
    integer :: nse
    integer :: nt

    if (jsflag(ns) .eq. 2) then
        ! Species is suppressed.
        si = -9999999.
        af = -9999999.
        go to 999
    end if

    if (jflag(ns) .eq. 30) then
        ! Species is active, but is required to satisfy equilibrium.
        si = 0.
        af = 0.
        go to 999
    end if

    nr1 = ndrsr(1,ns)
    nr2 = ndrsr(2,ns)
    nt = nr2 - nr1 + 1

    if (nt .lt. 2) then
        ! Species has no reaction. Take it to be in equilibrium with
        ! itself.
        si = 0.
        af = 0.
        go to 999
    end if

    if (xlks(ns) .ge. 9999999.) then
        ! The equilibrium constant is unknown (set to pseudo-infinity).
        ! Set the SI and affinity to negative pseudo-infinity.
        si = -9999999.
        af = -9999999.
        go to 999
    end if

    ! Calculate the SI and affinity in the normal way.
    si = -xlks(ns)

    do n = nr1,nr2
        nse = ndrs(n)

        if (nse .eq. 0) then
            go to 200
        end if

        si = si + cdrs(n)*actlg(nse)
    end do

    af = afcnst*si
    go to 999

200 continue

    ! The species has no valid reaction. It is a detached auxiliary
    ! basis species. Take it to be in equilibrium with itself.
    si = 0.
    af = 0.

999 continue
end subroutine afcalc