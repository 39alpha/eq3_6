subroutine fbassw(jcsort,jflag,mosp,narn1,narn2,nse,nsi,nsj,nstmax,weight,wsi)
    !! This subroutine attempts to find a candidate species to switch
    !! switch into the active basis set. The function here is similar
    !! to that EQLIB/fdomsp.f, which finds the species that dominates
    !! a mass balance.
    !! As presently contructed, this routine considers switches only
    !! among the set of aqueous species. This condition may be relaxed
    !! in the future.
    !! This subroutine is called by:
    !!   EQLIB/abswpk.f
    !!   EQ6/absswb.f
    !! Principal input:
    !!   jcsort = array of aqueous species indices in order of
    !!              increasing concentration
    !!   jflag  = jflag array defining constraint types imposed on
    !!              aqueous species
    !!   mosp   = array of numbers of moles of species
    !!   narn1  = index of the first aqueous species
    !!   narn2  = index of the last aqueous species
    !!   nse    = the data file basis species defining the current
    !!            mass balance
    !!   nsj    = the active basis species for the same mass balance;
    !!            nsj may be nse
    !!   weight = stoichiometric weighting factor
    !! Principal output:
    !!   nsi    = index of the aqeuous species that would represent the
    !!              optimal basis switch
    !!   wsi    = weight(nsi)
    implicit none

    ! Calling sequence variable declarations.
    integer :: nstmax

    integer :: jcsort(nstmax)
    integer :: jflag(nstmax)

    integer :: narn1
    integer :: narn2
    integer :: nse
    integer :: nsi
    integer :: nsj

    real(kind=8) :: mosp(nstmax)
    real(kind=8) :: weight(nstmax)

    real(kind=8) :: wsi

    ! Local variable declarations.
    integer :: ns
    integer :: nsc
    integer :: nss

    real(kind=8) :: ap0
    real(kind=8) :: ap1
    real(kind=8) :: m0
    real(kind=8) :: m1
    real(kind=8) :: p0
    real(kind=8) :: p1
    real(kind=8) :: rx
    real(kind=8) :: w0
    real(kind=8) :: w1

    ! Note that nsi = 0 is equivalent to nsi = nsj (the current
    ! active basis species associated with the present mass balance).
    ! Returning no candidate means keeping the existing basis species.
    nsi = 0
    wsi = 0.

    ! Make sure that the nse-th species is an aqueous species.
    if (nse.lt.narn1 .or. nse.gt.narn2) then
        go to 999
    end if

    ! Start by testing against the current active basis species.
    m0 = mosp(nsj)
    w0 = weight(nsj)
    p0 = w0*m0
    ap0 = abs(p0)

    ! Loop over all aqueous species, testing in order of decreasing
    ! number of moles. The current active basis species is not tested
    ! against itself in this loop because jflag(nsj) is not 30.
    ! Therefore, the following code cannot return nsi equal to nsj. If
    ! the corresponding data file basis species is not the current
    ! active basis species for the present mass balance, then jflag(nse)
    ! will be 30, and this species will be tested as a potential
    ! candidate.
    do nsc = narn1,narn2
        nss = narn2 + narn1 - nsc
        ns = jcsort(nss)
        w1 = weight(ns)
        m1 = mosp(ns)

        if (m1 .le. 0.) then
            go to 110
        end if

        if (w1 .ne. 0.) then
            if (jflag(ns) .eq. 30) then
                p1 = w1*m1
                ap1 = abs(p1)

                if (ap1 .gt. ap0) then
                    ! Have found a new leading candidate.
                    m0 = m1
                    w0 = w1
                    p0 = p1
                    ap0 = ap1
                    nsi = ns
                    go to 100
                end if
            end if
        end if

        ! Stop the search if the mole number ratio is now so high that
        ! it becomes unreasonable to expect the ratio of the weighting
        ! factors to possibly overcome this so as to yield ap1 > ap0.
        ! Note that the possibility of m1 being zero has been tested
        ! above, so the ratio calculation below should be safe.
        rx = m0/m1

        if (rx .gt. 100.) then
            go to 110
        end if

100 continue
    end do

110 continue
    wsi = w0

999 continue
end subroutine fbassw