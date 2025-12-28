subroutine ncmpvh(acflg,act,actlg,cdrs,cgxj,jflag,jsflag,losp,mosp,mtxj,nbasp,nbt,nbtmax,ndrs,ndrsmx,ndrsr,nrr1,nrr2,nstmax,xbar,xbarlg,xlks)
    !! This subroutine calculates part of the "expansion" of the basis
    !! set description regarding a site of a generic ion exchange
    !! phase. It assists EQLIB/ncmpve.h. The part of the expansion
    !! performed here does not take into account the need for an
    !! iterative process to complete the actual expansion; rather, it
    !! calculates a tentative part of the expansion. The iterative
    !! process is required (at least in the general case, for example
    !! Na+ for Ca++ following the Vanselow exchange model) because
    !! mtxj (the sum of the number of moles of the species on the site)
    !! is known only approximately at the start of the expansion process.
    !! This subroutine uses a tentative value of mtxj as an input.
    !! This subroutine is called by:
    !!   EQLIB/ncmpvh.f
    !! Principal input:
    !!   acflg  = array of logarithms of activity coefficients
    !!   actlg  = array of logarithms of species activities
    !!              (input consists of valid results for aqueous
    !!              basis species; output consists of the same
    !!              for exchanger species belonging to Vanselow
    !!              exchanger phases)
    !!   cdrs   = array of reaction coefficients
    !!   mosp   = array of numbers of moles of species
    !!              (input contains valid results for basis species
    !!              and carry-over values for non-basis species)
    !!   mtxj   = the sum of the number of moles of exchanger species
    !!              in the current site of the current generic ion
    !!              exchange phase
    !!   ndrs   = array parallel to cdrs giving the index of the
    !!              corresponding species
    !!   ndrsr  = array giving the range in the cdrs/ndrs arrays
    !!              corresonding to the reaction associated with a
    !!              given species
    !!   nrr1   = the start of the species index range for the
    !!              current site of the current ion exchange phase
    !!   nrr2   = the end of the species index range for the
    !!              current site of the current ion exchange phase
    !!   xlks   = array of equilibrium constants
    !! Principal output:
    !!   act    = array of species activities
    !!              (output consists of the subset for exchanger
    !!              species belonging to Vanselow exchanger phases)
    !!   actlg  = array of logarithms of species activities
    !!   mosp   = array of numbers of moles of species
    !!              (output consists of valid results for non-basis
    !!              species; valid results for basis species were
    !!              input and are retained)
    !!   losp   = array of logarithms of numbers of moles of species
    !!   xbar   = array of mole fractions of species
    !!   xbarlg = array of logarithms of mole fractions of species
    implicit none

    ! Calling sequence variable declarations.
    integer :: nbtmax
    integer :: ndrsmx
    integer :: nstmax

    integer :: jflag(nstmax)
    integer :: jsflag(nstmax)
    integer :: nbasp(nbtmax)
    integer :: ndrs(ndrsmx)
    integer :: ndrsr(2,nstmax)

    integer :: nbt
    integer :: nrr1
    integer :: nrr2

    real(kind=8) :: acflg(nstmax)
    real(kind=8) :: act(nstmax)
    real(kind=8) :: actlg(nstmax)
    real(kind=8) :: cdrs(ndrsmx)
    real(kind=8) :: losp(nstmax)
    real(kind=8) :: mosp(nstmax)
    real(kind=8) :: xbar(nstmax)
    real(kind=8) :: xbarlg(nstmax)
    real(kind=8) :: xlks(nstmax)

    real(kind=8) :: cgxj
    real(kind=8) :: mtxj

    ! Local variable declarations.
    integer :: n
    integer :: nb
    integer :: nr1
    integer :: nr2
    integer :: ns
    integer :: nss

    integer :: nbasis

    real(kind=8) :: cxs
    real(kind=8) :: lax
    real(kind=8) :: lx
    real(kind=8) :: lxx
    real(kind=8) :: mx
    real(kind=8) :: xx

    real(kind=8) :: texp
    real(kind=8) :: tlg

    ! Estimate the mole fractions and activities of the
    ! basis species for the current site.
    do ns = nrr1,nrr2
        nb = nbasis(nbasp,nbt,nbtmax,ns)

        if (nb.gt.0 .and. jsflag(ns).le.0) then
            xx = mosp(ns)/mtxj
            xbar(ns) = xx
            lxx = tlg(xx)
            xbarlg(ns) = lxx
            lax = cgxj*(lxx + acflg(ns))
            actlg(ns) = lax
            act(ns) = texp(lax)
        end if
    end do

    ! Estimate the mole fractions and activities of the
    ! non-basis species for the current site.
    do ns = nrr1,nrr2
        if (jflag(ns).eq.30 .and. jsflag(ns).le.0) then
            ! Calculate the activities and mole fractions for
            ! the exchanger species not in the basis set.
            nr1 = ndrsr(1,ns)
            nr2 = ndrsr(2,ns)
            cxs = cdrs(nr1)
            lax = -xlks(ns)

            do n = nr1 + 1,nr2
                nss = ndrs(n)
                lax = lax + cdrs(n)*actlg(nss)
            end do

            lax = -lax/cxs
            actlg(ns) = lax
            act(ns) = texp(lax)
            lxx = (lax/cgxj) - acflg(ns)
            xbarlg(ns) = lxx
            xx = texp(lxx)
            xbar(ns) = xx

            ! Also calculate the corresponding numbers of moles.
            mx = xx*mtxj
            lx = tlg(mx)
            mosp(ns) = mx
            losp(ns) = lx
        end if
    end do
end subroutine ncmpvh