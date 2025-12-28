subroutine clress(csts,iindx1,ipndx1,jpflag,jsflag,kdim,kmax,km1,kmt,kx1,kxt,loph,losp,moph,mosp,mtb,mtbaq,nbt,nbtmax,nptmax,nstmax,nsts,nstsmx,nstsr,ufixf,uzvec1,zvec1,zvclg1)
    !! This subroutine clears equilibrium system (ES) solids. These
    !! solids include pure minerals and solid solutions, but not fictive
    !! fugacity-fixing minerals.
    !! This subroutine is called by:
    !!   EQ6/eq6.f
    !!   EQ6/path.f
    !! Principal input:
    !!   csts   = array of stoichiometric coefficients for mass
    !!              balance expressions
    !!   losp   = array of log number of moles of species variables
    !!   mosp   = array of number of moles of species
    !!   mtb    = array of total numbers of moles of basis species
    !!   mtbaq  = array of total numbers of moles of basis species
    !!              in the aqueous solution
    !!   uzvec1 = array of master variable names
    !!   zvec1  = array of master variables
    !!   zvclg1 = array of log master variables
    !! Principal output:
    !!   losp   = array of log number of moles of species variables
    !!              (modified)
    !!   mosp   = array of number of moles of species (modified)
    !!   mtb    = array of total numbers of moles of basis species
    !!              (modified)
    !!   zvec1  = array of master variables (modified)
    !!   zvclg1 = array of log master variables (modified)
    implicit none

    ! Calling sequence variable declarations.
    integer :: kmax
    integer :: nbtmax
    integer :: nptmax
    integer :: nstmax
    integer :: nstsmx

    integer :: iindx1(kmax)
    integer :: ipndx1(kmax)
    integer :: jpflag(nptmax)
    integer :: jsflag(nstmax)
    integer :: nsts(nstsmx)
    integer :: nstsr(2,nstmax)
    integer :: kdim
    integer :: km1
    integer :: kmt
    integer :: kx1
    integer :: kxt
    integer :: nbt

    character(len=48) :: uzvec1(kmax)
    character(len=8) :: ufixf

    real(kind=8) :: csts(nstsmx)
    real(kind=8) :: loph(nptmax)
    real(kind=8) :: losp(nstmax)
    real(kind=8) :: moph(nptmax)
    real(kind=8) :: mosp(nstmax)
    real(kind=8) :: mtb(nbtmax)
    real(kind=8) :: mtbaq(nbtmax)
    real(kind=8) :: zvec1(kmax)
    real(kind=8) :: zvclg1(kmax)

    ! Local variable declarations.
    integer :: k
    integer :: kcol
    integer :: n
    integer :: nb
    integer :: nr1
    integer :: nr2
    integer :: np
    integer :: ns

    ! Clear the phase and species status flags for the items to be
    ! cleared from the ES.
    do kcol = km1,kmt
        if (uzvec1(kcol)(1:5) .ne. ufixf(1:5)) then
            ! Do not have a fictive fugacity-fixing mineral.
            np = ipndx1(kcol)
            ns = iindx1(kcol)
            jpflag(np) = 0
            jsflag(ns) = 0
            moph(np) = 0.
            loph(np) = -99999.
            mosp(ns) = 0.
            losp(ns) = -99999.
        end if
    end do

    do kcol = kx1,kxt
        np = ipndx1(kcol)
        ns = iindx1(kcol)
        jpflag(np) = 0
        jsflag(ns) = 0
        moph(np) = 0.
        loph(np) = -99999.
        mosp(ns) = 0.
        losp(ns) = -99999.
    end do

    ! Recompute the mass balance totals to account for the cleared
    ! phases.
    do nb = 1,nbt
        mtb(nb) = mtbaq(nb)
    end do

    do kcol = km1,kmt
        if (uzvec1(kcol)(1:5) .eq. ufixf(1:5)) then
            ! Have a fictive fugacity-fixing mineral. Adjust the
            ! mass balances to retain it in the system.
            ! Retain it.
            ns = iindx1(kcol)
            nr1 = nstsr(1,ns)
            nr2 = nstsr(2,ns)

            do n = nr1,nr2
                nb = nsts(n)
                mtb(nb) = mtb(nb) + csts(n)*zvec1(kcol)
            end do
        end if
    end do

    k = km1 - 1

    do kcol = km1,kmt
        if (uzvec1(kcol)(1:5) .eq. ufixf(1:5)) then
            ! Have a fictive fugacity-fixing mineral. Retain it in
            ! the indexing scheme for the cleared system.
            k = k + 1
            zvec1(k) = zvec1(kcol)
            zvclg1(k) = zvclg1(kcol)
            ipndx1(k) = ipndx1(kcol)
            iindx1(k) = iindx1(kcol)
        end if
    end do

    kmt = k
    kx1 = kmt + 1

    do kcol = kx1,kxt
        ! Clear all other minerals, both pure minerals and solid
        ! solutions.
        zvec1(kcol) = 0.
        zvclg1(kcol) = -99999.
        iindx1(kcol) = 0
        ipndx1(kcol) = 0
    end do

    kxt = kmt
    kdim = kmt
end subroutine clress