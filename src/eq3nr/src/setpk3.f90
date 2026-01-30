subroutine setpk3(electr,iindx1,jflag,jflgi,kbt,kdim,kmax,kmt,kprs,kwater,kxt,mtb,mtbi,mtbaq,mtbaqi,narn1,narn2,nbasp,nbaspd,nbaspi,nbti,nbtmax,ndrsrd,nern1,nern2,nobswt,nstmax,ntitl,ntitl2,ntitmx,omeglg,press,pressi,scamas,sigzi,tempc,tempci,ubmtbi,uobsw,uspec,utitl,utitl2,uzveci,uzvec1,zvclgi,zvclg1)
    !! This subroutine sets up certain variables and arrays for writing
    !! on the pickup file. This subroutine must be called prior to
    !! writing a pickup file.
    !! This subroutine is called by:
    !!   EQ3NR/eq3nr.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: kmax
    integer :: nbtmax
    integer :: nstmax
    integer :: ntitmx

    integer :: iindx1(kmax)
    integer :: jflag(nstmax)
    integer :: jflgi(nbtmax)
    integer :: nbasp(nbtmax)
    integer :: nbaspd(nbtmax)
    integer :: nbaspi(nbtmax)
    integer :: ndrsrd(2,nstmax)

    integer :: kbt
    integer :: kdim
    integer :: kmt
    integer :: kprs
    integer :: kwater
    integer :: kxt
    integer :: narn1
    integer :: narn2
    integer :: nern1
    integer :: nern2
    integer :: nbti
    integer :: nobswt
    integer :: ntitl
    integer :: ntitl2

    character(len=80) :: utitl(ntitmx)
    character(len=80) :: utitl2(ntitmx)
    character(len=48) :: ubmtbi(nbtmax)
    character(len=48) :: uobsw(2,nbtmax)
    character(len=48) :: uspec(nstmax)
    character(len=48) :: uzveci(kmax)
    character(len=48) :: uzvec1(kmax)

    real(kind=8) :: mtb(nbtmax)
    real(kind=8) :: mtbi(nbtmax)
    real(kind=8) :: mtbaq(nbtmax)
    real(kind=8) :: mtbaqi(nbtmax)
    real(kind=8) :: zvclgi(kmax)
    real(kind=8) :: zvclg1(kmax)

    real(kind=8) :: electr
    real(kind=8) :: omeglg
    real(kind=8) :: press
    real(kind=8) :: pressi
    real(kind=8) :: scamas
    real(kind=8) :: sigzi
    real(kind=8) :: tempc
    real(kind=8) :: tempci

    ! Local variable declarations.
    integer :: kcol
    integer :: krow
    integer :: n
    integer :: nb
    integer :: nbi
    integer :: nerr
    integer :: nsi
    integer :: ns1
    integer :: ns2
    integer :: nt1

    real(kind=8) :: lscama

    real(kind=8) :: tlg

    nerr = 0

    ! Old title.
    ntitl2 = ntitl

    do n = 1,ntitl
        utitl2(n) = utitl(n)
    end do

    ! Original temperature.
    tempci = tempc

    ! Original pressure.
    pressi = press

    ! Index limits.
    kprs = 0

    kmt = kdim
    kxt = kdim
    nbti = kbt

    ! Species for which mass balances are defined.
    !   ubmtbi = names of the data file basis species for which mass
    !              balances are defined
    !   jflgi  = jflag input for basis species
    !      0 = Retain as an active basis species
    !     30 = Convert to a dependent species; fold the mass balance
    !            total for this species into the mass balance totals
    !            of basis species which remain active
    do krow = 1,kbt
        nb = iindx1(krow)
        ns1 = nbaspd(nb)
        ns2 = nbasp(nb)
        nbi = krow
        nsi = nbaspi(nb)
        ubmtbi(nbi) = uspec(nsi)

        if (ns1.ge.narn1 .and. ns1.le.narn2) then
            mtbi(nbi) = scamas*mtb(nb)
            mtbaqi(nbi) = scamas*mtbaq(nb)
            nt1 = ndrsrd(2,ns1) - ndrsrd(1,ns1) + 1

            if (nt1 .lt. 2) then
                ! Have a strict basis species associated with the current
                ! mass balance.
                jflgi(nbi) = jflag(ns2)
            else
                ! Have an auxiliary basis species associated with the current
                ! mass balance.
                jflgi(nbi) = 30
            end if
        else if (ns1.ge.nern1 .and. ns1.le.nern2) then
            mtbi(nbi) = mtb(nb)
            mtbaqi(nbi) = mtbaq(nb)
            jflgi(nbi) = jflag(ns2)
        else
            mtbi(nbi) = mtb(nb)
            mtbaqi(nbi) = mtbaq(nb)
            jflgi(nbi) = jflag(ns2)
        end if
    end do

    ! Electrical balance.
    electr = scamas*sigzi

    ! Ordinary basis switching directives.
    n = 0

    do kcol = 1,kbt
        nb = iindx1(kcol)
        ns1 = nbaspd(nb)
        ns2 = nbasp(nb)

        if (ns1 .ne. ns2) then
            n = n + 1
            uobsw(1,n) = uspec(ns1)
            uobsw(2,n) = uspec(ns2)
        end if
    end do

    nobswt = n

    ! Matrix column variables and corresponding values.
    lscama = tlg(scamas)
    zvclg1(kwater) = omeglg

    do kcol = 1,kdim
        nb = iindx1(kcol)
        ns2 = nbasp(nb)
        uzveci(kcol) = uzvec1(kcol)

        if (ns2.ge.narn1 .and. ns2.le.narn2) then
            zvclgi(kcol) = zvclg1(kcol) + lscama
        else if (ns2.ge.nern1 .and. ns2.le.nern2) then
            zvclgi(kcol) = zvclg1(kcol)
        else
            zvclgi(kcol) = zvclg1(kcol)
        end if
    end do

    if (nerr .gt. 0) then
        stop
    end if
end subroutine setpk3
