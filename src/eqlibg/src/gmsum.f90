subroutine gmsum(conc,musumw,nmut,nmutmx,nmux,nstmax,pmu,uspec)
    !! This subroutine computes the following triple sum used in Pitzer's
    !! equations:
    !!   SUM(ijk) mu(ijk)*m(i)*m(j)*m(k)
    !! This sum is used to calculate the activity coefficient
    !! of water.
    !! This subroutine is called by:
    !!   EQLIBG/gcoeff.f
    !! Principal input:
    !!   conc   = array of concentration values
    !!   nmu    = number of mu values
    !!   nmux   = array of triples of species indices
    !!   pmu    = array of mu values
    !! Principal output:
    !!   musumw =  the sum: SUM(ijk) mu(ijk)*m(i)*m(j)*m(k)
    implicit none

    ! Calling sequence variable declarations.
    integer :: nmutmx
    integer :: nstmax

    integer :: nmux(3,nmutmx)
    integer :: nmut

    character(len=48) :: uspec(nstmax)

    real(kind=8) :: conc(nstmax)
    real(kind=8) :: pmu(nmutmx)
    real(kind=8) :: musumw

    ! Local variable declarations.
    integer :: nmu
    integer :: ns1
    integer :: ns2
    integer :: ns3

    real(kind=8) :: cx

    musumw = 0.

    do nmu = 1,nmut
        ns1 = nmux(1,nmu)
        ns2 = nmux(2,nmu)
        ns3 = nmux(3,nmu)
        cx = 6.

        if (ns1.eq.ns2 .and. ns1.ne.ns3) then
            cx = 3.
        end if

        if (ns2.eq.ns3 .and. ns1.ne.ns3) then
            cx = 3.
        end if

        if (ns1.eq.ns2 .and. ns1.eq.ns3) then
            cx = 1.
        end if

        musumw = musumw + cx*pmu(nmu)*conc(ns1)*conc(ns2)*conc(ns3)
    end do
end subroutine gmsum