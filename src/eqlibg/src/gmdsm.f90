subroutine gmdsm(conc,musum,na,natmax,nmutmx,nmxi,nmxmax,nmxx,ns,nstmax,pmu,uspec)
    !! This subroutine computes the following second order sum used
    !! in Pitzer's equations:
    !!   SUM(jk) mu(ijk)*m(j)*m(k)
    !! where i is represented by the aqueous species index na and
    !! overall species index ns.
    !! This subroutine is called by:
    !!   EQLIBG/gcoeff.f
    !! Principal input:
    !!   conc   = array of concentration values
    !!   na     = aqueous species index (i refers to ns, where ns is
    !!              the corresponding species index)
    !!   nmxi   = range pointer array into the nmxx array:
    !!              nmxi(1,na) and nmxi(2,na) are the first and last
    !!              values of the second subscript (nmx) of the nmxx
    !!              array for entries pertaining to the na-th
    !!              aqueous species (ns-th species)
    !!   nmxx   = pointer array:
    !!              nmxx(1,nmx) = the species index of the second
    !!              species in the nmu-th triplet, nmxx(2,nmx) is the
    !!              species index of the third species, and
    !!              nmxx(3,nmx) = nmu
    !!   ns     = species index corresponding to the aqueous species
    !!              index na
    !!   pmu    = array of mu values
    !!   uspec  = array of species names
    !! Principal output:
    !!   musum  = the sum: SUM(jk) mu(ijk)*m(j)*m(k)
    implicit none

    ! Calling sequence variable declarations.
    integer :: natmax
    integer :: nmutmx
    integer :: nmxmax
    integer :: nstmax

    integer :: nmxi(2,natmax)
    integer :: nmxx(3,nmxmax)
    integer :: na
    integer :: ns

    character(len=48) :: uspec(nstmax)

    real(kind=8) :: conc(nstmax)
    real(kind=8) :: pmu(nmutmx)
    real(kind=8) :: musum

    ! Local variable declarations.
    integer :: ifx
    integer :: ilx
    integer :: nmu
    integer :: nmx
    integer :: ns2
    integer :: ns3

    real(kind=8) :: cx

    musum = 0.
    ifx = nmxi(1,na)
    ilx = nmxi(2,na)

    if (ilx .ge. ifx) then
        do nmx = ifx,ilx
            ns2 = nmxx(1,nmx)
            ns3 = nmxx(2,nmx)
            nmu = nmxx(3,nmx)
            cx = 2.

            if (ns2 .eq. ns3) then
                cx = 1.
            end if

            musum = musum + cx*pmu(nmu)*conc(ns2)*conc(ns3)
        end do
    end if
end subroutine gmdsm