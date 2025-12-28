subroutine gsgsm(conc,dpslm,na,natmax,nsltmx,nstmax,nsxi,nsxmax,nsxx,pslm,slsum,slsump,uspec)
    !! This subroutine computes the following first-order sums used
    !! in Pitzer's equations:
    !!   SUM(j) S-lambda(ij)*m(j)
    !!   SUM(j) S-lambda'(ij)*m(j)
    !! This subroutine is called by:
    !!   EQLIBG/gcoeff.f
    !! Principal input:
    !!   conc   = array of species concentrations
    !!   na     = aqueous species index (i refers to ns, where ns is
    !!              the corresponding species index)
    !!   nsxi   = range pointer array into the nsxx array:
    !!              nsxi(1,na) and nsxi(2,na) are the first and last
    !!              values of the second subscript (nsx) of the nsxx
    !!              array for entries pertaining to the na-th
    !!              aqueous species (ns-th species)
    !!   nsxx   = pointer array:
    !!              nsxx(1,nsx) = the species index of the second
    !!              species in the nsl-th pair, where
    !!              nsxx(2,nsx) = nsl
    !!   pslm   = array of S-lambda functions
    !!   dpslm  = array of ionic strength derivatives of the
    !!              S-lambda functions
    !! Principal input:
    !!   slsum  = the sum: SUM(j) S-lambda(ij)*m(j)
    !!   slsump = the sum: SUM(j) S-lambda'(ij)*m(j)
    implicit none

    ! Calling sequence variable declarations.
    integer :: natmax
    integer :: nsltmx
    integer :: nstmax
    integer :: nsxmax

    integer :: nsxi(2,natmax)
    integer :: nsxx(2,nsxmax)
    integer :: na

    character(len=48) :: uspec(nstmax)

    real(kind=8) :: conc(nstmax)
    real(kind=8) :: dpslm(2,nsltmx)
    real(kind=8) :: pslm(nsltmx)
    real(kind=8) :: slsum
    real(kind=8) :: slsump

    ! Local variable declarations.
    integer :: ixf
    integer :: ixl
    integer :: nsl
    integer :: nsx
    integer :: ns2

    slsum = 0.
    slsump = 0.
    ixf = nsxi(1,na)
    ixl = nsxi(2,na)

    do nsx = ixf,ixl
        ns2 = nsxx(1,nsx)
        nsl = nsxx(2,nsx)
        slsum = slsum + pslm(nsl)*conc(ns2)
        slsump = slsump + dpslm(1,nsl)*conc(ns2)
    end do
end subroutine gsgsm