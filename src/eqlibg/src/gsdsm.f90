subroutine gsdsm(conc,dpslm,nslt,nsltmx,nslx,nstmax,spsum,spsump,uspec)
    !! This subroutine computes the following second order sums used
    !! in Pitzer's equations:
    !!   spsum:  SUM(jk) S-lambda'(jk)*m(j)*m(k)
    !!   spsump: SUM(jk) S-lambda''(jk)*m(j)*m(k)
    !! Here jk always refers to a cation-anion pair. Note that:
    !!   spsum  = 2 * SUM(ca) B'(ca)m(a)m(c)
    !!   spsump = 2 * SUM(ca) B''(ca)m(a)m(c)
    !! This subroutine is called by:
    !!   EQLIBG/gcoeff.f
    !!  Principal input:
    !!   conc   = array of concentration values
    !!   nslt   = number of S-lambda functions
    !!   nslx   = array of S-lambda species index pairs
    !!   dpslm  = array of ionic strength derivatives of S-lambda
    !!              functions
    !! Principal output:
    !!   spsum  = the sum: SUM(jk) S-lambda'(jk)*m(j)*m(k)
    !!   spsump = the sum: SUM(jk) S-lambda''(jk)*m(j)*m(k)
    implicit none

    ! Calling sequence variable declarations.
    integer :: nsltmx
    integer :: nstmax

    integer :: nslx(2,nsltmx)
    integer :: nslt

    character(len=48) :: uspec(nstmax)

    real(kind=8) :: conc(nstmax)
    real(kind=8) :: dpslm(2,nsltmx)

    real(kind=8) :: spsum
    real(kind=8) :: spsump

    ! Local variable declarations.
    integer :: nsl
    integer :: ns2
    integer :: ns3

    real(kind=8) :: cp

    spsum = 0.
    spsump = 0.

    do nsl = 1,nslt
        ns2 = nslx(1,nsl)
        ns3 = nslx(2,nsl)
        cp = conc(ns2)*conc(ns3)
        spsum = spsum + dpslm(1,nsl)*cp
        spsump = spsump + dpslm(2,nsl)*cp
    end do

    spsum = 2.*spsum
    spsump = 2.*spsump
end subroutine gsdsm