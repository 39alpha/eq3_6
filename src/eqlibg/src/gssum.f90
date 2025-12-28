subroutine gssum(conc,dpslm,fxi,nslt,nsltmx,nslx,nstmax,pslm,ssumw,uspec)
    !! This subroutine computes the following second order sum used
    !! in Pitzer's equations:
    !!   SUM(ij) [ S-lambda(ij) + I*S-lambda'(ij) ]*m(i)*m(j)
    !! This sum is used to calculate the activity coefficient
    !! of water.
    !! This subroutine is called by:
    !!   EQLIBG/gcoeff.f
    !! Principal input:
    !!   conc   = array of species concentrations
    !!   fxi    = the ionic strength (the 2nd-order electrostatic
    !!              moment function I)
    !!   nslt   = number of S-lambda functions
    !!   nslx   = array of S-lambda species index pairs
    !!   pslm   = array of S-lambda functions
    !!   dpslm  = array of ionic strength derivatives of S-lambda
    !!            functions
    !! Principal output:
    !!   ssumw  = the sum:
    !!              SUM(ij) [ S-lambda(ij) + I*S-lambda'(ij) ]*m(i)*m(j)
    implicit none

    ! Calling sequence variable declarations.
    integer :: nsltmx
    integer :: nstmax

    integer :: nslx(2,nsltmx)
    integer :: nslt

    character(len=48) :: uspec(nstmax)

    real(kind=8) :: conc(nstmax)
    real(kind=8) :: dpslm(2,nsltmx)
    real(kind=8) :: pslm(nsltmx)
    real(kind=8) :: fxi
    real(kind=8) :: ssumw

    ! Local variable declarations.
    integer :: nsl
    integer :: ns1
    integer :: ns2

    real(kind=8) :: cx

    ssumw = 0.

    do nsl = 1,nslt
        ns1 = nslx(1,nsl)
        ns2 = nslx(2,nsl)
        cx = 2.

        if (ns1 .eq. ns2) then
            cx = 1.
        end if

        ssumw = ssumw  + cx*(pslm(nsl) + fxi*dpslm(1,nsl))*conc(ns1)*conc(ns2)
    end do
end subroutine gssum