subroutine gdavie(acflgc,actwlc,adh,al10,fxi,narn1,narn2,nstmax,omega,sigmam,xbrwlc,zchsq2)
    !! This subroutine computes activity coefficients of aqueous species
    !! using the Davies (1961) equation. The activity of water is
    !! computed from an expression that was derived from the Davies
    !! equation using thermodynamic consistency relations.
    !! This subroutine is called by:
    !!    EQLIBG/gcoeff.f
    !! Principal input:
    !!   adh    = Debye-Huckel A(gamma) parameter
    !!   omega  = water constant; ~55.51.
    !!   sigmam = sum of solute molalities
    !!   xbrwlc = log mole fraction of water
    !!   fxi    = the ionic strength (the 2nd-order electrostatic
    !!              moment function I)
    !!   zchsq2 = one-half the charge squared array
    !! Principal output:
    !!   acflgc = array of log activity coefficients
    !!   actwlc = log activity of water
    implicit none

    ! Calling sequence variable declarations.
    integer :: nstmax

    integer :: narn1
    integer :: narn2

    real(kind=8) :: acflgc(nstmax)
    real(kind=8) :: zchsq2(nstmax)
    real(kind=8) :: actwlc
    real(kind=8) :: adh
    real(kind=8) :: al10
    real(kind=8) :: fxi
    real(kind=8) :: omega
    real(kind=8) :: sigmam
    real(kind=8) :: xbrwlc

    ! Local variable declarations.
    integer :: ns

    real(kind=8) :: factor
    real(kind=8) :: fxisqt
    real(kind=8) :: sga
    real(kind=8) :: sgx
    real(kind=8) :: xxp1

    fxisqt = sqrt(fxi)
    xxp1 = 1. + fxisqt

    ! Compute log lambda(w).
    sgx = (3./(fxisqt**3))*( xxp1 - (1./xxp1) - 2.*log(xxp1) )
    sga = sigmam/al10
    actwlc = ( -sga + 2.*adh*(fxi**1.5)*sgx/3. - 0.2*adh*fxi*fxi )/omega
    acflgc(narn1) = actwlc - xbrwlc

    ! Compute log gamma(i).
    factor = (fxisqt/xxp1) - 0.2*fxi
    factor = -2*adh*factor

    do ns = narn1 + 1,narn2
        acflgc(ns) = zchsq2(ns)*factor
    end do
end subroutine gdavie