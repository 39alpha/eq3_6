subroutine gfdho(aphi,bt,f,fp,fpp,fxi)
    !! This subroutine computes the Debye-Huckel function f and its
    !! derivatives with respect to ionic strength for the case of
    !! the Debye-Huckel-osmotic (DHO) model. This is the Debye-Huckel
    !! model in Pitzer's (1973, 1975) equations.
    !! This subroutine is called by:
    !!   EQLIBG/gcoeff.f
    !! Input:
    !!   aphi   = the Debye-Huckel A(phi) parameter
    !!   bt     = product of Debye-Huckel B and the hard core diameter,
    !!            usually taken as a constant (1.2)
    !!   fxi    = the ionic strength (the 2nd-order electrostatic
    !!              moment function I)
    !! Output:
    !!   f      = Debye-Huckel function f
    !!   fp     = f' (df/dI)
    !!   fpp    = f'' (d2f/dI2)
    implicit none

    ! Calling sequence variable declarations.
    real(kind=8) :: aphi
    real(kind=8) :: bt
    real(kind=8) :: f
    real(kind=8) :: fp
    real(kind=8) :: fpp
    real(kind=8) :: fxi

    ! Local variable declarations.
    real(kind=8) :: alvxp1
    real(kind=8) :: fc
    real(kind=8) :: fxisqt
    real(kind=8) :: vx
    real(kind=8) :: vxp1
    real(kind=8) :: vxp1sq

    fxisqt = sqrt(fxi)
    vx = bt*fxisqt
    vxp1 = vx + 1.
    vxp1sq = vxp1*vxp1
    alvxp1 = log(vxp1)

    fc = -4.*aphi/bt
    f = fc*fxi*alvxp1
    fp = fc*((vx / (2.*vxp1)) + alvxp1)
    fpp = (-aphi/fxisqt)*((1. + 2.*vxp1) / vxp1sq)
end subroutine gfdho