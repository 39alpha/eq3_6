subroutine gfdhc(adh,bt,fdhc,fdhcp,fdhcpp,fxi)
    !! This subroutine computes the Debye-Huckel function f and its
    !! derivatives with respect to ionic strength for the case of the
    !! Debye-Huckel-charging (DHC) model.
    !! This subroutine is called by:
    !!   EQLIBG/gcoeff.f
    !! Input:
    !!   adh    = the Debye-Huckel A(gamma) parameter
    !!   bt     = product of Debye-Huckel B and the hard core diameter
    !!   fxi    = the ionic strength (the 2nd-order electrostatic
    !!              moment function I)
    !! Output:
    !!   fdhc   = Debye-Huckel function f (DHC model)
    !!   fdhcp  = f' (df/dI)
    !!   fdhcpp = f'' (d2f/dI2)
    implicit none

    ! Calling sequence variable declarations.
    real(kind=8) :: adh
    real(kind=8) :: bt
    real(kind=8) :: fdhc
    real(kind=8) :: fdhcp
    real(kind=8) :: fdhcpp
    real(kind=8) :: fxi

    ! Local variable declarations.
    real(kind=8) :: fxisqt
    real(kind=8) :: vx
    real(kind=8) :: vxp1

    fxisqt = sqrt(fxi)
    vx = bt*fxisqt
    vxp1 = vx + 1.

    fdhc = -6.*adh*(vx*(vx-2.) + 2.*log(vxp1))/(bt**3)
    fdhcp = -6.*adh*fxisqt/vxp1
    fdhcpp = -3.*adh/(fxisqt*vxp1*vxp1)
end subroutine gfdhc