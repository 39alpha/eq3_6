subroutine gbdot(acflgc,actwlc,adh,al10,azero,bdh,bdot,cco2,fxi,insgf,narn1,narn2,natmax,nstmax,omega,sigmam,tempk,xbrwlc,zchar,zchsq2)
    !! This subroutine computes activity coefficients of aqueous species
    !! using the B-dot equation and a set of related approximations.
    !! This model is suitable only for relatively dilute solutions
    !! (ionic strength no greater than one molal, though demonstrable
    !! inaccuracy appears at much lower values).
    !! The B-dot equation is an extended Debye-Huckel equation (Helgeson,
    !! 1969). It is applied only to electrically charged species. The
    !! activity coefficients of electrically neutral solute species are
    !! treated according to the value of the insgf flag. If it is -1,
    !! the activity coefficient is computed from the Drummond (1981)
    !! equation. This is a fit of the log activity coefficient of
    !! CO2(aq) as a function of temperature and ionic strength (ignoring
    !! ion pairing) in pure aqueous sodium chloride solutions. This is
    !! an appropriate treatment for neutral species that are nonpolar.
    !! If insgf is 0, the log activity coefficient is set to zero.
    !! The activity of water is computed from an equation which was
    !! derived from the B-dot equation using thermodynamic consistency
    !! relations, though it is actually consistent with that equation
    !! only for the case of all ions having the same size and an
    !! electrical charge of +1 or -1.
    !! This subroutine is called by:
    !!   EQLIBG/gcoeff.f
    !! Principal input:
    !!   adh    = Debye-Huckel A(gamma) parameter
    !!   azero  = hard core diameter array
    !!   bdh    = Debye-Huckel B parameter
    !!   bdot   = B-dot parameter
    !!   cco2   = coefficients of the Drummond (1981) equation
    !!   insgf  = flag array for electrically neutral solute species
    !!   omega  = water constant; ~55.51.
    !!   sigmam = sum of solute molalities
    !!   xbrwlc = log mole fraction of water
    !!   fxi    = the ionic strength (the 2nd-order electrostatic
    !!              moment function I)
    !!   zchar  = electrical charge array
    !!   zchsq2 = one-half the charge squared array
    !!  Principal output:
    !!   acflgc = log activity coefficient array
    !!   actwlc = log activity of water
    implicit none

    ! Calling sequence variable declarations.
    integer :: natmax
    integer :: nstmax

    integer :: insgf(natmax)
    integer :: narn1
    integer :: narn2

    real(kind=8) :: acflgc(nstmax)
    real(kind=8) :: azero(natmax)
    real(kind=8) :: cco2(5)
    real(kind=8) :: zchar(nstmax)
    real(kind=8) :: zchsq2(nstmax)
    real(kind=8) :: actwlc
    real(kind=8) :: adh
    real(kind=8) :: al10
    real(kind=8) :: bdh
    real(kind=8) :: bdot
    real(kind=8) :: fxi
    real(kind=8) :: omega
    real(kind=8) :: sigmam
    real(kind=8) :: tempk
    real(kind=8) :: xbrwlc

    ! Local variable declarations.
    integer :: na
    integer :: ns

    real(kind=8) :: acfco2
    real(kind=8) :: art
    real(kind=8) :: bdtfxi
    real(kind=8) :: brt
    real(kind=8) :: fxisqt
    real(kind=8) :: sga
    real(kind=8) :: sgx
    real(kind=8) :: xx
    real(kind=8) :: xxp1

    fxisqt = sqrt(fxi)

    ! Compute log lambda(w).
    xx = 4.0*bdh*fxisqt
    xxp1 = xx + 1.
    sgx = (3./(xx**3)) * ( xxp1 - (1./xxp1) - 2.*log(xxp1) )
    sga = sigmam/al10
    actwlc = ( -sga + 2.*adh*(fxi**1.5)*sgx/3. - bdot*fxi*fxi )/omega
    acflgc(narn1) = actwlc - xbrwlc

    ! Compute log gamma(i).
    ! Compute the log activity coefficient of CO2(aq), using the
    ! equation proposed by Drummond (1981, thesis, Penn State
    ! University).
    acfco2 = ( cco2(1) + cco2(2)*tempk + cco2(3)/tempk )*fxi - ( cco2(4) + cco2(5)*tempk )*(fxi/(fxi + 1.))
    acfco2 = acfco2/al10

    art = 2*adh*fxisqt
    brt = bdh*fxisqt
    bdtfxi = bdot*fxi

    do ns = narn1 + 1,narn2
        na = ns - narn1 + 1

        if (zchar(ns) .ne. 0.) then
            ! Have an ion. Use the B-dot equation.
            acflgc(ns)  = -( (art*zchsq2(ns)) / (1.0 + brt*azero(na)) )    + bdtfxi
        else
            ! Have a neutral species.
            if (insgf(na) .le. -1) then
                ! For nonpolar species, using the result for CO2(aq) from
                ! the Drummond (1981) equation.
                acflgc(ns) = acfco2
            else
                ! For polar species, assign a log gamma value of zero.
                acflgc(ns) = 0.
            end if
        end if
    end do
end subroutine gbdot