subroutine gcoeff(abar,acflgc,actwlc,adh,adhh,adhv,al10,aphi,azero,a3bar,a3bars,bdh,bdhh,bdhv,bdot,bdoth,bdotv,cco2,conc,delam,dgpit,dpelm,dpslm,dselm,elam,fje,fxi,gpit,ielam,ifcphi1,ifcphi2,ifnnn,ifn2n,ifpsi1,ifpsi2,ifzeta,ilcphi1,ilcphi2,ilnnn,iln2n,ilpsi1,ilpsi2,ilzeta,insgf,iopg,ipbtmx,izmax,jcsort,nalpha,napmax,napt,narn1,narn2,natmax,nazmmx,nazpmx,nchlor,nhydr,nmut,nmutmx,nmux,nmxi,nmxmax,nmxx,nopgmx,noutpt,nslt,nsltmx,nslx,nstmax,nsxi,nsxmax,nsxx,nttyo,omega,palpha,pelm,pmu,press,pslamn,pslm,qhawep,qpit75,selm,sigmam,tempk,uspec,xbarwc,xbrwlc,zchar,zchsq2,zchcu6)
    !! This subroutine computes activity coefficients of aqueous species
    !! using various models. The model used is determined by the option
    !! flag iopg(1):
    !!    -1  Davies equation
    !!     0  B-dot equation plus others
    !!     1  Pitzer's equations. The use or non-use of higher order
    !!        electrostatic terms is determined by the E-theta flag
    !!        on the supporting data file.
    !!     2  HC + DHC
    !! The computed activity coefficients may be normalized to make
    !! them consistent with a given pH convention. This is determined
    !! by the option flag iopg(2):
    !!    -1  No normalization
    !!     0  Normalization to the NBS pH scale
    !!     1  Normalization to the Mesmer pH scale (pH is numerically
    !!        equal to -log m(H+))
    !! This subroutine is called by:
    !!   EQLIB/ngcadv.f
    !!   EQ3NR/arrset.f
    !!   EQ6/exivar.f
    !! Principal input:
    !!   adh    = Debye-Huckel A(gamma) parameter
    !!   azero  = array of hard core diameters of aqueous species
    !!   bdh    = Debye-Huckel B parameter
    !!   bdot   = B-dot parameter
    !!   cco2   = coefficients of the Drummond (1981) equation
    !!   conc   = array of species concentrations
    !!   fje    = the 3rd-order electrostatic moment function J
    !!   fxi    = the ionic strength (the 2nd-order electrostatic
    !!              moment function I)
    !!   ielam  = flag to not use (-1) or use (0) "higher order"
    !!              (higher than 2nd order) electrostatic terms in
    !!              those activity coefficient models which contain
    !!              provision for such terms
    !!   insgf  = array of activity coefficient flags for aqueous
    !!              species, used in the case of neutral solute
    !!              species when using the B-dot model
    !!   iopg   = array of activity coefficient option switches
    !!   izmax  = max norm of the electrical charges of the aqueous
    !!              species
    !!   jcsort = array of species indices, in order of increasing
    !!              concentration, but with sorting restricted to within
    !!              phase ranges
    !!   narn1  = index of the first aqueous species; this is also
    !!              the index of the solvent, water
    !!   narn2  = index of the last aqueous species
    !!   nchlor = index of the aqueous chloride ion
    !!   nhydr  = index of the aqueous hydrogen ion
    !!   omega  = water constant; ~55.51
    !!   press  = pressure, bars
    !!   sigmam = sum of solute molalities
    !!   tempk  = temperature, K
    !!   xbarwc = mole fraction of water (calculated)
    !!   xbrwlc = log mole fraction of water (calculated)
    !!   zchar  = array of electrical charge numbers for the various
    !!              species
    !!   zchsq2 = array of z**2/2 values
    !!   zchcu6 = array of z**3/6 values
    !! Principal Output:
    !!   acflgc = array of log activity coefficients (calculated)
    !!              of the various species
    !!   actwlc = log activity of water (calculated)
    implicit none

    ! Calling sequence variable declarations, for all variables except
    ! those intimately related to Pitzer's equations.
    integer :: ipbtmx
    integer :: natmax
    integer :: nopgmx
    integer :: nstmax

    integer :: insgf(natmax)
    integer :: iopg(nopgmx)
    integer :: jcsort(nstmax)

    integer :: ielam
    integer :: izmax
    integer :: narn1
    integer :: narn2
    integer :: nchlor
    integer :: nhydr
    integer :: noutpt
    integer :: nttyo

    integer :: ifcphi1
    integer :: ifcphi2
    integer :: ifnnn
    integer :: ifn2n
    integer :: ifpsi1
    integer :: ifpsi2
    integer :: ifzeta
    integer :: ilcphi1
    integer :: ilcphi2
    integer :: ilnnn
    integer :: iln2n
    integer :: ilpsi1
    integer :: ilpsi2
    integer :: ilzeta

    real(kind=8) :: acflgc(nstmax)
    real(kind=8) :: azero(natmax)
    real(kind=8) :: a3bars(natmax)
    real(kind=8) :: conc(nstmax)
    real(kind=8) :: cco2(5)
    real(kind=8) :: zchar(nstmax)
    real(kind=8) :: zchsq2(nstmax)
    real(kind=8) :: zchcu6(nstmax)

    character(len=48) :: uspec(nstmax)

    real(kind=8) :: adh
    real(kind=8) :: adhh
    real(kind=8) :: adhv
    real(kind=8) :: aphi
    real(kind=8) :: bdh
    real(kind=8) :: bdhh
    real(kind=8) :: bdhv
    real(kind=8) :: bdot
    real(kind=8) :: bdoth
    real(kind=8) :: bdotv

    real(kind=8) :: abar
    real(kind=8) :: actwlc
    real(kind=8) :: al10
    real(kind=8) :: a3bar
    real(kind=8) :: fje
    real(kind=8) :: fxi
    real(kind=8) :: omega
    real(kind=8) :: press
    real(kind=8) :: sigmam
    real(kind=8) :: tempk
    real(kind=8) :: xbarwc
    real(kind=8) :: xbrwlc

    ! Calling sequence variable declarations, for all variables
    ! intimately related to Pitzer's equations.
    ! Integer variables and arrays.
    integer :: napmax
    integer :: nazmmx
    integer :: nazpmx
    integer :: nmutmx
    integer :: nmxmax
    integer :: nsltmx
    integer :: nsxmax

    integer :: napt
    integer :: nmut
    integer :: nslt

    logical :: qhawep
    logical :: qpit75

    integer :: nalpha(nsltmx)
    integer :: nmux(3,nmutmx)
    integer :: nmxi(2,natmax)
    integer :: nmxx(3,nmxmax)
    integer :: nslx(2,nsltmx)
    integer :: nsxi(2,natmax)
    integer :: nsxx(2,nsxmax)

    ! Empirical interaction coefficients and related functions at
    ! temperature and pressure.
    real(kind=8) :: dgpit(2,ipbtmx,napmax)
    real(kind=8) :: dpslm(2,nsltmx)
    real(kind=8) :: gpit(ipbtmx,napmax)
    real(kind=8) :: palpha(ipbtmx,napmax)
    real(kind=8) :: pmu(nmutmx)
    real(kind=8) :: pslamn(0:ipbtmx,nsltmx)
    real(kind=8) :: pslm(nsltmx)

    ! Theoretical higher-order electrostatic interaction coefficients
    ! and related functions at temperature and pressure.
    real(kind=8) :: delam(2,nazpmx,nazpmx)
    real(kind=8) :: dpelm(2,nazpmx,nazpmx)
    real(kind=8) :: dselm(2,nazmmx:nazpmx)
    real(kind=8) :: elam(nazpmx,nazpmx)
    real(kind=8) :: pelm(nazpmx,nazpmx)
    real(kind=8) :: selm(nazmmx:nazpmx)

    ! Local variable declarations.
    integer :: iz
    integer :: j2
    integer :: n
    integer :: na
    integer :: nref
    integer :: ns
    integer :: nss
    integer :: ns1
    integer :: ns2
    integer :: ns3

    integer :: ilnobl

    character(len=8) :: ux8

    real(kind=8) :: acfnbs
    real(kind=8) :: afac
    real(kind=8) :: alm
    real(kind=8) :: bt
    real(kind=8) :: btp
    real(kind=8) :: chc
    real(kind=8) :: delacf
    real(kind=8) :: elsum
    real(kind=8) :: elsump
    real(kind=8) :: elsums
    real(kind=8) :: elsumw
    real(kind=8) :: f
    real(kind=8) :: fa
    real(kind=8) :: flsum
    real(kind=8) :: fp
    real(kind=8) :: fpp
    real(kind=8) :: lnaw
    real(kind=8) :: lngam
    real(kind=8) :: muterm
    real(kind=8) :: musum
    real(kind=8) :: musumw
    real(kind=8) :: slsum
    real(kind=8) :: slsump
    real(kind=8) :: spsum
    real(kind=8) :: spsump
    real(kind=8) :: ssumw
    real(kind=8) :: xx
    real(kind=8) :: xy
    real(kind=8) :: zchtot

    real(kind=8) :: csum
    real(kind=8) :: csumca
    real(kind=8) :: csum1
    real(kind=8) :: csum2
    real(kind=8) :: psum1
    real(kind=8) :: psum2
    real(kind=8) :: zesum

    real(kind=8) :: tlg

    ! Here btp is the effective product of Debye-Huckel B and the hard
    ! core diameter used in Pitzer's equations.
    data btp /1.2/

    ! Here chc is the hard core repulsion constant.
    data chc /1.2577e-3/

    ! Note: the following statements don't actually do anything except
    ! cause the compiler not to complain that fje, press, and zchcu6
    ! are not used.
    xx = fje
    xy = xx
    fje = xy

    xx = press
    xy = xx
    press = xy

    xx = zchcu6(1)
    xy = xx
    zchcu6(1) = xy

    ! Calculate the mole fraction of water from sigmam. Note that
    ! this value is subject to under-relaxation truncation constraints
    ! placed on sigmam.
    xbarwc = omega/(omega + sigmam)
    xbrwlc = tlg(xbarwc)

    ! Calculate the total charge per kg H2O (Z in Pitzer's equations).
    ! Note that this is a sorted summation. It is not necessary to
    ! exclude water, because it has zero charge.
    zchtot = 0.

    do nss = narn1,narn2
        ns = jcsort(nss)
        zchtot = zchtot + conc(ns)*abs(zchar(ns))
    end do

    if (iopg(1) .eq. -1) then
        ! Use the Davies equation.
        call gdavie(acflgc,actwlc,adh,al10,fxi,narn1,narn2,nstmax,omega,sigmam,xbrwlc,zchsq2)
    end if

    if (iopg(1) .eq. 0) then
        ! Use the B-dot equation and associated approximations.
        call gbdot(acflgc,actwlc,adh,al10,azero,bdh,bdot,cco2,fxi,insgf,narn1,narn2,natmax,nstmax,omega,sigmam,tempk,xbrwlc,zchar,zchsq2)
    end if

    if (iopg(1) .eq. 1) then
        ! Use Pitzer's equations.
        ! Compute the Debye-Huckel function f and its ionic strength
        ! derivatives.
        bt = btp
        call gfdho(aphi,bt,f,fp,fpp,fxi)

        elsumw = 0.
        elsums = 0.
        elsump = 0.

        if (ielam .ge. 0) then
            ! Get E-lambda (elam) and its derivatives (delam) for the
            ! various charge combinations.
            call gelam(aphi,delam,dpelm,elam,fxi,izmax,nazpmx,pelm,qpit75)

            ! Compute the following second order sums in E-lambda and its
            ! ionic strength derivatives:
            !   elsumw: SUM(ij) [E-lambda(ij) + I*E-lambda'(ij)]*m(i)*m(j)
            !   elsums: SUM(ij) E-lambda'(ij)*m(i)*m(j)
            !   elsump: SUM(ij) E-lambda''(ij)*m(i)*m(j) [Not used]
            ! Here ij refers to any species pair, but non-zero
            ! contributions only arise from cc' and aa' pairs. Note that:
            !   elsumw = 2*
            !         [SUM(c'>c) {E-theta(cc') + I*E-theta'(cc')}m(c)m(c')
            !        + SUM(a'>a) {E-theta(aa') + I*E-theta'(aa')}m(a)m(a')]
            !   elsums = 2*[SUM(c'>c) E-theta'(cc')m(c)m(c')
            !             + SUM(a'>a) E-theta'(aa')m(a)m(a')]
            !   elsump = 2*[SUM(c'>c) E-theta''(cc')m(c)m(c')
            !             + SUM(a'>a) E-theta''(aa')m(a)m(a')] [Not used]
            ! Note also that "flim" in EQ3/6 is 2 * F, as F is commonly
            ! defined in the Pitzer literature. This is because
            ! flim will be multipled by z(i)^2/2 instead of z(i)^2.
            call gesum(conc,delam,elam,elsump,elsums,elsumw,fxi,narn1,narn2,nazpmx,nstmax,zchar)

            ! Compute the following first order sum arrays in E-lambda
            ! and its ionic strength derivative:
            !   selm(i):  SUM(j) E-lambda(ij)*m(j)
            !   dselm(1,i): SUM(j) E-lambda'(ij)*m(j) [Not used]
            call gselm(conc,delam,dselm,elam,izmax,narn1,narn2,nazmmx,nazpmx,nstmax,selm,zchar)
        end if

        ! Compute the g(x) function (gpit) and its derivatives (dgpit)
        ! for pairs of alpha coefficients.
        call gdd(dgpit,fxi,gpit,ipbtmx,napmax,napt,palpha)

        ! Compute the S-lambda functions and their ionic strength
        ! derivatives.
        call gslam(dgpit,dpslm,gpit,ipbtmx,nalpha,napmax,nslt,nsltmx,pslamn,pslm)

        ! Compute the following second order sum in S-lambda and its
        ! ionic strength derivatives:
        !   ssumw: SUM(ij) [ S-lambda(ij) + I*S-lambda'(ij) ]*m(i)*m(j)
        call gssum(conc,dpslm,fxi,nslt,nsltmx,nslx,nstmax,pslm,ssumw,uspec)

        if (qhawep) then
            ! Use the form corresponding to C (from Cphi), psi, zeta,
            ! and any mu not derived from Cphi, psi, or zeta.
            ! Calculate the following sum that appears in the activity
            ! coefficients of both water and solutes:
            !   csumca:  SUM(ca) m(c)m(a)C(ca)
            ! This sum is technically a double sum, but it is evaluated
            ! here as an equivalent single sum over all ca combinations.
            csumca = 0.

            do n = ifcphi1,ilcphi1
                ns1 = nmux(1,n)
                ns3 = nmux(3,n)
                csumca = csumca + pmu(n)*conc(ns1)*conc(ns3)
            end do
        end if

        ! Compute the contributions to the activity coefficient of
        ! water from third-order terms.
        muterm = 0.

        if (qhawep) then
            ! Use the form corresponding to C (from Cphi), psi, zeta,
            ! and any mu not derived from Cphi, psi, or zeta.
            csum = zchtot*csumca

            psum1 = 0.

            ! This sum is technically a triple sum, but is evaluated
            ! as a single sum over all cc'a combinations.
            do n = ifpsi1,ilpsi1
                ns1 = nmux(1,n)
                ns2 = nmux(2,n)
                ns3 = nmux(3,n)
                psum1 = psum1 + pmu(n)*conc(ns1)*conc(ns2)*conc(ns3)
            end do

            psum2 = 0.

            ! This sum is technically a triple sum, but is evaluated
            ! as a single sum over all aa'c combinations.
            do n = ifpsi2,ilpsi2
                ns1 = nmux(1,n)
                ns2 = nmux(2,n)
                ns3 = nmux(3,n)
                psum2 = psum2 + pmu(n)*conc(ns1)*conc(ns2)*conc(ns3)
            end do

            zesum = 0.

            ! This sum is technically a triple sum, but is evaluated
            ! as a single sum over all nca combinations.
            do n = ifzeta,ilzeta
                ns1 = nmux(1,n)
                ns2 = nmux(2,n)
                ns3 = nmux(3,n)
                zesum = zesum + pmu(n)*conc(ns1)*conc(ns2)*conc(ns3)
            end do

            musum = 0.

            ! The following involves sums in mu(nnn) and mu(nnn').
            do n = ifnnn,ilnnn
                musum = musum + 2.*pmu(n)*conc(ns1)*conc(ns1)*conc(ns1)
            end do

            do n = ifn2n,iln2n
                ns1 = nmux(1,n)
                ns2 = nmux(2,n)
                ns3 = nmux(3,n)

                ! Note that it doesn't matter here whether ns1 = ns2 or
                ! ns2 = ns3.
                musum = musum + 6.*pmu(n)*conc(ns1)*conc(ns2)*conc(ns3)
            end do

            muterm = csum + psum1 + psum2 + zesum + musum

            ! Multiply muterm by 2 for use in calculating the activity
            ! of water directly.
            muterm = 2.*muterm
        else
            ! Use the original pure mu form.
            ! Compute the following third order sum in mu:
            !   musumw: SUM(ijk) mu(ijk)*m(i)*m(j)*m(k)
            call gmsum(conc,musumw,nmut,nmutmx,nmux,nstmax,pmu,uspec)
            muterm = 2.*musumw
        end if

        ! Compute the activity of water.
        actwlc = -(1./(al10*omega)) *  (sigmam + fxi*fp - f + elsumw + ssumw + muterm)

        ! Compute the activity coefficient of water. Note that this
        ! is mole-fraction based.
        acflgc(narn1) = actwlc - xbrwlc

        ! Compute the following second order sums in the first and
        ! second ionic strength derivatives of S-lambda:
        !   spsum:  SUM(jk) S-lambda'(jk)*m(j)*m(k)
        !   spsump: SUM(jk) S-lambda''(jk)*m(j)*m(k) [Not used]
        ! Here jk always refers to a cation-ion pair. Note that:
        !   spsum = 2 * SUM(ca) B'(ca)m(a)m(c)
        !   spsump = 2 * SUM(ca) B''(ca)m(a)m(c)
        ! Note also that "flim" in EQ3/6 is 2 * F, as F is commonly
        ! defined in the Pitzer literature. This is because
        ! flim will be multipled by z(i)^2/2 instead of z(i)^2.
        call gsdsm(conc,dpslm,nslt,nsltmx,nslx,nstmax,spsum,spsump,uspec)

        flsum = fp + elsums + spsum

        ! Compute the activity coefficients of the aqueous
        ! solute species.
        do ns = narn1 + 1,narn2
            na = ns - narn1 + 1

            ! Logically, na (or ns) denotes aqueous species i.
            acflgc(ns) = 0.
            iz = nint(zchar(ns))

            elsum = 0.

            if (ielam .ge. 0) then
                ! Get the E-lambda sum (this only depends on the
                ! magnitude of the charge number):
                !   selm(i): SUM(j) E-lambda(ij)*m(j)
                elsum = selm(iz)
            end if

            ! Compute the following first order sums in S-lambda
            ! and its ionic strength derivative:
            !   slsum(i):  SUM(j) S-lambda(ij)*m(j)
            !   slsump(i): SUM(j) S-lambda'(ij)*m(j)
            call gsgsm(conc,dpslm,na,natmax,nsltmx,nstmax,nsxi,nsxmax,nsxx,pslm,slsum,slsump,uspec)

            ! Compute the contributions from third-order terms.
            muterm = 0.

            if (qhawep) then
                ! Use the form corresponding to C (from Cphi), psi, zeta,
                ! and any mu not derived from Cphi, psi, or zeta.
                csum1 = 0.

                if (zchar(ns) .gt. 0.) then
                    ! Have a cation. Loop over C(ca).
                    do n = ifcphi1,ilcphi1
                        ns1 = nmux(1,n)

                        if (ns .eq. ns1) then
                            ns3 = nmux(3,n)
                            csum1 = csum1 + pmu(n)*conc(ns3)
                        end if
                    end do
                else if (zchar(ns) .lt. 0.) then
                    ! Have an anion. Loop over C(ac).
                    do n = ifcphi2,ilcphi2
                        ns1 = nmux(1,n)

                        if (ns .eq. ns1) then
                            ns3 = nmux(3,n)
                            csum1 = csum1 + pmu(n)*conc(ns3)
                        end if
                    end do
                end if

                csum1 = zchtot*csum1

                psum1 = 0.

                if (zchar(ns) .gt. 0.) then
                    ! Have a cation. Loop over psi(cc'a).
                    do n = ifpsi1,ilpsi1
                        ns1 = nmux(1,n)
                        ns2 = nmux(2,n)
                        ns3 = nmux(3,n)

                        if (ns .eq. ns1) then
                            psum1 = psum1 + pmu(n)*conc(ns2)*conc(ns3)
                        else if (ns .eq. ns2) then
                            psum1 = psum1 + pmu(n)*conc(ns1)*conc(ns3)
                        end if
                    end do
                else if (zchar(ns) .lt. 0.) then
                    ! Have an anion. Loop over psi(aa'c).
                    do n = ifpsi2,ilpsi2
                        ns1 = nmux(1,n)
                        ns2 = nmux(2,n)
                        ns3 = nmux(3,n)

                        if (ns .eq. ns1) then
                            psum1 = psum1 + pmu(n)*conc(ns2)*conc(ns3)
                        else if (ns .eq. ns2) then
                            psum1 = psum1 + pmu(n)*conc(ns1)*conc(ns3)
                        end if
                    end do
                end if

                psum2 = 0.

                ! This sum is technically a double sum, but is evaluated
                ! as a single sum over all distinct values of psi(aa'c)
                ! or psi(cc'a).
                if (zchar(ns) .gt. 0.) then
                    ! Have a cation. Loop over psi(aa'c).
                    do n = ifpsi2,ilpsi2
                        ns1 = nmux(1,n)
                        ns2 = nmux(2,n)
                        ns3 = nmux(3,n)

                        if (ns .eq. ns3) then
                            psum2 = psum2 + pmu(n)*conc(ns1)*conc(ns2)
                        end if
                    end do
                else if (zchar(ns) .lt. 0.) then
                    ! Have an anion. Loop over psi(cc'a).
                    do n = ifpsi1,ilpsi1
                        ns1 = nmux(1,n)
                        ns2 = nmux(2,n)
                        ns3 = nmux(3,n)

                        if (ns .eq. ns3) then
                            psum2 = psum2 + pmu(n)*conc(ns1)*conc(ns2)
                        end if
                    end do
                end if

                csum2 = abs(zchar(ns))*csumca

                zesum = 0.

                ! This sum is technically a double sum, but is evaluated
                ! as a single sum over all distinct values of zeta(nca).
                ! Have a cation, an anion, or a neutral. Loop over zeta(nca).
                do n = ifzeta,ilzeta
                    ns1 = nmux(1,n)
                    ns2 = nmux(2,n)
                    ns3 = nmux(3,n)

                    if (ns .eq. ns1) then
                        zesum = zesum + pmu(n)*conc(ns2)*conc(ns3)
                    else if (ns .eq. ns2) then
                        zesum = zesum + pmu(n)*conc(ns1)*conc(ns3)
                    else if (ns .eq. ns3) then
                        zesum = zesum + pmu(n)*conc(ns1)*conc(ns2)
                    end if
                end do

                musum = 0.

                if (iz .eq. 0) then
                    ! The following involves sums in mu(nnn) and mu(nnn').
                    do n = ifnnn,ilnnn
                        ns1 = nmux(1,n)

                        if (ns .eq. ns1) then
                            musum = musum + 3.*pmu(n)*conc(ns1)*conc(ns1)
                        end if
                    end do

                    do n = ifn2n,iln2n
                        ns1 = nmux(1,n)
                        ns2 = nmux(2,n)
                        ns3 = nmux(3,n)

                        if (ns .eq. ns1) then
                            if (ns1 .eq. ns2) then
                                musum = musum + 6.*pmu(n)*conc(ns1)*conc(ns3)
                            else
                                musum = musum + 3.*pmu(n)*conc(ns3)*conc(ns3)
                            end if
                        else if (ns .eq. ns3) then
                            if (ns2 .eq. ns3) then
                                musum = musum + 6.*pmu(n)*conc(ns1)*conc(ns3)
                            else
                                musum = musum + 3.*pmu(n)*conc(ns1)*conc(ns1)
                            end if
                        end if
                    end do
                end if

                muterm = csum1 + psum1 + psum2 + csum2 + zesum + musum
            else
                ! Use the original pure mu form.
                ! Compute the following second order sum in mu:
                !   musum(i): SUM(jk) mu(ijk)*m(j)*m(k)
                call gmdsm(conc,musum,na,natmax,nmutmx,nmxi,nmxmax,nmxx,ns,nstmax,pmu,uspec)

                muterm = 3.*musum
            end if

            ! Note that zchsq2(ns) is z(i)^2/2.
            acflgc(ns) = (zchsq2(ns)*flsum + 2.0*elsum + 2.0*slsum    + muterm) / al10
        end do
    end if

    if (iopg(1) .eq. 2) then
        ! Use the HC + DHC equations
        ! Get the average ion size (abar), the characteristic average
        ! cube of the distance of closest approach (a3bars), and the
        ! average cube of the distance of closest approach (a3bar).
        call cabar(abar,azero,conc,jcsort,fxi,narn1,narn2,natmax,nstmax,zchsq2)

        call ca3bar(azero,a3bar,a3bars,conc,jcsort,narn1,narn2,natmax,nstmax,sigmam)

        ! Get the Debye-Huckel function f and its ionic strength
        ! derivatives.
        bt = bdh*abar
        call gfdhc(adh,bt,f,fp,fpp,fxi)

        ! Compute the activity of water.
        xx = 1. + (sigmam/omega)
        alm = log(xx)
        lnaw = - alm - (fxi*fp - f)/omega  - (chc*a3bar*(sigmam**2))/(omega)
        actwlc = lnaw/al10
        acflgc(narn1) = actwlc - xbrwlc

        ! Compute the activity coefficients of the aqueous
        ! solute species.
        fa = 3.*fxi*fp - 2.*f

        do ns = narn1 + 1,narn2
            acflgc(ns) = 0.
            na = ns - narn1 + 1
            afac = (azero(na) - abar)/abar
            lngam = -alm + zchsq2(ns)*fp + (zchsq2(ns)*fa*afac)/fxi    + 2.*chc*sigmam**a3bars(na)
            acflgc(ns) = lngam/al10
        end do
    end if

    if (iopg(1).le.-2 .or. iopg(1).ge.3) then
        write (ux8,'(i5)') iopg(1)
        call lejust(ux8)
        j2 = ilnobl(ux8)
        write (noutpt,1000) ux8(1:j2)
        write (nttyo,1000) ux8(1:j2)
1000 format(/' * Error - (EQLIBG/gcoeff) Programming error trap:',' Have an',/7x,'illegal value of ',a,' for the option switch',' iopg(1), which determines',/7x,'the model to use for the',' activity coefficients of aqueous species.')

        stop
    end if

    ! Make activity coefficients consistent with a specified pH
    ! scale.
    !   iopg(2) = -1   Use activity coefficients as is
    !            = 0   Make consistent with NBS pH scale
    !            = 1   Make consistent with the Mesmer pH scale
    !                  (log gamma H+ = 0)
    if (iopg(2) .eq. 0)  then
        call nbsgam(acfnbs,adh,fxi,nchlor,noutpt,nttyo)
        delacf = acfnbs - acflgc(nchlor)
        nref = nchlor
        call gcscal(acflgc,delacf,narn1,narn2,nref,nstmax,zchar)
    end if

    if (iopg(2) .eq. 1) then
        delacf = -acflgc(nhydr)
        nref = nhydr
        call gcscal(acflgc,delacf,narn1,narn2,nref,nstmax,zchar)
    end if
end subroutine gcoeff
