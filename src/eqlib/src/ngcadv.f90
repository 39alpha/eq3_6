subroutine ngcadv(abar,acflg,acflgo,actwlc,adh,adhh,adhv,afcnst,al10,aphi,azero,a3bar,a3bars,bacfmx,bdh,bdhh,bdhv,bdot,bdoth,bdotv,bgamx,bpx,bsigmm,bfje,bfxi,cco2,cgexj,chfacf,chfsgm,conc,delam,dgpit,dpelm,dpslm,dselm,elam,eps100,fje,fjeo,fxi,fxio,gpit,ibpxt,ielam,ifcphi1,ifcphi2,ifnnn,ifn2n,ifpsi1,ifpsi2,ifzeta,ilcphi1,ilcphi2,ilnnn,iln2n,ilpsi1,ilpsi2,ilzeta,insgf,iopg,iter,ipndx1,ixrn1,ixrn2,izmax,jcsort,jern1,jern2,jgext,jsol,kx1,kxt,nalpha,napt,narn1,narn2,nchlor,ncmpr,net,nhydr,nmut,nmux,nmxi,nmxx,noutpt,nslt,nslx,nst,nsxi,nsxx,nttyo,omega,palpha,pelm,pmu,press,pslamn,pslm,qhawep,qpit75,qpracf,q6mode,rlxgam,selm,sigmam,sigmmo,tempk,ubacmx,ubgamx,uphase,uspec,wfac,xbar,xbarlg,xbarwc,xbrwlc,zchar,zchcu6,zchsq2)
    !! This subroutine recalculates the ionic strength, etc., and the
    !! activity of water and the molal activity coefficients of aqueous
    !! species. If necessary, it also recalculates the mole fraction
    !! activity coefficients of solid solution components. The
    !! values for the above parameters returned by this subroutine are
    !! subject to under-relaxation controls intended to aid convergence,
    !! particularly in the case of concentrated aqueous solutions.
    !! Associated residual functions such as bsigmm, bfxi, and bgamx
    !! are not affected by this under-relaxation.
    !! This subroutine no longer recalculates the concentrations and
    !! other properties of dependent species.
    !! This subroutine is called by:
    !!   EQLIB/newton.f
    !!   EQ3NR/arrset.f
    !!   EQ6/optmzr.f
    !! Principal input:
    !!   adh    = Debye-Huckel A(gamma) parameter
    !!   acflgo = array of old values of the activity coefficients of
    !!              the various species
    !!   aphi   = Debye-Huckel A(phi) parameter
    !!   azero  = array of hard core diameters of aqueous species
    !!   bdh    = Debye-Huckel B(gamma) parameter
    !!   bdot   = B-dot parameter
    !!   cco2   = coefficients of the Drummond (1981) equation
    !!   conc   = array of species concentrations
    !!   fjeo   = old value of the ionic asymmetry (J)
    !!   fxio   = old value of the ionic strength (I)
    !!   ibpxt  = array of the numbers of non-zero site-mixing
    !!              paramters for computing activity coefficients in
    !!              solid solutions
    !!   ielam  = flag to not use (-1) or use (0) "higher order"
    !!              (higher than 2nd order) electrostatic terms in
    !!              those activity coefficient models which contain
    !!              provision for such terms
    !!   insgf  = array of activity coefficient flags for aqueous
    !!              species, used in the case of neutral solute
    !!              species when using the B-dot model
    !!   iopg   = array of activity coefficient option switches
    !!   iter = Newton-Raphson iteration number
    !!   izmax  = max norm of the electrical charges of the aqueous
    !!              species
    !!   jcsort = array of species indices, in order of increasing
    !!              concentration, but with sorting restricted to within
    !!              phase ranges
    !!   jsol   = array of identifiers for solid solution activity
    !!              coefficient models
    !!   narn1  = index of the first aqueous species; this is also
    !!              the index of the solvent, water
    !!   narn2  = index of the last aqueous species
    !!   nchlor = index of the aqueous chloride ion
    !!   ncmpr  = array giving the range of species belonging to a given
    !!              phase: ncmpr(1,np) is the first such species for the
    !!              np-th phase, and ncmpr(2,np) is the last such species
    !!   nhydr  = index of the aqueous hydrogen ion
    !!   omega  = water constant; ~55.51
    !!   press  = pressure, bars
    !!   qpracf = debugging print flag, causes print of activity
    !!              coefficient update
    !!   q6mode = flag denoting usage for EQ3NR or EQ6:
    !!              .false. = EQ3NR
    !!              .true.  = EQ6NR
    !!   rlxgam = reported under-relaxation factor applied to corrections
    !!              in the activity coefficients of aqueous species
    !!   sigmmo = old value of sigmam
    !!   tempk  = temperature, K
    !!   uphase = array of phase names
    !!   uspec  = array of species names
    !!   wfac   = array of non-ideality parameters for solid solutions
    !!   xbarwc = mole fraction of water (calculated)
    !!   xbrwlc = log mole fraction of water (calculated)
    !!   zchar  = array of electrical charge numbers for the various
    !!              species
    !!   zchsq2 = array of (z**2)/2 values
    !!   zchcu6 = array of (z**3/)6 values
    !! Principal output:
    !!   abar     average ion size
    !!   acflg  = array of log activity coefficients of the various
    !!              species
    !!   actwlc = log activity of water (calculated)
    !!   a3bar    average cube of distance of closest apporach
    !!   a3bars   characteristic average cube of distance of closest
    !!              apporach for each solute species
    !!   fje    = the ionic asymmetry (the 3rd-order electrostatic
    !!              moment function J)
    !!   fxi    = the ionic strength (the 2nd-order electrostatic
    !!              moment function I)
    !!   sigmam = the sum of solute molalities
    implicit none

    include 'eqlib/eqldv.h'

    ! Calling sequence variable declarations.
    integer :: noutpt
    integer :: nttyo

    integer :: ibpxt(nxtmax)
    integer :: insgf(natmax)
    integer :: iopg(nopgmx)
    integer :: ipndx1(kmax)
    integer :: jcsort(nstmax)
    integer :: jern1(jetmax,netmax)
    integer :: jern2(jetmax,netmax)
    integer :: jgext(netmax)
    integer :: jsol(nxtmax)
    integer :: ncmpr(2,nptmax)

    integer :: nalpha(nsltmx)
    integer :: nmux(3,nmutmx)
    integer :: nmxi(2,natmax)
    integer :: nmxx(3,nmxmax)
    integer :: nslx(2,nsltmx)
    integer :: nsxi(2,natmax)
    integer :: nsxx(2,nsxmax)

    integer :: napt
    integer :: nmut
    integer :: nslt

    integer :: ielam
    integer :: iter
    integer :: ixrn1
    integer :: ixrn2
    integer :: izmax
    integer :: kx1
    integer :: kxt
    integer :: narn1
    integer :: narn2
    integer :: nchlor
    integer :: net
    integer :: nhydr
    integer :: nst

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

    logical :: qhawep
    logical :: qpit75
    logical :: qpracf
    logical :: q6mode

    character(len=48) :: uspec(nstmax)
    character(len=48) :: ubacmx
    character(len=48) :: ubgamx
    character(len=24) :: uphase(nptmax)

    real(kind=8) :: dgpit(2,ipbtmx,napmax)
    real(kind=8) :: dpslm(2,nsltmx)
    real(kind=8) :: gpit(ipbtmx,napmax)
    real(kind=8) :: palpha(ipbtmx,napmax)
    real(kind=8) :: pmu(nmutmx)
    real(kind=8) :: pslamn(0:ipbtmx,nsltmx)
    real(kind=8) :: pslm(nsltmx)

    real(kind=8) :: delam(2,nazpmx,nazpmx)
    real(kind=8) :: dpelm(2,nazpmx,nazpmx)
    real(kind=8) :: dselm(2,nazmmx:nazpmx)
    real(kind=8) :: elam(nazpmx,nazpmx)
    real(kind=8) :: pelm(nazpmx,nazpmx)
    real(kind=8) :: selm(nazmmx:nazpmx)

    real(kind=8) :: acflg(nstmax)
    real(kind=8) :: acflgo(nstmax)
    real(kind=8) :: azero(natmax)
    real(kind=8) :: a3bars(natmax)
    real(kind=8) :: bpx(ibpxmx,nxtmax)
    real(kind=8) :: cco2(5)
    real(kind=8) :: cgexj(jetmax,netmax)
    real(kind=8) :: conc(nstmax)
    real(kind=8) :: wfac(iktmax,nxtmax)
    real(kind=8) :: xbar(nstmax)
    real(kind=8) :: xbarlg(nstmax)
    real(kind=8) :: zchar(nstmax)
    real(kind=8) :: zchsq2(nstmax)
    real(kind=8) :: zchcu6(nstmax)

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
    real(kind=8) :: afcnst
    real(kind=8) :: al10
    real(kind=8) :: a3bar
    real(kind=8) :: bacfmx
    real(kind=8) :: bgamx
    real(kind=8) :: bsigmm
    real(kind=8) :: bfje
    real(kind=8) :: bfxi
    real(kind=8) :: chfacf
    real(kind=8) :: chfsgm
    real(kind=8) :: eps100
    real(kind=8) :: fje
    real(kind=8) :: fjeo
    real(kind=8) :: fxi
    real(kind=8) :: fxio
    real(kind=8) :: omega
    real(kind=8) :: press
    real(kind=8) :: rlxgam
    real(kind=8) :: sigmam
    real(kind=8) :: sigmmo
    real(kind=8) :: tempk
    real(kind=8) :: xbarwc
    real(kind=8) :: xbrwlc

    ! Local variable declarations.
    integer :: itgfix
    integer :: j2
    integer :: kcol
    integer :: np
    integer :: ns

    integer :: ilnobl

    logical :: qconc1
    logical :: qconc2
    logical :: qconc3

    real(kind=8) :: adfx
    real(kind=8) :: adfxmx
    real(kind=8) :: afjea
    real(kind=8) :: ax
    real(kind=8) :: bgamxo
    real(kind=8) :: chfsmi
    real(kind=8) :: dfje
    real(kind=8) :: dfx
    real(kind=8) :: fjec
    real(kind=8) :: fxic
    real(kind=8) :: rlxgcf
    real(kind=8) :: rlxgmo
    real(kind=8) :: rlxgml
    real(kind=8) :: sigmmc
    real(kind=8) :: sxl
    real(kind=8) :: sxu

    ! Compute the inverse change limit on "Sigma m" (the sum of solute
    ! molalities), the ionic strength, and the like.
    chfsmi = 1./chfsgm

    ! Recompute sigma m and compute the associated residual.
    call csigm(conc,jcsort,narn1,narn2,nstmax,sigmmc)
    bsigmm = 0.

    if (sigmmo .gt. 0.) then
        bsigmm = (sigmmc - sigmmo)/sigmmo
    end if

    sxu = chfsgm*sigmmo

    if (sigmmo .le. 0.) then
        sxu = sigmmc
    end if

    sigmam = min(sxu,sigmmc)
    sxl = chfsmi*sigmmo
    sigmam = max(sxl,sigmam)

    ! Recompute the ionic strength (the 2nd-order electrostatic
    ! moment function I) and compute the associated residual.
    call cfxi(conc,fxic,jcsort,narn1,narn2,nstmax,zchsq2)
    bfxi = 0.

    if (fxio .gt. 0.) then
        bfxi = (fxic - fxio)/fxio
    end if

    sxu = chfsgm*fxio

    if (fxio .le. 0.) then
        sxu = fxic
    end if

    fxi = min(sxu,fxic)
    sxl = chfsmi*fxio
    fxi = max(sxl,fxi)

    ! Recompute the 3rd-order electrostatic moment function J and
    ! the associated residual. Unlike sigma m and I, J can be zero
    ! or negative, so the treatment is slightly different.
    call cfje(conc,fjec,jcsort,narn1,narn2,nstmax,zchcu6)
    afjea = 0.5*(abs(fjeo) + abs(fjec))
    dfje = fjec - fjeo
    bfje = 0.

    if (afjea .gt. 0.) then
        bfje = dfje/afjea
    end if

    if (fjec.gt.0. .and. fjeo.gt.0.) then
        sxu = chfsgm*fjeo

        if (fjeo .le. 0.) then
            sxu = fjec
        end if

        fje = min(sxu,fjec)
        sxl = 0.
        fje = max(sxl,fje)
    else if (fjec.lt.0. .and. fjeo.lt.0.) then
        sxl = chfsgm*fjeo

        if (fjeo .le. 0.) then
            sxl = fjec
        end if

        fje = max(sxl,fjec)
        sxu = 0.
        fje = min(sxu,fje)
    else if (fjec.lt.0. .and. fjeo.gt.0.) then
        fje = 0.
    else if (fjec.gt.0. .and. fjeo.lt.0.) then
        fje = 0.
    else
        fje = fjec
    end if

    ! Compute the activity coefficients of aqueous species.
    ! Calling sequence substitutions:
    !   acflg for acflgc
    call gcoeff(abar,acflg,actwlc,adh,adhh,adhv,al10,aphi,azero,a3bar,a3bars,bdh,bdhh,bdhv,bdot,bdoth,bdotv,cco2,conc,delam,dgpit,dpelm,dpslm,dselm,elam,fje,fxi,gpit,ielam,ifcphi1,ifcphi2,ifnnn,ifn2n,ifpsi1,ifpsi2,ifzeta,ilcphi1,ilcphi2,ilnnn,iln2n,ilpsi1,ilpsi2,ilzeta,insgf,iopg,ipbtmx,izmax,jcsort,nalpha,napmax,napt,narn1,narn2,natmax,nazmmx,nazpmx,nchlor,nhydr,nmut,nmutmx,nmux,nmxi,nmxmax,nmxx,nopgmx,noutpt,nslt,nsltmx,nslx,nstmax,nsxi,nsxmax,nsxx,nttyo,omega,palpha,pelm,pmu,press,pslamn,pslm,qhawep,qpit75,selm,sigmam,tempk,uspec,xbarwc,xbrwlc,zchar,zchsq2,zchcu6)

    ! Calculate the activity coefficients of exchanger species.
    ! Calling sequence substitutions:
    !   acflg for acflgc
    call lamgex(acflg,cgexj,jern1,jern2,jetmax,jgext,net,netmax,nstmax,xbarlg)

    ! Compute activity coefficient residual norms.
    bgamxo = bgamx

    ! XX   The call to EQLIB/betacf.f is below. Reconcile the usage of
    ! XX   EQLIB/betgam.f and EQLIB/betacf.f when the activity coefficients
    ! XX   of species in non-aqueous phases are updated numerically the same
    ! XX   as those of aqueous species.
    call betgam(acflg,acflgo,bgamx,narn1,narn2,nstmax,ubgamx,uspec)

    ubacmx = ubgamx
    bacfmx = bgamx

    ! Algorithm for under-relaxing changes to the activity coefficients
    ! of aqueous species in concentrated solutions. This is designed to
    ! damp out oscillation in activity coefficient corrections. Here
    ! rlxgam is the under-relaxation factor. Note that activity
    ! coefficients may kept fixed for the first few iterations.
    ! Save the previous value of the under-relaxation factor.
    rlxgmo = rlxgam
    rlxgam = 1.0

    ! Compute a limit (rlxgml)on the under-relaxation factor which
    ! is intended to slow its growth from small values, especially
    ! in the case where the activity coefficients have been fixed
    ! (rlxgam = 0.)
    rlxgml = rlxgmo

    if (rlxgml .lt. 0.10) then
        rlxgml = 0.10
    else
        rlxgml = rlxgml + 0.10
    end if

    sxu = 1.0
    rlxgml = min(sxu,rlxgml)

    ! Determine the course of action according to how concentrated
    ! the solution is, based on the sum of solute molalities
    ! parameter (Sigma m).
    qconc1 = sigmam.ge.2.0 .or. sigmmo.ge.2.0
    qconc2 = sigmam.ge.8.0 .or. sigmmo.ge.8.0
    qconc3 = sigmam.ge.12.0 .or. sigmmo.ge.12.0

    ! Here itgfix is the number of iterations to hold the activity
    ! coefficients constant.
    itgfix = 0

    if (qconc1) then
        itgfix = 6
    end if

    if (qconc2) then
        itgfix = 8
    end if

    if (qconc3) then
        itgfix = 10
    end if

    if (iter .le. itgfix) then
        ! Fix the activity coefficients.
        rlxgam = 0.
    else
        ! Apply other under-relaxation controls if called for.
        if (bgamxo .gt. 0.) then
            ! Apply the limit based on the analysis of residual
            ! function behavior.
            ax = -bgamx/bgamxo

            if (ax .ge. -0.5) then
                rlxgam = rlxgmo*(1./(1. + ax))
                sxl = 0.05
                rlxgam = max(sxl,rlxgam)
            end if
        end if

        ! Apply the growth from small values limit.
        rlxgam = min(rlxgml,rlxgam)

        ! If the under-relaxation factor is close to one, make it one.
        if (abs(rlxgam - 1.0) .le. 0.05) then
            rlxgam = 1.0
        end if

        if (qconc2) then
            ! Force some minimum under-relaxation for more concentrated
            ! solutions.
            if (iter .le. 30) then
                sxu = 0.50
                rlxgam = min(sxu,rlxgam)
            end if
        end if

        if (qconc3) then
            ! Force some minimum under-relaxation for highly concentrated
            ! solutions.
            if (iter .le. 60) then
                sxu = 0.25
                rlxgam = min(sxu,rlxgam)
            end if
        end if
    end if

    if (abs(rlxgam) .gt. eps100) then
        ! Apply a limit to the change in the activity coefficients.
        ! This is yet another form of under-relaxation.
        adfx = abs(acflg(narn1) - acflgo(narn1))
        adfxmx = adfx

        do ns = narn1 + 1,narn2
            adfx = abs(acflg(ns) - acflgo(ns))
            adfxmx = max(adfxmx,adfx)
        end do

        if (adfxmx .gt. chfacf) then
            rlxgcf = chfacf/adfxmx
            rlxgam = min(rlxgcf,rlxgam)
        end if
    end if

    ! Apply under-relaxation, if any, to the activity coefficients.
    if (abs(rlxgam) .le. eps100) then
        ! Keep activity coefficients fixed.
        do ns = narn1,narn2
            acflg(ns) = acflgo(ns)
        end do
    else if (abs(rlxgam - 1.0) .gt. eps100) then
        ! Apply under-relaxation to activity coefficients.
        do ns = narn1,narn2
            dfx = acflg(ns) - acflgo(ns)
            acflg(ns) = acflgo(ns) + rlxgam*dfx
        end do
    end if

    if (qpracf) then
        write (noutpt,2300) sigmam,fxi,fje,xbrwlc,xbarwc
2300 format(/3x,'sigmam= ',1pe12.5,/6x,'fxi= ',e12.5,/6x,'fje= ',e12.5,/3x,'xbrwlc= ',0pf10.5,/3x,'xbarwc= ',1pe12.5)

        write (noutpt,1005)
1005 format(//11x,'Activity Coefficients of Aqueous Species',//7x,'Species',23x,'New',11x,'Old',/)

        do ns = narn1,narn2
            write (noutpt,1010) ns,uspec(ns),acflg(ns),acflgo(ns)
1010 format(1x,i3,2x,a24,2x,1pe12.5,2x,e12.5)
        end do
    end if

    if (q6mode) then
        ! Recalculate the activity coefficients of any solid solution
        ! components present in the ES.
        do kcol = kx1,kxt
            np = ipndx1(kcol)

            ! Calling sequence substitutions:
            !   acflg for acflgc
            call lambda(acflg,afcnst,bpx,ibpxmx,ibpxt,iktmax,ixrn1,ixrn2,jsol,ncmpr,noutpt,np,nptmax,nstmax,nttyo,nxtmax,wfac,xbar,xbarlg,uphase,uspec)

            if (qpracf) then
                j2 = ilnobl(uphase(np))
                write (noutpt,1020) uphase(np)(1:j2)
1020 format(/6x,'Activity Coefficients of Species in ',a,//7x,'Species',23x,'New',11x,'Old',/)

                do ns = narn1,narn2
                    write (noutpt,1010) ns,uspec(ns),acflg(ns),acflgo(ns)
                end do
            end if
        end do

        call betacf(acflg,acflgo,bacfmx,nst,nstmax,ubacmx,uspec)
    end if

999 continue
end subroutine ngcadv