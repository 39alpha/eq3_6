subroutine exivar(abar,acflg,acflgo,act,actlg,actwlc,adh,adhh,adhv,al10,aphi,azero,a3bar,a3bars,bdh,bdhh,bdhv,bdot,bdoth,bdotv,cco2,cdrs,cegexs,cgexj,conc,conclg,cpgexs,egexjc,egexjf,egexs,eps100,fje,fjeo,fo2,fo2lg,fsort,fugac,fugalg,fxi,fxio,ielam,iern1,iern2,ietmax,ifcphi1,ifcphi2,ifnnn,ifn2n,ifpsi1,ifpsi2,ifrn1,ifrn2,ifzeta,igas,igstak,iindx1,ilcphi1,ilcphi2,ilnnn,iln2n,ilpsi1,ilpsi2,ilrn1,ilrn2,ilzeta,imrn1,imrn2,insgf,iodb,iopg,ipbtmx,istack,ixrn1,ixrn2,izmax,jcsort,jern1,jern2,jetmax,jflag,jgext,jgsort,jgstak,jjsort,jpflag,jsflag,jsitex,jssort,jstack,kbt,kdim,kelect,kmax,km1,ko2gaq,kwater,kxt,loph,losp,lsort,mgext,moph,mosp,mrgexs,mtb,napmax,narn1,narn2,natmax,nazmmx,nazpmx,nbasp,nbaspd,nbt,nbtmax,nchlor,ncmpr,ndrs,ndrsmx,ndrsr,nelect,nern1,nern2,net,netmax,ngexsa,ngext,ngrn1,ngrn2,ngt,ngtmax,nhydr,nmutmx,nmxmax,nodbmx,nopgmx,noutpt,no2gaq,nphasx,npt,nptmax,nsltmx,nst,nstmax,nsxmax,nttyo,omega,omeglg,press,qhawep,qpit75,q6mode,sigmam,sigmmo,tempk,ugexj,ugexmo,uphase,uspec,xbar,xbarlg,xbarw,xbarwc,xbrwlc,xbrwlg,xlks,zchar,zchcu6,zchsq2,zgexj,zvclg1,zvec1)
    !! This subroutine expands the system description from the data
    !! read from the input file. This includes estimating the numbers
    !! of moles of all phases and species present, the concentrations,
    !! activity coefficients, and activities of all the species, the
    !! ionic strength, and the sum of the molalities of all aqueous
    !! aqueous solute species.
    !! Only a minimum expansion is carried out by the present
    !! subroutine. Further refinements, if any, are made by EQ6/optmzr.
    !! The iodb(3) debug print option switch is used in this subroutine
    !! to control certain prints that are analogous to those in
    !! EQ6/optmzr.
    !! This subroutine is called by:
    !!   EQ6/path.f
    !! Principal input:
    !!   adh    = Debye-Huckel A(gamma) parameter
    !!   aphi   = Debye-Huckel A(phi) parameter
    !!   azero  = array of hard core diameters of aqueous species
    !!   bdh    = Debye-Huckel B(gamma) parameter
    !!   bdot   = B-dot parameter
    !!   cco2   = coefficients of the Drummond (1981) equation
    !!   cdrs   = array of reaction coefficients
    !!   ielam  = flag to not use (-1) or use (0) higher order
    !!              electrostatic terms in relevant activity coefficient
    !!              models
    !!   insgf  = array of activity coefficient flags for aqueous
    !!              species, used in the case of neutral solute
    !!              species when using the B-dot model
    !!   iopg   = array of activity coefficient option switches
    !!   izmax  = max norm of the electrical charges of the aqueous
    !!              species
    !!   nbt    = the number of basis species
    !!   nchlor = index of the aqueous chloride ion
    !!   ncmpr  = array giving the range in arrays corresponding to
    !!              species of those species which belong to a given
    !!              phase
    !!   ndrs   = array parallel to cdrs giving the index of the
    !!              corresponding species
    !!   ndrsr  = array giving the range in the cdrs/ndrs arrays
    !!              corresonding to the reaction associated with a
    !!              given species
    !!   nhydr  = index of the aqueous hydrogen ion
    !!   press  = pressure, bars
    !!   tempk  = temperature, K
    !!   uphase = array of phase names
    !!   uspec  = array of species names
    !!   xlks   = array of equilibrium constants
    !!   zchar  = array of electrical charge numbers for the various
    !!              species
    !!   zchsq2 = array of (z**2)/2 values
    !!   zchcu6 = array of (z**3/)6 values
    !!   zvec1  = array of master variables
    !!   zvclg1 = array of logarithms of master variables
    !! Principal output:
    !!   abar   = average ion size
    !!   acflg  = array of logarithms of activity coefficients
    !!   acflgo = array of old values of the activity coefficients of
    !!              the various species
    !!   act    = array of species activities
    !!   actlg  = array of logarithms of species activities
    !!   actwlc = log activity of water (calculated)
    !!   a3bar    average cube of distance of closest apporach
    !!   a3bars   characteristic average cube of distance of closest
    !!              apporach for each solute species
    !!   conc   = array of species concentrations
    !!   conclg = array of logarithms of species concentrations
    !!   fo2    = the oxygen fugacity
    !!   fo2lg  = logarithm of the oxygen fugacity
    !!   fugac  = array of gas fugacities
    !!   fugalg = array of logarithms of gas fugacities
    !!   fxi    = the ionic strength
    !!   fxio   = old value of the ionic strength
    !!   moph   = array of numbers of moles of phases
    !!   loph   = array of logarithms of numbers of moles of phases
    !!   mosp   = array of numbers of moles of species
    !!   losp   = array of logarithms of numbers of moles of species
    !!   sigmam = the sum of solute molalities
    !!   sigmmo = old value of sigmam
    !!   xbar   = array of mole fractions of species
    !!   xbarlg = array of logarithms of mole fractions of species
    !!   xbarw  = mole fraction of water
    !!   xbarwc = mole fraction of water (calculated)
    !!   xbrwlc = log mole fraction of water (calculated)
    !!   xbrwlg = log mole fraction of water
    !! Modules.
    !! The module mod6pt contains data required to evaluate Pitzer's
    !! equations. Only a subset of these data is required here.
    use mod6pt

    implicit none

    ! Calling sequence variable declarations.
    integer :: ietmax
    integer :: ipbtmx
    integer :: jetmax
    integer :: kmax
    integer :: napmax
    integer :: natmax
    integer :: nazmmx
    integer :: nazpmx
    integer :: nbtmax
    integer :: ndrsmx
    integer :: netmax
    integer :: ngtmax
    integer :: nmutmx
    integer :: nmxmax
    integer :: nodbmx
    integer :: nopgmx
    integer :: nptmax
    integer :: nsltmx
    integer :: nstmax
    integer :: nsxmax

    integer :: noutpt
    integer :: nttyo

    integer :: igstak(ngtmax)
    integer :: iindx1(kmax)
    integer :: insgf(natmax)
    integer :: iodb(nodbmx)
    integer :: iopg(nopgmx)
    integer :: istack(nstmax)
    integer :: jcsort(nstmax)
    integer :: jern1(jetmax,netmax)
    integer :: jern2(jetmax,netmax)
    integer :: jflag(nstmax)
    integer :: jgext(netmax)
    integer :: jgsort(ngtmax)
    integer :: jgstak(ngtmax)
    integer :: jjsort(nstmax)
    integer :: jpflag(nptmax)
    integer :: jsflag(nstmax)
    integer :: jsitex(nstmax)
    integer :: jssort(nstmax)
    integer :: jstack(nstmax)
    integer :: nbasp(nbtmax)
    integer :: nbaspd(nbtmax)
    integer :: ncmpr(2,nptmax)
    integer :: ndrs(ndrsmx)
    integer :: ndrsr(2,nstmax)
    integer :: ngexsa(ietmax,jetmax,netmax)
    integer :: ngext(jetmax,netmax)
    integer :: nphasx(nstmax)

    integer :: ielam
    integer :: iern1
    integer :: iern2
    integer :: ifrn1
    integer :: ifrn2
    integer :: igas
    integer :: ilrn1
    integer :: ilrn2
    integer :: imrn1
    integer :: imrn2
    integer :: ixrn1
    integer :: ixrn2
    integer :: izmax
    integer :: kbt
    integer :: kdim
    integer :: kelect
    integer :: km1
    integer :: ko2gaq
    integer :: kwater
    integer :: kxt
    integer :: narn1
    integer :: narn2
    integer :: nbt
    integer :: nchlor
    integer :: nelect
    integer :: nern1
    integer :: nern2
    integer :: net
    integer :: ngrn1
    integer :: ngrn2
    integer :: ngt
    integer :: nhydr
    integer :: no2gaq
    integer :: npt
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
    logical :: q6mode

    character(len=48) :: uspec(nstmax)
    character(len=24) :: ugexmo(netmax)
    character(len=24) :: uphase(nptmax)
    character(len=8) :: ugexj(jetmax,netmax)

    real(kind=8) :: acflg(nstmax)
    real(kind=8) :: acflgo(nstmax)
    real(kind=8) :: act(nstmax)
    real(kind=8) :: actlg(nstmax)
    real(kind=8) :: azero(natmax)
    real(kind=8) :: a3bars(natmax)
    real(kind=8) :: cco2(5)
    real(kind=8) :: cdrs(ndrsmx)
    real(kind=8) :: cegexs(ietmax,jetmax,netmax)
    real(kind=8) :: cgexj(jetmax,netmax)
    real(kind=8) :: conc(nstmax)
    real(kind=8) :: conclg(nstmax)
    real(kind=8) :: cpgexs(ietmax,jetmax,netmax)
    real(kind=8) :: egexjc(jetmax,netmax)
    real(kind=8) :: egexjf(jetmax,netmax)
    real(kind=8) :: egexs(ietmax,jetmax,netmax)
    real(kind=8) :: fsort(ngtmax)
    real(kind=8) :: fugac(ngtmax)
    real(kind=8) :: fugalg(ngtmax)
    real(kind=8) :: loph(nptmax)
    real(kind=8) :: losp(nstmax)
    real(kind=8) :: lsort(nstmax)
    real(kind=8) :: mgext(jetmax,netmax)
    real(kind=8) :: moph(nptmax)
    real(kind=8) :: mosp(nstmax)
    real(kind=8) :: mrgexs(ietmax,jetmax,netmax)
    real(kind=8) :: mtb(nbtmax)
    real(kind=8) :: xbar(nstmax)
    real(kind=8) :: xbarlg(nstmax)
    real(kind=8) :: xlks(nstmax)
    real(kind=8) :: zchar(nstmax)
    real(kind=8) :: zchcu6(nstmax)
    real(kind=8) :: zchsq2(nstmax)
    real(kind=8) :: zgexj(jetmax,netmax)
    real(kind=8) :: zvclg1(kmax)
    real(kind=8) :: zvec1(kmax)

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
    real(kind=8) :: eps100
    real(kind=8) :: fje
    real(kind=8) :: fjeo
    real(kind=8) :: fo2
    real(kind=8) :: fo2lg
    real(kind=8) :: fxi
    real(kind=8) :: fxio
    real(kind=8) :: omega
    real(kind=8) :: omeglg
    real(kind=8) :: press
    real(kind=8) :: sigmam
    real(kind=8) :: sigmmo
    real(kind=8) :: tempk
    real(kind=8) :: xbarw
    real(kind=8) :: xbarwc
    real(kind=8) :: xbrwlc
    real(kind=8) :: xbrwlg

    ! Local variable declarations.
    integer :: jlen
    integer :: kcol
    integer :: nb
    integer :: ns

    logical :: qxbarw

    character(len=56) :: uspn56

    real(kind=8) :: zx1
    real(kind=8) :: zx2

    real(kind=8) :: texp
    real(kind=8) :: tlg

    data qxbarw/.false./

    if (iodb(3) .ge. 1) then
        write (noutpt,1000)
1000 format(/11x,' --- Starting initializing calculations ---',/)
    end if

    ! Initialize log mole fraction of water in aqueous solution,
    ! using a plausible value. This will later be corrected as needed.
    xbarwc = 1.0
    xbrwlc = 0.
    xbar(narn1) = xbarwc
    xbarlg(narn1) = xbrwlc

    ! Calculate the concentrations, etc., of basis and non-basis
    ! species. Here the activity coefficients are all zero.
    call ncmpex(acflg,act,actlg,cdrs,cegexs,cgexj,conc,conclg,cpgexs,egexjc,egexjf,egexs,eps100,fo2,fo2lg,fsort,fugac,fugalg,iern1,iern2,ietmax,ifrn1,ifrn2,igas,igstak,iindx1,ilrn1,ilrn2,imrn1,imrn2,istack,ixrn1,ixrn2,jcsort,jern1,jern2,jetmax,jflag,jgext,jgsort,jgstak,jjsort,jpflag,jsflag,jsitex,jssort,jstack,kbt,kdim,kelect,kmax,km1,ko2gaq,kwater,kxt,loph,losp,lsort,mgext,mrgexs,mtb,moph,mosp,narn1,narn2,nbasp,nbt,nbtmax,ncmpr,ndrs,ndrsmx,ndrsr,nelect,nern1,nern2,netmax,ngexsa,ngext,ngrn1,ngrn2,ngt,ngtmax,noutpt,no2gaq,nphasx,npt,nptmax,nst,nstmax,nttyo,omega,omeglg,press,qxbarw,q6mode,ugexj,ugexmo,uphase,uspec,xbar,xbarlg,xbarw,xbarwc,xbrwlc,xbrwlg,xlks,zchar,zgexj,zvclg1,zvec1)

    ! Clamp estimates of concentrations and numbers of moles of
    ! dependent aqueous species.
    do ns = narn1,narn2
        if (jflag(ns) .ge. 30) then
            conc(ns) = 0.
            mosp(ns) = 0.
            conclg(ns) = -99999.
            losp(ns) = -99999.
        end if
    end do

    ! Make the first estimates of SIGMA m, ionic strength, and J taking
    ! into account the concentrations of only data file basis species.
    ! This generally provides a set of safe values, as the computed
    ! concentrations of non-basis species at this point could be
    ! extremely large.
    sigmam = 0.
    fxi = 0.
    fje = 0.

    do nb = 1,nbt
        ns = nbaspd(nb)

        if (ns.gt.narn1 .and. ns.le.narn2) then
            sigmam = sigmam + conc(ns)
            fxi = fxi + conc(ns)*zchsq2(ns)
            fje = fje + conc(ns)*zchcu6(ns)
        end if
    end do

    ! Make the first real estimate of the mole fraction of water in
    ! aqueous solution.
    xbarwc = omega/(omega + sigmam)
    xbrwlc = tlg(xbarwc)
    xbarw = xbarwc
    xbrwlg = xbrwlc
    xbar(narn1) = xbarwc
    xbarlg(narn1) = xbrwlc

    ! Make the first estimates of the aqueous species activity
    ! coefficients.
    ! Calling sequence substitutions:
    !   acflg for acflgc
    call gcoeff(abar,acflg,actwlc,adh,adhh,adhv,al10,aphi,azero,a3bar,a3bars,bdh,bdhh,bdhv,bdot,bdoth,bdotv,cco2,conc,delam,dgpit,dpelm,dpslm,dselm,elam,fje,fxi,gpit,ielam,ifcphi1,ifcphi2,ifnnn,ifn2n,ifpsi1,ifpsi2,ifzeta,ilcphi1,ilcphi2,ilnnn,iln2n,ilpsi1,ilpsi2,ilzeta,insgf,iopg,ipbtmx,izmax,jcsort,nalpha,napmax,napt,narn1,narn2,natmax,nazmmx,nazpmx,nchlor,nhydr,nmut,nmutmx,nmux,nmxi,nmxmax,nmxx,nopgmx,noutpt,nslt,nsltmx,nslx,nstmax,nsxi,nsxmax,nsxx,nttyo,omega,palpha,pelm,pmu,press,pslamn,pslm,qhawep,qpit75,selm,sigmam,tempk,uspec,xbarwc,xbrwlc,zchar,zchsq2,zchcu6)

    ! Make the first estimates of the exchanger species activity
    ! coefficients.
    ! Calling sequence substitutions:
    !   acflg for acflgc
    call lamgex(acflg,cgexj,jern1,jern2,jetmax,jgext,net,netmax,nstmax,xbarlg)

    if (iodb(3) .ge. 1) then
        ! Print the attempted phase assemblage.
        write (noutpt,1010)
1010 format(/' Initial phase assemblage:',/)

        do kcol = 1,kbt
            nb = iindx1(kcol)
            ns = nbaspd(nb)

            ! Calling sequence substitutions:
            !   uspec(ns) for unam48
            call fmspnx(jlen,uspec(ns),uspn56)
            write (noutpt,1020) kcol,uspn56(1:jlen)
1020 format(2x,i3,2x,a)
        end do

        do kcol = km1,kxt
            ns = iindx1(kcol)

            ! Calling sequence substitutions:
            !   uspec(ns) for unam48
            call fmspnm(jlen,uspec(ns),uspn56)
            write (noutpt,1020) kcol,uspn56(1:jlen)
        end do

        write (noutpt,1030)
1030 format(1x)
    end if

    if (iodb(3) .ge. 2) then
        write (noutpt,1050)
1050 format(/16x,'--- Initialization Summary ---',//2x,'kcol   Name',32x,'zvclg1      zvec1',/)

        do kcol = 1,kbt
            nb = iindx1(kcol)
            ns = nbasp(nb)
            zx1 = zvclg1(kcol)
            zx2 = texp(zx1)

            ! Calling sequence substitutions:
            !   uspec(ns) for unam48
            call fmspnx(jlen,uspec(ns),uspn56)
            write (noutpt,1060) kcol,uspn56,zx1,zx2
1060 format(1x,i4,2x,a32,2x,f10.4,2x,1pe12.5)
        end do

        do kcol = km1,kxt
            ns = iindx1(kcol)
            zx1 = zvclg1(kcol)
            zx2 = texp(zx1)

            ! Calling sequence substitutions:
            !   uspec(ns) for unam48
            call fmspnm(jlen,uspec(ns),uspn56)
            write (noutpt,1060) kcol,uspn56,zx1,zx2
        end do

        write (noutpt,1030)
    end if

    if (iodb(3) .ge. 1) then
        write (noutpt,1120) sigmam,fxi,fje
1120 format(/8x,'sigmam= ',1pe12.5,/11x,'fxi= ',e12.5,/11x,'fje= ',e12.5,/)
    end if

    if (iodb(3) .ge. 4) then
        write (noutpt,1130)
1130 format(/7x,'Species',20x,'gamma',/)

        do ns = narn1,narn2
            write (noutpt,1140) uspec(ns),acflg(ns)
1140 format(5x,a24,3x,1pe12.5)
        end do

        write (noutpt,1030)
    end if

    sigmmo = sigmam
    fxio = fxi
    fjeo = fje

    do ns = narn1,narn2
        acflgo(ns) = acflg(ns)
    end do

    if (iodb(3) .ge. 1) then
        write (noutpt,1200)
1200 format(/11x,' --- Finished initializing calculations ---',/)
    end if

999 continue
end subroutine exivar