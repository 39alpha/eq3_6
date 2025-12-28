subroutine arrset(aamatr,abar,acflg,acflgo,act,actlg,adh,adhh,adhv,adhfs,adhfsx,advfs,advfsx,afcnst,alpha,al10,amtb,aphi,avcnst,azero,a3bar,a3bars,axhfs,axhfsx,axlks,axlksx,axvfs,axvfsx,bacfmx,bbig,bdh,bdhh,bdhv,bdot,bdoth,bdotv,beta,betamx,bfac,bgamx,bneg,bpx,bsigmm,bfje,bfxi,cco2,cdrs,cdrsx,cdrtw,cdrw,cegexs,cgexj,cjbasp,cnufac,conc,conclg,coval,cpgexs,csts,delam,delvec,dgpit,dhfs,dlogxw,dpelm,dpslm,dselm,dvfs,efac,egexjc,egexjf,egexs,eh,ehfac,elam,eps100,fje,fjeo,fo2,fo2lg,fsort,fugac,fugalg,fxi,fxio,gmmatr,gpit,ibetmx,ibpxt,ibswx,iction,iebal,ielam,iern1,iern2,ifcphi1,ifcphi2,ifnnn,ifn2n,ifpsi1,ifpsi2,ifrn1,ifrn2,ifzeta,igas,igstak,iindx1,ilcphi1,ilcphi2,ilnnn,iln2n,ilpsi1,ilpsi2,ilrn1,ilrn2,ilzeta,imrn1,imrn2,insgf,iopg,iodb,iopt,ipch,ipcv,ipivot,ipndx1,irdxc3,istack,ixbasp,ixrn1,ixrn2,izmax,jcsort,jern1,jern2,jflag,jgext,jgsort,jgstak,jjndex,jjsort,jpflag,jsflag,jsitex,jsol,jssort,jstack,ka1,kat,kbt,kct,kction,kdim,kebal,kelect,ker,ke1,ket,khydr,kkndex,km1,ko2gaq,kwater,kx1,kxt,loph,losp,lsort,mgext,moph,mosp,mrgexs,mtb,nalpha,napt,narn1,narn2,narxt,nbasp,nbaspd,nbaspx,nbt,nbtd,nbti,nbw,nchlor,ncmpr,ncosp,nct,ndecsp,ndrs,ndrsx,ndrsr,ndrsrd,ndrsrx,nelect,nern1,nern2,net,nfac,ngexsa,ngext,ngrn1,ngrn2,ngt,nhydr,nhydx,nmut,nmux,nmxi,nmxx,noutpt,no2gaq,nphasx,npt,nredox,nslt,nslx,nst,nsts,nstsr,nsxi,nsxx,ntfx,ntfxt,ntpr,nttyo,omega,omeglg,palpha,pe,pelm,pmu,presg,press,pslamn,pslm,qbassw,qchlor,qhawep,qpit75,qredox,q6mode,rhsvec,selm,sigmam,sigmmo,smp100,tempc,tempk,tfx,ubacmx,ubbig,ubgamx,ubneg,ubetmx,ucospi,ugexj,ugexmo,ujflls,uphase,uspec,uzvec1,weight,wfac,xbar,xbarlg,xbarw,xbarwc,xbrwlc,xbrwlg,xhfs,xlke,xlks,xvfs,zchar,zchsq2,zchcu6,zgexj,zvclg1,zvec1)
    !! This subroutine builds the iindx1 array. It sets up the matrix
    !! structure for hybrid Newton-Raphson iteration and computes
    !! initial values for the iteration variables. It optimizes these
    !! in preparation for hybrid Newton-Raphson iteration. The variable
    !! ker is returned as 0 if all went well, as 1 if the input
    !! constraints look suspiciously poor, and as 2 if they look
    !! really bad.
    !! This subroutine is somewhat analogous to EQ6/optmzr.f, which also
    !! performs optimization of iteration variables prior to hybrid
    !! Newton-Raphson iteration.
    !! This subroutine is called by:
    !!   EQ3NR/eq3nr.f
    !! Principal input:
    !!   narn1  = start of the range of aqueous species; this is
    !!              also the index of the solvent, water
    !!   narn2  = end of the range of aqueous species
    !! Principal output:
    !!   conc   = array of species concentrations
    !!   iindx1 = array of indices of components appearing as
    !!              matrix variables
    implicit none

    include 'eqlib/eqldv.h'

    ! Calling sequence variable declarations.
    integer :: noutpt
    integer :: nttyo

    integer :: ibpxt(nxtmax)
    integer :: ibswx(nbtmax)
    integer :: iction(nbtmax)
    integer :: igstak(ngtmax)
    integer :: iindx1(kmax)
    integer :: ipndx1(kmax)
    integer :: insgf(natmax)
    integer :: iodb(nodbmx)
    integer :: iopg(nopgmx)
    integer :: iopt(noptmx)
    integer :: ipivot(kmax)
    integer :: istack(nstmax)
    integer :: ixbasp(nbtmax)
    integer :: jcsort(nstmax)
    integer :: jern1(jetmax,netmax)
    integer :: jern2(jetmax,netmax)
    integer :: jflag(nstmax)
    integer :: jgext(netmax)
    integer :: jgsort(ngtmax)
    integer :: jgstak(ngtmax)
    integer :: jjndex(nbtmax)
    integer :: jjsort(nstmax)
    integer :: jsflag(nstmax)
    integer :: jsitex(nstmax)
    integer :: jsol(nxtmax)
    integer :: jpflag(nptmax)
    integer :: jssort(nstmax)
    integer :: jstack(nstmax)
    integer :: kction(nbtmax)
    integer :: kkndex(nbtmax)

    integer :: narxt(ntprmx)
    integer :: nbasp(nbtmax)
    integer :: nbaspd(nbtmax)
    integer :: nbaspx(nbtmax)
    integer :: ncmpr(2,nptmax)
    integer :: ncosp(nbtmax)
    integer :: ndecsp(nbtmax)
    integer :: ndrs(ndrsmx)
    integer :: ndrsx(ndrsmx)
    integer :: ndrsr(2,nstmax)
    integer :: ndrsrd(2,nstmax)
    integer :: ndrsrx(2,nstmax)
    integer :: nfac(nbtmax)
    integer :: ngexsa(ietmax,jetmax,netmax)
    integer :: ngext(jetmax,netmax)
    integer :: nphasx(nstmax)
    integer :: nsts(nstsmx)
    integer :: nstsr(2,nstmax)
    integer :: ntfx(ntfxmx)

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

    integer :: ibetmx
    integer :: iebal
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
    integer :: ipch
    integer :: ipcv
    integer :: irdxc3
    integer :: ixrn1
    integer :: ixrn2
    integer :: izmax
    integer :: ka1
    integer :: kat
    integer :: kbt
    integer :: kct
    integer :: kdim
    integer :: kebal
    integer :: kelect
    integer :: ker
    integer :: ke1
    integer :: ket
    integer :: khydr
    integer :: km1
    integer :: ko2gaq
    integer :: kwater
    integer :: kx1
    integer :: kxt
    integer :: narn1
    integer :: narn2
    integer :: nbt
    integer :: nbtd
    integer :: nbti
    integer :: nbw
    integer :: nchlor
    integer :: nct
    integer :: nelect
    integer :: nern1
    integer :: nern2
    integer :: net
    integer :: ngrn1
    integer :: ngrn2
    integer :: ngt
    integer :: nhydr
    integer :: nhydx
    integer :: no2gaq
    integer :: npt
    integer :: nredox
    integer :: nst
    integer :: ntfxt
    integer :: ntpr

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

    logical :: qbassw
    logical :: qchlor
    logical :: qhawep
    logical :: qpit75
    logical :: qredox
    logical :: q6mode

    character(len=48) :: ucospi(nbtmax)
    character(len=48) :: uspec(nstmax)
    character(len=48) :: uzvec1(kmax)
    character(len=48) :: ubacmx
    character(len=48) :: ubbig
    character(len=48) :: ubetmx
    character(len=48) :: ubgamx
    character(len=48) :: ubneg
    character(len=32) :: ujflls(0:njfmax)
    character(len=24) :: ugexmo(netmax)
    character(len=24) :: uphase(nptmax)
    character(len=8) :: ugexj(jetmax,netmax)

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

    real(kind=8) :: xhfs(nstmax)
    real(kind=8) :: xlks(nstmax)
    real(kind=8) :: xvfs(nstmax)

    real(kind=8) :: dhfs(ipchmx,nstmax)
    real(kind=8) :: dvfs(ipcvmx,nstmax)

    real(kind=8) :: axhfs(narxmx,ntprmx,nstmax)
    real(kind=8) :: axhfsx(narxmx,ntprmx,nstmax)
    real(kind=8) :: axlks(narxmx,ntprmx,nstmax)
    real(kind=8) :: axlksx(narxmx,ntprmx,nstmax)
    real(kind=8) :: axvfs(narxmx,ntprmx,nstmax)
    real(kind=8) :: axvfsx(narxmx,ntprmx,nstmax)

    real(kind=8) :: adhfs(narxmx,ntprmx,ipchmx,nstmax)
    real(kind=8) :: adhfsx(narxmx,ntprmx,ipchmx,nstmax)
    real(kind=8) :: advfs(narxmx,ntprmx,ipcvmx,nstmax)
    real(kind=8) :: advfsx(narxmx,ntprmx,ipcvmx,nstmax)

    real(kind=8) :: aamatr(kmax,kmax)
    real(kind=8) :: acflg(nstmax)
    real(kind=8) :: acflgo(nstmax)
    real(kind=8) :: act(nstmax)
    real(kind=8) :: actlg(nstmax)
    real(kind=8) :: alpha(kmax)
    real(kind=8) :: amtb(nbtmax)
    real(kind=8) :: azero(natmax)
    real(kind=8) :: a3bars(natmax)
    real(kind=8) :: beta(kmax)
    real(kind=8) :: bfac(nbtmax)
    real(kind=8) :: bpx(ibpxmx,nxtmax)
    real(kind=8) :: cco2(5)
    real(kind=8) :: cdrs(ndrsmx)
    real(kind=8) :: cdrsx(ndrsmx)
    real(kind=8) :: cdrtw(nstmax)
    real(kind=8) :: cdrw(nstmax)
    real(kind=8) :: cegexs(ietmax,jetmax,netmax)
    real(kind=8) :: cgexj(jetmax,netmax)
    real(kind=8) :: cjbasp(nbtmax)
    real(kind=8) :: cnufac(nstmax)
    real(kind=8) :: conc(nstmax)
    real(kind=8) :: conclg(nstmax)
    real(kind=8) :: coval(nbtmax)
    real(kind=8) :: cpgexs(ietmax,jetmax,netmax)
    real(kind=8) :: csts(nstsmx)
    real(kind=8) :: delvec(kmax)
    real(kind=8) :: dlogxw(nbtmax)
    real(kind=8) :: efac(nbtmax)
    real(kind=8) :: egexjc(jetmax,netmax)
    real(kind=8) :: egexjf(jetmax,netmax)
    real(kind=8) :: egexs(ietmax,jetmax,netmax)
    real(kind=8) :: fsort(ngtmax)
    real(kind=8) :: fugac(ngtmax)
    real(kind=8) :: fugalg(ngtmax)
    real(kind=8) :: gmmatr(kmax,kmax)
    real(kind=8) :: loph(nptmax)
    real(kind=8) :: losp(nstmax)
    real(kind=8) :: lsort(nstmax)

    real(kind=8) :: mgext(jetmax,netmax)
    real(kind=8) :: moph(nptmax)
    real(kind=8) :: mosp(nstmax)
    real(kind=8) :: mrgexs(ietmax,jetmax,netmax)
    real(kind=8) :: mtb(nbtmax)
    real(kind=8) :: rhsvec(kmax)
    real(kind=8) :: tfx(ntfxmx)
    real(kind=8) :: weight(nstmax)
    real(kind=8) :: wfac(iktmax,nxtmax)
    real(kind=8) :: xbar(nstmax)
    real(kind=8) :: xbarlg(nstmax)
    real(kind=8) :: zchar(nstmax)
    real(kind=8) :: zchsq2(nstmax)
    real(kind=8) :: zchcu6(nstmax)
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
    real(kind=8) :: afcnst
    real(kind=8) :: al10
    real(kind=8) :: avcnst
    real(kind=8) :: a3bar
    real(kind=8) :: bacfmx
    real(kind=8) :: bbig
    real(kind=8) :: betamx
    real(kind=8) :: bfje
    real(kind=8) :: bfxi
    real(kind=8) :: bgamx
    real(kind=8) :: bneg
    real(kind=8) :: bsigmm
    real(kind=8) :: btmxoe
    real(kind=8) :: btmxoo
    real(kind=8) :: eh
    real(kind=8) :: ehfac
    real(kind=8) :: eps100
    real(kind=8) :: fje
    real(kind=8) :: fjeo
    real(kind=8) :: fo2
    real(kind=8) :: fo2lg
    real(kind=8) :: fxi
    real(kind=8) :: fxio
    real(kind=8) :: omega
    real(kind=8) :: omeglg
    real(kind=8) :: pe
    real(kind=8) :: press
    real(kind=8) :: sigmam
    real(kind=8) :: sigmmo
    real(kind=8) :: smp100
    real(kind=8) :: tempc
    real(kind=8) :: tempk
    real(kind=8) :: xbarw
    real(kind=8) :: xbarwc
    real(kind=8) :: xbrwlc
    real(kind=8) :: xbrwlg

    ! Local variable declarations.
    integer :: iter
    integer :: jfl
    integer :: jlen
    integer :: jlen1
    integer :: jlen2
    integer :: j2
    integer :: kb
    integer :: kc
    integer :: kcol
    integer :: kount
    integer :: krow
    integer :: nb
    integer :: nb2
    integer :: ncycle
    integer :: ncylim
    integer :: negbfc
    integer :: nloop
    integer :: nlopmx
    integer :: npass
    integer :: nplim
    integer :: nr1
    integer :: ns
    integer :: nse
    integer :: nswtch
    integer :: ns1
    integer :: ns2

    integer :: ilnobl

    logical :: qabsw
    logical :: qawfix
    logical :: qbswx
    logical :: qcfxi
    logical :: qcgam
    logical :: qcsigm
    logical :: qloop
    logical :: qpracf
    logical :: qtestc
    logical :: qtestp
    logical :: qxbarw

    character(len=56) :: uspn56
    character(len=56) :: usp156
    character(len=56) :: usp256
    character(len=32) :: ujtp
    character(len=24) :: ux24

    real(kind=8) :: presg
    real(kind=8) :: xlke

    real(kind=8) :: av
    real(kind=8) :: azdel
    real(kind=8) :: betfnc
    real(kind=8) :: bxecor
    real(kind=8) :: bxp1
    real(kind=8) :: cecorr
    real(kind=8) :: chfacf
    real(kind=8) :: chfsgm
    real(kind=8) :: cx
    real(kind=8) :: cxl
    real(kind=8) :: cxn
    real(kind=8) :: cxo
    real(kind=8) :: dx
    real(kind=8) :: fjec
    real(kind=8) :: fxic
    real(kind=8) :: lx
    real(kind=8) :: rlxgam
    real(kind=8) :: sigmmc
    real(kind=8) :: sigza
    real(kind=8) :: sigzc
    real(kind=8) :: sigzi
    real(kind=8) :: sigzm
    real(kind=8) :: stx
    real(kind=8) :: tfxc
    real(kind=8) :: tolbig
    real(kind=8) :: tolbtf
    real(kind=8) :: tolgpt
    real(kind=8) :: tolneg
    real(kind=8) :: tolxpt
    real(kind=8) :: tolzpt
    real(kind=8) :: xecorr
    real(kind=8) :: zdel
    real(kind=8) :: zx1
    real(kind=8) :: zx2

    real(kind=8) :: coefdr
    real(kind=8) :: texp
    real(kind=8) :: tlg

    ! The following are iteration limits:
    !   nlopmx = the maximum number of auto basis switching loops
    !   nplim  = the maximum number of passes
    !   ncylim = the maximum number of cycles
    ! Passes refine estimates of the ionic strength, etc., the
    ! activity of water, and activity coefficients of aqueous species.
    ! Cycles are embedded in passes. They refine estimates of species
    ! concentrations before new estimates of ionic strength, etc.,
    ! are made.
    data nlopmx /12/,nplim  /7/,ncylim /15/

    ! The following are tolerance parameters for the optimization:
    data tolbtf /0.1/
    data tolzpt /0.25/
    data tolbig /0.5/,tolneg /-0.1/
    data tolxpt /0.5/,tolgpt /0.1/

    ! The following is needed for ncmpex.f, but is only relevant to EQ6.
    data qxbarw/.false./

    write (noutpt,1000)
    write (nttyo,1000)
1000 format(/' Starting Pre-Newton-Raphson Optimization.',/)

    qbswx = .false.

    do nb = 1,nbtmax
        kkndex(nb) = 0
        kction(nb) = 0
    end do

    bbig = 0.
    bneg = 0.
    ubbig = 'None'
    ubneg = 'None'

    bgamx = 0.
    ubgamx = 'None'

    ! Set up the structure of the iteration matrix. Its contents may be
    ! altered subsequently by automatic basis switching.
    ! Build the iindx1 array.
    kebal = 0
    kwater = 0
    khydr = 0
    kelect = 0
    ko2gaq = 0

    ka1 = 0
    kat = 0
    ke1 = 0
    ket = 0

    kb = 0
    kc = 0

    do nb = 1,nbt
        ns = nbasp(nb)
        jfl = jflag(ns)

        if (jfl.ne.-1 .and. jfl.ne.30) then
            kb = kb + 1
            iindx1(kb) = nb
            uzvec1(kb) = uspec(ns)

            if (nb .eq. iebal) then
                kebal = kb
            end if

            if (ns .eq. narn1) then
                kwater = kb
            end if

            if (ns .eq. nhydr) then
                khydr = kb
            end if

            if (ns .eq. nelect) then
                kelect = kb
            end if

            if (ns .eq. no2gaq) then
                ko2gaq = kb
            end if
        end if
    end do

    kbt = kb
    kdim = kbt

    ! XXX
    if (qchlor) then
        kct = nct -1
    else
        kct = nct
    end if

    ! Fill the kction array. This array is used to mark the columns
    ! belonging to basis species which are involved in constraints
    ! placed on other basis species.
    do nb = 1,nbt
        ns1 = ncosp(nb)

        do kcol = 1,kbt
            nb2 = iindx1(kcol)
            ns2 = nbasp(nb2)

            if (ns2 .eq. ns1) then
                kction(nb) = kcol
                go to 100
            end if
        end do

100 continue
    end do

    if (iodb(3) .ge. 2) then
        ! Print the active data file basis set.
        write (noutpt,1010)
1010 format(/16x,'--- Active Data File Basis Set ---',//2x,'krow   Name',30x,'Constraint',/)

        do krow = 1,kdim
            nb = iindx1(krow)
            ns = nbaspd(nb)
            ujtp = 'Defining equation'

            if (nb .eq. iebal) then
                ujtp = 'Electrical balance'
            else if (ns .eq. no2gaq) then
                if ((irdxc3 .eq. -1) .or. (irdxc3 .eq. -2)) then
                    ujtp = 'Eh'
                else if (irdxc3 .eq. 1) then
                    ujtp = 'Aqueous redox reaction'
                end if
            else
                jfl = jflag(ns)
                ujtp = ujflls(jfl)
            end if

            j2 = ilnobl(ujtp)

            ! Calling sequence substitutions:
            !   uspec(ns) for unam48
            call fmspnx(jlen,uzvec1(krow),uspn56)
            write (noutpt,1020) krow,uspn56,ujtp(1:j2)
1020 format(1x,i4,2x,a32,2x,a)
        end do

        write (noutpt,1030)
1030 format(/1x)
    end if

    ! The label below is a return point after an automatic basis switch.
    ! Here nloop is the loop counter for auto basis switching.
    nloop = -1

    qloop = .true.
200 continue
    nloop = nloop + 1

    if (iodb(3) .ge. 1) then
        write (noutpt,1040) nloop
    end if

1040 format(6x,'nloop= ',i2)

    if (iodb(3) .ge. 2) then
        if (qbassw) then
            ! Print the computational basis set.
            write (noutpt,1050)
1050 format(16x,'--- Basis Set Changes ---',//2x,'krow   Data File',25x,'Current',/)

            do krow = 1,kbt
                nb = iindx1(krow)
                ns1 = nbaspd(nb)
                ns2 = nbasp(nb)

                if (ns1 .ne. ns2) then
                    ! Calling sequence substitutions:
                    !   jlen1 for jlen
                    !   uspec(ns1) for unam48
                    !   usp156 for uspn56
                    call fmspnx(jlen1,uspec(ns1),usp156)

                    ! Calling sequence substitutions:
                    !   jlen2 for jlen
                    !   uspec(ns2) for unam48
                    !   usp256 for uspn56
                    call fmspnx(jlen2,uspec(ns2),usp256)
                    jlen2 = min(jlen2,32)
                    write (noutpt,1060) krow,usp156,usp256(1:jlen2)
1060 format(1x,i4,2x,a32,2x,a)
                end if
            end do

            write (noutpt,1030)
        end if
    end if

    ! Initialize concentrations and masses of the active basis species.
    do kcol = 1,kbt
        nb = iindx1(kcol)
        nse = nbasp(nb)
        jfl = jflag(nse)

        if (jfl .eq. -1) then
            conc(nse) = 0.
        else if (jfl.ge.0 .and. jfl.le.3) then
            ! Concentrations.
            ns = nbaspd(nb)

            if (nse .eq. ns) then
                conc(nse) = coval(nb)
            else
                cx = coefdr(cdrs,ndrs,ndrsmx,ndrsr,nse,ns,nstmax)
                nr1 = ndrsr(1,ns)
                stx = -cx/cdrs(nr1)
                conc(nse) = stx*coval(nb)
            end if
        else if (jfl.ge.7 .and. jfl.le.11) then
            ! Alkalinity.
            ux24 = uspec(nse)(1:24)
            tfxc = 1.

            if (ux24(1:6) .eq. 'HCO3- ') then
                tfxc = 1.
            end if

            if (ux24(1:6) .eq. 'CO3-- ') then
                tfxc = 2.
            end if

            conc(nse) = coval(nb)/tfxc
        else if (jfl .eq. 16) then
            ! Activity.
            cx = coval(nb)
            conc(nse) = texp(cx)
        else if (jfl.eq.19 .or. jfl.eq.20) then
            ! pX and pH.
            cx = -coval(nb)
            conc(nse) = texp(cx)
        else if (jfl.eq.22 .or. jfl.eq.23) then
            ! pmX and pmH.
            cx = -coval(nb)
            conc(nse) = texp(cx)
        else
            ! All other cases.
            conc(nse) = 1.e-7
        end if
    end do

    conc(narn1) = 0.

    if (nelect .gt. 0) then
        conc(nelect) = 0.
    end if

    if (no2gaq .gt. 0) then
        conc(no2gaq) = 0.
    end if

    ! Calculate the charge imbalance.
    zdel = 0.

    do nb = 1,nbt
        ns = nbasp(nb)

        if (ns.ge.narn1 .and. ns.le.narn2) then
            zdel = zdel + zchar(ns)*conc(ns)
        end if
    end do

    azdel = abs(zdel)

    ! Calculate a starting value for the SUM(i) m(i) function
    ! (sigmam). Treat the calculated charge imbalance among the
    ! basis species as the equivalent of a monovalent ion.
    call csigm(conc,jcsort,narn1,narn2,nstmax,sigmmc)
    sigmam = sigmmc + azdel

    ! Calculate a starting value for the ionic strength (fxi).
    ! Treat the calculated charge imbalance among the basis species
    ! as the equivalent of a monovalent ion.
    call cfxi(conc,fxic,jcsort,narn1,narn2,nstmax,zchsq2)
    fxi = fxic + 0.5*azdel

    ! Calculate a starting value for the J electrostatic moment
    ! function (fje). Treat the calculated charge imbalance among
    ! the basis species as the equivalent of a monovalent ion.
    call cfje(conc,fjec,jcsort,narn1,narn2,nstmax,zchcu6)
    fje = fjec + (-zdel/6.)

    ! Calculate the activity coefficients of aqueous species.
    ! Note that this also gets the starting value for the mole
    ! fraction of water.
    ! Calling sequence substitutions:
    !   acflg for acflgc
    call gcoeff(abar,acflg,actwlc,adh,adhh,adhv,al10,aphi,azero,a3bar,a3bars,bdh,bdhh,bdhv,bdot,bdoth,bdotv,cco2,conc,delam,dgpit,dpelm,dpslm,dselm,elam,fje,fxi,gpit,ielam,ifcphi1,ifcphi2,ifnnn,ifn2n,ifpsi1,ifpsi2,ifzeta,ilcphi1,ilcphi2,ilnnn,iln2n,ilpsi1,ilpsi2,ilzeta,insgf,iopg,ipbtmx,izmax,jcsort,nalpha,napmax,napt,narn1,narn2,natmax,nazmmx,nazpmx,nchlor,nhydr,nmut,nmutmx,nmux,nmxi,nmxmax,nmxx,nopgmx,noutpt,nslt,nsltmx,nslx,nstmax,nsxi,nsxmax,nsxx,nttyo,omega,palpha,pelm,pmu,press,pslamn,pslm,qhawep,qpit75,selm,sigmam,tempk,uspec,xbarwc,xbrwlc,zchar,zchsq2,zchcu6)

    ! Calculate the activity coefficients of exchanger species.
    ! Calling sequence substitutions:
    !   acflg for acflgc
    call lamgex(acflg,cgexj,jern1,jern2,jetmax,jgext,net,netmax,nstmax,xbarlg)

    ! Initialize the mole fraction of water.
    xbrwlg = xbrwlc
    xbarw = xbarwc

    ! Copy the ionic strength, etc., and the activity coefficients.
    sigmmo = sigmam
    fxio = fxi
    fjeo = fje

    do ns = narn1,narn2
        acflgo(ns) = acflg(ns)
    end do

    ! Load the entries of the conclg array corresponding
    ! to active basis species.
    do kcol = 1,kbt
        nb = iindx1(kcol)
        ns = nbasp(nb)
        cx = conc(ns)
        conclg(ns) = tlg(cx)
    end do

    ! Determine whether the constraints fix the activity of water.
    ! If so, the variable qawfix is set to .true.
    call dawfix(aamatr,cdrs,eps100,gmmatr,iindx1,iodb,irdxc3,jflag,jjndex,kbt,kkndex,kmax,narn1,nbasp,nbtmax,ncosp,ndrs,ndrsmx,ndrsr,nelect,nhydr,nodbmx,no2gaq,noutpt,nstmax,qawfix,uspec)

    ! XX   Need new coding to deal with phases assemblages that fix a(w).
    ! XX   Such assemblages are now merely trapped.
    !      Coding to deal with phase assemablages that fix the activity of
    !      water has not yet been implemented. Stop if this condition has
    !      been detected. The needed new coding is in EQ3NR/arrset.f and
    !      EQ3NR/ arrsim.f. If the following trap were not in place, the
    !      matrix constructed below in this subroutine would be singular.
    !      There should be no problem, however, in the coding for the hybrid
    !      Newton-Raphson method.
    if (qawfix) then
        write (noutpt,1100)
        write (nttyo,1100)
1100 format(/' * Error - (EQ3NR/arrset) The phase assemblage',/7x,'corresponding to the specified solubility constraints',/7x,'fixes the activity of water. This code is presently',/7x,'unable to solve problems of this type.')

        stop
    end if

    ! Here npass is the pass counter.
    npass = 0

    ! The label below is a return point for subsequent passes. A pass
    ! is an adjustment for the ionic strength, etc., the activity of
    ! water, and the activity coefficients of the solute species.
210 continue
    npass = npass + 1

    ! Note:
    !   betfnc = convergence function
    !   negbfc = the number of successive iterations that the
    !            convergence function betfnc has been zero or negative
    betfnc = 0.
    negbfc = 0
    btmxoe = 0.
    btmxoo = 0.

    if (iodb(3) .ge. 1) then
        write (noutpt,1110) npass
1110 format(/11x,'npass= ',i2)

        write (noutpt,1120) sigmam,fxi,fje,xbrwlc,xbarwc
1120 format(/13x,'sigmam= ',1pe12.5,/13x,'fxi= ',1pe12.5,/13x,'fje= ',1pe12.5,//13x,'xbrwlc= ',0pf9.5,/13x,'xbarwc= ',1pe12.5,/)
    end if

    ! Here ncycle is the cycle counter.
    ncycle = 0

    ! The label below is a return point for beginning a new cycle.
    ! A cycle is an structure within a pass in which the concentrations
    ! of the basis species are adjusted, while the ionic strength, etc.,
    ! and the activity coefficients are held constant.
220 continue
    ncycle = ncycle + 1

    if (iodb(3) .ge. 1) then
        write (noutpt,1200) ncycle
    end if

1200 format(16x,'ncycle= ',i2)

    ker = 0
    qabsw = iopt(11).ge.1 .and. npass.eq.1 .and. ncycle.eq.1

    ! Set up the zvclg1 and actlg array entries for the basis species
    ! whose concentrations do not have to be estimated simultaneously.
    do kcol = 1,kbt
        nb = iindx1(kcol)
        ns = nbasp(nb)

        if (ns.ge.narn1 .and. ns.le.narn2) then
            jfl = jflag(ns)

            if (kcol .eq. kwater) then
                zvclg1(kcol) = xbrwlg
                actlg(narn1) = xbrwlg + acflg(narn1)
            else if (jfl .le. 15) then
                zvclg1(kcol) = conclg(ns)
                actlg(ns) = conclg(ns) + acflg(ns)
            else if (jfl .eq. 16) then
                zvclg1(kcol) = coval(ns) - acflg(ns)
                conclg(ns) = zvclg1(kcol)
                actlg(ns) = conclg(ns) + acflg(ns)
            else if (jfl.eq.19 .or. jfl.eq.20) then
                zvclg1(kcol) = -coval(ns) - acflg(ns)
                conclg(ns) = zvclg1(kcol)
                actlg(ns) = conclg(ns) + acflg(ns)
            else if (jfl.eq.22 .or. jfl.eq.23) then
                zvclg1(kcol) = -coval(ns)
                conclg(ns) = zvclg1(kcol)
                actlg(ns) = conclg(ns) + acflg(ns)
            end if
        else if (ns.ge.nern1 .and. ns.le.nern2) then
            zvclg1(kcol) = conclg(ns)
            actlg(ns) = cjbasp(nb)*(xbarlg(ns) + acflg(ns))
        else
            zvclg1(kcol) = tlg(coval(nb)) + xbarlg(ns)
            actlg(ns) = xbarlg(ns) + acflg(ns)
        end if
    end do

    if (irdxc3 .eq. 0) then
        if (ko2gaq .gt. 0) then
            zvclg1(ko2gaq) = fo2lg
            actlg(no2gaq) = fo2lg
        end if

        if (kelect .gt. 0) then
            zvclg1(kelect) = -pe
            actlg(nelect) = -pe
        end if
    end if

    ! Make starting estimates that must be evaluated simultaneously.
    ! These include all cases of equilibrium constraints and compount
    ! activity constraints (e.g. pHCl) and any case in which log fO2
    ! log fO2 is constrained by Eh, pe-, or a redox couple.
    call arrsim(aamatr,acflg,actlg,bbig,cdrs,cjbasp,cnufac,conc,conclg,coval,delvec,dlogxw,eh,ehfac,eps100,gmmatr,iction,iindx1,iodb,ipivot,irdxc3,ixbasp,jcsort,jflag,jjndex,kbt,ker,khydr,kkndex,kmax,kwater,narn1,narn2,nbasp,nbt,nbti,nbtmax,nbw,ncosp,ndecsp,ndrs,ndrsmx,ndrsr,nelect,nern1,nern2,nhydr,nodbmx,no2gaq,noutpt,npass,nredox,nstmax,nttyo,omega,qawfix,rhsvec,ucospi,uspec,xbar,xbarlg,xbarw,xbrwlg,xlke,xlks,zchar,zvclg1)

    ! Recalculate the concentrations, etc., of dependent species.
    call ncmpex(acflg,act,actlg,cdrs,cegexs,cgexj,conc,conclg,cpgexs,egexjc,egexjf,egexs,eps100,fo2,fo2lg,fsort,fugac,fugalg,iern1,iern2,ietmax,ifrn1,ifrn2,igas,igstak,iindx1,ilrn1,ilrn2,imrn1,imrn2,istack,ixrn1,ixrn2,jcsort,jern1,jern2,jetmax,jflag,jgext,jgsort,jgstak,jjsort,jpflag,jsflag,jsitex,jssort,jstack,kbt,kdim,kelect,kmax,km1,ko2gaq,kwater,kxt,loph,losp,lsort,mgext,mrgexs,mtb,moph,mosp,narn1,narn2,nbasp,nbt,nbtmax,ncmpr,ndrs,ndrsmx,ndrsr,nelect,nern1,nern2,netmax,ngexsa,ngext,ngrn1,ngrn2,ngt,ngtmax,noutpt,no2gaq,nphasx,npt,nptmax,nst,nstmax,nttyo,omega,omeglg,press,qxbarw,q6mode,ugexj,ugexmo,uphase,uspec,xbar,xbarlg,xbarw,xbarwc,xbrwlc,xbrwlg,xlks,zchar,zgexj,zvclg1,zvec1)

    xbarw = xbar(narn1)
    xbrwlg = xbarlg(narn1)

    ! Compute the residuals.
    call betas(acflg,actlg,afcnst,alpha,amtb,bbig,beta,betamx,bneg,cdrs,conc,conclg,coval,csts,eh,ehfac,fo2lg,ibetmx,iebal,iindx1,irdxc3,jcsort,jflag,jsflag,jssort,kbt,kdim,kelect,khydr,kmax,km1,ko2gaq,kwater,kxt,mtb,mosp,narn1,narn2,nbasp,nbtmax,ncosp,ndrs,ndrsmx,ndrsr,nelect,nern1,nern2,nhydr,noutpt,no2gaq,nredox,nst,nstmax,nsts,nstsmx,nstsr,ntfx,ntfxmx,ntfxt,nttyo,omega,qredox,q6mode,tfx,ubbig,ubneg,ubetmx,uspec,uzvec1,weight,xbrwlg,xlke,xlks,zchar)

    ! Calculate the beta convergence function.
    betfnc = 0.

    if (mod(ncycle,2) .eq. 0) then
        if (btmxoe .ge. smp100) then
            betfnc = (btmxoe - betamx)/btmxoe
        end if

        btmxoe = betamx
    else
        if (btmxoo .ge. smp100) then
            betfnc = (btmxoo - betamx)/btmxoo
        end if

        btmxoo = betamx
    end if

    ! Print values of master iteration variables.
    if (iodb(3) .ge. 2) then
        write (noutpt,1230)
1230 format(//10x,'--- Pre-Newton-Raphson Optimization Summary ---',//2x,'kcol   Name',32x,'zvclg1      zvec1',/)

        do kcol = 1,kdim
            nb = iindx1(kcol)
            ns = nbasp(nb)
            zx1 = zvclg1(kcol)
            zx2 = texp(zx1)

            ! Calling sequence substitutions:
            !   uspec(ns) for unam48
            call fmspnx(jlen,uspec(ns),uspn56)
            write (noutpt,1240) kcol,uspn56,zx1,zx2
1240 format(1x,i4,2x,a32,2x,f10.4,2x,1pe12.5)
        end do

        write (noutpt,1250)
1250 format(/2x,'krow   Name',32x,'Beta',/)

        do krow = 1,kdim
            nb = iindx1(krow)
            ns = nbaspd(nb)

            ! Calling sequence substitutions:
            !   uspec(ns) for unam48
            call fmspnx(jlen,uspec(ns),uspn56)
            write (noutpt,1260) krow,uspn56,beta(krow)
1260 format(1x,i4,2x,a32,2x,1pe12.5)
        end do

        write (noutpt,1030)
    end if

    ! Identify the dominant species in each mass balance and
    ! compute the corresponding exponent for a continued
    ! fraction correction.
    call cfracf(cdrs,csts,efac,jcsort,jflag,jssort,kmax,mosp,narn1,narn2,nbasp,nbaspd,nbt,nbtmax,ndrs,ndrsmx,ndrsr,nern1,nern2,nfac,nst,nstmax,nsts,nstsmx,nstsr,q6mode,weight)

    ! Filter the data obtained from EQLIB/cfracf. The following
    ! filters are specific to EQ3NR.
    if (iebal .gt. 0) then
        ! For this species, a charge balance constraint is used instead
        ! of a mass balance constraint.
        nfac(iebal) = 0
        efac(iebal) = 1.0
    end if

    if (nbw .gt. 0) then
        ! For water, a mole fraction equation is used instead of a mass
        ! balance constraint.
        nfac(nbw) = 0
        efac(nbw) = 1.0
    end if

    do nb = 1,nbt
        ns = nbasp(nb)

        if (ns.eq.nhydr .or. ns.eq.nhydx) then
            ! For H+ or OH-, mass balance constraints are rarely used.
            if (jflag(ns) .gt. 15) then
                ! In the present case, a mass balance constraint is not
                ! being used.
                nfac(nb) = 0
                efac(nb) = 1.0
            end if
        else if (ns.eq.no2gaq .or. ns.eq.nelect) then
            ! For O2(g,aq) or e-, mass balance contraints are not used.
            nfac(nb) = 0
            efac(nb) = 1.0
        end if
    end do

    if (iodb(3) .ge. 3) then
        ! Write a table containing the preliminary results.
        kount = 0

        do kcol = 1,kbt
            nb = iindx1(kcol)
            ns = nbasp(nb)
            ns2 = nfac(nb)

            if (ns2.ne.0 .and. ns2.ne.ns) then
                kount = kount + 1
            end if
        end do

        if (kount .gt. 0) then
            write (noutpt,1300)
1300 format(16x,'--- Mass Balance Dominants ---',/)

            write (noutpt,1310)
1310 format(4x,'Master Species',21x,'Dominant Species',/)

            do kcol = 1,kbt
                nb = iindx1(kcol)
                ns = nbasp(nb)
                ns2 = nfac(nb)

                if (ns2.ne.0 .and. ns2.ne.ns) then
                    ! Calling sequence substitutions:
                    !   jlen1 for jlen
                    !   uspec(ns) for unam48
                    !   usp156 for uspn56
                    call fmspnx(jlen1,uspec(ns),usp156)

                    ! Calling sequence substitutions:
                    !   jlen2 for jlen
                    !   uspec(ns2) for unam48
                    !   usp256 for uspn56
                    call fmspnx(jlen2,uspec(ns2),usp256)
                    jlen2 = min(jlen2,32)
                    write (noutpt,1320) usp156,usp256(1:jlen2)
1320 format(2x,a32,3x,a)
                end if
            end do

            write (noutpt,1030)
        end if
    end if

    ! Set up the bfac correction factor array. In the continued fraction
    ! method, m(new) = m(old)/bfac, where bfac = (beta + 1)**efac.
    ! If the same species dominates more than one mass balance,
    ! then this algorithm can be applied to only one of the associated
    ! basis species. Otherwise, oscillatory behavior will occur. In each
    ! set of mass balances with a common dominating species, find the
    ! mass balance with the greatest bfac factor. This is usually nearly
    ! equivalent to finding the mass balance with the greater beta
    ! residual, as efac often has a value of unity.
    call gbfac(beta,bfac,efac,iindx1,kbt,kmax,nbt,nbtmax,nfac)

    if (iodb(3).ge.3 .and. .not.qbswx) then
        ! Write a table containing the modified results.
        write (noutpt,1330)
1330 format(/16x,'--- Factors for Continued Fraction',' Corrections ---',/)

        write (noutpt,1340)
1340 format(4x,'Master Species',22x,'bfac',10x,'efac',/)

        do nb = 1,nbt
            if (bfac(nb) .gt. 0.) then
                ns = nbasp(nb)

                ! Calling sequence substitutions:
                !   uspec(ns) for unam48
                call fmspnx(jlen,uspec(ns),uspn56)
                write (noutpt,1350) uspn56,bfac(nb),efac(nb)
1350 format(2x,a32,3x,1pe12.5,3x,1pe12.5)
            end if
        end do

        write (noutpt,1030)
    end if

    if (qabsw .and. qloop .and. nloop.lt.nlopmx) then
        ! In automatic basis switching mode (iopt(11) .ge. 1), try to
        ! first reduce the magntiude of large positive mass balance
        ! residuals by making one or more basis switches.
        call absswa(adhfs,adhfsx,advfs,advfsx,avcnst,axhfs,axhfsx,axlks,axlksx,axvfs,axvfsx,beta,cdrs,cdrsx,cdrtw,cdrw,csts,dhfs,dvfs,efac,eps100,ibswx,iebal,iindx1,iodb,ipch,ipchmx,ipcv,ipcvmx,jcsort,jflag,jsflag,jssort,kbt,kmax,mosp,narn1,narn2,narxmx,narxt,nbasp,nbaspd,nbaspx,nbt,nbtmax,nbw,ncosp,ndrs,ndrsmx,ndrsr,ndrsrx,ndrsx,nelect,nhydr,nodbmx,no2gaq,noutpt,nst,nstmax,nsts,nstsmx,nstsr,nswtch,ntpr,ntprmx,nttyo,presg,press,qbassw,qbswx,q6mode,tempc,uspec,uzvec1,weight,xvfs,xlks,xhfs)

        if (nswtch .le. 0) then
            ! No switches were made.
            qloop = .false.
            go to 250
        end if

        ! Reset the ixbasp and cjbasp arrays. The former is a flag
        ! array, each member of which denotes whether the
        ! thermodynamic activity of the corresponding basis species
        ! is defined in terms of molality (= 0) or mole fraction (= 1).
        ! The cjbasp array contains any site stoichiometric factors
        ! associated with the operational basis species.
        call gibasp(cgexj,cjbasp,iern1,ixbasp,jern1,jern2,jetmax,jgext,narn1,narn2,nbasp,nbt,nbtmax,nern1,nern2,netmax,nphasx,nstmax)

        ! Null some arrays.
        do ns = 1,nstmax
            conc(ns) = 0.
        end do

        av = -99999.
        call initav(conclg,nstmax,av)

        write (noutpt,1370) nloop,nswtch
        write (nttyo,1370) nloop,nswtch
1370 format(8x,'Completed loop ',i3,' after making ',i3,' basis switches.')

        ! Go back for another loop.
        go to 200
    end if

250 continue

    ! Calculate the electrical balance residual.
    bxecor = 0.
    xecorr = 0.

    if (kebal .gt. 0) then
        ns = nbaspd(iebal)
        cecorr = -alpha(kebal)/zchar(ns)
        xecorr = zchsq2(ns)*cecorr
        call gszm(conc,jcsort,narn1,narn2,nstmax,sigza,sigzc,sigzi,sigzm,zchar)
        bxecor = abs(xecorr)/sigzm
    end if

    if (iodb(3) .ge. 1) then
        ! Calling sequence substitutions:
        !   jlen1 for jlen
        !   ubbig for unam48
        !   usp156 for uspn56
        call fmspnx(jlen1,ubbig,usp156)

        ! Calling sequence substitutions:
        !   jlen2 for jlen
        !   ubneg for unam48
        !   usp256 for uspn56
        call fmspnx(jlen2,ubneg,usp256)

        write (noutpt,1500) betamx,betfnc,bbig,usp156(1:jlen1),bneg,usp256(1:jlen2)
1500 format(18x,'betamx= ',1pe12.5,', betfnc= ',1pe12.5,/18x,'  bbig= ',1pe12.5,', ubbig= ',a,/18x,'  bneg= ',1pe12.5,', ubneg= ',a,/)
    end if

    ! Test the balance residuals for mass and alkalinity to see if
    ! another cycle should be made before attempting to make an improved
    ! estimate of the ionic strength.
    qtestc = bbig.le.tolbig .and. bneg.ge.tolneg .and.  betamx.le.tolxpt .and. bxecor.le.tolzpt

    ! Quit doing cycles if:
    !   1. The cycle convergence criteria are met.
    !   2. The maximum number of cycles have been done.
    !   3. The convergence function betfnc indicates that
    !      the cycles are not converging.
    if (qtestc) then
        go to 300
    end if

    if (ncycle .ge. ncylim) then
        go to 300
    end if

    if (ncycle.gt.2 .and. betfnc.le.tolbtf) then
        negbfc = negbfc + 1

        if (negbfc .ge. 3) then
            go to 300
        end if
    else
        negbfc = 0
    end if

    ! Make improvements in concentration estimates and go back for
    ! another cycle. The algorithm employed here is the modified
    ! continued fraction method. Note that some constraints are placed
    ! on the use of the bfac correction factors. One is to avoid a
    ! blow-out in taking the base ten logarithm of such a factor, which
    ! would occur if the factor had a zero or negative value. Another
    ! limits the magnitude of the change in one step.
    do kcol = 1,kbt
        nb = iindx1(kcol)
        ns = nbasp(nb)

        if (nb.ne.iebal .and. ns.ne.narn1) then
            jfl = jflag(ns)
            bxp1 = beta(kcol) + 1.

            if (jfl.ge.0 .and. jfl.le.3) then
                dx = bfac(nb)
                lx = tlg(dx)

                if (lx .gt. 20.) then
                    lx = 20.
                end if

                conclg(ns) = conclg(ns) - lx
            else if (jfl.ge.7 .and. jfl.le.11) then
                if (bxp1 .le. 0.) then
                    bxp1 = 1.e-20
                end if

                lx = tlg(bxp1)

                if (lx .gt. 20.) then
                    lx = 20.
                end if

                conclg(ns) = conclg(ns) - lx
            end if
        end if
    end do

    if (iebal .gt. 0) then
        ! Electrical balance correction. Skip if mass balance
        ! residuals are way off.
        if (bbig.le.0.25 .and. bneg.ge.-0.25) then
            ns = nbasp(iebal)
            cxo = conclg(ns)
            cxn = conc(ns) + cecorr

            if (cecorr .ge. 0.) then
                cx = tlg(cxn)
                cxl = cxo + 2.0
                conclg(ns) = min(cx,cxl)
            else
                if (cxn .lt. 0.) then
                    cxn = 0.
                end if

                cx = tlg(cxn)
                cxl = cxo - 2.0
                conclg(ns) = max(cx,cxl)
            end if
        end if
    end if

    ! Go back for another cycle.
    go to 220

    ! The cycles for the current pass have been completed. Test to
    ! see if another pass should be made.
300 continue
    write (noutpt,1620) npass,ncycle
    write (nttyo,1620) npass,ncycle
1620 format(13x,'Completed pass ',i3,' in ',i3,' cycles.')

    ! Save the current values of the ionic strength, etc., and of the
    ! the ativity coefficients.
    sigmmo = sigmam
    fxio = fxi
    fjeo = fje

    do ns = narn1,narn2
        acflgo(ns)=acflg(ns)
    end do

    ! Determine the maximum allowed change factors for "Sigma m",
    ! etc. (chfsgm) and the log activity coefficients for aqueous
    ! species (chfacf).
    chfsgm = 1.3

    if (sigmam .le. 2.e-1) then
        chfsgm = 5.
    end if

    if (sigmam .le. 1.e-2) then
        chfsgm = 10.
    end if

    if (sigmam .le. 1.e-3) then
        chfsgm = 100.
    end if

    if (qtestc .and. chfsgm.lt.100.) then
        chfsgm = 100.
    end if

    chfacf = 0.05

    ! Note: setting iter = 0 would fix the activity coefficients
    ! in concentrated solutions. Setting it to a high value (here
    ! 99999) insures that the activity coefficients are recalculated
    ! as desired.
    iter = 99999
    rlxgam = 1.0
    qpracf = iodb(3) .ge. 4

    call ngcadv(abar,acflg,acflgo,actwlc,adh,adhh,adhv,afcnst,al10,aphi,azero,a3bar,a3bars,bacfmx,bdh,bdhh,bdhv,bdot,bdoth,bdotv,bgamx,bpx,bsigmm,bfje,bfxi,cco2,cgexj,chfacf,chfsgm,conc,delam,dgpit,dpelm,dpslm,dselm,elam,eps100,fje,fjeo,fxi,fxio,gpit,ibpxt,ielam,ifcphi1,ifcphi2,ifnnn,ifn2n,ifpsi1,ifpsi2,ifzeta,ilcphi1,ilcphi2,ilnnn,iln2n,ilpsi1,ilpsi2,ilzeta,insgf,iopg,iter,ipndx1,ixrn1,ixrn2,izmax,jcsort,jern1,jern2,jgext,jsol,kx1,kxt,nalpha,napt,narn1,narn2,nchlor,ncmpr,net,nhydr,nmut,nmux,nmxi,nmxx,noutpt,nslt,nslx,nst,nsxi,nsxx,nttyo,omega,palpha,pelm,pmu,press,pslamn,pslm,qhawep,qpit75,qpracf,q6mode,rlxgam,selm,sigmam,sigmmo,tempk,ubacmx,ubgamx,uphase,uspec,wfac,xbar,xbarlg,xbarwc,xbrwlc,zchar,zchcu6,zchsq2)

    ! Recalculate the concentrations, etc., of dependent species.
    call ncmpex(acflg,act,actlg,cdrs,cegexs,cgexj,conc,conclg,cpgexs,egexjc,egexjf,egexs,eps100,fo2,fo2lg,fsort,fugac,fugalg,iern1,iern2,ietmax,ifrn1,ifrn2,igas,igstak,iindx1,ilrn1,ilrn2,imrn1,imrn2,istack,ixrn1,ixrn2,jcsort,jern1,jern2,jetmax,jflag,jgext,jgsort,jgstak,jjsort,jpflag,jsflag,jsitex,jssort,jstack,kbt,kdim,kelect,kmax,km1,ko2gaq,kwater,kxt,loph,losp,lsort,mgext,mrgexs,mtb,moph,mosp,narn1,narn2,nbasp,nbt,nbtmax,ncmpr,ndrs,ndrsmx,ndrsr,nelect,nern1,nern2,netmax,ngexsa,ngext,ngrn1,ngrn2,ngt,ngtmax,noutpt,no2gaq,nphasx,npt,nptmax,nst,nstmax,nttyo,omega,omeglg,press,qxbarw,q6mode,ugexj,ugexmo,uphase,uspec,xbar,xbarlg,xbarw,xbarwc,xbrwlc,xbrwlg,xlks,zchar,zgexj,zvclg1,zvec1)

    xbarw = xbar(narn1)
    xbrwlg = xbarlg(narn1)

    if (iodb(3) .ge. 1) then
        write (noutpt,1510) bsigmm
1510 format(/13x,'bsigmm= ',1pe12.5)

        write (noutpt,1520) bfxi
1520 format(13x,'bfxi= ',1pe12.5)

        write (noutpt,1522) bfje
1522 format(13x,'bfje= ',1pe12.5)

        j2 = ilnobl(ubgamx(1:24))
        write (noutpt,1530) bgamx,ubgamx(1:j2)
1530 format(13x,'bgamx= ',1pe12.5,', ubgamx= ',a)
    end if

    qcsigm =abs(bsigmm) .le. tolgpt
    qcfxi = abs(bfxi) .le. tolgpt
    qcgam = bgamx .le. tolgpt

    qtestp = qcsigm .and. qcfxi .and. qcgam

    ! Are pass criteria satisfied?
    if (.not.qtestp) then
        ! Pass criteria are not satisfied. Test for maximum number
        ! of passes.
        if (npass .ge. nplim) then
            ! Quit. Optimization ended outside requested limits
            ! because the pass requirements were not satisfied.
            if (ker .lt. 2) then
                write (noutpt,1600)
                write (nttyo,1600)
1600 format(/'   Done. Optimization ended outside requested',' limits.',/)
            else
                write (noutpt,1610)
                write (nttyo,1610)
1610 format(/'   Done. Optimization ended outside allowable',' limits.',/)
            end if

            go to 999
        end if

        ! Do another pass.
        go to 210
    end if

    ! Are cycle criteria satisfied?
    if (qtestc) then
        ! Yes, optimization succeeded.
        write (noutpt,1630)
        write (nttyo,1630)
1630 format(/'   Done. Optimization ended within requested',' limits.',/)

        go to 999
    else if (npass .le. 2) then
        ! The pass convergence criteria are satisfied, but the cycle
        ! convergence criteria are not. Try another pass.
        go to 210
    else
        ! Quit. Optimization ended outside requested limits
        ! because cycle requirements were not met.
        if (ker .lt. 2) then
            write (noutpt,1600)
            write (nttyo,1600)
        else
            write (noutpt,1610)
            write (nttyo,1610)
        end if

        go to 999
    end if

999 continue
end subroutine arrset