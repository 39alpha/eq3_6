subroutine path(aadh,aadhh,aadhv,aaphi,abdh,abdhh,abdhv,abdot,abdoth,abdotv,adadhh,adadhv,adbdhh,adbdhv,adbdth,adbdtv,adh,adhh,adhv,afcnst,aftarg,al10,aphi,apx,atwt,avcnst,awmax,awmaxi,awmin,awmini,azero,bdh,bdhh,bdhv,bdot,bdoth,bdotv,bpx,cbsr,cbsri,cco2,cdac,cdrs,cdrsd,cdrsx,cegexs,cesr,cesri,cess,cpgexs,cscale,csigma,csts,dadhh,dadhv,dbdhh,dbdhv,dbdth,dbdtv,dlaplo,dlaprn,dleplo,dleprn,dlhplo,dlhprn,dloplo,dloprn,dltplo,dltpll,dltprl,dltprn,dlxdmp,dlxmax,dlxmin,dlxmx0,dlxplo,dlxpll,dlxprl,dlxprn,eact,egers,egersi,egexjf,ehfac,ehmax,ehmaxi,ehmin,ehmini,elecsr,electr,eps100,farad,fkrc,hact,iact,iapxt,iaqsln,ibpxt,ibsrti,ielam,iern1,iern2,iesrti,ifcphi1,ifcphi2,ifnnn,ifn2n,ifpsi1,ifpsi2,ifrn1,ifrn2,ifzeta,igas,igerti,iindx1,ilcphi1,ilcphi2,ilnnn,iln2n,ilpsi1,ilpsi2,ilrn1,ilrn2,ilzeta,imech,imrn1,imrn2,insgf,iodb,iopg,iopr,iopt,ipch,ipndx1,ipcv,irang,itermx,ixrn1,ixrn2,ixrti,izmax,jcode,jffg,jflag,jflagd,jflgi,jgerti,jpflag,jpress,jptffl,jreac,jsflag,jsitex,jsol,jtemp,kbt,kct,kdim,kelect,khydr,khydx,km1,kmt,ko2gaq,kprs,krdxsp,ksplmx,ksppmx,kstpmx,kwater,kxmod,kx1,kxt,loph,losp,modr,moffg,moph,morr,mosp,mprph,mprphi,mprsp,mprspi,mrgers,mrgexs,mtb,mtbi,mtbaq,mtbaqi,mte,mteaq,mwtrc,mwtsp,narn1,narn2,narxt,nat,nbasp,nbaspd,nbaspi,nbaspx,nbkupa,nbkupb,nbt,nbtd,nbti,nbw,nchlor,ncmpr,nct,ndac,ndact,nelect,nern1,nern2,ness,nessr,net,ndrs,ndrsd,ndrsx,ndrsr,ndrsrd,ndrsrx,nert,newin,nffg,nfrn1,nfrn2,ngrn1,ngrn2,ngt,nhydr,nhydx,nllnmx,nlrn1,nlrn2,nlt,nmrn1,nmrn2,nmrt,nmt,nobswt,noutpt,no2gaq,npchk,nphasx,nprob,nprpti,nprsti,npslmx,npt,nrct,nrdxsp,nrk,nrndex,nsbswt,nsk,nsrt,nsslmx,nst,nsts,nstsr,ntabx,ntf1,ntf1t,ntf2,ntf2t,ntitl1,ntitl2,ntitld,ntpr,ntprt,ntrymx,nttyo,nxmod,nxopex,nxopt,nxridx,nxrn1,nxrn2,nxrt,nxt,o2max,o2maxi,o2min,o2mini,phmax,phmaxi,phmin,phmini,prcinf,press,pressb,pressd,pressi,ptk,qcnpre,qcntmp,qdwipp,qecon,qgexsh,qhawep,qoptmz,qpit75,qredox,qscon,qtatxt,rconst,rcnstv,rk,rkb,rtcnst,rxbar,rxbari,sfcar,smp100,sscrew,ssfcar,tempc,tempcb,tempcd,tempci,tempcu,tempk,tf1,tf2,timemx,time1,timmxi,tistrt,tistti,trkb,ttk,tolaft,tolbt,toldl,tolsat,tolsst,tolxsf,tolxst,tolxsu,uaqsln,ubmtbi,ubsri,ucxri,udac,uelem,uesri,uffg,ufixf,ugerji,ugermo,ugersi,uinfor,ulinex,uobsw,uphase,uplatm,uprphi,uprspi,ureac,usbsw,uspec,usteq6,utitl1,utitl2,utitld,uveeq6,uxcat,uxmod,uxopex,uxopt,uzvec1,uzveci,vosp0,vreac,wfac,xgers,xgersi,ximax,ximaxi,xistrt,xistti,xi1,xlkffg,xlkmod,zchar,zchcu6,zchsq2,zklgmn,zklogl,zklogu,zvclgi,zvclg1,zvec1)
    !! This subroutine is the control subroutine for tracing a reaction
    !! path. It determines a variable step size, which is in terms of an
    !! overall reaction progress variable. The corresponding time step,
    !! if defined, is computed as a dependent variable. This subroutine
    !! also detects and responds to various events, such as saturation
    !! of irreversible reactions and exhaustion of reactants.
    !! The mathematical treatment of the algebraic relations is
    !! analogous to the predictor-corrector method of integrating
    !! differential equations. Finite difference expressions of order
    !! up to nordmx are used as predictor functions. These are
    !! manipulated into the form of truncated Taylor's series.
    !! EQ6/eqcalc.f is used as the "corrector." This subroutine does
    !! equilbrium calculations using a hybrid Newton-Raphson algorithm.
    !! Note that it corrects to satisfy algebraic equations, not
    !! differential equations.
    !! In kinetic mode (iopt(2) > 0), rate laws are integrated using an
    !! actual predictor-corrector algorithm. Algebraic equations continue
    !! to be dealt with as described above. Rate law integration is
    !! done using finite-difference expressions of the same order
    !! as those used to represent algebraic variables.
    !! Apart from no time/time, there are three distinct calculational
    !! modes:
    !!   (1) Normal mode- The step size is constrained to keep the
    !!     predictor functions fairly accurate. There is less burden
    !!     on the Newton-Raphson thermodynamic calculations.
    !!     Normal mode requires longer run times and is more
    !!     expensive than economy mode. Normal mode must be used
    !!     for some kinds of problems. Examples include the cases
    !!     iopt(2) > 0 (kinetic mode) and iopt(1) = 2, (flow-through
    !!     open system mode). Setting iopt(13) = 0 insures that
    !!     normal mode will be selected. It will automatically be
    !!     selected if economy or super-economy modes are not valid
    !!     options.
    !!   (2) Economy mode- The step size is allowed to become large
    !!     more quickly than in normal mode. The predictor functions
    !!     are limited to second order. There is no attempt to
    !!     constrain the step size in order to keep these functions
    !!     accurate. The burden on the equilibrium calculations is
    !!     heavier than in normal mode, because the initial values of
    !!     the unknowns will generally be farther off the mark.
    !!     Economy mode is intended to be faster and cheaper than
    !!     normal mode. The useful information density along the
    !!     reaction path with economy mode is essentially equivalent to
    !!     that obtained with normal mode. Use iopt(13) = 1 to select
    !!     economy mode.
    !!   (3) Super economy mode- This is a special from of economy
    !!     mode. The step size is typically large. It defaults to
    !!     dlxprn, the linear print interval. Phase boundaries are
    !!     ignored. The order of the finite differences is restricted
    !!     to zero. Super economy mode provides much less information
    !!     density along the reaction path than does economy mode or
    !!     normal mode. Use iopt(13) = 2 to select super economy mode.
    !! The step size is controlled by various constraints. One of these
    !! is the estimated accuracy of the predictor functions. If the
    !! system is rapidly changing, delxi will be small and the run time
    !! per unit reaction progress will be long. This is generally the
    !! case at the start of a reaction path simulation. When the system
    !! is changing slowly with respect to reaction progress, delxi may
    !! become large. Other controls on the step size include:
    !!    (1) The initial step size (for order zero)
    !!    (2) An arbitrary upper limit on step size
    !!    (3) The print interval requirements
    !!    (4) The plot interval requirements
    !!    (5) An arbitrary dump interval (flow-through model only)
    !!    (6) The exhaustion of any reactant
    !!    (7) The appearance or disappearance of any secondary
    !!        phase (unless iopt(3) .ge. 1)
    !!    (8) A maximum in the mass of any secondary mineral
    !!        (fluid-centered flow-through model only)
    !!    (9) Crossing 100 degrees C, if temperature is changing with
    !!        reaction progress
    !!   (10) The terminal value of reaction progress
    !!   (11) The terminal value of time, or infinite time
    !! In the flow-through open system model (iopt(1) = 2), masses of
    !! secondary minerals are shifted to the physically removed
    !! subsystem after any of (5) through (11) above.
    !! This subroutine is called by:
    !!   EQ6/eq6.f
    !! Principal input:
    !! Principal output:
    !! Modules.
    !! The module mod6pt contains data required to evaluate Pitzer's
    !! equations.
    use mod6pt

    ! The module mod6xf contains most of the standard-state
    ! thermodynamic data.
    use mod6xf

    implicit none

    include 'eqlib/eqlpar.h'
    include 'eqlib/eqldv.h'
    include 'eqlib/eqlge.h'
    include 'eqlib/eql1s.h'
    include 'eqlib/eqlwd.h'

    ! Calling sequence variable declarations.
    integer :: nbkupa
    integer :: nbkupb
    integer :: newin
    integer :: noutpt
    integer :: ntabx
    integer :: nttyo

    integer :: jcode(nrctmx)
    integer :: jreac(nrctmx)
    integer :: nrndex(nrctmx)
    integer :: nxridx(nrctmx)

    integer :: iact(imchmx,2,nrctmx)
    integer :: imech(2,nrctmx)
    integer :: ndac(ndctmx,imchmx,2,nrctmx)
    integer :: ndact(imchmx,2,nrctmx)
    integer :: nrk(2,nrctmx)
    integer :: nsk(nrctmx)

    integer :: ibsrti(nsrtmx)
    integer :: iesrti(nsrtmx)
    integer :: igerti(jetmax,nertmx)
    integer :: ixrti(nxrtmx)
    integer :: jflgi(nbtmax)
    integer :: jgerti(nertmx)

    integer :: iapxt(nxtmax)
    integer :: ibpxt(nxtmax)
    integer :: iindx1(kmax)
    integer :: ipndx1(kmax)
    integer :: insgf(natmax)
    integer :: iodb(nodbmx)
    integer :: iopg(nopgmx)
    integer :: iopr(noprmx)
    integer :: iopt(noptmx)
    integer :: jffg(nffgmx)
    integer :: jflag(nstmax)
    integer :: jflagd(nstmax)
    integer :: jpflag(nptmax)
    integer :: jsflag(nstmax)
    integer :: jsitex(nstmax)
    integer :: jsol(nxtmax)
    integer :: kxmod(nxmdmx)

    integer :: narxt(ntprmx)
    integer :: nbasp(nbtmax)
    integer :: nbaspd(nbtmax)
    integer :: nbaspi(nbtmax)
    integer :: nbaspx(nbtmax)
    integer :: ncmpr(2,nptmax)
    integer :: ness(nessmx)
    integer :: nessr(2,nstmax)
    integer :: ndrs(ndrsmx)
    integer :: ndrsd(ndrsmx)
    integer :: ndrsx(ndrsmx)
    integer :: ndrsr(2,nstmax)
    integer :: ndrsrd(2,nstmax)
    integer :: ndrsrx(2,nstmax)
    integer :: npchk(nptmax)
    integer :: nphasx(nstmax)
    integer :: nsts(nstsmx)
    integer :: nstsr(2,nstmax)
    integer :: ntf1(ntf1mx)
    integer :: ntf2(ntf2mx)

    integer :: iern1
    integer :: iern2
    integer :: ifrn1
    integer :: ifrn2
    integer :: ilrn1
    integer :: ilrn2
    integer :: imrn1
    integer :: imrn2
    integer :: ixrn1
    integer :: ixrn2

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

    integer :: nat
    integer :: nbt
    integer :: nct
    integer :: net
    integer :: ngt
    integer :: nlt
    integer :: nmt
    integer :: npt
    integer :: nst
    integer :: nxt

    integer :: narn1
    integer :: narn2
    integer :: nern1
    integer :: nern2
    integer :: nfrn1
    integer :: nfrn2
    integer :: ngrn1
    integer :: ngrn2
    integer :: nlrn1
    integer :: nlrn2
    integer :: nmrn1
    integer :: nmrn2
    integer :: nxrn1
    integer :: nxrn2

    integer :: kprs
    integer :: nbti
    integer :: nprpti
    integer :: nprsti

    integer :: iaqsln
    integer :: ielam
    integer :: igas
    integer :: ipch
    integer :: ipcv
    integer :: itermx
    integer :: irang
    integer :: izmax
    integer :: jpress
    integer :: jptffl
    integer :: jtemp
    integer :: kbt
    integer :: kct
    integer :: kdim
    integer :: kelect
    integer :: khydr
    integer :: khydx
    integer :: km1
    integer :: kmt
    integer :: ko2gaq
    integer :: krdxsp
    integer :: ksplmx
    integer :: ksppmx
    integer :: kstpmx
    integer :: kwater
    integer :: kx1
    integer :: kxt
    integer :: nbtd
    integer :: nbw
    integer :: nchlor
    integer :: nelect
    integer :: nert
    integer :: nffg
    integer :: nhydr
    integer :: nhydx
    integer :: nllnmx
    integer :: nmrt
    integer :: nobswt
    integer :: no2gaq
    integer :: nprob
    integer :: npslmx
    integer :: nrct
    integer :: nrdxsp
    integer :: nsbswt
    integer :: nsrt
    integer :: nsslmx
    integer :: ntf1t
    integer :: ntf2t
    integer :: ntitl1
    integer :: ntitl2
    integer :: ntitld
    integer :: ntpr
    integer :: ntprt
    integer :: ntrymx
    integer :: nxmod
    integer :: nxopex
    integer :: nxopt
    integer :: nxrt

    logical :: qcnpre
    logical :: qcntmp
    logical :: qdwipp
    logical :: qecon
    logical :: qgexsh
    logical :: qhawep
    logical :: qoptmz
    logical :: qpit75
    logical :: qredox
    logical :: qscon
    logical :: qstart
    logical :: qtatxt
    logical :: qvhfxi
    logical :: qvlsow

    character(len=nllnmx) :: ulinex
    character(len=80) :: utitl1(ntitmx)
    character(len=80) :: utitl2(ntitmx)
    character(len=80) :: utitld(ntidmx)
    character(len=48) :: ubmtbi(nbtmax)
    character(len=48) :: uprspi(nprsmx)
    character(len=48) :: uzveci(kmax)
    character(len=48) :: uobsw(2,nbtmax)
    character(len=48) :: usbsw(2,nbtmax)
    character(len=48) :: uspec(nstmax)
    character(len=48) :: uxmod(nxmdmx)
    character(len=48) :: uzvec1(kmax)
    character(len=24) :: ubsri(nbt1mx,nsrtmx)
    character(len=24) :: ucxri(iktmax,nxrtmx)
    character(len=24) :: ugersi(ietmax,jetmax,nertmx)
    character(len=24) :: uprphi(nprpmx)
    character(len=24) :: ugermo(nertmx)
    character(len=24) :: ureac(nrctmx)
    character(len=24) :: uffg(nffgmx)
    character(len=24) :: uphase(nptmax)
    character(len=24) :: uxcat(nxopmx)
    character(len=24) :: uxopex(nxpemx)
    character(len=24) :: udac(ndctmx,imchmx,2,nrctmx)
    character(len=24) :: uaqsln
    character(len=8) :: uesri(nctmax,nsrtmx)
    character(len=8) :: ugerji(jetmax,nertmx)
    character(len=8) :: uelem(nctmax)
    character(len=8) :: uxopt(nxopmx)
    character(len=8) :: ufixf
    character(len=8) :: uinfor
    character(len=8) :: uplatm
    character(len=8) :: usteq6
    character(len=8) :: uveeq6

    real(kind=8) :: aadh(narxmx,ntprmx)
    real(kind=8) :: aadhh(narxmx,ntprmx)
    real(kind=8) :: aadhv(narxmx,ntprmx)
    real(kind=8) :: aaphi(narxmx,ntprmx)
    real(kind=8) :: abdh(narxmx,ntprmx)
    real(kind=8) :: abdhh(narxmx,ntprmx)
    real(kind=8) :: abdhv(narxmx,ntprmx)
    real(kind=8) :: abdot(narxmx,ntprmx)
    real(kind=8) :: abdoth(narxmx,ntprmx)
    real(kind=8) :: abdotv(narxmx,ntprmx)

    real(kind=8) :: dadhh(ipchmx)
    real(kind=8) :: dadhv(ipcvmx)
    real(kind=8) :: dbdhh(ipchmx)
    real(kind=8) :: dbdhv(ipcvmx)
    real(kind=8) :: dbdth(ipchmx)
    real(kind=8) :: dbdtv(ipcvmx)

    real(kind=8) :: cbsri(nbt1mx,nsrtmx)
    real(kind=8) :: cesri(nctmax,nsrtmx)
    real(kind=8) :: egersi(ietmax,jetmax,nertmx)
    real(kind=8) :: mprphi(nprpmx)
    real(kind=8) :: mprspi(nprsmx)
    real(kind=8) :: mtbi(nbtmax)
    real(kind=8) :: mtbaqi(nbtmax)
    real(kind=8) :: rxbari(iktmax,nxrtmx)
    real(kind=8) :: xgersi(ietmax,jetmax,nertmx)
    real(kind=8) :: zvclgi(kmax)

    real(kind=8) :: adadhh(narxmx,ntprmx,ipchmx)
    real(kind=8) :: adadhv(narxmx,ntprmx,ipcvmx)
    real(kind=8) :: adbdhh(narxmx,ntprmx,ipchmx)
    real(kind=8) :: adbdhv(narxmx,ntprmx,ipcvmx)
    real(kind=8) :: adbdth(narxmx,ntprmx,ipchmx)
    real(kind=8) :: adbdtv(narxmx,ntprmx,ipcvmx)

    real(kind=8) :: cdac(ndctmx,imchmx,2,nrctmx)
    real(kind=8) :: csigma(imchmx,2,nrctmx)
    real(kind=8) :: eact(imchmx,2,nrctmx)
    real(kind=8) :: fkrc(nrctmx)
    real(kind=8) :: hact(imchmx,2,nrctmx)
    real(kind=8) :: rkb(imchmx,2,nrctmx)
    real(kind=8) :: trkb(imchmx,2,nrctmx)

    real(kind=8) :: cbsr(nbt1mx,nsrtmx)
    real(kind=8) :: cesr(nctmax,nsrtmx)
    real(kind=8) :: egers(ietmax,jetmax,nertmx)
    real(kind=8) :: elecsr(nsrtmx)
    real(kind=8) :: modr(nrctmx)
    real(kind=8) :: mrgers(ietmax,jetmax,nertmx)
    real(kind=8) :: morr(nrctmx)
    real(kind=8) :: mwtrc(nrctmx)
    real(kind=8) :: rxbar(iktmax,nxrtmx)
    real(kind=8) :: sfcar(nrctmx)
    real(kind=8) :: ssfcar(nrctmx)
    real(kind=8) :: vreac(nrctmx)
    real(kind=8) :: xgers(ietmax,jetmax,nertmx)

    real(kind=8) :: apx(iapxmx,nxtmax)
    real(kind=8) :: atwt(nctmax)
    real(kind=8) :: azero(natmax)
    real(kind=8) :: bpx(ibpxmx,nxtmax)
    real(kind=8) :: cco2(5)
    real(kind=8) :: cegexs(ietmax,jetmax,netmax)
    real(kind=8) :: cpgexs(ietmax,jetmax,netmax)
    real(kind=8) :: cess(nessmx)
    real(kind=8) :: cdrs(ndrsmx)
    real(kind=8) :: cdrsd(ndrsmx)
    real(kind=8) :: cdrsx(ndrsmx)
    real(kind=8) :: cscale(nstmax)
    real(kind=8) :: csts(nstsmx)
    real(kind=8) :: egexjf(jetmax,netmax)
    real(kind=8) :: loph(nptmax)
    real(kind=8) :: losp(nstmax)
    real(kind=8) :: moph(nptmax)
    real(kind=8) :: mosp(nstmax)
    real(kind=8) :: moffg(nffgmx)
    real(kind=8) :: mprph(nptmax)
    real(kind=8) :: mprsp(nstmax)
    real(kind=8) :: mrgexs(ietmax,jetmax,netmax)
    real(kind=8) :: mtb(nbtmax)
    real(kind=8) :: mtbaq(nbtmax)
    real(kind=8) :: mte(nctmax)
    real(kind=8) :: mteaq(nctmax)
    real(kind=8) :: mwtsp(nstmax)
    real(kind=8) :: ptk(nptkmx)
    real(kind=8) :: rk(imchmx,2,nrctmx)
    real(kind=8) :: sscrew(nsscmx)
    real(kind=8) :: tempcu(ntprmx)
    real(kind=8) :: tf1(ntf1mx)
    real(kind=8) :: tf2(ntf2mx)
    real(kind=8) :: ttk(nttkmx)
    real(kind=8) :: vosp0(nstmax)
    real(kind=8) :: wfac(iktmax,nxtmax)
    real(kind=8) :: xlkffg(nffgmx)
    real(kind=8) :: xlkmod(nxmdmx)
    real(kind=8) :: zchar(nstmax)
    real(kind=8) :: zchcu6(nstmax)
    real(kind=8) :: zchsq2(nstmax)
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

    real(kind=8) :: afcnst
    real(kind=8) :: aftarg
    real(kind=8) :: al10
    real(kind=8) :: avcnst
    real(kind=8) :: awmax
    real(kind=8) :: awmaxi
    real(kind=8) :: awmin
    real(kind=8) :: awmini
    real(kind=8) :: dlaplo
    real(kind=8) :: dlaprn
    real(kind=8) :: dleplo
    real(kind=8) :: dleprn
    real(kind=8) :: dlhplo
    real(kind=8) :: dlhprn
    real(kind=8) :: dloplo
    real(kind=8) :: dloprn
    real(kind=8) :: dltpll
    real(kind=8) :: dltplo
    real(kind=8) :: dltprl
    real(kind=8) :: dltprn
    real(kind=8) :: dlxdmp
    real(kind=8) :: dlxmax
    real(kind=8) :: dlxmin
    real(kind=8) :: dlxmx0
    real(kind=8) :: dlxpll
    real(kind=8) :: dlxplo
    real(kind=8) :: dlxprl
    real(kind=8) :: dlxprn
    real(kind=8) :: ehfac
    real(kind=8) :: ehmax
    real(kind=8) :: ehmaxi
    real(kind=8) :: ehmin
    real(kind=8) :: ehmini
    real(kind=8) :: electr
    real(kind=8) :: eps100
    real(kind=8) :: farad
    real(kind=8) :: o2max
    real(kind=8) :: o2maxi
    real(kind=8) :: o2min
    real(kind=8) :: o2mini
    real(kind=8) :: phmax
    real(kind=8) :: phmaxi
    real(kind=8) :: phmin
    real(kind=8) :: phmini
    real(kind=8) :: prcinf
    real(kind=8) :: press
    real(kind=8) :: pressb
    real(kind=8) :: pressd
    real(kind=8) :: pressi
    real(kind=8) :: rcnstv
    real(kind=8) :: rconst
    real(kind=8) :: rtcnst
    real(kind=8) :: smp100
    real(kind=8) :: tempc
    real(kind=8) :: tempcb
    real(kind=8) :: tempcd
    real(kind=8) :: tempci
    real(kind=8) :: tempk
    real(kind=8) :: timemx
    real(kind=8) :: time1
    real(kind=8) :: timmxi
    real(kind=8) :: tistrt
    real(kind=8) :: tistti
    real(kind=8) :: tolaft
    real(kind=8) :: tolbt
    real(kind=8) :: toldl
    real(kind=8) :: tolsat
    real(kind=8) :: tolsst
    real(kind=8) :: tolxsf
    real(kind=8) :: tolxst
    real(kind=8) :: tolxsu
    real(kind=8) :: ximax
    real(kind=8) :: ximaxi
    real(kind=8) :: xistrt
    real(kind=8) :: xistti
    real(kind=8) :: xi1
    real(kind=8) :: zklgmn
    real(kind=8) :: zklogl
    real(kind=8) :: zklogu

    ! Local variable declarations with static global dimensioning.
    character(len=32) :: uxtype(jso_par)

    real(kind=8) :: egexjc(jet_par,net_par)
    real(kind=8) :: egexpa(net_par)
    real(kind=8) :: egexpc(net_par)
    real(kind=8) :: egexs(iet_par,jet_par,net_par)
    real(kind=8) :: egexw(ket_par,net_par)
    real(kind=8) :: xgexw(ket_par,net_par)

    ! Local variable declarations with variable global dimensioning.
    integer, dimension(:), allocatable :: ibswx
    integer, dimension(:), allocatable :: idirec

    integer, dimension(:), allocatable :: igstak
    integer, dimension(:), allocatable :: ipivot
    integer, dimension(:), allocatable :: ipivtr
    integer, dimension(:), allocatable :: istack
    integer, dimension(:), allocatable :: ixbasp
    integer, dimension(:), allocatable :: jcsort
    integer, dimension(:), allocatable :: jgstak
    integer, dimension(:), allocatable :: jgsort
    integer, dimension(:), allocatable :: jjsort
    integer, dimension(:), allocatable :: jssort
    integer, dimension(:), allocatable :: jstack
    integer, dimension(:), allocatable :: kction

    integer, dimension(:), allocatable :: iemop
    integer, dimension(:), allocatable :: iemop0
    integer, dimension(:), allocatable :: iemos
    integer, dimension(:), allocatable :: iemos0
    integer, dimension(:), allocatable :: iexr
    integer, dimension(:), allocatable :: iindx0
    integer, dimension(:), allocatable :: ipndx0
    integer, dimension(:), allocatable :: jexr
    integer, dimension(:), allocatable :: jreac0
    integer, dimension(:), allocatable :: jsca
    integer, dimension(:), allocatable :: jscr

    integer, dimension(:,:), allocatable :: ncmpe
    integer, dimension(:,:), allocatable :: ncmpe0

    integer, dimension(:,:,:), allocatable :: ndactb

    integer, dimension(:,:,:,:), allocatable :: ndacb

    logical, dimension(:), allocatable :: qxknph

    character(len=48), dimension(:), allocatable :: uzvec0
    character(len=8), dimension(:), allocatable :: ulbeta
    character(len=8), dimension(:), allocatable :: uldel

    real(kind=8), dimension(:,:), allocatable :: akmat0
    real(kind=8), dimension(:,:), allocatable :: akmat1
    real(kind=8), dimension(:,:), allocatable :: daffp0
    real(kind=8), dimension(:,:), allocatable :: dafrc0
    real(kind=8), dimension(:,:), allocatable :: demop0
    real(kind=8), dimension(:,:), allocatable :: demos0
    real(kind=8), dimension(:,:), allocatable :: drer0
    real(kind=8), dimension(:,:), allocatable :: drer0s
    real(kind=8), dimension(:,:), allocatable :: dzvc0
    real(kind=8), dimension(:,:), allocatable :: dzvc0s
    real(kind=8), dimension(:,:), allocatable :: fdafm1
    real(kind=8), dimension(:,:), allocatable :: fdaf0
    real(kind=8), dimension(:,:), allocatable :: fdarm1
    real(kind=8), dimension(:,:), allocatable :: fdar0
    real(kind=8), dimension(:,:), allocatable :: fdpem1
    real(kind=8), dimension(:,:), allocatable :: fdpe0
    real(kind=8), dimension(:,:), allocatable :: fdrem1
    real(kind=8), dimension(:,:), allocatable :: fdre0
    real(kind=8), dimension(:,:), allocatable :: fdre1
    real(kind=8), dimension(:,:), allocatable :: fdrrm1
    real(kind=8), dimension(:,:), allocatable :: fdrr0
    real(kind=8), dimension(:,:), allocatable :: fdrr1
    real(kind=8), dimension(:,:), allocatable :: fdsem1
    real(kind=8), dimension(:,:), allocatable :: fdse0
    real(kind=8), dimension(:,:), allocatable :: fdzvm1
    real(kind=8), dimension(:,:), allocatable :: fdzv0

    real(kind=8), dimension(:), allocatable :: affp
    real(kind=8), dimension(:), allocatable :: affp0
    real(kind=8), dimension(:), allocatable :: affs
    real(kind=8), dimension(:), allocatable :: afrc0
    real(kind=8), dimension(:), allocatable :: afrc1
    real(kind=8), dimension(:), allocatable :: afrcp
    real(kind=8), dimension(:), allocatable :: cjbasp
    real(kind=8), dimension(:), allocatable :: cnufac
    real(kind=8), dimension(:), allocatable :: ctb
    real(kind=8), dimension(:), allocatable :: daw0
    real(kind=8), dimension(:), allocatable :: deh0
    real(kind=8), dimension(:), allocatable :: do20
    real(kind=8), dimension(:), allocatable :: dph0
    real(kind=8), dimension(:), allocatable :: drir0
    real(kind=8), dimension(:), allocatable :: drir0s
    real(kind=8), dimension(:), allocatable :: dxsm00
    real(kind=8), dimension(:), allocatable :: dxsm10
    real(kind=8), dimension(:), allocatable :: dxsm11
    real(kind=8), dimension(:), allocatable :: dxval0
    real(kind=8), dimension(:), allocatable :: d1zvc1
    real(kind=8), dimension(:), allocatable :: d1emp1
    real(kind=8), dimension(:), allocatable :: d2emp1
    real(kind=8), dimension(:), allocatable :: ehrc
    real(kind=8), dimension(:), allocatable :: emop
    real(kind=8), dimension(:), allocatable :: emop0
    real(kind=8), dimension(:), allocatable :: emos
    real(kind=8), dimension(:), allocatable :: emos0
    real(kind=8), dimension(:), allocatable :: fdaw0
    real(kind=8), dimension(:), allocatable :: fdawm1
    real(kind=8), dimension(:), allocatable :: fdeh0
    real(kind=8), dimension(:), allocatable :: fdehm1
    real(kind=8), dimension(:), allocatable :: fdo20
    real(kind=8), dimension(:), allocatable :: fdo2m1
    real(kind=8), dimension(:), allocatable :: fdph0
    real(kind=8), dimension(:), allocatable :: fdphm1
    real(kind=8), dimension(:), allocatable :: fdrim1
    real(kind=8), dimension(:), allocatable :: fdri0
    real(kind=8), dimension(:), allocatable :: fdri1
    real(kind=8), dimension(:), allocatable :: fo2lrc

    real(kind=8), dimension(:), allocatable :: acflg
    real(kind=8), dimension(:), allocatable :: acflgo
    real(kind=8), dimension(:), allocatable :: acflg0
    real(kind=8), dimension(:), allocatable :: act
    real(kind=8), dimension(:), allocatable :: actlg
    real(kind=8), dimension(:), allocatable :: ahrc
    real(kind=8), dimension(:), allocatable :: affpd
    real(kind=8), dimension(:), allocatable :: affsd
    real(kind=8), dimension(:), allocatable :: alpha
    real(kind=8), dimension(:), allocatable :: amtb
    real(kind=8), dimension(:), allocatable :: a3bars
    real(kind=8), dimension(:), allocatable :: beta
    real(kind=8), dimension(:), allocatable :: betao
    real(kind=8), dimension(:), allocatable :: cdrtw
    real(kind=8), dimension(:), allocatable :: cdrw
    real(kind=8), dimension(:), allocatable :: conc
    real(kind=8), dimension(:), allocatable :: conclg
    real(kind=8), dimension(:), allocatable :: cteaq
    real(kind=8), dimension(:), allocatable :: delvco
    real(kind=8), dimension(:), allocatable :: delvec
    real(kind=8), dimension(:), allocatable :: dlogxw
    real(kind=8), dimension(:), allocatable :: fsort
    real(kind=8), dimension(:), allocatable :: fugac
    real(kind=8), dimension(:), allocatable :: fugalg
    real(kind=8), dimension(:), allocatable :: lsort
    real(kind=8), dimension(:), allocatable :: ppmwb
    real(kind=8), dimension(:), allocatable :: ppmwe
    real(kind=8), dimension(:), allocatable :: rhsvec
    real(kind=8), dimension(:), allocatable :: rreacn
    real(kind=8), dimension(:), allocatable :: rreac0
    real(kind=8), dimension(:), allocatable :: rreac1
    real(kind=8), dimension(:), allocatable :: rrelrp
    real(kind=8), dimension(:), allocatable :: rrelr0
    real(kind=8), dimension(:), allocatable :: rrelr1
    real(kind=8), dimension(:), allocatable :: sfcar0
    real(kind=8), dimension(:), allocatable :: sidrph
    real(kind=8), dimension(:), allocatable :: sidrsp
    real(kind=8), dimension(:), allocatable :: wodr
    real(kind=8), dimension(:), allocatable :: worr
    real(kind=8), dimension(:), allocatable :: xbar
    real(kind=8), dimension(:), allocatable :: xbarlg

    real(kind=8), dimension(:), allocatable :: modr0
    real(kind=8), dimension(:), allocatable :: mophg
    real(kind=8), dimension(:), allocatable :: mophj
    real(kind=8), dimension(:), allocatable :: mopht
    real(kind=8), dimension(:), allocatable :: moph0
    real(kind=8), dimension(:), allocatable :: morr0
    real(kind=8), dimension(:), allocatable :: mospg
    real(kind=8), dimension(:), allocatable :: mospj
    real(kind=8), dimension(:), allocatable :: mospt
    real(kind=8), dimension(:), allocatable :: mosp0
    real(kind=8), dimension(:), allocatable :: mtb0
    real(kind=8), dimension(:), allocatable :: perc
    real(kind=8), dimension(:), allocatable :: voph
    real(kind=8), dimension(:), allocatable :: vophg
    real(kind=8), dimension(:), allocatable :: vophj
    real(kind=8), dimension(:), allocatable :: vopht
    real(kind=8), dimension(:), allocatable :: vosp
    real(kind=8), dimension(:), allocatable :: vospg
    real(kind=8), dimension(:), allocatable :: vospj
    real(kind=8), dimension(:), allocatable :: vospt
    real(kind=8), dimension(:), allocatable :: weight
    real(kind=8), dimension(:), allocatable :: woph
    real(kind=8), dimension(:), allocatable :: wophg
    real(kind=8), dimension(:), allocatable :: wophj
    real(kind=8), dimension(:), allocatable :: wopht
    real(kind=8), dimension(:), allocatable :: wosp
    real(kind=8), dimension(:), allocatable :: wospg
    real(kind=8), dimension(:), allocatable :: wospj
    real(kind=8), dimension(:), allocatable :: wospt
    real(kind=8), dimension(:), allocatable :: xirct
    real(kind=8), dimension(:), allocatable :: xirct0
    real(kind=8), dimension(:), allocatable :: zvclg0
    real(kind=8), dimension(:), allocatable :: zvec0

    real(kind=8), dimension(:), allocatable :: alphar
    real(kind=8), dimension(:), allocatable :: betar
    real(kind=8), dimension(:), allocatable :: delvcr
    real(kind=8), dimension(:), allocatable :: rhsvcr
    real(kind=8), dimension(:), allocatable :: dvjdte

    real(kind=8), dimension(:), allocatable :: hhcvec
    real(kind=8), dimension(:), allocatable :: xhcvec

    real(kind=8), dimension(:,:), allocatable :: aamatr
    real(kind=8), dimension(:,:), allocatable :: gmmatr
    real(kind=8), dimension(:,:), allocatable :: armatr
    real(kind=8), dimension(:,:), allocatable :: grmatr
    real(kind=8), dimension(:,:), allocatable :: aimatr
    real(kind=8), dimension(:,:), allocatable :: mmmatr
    real(kind=8), dimension(:,:), allocatable :: sgmatr
    real(kind=8), dimension(:,:), allocatable :: xxmatr
    real(kind=8), dimension(:,:), allocatable :: xymatr
    real(kind=8), dimension(:,:), allocatable :: rrxfi1

    real(kind=8), dimension(:,:,:,:), allocatable :: cdacb

    ! Local variable declarations.
    integer :: nbkupn

    integer :: nrct1

    integer :: i
    integer :: iavkdm
    integer :: ibtrmx
    integer :: idlrmx
    integer :: ier
    integer :: iexrt
    integer :: inmax
    integer :: iter
    integer :: iwdh2o
    integer :: iwnffg
    integer :: j
    integer :: jcut
    integer :: jexrt
    integer :: jexrtx
    integer :: jordlm
    integer :: jsawth
    integer :: jscat
    integer :: jscatx
    integer :: jscrt
    integer :: jscrtx
    integer :: j2
    integer :: j3
    integer :: kaft1
    integer :: kcol
    integer :: kdim0
    integer :: kly
    integer :: km10
    integer :: kmt0
    integer :: kord
    integer :: kordlm
    integer :: kordp1
    integer :: kordsv
    integer :: kpsat
    integer :: kpsst
    integer :: kstep
    integer :: kstpab
    integer :: kstppl
    integer :: kstppr
    integer :: kstpmn
    integer :: kstpze
    integer :: kx10
    integer :: kxt0
    integer :: kzmax
    integer :: n
    integer :: naft1
    integer :: nb
    integer :: nc
    integer :: ncorr
    integer :: ncut
    integer :: ndelay
    integer :: ne
    integer :: ng
    integer :: nlwffg
    integer :: nmax
    integer :: nord
    integer :: nordr
    integer :: nordrs
    integer :: nordsv
    integer :: nordz
    integer :: nordzs
    integer :: np
    integer :: npe
    integer :: npet
    integer :: npet0
    integer :: npts
    integer :: nrc
    integer :: nr1
    integer :: nr2
    integer :: ns
    integer :: nsawth
    integer :: nse
    integer :: nset
    integer :: nset0
    integer :: nswtch
    integer :: ntpr0
    integer :: nwdh2o
    integer :: nweope
    integer :: nwndpc
    integer :: nwnffg

    integer :: iarmxn
    integer :: ilnobl

    logical :: qabswx
    logical :: qadjdx
    logical :: qaflip
    logical :: qaft1
    logical :: qbassw
    logical :: qbseqc
    logical :: qbswx
    logical :: qbye
    logical :: qconst
    logical :: qodeok
    logical :: qftpr2
    logical :: qdmpr1
    logical :: qdmpr2
    logical :: qdump
    logical :: qhcon
    logical :: qmin
    logical :: qmod
    logical :: qmod1
    logical :: qmod2
    logical :: qphcl
    logical :: qplaw0
    logical :: qplaw1
    logical :: qpleh0
    logical :: qpleh1
    logical :: qplolt
    logical :: qplolx
    logical :: qplott
    logical :: qplotx
    logical :: qplo20
    logical :: qplo21
    logical :: qplph0
    logical :: qplph1
    logical :: qpr
    logical :: qpraw0
    logical :: qpraw1
    logical :: qpreh0
    logical :: qpreh1
    logical :: qprnlt
    logical :: qprnlx
    logical :: qprntt
    logical :: qprntx
    logical :: qpro20
    logical :: qpro21
    logical :: qprph0
    logical :: qprph1
    logical :: qrapch
    logical :: qreax
    logical :: qreq
    logical :: qrho
    logical :: qriinf
    logical :: qrpcfl
    logical :: qsawth
    logical :: qshoot
    logical :: qskip
    logical :: qstabl
    logical :: qstabr
    logical :: qstabz
    logical :: qstop
    logical :: qstopx
    logical :: qtplo
    logical :: qtprn
    logical :: qtrch
    logical :: qtvchk
    logical :: qwhcfa
    logical :: qx
    logical :: qxbarw
    logical :: qzdump
    logical :: qzplot
    logical :: qzprnt
    logical :: qztayl
    logical :: q6mode

    character(len=48) :: ubacmx
    character(len=48) :: ubgamx
    character(len=24) :: ustr24
    character(len=24) :: ux24
    character(len=16) :: ux16
    character(len=16) :: ux16a
    character(len=16) :: ux16b
    character(len=8) :: upkfor

    real(kind=8) :: abar
    real(kind=8) :: acfw
    real(kind=8) :: acfwlg
    real(kind=8) :: actw
    real(kind=8) :: actwlc
    real(kind=8) :: actwlg
    real(kind=8) :: adel
    real(kind=8) :: adlzlg
    real(kind=8) :: adx
    real(kind=8) :: aft1
    real(kind=8) :: aft0
    real(kind=8) :: aftm1
    real(kind=8) :: ah
    real(kind=8) :: ahmes
    real(kind=8) :: ahnbs
    real(kind=8) :: alk
    real(kind=8) :: alki
    real(kind=8) :: alk1
    real(kind=8) :: alk2
    real(kind=8) :: av
    real(kind=8) :: avkdim
    real(kind=8) :: avdlxi
    real(kind=8) :: awstrt
    real(kind=8) :: aw0
    real(kind=8) :: aw0plo
    real(kind=8) :: aw0prn
    real(kind=8) :: aw1
    real(kind=8) :: aw1plo
    real(kind=8) :: aw1prn
    real(kind=8) :: a3bar
    real(kind=8) :: bacfmx
    real(kind=8) :: bbig
    real(kind=8) :: betamx
    real(kind=8) :: bgamx
    real(kind=8) :: bneg
    real(kind=8) :: deltim
    real(kind=8) :: delxi
    real(kind=8) :: delxia
    real(kind=8) :: dlxilm
    real(kind=8) :: dlximx
    real(kind=8) :: dlxipl
    real(kind=8) :: dlxipr
    real(kind=8) :: dlxisv
    real(kind=8) :: dlxis2
    real(kind=8) :: dlxlim
    real(kind=8) :: dlxode
    real(kind=8) :: dlxtmx
    real(kind=8) :: dlxtpr
    real(kind=8) :: dlxtpl
    real(kind=8) :: dvoso
    real(kind=8) :: dwoso
    real(kind=8) :: dx
    real(kind=8) :: dxdmp
    real(kind=8) :: dxsave
    real(kind=8) :: dx1
    real(kind=8) :: dx2
    real(kind=8) :: eh
    real(kind=8) :: ehmes
    real(kind=8) :: ehnbs
    real(kind=8) :: ehstrt
    real(kind=8) :: eh0
    real(kind=8) :: eh0plo
    real(kind=8) :: eh0prn
    real(kind=8) :: eh1
    real(kind=8) :: eh1plo
    real(kind=8) :: eh1prn
    real(kind=8) :: fdlim
    real(kind=8) :: fdx
    real(kind=8) :: fje
    real(kind=8) :: fjeo
    real(kind=8) :: fje0
    real(kind=8) :: fjest
    real(kind=8) :: fo2
    real(kind=8) :: fo2lg
    real(kind=8) :: fo2lg0
    real(kind=8) :: fo2lg1
    real(kind=8) :: fxi
    real(kind=8) :: fxio
    real(kind=8) :: fxi0
    real(kind=8) :: fxist
    real(kind=8) :: fxprpl
    real(kind=8) :: morrw1
    real(kind=8) :: mx
    real(kind=8) :: lprcin
    real(kind=8) :: lx
    real(kind=8) :: osc
    real(kind=8) :: oscst
    real(kind=8) :: omega
    real(kind=8) :: omeglg
    real(kind=8) :: o2strt
    real(kind=8) :: o20plo
    real(kind=8) :: o20prn
    real(kind=8) :: o21plo
    real(kind=8) :: o21prn

    real(kind=8) :: pch
    real(kind=8) :: pe
    real(kind=8) :: pemes
    real(kind=8) :: penbs
    real(kind=8) :: ph
    real(kind=8) :: phcl
    real(kind=8) :: phmes
    real(kind=8) :: phnbs
    real(kind=8) :: phstrt
    real(kind=8) :: ph0
    real(kind=8) :: ph0plo
    real(kind=8) :: ph0prn
    real(kind=8) :: ph1
    real(kind=8) :: ph1plo
    real(kind=8) :: ph1prn
    real(kind=8) :: prminf
    real(kind=8) :: rho
    real(kind=8) :: rirecp
    real(kind=8) :: rirec0
    real(kind=8) :: rirec1
    real(kind=8) :: rx0
    real(kind=8) :: rx1
    real(kind=8) :: rxx
    real(kind=8) :: scale
    real(kind=8) :: scalim
    real(kind=8) :: scfcr
    real(kind=8) :: scfcrs
    real(kind=8) :: scfcz
    real(kind=8) :: scfczs
    real(kind=8) :: scnsti
    real(kind=8) :: scnstd
    real(kind=8) :: screwd
    real(kind=8) :: sigmam
    real(kind=8) :: sigmmo
    real(kind=8) :: sigmm0
    real(kind=8) :: sigmst
    real(kind=8) :: tempc0
    real(kind=8) :: thours
    real(kind=8) :: time0
    real(kind=8) :: tiplol
    real(kind=8) :: tiplot
    real(kind=8) :: tiplxx
    real(kind=8) :: tiprnl
    real(kind=8) :: tiprnt
    real(kind=8) :: tiprxx
    real(kind=8) :: tistrd
    real(kind=8) :: tistry
    real(kind=8) :: tistsv
    real(kind=8) :: tmins
    real(kind=8) :: tolsar
    real(kind=8) :: tolsrr
    real(kind=8) :: tyears
    real(kind=8) :: vodrt
    real(kind=8) :: tdays
    real(kind=8) :: vosoct
    real(kind=8) :: vosol
    real(kind=8) :: wfh2o
    real(kind=8) :: wftds
    real(kind=8) :: whcfac
    real(kind=8) :: wkgh2o
    real(kind=8) :: wwstrt
    real(kind=8) :: wkgsol
    real(kind=8) :: wkgwi
    real(kind=8) :: wodrt
    real(kind=8) :: woh2o
    real(kind=8) :: worrt
    real(kind=8) :: wosoct
    real(kind=8) :: wosol
    real(kind=8) :: wotds
    real(kind=8) :: wx
    real(kind=8) :: xbarw
    real(kind=8) :: xbarwc
    real(kind=8) :: xbrwlc
    real(kind=8) :: xbrwlg
    real(kind=8) :: xf
    real(kind=8) :: xidump
    real(kind=8) :: xilim
    real(kind=8) :: xim1
    real(kind=8) :: xiplol
    real(kind=8) :: xiplot
    real(kind=8) :: xiprnl
    real(kind=8) :: xiprnt
    real(kind=8) :: xistsv
    real(kind=8) :: xi0
    real(kind=8) :: xlf
    real(kind=8) :: xval0
    real(kind=8) :: xx

    real(kind=8) :: dxe0mx
    real(kind=8) :: dxe1mx
    real(kind=8) :: dxe0pl
    real(kind=8) :: dxe1pl
    real(kind=8) :: dxe0pr
    real(kind=8) :: dxe1pr
    real(kind=8) :: dxh0mx
    real(kind=8) :: dxh1mx
    real(kind=8) :: dxh0pl
    real(kind=8) :: dxh1pl
    real(kind=8) :: dxh0pr
    real(kind=8) :: dxh1pr
    real(kind=8) :: dxo0mx
    real(kind=8) :: dxo1mx
    real(kind=8) :: dxo0pl
    real(kind=8) :: dxo1pl
    real(kind=8) :: dxo0pr
    real(kind=8) :: dxo1pr
    real(kind=8) :: dxw0mx
    real(kind=8) :: dxw1mx
    real(kind=8) :: dxw0pl
    real(kind=8) :: dxw1pl
    real(kind=8) :: dxw0pr
    real(kind=8) :: dxw1pr

    real(kind=8) :: btrfnc
    real(kind=8) :: btrmax
    real(kind=8) :: btrmxo
    real(kind=8) :: dlrfnc
    real(kind=8) :: dlrmax
    real(kind=8) :: dlrmxo

    real(kind=8) :: mlmrra
    real(kind=8) :: mrmlra
    real(kind=8) :: rhoc
    real(kind=8) :: rhowc
    real(kind=8) :: tdsgks
    real(kind=8) :: tdsglw
    real(kind=8) :: tdspkc
    real(kind=8) :: tdsplc

    real(kind=8) :: fctrl
    real(kind=8) :: texp
    real(kind=8) :: tlg

    ! Variable declarations: Local data needed to assist in writing
    ! an EQ6 pickup file. These variables do not appear on that file
    ! itself.
    integer :: ibsrt1
    integer :: iesrt1

    character(len=24), dimension(:), allocatable :: ubsr1
    character(len=8), dimension(:), allocatable :: uesr1

    character(len=24) :: ureac1

    real(kind=8), dimension(:), allocatable :: cbsr1
    real(kind=8), dimension(:), allocatable :: cesr1

    ! Names of solid solution models.
    data uxtype(1) /'Ideal solution                  '/
    data uxtype(2) /'Binary, third-order Maclaurin   '/
    data uxtype(3) /'Binary, parabolic Maclaurin     '/
    data uxtype(4) /'Binary, cubic Maclaurin (P,T)   '/
    data uxtype(5) /'Binary, Guggenheim (T)          '/
    data uxtype(6) /'Ternary, regular                '/
    data uxtype(7) /'Plagioclase (Newton et al. 1980)'/

    ! Limit on the number of points of reaction progress beyond the
    ! starting point at which warnings are issued regarding fugacities
    ! not being at specified values.
    data nlwffg / 3 /

    ! The variable "alk" is not used in this subroutine except to
    ! satisfy certain EQLIB subroutine calls.
    data alk /0./

    data fxprpl /0.2/

    data qxbarw/.false./

    ! Allocate additional arrays needed to compute the reaction path.
    ALLOCATE(ctb(nbtmax))
    ALLOCATE(ppmwb(nbtmax))
    ALLOCATE(ehrc(nbtmax))
    ALLOCATE(cjbasp(nbtmax))
    ALLOCATE(ibswx(nbtmax))
    ALLOCATE(ixbasp(nbtmax))
    ALLOCATE(dlogxw(nbtmax))

    ALLOCATE(aamatr(kmax,kmax))
    ALLOCATE(gmmatr(kmax,kmax))

    ALLOCATE(iindx0(kmax))
    ALLOCATE(ipivot(kmax))
    ALLOCATE(ipndx0(kmax))

    ALLOCATE(alpha(kmax))
    ALLOCATE(beta(kmax))
    ALLOCATE(betao(kmax))
    ALLOCATE(delvco(kmax))
    ALLOCATE(delvec(kmax))
    ALLOCATE(d1zvc1(kmax))
    ALLOCATE(rhsvec(kmax))
    ALLOCATE(zvclg0(kmax))
    ALLOCATE(zvec0(kmax))

    ALLOCATE(uzvec0(kmax))
    ALLOCATE(uldel(kmax))
    ALLOCATE(ulbeta(kmax))

    ALLOCATE(a3bars(natmax))

    ALLOCATE(kction(nbtmax))

    ALLOCATE(ahrc(nbtmax))
    ALLOCATE(amtb(nbtmax))
    ALLOCATE(fo2lrc(nbtmax))
    ALLOCATE(mtb0(nbtmax))
    ALLOCATE(perc(nbtmax))

    ALLOCATE(igstak(ngtmax))
    ALLOCATE(jgstak(ngtmax))
    ALLOCATE(jgsort(ngtmax))

    ALLOCATE(cteaq(nctmax))
    ALLOCATE(ppmwe(nctmax))

    ALLOCATE(fsort(ngtmax))
    ALLOCATE(fugac(ngtmax))
    ALLOCATE(fugalg(ngtmax))

    ALLOCATE(ncmpe(2,npetmx))
    ALLOCATE(ncmpe0(2,npetmx))

    ALLOCATE(dzvc0(nrd1mx,kmax))
    ALLOCATE(dzvc0s(nrd1mx,kmax))
    ALLOCATE(fdzvm1(nrd1mx,kmax))
    ALLOCATE(fdzv0(nrd1mx,kmax))

    ALLOCATE(akmat0(nrd1mx,nrd1mx))
    ALLOCATE(akmat1(nrd1mx,nrd1mx))

    ALLOCATE(daffp0(nordmx,nptmax))
    ALLOCATE(fdafm1(nordmx,nptmax))
    ALLOCATE(fdaf0(nordmx,nptmax))

    ALLOCATE(dafrc0(nordmx,nrctmx))
    ALLOCATE(drer0(nrd1mx,nrctmx))
    ALLOCATE(drer0s(nrd1mx,nrctmx))
    ALLOCATE(fdarm1(nordmx,nrctmx))
    ALLOCATE(fdar0(nordmx,nrctmx))
    ALLOCATE(fdrem1(nordmx,nrctmx))
    ALLOCATE(fdre0(nordmx,nrctmx))
    ALLOCATE(fdre1(nordmx,nrctmx))
    ALLOCATE(fdrrm1(nrd1mx,nrctmx))
    ALLOCATE(fdrr0(nrd1mx,nrctmx))
    ALLOCATE(fdrr1(nrd1mx,nrctmx))

    ALLOCATE(demop0(nordmx,npetmx))
    ALLOCATE(fdpem1(nordmx,npetmx))
    ALLOCATE(fdpe0(nordmx,npetmx))

    ALLOCATE(demos0(nordmx,nsetmx))
    ALLOCATE(fdsem1(nordmx,nsetmx))
    ALLOCATE(fdse0(nordmx,nsetmx))

    ALLOCATE(daw0(nordmx))
    ALLOCATE(fdawm1(nordmx))
    ALLOCATE(fdaw0(nordmx))
    ALLOCATE(deh0(nordmx))
    ALLOCATE(fdehm1(nordmx))
    ALLOCATE(fdeh0(nordmx))
    ALLOCATE(do20(nordmx))
    ALLOCATE(fdo2m1(nordmx))
    ALLOCATE(fdo20(nordmx))
    ALLOCATE(dph0(nordmx))
    ALLOCATE(fdphm1(nordmx))
    ALLOCATE(fdph0(nordmx))

    ALLOCATE(dxsm00(nrd1mx))
    ALLOCATE(dxsm10(nrd1mx))
    ALLOCATE(dxsm11(nrd1mx))

    ALLOCATE(drir0(nrd1mx))
    ALLOCATE(drir0s(nrd1mx))

    ALLOCATE(fdrim1(nrd1mx))
    ALLOCATE(fdri0(nrd1mx))
    ALLOCATE(fdri1(nrd1mx))

    ALLOCATE(dxval0(nrd1mx))

    ALLOCATE(iemop(npetmx))
    ALLOCATE(iemop0(npetmx))

    ALLOCATE(d1emp1(npetmx))
    ALLOCATE(d2emp1(npetmx))
    ALLOCATE(emop(npetmx))
    ALLOCATE(emop0(npetmx))

    ALLOCATE(qxknph(nptmax))

    ALLOCATE(affp(nptmax))
    ALLOCATE(affp0(nptmax))
    ALLOCATE(affpd(nptmax))
    ALLOCATE(mophg(nptmax))
    ALLOCATE(mophj(nptmax))
    ALLOCATE(mopht(nptmax))
    ALLOCATE(moph0(nptmax))
    ALLOCATE(sidrph(nptmax))
    ALLOCATE(voph(nptmax))
    ALLOCATE(vophj(nptmax))
    ALLOCATE(vophg(nptmax))
    ALLOCATE(vopht(nptmax))
    ALLOCATE(woph(nptmax))
    ALLOCATE(wophg(nptmax))
    ALLOCATE(wophj(nptmax))
    ALLOCATE(wopht(nptmax))

    ALLOCATE(iexr(nrctmx))
    ALLOCATE(jexr(nrctmx))
    ALLOCATE(jreac0(nrctmx))
    ALLOCATE(jsca(nrctmx))
    ALLOCATE(jscr(nrctmx))

    ALLOCATE(afrc0(nrctmx))
    ALLOCATE(afrc1(nrctmx))
    ALLOCATE(afrcp(nrctmx))
    ALLOCATE(idirec(nrctmx))
    ALLOCATE(sfcar0(nrctmx))
    ALLOCATE(modr0(nrctmx))
    ALLOCATE(morr0(nrctmx))
    ALLOCATE(rreacn(nrctmx))
    ALLOCATE(rreac0(nrctmx))
    ALLOCATE(rreac1(nrctmx))
    ALLOCATE(rrelrp(nrctmx))
    ALLOCATE(rrelr0(nrctmx))
    ALLOCATE(rrelr1(nrctmx))
    ALLOCATE(wodr(nrctmx))
    ALLOCATE(worr(nrctmx))
    ALLOCATE(xirct(nrctmx))
    ALLOCATE(xirct0(nrctmx))

    ALLOCATE(rrxfi1(imchmx,nrctmx))

    ALLOCATE(iemos(nsetmx))
    ALLOCATE(iemos0(nsetmx))

    ALLOCATE(emos(nsetmx))
    ALLOCATE(emos0(nsetmx))

    ALLOCATE(istack(nstmax))
    ALLOCATE(jstack(nstmax))
    ALLOCATE(jcsort(nstmax))
    ALLOCATE(jjsort(nstmax))
    ALLOCATE(jssort(nstmax))

    ALLOCATE(acflg(nstmax))
    ALLOCATE(acflgo(nstmax))
    ALLOCATE(acflg0(nstmax))
    ALLOCATE(act(nstmax))
    ALLOCATE(actlg(nstmax))
    ALLOCATE(affs(nstmax))
    ALLOCATE(affsd(nstmax))
    ALLOCATE(cdrtw(nstmax))
    ALLOCATE(cdrw(nstmax))
    ALLOCATE(cnufac(nstmax))
    ALLOCATE(conc(nstmax))
    ALLOCATE(conclg(nstmax))
    ALLOCATE(lsort(nstmax))
    ALLOCATE(mospg(nstmax))
    ALLOCATE(mospj(nstmax))
    ALLOCATE(mospt(nstmax))
    ALLOCATE(mosp0(nstmax))
    ALLOCATE(sidrsp(nstmax))
    ALLOCATE(vosp(nstmax))
    ALLOCATE(vospg(nstmax))
    ALLOCATE(vospj(nstmax))
    ALLOCATE(vospt(nstmax))
    ALLOCATE(weight(nstmax))
    ALLOCATE(wosp(nstmax))
    ALLOCATE(wospg(nstmax))
    ALLOCATE(wospj(nstmax))
    ALLOCATE(wospt(nstmax))
    ALLOCATE(xbar(nstmax))
    ALLOCATE(xbarlg(nstmax))

    ! Allocate some additional arrays associated with the higher-order
    ! (stiff) ODE integrator.
    nrct1 = nrct + 1

    ALLOCATE(ipivtr(nrct1))

    ALLOCATE(alphar(nrct1))
    ALLOCATE(betar(nrct1))
    ALLOCATE(delvcr(nrct1))
    ALLOCATE(rhsvcr(nrct1))
    ALLOCATE(dvjdte(nrct))

    ALLOCATE(armatr(nrct1,nrct1))
    ALLOCATE(grmatr(nrct1,nrct1))

    ALLOCATE(aimatr(kmax,kmax))
    ALLOCATE(mmmatr(nrct,kmax))
    ALLOCATE(sgmatr(nrct,kmax))
    ALLOCATE(xxmatr(kmax,nrct))
    ALLOCATE(xymatr(nrct,nrct1))

    ALLOCATE(cdacb(nbt,imchmx,2,nrct))
    ALLOCATE(ndacb(nbt,imchmx,2,nrct))
    ALLOCATE(ndactb(imchmx,2,nrct))

    ALLOCATE(hhcvec(nordmx))
    ALLOCATE(xhcvec(nordmx))

    ! Allocate some arrays needed to assist in writing an EQ6 pickup
    ! file. These arrays do not appear on that file itself.
    ALLOCATE(uesr1(nctmax))
    ALLOCATE(cesr1(nctmax))
    ALLOCATE(ubsr1(nbt1mx))
    ALLOCATE(cbsr1(nbt1mx))

    ! Set the format ("W" or "D") for the pickup and backup files.
    if (iopr(17) .le. 0) then
        upkfor = uinfor
    else if (iopr(17) .eq. 1) then
        upkfor = 'W'
    else if (iopr(17) .ge. 2) then
        upkfor = 'D'
    end if

    ! Set the EQ6 calculational mode flag (.true. in EQ6).
    q6mode = .true.

    ! Set counters for warnings of problems associated with pressure
    ! corrections.
    !   nwndpc = number of warnings of no data to support pressure
    !              corrections
    !   nweope = number of warnings of excursions outside the
    !              recommended pressure envelope
    nwndpc = 0
    nweope = 0

    ! Set counter for the number of points of reaction progress at
    ! which warnings are issued regarding fugacities not being fixed at
    ! specified values.
    nwnffg = 0

    ! Set counters which control the writing of warnings that little
    ! solvent water remains in the equilibrium system.
    iwdh2o = 0
    nwdh2o = 0

    ! Set some constants that depend on the molecular weight of
    ! water.
    omega = 1000./mwtsp(narn1)
    omeglg = log10(omega)

    ! Initialize a value for the density of the aqueous solution.
    qrho = .false.
    rho = 0.0

    ! Reset the jflgi array so it contains the jflag values for
    ! the basis species.
    do nb = 1,nbt
        ns = nbasp(nb)
        jflgi(nb) = jflag(ns)
    end do

    ! Set up the cdrw array. This provides a fast way to get the
    ! reaction coefficient of H2O (Aqueous solution) in any reaction.
    call gcdrw(cdrs,cdrw,narn1,ndrs,ndrsmx,ndrsr,nst,nstmax)

    ! Set up the cdrtw array. This contains a stoichiometric sum of the
    ! number of times H2O (Aqueous solution) is implied as a solvent
    ! in any reaction.
    call gcdrtw(cdrs,cdrtw,narn1,narn2,ndrs,ndrsmx,ndrsr,nelect,no2gaq,nst,nstmax)

    ! Initialize concentration, activity coefficient, activity, etc.,
    ! arrays.
    do ns = 1,nstmax
        acflg(ns) = 0.
        acflgo(ns) = 0.
        act(ns) = 0.
        conc(ns) = 0.
        mosp(ns) = 0.
        xbar(ns) = 0.
        affs(ns) = 0.
    end do

    av = -99999.
    call initav(actlg,nstmax,av)
    call initav(conclg,nstmax,av)
    call initav(xbarlg,nstmax,av)
    call initav(loph,nptmax,av)

    nmax = nptmax
    call initaz(moph,nmax)
    call initaz(affp,nmax)

    nmax = nbtmax
    call initaz(amtb,nmax)

    nmax = 2*npetmx
    call initiz(ncmpe,nmax)

    ! Initialize (null) some other arrays.
    do ne = 1,netmax
        egexpa(ne) = 0.
        egexpc(ne) = 0.
    end do

    nmax = jetmax*netmax
    call initaz(egexjc,nmax)

    nmax = ietmax*jetmax*netmax
    call initaz(egexs,nmax)

    nmax = ketmax*netmax
    call initaz(egexw,nmax)

    ! Initialize the mole fractions and activities of pure solids and
    ! liquids.
    do ns = nlrn1,nlrn2
        act(ns) = 1.
        actlg(ns) = 0.
        xbar(ns) = 1.
        xbarlg(ns) = 0.
    end do

    do ns = nmrn1,nmrn2
        act(ns) = 1.
        actlg(ns) = 0.
        xbar(ns) = 1.
        xbarlg(ns) = 0.
    end do

    do ns = nfrn1,nfrn2
        act(ns) = 1.
        actlg(ns) = 0.
        xbar(ns) = 1.
        xbarlg(ns) = 0.
    end do

    ! Initialize the flag array which marks if compositions are known
    ! which maximize the affinities of the phases. Here assume that
    ! this is the case for all phases except solid solutions.
    do np = 1,npt
        qxknph(np) = .true.
    end do

    do np = ixrn1,ixrn2
        qxknph(np) = .false.
    end do

    if (net .gt. 0) then
        ! Initialize the numbers of moles of exchanger phases in the
        ! moph and loph arrays.
        do nb = 1,nbt
            ns = nbaspd(nb)

            if (ns.ge.nern1 .and. ns.le.nern2) then
                np = nphasx(ns)
                mx = mtb(nb)
                lx = tlg(mx)
                moph(np) = mx
                loph(np) = lx
            end if
        end do
    end if

    ! Initialize the jjsort, jssort, jcsort, and jgsort arrays by
    ! setting each element equal to its index.
    call initii(jcsort,nst)
    call initii(jjsort,nst)
    call initii(jssort,nst)
    call initii(jgsort,ngt)

    ! Expand the system description from the data read from the input
    ! file. This includes estimating the numbers of moles of all phases
    ! and species present, the concentrations, activity coefficients,
    ! and activities of all the species, the ionic strength, and the
    ! sum of the molalities of all aqueous aqueous solute species.
    call exivar(abar,acflg,acflgo,act,actlg,actwlc,adh,adhh,adhv,al10,aphi,azero,a3bar,a3bars,bdh,bdhh,bdhv,bdot,bdoth,bdotv,cco2,cdrs,cegexs,cgexj,conc,conclg,cpgexs,egexjc,egexjf,egexs,eps100,fje,fjeo,fo2,fo2lg,fsort,fugac,fugalg,fxi,fxio,ielam,iern1,iern2,ietmax,ifcphi1,ifcphi2,ifnnn,ifn2n,ifpsi1,ifpsi2,ifrn1,ifrn2,ifzeta,igas,igstak,iindx1,ilcphi1,ilcphi2,ilnnn,iln2n,ilpsi1,ilpsi2,ilrn1,ilrn2,ilzeta,imrn1,imrn2,insgf,iodb,iopg,ipbtmx,istack,ixrn1,ixrn2,izmax,jcsort,jern1,jern2,jetmax,jflag,jgext,jgsort,jgstak,jjsort,jpflag,jsflag,jsitex,jssort,jstack,kbt,kdim,kelect,kmax,km1,ko2gaq,kwater,kxt,loph,losp,lsort,mgext,moph,mosp,mrgexs,mtb,napmax,narn1,narn2,natmax,nazmmx,nazpmx,nbasp,nbaspd,nbt,nbtmax,nchlor,ncmpr,ndrs,ndrsmx,ndrsr,nelect,nern1,nern2,net,netmax,ngexsa,ngext,ngrn1,ngrn2,ngt,ngtmax,nhydr,nmutmx,nmxmax,nodbmx,nopgmx,noutpt,no2gaq,nphasx,npt,nptmax,nsltmx,nst,nstmax,nsxmax,nttyo,omega,omeglg,press,qhawep,qpit75,q6mode,sigmam,sigmmo,tempk,ugexj,ugexmo,uphase,uspec,xbar,xbarlg,xbarw,xbarwc,xbrwlc,xbrwlg,xlks,zchar,zchcu6,zchsq2,zgexj,zvclg1,zvec1)

    ! Set the composition known flag for solid solutions present
    ! in the ES.
    do np = ixrn1,ixrn2
        if (moph(np) .gt. 0.) then
            qxknph(np) = .true.
        end if
    end do

    ! Set up to make the initial equilibrium calculation.
    !   kord   = the largest possible order on a given step
    !   kstep  = step counter
    !   nord   = the actual order used
    !   npts   = the number of points available; restricts the
    !              value of kord
    kstep = 0
    kord = 0
    nord = 0
    npts = 1
    ncorr = 0

    km10 = km1
    kmt0 = kmt
    kx10 = kx1
    kxt0 = kxt
    kdim0 = kdim

    call copyia(iindx1,iindx0,kmax)
    call copyia(ipndx1,ipndx0,kmax)
    call copyca(uzvec1,uzvec0,kmax)
    call copyaa(zvclg1,zvclg0,kmax)
    call copyaa(zvec1,zvec0,kmax)

    qstart = .true.
    nordsv = 0

    xim1 = xistrt
    xi0 = xistrt
    xi1 = xistrt
    xistsv = xistrt

    xidump = 0.

    qriinf = .false.

    time0 = tistrt
    time1 = tistrt
    tistsv = tistrt

    delxi = 0.
    deltim = 0.

    rirec0 = 0.
    rirec1 = 0.

    do n = 1,nrctmx
        afrc0(n) = 0.
        afrc1(n) = 0.
        rreac0(n) = 0.
        rreac1(n) = 0.
        rrelr0(n) = 0.
        rrelr1(n) = 0.
    end do

    call copyia(jreac,jreac0,nrctmx)
    call copyaa(morr,morr0,nrctmx)
    call copyaa(modr,modr0,nrctmx)
    call copyaa(mtb,mtb0,nbtmax)
    call copyaa(sfcar,sfcar0,nrctmx)

    call initaz(xirct,nrctmx)
    call initaz(xirct0,nrctmx)

    iexrt = 0
    jexrt = 0
    jscat = 0
    jscrt = 0

    betamx = 0.
    ubacmx = 'None'
    ubgamx = 'None'

    avkdim = real(kdim)
    prminf = 1./prcinf
    screwd = sscrew(5)

    tolsar = 1.e-8
    tolsrr = 1.e-8

    tiprnt = prcinf
    tiprnl = prcinf
    tiplot = prcinf
    tiplol = prcinf

    qbye = .false.
    qconst = .false.
    qdump = .false.
    qftpr2 = .false.
    qmod = .false.
    qrapch = .false.
    qriinf = .false.
    qtvchk = .false.
    qtrch = .false.
    qzprnt = .false.
    qstopx = .false.
    qshoot = .false.
    qvlsow = .false.
    qvhfxi = .false.

    ! Zero various arrays pertaining to finite differences and
    ! equivalent derivatives.
    do n = 1,nrd1mx
        dxsm00(n) = 0.
        dxsm10(n) = 0.
        dxsm11(n) = 0.
        fdrim1(n) = 0.
        fdri0(n) = 0.
        fdri1(n) = 0.
        drir0(n) = 0.
        drir0s(n) = 0.
    end do

    nmax = nrd1mx*kmax
    call initaz(fdzv0,nmax)
    call initaz(fdzvm1,nmax)
    call initaz(dzvc0,nmax)
    call initaz(dzvc0s,nmax)

    do n = 1,nordmx
        fdawm1(n) = 0.
        fdaw0(n) = 0.
        daw0(n) = 0.
        fdehm1(n) = 0.
        fdeh0(n) = 0.
        deh0(n) = 0.
        fdo2m1(n) = 0.
        fdo20(n) = 0.
        do20(n) = 0.
        fdphm1(n) = 0.
        fdph0(n) = 0.
        dph0(n) = 0.
    end do

    nmax = nordmx*nptmax
    call initaz(fdaf0,nmax)
    call initaz(fdafm1,nmax)
    call initaz(daffp0,nmax)

    nmax = nordmx*npetmx
    call initaz(fdpe0,nmax)
    call initaz(fdpem1,nmax)
    call initaz(demop0,nmax)

    nmax = 2*npetmx
    call initiz(ncmpe0,nmax)

    nmax = nordmx*nsetmx
    call initaz(fdse0,nmax)
    call initaz(fdsem1,nmax)
    call initaz(demos0,nmax)

    nmax = nordmx*nrctmx
    call initaz(fdarm1,nmax)
    call initaz(fdar0,nmax)
    call initaz(dafrc0,nmax)
    call initaz(fdrem1,nmax)
    call initaz(fdre0,nmax)
    call initaz(fdre1,nmax)

    nmax = nrd1mx*nrctmx
    call initaz(fdrrm1,nmax)
    call initaz(fdrr0,nmax)
    call initaz(fdrr1,nmax)
    call initaz(drer0,nmax)
    call initaz(drer0s,nmax)

    ! Set up the ixbasp and cjbasp arrays. The former is a flag
    ! array, each member of which denotes whether the
    ! thermodynamic activity of the corresponding basis species
    ! is defined in terms of molality (= 0) or mole fraction (= 1).
    ! The cjbasp array contains any site stoichiometric factors
    ! associated with the operational basis species.
    call gibasp(cgexj,cjbasp,iern1,ixbasp,jern1,jern2,jetmax,jgext,narn1,narn2,nbasp,nbt,nbtmax,nern1,nern2,netmax,nphasx,nstmax)

    ! Make the equilibrium calculation at the initial point of reaction
    ! progress.
    call eqshel(aadh,aadhh,aadhv,aamatr,aaphi,abar,abdh,abdhh,abdhv,abdot,abdoth,abdotv,acflg,acflgo,act,actlg,adadhh,adadhv,adbdhh,adbdhv,adbdth,adbdtv,adh,adhh,adhv,afcnst,affp,affs,alpha,al10,amtb,aphi,apx,avcnst,azero,a3bar,a3bars,bacfmx,bbig,bdh,bdhh,bdhv,bdot,bdoth,bdotv,beta,betamx,betao,bgamx,bneg,bpx,cbsr,cco2,cdac,cegexs,cesr,cess,cdrs,cdrsd,cdrsx,cdrtw,cdrw,cjbasp,cnufac,conc,conclg,cpgexs,cscale,csigma,csts,dadhh,dadhv,dbdhh,dbdhv,dbdth,dbdtv,deltim,delvco,delvec,delxi,dlogxw,dlxmin,drer0,drir0,dzvc0,d1zvc1,eact,egers,egexjc,egexjf,egexs,eh,ehfac,elecsr,eps100,farad,fje,fjeo,fkrc,fo2,fo2lg,fsort,fugac,fugalg,fxi,fxio,gmmatr,hact,iact,iapxt,ibpxt,ibswx,ielam,ier,iern1,iern2,ifcphi1,ifcphi2,ifnnn,ifn2n,ifpsi1,ifpsi2,ifrn1,ifrn2,ifzeta,igas,igstak,iindx0,iindx1,ilcphi1,ilcphi2,ilnnn,iln2n,ilpsi1,ilpsi2,ilrn1,ilrn2,ilzeta,imech,imrn1,imrn2,insgf,iodb,iopg,iopt,ipch,ipivot,ipndx1,ipcv,istack,iter,itermx,ixbasp,ixrn1,ixrn2,izmax,jcode,jcsort,jflag,jgsort,jgstak,jjsort,jpflag,jpress,jptffl,jreac,jsflag,jsitex,jsol,jssort,jstack,jtemp,kbt,kction,kdim,kelect,khydr,khydx,km10,km1,kmt,kmt0,ko2gaq,kpsat,kpsst,krdxsp,kwater,kx1,kx10,kxt,kxt0,loph,losp,lsort,modr,modr0,moph,morr,morr0,mosp,mrgers,mrgexs,mtb,mtbaq,mtb0,mte,mteaq,mwtrc,narn1,narn2,narxt,nat,nbasp,nbaspd,nbaspx,nbt,nbtd,nbw,nchlor,ncmpr,ncorr,nct,ndac,ndact,ndrs,ndrsd,ndrsx,ndrsr,ndrsrd,ndrsrx,nelect,nern1,nern2,ness,nessr,net,nfrn1,nfrn2,ngrn1,ngrn2,ngt,nhydr,nhydx,nlrn1,nlrn2,nlt,nmrn1,nmrn2,nmt,nord,noutpt,no2gaq,npchk,nphasx,npslmx,npt,nrct,nrdxsp,nrk,nrndex,nsk,nsslmx,nst,nsts,nstsr,ntpr,ntrymx,ntprt,nttyo,nxridx,nxrn1,nxrn2,nxt,nweope,nwndpc,omega,omeglg,prcinf,press,pressb,pressd,ptk,qbassw,qbseqc,qbye,qcnpre,qcntmp,qhawep,qmod,qoptmz,qpit75,qredox,qriinf,qscon,qshoot,qstart,qtrch,qtvchk,qxknph,q6mode,rcnstv,rconst,rkb,rhsvec,rirec0,rk,rreacn,rreac1,rrelr0,rrelr1,rtcnst,rxbar,screwd,sidrph,sidrsp,sigmam,sigmmo,smp100,tdays,tempc,tempcb,tempcd,tempcu,tempc0,tempk,timemx,time0,time1,tiplol,tiplot,tiprnl,tiprnt,tistsv,tolbt,toldl,tolsat,tolsst,tolxst,trkb,ttk,ubacmx,ubgamx,udac,ufixf,ugermo,ulbeta,uldel,uphase,ureac,uspec,uzvec0,uzvec1,vreac,weight,wfac,wodr,worr,xbar,xbarlg,xbarw,xbarwc,xbrwlc,xbrwlg,xgers,xirct,xirct0,xistsv,xi0,xi1,zchar,zchcu6,zchsq2,zklogu,zvclg0,zvclg1,zvec0,zvec1)

    ! Note on ier codes returned by EQ6/eqshel.f:
    !   =    0  Okay
    !   =   10  Go back and take a smaller step size to avoid exceeding
    !             the supersaturation tolerance (tolsst)
    !   =  170  Too much of a phase was destroyed under the flow-through
    !             open system model; go back and first move part of the
    !             mass of protected phases in the ES to the PRS
    !   =  180  One of a number of problems occurred which may be
    !             resolvable, at least partially, by going back and
    !             cutting the step size
    !   =  190  Need to slide over a region of computational
    !             instability, but sliding is inhibited; go back, but
    !             terminate work on the current problem
    ! Note: EQ6/eqshel.f can't return ier = 10 or 170 when called at
    ! the initial point of reaction progress.
    if (ier .eq. 180) then
        ! Normally this error code results in a reduction in the step
        ! size. This can't be done at the initial point of reaction
        ! progress.
        write (noutpt,1010)
        write (nttyo,1010)
1010 format(/' * Error- (EQ6/path) The equilibrium calculation',' failed at the initial value',/7x,'of reaction progress.'," Can't cut the step size to try to recover.",/7x,'See',' previous notes and warnings for more information.')

        stop
    end if

    if (ier .gt. 0) then
        ! Here ier should have a value of 190.
        write (noutpt,1020)
        write (nttyo,1020)
1020 format(/' * Error - (EQ6/path) The equilibrium calculation',' failed at the initial value',/7x,'of reaction progress.')

        stop
    end if

    ! Get the activity of water.
    actwlg = actlg(narn1)
    actw = texp(actwlg)
    awstrt = actw
    aw0 = actw
    aw1 = actw

    ! Get the weights (masses) of solvent, total dissolved solutes,
    ! and aqueous solution, and get the aqeuous solution density.
    call gwdenp(adwipp,bdwipp,jcsort,mlmrra,mosp,mrmlra,mwtsp,narn1,narn2,nstmax,qdwipp,rhoc,rhowc,tdsgks,tdsglw,tdspkc,tdsplc,tempc,vosol,wfh2o,wftds,wkgwi,woh2o,wosol,wotds)

    if (qdwipp) then
        qrho = .true.
        rho = rhoc
    else
        qrho = .false.
        rho = 0.0
    end if

    ! Compute pH, Eh, and pe-, all with reference to appropriate
    ! pH scales. Also compute the pHCl.
    call gpheh(acflg,actlg,actwlg,adh,ah,ahmes,ahnbs,conc,eh,ehfac,ehmes,ehnbs,farad,fo2lg,fxi,iopg,mrmlra,nchlor,nhydr,nopgmx,noutpt,nstmax,nttyo,pch,pe,pemes,penbs,ph,phcl,phmes,phnbs,qphcl,qredox,qrho,xlke)

    wkgh2o = 1.e-3*woh2o
    wkgsol = 1.e-3*wosol
    wwstrt = wkgh2o
    phstrt = ph
    ph0 = ph
    ph1 = ph
    ehstrt = eh
    eh0 = eh
    eh1 = eh
    o2strt = fo2lg
    fo2lg0 = fo2lg
    fo2lg1 = fo2lg

    ! Do automatic basis switching. Optimize the basis set so that
    ! the basis species for each mass balance tends to be the species
    ! which dominates that mass balance.
    kstpab = 0

    if (iopt(12) .gt. 0) then
        call absswb(adhfs,adhfsx,advfs,advfsx,avcnst,axhfs,axhfsx,axlks,axlksx,axvfs,axvfsx,cdrs,cdrtw,cdrsx,cdrw,csts,dhfs,dvfs,eps100,ibswx,iindx1,iodb,ipch,ipchmx,ipcv,ipcvmx,jcsort,jflag,jsflag,kbt,kmax,mosp,mtb,narn1,narn2,narxmx,narxt,nbasp,nbaspd,nbaspx,nbt,nbtmax,nbw,ndrs,ndrsmx,ndrsr,ndrsrx,ndrsx,nelect,nhydr,nodbmx,no2gaq,noutpt,nst,nstmax,nsts,nstsmx,nstsr,nswtch,ntpr,ntprmx,nttyo,presg,press,qbassw,qbswx,tempc,uspec,uzvec1,weight,xhfs,xvfs,xlks)

        if (nswtch .gt. 0) then
            ! It is not necessary here to set the qabswx flag, which
            ! indicate that automatic basis switching has just been done.
            write (noutpt,1030) nswtch
1030 format(/' ',i2,' basis switches were executed automatically',' after solving',/'at the initial value of reaction',' progress.')

            ! Reset the ixbasp and cjbasp arrays. The former is a flag
            ! array, each member of which denotes whether the
            ! thermodynamic activity of the corresponding basis species
            ! is defined in terms of molality (= 0) or mole fraction (= 1).
            ! The cjbasp array contains any site stoichiometric factors
            ! associated with the operational basis species.
            call gibasp(cgexj,cjbasp,iern1,ixbasp,jern1,jern2,jetmax,jgext,narn1,narn2,nbasp,nbt,nbtmax,nern1,nern2,netmax,nphasx,nstmax)
        end if
    end if

    write (noutpt,1100) kstep,iter
1100 format(' Steps completed= ',i5,', iter= ',i3)

    ! Save some data at the initial point, as it is after the
    ! calculations at this point are complete. Note that the set of
    ! matrix variables may have been changed due to mineral
    ! precipitation, etc.
    call copyia(iindx1,iindx0,kdim)
    call copyia(ipndx1,ipndx0,kdim)

    km10 = km1
    kmt0 = kmt
    kx10 = kx1
    kxt0 = kxt
    kdim0 = kdim

    ! Check gases for which the fugacity is supposed to be fixed.
    do n = 1,nffg
        ng = jffg(n)
        xlf = xlkffg(n)
        xf = texp(xlf)
        dx = fugalg(ng) - xlkffg(n)

        if (abs(dx) .gt. toldl) then
            j2 = ilnobl(uffg(n))
            write (noutpt,1150) uffg(n)(1:j2),fugac(ng),xf
            write (nttyo,1150) uffg(n)(1:j2),fugac(ng),xf
1150 format(/' * Warning - (EQ6/eq6) The fugacity of ',a,' is',/7x,1pg12.5,' bars at the start of the run. It is supposed',/7x,'to be fixed at ',e12.5,' bars. There is an insufficient',/7x,'mass of the gas component present in the system to',/7x,'"saturate" it at the desired fugacity. You may wish to',/7x,'restart this run, adding such mass using the moffg',/7x,'parameter on the input file. Add only about 0.5 to 1.0',/7x,'mole at a time.')
        end if
    end do

    if (nrct .gt. 0) then
        ! Calculate the initial affinities and rates of the irreversible
        ! reactions.
        call raff(acflg,actlg,afcnst,affp,afrc1,bpx,cdrs,cgexj,ibpxmx,ibpxt,iern1,ietmax,iktmax,ixrn1,ixrn2,jcode,jern1,jern2,jetmax,jflag,jgext,jpflag,jsflag,jsol,ncmpr,ndrs,ndrsmx,ndrsr,nertmx,net,netmax,ngext,noutpt,nptmax,nrct,nrctmx,nrndex,nstmax,nttyo,nxridx,nxrtmx,nxtmax,rxbar,uphase,uspec,wfac,xbar,xbarlg,xgers,xlks)

        ! For each reactant, check the status flag against the affinity
        ! and the number of moles present to ensure that any reactant
        ! tagged as available to react, having an affinity favoring
        ! dissolution, and having no moles remaining is retagged as
        ! exhausted.
        do nrc = 1,nrct
            if (afrc1(nrc) .gt. -tolsar) then
                if (jreac(nrc).eq.0 .or. jreac(nrc).eq.-1) then
                    if (morr(nrc) .le. 0.) then
                        ! The following exception is intended to facilitate the
                        ! simulation of evaporation by declaring H2O as a
                        ! reactant with a negative relative rate and no initial
                        ! moles remaining. However, any reactant may be treated
                        ! analogously if desired.
                        if (.not.(nrk(1,nrc).eq.1 .and. rk(1,1,nrc).lt.0.)) then
                            jreac(nrc) = 1
                        end if
                    end if
                end if
            end if
        end do

        ! Check the reactants for saturation.
        call rsatch(csts,egers,egexs,iern1,ietmax,iindx1,iktmax,iopt,ipndx1,jcode,jern1,jern2,jetmax,jgext,jpflag,jreac,kmax,km1,kmt,kx1,kxt,loph,losp,moph,morr,mosp,mrgers,mtb,mtb0,nbaspd,nbtmax,ncmpr,nern1,nern2,nert,nertmx,netmax,ngext,noptmx,noutpt,nptmax,nrct,nrctmx,nrk,nrndex,nstmax,nsts,nstsmx,nstsr,nttyo,nxridx,nxrt,nxrtmx,qreq,rxbar,tolxsf,uphase,ureac,uspec,xbar,xbarlg,zvclg1,zvec1)

        ! Calculate rates at the initial point by evaluating the
        ! rate laws.
        call rtcalc(act,afrc1,cdac,csigma,eps100,fkrc,idirec,imchmx,imech,iodb,iopt,jcode,jreac,morr,morr0,mwtrc,ndac,ndact,ndctmx,nodbmx,noptmx,nord,noutpt,nrk,nrct,nrctmx,nsk,nstmax,nttyo,prcinf,prminf,qriinf,rirec1,rk,rreac1,rrelr1,rtcnst,rrxfi1,sfcar,sfcar0,ssfcar,udac,ureac)
    end if

    ! Calculate the total affinity (aft1).
    call caft1(afrc1,aft1,nrct,nrctmx,rrelr1)
    aft0 = aft1
    aftm1 = aft1
    qaft1 = .false.

    ! Compute apparent "whole-phase" equivalent fractions and mole
    ! fractions of the exchange ions present in generic ion exchanger
    ! phases. Cations and anions are treated separately in these
    ! calculations.
    call gegexw(cegexs,egexpc,egexpa,egexw,iern1,iern2,ietmax,jern1,jetmax,jgext,kern1,kern2,ketmax,kgexsa,moph,mosp,netmax,ngexsa,ngext,noutpt,nptmax,nstmax,nttyo,xgexw,zchar)

    ! Initialize arrays associated with finite-difference description
    ! of the number of mole of phases and species in the Equilibrium
    ! System (ES). Here the iemop and emop arrays respectively contain
    ! the indices and numbers of moles of phases in the ES. The iemos
    ! and emos array are the analogs for the species of these phases
    ! (but for the aqueous solution phase, only the species H2O(l)
    ! is tracked by this mechanism). The ncmpe array is a species range
    ! pointer array for the phases, analogous to ncmpr. The fdpe0 and
    ! fdse0 arrays contain the finite differences for the phases and
    ! species, respectively. The demop and demos arrays contain the
    ! corresponding derivatives.
    ! Initialize the index arrays.
    call iiemop(iemop,iemos,iindx1,ipndx1,jsflag,kdim,kmax,ncmpe,ncmpr,noutpt,npet,npetmx,npt,nptmax,nset,nsetmx,nstmax,nttyo,uaqsln,uspec,uphase)

    ! Load the corresponding numbers of moles.
    do npe = 1,npet
        np = iemop(npe)
        emop(npe) = moph(np)
        nr1 = ncmpe(1,npe)
        nr2 = ncmpe(2,npe)

        do nse = nr1,nr2
            ns = iemos(nse)
            emos(nse) = mosp(ns)
        end do
    end do

    npet0 = 0

    ! Calculate various secondary data at the initial point.
    call cdappl(acflg,acfw,acfwlg,actlg,actw,actwlg,adwipp,afcnst,affpd,affsd,ah,ahrc,alk,alk1,alk2,alki,atwt,bdwipp,cdrsd,cess,conc,csts,ctb,cteaq,dvoso,dwoso,eh,ehfac,ehrc,eps100,farad,fdpe0,fdse0,fjest,fo2lg,fo2lrc,fxist,iaqsln,iemop0,iemos0,iern1,iern2,ifrn1,ifrn2,ilrn1,ilrn2,imrn1,imrn2,iopt,ixrn1,ixrn2,jcode,jcsort,jern1,jflag,jflagd,jgext,jpflag,jsflag,jssort,modr,moph,mophg,mophj,mopht,morr,mosp,mospg,mospj,mospt,mprph,mprsp,mrgers,mrmlra,mtb,mtbaq,mte,mteaq,mwtges,mwtrc,mwtsp,narn1,narn2,nat,nbasp,nbaspd,nbt,nchlor,ncmpe0,ncmpr,nct,ndrsd,ndrsrd,nelect,nern1,nern2,nert,ness,nessr,net,nfrn1,nfrn2,ngext,ngrn1,ngrn2,ngt,nhydr,nlrn1,nlrn2,nlt,nmrn1,nmrn2,nmrt,nmt,no2gaq,npchk,npet,npet0,npt,npts,nrct,nrndex,nst,nsts,nstsr,ntf1,ntf1t,ntf2,ntf2t,nxridx,nxrn1,nxrn2,nxrt,nxt,osc,oscst,omega,pe,perc,ph,phmes,ppmwb,ppmwe,qriinf,qxknph,rreacn,rreac1,rxbar,sfcar,sidrsp,sidrph,sigmam,sigmst,tdays,tempc,tf1,tf2,thours,time1,tmins,tyears,uphase,uspec,vodrt,voph,vophg,vophj,vopht,vosoct,vosp,vospg,vospj,vospt,vosp0,vreac,wfh2o,wkgwi,wodr,wodrt,woph,wophg,wophj,wopht,worr,worrt,wosoct,wosp,wospg,wospj,wospt,xbar,xlke,xlksd,zchcu6,zchsq2)

    ! Print results at the initial point.
    call scripz(abar,acflg,acfw,acfwlg,actlg,actw,actwlg,affpd,affsd,afrc1,aft1,ah,ahmes,ahnbs,ahrc,alki,alk1,alk2,awmax,awmin,a3bar,cbsr,cdrsd,cegexs,cesr,conc,conclg,csts,ctb,cteaq,dvoso,dwoso,egers,egexjc,egexjf,egexpa,egexpc,egexs,egexw,eh,ehmax,ehmes,ehmin,ehnbs,ehrc,elecsr,electr,fje,fjest,fo2,fo2lg,fo2lrc,fugac,fugalg,fxi,fxist,iaqsln,iemop,iemop0,iemos,iemos0,iern1,iern2,iexr,iexrt,ifrn1,ifrn2,ilrn1,ilrn2,imech,imrn1,imrn2,iopg,iopr,iopt,ipndx1,ixrn1,ixrn2,jcode,jcsort,jern1,jern2,jexr,jexrt,jflag,jflagd,jflgi,jgext,jgsort,jpflag,jreac,jsca,jscat,jscr,jscrt,jsflag,jsol,jssort,kbt,kern1,kern2,kgexsa,km1,kmt,kx1,kxt,kstep,kstpmx,loph,losp,mlmrra,modr,moph,mophg,mophj,mopht,morr,mosp,mospg,mospj,mospt,mprph,mprsp,mrgers,mrmlra,mwtrc,mwtsp,narn1,narn2,nat,nbasp,nbaspd,nbt,ncmpe,ncmpe0,ncmpr,nct,ndrsd,ndrsrd,nelect,nern1,nern2,nert,net,nfrn1,nfrn2,ngext,ngexsa,ngrn1,ngrn2,ngt,nhydr,nhydx,nlrn1,nlrn2,nlt,nmrn1,nmrn2,nmrt,nmt,noutpt,no2gaq,npet,npet0,npt,npts,nrct,nrdxsp,nrk,nrndex,nst,nsts,nstsr,ntf1t,ntf2t,nxridx,nxrn1,nxrn2,nxrt,nxt,osc,oscst,omega,o2max,o2min,pch,pe,pemes,penbs,perc,ph,phcl,phmax,phmes,phmin,phnbs,ppmwe,presg,press,qaft1,qftpr2,qmod,qphcl,qredox,qrho,qriinf,qstopx,qvhfxi,qvlsow,qzprnt,rho,rhoc,rhowc,rk,rreacn,rreac1,rrelr1,rxbar,sfcar,sidrph,sidrsp,sigmst,sigmam,ssfcar,tdays,tdsglw,tdspkc,tdsplc,tempc,thours,time1,timemx,tmins,tolsat,tolxsf,tolxst,tolxsu,tyears,uelem,ugermo,ugexj,ugexmo,uphase,ureac,uspec,uxtype,vodrt,voph,vophg,vophj,vopht,vosoct,vosol,vosp,vospg,vospj,vospt,vreac,wfh2o,wftds,wkgwi,woh2o,wodr,wodrt,woph,wophg,wophj,wopht,worr,worrt,wosoct,wosol,wosp,wospg,wospj,wospt,wotds,xbar,xbarlg,xbarw,xbrwlg,xgers,xgexw,xi1,xidump,ximax,xistsv,xirct,zchar)

    ! Write results at the initial point on the tabx file. This file
    ! is used to create the plot file tab.
    if (iopt(18) .ge. 0) then
        if (qtatxt) then
            ! The TAB file is an ordinary text file.
            call wrtabx(actlg,afrc1,aft1,alk,cteaq,dvoso,dwoso,eh,fo2lg,iindx1,iktmax,iopt,ipndx1,kmax,km1,kmt,kstep,kx1,kxt,loph,ncmpr,modr,mopht,narn1,mosp,nct,nctmax,noptmx,nptmax,nrct,nrctmx,nstmax,ntabx,ntidmx,ntitl2,ntitld,ntitmx,nxtmax,pe,ph,ppmwe,prcinf,press,prminf,qbye,qmod,qriinf,tempc,time1,uelem,uphase,uplatm,ureac,uspec,usteq6,utitl2,utitld,uveeq6,vodrt,vosoct,wodrt,woh2o,wosoct,xbar,xi1)
        else
            ! The TAB file is a .csv file.
            call wrtabc(acflg,actlg,actw,afrc1,aft1,alk,conclg,cteaq,ctb,dvoso,dwoso,eh,fje,fo2lg,fugac,fxi,iktmax,iopt,jflag,jsflag,kmax,kstep,kx1,kxt,mrmlra,modr,mosp,mospt,moph,mopht,narn1,narn2,nbasp,nbaspd,nbt,nbtmax,ncmpr,nct,nctmax,nelect,ngrn1,ngrn2,ngtmax,nhydr,nhydx,nllnmx,no2gaq,noptmx,noutpt,npt,nptmax,nrct,nrctmx,nstmax,ntabx,ntidmx,ntitl2,ntitld,ntitmx,nttyo,nxrn1,nxrn2,nxtmax,pe,ph,phmes,ppmwb,ppmwe,prcinf,press,prminf,qrho,qriinf,rho,rhowc,sidrph,sigmam,tdsgks,tdsglw,tempc,time1,uelem,ulinex,uphase,uplatm,ureac,uspec,usteq6,utitl2,utitld,uveeq6,vodrt,vosoct,wkgh2o,wodrt,wosoct,xbar,xbarlg,xi1)
        end if
    end if

    write (nttyo,1200)
1200 format(1x)

    ! Write entertainment for the user, showing the progress
    ! of the current run.
    call wrentu(actw,eh,fo2lg,iopg,iopt,kstep,nopgmx,noptmx,nttyo,ph,qredox,time1,xi1)

    if (iopt(6) .gt. 0) then
        ! Clear ES solids at the starting point. Fictive fugacity-fixing
        ! minerals are not cleared.
        write (noutpt,1230)
        write (nttyo,1230)
1230 format(/' * Note (EQ6/path) Clearing solid phases from the',' equilibrium system (ES)',/7x,'at the starting point of',' reaction progress.')

        call clress(csts,iindx1,ipndx1,jpflag,jsflag,kdim,kmax,km1,kmt,kx1,kxt,loph,losp,moph,mosp,mtb,mtbaq,nbt,nbtmax,nptmax,nstmax,nsts,nstsmx,nstsr,ufixf,uzvec1,zvec1,zvclg1)
    end if

    ! Make sure that there is something to define a reaction path
    ! (reactants, changing temperature, changing pressure).
    if (nrct.le.0 .and. qcntmp .and. qcnpre) then
        write (noutpt,1250)
        write (nttyo,1250)
1250 format(/' * Note - (EQ6/path) No reaction path has been',' defined by on the input file.',/7x,'There are no specified',' reactants (irreversible reactions), and',/7x,'the',' temperature and pressure are fixed. No steps will be taken.')

        go to 990
    end if

    ! Check stop conditions.
    call chkstc(actw,awmax,awmin,eh,ehmax,ehmin,fo2lg,iopt,jreac,kstep,kstpmx,noptmx,noutpt,nrct,nrctmx,nttyo,o2max,o2min,ph,phmax,phmin,prcinf,qaft1,qcnpre,qcntmp,qconst,qredox,qstop,qvhfxi,qvlsow,timemx,time1,tolxst,tolxsu,ximax,xi1)

    if (qstop) then
        go to 990
    end if

    nbkupn = nbkupa

    if (iopt(16) .ge. 0) then
        ! Prepare to write results at the initial point to the
        ! backup file.
        call setpk6(actwlg,awmax,awmaxi,awmin,awmini,eh,ehmax,ehmaxi,ehmin,ehmini,fo2lg,iindx1,jflag,jflgi,kbt,kdim,kmax,kprs,mprph,mprphi,mprsp,mprspi,mtb,mtbi,mtbaq,mtbaqi,nbasp,nbaspd,nbaspi,nbti,nbtmax,ncmpr,nobswt,noutpt,nprpmx,nprpti,nprsmx,nprsti,npt,nptmax,nttyo,nstmax,o2max,o2maxi,o2min,o2mini,ph,phmax,phmaxi,phmin,phmini,prcinf,press,pressi,tempc,tempci,time1,timemx,timmxi,tistti,ubmtbi,uobsw,uphase,uprphi,uprspi,uspec,uzveci,uzvec1,xi1,ximax,ximaxi,xistti,zvclgi,zvclg1)

        ! Rewind the current backup file (BAKUPA or BAKUPB).
        rewind (nbkupn)

        ! Write results at the initial point to the backup file.
        if (upkfor(1:1) .eq. 'W') then
            ! Compact (W) format.
            ! Calling sequence substitutions:
            !   nbkupn for newin
            call wr6pkw(awmaxi,awmini,cbsri,cdac,cesri,cgexj,csigma,dlaplo,dlaprn,dleplo,dleprn,dlhplo,dlhprn,dloplo,dloprn,dltpll,dltplo,dltprl,dltprn,dlxdmp,dlxmx0,dlxpll,dlxplo,dlxprl,dlxprn,eact,egersi,ehmaxi,ehmini,electr,fkrc,iact,ibsrti,iesrti,ietmax,igerti,iktmax,imchmx,imech,iodb,iopg,iopr,iopt,itermx,ixrti,jcode,jetmax,jflgi,jgerti,jgext,jpress,jreac,jtemp,kbt,kct,kdim,kmax,kmt,kprs,ksplmx,ksppmx,kstpmx,kxmod,kxt,hact,modr,moffg,morr,mprphi,mprspi,mtbaqi,mtbi,mwtges,nbti,nbtmax,nbt1mx,nctmax,ndact,ndctmx,nert,nertmx,net,netmax,nbkupn,nffg,nffgmx,ngexrt,nobswt,nodbmx,nopgmx,noprmx,noptmx,nordmx,nprpmx,nprpti,nprsmx,nprsti,nptkmx,nrct,nrctmx,nrk,nsbswt,nsk,nsrt,nsrtmx,ntitl1,ntitl2,ntitmx,ntrymx,nttkmx,nxmdmx,nxmod,nxopex,nxopmx,nxopt,nxpemx,nxrt,nxrtmx,o2maxi,o2mini,phmaxi,phmini,pressb,pressi,ptk,qgexsh,rkb,rxbari,sfcar,ssfcar,tempcb,tempci,tgexp,timmxi,tistti,tolbt,toldl,tolsat,tolxsf,trkb,ttk,ubmtbi,ubsri,ucxri,udac,uesri,uffg,ugerji,ugermo,ugersi,ugexj,ugexmo,ugexp,ugexr,uhfgex,uobsw,uprphi,uprspi,ureac,usbsw,utitl1,utitl2,uvfgex,uxcat,uxkgex,uxmod,uxopex,uxopt,uzveci,vreac,xgersi,xhfgex,ximaxi,xistti,xlkffg,xlkgex,xlkmod,xvfgex,zgexj,zvclgi)
        else
            ! Menu-style (D) format.
            ! Calling sequence substitutions:
            !   nbkupn for newin
            call wr6pkd(awmaxi,awmini,cbsri,cdac,cesri,cgexj,csigma,dlaplo,dlaprn,dleplo,dleprn,dlhplo,dlhprn,dloplo,dloprn,dltpll,dltplo,dltprl,dltprn,dlxdmp,dlxmx0,dlxpll,dlxplo,dlxprl,dlxprn,eact,egersi,ehmaxi,ehmini,electr,fkrc,iact,ibsrti,iesrti,ietmax,igerti,iktmax,imchmx,imech,iodb,iopg,iopr,iopt,itermx,ixrti,jcode,jetmax,jflgi,jgerti,jgext,jpress,jreac,jtemp,kbt,kct,kdim,kmax,kmt,kprs,ksplmx,ksppmx,kstpmx,kxmod,kxt,hact,modr,moffg,morr,mprphi,mprspi,mtbaqi,mtbi,mwtges,nbti,nbtmax,nbt1mx,nctmax,ndact,ndctmx,nert,nertmx,net,netmax,nbkupn,nffg,nffgmx,ngexrt,nobswt,nodbmx,nopgmx,noprmx,noptmx,nordmx,nprpmx,nprpti,nprsmx,nprsti,nptkmx,nrct,nrctmx,nrk,nsbswt,nsk,nsrt,nsrtmx,ntitl1,ntitl2,ntitmx,ntrymx,nttkmx,nxmdmx,nxmod,nxopex,nxopmx,nxopt,nxpemx,nxrt,nxrtmx,o2maxi,o2mini,phmaxi,phmini,pressb,pressi,ptk,qgexsh,rkb,rxbari,sfcar,ssfcar,tempcb,tempci,tgexp,timmxi,tistti,tolbt,toldl,tolsat,tolxsf,trkb,ttk,ubmtbi,ubsri,ucxri,udac,uesri,uffg,ugerji,ugermo,ugersi,ugexj,ugexmo,ugexp,ugexr,uhfgex,uobsw,uprphi,uprspi,ureac,usbsw,utitl1,utitl2,uvfgex,uxcat,uxkgex,uxmod,uxopex,uxopt,uzveci,vreac,xgersi,xhfgex,ximaxi,xistti,xlkffg,xlkgex,xlkmod,xvfgex,zgexj,zvclgi)
        end if

        if (iopt(16) .eq. 0) then
            ! Switch to write on the other backup file the next time.
            nbkupn = nbkupb
        end if
    end if

    ! Initialize logical flags used in computing the reaction path.
    !   qzdump = .true. if the last value of Xi was a point
    !              corresponding to the PRS transfer interval.
    !   qdump  = .true. if a PRS transfer has just taken place.
    qzdump = .false.
    qzplot = .false.
    qskip = .false.
    qreq = .false.
    qreax = .false.
    qrpcfl = .false.
    qmod1 = .false.
    qmod2 = .false.
    qdmpr1 = .false.
    qdmpr2 = .false.

    qprntx = .false.
    qprntt = .false.
    qprnlx = .false.
    qprnlt = .false.

    qplotx = .false.
    qplott = .false.
    qplolx = .false.
    qplolt = .false.

    qprph0 = .false.
    qprph1 = .false.
    qpreh0 = .false.
    qpreh1 = .false.
    qpro20 = .false.
    qpro21 = .false.
    qpraw0 = .false.
    qpraw1 = .false.

    qplph0 = .false.
    qplph1 = .false.
    qpleh0 = .false.
    qpleh1 = .false.
    qplo20 = .false.
    qplo21 = .false.
    qplaw0 = .false.
    qplaw1 = .false.

    qhcon = .false.

    if (iopt(14) .eq. 2) then
        qhcon = .true.
    end if

    ! Set or initialize other variables used in computing the
    ! reaction path.
    kstpze = 0
    kstpmn = 0
    lprcin = tlg(prcinf)

    dlxtmx = prcinf
    xilim = prcinf

    scalim = 1.0
    scnsti = texp(sscrew(6)) - 1.
    scnstd = texp(-sscrew(6)) - 1.

    i = irang - 8
    fdx = real(i)
    xx = 100.
    fdx = min(fdx,xx)
    fdlim = texp(fdx)

    if (iodb(1) .ge. 1) then
        write (noutpt,1300) fdlim
        write (nttyo,1300) fdlim
1300 format(/' * Note - (EQ6/path) The magnitude of the finite',' difference functions',/7x,'is bounded by ',1pe10.3,'.',/)
    end if

    ! Re-set the arrays associated with finite-difference description
    ! of the number of mole of phases and species in the Equilibrium
    ! System (ES).
    ! Re-set the index arrays.
    call iiemop(iemop,iemos,iindx1,ipndx1,jsflag,kdim,kmax,ncmpe,ncmpr,noutpt,npet,npetmx,npt,nptmax,nset,nsetmx,nstmax,nttyo,uaqsln,uspec,uphase)

    ! Re-load the corresponding numbers of moles.
    do npe = 1,npet
        np = iemop(npe)
        emop(npe) = moph(np)
        nr1 = ncmpe(1,npe)
        nr2 = ncmpe(2,npe)

        do nse = nr1,nr2
            ns = iemos(nse)
            emos(nse) = mosp(ns)
        end do
    end do

    ! Initialize the akmat0 (F * B) and akmat1 (F * C) matrices. Each
    ! relates finite differences to the equivalent derivatives in a
    ! truncated Taylor's series. The akmat0 matrix is part of the
    ! predictor function. the akmat1 matrix is part of the corrector.
    nmax = nrd1mx*nrd1mx
    call initaz(akmat0,nmax)
    call initaz(akmat1,nmax)

    do i = 1,nrd1mx
        akmat0(i,i) = fctrl(i)
        akmat1(i,i) = fctrl(i)
    end do

    ! Set remaining step counters and order parameters.
    !   kordlm = the largest allowed value of kord (largest order
    !              considered for a given step), hence the largest
    !              allowed value of nord (the actual order used
    !              for a given step)
    !   ndelay = the number of points that must be done before
    !              allowing the order to build up from zero
    !   kly    = a counter used in restricting the growth of
    !               the step size; if it is >= 0, the scale limit
    !               (scalim) is set to a smaller value than it would be
    !               otherwise
    !   kstppr = counter of steps since the last print point
    !   kstppl = counter of steps since the last plot point
    kordlm = nordmx

    if (qecon) then
        kordlm = min(2,nordmx)
    end if

    if (qscon) then
        kordlm = 0
    end if

    ndelay = 1

    kstppr = 0
    kstppl = 0

    kly = 0
    delxi = dlxmx0

    naft1 = 0
    kaft1 = 0

    nsawth = 0
    jsawth = 0
    qsawth = .false.

    ! Calculate the first print, plot and dump points.
    call sippdp(actw,aw0plo,aw0prn,aw1plo,aw1prn,dlaplo,dlaprn,dleplo,dleprn,dlhplo,dlhprn,dloplo,dloprn,dltpll,dltplo,dltprl,dltprn,dlxdmp,dlxmx0,dlxpll,dlxplo,dlxprl,dlxprn,eh,eh0plo,eh0prn,eh1plo,eh1prn,eps100,fo2lg,lprcin,o20plo,o20prn,o21plo,o21prn,ph,ph0plo,ph0prn,ph1plo,ph1prn,prcinf,qredox,tiplol,tiplot,tiprnl,tiprnt,tistsv,xidump,xiplol,xiplot,xiprnl,xiprnt,xistsv)

    qtplo = dltplo.lt.prcinf .or. dltpll.lt.prcinf
    qtprn = dltprn.lt.prcinf .or. dltprl.lt.prcinf

    ! Make a new step in reaction progress.
100 continue
    kstep = kstep + 1
    kstppr = kstppr + 1
    kstppl = kstppl + 1
    kstpab = kstpab + 1
    avkdim = (avkdim*kstep + real(kdim))/real(kstep + 1)

    kly = kly - 1
    jsawth = jsawth - 1
    qmod2 = qmod1
    qmod1 = qmod .or. qreax
    qrapch = .false.

    aftm1 = aft0
    aft0 = aft1
    ncorr = 0
    dlxode = prcinf

    ! Check for instability (sawtoothing) in the ODE integration.
    if (iopt(2) .gt. 0) then
        rx0 = fdri0(1)
        rx1 = rirec1 - rirec0
        rxx = rx0*rx1

        if (rxx .lt. 0.) then
            nsawth = nsawth + 1

            if (nsawth .ge. 8) then
                ! Have sawtoothing. Take special action the next jsawth steps.
                jsawth = 20
                write (noutpt,'(" Sawtoothing detected")')
                write (nttyo,'(" Sawtoothing detected")')
            end if
        else
            nsawth = 0
        end if
    end if

    ! Pick the maximum allowed order for the current step (jordlm).
    if (jsawth .le. 0) then
        jordlm = kordlm
    else
        jordlm = min(jordlm,2)
    end if

    jordlm = min(kord,jordlm)

    ! Save information at the current point of reaction progress
    ! and update the finite differences for variables that are
    ! tracked using such.
    call stepfd(acflg,acflg0,affp0,affp,afrc0,afrc1,aw0,aw1,delxi,dxsm00,eh0,eh1,emop,emop0,emos,emos0,fdafm1,fdaf0,fdarm1,fdar0,fdawm1,fdaw0,fdehm1,fdeh0,fdlim,fdo2m1,fdo20,fdpem1,fdpe0,fdphm1,fdph0,fdrem1,fdre0,fdrim1,fdri0,fdrrm1,fdrr0,fdsem1,fdse0,fdzvm1,fdzv0,fje,fje0,fo2lg0,fo2lg1,fxi,fxi0,iemop,iemop0,iemos,iemos0,iindx0,iindx1,iodb,iopt,ipndx0,ipndx1,jcode,jreac,jreac0,jpflag,kdim,kdim0,kmax,km1,km10,kmt,kmt0,kord,kx1,kx10,kxt,kxt0,modr,modr0,moph,moph0,morr,morr0,mosp,mosp0,mtb,mtb0,nbt,nbtmax,ncmpe,ncmpe0,nodbmx,noptmx,nordmx,noutpt,npet,npetmx,npet0,npt,nptmax,npts,nrct,nrctmx,nrd1mx,nrndex,nset,nsetmx,nset0,nstmax,ph0,ph1,qredox,rirec0,rirec1,rreac0,rreac1,rrelr0,rrelr1,sfcar,sfcar0,sigmam,sigmm0,tempc,tempc0,time0,time1,uzvec0,uzvec1,xi0,xi1,xim1,xirct,xirct0,zvclg0,zvclg1,zvec0,zvec1)

    if (iopt(1) .eq. 2) then
        ! Doing the fluid-centered flow-through open system model.
        if (qdump) then
            ! Make a general partial shift of mass of the minerals in
            ! the ES to the PRS. This is part of the fluid-centered
            ! flow-through open system model. A partial shift of a phase
            ! leaves some of its mass in the ES.
            call pshfta(csts,emop,emop0,emos,emos0,fdpe0,fdpem1,fdse0,fdsem1,iemop,iemos,iern1,iern2,ietmax,iindx1,imrn1,imrn2,iodb,ipndx1,ixrn1,ixrn2,jcsort,jern1,jetmax,jgext,jpflag,jsflag,kbt,km1,kmax,kmt,kx1,kxt,loph,losp,moph,mosp,mprph,mprsp,mrgexs,mtb,mtb0,nbasp,nbaspd,nbt,nbtmax,ncmpe,ncmpr,netmax,ngext,nodbmx,nordmx,noutpt,npet,npetmx,npt,nptmax,nsetmx,nstmax,nsts,nstsmx,nstsr,nttyo,uaqsln,ufixf,uphase,uspec,xbar,xbarlg,zklgmn,zklogl,zvclg0,zvclg1,zvec0,zvec1)
        end if
    end if

    if (iopt(1) .eq. 2) then
        ! Doing the fluid-centered flow-through open system model.
        qbye = .false.

        if (qmod2 .or. kstep.eq.2) then
            ! On the previous step, a reactant became saturated or
            ! exhausted, a new phase was added to the ES, or the previous
            ! step was the first step of the run.
            if (kmt.ge.km1 .or. kxt.ge.kx1 .or. net.gt.0) then
                ! Make a total PRS shift of any mineral which is disappearing
                ! and purge it from the ES. The flag variable qbye is set
                ! to true if any such shift is made here.
                call dumpdp(csts,demop0,demos0,emop,emop0,emos,emos0,fdpe0,fdpem1,fdse0,fdsem1,iemop,iemos,iern1,iern2,ietmax,iindx0,iindx1,imrn1,imrn2,iodb,ipndx0,ipndx1,ixrn1,ixrn2,jcsort,jern1,jetmax,jgext,jpflag,jsflag,kbt,kdim,kdim0,kmax,km1,km10,kmt,kmt0,kord,kstep,kx1,kx10,kxt,kxt0,loph,losp,moph,mosp,mprph,mprsp,mrgexs,mtb,mtb0,nbasp,nbaspd,nbt,nbtmax,ncmpe,ncmpr,ndelay,netmax,ngext,nodbmx,nordmx,noutpt,npet,npetmx,npet0,npt,nptmax,npts,nset,nsetmx,nset0,nstmax,nsts,nstsmx,nstsr,nttyo,qbye,uaqsln,ufixf,uspec,uphase,uzvec0,uzvec1,xbar,xbarlg,zklgmn,zklogl,zvclg0,zvclg1,zvec0,zvec1)

                ! Note: if such a shift has been made here, the maximum
                ! order of the finite-differences (kord) has been set
                ! to zero.
                qmod2 = .false.
            end if
        end if
    end if

    ! Compute the step size and the order for the next increment. The
    ! order may be dropped if it appears that the high-order derivative
    ! estimates are not numerically significant. In economy mode, this
    ! algorithm is used to choose the order, but another algorithm (see
    ! below) is used to choose the step size.
110 continue
    if (qzdump) then
        delxi = max(dlxis2,delxi)
    end if

    qzprnt = .false.
    qprntx = .false.
    qprntt = .false.
    qprnlx= .false.
    qprnlt= .false.

    qzplot = .false.
    qplotx = .false.
    qplott = .false.
    qplolx = .false.
    qplolt = .false.

    qprph0 = .false.
    qprph1 = .false.
    qpreh0 = .false.
    qpreh1 = .false.
    qpro20 = .false.
    qpro21 = .false.
    qpraw0 = .false.
    qpraw1 = .false.

    qplph0 = .false.
    qplph1 = .false.
    qpleh0 = .false.
    qpleh1 = .false.
    qplo20 = .false.
    qplo21 = .false.
    qplaw0 = .false.
    qplaw1 = .false.

    qzdump = .false.
    qstabl = .false.
    qsawth = .false.

    if (qskip) then
        go to 120
    end if

    if (qshoot) then
        go to 125
    end if

    if (qecon) then
        ! Choose the step size and order when operating in economy
        ! mode.
        call chdxec(delxi,dlxmx0,dzvc0,iodb,kdim,kmax,km1,kxt,nodbmx,nord,noutpt,nrd1mx,qmin,qscon,scale,scalim,scnstd,scnsti,uzvec1,zvec0)
    else
        ! Choose the step size and order when operating in normal mode.
        ! Here control is based on the Gear approach. An order greater
        ! than zero is evaluated on the basis of the contributions
        ! of terms of one higher order. A hidden order is carried
        ! in the finite differences and derivatives of the z and r
        ! vectors to facilitate this approach.
        ! Get choices based on finite differences for the z vector.
        call chdxgz(delxi,dlxmx0,fdlim,fdzv0,iodb,iopt,jordlm,kdim,kmax,km1,kord,kxt,nodbmx,noptmx,nordmx,nordz,noutpt,nrd1mx,nsscmx,scalim,scfcz,sscrew,qmin,smp100,uzvec1,zklogu,zvec0,zvclg0)

        ! Get choices based on finite differences for the r vector.
        call chdxgr(delxi,dlxmx0,fdri0,fdrr0,iodb,jordlm,jreac,kord,nodbmx,nordmx,nordr,noutpt,nrct,nrctmx,nrd1mx,nsscmx,scalim,scfcr,sscrew,qriinf,rirec0,rrelr0,ureac)

        qstabz = .false.
        qstabr = .false.

        if (iopt(2).gt.0 .and. qsawth) then
            if (kord .ge. 2) then
                ! Compute derivatives for the z and r vectors from the
                ! corresponding finite differences.
                kordp1 = kord + 1

                ! First compute the matrix (akmat0 = F * B) used to
                ! calculate derivatives for truncated Taylor's series,
                ! given predictor-based finite differences.
                ! Calling sequence substitutions:
                !   kordp1 for nord
                call gakmat(akmat0,dxsm00,kordp1,nrd1mx)

                ! Now compute the corresponding derivatives.
                ! Calling sequence substitutions:
                !   kordp1 for nord
                call zderiv(akmat0,dzvc0,fdzv0,kdim,kmax,kordp1,nrd1mx)

                ! Calling sequence substitutions:
                !   kordp1 for nord
                call rderiv(akmat0,drer0,drir0,fdri0,fdrr0,jreac,kordp1,nrct,nrctmx,nrd1mx)

                ! Now compute the corresponding average derivatives.
                delxia = dxsm00(2)

                ! Calling sequence substitutions:
                !   kordp1 for nord
                call gsmdez(delxia,dzvc0,dzvc0s,kdim,kmax,kordp1,nrd1mx)

                ! Calling sequence substitutions:
                !   kordp1 for nord
                call gsmder(delxia,drer0,drer0s,drir0,drir0s,jreac,kordp1,nrct,nrctmx,nrd1mx)

                ! Choose step sizes and orders.
                call chdxtz(delxi,dlxmx0,dzvc0,iodb,iopt,jordlm,kdim,kmax,km1,kord,kxt,nodbmx,noptmx,nordmx,nordzs,noutpt,nrd1mx,nsscmx,scalim,scfczs,sscrew,qmin,smp100,uzvec1,zklogu,zvec0,zvclg0)

                call chdxtr(delxi,dlxmx0,drer0,drir0,iodb,jordlm,jreac,kord,nodbmx,nordmx,nordrs,noutpt,nrct,nrctmx,nrd1mx,nsscmx,scalim,scfcrs,sscrew,qriinf,rirec0,rrelr0,ureac)

                ! Resolve regular vs. "stable" results.
                if (scfczs .gt. scfcz) then
                    scfcz = scfczs
                    nordz = nordzs
                    qstabz = .true.
                end if

                if (scfcrs .gt. scfcr) then
                    scfcr = scfcrs
                    nordr = nordrs
                    qstabr = .true.
                end if
            end if
        end if

        ! Resolve z vector vs. r vector results.
        scale = scfcz
        nord = nordz
        qstabl = qstabz

        if (scfcr .lt. scfcz) then
            scale = scfcr
            nord = nordr
            qstabl = qstabr
        end if
    end if

    ! Calculate derivatives from finite differences. These are
    ! consistent with the actual order to be used (nord).
    if (nord .gt. 0) then
        ! First compute the matrix (akmat0 = F * B) used to calculate
        ! derivatives for truncated Taylor's series, given predictor-
        ! based finite differences.
        call gakmat(akmat0,dxsm00,nord,nrd1mx)

        ! Now compute the derivatives.
        call zderiv(akmat0,dzvc0,fdzv0,kdim,kmax,nord,nrd1mx)

        if (iopt(2) .gt. 0) then
            call rderiv(akmat0,drer0,drir0,fdri0,fdrr0,jreac,nord,nrct,nrctmx,nrd1mx)
        end if

        call aderiv(akmat0,daffp0,fdaf0,nord,nordmx,npt,nptmax,nrd1mx)
        call bderiv(akmat0,dafrc0,fdar0,jreac,nord,nordmx,nrct,nrctmx,nrd1mx)
        call pderiv(akmat0,demop0,fdpe0,nord,nordmx,npet,npetmx,nrd1mx)
        call sderiv(akmat0,demos0,fdse0,nord,nordmx,nrd1mx,nset,nsetmx)

        ! Calling sequence substitutions:
        !   dph0 for dxx0
        !   fdph0 for fdxx0
        call xderiv(akmat0,dph0,fdph0,nord,nordmx,nrd1mx)

        if (qredox) then
            ! Calling sequence substitutions:
            !   deh0 for dxx0
            !   fdeh0 for fdxx0
            call xderiv(akmat0,deh0,fdeh0,nord,nordmx,nrd1mx)

            ! Calling sequence substitutions:
            !   do20 for dxx0
            !   fdo20 for fdxx0
            call xderiv(akmat0,do20,fdo20,nord,nordmx,nrd1mx)
        end if

        ! Calling sequence substitutions:
        !   daw0 for dxx0
        !   fdaw0 for fdxx0
        call xderiv(akmat0,daw0,fdaw0,nord,nordmx,nrd1mx)

        if (qstabl) then
            delxia = dxsm00(nord)

            ! Compute average derivatives.
            call gsmdez(delxia,dzvc0,dzvc0s,kdim,kmax,nord,nrd1mx)
            call gsmder(delxia,drer0,drer0s,drir0,drir0s,jreac,nord,nrct,nrctmx,nrd1mx)

            ! Replace the actual with the average derivatives.
            do kcol = 1,kdim
                do n = 1,nord
                    dzvc0(n,kcol) = dzvc0s(n,kcol)
                end do
            end do

            do n = 1,nord
                drir0(n) = drir0s(n)
            end do

            do nrc = 1,nrct
                do n = 1,nord
                    drer0(n,nrc) = drer0s(n,nrc)
                end do
            end do
        end if
    end if

    ! Compute the new value of delxi.
    delxi = delxi*scale

    ! If the step size was previously cut because a phase boundary
    ! was overstepped by an unacceptable amount, limit delxi so that
    ! the new value of Xi will be less than the value at the last
    ! known point of overstep.
    dlxlim = 0.75*(xilim - xi0)

    if (dlxlim .lt. dlxmin) then
        ! The calculation has has now gotten very close to the previously
        ! determined overstep location. After the present step is
        ! completed, the calculation will be within a minimum step size
        ! distance of that location. Failure to detect either the
        ! phase boundary or the overstep at this point means that the
        ! overstep location previously calculated is not actually valid.
        ! What this means is that as the calculation gets closer to the
        ! phase boundary, its calculated position is moving to a greater
        ! value of reaction progress. This can occur when the mass
        ! transfer associated with a step must be calculated from an
        ! approximate, not an exact integration. Such is the case for
        ! kinetic mode calculations. When the calculation goes back
        ! and cuts the step size, the resulting calculation of the
        ! position of the phase boundary then becomes more accurate.
        dlxlim = dlxmin
        xilim = prcinf
    end if

    delxi = min(delxi,dlxlim)

    ! Make sure that delxi is not less than the minimum normally
    ! allowed value.
    delxi = max(delxi,dlxmin)

    if (iodb(5) .ge. 1) then
        write (noutpt,1320) nord,delxi
    end if

1320 format(/' --- Selection: nord= ',i2,', delxi= ',1pe11.4,' ---')

120 continue
    qskip = .false.

    ! Set other simple limits on the step size.
    if (iopt(1) .eq. 2) then
        dlxis2 = delxi
        dxdmp = xidump - xi0

        if (delxi .gt. dxdmp) then
            dxdmp = max(dxdmp,dlxmin)
            delxi = dxdmp

            if (iodb(5) .ge. 2) then
                write (noutpt,1330) delxi
            end if

1330 format(' --- delxi cut to ',1pe11.4,' to match the Xi PRS',' transfer interval ---')
        end if
    end if

    dlxisv = delxi
    dx1 = xiprnt - xi0
    dx2 = xiprnl - xi0
    dlxipr = min(dx1,dx2)

    if (delxi .gt. dlxipr) then
        delxi = dlxipr

        if (iodb(5) .ge. 2) then
            write (noutpt,1340) delxi
        end if

1340 format(' --- delxi cut to ',1pe11.4,' to match the Xi print',' interval ---')
    end if

    dx1 = xiplot - xi0
    dx2 = xiplol - xi0
    dlxipl = min(dx1,dx2)

    if (delxi .gt. dlxipl) then
        delxi = dlxipl

        if (iodb(5) .ge. 2) then
            write (noutpt,1350) delxi
        end if

1350 format(' --- delxi cut to ',1pe11.4,' to match the Xi plot',' interval ---')
    end if

    dlximx = ximax - xi0

    if (dlximx .lt. delxi) then
        delxi = dlximx

        if (iodb(5) .ge. 2) then
            write (noutpt,1360) delxi
        end if

1360 format(' --- delxi cut to  ',1pe11.4,' to match the maximum',' value of Xi  ---')
    end if

    dlxilm = dlxmax

    if (nord .le. 0) then
        dlxilm = dlxmx0
    end if

    if (dlxilm .lt. delxi) then
        delxi = dlxilm

        if (iodb(5) .ge. 2) then
            write (noutpt,1370) delxi
        end if

1370 format(' --- delxi cut to ',1pe11.4,' to match the upper limit',' for any step ---')
    end if

125 continue

    if (nord.gt.0 .and. .not.qscon) then
        ! Limit delxi when the variables corresponding to certain basis
        ! species are changing rapidly. Such variables may include the pH
        ! or the log fO2. This limit allows some minimum of information
        ! density to be obtained in the neighborhood of the rapid
        ! change.
        call ldlxrc(al10,delxi,dlxmin,dzvc0,iodb,iindx1,kbt,kdim,kelect,khydr,khydx,km1,kmax,ko2gaq,krdxsp,kwater,kxt,nbasp,nbtmax,nodbmx,nord,noutpt,nrd1mx,nstmax,nttyo,qrapch,uspec,zklogu,zvclg0,zvclg1,zvec0,zvec1)
    end if

    if (nord.gt.0 .and. .not.qscon .and. iopt(3).lt.1  .and. kxt.ge.km1 .and. iopt(1).le.1) then
        ! Limit delxi when a phase is rapidly disapparing from the
        ! equilibrium system. This is mechanism is only designed to
        ! help other mechanisms efficiently locate a phase disappearance
        ! boundary. It is not designed to increase information density
        ! near the boundary.
        call ldlxrd(delxi,dlxmin,dzvc0,fdzv0,iodb,ipndx1,kdim,km1,kmax,kxt,loph,nodbmx,nord,noutpt,nptmax,nrd1mx,nttyo,uphase,zklogu,zvec0)
    end if

    qdump = .false.

    if (iopt(1).eq.2 .and. nord.ge.1 .and. npet.gt.0) then
        ! Limit delxi by approximate position of significant maxima
        ! of mineral masses. Because the cost of restricting delxi to
        ! keep high accuracy in the estimates of first derivatives is
        ! prohibitive, control can not be based strictly on derivatives,
        ! but must be based mainly on how much of a mineral mass would
        ! be destroyed on a given step. This permits avoiding maxima
        ! which are false (due to Taylor's series inaccuracy) or trivial
        ! (if the mass of a mineral is very small). If the Taylor's
        ! series fail to detect a phase whose mass is going over a maximum
        ! and iopt(1) = 2, then go back instructions in EQ6/eqcalc.f will
        ! prevent a skip over the maximum.
        call fpbflo(al10,delxi,demop0,dlxmin,dxval0,d1emp1,d2emp1,emop,emop0,eps100,fdpe0,iemop,ier,iodb,nodbmx,nord,nordmx,noutpt,npet,npetmx,nptmax,nrd1mx,nttyo,qdump,toldl,uaqsln,ufixf,uphase,xim1,xi0,xi1,xval0,zklogu)

        if (ier .gt. 0) then
            ! A search for a maximum in subroutine fpbflo failed. Go back
            ! to the base point, shift mass to the PRS, and drop back
            ! to order zero.
            qdump = .true.
            kord = 0
            nord = 0
            npts = 1
            delxi = 0.
            xi1 = xi0

            do npe = 1,npet
                emop(npe) = emop0(npe)
            end do

            do nse = 1,nset
                emos(nse) = emos0(nse)
            end do

            do kcol = 1,kxt
                zvclg1(kcol) = zvclg0(kcol)
                zvec1(kcol) = zvec0(kcol)
            end do

            call pshfta(csts,emop,emop0,emos,emos0,fdpe0,fdpem1,fdse0,fdsem1,iemop,iemos,iern1,iern2,ietmax,iindx1,imrn1,imrn2,iodb,ipndx1,ixrn1,ixrn2,jcsort,jern1,jetmax,jgext,jpflag,jsflag,kbt,km1,kmax,kmt,kx1,kxt,loph,losp,moph,mosp,mprph,mprsp,mrgexs,mtb,mtb0,nbasp,nbaspd,nbt,nbtmax,ncmpe,ncmpr,netmax,ngext,nodbmx,nordmx,noutpt,npet,npetmx,npt,nptmax,nsetmx,nstmax,nsts,nstsmx,nstsr,nttyo,uaqsln,ufixf,uphase,uspec,xbar,xbarlg,zklgmn,zklogl,zvclg0,zvclg1,zvec0,zvec1)

            delxi = dlxmx0
            go to 110
        end if
    end if

    if (nord.gt.0 .and. .not.qscon) then
        ! Find the phase boundary for a newly appearing phase. This is
        ! here defined as a point at which the affinity to precipitate
        ! has a target value (aftarg) that is slightly greater than zero.
        ! If the Taylor's series fail to detect the phase boundary at
        ! which a new phase appears, this condition will be subsequently
        ! detected by a test in EQLIB/newton.f. This subroutine will then
        ! issue a go back instruction, which will prevent the boundary
        ! from being skipped over.
        call fpbnpp(affp,affp0,aftarg,daffp0,delxi,dlxmin,dxval0,eps100,iodb,iopt,jpflag,nodbmx,noptmx,nord,nordmx,noutpt,npchk,npt,nptmax,nrd1mx,nttyo,tolaft,tolsat,uphase,xi0,xi1,xval0)

        if (iopt(1).ne.2 .and. npet.gt.0) then
            ! Find the phase boundary for a disappearing phase in the ES.
            call fpbdpp(delxi,demop0,dlxmin,dxval0,emop,emop0,eps100,iemop,iodb,iopt,nodbmx,noptmx,nord,nordmx,noutpt,npet,nrd1mx,npetmx,nptmax,nttyo,uphase,xi0,xi1,xval0)
        end if
    end if

    if (nrct.gt.0 .and. nord.gt.0) then
        if (aft0 .gt. 0.01) then
            ! Check the signs of the predicted reactant affinities. If
            ! necessary, reduce delxi so that none of the affinities
            ! change sign. For the purposes of determining a crossover,
            ! a reactant affinity is taken to be zero if its magnitude is
            ! less than or equal to tolsar.
            ! This check is made to help deal with rate expressions that
            ! do not approach zero as the corresponding equilibria are
            ! approached. For example, one might specify a constant rate
            ! of dissolution that applies to any value of the corres-
            ! ponding affinity as long as it implies undersaturation. It
            ! is important not to extrapolate such a rate to any other
            ! condition (in this example, saturation or supersaturation).
            ! If close to overall equilibrium, this check is best avoided.
            ! Here the rates of "badly behaving" reactions may tend to
            ! oscillate about their equilibrium points. As long as the
            ! magnitudes of the oscillations are small, it is better to
            ! not try to pin down the cross-overs.
            call chksar(afrc0,afrcp,dafrc0,delxi,dlxmin,dxval0,eps100,iodb,jreac,nodbmx,noutpt,nord,nordmx,nrct,nrctmx,nrd1mx,nrk,nttyo,tolsar,ureac,xi0,xi1,xval0)
        end if
    end if

    if (iopt(2) .gt. 0) then
        if (nrct.gt.0 .and. nord.gt.0) then
            if (aft0 .gt. 0.01) then
                ! Check the sign of the predicted inverse rate. If necessary,
                ! reduce delxi so that the inverse rate does not become
                ! less than eps100. Physically, it makes no sense that
                ! the inverse rate should approach zero or become negative,
                ! but errors in finite differences that could cause a
                ! predicted inverse rate to do this. If so, the step size
                ! needs to be reduced.
                call chksir(delxi,dlxmin,drir0,dxval0,eps100,iodb,nodbmx,nord,noutpt,nrd1mx,nttyo,rirec0,rirecp,xi0,xi1,xval0)
            end if
        end if
    end if

    if (nrct.gt.0 .and. nord.gt.0) then
        if (aft0 .gt. 0.01) then
            ! Check the signs of the predicted relative rates. If,
            ! necessary reduce delxi so that none of the relative rates
            ! change sign. For the purposes of determining a crossover,
            ! a relative rate is taken to be zero if its magnitude is
            ! less than tolsrr.
            ! This check is best avoided near overall equilibrium. See
            ! comments made above in regard to a similar check on sign
            ! changes of reactant affinities.
            call chksrr(delxi,dlxmin,drer0,dxval0,eps100,iodb,jreac,nodbmx,noutpt,nord,nrct,nrctmx,nrd1mx,nttyo,rrelr0,rrelrp,tolsrr,ureac,xi0,xi1,xval0)
        end if
    end if

    if (iopt(2) .gt. 0) then
        if (nrct .gt. 0) then
            ! Make absolutely sure that the time increment is positive.
            ! Errors in numerical integration have the potential to
            ! produce a negative or zero result.
            call chksti(akmat0,drer0,drir0,deltim,delxi,dlxmin,fdri0,fdrr0,iodb,jreac,kly,kmax,nodbmx,kord,nord,noutpt,npts,nrct,nrctmx,nrd1mx,nttyo,prcinf,qriinf,rirec0,smp100,time0,time1,xi0,xi1)
        end if
    end if

    if (nrct .gt. 0) then
        ! Limit delxi by exhaustion of reactants.
        call fpexrc(delxi,dlxmin,drer0,dxval0,eps100,iodb,jreac,morr,morr0,nodbmx,nord,noutpt,nrct,nrctmx,nrd1mx,nttyo,qdump,rrelr0,ureac,xi0,xi1,xval0)
    end if

    if (iopt(2) .gt. 0) then
        ! Limit delxi so that the next time-based print point is not
        ! exceeded.
        if (qtprn) then
            call chktpr(delxi,dlxmin,dlxtpr,drir0,dxval0,eps100,iodb,nodbmx,nord,noutpt,nrd1mx,nttyo,prcinf,qdump,qriinf,rirec0,tiprnl,tiprnt,time0,time1,tolxst,xi0,xi1,xval0)
        end if
    end if

    ! Limit delxi so that the next pH-based print point is not exceeded.
    if (dlhprn .lt. prcinf) then
        call ckphpr(delxi,dlxmin,dph0,dxh0pr,dxh1pr,dxval0,eps100,iodb,nodbmx,nord,noutpt,nrd1mx,nttyo,ph0,ph1,ph0prn,ph1prn,prcinf,qdump,tolxsu,xi0,xi1,xval0)
    end if

    if (qredox) then
        ! Limit delxi so that the next Eh-based print point is not
        ! exceeded.
        if (dlhprn .lt. prcinf) then
            call ckehpr(delxi,dlxmin,deh0,dxe0pr,dxe1pr,dxval0,eps100,iodb,nodbmx,nord,noutpt,nrd1mx,nttyo,eh0,eh1,eh0prn,eh1prn,prcinf,qdump,tolxsu,xi0,xi1,xval0)
        end if

        ! Limit delxi so that the next log fO2-based print point is not
        ! exceeded.
        if (dloprn .lt. prcinf) then
            call cko2pr(delxi,dlxmin,do20,dxo0pr,dxo1pr,dxval0,eps100,iodb,nodbmx,nord,noutpt,nrd1mx,nttyo,fo2lg0,fo2lg1,o20prn,o21prn,prcinf,qdump,tolxsu,xi0,xi1,xval0)
        end if
    end if

    ! Limit delxi so that the next aw-based print point is not exceeded.
    if (dlaprn .lt. prcinf) then
        call ckawpr(delxi,dlxmin,daw0,dxw0pr,dxw1pr,dxval0,eps100,iodb,nodbmx,nord,noutpt,nrd1mx,nttyo,aw0,aw1,aw0prn,aw1prn,prcinf,qdump,tolxsu,xi0,xi1,xval0)
    end if

    if (iopt(2) .gt. 0) then
        ! Limit delxi so that the next time-based plot point is not
        ! exceeded.
        if (qtplo) then
            call chktpl(delxi,dlxmin,dlxtpl,drir0,dxval0,eps100,iodb,nodbmx,nord,noutpt,nrd1mx,nttyo,prcinf,qdump,qriinf,rirec0,tiplol,tiplot,time0,time1,tolxst,xi0,xi1,xval0)
        end if
    end if

    ! Limit delxi so that the next pH-based plot point is not exceeded.
    if (dlhplo .lt. prcinf) then
        call ckphpl(delxi,dlxmin,dph0,dxh0pl,dxh1pl,dxval0,eps100,iodb,nodbmx,nord,noutpt,nrd1mx,nttyo,ph0,ph1,ph0plo,ph1plo,prcinf,qdump,tolxsu,xi0,xi1,xval0)
    end if

    if (qredox) then
        ! Limit delxi so that the next Eh-based plot point is not
        ! exceeded.
        if (dlhplo .lt. prcinf) then
            call ckehpl(delxi,dlxmin,deh0,dxe0pl,dxe1pl,dxval0,eps100,iodb,nodbmx,nord,noutpt,nrd1mx,nttyo,eh0,eh1,eh0plo,eh1plo,prcinf,qdump,tolxsu,xi0,xi1,xval0)
        end if

        ! Limit delxi so that the next log fO2-based plot point is not
        ! exceeded.
        if (dloplo .lt. prcinf) then
            call cko2pl(delxi,dlxmin,do20,dxo0pl,dxo1pl,dxval0,eps100,iodb,nodbmx,nord,noutpt,nrd1mx,nttyo,fo2lg0,fo2lg1,o20plo,o21plo,prcinf,qdump,tolxsu,xi0,xi1,xval0)
        end if
    end if

    ! Limit delxi so that the next aw-based plot point is not exceeded.
    if (dlaplo .lt. prcinf) then
        call ckawpl(delxi,dlxmin,daw0,dxw0pl,dxw1pl,dxval0,eps100,iodb,nodbmx,nord,noutpt,nrd1mx,nttyo,aw0,aw1,aw0plo,aw1plo,prcinf,qdump,tolxsu,xi0,xi1,xval0)
    end if

    if (iopt(2) .gt. 0) then
        ! Limit delxi so that the requested maximum value of time
        ! is not exceeded.
        if (timemx .lt. prcinf) then
            call chktmx(delxi,dlxmin,dlxtmx,drir0,dxval0,eps100,iodb,nodbmx,nord,noutpt,nrd1mx,nttyo,prcinf,qdump,qriinf,rirec0,timemx,time0,time1,tolxst,xi0,xi1,xval0)
        end if
    end if

    ! Limit delxi so that the requested minimum value of pH
    ! is not exceeded.
    if (phmin .gt. -prcinf) then
        call ckphmn(delxi,dlxmin,dph0,dxh0mx,dxval0,eps100,iodb,nodbmx,nord,noutpt,nrd1mx,nttyo,ph0,ph1,phmin,prcinf,qdump,tolxsu,xi0,xi1,xval0)
    end if

    ! Limit delxi so that the requested maximum value of pH
    ! is not exceeded.
    if (phmax .lt. prcinf) then
        call ckphmx(delxi,dlxmin,dph0,dxh1mx,dxval0,eps100,iodb,nodbmx,nord,noutpt,nrd1mx,nttyo,ph0,ph1,phmax,prcinf,qdump,tolxsu,xi0,xi1,xval0)
    end if

    if (qredox) then
        ! Limit delxi so that the requested minimum value of Eh
        ! is not exceeded.
        if (ehmin .gt. -prcinf) then
            call ckehmn(delxi,dlxmin,deh0,dxe0mx,dxval0,eps100,iodb,nodbmx,nord,noutpt,nrd1mx,nttyo,eh0,eh1,ehmin,prcinf,qdump,tolxsu,xi0,xi1,xval0)
        end if

        ! Limit delxi so that the requested maximum value of Eh
        ! is not exceeded.
        if (ehmax .lt. prcinf) then
            call ckehmx(delxi,dlxmin,deh0,dxe1mx,dxval0,eps100,iodb,nodbmx,nord,noutpt,nrd1mx,nttyo,eh0,eh1,ehmax,prcinf,qdump,tolxsu,xi0,xi1,xval0)
        end if

        ! Limit delxi so that the requested minimum value of log fO2
        ! is not exceeded.
        if (o2min .gt. -prcinf) then
            call cko2mn(delxi,dlxmin,do20,dxo0mx,dxval0,eps100,iodb,nodbmx,nord,noutpt,nrd1mx,nttyo,fo2lg0,fo2lg1,o2min,prcinf,qdump,tolxsu,xi0,xi1,xval0)
        end if

        ! Limit delxi so that the requested maximum value of log fO2
        ! is not exceeded.
        if (o2max .lt. prcinf) then
            call cko2mx(delxi,dlxmin,do20,dxo1mx,dxval0,eps100,iodb,nodbmx,nord,noutpt,nrd1mx,nttyo,fo2lg0,fo2lg1,o2max,prcinf,qdump,tolxsu,xi0,xi1,xval0)
        end if
    end if

    ! Limit delxi so that the requested minimum value of the activity
    ! of water is not exceeded.
    if (awmin .gt. -prcinf) then
        call ckawmn(delxi,dlxmin,daw0,dxw0mx,dxval0,eps100,iodb,nodbmx,nord,noutpt,nrd1mx,nttyo,aw0,aw1,awmin,prcinf,qdump,tolxsu,xi0,xi1,xval0)
    end if

    ! Limit delxi so that the requested maximum value of the activity
    ! of water is not exceeded.
    if (awmax .lt. prcinf) then
        call ckawmx(delxi,dlxmin,daw0,dxw1mx,dxval0,eps100,iodb,nodbmx,nord,noutpt,nrd1mx,nttyo,aw0,aw1,awmax,prcinf,qdump,tolxsu,xi0,xi1,xval0)
    end if

    ! Make sure that delxi is not less than dlxmin unless one of the
    ! following special conditions is satisfied.
    if (abs(delxi - dlxipr).gt.0. .and. abs(delxi - dlxipl).gt.0. .and. abs(delxi - dlximx).gt.0. .and. abs(delxi - dlxtmx).gt.0.) then
        delxi = max(delxi,dlxmin)
    end if

    ncut = 0
    jcut = 0

130 continue
    xi1 = xi0 + delxi

    if (.not.qshoot) then
        ! The following trap keeps the code from doing lengthy cycling
        ! with delxi equal to zero.
        if (delxi .le. 0.) then
            kstpze = kstpze + 1
            qx = kstpze .ge. 2

            if (qx) then
                write (noutpt,1380) kstpze
                write (nttyo,1380) kstpze
1380 format(/' * Error - (EQ6/path) The step size has been zero',' now ',i2,' times in a row.')

                stop
            end if
        else
            kstpze = 0
        end if

        ! The following trap keeps the code from doing lengthy cycling
        ! with  delxi equal to dlxmin.
        if (delxi .le. dlxmin) then
            kstpmn = kstpmn + 1
            qx = kstpmn .ge. 100

            if (qx) then
                write (noutpt,1390) kstpmn
                write (nttyo,1390) kstpmn
1390 format(/' * Error - (EQ6/path) The step size has not',' exceeded the minimum',/7x,'value now ',i3,' times in a row.')

                stop
            end if
        else
            kstpmn = 0
        end if
    end if

    ! The following is a return point for doing an ODE corrector
    ! iteration.
140 continue

    if (nrct .gt. 0) then
        ! Increment the irreversible reactions (those associated with
        ! the "reactants"); update the ES mass balance totals accordingly.
        call reacts(cbsr,csts,delxi,drer0,iern1,ietmax,iktmax,iodb,jcode,jetmax,jgext,jreac,modr,modr0,morr,morr0,mrgers,mtb,mtb0,nbaspd,nbt,nbtmax,nbt1mx,ncmpr,nern1,nern2,nertmx,netmax,ngext,nodbmx,nord,noutpt,nptmax,nrct,nrctmx,nrd1mx,nrndex,nsrtmx,nstmax,nsts,nstsmx,nstsr,nttyo,nxridx,nxrtmx,rrelr0,rxbar,ureac,xirct,xirct0)
    end if

    if (iopt(2) .gt. 0) then
        if (qshoot) then
            time1 = prcinf
            deltim = prcinf
        else
            ! Increment the time and related functions.
            call timeca(deltim,delxi,drir0,iodb,nodbmx,nord,noutpt,nrd1mx,nttyo,prcinf,qriinf,rirec0,time0,time1)
        end if

        if (delxi .le. dlxmin) then
            ! Make sure that the calculated time does not exceed any
            ! any specified limits such as the maximum time just because
            ! delxi is at the minimum value.
            call tivchk(deltim,delxi,qtvchk,time1,time0,timemx,tiplol,tiplot,tiprnl,tiprnt,tolxst)
        end if
    end if

    if (.not.qcntmp .or. .not.qcnpre) then
        ! Recompute the temperature and pressure. Then recompute the
        ! thermodynamic and kinetic quantities which depend these
        ! variables.
        ntpr0 = ntpr
        call tpadv(aadh,aadhh,aadhv,aaphi,abdh,abdhh,abdhv,abdoth,abdot,abdotv,adadhh,adadhv,adbdhh,adbdhv,adbdth,adbdtv,adh,adhfe,adhh,adhv,adhfs,adhfsd,advfe,advfs,advfsd,afcnst,al10,amu,aslm,aphi,aprehw,apresg,apresh,apx,avcnst,axhfe,axhfs,axhfsd,axlke,axlks,axlksd,axvfe,axvfs,axvfsd,bdh,bdhh,bdhv,bdot,bdoth,bdotv,dadhh,dadhv,dbdhh,dbdhv,dbdth,dbdtv,dhfe,dhfs,dhfsd,dvfe,dvfs,dvfsd,eact,ehfac,farad,hact,iact,iapxmx,iktmax,imchmx,imech,iopg,iopt,ipbtmx,ipch,ipchmx,ipcv,ipcvmx,ixrn1,ixrn2,jpfcmx,jpress,jptffl,jsol,jtemp,narxmx,narxt,narxth,nbasp,nbaspd,nbt,nbtd,nbtmax,ncmpr,ndrsr,ndrsrd,nmut,nmutmx,nopgmx,noptmx,noutpt,nptkmx,nptmax,nrct,nrctmx,nrk,nslt,nsltmx,nst,nstmax,ntpr,ntprmx,ntprt,nttkmx,nttyo,nweope,nwndpc,nxt,nxtmax,pmu,presg,presh,press,pressb,pressd,pslamn,ptk,rcnstv,rconst,rk,rkb,rtcnst,tempc,tempcb,tempcd,tempcu,tempk,time1,trkb,ttk,uphase,uspec,wfac,xhfe,xhfs,xhfsd,xi1,xlke,xlks,xlksd,xvfe,xvfs,xvfsd)

        qtrch = ntpr0 .ne. ntpr

        if (qtrch) then
            nord = 0
        end if
    end if

    if (delxi .le. dlxmin) then
        nord = 0
        kord = 0
    end if

    ! Make a Taylor's series expansion of the z vector, applying
    ! change limits. It is important to make a protected expansion
    ! here to clear out unrecoverable values in the zvclg1 and
    ! zvec1 vectors that may have been generated in previous
    ! order/step size adjustments for the current step. This is
    ! true even if the current order is zero, as the order may
    ! have been higher earlier in the current step.
    qztayl = .true.
    call ztaylr(delxi,dzvc0,kdim,kmax,km1,kxt,nord,nrd1mx,qztayl,zklogu,zvclg0,zvclg1,zvec0,zvec1)

    qstart = xi1 .le. xistsv

    ! Make the equilibrium calculation at the current point of reaction
    ! progress. Subroutine eqshel forms a shell around EQ6/eqphas.f.
    ! It may make small advances in the value of the reaction progress
    ! variable in order to step over a small region (commonly on the
    ! order of 1 x 10-8 mol) in which the equilibrium state can not be
    ! calculated satisfying the normal tolerances. One category of such
    ! regions is associated with phase boundaries (where a phase such
    ! as a mineral either appears or disappears). The problem is that
    ! a calculation assuming the phase is present crashes, indicating
    ! that the phase should be deleted. One the other hand, a
    ! calculation assuming it is not present converges, but gives a
    ! supersaturation in excess of the supersaturation tolerance (a
    ! condition which normally instructs the code to add the phase to
    ! the phase assemblage and try again). Another category of such
    ! regions is associated with a rapid change in the value of a
    ! variable associated with an aqueous basis species. Typically
    ! this is the oxygen fugacity or equivalent redox variable. Such
    ! a change is called a redox jump. A similar jump might occur for
    ! some non-redox variable, such as pH. However, this is not likely
    ! for any non-redox variable because of the general presence of
    ! significant buffering effects of one kind or another (true
    ! buffering, mass buffering).
    call eqshel(aadh,aadhh,aadhv,aamatr,aaphi,abar,abdh,abdhh,abdhv,abdot,abdoth,abdotv,acflg,acflgo,act,actlg,adadhh,adadhv,adbdhh,adbdhv,adbdth,adbdtv,adh,adhh,adhv,afcnst,affp,affs,alpha,al10,amtb,aphi,apx,avcnst,azero,a3bar,a3bars,bacfmx,bbig,bdh,bdhh,bdhv,bdot,bdoth,bdotv,beta,betamx,betao,bgamx,bneg,bpx,cbsr,cco2,cdac,cegexs,cesr,cess,cdrs,cdrsd,cdrsx,cdrtw,cdrw,cjbasp,cnufac,conc,conclg,cpgexs,cscale,csigma,csts,dadhh,dadhv,dbdhh,dbdhv,dbdth,dbdtv,deltim,delvco,delvec,delxi,dlogxw,dlxmin,drer0,drir0,dzvc0,d1zvc1,eact,egers,egexjc,egexjf,egexs,eh,ehfac,elecsr,eps100,farad,fje,fjeo,fkrc,fo2,fo2lg,fsort,fugac,fugalg,fxi,fxio,gmmatr,hact,iact,iapxt,ibpxt,ibswx,ielam,ier,iern1,iern2,ifcphi1,ifcphi2,ifnnn,ifn2n,ifpsi1,ifpsi2,ifrn1,ifrn2,ifzeta,igas,igstak,iindx0,iindx1,ilcphi1,ilcphi2,ilnnn,iln2n,ilpsi1,ilpsi2,ilrn1,ilrn2,ilzeta,imech,imrn1,imrn2,insgf,iodb,iopg,iopt,ipch,ipivot,ipndx1,ipcv,istack,iter,itermx,ixbasp,ixrn1,ixrn2,izmax,jcode,jcsort,jflag,jgsort,jgstak,jjsort,jpflag,jpress,jptffl,jreac,jsflag,jsitex,jsol,jssort,jstack,jtemp,kbt,kction,kdim,kelect,khydr,khydx,km10,km1,kmt,kmt0,ko2gaq,kpsat,kpsst,krdxsp,kwater,kx1,kx10,kxt,kxt0,loph,losp,lsort,modr,modr0,moph,morr,morr0,mosp,mrgers,mrgexs,mtb,mtbaq,mtb0,mte,mteaq,mwtrc,narn1,narn2,narxt,nat,nbasp,nbaspd,nbaspx,nbt,nbtd,nbw,nchlor,ncmpr,ncorr,nct,ndac,ndact,ndrs,ndrsd,ndrsx,ndrsr,ndrsrd,ndrsrx,nelect,nern1,nern2,ness,nessr,net,nfrn1,nfrn2,ngrn1,ngrn2,ngt,nhydr,nhydx,nlrn1,nlrn2,nlt,nmrn1,nmrn2,nmt,nord,noutpt,no2gaq,npchk,nphasx,npslmx,npt,nrct,nrdxsp,nrk,nrndex,nsk,nsslmx,nst,nsts,nstsr,ntpr,ntrymx,ntprt,nttyo,nxridx,nxrn1,nxrn2,nxt,nweope,nwndpc,omega,omeglg,prcinf,press,pressb,pressd,ptk,qbassw,qbseqc,qbye,qcnpre,qcntmp,qhawep,qmod,qoptmz,qpit75,qredox,qriinf,qscon,qshoot,qstart,qtrch,qtvchk,qxknph,q6mode,rcnstv,rconst,rkb,rhsvec,rirec0,rk,rreacn,rreac1,rrelr0,rrelr1,rtcnst,rxbar,screwd,sidrph,sidrsp,sigmam,sigmmo,smp100,tdays,tempc,tempcb,tempcd,tempcu,tempc0,tempk,timemx,time0,time1,tiplol,tiplot,tiprnl,tiprnt,tistsv,tolbt,toldl,tolsat,tolsst,tolxst,trkb,ttk,ubacmx,ubgamx,udac,ufixf,ugermo,ulbeta,uldel,uphase,ureac,uspec,uzvec0,uzvec1,vreac,weight,wfac,wodr,worr,xbar,xbarlg,xbarw,xbarwc,xbrwlc,xbrwlg,xgers,xirct,xirct0,xistsv,xi0,xi1,zchar,zchcu6,zchsq2,zklogu,zvclg0,zvclg1,zvec0,zvec1)

    ! Note on ier codes returned by EQ6/eqshel.f:
    !   =    0  Okay
    !   =   10  Go back and take a smaller step size to avoid exceeding
    !             the supersaturation tolerance (tolsst)
    !   =  170  Too much of a phase was destroyed under the flow-through
    !             open system model; go back and first move part of the
    !             mass of protected phases in the ES to the PRS
    !   =  180  One of a number of problems occurred which may be
    !             resolvable, at least partially, by going back and
    !             cutting the step size
    !   =  190  Need to slide over a region of computational
    !             instability, but sliding is inhibited; go back, but
    !             terminate work on the current problem
    if (ier .le. 0) then
        go to 200
    end if

    if (ier .eq. 10) then
        ! An incoming phase boundary was stepped over by too much.
        ! Go back and step forward again, using a smaller step size.
        ! Note: this error code is not returned if the step size is less
        ! than or equal to the minimum value.
        xilim = xi1
        kly = 6
        go to 150
    end if

    ! Set up to go back to the previous point. Try a smaller step
    ! size or zero order Taylor's expansions, or write a pickup file
    ! for the last good point and stop.
150 continue

    call goback(acflg,acflg0,emop,emop0,emos,emos0,fje,fje0,fxi,fxi0,iemop,iemop0,iemos,iemos0,iindx0,iindx1,ipndx0,ipndx1,jpflag,jsflag,jreac,jreac0,kdim,kdim0,kmax,km1,km10,kmt,kmt0,kx1,kx10,kxt,kxt0,loph,losp,moph,moph0,mosp,mosp0,ncmpe,ncmpe0,npet,npetmx,npet0,npt,nptmax,nrct,nrctmx,nset,nsetmx,nset0,nst,nstmax,qreq,qriinf,sigmam,sigmm0,uzvec0,uzvec1,xi0,xi1)

    if (ier.eq.180 .and. delxi.le.dlxmin) then
        ! Normally this error code results in a reduction in the step
        ! size. This can't be done if the step size is already at the
        ! minimum value. Write a pickup file at the last good point
        ! and terminate work on the current problem.
        write (ux16,'(g12.5)') xi1
        call lejust(ux16)
        j2 = ilnobl(ux16)
        write (noutpt,1400) ux16(1:j2)
        write (nttyo,1400) ux16(1:j2)
1400 format(/' * Error- (EQ6/path) The equilibrium calculation',' failed at Xi= ',a,'.',/7x,"Can't cut the step size to try",' to recover. See previous notes',/7x,'and warnings for',' suggestions.')

        go to 180
    end if

    ! If a slide forward in EQ6/eqshel.f failed to get over a region of
    ! computational instability, write a pickup file at the last good
    ! point and terminate work on the current problem.
    if (ier .eq. 190) then
        go to 180
    end if

    if (ier .eq. 170) then
        ! The equilibrium calculation was completed successfully, but too
        ! much of a phase in the equilibrium system (ES) was destroyed to
        ! satisfy the requirements of the fluid-centered flow-through open
        ! system model. Go back and first transfer some of the phase's
        ! mass to the physically removed system (PRS).
        km1 = km10
        kmt = kmt0
        kx1 = kx10
        kxt = kxt0
        kdim = kdim0

        call copyaa(zvclg0,zvclg1,kdim)
        call copyaa(zvec0,zvec1,kdim)

        do kcol = 1,kdim
            uzvec1(kcol) = uzvec0(kcol)
            iindx1(kcol) = iindx0(kcol)
            ipndx1(kcol) = ipndx0(kcol)
        end do

        dxsave = 0.25*delxi
        kly = 6
        delxi = 0.
        deltim = 0.
        time1 = time0

        ! Recompute the temperature and pressure. Then recompute the
        ! thermodynamic and kinetic quantities which depend these
        ! variables.
        call tpadv(aadh,aadhh,aadhv,aaphi,abdh,abdhh,abdhv,abdoth,abdot,abdotv,adadhh,adadhv,adbdhh,adbdhv,adbdth,adbdtv,adh,adhfe,adhh,adhv,adhfs,adhfsd,advfe,advfs,advfsd,afcnst,al10,amu,aslm,aphi,aprehw,apresg,apresh,apx,avcnst,axhfe,axhfs,axhfsd,axlke,axlks,axlksd,axvfe,axvfs,axvfsd,bdh,bdhh,bdhv,bdot,bdoth,bdotv,dadhh,dadhv,dbdhh,dbdhv,dbdth,dbdtv,dhfe,dhfs,dhfsd,dvfe,dvfs,dvfsd,eact,ehfac,farad,hact,iact,iapxmx,iktmax,imchmx,imech,iopg,iopt,ipbtmx,ipch,ipchmx,ipcv,ipcvmx,ixrn1,ixrn2,jpfcmx,jpress,jptffl,jsol,jtemp,narxmx,narxt,narxth,nbasp,nbaspd,nbt,nbtd,nbtmax,ncmpr,ndrsr,ndrsrd,nmut,nmutmx,nopgmx,noptmx,noutpt,nptkmx,nptmax,nrct,nrctmx,nrk,nslt,nsltmx,nst,nstmax,ntpr,ntprmx,ntprt,nttkmx,nttyo,nweope,nwndpc,nxt,nxtmax,pmu,presg,presh,press,pressb,pressd,pslamn,ptk,rcnstv,rconst,rk,rkb,rtcnst,tempc,tempcb,tempcd,tempcu,tempk,time1,trkb,ttk,uphase,uspec,wfac,xhfe,xhfs,xhfsd,xi1,xlke,xlks,xlksd,xvfe,xvfs,xvfsd)

        if (nrct .ge. 1) then
            ! Increment the irreversible reactions (those associated with
            ! the "reactants"); update the ES mass balance totals
            ! accordingly.
            call reacts(cbsr,csts,delxi,drer0,iern1,ietmax,iktmax,iodb,jcode,jetmax,jgext,jreac,modr,modr0,morr,morr0,mrgers,mtb,mtb0,nbaspd,nbt,nbtmax,nbt1mx,ncmpr,nern1,nern2,nertmx,netmax,ngext,nodbmx,nord,noutpt,nptmax,nrct,nrctmx,nrd1mx,nrndex,nsrtmx,nstmax,nsts,nstsmx,nstsr,nttyo,nxridx,nxrtmx,rrelr0,rxbar,ureac,xirct,xirct0)
        end if

        ! Expand the system description.
        call ncmpex(acflg,act,actlg,cdrs,cegexs,cgexj,conc,conclg,cpgexs,egexjc,egexjf,egexs,eps100,fo2,fo2lg,fsort,fugac,fugalg,iern1,iern2,ietmax,ifrn1,ifrn2,igas,igstak,iindx1,ilrn1,ilrn2,imrn1,imrn2,istack,ixrn1,ixrn2,jcsort,jern1,jern2,jetmax,jflag,jgext,jgsort,jgstak,jjsort,jpflag,jsflag,jsitex,jssort,jstack,kbt,kdim,kelect,kmax,km1,ko2gaq,kwater,kxt,loph,losp,lsort,mgext,mrgexs,mtb,moph,mosp,narn1,narn2,nbasp,nbt,nbtmax,ncmpr,ndrs,ndrsmx,ndrsr,nelect,nern1,nern2,netmax,ngexsa,ngext,ngrn1,ngrn2,ngt,ngtmax,noutpt,no2gaq,nphasx,npt,nptmax,nst,nstmax,nttyo,omega,omeglg,press,qxbarw,q6mode,ugexj,ugexmo,uphase,uspec,xbar,xbarlg,xbarw,xbarwc,xbrwlc,xbrwlg,xlks,zchar,zgexj,zvclg1,zvec1)

        call pshfta(csts,emop,emop0,emos,emos0,fdpe0,fdpem1,fdse0,fdsem1,iemop,iemos,iern1,iern2,ietmax,iindx1,imrn1,imrn2,iodb,ipndx1,ixrn1,ixrn2,jcsort,jern1,jetmax,jgext,jpflag,jsflag,kbt,km1,kmax,kmt,kx1,kxt,loph,losp,moph,mosp,mprph,mprsp,mrgexs,mtb,mtb0,nbasp,nbaspd,nbt,nbtmax,ncmpe,ncmpr,netmax,ngext,nodbmx,nordmx,noutpt,npet,npetmx,npt,nptmax,nsetmx,nstmax,nsts,nstsmx,nstsr,nttyo,uaqsln,ufixf,uphase,uspec,xbar,xbarlg,zklgmn,zklogl,zvclg0,zvclg1,zvec0,zvec1)

        delxi = max(dxsave,dlxmin)
        nord = min(nord,kord)
        jcut = jcut + 1

        if (jcut .le. 10) then
            go to 130
        end if

        go to 180
    end if

    ! Cut the step size or order and try again or terminate.
    iexrt = 0
    jexrt = 0
    jscat = 0
    jscrt = 0
    qdump = .false.
    kly = 6

    if (delxi .gt. dlxmin) then
        write (noutpt,1440)
1440 format(' --- Cutting the step size and trying again ---')

        ncut = ncut + 1
        delxi = 0.25*delxi
        delxi = max(delxi,dlxmin)
        go to 130
    end if

    if (nord .gt. 0) then
        nord = 0
        write (noutpt,1430)
1430 format(' --- Cutting to order zero and trying again ---')

        go to 130
    end if

    ! Reaction path tracing has failed. Write a pickup file for the last
    ! good point and stop.
180 continue
    if (.not.qrpcfl .and. iopt(17).ge.0) then
        write (ux16,'(g12.5)') xi1
        call lejust(ux16)
        j2 = ilnobl(ux16)
        write (ux16a,'(g12.5)') xi0
        call lejust(ux16a)
        j3 = ilnobl(ux16a)
        write (noutpt,1450) ux16(1:j2),ux16a(1:j3)
        write (nttyo,1450)  ux16(1:j2),ux16a(1:j3)
1450 format(/' * Error - (EQ6/path) Reaction path tracing has',' failed',/7x,'at Xi= ',a,'. Will try to go back to',' Xi= ',a,',',/7x,'the last point of reaction progress',' at which the calculations',/7x,'were successful. If this',' succeeds, will then stop and write a',/7x,'pickup file.',' Try using this to restart the calculation.')

        qrpcfl = .true.
        iexrt = 0
        jexrt = 0
        jscat = 0
        jscrt = 0
        delxi = 0.
        qstopx = .true.
        go to 130
    else
        write (noutpt,1460)
        write (nttyo,1460)
1460 format(/" * Error - (EQ6/path) Can't recover by going back",' to the last point',/7x,'of reaction progress at which',' calculations were successful.',/7x,"Can't write a pickup",' file describing the system at that point.',/7x,'Try',' re-running the problem, setting the maximum number',' of steps',/7x,'(kstpmx) to the value which corresponds to',' that point. Then try',/7x,'restarting the calculation',' from that point.')

        stop
    end if

    ! The equilibrium calculation at the current point succeeded.
200 continue

    ! Get the activity of water.
    actwlg = actlg(narn1)
    actw = texp(actwlg)
    aw1 = actw

    ! Get the weights (masses) of solvent, total dissolved solutes,
    ! and aqueous solution, and get the aqeuous solution density.
    call gwdenp(adwipp,bdwipp,jcsort,mlmrra,mosp,mrmlra,mwtsp,narn1,narn2,nstmax,qdwipp,rhoc,rhowc,tdsgks,tdsglw,tdspkc,tdsplc,tempc,vosol,wfh2o,wftds,wkgwi,woh2o,wosol,wotds)

    if (qdwipp) then
        qrho = .true.
        rho = rhoc
    else
        qrho = .false.
        rho = 0.0
    end if

    ! Compute pH, Eh, and pe-, all with reference to appropriate
    ! pH scales. Also compute the pHCl.
    call gpheh(acflg,actlg,actwlg,adh,ah,ahmes,ahnbs,conc,eh,ehfac,ehmes,ehnbs,farad,fo2lg,fxi,iopg,mrmlra,nchlor,nhydr,nopgmx,noutpt,nstmax,nttyo,pch,pe,pemes,penbs,ph,phcl,phmes,phnbs,qphcl,qredox,qrho,xlke)

    wkgh2o = 1.e-3*woh2o
    wkgsol = 1.e-3*wosol
    ph1 = ph
    eh1 = eh
    fo2lg1 = fo2lg

    ! Calculate the new affinities and rates of the irreversible
    ! reactions. Determine if the kinetic mode requires a corrector
    ! step or a cut in the step size.
    if (nrct .le. 0) then
        go to 300
    end if

    ! Update the affinities of the irreversible reactions (afrc1).
    call raff(acflg,actlg,afcnst,affp,afrc1,bpx,cdrs,cgexj,ibpxmx,ibpxt,iern1,ietmax,iktmax,ixrn1,ixrn2,jcode,jern1,jern2,jetmax,jflag,jgext,jpflag,jsflag,jsol,ncmpr,ndrs,ndrsmx,ndrsr,nertmx,net,netmax,ngext,noutpt,nptmax,nrct,nrctmx,nrndex,nstmax,nttyo,nxridx,nxrtmx,nxtmax,rxbar,uphase,uspec,wfac,xbar,xbarlg,xgers,xlks)

    ! Check the irreversible reactions for saturation.
    call rsatch(csts,egers,egexs,iern1,ietmax,iindx1,iktmax,iopt,ipndx1,jcode,jern1,jern2,jetmax,jgext,jpflag,jreac,kmax,km1,kmt,kx1,kxt,loph,losp,moph,morr,mosp,mrgers,mtb,mtb0,nbaspd,nbtmax,ncmpr,nern1,nern2,nert,nertmx,netmax,ngext,noptmx,noutpt,nptmax,nrct,nrctmx,nrk,nrndex,nstmax,nsts,nstsmx,nstsr,nttyo,nxridx,nxrt,nxrtmx,qreq,rxbar,tolxsf,uphase,ureac,uspec,xbar,xbarlg,zvclg1,zvec1)

    ! Calculate rates of irreversible reactions at the new point by
    ! evaluating the corresponding rate laws.
    call rtcalc(act,afrc1,cdac,csigma,eps100,fkrc,idirec,imchmx,imech,iodb,iopt,jcode,jreac,morr,morr0,mwtrc,ndac,ndact,ndctmx,nodbmx,noptmx,nord,noutpt,nrk,nrct,nrctmx,nsk,nstmax,nttyo,prcinf,prminf,qriinf,rirec1,rk,rreac1,rrelr1,rtcnst,rrxfi1,sfcar,sfcar0,ssfcar,udac,ureac)

    if (qshoot) then
        go to 220
    end if

    if (iopt(2) .le. 0) then
        go to 230
    end if

    if (qriinf) then
        ! Trap an infinite time step.
        if (delxi .le. dlxmin) then
            ! Provisionally, take an infinite time step.
            deltim = prcinf
            time1 = prcinf

            ! Make sure that the calculated time does not exceed any
            ! any specified limits such as the maximum time just because
            ! delxi is at the minimum value.
            call tivchk(deltim,delxi,qtvchk,time1,time0,timemx,tiplol,tiplot,tiprnl,tiprnt,tolxst)

            if (qtvchk) then
                qriinf = .false.
            else
                write (noutpt,1490)
                write (nttyo,1490)
1490 format(/' --- Infinite time step ---',/)
            end if

            go to 230
        else
            ! Go back and cut the step size before taking an infinite
            ! time step.
            ier = 0
            go to 150
        end if
    end if

    if (ncorr .le. 0) then
        ! If using the predictor function, get predicted values for the
        ! inverse rate and the relative rates to compare with values
        ! obtained from rate law expressions. If using a corrector
        ! function, these variables are already available.
        call rtaylr(delxi,drer0,drir0,jreac,nord,nrct,nrctmx,nrd1mx,rirec0,rirecp,rrelr0,rrelrp)

        ! Zero the "old" values of btrmax and dlrmax, which are used
        ! to estimate the convergence rate for the ODE integrator.
        btrmxo = 0.
        dlrmxo = 0.
    end if

    if (iopt(14).le.0 .and. .not.qhcon) then
        ! If allowing both predictors, turn on the higher-order (stiff)
        ! corrector if the simple predictor is currently on but not
        ! rapidly converging.
        if (ncorr .ge. 4) then
            qhcon = .true.

            if (iodb(2) .gt. 1) then
                write (noutpt,1492)
                write (nttyo,1492)
1492 format(' Switching to the higher-order (stiff) corrector')
            end if
        end if
    end if

    ! Compute residual functions for the ODE integrator.
    call betars(alphar,betar,btrfnc,btrmax,btrmxo,ibtrmx,nrct,nrct1,nrctmx,rirec1,rirecp,rrelr1,rrelrp)

    if (iodb(8) .ge. 2) then
        write (noutpt,1500)
1500 format(/5x,'ODE residual functions:',//9x,'Reactant',22x,'betar',11x,'alphar',/)

        do nrc = 1,nrct
            write (noutpt,1510) ureac(nrc),betar(nrc),alphar(nrc)
1510 format(7x,a24,5x,1pe11.4,5x,e11.4)
        end do

        ux24 = 'Inverse rate'
        write (noutpt,1510) ux24,betar(nrct1),alphar(nrct1)
        write (noutpt,1512)
1512 format(1x)
    end if

    ! Test the accuracy of rate law integration. If the accuracy
    ! is sufficient, qodeok will be returned with a value of .true.
    call tstari(afrc1,alphar,betar,delxi,iodb,modr,nodbmx,noutpt,nrct,nrctmx,nrct1,nsscmx,qodeok,rirec1,rirecp,rrelr1,rrelrp,sscrew,time1,tistrt,ureac)

    if (iodb(8) .ge. 1) then
        write(noutpt,1520) btrmax,btrfnc
1520 format(5x,'btrmax= ',1pe10.3,5x,'btrfnc= ',e10.3)
    end if

    ! Satisfying any of the following tests causes the current
    ! results to be considered acceptable. Corrector iteration,
    ! if begun, is terminated.
    ! Force at a couple of corrector iterations.
    if (ncorr .lt. 2) then
        if (btrmax .gt. eps100) then
            qodeok = .false.
        end if
    end if

    if (qodeok) then
        go to 220
    end if

    if (qmod .or. qbye) then
        go to 220
    end if

    if (nord .le. 0) then
        go to 220
    end if

    ! Check for affinity sign flips.
    qaflip = .false.

    do j = 1,nrct
        if (afrc1(j) .ge. 0. .and. afrc0(j) .lt. 0.) then
            qaflip = .true.

            if (iodb(2) .gt. 1) then
                j2 = ilnobl(ureac(j))
                write (noutpt,1522) ureac(j)(1:j2)
                write (nttyo,1522) ureac(j)(1:j2)
1522 format(' Affinity sign flip (',a,': - to 0/+)')
            end if
        else if (afrc1(j) .lt. 0. .and. afrc0(j) .ge. 0.) then
            qaflip = .true.

            if (iodb(2) .gt. 1) then
                j2 = ilnobl(ureac(j))
                write (noutpt,1524) ureac(j)(1:j2)
                write (nttyo,1524) ureac(j)(1:j2)
1524 format(' Affinity sign flip (',a,': 0/+ to -)')
            end if
        end if
    end do

    ! Satisfying any of the following tests causes ODE correction
    ! to be terminated or foregone in favor of cutting the step
    ! size.
    if (qaflip) then
        go to 210
    end if

    if (iopt(14) .ge. 3) then
        go to 210
    end if

    if (ncorr.ge.6 .and. delxi.gt.dlxmin) then
        go to 210
    end if

    if (ncorr .ge. 8) then
        go to 210
    end if

    ! Set up to make an ODE corrector iteration step.
    ncorr = ncorr + 1

    if (ncorr .eq. 1) then
        qwhcfa = .true.
    else if (delxi .le. 0.) then
        qwhcfa = .true.
    else
        adx = abs(dlxode - delxi)
        qwhcfa = (adx/delxi) .gt. eps100
    end if

    if (qwhcfa) then
        ! This if block is normally entered only when ncorr = 1. If
        ! the step size is cut during ODE correction, the calculations
        ! in this block must be redone. When ncorr = 1, dlxode should
        ! be set to prcinf (practical infinity). This value is set at
        ! the start of a new step.
        ! Save the current step size.
        dlxode = delxi

        ! Calculate the dxsm10 vector from the step size (delxi) and the
        ! dxsm00 vector. Note that dxsm10(1) = delXi(0,-1) = -delxi, but
        ! dxsm10(2) = delXi(0,1), dxsm10(3) = delXi(0,2), etc. More
        ! simply put, dxsm10 is the vector of cumulative step sizes to
        ! the new point, but centered at point 0. In constrast, dxsm00
        ! is the vector of cumulative step sizes to point 0, centered
        ! at point 0.
        do j = 1,kordlm
            dxsm10(j + 1) = dxsm00(j)
        end do

        dxsm10(1) = -delxi

        ! Now calclate dxsm11, the vector of cumulative step sizes to
        ! the new point, centered at the new point.
        do n = 1,kordlm
            j = kordlm - n + 1
            dxsm11(j + 1) = dxsm00(j) + delxi
        end do

        dxsm11(1) = delxi

        ! Compute the matrix (akmat1 = FC) used to calculate derivatives
        ! for truncated Taylor's series, given corrector-based finite
        ! differences.
        ! Calling sequence substitutions:
        !   akmat1 for akmat0
        !   dxsm10 for dxsm00
        call gakmat(akmat1,dxsm10,nord,nrd1mx)

        if (qhcon) then
            ! Calculate the w factor (whcfac) needed for higher-order
            ! (stiff) ODE corrections. This depends only on the recent
            ! step size history, including the current step size.
            call gwhcfa(akmat1,delxi,dxsm11,hhcvec,nord,nordmx,nrd1mx,whcfac,xhcvec)
        end if
    end if

    if (.not. qhcon) then
        ! Use the simple corrector. The rates (rrelr1 and rirec1)
        ! calculated from the governing rate equations  will be used
        ! below as the rates at the new point when constructing new
        ! finite differences for the corrector. Here just record the
        ! correction in the rates at the new point.
        do j = 1,nrct1
            delvcr(j) = alphar(j)
        end do
    end if

    if (qhcon) then
        ! Use the higher-order (stiff) ODE integrator. Adjust rrelr1
        ! and rirec1 using this scheme. These new values will be used
        ! to compute a new finite-difference-based Taylor's-series-
        ! formatted corrector function. Basically, the current adjustment
        ! represents an application of the Newton-Raphson method.
        ! Set up the right-hand-size vector.
        do j = 1,nrct1
            rhsvcr(j) = -alphar(j)
        end do

        !        Recompute the cdacb, ndacb, and ndactb arrays prior to
        !        calculating the Jacobian matrix J[r].
        ! XX - Note: this can be changed by basis switching. Am probably
        ! XX   recalculating these arrays many times unnecessarily.
        call gndacb(cdac,cdacb,cdrs,eps100,imech,imchmx,jflag,nbasp,nbt,nbtmax,ndac,ndacb,ndact,ndactb,ndctmx,ndrs,ndrsmx,ndrsr,nrct,nrctmx,nstmax)

        ! Compute the Jacobian matrix J[r] (armatr).
        call garmat(act,afrc1,aimatr,al10,armatr,cdac,cdacb,cdrs,csigma,csts,delvec,dlogxw,dvjdte,eact,eps100,fkrc,gmmatr,hact,iact,idirec,iindx1,iktmax,imchmx,imech,ipivot,jcode,jreac,jtemp,kbt,kdim,kmax,mmmatr,morr,mwtrc,nbasp,nbt,nbtmax,ncmpr,ndac,ndacb,ndact,ndctmx,ndrs,ndrsmx,ndrsr,noutpt,nptmax,nrct,nrctmx,nrct1,nrk,nrndex,nsk,nstmax,nsts,nstsmx,nstsr,nttkmx,nttyo,nxridx,nxrtmx,rirec1,rk,rkb,rreac1,rrelr1,rrxfi1,rtcnst,rxbar,sfcar,sgmatr,ssfcar,tempc,tempcb,tempk,ttk,ureac,whcfac,xi1,xlks,xxmatr,xymatr)

        ! Solve for the correction vector (delvcr). If EQLIBU/msolvr.f
        ! can't solve the matrix, it is because the matrix is either
        ! zero (ier = 1) or non-zero, but computationally singular
        ! (ier = 2).
        ! Calling sequence substitutions:
        !   armatr for aamatr
        !   delvcr for delvec
        !   grmatr for gmmatr
        !   ipivtr for ipivot
        !   nrct1  for kdim
        !   nrct1  for kmax
        !   rhsvcr for rhsvec
        qpr = .false.
        call msolvr(armatr,delvcr,grmatr,ier,ipivtr,nrct1,nrct1,noutpt,nttyo,qpr,rhsvcr)

        if (ier .le. 0) then
            ! Correct the predicted values. Note that the results are
            ! stored in rrelr1 and rirec1, not rrelrp and rirecp.
            do j = 1,nrct
                rrelr1(j) = rrelrp(j) + delvcr(j)
            end do

            rirec1 = rirecp + delvcr(nrct1)
        else
            ! Higher-order (stiff) ODE correction failed on this corrector
            ! iteration. Drop back to the simple corrector.
            do j = 1,nrct1
                delvcr(j) = alphar(j)
            end do
        end if
    end if

    if (iodb(8) .ge. 2) then
        ux24 = 'Simple'

        if (qhcon) then
            ux24 = 'Higher-order (stiff)'
        end if

        j2 = ilnobl(ux24)
        write (noutpt,1530) ux24(1:j2)
1530 format(/5x,a,' ODE corrector:',//9x,'Reactant',22x,'delvcr',7x,'Corrected r',/)

        do nrc = 1,nrct
            write (noutpt,1540) ureac(nrc),delvcr(nrc),rrelr1(nrc)
1540 format(7x,a24,5x,1pe11.4,5x,e11.4)
        end do

        ux24 = 'Inverse rate'
        write (noutpt,1540) ux24,delvcr(nrct1),rirec1
        write (noutpt,1512)
    end if

    ! Find the max norm of the correction vector (delvcr). In the limit
    ! of the solution, the (unrelaxed) correction term bounds the error
    ! in the variables being calculated.
    idlrmx = iarmxn(delvcr,nrct1)
    dlrmax = 0.

    if (idlrmx .gt. 0) then
        dlrmax = abs(delvcr(idlrmx))
    end if

    ! Calculate the delvcr improvement function (dlrfnc).
    dlrfnc = 0.

    if (dlrmxo .gt. 0.) then
        dlrfnc = (dlrmxo -dlrmax)/dlrmxo
    end if

    dlrmxo = dlrmax

    if (iodb(8) .ge. 2) then
        write(noutpt,1550) dlrmax,dlrfnc
1550 format(9x,'dlrmax= ',1pe10.3,5x,'dlrfnc= ',e10.3)
    end if

    ! Compute the new set of finite differences for a corrector
    ! step.
    call corrfd(delxi,dxsm11,fdlim,fdre0,fdre1,fdri0,fdri1,fdrr0,fdrr1,iodb,iopt,jreac,nodbmx,noptmx,nord,nordmx,noutpt,npts,nrct,nrctmx,nrd1mx,rirec0,rirec1,rreac0,rreac1,rrelr0,rrelr1)

    ! Now recompute the corresponding set of derivatives for the
    ! original base point.
    ! Calling sequence substitutions:
    !   akmat1 for akmat0
    !   fdri1 for fdri0
    !   fdrr1 for fdrr0
    call rderiv(akmat1,drer0,drir0,fdri1,fdrr1,jreac,nord,nrct,nrctmx,nrd1mx)

    ! Save the current corrected rate values. They will be the
    ! predicted rates for the next corrector iteration. This
    ! saves the effort of recalculating them from finite
    ! differences.
    do j = 1,nrct
        rrelrp(j) = rrelr1(j)
    end do

    rirecp = rirec1

    ! Go back to the base point and try again with the current
    ! corrector function.
    go to 140

210 continue

    ! The corrector algorithm has failed to satsify the specified
    ! criteria at the current order and step size.
    if (delxi .gt. dlxmin) then
        ! Go back and cut the step size to satisfy the ODE corrector
        ! tolerance. When the step size is cut to the minimum value
        ! (dlxmin), the order will be reduced to zero and the
        ! ODE integration will be considered sufficiently accurate.
        ncorr = 0
        ier = 0
        go to 150
    end if

220 continue

    ! The AE/ODE algorithms have succeeded for the current step.
    ! However, it may be necessary to go back and cut the step size
    ! if an event has not been sufficiently accurately located.
    ! Make sure that the currently targeted time-based print point
    ! has not been exceeded.
    tiprxx = min(tiprnl,tiprnt)

    if (((time1 - tiprxx)/tiprxx) .gt. tolxst) then
        if (delxi .gt. dlxmin) then
            ! Go back and cut the step size.
            ncorr = 0
            ier = 0
            go to 150
        else
            ! Set the current time equal to the targeted time-based print
            ! point
            time1 = tiprxx
            deltim = tiprxx - time0
        end if
    end if

    ! Make sure that the currently targeted time-based plot point
    ! has not been exceeded.
    tiplxx = min(tiplol,tiplot)

    if (((time1 - tiplxx)/tiplxx) .gt. tolxst) then
        if (delxi .gt. dlxmin) then
            ! Go back and cut the step size.
            ncorr = 0
            ier = 0
            go to 150
        else
            ! Set the current time equal to the targeted time-based plot
            ! point
            time1 = tiplxx
            deltim = tiplxx - time0
        end if
    end if

    ! Make sure that the maximum time has not been exceeded.
    if (((time1 - timemx)/timemx) .gt. tolxst) then
        if (delxi .gt. dlxmin) then
            ! Go back and cut the step size.
            ncorr = 0
            ier = 0
            go to 150
        else
            ! Set the current time equal to the maximum time.
            time1 = timemx
            deltim = timemx - time0
        end if
    end if

    if (iopt(14).le.0 .and. qhcon) then
        ! If allowing both predictors, turn on the simple corrector
        ! if the higher-order (stiff)  predictor is currently on
        ! and rapidly converging.
        if (ncorr .le. 2) then
            qhcon = .false.

            if (iodb(2) .gt. 1) then
                write (noutpt,1560)
                write (nttyo,1560)
1560 format(' Switching to the low-order (simple) corrector')
            end if
        end if
    end if

230 continue

    ! Check for reactants that are newly exhausted and for formerly
    ! exhausted reactants that are starting to precipitate. Note
    ! that iexrt is the number of newly exhausted reactants, not
    ! the total number of currently exhausted reactants. Here
    ! jexrt is the number of formerly exhausted reactants that
    ! have been reactivated.
    iexrt = 0
    jexrt = 0

    do nrc = 1,nrct
        if (jreac(nrc).eq.0 .or. jreac(nrc).eq.-1) then
            ! The reactant is active. See if there are any moles
            ! remaining. If not, the reactant is exhausted if the
            ! affinity is positive, the relative rate is positive,
            ! or a finite amount was present at the previous point
            ! of reaction progress (all three conditions should
            ! agree).
            if (morr(nrc) .le. 0.) then
                if (afrc1(nrc).gt.eps100 .or. rrelr1(nrc).gt.eps100        .or. morr0(nrc).gt.0.) then
                    jreac(nrc) = 1
                    iexrt = iexrt + 1
                    iexr(iexrt) = nrc
                end if
            end if
        else if (jreac(nrc) .eq. 1) then
            ! The reactant is exhausted (has no moles remaining). See if
            ! the affinity favors formation. If so, activate the reactant
            ! so that this may proceed. Note that reactivation is not
            ! appropriate for the case in which nrk(2,nrc) = 0. In that
            ! case, precipitation is governed by partial equilibrium,
            ! not a rate law.
            if (nrk(2,nrc) .ne. 0) then
                if (afrc0(nrc).gt.eps100 .and. afrc1(nrc).le.eps100) then
                    jreac(nrc) = 0
                    jexrt = jexrt + 1
                    jexr(jexrt) = nrc
                end if
            end if
        end if
    end do

    ! Check to see that the point of reactivation of a formerly
    ! exhausted reactant has been located sufficiently accurately.
    ! A finite-difference based search does not guarantee this.
    jexrtx = 0

    do j = 1,jexrt
        nrc = jexr(j)

        if (abs(afrc1(nrc)) .gt. eps100) then
            jexrtx = jexrtx + 1

            if (iodb(1) .gt. 0) then
                j2 = ilnobl(ureac(nrc))
                write (noutpt,1570) ureac(nrc)(1:j2),afrc1(nrc)
1570 format(/3x,'The affinity of formerly exhausted reactant',' ',a,/5x,'is ',1pe11.4,' kcal. This is outside the',' initially permitted',/5x,'range for reactivation to',' allow precipitation according to',/5x,'a rate law.',/)
            end if
        end if
    end do

    if (jexrtx .gt. 0) then
        ! Determine whether or not to reduce delxi and go back and try
        ! again to better locate an event in question.
        call dadjdx(delxi,dlxmin,iodb,nodbmx,noutpt,qadjdx)

        if (qadjdx) then
            ier = 0
            go to 150
        end if
    end if

    ! Check for reactants that have affinities or reaction rates which
    ! are about to change sign. Here jscat is the number of reactants
    ! with affinities about to undergo a sign change, and jscrt is
    ! the number with rates about to undergo a sign change. Note that
    ! a reactant with an affinity undergoing a sign change may not have
    ! reaction rate undergoing a sign change, as the rate may be
    ! truncated to zero for a given sign of the affinity. Also, even
    ! if theoretically both the affinity and the rate should undergo
    ! a sign change simultaneously, the testing tolerances may be such
    ! that a sign change for only one or the other condition is caught
    ! and reported.
    jscat = 0

    do nrc = 1,nrct
        if (jreac(nrc).eq.0 .or. jreac(nrc).eq.-1) then
            if (nrk(2,nrc) .ne. 0) then
                if (afrc0(nrc).gt.tolsar .and. afrc1(nrc).le.tolsar) then
                    jscat = jscat + 1
                    jsca(jscat) = nrc
                else if (afrc0(nrc).lt.-tolsar .and. afrc1(nrc).ge.-tolsar) then
                    jscat = jscat + 1
                    jsca(jscat) = nrc
                end if
            end if
        end if
    end do

    jscrt = 0

    do nrc = 1,nrct
        if (jreac(nrc).eq.0 .or. jreac(nrc).eq.-1) then
            if (nrk(2,nrc) .ne. 0) then
                if (rrelr0(nrc).gt.tolsrr .and. rrelr1(nrc).le.tolsrr) then
                    jscrt = jscrt + 1
                    jscr(jscrt) = nrc
                else if (rrelr0(nrc).lt.-tolsrr .and. rrelr1(nrc).ge.-tolsrr) then
                    jscrt = jscrt + 1
                    jscr(jscrt) = nrc
                end if
            end if
        end if
    end do

    ! Check for oversteps. The finite-difference based searches do
    ! not guarantee sufficient accuracy.
    jscatx = 0

    do j = 1,jscat
        nrc = jsca(j)

        if (abs(afrc1(nrc)) .gt. tolsar) then
            jscatx = jscatx + 1

            if (iodb(1) .gt. 0) then
                j2 = ilnobl(ureac(nrc))
                write (noutpt,1580) ureac(nrc)(1:j2),afrc1(nrc)
1580 format(/3x,'The affinity of reactant ',a,' is ',1pe11.4,' kcal. This is',/5x,'too much of an overstep past the',' point at which the affinity of this',/5x,'reactant',' changes sign.',/)
            end if
        end if
    end do

    jscrtx = 0

    do j = 1,jscrt
        nrc = jscr(j)

        if (abs(rrelr1(nrc)) .gt. tolsrr) then
            jscrtx = jscrtx + 1

            if (iodb(1) .gt. 0) then
                j2 = ilnobl(ureac(nrc))
                write (noutpt,1590) ureac(nrc)(1:j2),rrelr1(nrc)
1590 format(/3x,'The relative rate of reactant ',a,' is ',1pe11.4,'. This is',/5x,'too much of an overstep past the',' point at which the relative rate of',/5x,'this',' reactant changes sign.',/)
            end if
        end if
    end do

    if (jscatx.gt.0 .or. jscrtx.gt.0) then
        ! Determine whether or not to reduce delxi and go back and try
        ! again to better locate an event in question.
        call dadjdx(delxi,dlxmin,iodb,nodbmx,noutpt,qadjdx)

        if (qadjdx) then
            ier = 0
            go to 150
        end if
    end if

300 continue

    ! Check for oversteps of the pH, Eh, log fO2, or aw with respect to
    ! currently defined print point values.
    call cophpr(actw,aw0prn,aw1prn,delxi,dlxmin,eh,eh0prn,eh1prn,fo2lg,iodb,nodbmx,noutpt,o20prn,o21prn,ph,ph0prn,ph1prn,qadjdx,qredox,tolxsu)

    if (qadjdx) then
        ier = 0
        go to 150
    end if

    ! Check for oversteps of the pH, Eh, log fO2, or aw with respect to
    ! currently defined plot point values.
    call cophpl(actw,aw0plo,aw1plo,delxi,dlxmin,eh,eh0plo,eh1plo,fo2lg,iodb,nodbmx,noutpt,o20plo,o21plo,ph,ph0plo,ph1plo,qadjdx,qredox,tolxsu)

    if (qadjdx) then
        ier = 0
        go to 150
    end if

    ! Check for oversteps of the pH, Eh, log fO2, or aw with respect to
    ! requested minimum and maximum values.
    call cophlm(actw,awmax,awmin,delxi,dlxmin,eh,ehmax,ehmin,fo2lg,iodb,nodbmx,noutpt,o2max,o2min,ph,phmax,phmin,qadjdx,qredox,tolxsu)

    if (qadjdx) then
        ier = 0
        go to 150
    end if

    ! Reset the iemop, emop, iemos, emos, etc., arrays.
    if (qmod) then
        ! The ES phase assemblage has changed. Re-set the index arrays
        ! associated with finite-difference description of the numbers of
        ! mole of phases and species in the Equilibrium System (ES).
        call iiemop(iemop,iemos,iindx1,ipndx1,jsflag,kdim,kmax,ncmpe,ncmpr,noutpt,npet,npetmx,npt,nptmax,nset,nsetmx,nstmax,nttyo,uaqsln,uspec,uphase)
    end if

    ! Update the numbers of moles variables for phases present in
    ! the ES.
    do npe = 1,npet
        np = iemop(npe)
        emop(npe) = moph(np)
        nr1 = ncmpe(1,npe)
        nr2 = ncmpe(2,npe)

        do nse = nr1,nr2
            ns = iemos(nse)
            emos(nse) = mosp(ns)
        end do
    end do

    ! Recalculate the total affinity (aft1).
    call caft1(afrc1,aft1,nrct,nrctmx,rrelr1)

    ! Write entertainment for the user, showing the progress
    ! of the current run.
    call wrentu(actw,eh,fo2lg,iopg,iopt,kstep,nopgmx,noptmx,nttyo,ph,qredox,time1,xi1)

    qprntx = (xi1 + eps100) .ge. xiprnt
    qprntt = (time1 + eps100) .ge. tiprnt
    qprnlx = (xi1 + eps100) .ge. xiprnl
    qprnlt = (time1 + eps100) .ge. tiprnl

    qprph0 = abs(ph - ph0prn) .le. tolxsu
    qprph1 = abs(ph - ph1prn) .le. tolxsu
    qpreh0 = abs(eh - eh0prn) .le. tolxsu
    qpreh1 = abs(eh - eh1prn) .le. tolxsu
    qpro20 = abs(fo2lg - o20prn) .le. tolxsu
    qpro21 = abs(fo2lg - o21prn) .le. tolxsu
    qpraw0 = abs(actw - aw0prn) .le. tolxsu
    qpraw1 = abs(actw - aw1prn) .le. tolxsu

    qzprnt = qprntx .or. qprntt .or. qprnlx .or. qprnlt .or. qprph0 .or. qprph1 .or. qpreh0 .or. qpreh1 .or. qpro20 .or. qpro21 .or. qpraw0 .or. qpraw1 .or. kstppr.ge.ksppmx

    qplotx = (delxi + eps100) .ge. (xiplot - xi0)
    qplott = (deltim + eps100) .ge. (tiplot - time0)
    qplolx = (delxi + eps100) .ge. (xiplol - xi0)
    qplolt = (deltim + eps100) .ge. (tiplol - time0)

    qplph0 = abs(ph - ph0plo) .le. tolxsu
    qplph1 = abs(ph - ph1plo) .le. tolxsu
    qpleh0 = abs(eh - eh0plo) .le. tolxsu
    qpleh1 = abs(eh - eh1plo) .le. tolxsu
    qplo20 = abs(fo2lg - o20plo) .le. tolxsu
    qplo21 = abs(fo2lg - o21plo) .le. tolxsu
    qplaw0 = abs(actw - aw0plo) .le. tolxsu
    qplaw1 = abs(actw - aw1plo) .le. tolxsu

    qzplot = qplotx .or. qplott .or. qplolx .or. qplolt .or. qplph0 .or. qplph1 .or. qpleh0 .or. qpleh1 .or. qplo20 .or. qplo21 .or. qplaw0 .or. qplaw1 .or. kstppl.ge.ksplmx

    qzdump = (xi1 + eps100) .ge. xidump

    if (qzdump) then
        qdump = .true.
    end if

    if (qrapch .and. .not.qmod .and. .not.qbye) then
        kly = max(kly,6)
    end if

    ! Check to see which variable is changing most rapidly.
    if (qmod .or. qbye) then
        write (noutpt,1900) kstep,iter,ncorr
1900 format(' Steps completed= ',i5,', iter= ',i3,', ncorr= ',i1)
    else
        inmax = 0
        kzmax = 0
        adlzlg = 0.

        do kcol = 1,kdim
            adel = abs(zvclg1(kcol) - zvclg0(kcol))

            if (adel .gt. adlzlg) then
                inmax = iindx1(kcol)
                kzmax = kcol
                adlzlg = adel
            end if
        end do

        adlzlg = 0.

        if (kzmax .le. 0) then
            write (noutpt,1900) kstep,iter,ncorr
        else
            adlzlg = zvclg1(kzmax)
            ustr24 = 'Error'

            if (kzmax .le. kbt) then
                ns = nbasp(inmax)
                ustr24 = uspec(ns)(1:24)
            else if (kzmax .le. kxt) then
                ustr24 = uspec(inmax)(1:24)
            end if

            j2 = ilnobl(ustr24)
            write (noutpt,1910) kstep,iter,ncorr,ustr24(1:j2),adlzlg
1910 format(' Steps completed= ',i5,', iter= ',i3,', ncorr= ',i1,/' Most rapidly changing is zvclg1(',a,')= ',f11.4)
        end if
    end if

    ! Check gases for which the fugacity is supposed to be fixed.
    iwnffg = 0

    do n = 1,nffg
        ng = jffg(n)
        xlf = xlkffg(n)
        xf = texp(xlf)
        dx = fugalg(ng) - xlkffg(n)

        if (abs(dx) .gt. toldl) then
            if (iwnffg .le. 0) then
                nwnffg = nwnffg + 1
            end if

            iwnffg = iwnffg + 1

            if (nwnffg .le. nlwffg) then
                j2 = ilnobl(uffg(n))
                write (noutpt,2000) uffg(n)(1:j2),fugac(ng),xf,xi1
                write (nttyo,2000) uffg(n)(1:j2),fugac(ng),xf,xi1
2000 format(/' * Warning - (EQ6/path) The fugacity of ',a,/7x,'is now ',1pg12.5,' bars. It is supposed to be fixed',/7x,'at ',e12.5,' bars. At the present value of reaction',/7x,'progress (Xi= ',e12.5,'), there is no longer a',/7x,'sufficient mass of the gas component in the system',/7x,'to "saturate" it at the desired fugacity value. If',/7x,'you restart the run, you can add such mass using the',/7x,'moffg parameter on the input file. Add only about',/7x,'0.5 to 1.0 mole at a time. You can also increase',/7x,'the amount of gas component present as the run',/7x,'progresses by specifying the gas as a reactant.')
            end if
        end if
    end do

    if (nwnffg .eq. nlwffg) then
        write (noutpt,2010)
        write (nttyo,2010)
2010 format(/' * Warning - (EQ6/path) No more warnings will be',/7x,'issued regarding gas fugacities not being at desired',/7x,'fixed values.')
    end if

    ! Test for very little remaining solvent water.
    if (wkgh2o .lt. 1.e-10) then
        if (zvec1(1) .lt. zvec0(1)) then
            write (noutpt,2020) wkgh2o,wkgsol
            write (nttyo,2020) wkgh2o,wkgsol
2020 format(/' * Note - (EQ6/path) The amount of remaining',' solvent water',/7x,'is only ',1pe12.5,' kg and is',' decreasing with reaction progress.',/7x,'The amount of',' remaining aqueous solution is ',e12.5,' kg.',/7x,'This code is not designed to deal with fully dry',' systems. This',/7x,'run will now be terminated.')

            qvlsow = .true.
        end if
    end if

    if (wkgh2o .lt. 1.e-8) then
        if (iwdh2o .ge. 0) then
            nwdh2o = nwdh2o + 1

            if (nwdh2o .le. 8) then
                if (zvec1(1) .lt. zvec0(1)) then
                    write (noutpt,2030) wkgh2o,wkgsol
                    write (nttyo,2030) wkgh2o,wkgsol
2030 format(/' * Warning - (EQ6/path) The amount of remaining',' solvent water',/7x,'is only ',1pe12.5,' kg and is',' decreasing with reaction progress.',/7x,'The amount of',' remaining aqueous solution is ',e12.5,' kg.',/7x,'This code is not designed to deal with fully dry',' systems. This',/7x,"run can't continue much farther.")
                else
                    write (noutpt,2040) wkgh2o,wkgsol
                    write (nttyo,2040) wkgh2o,wkgsol
2040 format(/' * Warning - (EQ6/path) The amount of remaining',' solvent water',/7x,'is only ',1pe12.5,' kg. The amount',' of remaining aqueous solution',/7x,'is only ',e12.5,' kg. This code is not designed to deal with',/7x,' fully dry systems. This run is near the limit of the'," code's capability.")
                end if

                iwdh2o = -10
            end if
        else
            iwdh2o = iwdh2o + 1
        end if
    end if

    ! Test for very high ionic strength.
    if (fxi .gt. 100.) then
        write (ux16,'(f10.3)') fxi
        call lejust(ux16)
        j2 = ilnobl(ux16)
        write (noutpt,2050) ux16(1:j2)
        write (nttyo,2050) ux16(1:j2)
2050 format(/' * Note - (EQ6/path) The ionic strength is ',a,' molal.',/7x,'This run will now be terminated.')

        qvhfxi = .true.
    end if

    ! Do automatic basis switching. Optimize the basis set so that
    ! the basis species for each mass balance tends to be the species
    ! which dominates that mass balance.
    qabswx = .false.

    if (iopt(12).gt.0 .and. kstpab.ge.20) then
        ! Automatic basis switching after solving at the current point
        ! of reaction progress is turned on, and this is the 20th step
        ! since the last time the need for this was checked.
        call absswb(adhfs,adhfsx,advfs,advfsx,avcnst,axhfs,axhfsx,axlks,axlksx,axvfs,axvfsx,cdrs,cdrtw,cdrsx,cdrw,csts,dhfs,dvfs,eps100,ibswx,iindx1,iodb,ipch,ipchmx,ipcv,ipcvmx,jcsort,jflag,jsflag,kbt,kmax,mosp,mtb,narn1,narn2,narxmx,narxt,nbasp,nbaspd,nbaspx,nbt,nbtmax,nbw,ndrs,ndrsmx,ndrsr,ndrsrx,ndrsx,nelect,nhydr,nodbmx,no2gaq,noutpt,nst,nstmax,nsts,nstsmx,nstsr,nswtch,ntpr,ntprmx,nttyo,presg,press,qbassw,qbswx,tempc,uspec,uzvec1,weight,xhfs,xvfs,xlks)

        if (nswtch .gt. 0) then
            ! Set a flag noting that automatic basis switching has just
            ! been done.
            qabswx = .true.

            write (noutpt,2100) nswtch
2100 format(/' ',i2,' basis switches were executed automatically',' after solving',/'at the current value of reaction',' progress.')
        end if

        ! Reset the step counter which contains the number of steps since
        ! the last time a check was made of the need for automatic basis
        ! switching of this kind.
        kstpab = 0
    end if

    ! Save the order parameters.
    nordsv = nord
    kordsv = kord

    ! Check to see if it is necessary to drop the order of the
    ! Taylor's series to zero because a point has been reached
    ! at which the derivatives of the variables being described
    ! by finite differences are not continuous.
    if (qmod .or. qbye .or. qabswx .or. iexrt.gt.0 .or.  jexrt.gt.0 .or. jscat.gt.0 .or. jscrt.gt.0 .or.  qtrch .or. qreq) then
        ! Drop the order of the Taylor's series to zero. Zero all
        ! finite-difference and derivative data.
        if (iopt(1) .eq. 2) then
            ! Set flags for dump and to control special print for
            ! the case of the fluid-centered flow-through open system.
            qdump = .true.
            qdmpr1 = .true.
            qdmpr2 = .true.
        end if

        xilim = prcinf
        delxi = dlxmx0
        npts = 0

        call initaz(dxsm00,nrd1mx)
        call initaz(dxsm10,nrd1mx)
        call initaz(dxsm11,nrd1mx)

        nmax = nrd1mx*kmax
        call initaz(dzvc0,nmax)
        call initaz(dzvc0s,nmax)
        call initaz(fdzv0,nmax)
        call initaz(fdzvm1,nmax)

        call initaz(fdri0,nrd1mx)
        call initaz(fdrim1,nrd1mx)
        call initaz(drir0,nrd1mx)
        call initaz(drir0s,nrd1mx)

        nmax = nordmx*nrctmx
        call initaz(fdar0,nmax)
        call initaz(fdarm1,nmax)
        call initaz(dafrc0,nmax)

        nmax = nrd1mx*nrctmx
        call initaz(fdrr0,nmax)
        call initaz(fdrrm1,nmax)
        call initaz(drer0,nmax)
        call initaz(drer0s,nmax)
    end if

    npts = npts + 1
    kord = npts - 1 - ndelay
    kord = min(kord,kordlm)
    kord = max(0,kord)

    ! Check to see that the reaction path is not constant.
    qconst = .false.

    if (qcntmp .and. qcnpre) then
        if (nrct .gt. 0) then
            do nrc = 1,nrct
                if (rrelr1(nrc) .ne. 0.) then
                    go to 330
                end if
            end do
        end if

        qconst = .true.
    end if

330 continue

    if (iopt(2) .le. 0) then
        ! If the temperature and pressure are constant, check to see if
        ! the total affinity is repeatedly nearly zero.
        if (qcntmp .and. qcnpre) then
            if (aft1 .le. tolsar) then
                naft1 = naft1 + 1
                qaft1 = naft1 .ge. 20
            else
                naft1 = 0
            end if
        end if
    end if

    if (iopt(2) .gt. 0) then
        if (qshoot) then
            jscat = 0
            jscrt = 0
        else
            ! Check to see if the code is just making small oscillations
            ! about the final equilibrium point.
            if (nord .le. 0) then
                if (aft1.gt.aft0 .and. aftm1.gt.aft0) then
                    kaft1 = kaft1 + 1

                    if (aft1 .le. 0.001) then
                        jscat = 0
                        jscrt = 0
                    end if
                else if (aft1.lt.aft0 .and. aftm1.lt.aft0) then
                    kaft1 = kaft1 + 1

                    if (aft1 .le. 0.001) then
                        jscat = 0
                        jscrt = 0
                    end if
                else
                    kaft1 = 0
                end if
            else
                kaft1 = 0
            end if

            if (kaft1 .ge. 6) then
                if (aft1 .le. 0.001) then
                    ! The total affinity is oscillating about some finite
                    ! value close to zero.
                    time1 = min(timemx,tiprnt,tiprnl,tiplot,tiplol)
                    deltim = time1 - time0
                    qshoot = .true.
                    qprntt = (time1 + eps100) .ge. tiprnt
                    qprnlt = (time1 + eps100) .ge. tiprnl
                    qzprnt = qzprnt .or. qprntt .or. qprnlt
                    qplott = (time1 + eps100) .ge. tiplot
                    qplolt = (time1 + eps100) .ge. tiplol
                    qzplot = qzplot .or. qplott .or. qplolt
                end if
            end if
        end if
    end if

    ! If doing a fluid-centered flow-through system mode run, check
    ! to see if a "short" print description is required after a transfer
    ! to the PRS.
    if (iopt(1).eq.2 .and. nord.gt.0) then
        qftpr2 = .not.qdmpr1 .and. qdmpr2
    end if

    ! Set upper limits on the scale factor for the next step. Here
    ! kly is a delay factor that comes into play when the step size
    ! has been cut recently in connection with a search.
    scalim = 10.

    if (kly .gt. 0) then
        scalim = 2.0
    end if

    ! Determine if a print and/or plot of the system should
    ! be made at the current point of reaction progress.
    if (xi1 .ge. ximax) then
        go to 340
    end if

    if (iopt(2) .gt. 0) then
        if (time1 .ge. ((1. - tolxst)*timemx)) then
            go to 340
        end if
    end if

    if (ph .le. (phmin + tolxsu)) then
        go to 340
    end if

    if (ph .ge. (phmax - tolxsu)) then
        go to 340
    end if

    if (qredox) then
        if (eh .le. (ehmin + tolxsu)) then
            go to 340
        end if

        if (eh .ge. (ehmax - tolxsu)) then
            go to 340
        end if

        if (fo2lg .le. (o2min + tolxsu)) then
            go to 340
        end if

        if (fo2lg .ge. (o2max - tolxsu)) then
            go to 340
        end if
    end if

    if (actw .le. (awmin + tolxsu)) then
        go to 340
    end if

    if (actw .ge. (awmax - tolxsu)) then
        go to 340
    end if

    if (kstep .ge. kstpmx) then
        go to 340
    end if

    if (qstopx) then
        go to 340
    end if

    if (iexrt.gt.0 .or. jexrt.gt.0) then
        go to 340
    end if

    if (jscat.gt.0 .or. jscrt.gt.0) then
        go to 340
    end if

    if (qconst .or. qmod .or. qbye .or. qreq .or. qtrch) then
        go to 340
    end if

    if (qaft1) then
        go to 340
    end if

    if (qvlsow .or. qvhfxi) then
        go to 340
    end if

    if (qzprnt .or. qftpr2 .or. qdump) then
        go to 350
    end if

    if (qzplot) then
        go to 360
    end if

    go to 100

340 continue
    kly = 0

    if (qftpr2) then
        qdmpr2 = .false.
    end if

    qftpr2 = .false.

350 continue

    ! Compute apparent "whole-phase" equivalent fractions and mole
    ! fractions of the exchange ions present in generic ion exchanger
    ! phases. Cations and anions are treated separately in these
    ! calculations.
    call gegexw(cegexs,egexpc,egexpa,egexw,iern1,iern2,ietmax,jern1,jetmax,jgext,kern1,kern2,ketmax,kgexsa,moph,mosp,netmax,ngexsa,ngext,noutpt,nptmax,nstmax,nttyo,xgexw,zchar)

    ! Compute various secondary parameters at the current point of
    ! reaction progress.
    call cdappl(acflg,acfw,acfwlg,actlg,actw,actwlg,adwipp,afcnst,affpd,affsd,ah,ahrc,alk,alk1,alk2,alki,atwt,bdwipp,cdrsd,cess,conc,csts,ctb,cteaq,dvoso,dwoso,eh,ehfac,ehrc,eps100,farad,fdpe0,fdse0,fjest,fo2lg,fo2lrc,fxist,iaqsln,iemop0,iemos0,iern1,iern2,ifrn1,ifrn2,ilrn1,ilrn2,imrn1,imrn2,iopt,ixrn1,ixrn2,jcode,jcsort,jern1,jflag,jflagd,jgext,jpflag,jsflag,jssort,modr,moph,mophg,mophj,mopht,morr,mosp,mospg,mospj,mospt,mprph,mprsp,mrgers,mrmlra,mtb,mtbaq,mte,mteaq,mwtges,mwtrc,mwtsp,narn1,narn2,nat,nbasp,nbaspd,nbt,nchlor,ncmpe0,ncmpr,nct,ndrsd,ndrsrd,nelect,nern1,nern2,nert,ness,nessr,net,nfrn1,nfrn2,ngext,ngrn1,ngrn2,ngt,nhydr,nlrn1,nlrn2,nlt,nmrn1,nmrn2,nmrt,nmt,no2gaq,npchk,npet,npet0,npt,npts,nrct,nrndex,nst,nsts,nstsr,ntf1,ntf1t,ntf2,ntf2t,nxridx,nxrn1,nxrn2,nxrt,nxt,osc,oscst,omega,pe,perc,ph,phmes,ppmwb,ppmwe,qriinf,qxknph,rreacn,rreac1,rxbar,sfcar,sidrsp,sidrph,sigmam,sigmst,tdays,tempc,tf1,tf2,thours,time1,tmins,tyears,uphase,uspec,vodrt,voph,vophg,vophj,vopht,vosoct,vosp,vospg,vospj,vospt,vosp0,vreac,wfh2o,wkgwi,wodr,wodrt,woph,wophg,wophj,wopht,worr,worrt,wosoct,wosp,wospg,wospj,wospt,xbar,xlke,xlksd,zchcu6,zchsq2)

    ! Print a description of the system at the current point of
    ! reaction progress.
    call scripz(abar,acflg,acfw,acfwlg,actlg,actw,actwlg,affpd,affsd,afrc1,aft1,ah,ahmes,ahnbs,ahrc,alki,alk1,alk2,awmax,awmin,a3bar,cbsr,cdrsd,cegexs,cesr,conc,conclg,csts,ctb,cteaq,dvoso,dwoso,egers,egexjc,egexjf,egexpa,egexpc,egexs,egexw,eh,ehmax,ehmes,ehmin,ehnbs,ehrc,elecsr,electr,fje,fjest,fo2,fo2lg,fo2lrc,fugac,fugalg,fxi,fxist,iaqsln,iemop,iemop0,iemos,iemos0,iern1,iern2,iexr,iexrt,ifrn1,ifrn2,ilrn1,ilrn2,imech,imrn1,imrn2,iopg,iopr,iopt,ipndx1,ixrn1,ixrn2,jcode,jcsort,jern1,jern2,jexr,jexrt,jflag,jflagd,jflgi,jgext,jgsort,jpflag,jreac,jsca,jscat,jscr,jscrt,jsflag,jsol,jssort,kbt,kern1,kern2,kgexsa,km1,kmt,kx1,kxt,kstep,kstpmx,loph,losp,mlmrra,modr,moph,mophg,mophj,mopht,morr,mosp,mospg,mospj,mospt,mprph,mprsp,mrgers,mrmlra,mwtrc,mwtsp,narn1,narn2,nat,nbasp,nbaspd,nbt,ncmpe,ncmpe0,ncmpr,nct,ndrsd,ndrsrd,nelect,nern1,nern2,nert,net,nfrn1,nfrn2,ngext,ngexsa,ngrn1,ngrn2,ngt,nhydr,nhydx,nlrn1,nlrn2,nlt,nmrn1,nmrn2,nmrt,nmt,noutpt,no2gaq,npet,npet0,npt,npts,nrct,nrdxsp,nrk,nrndex,nst,nsts,nstsr,ntf1t,ntf2t,nxridx,nxrn1,nxrn2,nxrt,nxt,osc,oscst,omega,o2max,o2min,pch,pe,pemes,penbs,perc,ph,phcl,phmax,phmes,phmin,phnbs,ppmwe,presg,press,qaft1,qftpr2,qmod,qphcl,qredox,qrho,qriinf,qstopx,qvhfxi,qvlsow,qzprnt,rho,rhoc,rhowc,rk,rreacn,rreac1,rrelr1,rxbar,sfcar,sidrph,sidrsp,sigmst,sigmam,ssfcar,tdays,tdsglw,tdspkc,tdsplc,tempc,thours,time1,timemx,tmins,tolsat,tolxsf,tolxst,tolxsu,tyears,uelem,ugermo,ugexj,ugexmo,uphase,ureac,uspec,uxtype,vodrt,voph,vophg,vophj,vopht,vosoct,vosol,vosp,vospg,vospj,vospt,vreac,wfh2o,wftds,wkgwi,woh2o,wodr,wodrt,woph,wophg,wophj,wopht,worr,worrt,wosoct,wosol,wosp,wospg,wospj,wospt,wotds,xbar,xbarlg,xbarw,xbrwlg,xgers,xgexw,xi1,xidump,ximax,xistsv,xirct,zchar)

    ! Write results at the current point of reaction progress on the
    ! tabx file. This file is used to create the plot file tab.
    if (ximax.gt.xistsv .and. iopt(18).ge.0) then
        if (qtatxt) then
            ! The TAB file is an ordinary text file.
            call wrtabx(actlg,afrc1,aft1,alk,cteaq,dvoso,dwoso,eh,fo2lg,iindx1,iktmax,iopt,ipndx1,kmax,km1,kmt,kstep,kx1,kxt,loph,ncmpr,modr,mopht,narn1,mosp,nct,nctmax,noptmx,nptmax,nrct,nrctmx,nstmax,ntabx,ntidmx,ntitl2,ntitld,ntitmx,nxtmax,pe,ph,ppmwe,prcinf,press,prminf,qbye,qmod,qriinf,tempc,time1,uelem,uphase,uplatm,ureac,uspec,usteq6,utitl2,utitld,uveeq6,vodrt,vosoct,wodrt,woh2o,wosoct,xbar,xi1)
        else
            ! The TAB file is a .csv file.
            call wrtabc(acflg,actlg,actw,afrc1,aft1,alk,conclg,cteaq,ctb,dvoso,dwoso,eh,fje,fo2lg,fugac,fxi,iktmax,iopt,jflag,jsflag,kmax,kstep,kx1,kxt,mrmlra,modr,mosp,mospt,moph,mopht,narn1,narn2,nbasp,nbaspd,nbt,nbtmax,ncmpr,nct,nctmax,nelect,ngrn1,ngrn2,ngtmax,nhydr,nhydx,nllnmx,no2gaq,noptmx,noutpt,npt,nptmax,nrct,nrctmx,nstmax,ntabx,ntidmx,ntitl2,ntitld,ntitmx,nttyo,nxrn1,nxrn2,nxtmax,pe,ph,phmes,ppmwb,ppmwe,prcinf,press,prminf,qrho,qriinf,rho,rhowc,sidrph,sigmam,tdsgks,tdsglw,tempc,time1,uelem,ulinex,uphase,uplatm,ureac,uspec,usteq6,utitl2,utitld,uveeq6,vodrt,vosoct,wkgh2o,wodrt,wosoct,xbar,xbarlg,xi1)
        end if
    end if

    if (iopt(16) .ge. 0) then
        ! Prepare to write results at the current point to the
        ! backup file.
        call setpk6(actwlg,awmax,awmaxi,awmin,awmini,eh,ehmax,ehmaxi,ehmin,ehmini,fo2lg,iindx1,jflag,jflgi,kbt,kdim,kmax,kprs,mprph,mprphi,mprsp,mprspi,mtb,mtbi,mtbaq,mtbaqi,nbasp,nbaspd,nbaspi,nbti,nbtmax,ncmpr,nobswt,noutpt,nprpmx,nprpti,nprsmx,nprsti,npt,nptmax,nttyo,nstmax,o2max,o2maxi,o2min,o2mini,ph,phmax,phmaxi,phmin,phmini,prcinf,press,pressi,tempc,tempci,time1,timemx,timmxi,tistti,ubmtbi,uobsw,uphase,uprphi,uprspi,uspec,uzveci,uzvec1,xi1,ximax,ximaxi,xistti,zvclgi,zvclg1)

        ! Rewind the current backup file (BAKUPA or BAKUPB).
        if (iopt(16) .eq. 0) then
            rewind (nbkupn)
        end if

        ! Write results at the current point to the backup file.
        if (upkfor(1:1) .eq. 'W') then
            ! Compact (W) format.
            ! Calling sequence substitutions:
            !   nbkupn for newin
            call wr6pkw(awmaxi,awmini,cbsri,cdac,cesri,cgexj,csigma,dlaplo,dlaprn,dleplo,dleprn,dlhplo,dlhprn,dloplo,dloprn,dltpll,dltplo,dltprl,dltprn,dlxdmp,dlxmx0,dlxpll,dlxplo,dlxprl,dlxprn,eact,egersi,ehmaxi,ehmini,electr,fkrc,iact,ibsrti,iesrti,ietmax,igerti,iktmax,imchmx,imech,iodb,iopg,iopr,iopt,itermx,ixrti,jcode,jetmax,jflgi,jgerti,jgext,jpress,jreac,jtemp,kbt,kct,kdim,kmax,kmt,kprs,ksplmx,ksppmx,kstpmx,kxmod,kxt,hact,modr,moffg,morr,mprphi,mprspi,mtbaqi,mtbi,mwtges,nbti,nbtmax,nbt1mx,nctmax,ndact,ndctmx,nert,nertmx,net,netmax,nbkupn,nffg,nffgmx,ngexrt,nobswt,nodbmx,nopgmx,noprmx,noptmx,nordmx,nprpmx,nprpti,nprsmx,nprsti,nptkmx,nrct,nrctmx,nrk,nsbswt,nsk,nsrt,nsrtmx,ntitl1,ntitl2,ntitmx,ntrymx,nttkmx,nxmdmx,nxmod,nxopex,nxopmx,nxopt,nxpemx,nxrt,nxrtmx,o2maxi,o2mini,phmaxi,phmini,pressb,pressi,ptk,qgexsh,rkb,rxbari,sfcar,ssfcar,tempcb,tempci,tgexp,timmxi,tistti,tolbt,toldl,tolsat,tolxsf,trkb,ttk,ubmtbi,ubsri,ucxri,udac,uesri,uffg,ugerji,ugermo,ugersi,ugexj,ugexmo,ugexp,ugexr,uhfgex,uobsw,uprphi,uprspi,ureac,usbsw,utitl1,utitl2,uvfgex,uxcat,uxkgex,uxmod,uxopex,uxopt,uzveci,vreac,xgersi,xhfgex,ximaxi,xistti,xlkffg,xlkgex,xlkmod,xvfgex,zgexj,zvclgi)
        else
            ! Menu-style (D) format.
            ! Calling sequence substitutions:
            !   nbkupn for newin
            call wr6pkd(awmaxi,awmini,cbsri,cdac,cesri,cgexj,csigma,dlaplo,dlaprn,dleplo,dleprn,dlhplo,dlhprn,dloplo,dloprn,dltpll,dltplo,dltprl,dltprn,dlxdmp,dlxmx0,dlxpll,dlxplo,dlxprl,dlxprn,eact,egersi,ehmaxi,ehmini,electr,fkrc,iact,ibsrti,iesrti,ietmax,igerti,iktmax,imchmx,imech,iodb,iopg,iopr,iopt,itermx,ixrti,jcode,jetmax,jflgi,jgerti,jgext,jpress,jreac,jtemp,kbt,kct,kdim,kmax,kmt,kprs,ksplmx,ksppmx,kstpmx,kxmod,kxt,hact,modr,moffg,morr,mprphi,mprspi,mtbaqi,mtbi,mwtges,nbti,nbtmax,nbt1mx,nctmax,ndact,ndctmx,nert,nertmx,net,netmax,nbkupn,nffg,nffgmx,ngexrt,nobswt,nodbmx,nopgmx,noprmx,noptmx,nordmx,nprpmx,nprpti,nprsmx,nprsti,nptkmx,nrct,nrctmx,nrk,nsbswt,nsk,nsrt,nsrtmx,ntitl1,ntitl2,ntitmx,ntrymx,nttkmx,nxmdmx,nxmod,nxopex,nxopmx,nxopt,nxpemx,nxrt,nxrtmx,o2maxi,o2mini,phmaxi,phmini,pressb,pressi,ptk,qgexsh,rkb,rxbari,sfcar,ssfcar,tempcb,tempci,tgexp,timmxi,tistti,tolbt,toldl,tolsat,tolxsf,trkb,ttk,ubmtbi,ubsri,ucxri,udac,uesri,uffg,ugerji,ugermo,ugersi,ugexj,ugexmo,ugexp,ugexr,uhfgex,uobsw,uprphi,uprspi,ureac,usbsw,utitl1,utitl2,uvfgex,uxcat,uxkgex,uxmod,uxopex,uxopt,uzveci,vreac,xgersi,xhfgex,ximaxi,xistti,xlkffg,xlkgex,xlkmod,xvfgex,zgexj,zvclgi)
        end if

        if (iopt(16) .eq. 0) then
            ! Switch to write on the other backup file the next time.
            if (nbkupn .eq. nbkupa) then
                nbkupn = nbkupb
            else if (nbkupn .eq. nbkupb) then
                nbkupn = nbkupa
            end if
        end if
    end if

    qdmpr1 = .false.

    if (qftpr2) then
        qdmpr2 = .false.
    end if

    kstppr = 0

    ! Plot a description of the system at the current point of
    ! reaction progress.
360 continue

    kstppl = 0

    ! 360 if (iplot .ge. 1) then
    !       Note: grafz has been removed from EQ6. Plotting is now done
    !       via the tab file
    !       call grafz
    !       kstppl = 0
    !     endif
    !     Check stop conditions.
    call chkstc(actw,awmax,awmin,eh,ehmax,ehmin,fo2lg,iopt,jreac,kstep,kstpmx,noptmx,noutpt,nrct,nrctmx,nttyo,o2max,o2min,ph,phmax,phmin,prcinf,qaft1,qcnpre,qcntmp,qconst,qredox,qstop,qvhfxi,qvlsow,timemx,time1,tolxst,tolxsu,ximax,xi1)

    if (qstop .or. qstopx) then
        go to 990
    end if

    ! Reset print, plot, and PRS transfer points as necessary.
    call adprpl(actw,aw0plo,aw0prn,aw1plo,aw1prn,dlaplo,dlaprn,dleplo,dleprn,dlhplo,dlhprn,dloplo,dloprn,dltpll,dltplo,dltprl,dltprn,dlxpll,dlxplo,dlxprl,dlxprn,dlxdmp,eh,eh0plo,eh0prn,eh1plo,eh1prn,fo2lg,fxprpl,o20plo,o20prn,o21plo,o21prn,ph,ph0plo,ph0prn,ph1plo,ph1prn,qplaw0,qplaw1,qpleh0,qpleh1,qplolt,qplolx,qplott,qplotx,qplo20,qplo21,qplph0,qplph1,qpraw0,qpraw1,qpreh0,qpreh1,qprnlx,qprnlt,qprntx,qprntt,qpro20,qpro21,qprph0,qprph1,qredox,time1,tiplol,tiplot,tiprnl,tiprnt,xidump,xiplol,xiplot,xiprnl,xiprnt,xi1)

    ! If shooting to the end, go do the next step.
    if (qshoot) then
        nord = 0
        delxi = dlxmin
        go to 100
    end if

    ! Continue the reaction path simulation. If the current point
    ! was hit only because of a print or plot requirement, and delxi
    ! is very small relative to the preceding step size, the current
    ! point is dropped as the base point (the point about which the
    ! Taylor's series expansion is made) in favor of the preceding
    ! point. This is done to protect the integrity of the finite
    ! difference functions, by avoiding the case in which two points
    ! are nearly indistinguishable. Here qskip is a logical flag to
    ! skip recalculating the finite differences and derivatives.
    qreax = iexrt.gt.0 .or. jexrt.gt.0 .or. jscat.gt.0 .or. jscrt.gt.0 .or. qreq
    qtrch = .false.
    iexrt = 0
    jexrt = 0
    jscat = 0
    jscrt = 0
    qskip = .false.

    if (qzprnt .or. qzplot) then
        if (qdump) then
            go to 100
        end if

        if (kord .le. 0) then
            go to 100
        end if

        if (delxi .ge. 0.01*(xi0 - xim1)) then
            go to 100
        end if
    else
        go to 100
    end if

    qskip = .true.
    nord = nordsv
    kord = kordsv
    delxi = max(delxi,dlxisv)
    xi1 = xi0
    time1 = time0

    km1 = km10
    kmt = kmt0
    kx1 = kx10
    kxt = kxt0
    kdim = kdim0

    call copyaa(zvclg0,zvclg1,kdim)
    call copyaa(zvec0,zvec1,kdim)

    do kcol = 1,kdim
        uzvec1(kcol) = uzvec0(kcol)
        iindx1(kcol) = iindx0(kcol)
        ipndx1(kcol) = ipndx0(kcol)
    end do

    if (iopt(2) .gt. 0) then
        rirec1 = rirec0

        call copyaa(rreac0,rreac1,nrct)
        call copyaa(rrelr0,rrelr1,nrct)
        call copyaa(sfcar0,sfcar,nrct)
    end if

    call copyaa(afrc0,afrc1,nrct)

    go to 110

990 continue
    if (kstpmx .le. 0) then
        ! Normal termination at the initial point.
        write (noutpt,2200)
        write (nttyo,2200)
2200 format(/,' ---  The reaction path terminated normally at',' the initial point ---',//)

        write (noutpt,2210) xistrt
        write (nttyo,2210) xistrt
2210 format(7x,'Xi is ',1pe12.5)

        write (noutpt,2220) kdim
        write (nttyo,2220) kdim
2220 format(7x,'The matrix dimension is ',i4,/)

        if (iopt(2) .gt. 0) then
            tistrd = tistrt/86400.
            tistry = tistrd/365.25
            write (noutpt,2230) tistrt,tistrd,tistry
            write (nttyo,2230) tistrt,tistrd,tistry
2230 format(7x,'The time is:',/14x,1pe12.5,' seconds',/14x,1pe12.5,' days',/14x,1pe12.5,' years',/)

            if (qriinf) then
                write (noutpt,2240)
                write (nttyo,2240)
2240 format(7x,'The system is now indistinguishable from',/7x,'what it would be at time equals infinity.')
            end if
        end if
    else
        if (.not.qstopx) then
            ! Normal termination of a reaction path calculation.
            write (noutpt,2300)
            write (nttyo,2300)
2300 format(/,' ---  The reaction path has terminated',' normally ---',//)
        else
            ! Early termination of a reaction path calculation.
            write (noutpt,2310)
            write (nttyo,2310)
2310 format(/,' ---  The reaction path has terminated',' early ---',//)
        end if

        write (noutpt,2320) kstep,xistsv,xi1
        write (nttyo,2320) kstep,xistsv,xi1
2320 format(7x,i5,' steps were taken',/7x,'Xi increased from ',/14x,1pe12.5,' to ',1pe12.5)

        avdlxi = 0.

        if (kstep .gt. 0) then
            avdlxi = (xi1 - xistsv)/kstep
        end if

        iavkdm = nint(avkdim)
        write (noutpt,2330) avdlxi,iavkdm
        write (nttyo,2330) avdlxi,iavkdm
2330 format(7x,'The average value of delxi was ',1pe12.5,/7x,'The average matrix dimension was ',i4,/)

        if (iopt(2) .gt. 0) then
            tistrd = tistrt/86400.
            tistry = tistrd/365.25

            if (qriinf) then
                write (noutpt,2340) tistrt,tistrd,tistry
                write (nttyo,2340) tistrt,tistrd,tistry
2340 format(7x,'The time increased from ',/14x,1pe12.5,' seconds to infinity',/14x,1pe12.5,' days to infinity',/14x,1pe12.5,' years to infinity',/)
            else
                tdays = time1/86400.
                tyears = tdays/365.25
                write (noutpt,2350) tistrt,time1,tistrd,tdays,tistry,tyears
                write (nttyo,2350) tistrt,time1,tistrd,tdays,tistry,tyears
2350 format(7x,'The time increased from ',/14x,1pe12.5,' to ',1pe12.5,' seconds',/14x,1pe12.5,' to ',1pe12.5,' days',/14x,1pe12.5,' to ',1pe12.5,' years',/)
            end if
        end if
    end if

    write (ux16a,'(f9.4)') phstrt
    call lejust(ux16a)
    j2 = ilnobl(ux16a)
    write (ux16b,'(f9.4)') ph
    call lejust(ux16b)
    j3 = ilnobl(ux16b)

    if (ph .gt. phstrt) then
        write (noutpt,2360) ux16a(1:j2),ux16b(1:j3)
        write (nttyo,2360) ux16a(1:j2),ux16b(1:j3)
2360 format(7x,'The pH increased from ',a,' to ',a)
    else if (ph .lt. phstrt) then
        write (noutpt,2362) ux16a(1:j2),ux16b(1:j3)
        write (nttyo,2362) ux16a(1:j2),ux16b(1:j3)
2362 format(7x,'The pH decreased from ',a,' to ',a)
    end if

    if (qredox) then
        write (ux16a,'(f9.4)') ehstrt
        call lejust(ux16a)
        j2 = ilnobl(ux16a)
        write (ux16b,'(f9.4)') eh
        call lejust(ux16b)
        j3 = ilnobl(ux16b)

        if (eh .gt. ehstrt) then
            write (noutpt,2370) ux16a(1:j2),ux16b(1:j3)
            write (nttyo,2370) ux16a(1:j2),ux16b(1:j3)
2370 format(7x,'The Eh increased from ',a,' to ',a,' v')
        else if (eh .lt. ehstrt) then
            write (noutpt,2372) ux16a(1:j2),ux16b(1:j3)
            write (nttyo,2372) ux16a(1:j2),ux16b(1:j3)
2372 format(7x,'The Eh decreased from ',a,' to ',a,' v')
        end if

        write (ux16a,'(f9.4)') o2strt
        call lejust(ux16a)
        j2 = ilnobl(ux16a)
        write (ux16b,'(f9.4)') fo2lg
        call lejust(ux16b)
        j3 = ilnobl(ux16b)

        if (fo2lg .gt. o2strt) then
            write (noutpt,2380) ux16a(1:j2),ux16b(1:j3)
            write (nttyo,2380) ux16a(1:j2),ux16b(1:j3)
2380 format(7x,'The log fO2 increased from ',a,' to ',a)
        else if (eh .lt. ehstrt) then
            write (noutpt,2382) ux16a(1:j2),ux16b(1:j3)
            write (nttyo,2382) ux16a(1:j2),ux16b(1:j3)
2382 format(7x,'The log fO2 decreased from ',a,' to ',a)
        end if
    end if

    write (ux16a,'(f9.4)') awstrt
    call lejust(ux16a)
    j2 = ilnobl(ux16a)
    write (ux16b,'(f9.4)') actw
    call lejust(ux16b)
    j3 = ilnobl(ux16b)

    if (actw .gt. awstrt) then
        write (noutpt,2390) ux16a(1:j2),ux16b(1:j3)
        write (nttyo,2390) ux16a(1:j2),ux16b(1:j3)
2390 format(7x,'The aw increased from ',a,' to ',a)
    else if (actw .lt. awstrt) then
        write (noutpt,2392) ux16a(1:j2),ux16b(1:j3)
        write (nttyo,2392) ux16a(1:j2),ux16b(1:j3)
2392 format(7x,'The aw decreased from ',a,' to ',a)
    end if

    write (ux16a,'(1pg12.5)') wwstrt
    call lejust(ux16a)
    j2 = ilnobl(ux16a)
    write (ux16b,'(1pg12.5)') wkgh2o
    call lejust(ux16b)
    j3 = ilnobl(ux16b)

    if (wkgh2o .gt. wwstrt) then
        write (noutpt,2394) ux16a(1:j2),ux16b(1:j3)
        write (nttyo,2394) ux16a(1:j2),ux16b(1:j3)
2394 format(7x,'The mass of solvent water increased from ',a,' to ',a,' kg')
    else if (wkgh2o .lt. wwstrt) then
        write (noutpt,2396) ux16a(1:j2),ux16b(1:j3)
        write (nttyo,2396) ux16a(1:j2),ux16b(1:j3)
2396 format(7x,'The mass of solvent water decreased from ',a,' to ',a,' kg')
    end if

    if (.not.qstopx) then
        if (iopt(7).gt.0 .and. iopt(1).ne.2) then
            ! Clear ES solids at the end of the run, unless the run
            ! terminated early due to calculational problems.
            ! Fictive fugacity-fixing minerals are not cleared.
            write (noutpt,2400)
            write (nttyo,2400)
2400 format(/' * Note (EQ6/path) Clearing equilibrium system',' (ES) solids',/7x,'at the end of the run.')

            call clress(csts,iindx1,ipndx1,jpflag,jsflag,kdim,kmax,km1,kmt,kx1,kxt,loph,losp,moph,mosp,mtb,mtbaq,nbt,nbtmax,nptmax,nstmax,nsts,nstsmx,nstsr,ufixf,uzvec1,zvec1,zvclg1)
        end if

        if (iopt(10) .gt. 0) then
            ! Clear PRS solids at the end of the run, unless the run
            ! terminated early due to calculational problems.
            write (noutpt,2410)
            write (nttyo,2410)
2410 format(/' * Note (EQ6/path) Clearing physically removed',' system (PRS) solids',/7x,'at the end of the run.')

            call initaz(mprsp,nstmax)
            call initaz(mprph,nptmax)
        end if
    end if

    if (nprob .le. 1) then
        ! Calculate data to describe the aqueous solution in the first
        ! problem on the input file as a special reactant. This defines
        ! the so-called "Fluid 2" to be described on the EQ6 pickup
        ! file under the iopt(20) = 1 option.
        ureac1 = 'Fluid 2'

        ! Scale the composition and reaction so that one "mole" of the
        ! fluid contains 1 kg of solvent water.
        morrw1 = mosp(narn1)/omega
        wx = 1./morrw1

        do nc = 1,nct
            n = nc
            uesr1(n) = uelem(nc)
            cesr1(n) = wx*mteaq(nc)
        end do

        iesrt1 = nct

        n = 1
        ubsr1(n) = ureac1
        cbsr1(n) = -1.0

        do nb = 1,nbt
            ns = nbaspd(nb)

            if (jflag(ns) .lt. 30) then
                n = n + 1
                ubsr1(n) = uspec(ns)(1:24)
                cbsr1(n) = wx*mtbaq(nb)
            end if
        end do

        ibsrt1 = n
    end if

    if (iopt(17) .ge. 0) then
        ! Prepare to write results at the final point to the
        ! pickup file.
        call setpk6(actwlg,awmax,awmaxi,awmin,awmini,eh,ehmax,ehmaxi,ehmin,ehmini,fo2lg,iindx1,jflag,jflgi,kbt,kdim,kmax,kprs,mprph,mprphi,mprsp,mprspi,mtb,mtbi,mtbaq,mtbaqi,nbasp,nbaspd,nbaspi,nbti,nbtmax,ncmpr,nobswt,noutpt,nprpmx,nprpti,nprsmx,nprsti,npt,nptmax,nttyo,nstmax,o2max,o2maxi,o2min,o2mini,ph,phmax,phmaxi,phmin,phmini,prcinf,press,pressi,tempc,tempci,time1,timemx,timmxi,tistti,ubmtbi,uobsw,uphase,uprphi,uprspi,uspec,uzveci,uzvec1,xi1,ximax,ximaxi,xistti,zvclgi,zvclg1)

        if (iopt(20) .gt. 0) then
            ! Prepare to write an advanced EQ6 pickup file. There is
            ! currently only one option. It is determined by the iopt(20)
            ! option switch. See comments in EQ6/stpkmd.f.
            call stpkmd(cbsri,cbsr1,cdac,cesri,cesr1,csigma,eact,fkrc,hact,iact,ibsrti,ibsrt1,iesrti,iesrt1,iktmax,imchmx,imech,iopt,ixrti,jcode,jreac,modr,morr,morrw1,nbt1mx,nctmax,ndact,ndctmx,noptmx,noutpt,nprob,nrct,nrctmx,nrk,nsk,nsrt,nsrtmx,ntitl1,ntitmx,nttyo,nxrt,nxrtmx,rkb,rxbari,sfcar,ssfcar,trkb,ubsri,ubsr1,ucxri,udac,uesri,uesr1,ureac,ureac1,utitl1,vreac)
        end if

        ! Write results at the final point to the pickup file.
        if (upkfor(1:1) .eq. 'W') then
            call wr6pkw(awmaxi,awmini,cbsri,cdac,cesri,cgexj,csigma,dlaplo,dlaprn,dleplo,dleprn,dlhplo,dlhprn,dloplo,dloprn,dltpll,dltplo,dltprl,dltprn,dlxdmp,dlxmx0,dlxpll,dlxplo,dlxprl,dlxprn,eact,egersi,ehmaxi,ehmini,electr,fkrc,iact,ibsrti,iesrti,ietmax,igerti,iktmax,imchmx,imech,iodb,iopg,iopr,iopt,itermx,ixrti,jcode,jetmax,jflgi,jgerti,jgext,jpress,jreac,jtemp,kbt,kct,kdim,kmax,kmt,kprs,ksplmx,ksppmx,kstpmx,kxmod,kxt,hact,modr,moffg,morr,mprphi,mprspi,mtbaqi,mtbi,mwtges,nbti,nbtmax,nbt1mx,nctmax,ndact,ndctmx,nert,nertmx,net,netmax,newin,nffg,nffgmx,ngexrt,nobswt,nodbmx,nopgmx,noprmx,noptmx,nordmx,nprpmx,nprpti,nprsmx,nprsti,nptkmx,nrct,nrctmx,nrk,nsbswt,nsk,nsrt,nsrtmx,ntitl1,ntitl2,ntitmx,ntrymx,nttkmx,nxmdmx,nxmod,nxopex,nxopmx,nxopt,nxpemx,nxrt,nxrtmx,o2maxi,o2mini,phmaxi,phmini,pressb,pressi,ptk,qgexsh,rkb,rxbari,sfcar,ssfcar,tempcb,tempci,tgexp,timmxi,tistti,tolbt,toldl,tolsat,tolxsf,trkb,ttk,ubmtbi,ubsri,ucxri,udac,uesri,uffg,ugerji,ugermo,ugersi,ugexj,ugexmo,ugexp,ugexr,uhfgex,uobsw,uprphi,uprspi,ureac,usbsw,utitl1,utitl2,uvfgex,uxcat,uxkgex,uxmod,uxopex,uxopt,uzveci,vreac,xgersi,xhfgex,ximaxi,xistti,xlkffg,xlkgex,xlkmod,xvfgex,zgexj,zvclgi)
        else
            ! Menu-style (D) format.
            call wr6pkd(awmaxi,awmini,cbsri,cdac,cesri,cgexj,csigma,dlaplo,dlaprn,dleplo,dleprn,dlhplo,dlhprn,dloplo,dloprn,dltpll,dltplo,dltprl,dltprn,dlxdmp,dlxmx0,dlxpll,dlxplo,dlxprl,dlxprn,eact,egersi,ehmaxi,ehmini,electr,fkrc,iact,ibsrti,iesrti,ietmax,igerti,iktmax,imchmx,imech,iodb,iopg,iopr,iopt,itermx,ixrti,jcode,jetmax,jflgi,jgerti,jgext,jpress,jreac,jtemp,kbt,kct,kdim,kmax,kmt,kprs,ksplmx,ksppmx,kstpmx,kxmod,kxt,hact,modr,moffg,morr,mprphi,mprspi,mtbaqi,mtbi,mwtges,nbti,nbtmax,nbt1mx,nctmax,ndact,ndctmx,nert,nertmx,net,netmax,newin,nffg,nffgmx,ngexrt,nobswt,nodbmx,nopgmx,noprmx,noptmx,nordmx,nprpmx,nprpti,nprsmx,nprsti,nptkmx,nrct,nrctmx,nrk,nsbswt,nsk,nsrt,nsrtmx,ntitl1,ntitl2,ntitmx,ntrymx,nttkmx,nxmdmx,nxmod,nxopex,nxopmx,nxopt,nxpemx,nxrt,nxrtmx,o2maxi,o2mini,phmaxi,phmini,pressb,pressi,ptk,qgexsh,rkb,rxbari,sfcar,ssfcar,tempcb,tempci,tgexp,timmxi,tistti,tolbt,toldl,tolsat,tolxsf,trkb,ttk,ubmtbi,ubsri,ucxri,udac,uesri,uffg,ugerji,ugermo,ugersi,ugexj,ugexmo,ugexp,ugexr,uhfgex,uobsw,uprphi,uprspi,ureac,usbsw,utitl1,utitl2,uvfgex,uxcat,uxkgex,uxmod,uxopex,uxopt,uzveci,vreac,xgersi,xhfgex,ximaxi,xistti,xlkffg,xlkgex,xlkmod,xvfgex,zgexj,zvclgi)
        end if
    end if
end subroutine path
