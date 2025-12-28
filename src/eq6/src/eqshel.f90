subroutine eqshel(aadh,aadhh,aadhv,aamatr,aaphi,abar,abdh,abdhh,abdhv,abdot,abdoth,abdotv,acflg,acflgo,act,actlg,adadhh,adadhv,adbdhh,adbdhv,adbdth,adbdtv,adh,adhh,adhv,afcnst,affp,affs,alpha,al10,amtb,aphi,apx,avcnst,azero,a3bar,a3bars,bacfmx,bbig,bdh,bdhh,bdhv,bdot,bdoth,bdotv,beta,betamx,betao,bgamx,bneg,bpx,cbsr,cco2,cdac,cegexs,cesr,cess,cdrs,cdrsd,cdrsx,cdrtw,cdrw,cjbasp,cnufac,conc,conclg,cpgexs,cscale,csigma,csts,dadhh,dadhv,dbdhh,dbdhv,dbdth,dbdtv,deltim,delvco,delvec,delxi,dlogxw,dlxmin,drer0,drir0,dzvc0,d1zvc1,eact,egers,egexjc,egexjf,egexs,eh,ehfac,elecsr,eps100,farad,fje,fjeo,fkrc,fo2,fo2lg,fsort,fugac,fugalg,fxi,fxio,gmmatr,hact,iact,iapxt,ibpxt,ibswx,ielam,ier,iern1,iern2,ifcphi1,ifcphi2,ifnnn,ifn2n,ifpsi1,ifpsi2,ifrn1,ifrn2,ifzeta,igas,igstak,iindx0,iindx1,ilcphi1,ilcphi2,ilnnn,iln2n,ilpsi1,ilpsi2,ilrn1,ilrn2,ilzeta,imech,imrn1,imrn2,insgf,iodb,iopg,iopt,ipch,ipivot,ipndx1,ipcv,istack,iter,itermx,ixbasp,ixrn1,ixrn2,izmax,jcode,jcsort,jflag,jgsort,jgstak,jjsort,jpflag,jpress,jptffl,jreac,jsflag,jsitex,jsol,jssort,jstack,jtemp,kbt,kction,kdim,kelect,khydr,khydx,km10,km1,kmt,kmt0,ko2gaq,kpsat,kpsst,krdxsp,kwater,kx1,kx10,kxt,kxt0,loph,losp,lsort,modr,modr0,moph,morr,morr0,mosp,mrgers,mrgexs,mtb,mtbaq,mtb0,mte,mteaq,mwtrc,narn1,narn2,narxt,nat,nbasp,nbaspd,nbaspx,nbt,nbtd,nbw,nchlor,ncmpr,ncorr,nct,ndac,ndact,ndrs,ndrsd,ndrsx,ndrsr,ndrsrd,ndrsrx,nelect,nern1,nern2,ness,nessr,net,nfrn1,nfrn2,ngrn1,ngrn2,ngt,nhydr,nhydx,nlrn1,nlrn2,nlt,nmrn1,nmrn2,nmt,nord,noutpt,no2gaq,npchk,nphasx,npslmx,npt,nrct,nrdxsp,nrk,nrndex,nsk,nsslmx,nst,nsts,nstsr,ntpr,ntrymx,ntprt,nttyo,nxridx,nxrn1,nxrn2,nxt,nweope,nwndpc,omega,omeglg,prcinf,press,pressb,pressd,ptk,qbassw,qbseqc,qbye,qcnpre,qcntmp,qhawep,qmod,qoptmz,qpit75,qredox,qriinf,qscon,qshoot,qstart,qtrch,qtvchk,qxknph,q6mode,rcnstv,rconst,rkb,rhsvec,rirec0,rk,rreacn,rreac1,rrelr0,rrelr1,rtcnst,rxbar,screwd,sidrph,sidrsp,sigmam,sigmmo,smp100,tdays,tempc,tempcb,tempcd,tempcu,tempc0,tempk,timemx,time0,time1,tiplol,tiplot,tiprnl,tiprnt,tistsv,tolbt,toldl,tolsat,tolsst,tolxst,trkb,ttk,ubacmx,ubgamx,udac,ufixf,ugermo,ulbeta,uldel,uphase,ureac,uspec,uzvec0,uzvec1,vreac,weight,wfac,wodr,worr,xbar,xbarlg,xbarw,xbarwc,xbrwlc,xbrwlg,xgers,xirct,xirct0,xistsv,xi0,xi1,zchar,zchcu6,zchsq2,zklogu,zvclg0,zvclg1,zvec0,zvec1)
    !! This subroutine provides a shell around EQ6/eqcalc.f that allows
    !! fine adjustments to the value of reaction progress in order to
    !! step over fine singularities, such as may exist in the middle of
    !! a major jump in the oxygen fugacity.
    !! This subroutine is called by:
    !!   EQ6/path.f
    !! Principal input:
    !! Principal output:
    !!   qbye comes in .true. if EQ6/path.f just changed the phase
    !!     assemblage by deleting one or more minerals. This
    !!     instructs this subroutine to print the index (iindx1)
    !!     structure of the system of equations as it does on the
    !!     starting call and whenever this subroutine itself adds or
    !!     deletes a phase.
    !!   iter is the number of Newton-Raphson iterations that
    !!     EQLIB/newton.f performed.
    !!   ier is an error parameter which is returned as 0 if the
    !!     calculations in this subroutine converged successfully and
    !!     there were no violations of possible open system constraints.
    !!   ier    = error flag:
    !!              Values returned from EQ6/eqphas.f:
    !!                =    0  Okay
    !!                =   10  Go back and take a smaller step size to
    !!                          avoid exceeding the supersaturation
    !!                          tolerance (tolsst). This is done only
    !!                          when an appropriate set of conditions
    !!                          is satisified. It isn't necessary for
    !!                          EQ6/path.f to analyze the situation
    !!                          when this value is returned to it by
    !!                          EQ6/eqshel.f.
    !!                =   20  Hit the maximum number of tries to find
    !!                          the correct phase assemblage
    !!                =   30  Caught in a region of computational
    !!                          instability about the phase boundary for
    !!                          an appearing phase. Iteration fails when
    !!                          the phase is added to the phase
    !!                          assemblage.
    !!                =   40  Caught in a region of computational
    !!                          instability about the phase boundary for
    !!                          a disappearing phase. The system is
    !!                          too supersaturated with respect to a
    !!                          phase if that phase is deleted from the
    !!                          phase assemblage.
    !!                =   50  Caught in a region of computational
    !!                          instability associated with solvent
    !!                          water. The amount of water in the
    !!                          system is probably very low.
    !!                =   60  Caught in a region of computational
    !!                          instability associated with the
    !!                          master redox variable. This is probably
    !!                          associated with a redox jump.
    !!                =   90  Caught in a region of computational
    !!                          instability associated with the
    !!                          aqueous activity coefficient model.
    !!                =  100  Detected out-of-range values for variables
    !!                          associated with the basis species before
    !!                          starting iteration.
    !!                =  110  Detected out-of-range values for variables
    !!                          associated with the basis species after
    !!                          an iteration crash.
    !!                =  150  Calculation failed, no diagnostics were
    !!                          generated.
    !!              Values returned by the present subroutine:
    !!                =    0  Okay
    !!                =   10  Go back and take a smaller step size to
    !!                          avoid exceeding the supersaturation
    !!                          tolerance (tolsst)
    !!                =  170  Too much of a phase was destroyed under
    !!                          the flow-through open system model;
    !!                          go back and first move part of the
    !!                          mass of protected phases in the ES to
    !!                          the PRS
    !!                =  180  One of a number of problems occurred which
    !!                          may be resolvable, at least partially,
    !!                          by going back and cutting the step size
    !!                =  190  Need to slide over a region of
    !!                          computational instability, but sliding
    !!                          is inhibited; go back, but terminate
    !!                          work on the current problem
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

    ! Calling sequence variable declarations.
    integer :: noutpt
    integer :: nttyo

    integer :: iact(imchmx,2,nrctmx)
    integer :: imech(2,nrctmx)
    integer :: ixbasp(nbtmax)
    integer :: ndac(ndctmx,imchmx,2,nrctmx)
    integer :: ndact(imchmx,2,nrctmx)
    integer :: nrk(2,nrctmx)
    integer :: nsk(nrctmx)

    integer :: jcode(nrctmx)
    integer :: jreac(nrctmx)
    integer :: nrndex(nrctmx)
    integer :: nxridx(nrctmx)

    integer :: iapxt(nxtmax)
    integer :: ibpxt(nxtmax)
    integer :: ibswx(nbtmax)
    integer :: igstak(ngtmax)
    integer :: iindx0(kmax)
    integer :: iindx1(kmax)
    integer :: ipivot(kmax)
    integer :: ipndx1(kmax)
    integer :: insgf(natmax)
    integer :: iodb(nodbmx)
    integer :: iopg(nopgmx)
    integer :: iopt(noptmx)
    integer :: istack(nstmax)
    integer :: jcsort(nstmax)
    integer :: jflag(nstmax)
    integer :: jgsort(ngtmax)
    integer :: jgstak(ngtmax)
    integer :: jjsort(nstmax)
    integer :: jpflag(nptmax)
    integer :: jsflag(nstmax)
    integer :: jsitex(nstmax)
    integer :: jsol(nxtmax)
    integer :: jssort(nstmax)
    integer :: jstack(nstmax)
    integer :: kction(nbtmax)

    integer :: narxt(ntprmx)
    integer :: nbasp(nbtmax)
    integer :: nbaspd(nbtmax)
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

    integer :: ielam
    integer :: ier
    integer :: igas
    integer :: ipch
    integer :: ipcv
    integer :: iter
    integer :: itermx
    integer :: izmax
    integer :: jpress
    integer :: jptffl
    integer :: jtemp
    integer :: kbt
    integer :: kdim
    integer :: kelect
    integer :: khydr
    integer :: khydx
    integer :: km1
    integer :: km10
    integer :: kmt
    integer :: kmt0
    integer :: ko2gaq
    integer :: kpsat
    integer :: kpsst
    integer :: krdxsp
    integer :: kwater
    integer :: kx1
    integer :: kx10
    integer :: kxt
    integer :: kxt0
    integer :: nbtd
    integer :: nbw
    integer :: nchlor
    integer :: ncorr
    integer :: nelect
    integer :: nhydr
    integer :: nhydx
    integer :: nord
    integer :: no2gaq
    integer :: npslmx
    integer :: nrct
    integer :: nrdxsp
    integer :: nsslmx
    integer :: ntpr
    integer :: ntprt
    integer :: ntrymx
    integer :: nweope
    integer :: nwndpc

    logical :: qxknph(nptmax)

    logical :: qbassw
    logical :: qbseqc
    logical :: qbye
    logical :: qcnpre
    logical :: qcntmp
    logical :: qhawep
    logical :: qmod
    logical :: qoptmz
    logical :: qpit75
    logical :: qredox
    logical :: qriinf
    logical :: qshoot
    logical :: qscon
    logical :: qstart
    logical :: qtrch
    logical :: qtvchk
    logical :: q6mode

    character(len=48) :: uspec(nstmax)
    character(len=48) :: uzvec0(kmax)
    character(len=48) :: uzvec1(kmax)
    character(len=48) :: ubacmx
    character(len=48) :: ubgamx
    character(len=24) :: ugermo(nertmx)
    character(len=24) :: ureac(nrctmx)
    character(len=24) :: udac(ndctmx,imchmx,2,nrctmx)
    character(len=24) :: uphase(nptmax)
    character(len=8) :: ulbeta(kmax)
    character(len=8) :: uldel(kmax)
    character(len=8) :: ufixf

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
    real(kind=8) :: rreacn(nrctmx)
    real(kind=8) :: rreac1(nrctmx)
    real(kind=8) :: rrelr1(nrctmx)
    real(kind=8) :: rxbar(iktmax,nxrtmx)
    real(kind=8) :: vreac(nrctmx)
    real(kind=8) :: wodr(nrctmx)
    real(kind=8) :: worr(nrctmx)
    real(kind=8) :: xgers(ietmax,jetmax,nertmx)

    real(kind=8) :: aamatr(kmax,kmax)
    real(kind=8) :: acflg(nstmax)
    real(kind=8) :: acflgo(nstmax)
    real(kind=8) :: act(nstmax)
    real(kind=8) :: actlg(nstmax)
    real(kind=8) :: affp(nptmax)
    real(kind=8) :: affs(nstmax)
    real(kind=8) :: alpha(kmax)
    real(kind=8) :: amtb(nbtmax)
    real(kind=8) :: apx(iapxmx,nxtmax)
    real(kind=8) :: azero(natmax)
    real(kind=8) :: a3bars(natmax)
    real(kind=8) :: beta(kmax)
    real(kind=8) :: betao(kmax)
    real(kind=8) :: bpx(ibpxmx,nxtmax)
    real(kind=8) :: cco2(5)
    real(kind=8) :: cegexs(ietmax,jetmax,netmax)
    real(kind=8) :: cess(nessmx)
    real(kind=8) :: cdrs(ndrsmx)
    real(kind=8) :: cdrsd(ndrsmx)
    real(kind=8) :: cdrsx(ndrsmx)
    real(kind=8) :: cdrtw(nstmax)
    real(kind=8) :: cdrw(nstmax)
    real(kind=8) :: cjbasp(nbtmax)
    real(kind=8) :: cnufac(nstmax)
    real(kind=8) :: conc(nstmax)
    real(kind=8) :: conclg(nstmax)
    real(kind=8) :: cpgexs(ietmax,jetmax,netmax)
    real(kind=8) :: cscale(nstmax)
    real(kind=8) :: csts(nstsmx)
    real(kind=8) :: drer0(nrd1mx,nrctmx)
    real(kind=8) :: delvco(kmax)
    real(kind=8) :: delvec(kmax)
    real(kind=8) :: dlogxw(nbtmax)
    real(kind=8) :: drir0(nrd1mx)
    real(kind=8) :: dzvc0(nrd1mx,kmax)
    real(kind=8) :: d1zvc1(kmax)
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
    real(kind=8) :: modr0(nrctmx)
    real(kind=8) :: moph(nptmax)
    real(kind=8) :: morr0(nrctmx)
    real(kind=8) :: mosp(nstmax)
    real(kind=8) :: mrgexs(ietmax,jetmax,netmax)
    real(kind=8) :: mtb(nbtmax)
    real(kind=8) :: mtbaq(nbtmax)
    real(kind=8) :: mtb0(nbtmax)
    real(kind=8) :: mte(nctmax)
    real(kind=8) :: mteaq(nctmax)
    real(kind=8) :: ptk(nttkmx)
    real(kind=8) :: rhsvec(kmax)
    real(kind=8) :: rk(imchmx,2,nrctmx)
    real(kind=8) :: rrelr0(nrctmx)
    real(kind=8) :: sidrph(nptmax)
    real(kind=8) :: sidrsp(nstmax)
    real(kind=8) :: tempcu(ntprmx)
    real(kind=8) :: ttk(nttkmx)
    real(kind=8) :: weight(nstmax)
    real(kind=8) :: wfac(iktmax,nxtmax)
    real(kind=8) :: xbar(nstmax)
    real(kind=8) :: xbarlg(nstmax)
    real(kind=8) :: xirct(nrctmx)
    real(kind=8) :: xirct0(nrctmx)
    real(kind=8) :: zchar(nstmax)
    real(kind=8) :: zchcu6(nstmax)
    real(kind=8) :: zchsq2(nstmax)
    real(kind=8) :: zvclg0(kmax)
    real(kind=8) :: zvclg1(kmax)
    real(kind=8) :: zvec0(kmax)
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
    real(kind=8) :: afcnst
    real(kind=8) :: al10
    real(kind=8) :: avcnst
    real(kind=8) :: a3bar
    real(kind=8) :: bacfmx
    real(kind=8) :: bbig
    real(kind=8) :: betamx
    real(kind=8) :: bgamx
    real(kind=8) :: bneg
    real(kind=8) :: deltim
    real(kind=8) :: delxi
    real(kind=8) :: dlxmin
    real(kind=8) :: eh
    real(kind=8) :: ehfac
    real(kind=8) :: eps100
    real(kind=8) :: farad
    real(kind=8) :: fje
    real(kind=8) :: fjeo
    real(kind=8) :: fo2
    real(kind=8) :: fo2lg
    real(kind=8) :: fxi
    real(kind=8) :: fxio
    real(kind=8) :: omega
    real(kind=8) :: omeglg
    real(kind=8) :: prcinf
    real(kind=8) :: press
    real(kind=8) :: pressb
    real(kind=8) :: pressd
    real(kind=8) :: rcnstv
    real(kind=8) :: rconst
    real(kind=8) :: rirec0
    real(kind=8) :: rtcnst
    real(kind=8) :: screwd
    real(kind=8) :: sigmam
    real(kind=8) :: sigmmo
    real(kind=8) :: smp100
    real(kind=8) :: tdays
    real(kind=8) :: tempc
    real(kind=8) :: tempcb
    real(kind=8) :: tempcd
    real(kind=8) :: tempk
    real(kind=8) :: tempc0
    real(kind=8) :: timemx
    real(kind=8) :: time0
    real(kind=8) :: time1
    real(kind=8) :: tiplol
    real(kind=8) :: tiplot
    real(kind=8) :: tiprnl
    real(kind=8) :: tiprnt
    real(kind=8) :: tistsv
    real(kind=8) :: tolbt
    real(kind=8) :: toldl
    real(kind=8) :: tolsat
    real(kind=8) :: tolsst
    real(kind=8) :: tolxst
    real(kind=8) :: xbarw
    real(kind=8) :: xbarwc
    real(kind=8) :: xbrwlc
    real(kind=8) :: xbrwlg
    real(kind=8) :: xistsv
    real(kind=8) :: xi0
    real(kind=8) :: xi1
    real(kind=8) :: zklogu

    ! Local variable declarations with global dimensioning.
    ! Variables for restoring the entering configuration, including
    ! the entering phase assemblage.
    integer :: isv_kmax
    integer :: isv_nbtmax
    integer :: isv_nstmax

    SAVE isv_kmax,isv_nbtmax,isv_nstmax

    integer, dimension(:), allocatable :: iindxs
    integer, dimension(:), allocatable :: ipndxs
    integer, dimension(:), allocatable :: nbasps

    SAVE iindxs,ipndxs,nbasps

    real(kind=8), dimension(:), allocatable :: acflgs
    real(kind=8), dimension(:), allocatable :: zvclgs

    SAVE acflgs,zvclgs

    ! The following do not need to be SAVEd.
    integer :: kdims
    integer :: km1s
    integer :: kmts
    integer :: kx1s
    integer :: kxts

    real(kind=8) :: xbarws
    real(kind=8) :: xbrwls

    ! Local variable declarations.
    integer :: jlen
    integer :: j2
    integer :: k
    integer :: kcol
    integer :: krow
    integer :: n
    integer :: nb
    integer :: nords
    integer :: np
    integer :: nphasl
    integer :: nplast
    integer :: np1
    integer :: nrdxsl
    integer :: ns
    integer :: ns2
    integer :: ntpr0

    integer :: ilnobl

    logical :: qbswok
    logical :: qslmod
    logical :: qsspgb
    logical :: qztayl

    character(len=56) :: uspn56

    real(kind=8) :: dlmoph
    real(kind=8) :: dltimd
    real(kind=8) :: lold
    real(kind=8) :: lxx
    real(kind=8) :: mold
    real(kind=8) :: mxx0

    real(kind=8) :: texp
    real(kind=8) :: tlg

    ! Allocate or reallocate local work arrays as needed.
    if (.not.ALLOCATED(iindxs)) then
        ! Local work arrays are not allocated. Zero the saved
        ! array size variables. Note that only one array is tested
        ! to see if it is allocated. It is assumed that all local
        ! work arrays are either allocated or not.
        isv_kmax = 0
        isv_nbtmax = 0
        isv_nstmax = 0
    else
        ! Local work arrays are allocated. Check to see if any of the
        ! array size variables have changed. If so, deallocate
        ! the corresponding local work arrays and zero the corresponding
        ! saved size variables.
        if (kmax .ne. isv_kmax) then
            DEALLOCATE(iindxs,ipndxs)
            DEALLOCATE(zvclgs)
            isv_kmax = 0
        end if

        if (nbtmax .ne. isv_nbtmax) then
            DEALLOCATE(nbasps)
            isv_nbtmax = 0
        end if

        if (nstmax .ne. isv_nstmax) then
            DEALLOCATE(acflgs)
            isv_nstmax = 0
        end if
    end if

    ! At this point, the saved array size values are zero if the
    ! corresponding arrays need to be allocated.
    if (isv_kmax .eq. 0) then
        ALLOCATE(iindxs(kmax),ipndxs(kmax))
        ALLOCATE(zvclgs(kmax))
        isv_kmax = kmax
    end if

    if (isv_nbtmax .eq. 0) then
        ALLOCATE(nbasps(nbtmax))
        isv_nbtmax = nbtmax
    end if

    if (isv_nstmax .eq. 0) then
        ALLOCATE(acflgs(nstmax))
        isv_nstmax = nstmax
    end if

    ! Zero the contents of the local work arrays.
    do k = 1,kmax
        iindxs(k) = 0
        ipndxs(k) = 0
    end do

    do k = 1,kmax
        zvclgs(k) = -99999.
    end do

    do n = 1,nbtmax
        nbasps(n) = 0
    end do

    do n = 1,nstmax
        acflgs(n) = 0.
    end do

    ! Initialize some parameters.
    qslmod = .false.
    qmod = .false.
    qbseqc = .false.

    ier = 0
    nrdxsl = 0
    nphasl = 0

    ! Save the "kernel" description for the current equilibrium system.
    ! This will be used to recover if the equilibrium calculation
    ! overseen by EQ6/eqphas.f fails for any reason. This kernel
    ! contains the minimum of information required to calculate the
    ! fully expanded description using EQLIB/ncmpex.f. Note that
    ! EQ6/eqphas.f (which the present subroutine calls) and EQ6/path.f
    ! (which calls the present subroutine) have their own distinct
    ! backup kernels.
    km1s = km1
    kmts = kmt
    kx1s = kx1
    kxts = kxt
    kdims = kdim

    call copyia(nbasp,nbasps,nbt)
    call copyia(iindx1,iindxs,kdim)
    call copyia(ipndx1,ipndxs,kdim)

    call copyaa(zvclg1,zvclgs,kdim)
    call copyaa(acflg,acflgs,nst)

    xbarws = xbarwc
    xbrwls = xbrwlc

    ! Save the order of the finite differences.
    nords = nord

100 continue

    ! This is a return point if the step size is changed within
    ! the present subroutine. There are two cases:
    !    1. Advance the step size through a region of critical
    !       redox instability (redox jump) from the last nonsingular
    !       point of reaction progress. Move the simulation over
    !       the region of singularity by using only predictor functions
    !       expanded from that last good point.
    !    2. Advance the step size across a region of singularity
    !       that is associated with the presence or absence of a phase.
    !       There are two cases. first, the calculations may fail to
    !       converge upon adding a newly supersaturated phase because
    !       the number of moles of the phase that would be produced is
    !       so small that the Jacobian matrix becomes effectively
    !       singular. Second, when a phase is disappearing, its number
    !       of moles may become so small that the above singlularity
    !       occurs. This is a normal mode for dropping a phase.
    !       Sometimes, the affinity of the phase exceeds the
    !       supersaturation tolerance when the system is solved for
    !       the case of the phase not present. In each case, the step
    !       size must be advanced across the region of apparent
    !       singularity.
    if (ncorr .le. 0) then
        write (noutpt,1000) xi1,delxi,nord
1000 format(/' Stepping to Xi= ',1pe11.4,', delxi= ',1pe11.4,', nord= ',i1)
    end if

    if (iopt(2) .gt. 0) then
        if (xi1 .gt. xistsv) then
            if (.not.qshoot) then
                call timeca(deltim,delxi,drir0,iodb,nodbmx,nord,noutpt,nrd1mx,nttyo,prcinf,qriinf,rirec0,time0,time1)

                if (delxi .le. dlxmin) then
                    ! Make sure that the calculated time does not exceed any
                    ! any specified limits such as the maximum time just because
                    ! delxi is at the minimum value.
                    call tivchk(deltim,delxi,qtvchk,time1,time0,timemx,tiplol,tiplot,tiprnl,tiprnt,tolxst)
                end if
            end if
        else
            time1 = tistsv
            deltim = 0.
        end if

        tdays = time1/86400.
        dltimd = deltim/86400.

        write (noutpt,1010) ncorr,tdays,dltimd
1010 format(3x,'ncorr= ',i2,', time= ',1pe11.4,' d, deltim= ',e11.4,' d')
    end if

    if (nord .gt. 0) then
        ! Make a Taylor's series expansion of the dz/d(xi) vector.
        ! This information is used to track how the reacting system
        ! is changing.
        call d1ztay(delxi,dzvc0,d1zvc1,kdim,kmax,nord,nrd1mx)
    end if

    if (nrct .ge. 1) then
        ! Increment the irreversible reactions (those associated with
        ! the "reactants"); update the ES mass balance totals accordingly.
        call reacts(cbsr,csts,delxi,drer0,iern1,ietmax,iktmax,iodb,jcode,jetmax,jgext,jreac,modr,modr0,morr,morr0,mrgers,mtb,mtb0,nbaspd,nbt,nbtmax,nbt1mx,ncmpr,nern1,nern2,nertmx,netmax,ngext,nodbmx,nord,noutpt,nptmax,nrct,nrctmx,nrd1mx,nrndex,nsrtmx,nstmax,nsts,nstsmx,nstsr,nttyo,nxridx,nxrtmx,rrelr0,rxbar,ureac,xirct,xirct0)
    end if

    ! Set a flag for a go back followed by a step size reduction when
    ! the equilibrium calculation indicates a supersaturation that
    ! exceeds the tolerance value used in locating the phase boundary
    ! for a newly appearing phase.
    qsspgb = iopt(3).le.1 .and. .not.qslmod .and.(.not.qscon) .and. delxi.gt.dlxmin .and. .not.qstart

    ! Make the equilibrium calculation. Basis switching may occur.
    ! Also, the phase assemblage may be changed.
    call eqphas(aamatr,abar,acflg,acflgo,act,actlg,adh,adhh,adhv,afcnst,affp,affs,alpha,al10,amtb,aphi,apx,avcnst,azero,a3bar,a3bars,bacfmx,bbig,bdh,bdhh,bdhv,bdot,bdoth,bdotv,beta,betamx,betao,bgamx,bneg,bpx,cco2,cegexs,cess,cdrs,cdrsd,cdrsx,cdrtw,cdrw,cjbasp,cnufac,conc,conclg,cpgexs,cscale,csts,delvco,delvec,d1zvc1,dlogxw,egexjc,egexjf,egexs,eh,ehfac,eps100,farad,fje,fjeo,fo2,fo2lg,fsort,fugac,fugalg,fxi,fxio,gmmatr,iapxt,ibpxt,ibswx,ielam,ier,iern1,iern2,ifcphi1,ifcphi2,ifnnn,ifn2n,ifpsi1,ifpsi2,ifrn1,ifrn2,ifzeta,igas,igstak,iindx0,iindx1,ilcphi1,ilcphi2,ilnnn,iln2n,ilpsi1,ilpsi2,ilrn1,ilrn2,ilzeta,imrn1,imrn2,insgf,iodb,iopg,iopt,ipch,ipivot,ipndx1,ipcv,istack,iter,itermx,ixbasp,ixrn1,ixrn2,izmax,jcsort,jflag,jgsort,jgstak,jjsort,jpflag,jsflag,jsitex,jsol,jssort,jstack,kbt,kction,kdim,kelect,khydr,khydx,km1,km10,kmt,kmt0,ko2gaq,kpsat,kpsst,krdxsp,kwater,kx1,kx10,kxt0,kxt,loph,losp,lsort,moph,mosp,mrgexs,mtb,mtbaq,narn1,narn2,narxt,nat,nbasp,nbaspd,nbaspx,nbt,nbtd,nbw,nchlor,ncmpr,nct,ndrs,ndrsd,ndrsx,ndrsr,ndrsrd,ndrsrx,nelect,nern1,nern2,ness,nessr,net,nfrn1,nfrn2,ngrn1,ngrn2,ngt,nhydr,nhydx,nlrn1,nlrn2,nlt,nmrn1,nmrn2,nmt,nord,no2gaq,noutpt,npchk,nphasx,npt,nrdxsp,nst,nsts,nstsr,ntpr,ntrymx,nttyo,nxrn1,nxrn2,nxt,omega,omeglg,prcinf,press,qbassw,qbseqc,qbye,qcnpre,qcntmp,qhawep,qmod,qoptmz,qpit75,qredox,qsspgb,qstart,qxknph,q6mode,rcnstv,rconst,rhsvec,rtcnst,screwd,sidrph,sidrsp,sigmam,sigmmo,smp100,tempc,tempk,tolbt,toldl,tolsat,tolsst,ubacmx,ubgamx,ulbeta,uldel,uphase,uspec,uzvec1,weight,wfac,xbar,xbarlg,xbarw,xbarwc,xbrwlc,xbrwlg,zchar,zchcu6,zchsq2,zvclg1,zvec1)

    if (ier .eq. 8) then
        ! Go back and reduce the step size to avoid exceeding the
        ! dimensioned limit on the iindx1 array (this also corresponds to
        ! the dimensions of the Jacobian matrix). EQ6/path.f will go back
        ! and cut the step size, if this is possible.
        ier = 180
        go to 999
    end if

    if (ier .eq. 10) then
        ! Go back and reduce the step size to avoid exceeding the
        ! supersaturation tolerance (tolsst). This is only done under
        ! an appropriate set of conditions, so no further analysis
        ! is required. EQ6/path.f will go back and cut the step size.
        go to 999
    end if

    if (ier .eq. 20) then
        ! EQ6/eqphas.f hit the maximum number of tries to find the correct
        ! phase without succeeding. Go back and let EQ6/path.f deal with
        ! the situation. It may be able to resolve the problem by reducing
        ! the step size.
        ier = 180
        go to 999
    end if

    if (ier.eq.30 .or. ier.eq.40) then
        ! The calculation appears to be in a region of critical
        ! instability in the ES phase assemblage (can't solve to within
        ! specified tolerances with a certain phase or phases in or out
        ! of the assemblage). Try to slide over it.
        if (nphasl .le. 0) then
            ! Am not currently engaged in such a slide. Need to start
            ! the process.
            if (npslmx .le. 0) then
                ! This kind of sliding is inhibited.
                write (noutpt,1100)
                write (nttyo,1100)
1100 format(/' * Warning - (EQ6/eqshel) Need to slide Xi',' forward to get over a region of',/7x,'critical',' instability in the ES phase assemblage. However, the',' maximum',/7x,'number of steps for this process (npslmx)',' is zero. Sliding of this kind',/7x,'is therefore',' inhibited.')

                ier = 190
                go to 999
            else
                ! Set up to start sliding.
                write (noutpt,1110)
1110 format(/' Setting up to slide Xi forward to get over a',' region of critical'/3x,'instability in the ES phase',' assemblage.')

                qslmod = .true.
            end if
        end if

        if (nphasl .ge. npslmx) then
            ! Check for the maximum number of tries sliding forward.
            write (noutpt,1120)
            write (nttyo,1120)
1120 format(/' * Warning - (EQ6/eqshel) Have done the maximum',' number of tries to',/7x,'slide over a region of critical',' instability in the ES phase assemblage.')

            ier = 190
            go to 999
        end if

        ! Slide Xi forward one increment. The increment gets progressively
        ! larger with the number of tries.
        nphasl = nphasl + 1
        delxi = delxi + 3.**(nphasl - 1)*dlxmin
        write (noutpt,1130)
1130 format(/' Trying to slide over a region of critical',' instability',/7x,'in the ES phase assemblage:')

        write (noutpt,1140) nphasl,delxi
1140 format(/3x,'Try ',i2,': delxi= ',1pe12.5,/)

        go to 300
    end if

    if (ier .eq. 50) then
        ! The calculation appears to be in a region of critical redox
        ! instability. Try to slide over it.
        if (nrdxsl .le. 0) then
            ! Am not currently engaged in such a slide. Need to start
            ! the process.
            if (nsslmx .le. 0) then
                ! This kind of sliding is inhibited.
                write (noutpt,1200)
                write (nttyo,1200)
1200 format(/' * Warning - (EQ6/eqshel) Need to slide Xi',' forward to get over a region of',/7x,'critical redox',' instability. However, the maximum number of steps',/7x,'for this process (nsslmx) is zero. Sliding of this',' kind is',/7x,'therefore inhibited.')

                ier = 190
                go to 999
            else
                ! Set up to start sliding.
                write (noutpt,1210)
1210 format(/' Setting up to slide Xi forward to get over a',' region of critical'/3x,'redox instability.')

                qslmod = .true.
            end if
        end if

        if (nrdxsl .ge. nsslmx) then
            ! Check for the maximum number of tries sliding forward.
            write (noutpt,1220)
            write (nttyo,1220)
1220 format(/' * Warning - (EQ6/eqshel) Have done the maximum',' number of tries',/7x,'to slide over a region of critical',' redox instability.')

            ier = 190
            go to 999
        end if

        ! Slide Xi forward one increment. The increment gets progressively
        ! larger with the number of tries.
        nrdxsl = nrdxsl + 1
        delxi = delxi + 3.**(nrdxsl - 1)*dlxmin
        write (noutpt,1230)
1230 format(/' Trying to slide over a region of critical redox',' instability:')

        write (noutpt,1240) nrdxsl,delxi
1240 format(/3x,'Try ',i2,': delxi= ',1pe12.5,/)

        go to 300
    end if

    if (ier.eq.60 .or. ier.eq.70) then
        ! The calculation appears to be in a region of critical
        ! instability associated with solvent water. Go back and let
        ! EQ6/path.f deal with the situation. It may be able to resolve
        ! the problem by reducing the step size.
        ier = 180
        go to 999
    end if

    if (ier .eq. 80) then
        ! EQ6/eqphas.f detected critical instability associated with the
        ! aqueous activity coefficient model. Go back and let EQ6/path.f
        ! deal with the situation. It may be able to at least partially
        ! resolve the problem by reducing the step size. However, it is
        ! likely that the reaction path is causing the activity
        ! coefficient model to be extrapolated outside its range of
        ! validity.
        ier = 180
        go to 999
    end if

    if (ier .eq. 100) then
        ! The calculation was about to start with out-of-range values
        ! for the iteration variables associated with the aqueous basis
        ! species. Go back and let EQ6/path.f deal with the situation.
        ! It may be able to resolve the problem by reducing the step size.
        ! Alternatively, there may be some kind of programming error.
        ier = 180
        go to 999
    end if

    if (ier .eq. 110) then
        ! The iteration ended with out-of-range values for the iteration
        ! variables associated with the aqueous basis species. Go back and
        ! let EQ6/path.f deal with the situation. It may be able to
        ! resolve the problem by reducing the step size. Alternatively,
        ! there may be some kind of programming error.
        ier = 180
        go to 999
    end if

    if (ier .eq. 150) then
        ! The iteration ended without any useful diagnostics being
        ! produced. Go back and let EQ6/path.f deal with the situation.
        ! It may be able to resolve the problem by reducing the step size.
        ! Alternatively, there may be some kind of programming error.
        ier = 180
        go to 999
    end if

    if (ier .eq. 0) then
        if (iopt(1) .eq. 2) then
            ! Check to see that no significant numbers of moles of solids
            ! were destroyed unexpectedly if computing a fluid-centered
            ! flow-through open system model.
            do kcol = km1s,kmts
                if ( uzvec0(kcol)(1:5) .ne. ufixf(1:5) ) then
                    if (zvclg0(kcol) .le. zklogu) then
                        ns = iindxs(kcol)
                        lxx = losp(ns)
                        dlmoph = zvec0(kcol) - texp(lxx)

                        if (dlmoph .gt. 0.) then
                            if (tlg(dlmoph) .gt. zklogu) then
                                ier = 170
                                mxx0 = zvec0(kcol)

                                ! Calling sequence substitutions:
                                !   uspec(ns) for unam48
                                call fmspnm(jlen,uspec(ns),uspn56)
                                write (noutpt,1300) uspn56(1:jlen),mxx0,mosp(ns),dlmoph
1300 format(' Some of ',a,' was unexpectedly destroyed.',/5x,'Previous number of moles was ',1pe12.5,/5x,'Current number of moles is   ',e12.5,/5x,'Amount destroyed was ',e12.5,/21x,/3x,"That's too much. Will go back and first",' transfer some of the current',/3x,'amount to the',' physically removed system (PRS).')
                            end if
                        end if
                    end if
                end if
            end do

            nplast = 0

            do kcol = kx1s,kxts
                ns = iindxs(kcol)
                np = ipndxs(kcol)

                if (np .ne. nplast) then
                    nplast = np
                    mold = 0.

                    do krow = kcol,kxts
                        np1 = ipndxs(krow)

                        if (np1 .ne. np) then
                            go to 200
                        end if

                        mxx0 = zvec0(krow)
                        mold = mold + mxx0
                    end do

200 continue

                    lold = tlg(mold)

                    if (lold .gt. zklogu) then
                        dlmoph = mold - moph(np)

                        if (dlmoph .gt. 0.) then
                            if (tlg(dlmoph) .gt. zklogu) then
                                ier = 170
                                j2 = ilnobl(uphase(np))
                                write (noutpt,1300) uphase(np)(1:j2),mold,moph(np),dlmoph
                            end if
                        end if
                    end if
                end if
            end do

            go to 999
        end if
    end if

    ! The calculation finished with no problems.
    go to 999

300 continue

    ! Restore the kernel for the equilibrium sytem from the backup.
    ! Undo any basis switches that may have been made by EQ6/optmzr.f.
    ! Note: some aspects of the backup kernel are not needed here.
    ! For example, the zvclg1 array will be recomputed further below
    ! by the call to EQ6/ztaylr.f, overwriting the values restored from
    ! the zvclgs array. However, a full restoration is done here for the
    ! sake of completeness.
    if (qbseqc) then
        do kcol = 1,kbt
            nb = iindx1(kcol)
            ns = nbasp(nb)
            ns2 = nbasps(nb)

            if (ns2 .ne. ns) then
                call switch(adhfs,adhfsx,advfs,advfsx,axhfs,axhfsx,axlks,axlksx,axvfs,axvfsx,cdrs,cdrsx,eps100,ipch,ipchmx,ipcv,ipcvmx,jflag,jsflag,narn1,narxmx,nbasp,nbaspd,nbaspx,nb,nbt,nbtmax,nbw,ndrs,ndrsmx,ndrsx,ndrsr,ndrsrx,noutpt,ns2,nst,nstmax,ntprmx,nttyo,qbassw,qbswok,uspec)
            end if
        end do

        qbseqc = .false.
        nord = nords
    end if

    km1 = km1s
    kmt = kmts
    kx1 = kx1s
    kxt = kxts
    kdim = kdims

    call copyia(nbasps,nbasp,nbt)
    call copyia(iindxs,iindx1,kdim)
    call copyia(ipndxs,ipndx1,kdim)

    call copyaa(zvclgs,zvclg1,kdim)
    call copyaa(acflgs,acflg,nst)

    xbarwc = xbarws
    xbrwlc = xbrwls

    ! Reset the ixbasp and cjbasp arrays. The former is a flag
    ! array, each member of which denotes whether the
    ! thermodynamic activity of the corresponding basis species
    ! is defined in terms of molality (= 0) or mole fraction (= 1).
    ! The cjbasp array contains any site stoichiometric factors
    ! associated with the operational basis species.
    call gibasp(cgexj,cjbasp,iern1,ixbasp,jern1,jern2,jetmax,jgext,narn1,narn2,nbasp,nbt,nbtmax,nern1,nern2,netmax,nphasx,nstmax)

    ! Slightly increase the value of reaction progress to try to slide
    ! over a region of critical instability.
    xi1 = xi0 + delxi

    ! Make a Taylor's series expansion of the z vector, applying
    ! change limits. This provides a protected set of values for
    ! the AE solver; i.e., a set not containing any values that
    ! will cause the solver to fail in an unrecoverable fashion.
    qztayl = .true.
    call ztaylr(delxi,dzvc0,kdim,kmax,km1,kxt,nord,nrd1mx,qztayl,zklogu,zvclg0,zvclg1,zvec0,zvec1)

    ! Save the new z vector expansion.
    call copyaa(zvclg1,zvclgs,kdim)

    if (.not.qcntmp .or. .not.qcnpre) then
        ! Recompute the temperature and pressure. Then recompute the
        ! thermodynamic and kinetic quantities which depend these
        ! variables.
        ntpr0 = ntpr
        call tpadv(aadh,aadhh,aadhv,aaphi,abdh,abdhh,abdhv,abdoth,abdot,abdotv,adadhh,adadhv,adbdhh,adbdhv,adbdth,adbdtv,adh,adhfe,adhh,adhv,adhfs,adhfsd,advfe,advfs,advfsd,afcnst,al10,amu,aslm,aphi,aprehw,apresg,apresh,apx,avcnst,axhfe,axhfs,axhfsd,axlke,axlks,axlksd,axvfe,axvfs,axvfsd,bdh,bdhh,bdhv,bdot,bdoth,bdotv,dadhh,dadhv,dbdhh,dbdhv,dbdth,dbdtv,dhfe,dhfs,dhfsd,dvfe,dvfs,dvfsd,eact,ehfac,farad,hact,iact,iapxmx,iktmax,imchmx,imech,iopg,iopt,ipbtmx,ipch,ipchmx,ipcv,ipcvmx,ixrn1,ixrn2,jpfcmx,jpress,jptffl,jsol,jtemp,narxmx,narxt,narxth,nbasp,nbaspd,nbt,nbtd,nbtmax,ncmpr,ndrsr,ndrsrd,nmut,nmutmx,nopgmx,noptmx,noutpt,nptkmx,nptmax,nrct,nrctmx,nrk,nslt,nsltmx,nst,nstmax,ntpr,ntprmx,ntprt,nttkmx,nttyo,nweope,nwndpc,nxt,nxtmax,pmu,presg,presh,press,pressb,pressd,pslamn,ptk,rcnstv,rconst,rk,rkb,rtcnst,tempc,tempcb,tempcd,tempcu,tempk,time1,trkb,ttk,uphase,uspec,wfac,xhfe,xhfs,xhfsd,xi1,xlke,xlks,xlksd,xvfe,xvfs,xvfsd)

        qtrch = ntpr .ne. ntpr0
    end if

    go to 100

999 continue
end subroutine eqshel