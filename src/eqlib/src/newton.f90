subroutine newton(aamatr,abar,acflg,acflgo,act,actlg,actwlc,adh,adhh,adhv,afcnst,alpha,al10,amtb,aphi,azero,a3bar,a3bars,bacfmx,bbig,beta,betamx,betao,bdh,bdhh,bdhv,bdot,bdoth,bdotv,bfje,bfxi,bgamx,bneg,bpx,bsigmm,cco2,cdrs,cdrtw,cdrw,cegexs,cgexj,cjbasp,cnufac,conc,conclg,coval,cpgexs,csts,delam,delmax,delvco,delvec,dgpit,dlogxw,dpelm,dpslm,dselm,egexjc,egexjf,egexs,eh,ehfac,elam,eps100,fje,fjeo,fo2,fo2lg,fsort,fugac,fugalg,fxi,fxio,gmmatr,gpit,ibpxt,idelmx,iebal,ielam,ier,iern1,iern2,ifcphi1,ifcphi2,ifnnn,ifn2n,ifpsi1,ifpsi2,ifrn1,ifrn2,ifzeta,igas,igstak,iindx1,ilcphi1,ilcphi2,ilnnn,iln2n,ilpsi1,ilpsi2,ilrn1,ilrn2,ilzeta,imrn1,imrn2,insgf,iodb,iopg,ipivot,ipndx1,irdxc3,istack,iter,itermx,ixbasp,ixrn1,ixrn2,izmax,jcsort,jern1,jern2,jflag,jgext,jgsort,jgstak,jjsort,jpflag,jsflag,jsitex,jsol,jssort,jstack,kbt,kction,kdim,kelect,khydr,km1,kmt,ko2gaq,kwater,kx1,kxt,loph,losp,lsort,mgext,moph,mosp,mrgexs,mtb,nalpha,napt,narn1,narn2,nbasp,nbt,nbw,nchlor,ncmpr,ncosp,ndrs,ndrsr,nelect,nern1,nern2,net,ngexsa,ngext,ngrn1,ngrn2,ngt,nhydr,nmut,nmux,nmxi,nmxx,noutpt,no2gaq,nphasx,npt,nredox,nslt,nslx,nst,nsts,nstsr,nsxi,nsxx,ntfx,ntfxt,nttyo,omega,omeglg,palpha,pelm,pmu,press,pslamn,pslm,qhawep,qpit75,qredox,q6mode,rhsvec,screwd,screwn,selm,sigmam,sigmmo,tempk,tfx,tolbt,toldl,ubacmx,ubbig,ubetmx,ubgamx,ubneg,ugexj,ugexmo,ulbeta,uldel,uphase,uspec,uzvec1,weight,wfac,xbar,xbarlg,xbarw,xbarwc,xbrwlc,xbrwlg,xlke,xlks,zchar,zchsq2,zchcu6,zgexj,zvclg1,zvec1)
    !! This subroutine performs hybrid Newton-Raphson iteration to solve
    !! for the equilibrium state of a chemical system.
    !! This subroutine is called by:
    !!   EQ3NR/eq3nr.f
    !!   EQ6/eqcalc.f
    !! Principal input:
    !!   ielam  = flag to not use (-1) or use (0) "higher order"
    !!              (higher than 2nd order) electrostatic terms in
    !!              those activity coefficient models which contain
    !!              provision for such terms
    !!   screwd = under-relaxation control parameter. It is used to
    !!              reduce the magnitude of the del vector, if necessary
    !!              so that the magnitude of the largest element of that
    !!              vector does not exceed screwd.
    !!   screwn = under-relaxation control parameter
    !!   tolbt  = convergence tolerance on betamx
    !!   toldl  = convergence tolerance on delmax
    !!   ulbeta = label for the quantity type for beta vector;
    !!              e.g., 'conc' or 'moles'
    !!   uldel  = label for the quantity type for del corrections;
    !!             e.g., 'conc' or 'moles'
    !!   uspec  = aqueous species name array
    !!   uzvec1 = name array corresponding to zvclg1
    !!   iodb   = array of debugging print options
    !!   itermx = maximum number of iterations
    !!   kdim   = dimension of aamatr
    !!   narn1  = start of species range for aqueous solution
    !!   narn2  = end of species range for aqueous solution
    !!   nbasp  = array of basis species
    !!   nbt    = number of basis species
    !!   ncmpr  = species range array for phases
    !!   nern1  = start of species range for generic ion exchangers
    !!   nern2  = end of species range for generic ion exchangers
    !!   nst    = number of aqeuous species
    !!   q6mode = flag denoting usage for EQ3NR or EQ6:
    !!              .false. = EQ3NR
    !!              .true.  = EQ6NR
    !!   wfac   = array of solid solution non-ideality parameters
    !! Principal input/output:
    !!   zvclg1 = the 'log z' array, the array corrected by
    !!              Newton-Raphson iteration
    !!   conc   = concentration array
    !!   acflg  = activity coefficient array
    !!   fje    = the ionic asymmetry (the 3rd-order electrostatic
    !!              moment function J)
    !!   fxi    = the ionic strength (the 2nd-order electrostatic
    !!              moment function I)
    !! Principal work space/output:
    !!   aamatr = Jacobian matrix
    !!   abar   = average ion size
    !!   a3bar  = average cube of the distance of closest approach
    !!              for all solute species
    !!   a3bars = characteristic average cube of the distance of closest
    !!              approach for a solute species
    !!   gmmatr = copy of aamatr
    !!   delvec = correction array
    !!   rhsvec = right hand side vector
    !!   beta   = normalized residual function array
    !!   alpha  = residual function array
    !!   acflgo = copy of activity coefficient array
    !!   betao  = old beta array
    !!   delvco = old del array
    !!   betamx = norm (largest magnitude) of the beta array
    !!   bbig   = largest magnitude positive mass balance residual
    !!   bneg   = largest magnitude negative mass balance residual
    !!   bacfmx = norm (largest magnitude) activity coefficient residual
    !!   ubbig  = name of species corresponding to bbig
    !!   ubneg  = name of species corresponding to bneg
    !!   ubacmx = name of species corresponding to bacfmx
    !!   ipivot = the pivot vector
    !!   iter   = Newton-Raphson iteration number
    !!   idelmx = kcol index corresponding to delmax
    !!   ier    = error flag (mostly returned from EQLIB/nrstep.f):
    !!              =  0  Okay
    !!              =  1  Encountered a zero matrix
    !!              =  2  Encountered a non-zero, computationally
    !!                      singular matrix
    !!              =  3  Iteration was diverging
    !!              =  4  Hit the maximum number of iterations (itermx)
    implicit none

    include 'eqlib/eqldv.h'

    ! Calling sequence variable declarations.
    integer :: noutpt
    integer :: nttyo

    integer :: ibpxt(nxtmax)
    integer :: igstak(ngtmax)
    integer :: iindx1(kmax)
    integer :: insgf(natmax)
    integer :: iodb(nodbmx)
    integer :: iopg(nopgmx)
    integer :: ipivot(kmax)
    integer :: ipndx1(kmax)
    integer :: istack(nstmax)
    integer :: ixbasp(nbtmax)
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
    integer :: jsol(nxtmax)
    integer :: jssort(nstmax)
    integer :: jstack(nstmax)
    integer :: kction(nbtmax)
    integer :: nbasp(nbtmax)
    integer :: ncmpr(2,nptmax)
    integer :: ncosp(nbtmax)
    integer :: ndrs(ndrsmx)
    integer :: ndrsr(2,nstmax)
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

    integer :: idelmx
    integer :: iebal
    integer :: ielam
    integer :: ier
    integer :: iern1
    integer :: iern2
    integer :: ifrn1
    integer :: ifrn2
    integer :: igas
    integer :: ilrn1
    integer :: ilrn2
    integer :: imrn1
    integer :: imrn2
    integer :: iter
    integer :: itermx
    integer :: irdxc3
    integer :: ixrn1
    integer :: ixrn2
    integer :: izmax
    integer :: kbt
    integer :: kdim
    integer :: kelect
    integer :: khydr
    integer :: km1
    integer :: kmt
    integer :: ko2gaq
    integer :: kwater
    integer :: kx1
    integer :: kxt
    integer :: narn1
    integer :: narn2
    integer :: nbt
    integer :: nbw
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
    integer :: nredox
    integer :: nst
    integer :: ntfxt

    logical :: qhawep
    logical :: qpit75
    logical :: qredox
    logical :: q6mode

    character(len=48) :: uspec(nstmax)
    character(len=48) :: uzvec1(kmax)
    character(len=48) :: ubbig
    character(len=48) :: ubetmx
    character(len=48) :: ubgamx
    character(len=48) :: ubneg
    character(len=48) :: ubacmx
    character(len=24) :: ugexmo(netmax)
    character(len=24) :: uphase(nptmax)
    character(len=8) :: ugexj(jetmax,netmax)
    character(len=8) :: ulbeta(kmax)
    character(len=8) :: uldel(kmax)

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
    real(kind=8) :: betao(kmax)
    real(kind=8) :: bpx(ibpxmx,nxtmax)
    real(kind=8) :: cco2(5)
    real(kind=8) :: cdrs(ndrsmx)
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
    real(kind=8) :: delvco(kmax)
    real(kind=8) :: delvec(kmax)
    real(kind=8) :: dlogxw(nbtmax)
    real(kind=8) :: egexjc(jetmax,netmax)
    real(kind=8) :: egexs(ietmax,jetmax,netmax)
    real(kind=8) :: egexjf(jetmax,netmax)

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
    real(kind=8) :: xlks(nstmax)
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
    real(kind=8) :: a3bar
    real(kind=8) :: bacfmx
    real(kind=8) :: bbig
    real(kind=8) :: betamx
    real(kind=8) :: bfje
    real(kind=8) :: bfxi
    real(kind=8) :: bgamx
    real(kind=8) :: bneg
    real(kind=8) :: bsigmm
    real(kind=8) :: delmax
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
    real(kind=8) :: press
    real(kind=8) :: screwd
    real(kind=8) :: screwn
    real(kind=8) :: sigmam
    real(kind=8) :: sigmmo
    real(kind=8) :: tempk
    real(kind=8) :: tolbt
    real(kind=8) :: toldl
    real(kind=8) :: xbarw
    real(kind=8) :: xbarwc
    real(kind=8) :: xbrwlc
    real(kind=8) :: xbrwlg
    real(kind=8) :: xlke

    ! Local variable declarations.
    integer :: ibetmx
    integer :: jlen
    integer :: j2
    integer :: kcol
    integer :: negbfc
    integer :: negdfc
    integer :: negxfc
    integer :: noibfc
    integer :: npconv
    integer :: npobfc
    integer :: npodfc
    integer :: ns

    integer :: ilnobl

    logical :: qcacf
    logical :: qcbeta
    logical :: qcdel
    logical :: qcgam
    logical :: qconv
    logical :: qcsigm
    logical :: qcfxi
    logical :: qpracf
    logical :: qxbarw

    character(len=56) :: uspn56
    character(len=8) :: udelmx
    character(len=8) :: utb
    character(len=8) :: utd

    real(kind=8) :: betfnc
    real(kind=8) :: betmxo
    real(kind=8) :: btfcnr
    real(kind=8) :: bx
    real(kind=8) :: dx
    real(kind=8) :: chfacf
    real(kind=8) :: chfsgm
    real(kind=8) :: delfnc
    real(kind=8) :: rdx
    real(kind=8) :: rlxfac
    real(kind=8) :: rlxgam

    iter = 0
    ier = 0
    ibetmx = 0
    idelmx = 0
    npconv = 0

    delmax = 0.
    betmxo = 0.
    betfnc = 0.
    delfnc = 0.
    btfcnr = 0.
    rlxfac = 1.
    rlxgam = 1.
    negdfc = 0
    negbfc = 0
    negxfc = 0
    noibfc = 0
    npodfc = 0
    npobfc = 0

    do kcol = 1,kdim
        delvec(kcol) = 0.
    end do

    qpracf = iodb(4) .ge. 4

    if (iodb(4) .ge. 2) then
        ! Debugging prints- state of the system prior to hybrid
        ! Newton-Raphson iteration.
        write (noutpt,1000)
1000 format(/8x,'Name',36x,'zvclg1',/)

        do kcol = 1,kdim
            ! Calling sequence substitutions:
            !   uzvec1(kcol) for unam48
            call fmspnx(jlen,uzvec1(kcol),uspn56)
            jlen = min(jlen,38)
            write (noutpt,1010) kcol,uspn56(1:jlen),zvclg1(kcol)
1010 format(1x,i3,2x,a,t46,1pe12.5)
        end do

        write (noutpt,1015)
1015 format(1x)
    end if

    if (qpracf) then
        write (noutpt,1020) sigmam,fxi,fje,xbrwlc,xbarwc
1020 format(/3x,'sigmam= ',1pe12.5,/6x,'fxi= ',e12.5,/6x,'fje= ',e12.5,//6x,'xbrwlc= ',0pf10.5,/6x,'xbarwc= ',1pe12.5)

        write (noutpt,1030)
1030 format(/9x,'Species',21x,'gamma',/)

        do ns = narn1,narn2
            write (noutpt,1040) ns,uspec(ns),acflg(ns)
1040 format(1x,i4,2x,a24,2x,1pe12.5)
        end do
    end if

    ! Compute the Newton-Raphson residual functions.
    call betas(acflg,actlg,afcnst,alpha,amtb,bbig,beta,betamx,bneg,cdrs,conc,conclg,coval,csts,eh,ehfac,fo2lg,ibetmx,iebal,iindx1,irdxc3,jcsort,jflag,jsflag,jssort,kbt,kdim,kelect,khydr,kmax,km1,ko2gaq,kwater,kxt,mtb,mosp,narn1,narn2,nbasp,nbtmax,ncosp,ndrs,ndrsmx,ndrsr,nelect,nern1,nern2,nhydr,noutpt,no2gaq,nredox,nst,nstmax,nsts,nstsmx,nstsr,ntfx,ntfxmx,ntfxt,nttyo,omega,qredox,q6mode,tfx,ubbig,ubneg,ubetmx,uspec,uzvec1,weight,xbrwlg,xlke,xlks,zchar)

    ! The label below is the return point for iter >= 1.
100 continue
    sigmmo = sigmam
    fxio = fxi
    fjeo = fje

    ! Save the activity coefficients.
    call copyaa(acflg,acflgo,nst)

    bx = 0.
    dx = 0.

    if (ibetmx .gt. 0) then
        bx = beta(ibetmx)
    end if

    if (idelmx .gt. 0) then
        dx = delvec(idelmx)
    end if

    ! Calculate the beta convergence function.
    betfnc = 0.

    if (betmxo .gt. 0.) then
        betfnc = (betmxo - betamx)/betmxo
    end if

    betmxo = betamx

    if (iodb(4) .ge. 2) then
        write (noutpt,1100)
1100 format(/8x,'Name',37x,'Beta',/)

        do kcol = 1,kdim
            ! Calling sequence substitutions:
            !   uzvec1(kcol) for unam48
            call fmspnx(jlen,uzvec1(kcol),uspn56)
            jlen = min(jlen,38)
            write (noutpt,1110) kcol,uspn56(1:jlen),beta(kcol)
1110 format(1x,i3,2x,a,t46,1pe12.5)
        end do

        write (noutpt,1120)
1120 format(/1x)
    end if

    if (iodb(4) .ge. 1) then
        ubetmx = ' '
        utb = ' '
        udelmx = ' '
        utd = ' '

        if (ibetmx .gt. 0) then
            utb = ulbeta(ibetmx)
            ubetmx = uzvec1(ibetmx)
        end if

        if (idelmx .gt. 0) then
            utd = uldel(idelmx)
            udelmx = uzvec1(idelmx)(1:8)
        end if

        write (noutpt,1200) iter
1200 format(' iter= ',i3)

        if (abs(rlxgam) .le. eps100) then
            write (noutpt,1210)
1210 format(7x,'Gammas are fixed')
        else if (abs(rlxgam - 1.0) .gt. eps100) then
            write (noutpt,1220) rlxgam
1220 format(7x,'Gamma relaxation factor= ',g12.5)
        end if

        write (noutpt,1230) utd,udelmx,dx,delfnc
1230 format(5x,'delvec(',2a8,')= ',1pe12.5,', delfnc= ',e12.5)

        if (abs(rlxfac - 1.0) .gt. eps100) then
            rdx = rlxfac*dx
            write (noutpt,1240) rlxfac,utd,udelmx,rdx
1240 format(7x,'Relaxation factor= ',g12.5,/5x,'Relaxed delvec(',2a8,')= ',1pe12.5)
        end if

        write (noutpt,1250) utb,ubetmx,bx,betfnc
1250 format(7x,'beta(',2a8,')= ',1pe12.5,', betfnc= ',e12.5)

        ! Calling sequence substitutions:
        !   ubbig for unam48
        call fmspnx(jlen,ubbig,uspn56)
        write (noutpt,1260) bbig,uspn56(1:jlen)
1260 format(9x,'bbig= ',1pe12.5,',   ubbig= ',a)

        ! Calling sequence substitutions:
        !   ubneg for unam48
        call fmspnx(jlen,ubneg,uspn56)
        write (noutpt,1270) bneg,uspn56(1:jlen)
1270 format(9x,'bneg= ',1pe12.5,',   ubneg= ',a)

        j2 = ilnobl(ubgamx(1:24))
        write (noutpt,1280) bgamx,ubgamx(1:j2)
1280 format(8x,'bgamx= ',1pe12.5,',  ubgamx= ',a)

        write (noutpt,1290) bsigmm,bfxi,bfje
1290 format(7x,'bsigmm= ',1pe12.5,/9x,'bfxi= ',e12.5,/9x,'bfje= ',1pe12.5)

        write (noutpt,1300) btfcnr
1300 format(7x,'btfcnr= ',1pe12.5,/)
    end if

    ! Check to see if the iteration satisfies the convergence criteria.
    ! Both residual functions and correction terms are tested.
    ! Iteration may terminate acceptably without satisfying the
    ! constraint on the correction terms (see below).
    qcbeta = betamx .le. tolbt
    qcdel = delmax .le. toldl
    qconv = qcbeta .and. qcdel
    qcsigm = abs(bsigmm) .le. tolbt
    qcfxi = abs(bfxi) .le. tolbt
    qcgam = bgamx .le. tolbt
    qcacf = qcsigm .and. qcfxi .and. qcgam

    ! XX   Redefine qcacf when the activity coefficients of species in
    ! XX   non-aqueous phases are updated numerically the same as those
    ! XX   of aqueous species.
    !      qcacf = bacfmx .le. tolbt
    !      Force to run at least one iteration
    if (iter .ge. 1) then
        ! Test for convergence
        if (qconv .and. qcacf) then
            go to 999
        end if
    end if

    qxbarw = .false.

    ! XX   qxbarw = abs(bsigmm).le.0.05 .and. bgamx.le.0.05 .and.
    ! XX  $ abs(rlxgam - 1.0).le.eps100
    !      Do a Newton-Raphson step.
    call nrstep(aamatr,acflg,act,actlg,afcnst,alpha,al10,amtb,bbig,beta,betamx,betao,betfnc,betmxo,bneg,bpx,btfcnr,cdrs,cdrtw,cdrw,cegexs,cgexj,cjbasp,cnufac,conc,conclg,cpgexs,csts,coval,delfnc,delmax,delvco,delvec,dlogxw,egexjc,egexjf,egexs,eh,ehfac,eps100,fo2,fo2lg,fsort,fugac,fugalg,gmmatr,ibetmx,ibpxmx,idelmx,iebal,ier,iern1,iern2,ietmax,ifrn1,ifrn2,igas,igstak,iindx1,ilrn1,ilrn2,imrn1,imrn2,iodb,ipivot,ipndx1,irdxc3,istack,iter,itermx,ixbasp,ixrn1,ixrn2,jcsort,jern1,jern2,jetmax,jflag,jgext,jgsort,jgstak,jjsort,jpflag,jsflag,jsitex,jssort,jstack,kbt,kction,kdim,kelect,khydr,kmax,kmt,km1,ko2gaq,kwater,kxt,kx1,loph,losp,lsort,mgext,moph,mosp,mrgexs,mtb,narn1,narn2,nbasp,nbt,nbtmax,nbw,ncmpr,ncosp,ndrs,ndrsmx,ndrsr,negbfc,negdfc,negxfc,nelect,nern1,nern2,netmax,ngexsa,ngext,ngrn1,ngrn2,ngt,ngtmax,nhydr,nodbmx,noibfc,noutpt,no2gaq,npconv,nphasx,npobfc,npodfc,npt,nptmax,nredox,nst,nstmax,nsts,nstsmx,nstsr,ntfx,ntfxmx,ntfxt,nttyo,nxtmax,omega,omeglg,press,qcacf,qcbeta,qredox,qxbarw,q6mode,rhsvec,rlxfac,screwd,screwn,sigmam,sigmmo,tfx,ubetmx,ubbig,ubneg,ugexj,ugexmo,uphase,uspec,uzvec1,weight,xbar,xbarlg,xbarw,xbarwc,xbrwlc,xbrwlg,xlke,xlks,zchar,zgexj,zvclg1,zvec1)

    if (ier .gt. 0) then
        ! Iteration has been stopped.
        if (qcbeta .and. qcgam) then
            ! Have pseudo-convergence. Reset ier to 0.
            ier = 0

            if (iodb(4) .ge. 1) then
                write (noutpt,1310)
1310 format('  Hybrid Newton-Raphson iteration has',' pseudo-converged.')

                do kcol = 1,kdim
                    if (abs(delvec(kcol)) .gt. toldl) then
                        ! Calling sequence substitutions:
                        !   uzvec1(kcol) for unam48
                        call fmspnx(jlen,uzvec1(kcol),uspn56)
                        write (noutpt,1320) uspn56(1:jlen),delvec(kcol)
1320 format(11x,'delvec(',a,') = ',g12.5)
                    end if
                end do
            end if

            go to 999
        else
            ! Iteration has failed.
            write (noutpt,1330) iter
1330 format('  Hybrid Newton-Raphson iteration has gone sour',' (iter= ',i3,').')

            go to 999
        end if
    end if

    ! Recompute sigmam, fxi, and fje. Recompute the activity
    ! coefficients. Then recompute the concentrations of dependent
    ! species. After that, recalculate the Newton-Raphson residual
    ! functions. Finally, go back and see whether or not to do another
    ! iteration.
    if (iodb(4) .ge. 2) then
        write (noutpt,1340)
    end if

1340 format(/3x,'--- Post-Newton-Raphson update of activity',' coefficients ---',/)

    ! Determine the maximum allowed change factors for "Sigma m",
    ! etc. (chfsgm) and the log activity coefficients for aqueous
    ! species (chfacf).
    chfsgm = 1.5

    if (sigmam .le. 2.e-1) then
        chfsgm = 5.
    end if

    if (sigmam .le. 1.e-2) then
        chfsgm = 10.
    end if

    chfacf = 0.50

    call ngcadv(abar,acflg,acflgo,actwlc,adh,adhh,adhv,afcnst,al10,aphi,azero,a3bar,a3bars,bacfmx,bdh,bdhh,bdhv,bdot,bdoth,bdotv,bgamx,bpx,bsigmm,bfje,bfxi,cco2,cgexj,chfacf,chfsgm,conc,delam,dgpit,dpelm,dpslm,dselm,elam,eps100,fje,fjeo,fxi,fxio,gpit,ibpxt,ielam,ifcphi1,ifcphi2,ifnnn,ifn2n,ifpsi1,ifpsi2,ifzeta,ilcphi1,ilcphi2,ilnnn,iln2n,ilpsi1,ilpsi2,ilzeta,insgf,iopg,iter,ipndx1,ixrn1,ixrn2,izmax,jcsort,jern1,jern2,jgext,jsol,kx1,kxt,nalpha,napt,narn1,narn2,nchlor,ncmpr,net,nhydr,nmut,nmux,nmxi,nmxx,noutpt,nslt,nslx,nst,nsxi,nsxx,nttyo,omega,palpha,pelm,pmu,press,pslamn,pslm,qhawep,qpit75,qpracf,q6mode,rlxgam,selm,sigmam,sigmmo,tempk,ubacmx,ubgamx,uphase,uspec,wfac,xbar,xbarlg,xbarwc,xbrwlc,zchar,zchcu6,zchsq2)

    ! Recalculate the concentrations, etc., of dependent species.
    call ncmpex(acflg,act,actlg,cdrs,cegexs,cgexj,conc,conclg,cpgexs,egexjc,egexjf,egexs,eps100,fo2,fo2lg,fsort,fugac,fugalg,iern1,iern2,ietmax,ifrn1,ifrn2,igas,igstak,iindx1,ilrn1,ilrn2,imrn1,imrn2,istack,ixrn1,ixrn2,jcsort,jern1,jern2,jetmax,jflag,jgext,jgsort,jgstak,jjsort,jpflag,jsflag,jsitex,jssort,jstack,kbt,kdim,kelect,kmax,km1,ko2gaq,kwater,kxt,loph,losp,lsort,mgext,mrgexs,mtb,moph,mosp,narn1,narn2,nbasp,nbt,nbtmax,ncmpr,ndrs,ndrsmx,ndrsr,nelect,nern1,nern2,netmax,ngexsa,ngext,ngrn1,ngrn2,ngt,ngtmax,noutpt,no2gaq,nphasx,npt,nptmax,nst,nstmax,nttyo,omega,omeglg,press,qxbarw,q6mode,ugexj,ugexmo,uphase,uspec,xbar,xbarlg,xbarw,xbarwc,xbrwlc,xbrwlg,xlks,zchar,zgexj,zvclg1,zvec1)

    ! Recalculate the Newton-Raphson residual functions.
    call betas(acflg,actlg,afcnst,alpha,amtb,bbig,beta,betamx,bneg,cdrs,conc,conclg,coval,csts,eh,ehfac,fo2lg,ibetmx,iebal,iindx1,irdxc3,jcsort,jflag,jsflag,jssort,kbt,kdim,kelect,khydr,kmax,km1,ko2gaq,kwater,kxt,mtb,mosp,narn1,narn2,nbasp,nbtmax,ncosp,ndrs,ndrsmx,ndrsr,nelect,nern1,nern2,nhydr,noutpt,no2gaq,nredox,nst,nstmax,nsts,nstsmx,nstsr,ntfx,ntfxmx,ntfxt,nttyo,omega,qredox,q6mode,tfx,ubbig,ubneg,ubetmx,uspec,uzvec1,weight,xbrwlg,xlke,xlks,zchar)

    go to 100

999 continue
end subroutine newton