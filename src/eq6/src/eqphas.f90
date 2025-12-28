subroutine eqphas(aamatr,abar,acflg,acflgo,act,actlg,adh,adhh,adhv,afcnst,affp,affs,alpha,al10,amtb,aphi,apx,avcnst,azero,a3bar,a3bars,bacfmx,bbig,bdh,bdhh,bdhv,bdot,bdoth,bdotv,beta,betamx,betao,bgamx,bneg,bpx,cco2,cegexs,cess,cdrs,cdrsd,cdrsx,cdrtw,cdrw,cjbasp,cnufac,conc,conclg,cpgexs,cscale,csts,delvco,delvec,d1zvc1,dlogxw,egexjc,egexjf,egexs,eh,ehfac,eps100,farad,fje,fjeo,fo2,fo2lg,fsort,fugac,fugalg,fxi,fxio,gmmatr,iapxt,ibpxt,ibswx,ielam,ier,iern1,iern2,ifcphi1,ifcphi2,ifnnn,ifn2n,ifpsi1,ifpsi2,ifrn1,ifrn2,ifzeta,igas,igstak,iindx0,iindx1,ilcphi1,ilcphi2,ilnnn,iln2n,ilpsi1,ilpsi2,ilrn1,ilrn2,ilzeta,imrn1,imrn2,insgf,iodb,iopg,iopt,ipch,ipivot,ipndx1,ipcv,istack,iter,itermx,ixbasp,ixrn1,ixrn2,izmax,jcsort,jflag,jgsort,jgstak,jjsort,jpflag,jsflag,jsitex,jsol,jssort,jstack,kbt,kction,kdim,kelect,khydr,khydx,km1,km10,kmt,kmt0,ko2gaq,kpsat,kpsst,krdxsp,kwater,kx1,kx10,kxt0,kxt,loph,losp,lsort,moph,mosp,mrgexs,mtb,mtbaq,narn1,narn2,narxt,nat,nbasp,nbaspd,nbaspx,nbt,nbtd,nbw,nchlor,ncmpr,nct,ndrs,ndrsd,ndrsx,ndrsr,ndrsrd,ndrsrx,nelect,nern1,nern2,ness,nessr,net,nfrn1,nfrn2,ngrn1,ngrn2,ngt,nhydr,nhydx,nlrn1,nlrn2,nlt,nmrn1,nmrn2,nmt,nord,no2gaq,noutpt,npchk,nphasx,npt,nrdxsp,nst,nsts,nstsr,ntpr,ntrymx,nttyo,nxrn1,nxrn2,nxt,omega,omeglg,prcinf,press,qbassw,qbseqc,qbye,qcnpre,qcntmp,qhawep,qmod,qoptmz,qpit75,qredox,qsspgb,qstart,qxknph,q6mode,rcnstv,rconst,rhsvec,rtcnst,screwd,sidrph,sidrsp,sigmam,sigmmo,smp100,tempc,tempk,tolbt,toldl,tolsat,tolsst,ubacmx,ubgamx,ulbeta,uldel,uphase,uspec,uzvec1,weight,wfac,xbar,xbarlg,xbarw,xbarwc,xbrwlc,xbrwlg,zchar,zchcu6,zchsq2,zvclg1,zvec1)
    !! This subroutine attempts to calculate the equilibrium state of
    !! the equilibrium system. As necessary, it makes repeated calls
    !! to EQ6/eqcalc.f, which attempts to make this calculation for a
    !! given specified phase assemblage. The present subroutine adjusts
    !! this assemblage as needed.
    !! This subroutine is called by:
    !!   EQ6/eqshel.f
    !! Principal input:
    !! Principal output:
    !!   qbye   = flag which is .true. if EQ6/path.f just changed the
    !!              phase assemblage by deleting one or more phases.
    !!              This condition instructs the present subroutine to
    !!              print the index structure of the system of equations
    !!              as it does on the starting call and whenever this
    !!              subroutine itself adds or deletes a phase.
    !!   iter   = the number of hydrbid Newton-Raphson iterations
    !!              done by EQLIB/newton.f
    !!   ier    = error flag:
    !!              Values returned from EQ6/eqcalc.f:
    !!                =    0  Okay
    !!                =    1  Encountered a zero matrix
    !!                =    2  Encountered a non-zero, computationally
    !!                          singular matrix
    !!                =    3  Iteration was diverging
    !!                =    4  Hit the maximum number of iterations
    !!                          (itermx)
    !!              Values returned from EQ6/miidxz.f:
    !!                =    0  Okay
    !!                =    8  Exceeded the dimension of the iindx1
    !!                          array
    !!              Values returned by the present subroutine:
    !!                =    0  Okay
    !!                =    8  Exceeded the dimension of the iindx1
    !!                          array; go back and cut the step size
    !!                          if possible
    !!                =   10  Go back and take a smaller step size to
    !!                          avoid exceeding the supersaturation
    !!                          tolerance (tolsst). This is done only
    !!                          when an appropriate set of conditions
    !!                          is satisified. It isn't necessary for
    !!                          EQ6/path.f to analyze the situation
    !!                          when this value is returned to it by
    !!                          EQ6/ eqshel.f.
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
    !!                          instability associated with the
    !!                          master redox variable. This is probably
    !!                          associated with a redox jump.
    !!                =   60  Caught in a region of computational
    !!                          instability associated with solvent
    !!                          water. The amount of water in the
    !!                          system is probably very low.
    !!                =   70  Caught in a region of computational
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

    ! Calling sequence variable declarations.
    integer :: noutpt
    integer :: nttyo

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
    integer :: ixbasp(nbtmax)
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
    integer :: nfac(nbtmax)
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
    integer :: nelect
    integer :: nhydr
    integer :: nhydx
    integer :: nord
    integer :: no2gaq
    integer :: nrdxsp
    integer :: ntpr
    integer :: ntrymx

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
    logical :: qsspgb
    logical :: qstart
    logical :: q6mode

    character(len=48) :: uspec(nstmax)
    character(len=48) :: uzvec1(kmax)
    character(len=48) :: ubacmx
    character(len=48) :: ubgamx
    character(len=24) :: uphase(nptmax)
    character(len=8) :: ulbeta(kmax)
    character(len=8) :: uldel(kmax)

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
    real(kind=8) :: delvco(kmax)
    real(kind=8) :: delvec(kmax)
    real(kind=8) :: dlogxw(nbtmax)
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
    real(kind=8) :: moph(nptmax)
    real(kind=8) :: mosp(nstmax)
    real(kind=8) :: mrgexs(ietmax,jetmax,netmax)
    real(kind=8) :: mtb(nbtmax)
    real(kind=8) :: mtbaq(nbtmax)
    real(kind=8) :: rhsvec(kmax)
    real(kind=8) :: sidrph(nptmax)
    real(kind=8) :: sidrsp(nstmax)
    real(kind=8) :: weight(nstmax)
    real(kind=8) :: wfac(iktmax,nxtmax)
    real(kind=8) :: xbar(nstmax)
    real(kind=8) :: xbarlg(nstmax)
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
    real(kind=8) :: rcnstv
    real(kind=8) :: rconst
    real(kind=8) :: rtcnst
    real(kind=8) :: screwd
    real(kind=8) :: sigmam
    real(kind=8) :: sigmmo
    real(kind=8) :: smp100
    real(kind=8) :: tempc
    real(kind=8) :: tempk
    real(kind=8) :: tolbt
    real(kind=8) :: toldl
    real(kind=8) :: tolsat
    real(kind=8) :: tolsst
    real(kind=8) :: xbarw
    real(kind=8) :: xbarwc
    real(kind=8) :: xbrwlc
    real(kind=8) :: xbrwlg

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

    ! Local variable declarations with special dimensioning.
    ! Data for the eight most supersaturated minerals.
    integer :: nssppa
    parameter (nssppa = 8)

    integer :: nsspmx

    integer :: nssp(nssppa)

    character(len=24) :: ussp(nssppa)

    real(kind=8) :: afssp(nssppa)
    real(kind=8) :: afssps(nssppa)
    real(kind=8) :: msspmx(nssppa)

    ! Local variable declarations.
    integer :: ibetmx
    integer :: idelmx
    integer :: ifail
    integer :: inext
    integer :: irow1
    integer :: irow2
    integer :: isave
    integer :: isspt
    integer :: j
    integer :: jcol1
    integer :: jcol2
    integer :: jj
    integer :: jkl
    integer :: jlen
    integer :: j2
    integer :: j3
    integer :: k
    integer :: kcol
    integer :: n
    integer :: nb
    integer :: nerr
    integer :: nords
    integer :: np
    integer :: npadd
    integer :: npaddi
    integer :: npdel
    integer :: npdeli
    integer :: nrn1
    integer :: nrn2
    integer :: nr1
    integer :: nr2
    integer :: ns
    integer :: nss
    integer :: nsspt
    integer :: nswtch
    integer :: ns2
    integer :: nt
    integer :: ntry
    integer :: ntryai
    integer :: ntrydi
    integer :: nxx

    integer :: ilnobl

    logical :: qadd
    logical :: qbswok
    logical :: qphruv

    character(len=56) :: uspn56
    character(len=8) :: ux8

    real(kind=8) :: afscal
    real(kind=8) :: bfje
    real(kind=8) :: bfxi
    real(kind=8) :: bsigmm
    real(kind=8) :: cx
    real(kind=8) :: cxx
    real(kind=8) :: delmax
    real(kind=8) :: mophmx
    real(kind=8) :: mophn
    real(kind=8) :: mx
    real(kind=8) :: mxmadd
    real(kind=8) :: mpmxaq
    real(kind=8) :: mpmxes
    real(kind=8) :: msmxaq
    real(kind=8) :: msmxes

    real(kind=8) :: tlg

    ! Local dimensioning variables.
    nsspmx = nssppa

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
    ier = 0

    qadd = .false.
    qmod = .false.
    qbseqc = .false.

    ifail = 0

    npadd = 0
    npaddi = 0

    npdel = 0
    npdeli = 0

    ntry = 1
    ntryai = 0
    ntrydi = 0

    kpsst = 0
    nsspt = 0
    call initiz(nssp,nsspmx)

    ! Here is a return point if the phase assemblage has been modified.
100 continue
    bgamx = 0.
    bsigmm = 0.
    bfxi = 0.
    bfje = 0.

    if (ntry.gt.1 .or. qstart .or. qbye .or. iodb(1).ge.3) then
        ! Print the current phase assemblage.
        write (noutpt,1000) ntry
1000 format(/' Attempted phase assemblage number ',i3,/)

        do kcol = 1,kbt
            nb = iindx1(kcol)
            ns = nbaspd(nb)

            ! Calling sequence substitutions:
            !   uspec(ns) for unam48
            call fmspnx(jlen,uspec(ns),uspn56)
            write (noutpt,1010) kcol,uspn56(1:jlen)
1010 format(2x,i3,2x,a)
        end do

        do kcol = km1,kxt
            ns = iindx1(kcol)

            ! Calling sequence substitutions:
            !   uspec(ns) for unam48
            call fmspnm(jlen,uspec(ns),uspn56)
            write (noutpt,1010) kcol,uspn56(1:jlen)
        end do

        write (noutpt,1020)
1020 format(1x)
    end if

    ! Save the "kernel" description for the current equilibrium system.
    ! This will be used to recover if an equilibrium calculation made
    ! by EQ6/eqcalc.f fails for any reason. This kernel contains the
    ! minimum of information required to calculate the fully expanded
    ! description using EQLIB/ncmpex.f. Note that EQ6/eqshel.f, which
    ! calls the present subroutine, has its own backup kernel.
    ! EQ6/path.f, which calls EQ6/eqshel.f, also has its own backup
    ! kernel, which is used to restore the equilibrium system at the
    ! previous point of reaction progress.
    ! Note that the present backup kernel is not a compact description
    ! of the equilibrium system from the last call to EQ6/eqcalc.f.
    ! Upon entry to the present subroutine, the backup kernel
    ! corresponds to a description of the system either as read from
    ! the input file (at the starting value of reaction progress) or
    ! as obtained by finite-difference prediction in stepping to the
    ! current point of reaction progress. Otherwise, it corresponds to
    ! what resulted from the last successful call to EQ6/eqcalc.f from
    ! the present subroutine, modified by the addition, deletion, or
    ! replacement of one phase.
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

    ! Check for aqueous basis species (or any other species switched
    ! with such) for cases in which the log number of moles variable
    ! is -99999. This implies a number of moles values which is zero,
    ! which in turn implies that the species is not actually present in
    ! the equilibrium system. Note that in this code -99999. is the
    ! conventional value for log10(0). Conversely, 10**(-99999.) returns
    ! a value of zero.
    nerr = 0

    do kcol = 1,kbt
        if (zvclg1(kcol) .le. -99999.) then
            ! Calling sequence substitutions:
            !   uzvec1 for unam48
            call fmspnx(jlen,uzvec1(kcol),uspn56)

            if (kcol .ne. krdxsp) then
                write (noutpt,1050) uspn56(1:jlen),zvclg1(kcol)
                write (nttyo,1050) uspn56(1:jlen),zvclg1(kcol)
1050 format(/' * Warning - (EQ6/eqphas) The log number of moles',' variable for basis species',/7x,a,' has the out-of-range',' value ',1pe12.5,'.')
            else if (qredox) then
                if (krdxsp .eq. ko2gaq) then
                    write (noutpt,1060) zvclg1(kcol)
                    write (nttyo,1060) zvclg1(kcol)
1060 format(/' * Warning - (EQ6/eqphas) The log oxygen',' fugacity variable has the out-of-range',/7x,'value ',1pe12.5,'.')
                else if (krdxsp .eq. kelect) then
                    write (noutpt,1070) zvclg1(kcol)
                    write (nttyo,1070) zvclg1(kcol)
1070 format(/' * Warning - (EQ6/eqphas) The log electron',' activity variable has the',/7x,'out-of-range value',' ',1pe12.5,'.')
                else
                    write (noutpt,1050) uspn56(1:jlen),zvclg1(kcol)
                    write (nttyo,1050) uspn56(1:jlen),zvclg1(kcol)
                end if
            end if

            nerr = nerr + 1
        end if
    end do

    if (nerr .gt. 0) then
        ier = 100
        go to 999
    end if

    ! Check the current phase assemblage for mineral species for which
    ! the log number of moles variable is -99999. This implies that the
    ! corresponding number of moles variable is zero, and in turn that
    ! the species is not actually present in the equilibrium system.
    ! In such a case, the corresponding phase must be deleted from the
    ! current phase assemblage.
    do kcol = km1,kxt
        if (zvclg1(kcol) .le. -99999.) then
            npdel = ipndx1(kcol)
            go to 700
        end if
    end do

    call eqcalc(aamatr,abar,acflg,acflgo,act,actlg,adh,adhh,adhv,afcnst,alpha,al10,amtb,aphi,apx,avcnst,azero,a3bar,a3bars,bacfmx,bbig,bdh,bdhh,bdhv,bdot,bdoth,bdotv,beta,betamx,betao,bfje,bfxi,bgamx,bneg,bpx,bsigmm,cco2,cegexs,cess,cdrs,cdrsd,cdrsx,cdrtw,cdrw,cjbasp,cnufac,conc,conclg,cpgexs,cscale,csts,delmax,delvco,delvec,dlogxw,egexjc,egexjf,egexs,eh,ehfac,eps100,farad,fje,fjeo,fo2,fo2lg,fsort,fugac,fugalg,fxi,fxio,gmmatr,iapxt,ibetmx,ibpxt,ibswx,idelmx,ielam,ier,iern1,iern2,ifcphi1,ifcphi2,ifnnn,ifn2n,ifpsi1,ifpsi2,ifrn1,ifrn2,ifzeta,igas,igstak,iindx1,ilcphi1,ilcphi2,ilnnn,iln2n,ilpsi1,ilpsi2,ilrn1,ilrn2,ilzeta,imrn1,imrn2,insgf,iodb,iopg,iopt,ipch,ipcv,ipivot,ipndx1,istack,iter,itermx,ixbasp,ixrn1,ixrn2,izmax,jcsort,jflag,jgsort,jgstak,jjsort,jpflag,jsflag,jsitex,jsol,jssort,jstack,kbt,kction,kdim,kelect,khydr,khydx,km1,kmt,ko2gaq,krdxsp,kwater,kx1,kxt,loph,losp,lsort,moph,mosp,mrgexs,mtb,narn1,narn2,narxt,nat,nbasp,nbaspd,nbaspx,nbt,nbtd,nbw,nchlor,ncmpr,nct,ndrs,ndrsd,ndrsx,ndrsr,ndrsrd,ndrsrx,nelect,nern1,nern2,ness,nessr,net,nfac,nfrn1,nfrn2,ngrn1,ngrn2,ngt,nhydr,nhydx,nlrn1,nlrn2,nlt,nmrn1,nmrn2,nmt,noutpt,no2gaq,nphasx,npt,nrdxsp,nst,nsts,nstsr,ntpr,nttyo,nxrn1,nxrn2,nxt,omega,omeglg,press,qbassw,qhawep,qoptmz,qpit75,qredox,q6mode,rhsvec,screwd,sigmam,sigmmo,smp100,tempc,tempk,tolbt,toldl,ubacmx,ubgamx,ulbeta,uldel,uphase,uspec,uzvec1,weight,wfac,xbar,xbarlg,xbarw,xbarwc,xbrwlc,xbrwlg,zchar,zchcu6,zchsq2,zvclg1,zvec1)

    ! Were any basis switches made by EQ6/optmzr.f?
    nswtch = 0

    do kcol = 1,kbt
        nb = iindx1(kcol)
        ns = nbasp(nb)
        ns2 = nbasps(nb)

        if (ns2 .ne. ns) then
            nswtch = nswtch + 1
        end if
    end do

    qbseqc = nswtch .gt. 0

    if (qbseqc) then
        ! Drop the order to zero. The finite difference data are tied
        ! to the old basis set. The original order and the use of these
        ! data can be restored later in the present subroutine if the
        ! switches just made are undone.
        nord = nords
    end if

    if (ier .gt. 0) then
        go to 500
    end if

    ! Hybrid Newton-Raphson iteration has converged. Check for super-
    ! saturations (other than those which are permitted to exist
    ! metastably).
    ifail = 0
    npadd = 0
    nsspt = 0
    call initiz(nssp,nsspmx)

    ! Calculate the total number of moles of each basis species
    ! present in the aqueous phase.
    call initaz(mtbaq,nbt)

    do nss = narn1,narn2
        ns = jcsort(nss)

        if (mosp(ns) .ne. 0.) then
            nr1 = nstsr(1,ns)
            nr2 = nstsr(2,ns)

            do n = nr1,nr2
                nb = nsts(n)
                mtbaq(nb) = mtbaq(nb) + csts(n)*mosp(ns)
            end do
        end if
    end do

    ! Check for supersaturations.
    call satchk(acflg,act,actlg,afcnst,affp,affs,apx,bpx,cdrs,eps100,iindx1,iodb,iopt,iapxmx,ibpxmx,iktmax,ixrn1,jflag,jpflag,jsflag,jsol,kmax,km1,kpsat,kpsst,kxt,nbasp,nbt,nbtmax,ncmpr,ndrs,ndrsmx,ndrsr,nodbmx,noptmx,noutpt,npchk,npt,nptmax,nstmax,nttyo,nxrn1,nxrn2,nxtmax,qxknph,sidrph,sidrsp,tolsat,uphase,uspec,wfac,xbar,xbarlg,xlks)

    if (kpsst .gt. 0) then
        if (iodb(4) .le. 0) then
            ux8 = ' '
            write (ux8,'(i7)') iter
            call lejust(ux8)
            j2 = ilnobl(ux8)
            write (noutpt,1090) ux8
1090 format(' iter= ',a)
        end if

        ux8 = ' '
        write (ux8,'(i7)') kpsst
        call lejust(ux8)
        j2 = ilnobl(ux8)
        write (noutpt,1100) ux8(1:j2)
1100 format(/'   Have ',a,' supersaturated phases.')
    end if

    if (kpsst .le. 0) then
        ! The last equilibrium calculation converged and there are no
        ! supersaturations to eliminate. The calculation overseen by
        ! the present subroutine has completed successfully.
        go to 990
    end if

    ! Have one or more supersaturations. Find a phase to add to the
    ! current phase assemblage. To do this, scale the affinities of the
    ! supersaturated phases according to the number of ions produced
    ! and destroyed in the respective dissolution reactions. The
    ! scaled affinity is arbitrarily reduced if the maximum possible
    ! number of moles of the new phase as estimated from the aqueous
    ! solution chemistry is very small. The supersaturated phases are
    ! arranged in decreasing order of scaled affinity. The phase with
    ! the greatest scaled affinity is added to the phase assemblage.
    ! In general, this algorithm picks the correct phase about four
    ! out of five times.
    call initcb(ussp,nsspmx)
    call initaz(afssp,nsspmx)
    call initaz(afssps,nsspmx)
    call initaz(msspmx,nsspmx)

    isspt = 0

    do np = 1,npt
        if (jpflag(np) .eq. -2) then
            jpflag(np) = 0
            isspt = isspt + 1

            cxx = 0.
            nxx = 0
            nrn1 = ncmpr(1,np)
            nrn2 = ncmpr(2,np)
            nt = nrn2 - nrn1 + 1
            mpmxaq = 0.
            mpmxes = 0.

            do ns = nrn1,nrn2
                if (jsflag(ns) .le. 0) then
                    cxx = cxx + cscale(ns)
                    nxx = nxx + 1
                    msmxaq = 1.e+38
                    msmxes = 1.e+38
                    nr1 = nstsr(1,ns)
                    nr2 = nstsr(2,ns)

                    do n = nr1,nr2
                        cx = csts(n)

                        if (cx .ne. 0.) then
                            nb = nsts(n)
                            nss = nbaspd(nb)

                            if (nss.ne.nhydr .and. nss.ne.nhydx .and.              nss.ne.no2gaq .and. nss.ne.nelect) then
                                ! Note: the total number of moles of a basis species
                                ! can not be used to bound the number of moles of
                                ! the species currently being examined (and hence
                                ! also the number of moles of the corresponding phase)
                                ! if the current basis species is H+, OH-, O2(g,aq),
                                ! or e-. The latter two species are fictive. More to
                                ! the point, the total number of moles of such
                                ! basis species may be zero or even a negative number.
                                if (cx .gt. 0.) then
                                    if (mtbaq(nb) .gt. 0.) then
                                        mx = mtbaq(nb)/cx
                                        msmxaq = min(msmxaq,mx)
                                    end if

                                    if (mtb(nb) .gt. 0.) then
                                        mx = mtb(nb)/cx
                                        msmxes = min(msmxes,mx)
                                    end if
                                end if
                            end if
                        end if
                    end do

                    if (msmxaq .ge. 1.e+38) then
                        msmxaq = 1.0
                    end if

                    if (msmxes .ge. 1.e+38) then
                        msmxes = 1.0
                    end if

                    mpmxaq = max(mpmxaq,msmxaq)
                    mpmxes = max(mpmxes,msmxes)
                end if
            end do

            ! Pick a maximum number of moles for the current phase based
            ! on the moles of components in the aqueous solution and in
            ! the equilibrium system.
            mophmx = max(mpmxaq,(0.05*mpmxes))

            ! If the phase consists of more than one species, compute and
            ! use the average affinity scaling factor.
            cxx = cxx/nxx
            afscal = affp(np)/cxx

            ! Reduce the scaled affinity if the maximum precipitable
            ! number of moles is small.
            if (mpmxaq .le. 5.e-6) then
                afscal = 0.01*afscal
            end if

            do jj = 1,nsspmx
                if (afscal .ge. afssps(jj)) then
                    jkl = nsspmx - jj

                    do j = 1,jkl
                        k = nsspmx - j
                        nssp(k + 1) = nssp(k)
                        ussp(k + 1) = ussp(k)
                        afssps(k + 1) = afssps(k)
                        afssp(k + 1) = afssp(k)
                        msspmx(k + 1) = msspmx(k)
                    end do

                    nssp(jj) = np
                    ussp(jj) = uphase(np)
                    afssps(jj) = afscal
                    afssp(jj) = affp(np)
                    msspmx(jj) = mophmx
                    go to 200
                end if
            end do

200 continue
        end if
    end do

    nsspt = min(isspt,nsspmx)

    if (iodb(1) .ge. 2) then
        write (noutpt,1110)
1110 format(//5x,'The most supersaturated phases:',//31x,'Scaled',6x,'Affinity,',5x,'Maximum',/3x,'Name',23x,'Affinity',6x,'kcal/mol',4x,'Precipitable',/59x,'Moles',/)

        do n = 1,nsspt
            write (noutpt,1120) ussp(n),afssps(n),afssp(n),msspmx(n)
1120 format(1x,a24,2x,f12.7,2x,f12.7,2x,1pe12.5)
        end do
    end if

    if (afssp(1) .gt. tolsst) then
        if (qsspgb .and. ntry.le.1) then
            ! Return a special error flag so that subroutine path will cut
            ! the step size in order to be able to home in on the value of
            ! reaction progress which corresponds to the target
            ! supersaturation.
            write (noutpt,1150)
1150 format(/' --- The extent of supersaturation exceeds the',' normal tolerance ---')

            ier = 10
            go to 999
        end if
    end if

    npadd = nssp(1)
    mxmadd = msspmx(1)

    if (npadd.eq.npaddi .and. ntry.eq.(ntryai + 2)) then
        ! The phase chosen to add was added two steps previously and
        ! removed on the previous step.
        if (iodb(1) .gt. 0) then
            write (noutpt,1160)
            write (nttyo,1160)
1160 format(/' * Note - (EQ6/eqphas) Caught in a region of',/7x,'critical instability in the ES phase assemblage. A',/7x,"supersaturated phase can't be precipitated due to",/7x,'lack of convergence when it is included in the ES',/7x,'phase assemblage.')
        end if

        ier = 30
        go to 999
    end if

    ! The following is a return point to add one new member to the
    ! current phase assemblage and try again.
300 continue
    j2 = ilnobl(uphase(npadd))
    write (ux8,'(i7)') npadd
    call lejust(ux8)
    j3 = ilnobl(ux8)
    write (noutpt,1180) uphase(npadd)(1:j2),ux8(1:j3)
1180 format(/'   The phase to be added is ',a,1x,'(',a,')')

    qadd = .true.
    npaddi = npadd
    ntryai = ntry

    do kcol = 1,kmax
        delvec(kcol) = 0.
    end do

    ! The starting number of moles is set to some fraction of the
    ! probable maximum value. Here it is set to 100% of that.
    mophn = mxmadd
    np = npadd
    jpflag(np) = -1
    loph(np) = tlg(mophn)
    nr1 = ncmpr(1,np)
    nr2 = ncmpr(2,np)
    nt = nr2 - nr1 + 1

    if (nt .eq. 1) then
        jsflag(nr1) = -1
        losp(nr1) = loph(np)
    else
        do ns = nr1,nr2
            if (jsflag(ns) .le. 0) then
                jsflag(ns) = -1
                losp(ns) = loph(np) + xbarlg(ns)
            end if
        end do
    end if

    ! Modify the current phase assemblage.
400 continue
    ntry = ntry + 1

    if (ntry .gt. ntrymx) then
        write (noutpt,1200) ntrymx
        write (nttyo,1200) ntrymx
1200 format(/' * Note - (EQ6/eqphas) Have done the maximum ',i3,' tries to find the',/7x,'correct phase assemblage without',' succeeding. Suggest increasing the',/7x,'value of the ntrymx',' variable on the input file.')

        ier = 20
        go to 999
    end if

    qmod = .true.
    kpsst = 0

    ! Modify the indexing of the Jacobian.
    call miidxz(ier,iindx1,ipndx1,jpflag,jsflag,kbt,kdim,kmax,km1,kmt,kx1,kxt,losp,ncmpr,noutpt,npt,nptmax,nstmax,nttyo,uspec,uzvec1,zvclg1,zvec1)

    if (ier .gt. 0) then
        ier = 8
        go to 999
    end if

    write (noutpt,1020)

    ! Go back and try to make the equilibrium calculation with the
    ! new phase assemblage.
    go to 100

500 continue

    ! Iteration has terminated without achieving convergence. Attempt to
    ! diagnose the cause. It may be that a phase must be removed from
    ! the current phase assemblage.
    ifail = ifail + 1

    ! If there are no mineral phases in the equilibrium system phase
    ! assemblage to consider deleting, look for another cause for
    ! the failure to converge.
    if (kxt .le. kbt) then
        go to 980
    end if

    ! Try to find a phase do delete from the phase assemblage.
    npdel = 0

    ! Trap out-of-range values for the variables associated with the
    ! mineral species. An out-of-range value ("-99999.") is likely to
    ! be encountered only as a result of having read it from the input
    ! file.
    do kcol = km1,kxt
        ! If a value of -99999. is found for the log number of moles
        ! variable of any mineral species, go remove the associated
        ! phase from the current phase assemblage and try again.
        if (zvclg1(kcol) .le. -99999.) then
            npdel = ipndx1(kcol)
            go to 700
        end if
    end do

    if ((qadd .or. qstart) .and. (kxt - km1 + 1).ge.2) then
        ! Check for the possibility of a mineralogic phase rule violation.
        ! The mineralogic phase rule is 'f = c - p', whereas the phase
        ! rule of thermodynamics is 'f = c - p + 2'. The mineralogic
        ! phase rule applies here because the present computations
        ! are made for fixed temperature and pressure. A violation of
        ! the mineralogic phase rule is equivalent to a condition of
        ! linear dependence among the mineral mass action rows
        ! (krow = km1,kxt).
        irow1 = km1
        irow2 = kmt

        ! XX     irow2 = kxt
        jcol1 = 1
        jcol2 = kmt

        ! XX     jcol2 = kxt
        !        Calling sequence substitutions:
        !          qphruv for qldep
        call lindep(aamatr,eps100,irow1,irow2,jcol1,jcol2,kmax,qphruv)

        if (qphruv) then
            write (noutpt,1400)
            write (nttyo,1400)
1400 format(/' * Note - (EQ6/eqphas) The current phase assemblage',' violates the',/7x,'mineralogic phase rule.')

            if (iopt(1).eq.2 .and. .not.qstart .and.      .not.(qcntmp.and.qcnpre)) then
                ! A mineralogic phase rule violation occurred while running
                ! a non-isothermal fluid-centered flow-through open system
                ! simulation. The violation may be due to crossing a
                ! univariant curve.
                write (noutpt,1500)
1500 format('   --- Possibly attempting to cross a univariant ','curve ---',/)
            end if

            ! Identify a phase which must be deleted in order to avoid the
            ! violation. Ordinarily, there would only be one.
            call jgibbs(aamatr,afcnst,affp,cdrs,cscale,csts,delvec,eps100,gmmatr,iindx1,iodb,ipivot,ipndx1,jpflag,kbt,kdim,kmax,km1,kmt,kx1,kxt,mtb,nbasp,nbtmax,ndrs,ndrsmx,ndrsr,nodbmx,noutpt,npadd,npdel,nptmax,nstmax,nsts,nstsmx,nstsr,nttyo,rhsvec,uphase,uspec,xlks)

            if (npdel .gt. 0) then
                ! Go delete the offending phase and try again.
                go to 700
            end if
        end if
    end if

    ! Try to find a phase to delete using the following two empirical
    ! criteria:
    !    1. The log number of moles of a phase was becoming very
    !       negative (or was already very negative if iter .le. 1).
    !    2. Its derivative with respect to Xi is very negative
    !       (must have nord.ge.1 in order to apply this).
    call phsdrp(d1zvc1,iindx0,iindx1,iodb,ipndx1,iter,kmax,km1,km1s,kxt,kxts,nodbmx,nord,noutpt,npadd,npdel,nptmax,ntry,uphase,zvclgs,zvclg1)

    if (npdel .le. 0) then
        ! Couldn't find a phase to delete.
        go to 980
    end if

700 continue

    ! Have found a phase to delete from the phase assemblage.
    if (npdel.eq.npdeli .and. ntry.eq.(ntrydi + 2)) then
        ! The phase chosen to delete was deleted two steps previously and
        ! then added on the previous step.
        if (iodb(1) .gt. 0) then
            write (noutpt,1510)
            write (nttyo,1510)
1510 format(/' * Note - (EQ6/eqphas) Caught in a region of',/7x,'critical instability in the ES phase assemblage. A',/7x,"phase can't be deleted because when it is not included",/7x,'in the ES phase assemblage the system is supersaturated',/7x,'beyond the allowed tolerance.')
        end if

        ier = 40
        go to 999
    end if

    ! Delete a phase from the phase assemblage and try again.
    ier = 0
    npdeli = npdel
    ntrydi = ntry

    j2 = ilnobl(uphase(npdel))
    write (ux8,'(i7)') npdel
    call lejust(ux8)
    j3 = ilnobl(ux8)
    write (noutpt,1520) uphase(npdel)(1:j2),ux8(1:j3)
1520 format(/'   The phase to be dropped is ',a,1x,'(',a,')')

    if (qbseqc) then
        ! If any basis switches were done by EQ6/optmzr.f, undo them.
        ! Restore the kernel for the equilibrium sytem from the backup.
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
    end if

    ! Reset the log number of moles variables (losp) for all mineral
    ! species.
    do kcol = km1,kxt
        ns = iindx1(kcol)
        losp(ns) = zvclg1(kcol)
    end do

    ! Delete the desired phase.
    np = npdel
    jpflag(np) = 0
    loph(np) = -99999.
    moph(np) = 0.
    nr1 = ncmpr(1,np)
    nr2 = ncmpr(2,np)

    do ns = nr1,nr2
        if (jsflag(ns) .le. 0) then
            jsflag(ns) = 0
            losp(ns) = -99999.
            mosp(ns) = 0.
        end if
    end do

    ! If the phase to be deleted was the last one added, it is replaced
    ! by the phase which had the next highest scaled affinity (as long
    ! as such can be found in the limited list maintained for this
    ! purpose). Otherwise, the phase to be deleted is just removed.
    qadd = .false.

    if (npdel .eq. npaddi) then
        inext = ifail + 1

        if (inext .le. nsspt) then
            npadd = nssp(inext)

            ! Replace the phase to be deleted with another.
            go to 300
        else
            ! Just delete the phase.
            ! Note: from here the calculation needs to slide forward
            ! in delxi in order to get past a region of computational
            ! instability associated with a phase appearance boundary.
            go to 400
        end if
    else
        ! Just delete the phase.
        go to 400
    end if

980 continue

    ! Test for critical redox instability.
    if (qredox) then
        if (idelmx.eq.krdxsp .or. abs(delvec(krdxsp)).ge.2.0) then
            write (noutpt,2000)
            write (nttyo,2000)
2000 format(/' * Note - (EQ6/eqphas) Caught in a region of',' critical redox instability.',/7x,'The value of the redox',' master variable (log fO2, Eh, pe-, etc.)',/7x,'is',' essentially indifferent to the constraints defining',' the state',/7x,'of the ES.')

            ier = 50
            go to 999
        end if
    end if

    ! Couldn't find a phase to delete from the current phase assemblage.
    ! Look for other explanations for the failure to converge.
    ! Test for critical instability associated with solvent water.
    if (kwater .gt. 0) then
        if (idelmx.eq.kwater .or. abs(delvec(kwater)).ge.2.0) then
            write (noutpt,2010)
            write (nttyo,2010)
2010 format(/' * Note - (EQ6/eqphas) Caught in a region of',' critical instability associated',/7x,'with solvent water.')

            if (zvclg1(kwater) .le. -14.0) then
                write (noutpt,2020)
                write (nttyo,2020)
2020 format(/7x,'There is almost no solvent water present.')
            end if

            ier = 60
            go to 999
        end if
    end if

    ! Test for instability that appears to be associated with the
    ! aqueous activity coefficient model.
    if (bgamx.ge.1.0 .or. bfxi.ge.1.0 .or. bsigmm.ge.1.0) then
        write (noutpt,2050) bgamx,fxi,fxio,sigmam,sigmmo
        write (nttyo,2050) bgamx,fxi,fxio,sigmam,sigmmo
2050 format(/' * Note - (EQ6/eqphas) Caught in a region of',' instability that may be',/7x,'associated with the',' aqueous activity coefficient model. The max',/7x,'norm of',' the activity coefficients for the last change was ',1pe12.5,/7x,'The last value of the ionic strength in the iteration',' was',/7x,e12.5,' molal, the immediately preceding value',' was',/7x,e12.5,'. The last value of the sum of the solute',' molalities in',/7x,'the iteration was ',e12.5,', the',' immediately preceding value was',/7x,e12.5,'.')

        if (zvclg1(kwater) .le. -14.0) then
            write (noutpt,2060)
            write (nttyo,2060)
2060 format(/7x,'There is almost no solvent water present.')
        end if

        ier = 80
        go to 999
    end if

    ! Trap out-of-range values for the variables associated with the
    ! basis species.
    nerr = 0

    do kcol = 1,kbt
        if (zvclg1(kcol) .le. -99999.) then
            ! Calling sequence substitutions:
            !   uzvec1 for unam48
            call fmspnx(jlen,uzvec1(kcol),uspn56)

            if (kcol .ne. krdxsp) then
                write (noutpt,2100) uspn56(1:jlen),zvclg1(kcol)
                write (nttyo,2100) uspn56(1:jlen),zvclg1(kcol)
2100 format(/' * Warning - (EQ6/eqphas) The log number of moles',' variable for basis species',/7x,a,' has the out-of-range',' value ',1pe12.5,/7x,'after an iteration crash.')
            else if (qredox) then
                if (krdxsp .eq. ko2gaq) then
                    write (noutpt,2110) zvclg1(kcol)
                    write (nttyo,2110) zvclg1(kcol)
2110 format(/' * Warning - (EQ6/eqphas) The log oxygen',' fugacity variable has the out-of-range',/7x,'value ',1pe12.5,' after an iteration crash.')
                else if (krdxsp .eq. kelect) then
                    write (noutpt,2120) zvclg1(kcol)
                    write (nttyo,2120) zvclg1(kcol)
2120 format(/' * Warning - (EQ6/eqphas) The log electron',' activity variable has the',/7x,'out-of-range value',' ',1pe12.5,' after an iteration crash.')
                else
                    write (noutpt,2100) uspn56(1:jlen),zvclg1(kcol)
                    write (nttyo,2100) uspn56(1:jlen),zvclg1(kcol)
                end if
            end if

            nerr = nerr + 1
        end if
    end do

    if (nerr .gt. 0) then
        ier = 110
        go to 999
    end if

    ! Test simply for almost no solvent water present. This can
    ! sometimes be the problem even if some of the other tests
    ! applied above fail to detect it.
    if (kwater .gt. 0) then
        if (zvclg1(kwater) .le. -14.0) then
            write (noutpt,2030)
            write (nttyo,2030)
2030 format(/' * Note - (EQ6/eqphas) There appears to be almost',' no solvent water present.')

            ier = 70
            go to 999
        end if
    end if

    ! Okay, can't diagnose what caused the iteration to crash.
    write (noutpt,2150)
    write (nttyo,2150)
2150 format(/' * Note - (EQ6/eqphas) An equilibrium calculation',' failed for reasons',/7x,"that couldn't be diagnosed. Will",' try to recover.')

    ier = 150
    go to 999

990 continue

    ! Check to see if the new phase assemblage is the same as the old
    ! one. If qmod is true at this point, it only means that the
    ! putative phase assemblage has changed in the course of the
    ! equilibrium calculations. It may, however, have ended up being
    ! the original assemblage.
    if (qmod) then
        if (km1.eq.km10 .and. kmt.eq.kmt0) then
            if (kx1.eq.kx10 .and. kxt.eq.kxt0) then
                do kcol = km1,kxt
                    if (iindx1(kcol) .ne. iindx0(kcol)) then
                        go to 995
                    end if
                end do

                qmod = .false.
            end if
        end if
    end if

995 continue

999 continue
end subroutine eqphas