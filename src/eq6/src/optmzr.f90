subroutine optmzr(aamatr,abar,acflg,acflgo,act,actlg,adh,adhh,adhv,afcnst,al10,alpha,amtb,aphi,avcnst,azero,a3bar,a3bars,bacfmx,bbig,bdh,bdhh,bdhv,bdot,bdoth,bdotv,beta,betamx,bgamx,bneg,bpx,cco2,cdrs,cdrsx,cdrtw,cdrw,cegexs,cgexj,cjbasp,cnufac,conc,conclg,cpgexs,cscale,csts,coval,delvec,dlogxw,egexjc,egexjf,egexs,ehfac,eps100,fje,fjeo,fo2,fo2lg,fsort,fugac,fugalg,fxi,fxio,gmmatr,ibetmx,ibpxt,ibswx,ielam,iern1,iern2,ifcphi1,ifcphi2,ifnnn,ifn2n,ifpsi1,ifpsi2,ifrn1,ifrn2,ifzeta,igas,igstak,iindx1,ilcphi1,ilcphi2,ilnnn,iln2n,ilpsi1,ilpsi2,ilrn1,ilrn2,ilzeta,imrn1,imrn2,insgf,iodb,iopg,iopt,ipch,ipcv,ipivot,ipndx1,irdxc3,istack,ixbasp,ixrn1,ixrn2,izmax,jcsort,jern1,jern2,jflag,jgext,jgsort,jgstak,jjsort,jpflag,jsflag,jsitex,jsol,jssort,jstack,kbt,kction,kdim,kelect,khydr,km1,kmt,ko2gaq,kwater,kx1,kxt,loph,losp,lsort,mgext,moph,mosp,mrgexs,mtb,narn1,narn2,narxt,nat,nbasp,nbaspd,nbaspx,nbw,nbt,nbtd,nchlor,ncmpr,ncosp,ndrs,ndrsx,ndrsr,ndrsrd,ndrsrx,nelect,nern1,nern2,net,ngexsa,nfac,ngext,nhydr,nhydx,ngrn1,ngrn2,ngt,noutpt,no2gaq,nphasx,npt,nst,nsts,nstsr,ntfx,ntfxt,ntpr,nttyo,omega,omeglg,press,qbassw,qblamx,qhawep,qpit75,qredox,q6mode,rhsvec,sigmam,sigmmo,smp100,tempc,tempk,tfx,tolbig,tolneg,tolxpt,ubacmx,ubbig,ubetmx,ubgamx,ubneg,ugexj,ugexmo,uphase,uspec,uzvec1,weight,wfac,xbar,xbarlg,xbarw,xbarwc,xbrwlc,xbrwlg,zchar,zchcu6,zchsq2,zgexj,zvclg1,zvec1)
    !! This subroutine optimizes the starting values prior to hybrid
    !! Newton-Raphson iteration. The structure is somewhat similar to
    !! that in EQ3NR/arrset.f. Automatic basis switching is done in
    !! "loops." Inside loops are "passes," in which the activity
    !! coefficients are readjusted. Inside passes are "cycles," in which
    !! the chief iteration variables are adjusted.
    !! In the present subroutine, the cycle algorithm is based on
    !! minimizing the function:
    !!   aleph = Sum(i) zeta(i)**2
    !! where i spans the set of primary iteration variables and zeta
    !! is a residual vector derived from the beta residual vectors
    !! normally used in EQ3/6. In EQ3NR/arrset.f, the cycle algorithm
    !! is the continued-fraction algorithm. That algorithm is not
    !! well suited to dealing with nonpphysical mass balances,
    !! which must be dealt with for basis species like H2O(l), H+,
    !! O2(g,aq), and e-. Nor is it well-suited for dealing with
    !! non-aqueous species (e.g., pure minerals) when they are not
    !! switched into the active basis set but are instead treated as
    !! "extended" basis species. Hence the change of algorithm.
    !! This subroutine does not change the presumed phase assemblage.
    !! This subroutine is called by:
    !!   EQ6/eqcalc.f
    !! Principal input:
    !! Principal output:
    !!   nhydr  = species index of the hydrogen ion
    !!   nhydx  = species index of the hydroxide ion
    !!   nchlor = species index of the chloride ion
    !!   qblamx = .true. if solid solutions are in the matrix
    !! Modules.
    !! The module mod6pt contains data required to evaluate Pitzer's
    !! equations.
    use mod6pt

    ! The module mod6xf contains most of the standard-state
    ! thermodynamic data.
    use mod6xf

    implicit none

    include 'eqlib/eqldv.h'

    ! Calling sequence variable declarations.
    integer :: noutpt
    integer :: nttyo

    integer :: ibpxt(nxtmax)
    integer :: ibswx(nbtmax)
    integer :: igstak(ngtmax)
    integer :: iindx1(kmax)
    integer :: insgf(natmax)
    integer :: iodb(nodbmx)
    integer :: iopg(nopgmx)
    integer :: iopt(noptmx)
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
    integer :: narxt(ntprmx)
    integer :: nbasp(nbtmax)
    integer :: nbaspd(nbtmax)
    integer :: nbaspx(nbtmax)
    integer :: ncmpr(2,nptmax)
    integer :: ncosp(nbtmax)
    integer :: ndrs(ndrsmx)
    integer :: ndrsx(ndrsmx)
    integer :: ndrsr(2,nstmax)
    integer :: ndrsrd(2,nstmax)
    integer :: ndrsrx(2,nstmax)
    integer :: ngexsa(ietmax,jetmax,netmax)
    integer :: ngext(jetmax,netmax)
    integer :: nphasx(nstmax)
    integer :: nsts(nstsmx)
    integer :: nstsr(2,nstmax)
    integer :: ntfx(ntfxmx)

    integer :: nfac(nbtmax)

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

    integer :: ibetmx
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
    integer :: nat
    integer :: nbw
    integer :: nbt
    integer :: nbtd
    integer :: nchlor
    integer :: nelect
    integer :: nern1
    integer :: nern2
    integer :: net
    integer :: nhydr
    integer :: nhydx
    integer :: ngrn1
    integer :: ngrn2
    integer :: ngt
    integer :: no2gaq
    integer :: npt
    integer :: nst
    integer :: ntpr
    integer :: ntfxt

    logical :: qabsw
    logical :: qbassw
    logical :: qblamx
    logical :: qhawep
    logical :: qpit75
    logical :: qredox
    logical :: q6mode

    character(len=48) :: uspec(nstmax)
    character(len=48) :: uzvec1(kmax)
    character(len=48) :: ubacmx
    character(len=48) :: ubbig
    character(len=48) :: ubetmx
    character(len=48) :: ubgamx
    character(len=48) :: ubneg
    character(len=24) :: ugexmo(netmax)
    character(len=24) :: uphase(nptmax)
    character(len=8) :: ugexj(jetmax,netmax)

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
    real(kind=8) :: cscale(nstmax)
    real(kind=8) :: csts(nstsmx)

    real(kind=8) :: delvec(kmax)
    real(kind=8) :: dlogxw(nbtmax)
    real(kind=8) :: egexjc(jetmax,netmax)
    real(kind=8) :: egexjf(jetmax,netmax)
    real(kind=8) :: egexs(ietmax,jetmax,netmax)
    real(kind=8) :: fsort(ngtmax)
    real(kind=8) :: fugac(ngtmax)
    real(kind=8) :: fugalg(ngtmax)
    real(kind=8) :: gmmatr(kmax,kmax)
    real(kind=8) :: mgext(jetmax,netmax)
    real(kind=8) :: mrgexs(ietmax,jetmax,netmax)
    real(kind=8) :: mtb(nbtmax)
    real(kind=8) :: loph(nptmax)
    real(kind=8) :: losp(nstmax)
    real(kind=8) :: lsort(nstmax)
    real(kind=8) :: moph(nptmax)
    real(kind=8) :: mosp(nstmax)
    real(kind=8) :: rhsvec(kmax)
    real(kind=8) :: tfx(ntfxmx)
    real(kind=8) :: weight(nstmax)
    real(kind=8) :: wfac(iktmax,nxtmax)
    real(kind=8) :: xbar(nstmax)
    real(kind=8) :: xbarlg(nstmax)
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
    real(kind=8) :: afcnst
    real(kind=8) :: al10
    real(kind=8) :: avcnst
    real(kind=8) :: a3bar
    real(kind=8) :: bacfmx
    real(kind=8) :: bbig
    real(kind=8) :: betamx
    real(kind=8) :: bgamx
    real(kind=8) :: bneg
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
    real(kind=8) :: sigmam
    real(kind=8) :: sigmmo
    real(kind=8) :: smp100
    real(kind=8) :: tempc
    real(kind=8) :: tempk
    real(kind=8) :: tolbig
    real(kind=8) :: tolneg
    real(kind=8) :: tolxpt
    real(kind=8) :: xbarw
    real(kind=8) :: xbarwc
    real(kind=8) :: xbrwlc
    real(kind=8) :: xbrwlg
    real(kind=8) :: zetsqm
    real(kind=8) :: zetsqo

    ! Local variable declarations with global dimensioning.
    integer, dimension(:), allocatable :: iastak
    integer, dimension(:), allocatable :: jasort
    integer, dimension(:), allocatable :: jastak
    integer, dimension(:), allocatable :: kmvar
    integer, dimension(:), allocatable :: knflag

    real(kind=8), dimension(:), allocatable :: daleph
    real(kind=8), dimension(:), allocatable :: delprc
    real(kind=8), dimension(:), allocatable :: dijmaj
    real(kind=8), dimension(:), allocatable :: zesort
    real(kind=8), dimension(:), allocatable :: zeta
    real(kind=8), dimension(:), allocatable :: zetasq

    real(kind=8), dimension(:), allocatable :: efac

    integer :: isv_nbtmax
    integer :: isv_kmax

    SAVE isv_nbtmax,isv_kmax

    SAVE daleph,delprc,dijmaj,efac,iastak,jasort,jastak,kmvar,knflag,zesort,zeta,zetasq

    ! Local variable declarations.
    integer :: icorr
    integer :: kdirec
    integer :: iebal
    integer :: iter
    integer :: j
    integer :: jdim
    integer :: jlen
    integer :: jlen1
    integer :: jlen2
    integer :: j2
    integer :: k
    integer :: kbig
    integer :: kcol
    integer :: kcorr
    integer :: kcscan
    integer :: ker
    integer :: kk
    integer :: kkdim
    integer :: kount
    integer :: krow
    integer :: krscan
    integer :: n
    integer :: nb
    integer :: nbigger
    integer :: nchange
    integer :: ncycle
    integer :: ncylim
    integer :: negafc
    integer :: nloop
    integer :: nlopmx
    integer :: npass
    integer :: nplim
    integer :: nredox
    integer :: ns
    integer :: nscan
    integer :: nsd
    integer :: nswtch
    integer :: ns1
    integer :: ns2

    integer :: ilnobl

    logical :: qbswx
    logical :: qcgam
    logical :: qcfxi
    logical :: qcsigm
    logical :: qcycnc
    logical :: qpracf
    logical :: qscan2
    logical :: qsclim
    logical :: qscosc
    logical :: qstops
    logical :: qtestc
    logical :: qtestp
    logical :: qxbarw

    character(len=56) :: uspn56
    character(len=56) :: usp156
    character(len=56) :: usp256
    character(len=48) :: uzebig
    character(len=48) :: uzeneg
    character(len=24) :: ujtp

    real(kind=8) :: abig
    real(kind=8) :: actwlc
    real(kind=8) :: adij
    real(kind=8) :: alefnc
    real(kind=8) :: alepoe
    real(kind=8) :: alepoo
    real(kind=8) :: aleph
    real(kind=8) :: alephm
    real(kind=8) :: av
    real(kind=8) :: ax
    real(kind=8) :: azex
    real(kind=8) :: betfnc
    real(kind=8) :: bfje
    real(kind=8) :: bfxi
    real(kind=8) :: bsigmm
    real(kind=8) :: btmxoe
    real(kind=8) :: btmxoo
    real(kind=8) :: bx
    real(kind=8) :: chfacf
    real(kind=8) :: chfsgm
    real(kind=8) :: cx
    real(kind=8) :: dij
    real(kind=8) :: dscan
    real(kind=8) :: dscmin
    real(kind=8) :: dx
    real(kind=8) :: dzx
    real(kind=8) :: dzxclm
    real(kind=8) :: eh
    real(kind=8) :: eps1hi
    real(kind=8) :: rlxgam
    real(kind=8) :: tolatf
    real(kind=8) :: tolgpt
    real(kind=8) :: zebig
    real(kind=8) :: zeneg
    real(kind=8) :: zetamx
    real(kind=8) :: zetsum
    real(kind=8) :: zex
    real(kind=8) :: zvclim
    real(kind=8) :: zvclgi
    real(kind=8) :: zvclgo
    real(kind=8) :: zvcnew
    real(kind=8) :: zx
    real(kind=8) :: zx1
    real(kind=8) :: zx2

    real(kind=8) :: texp

    ! The following are iteration limits:
    !   nlopmx = the maximum number of auto basis switching loops
    !   nplim  = the maximum number of passes
    !   ncylim = the maximum number of cycles
    ! Passes refine estimates of the ionic strength, etc., the
    ! activity of water, and activity coefficients of aqueous species.
    ! Cycles are embedded in passes. They refine estimates of species
    ! concentrations before new estimates of ionic strength, etc.,
    data nlopmx/12/,nplim  /7/,ncylim /15/

    ! The following are local tolerance parameters for the
    ! optimization:
    data tolatf /0.005/,tolgpt /0.1/

    ! The following is the minimum scan increment.
    data dscmin /0.05/

    ! The following are not used by EQ6, but are needed to satisfy
    ! some subroutine calls.
    data iebal  /0/,nredox /0/
    data eh     /0./

    ! The following is set to .false. to turn off iterative improvement
    ! of xbarw in ncmpex.f.
    data qxbarw/.false./

    ! Allocate or reallocate local work arrays as needed.
    if (.not.ALLOCATED(iastak)) then
        ! Local work arrays are not allocated. Zero the saved
        ! array size variables. Note that only one array is tested
        ! to see if it is allocated. It is assumed that all local
        ! work arrays are either allocated or not.
        isv_kmax = 0
        isv_nbtmax = 0
    else
        ! Local work arrays are allocated. Check to see if any of the
        ! array size variables have changed. If so, deallocate
        ! the corresponding local work arrays and zero the corresponding
        ! saved size variables.
        if (kmax .ne. isv_kmax) then
            DEALLOCATE(iastak,jasort,jastak)
            DEALLOCATE(kmvar,knflag)
            DEALLOCATE(daleph,delprc,dijmaj)
            DEALLOCATE(zesort,zeta,zetasq)
            isv_kmax = 0
        end if

        if (nbtmax .ne. isv_nbtmax) then
            DEALLOCATE(efac)
            isv_nbtmax = 0
        end if
    end if

    ! At this point, the saved array size values are zero if the
    ! corresponding arrays need to be allocated.
    if (isv_kmax .eq. 0) then
        ALLOCATE(iastak(kmax),jasort(kmax),jastak(kmax))
        ALLOCATE(kmvar(kmax),knflag(kmax))
        ALLOCATE(daleph(kmax),delprc(kmax),dijmaj(kmax))
        ALLOCATE(zesort(kmax),zeta(kmax),zetasq(kmax))
        isv_kmax = kmax
    end if

    if (isv_nbtmax .eq. 0) then
        ALLOCATE(efac(nbtmax))
        isv_nbtmax = nbtmax
    end if

    ! Zero the contents of the local work arrays.
    do k = 1,kmax
        iastak(k) = 0
        jasort(k) = 0
        jastak(k) = 0
        kmvar(k) = 0
        knflag(k) = 0
    end do

    do k = 1,kmax
        daleph(k) = 0.
        delprc(k) = 0.
        dijmaj(k) = 0.
        zesort(k) = 0.
        zeta(k) = 0.
        zetasq(k) = 0.
    end do

    do n = 1,nbtmax
        efac(n) = 0.
    end do

    if (iodb(3) .ge. 1) then
        write (noutpt,1000)
1000 format(/' Starting Pre-Newton-Raphson Optimization.',/)
    end if

    qbswx = .false.
    eps1hi = 1.0/eps100

    bsigmm = 0.
    bfxi = 0.
    bgamx = 0.

    if (iodb(3) .ge. 3) then
        ! Print the active data file basis set.
        write (noutpt,1010)
1010 format(/16x,'--- Active Data File Basis Set ---',//2x,'krow   Name',30x,'Constraint',/)

        do krow = 1,kdim
            if (krow .le. kbt) then
                ujtp = 'Mass balance, moles'
            else if (krow.ge.km1 .and. krow.le.kxt) then
                ujtp = 'Mass action'
            else
                ujtp = 'Defining equation'
            end if

            j2 = ilnobl(ujtp)

            if (krow .le. kbt) then
                nb = iindx1(krow)
                ns = nbaspd(nb)
            else
                ns = iindx1(krow)
            end if

            ! Calling sequence substitutions:
            !   uspec(ns) for unam48
            call fmspnx(jlen,uspec(ns),uspn56)

            write (noutpt,1020) krow,uspn56,ujtp(1:j2)
1020 format(1x,i4,2x,a32,2x,a)
        end do

        write (noutpt,1030)
1030 format(/1x)
    end if

    ! The label below is a return point after an automatic basis switch.
    ! Here nloop is the loop counter for auto basis switching.
    nloop = -1

200 continue
    nloop = nloop + 1

    if (iodb(3) .ge. 1) then
        write (noutpt,1040) nloop
    end if

1040 format(6x,'nloop= ',i2)

    if (iodb(3) .ge. 3) then
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

            if (kxt .gt. kbt) then
                write (noutpt,1070)
            end if

1070 format(16x,'--- Extended Basis Species ---',//2x,'krow   Species',/)

            do krow = km1,kxt
                ns = iindx1(krow)
                call fmspnx(jlen,uspec(ns),uspn56)
                jlen = min(jlen,32)
                write (noutpt,1080) krow,uspn56(1:jlen)
1080 format(1x,i4,2x,a)
            end do

            write (noutpt,1030)
        end if
    end if

    ! Note: At this point, do not recalculate values of the following:
    !   The SUM(i) m(i) function (sigmam)
    !   The ionic strength (fxi)
    !   The J electrostatic moment function (fje)
    !   The activity coefficients of aqueous species
    !   The mole fraction of water
    !   The activity coefficients of exchanger species
    !   The activity coefficients of solid solution components
    !     present in the equilibrium system.
    ! This subroutine is always called with initial values for all of
    ! these parameters. At the start of an EQ6 run, all these parameters
    ! are initialized by EQ6/exivar.f. After that, new values are
    ! subsequently generated by normal stepping procedures in
    ! EQ6/path.f.
    ! Here npass is the pass counter.
    npass = -1

    ! The label below is a return point for subsequent passes. A pass
    ! is an adjustment for the ionic strength, etc., the activity of
    ! water, and the activity coefficients of the solute species.
210 continue
    npass = npass + 1

    ! Note:
    !   alefnc = aleph convergence function
    !   negafc = the number of successive iterations that the
    !            convergence function alefnc
    !   betfnc = beta convergence function
    alefnc = 0.
    alepoe = 0.
    alepoo = 0.
    negafc = 0

    betfnc = 0.
    btmxoe = 0.
    btmxoo = 0.

    if (iodb(3) .ge. 1) then
        write (noutpt,1110) npass
1110 format(/11x,'npass= ',i2)

        write (noutpt,1120) sigmam,fxi,fje,xbrwlc,xbarwc
1120 format(/13x,'sigmam= ',1pe12.5,/13x,'fxi= ',1pe12.5,/13x,'fje= ',1pe12.5,//13x,'xbrwlc= ',0pf9.5,/13x,'xbarwc= ',1pe12.5,/)
    end if

    ! Here ncycle is the cycle counter.
    ncycle = -1

    qcsigm =abs(bsigmm) .le. tolgpt
    qcfxi = abs(bfxi) .le. tolgpt
    qcgam = bgamx .le. tolgpt

    qtestp = qcsigm .and. qcfxi .and. qcgam

    do krow = 1,kdim
        delprc(krow) = 0.
    end do

    if (iodb(3) .ge. 1) then
        write (noutpt,1190)
    end if

1190 format(/14x,'Beginning cycle corrections',/)

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

    ! Expand the description of the equilibrium system from the
    ! extended basis variables.
    call ncmpex(acflg,act,actlg,cdrs,cegexs,cgexj,conc,conclg,cpgexs,egexjc,egexjf,egexs,eps100,fo2,fo2lg,fsort,fugac,fugalg,iern1,iern2,ietmax,ifrn1,ifrn2,igas,igstak,iindx1,ilrn1,ilrn2,imrn1,imrn2,istack,ixrn1,ixrn2,jcsort,jern1,jern2,jetmax,jflag,jgext,jgsort,jgstak,jjsort,jpflag,jsflag,jsitex,jssort,jstack,kbt,kdim,kelect,kmax,km1,ko2gaq,kwater,kxt,loph,losp,lsort,mgext,mrgexs,mtb,moph,mosp,narn1,narn2,nbasp,nbt,nbtmax,ncmpr,ndrs,ndrsmx,ndrsr,nelect,nern1,nern2,netmax,ngexsa,ngext,ngrn1,ngrn2,ngt,ngtmax,noutpt,no2gaq,nphasx,npt,nptmax,nst,nstmax,nttyo,omega,omeglg,press,qxbarw,q6mode,ugexj,ugexmo,uphase,uspec,xbar,xbarlg,xbarw,xbarwc,xbrwlc,xbrwlg,xlks,zchar,zgexj,zvclg1,zvec1)

    xbarw = xbar(narn1)
    xbrwlg = xbarlg(narn1)

    ! Compute the usual residuals (the alpha vector and beta vectors).
    call betas(acflg,actlg,afcnst,alpha,amtb,bbig,beta,betamx,bneg,cdrs,conc,conclg,coval,csts,eh,ehfac,fo2lg,ibetmx,iebal,iindx1,irdxc3,jcsort,jflag,jsflag,jssort,kbt,kdim,kelect,khydr,kmax,km1,ko2gaq,kwater,kxt,mtb,mosp,narn1,narn2,nbasp,nbtmax,ncosp,ndrs,ndrsmx,ndrsr,nelect,nern1,nern2,nhydr,noutpt,no2gaq,nredox,nst,nstmax,nsts,nstsmx,nstsr,ntfx,ntfxmx,ntfxt,nttyo,omega,qredox,q6mode,tfx,ubbig,ubneg,ubetmx,uspec,uzvec1,weight,xbrwlg,xlke,xlks,zchar)

    do krow = 1,kdim
        bx = beta(krow)

        if (krow .le. kbt) then
            ! Mass balance element.
            if (bx .ge. 0.0) then
                ! Have a non-negative beta(i) value. Set zeta(i) = beta(i).
                zeta(krow) = bx
            else
                ! Have a negative beta(i) value.
                dx = 1.0 + bx

                if (dx .gt. eps100) then
                    ! Have (-1.0 + eps100) < beta(i) < 0.0.
                    ! Set zeta(i) = 1.0 - 1.0/(1.0 + beta(i)).
                    zeta(krow) = 1.0 - (1.0/dx)
                else
                    ! Have beta(i) < (-1.0 + eps100).
                    ! Set zeta(i) = 1.0 - (1.0/eps100).
                    zeta(krow) = 1.0 - eps1hi
                end if
            end if
        else
            ! Mass action element.
            ! zeta(i) = 0.001*beta(i)/cscale(i).
            ns = iindx1(krow)
            zeta(krow) = 0.001*bx/cscale(ns)
        end if
    end do

    ! Calculate the aleph and zeta-squared functions.
    aleph = 0.

    do krow = 1,kdim
        zetasq(krow) = zeta(krow)*zeta(krow)
        aleph = aleph + zetasq(krow)
    end do

    ! Calculate the aleph convergence function.
    alefnc = 0.

    if (mod(ncycle,2) .eq. 0) then
        if (alepoe .ge. smp100) then
            alefnc = (alepoe - aleph)/alepoe
        end if

        alepoe = aleph
    else
        if (alepoo .ge. smp100) then
            alefnc = (alepoo - aleph)/alepoo
        end if

        alepoo = aleph
    end if

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
    if (iodb(3) .ge. 3) then
        write (noutpt,1230)
1230 format(//10x,'--- Pre-Newton-Raphson Optimization Summary ---',//2x,'kcol   Name',32x,'zvclg1      zvec1',/)

        do kcol = 1,kbt
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

        do kcol = km1,kxt
            ns = iindx1(kcol)
            zx1 = zvclg1(kcol)
            zx2 = texp(zx1)

            ! Calling sequence substitutions:
            !   uspec(ns) for unam48
            call fmspnx(jlen,uspec(ns),uspn56)
            write (noutpt,1240) kcol,uspn56,zx1,zx2
        end do

        write (noutpt,1250)
1250 format(/1x,'krow   Name',28x,'Alpha',9x,'Beta',9x,'Zeta',/)

        do krow = 1,kbt
            nb = iindx1(krow)
            ns = nbaspd(nb)

            ! Calling sequence substitutions:
            !   uspec(ns) for unam48
            call fmspnx(jlen,uspec(ns),uspn56)
            write (noutpt,1260) krow,uspn56,alpha(krow),beta(krow),zeta(krow)
1260 format(1x,i3,2x,a28,2x,1pe12.5,2x,e12.5,2x,e12.5)
        end do

        do krow = km1,kxt
            ns = iindx1(krow)

            ! Calling sequence substitutions:
            !   uspec(ns) for unam48
            call fmspnx(jlen,uspec(ns),uspn56)
            write (noutpt,1260) krow,uspn56,alpha(krow),beta(krow),zeta(krow)
        end do

        write (noutpt,1030)
    end if

    if (iodb(3) .ge. 1) then
        ! Calculate some data on the state of the zeta residual
        ! functions.
        zetamx = 0.
        zebig = 0
        uzebig = 'None'
        zeneg = 0
        uzeneg = 'None'

        do kcol = 1,kdim
            zex = zeta(kcol)

            if (zex .gt. zebig) then
                zebig = zex
                uzebig = uzvec1(kcol)
            end if

            if (zex .lt. zeneg) then
                zeneg = zex
                uzeneg = uzvec1(kcol)
            end if

            azex = abs(zex)

            if (azex .gt. zetamx) then
                zetamx = azex
            end if
        end do

        ! Write some data on the state of the aleph and beta residual
        ! functions.
        write (noutpt,1270) aleph,alefnc
1270 format(/19x,'aleph= ',1pe12.5,', alefnc= ',e12.5)

        ! Write some data on the state of the beta residual functions.
        ! Calling sequence substitutions:
        !   jlen1 for jlen
        !   ubbig for unam48
        !   usp156 for uspn56
        call fmspnx(jlen1,ubbig,usp156)

        ! Calling sequence substitutions:
        !   jlen2 for jlen
        !   ubneg for unam48
        !   usp246 for uspn56
        call fmspnx(jlen2,ubneg,usp256)

        write (noutpt,1280) betamx,betfnc,bbig,usp156(1:jlen1),bneg,usp256(1:jlen2)
1280 format(/18x,'betamx= ',1pe12.5,', betfnc= ',e12.5,/18x,'  bbig= ',1pe12.5,', ubbig= ',a,/18x,'  bneg= ',1pe12.5,', ubneg= ',a)

        ! Write some data on the state of the zeta residual functions.
        ! Calling sequence substitutions:
        !   jlen1 for jlen
        !   uzebig for unam48
        !   usp156 for uspn56
        call fmspnx(jlen1,uzebig,usp156)

        ! Calling sequence substitutions:
        !   jlen2 for jlen
        !   uzeneg for unam48
        !   usp256 for uspn56
        call fmspnx(jlen2,uzeneg,usp256)

        write (noutpt,1290) zetamx,zebig,usp156(1:jlen1),zeneg,usp256(1:jlen2)
1290 format(/18x,'zetamx= ',1pe12.5,/18x,' zebig= ',1pe12.5,', uzebig= ',a,/18x,' zeneg= ',1pe12.5,', uzeneg= ',a,/)
    end if

    ! Test the mass balance residuals see if another cycle should be
    ! made before attempting to make an improved estimate of the
    ! ionic strength.
    qtestc = betamx.le.tolxpt .and. bbig.le.tolbig .and.  bneg.ge.tolneg

    if (ncycle.gt.2 .and. alefnc.le.tolatf) then
        negafc = negafc + 1
    else
        negafc = 0
    end if

    qcycnc = negafc .ge. 4

    ! Quit doing cycles if:
    !   1. The cycle convergence criteria are met.
    !   2. The maximum number of cycles have been done.
    !   3. The convergence function alefnc indicates that the cycles
    !        are not converging.
    if (qtestc) then
        go to 400
    end if

    if (ncycle .ge. ncylim) then
        go to 400
    end if

    if (qcycnc) then
        if (iodb(3) .ge. 1) then
            write (noutpt,1292)
1292 format(11x,'The cycles are not converging.',/)
        end if

        go to 400
    end if

    ! Identify the dominant species in each mass balance and
    ! compute the corresponding exponent for a continued
    ! fraction correction.
    ! Presently the species determined as the dominant can only be an
    ! aqueous species or an exchanger species. Minerals and solid
    ! solution component species are ignored.
    call cfracf(cdrs,csts,efac,jcsort,jflag,jssort,kmax,mosp,narn1,narn2,nbasp,nbaspd,nbt,nbtmax,ndrs,ndrsmx,ndrsr,nern1,nern2,nfac,nst,nstmax,nsts,nstsmx,nstsr,q6mode,weight)

    if (iodb(3) .ge. 3) then
        ! Write a table containing the preliminary results.
        kount = 0

        do kcol = 1,kbt
            nb = iindx1(kcol)
            ns = nbasp(nb)
            nsd = nfac(nb)

            if (nsd.ne.0 .and. nsd.ne.ns) then
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
                nsd = nfac(nb)

                if (nsd.ne.0 .and. nsd.ne.ns) then
                    ! Calling sequence substitutions:
                    !   jlen1 for jlen
                    !   uspec(ns) for unam48
                    !   usp156 for uspn56
                    call fmspnx(jlen1,uspec(ns),usp156)

                    ! Calling sequence substitutions:
                    !   jlen2 for jlen
                    !   uspec(nsd) for unam48
                    !   usp256 for uspn56
                    call fmspnx(jlen2,uspec(nsd),usp256)
                    jlen2 = min(jlen2,32)
                    write (noutpt,1320) usp156,usp256(1:jlen2)
1320 format(2x,a32,3x,a)
                end if
            end do

            write (noutpt,1030)
        end if
    end if

    if (qabsw .and. nloop.lt.nlopmx) then
        ! In automatic basis switching mode (iopt(11) .ge. 1), try to
        ! first reduce the magntiude of large positive mass balance
        ! residuals by making one or more basis switches.
        ! Note that this block utilizes the efac array calculated above.
        ! This array is primarily thought of as being connected with
        ! the continued fraction (cycle) algorithm, but it has its use
        ! in automatic basis switching mode.
        call absswa(adhfs,adhfsx,advfs,advfsx,avcnst,axhfs,axhfsx,axlks,axlksx,axvfs,axvfsx,beta,cdrs,cdrsx,cdrtw,cdrw,csts,dhfs,dvfs,efac,eps100,ibswx,iebal,iindx1,iodb,ipch,ipchmx,ipcv,ipcvmx,jcsort,jflag,jsflag,jssort,kbt,kmax,mosp,narn1,narn2,narxmx,narxt,nbasp,nbaspd,nbaspx,nbt,nbtmax,nbw,ncosp,ndrs,ndrsmx,ndrsr,ndrsrx,ndrsx,nelect,nhydr,nodbmx,no2gaq,noutpt,nst,nstmax,nsts,nstsmx,nstsr,nswtch,ntpr,ntprmx,nttyo,presg,press,qbassw,qbswx,q6mode,tempc,uspec,uzvec1,weight,xvfs,xlks,xhfs)

        if (nswtch .gt. 0) then
            do ns = 1,nstmax
                mosp(ns) = 0.
            end do

            av = -99999.
            call initav(losp,nstmax,av)

            ! Reset the ixbasp and cjbasp arrays. The former is a flag
            ! array, each member of which denotes whether the
            ! thermodynamic activity of the corresponding basis species
            ! is defined in terms of molality (= 0) or mole fraction (= 1).
            ! The cjbasp array contains any site stoichiometric factors
            ! associated with the operational basis species.
            call gibasp(cgexj,cjbasp,iern1,ixbasp,jern1,jern2,jetmax,jgext,narn1,narn2,nbasp,nbt,nbtmax,nern1,nern2,netmax,nphasx,nstmax)

            ! Go back for another loop.
            go to 200
        end if
    end if

    ! Sort the zeta-squared terms.
    ! Calling sequence substitutions:
    !   zesort for asort
    !   zetasq for aval
    !   jasort for jsort
    !   iastak for istack
    !   jastak for jstack
    !   kmax for nmax
    !   kdim for nval
    ! Insure that each sort made by calling EQLIBU/qsortw.f is made
    ! starting from scratch. That is, the jasort array left over
    ! from a previous call to EQLIBU/qsortw is not used as a starting
    ! point for the current sort.
    do k = 1,kdim
        jasort(k) = k
    end do

    call qsortw(zesort,zetasq,iastak,jasort,jastak,kmax,noutpt,nttyo,kdim)

    ! Reverse the sorting order, so that the values proceed from
    ! greatest to smallest.
    do k = 1,kdim
        jastak(k) = jasort(kdim + 1 - k)
    end do

    do k = 1,kdim
        krow = jastak(k)
        jasort(k) = krow
        zesort(k) = zetasq(krow)
    end do

    ! Find the set of zeta-squared terms to work on. Consider only
    ! the kkdim largest ones, those that account for 99% of aleph.
    ! If there are more than three, truncate to three.
    ax = 0.
    kkdim = 0
    zetsum = 0.

    do k = 1,kdim
        kkdim = kkdim + 1
        zetsum = zetsum + zesort(k)
        ax = zetsum/aleph

        if (ax .ge. 0.99) then
            go to 270
        end if
    end do

270 continue

    if (kkdim .gt. 3) then
        kkdim = 3
    end if

    ! Compute the vector of partial derivatives d aleph/d log z(j),
    ! where aleph is the sum of squares of the zeta residuals and
    ! log z(j) is j-th primary iteration variable. At present, z(j)
    ! is always the number of moles of the j-th member of the
    ! "extended" basis set. Note that:
    !   d aleph/d log z(j) = 2 * Sum(i) zeta(i)*(d zeta(i)/d log z(j))
    ! where:
    !   d zeta(i)/d log z(j) = d alpha(i)/d log z(j)
    !        for zeta(i) = alpha(i);
    !   d zeta(i)/d log z(j) = d beta(i)/d log z(j)
    !     = [1./amtb(i)] * d alpha(i)/d log z(j)
    !        for zeta(i) = beta(i)
    !   where amtb(i) is treated as a constant, which it really isn't
    !   when i refers to a species like H2O, H+, OH-, O2(g,aq), or e-;
    !   and
    !   d zeta(i)/d log z(j) = [0.001/cscale(i)] * d alpha(i)/d log z(j)
    !        for zeta(i) = [0.001/cscale(i)] * alpha(i)
    ! Note that  d alpha(i)/d log z(j) = Jacob(i,j), where [Jacob] is
    ! the Jacobian matrix used in Newton-Raphson iteration.
    ! Get the Jacobian.
    call matrix(aamatr,al10,bpx,cdrs,cdrtw,cdrw,cjbasp,cnufac,conc,csts,dlogxw,eps100,ibpxmx,iebal,iern1,ietmax,iindx1,ipndx1,irdxc3,ixbasp,ixrn1,jcsort,jern1,jern2,jetmax,jflag,jjsort,jsitex,kbt,kction,kdim,kelect,khydr,kmax,kmt,km1,ko2gaq,kwater,kxt,kx1,mosp,narn1,narn2,nbasp,nbt,nbtmax,nbw,ncosp,ndrs,ndrsmx,ndrsr,nelect,nern1,nern2,netmax,noutpt,no2gaq,nphasx,nredox,nst,nstmax,nsts,nstsmx,nstsr,ntfx,ntfxmx,ntfxt,nttyo,nxtmax,omega,qredox,q6mode,tfx,ugexmo,uspec,weight,xbar,xbarw,xbarwc,zchar)

    ! Compute the d aleph/d log z(j) vector.
    do kcol = 1,kdim
        daleph(kcol) = 0.
    end do

    do krow = 1,kbt
        nb = iindx1(krow)
        bx = beta(krow)

        if (bx .ge. 0.0) then
            ! Have a non-negative beta(i) value: zeta(i) = beta(i).
            if (amtb(nb) .gt. eps100) then
                cx = 2.*zeta(krow)/amtb(nb)
            else
                cx = 2.*zeta(krow)*eps1hi
            end if
        else
            ! Have a negative beta(i) value.
            dx = 1.0 + bx

            if (dx .gt. eps100) then
                ! Have (-1.0 + eps100) < beta(i) < 0.0:
                ! zeta(i) = 1.0 - 1.0/(1.0 + beta(i)).
                cx = 2.*zeta(krow)**3/amtb(nb)
            else
                ! Have beta(i) < (-1.0 + eps100):
                ! zeta(i) = 1.0 - (1.0/eps100).
                if (amtb(nb) .gt. eps100) then
                    cx = 2.*eps1hi**3/amtb(nb)
                else
                    cx = 2.*eps1hi**4
                end if
            end if
        end if

        do kcol = 1,kdim
            daleph(kcol) = daleph(kcol) + cx*aamatr(krow,kcol)
        end do
    end do

    do krow = km1,kxt
        ns = iindx1(krow)
        cx = 0.001/cscale(ns)

        do kcol = 1,kdim
            daleph(kcol) = daleph(kcol) + cx*aamatr(krow,kcol)
        end do
    end do

    if (iodb(3) .ge. 4) then
        write (noutpt,1330)
1330 format(/7x,'--- Partial Derivatives of Aleph ---','  krow    Variable',25x,'daleph(krow)',/)

        do krow = 1,kdim
            ! Calling sequence substitutions:
            !   uzvec1(krow) for unam48
            call fmspnx(jlen,uzvec1(krow),uspn56)

            write (noutpt,1340) krow,uspn56,daleph(krow)
1340 format(2x,i3,3x,a32,3x,1pe12.5)
        end do

        write (noutpt,1030)
    end if

    ! Now find the primary iteration variable that is the principal
    ! determining variable for each of the most significant zeta-
    ! squared variables.
    do k = 1,kkdim
        kmvar(k) = 0
    end do

    do k = 1,kkdim
        krow = jasort(k)

        if (krow .le. kbt) then
            nb = iindx1(krow)
            bx = beta(krow)

            if (bx .ge. 0.0) then
                ! Have a non-negative beta(i) value: zeta(i) = beta(i).
                if (amtb(nb) .gt. eps100) then
                    cx = 2.*zeta(krow)/amtb(nb)
                else
                    cx = 2.*zeta(krow)*eps1hi
                end if
            else
                ! Have a negative beta(i) value.
                dx = 1.0 + bx

                if (dx .gt. eps100) then
                    ! Have (-1.0 + eps100) < beta(i) < 0.0:
                    ! zeta(i) = 1.0 - 1.0/(1.0 + beta(i)).
                    if (amtb(nb) .gt. eps100) then
                        cx = 2.*zeta(krow)**3/amtb(nb)
                    else
                        cx = 2.*zeta(krow)**4
                    end if
                else
                    ! Have beta(i) < (-1.0 + eps100):
                    ! zeta(i) = 1.0 - (1.0/eps1000.
                    if (amtb(nb) .gt. eps100) then
                        cx = 2.*eps1hi**3/amtb(nb)
                    else
                        cx = 2.*eps1hi**4
                    end if
                end if
            end if
        else
            ns = iindx1(krow)
            cx = 0.002*zeta(krow)/cscale(ns)
        end if

        kbig = 0
        abig = 0.

        do kcol = 1,kdim
            dij = cx*aamatr(krow,kcol)
            adij = abs(dij)

            if (adij .gt. abig) then
                if (krow.gt.kbt .or. kcol.ne.kwater .or. krow.eq.kwater)        then
                    ! For a mass balance row other than that for H2O, disallow
                    ! H2O as the principal determining species.
                    kbig = kcol
                    abig = adij
                    dijmaj(krow) = dij
                end if
            end if
        end do

        kmvar(k) = kbig
    end do

    ! Special restrictions apply to H2O, H+, O2(g,aq), and e-.
    do k = 1,kkdim
        krow = jasort(k)
        kcol = kmvar(k)

        if (krow.eq.kwater .or. krow.eq.khydr .or.    krow.eq.ko2gaq .or. krow.eq.kelect) then
            ! Only H2O may be adjusted for the H2O row, only H+ for the H+
            ! row, only O2(g,aq) for the O2 row, and e- for the e- row.
            if (kcol .ne. krow) then
                kmvar(k) = krow
            end if

            go to 290
        end if

        if (krow .le. kbt) then
            kcol = kmvar(k)

            if (kcol.eq.kwater .or. kcol.eq.khydr .or.      kcol.eq.ko2gaq .or. kcol.eq.kelect) then
                ! H2O may be adjusted only for the H2O row, among mass
                ! balance rows, H+ for the H+ row, O2(g,aq) for the O2 row,
                ! and e- for the e- row.
                if (kcol .ne. krow) then
                    kmvar(k) = krow
                end if
            end if
        end if

290 continue
    end do

    ! The following is a return point if no cycle corrections were
    ! generated and an alternate set of variables to correct is to be
    ! tried.
350 continue

    if (iodb(3) .ge. 2) then
        write (noutpt,1400)
1400 format(/11x,'--- Principal Dependencies ---',//4x,'Row Variable (i)',19x,'Column Variable (j)',/)

        do k = 1, kkdim
            krow = jasort(k)

            if (krow .le. kbt) then
                nb = iindx1(krow)
                ns = nbaspd(nb)
            else
                ns = iindx1(krow)
            end if

            kcol = kmvar(k)

            ! Calling sequence substitutions:
            !   jlen1 for jlen
            !   uspec(ns) for unam48
            !   usp156 for uspn56
            call fmspnx(jlen1,uspec(ns),usp156)

            if (kcol .gt. 0) then
                ! Calling sequence substitutions:
                !   jlen2 for jlen
                !   uzvec1(kcol) for unam48
                !   usp256 for uspn56
                call fmspnx(jlen2,uzvec1(kcol),usp256)
            else
                usp256 = 'None'
                jlen2 = 4
            end if

            write (noutpt,1410) usp156,usp256(1:jlen2)
1410 format(2x,a32,3x,a)
        end do
    end if

    if (iodb(3) .ge. 2) then
        write (noutpt,1030)
    end if

    do k = 1,kkdim
        if (kmvar(k) .eq. 0) then
            krow = jasort(k)

            if (krow .le. kbt) then
                nb = iindx1(krow)
                ns = nbaspd(nb)
            else if (krow.ge.km1 .and. krow.le.kxt) then
                ns = iindx1(krow)
            end if

            kcol = krow
            kmvar(k) = kcol

            ! Calling sequence substitutions:
            !   jlen1 for jlen
            !   uspec(ns) for unam48
            !   usp156 for uspn56
            call fmspnx(jlen1,uspec(ns),usp156)

            ! Calling sequence substitutions:
            !   jlen2 for jlen
            !   uzvec1(kcol) for unam48
            !   usp256 for uspn56
            call fmspnx(jlen2,uzvec1(kcol),usp256)

            write (noutpt,1440) usp156(1:jlen1),usp256(1:jlen2)
            write (nttyo,1440) usp156(1:jlen1),usp256(1:jlen2)
1440 format(/" * Warning- (EQ6/optmzr) Couldn't find a variable",' to adjust to reduce',/7x,'the aleph residual for ',a,'. Will use the variable associated',/7x,'with ',a,'.')
        end if
    end do

    ! Find and flag all instances in which a principal determining
    ! variable duplicates a prior entry in the list. Set knflag(k) = 0
    ! for such duplications, otherwise set knflag(k) = 1.
    do k = 1,kkdim
        knflag(k) = 1
    end do

    do k = 1,kkdim - 1
        do kk = k + 1,kkdim
            if (kmvar(kk) .eq. kmvar(k)) then
                knflag(kk) = 0
            end if
        end do
    end do

    ! Find the number of variables to be corrected in the current
    ! cycle.
    jdim = 0

    do k = 1,kkdim
        if (knflag(k) .eq. 1) then
            jdim = jdim + 1
        end if
    end do

    if (iodb(3) .ge. 3) then
        ! Write a table containing data on the structure of the
        ! corrections to be made for the current step.
        write (noutpt,1530)
1530 format(/6x,'--- Structure for Current Cycle Corrections',' ---',/)

        write (noutpt,1540)
1540 format(4x,'Row Variable (i)',19x,'zeta(i)**2',/)

        do k = 1,kkdim
            krow = jasort(k)

            if (krow .le. kbt) then
                nb = iindx1(krow)
                ns = nbaspd(nb)
            else
                ns = iindx1(krow)
            end if

            ! Calling sequence substitutions:
            !   jlen1 for jlen
            !   uspec(ns) for unam48
            !   usp156 for uspn56
            call fmspnx(jlen1,uspec(ns),usp156)

            write (noutpt,1550) usp156,zesort(k)
1550 format(2x,a32,3x,1pe12.5)
        end do

        write (noutpt,1560)
1560 format(//4x,'Column Variable (j)',9x,'d zeta(i)**2/d log z(j)',/)

        do k = 1,kkdim
            ! Note: write all duplications, so don't test to see if
            ! knflag(k) = 1.
            krow = jasort(k)
            kcol = kmvar(k)

            ! Calling sequence substitutions:
            !   jlen2 for jlen
            !   uzvec1(kcol) for unam48
            !   usp256 for uspn56
            call fmspnx(jlen,uzvec1(kcol),usp256)

            write (noutpt,1550) usp256,dijmaj(krow)
        end do

        write (noutpt,1030)
    end if

    if (jdim .gt. 1) then
        ! If the magnitude of aleph is quite large, correct only one
        ! variable per cycle.
        if (aleph .gt. 1.e+10) then
            if (iodb(3) .ge. 3) then
                write (noutpt,1570)
1570 format(' Because of the large magnitude of aleph, only',' one variable will be',/' corrected in this cycle.',/)
            end if

            jdim = 1

            do k = 2,kkdim
                knflag(k) = 0
            end do
        end if
    end if

    ! Execute the cycle algorithm.
    ! Testing has shown that the alpha residual for the O2(g,aq) mass
    ! balance taken as a function of the log fO2 variable can have
    ! a shape like the following:
    !              |    x                              |
    !              |    x                              |
    !              |    x                              |
    !              |    x                              |
    !              |    x                              |
    !        ^     |    x                              |
    !        |     |     x                             |
    !              |      x                            |
    !      alpha   |       xxxxxxxxxxxxxx              |
    !            0 |-----------------------x-----------|
    !              |                          x        |
    !              |                           x       |
    !              |                            x      |
    !              |                            x      |
    !              |                            x      |
    !              |                            x      |
    !              |                            x      |
    !              |                            x      |
    !              |                            x      |
    !                           log fO2 ->
    ! For this variable, the zeta residual is the alpha residual.
    ! The function is nearly flat over about an 8-10 unit range of
    ! log fO2, and very steep on the sides. It looks like a step
    ! function rotated 90 degrees. The corresponding plot of aleph
    ! versus log fO2 looks like a square well.
    ! On the limbs, the slopes are so steep that attempting to
    ! extrapolate to a zero residual using first-order partial
    ! derivatives produces corrections that are practically
    ! infinitesimal. Thus, any method that generates the values of
    ! corrections from first-order partial derivatives is pretty much
    ! guaranteed to fail. This rules out using whole classes of
    ! methods.
    ! Within the "flat" region, the Newton-Raphson method is
    ! pretty much guaranteed to succeed. In the present subroutine
    ! what we need is something that will get us in or near the flat
    ! region when we are literally "out on a limb."
    ! One possibility is to use a scanning procedure. Another is to
    ! use a non-scanning method that does not use derivatives. An
    ! example of the latter would be the cycle algorithm that
    ! is employed for pre-Newton-Raphson optimiztion by EQ3NR. However,
    ! that algorithm does not extend well to EQ6 problems. In EQ6, mass
    ! balances must be considered for basis species such as H2O(l), H+,
    ! and O2(g,aq) (or alternatives such as OH- and e-). These mass
    ! balances are distinct in that some species can make negative
    ! contributions. Also, EQ6 problems require consideration of some
    ! heterogeneous mass action equations when the active basis set is
    ! "extended" to include non-aqueous species (e.g., minerals).
    ! In the present version of this subroutine, we will use a simple
    ! scannning approach. Although the first-order partial derivatives
    ! appear to be not useful in estimating the actual magnitude of
    ! needed corrections, they do appear to identify the primary
    ! iteration variables that are most in need of correction. They
    ! also appear to define the direction in which one should be
    ! scanning.
    kcorr = 0

    do krow = 1,kdim
        delvec(krow) = 0.
    end do

    j = 0

    do k = 1,kkdim
        if (knflag(k) .eq. 1) then
            j = j + 1
            krscan = jasort(k)
            kcscan = kmvar(k)

            if (iodb(3) .ge. 2) then
                ! Calling sequence substitutions:
                !   uzvec1(kcscan) for unam48
                call fmspnx(jlen,uzvec1(kcscan),uspn56)

                write (noutpt,1610) uspn56(1:jlen)
1610 format(/'   --- Scanning on ',a,' ---','  zvclg1(kcscan)  zetasq(krscan)',3x,'aleph',10x,'dscan',/)

                dscan = 0.
                write (noutpt,1620) zvclg1(kcscan),zetasq(krscan),aleph,dscan
1620 format(5x,f9.4,3x,1pe12.5,3x,e12.5,3x,e12.5)
            end if

            ! Get the direction in which to start scanning. Look at the
            ! sign of the appropriate partial derivative to get this.
            ! The other direction will be scanned if no improvement can
            ! be made in the initial direction.
            ! Note: could try dijmaj(krscan) in place of dalpeh(kcscan)
            ! below.
            if (daleph(kcscan) .gt. 0.) then
                kdirec = -1
            else
                kdirec = 1
            end if

            if (kcscan .eq. krscan) then
                if (zeta(krscan) .gt. 0.) then
                    kdirec = -1
                else
                    kdirec = 1
                end if
            end if

            qscan2 = .false.

            ! The label below is a return point if the scan is being
            ! continued in the opposite direction.
300 continue

            ! Set the initial value of the scan increment.
            if (kcscan .eq. kwater) then
                dscan = 0.0625

                if (abs(zeta(krscan)) .gt. 10.) then
                    dscan = 0.125
                end if

                if (abs(zeta(krscan)) .gt. 1.e+2) then
                    dscan = 0.25
                end if

                if (abs(zeta(krscan)) .gt. 1.e+4) then
                    dscan = 0.5
                end if

                if (abs(zeta(krscan)) .gt. 1.e+8) then
                    dscan = 1.0
                end if
            else if (kcscan .eq. khydr) then
                dscan = 0.0625

                if (abs(zeta(krscan)) .gt. 10.) then
                    dscan = 0.125
                end if

                if (abs(zeta(krscan)) .gt. 1.e+2) then
                    dscan = 0.25
                end if

                if (abs(zeta(krscan)) .gt. 1.e+4) then
                    dscan = 0.5
                end if

                if (abs(zeta(krscan)) .gt. 1.e+8) then
                    dscan = 1.0
                end if

                if (abs(zeta(krscan)) .gt. 1.e+16) then
                    dscan = 2.0
                end if
            else if (kcscan .eq. ko2gaq) then
                dscan = 0.25

                if (abs(zeta(krscan)) .gt. 10.) then
                    dscan = 0.5
                end if

                if (abs(zeta(krscan)) .gt. 1.e+2) then
                    dscan = 1.0
                end if

                if (abs(zeta(krscan)) .gt. 1.e+4) then
                    dscan = 2.0
                end if

                if (abs(zeta(krscan)) .gt. 1.e+8) then
                    dscan = 4.0
                end if

                if (abs(zeta(krscan)) .gt. 1.e+12) then
                    dscan = 6.0
                end if

                if (abs(zeta(krscan)) .gt. 1.e+16) then
                    dscan = 8.0
                end if

                if (abs(zeta(krscan)) .gt. 1.e+24) then
                    dscan = 12.0
                end if

                if (abs(zeta(krscan)) .gt. 1.e+32) then
                    dscan = 16.0
                end if
            else if (kcscan .eq. kelect) then
                dscan = 0.125

                if (abs(zeta(krscan)) .gt. 10.) then
                    dscan = 0.5
                end if

                if (abs(zeta(krscan)) .gt. 1.e+2) then
                    dscan = 1.0
                end if

                if (abs(zeta(krscan)) .gt. 1.e+4) then
                    dscan = 2.0
                end if

                if (abs(zeta(krscan)) .gt. 1.e+8) then
                    dscan = 4.0
                end if

                if (abs(zeta(krscan)) .gt. 1.e+12) then
                    dscan = 6.0
                end if

                if (abs(zeta(krscan)) .gt. 1.e+16) then
                    dscan = 8.0
                end if

                if (abs(zeta(krscan)) .gt. 1.e+24) then
                    dscan = 12.0
                end if

                if (abs(zeta(krscan)) .gt. 1.e+32) then
                    dscan = 16.0
                end if
            else
                dscan = 0.125

                if (abs(zeta(krscan)) .gt. 10.) then
                    dscan = 0.5
                end if

                if (abs(zeta(krscan)) .gt. 1.e+2) then
                    dscan = 1.0
                end if

                if (abs(zeta(krscan)) .gt. 1.e+4) then
                    dscan = 2.0
                end if

                if (abs(zeta(krscan)) .gt. 1.e+8) then
                    dscan = 4.0
                end if

                if (abs(zeta(krscan)) .gt. 1.e+12) then
                    dscan = 6.0
                end if

                if (abs(zeta(krscan)) .gt. 1.e+16) then
                    dscan = 8.0
                end if

                if (abs(zeta(krscan)) .gt. 1.e+24) then
                    dscan = 12.0
                end if

                if (abs(zeta(krscan)) .gt. 1.e+32) then
                    dscan = 16.0
                end if
            end if

            icorr = 0
            zvclgi = zvclg1(kcscan)

            if (kdirec .gt. 0) then
                zvclim = 99999.
            else
                zvclim = -99999.
            end if

            nscan = 0
            qstops = .false.

            zvclgo = zvclg1(kcscan)
            zetsqo = zetasq(krscan)
            zetsqm = zetsqo
            alephm = aleph
            nbigger = 0

            ! The following is a return point for continuing the scan.
410 continue
            nscan = nscan + 1

            ! Increment the variable being scanned.
            zvclg1(kcscan) = zvclg1(kcscan) + kdirec*dscan

            ! Apply a total change limit for a given cycle.
            dzxclm = 40.0

            if (kcscan .eq. kwater) then
                dzxclm = 2.0
            end if

            if (kcscan .eq. khydr) then
                dzxclm = 4.0
            end if

            if (kcscan.ge.km1 .and.kcscan.le.kxt) then
                dzxclm = 0.25
            end if

            qsclim = .false.
            dzx = zvclg1(kcscan) - zvclgi

            if (abs(dzx) .gt. dzxclm) then
                qsclim = .true.

                if (kdirec .gt. 0) then
                    zvclg1(kcscan) = zvclgi + dzxclm
                else
                    zvclg1(kcscan) = zvclgi - dzxclm
                end if
            end if

            ! Expand the description of the equilibrium system using the
            ! incremented variable.
            call ncmpex(acflg,act,actlg,cdrs,cegexs,cgexj,conc,conclg,cpgexs,egexjc,egexjf,egexs,eps100,fo2,fo2lg,fsort,fugac,fugalg,iern1,iern2,ietmax,ifrn1,ifrn2,igas,igstak,iindx1,ilrn1,ilrn2,imrn1,imrn2,istack,ixrn1,ixrn2,jcsort,jern1,jern2,jetmax,jflag,jgext,jgsort,jgstak,jjsort,jpflag,jsflag,jsitex,jssort,jstack,kbt,kdim,kelect,kmax,km1,ko2gaq,kwater,kxt,loph,losp,lsort,mgext,mrgexs,mtb,moph,mosp,narn1,narn2,nbasp,nbt,nbtmax,ncmpr,ndrs,ndrsmx,ndrsr,nelect,nern1,nern2,netmax,ngexsa,ngext,ngrn1,ngrn2,ngt,ngtmax,noutpt,no2gaq,nphasx,npt,nptmax,nst,nstmax,nttyo,omega,omeglg,press,qxbarw,q6mode,ugexj,ugexmo,uphase,uspec,xbar,xbarlg,xbarw,xbarwc,xbrwlc,xbrwlg,xlks,zchar,zgexj,zvclg1,zvec1)

            xbarw = xbar(narn1)
            xbrwlg = xbarlg(narn1)

            ! Recompute the residuals.
            call betas(acflg,actlg,afcnst,alpha,amtb,bbig,beta,betamx,bneg,cdrs,conc,conclg,coval,csts,eh,ehfac,fo2lg,ibetmx,iebal,iindx1,irdxc3,jcsort,jflag,jsflag,jssort,kbt,kdim,kelect,khydr,kmax,km1,ko2gaq,kwater,kxt,mtb,mosp,narn1,narn2,nbasp,nbtmax,ncosp,ndrs,ndrsmx,ndrsr,nelect,nern1,nern2,nhydr,noutpt,no2gaq,nredox,nst,nstmax,nsts,nstsmx,nstsr,ntfx,ntfxmx,ntfxt,nttyo,omega,qredox,q6mode,tfx,ubbig,ubneg,ubetmx,uspec,uzvec1,weight,xbrwlg,xlke,xlks,zchar)

            ! Recompute the zeta residuals and the aleph residual.
            do krow = 1,kdim
                bx = beta(krow)

                if (krow .le. kbt) then
                    ! Mass balance element.
                    if (bx .ge. 0.0) then
                        ! Have a non-negative beta(i) value. Set zeta(i) = beta(i).
                        zeta(krow) = bx
                    else
                        ! Have a negative beta(i) value.
                        dx = 1.0 + bx

                        if (dx .gt. eps100) then
                            ! Have (-1.0 + eps100) < beta(i) < 0.0.
                            ! Set zeta(i) = 1.0 - 1.0/(1.0 + beta(i)).
                            zeta(krow) = 1.0 - (1.0/dx)
                        else
                            ! Have beta(i) < (-1.0 + eps100).
                            ! Set zeta(i) = 1.0 - (1.0/eps100).
                            zeta(krow) = 1.0 - eps1hi
                        end if
                    end if
                else
                    ! Mass action element.
                    ! zeta(i) = 0.001*beta(i)/cscale(i).
                    ns = iindx1(krow)
                    zeta(krow) = 0.001*bx/cscale(ns)
                end if
            end do

            ! Calculate the aleph and zeta-squared functions.
            aleph = 0.

            do krow = 1,kdim
                zetasq(krow) = zeta(krow)*zeta(krow)
                aleph = aleph + zetasq(krow)
            end do

            if (iodb(3) .ge. 2) then
                write (noutpt,1620) zvclg1(kcscan),zetasq(krscan),aleph,dscan

                if (nscan .le. 0) then
                    ! The best point for this scan within the limits of the
                    ! minimum non-zero scan increment was the initial point.
                    write (noutpt,1630)
1630 format(/5x,'No correction was made.')
                end if
            end if

            ! If the scan has stepped back to a previous point because
            ! it gave a smaller residual, exit the scan now that all
            ! residual functions have been recalculated.
            if (qstops) then
                go to 420
            end if

            ! Now check the results at the current point of the scan.
            ! Is the new point better or not? If either the current
            ! zeta-squared function or the aleph function has increased,
            ! then the scan has over-corrected.
            if (aleph .gt. alephm) then
                nbigger = nbigger + 1
            end if

            if (zetasq(krscan).gt.zetsqm .or.      (aleph.gt.alephm .and. nbigger.ge.2)) then
                ! Have over-corrected. Back up.
                nscan = nscan - 1
                nbigger = nbigger - 1
                zvclim = zvclg1(kcscan)
                zvclg1(kcscan) = zvclgo

                ! Note: not everything has been reset here. It is still
                ! necessary to go back and re-expand the description of the
                ! equilibrium system.
                if (dscan .gt. dscmin) then
                    ! Reduce the scan increment.
                    dscan = 0.5*dscan

                    ! Reduce the scan increment again if the new value of the
                    ! iteration variable would match the previously tried
                    ! extreme value.
                    zvcnew = zvclg1(kcscan) + kdirec*dscan

                    if (kdirec .gt. 0) then
                        if (zvcnew .gt. (zvclim - eps100)) then
                            dscan = 0.5*dscan
                        end if
                    else
                        if (zvcnew .lt. (zvclim + eps100)) then
                            dscan = 0.5*dscan
                        end if
                    end if

                    if (dscan .lt. dscmin) then
                        dscan = dscmin
                    end if

                    ! Go back and try the new, reduced scan increment.
                else
                    ! The previous value for this scan was better than any
                    ! subsequent adjusted value within the limitations imposed
                    ! by the minimum scan increment.
                    dscan = 0.
                    nscan = nscan  -1

                    ! Set up to stop the current scan.
                    qstops = .true.
                end if

                ! Go back and re-expand the system and recalculate all
                ! residual functions. If qstops is .true., the scan increment
                ! is zero and the only purpose in going back is to ensure
                ! that all the residual functions are consistent with the
                ! last point of the scan before continuing with the next
                ! optimization step (e.g., next scan, cycle, pass), if any.
                qsclim = .false.
                go to 410
            end if

            ! At this point, a non-zero correction has resulted in an
            ! improvement, at least in the current zeta-squared function.
            icorr = icorr + 1
            zvclgo = zvclg1(kcscan)
            zetsqo = zetasq(krscan)
            zetsqm = min(zetsqm,zetsqo)

            if (aleph.lt. alephm) then
                ! The aleph function was also improved.
                alephm = aleph
                nbigger = 0
            end if

            if (.not.qsclim .and. nscan .lt. 10) then
                zx = zvclg1(kcscan) - zvclgi

                if (abs(zx) .lt. dzxclm) then
                    ! Continue the current scan, using the same scan increment.
                    go to 410
                end if
            end if

420 continue
            if (iodb(3) .ge. 2) then
                write (noutpt,1030)
            end if

            ! At this point, a scan in one direction is complete.
            if (icorr .eq. 0) then
                ! No correction was made by the scan in the current
                ! direction.
                if (.not.qscan2) then
                    ! Go back and scan in the other direction.
                    qscan2 = .true.
                    kdirec = -kdirec

                    if (iodb(3) .ge. 2) then
                        write (noutpt,1640)
1640 format(5x,'Scanning in the opposite direction',/)
                    end if

                    go to 300
                end if
            end if

            ! End corrections for the current variable on this cycle.
            ! Calculate the total correction for this scan.
            delvec(kcscan) = zvclg1(kcscan) - zvclgi

            ! Count the number of variables which were actually corrected.
            if (icorr .gt. 0) then
                kcorr = kcorr + 1
            end if
        end if
    end do

    if (kcorr .le. 0) then
        if (iodb(3) .ge. 2) then
            write (noutpt,1650)
1650 format(11x,'No corrections were generated.')
        end if

        ! See if an alternate set of variables to correct should be
        ! tried. For example, if H+ was to be corrected to reduce the
        ! residual for Al+++, try correcting Al+++.
        nchange = 0

        do k = 1,kkdim
            krow = jasort(k)
            kcol = kmvar(k)

            if (kcol .ne. krow) then
                kmvar(k) = krow
                nchange = nchange + 1
            end if
        end do

        if (nchange .gt. 0) then
            if (iodb(3) .ge. 2) then
                write (noutpt,1660)
1660 format(11x,'Changing the set of variables to be',' corrected',/11x,'and trying again.',/)
            end if

            go to 350
        end if

        write (noutpt,1662)
1662 format(1x)
    end if

    ! Check for an oscillating scan.
    qscosc = .false.

    do kcol = 1,kdim
        if (abs(delvec(kcol) + delprc(kcol)) .gt. eps100) then
            go to 370
        end if
    end do

    ! The scan is oscillating. See if an alternate set of variables
    ! to correct should be tried. For example, if H+ was to be
    ! corrected to reduce the residual for Al+++, try correcting Al+++.
    if (iodb(3) .ge. 2) then
        write (noutpt,1670)
1670 format(/11x,'The corrections are oscillating.')
    end if

    nchange = 0

    do k = 1,kkdim
        krow = jasort(k)
        kcol = kmvar(k)

        if (kcol .ne. krow) then
            kmvar(k) = krow
            nchange = nchange + 1
        end if
    end do

    if (nchange .gt. 0) then
        if (iodb(3) .ge. 2) then
            write (noutpt,1660)
        end if

        go to 350
    end if

    write (noutpt,1662)

    qscosc = .true.

370 continue

    do kcol = 1,kdim
        delprc(kcol) = delvec(kcol)
    end do

    if (iodb(3).ge.1 .and. kcorr.gt.0) then
        ! Summarize the scan corrections for the current cycle.
        write (noutpt,1680)
1680 format(/16x,'--- Cycle Corrections ---',//2x,'krow   Name',32x,'zvclgi',6x,'zvclg1',6x,'delvec',/)

        do k = 1,kkdim
            if (knflag(k) .eq. 1) then
                krow = jasort(k)
                kcol = kmvar(k)

                if (delvec(kcol) .ne. 0.) then
                    ! Calling sequence substitutions:
                    !   uzvec1(kcol) for unam48
                    call fmspnx(jlen,uzvec1(kcol),uspn56)
                    zvclgi = zvclg1(kcol) - delvec(kcol)
                    write (noutpt,1690) krow,uspn56,zvclgi,zvclg1(kcol),delvec(kcol)
1690 format(1x,i4,2x,a32,2x,f10.4,2x,f10.4,2x,f10.4)
                end if
            end if
        end do

        write (noutpt,1030)
    end if

    if (kcorr.le.0 .or. qscosc .or. qsclim) then
        go to 400
    end if

    ! Go back for another cycle.
    go to 220

    ! The cycles for the current pass have been completed. Test to
    ! see if another pass should be made.
400 continue
    if (iodb(3) .ge. 1) then
        write (noutpt,1700) npass,ncycle
1700 format(13x,'Completed pass ',i3,' in ',i3,' cycles.')
    end if

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
        write (noutpt,1720) bsigmm
1720 format(/13x,'bsigmm= ',1pe12.5)

        write (noutpt,1730) bfxi
1730 format(13x,'bfxi= ',1pe12.5)

        write (noutpt,1740) bfje
1740 format(13x,'bfje= ',1pe12.5)

        j2 = ilnobl(ubgamx(1:24))
        write (noutpt,1750) bgamx,ubgamx(1:j2)
1750 format(13x,'bgamx= ',1pe12.5,', ubgamx= ',a)
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
            if (iodb(3) .ge. 1) then
                if (ker .lt. 2) then
                    write (noutpt,1800)
1800 format(/'   Done. Optimization ended outside requested',' limits.',/)
                else
                    write (noutpt,1810)
1810 format(/'   Done. Optimization ended outside allowable',' limits.',/)
                end if
            end if

            go to 999
        end if

        ! Do another pass.
        go to 210
    end if

    ! Are cycle criteria satisfied?
    if (qtestc) then
        ! Yes, optimization succeeded.
        if (iodb(3) .ge. 1) then
            write (noutpt,1830)
1830 format(/'   Done. Optimization ended within requested',' limits.',/)
        end if

        go to 999
    else if (ncycle.eq.0 .and. kcorr .le. 0) then
        ! No cycle corrections were generated after starting a new
        ! pass. Quit.
        if (iodb(3) .ge. 1) then
            if (ker .lt. 2) then
                write (noutpt,1800)
            else
                write (noutpt,1810)
            end if
        end if

        go to 999
    else if (npass .le. 2) then
        ! The pass convergence criteria are satisfied, but the cycle
        ! convergence criteria are not.
        ! Try another pass.
        go to 210
    else
        ! Quit. Optimization ended outside requested limits
        ! because cycle requirements were not met.
        if (iodb(3) .ge. 1) then
            if (ker .lt. 2) then
                write (noutpt,1800)
            else
                write (noutpt,1810)
            end if
        end if

        go to 999
    end if

999 continue
end subroutine optmzr