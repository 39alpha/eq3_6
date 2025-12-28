subroutine eqcalc(aamatr,abar,acflg,acflgo,act,actlg,adh,adhh,adhv,afcnst,alpha,al10,amtb,aphi,apx,avcnst,azero,a3bar,a3bars,bacfmx,bbig,bdh,bdhh,bdhv,bdot,bdoth,bdotv,beta,betamx,betao,bfje,bfxi,bgamx,bneg,bpx,bsigmm,cco2,cegexs,cess,cdrs,cdrsd,cdrsx,cdrtw,cdrw,cjbasp,cnufac,conc,conclg,cpgexs,cscale,csts,delmax,delvco,delvec,dlogxw,egexjc,egexjf,egexs,eh,ehfac,eps100,farad,fje,fjeo,fo2,fo2lg,fsort,fugac,fugalg,fxi,fxio,gmmatr,iapxt,ibetmx,ibpxt,ibswx,idelmx,ielam,ier,iern1,iern2,ifcphi1,ifcphi2,ifnnn,ifn2n,ifpsi1,ifpsi2,ifrn1,ifrn2,ifzeta,igas,igstak,iindx1,ilcphi1,ilcphi2,ilnnn,iln2n,ilpsi1,ilpsi2,ilrn1,ilrn2,ilzeta,imrn1,imrn2,insgf,iodb,iopg,iopt,ipch,ipcv,ipivot,ipndx1,istack,iter,itermx,ixbasp,ixrn1,ixrn2,izmax,jcsort,jflag,jgsort,jgstak,jjsort,jpflag,jsflag,jsitex,jsol,jssort,jstack,kbt,kction,kdim,kelect,khydr,khydx,km1,kmt,ko2gaq,krdxsp,kwater,kx1,kxt,loph,losp,lsort,moph,mosp,mrgexs,mtb,narn1,narn2,narxt,nat,nbasp,nbaspd,nbaspx,nbt,nbtd,nbw,nchlor,ncmpr,nct,ndrs,ndrsd,ndrsx,ndrsr,ndrsrd,ndrsrx,nelect,nern1,nern2,ness,nessr,net,nfac,nfrn1,nfrn2,ngrn1,ngrn2,ngt,nhydr,nhydx,nlrn1,nlrn2,nlt,nmrn1,nmrn2,nmt,noutpt,no2gaq,nphasx,npt,nrdxsp,nst,nsts,nstsr,ntpr,nttyo,nxrn1,nxrn2,nxt,omega,omeglg,press,qbassw,qhawep,qoptmz,qpit75,qredox,q6mode,rhsvec,screwd,sigmam,sigmmo,smp100,tempc,tempk,tolbt,toldl,ubacmx,ubgamx,ulbeta,uldel,uphase,uspec,uzvec1,weight,wfac,xbar,xbarlg,xbarw,xbarwc,xbrwlc,xbrwlg,zchar,zchcu6,zchsq2,zvclg1,zvec1)
    !! This subroutine attempts to calculate the chemical equilibrium
    !! state of a system for a pre-defined phase assemblage. Unlike in
    !! previous versions of EQ6, this subroutine no longer has the
    !! function of determining the actual or most stable phase
    !! assemblage. Adjustments to the phase assemblage are now carried
    !! out by EQ6/eqphas.f, which as necessary makes repeated calls to
    !! the present subroutine. The present subroutine drives
    !! EQ6/optmzr.f and EQLIB/newton.f, which respectively carry out
    !! pre-Newton-Raphson optimization of starting estimates and hybrid
    !! Newton-Raphson iteration.
    !! Note that there are no attempts in this subroutine to diagnose
    !! the cause of any failure to converge. Such a failure might occur,
    !! for example, because a phase needs to be deleted from the phase
    !! assemblage. It is intended that this subroutine normally be
    !! called by EQ6/eqphas.f, which does perform such diagnostics.
    !! This subroutine is called by:
    !!   EQ6/eqphas.f
    !! Principal input:
    !! Principal output:
    !!   iter   = the number of hydrbid Newton-Raphson iterations
    !!              done by EQLIB/newton.f
    !!   ier    = error flag (returned from EQLIB/newton.f):
    !!              =  0  Okay
    !!              =  1  Encountered a zero matrix
    !!              =  2  Encountered a non-zero, computationally
    !!                      singular matrix
    !!              =  3  Iteration was diverging
    !!              =  4  Hit the maximum number of iterations (itermx)
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

    integer :: nfac(nbtmax)

    integer :: iapxt(nxtmax)
    integer :: ibpxt(nxtmax)
    integer :: ibswx(nbtmax)
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
    integer :: ndrs(ndrsmx)
    integer :: ndrsd(ndrsmx)
    integer :: ndrsx(ndrsmx)
    integer :: ndrsr(2,nstmax)
    integer :: ndrsrd(2,nstmax)
    integer :: ndrsrx(2,nstmax)
    integer :: ness(nessmx)
    integer :: nessr(2,nstmax)
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

    integer :: ibetmx
    integer :: idelmx
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
    integer :: kmt
    integer :: ko2gaq
    integer :: krdxsp
    integer :: kwater
    integer :: kx1
    integer :: kxt
    integer :: nbtd
    integer :: nbw
    integer :: nchlor
    integer :: nelect
    integer :: nhydr
    integer :: nhydx
    integer :: no2gaq
    integer :: nrdxsp
    integer :: ntpr

    logical :: qbassw
    logical :: qhawep
    logical :: qoptmz
    logical :: qpit75
    logical :: qredox
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
    real(kind=8) :: rhsvec(kmax)
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
    real(kind=8) :: bfje
    real(kind=8) :: bfxi
    real(kind=8) :: bgamx
    real(kind=8) :: bneg
    real(kind=8) :: bsigmm
    real(kind=8) :: delmax
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
    real(kind=8) :: press
    real(kind=8) :: screwd
    real(kind=8) :: sigmam
    real(kind=8) :: sigmmo
    real(kind=8) :: smp100
    real(kind=8) :: tempc
    real(kind=8) :: tempk
    real(kind=8) :: tolbt
    real(kind=8) :: toldl
    real(kind=8) :: xbarw
    real(kind=8) :: xbarwc
    real(kind=8) :: xbrwlc
    real(kind=8) :: xbrwlg

    ! Local variable declarations with global dimensioning.
    ! Note: the following variables are not used by EQ6 except to
    ! satisfy calls to the following subroutines shared with EQ3NR:
    !   EQLIB/betas.f
    !   EQLIB/newton.f
    ! These subroutines are called directly in this subroutine, and
    ! indirectly through EQ6/optmzr.f. Some of these variables are
    ! not arrays.
    integer :: isv_nbtmax
    integer :: isv_ntfxmx

    SAVE isv_nbtmax,isv_ntfxmx

    integer, dimension(:), allocatable :: ncosp
    integer, dimension(:), allocatable :: ntfx

    SAVE ncosp,ntfx

    real(kind=8), dimension(:), allocatable :: coval
    real(kind=8), dimension(:), allocatable :: tfx

    SAVE coval,tfx

    ! Local variable declarations.
    integer :: j2
    integer :: j3
    integer :: kcol
    integer :: n

    integer :: iebal
    integer :: irdxc3
    integer :: nredox
    integer :: ntfxt

    integer :: ilnobl

    logical :: qblamx
    logical :: qxbarw

    character(len=48) :: ubbig
    character(len=48) :: ubetmx
    character(len=48) :: ubneg
    character(len=8) :: uaff
    character(len=8) :: umoles
    character(len=8) :: ux8
    character(len=8) :: ux8a
    character(len=8) :: ux8b

    real(kind=8) :: actwlc
    real(kind=8) :: screwn
    real(kind=8) :: tolxpt
    real(kind=8) :: tolbig
    real(kind=8) :: tolneg

    data iebal  /0/
    data irdxc3 /0/
    data nredox /0/
    data ntfxt  /0/

    data qxbarw/.false./

    ! Set the screwn parameter.
    data screwn /0.10/

    ! The following are tolerance parameters for the optimization
    ! carried out by EQ6/optmzr.f. These are convergence tolerances
    ! within that subroutine; however, they are used in the present
    ! subroutine to determine if pre-Newton-Raphson iteration is
    ! necessary.
    data tolbig /10.0/,tolneg /-0.9/
    data tolxpt /10.0/

    data uaff /'aff     '/,umoles /'moles   '/

    ! Allocate or reallocate local work arrays as needed.
    if (.not.ALLOCATED(ncosp)) then
        ! Local work arrays are not allocated. Zero the saved
        ! array size variables. Note that only one array is tested
        ! to see if it is allocated. It is assumed that all local
        ! work arrays are either allocated or not.
        isv_nbtmax = 0
        isv_ntfxmx = 0
    else
        ! Local work arrays are allocated. Check to see if any of the
        ! array size variables have changed. If so, deallocate
        ! the corresponding local work arrays and zero the corresponding
        ! saved size variables.
        if (nbtmax .ne. isv_nbtmax) then
            DEALLOCATE(ncosp)
            DEALLOCATE(coval)
            isv_nbtmax = 0
        end if

        if (ntfxmx .ne. isv_ntfxmx) then
            DEALLOCATE(ntfx)
            DEALLOCATE(tfx)
            isv_ntfxmx = 0
        end if
    end if

    ! At this point, the saved array size values are zero if the
    ! corresponding arrays need to be allocated.
    if (isv_nbtmax .eq. 0) then
        ALLOCATE(ncosp(nbtmax))
        ALLOCATE(coval(nbtmax))
        isv_nbtmax = nbtmax
    end if

    if (isv_ntfxmx .eq. 0) then
        ALLOCATE(ntfx(ntfxmx))
        ALLOCATE(tfx(ntfxmx))
        isv_ntfxmx = ntfxmx
    end if

    ! Zero the contents of the local work arrays.
    do n = 1,nbtmax
        ncosp(n) = 0
    end do

    do n = 1,nbtmax
        dlogxw(n) = 0.
        coval(n) = 0.
    end do

    do n = 1,ntfxmx
        ntfx(n) = 0
    end do

    do n = 1,ntfxmx
        tfx(n) = 0.
    end do

    ! Set up the ulbeta and uldel arrays.
    call initcv(ulbeta,kbt,umoles)
    call initcv(uldel,kbt,umoles)

    do kcol = km1,kxt
        ulbeta(kcol) = uaff
        uldel(kcol) = umoles
    end do

    ! Set the qblamx flag, which indicates the presence of solid
    ! solutions in the currently specified phase assemblage.
    qblamx = kxt .ge. kx1

    ! Expand the system description.
    call ncmpex(acflg,act,actlg,cdrs,cegexs,cgexj,conc,conclg,cpgexs,egexjc,egexjf,egexs,eps100,fo2,fo2lg,fsort,fugac,fugalg,iern1,iern2,ietmax,ifrn1,ifrn2,igas,igstak,iindx1,ilrn1,ilrn2,imrn1,imrn2,istack,ixrn1,ixrn2,jcsort,jern1,jern2,jetmax,jflag,jgext,jgsort,jgstak,jjsort,jpflag,jsflag,jsitex,jssort,jstack,kbt,kdim,kelect,kmax,km1,ko2gaq,kwater,kxt,loph,losp,lsort,mgext,mrgexs,mtb,moph,mosp,narn1,narn2,nbasp,nbt,nbtmax,ncmpr,ndrs,ndrsmx,ndrsr,nelect,nern1,nern2,netmax,ngexsa,ngext,ngrn1,ngrn2,ngt,ngtmax,noutpt,no2gaq,nphasx,npt,nptmax,nst,nstmax,nttyo,omega,omeglg,press,qxbarw,q6mode,ugexj,ugexmo,uphase,uspec,xbar,xbarlg,xbarw,xbarwc,xbrwlc,xbrwlg,xlks,zchar,zgexj,zvclg1,zvec1)

    ! Compute the residual functions prior to iteration.
    call betas(acflg,actlg,afcnst,alpha,amtb,bbig,beta,betamx,bneg,cdrs,conc,conclg,coval,csts,eh,ehfac,fo2lg,ibetmx,iebal,iindx1,irdxc3,jcsort,jflag,jsflag,jssort,kbt,kdim,kelect,khydr,kmax,km1,ko2gaq,kwater,kxt,mtb,mosp,narn1,narn2,nbasp,nbtmax,ncosp,ndrs,ndrsmx,ndrsr,nelect,nern1,nern2,nhydr,noutpt,no2gaq,nredox,nst,nstmax,nsts,nstsmx,nstsr,ntfx,ntfxmx,ntfxt,nttyo,omega,qredox,q6mode,tfx,ubbig,ubneg,ubetmx,uspec,uzvec1,weight,xbrwlg,xlke,xlks,zchar)

    if (qoptmz) then
        if (betamx.gt.tolxpt .or. bbig.gt.tolbig .or. bneg.lt.tolneg)    then
            ! XXX
            !          write (noutpt,3000) betamx,tolxpt
            !          write (nttyo,3000) betamx,tolxpt
            ! 3000     format(//' ***** betamx= ',1pe11.4,' and tolxpt= ',1pe11.4,
            !     $    ' ****')
            !          j2 = ilnobl(ubetmx)
            !          write (noutpt,3005) ubetmx(1:j2)
            !          write (nttyo,3005) ubetmx(1:j2)
            ! 3005     format('   ***** ubetmx= ',a,' *****',/)
            !          write (noutpt,3010) bbig,tolbig
            !          write (nttyo,3010) bbig,tolbig
            ! 3010     format(/' ***** bbig= ',1pe11.4,' and tolbig= ',1pe11.4,
            !     $    ' ****')
            !          j2 = ilnobl(ubbig)
            !          write (noutpt,3015) ubbig(1:j2)
            !          write (nttyo,3015) ubbig(1:j2)
            ! 3015     format('   ***** ubbig= ',a,' *****',/)
            !          write (noutpt,3020) bneg,tolneg
            !          write (nttyo,3020) bneg,tolneg
            ! 3020     format(/' ***** bneg= ',1pe11.4,' and tolneg= ',1pe11.4,
            !     $    ' ****')
            !          j2 = ilnobl(ubneg)
            !          write (noutpt,3025) ubneg(1:j2)
            !          write (nttyo,3025) ubneg(1:j2)
            ! 3025     format('   ***** ubbneg= ',a,' *****',/)
            ! XXX
            !          Optimize the iteration variables before starting hybrid
            !          Newton-Raphson iteration. The phase assemblage is fixed
            !          in this process.
            call optmzr(aamatr,abar,acflg,acflgo,act,actlg,adh,adhh,adhv,afcnst,al10,alpha,amtb,aphi,avcnst,azero,a3bar,a3bars,bacfmx,bbig,bdh,bdhh,bdhv,bdot,bdoth,bdotv,beta,betamx,bgamx,bneg,bpx,cco2,cdrs,cdrsx,cdrtw,cdrw,cegexs,cgexj,cjbasp,cnufac,conc,conclg,cpgexs,cscale,csts,coval,delvec,dlogxw,egexjc,egexjf,egexs,ehfac,eps100,fje,fjeo,fo2,fo2lg,fsort,fugac,fugalg,fxi,fxio,gmmatr,ibetmx,ibpxt,ibswx,ielam,iern1,iern2,ifcphi1,ifcphi2,ifnnn,ifn2n,ifpsi1,ifpsi2,ifrn1,ifrn2,ifzeta,igas,igstak,iindx1,ilcphi1,ilcphi2,ilnnn,iln2n,ilpsi1,ilpsi2,ilrn1,ilrn2,ilzeta,imrn1,imrn2,insgf,iodb,iopg,iopt,ipch,ipcv,ipivot,ipndx1,irdxc3,istack,ixbasp,ixrn1,ixrn2,izmax,jcsort,jern1,jern2,jflag,jgext,jgsort,jgstak,jjsort,jpflag,jsflag,jsitex,jsol,jssort,jstack,kbt,kction,kdim,kelect,khydr,km1,kmt,ko2gaq,kwater,kx1,kxt,loph,losp,lsort,mgext,moph,mosp,mrgexs,mtb,narn1,narn2,narxt,nat,nbasp,nbaspd,nbaspx,nbw,nbt,nbtd,nchlor,ncmpr,ncosp,ndrs,ndrsx,ndrsr,ndrsrd,ndrsrx,nelect,nern1,nern2,net,ngexsa,nfac,ngext,nhydr,nhydx,ngrn1,ngrn2,ngt,noutpt,no2gaq,nphasx,npt,nst,nsts,nstsr,ntfx,ntfxt,ntpr,nttyo,omega,omeglg,press,qbassw,qblamx,qhawep,qpit75,qredox,q6mode,rhsvec,sigmam,sigmmo,smp100,tempc,tempk,tfx,tolbig,tolneg,tolxpt,ubacmx,ubbig,ubetmx,ubgamx,ubneg,ugexj,ugexmo,uphase,uspec,uzvec1,weight,wfac,xbar,xbarlg,xbarw,xbarwc,xbrwlc,xbrwlg,zchar,zchcu6,zchsq2,zgexj,zvclg1,zvec1)
        end if
    end if

    ! Do hybrid Newton-Raphson iteration for a given phase assemblage.
    if (iodb(4) .ge. 1) then
        write (noutpt,1600)
1600 format(/,' Starting hybrid Newton-Raphson iteration.',/)
    end if

    call newton(aamatr,abar,acflg,acflgo,act,actlg,actwlc,adh,adhh,adhv,afcnst,alpha,al10,amtb,aphi,azero,a3bar,a3bars,bacfmx,bbig,beta,betamx,betao,bdh,bdhh,bdhv,bdot,bdoth,bdotv,bfje,bfxi,bgamx,bneg,bpx,bsigmm,cco2,cdrs,cdrtw,cdrw,cegexs,cgexj,cjbasp,cnufac,conc,conclg,coval,cpgexs,csts,delam,delmax,delvco,delvec,dgpit,dlogxw,dpelm,dpslm,dselm,egexjc,egexjf,egexs,eh,ehfac,elam,eps100,fje,fjeo,fo2,fo2lg,fsort,fugac,fugalg,fxi,fxio,gmmatr,gpit,ibpxt,idelmx,iebal,ielam,ier,iern1,iern2,ifcphi1,ifcphi2,ifnnn,ifn2n,ifpsi1,ifpsi2,ifrn1,ifrn2,ifzeta,igas,igstak,iindx1,ilcphi1,ilcphi2,ilnnn,iln2n,ilpsi1,ilpsi2,ilrn1,ilrn2,ilzeta,imrn1,imrn2,insgf,iodb,iopg,ipivot,ipndx1,irdxc3,istack,iter,itermx,ixbasp,ixrn1,ixrn2,izmax,jcsort,jern1,jern2,jflag,jgext,jgsort,jgstak,jjsort,jpflag,jsflag,jsitex,jsol,jssort,jstack,kbt,kction,kdim,kelect,khydr,km1,kmt,ko2gaq,kwater,kx1,kxt,loph,losp,lsort,mgext,moph,mosp,mrgexs,mtb,nalpha,napt,narn1,narn2,nbasp,nbt,nbw,nchlor,ncmpr,ncosp,ndrs,ndrsr,nelect,nern1,nern2,net,ngexsa,ngext,ngrn1,ngrn2,ngt,nhydr,nmut,nmux,nmxi,nmxx,noutpt,no2gaq,nphasx,npt,nredox,nslt,nslx,nst,nsts,nstsr,nsxi,nsxx,ntfx,ntfxt,nttyo,omega,omeglg,palpha,pelm,pmu,press,pslamn,pslm,qhawep,qpit75,qredox,q6mode,rhsvec,screwd,screwn,selm,sigmam,sigmmo,tempk,tfx,tolbt,toldl,ubacmx,ubbig,ubetmx,ubgamx,ubneg,ugexj,ugexmo,ulbeta,uldel,uphase,uspec,uzvec1,weight,wfac,xbar,xbarlg,xbarw,xbarwc,xbrwlc,xbrwlg,xlke,xlks,zchar,zchsq2,zchcu6,zgexj,zvclg1,zvec1)

    if (ier .gt. 0) then
        if (iodb(4) .gt. 0) then
            write (noutpt,1610)
1610 format(/' * Note - (EQ6/eqcalc) Hybrid Newton-Raphson',' iteration failed.')

            ux8a = ' '
            write (ux8a,'(i7)') iter
            call lejust(ux8a)
            j2 = ilnobl(ux8a)

            if (ier .eq. 1) then
                write (noutpt,1630) ux8a(1:j2)
1630 format(7x,'after ',a,' iterations because a zero matrix',' was encountered. This is',/7x,'probably due to a',' programming error.')
            else if (ier .eq. 2) then
                write (noutpt,1640) ux8a(1:j2)
1640 format(7x,'after ',a,' iterations because a non-zero,',' computationally singular',/7x,'matrix was encountered.')
            else if (ier .eq. 3) then
                write (noutpt,1650) ux8a(1:j2)
1650 format(7x,'after ',a,' iterations because the code',' detected that',/7x,'iteration was diverging.')
            else if (ier .eq. 4) then
                write (noutpt,1660) ux8a(1:j2)
1660 format(7x,'after ',a,' iterations because the maximum',' number of iterations',/7x,'was done.')
            else
                ux8b = ' '
                write (ux8b,'(i7)') ier
                call lejust(ux8b)
                j3 = ilnobl(ux8b)
                write (noutpt,1670) ux8a(1:j2),ux8b(1:j3)
1670 format(7x,'after ',a,' iterations because an unknown event',' occurred. The ier',/7x,'error code has the unknown value',' ',a,'. This condition is a',/7x,'programming error.')
            end if
        end if

        go to 999
    end if

    if (iodb(4) .ge. 1) then
        ux8 = ' '
        write (ux8,'(i7)') iter
        call lejust(ux8)
        j2 = ilnobl(ux8)
        write (noutpt,1690) ux8(1:j2)
1690 format('   Done. Hybrid Newton-Raphson iteration converged in ',a,' iterations.',/)
    end if

999 continue
end subroutine eqcalc