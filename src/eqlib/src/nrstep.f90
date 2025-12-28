subroutine nrstep(aamatr,acflg,act,actlg,afcnst,alpha,al10,amtb,bbig,beta,betamx,betao,betfnc,betmxo,bneg,bpx,btfcnr,cdrs,cdrtw,cdrw,cegexs,cgexj,cjbasp,cnufac,conc,conclg,cpgexs,csts,coval,delfnc,delmax,delvco,delvec,dlogxw,egexjc,egexjf,egexs,eh,ehfac,eps100,fo2,fo2lg,fsort,fugac,fugalg,gmmatr,ibetmx,ibpxmx,idelmx,iebal,ier,iern1,iern2,ietmax,ifrn1,ifrn2,igas,igstak,iindx1,ilrn1,ilrn2,imrn1,imrn2,iodb,ipivot,ipndx1,irdxc3,istack,iter,itermx,ixbasp,ixrn1,ixrn2,jcsort,jern1,jern2,jetmax,jflag,jgext,jgsort,jgstak,jjsort,jpflag,jsflag,jsitex,jssort,jstack,kbt,kction,kdim,kelect,khydr,kmax,kmt,km1,ko2gaq,kwater,kxt,kx1,loph,losp,lsort,mgext,moph,mosp,mrgexs,mtb,narn1,narn2,nbasp,nbt,nbtmax,nbw,ncmpr,ncosp,ndrs,ndrsmx,ndrsr,negbfc,negdfc,negxfc,nelect,nern1,nern2,netmax,ngexsa,ngext,ngrn1,ngrn2,ngt,ngtmax,nhydr,nodbmx,noibfc,noutpt,no2gaq,npconv,nphasx,npobfc,npodfc,npt,nptmax,nredox,nst,nstmax,nsts,nstsmx,nstsr,ntfx,ntfxmx,ntfxt,nttyo,nxtmax,omega,omeglg,press,qcacf,qcbeta,qredox,qxbarw,q6mode,rhsvec,rlxfac,screwd,screwn,sigmam,sigmmo,tfx,ubetmx,ubbig,ubneg,ugexj,ugexmo,uphase,uspec,uzvec1,weight,xbar,xbarlg,xbarw,xbarwc,xbrwlc,xbrwlg,xlke,xlks,zchar,zgexj,zvclg1,zvec1)
    !! This subroutine performs one Newton-Raphson step.
    !! This subroutine is called by:
    !!   EQLIB/newton.f
    !! Principal input:
    !!   iodb   = array of debugging print options
    !!   ipivot = the pivot vector
    !!   screwd = under-relaxation control parameter. It is used to
    !!              reduce the magnitude of the delvec vector, if
    !!              necessary, so that the magnitude of the largest
    !!              element of that vector does not exceed screwd
    !!   screwn = under-relaxation control parameter
    !!   itermx = maximum number of iterations
    !!   kdim   = dimension of aamatr
    !!   kmax   = maximum dimension of aamatr
    !!   qcbeta = true if the convergence criterion on betamx was met
    !!   qcacf  = true if the convergence criterion on activity
    !!              coefficients was met
    !! Principal output:
    !!   aamatr = Jacobian matrix
    !!   gmmatr = copy of aamatr
    !!   rhsvec = right hand side vector
    !!   alpha  = residual function array
    !!   idelmx = kcol index corresponding to delmax
    !!   delfnc = convergence function, defined by reference to delmax
    !!   betfnc = convergence function, defined by reference to betamx
    !!   uzvec1 = name array corresponding to zvclg1
    !!   ier    = error flag:
    !!              =  0  Okay
    !!              =  1  Encountered a zero matrix
    !!              =  2  Encountered a non-zero, computationally
    !!                      singular matrix
    !!              =  3  Iteration was diverging
    !!              =  4  Hit the maximum number of iterations (itermx)
    !! Principal input/output:
    !!   iter   = the number of Newton-Raphson iterations
    !!   negdfc = the number of successive iterations that delmax
    !!              has been greater than zero and the corresponding
    !!              convergence function delfnc has been less than
    !!              or equal to zero
    !!   negbfc = the number of successive iterations that betamx
    !!              has been greater than zero and the corresponding
    !!              convergence function betfnc has been less than
    !!              or equal to zero
    !!   negxfc = the number of successive iterations that either
    !!              delmax has been greater than zero while delfnc
    !!              has been less than or equal to zero, or betamx
    !!              has been greater than zero while betfnc has
    !!              been less than or equal to zero
    !!   noibfc = the number of successive iterations that betamx
    !!              has been greater than zero and the corresponding
    !!              convergence function betfnc has been less than
    !!              1.e-6
    !!   npconv = number of successive steps in which the convergence
    !!              criteria on the residual norm has been satisfied
    !!   npodfc = the number of successive iterations that the
    !!              convergence function delfnc has been greater
    !!              than zero
    !!   npobfc = the number of successive iterations that the
    !!              convergence function betfnc has been greater
    !!              than zero
    !!   zvclg1 = the 'log z' array, the array corrected by
    !!              Newton-Raphson iteration
    !!   beta   = normalized residual function array
    !!   delvec = correction array
    !!   betao  = old beta array
    !!   delvco = old delvec array
    !!   betamx = max norm of the beta array
    !!   delmax = max norm of the delvec array
    !!   rlxfac = under-relaxation factor
    implicit none

    ! Calling sequence variable declarations.
    integer :: ibpxmx
    integer :: ietmax
    integer :: jetmax
    integer :: kmax
    integer :: nbtmax
    integer :: ndrsmx
    integer :: netmax
    integer :: ngtmax
    integer :: nodbmx
    integer :: nptmax
    integer :: nstmax
    integer :: nstsmx
    integer :: ntfxmx
    integer :: nxtmax

    integer :: noutpt
    integer :: nttyo

    integer :: igstak(ngtmax)
    integer :: iindx1(kmax)
    integer :: iodb(nodbmx)
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

    integer :: ibetmx
    integer :: idelmx
    integer :: iebal
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
    integer :: irdxc3
    integer :: iter
    integer :: itermx
    integer :: ixrn1
    integer :: ixrn2
    integer :: kbt
    integer :: kdim
    integer :: kelect
    integer :: khydr
    integer :: kmt
    integer :: km1
    integer :: ko2gaq
    integer :: kwater
    integer :: kxt
    integer :: kx1
    integer :: narn1
    integer :: narn2
    integer :: nbt
    integer :: nbw
    integer :: negbfc
    integer :: negdfc
    integer :: negxfc
    integer :: nelect
    integer :: nern1
    integer :: nern2
    integer :: ngrn1
    integer :: ngrn2
    integer :: ngt
    integer :: nhydr
    integer :: noibfc
    integer :: no2gaq
    integer :: npconv
    integer :: npobfc
    integer :: npodfc
    integer :: npt
    integer :: nredox
    integer :: nst
    integer :: ntfxt

    logical :: qcacf
    logical :: qcbeta
    logical :: qredox
    logical :: qxbarw
    logical :: q6mode

    character(len=48) :: uspec(nstmax)
    character(len=48) :: uzvec1(kmax)
    character(len=48) :: ubbig
    character(len=48) :: ubneg
    character(len=48) :: ubetmx
    character(len=24) :: ugexmo(netmax)
    character(len=24) :: uphase(nptmax)
    character(len=8) :: ugexj(jetmax,netmax)

    real(kind=8) :: aamatr(kmax,kmax)
    real(kind=8) :: acflg(nstmax)
    real(kind=8) :: act(nstmax)
    real(kind=8) :: actlg(nstmax)
    real(kind=8) :: alpha(kmax)
    real(kind=8) :: amtb(nbtmax)
    real(kind=8) :: beta(kmax)
    real(kind=8) :: betao(kmax)
    real(kind=8) :: bpx(ibpxmx,nxtmax)
    real(kind=8) :: cdrs(ndrsmx)
    real(kind=8) :: cdrtw(nstmax)
    real(kind=8) :: cdrw(nstmax)
    real(kind=8) :: cegexs(ietmax,jetmax,netmax)
    real(kind=8) :: cgexj(jetmax,netmax)
    real(kind=8) :: cjbasp(nbtmax)
    real(kind=8) :: conc(nstmax)
    real(kind=8) :: conclg(nstmax)
    real(kind=8) :: cnufac(nstmax)
    real(kind=8) :: coval(nbtmax)
    real(kind=8) :: cpgexs(ietmax,jetmax,netmax)
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
    real(kind=8) :: mgext(jetmax,netmax)
    real(kind=8) :: moph(nptmax)
    real(kind=8) :: mosp(nstmax)
    real(kind=8) :: mrgexs(ietmax,jetmax,netmax)
    real(kind=8) :: mtb(nbtmax)
    real(kind=8) :: rhsvec(kmax)
    real(kind=8) :: tfx(ntfxmx)
    real(kind=8) :: weight(nstmax)
    real(kind=8) :: xbar(nstmax)
    real(kind=8) :: xbarlg(nstmax)
    real(kind=8) :: xlks(nstmax)
    real(kind=8) :: zchar(nstmax)
    real(kind=8) :: zgexj(jetmax,netmax)
    real(kind=8) :: zvclg1(kmax)
    real(kind=8) :: zvec1(kmax)

    real(kind=8) :: afcnst
    real(kind=8) :: al10
    real(kind=8) :: bbig
    real(kind=8) :: betamx
    real(kind=8) :: betfnc
    real(kind=8) :: betmxo
    real(kind=8) :: bneg
    real(kind=8) :: btfcnr
    real(kind=8) :: delfnc
    real(kind=8) :: delmax
    real(kind=8) :: eh
    real(kind=8) :: ehfac
    real(kind=8) :: eps100
    real(kind=8) :: fo2
    real(kind=8) :: fo2lg
    real(kind=8) :: omega
    real(kind=8) :: omeglg
    real(kind=8) :: press
    real(kind=8) :: rlxfac
    real(kind=8) :: screwd
    real(kind=8) :: screwn
    real(kind=8) :: sigmam
    real(kind=8) :: sigmmo
    real(kind=8) :: xbarw
    real(kind=8) :: xbarwc
    real(kind=8) :: xbrwlc
    real(kind=8) :: xbrwlg
    real(kind=8) :: xlke

    ! Local variable declarations.
    integer :: i
    integer :: jlen
    integer :: kcol
    integer :: krow
    integer :: ncut

    integer :: iarmxn

    logical :: qconc1
    logical :: qconc2
    logical :: qconc3
    logical :: qpr

    character(len=56) :: uspn56

    real(kind=8) :: delmxo
    real(kind=8) :: divfmx
    real(kind=8) :: divfnc
    real(kind=8) :: fx
    real(kind=8) :: rdx
    real(kind=8) :: sxu

    data qpr    /.false./

    ! Save the current values of the beta and delvec arrays for use in
    ! under-relaxation schemes.
    call copyaa(beta,betao,kdim)
    call copyaa(delvec,delvco,kdim)
    delmxo = delmax

    ! Check to see if iteration should be stopped because there is
    ! little probability of further improvement. If any of the
    ! following conditions cause iteration to stop, the result is
    ! considered acceptable if it satisfies the convergence criteria
    ! on residual functions, but not on the delvec vector.
    ! Note: betfnc is the convergence function based on the beta
    ! residual vector, and delfnc is that based on the delvec
    ! correction vector. If convergence is occurring, each function
    ! should approach unity (the maximum possible value).
    ! Update counters. These all count the number of times in a row
    ! that a condition has been satisfied.
    if (qcbeta .and. qcacf) then
        npconv = npconv + 1
    else
        npconv = 0
    end if

    if (delmax.gt.0. .and. delfnc.le.0.) then
        negdfc = negdfc + 1
    else
        negdfc = 0
    end if

    if (betamx.gt.0. .and. betfnc.le.0.) then
        negbfc = negbfc + 1
    else
        negbfc = 0
    end if

    if ((betamx.gt.0. .and. betfnc.le.0.) .or.  (delmax.gt.0. .and. delfnc.le.0.)) then
        negxfc = negxfc + 1
    else
        negxfc = 0
    end if

    if (delfnc .gt. 0.) then
        npodfc = npodfc + 1
    else
        npodfc = 0
    end if

    if (betfnc .gt. 0.) then
        npobfc = npobfc + 1
    else
        npobfc = 0
    end if

    if (betamx.gt.0. .and. betfnc.le.1.e-6) then
        noibfc = noibfc + 1
    else
        noibfc = 0
    end if

    ! Stop iteration if the convergence criteria on residual functions
    ! (derived from the beta vector) have been satisfied for four
    ! consecutive iterations.
    if (npconv .ge. 4) then
        go to 300
    end if

    ! Stop iteration if the behavior of delfnc indicates that the
    ! iteration is diverging, and this indication is not contradicted
    ! by the behavior of betfnc. Technically, delfnc is negative or
    ! zero for twelve consecutive iterations and betfnc has not been
    ! positive for three or more consecutive iterations.
    if (negdfc.ge.12 .and. npobfc.lt.3) then
        go to 300
    end if

    ! Stop iteration if the behavior of betfnc indicates that the
    ! iteration is diverging, and this indication is not contradicted
    ! by the behavior of delfnc. Technically, betfnc is negative or
    ! zero for twelve consecutive iterations and delfnc has not been
    ! positive for three or more consecutive iterations.
    if (negbfc.ge.12 .and. npodfc.lt.3) then
        go to 300
    end if

    ! Stop iteration if the behavior of delfnc and betfnc together
    ! suggests that the iteration is diverging. Technically, either
    ! delfnc or betfnc is negative or zero for twelve consecutive
    ! iterations, and iter is greater than or equal to 40.
    if (negxfc.ge.12 .and. iter.ge.40) then
        go to 300
    end if

    ! Stop iteration if the behavior of betfnc indicates that no
    ! significant improvement is occurring. Technically, betfnc is
    ! less than or equal to 1.e-6 for twelve consecutive iterations,
    ! and iter is greater than or equal to 20.
    if (noibfc.ge.12 .and. iter.ge.20) then
        go to 300
    end if

    ! Stop iteration if the maximum number of iterations have
    ! been done.
    if (iter .ge. itermx) then
        write (noutpt,1000) itermx
        write (nttyo,1000) itermx
1000 format(/' * Note - (EQLIB/nrstep) Have completed ',i4,' hybrid',/7x,'Newton-Raphson iterations. This is the maximum number',' permitted.')

        ier = 4
        go to 999
    end if

    ! Do another Newton-Raphson iteration.
    iter = iter + 1

    ! Calculate the Jacobian matrix (aamatr).
    call matrix(aamatr,al10,bpx,cdrs,cdrtw,cdrw,cjbasp,cnufac,conc,csts,dlogxw,eps100,ibpxmx,iebal,iern1,ietmax,iindx1,ipndx1,irdxc3,ixbasp,ixrn1,jcsort,jern1,jern2,jetmax,jflag,jjsort,jsitex,kbt,kction,kdim,kelect,khydr,kmax,kmt,km1,ko2gaq,kwater,kxt,kx1,mosp,narn1,narn2,nbasp,nbt,nbtmax,nbw,ncosp,ndrs,ndrsmx,ndrsr,nelect,nern1,nern2,netmax,noutpt,no2gaq,nphasx,nredox,nst,nstmax,nsts,nstsmx,nstsr,ntfx,ntfxmx,ntfxt,nttyo,nxtmax,omega,qredox,q6mode,tfx,ugexmo,uspec,weight,xbar,xbarw,xbarwc,zchar)

    ! Get the right-hand-side vector (rhsvec = -alpha).
    do krow = 1,kdim
        rhsvec(krow) = -alpha(krow)
    end do

    if (iodb(4) .ge. 3) then
        write (noutpt,1010)
1010 format(/16x,'--- The aamatr array, with rhsvec ---',/)

        do krow = 1,kdim
            ! Calling sequence substitutions:
            !   uzvec1(krow) for unam48
            call fmspnx(jlen,uzvec1(krow),uspn56)
            write (noutpt,1020) krow,uspn56(1:jlen)
1020 format(1x,i3,2x,a)

            write (noutpt,1030) (aamatr(krow,i),i = 1,kdim),rhsvec(krow)
1030 format(( 6(2x,g10.3)) )

            write(noutpt,1040)
1040 format(' ')
        end do
    end if

    ! Solve the Jacobian system. If EQLIBU/msolvr.f can't solve the
    ! matrix, it is because the matrix is either zero (ier = 1) or
    ! non-zero, but computationally singular (ier = 2).
    call msolvr(aamatr,delvec,gmmatr,ier,ipivot,kdim,kmax,noutpt,nttyo,qpr,rhsvec)

    if (ier .gt. 0) then
        go to 999
    end if

    ! Find the max norm of the correction vector (delvec).
    ! In the limit of the solution, the (unrelaxed) correction
    ! term bounds the error in the variables being calculated.
    idelmx = iarmxn(delvec,kdim)
    delmax = 0.

    if (idelmx .gt. 0) then
        delmax = abs(delvec(idelmx))
    end if

    ! Calculate the delvec improvement function (delfnc).
    delfnc = 0.

    if (delmxo .gt. 0.) then
        delfnc = (delmxo -delmax)/delmxo
    end if

    ! Apply under-relaxation according to some simple criteria.
    ! The under-relaxation at this point is done before recomputing
    ! the residual functions. A subsequent under-relaxation
    ! will be done to insure that the residual functions do not
    ! increase too much.
    rlxfac = 1.0

    ! Limit the magnitude of change per iteration.
    if (delmax .gt. screwd) then
        rlxfac = screwd/delmax
    end if

    ! If within the first eight iterations and the residual norm
    ! (betamx) is large, go slow.
    if (iter.le.8 .and. betamx.gt.screwn) then
        rlxfac = 0.5*rlxfac
    end if

    qconc1 = sigmam.ge.2.0 .or. sigmmo.ge.2.0
    qconc2 = sigmam.ge.8.0 .or. sigmmo.ge.8.0
    qconc3 = sigmam.ge.12.0 .or. sigmmo.ge.12.0

    if (qconc1) then
        ! Force some minimum under-relaxation for modestly concentrated
        ! solutions.
        if (iter .le. 15) then
            sxu = 0.50
            rlxfac = min(sxu,rlxfac)
        end if
    end if

    if (qconc2) then
        ! Force some minimum under-relaxation for more concentrated
        ! solutions.
        if (iter .le. 15) then
            sxu = 0.35
            rlxfac = min(sxu,rlxfac)
        end if
    end if

    if (qconc3) then
        ! Force some minimum under-relaxation for highly concentrated
        ! solutions.
        if (iter .le. 15) then
            sxu = 0.25
            rlxfac = min(sxu,rlxfac)
        end if
    end if

    ! Apply correction terms. The label below is a return point
    ! under-relaxation which restricts the growth of the residual
    ! functions.
    ncut = 0
140 continue

    if (iodb(4) .ge. 2) then
        write (noutpt,1050) ncut
    end if

1050 format(/3x,'ncut= ',i3,/)

    do kcol = 1,kdim
        zvclg1(kcol) = zvclg1(kcol) + rlxfac*delvec(kcol)
    end do

    if (iodb(4) .ge. 2) then
        write (noutpt,1060)
1060 format(/8x,'Name',38x,'del',9x,'zvclg1'/)

        do kcol = 1,kdim
            ! Calling sequence substitutions:
            !   uzvec1(kcol) for unam48
            call fmspnx(jlen,uzvec1(kcol),uspn56)
            jlen = min(jlen,38)
            rdx = rlxfac*delvec(kcol)
            write (noutpt,1062) kcol,uspn56(1:jlen),rdx,zvclg1(kcol)
1062 format(1x,i3,2x,a,t46,1pe12.5,2x,e12.5)
        end do

        write (noutpt,1064)
1064 format(1x)
    end if

    ! Expand the system and recalculate the Newton-Raphson residual
    ! functions. It is important here to keep the activity coefficients
    ! constant,that is, not to update them until all under-relaxation
    ! has been completed.
    ! Recalculate the concentrations, etc., of dependent species.
    call ncmpex(acflg,act,actlg,cdrs,cegexs,cgexj,conc,conclg,cpgexs,egexjc,egexjf,egexs,eps100,fo2,fo2lg,fsort,fugac,fugalg,iern1,iern2,ietmax,ifrn1,ifrn2,igas,igstak,iindx1,ilrn1,ilrn2,imrn1,imrn2,istack,ixrn1,ixrn2,jcsort,jern1,jern2,jetmax,jflag,jgext,jgsort,jgstak,jjsort,jpflag,jsflag,jsitex,jssort,jstack,kbt,kdim,kelect,kmax,km1,ko2gaq,kwater,kxt,loph,losp,lsort,mgext,mrgexs,mtb,moph,mosp,narn1,narn2,nbasp,nbt,nbtmax,ncmpr,ndrs,ndrsmx,ndrsr,nelect,nern1,nern2,netmax,ngexsa,ngext,ngrn1,ngrn2,ngt,ngtmax,noutpt,no2gaq,nphasx,npt,nptmax,nst,nstmax,nttyo,omega,omeglg,press,qxbarw,q6mode,ugexj,ugexmo,uphase,uspec,xbar,xbarlg,xbarw,xbarwc,xbrwlc,xbrwlg,xlks,zchar,zgexj,zvclg1,zvec1)

    call betas(acflg,actlg,afcnst,alpha,amtb,bbig,beta,betamx,bneg,cdrs,conc,conclg,coval,csts,eh,ehfac,fo2lg,ibetmx,iebal,iindx1,irdxc3,jcsort,jflag,jsflag,jssort,kbt,kdim,kelect,khydr,kmax,km1,ko2gaq,kwater,kxt,mtb,mosp,narn1,narn2,nbasp,nbtmax,ncosp,ndrs,ndrsmx,ndrsr,nelect,nern1,nern2,nhydr,noutpt,no2gaq,nredox,nst,nstmax,nsts,nstsmx,nstsr,ntfx,ntfxmx,ntfxt,nttyo,omega,qredox,q6mode,tfx,ubbig,ubneg,ubetmx,uspec,uzvec1,weight,xbrwlg,xlke,xlks,zchar)

    if (iodb(4) .ge. 2) then
        write (noutpt,1080)
1080 format(/8x,'Name',37x,'Beta',/)

        do kcol = 1,kdim
            ! Calling sequence substitutions:
            !   uzvec1(kcol) for unam48
            call fmspnx(jlen,uzvec1(kcol),uspn56)
            jlen = min(jlen,38)
            write (noutpt,1085) kcol,uspn56(1:jlen),beta(kcol)
1085 format(1x,i3,2x,a,t46,1pe12.5)
        end do

        write (noutpt,1090)
1090 format(/1x)
    end if

    btfcnr = 0.

    if (betmxo .gt. 0.) then
        btfcnr = (betmxo - betamx)/betmxo
    end if

    ! Use under-relaxation to control residual growth, if any.
    ! Terminate iteration if residual tracking indicates that the
    ! iteration is diverging.
    if (iter .le. 12) then
        divfmx = 0.

        do kcol = 1,kdim
            divfnc = abs(beta(kcol)) - betao(kcol)

            if (divfmx .gt. divfnc) then
                divfmx = divfnc
            end if
        end do

        if (divfmx .gt. 0.10) then
            if (ncut .le. 7) then
                ncut = ncut + 1

                if (iodb(4) .ge. 2) then
                    write (noutpt,1100)
                end if

1100 format('  Residual tracking requires under-relaxation.')

                fx = 1./rlxfac

                do kcol = 1,kdim
                    zvclg1(kcol) = zvclg1(kcol) - fx*delvec(kcol)
                end do

                rlxfac = 0.25*rlxfac
                go to 140
            else if (iter .ge. 8) then
                write (noutpt,1110)
1110 format('  Residual tracking indicates that iteration',' is diverging.')

                go to 300
            end if
        end if
    end if

    if (iodb(4) .ge. 2) then
        write (noutpt,1120)
    end if

1120 format(' A Newton-Raphson correction has been completed.',/)

    go to 999

300 continue

    ! The code detected that iteration was diverging.
    ier = 3

999 continue
end subroutine nrstep