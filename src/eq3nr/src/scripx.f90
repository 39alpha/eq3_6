subroutine scripx(abar,acflg,act,actlg,adh,afcnst,affpd,affsd,ahrc,alki,apx,atwt,a3bar,a3bars,bpx,cdrsd,cegexs,cess,conc,conclg,coval,csts,ctb,cteaq,egexjc,egexjf,egexpa,egexpc,egexs,egexw,ehfac,ehrc,eps100,eh,farad,fje,fo2,fo2lg,fo2lrc,fugac,fugalg,fxi,iapxmx,ibpxmx,iebal,iern1,iern2,ietmax,igas,iktmax,iopg,iopr,iopt,ilrn1,ilrn2,imrn1,imrn2,ixrn1,ixrn2,jcsort,jern1,jern2,jetmax,jflag,jflagd,jflgi,jfleba,jgext,jgsort,jpflag,jsflag,jsol,jsomax,kern1,kern2,ketmax,kgexsa,mlmrra,mrmlra,moph,mosp,mte,mteaq,mwtsp,narn1,narn2,natmax,nbasp,nbaspd,nbt,nbtmax,nchlor,ncmpr,nct,nctmax,ndrsd,ndrsmx,ndrsrd,nelect,ness,nessmx,nessr,net,neti,netmax,ngexpi,ngexsa,ngext,ngrn1,ngrn2,ngt,ngtmax,nhydr,nhydx,nopgmx,noprmx,noptmx,noutpt,no2gaq,npnxp,npt,nptmax,nrdxsp,nst,nstmax,nsts,nstsmx,nstsr,ntf1,ntf1mx,ntf1t,ntf2,ntf2mx,ntf2t,nttyo,nxrn1,nxrn2,nxt,nxti,nxtimx,nxtmax,omega,pe,perc,ppmwe,qrho,qxknph,rho,rhoc,rhowc,sidrph,sidrsp,sigmam,sigzi,tdsglw,tdspkc,tdspkg,tdspl,tdsplc,tempc,tf1,tf2,tolspf,uelem,ugexj,ugexmo,uphase,uspec,uxtype,vosol,wfac,wfh2o,wftds,wkgwi,woh2o,wosol,wotds,xbar,xbarlg,xbarw,xbrwlg,xgexw,xlke,xlksd,zchar,zchcu6,zchsq2)
    !! This subroutine writes a description of the computed aqueous
    !! solution model on the output file.
    !! This subroutine is called by:
    !!   EQ3NR/eq3nr.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: iapxmx
    integer :: ibpxmx
    integer :: ietmax
    integer :: iktmax
    integer :: jetmax
    integer :: jsomax
    integer :: ketmax
    integer :: natmax
    integer :: nbtmax
    integer :: nctmax
    integer :: ndrsmx
    integer :: nessmx
    integer :: netmax
    integer :: ngtmax
    integer :: nopgmx
    integer :: noprmx
    integer :: noptmx
    integer :: nptmax
    integer :: nstmax
    integer :: nstsmx
    integer :: ntf1mx
    integer :: ntf2mx
    integer :: nxtimx
    integer :: nxtmax

    integer :: noutpt
    integer :: nttyo

    integer :: iopg(nopgmx)
    integer :: iopr(noprmx)
    integer :: iopt(noptmx)
    integer :: jcsort(nstmax)
    integer :: jern1(jetmax,netmax)
    integer :: jern2(jetmax,netmax)
    integer :: jgsort(ngtmax)
    integer :: jflag(nstmax)
    integer :: jflagd(nstmax)
    integer :: jflgi(nbtmax)
    integer :: jgext(netmax)
    integer :: jpflag(nptmax)
    integer :: jsflag(nstmax)
    integer :: jsol(nxtmax)
    integer :: kern1(netmax)
    integer :: kern2(netmax)
    integer :: kgexsa(ketmax,netmax)
    integer :: nbasp(nbtmax)
    integer :: nbaspd(nbtmax)
    integer :: ncmpr(2,nptmax)
    integer :: ndrsd(ndrsmx)
    integer :: ndrsrd(2,nstmax)
    integer :: ness(nessmx)
    integer :: nessr(2,nstmax)
    integer :: ngexpi(netmax)
    integer :: ngexsa(ietmax,jetmax,netmax)
    integer :: ngext(jetmax,netmax)
    integer :: npnxp(nxtimx)
    integer :: nsts(nstsmx)
    integer :: nstsr(2,nstmax)
    integer :: ntf1(ntf1mx)
    integer :: ntf2(ntf2mx)

    integer :: iebal
    integer :: iern1
    integer :: iern2
    integer :: igas
    integer :: ilrn1
    integer :: ilrn2
    integer :: imrn1
    integer :: imrn2
    integer :: ixrn1
    integer :: ixrn2
    integer :: jfleba
    integer :: narn1
    integer :: narn2
    integer :: nbt
    integer :: nchlor
    integer :: nct
    integer :: nelect
    integer :: net
    integer :: neti
    integer :: ngrn1
    integer :: ngrn2
    integer :: ngt
    integer :: nhydr
    integer :: nhydx
    integer :: no2gaq
    integer :: npt
    integer :: nrdxsp
    integer :: nst
    integer :: ntf1t
    integer :: ntf2t
    integer :: nxrn1
    integer :: nxrn2
    integer :: nxt
    integer :: nxti

    logical :: qxknph(nptmax)

    logical :: qrho

    character(len=48) :: uspec(nstmax)
    character(len=32) :: uxtype(jsomax)
    character(len=24) :: ugexmo(netmax)
    character(len=24) :: uphase(nptmax)
    character(len=8) :: uelem(nctmax)
    character(len=8) :: ugexj(jetmax,netmax)

    real(kind=8) :: acflg(nstmax)
    real(kind=8) :: act(nstmax)
    real(kind=8) :: actlg(nstmax)
    real(kind=8) :: ahrc(nbtmax)
    real(kind=8) :: atwt(nctmax)
    real(kind=8) :: a3bars(natmax)
    real(kind=8) :: affpd(nptmax)
    real(kind=8) :: affsd(nstmax)
    real(kind=8) :: apx(iapxmx,nxtmax)
    real(kind=8) :: bpx(ibpxmx,nxtmax)
    real(kind=8) :: cdrsd(ndrsmx)
    real(kind=8) :: cegexs(ietmax,jetmax,netmax)
    real(kind=8) :: cess(nessmx)
    real(kind=8) :: conc(nstmax)
    real(kind=8) :: conclg(nstmax)
    real(kind=8) :: coval(nbtmax)
    real(kind=8) :: csts(nstsmx)
    real(kind=8) :: ctb(nbtmax)
    real(kind=8) :: cteaq(nctmax)
    real(kind=8) :: egexjc(jetmax,netmax)
    real(kind=8) :: egexjf(jetmax,netmax)
    real(kind=8) :: egexpa(netmax)
    real(kind=8) :: egexpc(netmax)
    real(kind=8) :: egexs(ietmax,jetmax,netmax)
    real(kind=8) :: egexw(ketmax,netmax)
    real(kind=8) :: ehrc(nbtmax)
    real(kind=8) :: fo2lrc(nbtmax)
    real(kind=8) :: fugac(ngtmax)
    real(kind=8) :: fugalg(ngtmax)

    real(kind=8) :: moph(nptmax)
    real(kind=8) :: mosp(nstmax)
    real(kind=8) :: mte(nctmax)
    real(kind=8) :: mteaq(nctmax)
    real(kind=8) :: mwtsp(nstmax)
    real(kind=8) :: perc(nbtmax)
    real(kind=8) :: ppmwe(nctmax)
    real(kind=8) :: sidrph(nptmax)
    real(kind=8) :: sidrsp(nstmax)
    real(kind=8) :: tf1(ntf1mx)
    real(kind=8) :: tf2(ntf2mx)
    real(kind=8) :: wfac(iktmax,nxtmax)
    real(kind=8) :: xbar(nstmax)
    real(kind=8) :: xbarlg(nstmax)
    real(kind=8) :: xgexw(ketmax,netmax)
    real(kind=8) :: xlksd(nstmax)
    real(kind=8) :: zchar(nstmax)
    real(kind=8) :: zchcu6(nstmax)
    real(kind=8) :: zchsq2(nstmax)

    real(kind=8) :: abar
    real(kind=8) :: adh
    real(kind=8) :: afcnst
    real(kind=8) :: alki
    real(kind=8) :: a3bar
    real(kind=8) :: ehfac
    real(kind=8) :: eps100
    real(kind=8) :: eh
    real(kind=8) :: farad
    real(kind=8) :: fje
    real(kind=8) :: fo2
    real(kind=8) :: fo2lg
    real(kind=8) :: fxi
    real(kind=8) :: omega
    real(kind=8) :: pe
    real(kind=8) :: rho
    real(kind=8) :: sigmam
    real(kind=8) :: sigzi
    real(kind=8) :: tdspkg
    real(kind=8) :: tdspl
    real(kind=8) :: tempc
    real(kind=8) :: tolspf
    real(kind=8) :: vosol
    real(kind=8) :: wfh2o
    real(kind=8) :: wftds
    real(kind=8) :: wkgwi
    real(kind=8) :: woh2o
    real(kind=8) :: wosol
    real(kind=8) :: wotds
    real(kind=8) :: xbarw
    real(kind=8) :: xbrwlg
    real(kind=8) :: xlke

    real(kind=8) :: mlmrra
    real(kind=8) :: mrmlra
    real(kind=8) :: rhoc
    real(kind=8) :: rhowc
    real(kind=8) :: tdsglw
    real(kind=8) :: tdspkc
    real(kind=8) :: tdsplc

    ! Local variable declarations.
    integer :: ier
    integer :: j2
    integer :: nb
    integer :: ncount
    integer :: nei
    integer :: np
    integer :: nr1
    integer :: nr2
    integer :: ns
    integer :: nss
    integer :: ns1
    integer :: ns2
    integer :: nt
    integer :: nx
    integer :: nxi

    integer :: ilnobl

    logical :: qphcl
    logical :: qredox

    character(len=8) :: uadj
    character(len=8) :: ufinal
    character(len=8) :: uinput
    character(len=8) :: ux1
    character(len=8) :: ux2
    character(len=8) :: ux3

    real(kind=8) :: acfw
    real(kind=8) :: acfwlg
    real(kind=8) :: alk1
    real(kind=8) :: alk2
    real(kind=8) :: actw
    real(kind=8) :: actwlg
    real(kind=8) :: ah
    real(kind=8) :: ahmes
    real(kind=8) :: ahnbs
    real(kind=8) :: axc
    real(kind=8) :: axd
    real(kind=8) :: axx
    real(kind=8) :: cte
    real(kind=8) :: ctebal
    real(kind=8) :: cx
    real(kind=8) :: cxv
    real(kind=8) :: cxw
    real(kind=8) :: ehmes
    real(kind=8) :: ehnbs
    real(kind=8) :: fjest
    real(kind=8) :: fxist
    real(kind=8) :: msigzm
    real(kind=8) :: osc
    real(kind=8) :: oscst
    real(kind=8) :: osfac
    real(kind=8) :: pch
    real(kind=8) :: pemes
    real(kind=8) :: penbs
    real(kind=8) :: ph
    real(kind=8) :: phc
    real(kind=8) :: phcl
    real(kind=8) :: phd
    real(kind=8) :: phnbs
    real(kind=8) :: phmes
    real(kind=8) :: phx
    real(kind=8) :: ppmv
    real(kind=8) :: ppmw
    real(kind=8) :: pxc
    real(kind=8) :: pxd
    real(kind=8) :: pxx
    real(kind=8) :: pzmean
    real(kind=8) :: pztot
    real(kind=8) :: sanion
    real(kind=8) :: sigmst
    real(kind=8) :: sigza
    real(kind=8) :: sigzc
    real(kind=8) :: sigzm
    real(kind=8) :: sx

    real(kind=8) :: coefst
    real(kind=8) :: texp
    real(kind=8) :: tlg

    data uinput /'Input   '/,ufinal /'Final   '/,uadj   /'Adj     '/

    write (noutpt,1000)
1000 format(/' - - - - - - - - - - - - - - - - - - - - - - - - - - -',' - - - - - - - - - - - -',/)

    ! XX   Begin jsol values with 0 or 1? Need to work something out with
    ! XX   the data base group as to which. Wait until solid solutions are
    ! XX   redone.
    do nx = 1,nxt
        if (jsol(nx) .eq. 0) then
            jsol(nx) = 1
        end if
    end do

    ! Print a table of the elemental composition of the aqueous
    ! solution.
    call prteca(cteaq,mrmlra,nct,nctmax,noutpt,ppmwe,qrho,rho,uelem)

    ! Compute and save the calculated total concentration of the
    ! species adjusted for electrical balance.
    if (iebal .gt. 0) then
        ctebal = 0.
        nb = iebal

        do nss = narn1,narn2
            ns = jcsort(nss)

            if (jsflag(ns) .le. 0) then
                cx = coefst(csts,nsts,nstsmx,nstsr,nb,ns,nstmax)
                ctebal = ctebal + cx*conc(ns)
            end if
        end do
    end if

    ! Print a table of the numerical composition of the aqueous
    ! solution.
    call prtnca(ctb,jflag,jsflag,mrmlra,mwtsp,narn1,narn2,nbasp,nbaspd,nbt,nbtmax,noutpt,nstmax,qrho,rho,uspec,wfh2o)

    ! Compute and print a table of the sensible composition of the
    ! aqueous solution.
    call prtsca(ctb,jflag,jsflag,mrmlra,mwtsp,narn1,narn2,nbasp,nbaspd,nbt,nbtmax,nelect,nhydr,nhydx,noutpt,no2gaq,nstmax,qrho,rho,uspec,wfh2o)

    ! Compute and print various aqueous solution parameters.
    actwlg = actlg(narn1)
    actw = texp(actwlg)

    acfwlg = acflg(narn1)
    acfw = texp(acfwlg)

    ! Compute the sum of the stoichiometric molalities.
    sigmst = 0.

    do nb = 1,nbt
        ns = nbaspd(nb)

        ! The 'abs' comes into play below only when ns = nhydr or nhydx
        ! and the value is negative (e.g., a negative ctb for H+ is
        ! interpreted as a positive value for OH-).
        if (ns.ne.no2gaq .and. ns.ne.narn1) then
            sigmst = sigmst + abs(ctb(nb))
        end if
    end do

    ! Compute the true and stoichiometric osmotic coefficients.
    osfac = -omega * log(actw)
    osc = 0.

    if (sigmam .gt. eps100) then
        osc = osfac/sigmam
    end if

    oscst = 0.

    if (sigmst .gt. eps100) then
        oscst = osfac/sigmst
    end if

    ! Compute the stoichiometric ionic strength.
    ! Calling sequence substitution:
    !   fxist for fxistc
    call cfxist(ctb,fxist,nbaspd,nbt,nbtmax,nstmax,zchsq2)

    ! Compute the stoichiometric ionic asymmetry.
    ! Calling sequence substitution:
    !   fjest for fjestc
    call cfjest(ctb,fjest,nbaspd,nbt,nbtmax,nstmax,zchcu6)

    call prtvpa(abar,acfw,acfwlg,actw,actwlg,a3bar,fje,fjest,fo2,fo2lg,fxi,fxist,iopg,mlmrra,mrmlra,nopgmx,noutpt,osc,oscst,qrho,rhoc,rhowc,sigmam,sigmst,tdsglw,tdspkc,tdsplc,vosol,wfh2o,wftds,woh2o,wosol,wotds,xbarw,xbrwlg)

    ! Compute and print pH, Eh, and pe-, all with reference to
    ! appropriate pH scales. Also compute and print the pHCl.
    qredox = .true.

    call gpheh(acflg,actlg,actwlg,adh,ah,ahmes,ahnbs,conc,eh,ehfac,ehmes,ehnbs,farad,fo2lg,fxi,iopg,mrmlra,nchlor,nhydr,nopgmx,noutpt,nstmax,nttyo,pch,pe,pemes,penbs,ph,phcl,phmes,phnbs,qphcl,qredox,qrho,xlke)

    call prpheh(ah,ahmes,ahnbs,eh,ehmes,ehnbs,iopg,nopgmx,noutpt,pch,pe,pemes,penbs,ph,phcl,phmes,phnbs,qphcl,qredox,qrho)

    ! Print a table of computed alkalinity parameters.
    ! Compute the HCO3-CO3-OH total alkalinity.
    ! Calling sequence substitutions:
    !   alk1 for alkc
    !   ntf1 for ntfx
    !   ntf1mx for ntfxmx
    !   ntf1t for ntfxt
    !   tf1 for tfx
    call calk(alk1,conc,nstmax,ntf1,ntf1mx,ntf1t,tf1)

    ! Compute the extended total alkalinity.
    ! Calling sequence substitutions:
    !   alk2 for alkc
    !   ntf2 for ntfx
    !   ntf2mx for ntfxmx
    !   ntf2t for ntfxt
    !   tf2 for tfx
    call calk(alk2,conc,nstmax,ntf2,ntf2mx,ntf2t,tf2)

    call prtalk(alki,alk1,alk2,mrmlra,noutpt,ntf1t,ntf2t,qrho,rho,tempc,wfh2o)

    ! Calculate and print the electrical balance and the cation and
    ! anion contributions.
    call gszm(conc,jcsort,narn1,narn2,nstmax,sigza,sigzc,sigzi,sigzm,zchar)

    msigzm = 0.5*sigzm
    sanion = -sigza

    write (noutpt,1250)
1250 format(/11x,'--- Electrical Balance Totals ---',// 34x,'eq/kg.H2O',/)

    write (noutpt,1260) sigzc,sanion,sigzm,msigzm,sigzi
1260 format(8x,'Sigma(mz) cations= ',3x,1pe17.10,/9x,'Sigma(mz) anions= ',3x,1pe17.10,/13x,'Total charge= ',3x,1pe17.10,/14x,'Mean charge= ',3x,1pe17.10,/9x,'Charge imbalance= ',3x,1pe17.10)

    sx = 100.*sigzi
    pztot = sx/sigzm
    pzmean = sx/msigzm
    write (noutpt,1270) pztot,pzmean
1270 format(//9x,'The electrical imbalance is:',//11x,f9.4,' per cent of the total charge',/11x,f9.4,' per cent of the mean charge',/)

    ! Print results of electrical balancing.
    if (iebal .gt. 0) then
        ns = nbaspd(iebal)
        j2 = ilnobl(uspec(ns)(1:24))
        write (noutpt,1280) uspec(ns)(1:j2)
1280 format(/6x,'--- Electrical Balancing on ',a,' ---')

        ux1 = uinput
        ux2 = ufinal
        ux3 = uadj

        if (jfleba .eq. 0) then
            write (noutpt,1290)
1290 format(/17x,'mg/L',10x,'mg/kg.sol',7x,'Molality',/)

            cte = coval(iebal)
            cxw = 1000.*cte*mwtsp(ns)*wfh2o
            cxv = cxw*rho
            write (noutpt,1300) ux1,cxv,cxw,cte
1300 format(5x,a6,f13.4,3x,f13.4,3x,1pe17.10)

            cte = ctebal
            cxw = 1000.*cte*mwtsp(ns)*wfh2o
            cxv = cxw*rho
            write (noutpt,1300) ux2,cxv,cxw,cte
            cte = ctebal - coval(iebal)
            cxw = 1000.*cte*mwtsp(ns)*wfh2o
            cxv = cxw*rho
            write (noutpt,1300) ux3,cxv,cxw,cte
        else if (jfleba .eq. 16) then
            write (noutpt,1310)
1310 format(/3x,'        Log Activity',/)

            axx = coval(iebal)
            axc = actlg(ns)
            write (noutpt,1320) ux1,axx
            write (noutpt,1320) ux2,axc
1320 format(5x,a6,f13.4)

            axd = axc - axx
            write (noutpt,1320) ux3,axd
        else if (jfleba .eq. 19) then
            write (noutpt,1312)
1312 format(/3x,'            pX',/)

            pxx = coval(iebal)
            pxc = -actlg(ns)
            write (noutpt,1320) ux1,pxx
            write (noutpt,1320) ux2,pxc
            pxd = pxc - pxx
            write (noutpt,1320) ux3,pxd
        else if (jfleba .eq. 20) then
            write (noutpt,1314)
1314 format(/3x,'            pH',/)

            phx = coval(iebal)
            phc = -actlg(ns)
            write (noutpt,1320) ux1,phx
            write (noutpt,1320) ux2,phc
            phd = phc - phx
            write (noutpt,1320) ux3,phd
        else if (jfleba .eq. 21) then
            write (noutpt,1315)
1315 format(/3x,'            pHCl',/)

            phx = coval(iebal)
            phc = -actlg(nhydr)-actlg(nchlor)
            write (noutpt,1320) ux1,phx
            write (noutpt,1320) ux2,phc
            phd = phc - phx
            write (noutpt,1320) ux3,phd
        else if (jfleba .eq. 22) then
            write (noutpt,1316)
1316 format(/3x,'            pmH',/)

            phx = coval(iebal)
            phc = -conclg(ns)
            write (noutpt,1320) ux1,phx
            write (noutpt,1320) ux2,phc
            phd = phc - phx
            write (noutpt,1320) ux3,phd
        else if (jfleba .eq. 23) then
            write (noutpt,1318)
1318 format(/3x,'            pxH',/)

            pxx = coval(iebal)
            pxc = -conclg(ns)
            write (noutpt,1320) ux1,pxx
            write (noutpt,1320) ux2,pxc
            pxd = pxc - pxx
            write (noutpt,1320) ux3,pxd
        end if

        write (noutpt,1322)
1322 format(1x)
    end if

    if (iopr(5) .gt. 0) then
        ! Compute and print activity ratios of aqueous species.
        call prtacr(actlg,iopr,jsflag,nbaspd,nbt,nbtmax,nelect,nhydr,noprmx,no2gaq,noutpt,nstmax,uspec,zchar)
    end if

    ! Print the aqueous species distribution.
    call prtaqs(acflg,actlg,conc,conclg,iopr,jcsort,narn1,narn2,noprmx,noutpt,nstmax,uspec)

    if (iopr(9) .gt. 0) then
        ! Compute and print the mean ionic activities and activity
        ! coefficients.
        call prtmip(acflg,actlg,conclg,ctb,nbaspd,nbt,nbtmax,nelect,nhydr,nhydx,noutpt,nstmax,uspec,zchar)
    end if

    ! Print the contributions of aqueous species to each mass balance.
    call prtpct(conc,csts,ctb,iopr,jcsort,jflag,narn1,narn2,nbasp,nbaspd,nbt,nbtmax,nelect,nhydr,noprmx,no2gaq,noutpt,nstmax,nsts,nstsmx,nstsr,uspec)

    ! Compute and print the state of redox reactions that are not
    ! constrained to be at equilibrium.
    call cdardx(actlg,actwlg,ah,ahrc,cdrsd,eh,ehfac,ehrc,farad,fo2lg,fo2lrc,jsflag,mosp,nbasp,nbaspd,nbt,nbtmax,ndrsd,ndrsmx,ndrsrd,no2gaq,nstmax,pe,perc,ph,xlke,xlksd)

    call prtrdx(ah,ahrc,cdrsd,eh,ehrc,fo2lg,fo2lrc,jflgi,jsflag,narn1,nbasp,nbaspd,nbt,nbtmax,ndrsd,ndrsmx,ndrsrd,nelect,nhydr,no2gaq,noutpt,nstmax,pe,perc,uspec)

    ! Compute the saturation index for the reaction associated to
    ! each species. Do this for the reactions as they were written
    ! on the data file.
    if (nxrn1 .gt. 0) then
        do ns = nxrn1,nxrn2
            if (actlg(ns) .le. -99999.) then
                actlg(ns) = 0.0
            end if
        end do
    end if

    ! Calculate affinities and saturation indices using the 'd' set
    ! of reactions.
    call gaffsd(actlg,afcnst,affpd,affsd,cdrsd,jflagd,jpflag,ncmpr,ndrsd,ndrsmx,ndrsrd,npt,nptmax,nst,nstmax,qxknph,sidrph,sidrsp,uphase,uspec,xbar,xlksd)

    ! Compute and print saturation index and affinity tables for
    ! reactions in the aqueous phase not constrained to be at
    ! equilibrium.
    call prtsia(affsd,jflagd,jflgi,jsflag,narn1,narn2,nbasp,nbaspd,nbt,nbtmax,ndrsd,ndrsmx,ndrsrd,nhydr,noutpt,nrdxsp,nstmax,sidrsp,uspec)

    if (iopr(7) .ge. 0) then
        ! Print saturation index and affinity tables for the
        ! various non-aqueous phases.
        call prtsat(affpd,iern1,iern2,ilrn1,ilrn2,imrn1,imrn2,iopr,iopt,ixrn1,ixrn2,jpflag,noutpt,noprmx,noptmx,nptmax,sidrph,tolspf,uphase)
    end if

    if (iopt(4).ge.1 .and. nxti.gt.0) then
        ! Print affinities and saturation indices for solid solutions
        ! for which compositions were entered on the input file.
        write (noutpt,1600)
1600 format(//16x,'--- Saturation States of Input Solid Solutions',' ---',/)

        ncount = 0

        do nxi = 1,nxti
            np = npnxp(nxi)
            ncount = ncount + 1
            call prtsso(acflg,actlg,affpd,affsd,ixrn1,ixrn2,jsol,jsomax,ncmpr,noutpt,np,nptmax,nstmax,nxtmax,sidrph,sidrsp,tolspf,uspec,uphase,uxtype,xbar,xbarlg)
        end do

        if (ncount .le. 0) then
            write (noutpt,1710)
        end if

1710 format(1x,'None')
    end if

    if (iopt(4).ge.1 .and. ixrn2.ge.ixrn1) then
        ! Compute hypothetical saturation states for solid solutions.
        write (noutpt,1720)
1720 format(/16x,'--- Saturation States of Hypothetical Solid',' Solutions ---',/)

        ncount = 0

        do np = ixrn1,ixrn2
            nr1 = ncmpr(1,np)
            nr2 = ncmpr(2,np)
            nt = nr2 - nr1 + 1

            if (nt .gt. 1) then
                ier = 0
                ncount = ncount + 1

                ! Calling sequence substitutions:
                !   affpd for affp
                !   affsd for affs
                !   cdrsd for cdrs
                !   ndrsd for ndrs
                !   ndrsrd for ndrsr
                !   xlksd for xlks
                call hpsat(acflg,act,actlg,afcnst,affpd,affsd,apx,bpx,cdrsd,eps100,iapxmx,ibpxmx,ier,iktmax,ixrn1,jflag,jpflag,jsflag,jsol,ncmpr,ndrsd,ndrsmx,ndrsrd,noutpt,np,nptmax,nstmax,nttyo,nxrn1,nxrn2,nxtmax,sidrsp,sidrph,uphase,uspec,wfac,xbar,xbarlg,xlksd)

                ! Check to see if the hypothetical affinity calculation
                ! converged.
                if (ier .le. 0) then
                    qxknph(np) = .true.
                else
                    qxknph(np) = .false.
                    j2 = ilnobl(uphase(np))
                    write (noutpt,1730) uphase(np)(1:j2)
                    write (nttyo,1730) uphase(np)(1:j2)
1730 format(/' * Warning - (EQ3NR/scripz) The hypothetical',/7x,'affinity calculation failed for ',a,'.')

                    go to 310
                end if

                call prtsso(acflg,actlg,affpd,affsd,ixrn1,ixrn2,jsol,jsomax,ncmpr,noutpt,np,nptmax,nstmax,nxtmax,sidrph,sidrsp,tolspf,uspec,uphase,uxtype,xbar,xbarlg)
            end if

310 continue
        end do

        if (ncount .le. 0) then
            write (noutpt,1740)
        end if

1740 format(1x,'None')
    end if

    ! XX   Temporarily suppress these tables by setting neti = 0.
    ! XX   Maybe reactivate them in the future.
    neti = 0

    ! XX
    if (neti .gt. 0) then
        ! Print affinities and saturation indices for generic ion
        ! exchangers for which compositions were entered on the input
        ! file.
        write (noutpt,1750)
1750 format(//11x,'--- Saturation States of Input Generic Ion',' Exchangers ---',/)

        ncount = 0

        do nei = 1,neti
            np = ngexpi(nei)
            ncount = ncount + 1

            call prtgex(acflg,actlg,affpd,affsd,cegexs,conc,egexjc,egexjf,egexpa,egexpc,egexs,egexw,iern1,iern2,ietmax,jern1,jern2,jetmax,jgext,kern1,kern2,ketmax,kgexsa,moph,mosp,netmax,ngexsa,ngext,noutpt,np,nptmax,nstmax,sidrph,sidrsp,tolspf,ugexj,ugexmo,uspec,uphase,xbar,xbarlg,xgexw,wkgwi)
        end do

        if (ncount .le. 0) then
            write (noutpt,1760)
        end if

1760 format(1x,'None')
    end if

    if (net .gt. 0) then
        ! Print compositions of generic ion exchangers.
        write (noutpt,1770)
1770 format(/11x,'--- Generic Ion Exchangers ---',/)

        ncount = 0

        do np = iern1,iern2
            nr1 = ncmpr(1,np)
            nr2 = ncmpr(2,np)
            nt = nr2 - nr1 + 1

            if (nt .gt. 1) then
                ier = 0
                ncount = ncount + 1

                call prtgex(acflg,actlg,affpd,affsd,cegexs,conc,egexjc,egexjf,egexpa,egexpc,egexs,egexw,iern1,iern2,ietmax,jern1,jern2,jetmax,jgext,kern1,kern2,ketmax,kgexsa,moph,mosp,netmax,ngexsa,ngext,noutpt,np,nptmax,nstmax,sidrph,sidrsp,tolspf,ugexj,ugexmo,uspec,uphase,xbar,xbarlg,xgexw,wkgwi)
            end if
        end do

        if (ncount .le. 0) then
            write (noutpt,1790)
        end if

1790 format(1x,'None')
    end if

    if (iopr(8) .ge. 0) then
        ! Print a table of equilibrium fugacities.
        call prtfug(jgsort,fugac,fugalg,jsflag,ngrn1,ngt,ngtmax,noutpt,nstmax,uspec)
    end if

    write (noutpt,2000)
2000 format(/' - - - - - - - - - - - - - - - - - - - - - - - - - - -',' - - - - - - - - - - - -',/)
end subroutine scripx