subroutine cdappl(acflg,acfw,acfwlg,actlg,actw,actwlg,adwipp,afcnst,affpd,affsd,ah,ahrc,alk,alk1,alk2,alki,atwt,bdwipp,cdrsd,cess,conc,csts,ctb,cteaq,dvoso,dwoso,eh,ehfac,ehrc,eps100,farad,fdpe0,fdse0,fjest,fo2lg,fo2lrc,fxist,iaqsln,iemop0,iemos0,iern1,iern2,ifrn1,ifrn2,ilrn1,ilrn2,imrn1,imrn2,iopt,ixrn1,ixrn2,jcode,jcsort,jern1,jflag,jflagd,jgext,jpflag,jsflag,jssort,modr,moph,mophg,mophj,mopht,morr,mosp,mospg,mospj,mospt,mprph,mprsp,mrgers,mrmlra,mtb,mtbaq,mte,mteaq,mwtges,mwtrc,mwtsp,narn1,narn2,nat,nbasp,nbaspd,nbt,nchlor,ncmpe0,ncmpr,nct,ndrsd,ndrsrd,nelect,nern1,nern2,nert,ness,nessr,net,nfrn1,nfrn2,ngext,ngrn1,ngrn2,ngt,nhydr,nlrn1,nlrn2,nlt,nmrn1,nmrn2,nmrt,nmt,no2gaq,npchk,npet,npet0,npt,npts,nrct,nrndex,nst,nsts,nstsr,ntf1,ntf1t,ntf2,ntf2t,nxridx,nxrn1,nxrn2,nxrt,nxt,osc,oscst,omega,pe,perc,ph,phmes,ppmwb,ppmwe,qriinf,qxknph,rreacn,rreac1,rxbar,sfcar,sidrsp,sidrph,sigmam,sigmst,tdays,tempc,tf1,tf2,thours,time1,tmins,tyears,uphase,uspec,vodrt,voph,vophg,vophj,vopht,vosoct,vosp,vospg,vospj,vospt,vosp0,vreac,wfh2o,wkgwi,wodr,wodrt,woph,wophg,wophj,wopht,worr,worrt,wosoct,wosp,wospg,wospj,wospt,xbar,xlke,xlksd,zchcu6,zchsq2)
    !! This subroutine computes various data pertaining to the modeled
    !! system for printing and plotting at the current point of
    !! reaction progress (xi1).
    !! This subroutine is called by:
    !!   EQ6/path.f
    !! Principal input:
    !! Principal output:
    implicit none

    include 'eqlib/eqlpar.h'
    include 'eqlib/eqldv.h'

    ! Calling sequence variable declarations.
    integer :: iemop0(npetmx)
    integer :: iemos0(nsetmx)
    integer :: iopt(noptmx)
    integer :: jcode(nrctmx)
    integer :: jcsort(nstmax)
    integer :: jern1(jetmax,netmax)
    integer :: jflag(nstmax)
    integer :: jflagd(nstmax)
    integer :: jgext(netmax)
    integer :: jpflag(nptmax)
    integer :: jsflag(nstmax)
    integer :: jssort(nstmax)
    integer :: nbasp(nbtmax)
    integer :: nbaspd(nbtmax)
    integer :: ncmpe0(2,npetmx)
    integer :: ncmpr(2,nptmax)
    integer :: ndrsd(ndrsmx)
    integer :: ndrsrd(2,nstmax)
    integer :: ness(nessmx)
    integer :: nessr(2,nstmax)
    integer :: ngext(jetmax,netmax)
    integer :: npchk(nptmax)
    integer :: nrndex(nrctmx)
    integer :: nsts(nstsmx)
    integer :: nstsr(2,nstmax)
    integer :: ntf1(ntf1mx)
    integer :: ntf2(ntf2mx)
    integer :: nxridx(nrctmx)

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

    integer :: nat
    integer :: nbt
    integer :: nct
    integer :: net
    integer :: nlt
    integer :: nmt
    integer :: nst
    integer :: npt
    integer :: ngt
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

    integer :: iaqsln
    integer :: nchlor
    integer :: nelect
    integer :: nert
    integer :: nhydr
    integer :: nmrt
    integer :: no2gaq
    integer :: npet
    integer :: npet0
    integer :: npts
    integer :: nrct
    integer :: ntf1t
    integer :: ntf2t
    integer :: nxrt

    logical :: qxknph(nptmax)

    logical :: qriinf

    character(len=48) :: uspec(nstmax)
    character(len=24) :: uphase(nptmax)

    real(kind=8) :: acflg(nstmax)
    real(kind=8) :: actlg(nstmax)
    real(kind=8) :: affpd(nptmax)
    real(kind=8) :: affsd(nstmax)
    real(kind=8) :: ahrc(nbtmax)
    real(kind=8) :: atwt(nctmax)
    real(kind=8) :: cdrsd(ndrsmx)
    real(kind=8) :: cess(nessmx)
    real(kind=8) :: conc(nstmax)
    real(kind=8) :: csts(nstsmx)
    real(kind=8) :: ctb(nbtmax)
    real(kind=8) :: cteaq(nctmax)
    real(kind=8) :: ehrc(nbtmax)
    real(kind=8) :: fdpe0(nordmx,npetmx)
    real(kind=8) :: fdse0(nordmx,nsetmx)
    real(kind=8) :: fo2lrc(nbtmax)

    real(kind=8) :: modr(nrctmx)
    real(kind=8) :: moph(nptmax)
    real(kind=8) :: mophg(nptmax)
    real(kind=8) :: mophj(nptmax)
    real(kind=8) :: mopht(nptmax)
    real(kind=8) :: morr(nrctmx)
    real(kind=8) :: mosp(nstmax)
    real(kind=8) :: mospg(nstmax)
    real(kind=8) :: mospj(nstmax)
    real(kind=8) :: mospt(nstmax)
    real(kind=8) :: mprph(nptmax)
    real(kind=8) :: mprsp(nstmax)
    real(kind=8) :: mrgers(ietmax,jetmax,nertmx)
    real(kind=8) :: mtb(nbtmax)
    real(kind=8) :: mtbaq(nbtmax)
    real(kind=8) :: mte(nctmax)
    real(kind=8) :: mteaq(nctmax)
    real(kind=8) :: mwtges(netmax)
    real(kind=8) :: mwtrc(nrctmx)
    real(kind=8) :: mwtsp(nstmax)

    real(kind=8) :: perc(nbtmax)
    real(kind=8) :: ppmwb(nbtmax)
    real(kind=8) :: ppmwe(nctmax)
    real(kind=8) :: rreacn(nrctmx)
    real(kind=8) :: rreac1(nrctmx)
    real(kind=8) :: rxbar(iktmax,nxrtmx)
    real(kind=8) :: sfcar(nrctmx)
    real(kind=8) :: sidrph(nptmax)
    real(kind=8) :: sidrsp(nstmax)
    real(kind=8) :: tf1(ntf1mx)
    real(kind=8) :: tf2(ntf2mx)
    real(kind=8) :: voph(nptmax)
    real(kind=8) :: vophg(nptmax)
    real(kind=8) :: vophj(nptmax)
    real(kind=8) :: vopht(nptmax)
    real(kind=8) :: vosp(nstmax)
    real(kind=8) :: vospg(nstmax)
    real(kind=8) :: vospj(nstmax)
    real(kind=8) :: vospt(nstmax)
    real(kind=8) :: vosp0(nstmax)
    real(kind=8) :: vreac(nrctmx)
    real(kind=8) :: wodr(nrctmx)
    real(kind=8) :: woph(nptmax)
    real(kind=8) :: wophg(nptmax)
    real(kind=8) :: wophj(nptmax)
    real(kind=8) :: wopht(nptmax)
    real(kind=8) :: worr(nrctmx)
    real(kind=8) :: wosp(nstmax)
    real(kind=8) :: wospg(nstmax)
    real(kind=8) :: wospj(nstmax)
    real(kind=8) :: wospt(nstmax)
    real(kind=8) :: xbar(nstmax)
    real(kind=8) :: xlksd(nstmax)
    real(kind=8) :: zchcu6(nstmax)
    real(kind=8) :: zchsq2(nstmax)

    real(kind=8) :: acfw
    real(kind=8) :: acfwlg
    real(kind=8) :: actw
    real(kind=8) :: actwlg
    real(kind=8) :: adwipp
    real(kind=8) :: afcnst
    real(kind=8) :: ah
    real(kind=8) :: alk
    real(kind=8) :: alki
    real(kind=8) :: alk1
    real(kind=8) :: alk2
    real(kind=8) :: bdwipp
    real(kind=8) :: dvoso
    real(kind=8) :: dwoso
    real(kind=8) :: eh
    real(kind=8) :: ehfac
    real(kind=8) :: eps100
    real(kind=8) :: farad
    real(kind=8) :: fjest
    real(kind=8) :: fo2lg
    real(kind=8) :: fxist
    real(kind=8) :: mrmlra
    real(kind=8) :: osc
    real(kind=8) :: oscst
    real(kind=8) :: omega
    real(kind=8) :: pe
    real(kind=8) :: ph
    real(kind=8) :: phmes
    real(kind=8) :: sigmam
    real(kind=8) :: sigmst
    real(kind=8) :: tdays
    real(kind=8) :: tempc
    real(kind=8) :: thours
    real(kind=8) :: time1
    real(kind=8) :: tmins
    real(kind=8) :: tyears
    real(kind=8) :: vosoct
    real(kind=8) :: vodrt
    real(kind=8) :: wfh2o
    real(kind=8) :: wkgwi
    real(kind=8) :: wosoct
    real(kind=8) :: wodrt
    real(kind=8) :: worrt
    real(kind=8) :: xlke

    ! Local variable declarations.
    integer :: ie
    integer :: je
    integer :: ik
    integer :: n
    integer :: nb
    integer :: nc
    integer :: ne
    integer :: ner
    integer :: np
    integer :: npe
    integer :: np1
    integer :: ns
    integer :: nse
    integer :: nrc
    integer :: nr1
    integer :: nr2
    integer :: nss
    integer :: nxr

    logical :: qrct

    real(kind=8) :: mophx
    real(kind=8) :: mospx
    real(kind=8) :: osfac
    real(kind=8) :: vophx
    real(kind=8) :: vospx
    real(kind=8) :: wophx
    real(kind=8) :: wospx

    real(kind=8) :: texp

    ! Compute the model time in various units.
    if (iopt(2) .ge. 1) then
        if (qriinf) then
            tmins = 1.e+38
            thours = 1.e+38
            tdays = 1.e+38
            tyears = 1.e+38
        else
            tmins = time1/60.
            thours = tmins/60.
            tdays = thours/24.
            tyears = tdays/365.25
        end if
    end if

    ! Compute some data for reactants and corresponding reactions.
    vodrt = 0.
    wodrt = 0.
    worrt = 0.

    do nrc = 1,nrct
        worr(nrc) = mwtrc(nrc)*morr(nrc)
        wodr(nrc) = mwtrc(nrc)*modr(nrc)
        worrt = worrt + worr(nrc)
        wodrt = wodrt + wodr(nrc)
        vodrt = vodrt + modr(nrc)*vreac(nrc)
    end do

    if (iopt(2) .ge. 1) then
        do nrc = 1,nrct
            rreacn(nrc) = 0.

            if (sfcar(nrc) .gt. 0.) then
                rreacn(nrc) = rreac1(nrc)/sfcar(nrc)
            end if
        end do
    end if

    ! Compute data for a summary of the elemental composition of the
    ! aqueous phase.
    call initaz(mte,nctmax)
    call initaz(mteaq,nctmax)

    do nss = 1,nst
        ns = jssort(nss)

        if (mosp(ns) .gt. 0.) then
            nr1 = nessr(1,ns)
            nr2 = nessr(2,ns)

            do n = nr1,nr2
                nc = ness(n)

                if (nc .gt. 0) then
                    mte(nc) = mte(nc) + cess(n)*mosp(ns)
                end if
            end do
        end if
    end do

    do nss = narn1,narn2
        ns = jcsort(nss)

        if (ns .ne. no2gaq) then
            nr1 = nessr(1,ns)
            nr2 = nessr(2,ns)

            do n = nr1,nr2
                nc = ness(n)

                if (nc .gt. 0) then
                    mteaq(nc) = mteaq(nc) + cess(n)*mosp(ns)
                end if
            end do
        end if
    end do

    ! Compute the elemental composition of the aqueous solution
    ! (molalities and ppm: mg/kg.sol).
    call initaz(cteaq,nctmax)
    call initaz(ppmwe,nctmax)

    do nc = 1,nct
        cteaq(nc) = wkgwi*mteaq(nc)
        ppmwe(nc) = 1000.*cteaq(nc)*atwt(nc)*wfh2o
    end do

    ! Compute the mass and concentration totals in the aqueous phase
    ! in terms of data file basis species.
    call initaz(mtbaq,nbtmax)
    call initaz(ctb,nbtmax)
    call initaz(ppmwb,nbtmax)

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

    do nb = 1,nbt
        ctb(nb) = wkgwi*mtbaq(nb)
        ns = nbaspd(nb)
        ppmwb(nb) = 1000.*ctb(nb)*mwtsp(ns)*wfh2o
    end do

    ! Compute various aqueous phase parameters.
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

    ! Compute alkalinity parameters.
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

    alki = 0.
    alk = alk2

    ! Compute data for redox reactions that are not constrained to be
    ! at equilibrium.
    call cdardx(actlg,actwlg,ah,ahrc,cdrsd,eh,ehfac,ehrc,farad,fo2lg,fo2lrc,jsflag,mosp,nbasp,nbaspd,nbt,nbtmax,ndrsd,ndrsmx,ndrsrd,no2gaq,nstmax,pe,perc,ph,xlke,xlksd)

    ! Calculate affinities and saturation indices using the 'd' set
    ! of reactions.
    call gaffsd(actlg,afcnst,affpd,affsd,cdrsd,jflagd,jpflag,ncmpr,ndrsd,ndrsmx,ndrsrd,npt,nptmax,nst,nstmax,qxknph,sidrph,sidrsp,uphase,uspec,xbar,xlksd)

    vosoct = 0.
    wosoct = 0.

    call initaz(voph,nptmax)
    call initaz(woph,nptmax)
    call initaz(vosp,nstmax)
    call initaz(wosp,nstmax)

    call initaz(mopht,nptmax)
    call initaz(vopht,nptmax)
    call initaz(wopht,nptmax)
    call initaz(mospt,nstmax)
    call initaz(vospt,nstmax)
    call initaz(wospt,nstmax)

    call initaz(mophg,nptmax)
    call initaz(vophg,nptmax)
    call initaz(wophg,nptmax)
    call initaz(mospg,nstmax)
    call initaz(vospg,nstmax)
    call initaz(wospg,nstmax)

    if (iopt(1).eq.0 .or. iopt(1).eq.1) then
        ! Compute data which describe the solid phases in the ES.
        if (npet .gt. 1) then
            do np = 1,npt
                if (np .ne. iaqsln) then
                    if (npchk(np).eq.0 .and. jpflag(np).eq.-1) then
                        nr1 = ncmpr(1,np)
                        nr2 = ncmpr(2,np)
                        vophx = 0.
                        wophx = 0.

                        do ns = nr1,nr2
                            wospx = mosp(ns)*mwtsp(ns)
                            wosp(ns) = wospx
                            wophx = wophx + wospx
                            vospx = mosp(ns)*vosp0(ns)
                            vosp(ns) = vospx
                            vophx = vophx + vospx
                        end do

                        if (np.ge.iern1 .and. np.le.iern2) then
                            ne = np - iern1 + 1
                            wophx = wophx + moph(np)*mwtges(ne)
                        end if

                        woph(np) = wophx
                        voph(np) = vophx
                        wosoct = wosoct + wophx
                        vosoct = vosoct + vophx
                    end if
                end if
            end do
        end if
    end if

    if (iopt(1) .eq. 2) then
        ! Compute the product mineral assemblage which is instantaneous
        ! with respect to Xi (for the fluid-centered flow-through system
        ! only).
        call initaz(mophj,nptmax)
        call initaz(vophj,nptmax)
        call initaz(wophj,nptmax)
        call initaz(mospj,nstmax)
        call initaz(vospj,nstmax)
        call initaz(wospj,nstmax)

        if (npts .gt. 0) then
            ! Use the two-point derivative.
            do npe = 1,npet0
                np = iemop0(npe)
                mophj(np) = fdpe0(1,npe)

                if (mophj(np) .lt. 0.) then
                    mophj(np) = 0.
                end if

                wophx = 0.
                vophx = 0.

                nr1 = ncmpe0(1,npe)
                nr2 = ncmpe0(2,npe)

                do nse = nr1,nr2
                    ns = iemos0(nse)
                    mospx = fdse0(1,nse)
                    mospj(ns) = mospx
                    wospx = mospx*mwtsp(ns)
                    wospj(ns) = wospx
                    wophx = wophx + wospx
                    vospx = mospx*vosp0(ns)
                    vospj(ns) = vospx
                    vophx = vophx + vospx
                end do

                if (np.ge.iern1 .and. np.le.iern2) then
                    ne = np - iern1 + 1
                    wophx = wophx + mophj(np)*mwtges(ne)
                end if

                wophj(np) = wophx
                vophj(np) = vophx
            end do
        end if
    end if

    if (iopt(1) .eq. 2) then
        ! Compute data which describe the solid phases in the system
        ! ES + PRS.
        do np = 1,npt
            if (np .ne. iaqsln) then
                if (jpflag(np) .lt. 2) then
                    mopht(np) = moph(np) + mprph(np)
                    nr1 = ncmpr(1,np)
                    nr2 = ncmpr(2,np)
                    vophx = 0.
                    wophx = 0.

                    do ns = nr1,nr2
                        mospx = mosp(ns) + mprsp(ns)
                        mospt(ns) = mospx
                        wospx = mospx*mwtsp(ns)
                        wospt(ns) = wospx
                        wophx = wophx + wospx
                        vospx = mospx*vosp0(ns)
                        vospt(ns) = vospx
                        vophx = vophx + vospx
                    end do

                    if (np.ge.iern1 .and. np.le.iern2) then
                        ne = np - iern1 + 1
                        wophx = wophx + mopht(np)*mwtges(ne)
                    end if

                    wopht(np) = wophx
                    vopht(np) = vophx
                    wosoct = wosoct + wophx
                    vosoct = vosoct + vophx
                end if
            end if
        end do
    end if

    if ((nmrt + nxrt + nert) .gt. 0) then
        ! Print a grand summary of the solid phases in the system
        ! ES + PRS + reactants.
        do np = imrn1,imrn2
            mophx = moph(np) + mprph(np)
            qrct = .false.

            if (nmrt .gt. 0) then
                do nrc = 1,nrct
                    if (jcode(nrc) .eq. 0) then
                        np1 = nrndex(nrc)

                        if (np1 .eq. np) then
                            mophx = mophx + morr(nrc)
                            qrct = .true.
                            go to 607
                        end if
                    end if
                end do
            end if

607 continue
            if (mophx.ne.0. .or. qrct) then
                mophg(np) = mophx

                ns = ncmpr(1,np)
                mospx = mophx
                mospg(ns) = mospx
                wospx = mospx*mwtsp(ns)
                wospg(ns) = wospx
                vospx = mospx*vosp0(ns)
                vospg(ns) = vospx

                wophx = wospx
                vophx = vospx
                wophg(np) = wophx
                vophg(np) = vophx
            end if
        end do

        do np = ifrn1,ifrn2
            mophx = moph(np) + mprph(np)
            qrct = .false.

            if (nmrt .gt. 0) then
                do nrc = 1,nrct
                    if (jcode(nrc) .eq. 0) then
                        np1 = nrndex(nrc)

                        if (np1 .eq. np) then
                            mophx = mophx + morr(nrc)
                            qrct = .true.
                            go to 609
                        end if
                    end if
                end do
            end if

609 continue
            if (mophx.ne.0. .or. qrct) then
                mophg(np) = mophx

                ns = ncmpr(1,np)
                mospx = mophx
                mospg(ns) = mospx
                wospx = mospx*mwtsp(ns)
                wospg(ns) = wospx
                vospx = mospx*vosp0(ns)
                vospg(ns) = vospx

                wophx = wospx
                vophx = vospx
                wophg(np) = wophx
                vophg(np) = vophx
            end if
        end do

        if (iopt(4) .gt. 0) then
            do np = ixrn1,ixrn2
                if (jpflag(np) .lt. 2) then
                    mophx = moph(np) + mprph(np)
                    nr1 = ncmpr(1,np)
                    nr2 = ncmpr(2,np)
                    qrct = .false.
                    nrc = 0

                    if (nxrt .gt. 0) then
                        do nrc = 1,nrct
                            if (jcode(nrc) .eq. 1) then
                                np1 = nrndex(nrc)

                                if (np1 .eq. np) then
                                    mophx = mophx + morr(nrc)
                                    qrct = .true.
                                    go to 612
                                end if
                            end if
                        end do
                    end if

612 continue
                    if (mophx.ne.0. .or. qrct) then
                        mophg(np) = mophx

                        if (qrct) then
                            nxr = nxridx(nrc)
                        end if

                        ik = 0
                        vophx = 0.
                        wophx = 0.

                        do ns = nr1,nr2
                            mospx = mosp(ns) + mprsp(ns)
                            ik = ik + 1

                            if (qrct) then
                                mospx = mospx + morr(nrc)*rxbar(ik,nxr)
                            end if

                            mospg(ns) = mospx
                            wospx = mospx*mwtsp(ns)
                            wospg(ns) = wospx
                            wophx = wophx + wospx
                            vospx = mospx*vosp0(ns)
                            vospg(ns) = vospx
                            vophx = vophx + vospx
                        end do

                        vophg(np) = vophx
                        wophg(np) = wophx
                    end if
                end if
            end do
        end if

        do np = iern1,iern2
            if (jpflag(np) .lt. 2) then
                mophx = moph(np) + mprph(np)
                qrct = .false.
                nrc = 0

                if (nert .gt. 0) then
                    do nrc = 1,nrct
                        if (jcode(nrc) .eq. 5) then
                            np1 = nrndex(nrc)

                            if (np1 .eq. np) then
                                mophx = mophx + morr(nrc)
                                qrct = .true.
                                go to 614
                            end if
                        end if
                    end do
                end if

614 continue
                if (mophx.ne.0. .or. qrct) then
                    mophg(np) = mophx

                    if (qrct) then
                        ner = nxridx(nrc)
                    end if

                    vophx = 0.
                    wophx = 0.
                    ne = np - iern1 + 1

                    do je = 1,jgext(ne)
                        ns = jern1(je,ne) - 1

                        do ie = 1,ngext(je,ne)
                            ns = ns + 1
                            mospx = mosp(ns) + mprsp(ns)

                            if (qrct) then
                                mospx = mospx + morr(nrc)*mrgers(ie,je,ner)
                            end if

                            mospg(ns) = mospx
                            wospx = mospx*mwtsp(ns)
                            wospg(ns) = wospx
                            wophx = wophx + wospx
                            vospx = mospx*vosp0(ns)
                            vospg(ns) = vospx
                            vophx = vophx + vospx
                        end do
                    end do

                    wophx = wophx + mophg(np)*mwtges(ne)

                    vophg(np) = vophx
                    wophg(np) = wophx
                end if
            end if
        end do
    end if

    dvoso = vosoct - vodrt
    dwoso = wosoct - wodrt

999 continue
end subroutine cdappl
