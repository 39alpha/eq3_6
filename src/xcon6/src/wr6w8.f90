subroutine wr6w8(awmaxi,awmini,cbsri,cdac,cesri,cgexj,csigma,dlaplo,dlaprn,dleplo,dleprn,dlhplo,dlhprn,dloplo,dloprn,dltpll,dltplo,dltprl,dltprn,dlxdmp,dlxmx0,dlxpll,dlxplo,dlxprl,dlxprn,eact,egersi,ehmaxi,ehmini,electr,fkrc,iact,ibsrti,iesrti,ietmax,igerti,iktmax,imchmx,imech,iodb,iopg,iopr,iopt,itermx,ixrti,jcode,jetmax,jflgi,jgerti,jgext,jpress,jreac,jtemp,kbt,kct,kdim,kmax,kmt,kprs,ksplmx,ksppmx,kstpmx,kxmod,kxt,hact,modr,moffg,morr,mprphi,mprspi,mtbaqi,mtbi,mwtges,nbti,nbtmax,nbt1mx,nctmax,ndact,ndctmx,nert,nertmx,net,netmax,newin,nffg,nffgmx,ngexrt,nobswt,nodbmx,nopgmx,noprmx,noptmx,nordmx,nprpmx,nprpti,nprsmx,nprsti,nptkmx,nrct,nrctmx,nrk,nsbswt,nsk,nsrt,nsrtmx,ntitl1,ntitl2,ntitmx,ntrymx,nttkmx,nxmdmx,nxmod,nxopex,nxopmx,nxopt,nxpemx,nxrt,nxrtmx,o2maxi,o2mini,phmaxi,phmini,pressb,pressi,ptk,qgexsh,rkb,rxbari,sfcar,ssfcar,tempcb,tempci,tgexp,timmxi,tistti,tolbt,toldl,tolsat,tolxsf,trkb,ttk,ubmtbi,ubsri,ucxri,udac,uesri,uffg,ugerji,ugermo,ugersi,ugexj,ugexmo,ugexp,ugexr,uhfgex,uobsw,uprphi,uprspi,ureac,usbsw,utitl1,utitl2,uvfgex,uxcat,uxkgex,uxmod,uxopex,uxopt,uzveci,vreac,xgersi,xhfgex,ximaxi,xistti,xlkffg,xlkgex,xlkmod,xvfgex,zgexj,zvclgi)
    !! This subroutine writes the EQ6 INPUT file in compact ("W") format
    !! for version 8.0.
    !! This subroutine is a near-clone of EQLIB/wr6pkw.f.
    !! The calling sequence of this subroutine is identical to that of
    !! XCON6/wr6d8.f, EQLIB/wr6pkw.f, and EQLIB/wr6pkd.f.
    !! The calling sequence of this subroutine is identical to that of
    !! EQ6/rd6inw.f, EQ6/rd6ind.f, XCON6/rd6w8.f, and XCON6/rd6d8.f,
    !! except that newin is added and ninpts, nprob, noutpt, qend,
    !! and qrderr are deleted.
    !! This subroutine is called by:
    !!   XCON6/xcon6.f
    !! Principal input:
    !! Principal output:
    !!   None
    implicit none

    ! Calling sequence variable declarations.
    integer :: ietmax
    integer :: iktmax
    integer :: imchmx
    integer :: jetmax
    integer :: kmax
    integer :: nbtmax
    integer :: nbt1mx
    integer :: nctmax
    integer :: ndctmx
    integer :: nertmx
    integer :: netmax
    integer :: nffgmx
    integer :: nodbmx
    integer :: nopgmx
    integer :: noprmx
    integer :: noptmx
    integer :: nordmx
    integer :: nprpmx
    integer :: nprsmx
    integer :: nptkmx
    integer :: nrctmx
    integer :: nsrtmx
    integer :: ntitmx
    integer :: nttkmx
    integer :: nxmdmx
    integer :: nxopmx
    integer :: nxpemx
    integer :: nxrtmx

    integer :: newin

    integer :: iact(imchmx,2,nrctmx)
    integer :: ibsrti(nsrtmx)
    integer :: iesrti(nsrtmx)
    integer :: igerti(jetmax,nertmx)
    integer :: imech(2,nrctmx)
    integer :: iodb(nodbmx)
    integer :: iopg(nopgmx)
    integer :: iopr(noprmx)
    integer :: iopt(noptmx)
    integer :: ixrti(nxrtmx)
    integer :: jcode(nrctmx)
    integer :: jflgi(nbtmax)
    integer :: jgerti(nertmx)
    integer :: jgext(netmax)
    integer :: jreac(nrctmx)
    integer :: kxmod(nxmdmx)
    integer :: ndact(imchmx,2,nrctmx)
    integer :: ngexrt(jetmax,netmax)
    integer :: nrk(2,nrctmx)
    integer :: nsk(nrctmx)

    integer :: itermx
    integer :: jpress
    integer :: jtemp
    integer :: kbt
    integer :: kct
    integer :: kdim
    integer :: kmt
    integer :: kprs
    integer :: ksplmx
    integer :: ksppmx
    integer :: kstpmx
    integer :: kxt
    integer :: nbti
    integer :: nffg
    integer :: nert
    integer :: net
    integer :: nobswt
    integer :: nprpti
    integer :: nprsti
    integer :: nrct
    integer :: nsbswt
    integer :: nsrt
    integer :: ntitl1
    integer :: ntitl2
    integer :: ntrymx
    integer :: nxmod
    integer :: nxopex
    integer :: nxopt
    integer :: nxrt

    logical :: qgexsh

    character(len=80) :: utitl1(ntitmx)
    character(len=80) :: utitl2(ntitmx)
    character(len=56) :: ugexr(ietmax,jetmax,netmax)
    character(len=48) :: ubmtbi(nbtmax)
    character(len=48) :: uobsw(2,nbtmax)
    character(len=48) :: uprspi(nprsmx)
    character(len=48) :: usbsw(2,nbtmax)
    character(len=48) :: uxmod(nxmdmx)
    character(len=48) :: uzveci(kmax)
    character(len=24) :: ubsri(nbt1mx,nsrtmx)
    character(len=24) :: ucxri(iktmax,nxrtmx)
    character(len=24) :: udac(ndctmx,imchmx,2,nrctmx)
    character(len=24) :: uffg(nffgmx)
    character(len=24) :: ugermo(nertmx)
    character(len=24) :: ugersi(ietmax,jetmax,nertmx)
    character(len=24) :: ugexmo(netmax)
    character(len=24) :: ugexp(netmax)
    character(len=24) :: uprphi(nprpmx)
    character(len=24) :: ureac(nrctmx)
    character(len=24) :: uxcat(nxopmx)
    character(len=24) :: uxopex(nxpemx)
    character(len=8) :: uesri(nctmax,nsrtmx)
    character(len=8) :: ugerji(jetmax,nertmx)
    character(len=8) :: ugexj(jetmax,netmax)
    character(len=8) :: uhfgex(ietmax,jetmax,netmax)
    character(len=8) :: uvfgex(ietmax,jetmax,netmax)
    character(len=8) :: uxkgex(ietmax,jetmax,netmax)
    character(len=8) :: uxopt(nxopmx)

    real(kind=8) :: cbsri(nbt1mx,nsrtmx)
    real(kind=8) :: cdac(ndctmx,imchmx,2,nrctmx)
    real(kind=8) :: cesri(nctmax,nsrtmx)
    real(kind=8) :: cgexj(jetmax,netmax)
    real(kind=8) :: csigma(imchmx,2,nrctmx)
    real(kind=8) :: eact(imchmx,2,nrctmx)
    real(kind=8) :: egersi(ietmax,jetmax,nertmx)
    real(kind=8) :: fkrc(nrctmx)
    real(kind=8) :: hact(imchmx,2,nrctmx)
    real(kind=8) :: modr(nrctmx)
    real(kind=8) :: moffg(nffgmx)
    real(kind=8) :: morr(nrctmx)
    real(kind=8) :: mprphi(nprpmx)
    real(kind=8) :: mprspi(nprsmx)
    real(kind=8) :: mtbaqi(nbtmax)
    real(kind=8) :: mtbi(nbtmax)
    real(kind=8) :: mwtges(netmax)
    real(kind=8) :: ptk(nptkmx)
    real(kind=8) :: rkb(imchmx,2,nrctmx)
    real(kind=8) :: rxbari(iktmax,nxrtmx)
    real(kind=8) :: sfcar(nrctmx)
    real(kind=8) :: ssfcar(nrctmx)
    real(kind=8) :: tgexp(netmax)
    real(kind=8) :: trkb(imchmx,2,nrctmx)
    real(kind=8) :: ttk(nttkmx)
    real(kind=8) :: vreac(nrctmx)
    real(kind=8) :: xgersi(ietmax,jetmax,nertmx)
    real(kind=8) :: xhfgex(ietmax,jetmax,netmax)
    real(kind=8) :: xlkffg(nffgmx)
    real(kind=8) :: xlkgex(ietmax,jetmax,netmax)
    real(kind=8) :: xlkmod(nxmdmx)
    real(kind=8) :: xvfgex(ietmax,jetmax,netmax)
    real(kind=8) :: zvclgi(kmax)
    real(kind=8) :: zgexj(jetmax,netmax)

    real(kind=8) :: awmaxi
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
    real(kind=8) :: dlxmx0
    real(kind=8) :: dlxpll
    real(kind=8) :: dlxplo
    real(kind=8) :: dlxprl
    real(kind=8) :: dlxprn
    real(kind=8) :: ehmaxi
    real(kind=8) :: ehmini
    real(kind=8) :: electr
    real(kind=8) :: o2maxi
    real(kind=8) :: o2mini
    real(kind=8) :: phmaxi
    real(kind=8) :: phmini
    real(kind=8) :: pressb
    real(kind=8) :: pressi
    real(kind=8) :: tempcb
    real(kind=8) :: tempci
    real(kind=8) :: timmxi
    real(kind=8) :: tistti
    real(kind=8) :: tolbt
    real(kind=8) :: toldl
    real(kind=8) :: tolsat
    real(kind=8) :: tolxsf
    real(kind=8) :: ximaxi
    real(kind=8) :: xistti

    ! Local variable declarations.
    integer :: i
    integer :: iei
    integer :: iki
    integer :: je
    integer :: jei
    integer :: j2
    integer :: j3
    integer :: j4
    integer :: kcol
    integer :: krow
    integer :: n
    integer :: nbi
    integer :: nci
    integer :: ne
    integer :: ner
    integer :: npi
    integer :: nrc
    integer :: nsi
    integer :: nsr
    integer :: nxr

    integer :: ilnobl

    logical :: qx

    character(len=8) :: uendit

    data uendit /'endit.  '/

    ! Main title.
    do n = 1,ntitl1
        j2 = ilnobl(utitl1(n))
        write (newin,1020) utitl1(n)(1:j2)
1020 format(a)
    end do

    j3 = ilnobl(uendit)

    if (ntitl1 .lt. ntitmx) then
        write (newin,1020) uendit(1:j3)
    end if

    ! Temperature parameters.
    ! Note: ttk(1) = ttk1, etc.
    write (newin,1130) jtemp,tempcb,ttk(1),ttk(2)
1130 format(5x,'jtemp= ',i2,/4x,'tempcb= ',1pe12.5,/6x,'ttk1= ',1pe12.5,6x,'ttk2= ',e12.5)

    ! Pressure parameters.
    ! Note: ptk(1) = ptk1, etc.
    write (newin,1140) jpress,pressb,ptk(1),ptk(2)
1140 format(4x,'jpress= ',i2,/4x,'pressb= ',1pe12.5,/6x,'ptk1= ',1pe12.5,6x,'ptk2= ',e12.5)

    ! Number of reactants.
    write (newin,1490) nrct
1490 format(6x,'nrct= ',i2)

    write (newin,1500)
1500 format('*-------------------------------------------------','----------------------------')

    ! Reactants.
    nxr = 0
    nsr = 0
    ner = 0

    do nrc = 1,nrct
        ! Name, flags, and masses.
        j2 = ilnobl(ureac(nrc))
        write (newin,1540) ureac(nrc)(1:j2),jcode(nrc),jreac(nrc),morr(nrc),modr(nrc)
1540 format(2x,'reactant= ',a,/5x,'jcode= ',i2,15x,'jreac= ',i2,/6x,'morr= ',1pe12.5,6x,'modr= ',e12.5)

        if (jcode(nrc) .eq. 1) then
            ! Solid solution compositions.
            nxr = nxr + 1

            do iki = 1,ixrti(nxr)
                write (newin,1580) ucxri(iki,nxr),rxbari(iki,nxr)
1580 format(3x,a24,3x,1pe12.5)
            end do

            j3 = ilnobl(uendit)
            write (newin,1570) uendit(1:j3)
1570 format(3x,a)
        else if (jcode(nrc) .eq. 2) then
            ! Special reactant compositions.
            nsr = nsr + 1
            write (newin,1620) vreac(nrc)
1620 format(5x,'vreac= ',e12.5)

            write (newin,1640)
1640 format ('* Elemental composition')

            do nci = 1,iesrti(nsr)
                write (newin,1650) uesri(nci,nsr),cesri(nci,nsr)
1650 format(3x,a8,3x,1pe22.15)
            end do

            j3 = ilnobl(uendit)
            write (newin,1570) uendit(1:j3)

            write (newin,1680)
1680 format ('* Reaction')

            do nbi = 1,ibsrti(nsr)
                write (newin,1690) ubsri(nbi,nsr)(1:24),cbsri(nbi,nsr)
1690 format(3x,a24,3x,1pe22.15)
            end do

            write (newin,1570) uendit(1:j3)
        else if (jcode(nrc) .eq. 5) then
            ! Generic ion exchanger phase compositions.
            ner = ner + 1

            j3 = ilnobl(ugermo(ner))
            write (newin,1692) ugermo(ner)(1:j3)
1692 format(4x,'ugermo= ',a)

            j4 = ilnobl(uendit)

            do jei = 1,jgerti(ner)
                j3 = ilnobl(ugerji(jei,ner))
                write (newin,1700) ugerji(jei,ner)(1:j3)
1700 format(3x,a)

                do iei = 1,igerti(jei,ner)
                    write (newin,1710) ugersi(iei,jei,ner),egersi(iei,jei,ner)
1710 format(6x,a24,3x,1pe12.5)
                end do

                write (newin,1720) uendit(1:j4)
1720 format(6x,a)
            end do

            write (newin,1570) uendit(1:j4)
        end if

        ! Surface area parameters.
        write (newin,1750) nsk(nrc),sfcar(nrc),ssfcar(nrc),fkrc(nrc)
1750 format(7x,'nsk= ',i2,15x,'sfcar= ',1pe12.5,4x,'ssfcar= ',e12.5,/6x,'fkrc= ',e12.5)

        ! Kinetic rate laws.
        write (newin,1770) nrk(1,nrc),nrk(2,nrc)
1770 format(6x,'nrk1= ',i2,16x,'nrk2= ',i2)

        ! Rate law parameters, forward direction (destruction).
        if (nrk(1,nrc) .eq. 1) then
            ! Arbitrary kinetics.
            write (newin,1800) (rkb(i,1,nrc), i = 1,3)
1800 format(6x,'rkb1= ',1pe12.5,6x,'rkb2= ',e12.5,6x,'rkb3= ',e12.5)
        else if (nrk(1,nrc) .eq. 2) then
            ! Transition state theory.
            write (newin,1820) imech(1,nrc)
1820 format(5x,'imech= ',i2)

            do i = 1,imech(1,nrc)
                write (newin,1840) rkb(i,1,nrc),trkb(i,1,nrc),iact(i,1,nrc)
1840 format(7x,'rkb= ',1pe12.5,6x,'trkb= ',1pe12.5,6x,'iact= ',i2)

                write (newin,1860) eact(i,1,nrc),hact(i,1,nrc)
1860 format(6x,'eact= ',1pe12.5,6x,'hact= ',e12.5)

                write (newin,1880) ndact(i,1,nrc),csigma(i,1,nrc)
1880 format(5x,'ndact= ',i2,14x,'csigma= ',1pe12.5)

                do n = 1,ndact(i,1,nrc)
                    write (newin,1920) udac(n,i,1,nrc),cdac(n,i,1,nrc)
1920 format(6x,'udac= ',a24,6x,'cdac= ',1pe12.5)
                end do
            end do
        else if (nrk(1,nrc) .eq. 3) then
            ! Linear rate law.
            i = imech(1,nrc)
            write (newin,1840) rkb(i,1,nrc),trkb(i,1,nrc),iact(i,1,nrc)
            write (newin,1860) eact(i,1,nrc),hact(i,1,nrc)
        end if

        ! Rate law parameters, backward direction (formation).
        if (nrk(2,nrc) .eq. 1) then
            ! Arbitrary kinetics.
            write (newin,1800) (rkb(i,2,nrc), i = 1,3)
        else if (nrk(2,nrc) .eq. 2) then
            ! Transition state theory.
            write (newin,1820) imech(2,nrc)

            do i = 1,imech(2,nrc)
                write (newin,1840) rkb(i,2,nrc),trkb(i,2,nrc),iact(i,2,nrc)

                write (newin,1860) eact(i,2,nrc),hact(i,2,nrc)

                write (newin,1880) ndact(i,2,nrc),csigma(i,2,nrc)

                do n = 1,ndact(i,2,nrc)
                    write (newin,1920) udac(n,i,2,nrc),cdac(n,i,2,nrc)
                end do
            end do
        else if (nrk(2,nrc) .eq. 3) then
            ! Linear rate law.
            i = imech(2,nrc)
            write (newin,1840) rkb(i,2,nrc),trkb(i,2,nrc),iact(i,1,nrc)
            write (newin,1860) eact(i,2,nrc),hact(i,2,nrc)
        end if

        write (newin,1500)
    end do

    ! Starting, minimum, and maximum values of key run parameters.
    write (newin,1150) xistti,ximaxi,tistti,timmxi
1150 format(4x,'xistti= ',1pe12.5,4x,'ximaxi= ',e12.5,/4x,'tistti= ',1pe12.5,4x,'timmxi= ',e12.5)

    write (newin,1160) phmini,phmaxi,ehmini,ehmaxi,o2mini,o2maxi,awmini,awmaxi,kstpmx
1160 format(4x,'phmini= ',1pe12.5,4x,'phmaxi= ',e12.5,/4x,'ehmini= ',1pe12.5,4x,'ehmaxi= ',e12.5,/4x,'o2mini= ',1pe12.5,4x,'o2maxi= ',e12.5,/4x,'awmini= ',1pe12.5,4x,'awmaxi= ',e12.5,/4x,'kstpmx= ',i12)

    ! Print interval parameters.
    write (newin,1170) dlxprn,dlxprl,dltprn,dltprl
1170 format(4x,'dlxprn= ',1pe12.5,4x,'dlxprl= ',e12.5,/4x,'dltprn= ',e12.5,4x,'dltprl= ',e12.5)

    write (newin,1180) dlhprn,dleprn,dloprn,dlaprn,ksppmx
1180 format(4x,'dlhprn= ',1pe12.5,4x,'dleprn= ',e12.5,/4x,'dloprn= ',e12.5,4x,'dlaprn= ',e12.5,/4x,'ksppmx= ',i12)

    ! Plot interval parameters.
    write (newin,1190) dlxplo,dlxpll,dltplo,dltpll
1190 format(4x,'dlxplo= ',1pe12.5,4x,'dlxpll= ',e12.5,/4x,'dltplo= ',e12.5,4x,'dltpll= ',e12.5)

    write (newin,1195) dlhplo,dleplo,dloplo,dlaplo,ksplmx
1195 format(4x,'dlhplo= ',1pe12.5,4x,'dleplo= ',e12.5,/4x,'dloplo= ',e12.5,4x,'dlaplo= ',e12.5,/4x,'ksplmx= ',i12)

    ! Comment header for iopt, iopr, iodb, and iopg option switches.
    write (newin,1200)
1200 format('*',15x,'1    2    3    4    5    6    7    8    9   10')

    ! Iopt option switches.
    ! Note: iopt(1) = iopt1, etc.
    write (newin,1220) (iopt(i), i = 1,20)
1220 format(2x,'iopt1-10= ',10i5,/1x,'iopt11-20= ',10i5)

    ! Iopr option switches.
    ! Note: iopr(1) = iopr1, etc.
    write (newin,1240) (iopr(i), i = 1,20)
1240 format(2x,'iopr1-10= ',10i5,/1x,'iopr11-20= ',10i5)

    ! Iodb option switches.
    ! Note: iodb(1) = iodb1, etc.
    write (newin,1260) (iodb(i), i = 1,20)
1260 format(2x,'iodb1-10= ',10i5,/1x,'iodb11-20= ',10i5)

    ! Number of nxopt options.
    write (newin,1310) nxopt
1310 format(5x,'nxopt= ',i2)

    ! Nxopt options.
    if (nxopt .gt. 0) then
        do n = 1,nxopt
            j2 = ilnobl(uxcat(n))
            write (newin,1340) uxopt(n),uxcat(n)(1:j2)
1340 format(4x,'option= ',a6,1x,a)
        end do

        write (newin,1360) nxopex
1360 format(4x,'nxopex= ',i2)

        if (nxopex .gt. 0) then
            do n = 1,nxopex
                j2 = ilnobl(uxopex(n))
                write (newin,1390) uxopex(n)(1:j2)
1390 format(' exception= ',a)
            end do
        end if
    end if

    ! Number of nffg options.
    write (newin,1400) nffg
1400 format(6x,'nffg= ',i2)

    ! Nffg options.
    if (nffg .gt. 0) then
        do n = 1,nffg
            j2 = ilnobl(uffg(n))
            write (newin,1440) uffg(n)(1:j2),moffg(n),xlkffg(n)
1440 format (3x,'species= ',a,/5x,'moffg= ',1pe12.5,4x,'xlkffg= ',e12.5)
        end do
    end if

    ! Maximum finite-difference order.
    write (newin,2020) nordmx
2020 format(4x,'nordmx= ',i3)

    ! Newton-Raphson convergence tolerances.
    write (newin,2030) tolbt,toldl
2030 format(5x,'tolbt= ',1pe12.5,5x,'toldl= ',e12.5)

    ! Maximum number of Newton-Raphson iterations.
    write (newin,2110) itermx
2110 format(4x,'itermx= ',i3)

    ! Search/find tolerance.
    write (newin,2032) tolxsf
2032 format(4x,'tolxsf= ',1pe12.5)

    ! Saturation tolerance.
    write (newin,2034) tolsat
2034 format(4x,'tolsat= ',1pe12.5)

    ! Maximum number of phase assemblage tries.
    write (newin,2112) ntrymx
2112 format(4x,'ntrymx= ',i3)

    ! Zero-order step size (in Xi).
    write (newin,2090) dlxmx0
2090 format(4x,'dlxmx0= ',1pe12.5)

    ! PRS transfer interval in Xi.
    write (newin,2010) dlxdmp
2010 format(4x,'dlxdmp= ',1pe12.5)

    write (newin,1500)

    ! Write the bottom half of the PICKUP file.
    ! Old title.
    do n = 1,ntitl2
        j2 = ilnobl(utitl2(n))
        write (newin,1020) utitl2(n)(1:j2)
    end do

    j3 = ilnobl(uendit)

    if (ntitl1 .lt. ntitmx) then
        write (newin,1020) uendit(1:j3)
    end if

    ! Special basis switches.
    write (newin,2610)
2610 format('*   Special basis switches')

    write (newin,2630) nsbswt
2630 format(4x,'nsbswt= ',i3)

    do n = 1,nsbswt
        j2 = ilnobl(usbsw(1,n))
        write (newin,2650) usbsw(1,n)(1:j2)
2650 format('species= ',a)

        j2 = ilnobl(usbsw(2,n))
        write (newin,2670) usbsw(2,n)(1:j2)
2670 format(2x,'switch with= ',a)
    end do

    ! Original temperature.
    write (newin,2690) tempci
2690 format(4x,'tempci= ',1pe12.5)

    ! Original pressure.
    write (newin,2710) pressi
2710 format(4x,'pressi= ',1pe12.5)

    ! Ion exchanger creation.
    write (newin,2790)
2790 format('* Ion exchanger creation')

    write (newin,2795) qgexsh
2795 format(4x,'qgexsh= ',l8)

    write (newin,2810) net
2810 format(7x,'net= ',i3)

    do ne = 1,net
        j2 = ilnobl(ugexp(ne))
        write (newin,2840) ugexp(ne)(1:j2)
2840 format(5x,'ugexp= ',a)

        write (newin,2860) mwtges(ne)
2860 format(4x,'mwtges= ',1pe12.5)

        j3 = ilnobl(ugexmo(ne))
        write (newin,2870) ugexmo(ne)(1:j3)
2870 format(4x,'ugexmo= ',a)

        write (newin,2880) tgexp(ne)
2880 format(5x,'tgexp= ',1pe12.5)

        write (newin,2900) jgext(ne)
2900 format(5x,'jgext= ',i3)

        do je = 1,jgext(ne)
            j3 = ilnobl(ugexj(je,ne))
            write (newin,2940) ugexj(je,ne)(1:j3)
2940 format(5x,'ugexj= ',a)

            write (newin,2960) cgexj(je,ne),zgexj(je,ne)
2960 format(5x,'cgexj= ',1pe12.5,5x,'zgexj= ',e12.5)

            write (newin,3010) ngexrt(je,ne)
3010 format(4x,'ngexrt= ',i3)

            do n = 1,ngexrt(je,ne)
                j2 = ilnobl(ugexr(n,je,ne))
                write (newin,3040) ugexr(n,je,ne)(1:j2)
3040 format(5x,'ugexr= ',a)

                j4 = ilnobl(uxkgex(n,je,ne))
                write (newin,3060) xlkgex(n,je,ne),uxkgex(n,je,ne)(1:j4)
3060 format(4x,'xlkgex= ',1pe12.5,5x,'units= ',a)

                j4 = ilnobl(uhfgex(n,je,ne))
                write (newin,3070) xhfgex(n,je,ne),uhfgex(n,je,ne)(1:j4)
3070 format(4x,'xhfgex= ',1pe12.5,5x,'units= ',a)

                j4 = ilnobl(uvfgex(n,je,ne))
                write (newin,3080) xvfgex(n,je,ne),uvfgex(n,je,ne)(1:j4)
3080 format(4x,'xvfgex= ',1pe12.5,5x,'units= ',a)
            end do
        end do
    end do

    ! Number of nxmod options.
    write (newin,3130) nxmod
3130 format(5x,'nxmod= ',i2)

    ! Nxmod options.
    do n = 1,nxmod
        j2 = ilnobl(uxmod(n))
        write (newin,3160) uxmod(n)(1:j2),kxmod(n),xlkmod(n)
3160 format(3x,'species= ',a,/4x,'option= ',i2,14x,'xlkmod= ',1pe12.5)
    end do

    ! Iopg options.
    ! Note: iopg(1) = iopg1, etc.
    write (newin,1200)
    write (newin,3180) (iopg(i), i = 1,20)
3180 format(2x,'iopg1-10= ',10i5,/1x,'iopg11-20= ',10i5)

    ! Index limits.
    write (newin,3230) kct,kbt,kmt,kxt,kdim,kprs
3230 format(7x,'kct= ',i2,17x,'kbt= ',i2,17x,'kmt= ',i2,/ 7x,'kxt= ',i2,16x,'kdim= ',i2,16x,'kprs= ',i2)

    ! Species for which mass balances are defined.
    write (newin,3330)
3330 format('* Data file basis species and jflgi values')

    do nbi = 1,nbti
        write (newin,3350) ubmtbi(nbi),jflgi(nbi)
3350 format(3x,a48,3x,i2)
    end do

    ! Mass balance totals.
    write (newin,3400)
3400 format('*',7x,'Mass balance totals     Aqueous mass balance totals')

    do nbi = 1,nbti
        write (newin,3420) mtbi(nbi),mtbaqi(nbi)
3420 format(6x,1pe22.15,6x,e22.15)
    end do

    write (newin,3440) electr
3440 format(9x,'Electrical imbalance= ',3x,1pe22.15)

    ! Ordinary basis switches.
    write (newin,3450)
3450 format('*   Ordinary basis switches')

    write (newin,3460) nobswt
3460 format(4x,'nobswt= ',i3)

    do n = 1,nobswt
        j2 = ilnobl(uobsw(1,n))
        write (newin,2650) uobsw(1,n)(1:j2)
        j2 = ilnobl(uobsw(2,n))
        write (newin,2670) uobsw(2,n)(1:j2)
    end do

    ! Matrix column variables.
    write (newin,3470)
3470 format('* Matrix species or entities')

    do krow = 1,kdim
        j2 = ilnobl(uzveci(krow))
        write (newin,3490) uzveci(krow)(1:j2)
3490 format(3x,a)
    end do

    ! Matrix variable values.
    write (newin,3510)
3510 format('*   Values of matrix variables')

    do kcol = 1,kdim
        write (newin,3530) zvclgi(kcol)
3530 format(3x,1pe22.15)
    end do

    qx = kprs.gt.0 .and. nprpti.gt.0 .and. nprsti.gt.0

    if (qx) then
        write (newin,3600)
3600 format('*  Phases and species in the PRS')

        nsi = 0

        do npi = 1,nprpti
            write (newin,3630) uprphi(npi),mprphi(npi)
3630 format(1x,a24,6x,1pe22.15)

            do iki = 1,iktmax
                if (uprspi(nsi + 1)(25:48) .eq. uprphi(npi)(1:24)) then
                    nsi = nsi + 1
                    write (newin,3670) uprspi(nsi)(1:24),mprspi(nsi)
3670 format(3x,a24,6x,1pe22.15)
                else
                    j3 = ilnobl(uendit)
                    write (newin,3660) uendit(1:j3)
3660 format(3x,a)

                    go to 300
                end if
            end do

300 continue
        end do

        j3 = ilnobl(uendit)
        write (newin,3620) uendit(1:j3)
3620 format(1x,a)
    end if
end subroutine wr6w8
