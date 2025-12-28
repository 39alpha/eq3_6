subroutine scripz(abar,acflg,acfw,acfwlg,actlg,actw,actwlg,affpd,affsd,afrc1,aft1,ah,ahmes,ahnbs,ahrc,alki,alk1,alk2,awmax,awmin,a3bar,cbsr,cdrsd,cegexs,cesr,conc,conclg,csts,ctb,cteaq,dvoso,dwoso,egers,egexjc,egexjf,egexpa,egexpc,egexs,egexw,eh,ehmax,ehmes,ehmin,ehnbs,ehrc,elecsr,electr,fje,fjest,fo2,fo2lg,fo2lrc,fugac,fugalg,fxi,fxist,iaqsln,iemop,iemop0,iemos,iemos0,iern1,iern2,iexr,iexrt,ifrn1,ifrn2,ilrn1,ilrn2,imech,imrn1,imrn2,iopg,iopr,iopt,ipndx1,ixrn1,ixrn2,jcode,jcsort,jern1,jern2,jexr,jexrt,jflag,jflagd,jflgi,jgext,jgsort,jpflag,jreac,jsca,jscat,jscr,jscrt,jsflag,jsol,jssort,kbt,kern1,kern2,kgexsa,km1,kmt,kx1,kxt,kstep,kstpmx,loph,losp,mlmrra,modr,moph,mophg,mophj,mopht,morr,mosp,mospg,mospj,mospt,mprph,mprsp,mrgers,mrmlra,mwtrc,mwtsp,narn1,narn2,nat,nbasp,nbaspd,nbt,ncmpe,ncmpe0,ncmpr,nct,ndrsd,ndrsrd,nelect,nern1,nern2,nert,net,nfrn1,nfrn2,ngext,ngexsa,ngrn1,ngrn2,ngt,nhydr,nhydx,nlrn1,nlrn2,nlt,nmrn1,nmrn2,nmrt,nmt,noutpt,no2gaq,npet,npet0,npt,npts,nrct,nrdxsp,nrk,nrndex,nst,nsts,nstsr,ntf1t,ntf2t,nxridx,nxrn1,nxrn2,nxrt,nxt,osc,oscst,omega,o2max,o2min,pch,pe,pemes,penbs,perc,ph,phcl,phmax,phmes,phmin,phnbs,ppmwe,presg,press,qaft1,qftpr2,qmod,qphcl,qredox,qrho,qriinf,qstopx,qvhfxi,qvlsow,qzprnt,rho,rhoc,rhowc,rk,rreacn,rreac1,rrelr1,rxbar,sfcar,sidrph,sidrsp,sigmst,sigmam,ssfcar,tdays,tdsglw,tdspkc,tdsplc,tempc,thours,time1,timemx,tmins,tolsat,tolxsf,tolxst,tolxsu,tyears,uelem,ugermo,ugexj,ugexmo,uphase,ureac,uspec,uxtype,vodrt,voph,vophg,vophj,vopht,vosoct,vosol,vosp,vospg,vospj,vospt,vreac,wfh2o,wftds,wkgwi,woh2o,wodr,wodrt,woph,wophg,wophj,wopht,worr,worrt,wosoct,wosol,wosp,wospg,wospj,wospt,wotds,xbar,xbarlg,xbarw,xbrwlg,xgers,xgexw,xi1,xidump,ximax,xistsv,xirct,zchar)
    !! This subroutine writes a detailed description on the output file
    !! of the modeled system at the current value of reaction progress
    !! (xi1).
    !! NOTE: This subroutine is designed to carry out a largely pure
    !! write function. Very few things to be printed by this subroutine
    !! should be computed within it. Unless the data are already
    !! available in EQ6/path.f, they should be calculated in
    !! EQ6/cdappl.f, which is called by EQ6/path.f prior to printing or
    !! plotting. This practice should be followed so that data needed for
    !! printing are also available for plotting. It is acceptable to make
    !! certain simple calculations in this subroutine, such as taking
    !! logarithms or converting units. If need be, these can be repeated
    !! elsewhere for plotting.
    !! Note that qftpr2 is a logical flag denoting a print being made
    !! only for a special fluid-centered flow-through open system just
    !! after a discontinuity in the derivatives along the reaction
    !! path.
    !! This subroutine is called by:
    !!   EQ6/path.f
    !! Principal input:
    !! Principal output:
    implicit none

    include 'eqlib/eqlpar.h'
    include 'eqlib/eqldv.h'

    ! Calling sequence variable declarations.
    integer :: noutpt

    integer :: jcode(nrctmx)
    integer :: jreac(nrctmx)
    integer :: nrndex(nrctmx)
    integer :: nxridx(nrctmx)

    integer :: iemop(npetmx)
    integer :: iemop0(npetmx)
    integer :: iemos(nsetmx)
    integer :: iemos0(nsetmx)
    integer :: iexr(nrctmx)
    integer :: imech(2,nrctmx)
    integer :: iopg(nopgmx)
    integer :: iopr(noprmx)
    integer :: iopt(noptmx)
    integer :: ipndx1(kmax)
    integer :: jcsort(nstmax)
    integer :: jern1(jetmax,netmax)
    integer :: jern2(jetmax,netmax)
    integer :: jexr(nrctmx)
    integer :: jflag(nstmax)
    integer :: jflagd(nstmax)
    integer :: jflgi(nbtmax)
    integer :: jgext(netmax)
    integer :: jgsort(ngtmax)
    integer :: jpflag(nptmax)
    integer :: jsca(nrctmx)
    integer :: jscr(nrctmx)
    integer :: jsflag(nstmax)
    integer :: jsol(nxtmax)
    integer :: jssort(nstmax)
    integer :: kern1(netmax)
    integer :: kern2(netmax)
    integer :: kgexsa(ketmax,netmax)
    integer :: nbasp(nbtmax)
    integer :: nbaspd(nbtmax)
    integer :: ncmpe(2,npetmx)
    integer :: ncmpe0(2,npetmx)
    integer :: ncmpr(2,nptmax)
    integer :: ndrsd(ndrsmx)
    integer :: ndrsrd(2,nstmax)
    integer :: ngexsa(ietmax,jetmax,netmax)
    integer :: ngext(jetmax,netmax)
    integer :: nrk(2,nrctmx)
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

    integer :: iaqsln
    integer :: iexrt
    integer :: jexrt
    integer :: jscat
    integer :: jscrt
    integer :: kbt
    integer :: km1
    integer :: kmt
    integer :: kx1
    integer :: kxt
    integer :: kstep
    integer :: kstpmx
    integer :: nelect
    integer :: nert
    integer :: nhydr
    integer :: nhydx
    integer :: nmrt
    integer :: no2gaq
    integer :: npet
    integer :: npet0
    integer :: npts
    integer :: nrct
    integer :: nrdxsp
    integer :: ntf1t
    integer :: ntf2t
    integer :: nxrt

    logical :: qaft1
    logical :: qftpr2
    logical :: qmod
    logical :: qphcl
    logical :: qredox
    logical :: qrho
    logical :: qriinf
    logical :: qstopx
    logical :: qvhfxi
    logical :: qvlsow
    logical :: qzprnt

    character(len=48) :: uspec(nstmax)
    character(len=32) :: uxtype(jsomax)
    character(len=24) :: ugermo(nertmx)
    character(len=24) :: ureac(nrctmx)
    character(len=24) :: ugexmo(netmax)
    character(len=24) :: uphase(nptmax)
    character(len=8) :: uelem(nctmax)
    character(len=8) :: ugexj(jetmax,netmax)

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
    real(kind=8) :: sfcar(nrctmx)
    real(kind=8) :: ssfcar(nrctmx)
    real(kind=8) :: vreac(nrctmx)
    real(kind=8) :: wodr(nrctmx)
    real(kind=8) :: worr(nrctmx)
    real(kind=8) :: xgers(ietmax,jetmax,nertmx)

    real(kind=8) :: acflg(nstmax)
    real(kind=8) :: actlg(nstmax)
    real(kind=8) :: affpd(nptmax)
    real(kind=8) :: affsd(nstmax)
    real(kind=8) :: afrc1(nrctmx)
    real(kind=8) :: ahrc(nbtmax)
    real(kind=8) :: cdrsd(ndrsmx)
    real(kind=8) :: cegexs(ietmax,jetmax,netmax)
    real(kind=8) :: conc(nstmax)
    real(kind=8) :: conclg(nstmax)
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
    real(kind=8) :: loph(nptmax)
    real(kind=8) :: losp(nstmax)
    real(kind=8) :: moph(nptmax)
    real(kind=8) :: mophg(nptmax)
    real(kind=8) :: mophj(nptmax)
    real(kind=8) :: mopht(nptmax)
    real(kind=8) :: mosp(nstmax)
    real(kind=8) :: mospg(nstmax)
    real(kind=8) :: mospj(nstmax)
    real(kind=8) :: mospt(nstmax)
    real(kind=8) :: mprph(nptmax)
    real(kind=8) :: mprsp(nstmax)
    real(kind=8) :: mwtsp(nstmax)

    real(kind=8) :: perc(nbtmax)
    real(kind=8) :: ppmwe(nctmax)
    real(kind=8) :: rk(imchmx,2,nrctmx)
    real(kind=8) :: sidrph(nptmax)
    real(kind=8) :: sidrsp(nstmax)
    real(kind=8) :: voph(nptmax)
    real(kind=8) :: vophg(nptmax)
    real(kind=8) :: vophj(nptmax)
    real(kind=8) :: vopht(nptmax)
    real(kind=8) :: vosp(nstmax)
    real(kind=8) :: vospg(nstmax)
    real(kind=8) :: vospj(nstmax)
    real(kind=8) :: vospt(nstmax)
    real(kind=8) :: woph(nptmax)
    real(kind=8) :: wophg(nptmax)
    real(kind=8) :: wophj(nptmax)
    real(kind=8) :: wopht(nptmax)
    real(kind=8) :: wosp(nstmax)
    real(kind=8) :: wospg(nstmax)
    real(kind=8) :: wospj(nstmax)
    real(kind=8) :: wospt(nstmax)
    real(kind=8) :: xbar(nstmax)
    real(kind=8) :: xbarlg(nstmax)
    real(kind=8) :: xgexw(ketmax,netmax)
    real(kind=8) :: xirct(nrctmx)
    real(kind=8) :: zchar(nstmax)

    real(kind=8) :: abar
    real(kind=8) :: acfw
    real(kind=8) :: acfwlg
    real(kind=8) :: actw
    real(kind=8) :: actwlg
    real(kind=8) :: aft1
    real(kind=8) :: ah
    real(kind=8) :: ahmes
    real(kind=8) :: ahnbs
    real(kind=8) :: alki
    real(kind=8) :: alk1
    real(kind=8) :: alk2
    real(kind=8) :: awmax
    real(kind=8) :: awmin
    real(kind=8) :: a3bar
    real(kind=8) :: dvoso
    real(kind=8) :: dwoso
    real(kind=8) :: eh
    real(kind=8) :: ehmax
    real(kind=8) :: ehmin
    real(kind=8) :: ehmes
    real(kind=8) :: ehnbs
    real(kind=8) :: electr
    real(kind=8) :: fje
    real(kind=8) :: fjest
    real(kind=8) :: fo2
    real(kind=8) :: fo2lg
    real(kind=8) :: fxi
    real(kind=8) :: fxist
    real(kind=8) :: osc
    real(kind=8) :: oscst
    real(kind=8) :: omega
    real(kind=8) :: o2max
    real(kind=8) :: o2min
    real(kind=8) :: pch
    real(kind=8) :: pe
    real(kind=8) :: penbs
    real(kind=8) :: pemes
    real(kind=8) :: ph
    real(kind=8) :: phcl
    real(kind=8) :: phmax
    real(kind=8) :: phmes
    real(kind=8) :: phmin
    real(kind=8) :: phnbs
    real(kind=8) :: presg
    real(kind=8) :: press
    real(kind=8) :: rho
    real(kind=8) :: sigmam
    real(kind=8) :: sigmst
    real(kind=8) :: tdays
    real(kind=8) :: tempc
    real(kind=8) :: thours
    real(kind=8) :: time1
    real(kind=8) :: timemx
    real(kind=8) :: tmins
    real(kind=8) :: tolsat
    real(kind=8) :: tolxsf
    real(kind=8) :: tolxst
    real(kind=8) :: tolxsu
    real(kind=8) :: tyears
    real(kind=8) :: vodrt
    real(kind=8) :: vosoct
    real(kind=8) :: vosol
    real(kind=8) :: wfh2o
    real(kind=8) :: wftds
    real(kind=8) :: wkgwi
    real(kind=8) :: wodrt
    real(kind=8) :: woh2o
    real(kind=8) :: worrt
    real(kind=8) :: wosoct
    real(kind=8) :: wosol
    real(kind=8) :: wotds
    real(kind=8) :: xbarw
    real(kind=8) :: xbrwlg
    real(kind=8) :: xi1
    real(kind=8) :: xidump
    real(kind=8) :: ximax
    real(kind=8) :: xistsv

    real(kind=8) :: mlmrra
    real(kind=8) :: mrmlra
    real(kind=8) :: rhoc
    real(kind=8) :: rhowc
    real(kind=8) :: tdsglw
    real(kind=8) :: tdspkc
    real(kind=8) :: tdsplc

    ! Local variable declarations.
    integer :: i
    integer :: ie
    integer :: j
    integer :: je
    integer :: j1
    integer :: j2
    integer :: j3
    integer :: kcol
    integer :: ncount
    integer :: ne
    integer :: np
    integer :: npe
    integer :: nplast
    integer :: np1
    integer :: ns
    integer :: nsr
    integer :: nrc
    integer :: nr1
    integer :: nr2
    integer :: nse
    integer :: nt

    integer :: ilnobl

    logical :: qrct
    logical :: qtimmx

    character(len=16) :: ux16a
    character(len=16) :: ux16b
    character(len=16) :: ux16c

    real(kind=8) :: betchb
    real(kind=8) :: chdrft
    real(kind=8) :: chdrfs
    real(kind=8) :: dp
    real(kind=8) :: elects
    real(kind=8) :: lxx
    real(kind=8) :: mxx
    real(kind=8) :: sigza
    real(kind=8) :: sigzc
    real(kind=8) :: sigzi
    real(kind=8) :: sigzia
    real(kind=8) :: sigzis
    real(kind=8) :: sigzm
    real(kind=8) :: sigzma
    real(kind=8) :: sigzms
    real(kind=8) :: tolspf
    real(kind=8) :: wconst
    real(kind=8) :: xilg

    real(kind=8) :: tlg

    wconst = omega/mosp(narn1)

    write (noutpt,1000)
1000 format(/' - - - - - - - - - - - - - - - - - - - - - - - - -',' - - - - - - -',/)

    xilg = -99999.

    if (xi1 .gt. 0.) then
        xilg = tlg(xi1)
    end if

    write (noutpt,1010) xi1,xilg
1010 format(/20x,'Xi= ',1pe12.5,/16x,'Log Xi= ',0pf12.5,/)

    if (iopt(2) .ge. 1) then
        if (qriinf) then
            write (noutpt,1020)
1020 format(/12x,'Time= infinity',/)
        else
            write (noutpt,1030) time1,tmins,thours,tdays,tyears
1030 format(/12x,'Time= ',1pe10.3,' seconds',/16x,'= ',1pe10.3,' minutes',/16x,'= ',1pe10.3,' hours',/16x,'= ',1pe10.3,' days',/16x,'= ',1pe10.3,' years',/)
        end if
    end if

    write (noutpt,1050) tempc
1050 format(/' Temperature= ',f6.2,' C')

    dp = press - presg

    if (abs(dp) .le. 1.e-4) then
        write (noutpt,1060) press
1060 format(/' Pressure= ',1pg12.5,' bars',/)
    else
        write (noutpt,1070) press,presg,dp
1070 format(/' Pressure= ',1pg12.5,' bars',/' Data file reference curve pressure= ',g12.5,' bars',/' Pressure difference= ',g12.5,' bars',/)
    end if

    if (qftpr2) then
        write(noutpt,1090)
1090 format(/' Fluid-centered flow-through open system print',/' following a discontinuity.')

        ! If this path is taken, only a partial description of
        ! the system is printed.
        go to 200
    end if

    if (qzprnt) then
        write (noutpt,1100)
    end if

1100 format(/' Required print point.')

    if (xi1.ge.xidump .and. xi1.gt.xistsv) then
        write (noutpt,1110)
    end if

1110 format(/' Required PRS shift.')

    if (qaft1) then
        write (noutpt,1120)
    end if

1120 format(/' The total affinity is repeatedly nearly zero.')

    if (iexrt .gt. 0) then
        write (noutpt,1130)
1130 format(/' Exhaustion of the following reactant(s):')

        do i = 1,iexrt
            nrc = iexr(i)
            j2 = ilnobl(ureac(nrc))
            write (noutpt,1140) ureac(nrc)(1:j2)
1140 format(5x,a)
        end do
    end if

    if (jexrt .gt. 0) then
        write (noutpt,1150)
1150 format(/' Step size is limited by the reactivation of',' the following formerly',/1x,'exhausted reactant(s):',/)

        do j = 1,jexrt
            nrc = jexr(j)
            j2 = ilnobl(ureac(nrc))
            write (noutpt,1140) ureac(nrc)
        end do
    end if

    if (jscat .gt. 0) then
        write (noutpt,1160)
1160 format(/' Step size is limited by the sign change of',' the affinity of the',/1x,'following reactant(s):',/)

        do j = 1,jscat
            nrc = jsca(j)
            j2 = ilnobl(ureac(nrc))
            write (noutpt,1140) ureac(nrc)
        end do
    end if

    if (jscrt .gt. 0) then
        write (noutpt,1170)
1170 format(/,' Step size is limited by the sign change of',' the relative rate of the',/1x,'following reactant(s):',/)

        do j = 1,jscrt
            nrc = jscr(j)
            j2 = ilnobl(ureac(nrc))
            write (noutpt,1140) ureac(nrc)
        end do
    end if

    if (xi1 .le. xistsv) then
        write (noutpt,1200)
    end if

1200 format(/' Start or restart of the run.')

    if (qmod) then
        write (noutpt,1210)
    end if

1210 format(/' Have a change in the ES phase assemblage.')

    if (xi1 .ge. ximax) then
        write (noutpt,1220)
    end if

1220 format(/' Maximum value of Xi.')

    qtimmx = iopt(2) .ge. 1 .and. time1.ge.( (1. - tolxst)*timemx )

    if (qtimmx) then
        write (noutpt,1230)
    end if

1230 format(/' Maximum value of time.')

    if (abs(ph - phmin) .le. tolxsu) then
        write (noutpt,1240)
1240 format(/' Minimum value of pH.')
    else if (ph .lt. phmin) then
        write (noutpt,1250)
1250 format(/' pH is less than the minimum value.')
    end if

    if (abs(ph - phmax) .le. tolxsu) then
        write (noutpt,1260)
1260 format(/' Maximum value of pH.')
    else if (ph .gt. phmax) then
        write (noutpt,1270)
1270 format(/' pH is greater than the maximum value.')
    end if

    if (qredox) then
        if (abs(eh - ehmin) .le. tolxsu) then
            write (noutpt,1280)
1280 format(/' Minimum value of Eh (v).')
        else if (eh .lt. ehmin) then
            write (noutpt,1290)
1290 format(/' Eh (v) is less than the minimum value.')
        end if

        if (abs(eh - ehmax) .le. tolxsu) then
            write (noutpt,1300)
1300 format(/' Maximum value of Eh (v).')
        else if (eh .gt. ehmax) then
            write (noutpt,1310)
1310 format(/' Eh (v) is greater than the maximum value.')
        end if
    end if

    if (qredox) then
        if (abs(fo2lg - o2min) .le. tolxsu) then
            write (noutpt,1320)
1320 format(/' Minimum value of log fO2.')
        else if (fo2lg .lt. o2min) then
            write (noutpt,1330)
1330 format(/' Log fO2 is less than the minimum value.')
        end if

        if (abs(fo2lg - o2max) .le. tolxsu) then
            write (noutpt,1340)
1340 format(/' Maximum value of log fO2.')
        else if (fo2lg .gt. o2max) then
            write (noutpt,1350)
1350 format(/' log fO2 is greater than the maximum value.')
        end if
    end if

    if (abs(actw - awmin) .le. tolxsu) then
        write (noutpt,1360)
1360 format(/' Minimum value of the activity of water.')
    else if (actw .lt. awmin) then
        write (noutpt,1370)
1370 format(/' The activity of water is less than the minimum',' value.')
    end if

    if (abs(actw - awmax) .le. tolxsu) then
        write (noutpt,1380)
1380 format(/' Maximum value of the activity of water.')
    else if (actw .gt. awmax) then
        write (noutpt,1390)
1390 format(/' The activity of water is greater than the maximum',' value.')
    end if

    if (kstep .ge. kstpmx) then
        write (noutpt,1450)
    end if

1450 format(/' Maximum number of steps.')

    if (qvlsow) then
        write (noutpt,1460)
    end if

1460 format(/' Have very nearly exhausted solvent water.')

    if (qvlsow) then
        write (noutpt,1470)
    end if

1470 format(/' Have extremely high ionic strength.')

    if (qstopx) then
        write (noutpt,1480)
    end if

1480 format(/' Early termination.')

    ! Print tables of data for reactants and reaction rates.
    call prtrct(afrc1,aft1,imchmx,imech,iopt,modr,morr,noptmx,noutpt,nrct,nrctmx,nrk,rk,rreacn,rreac1,rrelr1,sfcar,ureac,wodr,wodrt,worr,worrt,xi1,xistsv)

    ! Print a table of the elemental composition of the aqueous
    ! solution.
    call prteca(cteaq,mrmlra,nct,nctmax,noutpt,ppmwe,qrho,rho,uelem)

    ! Compute and print a table of the numerical compostion of the
    ! aqueous solution.
    call prtnca(ctb,jflag,jsflag,mrmlra,mwtsp,narn1,narn2,nbasp,nbaspd,nbt,nbtmax,noutpt,nstmax,qrho,rho,uspec,wfh2o)

    ! Compute and print a table of the sensible compostion of the
    ! aqueous solution.
    call prtsca(ctb,jflag,jsflag,mrmlra,mwtsp,narn1,narn2,nbasp,nbaspd,nbt,nbtmax,nelect,nhydr,nhydx,noutpt,no2gaq,nstmax,qrho,rho,uspec,wfh2o)

    ! Print pH, Eh, and pe-, all with reference to appropriate
    ! pH scales. Also print the pHCl.
    call prpheh(ah,ahmes,ahnbs,eh,ehmes,ehnbs,iopg,nopgmx,noutpt,pch,pe,pemes,penbs,ph,phcl,phmes,phnbs,qphcl,qredox,qrho)

    ! Print various aqueous phase parameters.
    call prtvpa(abar,acfw,acfwlg,actw,actwlg,a3bar,fje,fjest,fo2,fo2lg,fxi,fxist,iopg,mlmrra,mrmlra,nopgmx,noutpt,osc,oscst,qrho,rhoc,rhowc,sigmam,sigmst,tdsglw,tdspkc,tdsplc,vosol,wfh2o,wftds,woh2o,wosol,wotds,xbarw,xbrwlg)

    ! Print more precise mass results for the aqueous phase.
    write (ux16a,'(1pg15.8)') woh2o
    call lejust(ux16a)
    j1 = ilnobl(ux16a)
    write (ux16b,'(1pg15.8)') wotds
    call lejust(ux16b)
    j2 = ilnobl(ux16b)
    write (ux16c,'(1pg15.8)') wosol
    call lejust(ux16c)
    j3 = ilnobl(ux16c)
    write (noutpt,1500) ux16a(1:j1),ux16b(1:j2),ux16c(1:j3)
1500 format(//11x,'--- More Precise Aqueous Phase Masses ---',//23x,'Solvent mass= ',a,' g',/17x,'Solutes (TDS) mass= ',a,' g',/14x,'Aqueous solution mass= ',a,' g',//)

    ! Print a table of computed alkalinity parameters.
    call prtalk(alki,alk1,alk2,mrmlra,noutpt,ntf1t,ntf2t,qrho,rho,tempc,wfh2o)

    ! Compute and print a table describing the aqueous phase charge
    ! balance.
    call gszm(conc,jcsort,narn1,narn2,nstmax,sigza,sigzc,sigzi,sigzm,zchar)

    sigzia = sigzi/wconst
    sigzma = sigzm/wconst

    chdrft = sigzia - electr

    do nrc = 1,nrct
        if (jcode(nrc) .eq. 2) then
            nsr = nrndex(nrc)
            chdrft = chdrft - elecsr(nsr)*xirct(nrc)
        end if
    end do

    betchb = chdrft/sigzma
    sigzis = wfh2o*sigzi
    elects = wfh2o*electr
    chdrfs = wfh2o*chdrft
    sigzms = wfh2o*sigzma

    write (noutpt,1510) sigzia,electr,chdrft,sigzma,sigzis,elects,chdrfs,sigzms,betchb
1510 format(//7x,' --- Aqueous Solution Charge Balance ---',//6x,'    Actual Charge imbalance= ',1pe11.4,' eq',/6x,'  Expected Charge imbalance= ',e11.4,' eq',/6x,'         Charge discrepancy= ',e11.4,' eq',/6x,'        Sigma |equivalents|= ',e11.4,' eq',//6x,'    Actual Charge imbalance= ',e11.4,' eq/kg.solu',/6x,'  Expected Charge imbalance= ',e11.4,' eq/kg.solu',/6x,'         Charge discrepancy= ',e11.4,' eq/kg.solu',/6x,'        Sigma |equivalents|= ',e11.4,' eq/kg.solu',//6x,'Relative charge discrepancy= ',e11.4)

    ! Print detailed listing of aqueous species.
    call prtaqs(acflg,actlg,conc,conclg,iopr,jcsort,narn1,narn2,noprmx,noutpt,nstmax,uspec)

    ! Print the contributions of species to each aqueous mass total.
    call prtpct(conc,csts,ctb,iopr,jcsort,jflag,narn1,narn2,nbasp,nbaspd,nbt,nbtmax,nelect,nhydr,noprmx,no2gaq,noutpt,nstmax,nsts,nstsmx,nstsr,uspec)

    if (iopr(9) .gt. 0) then
        ! Compute and print the mean ionic activities and activity
        ! coefficients.
        call prtmip(acflg,actlg,conclg,ctb,nbaspd,nbt,nbtmax,nelect,nhydr,nhydx,noutpt,nstmax,uspec,zchar)
    end if

    if (iopr(5) .gt. 0) then
        ! Compute and print activity ratios of aqueous species.
        call prtacr(actlg,iopr,jsflag,nbaspd,nbt,nbtmax,nelect,nhydr,noprmx,no2gaq,noutpt,nstmax,uspec,zchar)
    end if

    if (qredox) then
        ! Print data for redox reactions that are not constrained to be
        ! at equilibrium.
        call prtrdx(ah,ahrc,cdrsd,eh,ehrc,fo2lg,fo2lrc,jflag,jsflag,narn1,nbasp,nbaspd,nbt,nbtmax,ndrsd,ndrsmx,ndrsrd,nelect,nhydr,no2gaq,noutpt,nstmax,pe,perc,uspec)
    end if

    if (iopt(1).eq.0 .or. iopt(1).eq.1) then
        ! List the solid phases in the equilibrium system (ES).
        write (noutpt,1600)
1600 format(//21x,'--- Summary of Solid Phases (ES) ---')

        write (noutpt,1610)
1610 format(/'   Phase/End-member',12x,'Log moles',5x,'Moles',8x,'Grams',5x,'Volume, cm3',/)

        ncount = 0

        do kcol = km1,kmt
            ncount = ncount + 1
            np = ipndx1(kcol)
            write (noutpt,1620) uphase(np),loph(np),moph(np),woph(np),voph(np)
1620 format(1x,a24,4x,f11.4,3(2x,1pe11.4))
        end do

        nplast = 0

        do kcol = kx1,kxt
            np = ipndx1(kcol)

            if (np .ne. nplast) then
                if (ncount .gt. 0) then
                    write (noutpt,1630)
                end if

1630 format(1x)

                ncount = ncount + 1
                write (noutpt,1620) uphase(np),loph(np),moph(np),woph(np),voph(np)

                nr1 = ncmpr(1,np)
                nr2 = ncmpr(2,np)
                nt = nr2 - nr1 + 1

                if (nt .ge. 2) then
                    do ns = nr1,nr2
                        if (jsflag(ns) .le. 0) then
                            write (noutpt,1640) uspec(ns),losp(ns),mosp(ns),wosp(ns),vosp(ns)
1640 format(3x,a24,2x,f11.4,3(2x,1pe11.4))
                        end if
                    end do
                end if

                nplast = np
            end if
        end do

        nplast = 0

        do np = iern1,iern2
            if (moph(np) .gt. 0.) then
                if (np .ne. nplast) then
                    if (ncount .gt. 0) then
                        write (noutpt,1630)
                    end if

                    ncount = ncount + 1
                    write (noutpt,1620) uphase(np),loph(np),moph(np),woph(np),voph(np)

                    ne = np - iern1 + 1

                    do je = 1,jgext(ne)
                        if (je .gt. 1) then
                            write (noutpt,1642)
                        end if

1642 format(1x)

                        ns = jern1(je,ne) - 1

                        do ie = 1,ngext(je,ne)
                            ns = ns + 1

                            if (jsflag(ns) .le. 0) then
                                write (noutpt,1640) uspec(ns),losp(ns),mosp(ns),wosp(ns),vosp(ns)
                            end if
                        end do
                    end do

                    nplast = np
                end if
            end if
        end do

        if (ncount .le. 0) then
            write (noutpt,1650)
        end if

1650 format(1x,'None')
    end if

200 continue

    if (iopt(1) .eq. 2) then
        ! List the mineral assemblage which is forming instantaneously
        ! with respect to Xi (on a per unit Xi basis). When this is
        ! calculated at a phase boundary, the assemblage given is the one
        ! to the left of the boundary. When this is calculated just
        ! after the boundary, it is the one to the right of the boundary.
        if (npts .gt. 0) then
            write (noutpt,1700)
1700 format(//3x,'--- Summary of Solid Phases (ES, Instantaneous',' Basis) ---',/19x,'(defined by derivatives)')

            write (noutpt,1710)
1710 format(//3x,'Phase/End-member',11x,'d mol/d Xi',2x,'d grams/d Xi',2x,'d V(cm3)/d Xi',/)

            ncount = 0

            do npe = 2,npet0
                np = iemop0(npe)
                ncount = ncount + 1

                if (np.ge.imrn1 .and. np.le.imrn2) then
                    ! Have a pure mineral.
                    write (noutpt,1720) uphase(np),mophj(np),wophj(np),vophj(np)
1720 format(1x,a24,2x,2(2x,1pe11.4),2x,1pe12.5,2x,1pe11.4)
                end if
            end do

            do npe = 2,npet0
                np = iemop0(npe)
                ncount = ncount + 1

                if (np.ge.ixrn1 .and. np.le.ixrn2) then
                    ! Have a solid solution.
                    write (noutpt,1730) uphase(np),mophj(np),wophj(np),vophj(np)
1730 format(/1x,a24,2x,2(2x,1pe11.4),2x,1pe12.5,2x,1pe11.4)

                    nr1 = ncmpe0(1,npe)
                    nr2 = ncmpe0(2,npe)

                    do nse = nr1,nr2
                        ns = iemos0(nse)
                        write (noutpt,1740) uspec(ns),mospj(ns),wospj(ns),vospj(ns)
1740 format(3x,a24,2(2x,1pe11.4),2x,1pe12.5,2x,1pe11.4)
                    end do
                end if
            end do

            do npe = 2,npet0
                np = iemop0(npe)
                ncount = ncount + 1

                if (np.ge.iern1 .and. np.le.iern2) then
                    ! Have a generic ion exchange phase.
                    write (noutpt,1730) uphase(np),mophj(np),wophj(np),vophj(np)
                    nr1 = ncmpe0(1,npe)
                    nr2 = ncmpe0(2,npe)

                    do nse = nr1,nr2
                        ns = iemos0(nse)
                        write (noutpt,1740) uspec(ns),mospj(ns),wospj(ns),vospj(ns)
                    end do
                end if
            end do

            if (ncount .le. 0) then
                write (noutpt,1650)
            end if
        end if
    end if

    if (iopt(1) .eq. 2) then
        ! Print a summary of the solid phases in the ES and the PRS.
        write (noutpt,1800)
1800 format(//12x,'--- Summary of Solid Phases (ES + PRS) ---')

        write (noutpt,1610)

        ncount = 0

        do np = 1,npt
            if (np .ne. iaqsln) then
                if (jpflag(np) .lt. 2) then
                    if (mopht(np) .gt. 0.) then
                        nr1 = ncmpr(1,np)
                        nr2 = ncmpr(2,np)
                        nt = nr2 - nr1 + 1
                        mxx = mopht(np)
                        lxx = tlg(mxx)

                        if (nt .lt. 2) then
                            ncount = ncount + 1
                            write (noutpt,1620) uphase(np),lxx,mxx,wopht(np),vopht(np)
                        else
                            if (ncount .gt. 0) then
                                write (noutpt,1630)
                            end if

                            ncount = ncount + 1
                            write (noutpt,1620) uphase(np),lxx,mxx,wopht(np),vopht(np)

                            do ns = nr1,nr2
                                if (mospt(ns) .gt. 0.) then
                                    mxx = mospt(ns)
                                    lxx = tlg(mxx)
                                    write (noutpt,1640) uspec(ns),lxx,mxx,wospt(ns),vospt(ns)
                                end if
                            end do
                        end if
                    end if
                end if
            end if
        end do

        if (ncount .le. 0) then
            write (noutpt,1650)
        end if
    end if

    if ((nmrt + nxrt + nert) .gt. 0) then
        write (noutpt,1810)
1810 format(//12x,'--- Grand Summary of Solid Phases (ES + PRS',' + Reactants) ---')

        write (noutpt,1610)

        ncount = 0

        do np = imrn1,imrn2
            qrct = .false.

            if (nmrt .gt. 0) then
                do nrc = 1,nrct
                    if (jcode(nrc) .eq. 0) then
                        np1 = nrndex(nrc)

                        if (np1 .eq. np) then
                            qrct = .true.
                            go to 300
                        end if
                    end if
                end do
            end if

300 continue
            if (mophg(np).gt.0. .or. qrct) then
                ncount = ncount + 1
                ns = ncmpr(1,np)
                mxx = mophg(np)
                lxx = tlg(abs(mxx))
                write (noutpt,1620) uphase(np),lxx,mxx,wophg(np),vophg(np)
            end if
        end do

        do np = ifrn1,ifrn2
            qrct = .false.

            if (nmrt .gt. 0) then
                do nrc = 1,nrct
                    if (jcode(nrc) .eq. 0) then
                        np1 = nrndex(nrc)

                        if (np1 .eq. np) then
                            qrct = .true.
                            go to 310
                        end if
                    end if
                end do
            end if

310 continue
            if (mophg(np).gt.0. .or. qrct) then
                ncount = ncount + 1
                ns = ncmpr(1,np)
                mxx = mophg(np)
                lxx = tlg(abs(mxx))
                write (noutpt,1620) uphase(np),lxx,mxx,wophg(np),vophg(np)
            end if
        end do

        if (iopt(4) .gt. 0) then
            do np = ixrn1,ixrn2
                if (jpflag(np) .lt. 2) then
                    qrct = .false.

                    if (nxrt .gt. 0) then
                        do nrc = 1,nrct
                            if (jcode(nrc) .eq. 1) then
                                np1 = nrndex(nrc)

                                if (np1 .eq. np) then
                                    qrct = .true.
                                    go to 320
                                end if
                            end if
                        end do
                    end if

320 continue
                    if (mophg(np).gt.0. .or. qrct) then
                        if (ncount .gt. 0) then
                            write (noutpt,1630)
                        end if

                        ncount = ncount + 1
                        mxx = mophg(np)
                        lxx = tlg(abs(mxx))
                        write (noutpt,1620) uphase(np),lxx,mxx,wophg(np),vophg(np)
                        nr1 = ncmpr(1,np)
                        nr2 = ncmpr(2,np)

                        do ns = nr1,nr2
                            mxx = mospg(ns)
                            lxx = tlg(abs(mxx))
                            write (noutpt,1640) uspec(ns),lxx,mxx,wospg(ns),vospg(ns)
                        end do
                    end if
                end if
            end do
        end if

        do np = iern1,iern2
            if (jpflag(np) .lt. 2) then
                qrct = .false.

                if (nert .gt. 0) then
                    do nrc = 1,nrct
                        if (jcode(nrc) .eq. 5) then
                            np1 = nrndex(nrc)

                            if (np1 .eq. np) then
                                qrct = .true.
                                go to 330
                            end if
                        end if
                    end do
                end if

330 continue
                if (mophg(np).gt.0. .or. qrct) then
                    if (ncount .gt. 0) then
                        write (noutpt,1630)
                    end if

                    ncount = ncount + 1
                    mxx = mophg(np)
                    lxx = tlg(abs(mxx))
                    write (noutpt,1620) uphase(np),lxx,mxx,wophg(np),vophg(np)
                    nr1 = ncmpr(1,np)
                    nr2 = ncmpr(2,np)

                    do ns = nr1,nr2
                        mxx = mospg(ns)

                        if (mxx .gt. 0.) then
                            lxx = tlg(abs(mxx))
                            write (noutpt,1640) uspec(ns),lxx,mxx,wospg(ns),vospg(ns)
                        end if
                    end do
                end if
            end if
        end do
    end if

    write (noutpt,1820) wosoct,vosoct,wodrt,vodrt,dwoso,dvoso
1820 format(//26x,'Mass, grams',7x,'Volume, cm3',//11x,'Created  ',6x,1pe12.5,6x,e12.5,/11x,'Destroyed',6x,e12.5,6x,e12.5,/11x,'Net      ',6x,e12.5,6x,e12.5,//10x,'These volume totals may be incomplete because of missing',/10x,'partial molar volume data in the data base.',/)

    if (qftpr2) then
        write (noutpt,1000)
        go to 999
    end if

    ! Print saturation index and affinity tables for reactions in
    ! the aqueous phase not constrained to be at equilibrium. The
    ! data correspond to the reactions and thermodynamic data in
    ! the 'd' set.
    call prtsia(affsd,jflagd,jflgi,jsflag,narn1,narn2,nbasp,nbaspd,nbt,nbtmax,ndrsd,ndrsmx,ndrsrd,nhydr,noutpt,nrdxsp,nstmax,sidrsp,uspec)

    if (iopr(7) .ge. 0) then
        ! Print saturation index and affinity tables for the
        ! various non-aqueous phases.
        call prtsat(affpd,iern1,iern2,ilrn1,ilrn2,imrn1,imrn2,iopr,iopt,ixrn1,ixrn2,jpflag,noutpt,noprmx,noptmx,nptmax,sidrph,tolsat,uphase)
    end if

    if (iopt(4) .gt. 0) then
        ! Print tables of the compositions and affinities of solid
        ! solution phases present in the equilbrium system (ES).
        write (noutpt,1900)
1900 format(//11x,'--- Solid Solution Product Phases ---',/)

        ncount = 0

        do np = ixrn1,ixrn2
            if (jpflag(np) .eq. -1) then
                ncount = ncount + 1
                call prtsso(acflg,actlg,affpd,affsd,ixrn1,ixrn2,jsol,jsomax,ncmpr,noutpt,np,nptmax,nstmax,nxtmax,sidrph,sidrsp,tolsat,uspec,uphase,uxtype,xbar,xbarlg)
            end if
        end do

        if (ncount .le. 0) then
            write (noutpt,1910)
        end if

1910 format(1x,'None')
    end if

    if (net .gt. 0) then
        ! Print tables of the compositions and affinities of generic
        ! ion exchange phases present in the equilbrium system (ES).
        write (noutpt,1930)
1930 format(//11x,'--- Generic Ion Exchange Phases ---',/)

        tolspf = tolsat
        ncount = 0

        do np = iern1,iern2
            ncount = ncount + 1

            call prtgex(acflg,actlg,affpd,affsd,cegexs,conc,egexjc,egexjf,egexpa,egexpc,egexs,egexw,iern1,iern2,ietmax,jern1,jern2,jetmax,jgext,kern1,kern2,ketmax,kgexsa,moph,mosp,netmax,ngexsa,ngext,noutpt,np,nptmax,nstmax,sidrph,sidrsp,tolspf,ugexj,ugexmo,uspec,uphase,xbar,xbarlg,xgexw,wkgwi)
        end do

        if (ncount .le. 0) then
            write (noutpt,1940)
        end if

1940 format(1x,'None')
    end if

    if (iopr(8) .ge. 0) then
        ! Print table of equilibrium fugacities.
        call prtfug(jgsort,fugac,fugalg,jsflag,ngrn1,ngt,ngtmax,noutpt,nstmax,uspec)
    end if

    write (noutpt,2000)
2000 format(/' - - - - - - - - - - - - - - - - - - - - - - - - -',' - - - - - - -',//)

999 continue
end subroutine scripz