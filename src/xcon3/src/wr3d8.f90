subroutine wr3d8(cgexj,cgexpi,covali,ehi,egexsi,fo2lgi,iebal3,ietmax,iodb,iopg,iopr,iopt,irdxc3,itdsf3,itermx,jgext,jetmax,jflgi,jgexti,jpres3,kxmod,mwtges,nbti,nbtmax,ncmpri,net,neti,netmax,newin,ngexti,ngexrt,nobswt,nodbmx,nopgmx,noprmx,noptmx,nsbswt,ntitl,ntitmx,nttyo,nxicmx,nxmdmx,nxmod,nxti,nxtimx,pei,press,qgexsh,rho,scamas,tdspkg,tdspl,tempc,tgexp,tolbt,toldl,tolspf,ucospi,uebal,ugexj,ugexji,ugexmo,ugexp,ugexpi,ugexr,ugexsi,umemi,uobsw,uredox,usbsw,usoli,uspeci,utitl,uhfgex,uvfgex,uxkgex,uxmod,xbari,xgexsi,xhfgex,xlkgex,xvfgex,xlkmod,zgexj)
    !! This subroutine writes the EQ3NR input file in menu-style ("D")
    !! format for version 8.0.
    !! The calling sequence of this subroutine is identical to that of
    !! XCON3/wr3w8.f.
    !! The calling sequence of this subroutine is identical to that of
    !! EQ3NR/rd3inw.f, EQ3NR/rd3ind.f, XCON3/rd3w8.f, and XCON3/rd3d8.f,
    !! except that newin is added and ninpts, nprob, noutpt, qend, and
    !! qrderr are deleted.
    !! This subroutine is called by:
    !!   XCON3/xcon3.f
    implicit none

    ! Calling sequence variable declarations.
    integer :: ietmax
    integer :: jetmax
    integer :: nbtmax
    integer :: netmax
    integer :: nodbmx
    integer :: nopgmx
    integer :: noprmx
    integer :: noptmx
    integer :: ntitmx
    integer :: nxmdmx
    integer :: nxicmx
    integer :: nxtimx

    integer :: newin
    integer :: nttyo

    integer :: iodb(nodbmx)
    integer :: iopg(nopgmx)
    integer :: iopr(noprmx)
    integer :: iopt(noptmx)
    integer :: jgext(netmax)
    integer :: jgexti(netmax)
    integer :: jflgi(nbtmax)
    integer :: kxmod(nxmdmx)
    integer :: ncmpri(2,nxtimx)
    integer :: ngexrt(jetmax,netmax)
    integer :: ngexti(jetmax,netmax)

    integer :: iebal3
    integer :: irdxc3
    integer :: itdsf3
    integer :: itermx
    integer :: jpres3
    integer :: nbti
    integer :: net
    integer :: neti
    integer :: nobswt
    integer :: nsbswt
    integer :: ntitl
    integer :: nxmod
    integer :: nxti

    logical :: qgexsh

    character(len=80) :: utitl(ntitmx)
    character(len=56) :: ugexr(ietmax,jetmax,netmax)
    character(len=48) :: ucospi(nbtmax)
    character(len=48) :: uobsw(2,nbtmax)
    character(len=48) :: usbsw(2,nbtmax)
    character(len=48) :: uspeci(nbtmax)
    character(len=48) :: uxmod(nxmdmx)
    character(len=24) :: ugexmo(netmax)
    character(len=24) :: ugexp(netmax)
    character(len=24) :: ugexpi(netmax)
    character(len=24) :: ugexsi(ietmax,jetmax,netmax)
    character(len=24) :: umemi(nxicmx)
    character(len=24) :: usoli(nxtimx)
    character(len=24) :: uebal
    character(len=24) :: uredox
    character(len=8) :: ugexj(jetmax,netmax)
    character(len=8) :: ugexji(jetmax,netmax)
    character(len=8) :: uhfgex(ietmax,jetmax,netmax)
    character(len=8) :: uvfgex(ietmax,jetmax,netmax)
    character(len=8) :: uxkgex(ietmax,jetmax,netmax)

    real(kind=8) :: cgexj(jetmax,netmax)
    real(kind=8) :: cgexpi(netmax)
    real(kind=8) :: covali(nbtmax)
    real(kind=8) :: egexsi(ietmax,jetmax,netmax)
    real(kind=8) :: mwtges(netmax)
    real(kind=8) :: tgexp(netmax)
    real(kind=8) :: xbari(nxicmx)
    real(kind=8) :: xhfgex(ietmax,jetmax,netmax)
    real(kind=8) :: xlkgex(ietmax,jetmax,netmax)
    real(kind=8) :: xvfgex(ietmax,jetmax,netmax)
    real(kind=8) :: xlkmod(nxmdmx)
    real(kind=8) :: zgexj(jetmax,netmax)
    real(kind=8) :: xgexsi(ietmax,jetmax,netmax)

    real(kind=8) :: ehi
    real(kind=8) :: fo2lgi
    real(kind=8) :: pei
    real(kind=8) :: press
    real(kind=8) :: rho
    real(kind=8) :: scamas
    real(kind=8) :: tdspkg
    real(kind=8) :: tdspl
    real(kind=8) :: tempc
    real(kind=8) :: tolbt
    real(kind=8) :: toldl
    real(kind=8) :: tolspf

    include 'eqlib/eqlj8.h'
    include 'eqlib/eqlo8.h'

    ! Local variable declarations.
    integer :: i
    integer :: iei
    integer :: ij
    integer :: j
    integer :: je
    integer :: jei
    integer :: jfl
    integer :: jj
    integer :: j2
    integer :: j3
    integer :: j4
    integer :: j5
    integer :: j6
    integer :: k
    integer :: kxmd
    integer :: n
    integer :: nb
    integer :: ne
    integer :: nei
    integer :: nr1
    integer :: nr2
    integer :: nxi
    integer :: nxic

    integer :: ilnobl

    logical :: qgexbf
    logical :: qgexef
    logical :: qstop

    character(len=80) :: ux80
    character(len=48) :: ux48
    character(len=48) :: ux48a
    character(len=48) :: ux48b
    character(len=24) :: ux24
    character(len=16) :: ux16
    character(len=8) :: uxn
    character(len=8) :: uxo
    character(len=8) :: ux8
    character(len=1) :: ux1(5)
    character(len=1) :: uxs

    real(kind=8) :: xx

    ! Check some dimensioning parameters.
    qstop = .false.

    if (ndbxpa .ne. nodbmx) then
        write (nttyo,3000) ndbxpa,nodbmx
3000 format(/' * Error - (XCON3/wr3d8) The dimensioning parameter',' for the',/7x,'number of iodb debugging print option switches',' with string definitions',/7x,'(ndbxpa) has a value of ',i3,', but the dimensioning',/7x,'parameter for the number',' of such switches (nodbpa) has a',/7x,'value of ',i3,'.')

        qstop = .true.
    end if

    if (npgxpa .ne. nopgmx) then
        write (nttyo,3010) npgxpa,nopgmx
3010 format(/' * Error - (XCON3/wr3d8) The dimensioning parameter',' for the',/7x,'number of iopg activity coefficient option',' switches with string definitions',/7x,'(npgxpa) has a value',' of ' ,i3,', but the dimensioning',/7x,'parameter for the',' number of such switches (nopgpa) has a',/7x,'value of ',i3,'.')

        qstop = .true.
    end if

    if (nprxpa .ne. noprmx) then
        write (nttyo,3020) nprxpa,noprmx
3020 format(/' * Error - (XCON3/wr3d8) The dimensioning parameter',' for the',/7x,'number of iopr print option switches',' with string definitions',/7x,'(nprxpa) has a value of ',i3,', but the dimensioning',/7x,'parameter for the number',' of such switches (noprpa) has a',/7x,'value of ',i3,'.')

        qstop = .true.
    end if

    if (nptxpa .ne. noptmx) then
        write (nttyo,3030) nptxpa,noptmx
3030 format(/' * Error - (XCON3/wr3d8) The dimensioning parameter',' for the',/7x,'number of iopt model option switches',' with string definitions',/7x,'(nptxpa) has a value of ',i3,', but the dimensioning',/7x,'parameter for the number',' of such switches (noptpa) has a',/7x,'value of ',i3,'.')

        qstop = .true.
    end if

    if (qstop) then
        stop
    end if

    ! Title.
    write (newin,2050)
2050 format('|',78('-'),'|')

    write (newin,2040)
2040 format('| Title',18x,'| (utitl(n))',42x,'|')

    write (newin,2050)

    do n = 1,ntitl
        write (newin,2070) utitl(n)
2070 format('|',a78,'|')
    end do

    write (newin,2050)

    ! Special basis switches.
    write (newin,2100)
2100 format('|Special Basis Switches (for model definition only)',7x,'| (nsbswt)',11x,'|')

    write (newin,2050)

    do n = 1,nsbswt
        write (newin,2120) usbsw(1,n),usbsw(2,n)
2120 format('|Replace |',a48,'| (usbsw(1,n))',7x,'|',/'|   with |',a48,'| (usbsw(2,n))',7x,'|')

        write (newin,2050)
    end do

    if (nsbswt .le. 0) then
        ux48a = 'None'
        ux48b = 'None'
        write (newin,2120) ux48a,ux48b
        write (newin,2050)
    end if

    ! Temperature.
    write (newin,2140) tempc
2140 format('|Temperature (C)         |',1pe12.5,'| (tempc)',t80,'|')

    write (newin,2050)

    ! Pressure.
    ux1(1) = ' '
    ux1(2) = ' '
    ux1(3) = ' '

    if (jpres3 .le. 0) then
        press = 0.
        ux1(1) = 'x'
    else if (jpres3 .eq. 1) then
        press = 0.
        ux1(2) = 'x'
    else if (jpres3 .ge. 2) then
        ux1(3) = 'x'
    end if

    write (newin,2160) ux1(1),ux1(2),ux1(3),press
2160 format('|Pressure option (jpres3):',t80,'|',/'|  [',a1,'] ( 0) Data file reference curve value',t80,'|',/'|  [',a1,'] ( 1) 1.013-bar/steam-saturation curve value',t80,'|',/'|  [',a1,'] ( 2) Value (bars) |',1pe12.5,'| (press)',t80,'|')

    write (newin,2050)

    ! Density.
    write (newin,2180) rho
2180 format('|Density (g/cm3)         |',1pe12.5,'| (rho)',34x,'|')

    write (newin,2050)

    ! Total dissolved solutes (TDS).
    ux1(1) = ' '
    ux1(2) = ' '

    if (itdsf3 .le. 0) then
        ux1(1) = 'x'
    else if (itdsf3 .ge. 1) then
        ux1(2) = 'x'
    end if

    write (newin,2200) ux1(1),tdspkg,ux1(2),tdspl
2200 format('|Total dissolved solutes option (itdsf3):',t80,'|',/'|  [',a1,'] ( 0) Value (mg/kg.sol) |',1pe12.5,'| (tdspkg)',t80,'|',/'|  [',a1,'] ( 1) Value (mg/L)      |',e12.5,'| (tdspl)',t80,'|')

    write (newin,2050)

    ! Species to adjust for electrical balancing.
    ux1(1) = ' '
    ux1(2) = ' '

    if (iebal3 .le. 0) then
        ux1(1) = 'x'
    else
        ux1(2) = 'x'
    end if

    write (newin,2220) ux1(1),ux1(2),uebal
2220 format('|Electrical balancing option (iebal3):',t80,'|',/'|  [',a1,'] ( 0) No balancing is done',t80,'|',/'|  [',a1,'] ( 1) Balance on species |',a24,'| (uebal)',t80,'|')

    write (newin,2050)

    ! Default redox constraint.
    ux1(1) = ' '
    ux1(2) = ' '
    ux1(3) = ' '
    ux1(4) = ' '
    ux1(5) = ' '

    if (irdxc3 .le. -3) then
        ux1(1) = 'x'
    else if (irdxc3 .eq. -2) then
        ux1(2) = 'x'
    else if (irdxc3 .eq. -1) then
        ux1(3) = 'x'
    else if (irdxc3 .eq. 0) then
        ux1(4) = 'x'
    else if (irdxc3 .ge. 1) then
        ux1(5) = 'x'
    end if

    write (newin,2250) ux1(1),ux1(2),pei,ux1(3),ehi,ux1(4),fo2lgi,ux1(5),uredox
2250 format('|Default redox constraint (irdxc3):',t80,'|',/'|  [',a1,'] (-3) Use O2(g) line in the aqueous basis species',' block',t80,'|',/'|  [',a1,'] (-2) pe (pe units)      |',1pe12.5,'| (pei)',t80,'|',/'|  [',a1,'] (-1) Eh (volts)         |',e12.5,'| (ehi)',t80,'|',/'|  [',a1,'] ( 0) Log fO2 (log bars) |',e12.5,'| (fo2lgi)',t80,'|',/'|  [',a1,'] ( 1) Couple (aux. sp.)  |',a24,'| (uredox)',t80,'|')

    write (newin,2050)

    ! Aqueous basis species and associated constraints.
    write (newin,2300)
2300 format('|Aqueous Basis Species/Constraint Species',8x,'|Conc., etc. |Units/Constraint|',/'| (uspeci(n)/ucospi(n))',26x,'| (covali(n))|(ujf3(jflgi(n)))|')

    write (newin,2050)

    qstop = .false.

    do nb = 1,nbti
        jfl = jflgi(nb)

        if (jfl .gt. njfxpa) then
            j2 = ilnobl(uspeci(nb))
            write (nttyo,3050) jfl,uspeci(nb)(1:j2),njfxpa
3050 format(/' * Error - (XCON3/wr3d8) The jflgi value of ',i3,' for',/7x,'the species ',a,' exceeds the value of ',i3,' for the dimensioning',/7x,'parameter (njfxpa) of the',' associated string array.')

            qstop = .true.
        end if

        write (newin,1130) uspeci(nb),covali(nb),ujf3(jfl)
1130 format('|',a48,'|',1pe12.5,'|',a16,'|')

        if (jfl.eq.17 .or. jfl.eq.18 .or. jfl.eq.25) then
            write (newin,1160) ucospi(nb)
1160 format('|->|',a48,'| (ucospi(n))',14x,'|')
        end if
    end do

    if (qstop) then
        stop
    end if

    if (nbti .le. 0) then
        ux48 = 'None'
        ux16 = 'None'
        xx = 0.
        write (newin,1130) ux48,xx,ux16
    end if

    write (newin,2050)

    write (newin,1162)
1162 format('* Valid jflag strings (ujf3(jflgi(n))) are:',t80,'*')

    ux80(1:1) = '*'
    ux80(2:79) = ' '
    ux80(80:80) = '*'
    j3 = 6
    k = 0

    do jfl = -1,njfxpa
        if (ujf3(jfl)(1:6) .ne. 'ERROR ') then
            jj = j3 + 19
            ux80(j3:jj) = ujf3(jfl)
            j3 = jj + 1
            k = k + 1

            if (k .ge. 3) then
                write (newin,1164) ux80(1:80)
1164 format(a)

                ux80(2:79) = ' '
                j3 = 6
                k = 0
            end if
        end if
    end do

    if (k .gt. 0) then
        write (newin,1164) ux80(1:80)
    end if

    write (newin,2052)
2052 format('*',78('-'),'*')

    ! Ion exchanger creation.
    write (newin,1600)
1600 format('|Create Ion Exchangers',2x,'| (net)',48x,'|')

    write (newin,2050)

    ! Set the advisory for exchanger creation blocks to follow or not.
    qgexbf = qgexsh .or. net.gt.0

    if (qgexbf) then
        write (newin,1603)
1603 format('|Advisory: at least one exchanger creation block',' follows on this file.',t80,'|')
    else
        write (newin,1605)
1605 format('|Advisory: no exchanger creation blocks follow',' on this file.',t80,'|')
    end if

    ! Flag: on further processing, show at least one such block.
    ux1(1) = ' '

    if (qgexsh) then
        ux1(1) = 'x'
    end if

    write (newin,1607) ux1(1)
1607 format('|Option: on further processing (writing a PICKUP file',' or running XCON3 on the',t80,'|',/'|present file), force the',' inclusion of at least one such block (qgexsh):',t80,'|',/'|  [',a1,'] (.true.)',t80,'|')

    write (newin,2050)

    do ne = 1,net
        write (newin,1610) ugexp(ne)
1610 format('|Exchanger phase |',a24,'| (ugexp(n))',25x,'|')

        write (newin,2050)
        write (newin,1620) mwtges(ne)
1620 format('|->|Mol. wt. (Z)    |',1pg12.5,'| (mwtges(n))',33x,'|')

        write (newin,2050)
        write (newin,1630) ugexmo(ne)
1630 format('|->|Exchange model  |',a24,'| (ugexmo(n))',21x,'|')

        write (newin,2050)
        write (newin,1640) tgexp(ne)
1640 format('|->|Ref. temp. (C)  |',1pg12.5,'| (tgexp(n))',34x,'|')

        write (newin,2050)

        do je = 1,jgext(ne)
            write (newin,1650) ugexj(je,ne)
1650 format('|->|Exchange site   |',a8,'| (ugexj(j,n))',36x,'|')

            write (newin,2050)
            write (newin,1660) cgexj(je,ne)
1660 format('|--->|Stoich. number  |',1pg12.5,'| (cgexj(j,n))',30x,'|')

            write (newin,2050)
            write (newin,1670) zgexj(je,ne)
1670 format('|--->|Electr. charge  |',1pg12.5,'| (zgexj(j,n))',30x,'|')

            write (newin,2050)

            do n = 1,ngexrt(je,ne)
                write (newin,1680) ugexr(n,je,ne)
1680 format('|--->|Reaction        |',a24,'| (ugexr(i,j,n)',17x,'|')

                write (newin,2050)
                write (newin,1862)
1862 format('|----->|Parameter |Value',7x,'|Units   | (this',' is a table header)',13x,'|')

                write (newin,2050)
                write (newin,1684) xlkgex(n,je,ne),uxkgex(n,je,ne)
1684 format('|----->|K func.   |',1pe12.5,'|',a8,'| (xlkgex(i,j,n), uxkgex(i,j,n))',7x,'|')

                write (newin,1686) xhfgex(n,je,ne),uhfgex(n,je,ne)
1686 format('|----->|DelH0r    |',1pe12.5,'|',a8,'| (xhfgex(i,j,n), uhfgex(i,j,n))',7x,'|')

                write (newin,1688) xvfgex(n,je,ne),uvfgex(n,je,ne)
1688 format('|----->|DelV0r    |',1pe12.5,'|',a8,'| (xvfgex(i,j,n), uvfgex(i,j,n))',7x,'|')

                write (newin,2050)
            end do
        end do
    end do

    if (net.le.0 .and. qgexsh) then
        ux24 = 'None'
        ux8 = 'None'
        xx = 0.
        write (newin,1610) ux24
        write (newin,2050)
        write (newin,1620) xx
        write (newin,2050)
        write (newin,1630) ux24
        write (newin,2050)
        write (newin,1640) xx
        write (newin,2050)
        write (newin,1650) ux8
        write (newin,2050)
        write (newin,1660) xx
        write (newin,2050)
        write (newin,1670) xx
        write (newin,2050)
        write (newin,1680) ux24
        write (newin,2050)
        write (newin,1862)
        write (newin,2050)
        write (newin,1684) xx,ux8
        write (newin,1686) xx,ux8
        write (newin,1688) xx,ux8
        write (newin,2050)
    end if

    if (net.gt.0 .or. qgexsh) then
        write (newin,1694)
1694 format('* Valid exchange model strings (ugexmo(n)) are:',t80,'*')

        ux80(1:1) = '*'
        ux80(2:79) = ' '
        ux80(80:80) = '*'
        j3 = 6
        k = 0

        do n = 1,nvet_par
            if (ugexmv(n)(1:6) .ne. 'ERROR ') then
                jj = j3 + 27
                ux80(j3:jj) = ugexmv(n)
                j3 = jj + 1
                k = k + 1

                if (k .ge. 2) then
                    write (newin,1696) ux80(1:80)
1696 format(a)

                    ux80(2:79) = ' '
                    j3 = 6
                    k = 0
                end if
            end if
        end do

        if (k .gt. 0) then
            write (newin,1696) ux80(1:80)
        end if

        write (newin,2052)

        write (newin,1690)
1690 format('* Valid units strings (uxkgex(i,j,n)/uhfgex(i,j,n)/','uvfgex(i,j,n)) are:',t80,'*')

        ux80(1:1) = '*'
        ux80(2:79) = ' '
        ux80(80:80) = '*'
        j3 = 6
        k = 0

        do n = 1,4
            jj = j3 + 19
            ux80(j3:jj) = uxfuni(n)
            j3 = jj + 1
            k = k + 1

            if (k .ge. 3) then
                write (newin,1692) ux80(1:80)
1692 format(a)

                ux80(2:79) = ' '
                j3 = 6
                k = 0
            end if
        end do

        if (k .gt. 0) then
            write (newin,1692) ux80(1:80)
        end if

        write (newin,2052)
    end if

    ! Specified generic ion exchanger compositions.
    write (newin,1700)
1700 format('|Ion Exchanger Compositions      | (neti)',38x,'|')

    write (newin,2050)

    do nei = 1,neti
        ! Find the corresponding ne index.
        qgexef = .true.

        do ne = 1,net
            j2 = ilnobl(ugexp(ne))
            j3 = ilnobl(ugexpi(nei))

            if (ugexp(nei)(1:j3) .eq. ugexp(ne)(1:j2)) then
                j2 = ilnobl(ugexmo(ne))

                if (ugexmo(ne)(1:j2) .eq. 'Gapon' .or.        ugexmo(ne)(1:6) .eq. 'Gapon-' .or.        ugexmo(ne)(1:j2) .eq. 'Vanselow' .or.        ugexmo(ne)(1:9) .eq. 'Vanselow-') then
                    ! Input composition is described in terms of equivalent
                    ! fractions on the sites.
                    qgexef = .true.
                else
                    ! Input composition is described in terms of mole fractions
                    ! on the sites.
                    qgexef = .false.
                end if

                go to 240
            end if
        end do

240 continue

        write (newin,1710) ugexpi(nei)
1710 format('|Exchanger phase |',a24,'| (ugexpi(n))',24x,'|')

        write (newin,2050)
        write (newin,1720) cgexpi(nei)
1720 format('|->|Moles/kg.H2O    |',1pg12.5,'| (cgexpi(n))',33x,'|')

        write (newin,2050)

        do jei = 1,jgexti(nei)
            write (newin,1730) ugexji(jei,nei)
1730 format('|->|Exchange site   |',a8,'| (ugexji(j,n))',35x,'|')

            write (newin,2050)

            if (qgexef) then
                write (newin,1740)
1740 format('|--->|Exchange species',8x,'|Eq. frac.   | (this is a table header)',10x,'|')

                write (newin,2050)

                do iei = 1,ngexti(jei,nei)
                    write (newin,1750) ugexsi(iei,jei,nei),egexsi(iei,jei,nei)
1750 format('|--->|',a24,'|',1pe12.5,'| (ugexsi(i,j,n), egexsi(i,j,n))',4x,'|')
                end do
            else
                write (newin,1760)
1760 format('|--->|Exchange species',8x,'|Mole frac.  | (this is a table header)',10x,'|')

                write (newin,2050)

                do iei = 1,ngexti(jei,nei)
                    write (newin,1770) ugexsi(iei,jei,nei),xgexsi(iei,jei,nei)
1770 format('|--->|',a24,'|',1pe12.5,'| (ugexsi(i,j,n), xgexsi(i,j,n))',4x,'|')
                end do
            end if

            write (newin,2050)
        end do
    end do

    if (neti .le. 0) then
        ux24 = 'None'
        ux8 = 'None'
        xx = 0.
        write (newin,1710) ux24
        write (newin,2050)
        write (newin,1720) xx
        write (newin,2050)
        write (newin,1730) ux8
        write (newin,2050)
        write (newin,1740)
        write (newin,2050)
        write (newin,1750) ux24,xx
        write (newin,2050)
    end if

    ! Solid solutions.
    write (newin,1800)
1800 format('|Solid Solution Compositions     | (nxti)',38x,'|')

    write (newin,2050)

    do nxi = 1,nxti
        write (newin,1810) usoli(nxi)
1810 format('|Solid Solution',10x,'|',a24,'| (usoli(n))',17x,'|')

        write (newin,2050)
        nr1 = ncmpri(1,nxi)
        nr2 = ncmpri(2,nxi)
        write (newin,1812)
1812 format('|->|Component',15x,'|Mole frac.  | (this is a',' table header)',12x,'|')

        write (newin,2050)

        do nxic = nr1,nr2
            write (newin,1815) umemi(nxic),xbari(nxic)
1815 format('|->|',a24,'|',1pe12.5,'| (umemi(i,n), xbari(i,n))',12x,'|')
        end do

        write (newin,2050)
    end do

    if (nxti .le. 0) then
        ux24 = 'None'
        xx = 0.
        write (newin,1810) ux24
        write (newin,2050)
        write (newin,1812)
        write (newin,2050)
        write (newin,1815) ux24,xx
        write (newin,2050)
    end if

    ! Nxmod options.
    write (newin,3080)
3080 format('|Alter/Suppress Options  | (nxmod)',45x,'|')

    write (newin,2050)
    write (newin,3082)
3082 format('|Species',41x,'|Option',10x,'|Alter value |',/'| (uxmod(n))',37x,'|(ukxm(kxmod(n)))| (xlkmod(n))|')

    write (newin,2050)

    do n = 1,nxmod
        write (newin,3090) uxmod(n),ukxm(kxmod(n)),xlkmod(n)
3090 format('|',a48,'|',a16,'|',1pe12.5,'|')
    end do

    if (nxmod .le. 0) then
        ux48 = 'None'
        ux16 = 'None'
        xx = 0.
        write (newin,3090) ux48,ux16,xx
    end if

    write (newin,2050)

    write (newin,3092)
3092 format('* Valid alter/suppress strings (ukxm(kxmod(n))) are:',t80,'*')

    ux80(1:1) = '*'
    ux80(2:79) = ' '
    ux80(80:80) = '*'
    j3 = 6
    k = 0

    do kxmd = -1,2
        jj = j3 + 19
        ux80(j3:jj) = ukxm(kxmd)
        j3 = jj + 1
        k = k + 1

        if (k .ge. 3) then
            write (newin,3094) ux80(1:80)
3094 format(a)

            ux80(2:79) = ' '
            j3 = 6
            k = 0
        end if
    end do

    if (k .gt. 0) then
        write (newin,3094) ux80(1:80)
    end if

    write (newin,2052)

    ! Iopt model option switches.
    ! Note: iopt(1) = iopt1, etc.
    write (newin,1460)
1460 format('|Iopt Model Option Switches',' ("( 0)" marks default choices)',t80,'|')

    write (newin,2050)

    uxo = 'iopt('

    do n = 1,noptmx
        j6 = index(uoptcx(n),'EQ3NR')

        if (j6 .gt. 0) then
            uxn = ' '
            write (uxn,'(i2)') n
            call lejust(uxn)
            j5 = ilnobl(uxn)
            uxo(6:8) = uxn(1:j5)
            j4 = ilnobl(uxo) + 1
            uxo(j4:j4) = ')'
            j3 = ilnobl(uopttx(n))
            write (newin,1465) uxo(1:j4),uopttx(n)(1:j3)
1465 format('|',a,' - ',a,':',t80,'|')

            i = iopt(n)

            do j = 1,jptxpa
                ij = ioptox(j,n)
                uxs = ' '

                if (ij .eq. i) then
                    uxs = 'x'
                end if

                j2 = ilnobl(uoptox(j,n))

                if (j2 .gt. 0) then
                    write (newin,1470) uxs,ij,uoptox(j,n)(1:j2)
1470 format('|  [',a1,'] (',i2,') ',a,t80,'|')
                end if
            end do

            write (newin,2050)
        end if
    end do

    ! Iopg activity coefficient option switches.
    ! Note: iopg(1) = iopg1, etc.
    write (newin,1480)
1480 format('|Iopg Activity Coefficient Option Switches',' ("( 0)" marks default choices)',t80,'|')

    write (newin,2050)

    uxo = 'iopg('

    do n = 1,nopgmx
        j6 = index(uopgcx(n),'EQ3NR')

        if (j6 .gt. 0) then
            uxn = ' '
            write (uxn,'(i2)') n
            call lejust(uxn)
            j5 = ilnobl(uxn)
            uxo(6:8) = uxn(1:j5)
            j4 = ilnobl(uxo) + 1
            uxo(j4:j4) = ')'
            j3 = ilnobl(uopgtx(n))
            write (newin,1465) uxo(1:j4),uopgtx(n)(1:j3)
            i = iopg(n)

            do j = 1,jpgxpa
                ij = iopgox(j,n)
                uxs = ' '

                if (ij .eq. i) then
                    uxs = 'x'
                end if

                j2 = ilnobl(uopgox(j,n))

                if (j2 .gt. 0) then
                    write (newin,1470) uxs,ij,uopgox(j,n)(1:j2)
                end if
            end do

            write (newin,2050)
        end if
    end do

    ! Iopr print option switches.
    ! Note: iopr(1) = iopr1, etc.
    write (newin,1490)
1490 format('|Iopr Print Option Switches',' ("( 0)" marks default choices)',t80,'|')

    write (newin,2050)

    uxo = 'iopr('

    do n = 1,noprmx
        j6 = index(uoprcx(n),'EQ3NR')

        if (j6 .gt. 0) then
            uxn = ' '
            write (uxn,'(i2)') n
            call lejust(uxn)
            j5 = ilnobl(uxn)
            uxo(6:8) = uxn(1:j5)
            j4 = ilnobl(uxo) + 1
            uxo(j4:j4) = ')'
            j3 = ilnobl(uoprtx(n))
            write (newin,1465) uxo(1:j4),uoprtx(n)(1:j3)
            i = iopr(n)

            do j = 1,jprxpa
                ij = ioprox(j,n)
                uxs = ' '

                if (ij .eq. i) then
                    uxs = 'x'
                end if

                j2 = ilnobl(uoprox(j,n))

                if (j2 .gt. 0) then
                    write (newin,1470) uxs,ij,uoprox(j,n)(1:j2)
                end if
            end do

            write (newin,2050)
        end if
    end do

    ! Iodb debugging print option switches.
    ! Note: iodb(1) = iodb1, etc.
    write (newin,1495)
1495 format('|Iodb Debugging Print Option Switches',' ("( 0)" marks default choices)',t80,'|')

    write (newin,2050)

    uxo = 'iodb('

    do n = 1,nodbmx
        j6 = index(uodbcx(n),'EQ3NR')

        if (j6 .gt. 0) then
            uxn = ' '
            write (uxn,'(i2)') n
            call lejust(uxn)
            j5 = ilnobl(uxn)
            uxo(6:8) = uxn(1:j5)
            j4 = ilnobl(uxo) + 1
            uxo(j4:j4) = ')'
            j3 = ilnobl(uodbtx(n))
            write (newin,1465) uxo(1:j4),uodbtx(n)(1:j3)
            i = iodb(n)

            do j = 1,jdbxpa
                ij = iodbox(j,n)
                uxs = ' '

                if (ij .eq. i) then
                    uxs = 'x'
                end if

                j2 = ilnobl(uodbox(j,n))

                if (j2 .gt. 0) then
                    write (newin,1470) uxs,ij,uodbox(j,n)(1:j2)
                end if
            end do

            write (newin,2050)
        end if
    end do

    ! Numerical Parameters.
    write (newin,1500)
1500 format('|Numerical Parameters',t80,'|')

    write (newin,2050)

    write (newin,1510) tolbt
1510 format('| Beta convergence tolerance',t35,'|',1pe12.5,'| (tolbt)',t80,'|')

    write (newin,1520) toldl
1520 format('| Del convergence tolerance',t35,'|',1pe12.5,'| (toldl)',t80,'|')

    write (newin,1524) itermx
1524 format('| Max. Number of N-R Iterations',t35,'|',1x,i3,8x,'| (itermx)',t80,'|')

    write (newin,2050)

    ! Ordinary basis switches.
    write (newin,2700)
2700 format('|Ordinary Basis Switches (for numerical purposes only)',4x,'| (nobswt)',11x,'|')

    write (newin,2050)

    do n = 1,nobswt
        write (newin,2710) uobsw(1,n),uobsw(2,n)
2710 format('|Replace |',a48,'| (uobsw(1,n))',7x,'|',/'|   with |',a48,'| (uobsw(2,n))',7x,'|')

        write (newin,2050)
    end do

    if (nobswt .le. 0) then
        ux48a = 'None'
        ux48b = 'None'
        write (newin,2710) ux48a,ux48b
        write (newin,2050)
    end if

    ! Saturation flag tolerance.
    write (newin,2720) tolspf
2720 format('|Sat. flag tolerance     |',1pe12.5,'| (tolspf)',31x,'|')

    write (newin,2050)

    ! Scale factor for the mass of aqueous solution to write
    ! on the PICKUP file.
    write (newin,2790) scamas
2790 format('|Aq. Phase Scale Factor  |',1pe12.5,'| (scamas)',31x,'|')

    write (newin,2050)

    ! End of problem.
    write (newin,2900)
2900 format('|End of problem',64x,'|')

    write (newin,2050)

999 continue
end subroutine wr3d8