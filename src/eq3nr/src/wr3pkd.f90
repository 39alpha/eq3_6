subroutine wr3pkd(electr,cgexj,ietmax,iopg,jetmax,jflgi,jgext,kbt,kct,kdim,kmax,kmt,kxmod,kxt,mtbi,mtbaqi,mwtges,nbti,nbtmax,net,netmax,newin,ngexrt,nobswt,nopgmx,nsbswt,ntitl2,ntitmx,nxmdmx,nxmod,qgexsh,pressi,tempci,tgexp,ubmtbi,ugexj,ugexmo,ugexp,ugexr,uhfgex,uobsw,usbsw,utitl2,uvfgex,uxkgex,uxmod,uzveci,xhfgex,xlkgex,xlkmod,xvfgex,zgexj,zvclgi)
    !! This subroutine writes the pickup file in menu-style ("D") format.
    !! This file is used to communicate data from EQ3NR to EQ6 (it
    !! comprises the bottom half of an EQ6 input file).
    !! This subroutine is a near clone of parts of XCON6/wr6d8.f and
    !! EQ6/wr6pkd.f.
    !! The calling sequence of this subroutine is identical to that of
    !! EQ3NR/wr3pkw.f.
    !! This subroutine is called by:
    !!   EQ3NR/eq3nr.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: ietmax
    integer :: jetmax
    integer :: kmax
    integer :: nbtmax
    integer :: netmax
    integer :: nopgmx
    integer :: ntitmx
    integer :: nxmdmx

    integer :: newin

    integer :: iopg(nopgmx)
    integer :: jgext(netmax)
    integer :: jflgi(nbtmax)
    integer :: kxmod(nxmdmx)
    integer :: ngexrt(jetmax,netmax)

    integer :: kbt
    integer :: kct
    integer :: kdim
    integer :: kmt
    integer :: kxt
    integer :: nbti
    integer :: net
    integer :: nobswt
    integer :: nsbswt
    integer :: ntitl2
    integer :: nxmod

    logical :: qgexsh

    character(len=80) :: utitl2(ntitmx)
    character(len=56) :: ugexr(ietmax,jetmax,netmax)
    character(len=48) :: ubmtbi(nbtmax)
    character(len=48) :: uobsw(2,nbtmax)
    character(len=48) :: usbsw(2,nbtmax)
    character(len=48) :: uxmod(nxmdmx)
    character(len=48) :: uzveci(kmax)
    character(len=24) :: ugexmo(netmax)
    character(len=24) :: ugexp(netmax)
    character(len=8) :: ugexj(jetmax,netmax)
    character(len=8) :: uhfgex(ietmax,jetmax,netmax)
    character(len=8) :: uvfgex(ietmax,jetmax,netmax)
    character(len=8) :: uxkgex(ietmax,jetmax,netmax)

    real(kind=8) :: cgexj(jetmax,netmax)
    real(kind=8) :: mtbaqi(nbtmax)
    real(kind=8) :: mtbi(nbtmax)
    real(kind=8) :: mwtges(netmax)
    real(kind=8) :: tgexp(netmax)
    real(kind=8) :: xhfgex(ietmax,jetmax,netmax)
    real(kind=8) :: xlkgex(ietmax,jetmax,netmax)
    real(kind=8) :: xlkmod(nxmdmx)
    real(kind=8) :: xvfgex(ietmax,jetmax,netmax)
    real(kind=8) :: zgexj(jetmax,netmax)
    real(kind=8) :: zvclgi(kmax)

    real(kind=8) :: electr
    real(kind=8) :: pressi
    real(kind=8) :: tempci

    include 'eqlib/eqlk8.h'
    include 'eqlib/eqlo8.h'

    ! Local variable declarations.
    integer :: i
    integer :: ij
    integer :: j
    integer :: je
    integer :: jfl
    integer :: jfli
    integer :: jj
    integer :: j2
    integer :: j3
    integer :: j4
    integer :: j5
    integer :: j6
    integer :: k
    integer :: kprs
    integer :: krow
    integer :: kxmd
    integer :: n
    integer :: nbi
    integer :: ne

    integer :: ilnobl

    logical :: qgexbf
    logical :: qstop

    character(len=80) :: ux80
    character(len=48) :: ux48
    character(len=48) :: ux48a
    character(len=48) :: ux48b
    character(len=24) :: ux24
    character(len=16) :: ux16
    character(len=8) :: uxn
    character(len=8) :: uxo
    character(len=8) :: uxs
    character(len=8) :: ux8
    character(len=1) :: ux1(5)

    real(kind=8) :: xx

    ! Mark the start of the bottom half of the file.
    write (newin,7170)
7170 format('* Start of the bottom half of the input file',t80,'*')

    write (newin,1052)

    ! Secondary title.
    write (newin,7190)
7190 format('| Secondary Title',8x,'| (utitl2(n))',41x,'|')

    write (newin,1050)
1050 format('|',78('-'),'|')

    do n = 1,ntitl2
        write (newin,1090) utitl2(n)
1090 format('|',a78,'|')
    end do

    write (newin,1050)

    ! Special basis switches.
    write (newin,7210)
7210 format('|Special Basis Switches (for model definition only)',7x,'| (nsbswt)',11x,'|')

    write (newin,1050)

    do n = 1,nsbswt
        write (newin,7230) usbsw(1,n),usbsw(2,n)
7230 format('|Replace |',a48,'| (usbsw(1,n))',7x,'|',/'|   with |',a48,'| (usbsw(2,n))',7x,'|')

        write (newin,1050)
    end do

    if (nsbswt .le. 0) then
        ux48a = 'None'
        ux48b = 'None'
        write (newin,7230) ux48a,ux48b
        write (newin,1050)
    end if

    ! Original temperature.
    write (newin,7310) tempci
7310 format('|Original temperature (C) |',1pe12.5,'| (tempci)',t80,'|')

    write (newin,1050)

    ! Original pressure.
    write (newin,7330) pressi
7330 format('|Original pressure (bars) |',1pe12.5,'| (pressi)',t80,'|')

    write (newin,1050)

    ! Ion exchanger creation.
    write (newin,1600)
1600 format('|Create Ion Exchangers',2x,'| (net)',48x,'|')

    write (newin,1050)

    ! Set the advisory for exchanger creation blocks to follow or not.
    qgexbf = qgexsh .or. net.gt.0

    if (qgexbf) then
        write (newin,1603)
1603 format('|Advisory: at least one exchanger creation block',' follows on this file.',t80,'|')
    else
        write (newin,1605)
1605 format('|Advisory: no exchanger creation blocks follow on',' this file.',t80,'|')
    end if

    ! Flag: on further processing, show at least one such block.
    ux1(1) = ' '

    if (qgexsh) then
        ux1(1) = 'x'
    end if

    write (newin,1607) ux1(1)
1607 format('|Option: on further processing (writing a pickup file',' or running XCON6 on the',t80,'|',/'|present file), force the',' inclusion of at least one such block (qgexsh):',t80,'|',/'|  [',a1,'] (.true.)',t80,'|')

    write (newin,1050)

    do ne = 1,net
        write (newin,1610) ugexp(ne)
1610 format('|Exchanger phase |',a24,'| (ugexp(n))',25x,'|')

        write (newin,1050)
        write (newin,1620) mwtges(ne)
1620 format('|->|Mol. wt. (Z)    |',1pg12.5,'| (mwtges(n))',33x,'|')

        write (newin,1050)
        write (newin,1630) ugexmo(ne)
1630 format('|->|Exchange model  |',a24,'| (ugexmo(n))',21x,'|')

        write (newin,1050)
        write (newin,1640) tgexp(ne)
1640 format('|->|Ref. temp. (C)  |',1pg12.5,'| (tgexp(n))',34x,'|')

        write (newin,1050)

        do je = 1,jgext(ne)
            write (newin,1650) ugexj(je,ne)
1650 format('|->|Exchange site   |',a8,'| (ugexj(j,n))',36x,'|')

            write (newin,1050)
            write (newin,1660) cgexj(je,ne)
1660 format('|--->|Stoich. number  |',1pg12.5,'| (cgexj(j,n))',30x,'|')

            write (newin,1050)
            write (newin,1670) zgexj(je,ne)
1670 format('|--->|Electr. charge  |',1pg12.5,'| (zgexj(j,n))',30x,'|')

            write (newin,1050)

            do n = 1,ngexrt(je,ne)
                write (newin,1680) ugexr(n,je,ne)
1680 format('|--->|Reaction        |',a24,'| (ugexr(i,j,n)',17x,'|')

                write (newin,1050)
                write (newin,1682)
1682 format('|----->|Parameter |Value',7x,'|Units   | (this',' is a table header)',13x,'|')

                write (newin,1050)
                write (newin,1684) xlkgex(n,je,ne),uxkgex(n,je,ne)
1684 format('|----->|K func.   |',1pe12.5,'|',a8,'| (xlkgex(i,j,n), uxkgex(i,j,n))',7x,'|')

                write (newin,1686) xhfgex(n,je,ne),uhfgex(n,je,ne)
1686 format('|----->|DelH0r    |',1pe12.5,'|',a8,'| (xhfgex(i,j,n), uhfgex(i,j,n))',7x,'|')

                write (newin,1688) xvfgex(n,je,ne),uvfgex(n,je,ne)
1688 format('|----->|DelV0r    |',1pe12.5,'|',a8,'| (xvfgex(i,j,n), uvfgex(i,j,n))',7x,'|')

                write (newin,1050)
            end do
        end do
    end do

    if (net.le.0 .and. qgexsh) then
        ux24 = 'None'
        ux8 = 'None'
        xx = 0.
        write (newin,1610) ux24
        write (newin,1050)
        write (newin,1620) xx
        write (newin,1050)
        write (newin,1630) ux24
        write (newin,1050)
        write (newin,1640) xx
        write (newin,1050)
        write (newin,1650) ux8
        write (newin,1050)
        write (newin,1660) xx
        write (newin,1050)
        write (newin,1670) xx
        write (newin,1050)
        write (newin,1680) ux24
        write (newin,1050)
        write (newin,1682)
        write (newin,1050)
        write (newin,1684) xx,ux8
        write (newin,1686) xx,ux8
        write (newin,1688) xx,ux8
        write (newin,1050)
    end if

    if (net.gt.0 .or. qgexsh) then
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

        write (newin,1052)
1052 format('*',78('-'),'*')

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

        write (newin,1052)
    end if

    ! Nxmod options.
    write (newin,5080)
5080 format('|Alter/Suppress Options  | (nxmod)',45x,'|')

    write (newin,1050)
    write (newin,5082)
5082 format('|Species',41x,'|Option',10x,'|Alter value |',/'| (uxmod(n))',37x,'|(ukxm(kxmod(n)))| (xlkmod(n))|')

    write (newin,1050)

    do n = 1,nxmod
        write (newin,5090) uxmod(n),ukxm(kxmod(n)),xlkmod(n)
5090 format('|',a48,'|',a16,'|',1pe12.5,'|')
    end do

    if (nxmod .le. 0) then
        ux48 = 'None'
        ux16 = 'None'
        xx = 0.
        write (newin,5090) ux48,ux16,xx
    end if

    write (newin,1050)

    write (newin,5092)
5092 format('* Valid alter/suppress strings (ukxm(kxmod(n))) are:',t80,'*')

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
            write (newin,5094) ux80(1:80)
5094 format(a)

            ux80(2:79) = ' '
            j3 = 6
            k = 0
        end if
    end do

    if (k .gt. 0) then
        write (newin,5094) ux80(1:80)
    end if

    write (newin,1052)

    ! Iopg activity coefficient option switches.
    ! Note: iopg(1) = iopg1, etc.
    write (newin,7510)
7510 format('|Iopg Activity Coefficient Option Switches',' ("( 0)" marks default choices)',t80,'|')

    write (newin,1050)

    uxo = 'iopg('

    do n = 1,nopgmx
        j6 = index(uopgcx(n),'EQ6')

        if (j6 .gt. 0) then
            uxn = ' '
            write (uxn,'(i2)') n
            call lejust(uxn)
            j5 = ilnobl(uxn)
            uxo(6:8) = uxn(1:j5)
            j4 = ilnobl(uxo) + 1
            uxo(j4:j4) = ')'
            j3 = ilnobl(uopgtx(n))
            write (newin,7530) uxo(1:j4),uopgtx(n)(1:j3)
7530 format('|',a,' - ',a,':',t80,'|')

            i = iopg(n)

            do j = 1,jpgxpa
                ij = iopgox(j,n)
                uxs = ' '

                if (ij .eq. i) then
                    uxs = 'x'
                end if

                j2 = ilnobl(uopgox(j,n))

                if (j2 .gt. 0) then
                    write (newin,7550) uxs,ij,uopgox(j,n)(1:j2)
7550 format('|  [',a1,'] (',i2,') ',a,t80,'|')
                end if
            end do

            write (newin,1050)
        end if
    end do

    ! Index limits.
    write (newin,7610)
7610 format('|Matrix Index Limits',t80,'|')

    write (newin,1050)
    write (newin,7630) kct
7630 format('|No. of chem. elements   |',i5,'| (kct)',t80,'|')

    write (newin,7650) kbt
7650 format('|No. of basis species    |',i5,'| (kbt)',t80,'|')

    write (newin,7670) kmt
7670 format('|Index of last pure min. |',i5,'| (kmt)',t80,'|')

    write (newin,7690) kxt
7690 format('|Index of last sol-sol.  |',i5,'| (kxt)',t80,'|')

    write (newin,7710) kdim
7710 format('|Matrix size             |',i5,'| (kdim)',t80,'|')

    kprs = 0
    write (newin,7730) kprs
7730 format('|PRS data flag           |',i5,'| (kprs)',t80,'|')

    write (newin,1050)

    ! Species for which mass balances are defined.
    write (newin,7750)
7750 format('|Mass Balance Species (Matrix Row Variables)',5x,'|Units/Constraint| --',9x,'|',/'| (ubmtbi(n))',36x,'|(ujf6(jflgi(n)))| --',9x,'|')

    write (newin,1050)

    qstop = .false.

    do nbi = 1,nbti
        jfl = jflgi(nbi)
        jfli = 1

        if (jfl .eq. 30) then
            jfli = 2
        end if

        write (newin,7790) ubmtbi(nbi),ujf6(jfli)
7790 format('|',a48,'|',a16,'| --',9x,'|')
    end do

    if (qstop) then
        stop
    end if

    if (nbti .le. 0) then
        ux48 = 'None'
        ux16 = 'None'
        xx = 0.
        write (newin,7790) ux48,xx,ux16
    end if

    write (newin,1050)

    write (newin,7810)
7810 format('* Valid jflag strings (ujf6(jflgi(n))) are:',t80,'*')

    ux80(1:1) = '*'
    ux80(2:79) = ' '
    ux80(80:80) = '*'
    j3 = 6
    k = 0

    do jfli = 1,2
        jj = j3 + 19
        ux80(j3:jj) = ujf6(jfli)
        j3 = jj + 1
        k = k + 1

        if (k .ge. 3) then
            write (newin,7830) ux80(1:80)
7830 format(a)

            ux80(2:79) = ' '
            j3 = 6
            k = 0
        end if
    end do

    if (k .gt. 0) then
        write (newin,7830) ux80(1:80)
    end if

    write (newin,1052)

    ! Mass balance totals.
    write (newin,7850)
7850 format('|Mass Balance Totals (moles)',51x,'|')

    write (newin,1050)

    write (newin,7870)
7870 format('|Basis species (info. only)',6x,'|Equilibrium System',4x,'|Aqueous Solution',6x,'|',/'| (ubmtbi(n))',20x,'| (mtbi(n))',12x,'| (mtbaqi(n))',10x,'|')

    write (newin,1050)

    do nbi = 1,nbti
        write (newin,7890) ubmtbi(nbi),mtbi(nbi),mtbaqi(nbi)
7890 format('|',a32,'|',1pe22.15,'|',e22.15,'|')
    end do

    ux48 = 'Electrical imbalance'
    write (newin,7890) ux48,electr,electr
    write (newin,1050)

    ! Ordinary basis switches.
    write (newin,7970)
7970 format('|Ordinary Basis Switches (for numerical purposes only)',4x,'| (nobswt)',11x,'|')

    write (newin,1050)

    do n = 1,nobswt
        write (newin,7990) uobsw(1,n),uobsw(2,n)
7990 format('|Replace |',a48,'| (uobsw(1,n))',7x,'|',/'|   with |',a48,'| (uobsw(2,n))',7x,'|')

        write (newin,1050)
    end do

    if (nobswt .le. 0) then
        ux48a = 'None'
        ux48b = 'None'
        write (newin,7990) ux48a,ux48b
        write (newin,1050)
    end if

    ! Matrix column variables and values.
    write (newin,7910)
7910 format('|Matrix Column Variables and Values',44x,'|')

    write (newin,1050)

    write (newin,7930)
7930 format('|Basis species (uzveci(n))',23x,'|Log moles (zvclgi(n)) | --   |')

    write (newin,1050)

    do krow = 1,kdim
        write (newin,7950) uzveci(krow),zvclgi(krow)
7950 format('|',a48,'|',1pe22.15,'| --   |')
    end do

    write (newin,1050)

    ! Phases and species in the PRS (there are none in the context of
    ! EQ3NR models).
    write (newin,8110)
8110 format('|Phases and Species in the PRS',49x,'|')

    write (newin,1050)

    ux24 = 'None'
    xx = 0.
    write (newin,8130) ux24
8130 format('|Phase',11x,'|',a24,'| (uprphi(n))',t80,'|')

    write (newin,1050)
    write (newin,8150) xx
8150 format('|->|No. of Moles',4x,'|',1pe22.15,'| (mprphi(n))',t80,'|')

    write (newin,1050)
    write (newin,8170)
8170 format('|--->|Species',17x,'|No. of Moles',10x,'| --',t80,'|',/'|--->| (uprspi(i,n))',10x,'| (mprspi(i,n))',8x,'| --',t80,'|')

    write (newin,1050)
    write (newin,8190) ux24,xx
8190 format('|--->|',a24,'|',1pe22.15,'| --',t80,'|')

    write (newin,1050)

    ! End of problem.
    write (newin,8310)
8310 format('|End of problem',64x,'|')

    write (newin,1050)

999 continue
end subroutine wr3pkd