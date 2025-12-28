subroutine wr3d7(cspb,fep,iktmax,iodb,iopg,iopr,iopt,itermx,jflagb,jxmod,kxmod,ncompb,newin,nodbmx,nopgmx,noprmx,noptmx,nsq,nsqmax,ntitl,ntitmx,nttyo,nxmdmx,nxmod,nxtb,nxtmax,rho,tempc,tdspkg,tdspl,tolbt,toldl,tolsat,ubasis,uebal,umemb,uphas1,uphas2,uredox,usolb,uspecb,utitl,xbarb,uxmd24,xlkmod)
    !! This subroutine writes the EQ3NR input file in menu-style ("D")
    !! format for version 7.0-7.1.
    !! This subroutine is called by:
    !!   XCON3/xcon3.f
    implicit none

    ! Calling sequence variable declarations.
    integer :: iktmax
    integer :: nodbmx
    integer :: nopgmx
    integer :: noprmx
    integer :: noptmx
    integer :: nsqmax
    integer :: ntitmx
    integer :: nxmdmx
    integer :: nxtmax

    integer :: iodb(nodbmx)
    integer :: iopg(nopgmx)
    integer :: iopr(noprmx)
    integer :: iopt(noptmx)
    integer :: jflagb(nsqmax)
    integer :: jxmod(nxmdmx)
    integer :: kxmod(nxmdmx)
    integer :: ncompb(nxtmax)

    integer :: itermx
    integer :: newin
    integer :: nsq
    integer :: ntitl
    integer :: nttyo
    integer :: nxmod
    integer :: nxtb

    character(len=80) :: utitl(ntitmx)
    character(len=24) :: ubasis(nsqmax)
    character(len=24) :: umemb(iktmax,nxtmax)
    character(len=24) :: uphas1(nsqmax)
    character(len=24) :: uphas2(nsqmax)
    character(len=24) :: usolb(nxtmax)
    character(len=24) :: uspecb(nsqmax)
    character(len=24) :: uxmd24(nxmdmx)
    character(len=24) :: uebal
    character(len=24) :: uredox

    real(kind=8) :: cspb(nsqmax)
    real(kind=8) :: xbarb(iktmax,nxtmax)
    real(kind=8) :: xlkmod(nxmdmx)

    real(kind=8) :: fep
    real(kind=8) :: rho
    real(kind=8) :: tempc
    real(kind=8) :: tdspkg
    real(kind=8) :: tdspl
    real(kind=8) :: tolbt
    real(kind=8) :: toldl
    real(kind=8) :: tolsat

    include 'xcon3/x3op7.h'

    ! Local variable declarations.
    integer :: i
    integer :: idesc
    integer :: idescx
    integer :: ik
    integer :: iktb
    integer :: j
    integer :: jfl
    integer :: jj
    integer :: j2
    integer :: j3
    integer :: n
    integer :: ncount
    integer :: nredox
    integer :: ns
    integer :: nx
    integer :: ilnobl

    logical :: qlcase

    character(len=24) :: unone
    character(len=24) :: ux24
    character(len=16) :: udef(nto3pa)
    character(len=16) :: ujx(0:3)
    character(len=16) :: ux16
    character(len=8) :: usup(-1:2)
    character(len=1) :: ux1(3)

    real(kind=8) :: xx

    data unone  /'none                    '/

    data udef(1)  / '1.e-10          '/
    data udef(2)  / '1.e-10          '/
    data udef(3)  / '0.5             '/
    data udef(4)  / '30              '/

    ! data udef(1)  / '1.0e-10  tolbt  '/
    ! data udef(2)  / '1.0e-10  toldl  '/
    ! data udef(3)  / '0.5      tolsat '/
    ! data udef(4)  / '30       itermx '/
    data ujx(0)   / 'aqueous         '/
    data ujx(1)   / 'mineral         '/
    data ujx(2)   / 'gas             '/
    data ujx(3)   / 'solid solution  '/

    data usup(-1) / 'suppress'/
    data usup(0)  / 'replace '/
    data usup(1)  / 'augmentk'/
    data usup(2)  / 'augmentg'/

    ! Write the new input file.
    ! Title.
    write (newin,1090)
1090 format('|',70('-'),'|')

    do 110 n = 1,ntitl
        write (newin,1100) utitl(n)
1100 format('|',a70,'|')

110 continue

        write (newin,1090)

        ! Temperature and density.
        write (newin,1050) tempc,rho
1050 format('|Temperature (C)         |',f6.2,8x,'|Density(gm/cm3)|',f9.5,t72,'|')

        write (newin,1090)

        ! Total dissolved salts.
        ux1(1) = ' '
        ux1(2) = ' '
        ux1(3) = ' '

        if (tdspkg .gt. 0.) then
            ux1(1) = '*'
            xx = tdspkg
        else if (tdspl .gt. 0.) then
            ux1(2) = '*'
            xx = tdspl
        else
            ux1(3) = '*'
            xx = 0.
        end if

        if (xx .gt. 0.) then
            write (newin,1060) xx,(ux1(i), i = 1,3)
1060 format('|Total Dissolved Salts   | ',g12.5,t41,'|',a1,'mg/kg',t49,'|',a1,'mg/l',t57,'|',a1,'not used',t72,'|')
        else
            write (newin,1070) (ux1(i), i = 1,3)
1070 format('|Total Dissolved Salts   | ',t41,'|',a1,'mg/kg',t49,'|',a1,'mg/l',t57,'|',a1,'not used',t72,'|')
        end if

        write (newin,1090)

        ! Species for electrical balancing.
        ux1(1) = ' '
        ux1(2) = ' '

        if (uebal(1:8) .eq. 'pick1.  ') then
            ux1(1) = '*'
            uebal = ' '
        else if (uebal(1:8) .eq. '        ') then
            ux1(2) = '*'
        end if

        write (newin,1080) uebal(1:14),ux1(1),ux1(2)
1080 format('|Electrical Balancing on |',a14,'|',a1,'code selects  |',a1,'not performed|')

        write (newin,1090)

        ! Basis species and associated constraints.
        write (newin,1110)
1110 format('|SPECIES   | BASIS SWITCH/CONSTRAINT |',' CONCENTRATION| UNITS OR TYPE    |')

        write (newin,1090)

        if (nsq .le. 0) then
            write (newin,1120)
1120 format('| none',t12,'|',t38,'|',t51,'|',t72,'|')

            write (newin,1090)
        end if

        ux16 = ' '

        if (iopt(1) .eq. -2) then
            ux16 = 'pe'
        else if (iopt(1) .eq. -1) then
            ux16 = 'Eh'
        else if (iopt(1) .eq. 0) then
            ux16 = 'LogfO2'
        end if

        if (ux16(1:8) .ne. '        ') then
            write (newin,1130) fep,ux16
1130 format('|redox     |',24x,' | ',g10.4,3x,'|',a16,t72,'|')
        end if

        if (iopt(1) .eq. -3) then
            qlcase = .true.

            do 115 ns = 1, nsq
                ux24 = uspecb(ns)
                call locase(ux24)

                if (ux24(1:24) .ne. uspecb(ns)(1:24)) then
                    qlcase = .false.
                    go to 117
                end if

115 continue

117 continue
                if (qlcase) then
                    uredox = 'o2(g)'
                else
                    uredox = 'O2(g)'
                end if
            end if

            nredox = 0

            if (iopt(1).eq.1 .or. iopt(1).eq.-3) then
                write (newin,1140)
1140 format('|redox',t12,'| Defined by next species',t38,'|',t53,'|','Redox couple',t72,'|')

                do 120 ns = 1,nsq
                    if (uredox(1:24) .eq. uspecb(ns)(1:24)) then
                        nredox = ns
                        go to 130
                    end if

120 continue

                    j2 = ilnobl(uredox)
                    write (nttyo,1145) uredox(1:j2)
                    write (newin,1145) uredox(1:j2)
1145 format(" * Error - (XCON3/wr3d7) Can't find the redox",/7x,'couple-defining  species "',a,'" in the set of',/7x,'basis species.')

                    go to 999

130 continue
                    if (nredox .ne. 1) then
                        ux24 = uspecb(1)
                        uspecb(1) = uspecb(nredox)
                        uspecb(nredox) = ux24

                        xx = cspb(1)
                        cspb(1) = cspb(nredox)
                        cspb(nredox) = xx

                        jfl = jflagb(1)
                        jflagb(1) = jflagb(nredox)
                        jflagb(nredox) = jfl

                        ux24 = ubasis(1)
                        ubasis(1) = ubasis(nredox)
                        ubasis(nredox) = ux24

                        ux24 = uphas1(1)
                        uphas1(1) = uphas1(nredox)
                        uphas1(nredox) = ux24

                        ux24 = uphas2(1)
                        uphas2(1) = uphas2(nredox)
                        uphas2(nredox) = ux24

                        nredox = 1
                    end if
                end if

                do 180 ns = 1,nsq
                    if (uphas1(ns)(1:8) .ne. '        ') then
                        ux24 = uphas1(ns)
                    else
                        ux24 = ubasis(ns)
                    end if

                    jfl = jflagb(ns)
                    write (newin,1150) uspecb(ns),ux24,cspb(ns),ujflg7(jfl)
1150 format('|',a10,'|',a24,t38,'|',g11.5,t53,'|',a18,'|')

                    if (jfl .eq. 20) then
                        do 140 nx = 1,nxtb
                            if (uphas1(ns)(1:24) .eq. usolb(nx)(1:24)) then
                                go to 150
                            end if

140 continue

                            j2 = ilnobl(uphas1(ns))
                            write (nttyo,1160) uphas1(ns)(1:j2)
                            write (newin,1160) uphas1(ns)(1:j2)
1160 format(" * Error - (XCON3/wr3d7) Can't find the solid",/7x,'solution "',a,'" in the set of solid solutions',/7x,'for which compositions are defined.')

                            go to 999

150 continue

                            iktb = ncompb(nx)

                            do 160 ik = 1,iktb
                                if (uphas2(ns)(1:24) .eq. umemb(ik,nx)(1:24)) then
                                    go to 170
                                end if

160 continue

                                j3 = ilnobl(uphas2(ns))
                                write (nttyo,1170) uphas2(ns)(1:j3),uphas1(ns)(1:j2)
                                write (newin,1170) uphas2(ns)(1:j3),uphas1(ns)(1:j2)
1170 format(" * Error - (XCON3/wr3d7) Can't find the end-member",/7x,'"',a,'" in the set of components used to,'    /7x,'specify the composition of solid solution "',a,'".')

                                go to 999

170 continue
                                write (newin,1180) umemb(ik,nx),xbarb(ik,nx)
1180 format('|',10x,'|',a24,t38,'|',g10.4,t53,'|',t72,'|')

                                do 175 i = 1,iktb
                                    if (i .ne. ik) then
                                        write (newin,1180) umemb(i,nx),xbarb(i,nx)
                                    end if

175 continue
                                end if

180 continue

190 continue
                                write (newin,1090)

                                ! Input solid solutions, excluding those defined as part of
                                ! solubility equilibrium constraints.
                                write (newin,1800)
1800 format('|Input Solid Solutions',t72,'|')

                                write (newin,1090)

                                ncount = 0

                                if (nxtb .gt. 0) then
                                    do 210 nx = 1,nxtb
                                        do 200 ns = 1,nsq
                                            if (jflagb(ns) .eq. 20) then
                                                if (usolb(nx)(1:24) .ne. uphas1(ns)(1:24)) then
                                                    ncount = ncount + 1
                                                    write (newin,1810) usolb(nx)
1810 format('|',a10,'|',t38,'|',t53,'|',t72,'|')

                                                    iktb = ncompb(nx)

                                                    do 195 i = 1,iktb
                                                        write (newin,1180) umemb(i,nx),xbarb(i,nx)
195 continue
                                                    end if
                                                end if

200 continue

210 continue
                                            end if

                                            if (ncount .le. 0) then
                                                write (newin,1820)
1820 format('| none',5x,'|',t38,'|',t53,'|',t72,'|')
                                            end if

                                            write (newin,1090)

                                            ! Nxmod options.
                                            write (newin,3080)
3080 format('|SUPPRESSED SPECIES',3x,'(suppress,replace,augmentk,augmentg)    value',t72,'|')

                                            write (newin,1090)

                                            if (nxmod .gt. 0) then
                                                do 520 n = 1,nxmod
                                                    write (newin,3090) uxmd24(n),ujx(jxmod(n)),usup(kxmod(n)),xlkmod(n)
3090 format('| ',a24,t26,'| ',a14,' | ',a8,' | ',g12.5,t72,'|')

520 continue
                                                else
                                                    write(newin,3095) unone
3095 format('| ',a8,t26,'|',t43,'|',t54,'|',t72,'|')
                                                end if

                                                write (newin,1090)

                                                ! Options. Iopt, iopg, and iopr option switches, but not iodb
                                                ! option switches. Skip any of the above which happen to be
                                                ! development options.
                                                ! Note: iopt(1) = iopt1, etc.
                                                write (newin,1460)
1460 format('|OPTIONS',t72,'|')

                                                write (newin,1090)

                                                do 410 idescx = 1,nop3pa
                                                    if (uvar3(idescx)(1:4) .ne. 'iodb') then
                                                        j = ilnobl(uopt3(idescx))
                                                        write (newin,1465) uopt3(idescx)(1:j)
1465 format('| - ',a,' -',t72,'|')

                                                        do 400 idesc = 1,nod3pa
                                                            if (iopti3(idesc) .eq. idescx) then
                                                                if (uvar3(idescx)(1:4) .eq. 'iopt') then
                                                                    if (iopt(index3(idescx)) .eq. ivalu3(idesc)) then
                                                                        ux1(1) = '*'
                                                                    else
                                                                        ux1(1) = ' '
                                                                    end if
                                                                else if (uvar3(idescx)(1:4) .eq. 'iopg') then
                                                                    if (iopg(index3(idescx)) .eq. ivalu3(idesc)) then
                                                                        ux1(1) = '*'
                                                                    else
                                                                        ux1(1) = ' '
                                                                    end if
                                                                else if (uvar3(idescx)(1:4) .eq. 'iopr') then
                                                                    if (iopr(index3(idescx)) .eq. ivalu3(idesc)) then
                                                                        ux1(1) = '*'
                                                                    else
                                                                        ux1(1) = ' '
                                                                    end if
                                                                else
                                                                    ux1(1) = ' '
                                                                end if

                                                                j = ilnobl(udesc3(idesc))
                                                                write (newin,1470) ux1(1),udesc3(idesc)(1:j)
1470 format('|',t5,a1,t7,a,t72,'|')
                                                            end if

400 continue
                                                        end if

410 continue

                                                        write (newin,1090)

                                                        ! Iodb switches.
                                                        write (newin,1480)
1480 format('|DEBUGGING SWITCHES (o-off, 1,2-on, default is off)',t72,'|')

                                                        write (newin,1090)

                                                        do 420 idescx = 1,ndb3pa
                                                            i = idbugi(idescx)
                                                            j2 = ilnobl(udebug(idescx))
                                                            j = index(udebug(idescx),'generic debugging')

                                                            if (j .eq. 0) then
                                                                j = index(udebug(idescx),'pre-Newton-Raphson')
                                                            end if

                                                            if (j .eq. 0) then
                                                                j = index(udebug(idescx),'stoichiometric factors')

                                                                if (j .gt. 0) then
                                                                    jj = index(udebug(idescx),'factors calculation')

                                                                    if (jj .gt. 0) then
                                                                        j = 0
                                                                    end if
                                                                end if
                                                            end if

                                                            if (j .eq. 0) then
                                                                write (newin,1490) iodb(i),udebug(idescx)(1:j2)
1490 format('|',i1,2x,a,t72,'|')
                                                            else
                                                                write (newin,1495) iodb(i),udebug(idescx)(1:j2)
1495 format('|',i1,2x,a,t72,'|2')
                                                            end if

420 continue

                                                            write (newin,1090)

                                                            ! Development options (there are none).
                                                            write (newin,1700)
1700 format('|DEVELOPMENT OPTIONS  (used for code development)',t72,'|')

                                                            write (newin,1090)
                                                            write (newin,1710)
1710 format('| none',t72,'|')

                                                            write (newin,1090)

                                                            ! Tolerances.
                                                            write (newin,1500)
1500 format('|TOLERANCES',19x,'(desired values)',8x,'(defaults)',7x,'|')

                                                            write (newin,1090)

                                                            if (tolbt .eq. 0) then
                                                                write (newin,1510) utol3(1),udef(1)
1510 format('|',a27,t28,'|',t51,'|',a16,t72,'|')
                                                            else
                                                                write (newin,1520) utol3(1),tolbt,udef(1)
1520 format('|',a27,t28,'|',i10,t51,'|',a16,t72,'|')
                                                            end if

                                                            if (toldl .eq. 0) then
                                                                write (newin,1510) utol3(2),udef(2)
                                                            else
                                                                write (newin,1520) utol3(2),toldl,udef(2)
                                                            end if

                                                            if (tolsat .eq. 0) then
                                                                write (newin,1510) utol3(3),udef(3)
                                                            else
                                                                write (newin,1520) utol3(3),tolsat,udef(3)
                                                            end if

                                                            if (itermx .eq. 0) then
                                                                write (newin,1510) utol3(4),udef(4)
                                                            else
                                                                write (newin,1520) utol3(4),itermx,udef(4)
                                                            end if

                                                            write (newin,1090)

999 continue
                                                        end subroutine wr3d7