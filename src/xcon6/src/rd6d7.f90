subroutine rd6d7(cdac,cesrb,csigma,dlzmx1,dlzmx2,dlzidp,dzprlg,dzprnt,eact,electr,fk,iact,iktbt,iktmax,imchmx,imech,iodb,iopg,iopr,iopt,ioscan,itermx,jcode,jreac,jtemp,jxmod,kct,kdim,kmax,kmt,kprs,ksq,ksppmx,kstpmx,kxmod,kxt,hact,modr,moffg,morr,mprs,mteaqb,mteb,nctmax,ndact,ndctmx,nesrbt,nffg,nffgmx,ninpts,nmodl1,nmodl2,nodbmx,nopgmx,noprmx,noptmx,nordlm,npslmx,nprmn,nprmx,nprsmx,nrct,nrctmx,nrk,nsk,nsrtmx,nsscmx,nsslmx,ntitl1,ntitl2,ntitmx,ntrymx,nttkmx,nttyo,nxmdmx,nxmod,nxopex,nxopmx,nxopt,nxpemx,nxrtmx,qend,qrderr,rk0,rxbarb,sscrew,sk,tempci,tempc0,timemx,tolbt,toldl,tolsat,tolsst,tolx,trk0,tstrt,ttk,udac,uelemb,uendb,uesrb,uffg,undms,unrms,uprs,ureac,utitl1,utitl2,uxct16,uxmd24,uxopex,uxopt,vreac,xlkffg,xlkmod,zimax,zistrt,zkfac,zklogl,zklogu,zvclgi)
    !! This subroutine reads the EQ6 input file in menu-style ("D")
    !! format for versions 7.0-7.2. It thus encompasses two version
    !! levels, '7.0' and '7.2'. The input file in this format is not
    !! identical for these two version levels. However, the line
    !! parsing capability can handle the differences, which consist
    !! only of differences in field sizes.
    !! This subroutine is called by:
    !!   XCON6/xcon6.f
    implicit none

    ! Calling sequence variable declarations.
    integer :: iktmax
    integer :: imchmx
    integer :: kmax
    integer :: nctmax
    integer :: ndctmx
    integer :: nffgmx
    integer :: nodbmx
    integer :: nopgmx
    integer :: noprmx
    integer :: noptmx
    integer :: nprsmx
    integer :: nrctmx
    integer :: nsrtmx
    integer :: nsscmx
    integer :: ntitmx
    integer :: nttkmx
    integer :: nxmdmx
    integer :: nxopmx
    integer :: nxpemx
    integer :: nxrtmx

    integer :: iact(imchmx,2,nrctmx)
    integer :: iktbt(nxrtmx)
    integer :: imech(2,nrctmx)
    integer :: iodb(nodbmx)
    integer :: iopg(nopgmx)
    integer :: iopr(noprmx)
    integer :: iopt(noptmx)
    integer :: jcode(nrctmx)
    integer :: jreac(nrctmx)
    integer :: jxmod(nxmdmx)
    integer :: kxmod(nxmdmx)
    integer :: ndact(imchmx,2,nrctmx)
    integer :: nesrbt(nsrtmx)
    integer :: nrk(2,nrctmx)
    integer :: nsk(nrctmx)

    integer :: ioscan
    integer :: itermx
    integer :: jtemp
    integer :: kct
    integer :: kdim
    integer :: ksq
    integer :: kmt
    integer :: kprs
    integer :: ksppmx
    integer :: kstpmx
    integer :: kxt
    integer :: nffg
    integer :: ninpts
    integer :: nmodl1
    integer :: nmodl2
    integer :: nordlm
    integer :: npslmx
    integer :: nprmn
    integer :: nprmx
    integer :: nrct
    integer :: nsslmx
    integer :: ntitl1
    integer :: ntitl2
    integer :: ntrymx
    integer :: nttyo
    integer :: nxmod
    integer :: nxopex
    integer :: nxopt

    logical :: qend
    logical :: qrderr

    character(len=80) :: utitl1(ntitmx)
    character(len=80) :: utitl2(ntitmx)
    character(len=48) :: uprs(nprsmx)
    character(len=24) :: udac(ndctmx,imchmx,2,nrctmx)
    character(len=24) :: uendb(iktmax,nxrtmx)
    character(len=24) :: uffg(nffgmx)
    character(len=24) :: undms(kmax)
    character(len=24) :: unrms(kmax)
    character(len=24) :: ureac(nrctmx)
    character(len=24) :: uxmd24(nxmdmx)
    character(len=24) :: uxopex(nxpemx)
    character(len=16) :: uxct16(nxopmx)
    character(len=8) :: uelemb(nctmax)
    character(len=8) :: uesrb(nctmax,nsrtmx)
    character(len=8) :: uxopt(nxopmx)

    real(kind=8) :: cdac(ndctmx,imchmx,2,nrctmx)
    real(kind=8) :: cesrb(nctmax,nsrtmx)
    real(kind=8) :: csigma(imchmx,2,nrctmx)
    real(kind=8) :: eact(imchmx,2,nrctmx)
    real(kind=8) :: fk(nrctmx)
    real(kind=8) :: hact(imchmx,2,nrctmx)
    real(kind=8) :: modr(nrctmx)
    real(kind=8) :: moffg(nffgmx)
    real(kind=8) :: morr(nrctmx)
    real(kind=8) :: mprs(nprsmx)
    real(kind=8) :: mteaqb(nctmax)
    real(kind=8) :: mteb(nctmax)
    real(kind=8) :: rk0(imchmx,2,nrctmx)
    real(kind=8) :: rxbarb(iktmax,nxrtmx)
    real(kind=8) :: sscrew(nsscmx)
    real(kind=8) :: sk(nrctmx)
    real(kind=8) :: ttk(nttkmx)
    real(kind=8) :: trk0(imchmx,2,nrctmx)
    real(kind=8) :: vreac(nrctmx)
    real(kind=8) :: xlkffg(nffgmx)
    real(kind=8) :: xlkmod(nxmdmx)
    real(kind=8) :: zvclgi(kmax)

    real(kind=8) :: dlzmx1
    real(kind=8) :: dlzmx2
    real(kind=8) :: dlzidp
    real(kind=8) :: dzprlg
    real(kind=8) :: dzprnt
    real(kind=8) :: electr
    real(kind=8) :: tempci
    real(kind=8) :: tempc0
    real(kind=8) :: timemx
    real(kind=8) :: tolbt
    real(kind=8) :: toldl
    real(kind=8) :: tolsat
    real(kind=8) :: tolsst
    real(kind=8) :: tolx
    real(kind=8) :: tstrt
    real(kind=8) :: zkfac
    real(kind=8) :: zklogl
    real(kind=8) :: zklogu
    real(kind=8) :: zimax
    real(kind=8) :: zistrt

    include 'xcon6/x6op7.h'

    ! Local parameter declarations.
    !   nfldpa = maximum number of fields per line
    !   nlchpa = character length of a line
    integer :: nfldpa
    integer :: nlchpa

    parameter (nfldpa = 8,nlchpa = 80)

    ! Local variable declarations.
    integer :: i
    integer :: im
    integer :: idesc
    integer :: idescx
    integer :: iktb
    integer :: ivar
    integer :: j
    integer :: jdesc
    integer :: j1
    integer :: j2
    integer :: ksb
    integer :: n
    integer :: ncb
    integer :: nfldmx
    integer :: nfldt
    integer :: nfldtx
    integer :: nlchmx
    integer :: nmark
    integer :: nmarks
    integer :: nrc
    integer :: nsrt
    integer :: nxrt

    integer :: ilnobl

    character(len=nlchpa) :: ufield(nfldpa)
    character(len=nlchpa) :: uline1
    character(len=nlchpa) :: uline2
    character(len=nlchpa) :: ulscr
    character(len=nlchpa) :: uheadr
    character(len=nlchpa) :: uheadx

    character(len=80) :: ustr
    character(len=24) :: ux
    character(len=1) :: ux1

    real(kind=8) :: var

    nfldmx = nfldpa
    nlchmx = nlchpa

    qrderr = .false.

    ! Title.
    qend = .false.
    read (ninpts,1000,end=100,err=990) uline1
1000 format(a80)

    call parsln(nfldmx,nfldt,nlchmx,ufield,uline1,ulscr)
    ustr = ufield(1)

    if (ustr(1:8) .ne. '--------') then
        write (nttyo,1010) uline1
1010 format(/' * Error - (XCON6/rd6d7) The first line of a "D"',/7x,'format input file must begin with "|--------".',/7x,'The first line read from the old input file begins with-',/7x,'"',a70,'".')

        go to 990
    end if

    go to 105

100 continue
    qend = .true.
    go to 999

105 continue
    do 110 n = 1,ntitmx + 1
        read (ninpts,1000,err=990) uline1
        call parsln(nfldmx,nfldt,nlchmx,ufield,uline1,ulscr)
        ustr = ufield(1)

        if (ustr(1:8) .eq. '--------') then
            go to 120
        end if

        utitl1(n) = ufield(1)
110 continue

        write (nttyo,1015) ntitmx
1015 format(/' * Error - (XCON6/rd6d7) Have too many lines in the',/7x,'main title. The code is only dimensioned for ',i4,/7x,'lines. Reduce the size of the title or increase the',/7x,'dimensioning parameter ntitpa.')

        go to 990

120 continue
        ntitl1 = n - 1

        ! Nmodl option switches.
        nmodl2 = 0
        nmodl1 = 0

        uheadx = 'calculational mode'
        nfldtx = 4
        call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

        if (qrderr) then
            go to 999
        end if

        call gmarko(nfldmx,nfldt,nmark,nttyo,ufield,uline1)

        if (nmark .gt. 0) then
            nmodl2 = nmark - 2
        end if

        uheadx = 'model type'
        nfldtx = 4
        call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

        if (qrderr) then
            go to 999
        end if

        call gmarko(nfldmx,nfldt,nmark,nttyo,ufield,uline1)

        if (nmark .gt. 0) then
            nmodl1 = nmark - 1
        end if

        ! Temperature parameters.
        ! Note: ttk(1) = tk1, etc.
        jtemp = 0

        uheadx = 'temperature model'
        nfldtx = 3
        call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

        if (qrderr) then
            go to 999
        end if

        call gmarko(nfldmx,nfldt,nmark,nttyo,ufield,uline1)

        if (nmark .gt. 0) then
            jtemp = nmark - 2
        end if

        uheadx = 'tstart(c)'
        nfldtx = 8
        call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

        if (qrderr) then
            go to 999
        end if

        ustr = ufield(2)
        call chreal(nttyo,qrderr,ustr,var)

        if (qrderr) then
            go to 990
        end if

        tempc0 = var

        do 130 i = 1,3
            ustr = ufield(2*i + 2)
            call chreal(nttyo,qrderr,ustr,var)

            if (qrderr) then
                go to 990
            end if

            ttk(i) = var
130 continue

            ! Zi and time parameters.
            ! Note: cplim does not appear on the "D" format input file.
            uheadx = 'starting value of zi'
            nfldtx = 4
            call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

            if (qrderr) then
                go to 999
            end if

            ustr = ufield(2)
            call chreal(nttyo,qrderr,ustr,var)

            if (qrderr) then
                go to 990
            end if

            zistrt = var
            ustr = ufield(4)
            call chreal(nttyo,qrderr,ustr,var)

            if (qrderr) then
                go to 990
            end if

            zimax = var

            uheadx = 'starting time (sec)'
            nfldtx = 4
            call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

            if (qrderr) then
                go to 999
            end if

            ustr = ufield(2)
            call chreal(nttyo,qrderr,ustr,var)

            if (qrderr) then
                go to 990
            end if

            tstrt = var
            ustr = ufield(4)
            call chreal(nttyo,qrderr,ustr,var)

            if (qrderr) then
                go to 990
            end if

            timemx = var

            uheadx = 'max. steps'
            nfldtx = 4
            call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

            if (qrderr) then
                go to 999
            end if

            ustr = ufield(2)
            call chrint(ivar,nttyo,qrderr,ustr)

            if (qrderr) then
                go to 990
            end if

            kstpmx = ivar
            ustr = ufield(4)
            call chrint(ivar,nttyo,qrderr,ustr)

            if (qrderr) then
                go to 990
            end if

            ksppmx = ivar

            ! Print interval parameters.
            uheadx = 'linear print interval'
            nfldtx = 4
            call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

            if (qrderr) then
                go to 999
            end if

            ustr = ufield(2)
            call chreal(nttyo,qrderr,ustr,var)

            if (qrderr) then
                go to 990
            end if

            dzprnt = var
            ustr = ufield(4)
            call chreal(nttyo,qrderr,ustr,var)

            if (qrderr) then
                go to 990
            end if

            dzprlg = var

            ! Plot interval parameters.
            ! Note: dzplot, dzpllg, and ksplmx currently do not appear on the
            ! "D" format input file.
            ! Ifile.
            ! Note: ifile does not appear on the "D" format input file.
            ! Nxopt mineral subset selection suppression options.
            nxopt = 0

            uheadx = 'suppress mineral phases'
            nfldtx = 1
            call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

            if (qrderr) then
                go to 999
            end if

            uheadx = 'phases w/ elements'
            nfldtx = 3
            call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)

            if (qrderr) then
                go to 999
            end if

140 continue
            do 150 n = 2,3
                ustr = ufield(n)
                call locase(ustr)

                if (ustr(1:8) .ne. '        ') then
                    nxopt = nxopt + 1

                    if (nxopt .gt. nxopmx) then
                        write (nttyo,1025) nxopmx
1025 format(/' * Error - (XCON6/rd6d7) Have too many mineral',/7x,'subset-selection suppression options. The code is',/7x,'only dimensioned for ',i3,' such options. Reduce the',/7x,'number of options or increase the dimensioning',/7x,'parameter nxoppa.')

                        go to 990
                    end if

                    if (ustr(1:8) .eq. 'all     ') then
                        uxopt(nxopt) = 'all'
                        uxct16(nxopt) = '        '
                    else
                        uxopt(nxopt) = 'alwith'
                        uxct16(nxopt) = ustr
                    end if
                end if

150 continue

                nfldtx = 0
                call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

                if (qrderr) then
                    go to 999
                end if

                ustr = ufield(1)

                if (ustr(1:8) .eq. '--------') then
                    go to 175
                end if

                uheadr = ufield(1)
                call locase(uheadr)
                j2 = ilnobl(uheadx)
                j = index(uheadr,uheadx(1:j2))

                if (j .gt. 0) then
                    nfldtx = 3

                    if (nfldt .ne. nfldtx) then
                        write (nttyo,1030) nfldt,nfldtx,uline1
1030 format(/' * Warning - (XCON6/rd6d7) Found ',i2,' fields',/7x,'where ',i2,' were expected on the line beginning with-',/7x,'"',a70,'".')
                    end if

                    go to 140
                end if

                nxopex = 0

                uheadx = 'phases except'
                j2 = ilnobl(uheadx)
                j = index(uheadr,uheadx(1:j2))

                if (j .eq. 0) then
                    write (nttyo,1020) uheadx,uline1
1020 format(/' * Error - (XCON6/rd6d7) Was expecting to find the',/7x,'header beginning with-',/7x,'"',a70,'"',/7x,'on the line beginning with-',/7x,'"',a70,'".')

                    go to 990
                end if

160 continue
                do 170 n = 2,3
                    ustr = ufield(n)

                    if (ustr(1:8) .ne. '        ') then
                        nxopex = nxopex + 1

                        if (nxopex .gt. nxpemx) then
                            write (nttyo,1027) nxpemx
1027 format(/' * Error - (XCON6/rd6d7) Have too many',/7x,'exceptions specified to the mineral subset-selection',/7x,'suppression options. The code is only dimensioned',/7x,'for ',i3,'exceptions. Reduce the number of exceptions',/7x,'or increase the dimensioning parameter nxoppa.')

                            go to 990
                        end if

                        uxopex(nxopex) = ustr
                    end if

170 continue

                    nfldtx = 0
                    call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

                    if (qrderr) then
                        go to 999
                    end if

                    uheadr = ufield(1)
                    call locase(uheadr)
                    j2 = ilnobl(uheadx)
                    j = index(uheadr,uheadx(1:j2))

                    if (j .gt. 0) then
                        nfldtx = 3

                        if (nfldt .ne. nfldtx) then
                            write (nttyo,1030) nfldt,nfldtx,uline1
                        end if

                        go to 160
                    end if

                    ustr = ufield(1)

                    if (ustr(1:8) .ne. '--------') then
                        write (nttyo,1040) uline1
1040 format(/' * Error - (XCON6/rd6d7) Found the line beginning',' with-',/7x,'"',a70,'"',/7x,'where a dashed separator line was expected.')

                        go to 990
                    end if

175 continue

                    ! Nffg options.
                    nffg = 0

                    uheadx = 'fixed fugacity phases- species, '
                    nfldtx = 1
                    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

                    if (qrderr) then
                        go to 999
                    end if

                    uheadr = ufield(1)
                    call locase(uheadr)
                    j2 = ilnobl(uheadx)
                    j = index(uheadr,uheadx(1:j2))

                    if (j .eq. 0) then
                        write (nttyo,1020) uheadx,uline1
                        go to 990
                    end if

                    nfldtx = 0
                    call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

                    if (qrderr) then
                        go to 999
                    end if

                    ustr = ufield(1)
                    call locase(ustr)

                    if (ustr(1:4) .eq. 'none') then
                        go to 190
                    end if

                    nfldtx = 3

                    if (nfldt .ne. nfldtx) then
                        write (nttyo,1030) nfldt,nfldtx,uline1
                    end if

180 continue
                    nffg = nffg + 1

                    if (nffg .gt. nffgmx) then
                        write (nttyo,1035) nffgmx
1035 format(/' * Error - (XCON6/rd6d7) Have too many gases whose',/7x,'fugacities are to be fixed. The code is only dimensioned',/7x,'for ',i4,' such gases. Reduce the number of gases or',/7x,'increase the dimensioning parameter nffgpa.')

                        go to 990
                    end if

                    uffg(nffg) = ufield(1)
                    ustr = ufield(2)
                    call chreal(nttyo,qrderr,ustr,var)

                    if (qrderr) then
                        go to 990
                    end if

                    moffg(nffg) = var
                    ustr = ufield(3)
                    call chreal(nttyo,qrderr,ustr,var)

                    if (qrderr) then
                        go to 990
                    end if

                    xlkffg(nffg) = var

190 continue
                    nfldtx = 0
                    call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

                    if (qrderr) then
                        go to 999
                    end if

                    ustr = ufield(1)

                    if (ustr(1:8) .eq. '--------') then
                        go to 200
                    end if

                    nfldtx = 3

                    if (nfldt .ne. nfldtx) then
                        write (nttyo,1030) nfldt,nfldtx,uline1
                    end if

                    go to 180

                    ! Read a second dashed separator line.
200 continue
                    nfldtx = 0
                    call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

                    if (qrderr) then
                        go to 999
                    end if

                    ustr = ufield(1)

                    if (ustr(1:8) .ne. '--------') then
                        write (nttyo,1040) uline1
                        go to 990
                    end if

                    ! Reactants.
                    nrct = 0
                    nsrt = 0
                    nxrt = 0

                    uheadx = 'reactants'
                    nfldtx = 1
                    call rdd2l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,uline2,ulscr)

                    if (qrderr) then
                        go to 999
                    end if

                    uheadr = ufield(1)
                    call locase(uheadr)
                    j2 = ilnobl(uheadx)
                    j = index(uheadr,uheadx(1:j2))

                    if (j .eq. 0) then
                        write (nttyo,1020) uheadx,uline1
                        go to 990
                    end if

210 continue
                    nfldtx = 0
                    call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

                    if (qrderr) then
                        go to 999
                    end if

                    ustr = ufield(1)

                    if (ustr(1:8).eq.'--------') then
                        go to 330
                    end if

                    uheadr = ufield(1)
                    call locase(uheadr)

                    uheadx = 'REACTANT '
                    call locase(uheadx)
                    j2 = ilnobl(uheadx)
                    j = index(uheadr,uheadx(1:j2))

                    if (j .eq. 0) then
                        write (nttyo,1020) uheadx,uline1
                        go to 990
                    end if

                    nfldtx = 4

                    if (nfldt .ne. nfldtx) then
                        write (nttyo,1030) nfldt,nfldtx,uline1
                    end if

                    ustr = ufield(2)
                    call locase(ustr)

                    if (ustr(1:4) .eq. 'none') then
                        nfldtx = 1
                        call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

                        if (qrderr) then
                            go to 999
                        end if

                        ustr = ufield(1)

                        if (ustr(1:8) .ne. '--------') then
                            write (nttyo,1040) uline1
                            go to 990
                        end if

                        go to 330
                    end if

215 continue
                    nrc = nrc + 1

                    if (nrc .gt. nrctmx) then
                        write (nttyo,1037) nrctmx
1037 format(/' * Error - (XCON6/rd6d7) Have too many reactants',/7x,'The code is only dimensioned for ',i4,' reactants.',/7x,'Reduce the number of reactants or increase the',/7x,'dimensioning parameter nrctpa.')

                        go to 990
                    end if

                    nrct = nrc
                    ureac(nrc) = ufield(2)
                    ustr = ufield(4)
                    call chrint(ivar,nttyo,qrderr,ustr)

                    if (qrderr) then
                        go to 990
                    end if

                    jreac(nrc) = ivar

                    uheadx = 'moles remaining'
                    nfldtx = 4
                    call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)

                    if (qrderr) then
                        go to 999
                    end if

                    ustr = ufield(2)
                    call chreal(nttyo,qrderr,ustr,var)

                    if (qrderr) then
                        go to 990
                    end if

                    morr(nrc) = var
                    ustr = ufield(4)
                    call chreal(nttyo,qrderr,ustr,var)

                    if (qrderr) then
                        go to 990
                    end if

                    modr(nrc) = var

                    uheadx = 'reactant type'
                    nfldtx = 4
                    call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)

                    if (qrderr) then
                        go to 999
                    end if

                    ustr = ufield(2)

                    if (ustr(1:7) .eq. 'mineral') then
                        jcode(nrc) = 0
                    else if (ustr(1:14) .eq. 'solid solution') then
                        jcode(nrc) = 1
                    else if (ustr(1:7) .eq. 'special') then
                        jcode(nrc) = 2
                    else if (ustr(1:7) .eq. 'aqueous') then
                        jcode(nrc) = 3
                    else if (ustr(1:3) .eq. 'gas') then
                        jcode(nrc) = 4
                    else
                        jcode(nrc) = 0
                    end if

                    ustr = ufield(4)
                    call chreal(nttyo,qrderr,ustr,var)

                    if (qrderr) then
                        go to 990
                    end if

                    sk(nrc) = var

                    uheadx = 'surface type'
                    nfldtx = 4
                    call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)

                    if (qrderr) then
                        go to 999
                    end if

                    ustr = ufield(2)
                    call chrint(ivar,nttyo,qrderr,ustr)

                    if (qrderr) then
                        go to 990
                    end if

                    nsk(nrc) = ivar
                    ustr = ufield(4)
                    call chreal(nttyo,qrderr,ustr,var)

                    if (qrderr) then
                        go to 990
                    end if

                    fk(nrc) = var

                    if (jcode(nrc) .eq. 1) then
                        nxrt = nxrt + 1

                        if (nxrt .gt. nxrtmx) then
                            write (nttyo,1038) nxrtmx
1038 format(/' * Error - (XCON6/rd6d7) Have too many solid',/7x,'solution reactants. The code is only dimensioned',/7x,'for ',i4,' such reactants. Reduce the number of such',/7x,'reactants or increase the dimensioning parameter',' nxrtpa.')

                            go to 990
                        end if

                        uheadx = 'end-member'
                        nfldtx = 4
                        call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)

                        if (qrderr) then
                            go to 999
                        end if

                        iktb = 0
220 continue
                        iktb = iktb + 1

                        if (iktb .gt. iktmax) then
                            write (nttyo,1039) ureac(nrc),iktmax
1039 format(/' * Error - (XCON6/rd6d7) Have too many end-members',/7x,'in the solid solution reactant "',a24,'".',/7x,'The code is only dimensioned for ',i4,' end-members per',/7x,'solid solution. Reduce the number of end-members or',/7x,'increase the dimensioning parameter iktpar.')

                            go to 990
                        end if

                        uendb(iktb,nxrt) = ufield(2)
                        ustr = ufield(4)
                        call chreal(nttyo,qrderr,ustr,var)

                        if (qrderr) then
                            go to 990
                        end if

                        rxbarb(iktb,nxrt) = var

                        nfldtx = 4
                        call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

                        if (qrderr) then
                            go to 999
                        end if

                        uheadr = ufield(1)
                        call locase(uheadr)
                        j2 = ilnobl(uheadx)
                        j = index(uheadr,uheadx(1:j2))

                        if (j .gt. 0) then
                            go to 220
                        end if

                        iktbt(nxrt) = iktb
                    else
                        uheadx = 'end-member'
                        nfldtx = 4
230 continue
                        call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

                        if (qrderr) then
                            go to 999
                        end if

                        uheadr = ufield(1)
                        call locase(uheadr)
                        j2 = ilnobl(uheadx)
                        j = index(uheadr,uheadx(1:j2))

                        if (j .gt. 0) then
                            go to 230
                        end if
                    end if

                    if (jcode(nrc) .eq. 2) then
                        nsrt = nsrt + 1

                        if (nsrt .gt. nsrtmx) then
                            write (nttyo,1042) nsrtmx
1042 format(/' * Error - (XCON6/rd6d7) Have too many special',/7x,'reactants. The code is only dimensioned for ',i4,/7x,'such reactants. Reduce the number of such reactants',/7x,'or increase the dimensioning parameter nsrtpa.')

                            go to 990
                        end if

                        nrct = nrc
                        uheadx = 'volume'
                        j2 = ilnobl(uheadx)
                        j = index(uheadr,uheadx(1:j2))

                        if (j .eq. 0) then
                            write (nttyo,1020) uheadx,uline1
                            go to 990
                        end if

                        ustr = ufield(2)
                        call chreal(nttyo,qrderr,ustr,var)

                        if (qrderr) then
                            go to 990
                        end if

                        vreac(nsrt) = var

                        uheadx = 'element'
                        nfldtx = 4
                        call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)

                        if (qrderr) then
                            go to 999
                        end if

                        ncb = 0
240 continue
                        ncb = ncb + 1

                        if (ncb .gt. nctmax) then
                            write (nttyo,1044) ureac(nrc),nctmax
1044 format(/' * Error - (XCON6/rd6d7) Have too many chemical',/7x,'elements in the special  reactant "',a24,'".',/7x,'The code is only dimensioned for ',i4,' elements.',/7x,'Reduce the number of elements or increase the',/7x,'dimensioning parameter nctpar.')

                            go to 990
                        end if

                        uesrb(ncb,nsrt) = ufield(2)
                        ustr = ufield(4)
                        call chreal(nttyo,qrderr,ustr,var)

                        if (qrderr) then
                            go to 990
                        end if

                        cesrb(ncb,nsrt) = var

                        nfldtx = 4
                        call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

                        if (qrderr) then
                            go to 999
                        end if

                        uheadr = ufield(1)
                        call locase(uheadr)
                        j2 = ilnobl(uheadx)
                        j = index(uheadr,uheadx(1:j2))

                        if (j .gt. 0) then
                            go to 240
                        end if

                        nesrbt(nsrt) = ncb
                    else
                        uheadx = 'volume'
                        j2 = ilnobl(uheadx)
                        j = index(uheadr,uheadx(1:j2))

                        uheadx = 'element'
                        nfldtx = 4
250 continue
                        call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

                        if (qrderr) then
                            go to 999
                        end if

                        uheadr = ufield(1)
                        call locase(uheadr)
                        j2 = ilnobl(uheadx)
                        j = index(uheadr,uheadx(1:j2))

                        if (j .gt. 0) then
                            go to 250
                        end if
                    end if

                    uheadx = 'DISSOLUTION LAW'
                    call locase(uheadx)
                    j2 = ilnobl(uheadx)
                    j = index(uheadr,uheadx(1:j2))

                    if (j .eq. 0) then
                        write (nttyo,1020) uheadx,uline1
                        go to 990
                    end if

                    ustr = ufield(2)
                    call chrint(ivar,nttyo,qrderr,ustr)

                    if (qrderr) then
                        go to 990
                    end if

                    nrk(1,nrc) = ivar

                    ! Rate law parameters, forward direction.
                    nfldtx = 4
                    call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

                    if (qrderr) then
                        go to 999
                    end if

                    if (nrk(1,nrc) .le. 0) then
                        imech(1,nrc) = 0
                        go to 290
                    end if

                    im = 0
260 continue
                    im = im + 1

                    uheadr = ufield(1)
                    call locase(uheadr)

                    uheadx = 'PRECIPITATION LAW'
                    call locase(uheadx)
                    j2 = ilnobl(uheadx)
                    j = index(uheadr,uheadx(1:j2))

                    if (j .gt. 0) then
                        im = im - 1
                        go to 290
                    end if

                    if (nrk(1,nrc) .le. 0) then
                        write (nttyo,1020) uheadx,uline1
                        go to 990
                    end if

                    uheadx = 'rate constant rk'
                    write (ux1,'(i1)') im
                    uheadx(17:17) = ux1(1:1)
                    j2 = ilnobl(uheadx)
                    j = index(uheadr,uheadx(1:j2))

                    if (j .eq. 0) then
                        write (nttyo,1020) uheadx,uline1
                        go to 990
                    end if

                    if (im .gt. imchmx) then
                        write (nttyo,1045) ureac(nrc),imchmx
1045 format(/' * Error - (XCON6/rd6d7) Have too many rate',/7x,'constants in the forward direction rate law for',/7x,'reactant "',a24,'". The code is only',/7x,'dimensioned for ',i2,' rate constants per rate law.',/7x,'Reduce the number of rate constants or increase the',/7x,'dimensioning parameter imchpa.')

                        go to 990
                    end if

                    i = im
                    ustr = ufield(2)
                    call chreal(nttyo,qrderr,ustr,var)

                    if (qrderr) then
                        go to 990
                    end if

                    rk0(i,1,nrc) = var

                    if (nrk(1,nrc) .eq. 2) then
                        ustr = ufield(4)
                        call chreal(nttyo,qrderr,ustr,var)

                        if (qrderr) then
                            go to 990
                        end if

                        csigma(i,1,nrc) = var
                    end if

                    if (nrk(1,nrc) .eq. 1) then
                        go to 280
                    end if

                    n = 0
                    uheadx = 'aqueous species'
                    nfldtx = 4
270 continue
                    call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

                    if (qrderr) then
                        go to 999
                    end if

                    uheadr = ufield(1)
                    call locase(uheadr)
                    j2 = ilnobl(uheadx)
                    j = index(uheadr,uheadx(1:j2))

                    if (j .gt. 0) then
                        n = n + 1

                        if (n .gt. ndctmx) then
                            write (nttyo,1046) ureac(nrc),i,ndctmx
1046 format(/' * Error - (XCON6/rd6d7) Have too many species',/7x,'in the activity product in term ',i2,' of the',/7x,'forward direction rate law for reactant',/7x,'"',a24,'". The code is only dimensioned for ',i3,/7x,'such species. Reduce the number of such species',/7x,'or increase the dimensioning parameter ndctpa.')

                            go to 990
                        end if

                        udac(n,i,1,nrc) = ufield(2)
                        ustr = ufield(4)
                        call chreal(nttyo,qrderr,ustr,var)

                        if (qrderr) then
                            go to 990
                        end if

                        cdac(n,i,1,nrc) = var
                        go to 270
                    end if

                    ndact(i,1,nrc) = n

                    uheadx = 'PRECIPITATION LAW'
                    call locase(uheadx)
                    j2 = ilnobl(uheadx)
                    j = index(uheadr,uheadx(1:j2))

                    if (j .gt. 0) then
                        go to 290
                    end if

                    uheadx = 'temperature (c)'
                    j2 = ilnobl(uheadx)
                    j = index(uheadr,uheadx(1:j2))

                    if (j .eq. 0) then
                        write (nttyo,1020) uheadx,uline1
                        go to 990
                    end if

                    ustr = ufield(2)
                    call chreal(nttyo,qrderr,ustr,var)

                    if (qrderr) then
                        go to 990
                    end if

                    trk0(i,1,nrc) = var

                    iact(i,1,nrc) = 0

                    nfldtx = 4
                    call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

                    if (qrderr) then
                        go to 999
                    end if

                    uheadr = ufield(1)
                    call locase(uheadr)

                    uheadx = 'PRECIPITATION LAW'
                    call locase(uheadx)
                    j2 = ilnobl(uheadx)
                    j = index(uheadr,uheadx(1:j2))

                    if (j .gt. 0) then
                        go to 290
                    end if

                    uheadx = 'rate constant rk'
                    call locase(uheadx)
                    j2 = ilnobl(uheadx)
                    j = index(uheadr,uheadx(1:j2))

                    if (j .gt. 0) then
                        go to 260
                    end if

                    uheadx = 'act. energy-kcal'
                    j2 = ilnobl(uheadx)
                    j = index(uheadr,uheadx(1:j2))

                    if (j .eq. 0) then
                        write (nttyo,1020) uheadx,uline1
                        go to 990
                    end if

                    ustr = ufield(2)
                    call chreal(nttyo,qrderr,ustr,var)

                    if (qrderr) then
                        go to 990
                    end if

                    eact(i,1,nrc) = var

                    if (var .gt. 0.) then
                        iact(i,1,nrc) = 1
                    end if

                    ustr = ufield(4)
                    call chreal(nttyo,qrderr,ustr,var)

                    if (qrderr) then
                        go to 990
                    end if

                    hact(i,1,nrc) = var

                    if (var .gt. 0.) then
                        iact(i,1,nrc) = 2
                    end if

280 continue
                    nfldtx = 4
                    call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

                    if (qrderr) then
                        go to 999
                    end if

                    go to 260

290 continue
                    imech(1,nrc) = im

                    ! Rate law parameters, backward direction.
                    ustr = ufield(2)
                    call chrint(ivar,nttyo,qrderr,ustr)

                    if (qrderr) then
                        go to 990
                    end if

                    nrk(2,nrc) = ivar

                    if (nrk(2,nrc) .le. 0) then
                        imech(2,nrc) = 0
                        go to 210
                    end if

                    im = 0
300 continue
                    im = im + 1

                    nfldtx = 0
                    call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

                    if (qrderr) then
                        go to 999
                    end if

                    ustr = ufield(1)

                    if (ustr(1:8) .eq. '--------') then
                        imech(2,nrc) = 0
                        go to 330
                    end if

                    nfldtx = 4

                    if (nfldt .ne. nfldtx) then
                        write (nttyo,1030) nfldt,nfldtx,uline1
                    end if

                    uheadr = ufield(1)
                    call locase(uheadr)

                    uheadx = 'REACTANT'
                    call locase(uheadx)
                    j2 = ilnobl(uheadx)
                    j = index(uheadr,uheadx(1:j2))

                    if (j .gt. 0) then
                        im = im - 1
                        imech(2,nrc) = im
                        go to 215
                    end if

                    uheadx = 'rate constant rk'
                    write (ux1,'(i1)') im
                    uheadx(17:17) = ux1(1:1)
                    j2 = ilnobl(uheadx)
                    j = index(uheadr,uheadx(1:j2))

                    if (j .eq. 0) then
                        write (nttyo,1020) uheadx,uline1
                        go to 990
                    end if

                    if (im .gt. imchmx) then
                        write (nttyo,1047) ureac(nrc),imchmx
1047 format(/' * Error - (XCON6/rd6d7) Have too many rate',/7x,'constants in the backward direction rate law for',/7x,'reactant "',a24,'". The code is only',/7x,'dimensioned for ',i2,' rate constants per rate law.',/7x,'Reduce the number of rate constants or increase the',/7x,'dimensioning parameter imchpa.')

                        go to 990
                    end if

                    i = im
                    ustr = ufield(2)
                    call chreal(nttyo,qrderr,ustr,var)

                    if (qrderr) then
                        go to 990
                    end if

                    rk0(i,2,nrc) = var

                    if (nrk(2,nrc) .eq. 2) then
                        ustr = ufield(4)
                        call chreal(nttyo,qrderr,ustr,var)

                        if (qrderr) then
                            go to 990
                        end if

                        csigma(i,2,nrc) = var
                    end if

                    if (nrk(2,nrc) .eq. 1) then
                        go to 210
                    end if

                    n = 0
                    uheadx = 'aqueous species'
                    nfldtx = 4
310 continue
                    call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

                    if (qrderr) then
                        go to 999
                    end if

                    uheadr = ufield(1)
                    call locase(uheadr)
                    j2 = ilnobl(uheadx)
                    j = index(uheadr,uheadx(1:j2))

                    if (j .gt. 0) then
                        n = n + 1

                        if (n .gt. ndctmx) then
                            write (nttyo,1048) ureac(nrc),i,ndctmx
1048 format(/' * Error - (XCON6/rd6d7) Have too many species',/7x,'in the activity product in term ',i2,' of the',/7x,'forward direction rate law for reactant',/7x,'"',a24,'". The code is only dimensioned for ',i3,/7x,'such species. Reduce the number of such species',/7x,'or increase the dimensioning parameter ndctpa.')

                            go to 990
                        end if

                        udac(n,i,2,nrc) = ufield(2)
                        ustr = ufield(4)
                        call chreal(nttyo,qrderr,ustr,var)

                        if (qrderr) then
                            go to 990
                        end if

                        cdac(n,i,2,nrc) = var
                        go to 310
                    end if

                    ndact(i,2,nrc) = n

                    ustr = ufield(1)

                    if (ustr(1:8) .eq. '--------') then
                        imech(2,nrc) = im
                        go to 330
                    end if

                    uheadx = 'temperature (c)'
                    j2 = ilnobl(uheadx)
                    j = index(uheadr,uheadx(1:j2))

                    if (j .eq. 0) then
                        write (nttyo,1020) uheadx,uline1
                        go to 990
                    end if

                    nfldtx = 4

                    if (nfldt .ne. nfldtx) then
                        write (nttyo,1030) nfldt,nfldtx,uline1
                    end if

                    ustr = ufield(2)
                    call chreal(nttyo,qrderr,ustr,var)

                    if (qrderr) then
                        go to 990
                    end if

                    trk0(i,2,nrc) = var

                    iact(i,2,nrc) = 0

                    nfldtx = 0
                    call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

                    if (qrderr) then
                        go to 999
                    end if

                    ustr = ufield(1)

                    if (ustr(1:8) .eq. '--------') then
                        imech(2,nrc) = im
                        go to 330
                    end if

                    uheadr = ufield(1)
                    call locase(uheadr)

                    uheadx = 'REACTANT'
                    call locase(uheadx)
                    j2 = ilnobl(uheadx)
                    j = index(uheadr,uheadx(1:j2))

                    if (j .gt. 0) then
                        imech(2,nrc) = im
                        go to 215
                    end if

                    uheadx = 'rate constant rk'
                    call locase(uheadx)
                    j2 = ilnobl(uheadx)
                    j = index(uheadr,uheadx(1:j2))

                    if (j .gt. 0) then
                        go to 300
                    end if

                    uheadx = 'act. energy-kcal'
                    j2 = ilnobl(uheadx)
                    j = index(uheadr,uheadx(1:j2))

                    if (j .eq. 0) then
                        write (nttyo,1020) uheadx,uline1
                        go to 990
                    end if

                    nfldtx = 4

                    if (nfldt .ne. nfldtx) then
                        write (nttyo,1030) nfldt,nfldtx,uline1
                    end if

                    ustr = ufield(2)
                    call chreal(nttyo,qrderr,ustr,var)

                    if (qrderr) then
                        go to 990
                    end if

                    eact(i,2,nrc) = var

                    if (var .gt. 0.) then
                        iact(i,2,nrc) = 1
                    end if

                    ustr = ufield(4)
                    call chreal(nttyo,qrderr,ustr,var)

                    if (qrderr) then
                        go to 990
                    end if

                    hact(i,2,nrc) = var

                    if (var .gt. 0.) then
                        iact(i,2,nrc) = 2
                    end if

                    nfldtx = 0
                    call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

                    if (qrderr) then
                        go to 999
                    end if

                    go to 300

330 continue

                    ! Options. These are iopt, iopr, and iodb option switches,
                    ! skipping ones which may be classified as development options.
                    ! The development options are read below, in a somewhat different
                    ! manner. The iopg option switches are read still further below,
                    ! in a manner similar to that employed here.
                    ! Note: iopt(1) = iopt1, etc.
                    nfldtx = 1
                    uheadx = 'options'
                    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

                    if (qrderr) then
                        go to 999
                    end if

                    i = 0
                    nfldtx = 1
                    call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

                    if (qrderr) then
                        go to 999
                    end if

                    ustr = ufield(1)

                    if (ustr(1:8) .eq. '--------') then
                        go to 420
                    end if

                    ulscr = ufield(1)

                    if (ulscr(1:2) .ne. '- ') then
                        write (nttyo,1050) uline1
1050 format(/' * Error - (XCON6/rd6d7) Found the line beginning',' with-',/7x,'"',a70,'"',/7x,'where an iopt, iopg, iopr, or iodb option switch header',/7x,'was expected.')

                        go to 990
                    end if

                    ! Have found an option header.
340 continue
                    i = i + 1
                    nmark = 0
                    nmarks = 0

                    ! Put the option header string minus the preceding "- " string
                    ! into the variable uheadr.
                    uheadr = ulscr(3:nlchmx)
                    call locase(uheadr)

                    ! Remove the succeeding " -" string, if any.
                    j = ilnobl(uheadr)

                    if (uheadr(j - 1:j) .eq. ' -') then
                        uheadr(j:j) = ' '
                    end if

                    ! Identify the corresponding option.
                    do 370 idescx = 1,nop6pa
                        uheadx = uopt6(idescx)
                        call locase(uheadx)

                        if (uheadr(1:40) .eq. uheadx(1:40)) then
                            if (idescx .ne. i) then
                                write (nttyo,1060) uopt6(idescx),uvar6(idescx),index6(idescx),i,idescx
1060 format(/' * Warning - (XCON6/rd6d7) The input for the',/7x,'option whose identifying string which begins with-',/7x,'"',a70,'"',/7x,'(',a4,'(',i2,')) is out of order, in place ',i3,' instead of place ',i3,'.')
                            end if

                            go to 380
                        end if

370 continue

                        write (nttyo,1070) uheadr
1070 format(/" * Error - (XCON6/rd6d7) Can't identify the following",/7x,'iopt, iopg, iopr, or iodb option header-',/7x,'"',a70,'".')

                        go to 990

380 continue

                        ! Read the first line after an option  header.
                        jdesc = 0
390 continue
                        jdesc = jdesc + 1
                        nfldtx = 1
                        call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

                        if (qrderr) then
                            go to 999
                        end if

                        ulscr = ufield(1)

                        if (ulscr(1:8) .eq. '--------') then
                            if (nmarks .gt. 1) then
                                write (nttyo,1080) uopt6(idescx),uvar6(idescx),index6(idescx)
1080 format(/' * Error - (XCON6/rd6d7) More than one choice',/7x,'was marked for the following option-',/7x,'"',a70,'"',/9x,a4,'(',i2,').')

                                go to 990
                            end if

                            go to 420
                        end if

                        ! Is the current line a new option header?
                        if (ulscr(1:2) .eq. '- ') then
                            if (nmarks .gt. 1) then
                                write (nttyo,1080) uopt6(idescx),uvar6(idescx),index6(idescx)
                                go to 990
                            end if

                            go to 340
                        end if

                        j = 1
                        nmark = 0

                        if (ulscr(1:1) .eq. '*') then
                            j = 2
                            nmark = 1
                            nmarks = nmarks + 1
                        end if

                        uheadr = ulscr(j:nlchmx)
                        call locase(uheadr)
                        call lejust(uheadr)

                        do 400 idesc = 1,nod6pa
                            if (iopti6(idesc) .eq. idescx) then
                                uheadx = udesc6(idesc)
                                call locase(uheadx)

                                if (uheadr(1:40) .eq. uheadx(1:40)) then
                                    go to 410
                                end if
                            end if

400 continue

                            write (nttyo,1090) uheadr,uopt6(idescx),uvar6(idescx),index6(idescx)
1090 format(/" * Error - (XCON6/rd6d7) Can't identify the following",/7x,'option string-',/7x,'"',a70,'".',/7x,'It was given for the following option-',/7x,'"',a70,'"',/9x,a4,'(',i2,').')

                            go to 990

410 continue
                            if (nmark .eq. 1) then
                                ustr = uvar6(idescx)
                                call locase(ustr)

                                if (ustr(1:4) .eq. 'iopt') then
                                    iopt(index6(idescx)) = ivalu6(idesc)
                                else if (ustr(1:4) .eq. 'iopg') then
                                    iopg(index6(idescx)) = ivalu6(idesc)
                                else if (ustr(1:4) .eq. 'iopr') then
                                    iopr(index6(idescx)) = ivalu6(idesc)
                                else if (ustr(1:4) .eq. 'iodb') then
                                    iodb(index6(idescx)) = ivalu6(idesc)
                                else
                                    write (nttyo,1100) uheadr,uvar6(idescx)
1100 format(/" * Error - (XCON6/rd6d7) Don't recognize the",/7x,'option type string "',a4,'", which was given for',/7x,'the following option-',/7x,'"',a70,'".')

                                    go to 990
                                end if
                            end if

                            go to 390

420 continue

                            ! Get the time frame flag (iopt(1)).
                            iopt(1) = 0

                            if (nrct .gt. 0) then
                                do 430 nrc = 1,nrct
                                    if (nrk(1,nrc) .ge. 2) then
                                        iopt(1) = 1
                                        go to 435
                                    end if

                                    if (nrk(2,nrc) .ge. 2) then
                                        iopt(1) = 1
                                        go to 435
                                    end if

430 continue

435 continue
                                end if

                                ! Development options.
                                uheadx = 'development options'
                                nfldtx = 1
                                call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

                                if (qrderr) then
                                    go to 999
                                end if

                                uheadr = ufield(1)
                                call locase(uheadr)
                                j2 = ilnobl(uheadx)
                                j = index(uheadr,uheadx(1:j2))

                                if (j .eq. 0) then
                                    write (nttyo,1020) uheadx,uline1
                                    go to 990
                                end if

                                i = 0
440 continue
                                i = i + 1
                                nfldtx = 1
                                call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

                                if (qrderr) then
                                    go to 999
                                end if

                                ulscr = ufield(1)

                                if (ulscr(1:8) .eq. '--------') then
                                    go to 470
                                end if

                                j = index(ulscr,' ')
                                j2 = j - 1

                                if (j2 .eq. 0) then
                                    j2 = 1
                                end if

                                ustr = ulscr(1:j2)
                                j1 = j + 1

                                if (j1 .gt. nlchmx) then
                                    j1 = nlchmx
                                end if

                                uheadr = ulscr(j1:nlchmx)
                                call locase(uheadr)
                                call lejust(uheadr)

                                do 450 idescx = 1,ndv6pa
                                    uheadx = udevl6(idescx)
                                    call locase(uheadx)

                                    if (uheadr(1:40) .eq. uheadx(1:40)) then
                                        go to 460
                                    end if

450 continue

                                    write (nttyo,1110) uheadr
1110 format(/" * Error - (XCON6/rd6d7) Can't identify the following",/7x,'development option header-',/7x,'"',a70,'".')

                                    go to 990

460 continue
                                    call chrint(ivar,nttyo,qrderr,ustr)

                                    if (qrderr) then
                                        go to 990
                                    end if

                                    ustr = udv6vr(idescx)
                                    call locase(ustr)

                                    if (ustr(1:4) .eq. 'iopt') then
                                        iopt(index6(idescx)) = ivar
                                    else if (ustr(1:4) .eq. 'iopr') then
                                        iopr(index6(idescx)) = ivar
                                    else if (ustr(1:4) .eq. 'iodb') then
                                        iodb(index6(idescx)) = ivar
                                    else if (ustr(1:4) .eq. 'iopg') then
                                        iopg(index6(idescx)) = ivar
                                    else
                                        write (nttyo,1120) uheadr,udv6vr(idescx)
1120 format(/" * Error - (XCON6/rd6d7) Don't recognize the",/7x,'option type string "',a4,'", which was  given for',/7x,'the following option-',/7x,'"',a70,'".')

                                        go to 990
                                    end if

                                    go to 440

470 continue

                                    ! Tolerances.
                                    uheadx = 'tolerances'
                                    nfldtx = 1
                                    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

                                    if (qrderr) then
                                        go to 999
                                    end if

                                    uheadr = ufield(1)
                                    call locase(uheadr)
                                    j2 = ilnobl(uheadx)
                                    j = index(uheadr,uheadx(1:j2))

                                    if (j .eq. 0) then
                                        write (nttyo,1020) uheadx,uline1
                                        go to 990
                                    end if

                                    i = 0
510 continue
                                    i = i + 1
                                    nfldtx = 0
                                    call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

                                    if (qrderr) then
                                        go to 999
                                    end if

                                    ustr = ufield(1)

                                    if (ustr(1:8) .eq. '--------') then
                                        go to 520
                                    end if

                                    nfldtx = 3

                                    if (nfldt .ne. nfldtx) then
                                        write (nttyo,1030) nfldt,nfldtx,uline1
                                    end if

                                    uheadr = ufield(1)
                                    call locase(uheadr)

                                    do 515 idescx = 1,nto6pa
                                        uheadx = utol6(idescx)
                                        call locase(uheadx)

                                        if (uheadr(1:32) .eq. uheadx(1:32)) then
                                            go to 517
                                        end if

515 continue

                                        write (nttyo,1130) uheadr
1130 format(/" * Error - (XCON6/rd6d7) Can't identify the following",/7x,'tolerance parameter header-',/7x,'"',a70,'".')

                                        go to 990

517 continue
                                        ustr = ufield(2)

                                        if (idescx .eq. 1) then
                                            ! Itermx.
                                            call chrint(ivar,nttyo,qrderr,ustr)

                                            if (qrderr) then
                                                go to 999
                                            end if

                                            itermx = ivar
                                        else if (idescx .eq. 2) then
                                            ! Dlzidp.
                                            call chreal(nttyo,qrderr,ustr,var)

                                            if (qrderr) then
                                                go to 999
                                            end if

                                            dlzidp = var
                                        else if (idescx .eq. 3) then
                                            ! Tolbt.
                                            call chreal(nttyo,qrderr,ustr,var)

                                            if (qrderr) then
                                                go to 999
                                            end if

                                            tolbt = var
                                        else if (idescx .eq. 4) then
                                            ! Toldl.
                                            call chreal(nttyo,qrderr,ustr,var)

                                            if (qrderr) then
                                                go to 999
                                            end if

                                            toldl = var
                                        else if (idescx .eq. 5) then
                                            ! Tolx.
                                            call chreal(nttyo,qrderr,ustr,var)

                                            if (qrderr) then
                                                go to 999
                                            end if

                                            tolx = var
                                        else if (idescx .eq. 6) then
                                            ! Tolsat.
                                            call chreal(nttyo,qrderr,ustr,var)

                                            if (qrderr) then
                                                go to 999
                                            end if

                                            tolsat = var
                                        else if (idescx .eq. 7) then
                                            ! Tolsst.
                                            call chreal(nttyo,qrderr,ustr,var)

                                            if (qrderr) then
                                                go to 999
                                            end if

                                            tolsst = var
                                        else if (idescx .eq. 8) then
                                            ! Sscrew1.
                                            call chreal(nttyo,qrderr,ustr,var)

                                            if (qrderr) then
                                                go to 999
                                            end if

                                            sscrew(1) = var
                                        else if (idescx .eq. 9) then
                                            ! Sscrew2.
                                            call chreal(nttyo,qrderr,ustr,var)

                                            if (qrderr) then
                                                go to 999
                                            end if

                                            sscrew(2) = var
                                        else if (idescx .eq. 10) then
                                            ! Sscrew3.
                                            call chreal(nttyo,qrderr,ustr,var)

                                            if (qrderr) then
                                                go to 999
                                            end if

                                            sscrew(3) = var
                                        else if (idescx .eq. 11) then
                                            ! Sscrew4.
                                            call chreal(nttyo,qrderr,ustr,var)

                                            if (qrderr) then
                                                go to 999
                                            end if

                                            sscrew(4) = var
                                        else if (idescx .eq. 12) then
                                            ! Sscrew5.
                                            call chreal(nttyo,qrderr,ustr,var)

                                            if (qrderr) then
                                                go to 999
                                            end if

                                            sscrew(5) = var
                                        else if (idescx .eq. 13) then
                                            ! Sscrew6.
                                            call chreal(nttyo,qrderr,ustr,var)

                                            if (qrderr) then
                                                go to 999
                                            end if

                                            sscrew(6) = var
                                        else if (idescx .eq. 14) then
                                            ! Zklogu.
                                            call chreal(nttyo,qrderr,ustr,var)

                                            if (qrderr) then
                                                go to 999
                                            end if

                                            zklogu = var
                                        else if (idescx .eq. 15) then
                                            ! Zklogl.
                                            call chreal(nttyo,qrderr,ustr,var)

                                            if (qrderr) then
                                                go to 999
                                            end if

                                            zklogl = var
                                        else if (idescx .eq. 16) then
                                            ! Zkfac.
                                            call chreal(nttyo,qrderr,ustr,var)

                                            if (qrderr) then
                                                go to 999
                                            end if

                                            zkfac = var
                                        else if (idescx .eq. 17) then
                                            ! Dlzmx1.
                                            call chreal(nttyo,qrderr,ustr,var)

                                            if (qrderr) then
                                                go to 999
                                            end if

                                            dlzmx1 = var
                                        else if (idescx .eq. 18) then
                                            ! Dlzmx2.
                                            call chreal(nttyo,qrderr,ustr,var)

                                            if (qrderr) then
                                                go to 999
                                            end if

                                            dlzmx2 = var
                                        else if (idescx .eq. 19) then
                                            ! Nordlm.
                                            call chrint(ivar,nttyo,qrderr,ustr)

                                            if (qrderr) then
                                                go to 999
                                            end if

                                            nordlm = ivar
                                        else if (idescx .eq. 20) then
                                            ! Ntrymx.
                                            call chrint(ivar,nttyo,qrderr,ustr)

                                            if (qrderr) then
                                                go to 999
                                            end if

                                            ntrymx = ivar
                                        else if (idescx .eq. 21) then
                                            ! Npslmx.
                                            call chrint(ivar,nttyo,qrderr,ustr)

                                            if (qrderr) then
                                                go to 999
                                            end if

                                            npslmx = ivar
                                        else if (idescx .eq. 22) then
                                            ! Nsslmx.
                                            call chrint(ivar,nttyo,qrderr,ustr)

                                            if (qrderr) then
                                                go to 999
                                            end if

                                            nsslmx = ivar
                                        else if (idescx .eq. 23) then
                                            ! Ioscan.
                                            call chrint(ivar,nttyo,qrderr,ustr)

                                            if (qrderr) then
                                                go to 999
                                            end if

                                            ioscan = ivar
                                        end if

                                        go to 510

520 continue

                                        ! Process the bottom half of the current output file.
                                        ! Title.
                                        n = 1

                                        do 530 n = 1,ntitmx + 1
                                            read (ninpts,1000,err=990) uline1
                                            call parsln(nfldmx,nfldt,nlchmx,ufield,uline1,ulscr)
                                            ustr = ufield(1)

                                            if (ustr(1:8) .eq. '--------') then
                                                go to 540
                                            end if

                                            utitl2(n) = ufield(1)
530 continue

                                            write (nttyo,1115) ntitmx
1115 format(/' * Error - (XCON6/rd6d7) Have too many lines in the',/7x,'secondary title. The code is only dimensioned for ',i4,/7x,'lines. Reduce the size of the title or increase the',/7x,'dimensioning parameter ntitpa.')

                                            go to 990

540 continue
                                            ntitl2 = n - 1

                                            ! Original temperature.
                                            uheadx = 'temperature (c)'
                                            nfldtx = 2
                                            call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

                                            if (qrderr) then
                                                go to 999
                                            end if

                                            ustr = ufield(2)
                                            call chreal(nttyo,qrderr,ustr,var)

                                            if (qrderr) then
                                                go to 990
                                            end if

                                            tempci = var

                                            ! Electrical imbalance.
                                            uheadx = 'electrical imbalance'
                                            nfldtx = 2
                                            call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

                                            if (qrderr) then
                                                go to 999
                                            end if

                                            ustr = ufield(2)
                                            call chreal(nttyo,qrderr,ustr,var)

                                            if (qrderr) then
                                                go to 990
                                            end if

                                            electr = var

                                            ! Number of aqueous basis species.
                                            uheadx = 'number of aqueous master species'
                                            nfldtx = 2
                                            call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

                                            if (qrderr) then
                                                go to 999
                                            end if

                                            ustr = ufield(2)
                                            call chrint(ivar,nttyo,qrderr,ustr)

                                            if (qrderr) then
                                                go to 990
                                            end if

                                            ksb = ivar

                                            ksq = ksb
                                            kct = ksb - 1

                                            ! Postion of last pure mineral.
                                            uheadx = 'position of last pure mineral'
                                            nfldtx = 2
                                            call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

                                            if (qrderr) then
                                                go to 999
                                            end if

                                            ustr = ufield(2)
                                            call chrint(ivar,nttyo,qrderr,ustr)

                                            if (qrderr) then
                                                go to 990
                                            end if

                                            kmt = ivar

                                            ! Postion of last solid solution.
                                            uheadx = 'position of last solid solution'
                                            nfldtx = 2
                                            call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

                                            if (qrderr) then
                                                go to 999
                                            end if

                                            ustr = ufield(2)
                                            call chrint(ivar,nttyo,qrderr,ustr)

                                            if (qrderr) then
                                                go to 990
                                            end if

                                            kxt = ivar

                                            ! Nxmod options.
                                            uheadx = 'suppressed species'
                                            nfldtx = 1
                                            call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

                                            if (qrderr) then
                                                go to 999
                                            end if

                                            nfldtx = 0
                                            call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

                                            if (qrderr) then
                                                go to 999
                                            end if

                                            n = 0
                                            ustr = ufield(1)
                                            call locase(ustr)

                                            if (ustr(1:8) .eq. '--------') then
                                                go to 570
                                            end if

                                            if (ustr(1:8) .eq. 'none') then
                                                go to 560
                                            end if

                                            nfldtx = 4

                                            if (nfldt .ne. nfldtx) then
                                                write (nttyo,1030) nfldt,nfldtx,uline1
                                            end if

550 continue
                                            n = n + 1

                                            if (n .gt. nxmdmx) then
                                                write (nttyo,1132) nxmdmx
1132 format(/' * Error - (XCON6/rd6d7) Have too many nxmod',/7x,'alter/suppress options. The code is only dimensioned',/7x,'for ',i3,' such options. Reduce the number of such',' options',/7x,'or increase the dimensioning parameter',' nxmdpa.')

                                                go to 990
                                            end if

                                            uxmd24(n) = ufield(1)

                                            ustr = ufield(2)
                                            call locase(ustr)

                                            if (ustr(1:7) .eq. 'aqueous') then
                                                jxmod(n) = 0
                                            else if (ustr(1:7) .eq. 'mineral') then
                                                jxmod(n) = 1
                                            else if (ustr(1:3) .eq. 'gas') then
                                                jxmod(n) = 2
                                            else if (ustr(1:14) .eq. 'solid solution') then
                                                jxmod(n) = 3
                                            else
                                                j2 = ilnobl(ustr)
                                                write (nttyo,1135) ustr(1:j2)
1135 format(/" * Error - (XCON6/rd6d7) Can't identify the",/7x,'following alter/suppress species type string- "',a,'".',/7x,'This must be one of "aqueous", "mineral", "gas",',/7x,' or "solid solution".')

                                                go to 990
                                            end if

                                            ustr = ufield(3)
                                            call locase(ustr)

                                            if (ustr(1:8) .eq. 'suppress') then
                                                kxmod(n) = -1
                                            else if (ustr(1:7) .eq. 'replace') then
                                                kxmod(n) = 0
                                            else if (ustr(1:8) .eq. 'augmentk') then
                                                kxmod(n) = 1
                                            else if (ustr(1:8) .eq. 'augmentg') then
                                                kxmod(n) = 2
                                            else
                                                j2 = ilnobl(ustr)
                                                write (nttyo,1140) ustr(1:j2)
1140 format(/" * Error - (XCON6/rd6d7) Can't identify the",/7x,'following alter/suppress option string- "',a,'".',/7x,'This must be one of "suppress", "replace", "augmentk",',/7x,' or "augmentg".')

                                                go to 990
                                            end if

                                            ustr = ufield(4)
                                            call chreal(nttyo,qrderr,ustr,var)

                                            if (qrderr) then
                                                go to 990
                                            end if

                                            xlkmod(n) = var

560 continue
                                            nfldtx = 0
                                            call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

                                            if (qrderr) then
                                                go to 999
                                            end if

                                            ustr = ufield(1)
                                            call locase(ustr)

                                            if (ustr(1:8) .eq. '--------') then
                                                go to 570
                                            end if

                                            nfldtx = 4

                                            if (nfldt .ne. nfldtx) then
                                                write (nttyo,1030) nfldt,nfldtx,uline1
                                            end if

                                            go to 550

570 continue
                                            nxmod = n

                                            ! Iopg options.
                                            ! Note: iopg(1) = iopg1, etc.
                                            uheadx = 'iopg options'
                                            nfldtx = 1
                                            call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

                                            if (qrderr) then
                                                go to 999
                                            end if

                                            i = 0
                                            nfldtx = 1
                                            call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

                                            if (qrderr) then
                                                go to 999
                                            end if

                                            ulscr = ufield(1)

                                            if (ulscr(1:8) .eq. '--------') then
                                                go to 720
                                            end if

                                            if (ulscr(1:2) .ne. '- ') then
                                                write (nttyo,1050) uline1
                                                go to 990
                                            end if

                                            ! Have found an option header.
640 continue
                                            i = i + 1
                                            nmark = 0
                                            nmarks = 0

                                            ! Put the option header string minus the preceding "- " string
                                            ! into the variable uheadr.
                                            uheadr = ulscr(3:nlchmx)
                                            call locase(uheadr)

                                            ! Remove the succeeding " -" string, if any.
                                            j = ilnobl(uheadr)

                                            if (uheadr(j - 1:j) .eq. ' -') then
                                                uheadr(j:j) = ' '
                                            end if

                                            ! Identify the corresponding option.
                                            do 670 idescx = 1,nop6pa
                                                uheadx = uopt6(idescx)
                                                call locase(uheadx)

                                                if (uheadr(1:40) .eq. uheadx(1:40)) then
                                                    !      This warning test is commented out.  The warning doesn't
                                                    !      really apply as the version 7 options are structured
                                                    !      (the strings and such for the iopg options are mixed in
                                                    !      with those of the iopt, iopr, and iodb options).
                                                    !      if (idescx .ne. i) then
                                                    !        write (nttyo,1060) uopt6(idescx),uvar6(idescx),
                                                    ! $      index6(idescx),i,idescx
                                                    !      endif
                                                    go to 680
                                                end if

670 continue

                                                write (nttyo,1070) uheadr
                                                go to 990

680 continue

                                                ! Read the first line after an option  header.
                                                jdesc = 0
690 continue
                                                jdesc = jdesc + 1
                                                nfldtx = 1
                                                call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

                                                if (qrderr) then
                                                    go to 999
                                                end if

                                                ulscr = ufield(1)

                                                if (ulscr(1:8) .eq. '--------') then
                                                    if (nmarks .gt. 1) then
                                                        write (nttyo,1080) uopt6(idescx),uvar6(idescx),index6(idescx)
                                                    end if

                                                    go to 720
                                                end if

                                                ! Is the current line a new option header?
                                                if (ulscr(1:2) .eq. '- ') then
                                                    if (nmarks .gt. 1) then
                                                        write (nttyo,1080) uopt6(idescx),uvar6(idescx),index6(idescx)
                                                        go to 990
                                                    end if

                                                    go to 640
                                                end if

                                                j = 1
                                                nmark = 0

                                                if (ulscr(1:1) .eq. '*') then
                                                    j = 2
                                                    nmark = 1
                                                    nmarks = nmarks + 1
                                                end if

                                                uheadr = ulscr(j:nlchmx)
                                                call locase(uheadr)
                                                call lejust(uheadr)

                                                do 700 idesc = 1,nod6pa
                                                    if (iopti6(idesc) .eq. idescx) then
                                                        uheadx = udesc6(idesc)
                                                        call locase(uheadx)

                                                        if (uheadr(1:40) .eq. uheadx(1:40)) then
                                                            go to 710
                                                        end if
                                                    end if

700 continue

                                                    write (nttyo,1090) uheadr,uopt6(idescx),uvar6(idescx),index6(idescx)
                                                    go to 990

710 continue
                                                    if (nmark .eq. 1) then
                                                        ustr = uvar6(idescx)
                                                        call locase(ustr)

                                                        if (ustr(1:4) .eq. 'iopt') then
                                                            iopt(index6(idescx)) = ivalu6(idesc)
                                                        else if (ustr(1:4) .eq. 'iopg') then
                                                            iopg(index6(idescx)) = ivalu6(idesc)
                                                        else if (ustr(1:4) .eq. 'iopr') then
                                                            iopr(index6(idescx)) = ivalu6(idesc)
                                                        else if (ustr(1:4) .eq. 'iodb') then
                                                            iodb(index6(idescx)) = ivalu6(idesc)
                                                        else
                                                            write (nttyo,1100) uheadr,uvar6(idescx)
                                                            go to 990
                                                        end if
                                                    end if

                                                    go to 690

720 continue

                                                    ! Balance totals.
                                                    uheadx = 'elements, moles and moles aqueous'
                                                    nfldtx = 1
                                                    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

                                                    if (qrderr) then
                                                        go to 999
                                                    end if

                                                    i = 0
730 continue
                                                    nfldtx = 0
                                                    call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

                                                    if (qrderr) then
                                                        go to 999
                                                    end if

                                                    ulscr = ufield(1)

                                                    if (ulscr(1:8) .eq. '--------') then
                                                        go to 750
                                                    end if

                                                    nfldtx = 3

                                                    if (nfldt .ne. nfldtx) then
                                                        write (nttyo,1030) nfldt,nfldtx,uline1
                                                    end if

                                                    i = i + 1

                                                    if (i .gt. nctmax) then
                                                        write (nttyo,1134) nctmax
1134 format(/' * Error - (XCON6/rd6d7) Have too many chemical',/7x,'elements present. The code is only dimensioned',/7x,'for ',i3,' elements. Reduce the number of elements',/7x,'or increase the dimensioning parameter nctpar.')

                                                        go to 990
                                                    end if

                                                    uelemb(i) = ufield(1)
                                                    ustr = ufield(2)
                                                    call chreal(nttyo,qrderr,ustr,var)

                                                    if (qrderr) then
                                                        go to 990
                                                    end if

                                                    mteb(i) = var
                                                    ustr = ufield(3)
                                                    call chreal(nttyo,qrderr,ustr,var)

                                                    if (qrderr) then
                                                        go to 990
                                                    end if

                                                    mteaqb(i) = var
                                                    go to 730

750 continue

                                                    ! Basis variable data.
                                                    uheadx = 'master species and logarithmic basis variables'
                                                    nfldtx = 1
                                                    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

                                                    if (qrderr) then
                                                        go to 999
                                                    end if

                                                    i = 0
760 continue
                                                    nfldtx = 0
                                                    call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

                                                    if (qrderr) then
                                                        go to 999
                                                    end if

                                                    ulscr = ufield(1)

                                                    if (ulscr(1:8) .eq. '--------') then
                                                        go to 770
                                                    end if

                                                    nfldtx = 3

                                                    if (nfldt .ne. nfldtx) then
                                                        write (nttyo,1030) nfldt,nfldtx,uline1
                                                    end if

                                                    i = i + 1

                                                    if (i .gt. kmax) then
                                                        write (nttyo,1137) kmax
1137 format(/' * Error - (XCON6/rd6d7) Have too many master',/7x,'variables. The code is only dimensioned for ',i3,/7x,'master variables. Reduce the number of such variables',/7x,'or increase the dimensioning parameter kpar.')

                                                        go to 990
                                                    end if

                                                    unrms(i) = ufield(1)
                                                    undms(i) = ufield(2)
                                                    ustr = ufield(3)
                                                    call chreal(nttyo,qrderr,ustr,var)

                                                    if (qrderr) then
                                                        go to 990
                                                    end if

                                                    zvclgi(i) = var
                                                    go to 760

770 continue
                                                    kdim = i

                                                    ! Physically removed system data.
                                                    uheadx = 'physically removed subsystem'
                                                    nfldtx = 1
                                                    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

                                                    if (qrderr) then
                                                        go to 999
                                                    end if

                                                    nprmn = 0
                                                    nprmx = 0
                                                    kprs = 0

                                                    nfldtx = 0
                                                    call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

                                                    if (qrderr) then
                                                        go to 999
                                                    end if

                                                    ustr = ufield(1)

                                                    if (ustr(1:8) .eq. '--------') then
                                                        go to 820
                                                    end if

                                                    if (ustr(1:4) .eq. 'none') then
                                                        nfldtx = 1
                                                        call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

                                                        if (qrderr) then
                                                            go to 999
                                                        end if

                                                        ustr = ufield(1)

                                                        if (ustr(1:8) .ne. '--------') then
                                                            write (nttyo,1040) uline1
                                                            go to 990
                                                        end if

                                                        go to 820
                                                    end if

                                                    nfldtx = 3

                                                    if (nfldt .ne. nfldtx) then
                                                        write (nttyo,1030) nfldt,nfldtx,uline1
                                                    end if

                                                    ulscr = ufield(1)
                                                    call locase(ulscr)

                                                    ! Pure minerals.
780 continue
                                                    ustr = ufield(1)

                                                    if (ustr(1:8) .eq. '--------') then
                                                        go to 820
                                                    end if

                                                    nfldtx = 3

                                                    if (nfldt .ne. nfldtx) then
                                                        write (nttyo,1030) nfldt,nfldtx,uline1
                                                    end if

                                                    ulscr = ufield(1)

                                                    if (ulscr(1:8) .ne. '        ') then
                                                        go to 790
                                                    end if

                                                    nprmn = nprmn + 1

                                                    if (nprmn .gt. nprsmx) then
                                                        write (nttyo,1150) nprsmx
1150 format(/' * Error - (XCON6/rd6d7) Have too many mineral',/7x,'species in the physically removed system. The code is',/7x,'only dimensioned for ',i3,' such species. Reduce the',/7x,'number of such species or increase the dimensioning',/7x,'parameter nprspa.')

                                                        go to 990
                                                    end if

                                                    uprs(nprmn)(1:24) = ufield(2)
                                                    uprs(nprmn)(25:48) = ufield(2)
                                                    ustr = ufield(3)
                                                    call chreal(nttyo,qrderr,ustr,var)

                                                    if (qrderr) then
                                                        go to 990
                                                    end if

                                                    mprs(nprmn) = var

                                                    nfldtx = 0
                                                    call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

                                                    if (qrderr) then
                                                        go to 999
                                                    end if

                                                    go to 780

                                                    ! Solid solutions.
790 continue
                                                    nprmx = nprmn

800 continue
                                                    ux = ufield(1)

810 continue
                                                    nfldtx = 0
                                                    call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

                                                    if (qrderr) then
                                                        go to 999
                                                    end if

                                                    ustr = ufield(1)

                                                    if (ustr(1:8) .eq. '--------') then
                                                        go to 820
                                                    end if

                                                    nfldtx = 3

                                                    if (nfldt .ne. nfldtx) then
                                                        write (nttyo,1030) nfldt,nfldtx,uline1
                                                    end if

                                                    ulscr = ufield(1)

                                                    if (ulscr(1:8) .ne. '        ') then
                                                        go to 800
                                                    end if

                                                    nprmx = nprmx + 1

                                                    if (nprmn .gt. nprsmx) then
                                                        write (nttyo,1150) nprsmx
                                                        go to 990
                                                    end if

                                                    uprs(nprmx)(1:24) = ux
                                                    uprs(nprmx)(25:48) = ufield(2)
                                                    ustr = ufield(3)
                                                    call chreal(nttyo,qrderr,ustr,var)

                                                    if (qrderr) then
                                                        go to 990
                                                    end if

                                                    mprs(nprmx) = var
                                                    go to 810

820 continue
                                                    n = nprmn + nprmx

                                                    if (n .gt. 0) then
                                                        kprs = 1
                                                    end if

                                                    go to 999

990 continue
                                                    qrderr = .true.

999 continue
                                                end subroutine rd6d7