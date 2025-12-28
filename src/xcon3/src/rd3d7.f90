subroutine rd3d7(cspb,fep,iktmax,iodb,iopg,iopr,iopt,itermx,jflagb,jxmod,kxmod,ncompb,ninpts,nodbmx,nopgmx,noprmx,noptmx,nsq,nsqmax,ntitl,ntitmx,nttyo,nxmdmx,nxmod,nxtb,nxtmax,qend,qrderr,rho,tempc,tdspkg,tdspl,tolbt,toldl,tolsat,ubasis,uebal,umemb,uphas1,uphas2,uredox,usolb,uspecb,utitl,xbarb,uxmd24,xlkmod)
    !! This subroutine reads the EQ3NR input file in menu-style ("D")
    !! format for versions 7.0-7.2. It thus encompasses two version
    !! levels.
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
    integer :: ninpts
    integer :: nsq
    integer :: ntitl
    integer :: nttyo
    integer :: nxmod
    integer :: nxtb

    logical :: qend
    logical :: qrderr

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

    ! Local parameter declarations.
    !   nfldpa = maximum number of fields per line
    !   nlchpa = character length of a line
    integer :: nfldpa
    integer :: nlchpa

    parameter (nfldpa = 8,nlchpa = 80)

    ! Local variable declarations.
    integer :: i
    integer :: idesc
    integer :: idescx
    integer :: iktb
    integer :: ivar
    integer :: j
    integer :: jdesc
    integer :: j1
    integer :: j2
    integer :: j3
    integer :: n
    integer :: ncount
    integer :: nfldmx
    integer :: nfldt
    integer :: nfldtx
    integer :: nlchmx
    integer :: nmark
    integer :: nmarks
    integer :: ilnobl

    logical :: qrdxcp

    character(len=nlchpa) :: ufield(nfldpa)
    character(len=nlchpa) :: uline1
    character(len=nlchpa) :: uline2
    character(len=nlchpa) :: ulscr
    character(len=nlchpa) :: uheadr
    character(len=nlchpa) :: uheadx

    character(len=80) :: ustr

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
        j2 = ilnobl(uline1)
        write (nttyo,1010) uline1(1:j2)
1010 format(/' * Error - (XCON3/rd3d7) The first line of a "D"',/7x,'format input file must begin with "|--------".',/7x,'The first line read from the old input file begins with-',/7x,'"',a,'".')

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

        utitl(n) = ufield(1)
110 continue

        write (nttyo,1015) ntitmx
1015 format(/' * Error - (XCON3/rd3d7) Have too many lines in the',/7x,'main title. The code is only dimensioned for ',i4,/7x,'lines. Reduce the size of the title or increase the',/7x,'dimensioning parameter ntitpa.')

        go to 990

120 continue
        ntitl = n - 1

        ! Temperature and density.
        tempc = 0.
        rho = 0.

        uheadx = 'Temperature (C)'
        nfldtx = 4
        call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

        if (qrderr) then
            go to 999
        end if

        ustr = ufield(2)
        call chreal(nttyo,qrderr,ustr,var)

        if (qrderr) then
            go to 999
        end if

        tempc = var
        ustr = ufield(4)
        call chreal(nttyo,qrderr,ustr,var)

        if (qrderr) then
            go to 999
        end if

        rho = var

        ! Total dissolved salts.
        tdspkg = 0.
        tdspl = 0.

        uheadx = 'Total dissolved salts'
        nfldtx = 5
        call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

        if (qrderr) then
            go to 999
        end if

        nmark = 0
        ncount = 0

        do 125 n = 3,5
            j = index(ufield(n),'*')

            if (j .gt. 0) then
                ncount = ncount + 1

                if (ncount .eq. 1) then
                    nmark = n
                end if
            end if

125 continue

            if (ncount .eq. 0) then
                j2 = ilnobl(uline1)
                write (nttyo,1017) uline1(1:j2)
1017 format(/' * Warning - (XCON3/rd3d7) None of the options'  /7x,'is marked with an asterisk on the line beginning with',/7x,'"',a,'".')
            end if

            if (ncount .gt. 1) then
                j2 = ilnobl(uline1)
                write (nttyo,1018) uline1(1:j2)
1018 format(/' * Warning - (XCON3/rd3d7) More than one of the'  /7x,'options is marked with an asterisk on the line',' beginning with',/7x,'"',a,'".')
            end if

            ustr = ufield(2)
            call chreal(nttyo,qrderr,ustr,var)

            if (qrderr) then
                go to 999
            end if

            if (nmark .eq. 3) then
                tdspkg = var
            else if (nmark .eq. 4) then
                tdspl = var
            end if

            ! Electrical balancing species.
            uebal = ' '

            uheadx = 'Electrical balancing'
            nfldtx = 4
            call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

            if (qrderr) then
                go to 999
            end if

            nmark = 0
            ncount = 0

            do 127 n = 3,4
                j = index(ufield(n),'*')

                if (j .gt. 0) then
                    ncount = ncount + 1

                    if (ncount .eq. 1) then
                        nmark = n
                    end if
                end if

127 continue

                if (ncount .gt. 1) then
                    write (nttyo,1018) uline1
                end if

                uebal = ufield(2)

                if (nmark .eq. 3) then
                    uebal = 'pick1.'
                end if

                if (nmark .eq. 4) then
                    uebal = ' '
                end if

                ! Basis species and associated constraints. If a solid solution
                ! end-member is part of a constraint, read the solid solution
                ! composition here. Other solid solution compositions are read
                ! in the following block.
                uheadx = 'SPECIES'
                nfldtx = 4
                call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

                if (qrderr) then
                    go to 999
                end if

                nxtb = 0
                uredox = ' '
                fep = 0.
                qrdxcp = .false.

                nsq = 0
                nfldtx = 0
                call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

                if (qrderr) then
                    go to 999
                end if

130 continue
                ustr = ufield(1)

                if (ustr(1:8) .eq. '--------') then
                    go to 180
                end if

                nfldtx = 4

                if (nfldt .ne. nfldtx) then
                    j2 = ilnobl(uline1)
                    write (nttyo,1030) nfldt,nfldtx,uline1(1:j2)
1030 format(/' * Warning - (XCON3/rd3d7) Found ',i2,' fields',/7x,'where ',i2,' were expected on the line beginning with-',/7x,'"',a,'".')
                end if

                uheadr = ufield(1)
                call locase(uheadr)

                ustr = ufield(3)
                call chreal(nttyo,qrderr,ustr,var)

                if (qrderr) then
                    go to 999
                end if

                ustr = ufield(4)
                call locase(ustr)

                if (uheadr(1:8) .eq. 'redox') then
                    fep = var

                    if (ustr(1:6) .eq. 'logfo2') then
                        iopt(1) = 0
                    else if (ustr(1:2) .eq. 'eh') then
                        iopt(1) = -1
                    else if (ustr(1:2) .eq. 'pe') then
                        iopt(1) = -2
                    else if (ustr(1:12) .eq. 'redox couple') then
                        iopt(1) = 1
                        qrdxcp = .true.
                    else
                        j2 = ilnobl(ustr)
                        write (nttyo,1040) ustr(1:j2)
1040 format(/" * Error - (XCON3/rd3d7) Can't identify the",/7x,'following redox option string- "',a,'". This must',/7x,'be one of "LogfO2", "Eh", "pe", or "redox couple".')

                        go to 990
                    end if

                    nfldtx = 0
                    call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

                    if (qrderr) then
                        go to 999
                    end if

                    go to 130
                end if

                nsq = nsq + 1

                if (nsq .gt. nsqmax) then
                    write (nttyo,1050) nsqmax
1050 format(/' * Error - (XCON3/rd3d7) Have too many basis',/7x,'species. The code is only dimensioned for ',i3,/7x,'such species. Reduce the number of such species',/7x,'or increase the dimensioning parameter nsqpar.')

                    go to 990
                end if

                uspecb(nsq) = ufield(1)

                if (qrdxcp) then
                    uredox = ufield(1)
                end if

                qrdxcp = .false.

                if (uspecb(nsq)(1:5).eq.'o2(g)' .or.    uspecb(nsq)(1:5).eq.'O2(g)') then
                    iopt(1) = -3
                    uredox = ' '
                end if

                cspb(nsq) = var

                do 140 n = -1,njf7pa
                    uheadx = ujflg7(n)
                    call locase(uheadx)

                    if (ustr(1:16) .eq. uheadx(1:16)) then
                        jflagb(nsq) = n
                        go to 150
                    end if

140 continue

                    j2 = ilnobl(ustr)
                    write (nttyo,1060) ustr(1:j2)
1060 format(/" * Error - (XCON3/rd3d7) Can't identify the",/7x,'following jflag option string- "',a,'". This should',/7x,'be one of the strings defined in the ujflg7 array.')

                    go to 990

150 continue
                    ustr = ufield(2)

                    if (jflagb(nsq).ge.17 .and. jflagb(nsq).le.21) then
                        uphas1(nsq) = ustr
                    else
                        ubasis(nsq) = ustr
                    end if

                    if (jflagb(nsq) .eq. 20) then
                        nxtb = nxtb + 1

                        if (nxtb .gt. nxtmax) then
                            write (nttyo,1310) nxtmax
1310 format(/' * Error - (XCON3/rd3d7) Have too many solid',/7x,'solutions present. The code is only dimensioned',/7x,'for ',i3,' solid solutions. Reduce the number of such',/7x,'phases or increase the dimensioning parameter nxtpar.')

                            go to 990
                        end if

                        usolb(nxtb) = ustr
                        iktb = 0
160 continue
                        nfldtx = 0
                        call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

                        if (qrderr) then
                            go to 999
                        end if

                        if (ufield(1)(1:8) .eq. '        ') then
                            nfldtx = 4

                            if (nfldt .ne. nfldtx) then
                                write (nttyo,1030) nfldt,nfldtx,uline1
                            end if

                            iktb = iktb + 1

                            if (iktb .eq. 1) then
                                uphas2(nsq) = ufield(2)
                            end if

                            umemb(iktb,nxtb) = ufield(2)
                            ustr = ufield(3)
                            call chreal(nttyo,qrderr,ustr,var)

                            if (qrderr) then
                                go to 999
                            end if

                            xbarb(iktb,nxtb) = var
                            go to 160
                        else
                            ncompb(nxtb) = iktb
                            go to 130
                        end if
                    end if

                    nfldtx = 0
                    call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

                    if (qrderr) then
                        go to 999
                    end if

                    go to 130

180 continue

                    ! Mole fractions of solid solutions.
                    uheadx = 'input solid solutions'
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

                    ustr = ufield(1)
                    call locase(ustr)

                    if (ustr(1:8) .eq. '--------') then
                        go to 260
                    end if

                    nfldtx = 4

                    if (nfldt .ne. nfldtx) then
                        write (nttyo,1030) nfldt,nfldtx,uline1
                    end if

                    if (ustr(1:4) .eq. 'none') then
                        nfldtx = 1
                        call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

                        if (qrderr) then
                            go to 999
                        end if

                        ustr = ufield(1)
                        call locase(ustr)

                        if (ustr(1:8) .ne. '--------') then
                            j2 = ilnobl(uline1)
                            write (nttyo,1305) uline1(1:j2)
1305 format(/' * Error - (XCON3/rd3d7) Found the line beginning',' with-',/7x,'"',a,'"',/7x,'where a dashed separator line was expected.')

                            go to 990
                        end if

                        go to 260
                    end if

220 continue
                    nxtb = nxtb + 1

                    if (nxtb .gt. nxtmax) then
                        write (nttyo,1310) nxtmax
                        go to 990
                    end if

                    usolb(nxtb) = ufield(1)
                    iktb = 0

230 continue
                    nfldtx = 0
                    call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

                    if (qrderr) then
                        go to 999
                    end if

                    ustr = ufield(1)
                    call locase(ustr)

                    if (ustr(1:8) .ne. '        ') then
                        ncompb(nxtb) = iktb

                        if (ustr(1:8) .eq. '--------') then
                            go to 260
                        end if

                        go to 220
                    end if

                    nfldtx = 4

                    if (nfldt .ne. nfldtx) then
                        write (nttyo,1030) nfldt,nfldtx,uline1
                    end if

                    ustr = ufield(1)

                    if (ustr(1:8) .ne. '        ') then
                        ncompb(nxtb) = iktb
                        go to 220
                    end if

                    iktb = iktb + 1

                    if (iktb .gt. iktmax) then
                        j2 = ilnobl(usolb(nxtb))
                        write (nttyo,1330) usolb(nxtb)(1:j2),iktmax
1330 format(/' * Error - (XCON3/rd3d7) Solid solution',/7x,'"',a,'" has too many end-members present.',/7x,'This code is only dimensioned for ',i3,' end-members',/7x,'per solid solution. Reduce the number of end-members',/7x,'or increase the dimensioning parameter iktpar.')

                        go to 990
                    end if

                    umemb(iktb,nxtb) = ufield(2)
                    ustr = ufield(3)
                    call chreal(nttyo,qrderr,ustr,var)

                    if (qrderr) then
                        go to 999
                    end if

                    xbarb(iktb,nxtb) = var

                    go to 230

260 continue

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

                    if (ustr(1:4) .eq. 'none') then
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
1132 format(/' * Error - (XCON3/rd3d7) Have too many nxmod',/7x,'alter/suppress options. The code is only dimensioned',/7x,'for ',i3,' such options. Reduce the number of such',' options',/7x,'or increase the dimensioning parameter',' nxmdpa.')

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
1135 format(/" * Error - (XCON3/rd3d7) Can't identify the",/7x,'following alter/suppress species type string- "',a,'".',/7x,'This must be one of "aqueous", "mineral", "gas",',/7x,' or "solid solution".')

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
1140 format(/" * Error - (XCON3/rd3d7) Can't identify the",/7x,'following alter/suppress option string- "',a,'".',/7x,'This must be one of "suppress", "replace", "augmentk",',/7x,' or "augmentg".')

                        go to 990
                    end if

                    ustr = ufield(4)
                    call chreal(nttyo,qrderr,ustr,var)

                    if (qrderr) then
                        go to 999
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

                    ! Options. These are iopt, iopg, and iopr option switches,
                    ! skipping ones which may be classified as development options.
                    ! The development options are read below, in a somewhat different
                    ! manner. The iodb option switches are read still further below,
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
                        j = ilnobl(uline1(1:70))
                        write (nttyo,1250) uline1(1:j)
1250 format(/' * Error - (XCON3/rd3d7) Found the line beginning',' with-',/7x,'"',a,'"',/7x,'where an iopt, iopg, iopr, or iodb option switch header',/7x,'was expected.')

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
                    do 370 idescx = 1,nop3pa
                        uheadx = uopt3(idescx)
                        call locase(uheadx)

                        if (uheadr(1:40) .eq. uheadx(1:40)) then
                            if (idescx .ne. i) then
                                j2 = ilnobl(uopt3(idescx))
                                j3 = ilnobl(uvar3(idescx))
                                write (nttyo,1260) uopt3(idescx)(1:j2),uvar3(idescx)(1:j3),index3(idescx),i,idescx
1260 format(/' * Warning - (XCON3/rd3d7) The input for the',/7x,'option whose identifying string which begins with-',/7x,'"',a,'"',/7x,'(',a,'(',i2,')) is out of order,',' in place ',i3,' instead of place ',i3,'.')
                            end if

                            go to 380
                        end if

370 continue

                        j2 = ilnobl(uheadr)
                        write (nttyo,1270) uheadr(1:j2)
1270 format(/" * Error - (XCON3/rd3d7) Can't identify the following",/7x,'iopt, iopg, iopr, or iodb option header-',/7x,'"',a,'".')

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
                                j2 = ilnobl(uopt3(idescx))
                                j3 = ilnobl(uvar3(idescx))
                                write (nttyo,1280) uopt3(idescx)(1:j2),uvar3(idescx)(1:j3),index3(idescx)
1280 format(/' * Error - (XCON3/rd3d7) More than one choice',/7x,'was marked for the following option-',/7x,'"',a,'"',/9x,a,'(',i2,').')

                                go to 990
                            end if

                            go to 420
                        end if

                        ! Is the current line a new option header?
                        if (ulscr(1:2) .eq. '- ') then
                            if (nmarks .gt. 1) then
                                j2 = ilnobl(uopt3(idescx))
                                j3 = ilnobl(uvar3(idescx))
                                write (nttyo,1280) uopt3(idescx)(1:j2),uvar3(idescx)(1:j3),index3(idescx)
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

                        do 400 idesc = 1,nod3pa
                            if (iopti3(idesc) .eq. idescx) then
                                uheadx = udesc3(idesc)
                                call locase(uheadx)

                                if (uheadr(1:40) .eq. uheadx(1:40)) then
                                    go to 410
                                end if
                            end if

400 continue

                            j = ilnobl(uheadr)
                            j2 = ilnobl(uopt3(idescx))
                            j3 = ilnobl(uvar3(idescx))
                            write (nttyo,1290) uheadr(1:j),uopt3(idescx)(1:j2),uvar3(idescx)(1:j3),index3(idescx)
1290 format(/" * Error - (XCON3/rd3d7) Can't identify the following",/7x,'option string-',/7x,'"',a,'".',/7x,'It was given for the following option-',/7x,'"',a,'"',/9x,a,'(',i2,').')

                            go to 990

410 continue
                            if (nmark .eq. 1) then
                                ustr = uvar3(idescx)
                                call locase(ustr)

                                if (ustr(1:4) .eq. 'iopt') then
                                    iopt(index3(idescx)) = ivalu3(idesc)
                                else if (ustr(1:4) .eq. 'iopg') then
                                    iopg(index3(idescx)) = ivalu3(idesc)
                                else if (ustr(1:4) .eq. 'iopr') then
                                    iopr(index3(idescx)) = ivalu3(idesc)
                                else if (ustr(1:4) .eq. 'iodb') then
                                    iodb(index3(idescx)) = ivalu3(idesc)
                                else
                                    j = ilnobl(uheadr)
                                    j3 = ilnobl(uvar3(idescx))
                                    write (nttyo,1300) uvar3(idescx)(1:j3),uheadr(1:j)
1300 format(/" * Error - (XCON3/rd3d7) Don't recognize the",/7x,'option type string "',a,'", which was given for',/7x,'the following option-',/7x,'"',a,'".')

                                    go to 990
                                end if
                            end if

                            go to 390

420 continue

                            ! Iodb options.
                            ! Note: iodb(1) = iodb, etc.
                            nfldtx = 1
                            uheadx = 'debugging switches'
                            call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

                            if (qrderr) then
                                go to 999
                            end if

                            uheadr = ufield(1)
                            call locase(uheadr)
                            j2 = ilnobl(uheadx)
                            j = index(uheadr,uheadx(1:j2))

                            if (j .eq. 0) then
                                j2 = ilnobl(uheadx)
                                j3 = ilnobl(uline1)
                                write (nttyo,1020) uheadx(1:j2),uline1(1:j3)
1020 format(/' * Error - (XCON3/rd3d7) Was expecting to find the',/7x,'header beginning with-',/7x,'"',a,'"',/7x,'on the line beginning with-',/7x,'"',a,'".')

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

                            do 450 idescx = 1,ndb3pa
                                uheadx = udebug(idescx)
                                call locase(uheadx)

                                if (uheadr(1:40) .eq. uheadx(1:40)) then
                                    go to 460
                                end if

450 continue

                                j2 = ilnobl(uheadr)
                                write (nttyo,1110) uheadr(1:j2)
1110 format(/" * Error - (XCON3/rd3d7) Can't identify the following",/7x,'iodb option header-',/7x,'"',a,'".')

                                go to 990

460 continue
                                call chrint(ivar,nttyo,qrderr,ustr)

                                if (qrderr) then
                                    go to 999
                                end if

                                iodb(idbugi(idescx)) = ivar
                                go to 440

470 continue

                                ! Development options. There are none.
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

                                nfldtx = 1
                                uheadx = 'none'
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

                                do 515 idescx = 1,nto3pa
                                    uheadx = utol3(idescx)
                                    call lejust(uheadx)
                                    call locase(uheadx)

                                    if (uheadr(1:32) .eq. uheadx(1:32)) then
                                        go to 517
                                    end if

515 continue

                                    j2 = ilnobl(uheadr)
                                    write (nttyo,1130) uheadr(1:j2)
1130 format(/" * Error - (XCON3/rd3d7) Can't identify the following",/7x,'tolerance parameter header-',/7x,'"',a,'".')

                                    go to 990

517 continue
                                    ustr = ufield(2)

                                    if (idescx .eq. 1) then
                                        ! Tolbt.
                                        call chreal(nttyo,qrderr,ustr,var)

                                        if (qrderr) then
                                            go to 999
                                        end if

                                        tolbt = var
                                    else if (idescx .eq. 2) then
                                        ! Toldl.
                                        call chreal(nttyo,qrderr,ustr,var)

                                        if (qrderr) then
                                            go to 999
                                        end if

                                        toldl = var
                                    else if (idescx .eq. 3) then
                                        ! Tolsat.
                                        call chreal(nttyo,qrderr,ustr,var)

                                        if (qrderr) then
                                            go to 999
                                        end if

                                        tolsat = var
                                    else if (idescx .eq. 4) then
                                        ! Itermx.
                                        call chrint(ivar,nttyo,qrderr,ustr)

                                        if (qrderr) then
                                            go to 999
                                        end if

                                        itermx = ivar
                                    end if

                                    go to 510

520 continue

                                    go to 999

990 continue
                                    qrderr = .true.

999 continue
                                end subroutine rd3d7