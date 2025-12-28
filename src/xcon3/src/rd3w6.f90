subroutine rd3w6(cspb,fep,iktmax,iodb,iopg,iopr,iopt,itermx,jflagb,jxmod,kxmod,ncompb,ninpts,nodbmx,nopgmx,noprmx,noptmx,nsq,nsqmax,ntitl,ntitmx,nttyo,nxmdmx,nxmod,nxtb,nxtmax,qend,qrderr,rho,tempc,tdspkg,tdspl,tolbt,toldl,tolsat,uacion,ubasis,uebal,umemb,uphas1,uphas2,uredox,usolb,uspecb,utitl,xbarb,uxmd24,xlkmod)
    !! This subroutine reads the EQ3NR input file in compact ("W")
    !! format for versions 6.0-6.1.
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
    character(len=24) :: uacion
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

    ! Local variable declarations.
    integer :: i
    integer :: iktb
    integer :: j
    integer :: j2
    integer :: n
    integer :: ilnobl

    character(len=80) :: uline
    character(len=24) :: ux24
    character(len=8) :: uendit
    character(len=8) :: ux8

    real(kind=8) :: xx

    data uendit /'endit.  '/

    qrderr = .false.

    ! Main title.
    ! Note: if there are exactly ntitpa lines in the title, the title
    ! need not be terminated by an 'endit.'. The 'endit.' if present
    ! is here not considered to be part of the title. The secondary
    ! title is treated in the same manner.
    qend = .false.
    read (ninpts,1000,end=100,err=990) uline
1000 format(a80)

    go to 105

100 continue
    qend = .true.
    go to 999

105 continue
    n = 1

    if (uline(1:8) .eq. uendit(1:8)) then
        go to 120
    end if

    utitl(1) = uline

    do 110 n = 2,ntitmx
        read (ninpts,1000,err=990) uline

        if (uline(1:8) .eq. uendit(1:8)) then
            go to 120
        end if

        utitl(n) = uline
110 continue

        n = n + 1

        read (ninpts,1000,err=990) uline
        call locase(uline)
        j = index(uline,'tempc=')

        if (j .gt. 0) then
            backspace(ninpts)
        else
            write (nttyo,1030) ntitmx
1030 format(/' * Error - (XCON3/rd3w6) Have too many lines in the',/7x,'secondary title. The code is only dimensioned for ',i4,/7x,'lines. Reduce the size of the title or increase the',/7x,'dimensioning parameter ntitpa.')

            go to 990
        end if

120 continue
        ntitl = n - 1

        ! Temperature.
        read (ninpts,1040,err=990) tempc
1040 format(12x,e12.5)

        ! Density, total dissolved salts (per kg), and total dissolved
        ! salts (per liter).
        read (ninpts,1050,err=990) rho,tdspkg,tdspl
1050 format (3(12x,e12.5))

        ! Redox parameter and name of redox species defining a controlling
        ! couple, if any.
        read (ninpts,1060,err=990) fep,uredox
1060 format (12x,e12.5,12x,a24)

        ! Convergence tolerances (tolbt and toldl) and saturation
        ! tolerance (tolsat).
        read (ninpts,1070,err=990) tolbt,toldl,tolsat
1070 format (3(12x,e12.5))

        ! Maximum number of iterations.
        read (ninpts,1080,err=990) itermx
1080 format (12x,i2)

        ! Iopt option switches.
        ! Note: iopt(1) = iopt1, etc.
        read (ninpts,1100,err=990) (iopt(i), i = 1,10)
1100 format(12x,10i5)

        ! Iopg option switches.
        ! Note: iopg(1) = iopg1, etc.
        read (ninpts,1110,err=990) (iopg(i), i = 1,10)
1110 format(12x,10i5)

        ! Iopr option switches.
        ! Note: iopr(1) = iopr1, etc.
        read (ninpts,1120,err=990) (iopr(i), i = 1,20)
1120 format(12x,10i5)

        ! Iodb option switches.
        ! Note: iodb(1) = iodb1, etc.
        read (ninpts,1130,err=990) (iodb(i), i = 1,10)
1130 format(12x,10i5)

        ! Species for electrical balancing.
        read (ninpts,1150,err=990) uebal
1150 format(12x,a24)

        ! Species for defining equivalent stoichiometric ionic strength.
        read (ninpts,1160,err=990) uacion
1160 format(12x,a24)

        ! Number of nxmod options.
        read (ninpts,1200,err=990) nxmod
1200 format(12x,i2)

        if (nxmod .gt. nxmdmx) then
            write (nttyo,1210) nxmdmx
1210 format(/' * Error - (XCON3/rd3w6) Have too many nxmod',/7x,'alter/suppress options. The code is only dimensioned',/7x,'for ',i3,' such options. Reduce the number of such',' options',/7x,'or increase the dimensioning parameter',' nxmdpa.')

            go to 990
        end if

        ! Nxmod options.
        if (nxmod .gt. 0) then
            do 420 n = 1,nxmod
                read (ninpts,1220,err=990) uxmd24(n),jxmod(n),kxmod(n),xlkmod(n)
1220 format(12x,a24,/12x,i2,22x,i2,22x,e12.5)

420 continue
            end if

            ! Basis species and associated constraints.
            nsq = 0

430 continue
            read (ninpts,1250,err=990) ux8,ux24
1250 format(a8,18x,a24)

            if (ux8(1:8) .eq. uendit(1:8)) then
                go to 440
            end if

            nsq = nsq + 1

            if (nsq .gt. nsqmax) then
                write (nttyo,1260) nsqmax
1260 format(/' * Error - (XCON3/rd3w6) Have too many basis',/7x,'species present. The code is only dimensioned',/7x,'for ',i3,' basis species. Reduce the number of such',/7x,'species or increase the dimensioning parameter nsqpar.')

                go to 990
            end if

            uspecb(nsq) = ux24

            read (ninpts,1270,err=990) ubasis(nsq)
1270 format(24x,a24)

            read (ninpts,1280,err=990) jflagb(nsq),cspb(nsq)
1280 format(10x,i2,8x,e12.5)

            if (jflagb(nsq).ge.17 .and. jflagb(nsq).le.21) then
                read (ninpts,1290,err=990) uphas1(nsq),uphas2(nsq)
1290 format(10x,a24,11x,a24)
            end if

            go to 430

440 continue

            ! Mole fractions of solid solutions.
            nxtb = 0

            if (iopt(4) .ge. 2) then
450 continue
                read (ninpts,1300,err=990) ux24
1300 format(3x,a24)

                if (ux24(1:8) .eq. uendit(1:8)) then
                    go to 470
                end if

                nxtb = nxtb + 1

                if (nxtb .gt. nxtmax) then
                    write (nttyo,1310) nxtmax
1310 format(/' * Error - (XCON3/rd3w6) Have too many solid',/7x,'solutions present. The code is only dimensioned',/7x,'for ',i3,' solid solutions. Reduce the number of such',/7x,'phases or increase the dimensioning parameter nxtpar.')

                    go to 990
                end if

                usolb(nxtb) = ux24
                iktb = 0

460 continue
                read (ninpts,1320,err=990) ux24,xx
1320 format(6x,a24,3x,f10.4)

                if (ux24(1:8) .eq. uendit(1:8)) then
                    ncompb(nxtb) = iktb
                    go to 450
                end if

                iktb = iktb + 1

                if (iktb .gt. iktmax) then
                    j2 = ilnobl(usolb(nxtb))
                    write (nttyo,1330) usolb(nxtb)(1:j2),iktmax
1330 format(/' * Error - (XCON3/rd3w6) Solid solution',/7x,'"',a,'" has too many end-members present.',/7x,'This code is only dimensioned for ',i3,' end-members',/7x,'per solid solution. Reduce the number of end-members',/7x,'or increase the dimensioning parameter iktpar.')

                    go to 990
                end if

                umemb(iktb,nxtb) = ux24
                xbarb(iktb,nxtb) = xx
                go to 460
            end if

470 continue

            go to 999

990 continue
            qrderr = .true.

999 continue
        end subroutine rd3w6