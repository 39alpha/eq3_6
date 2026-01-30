subroutine rd3ind(cgexj,cgexpi,covali,ehi,egexsi,fo2lgi,iebal3,ietmax,iodb,iopg,iopr,iopt,irdxc3,itdsf3,itermx,jgext,jetmax,jflgi,jgexti,jpres3,kxmod,mwtges,nbti,nbtmax,ncmpri,net,neti,netmax,ngexti,ninpts,ngexrt,nobswt,nodbmx,nopgmx,noprmx,noptmx,noutpt,nprob,nsbswt,ntitl,ntitmx,nttyo,nxicmx,nxmdmx,nxmod,nxti,nxtimx,pei,press,qend,qgexsh,qrderr,rho,scamas,tdspkg,tdspl,tempc,tgexp,tolbt,toldl,tolspf,ucospi,uebal,ugexj,ugexji,ugexmo,ugexp,ugexpi,ugexr,ugexsi,umemi,uobsw,uredox,usbsw,usoli,uspeci,utitl,uhfgex,uvfgex,uxkgex,uxmod,xbari,xgexsi,xhfgex,xlkgex,xvfgex,xlkmod,zgexj)
    !! This subroutine reads the EQ3NR input file in menu-style ("D")
    !! format.
    !! This subroutine is a near-clone of XCON3/rd3d8.f. The present
    !! subroutine differs from the latter, in that in addition to the
    !! pure read function, it writes an instant echo of what is read to
    !! the output file, and sandwiches this between prefatory and
    !! ending messages. To maintain close consistency with XCON3/rd3d8.f,
    !! this subroutine assigns no default values and only performs such
    !! checking of what is read to ensure that what follows is readable.
    !! The calling sequence of this subroutine is identical to that of
    !! EQ3NR/rd3inw.f, XCON3/rd3d8.f, and XCON3/rd3w8.f.
    !! This subroutine is called by:
    !!   EQ3NR/eq3nr.f
    !! Principal input:
    !!   ninpts = unit number of the stripped input file
    !!   noutpt = unit number of the output file
    !!   nttyo  = unit number of the screen file
    !! Principal output:
    !!   qrderr = flag denoting a read error
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

    integer :: ninpts
    integer :: noutpt
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
    integer :: nprob
    integer :: nsbswt
    integer :: ntitl
    integer :: nxmod
    integer :: nxti

    logical :: qend
    logical :: qgexsh
    logical :: qrderr

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

    ! Local parameter declarations.
    !   nfldpa = maximum number of fields per line
    !   nlchpa = character length of a line
    integer :: nfldpa
    integer :: nlchpa

    parameter (nfldpa = 8,nlchpa = 80)

    ! Local variable declarations.
    integer :: icount
    integer :: iei
    integer :: ii
    integer :: ival
    integer :: ivar
    integer :: j
    integer :: je
    integer :: jei
    integer :: jj
    integer :: jlast
    integer :: j2
    integer :: j3
    integer :: j4
    integer :: j5
    integer :: k1
    integer :: k2
    integer :: k3
    integer :: k4
    integer :: kxmd
    integer :: n
    integer :: nbi
    integer :: ne
    integer :: nei
    integer :: nfldmx
    integer :: nfldt
    integer :: nfldtx
    integer :: nlchmx
    integer :: nmark
    integer :: nn
    integer :: nnn
    integer :: nobsw
    integer :: nsbsw
    integer :: nxi
    integer :: nxic

    integer :: i
    integer :: jfl
    integer :: k

    integer :: ilnobl

    logical :: qgexbf
    logical :: qgexbs
    logical :: qgexef
    logical :: qgexrd
    logical :: qmark
    logical :: qnone
    logical :: qnonei
    logical :: qnonep
    logical :: qnones
    logical :: qnoner
    logical :: qokay
    logical :: qstop

    character(len=nlchpa) :: ufield(nfldpa)
    character(len=nlchpa) :: uline1
    character(len=nlchpa) :: uline2
    character(len=nlchpa) :: ulscr
    character(len=nlchpa) :: uheadr
    character(len=nlchpa) :: uheadx

    character(len=80) :: ustr
    character(len=80) :: ux80
    character(len=48) :: ux48
    character(len=24) :: ustr24
    character(len=24) :: ustrn
    character(len=1) :: ux1

    real(kind=8) :: var

    nfldmx = nfldpa
    nlchmx = nlchpa

    qrderr = .false.

    ! Check some dimensioning parameters.
    qstop = .false.

    if (ndbxpa .ne. nodbmx) then
        write (noutpt,3000) ndbxpa,nodbmx
        write (nttyo,3000) ndbxpa,nodbmx
3000 format(/' * Error - (EQ3NR/rd3ind) The dimensioning parameter',' for the',/7x,'number of iodb debugging print option switches',' with string definitions',/7x,'(ndbxpa) has a value of ',i3,', but the dimensioning',/7x,'parameter for the number',' of such switches (nodbpa) has a',/7x,'value of ',i3,'.')

        qstop = .true.
    end if

    if (npgxpa .ne. nopgmx) then
        write (noutpt,3010) npgxpa,nopgmx
        write (nttyo,3010) npgxpa,nopgmx
3010 format(/' * Error - (EQ3NR/rd3ind) The dimensioning parameter',' for the',/7x,'number of iopg activity coefficient option',' switches with string definitions',/7x,'(npgxpa) has a value',' of ' ,i3,', but the dimensioning',/7x,'parameter for the',' number of such switches (nopgpa) has a',/7x,'value of ',i3,'.')

        qstop = .true.
    end if

    if (nprxpa .ne. noprmx) then
        write (noutpt,3020) nprxpa,noprmx
        write (nttyo,3020) nprxpa,noprmx
3020 format(/' * Error - (EQ3NR/rd3ind) The dimensioning parameter',' for the',/7x,'number of iopr print option switches',' with string definitions',/7x,'(nprxpa) has a value of ',i3,', but the dimensioning',/7x,'parameter for the number',' of such switches (noprpa) has a',/7x,'value of ',i3,'.')

        qstop = .true.
    end if

    if (nptxpa .ne. noptmx) then
        write (noutpt,3030) nptxpa,noptmx
        write (nttyo,3030) nptxpa,noptmx
3030 format(/' * Error - (EQ3NR/rd3ind) The dimensioning parameter',' for the',/7x,'number of iopt model option switches',' with string definitions',/7x,'(nptxpa) has a value of ',i3,', but the dimensioning',/7x,'parameter for the number',' of such switches (noptpa) has a',/7x,'value of ',i3,'.')

        qstop = .true.
    end if

    if (qstop) then
        stop
    end if

    ! Title.
    qend = .false.

    ! Read the first line of the file. This must be a separator line
    ! ("|----- ... -----|").
    read (ninpts,1000,end=100,err=990) uline1
1000 format(a80)

    if (uline1(1:8) .ne. '|-------') then
        j2 = ilnobl(uline1)
        j2 = min(j2,70)
        write (noutpt,1010) uline1(1:j2)
        write (nttyo,1010) uline1(1:j2)
1010 format(/' * Error - (EQ3NR/rd3ind) The first line of a "D"',' format input file',/7x,'should be a separator line;',' therefore, it should begin with',/7x,'"|-------". The first',' line begins instead with',/7x,'"',a,'".')

        go to 990
    end if

    go to 105

100 continue
    qend = .true.
    go to 999

105 continue
    write (noutpt,1012) nprob
    write (nttyo,1012) nprob
1012 format(//' Reading problem ',i3,' from the input file ...',/)

    write (noutpt,1014) uline1
1014 format(a80)

    ! Read the block title ("Title") from a two-line header.
    uheadx = 'Title'
    nfldtx = 2
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    write (noutpt,1016) uline1,uline2
1016 format(a80,/a80)

    ! Now read the title itself.
    n = 0

    do nn = 1,ntitmx + 1
        read (ninpts,1000,err=990) uline1
        write (noutpt,1014) uline1
        call parsln(nfldmx,nfldt,nlchmx,ufield,uline1,ulscr)
        ustr = ufield(1)

        ! A separator line terminates the this block. It is not part
        ! of the title itself.
        if (ustr(1:8) .eq. '--------') then
            go to 120
        end if

        n = n + 1

        if (n .gt. ntitmx) then
            write (noutpt,1015) ntitmx
            write (nttyo,1015) ntitmx
1015 format(/' * Error - (EQ3NR/rd3ind) Have too many lines in',/7x,'the title. The code is only dimensioned for ',i4,/7x,'lines. Reduce the size of the title or increase the',/7x,'dimensioning parameter ntitpa.')

            go to 990
        end if

        utitl(n) = ufield(1)
    end do

120 continue
    ntitl = n

    if (ntitl .le. 0) then
        write (noutpt,1032)
        write (nttyo,1032)
1032 format(3x,'There is no input problem title.')
    end if

    if (ntitl .gt. 0) then
        ! Write the first 5 lines of the input problem title to the
        ! screen file.
        write (nttyo,1040)
1040 format(3x,'The input problem title is (first 5 lines',' maximum):',/)

        i = min(ntitl,5)

        do n = 1,i
            j2 = ilnobl(utitl(n))
            j2 = min(j2,74)
            write (nttyo,1052) utitl(n)(1:j2)
1052 format(5x,a)
        end do
    end if

    write (nttyo,1062)
1062 format(/3x,'Continuing to read the problem input ...')

    ! Special basis switches.
    ! Read the block title from a two-line header.
    uheadx = 'Special Basis Switches'
    nfldtx = 2
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    write (noutpt,1016) uline1,uline2

    nsbsw = 0

    do n = 1,nbtmax + 1
        ! Read a line. If the block has not been completely read,
        ! this contains the name of a species to "Replace", and a
        ! sub-block for that species follows. Otherwise, this line is
        ! the first line of the next block.
        nfldtx = 0
        call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

        if (qrderr) then
            go to 999
        end if

        ustr = ufield(1)
        uheadx = 'Replace'
        call locase(ustr)
        call locase(uheadx)
        j2 = ilnobl(ustr)
        j3 = ilnobl(uheadx)

        if (ustr(1:j2) .ne. uheadx(1:j3)) then
            ! Back up.
            backspace ninpts
            go to 125
        end if

        write (noutpt,1014) uline1
        ustr = ufield(2)
        ustrn = ustr(1:24)
        call locase(ustrn)

        if (ustrn(1:5).eq.'none ' .or. ustr(1:1).eq.' ') then
            ustr = 'None'
        end if

        qnone = ustr(1:5).eq.'None '

        if (.not.qnone) then
            nsbsw = nsbsw + 1
            usbsw(1,nsbsw) = ufield(2)(1:48)
        end if

        ! Read the name of the "with" species from a two-line header.
        uheadx = 'with'
        nfldtx = 3
        call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

        if (qrderr) then
            go to 999
        end if

        write (noutpt,1016) uline1,uline2

        if (.not.qnone) then
            usbsw(2,nsbsw) = ufield(2)(1:48)
        end if
    end do

125 continue
    nsbswt = nsbsw

    ! Temperature.
    tempc = 0.

    ! Read the data from a two-line header.
    uheadx = 'Temperature (C)'
    nfldtx = 3
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    write (noutpt,1016) uline1,uline2
    ustr = ufield(2)
    call chreal(nttyo,qrderr,ustr,var)

    if (qrderr) then
        go to 999
    end if

    tempc = var

    ! Pressure.
    jpres3 = 0
    press = 0.
    icount = 0

    ! Read a one-line header.
    uheadx = 'Pressure option (jpres3):'
    nfldtx = 1
    call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)

    if (qrderr) then
        go to 999
    end if

    write (noutpt,1014) uline1

    ! Read the first option from a data line.
    nfldtx = 1
    call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

    if (qrderr) then
        go to 999
    end if

    write (noutpt,1014) uline1
    ustr = ufield(1)
    ustr24 = ustr(5:29)
    uheadx = '( 0) Data file reference'
    call locase(ustr24)
    call locase(uheadx)
    j2 = ilnobl(ustr24)
    j3 = ilnobl(uheadx)

    if (ustr24(1:j2) .ne. uheadx(1:j3)) then
        write (noutpt,1020) uheadx(1:j3),ustr24(1:j2)
        write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
1020 format(/' * Error - (EQ3NR/rd3ind) Was expecting to find a',' line beginning with',/7x,'"',a,'", instead found one',' beginning with',/7x,'"',a,'".')

        qrderr = .true.
        go to 999
    end if

    if (ustr(1:1).ne.'[' .or. ustr(3:3).ne.']') then
        j2 = ilnobl(ustr)
        j2 = min(j2,70)
        write (noutpt,1022) ustr(1:j2)
        write (nttyo,1022) ustr(1:j2)
1022 format(/' * Error - (EQ3NR/rd3ind) Was expecting to find an',/7x,'option check box "[ ]" in the line beginning with',/7x,'"',a,'".')

        qrderr = .true.
        go to 999
    end if

    ux1 = ustr(2:2)

    if (ux1.eq.'*' .or. ux1.eq.'x' .or. ux1.eq.'X') then
        jpres3 = 0
        press = 0.
        icount = icount + 1
    end if

    ! Read the second option from a data line.
    nfldtx = 1
    call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

    if (qrderr) then
        go to 999
    end if

    write (noutpt,1014) uline1
    ustr = ufield(1)
    ustr24 = ustr(5:29)
    uheadx = '( 1) 1.013-bar/steam-sat'
    call locase(ustr24)
    call locase(uheadx)
    j2 = ilnobl(ustr24)
    j3 = ilnobl(uheadx)

    if (ustr24(1:j2) .ne. uheadx(1:j3)) then
        write (noutpt,1020) uheadx(1:j3),ustr24(1:j2)
        write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
        qrderr = .true.
        go to 999
    end if

    if (ustr(1:1).ne.'[' .or. ustr(3:3).ne.']') then
        j2 = ilnobl(ustr)
        j2 = min(j2,70)
        write (noutpt,1022) ustr(1:j2)
        write (nttyo,1022) ustr(1:j2)
        qrderr = .true.
        go to 999
    end if

    ux1 = ustr(2:2)

    if (ux1.eq.'*' .or. ux1.eq.'x' .or. ux1.eq.'X') then
        jpres3 = 1
        press = 0.
        icount = icount + 1
    end if

    ! Read the third (last) option from a two-line combination
    ! (a data line plus a separator line).
    nfldtx = 3
    call rdd2l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    write (noutpt,1016) uline1,uline2
    ustr = ufield(1)
    ustr24 = ustr(5:29)
    uheadx = '( 2) Value (bars)'
    call locase(ustr24)
    call locase(uheadx)
    j2 = ilnobl(ustr24)
    j3 = ilnobl(uheadx)

    if (ustr24(1:j2) .ne. uheadx(1:j3)) then
        write (noutpt,1020) uheadx(1:j3),ustr24(1:j2)
        write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
        go to 999
    end if

    if (ustr(1:1).ne.'[' .or. ustr(3:3).ne.']') then
        j2 = ilnobl(ustr)
        j2 = min(j2,70)
        write (noutpt,1022) ustr(1:j2)
        write (nttyo,1022) ustr(1:j2)
        qrderr = .true.
        go to 999
    end if

    ux1 = ustr(2:2)

    if (ux1.eq.'*' .or. ux1.eq.'x' .or. ux1.eq.'X') then
        jpres3 = 2
        ustr = ufield(2)
        call chreal(nttyo,qrderr,ustr,var)

        if (qrderr) then
            go to 999
        end if

        press = var
        icount = icount + 1
    end if

    if (icount .eq. 0) then
        write (noutpt,1024)
        write (nttyo,1024)
1024 format(/' * Note - (EQ3NR/rd3ind) No option was selected for',' the pressure.',/7x,'The pressure will be set to be in',' accord with the',/7x,'data file reference pressure curve.')

        jpres3 = 0
        press = 0.
    else if (icount .gt. 1) then
        write (noutpt,1026)
        write (nttyo,1026)
1026 format(/' * Warning - (EQ3NR/rd3ind) Multiple options were',' selected for',/7x,'the pressure. The pressure will be set',' to be in accord with the',/7x,'data file reference pressure',' curve.')

        jpres3 = 0
        press = 0.
    end if

    ! Density.
    rho = 0.

    ! Read the data from a two-line header.
    uheadx = 'Density'
    nfldtx = 3
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    write (noutpt,1016) uline1,uline2
    ustr = ufield(2)
    call chreal(nttyo,qrderr,ustr,var)

    if (qrderr) then
        go to 999
    end if

    rho = var

    ! Total dissolved solutes.
    itdsf3 = 0
    tdspkg = 0.
    tdspl = 0.
    icount = 0

    ! Read the data from a two-line header.
    uheadx = 'Total dissolved solutes option (itdsf3):'
    nfldtx = 1
    call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)

    if (qrderr) then
        go to 999
    end if

    write (noutpt,1014) uline1

    ! Read the first option from a data line.
    nfldtx = 3
    call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

    if (qrderr) then
        go to 999
    end if

    write (noutpt,1014) uline1
    ustr = ufield(1)
    ustr24 = ustr(5:29)
    uheadx = '( 0) Value (mg/kg.sol)'
    call locase(ustr24)
    call locase(uheadx)
    j2 = ilnobl(ustr24)
    j3 = ilnobl(uheadx)

    if (ustr24(1:j2) .ne. uheadx(1:j3)) then
        write (noutpt,1020) uheadx(1:j3),ustr24(1:j2)
        write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
        qrderr = .true.
        go to 999
    end if

    if (ustr(1:1).ne.'[' .or. ustr(3:3).ne.']') then
        j2 = ilnobl(ustr)
        j2 = min(j2,70)
        write (noutpt,1022) ustr(1:j2)
        write (nttyo,1022) ustr(1:j2)
        qrderr = .true.
        go to 999
    end if

    ux1 = ustr(2:2)

    if (ux1.eq.'*' .or. ux1.eq.'x' .or. ux1.eq.'X') then
        itdsf3 = 0
        ustr = ufield(2)
        call chreal(nttyo,qrderr,ustr,var)

        if (qrderr) then
            go to 999
        end if

        tdspkg = var
        icount = icount + 1
    end if

    ! Read the second (last) option from a two-line combination
    ! (a data line plus a separator line).
    nfldtx = 3
    call rdd2l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    write (noutpt,1016) uline1,uline2
    ustr = ufield(1)
    ustr24 = ustr(5:29)
    uheadx = '( 1) Value (mg/L)'
    call locase(ustr24)
    call locase(uheadx)
    j2 = ilnobl(ustr24)
    j3 = ilnobl(uheadx)

    if (ustr24(1:j2) .ne. uheadx(1:j3)) then
        write (noutpt,1020) uheadx(1:j3),ustr24(1:j2)
        write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
        go to 999
    end if

    if (ustr(1:1).ne.'[' .or. ustr(3:3).ne.']') then
        j2 = ilnobl(ustr)
        j2 = min(j2,70)
        write (noutpt,1022) ustr(1:j2)
        write (nttyo,1022) ustr(1:j2)
        qrderr = .true.
        go to 999
    end if

    ux1 = ustr(2:2)

    if (ux1.eq.'*' .or. ux1.eq.'x' .or. ux1.eq.'X') then
        itdsf3 = 1
        ustr = ufield(2)
        call chreal(nttyo,qrderr,ustr,var)

        if (qrderr) then
            go to 999
        end if

        tdspl = var
        icount = icount + 1
    end if

    if (icount .eq. 0) then
        write (noutpt,1034)
        write (nttyo,1034)
1034 format(/' * Note - (EQ3NR/rd3ind) No option was selected',' for the total',/7x,'dissolved solutes (TDS). The TDS will',' be set to 0 mg/kg.sol.')

        itdsf3 = 0
        tdspl = 0.
        tdspkg = 0.
    else if (icount .gt. 1) then
        write (noutpt,1036) tdspkg
        write (nttyo,1036) tdspkg
1036 format(/' * Warning - (EQ3NR/rd3ind) Multiple options were',' selected for the',/7x,'total dissolved solutes (TDS).',' The TDS will be set to',/7x,g12.5,' mg/kg.sol.')

        itdsf3 = 0
    end if

    ! Electrical balancing species.
    iebal3 = 0
    uebal = ' '
    icount = 0

    ! Read a one-line header.
    uheadx = 'Electrical balancing option (iebal3):'
    nfldtx = 1
    call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)

    if (qrderr) then
        go to 999
    end if

    write (noutpt,1014) uline1

    ! Read the first option from a data line.
    nfldtx = 1
    call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

    if (qrderr) then
        go to 999
    end if

    write (noutpt,1014) uline1
    ustr = ufield(1)
    ustr24 = ustr(5:29)
    uheadx = '( 0) No balancing is don'
    call locase(ustr24)
    call locase(uheadx)
    j2 = ilnobl(ustr24)
    j3 = ilnobl(uheadx)

    if (ustr24(1:j2) .ne. uheadx(1:j3)) then
        write (noutpt,1020) uheadx(1:j3),ustr24(1:j2)
        write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
        qrderr = .true.
        go to 999
    end if

    if (ustr(1:1).ne.'[' .or. ustr(3:3).ne.']') then
        j2 = ilnobl(ustr)
        j2 = min(j2,70)
        write (noutpt,1022) ustr(1:j2)
        write (nttyo,1022) ustr(1:j2)
        qrderr = .true.
        go to 999
    end if

    ux1 = ustr(2:2)

    if (ux1.eq.'*' .or. ux1.eq.'x' .or. ux1.eq.'X') then
        iebal3 = 0
        uebal = 'None'
        icount = icount + 1
    end if

    ! Read the second (last) option from a two-line combination
    ! (a data line plus a separator line).
    nfldtx = 3
    call rdd2l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    write (noutpt,1016) uline1,uline2
    ustr = ufield(1)
    ustr24 = ustr(5:29)
    uheadx = '( 1) Balance on species'
    call locase(ustr24)
    call locase(uheadx)
    j2 = ilnobl(ustr24)
    j3 = ilnobl(uheadx)

    if (ustr24(1:j2) .ne. uheadx(1:j3)) then
        write (noutpt,1020) uheadx(1:j3),ustr24(1:j2)
        write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
        qrderr = .true.
        go to 999
    end if

    if (ustr(1:1).ne.'[' .or. ustr(3:3).ne.']') then
        j2 = ilnobl(ustr)
        j2 = min(j2,70)
        write (noutpt,1022) ustr(1:j2)
        write (nttyo,1022) ustr(1:j2)
        qrderr = .true.
        go to 999
    end if

    ux1 = ustr(2:2)

    if (ux1.eq.'*' .or. ux1.eq.'x' .or. ux1.eq.'X') then
        iebal3 = 1
        ustr = ufield(2)
        uebal = ustr(1:24)
        icount = icount + 1
    end if

    if (icount .eq. 0) then
        write (noutpt,1044)
        write (nttyo,1044)
1044 format(/' * Note - (EQ3NR/rd3ind) No electrical balancing',' option was',/7x,'selected. No balancing will be done.')

        iebal3 = 0
        uebal = 'None'
    else if (icount .gt. 1) then
        write (noutpt,1046)
        write (nttyo,1046)
1046 format(/' * Warning - (EQ3NR/rd3ind) Multiple electrical',' balancing options.',/7x,'were selected. No balancing',' will be done.')

        iebal3 = 0
        uebal = 'None'
    end if

    ! Default redox constraint.
    irdxc3 = 0
    fo2lgi = 0.
    ehi = 0.
    pei = 0.
    uredox = ' '
    icount = 0

    ! Read a one-line header.
    uheadx = 'Default redox constraint (irdxc3):'
    nfldtx = 1
    call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)

    if (qrderr) then
        go to 999
    end if

    write (noutpt,1014) uline1

    ! Read the first option from a data line.
    nfldtx = 1
    call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

    if (qrderr) then
        go to 999
    end if

    write (noutpt,1014) uline1
    ustr = ufield(1)
    ustr24 = ustr(5:29)
    uheadx = '(-3) Use O2(g) line in t'
    call locase(ustr24)
    call locase(uheadx)
    j2 = ilnobl(ustr24)
    j3 = ilnobl(uheadx)

    if (ustr24(1:j2) .ne. uheadx(1:j3)) then
        write (noutpt,1020) uheadx(1:j3),ustr24(1:j2)
        write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
        qrderr = .true.
        go to 999
    end if

    if (ustr(1:1).ne.'[' .or. ustr(3:3).ne.']') then
        j2 = ilnobl(ustr)
        j2 = min(j2,70)
        write (noutpt,1022) ustr(1:j2)
        write (nttyo,1022) ustr(1:j2)
        qrderr = .true.
        go to 999
    end if

    ux1 = ustr(2:2)

    if (ux1.eq.'*' .or. ux1.eq.'x' .or. ux1.eq.'X') then
        irdxc3 = -3
        icount = icount + 1
    end if

    ! Read the second option from a data line.
    nfldtx = 3
    call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

    if (qrderr) then
        go to 999
    end if

    write (noutpt,1014) uline1
    ustr = ufield(1)
    ustr24 = ustr(5:29)
    uheadx = '(-2) pe (pe units)'
    call locase(ustr24)
    call locase(uheadx)
    j2 = ilnobl(ustr24)
    j3 = ilnobl(uheadx)

    if (ustr24(1:j2) .ne. uheadx(1:j3)) then
        write (noutpt,1020) uheadx(1:j3),ustr24(1:j2)
        write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
        qrderr = .true.
        go to 999
    end if

    if (ustr(1:1).ne.'[' .or. ustr(3:3).ne.']') then
        j2 = ilnobl(ustr)
        j2 = min(j2,70)
        write (noutpt,1022) ustr(1:j2)
        write (nttyo,1022) ustr(1:j2)
        qrderr = .true.
        go to 999
    end if

    ux1 = ustr(2:2)

    if (ux1.eq.'*' .or. ux1.eq.'x' .or. ux1.eq.'X') then
        irdxc3 = -2
        ustr = ufield(2)
        call chreal(nttyo,qrderr,ustr,var)

        if (qrderr) then
            go to 999
        end if

        pei = var
        icount = icount + 1
    end if

    ! Read the third option from a data line.
    nfldtx = 3
    call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

    if (qrderr) then
        go to 999
    end if

    write (noutpt,1014) uline1
    ustr = ufield(1)
    ustr24 = ustr(5:29)
    uheadx = '(-1) Eh (volts)'
    call locase(ustr24)
    call locase(uheadx)
    j2 = ilnobl(ustr24)
    j3 = ilnobl(uheadx)

    if (ustr24(1:j2) .ne. uheadx(1:j3)) then
        write (noutpt,1020) uheadx(1:j3),ustr24(1:j2)
        write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
        qrderr = .true.
        go to 999
    end if

    if (ustr(1:1).ne.'[' .or. ustr(3:3).ne.']') then
        j2 = ilnobl(ustr)
        j2 = min(j2,70)
        write (noutpt,1022) ustr(1:j2)
        write (nttyo,1022) ustr(1:j2)
        qrderr = .true.
        go to 999
    end if

    ux1 = ustr(2:2)

    if (ux1.eq.'*' .or. ux1.eq.'x' .or. ux1.eq.'X') then
        irdxc3 = -1
        ustr = ufield(2)
        call chreal(nttyo,qrderr,ustr,var)

        if (qrderr) then
            go to 999
        end if

        ehi = var
        icount = icount + 1
    end if

    ! Read the fourth option from a data line.
    nfldtx = 3
    call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

    if (qrderr) then
        go to 999
    end if

    write (noutpt,1014) uline1
    ustr = ufield(1)
    ustr24 = ustr(5:29)
    uheadx = '( 0) Log fO2 (log bars)'
    call locase(ustr24)
    call locase(uheadx)
    j2 = ilnobl(ustr24)
    j3 = ilnobl(uheadx)

    if (ustr24(1:j2) .ne. uheadx(1:j3)) then
        write (noutpt,1020) uheadx(1:j3),ustr24(1:j2)
        write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
        qrderr = .true.
        go to 999
    end if

    if (ustr(1:1).ne.'[' .or. ustr(3:3).ne.']') then
        j2 = ilnobl(ustr)
        j2 = min(j2,70)
        write (noutpt,1022) ustr(1:j2)
        write (nttyo,1022) ustr(1:j2)
        qrderr = .true.
        go to 999
    end if

    ux1 = ustr(2:2)

    if (ux1.eq.'*' .or. ux1.eq.'x' .or. ux1.eq.'X') then
        irdxc3 = 0
        ustr = ufield(2)
        call chreal(nttyo,qrderr,ustr,var)

        if (qrderr) then
            go to 999
        end if

        fo2lgi = var
        icount = icount + 1
    end if

    ! Read the fifth (last) option from a two-line combination
    ! (a data line plus a separator line).
    nfldtx = 3
    call rdd2l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    write (noutpt,1016) uline1,uline2
    ustr = ufield(1)
    ustr24 = ustr(5:29)
    uheadx = '( 1) Couple (aux. sp.)'
    call locase(ustr24)
    call locase(uheadx)
    j2 = ilnobl(ustr24)
    j3 = ilnobl(uheadx)

    if (ustr24(1:j2) .ne. uheadx(1:j3)) then
        write (noutpt,1020) uheadx(1:j3),ustr24(1:j2)
        write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
        qrderr = .true.
        go to 999
    end if

    if (ustr(1:1).ne.'[' .or. ustr(3:3).ne.']') then
        j2 = ilnobl(ustr)
        j2 = min(j2,70)
        write (noutpt,1022) ustr(1:j2)
        write (nttyo,1022) ustr(1:j2)
        qrderr = .true.
        go to 999
    end if

    ux1 = ustr(2:2)

    if (ux1.eq.'*' .or. ux1.eq.'x' .or. ux1.eq.'X') then
        irdxc3 = 1
        ustr = ufield(2)
        uredox = ustr(1:24)
        icount = icount + 1
    end if

    if (icount .eq. 0) then
        write (noutpt,1047) fo2lgi
        write (nttyo,1047) fo2lgi
1047 format(/' * Note - (EQ3NR/rd3ind) No default redox constraint',' option was',/7x,'selected. A default condition of log fO2= ',g12.5,' (log bars) will be used.')

        irdxc3 = 0
    else if (icount .gt. 1) then
        write (noutpt,1048)
        write (nttyo,1048)
1048 format(/' * Warning - (EQ3NR/rd3ind) Multiple default redox',' constraint options.',/7x,'were selected. The last option',' that was selected will be used.')
    end if

    ! Aqueous basis species and associated constraints.
    ! Read the first part of the block title from a one-line header.
    uheadx = 'Aqueous Basis Species/Constraint Species'
    nfldtx = 3
    call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)

    if (qrderr) then
        go to 999
    end if

    write (noutpt,1014) uline1

    ! Read the second part of the block title from a two-line header.
    uheadx = '(uspeci(n)/ucospi(n))'
    nfldtx = 3
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    write (noutpt,1016) uline1,uline2

    nbi = 0

    ! Read the first data line.
    nfldtx = 0
    call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

    if (qrderr) then
        go to 999
    end if

    write (noutpt,1014) uline1

    ! The label below marks a return point for processing subsequent
    ! lines in the current data block.
130 continue
    ustr = ufield(1)

    ! There are no data remaining in the current block if a
    ! separator line has been encountered.
    if (ustr(1:8) .eq. '--------') then
        go to 180
    end if

    nfldtx = 3

    if (nfldt .ne. nfldtx) then
        j2 = ilnobl(uline1)
        j2 = min(j2,70)
        write (noutpt,1030) nfldt,nfldtx,uline1(1:j2)
        write (nttyo,1030) nfldt,nfldtx,uline1(1:j2)
1030 format(/' * Warning - (EQ3NR/rd3ind) Found ',i2,' fields',' where ',i2,/7x,'were expected on the line which begins with',/7x,'"',a,'".')
    end if

    ux48 = ufield(1)(1:48)

    nbi = nbi + 1

    if (nbi .gt. nbtmax) then
        j2 = ilnobl(ux48)
        write (noutpt,1050) nbtmax,ux48(1:j2)
        write (nttyo,1050) nbtmax,ux48(1:j2)
1050 format(/' * Error - (EQ3NR/rd3ind) The number of basis',' species read',/7x,'from the input file exceeded the',' maximum of ',i3,' while',/7x,'trying to read data for',' the species ',a,'.',/7x,'Increase the dimensioning',' parameter nbtpar.')

        go to 990
    end if

    uspeci(nbi) = ux48

    ustr = ufield(2)
    call chreal(nttyo,qrderr,ustr,var)

    if (qrderr) then
        go to 999
    end if

    covali(nbi) = var

    ustr = ufield(3)
    call locase(ustr)

    do n = -1,njfxpa
        uheadx = ujf3(n)
        call locase(uheadx)

        if (ustr(1:16) .eq. uheadx(1:16)) then
            jflgi(nbi) = n
            go to 150
        end if
    end do

    j2 = ilnobl(ustr)
    write (noutpt,1060) ustr(1:j2)
    write (nttyo,1060) ustr(1:j2)
1060 format(/" * Error - (EQ3NR/rd3ind) Don't recognize the",' jflag option string',/7x,' "',a,'". This should',' be one of the strings defined in the',/7x,'ujf3 array.')

    go to 990

150 continue
    if (jflgi(nbi).eq.17 .or. jflgi(nbi).eq.18 .or. jflgi(nbi).eq.25) then
        ! Have an option that requires a second line of data to
        ! complete the constraint.
        nfldtx = 3
        call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

        if (qrderr) then
            go to 999
        end if

        write (noutpt,1014) uline1
        ustr24 = ufield(1)(1:24)
        uheadx = '->'
        j2 = ilnobl(ustr24)
        j3 = ilnobl(uheadx)

        if (ustr24(1:j2) .ne. uheadx(1:j3)) then
            write (noutpt,1020) uheadx(1:j3),ustr24(1:j2)
            write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
            qrderr = .true.
            go to 999
        end if

        ustr = ufield(2)
        ucospi(nbi)(1:48) = ustr(1:48)
    end if

    ! Read the next data line. Go back to process it.
    nfldtx = 0
    call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

    if (qrderr) then
        go to 999
    end if

    write (noutpt,1014) uline1
    go to 130

180 continue

    nbti = nbi

    ! Add comment lines listing valid option strings.
    write (noutpt,1162)
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
                write (noutpt,1164) ux80(1:80)
1164 format(a)

                ux80(2:79) = ' '
                j3 = 6
                k = 0
            end if
        end if
    end do

    if (k .gt. 0) then
        write (noutpt,1164) ux80(1:80)
    end if

    write (noutpt,2052)
2052 format('*',78('-'),'*')

    ! Ion exchanger creation.
    net = 0
    qnonep = .false.

    ! Read a two-line header for the block.
    uheadx = 'Create Ion Exchangers'
    nfldtx = 2
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    write (noutpt,1016) uline1,uline2

    ! Read the advisory from one line: one or more exchanger creation
    ! blocks follow or not. This advisory will be confirmed by
    ! examining the input that actually follows. Here qgexbf is
    ! the advisory flag. It is .true. if the advisory indicates
    ! that one or more exchanger creation blocks follow. Here
    ! also qgexbs is .true. if qgexbf has been set.
    qgexbs = .false.
    nfldtx = 1
    call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

    if (qrderr) then
        go to 999
    end if

    write (noutpt,1014) uline1

    ustr = ufield(1)
    call lejust(ustr)
    uheadx = 'Advisory:'
    call locase(ustr)
    call locase(uheadx)

    if (ustr(1:9) .ne. uheadx(1:9)) then
        qrderr = .true.
        j3 = ilnobl(uline1(1:70))
        write (nttyo,1210) uheadx(1:9),uline1(1:j3)
        write (noutpt,1210) uheadx(1:9),uline1(1:j3)
1210 format(/' * Error - (EQ3NR/rd3ind) Was expecting to find the',/7x,'string "',a,'" at the start of the line beginning with',/7x,'"',a,'".')

        go to 999
    end if

    uheadx = 'Advisory: at least one exchanger creation block follows on this file.'
    call locase(ustr)
    call locase(uheadx)
    j2 = ilnobl(ustr)
    j3 = ilnobl(uheadx)

    if (ustr(1:j2) .eq. uheadx(1:j3)) then
        qgexbf = .true.
        qgexbs = .true.
    else
        uheadx =  'Advisory: no exchanger creation blocks follow on this file.'
        call locase(uheadx)
        j3 = ilnobl(uheadx)

        if (ustr(1:j2) .eq. uheadx(1:j3)) then
            qgexbf = .false.
            qgexbs = .true.
        else
            qgexbf = .false.
            qgexbs = .false.
            j3 = ilnobl(uline1(1:70))
            write (nttyo,1220) uline1(1:j3)
            write (noutpt,1220) uline1(1:j3)
1220 format(/' * Warning - (EQ3NR/rd3ind) Could not interpret the',' advisory on whether',/7x,'or not one or more exchanger',' blocks follow. The problem is with',/7x,'the line',' beginning with',/7x,'"',a,'".',/7x,'The presence of such',' blocks will be determined directly.')
        end if
    end if

    ! Read the first line of the qgexsh option (on processing, show
    ! at least one exchanger creation block on a "D" format input or
    ! pickup file).
    uheadx = 'Option: on further processing (writing a pickup file or running XCON3 on the'
    nfldtx = 1
    call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)

    if (qrderr) then
        go to 999
    end if

    write (noutpt,1014) uline1

    ! Read the second line of the qgexsh option.
    uheadx = 'present file), force the inclusion of at least one such block (qgexsh):'
    nfldtx = 1
    call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)

    if (qrderr) then
        go to 999
    end if

    write (noutpt,1014) uline1

    ! Read the third and last line of the qgexsh option, plus the
    ! following separator line.
    nfldtx = 1
    call rdd2l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    write (noutpt,1016) uline1,uline2

    k1 = index(uline1,'[')
    k2 = index(uline1,']')
    j3 = ilnobl(uline1)

    if (k1.le.0 .or. k2.le.0 .or. k2.lt.(k1 + 2)) then
        write (nttyo,1230) uline1(1:j3)
        write (noutpt,1230) uline1(1:j3)
1230 format(/' * Error - (EQ3NR/rd3ind) Was expecting to find an',/7x,'option check box "[ ]" in the line beginning with',/7x,'"',a,'"',/7x,'The expected line contains the box for',' the qgexsh flag.')

        go to 990
    end if

    ! Check the full option choice string.
    ustr = '(.true.)'

    if (index(uline1(1:j3),ustr(1:8)) .le. 0) then
        write (nttyo,1240) uline1(1:j2),ustr(1:8)
        write (noutpt,1240) uline1(1:j2),ustr(1:8)
1240 format(/' * Error - (EQ3NR/rd3ind) The qgexsh flag line',/7x,'"',a,'"',/7x,"read from the input file doesn't contain",' the expected string',/7x,'"',a,'".')

        go to 990
    end if

    ustr24 = uline1(k1 + 1:k2 -1)
    call lejust(ustr24)
    j2 = ilnobl(ustr24)
    qmark = .false.

    if (j2 .gt. 0) then
        if (index(ustr24,'*') .ge. 1) then
            qmark = .true.
        else if (index(ustr24,'x') .ge. 1) then
            qmark = .true.
        else if (index(ustr24,'X') .ge. 1) then
            qmark = .true.
        else
            write (nttyo,1250) ustr24(1:j2),uline1(1:j3)
1250 format(/" * Error - (EQ3NR/rd3ind) Don't recognize the",' string "',a,'"',/7x,'that appears on the qgexsh',' option switch choice line',/7x,'"',a,'"',/7x,'read from the input file. An option choice should',' be chosen by',/7x,'placing a "*", "x", or "X" in the',' checkbox ("[ ]").')

            go to 990
        end if
    end if

    qgexsh = qmark

    ! Check to see if an exchanger block actually follows. This is a
    ! test of the qgexbf flag. Read a line. If it contains the string
    ! 'Exchanger phase' in the first field, an exchanger block is
    ! present.
    nfldtx = 0
    call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

    if (qrderr) then
        go to 999
    end if

    ustr = ufield(1)
    uheadx = 'Exchanger phase'
    call locase(ustr)
    call locase(uheadx)
    j2 = ilnobl(ustr)
    j3 = ilnobl(uheadx)

    if (ustr(1:j2) .eq. uheadx(1:j3)) then
        ! An exchanger block is present.
        qgexrd = .true.

        if (qgexbs .and. .not.qgexbf) then
            write (nttyo,1260)
            write (noutpt,1260)
1260 format(/' * Note - (EQ3NR/rd3ind) The advisory on the',' presence of',/7x,'ion exchanger blocks was negative,',' but one or more blocks is',/7x,'actually present.')
        end if
    else
        ! An exchanger block is not present.
        qgexrd = .false.

        if (qgexbs .and. qgexbf) then
            write (nttyo,1270)
            write (noutpt,1270)
1270 format(/' * Note - (EQ3NR/rd3ind) The advisory on the',' presence of',/7x,'ion exchanger blocks was positive,',' but no blocks are actually',/7x,'present.')
        end if
    end if

    ! Back up.
    backspace ninpts

    if (.not.qgexrd) then
        go to 575
    end if

    ! Loop on exchanger phases.
    ne = 0

    do nn = 1,netmax + 1
        ! Read a line. If the block has not been completely read,
        ! this contains the name of an exchanger phase, and a sub-block
        ! for that phase follows. Otherwise, this line is the first line
        ! of the next block.
        nfldtx = 0
        call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

        if (qrderr) then
            go to 999
        end if

        ustr = ufield(1)
        uheadx = 'Exchanger phase'
        call locase(ustr)
        call locase(uheadx)
        j2 = ilnobl(ustr)
        j3 = ilnobl(uheadx)

        if (ustr(1:j2) .ne. uheadx(1:j3)) then
            ! Back up.
            backspace ninpts
            go to 240
        end if

        write (noutpt,1014) uline1
        ustr = ufield(2)
        ustrn = ustr(1:24)
        call locase(ustrn)

        if (ustrn(1:5).eq.'none ' .or. ustr(1:1).eq.' ') then
            ustr = 'None'
        end if

        qnonep = ustr(1:5).eq.'None '

        if (.not.qnonep) then
            net = net + 1
            ne = ne + 1

            if (ne .gt. netmax) then
                write (noutpt,1620) netmax,ne
                write (nttyo,1620) netmax,ne
1620 format(/' * Error - (EQ3NR/rd3ind) Have exceeded the',' maximum number of ',i3,/7x,'generic ion exchange phases',' while reading the data to create',/7x,'such phases.',' Increase the dimensioning parameter netpar',/7x,'to at',' least ',i3,'.')

                go to 990
            end if

            ugexp(ne) = ustr(1:24)
            jgext(ne) = 0
        end if

        ! Read the separator line following the data line containing
        ! the name of an exchanger phase.
        nfldtx = 1
        call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

        if (qrderr) then
            go to 999
        end if

        write (noutpt,1014) uline1
        ustr24 = ufield(1)(1:24)
        uheadx = '--------'

        if (ustr24(1:8) .ne. uheadx(1:8)) then
            j2 = ilnobl(uline1)
            j2 = min(j2,70)
            write (noutpt,1070) uline1(1:j2)
            write (nttyo,1070) uline1(1:j2)
1070 format(/' * Error - (EQ3NR/rd3ind) Read the following line',' where a separator line',/7x,'("|------- ... -------|")',' was expected:',/7x,'"',a,'".')

            qrderr = .true.
            go to 999
        end if

        ! Read the molecular weight of the bare exchange ligand (Z)
        ! from a two-line header.
        uheadx = '->'
        nfldtx = 4
        call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

        if (qrderr) then
            go to 999
        end if

        write (noutpt,1016) uline1,uline2
        ustr24 = ufield(2)(1:24)
        uheadx = 'Mol. Wt. (Z)'
        call locase(ustr24)
        call locase(uheadx)
        j2 = ilnobl(ustr24)
        j3 = ilnobl(uheadx)

        if (ustr24(1:j2) .ne. uheadx(1:j3)) then
            write (noutpt,1020) uheadx(1:j3),ustr24(1:j2)
            write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
            qrderr = .true.
            go to 999
        end if

        ustr = ufield(3)
        call chreal(nttyo,qrderr,ustr,var)

        if (qrderr) then
            go to 999
        end if

        if (.not.qnonep) then
            mwtges(ne) = var
        end if

        ! Read the exchange model name from a two-line header.
        uheadx = '->'
        nfldtx = 4
        call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

        if (qrderr) then
            go to 999
        end if

        write (noutpt,1016) uline1,uline2
        ustr24 = ufield(2)(1:24)
        uheadx = 'Exchange model'
        call locase(ustr24)
        call locase(uheadx)
        j2 = ilnobl(ustr24)
        j3 = ilnobl(uheadx)

        if (ustr24(1:j2) .ne. uheadx(1:j3)) then
            write (noutpt,1020) uheadx(1:j3),ustr24(1:j2)
            write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
            qrderr = .true.
            go to 999
        end if

        if (.not.qnonep) then
            ugexmo(ne) = ufield(3)(1:24)
        end if

        ! Read the reference temperature (C) for the thermodynamic
        ! data from a two-line header.
        uheadx = '->'
        nfldtx = 4
        call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

        if (qrderr) then
            go to 999
        end if

        write (noutpt,1016) uline1,uline2
        ustr = ufield(3)
        call chreal(nttyo,qrderr,ustr,var)

        if (qrderr) then
            go to 999
        end if

        ustr24 = ufield(2)(1:24)
        uheadx = 'Ref. Temp. (C)'
        call locase(ustr24)
        call locase(uheadx)
        j2 = ilnobl(ustr24)
        j3 = ilnobl(uheadx)

        if (ustr24(1:j2) .ne. uheadx(1:j3)) then
            write (noutpt,1020) uheadx(1:j3),ustr24(1:j2)
            write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
            qrderr = .true.
            go to 999
        end if

        if (.not.qnonep) then
            tgexp(ne) = var
        end if

        ! Loop on exchange sites.
        je = 0

        do jj = 1,jetmax + 1
            ! Read a line. If the sub-block (for the current exchanger
            ! phase) has not been completely read, this contains the name
            ! of an exchange site, and a sub-sub-block for that site
            ! follows. Otherwise, this line is the first line of the next
            ! sub-block (for the next exchanger phase).
            nfldtx = 0
            call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

            if (qrderr) then
                go to 999
            end if

            ustr = ufield(1)
            uheadx = '->'
            call locase(ustr)
            call locase(uheadx)
            j2 = ilnobl(ustr)
            j3 = ilnobl(uheadx)

            if (ustr(1:j2) .eq. uheadx(1:j3)) then
                write (noutpt,1014) uline1
                ustr24 = ufield(2)(1:24)
                uheadx = 'Exchange site'
                call locase(ustr24)
                call locase(uheadx)
                j2 = ilnobl(ustr24)
                j3 = ilnobl(uheadx)

                if (ustr24(1:j2) .ne. uheadx(1:j3)) then
                    write (noutpt,1020) uheadx(1:j3),ustr24(1:j2)
                    write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
                    qrderr = .true.
                    go to 999
                end if

                ustr = ufield(3)
                ustrn = ustr(1:24)
                call locase(ustrn)

                if (ustrn(1:5).eq.'none ' .or. ustr(1:1).eq.' ') then
                    ustr = 'None'
                end if

                qnones = ustr(1:5).eq.'None '

                qokay = .not.qnones .and. .not.qnonep

                if (qokay) then
                    ! Have found another exchange site sub-sub-block.
                    jgext(ne) = jgext(ne) + 1
                    je = je + 1

                    if (je .gt. jetmax) then
                        j2 = ilnobl(ugexp(ne))
                        write (noutpt,1710) jetmax,ugexp(ne)(1:j2),je
                        write (nttyo,1710) jetmax,ugexp(ne)(1:j2),je
1710 format(/' * Error - (EQ3NR/rd3ind) Have exceeded the',' maximum number of ',i3,/7x,'exchange sites on a',' generic ion exchange phase while reading',/7x,'the',' data to create ',a,'. Increase the',/7x,'dimensioning parameter jetpar to at least ',i3,'.')

                        go to 990
                    end if

                    ugexj(je,ne) = ufield(3)(1:8)
                    ngexrt(je,ne) = 0
                end if
            else
                ! Back up.
                backspace ninpts

                ustr = ufield(1)
                uheadx = 'Exchanger phase'
                call locase(ustr)
                call locase(uheadx)
                j2 = ilnobl(ustr)
                j3 = ilnobl(uheadx)

                if (ustr(1:j2) .eq. uheadx(1:j3)) then
                    ! Have found another exchanger phase.
                    go to 230
                end if

                ! Have found the end of the current block.
                go to 240
            end if

            ! Read the separator line following the data line containing
            ! the name of an exchange site.
            nfldtx = 1
            call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

            if (qrderr) then
                go to 999
            end if

            write (noutpt,1014) uline1
            ustr24 = ufield(1)(1:24)
            uheadx = '--------'

            if (ustr24(1:8) .ne. uheadx(1:8)) then
                j2 = ilnobl(uline1)
                j2 = min(j2,70)
                write (noutpt,1070) uline1(1:j2)
                write (nttyo,1070) uline1(1:j2)
                qrderr = .true.
                go to 999
            end if

            ! Read the site stoichiometric number (moles of site per mole
            ! of bare exchange ligand) from a two-line header.
            uheadx = '--->'
            nfldtx = 4
            call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

            if (qrderr) then
                go to 999
            end if

            write (noutpt,1016) uline1,uline2
            ustr24 = ufield(2)(1:24)
            uheadx = 'Stoich. number'
            call locase(ustr24)
            call locase(uheadx)
            j2 = ilnobl(ustr24)
            j3 = ilnobl(uheadx)

            if (ustr24(1:j2) .ne. uheadx(1:j3)) then
                write (noutpt,1020) uheadx(1:j3),ustr24(1:j2)
                write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
                qrderr = .true.
                go to 999
            end if

            ustr = ufield(3)
            call chreal(nttyo,qrderr,ustr,var)

            if (qrderr) then
                go to 999
            end if

            if (qokay) then
                cgexj(je,ne) = var
            end if

            ! Read the intrinsic electrical charge of the site (charge
            ! number for one mole of site) from a two-line header.
            uheadx = '--->'
            nfldtx = 4
            call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

            if (qrderr) then
                go to 999
            end if

            write (noutpt,1016) uline1,uline2
            ustr24 = ufield(2)(1:24)
            uheadx = 'Electr. charge'
            call locase(ustr24)
            call locase(uheadx)
            j2 = ilnobl(ustr24)
            j3 = ilnobl(uheadx)

            if (ustr24(1:j2) .ne. uheadx(1:j3)) then
                write (noutpt,1020) uheadx(1:j3),ustr24(1:j2)
                write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
                qrderr = .true.
                go to 999
            end if

            ustr = ufield(3)
            call chreal(nttyo,qrderr,ustr,var)

            if (qrderr) then
                go to 999
            end if

            if (qokay) then
                zgexj(je,ne) = var
            end if

            ! Loop on exchange reactions.
            n = 0

            do nnn = 1,ietmax + 1
                ! Read a line. If the sub-sub-block (for the current
                ! exchange site) has not been completely read, this contains
                ! a string denoting an exchange reaction in condensed format,
                ! and a sub-sub-sub-block for that reaction follows.
                ! Otherwise, this line is the first line of the next
                ! sub-sub-block (for the next site) or the first line of the
                ! next sub-block for the next exchanger phase).
                nfldtx = 0
                call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

                if (qrderr) then
                    go to 999
                end if

                ustr = ufield(1)
                uheadx = '--->'
                call locase(ustr)
                call locase(uheadx)
                j2 = ilnobl(ustr)
                j3 = ilnobl(uheadx)

                if (ustr(1:j2) .eq. uheadx(1:j3)) then
                    write (noutpt,1014) uline1
                    ustr24 = ufield(2)(1:24)
                    uheadx = 'Reaction'
                    call locase(ustr24)
                    call locase(uheadx)
                    j2 = ilnobl(ustr24)
                    j3 = ilnobl(uheadx)

                    if (ustr24(1:j2) .ne. uheadx(1:j3)) then
                        write (noutpt,1020) uheadx(1:j3),ustr24(1:j2)
                        write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
                        qrderr = .true.
                        go to 999
                    end if

                    ustr = ufield(3)
                    ustrn = ustr(1:24)
                    call locase(ustrn)

                    if (ustrn(1:5).eq.'none ' .or. ustr(1:1).eq.' ') then
                        ustr = 'None'
                    end if

                    qnoner = ustr(1:5).eq.'None '

                    qokay = .not.qnoner .and. .not.qnones .and. .not.qnonep

                    if (qokay) then
                        ! Have found another exchange reaction sub-sub-sub-block.
                        ngexrt(je,ne) = ngexrt(je,ne) + 1
                        n = n + 1

                        if (n .gt. ietmax) then
                            j2 = ilnobl(ugexp(ne))
                            j3 = ilnobl(ugexj(ne,ne))
                            write (noutpt,1820) netmax,ugexj(je,ne)(1:j3),ugexp(ne)(1:j2),n
                            write (nttyo,1820) netmax,ugexj(je,ne)(1:j3),ugexp(ne)(1:j2),n
1820 format(/' * Error - (EQ3NR/rd3ind) Have exceeded the',' maximum number of ',i3,/7x,'reactions for a site',' belonging to a generic ion exchange',/7x,'phase',' while reading the data for site ',a,' of exchange',' phase',/7x,a,'. Increase the dimensioning',' parameter',/7x,'ietpar to at least ',i3,'.')

                            go to 990
                        end if

                        ugexr(n,je,ne) = ufield(3)(1:56)
                    end if
                else
                    ! Back up.
                    backspace ninpts

                    ustr = ufield(1)
                    uheadx = '->'
                    call locase(ustr)
                    call locase(uheadx)
                    j2 = ilnobl(ustr)
                    j3 = ilnobl(uheadx)

                    if (ustr(1:j2) .eq. uheadx(1:j3)) then
                        ustr = ufield(2)
                        uheadx = 'Exchange site'
                        call locase(ustr)
                        call locase(uheadx)
                        j2 = ilnobl(ustr)
                        j3 = ilnobl(uheadx)

                        if (ustr(1:j2) .eq. uheadx(1:j3)) then
                            ! Have found another exchanger site.
                            go to 220
                        end if
                    end if

                    ustr = ufield(1)
                    uheadx = 'Exchanger phase'
                    call locase(ustr)
                    call locase(uheadx)
                    j2 = ilnobl(ustr)
                    j3 = ilnobl(uheadx)

                    if (ustr(1:j2) .eq. uheadx(1:j3)) then
                        ! Have found another exchanger phase.
                        go to 230
                    end if

                    ! Have found the end of the current block.
                    go to 240
                end if

                ! Read the separator line following the data line containing
                ! the string containing an exchange reaction in condensed
                ! format.
                nfldtx = 1
                call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

                if (qrderr) then
                    go to 999
                end if

                write (noutpt,1014) uline1
                ustr24 = ufield(1)(1:24)
                uheadx = '--------'

                if (ustr24(1:8) .ne. uheadx(1:8)) then
                    j2 = ilnobl(uline1)
                    j2 = min(j2,70)
                    write (noutpt,1070) uline1(1:j2)
                    write (nttyo,1070) uline1(1:j2)
                    qrderr = .true.
                    go to 999
                end if

                ! Read the table header for the thermodynamic data for the
                ! current exchange reaction from a two-line header.
                uheadx = '----->'
                nfldtx = 5
                call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

                if (qrderr) then
                    go to 999
                end if

                write (noutpt,1016) uline1,uline2
                ustr24 = ufield(2)(1:24)
                uheadx = 'Parameter'
                call locase(ustr24)
                call locase(uheadx)
                j2 = ilnobl(ustr24)
                j3 = ilnobl(uheadx)

                if (ustr24(1:j2) .ne. uheadx(1:j3)) then
                    write (noutpt,1020) uheadx(1:j3),ustr24(1:j2)
                    write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
                    qrderr = .true.
                    go to 999
                end if

                ! Read the log K/Delta G data for the reaction from a
                ! one-line header.
                uheadx = '----->'
                nfldtx = 5
                call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

                if (qrderr) then
                    go to 999
                end if

                write (noutpt,1014) uline1
                ustr24 = ufield(2)(1:24)
                uheadx = 'K func.'
                call locase(ustr24)
                call locase(uheadx)
                j2 = ilnobl(ustr24)
                j3 = ilnobl(uheadx)

                if (ustr24(1:j2) .ne. uheadx(1:j3)) then
                    write (noutpt,1020) uheadx(1:j3),ustr24(1:j2)
                    write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
                    qrderr = .true.
                    go to 999
                end if

                ustr = ufield(3)
                call chreal(nttyo,qrderr,ustr,var)

                if (qrderr) then
                    go to 999
                end if

                if (qokay) then
                    xlkgex(n,je,ne) = var
                    uxkgex(n,je,ne) = ufield(4)(1:8)
                end if

                ! Read the Delta H data for the reaction from a one-line
                ! header.
                uheadx = '----->'
                nfldtx = 5
                call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

                if (qrderr) then
                    go to 999
                end if

                write (noutpt,1014) uline1
                ustr24 = ufield(2)(1:24)
                uheadx = 'DelH0r'
                call locase(ustr24)
                call locase(uheadx)
                j2 = ilnobl(ustr24)
                j3 = ilnobl(uheadx)

                if (ustr24(1:j2) .ne. uheadx(1:j3)) then
                    write (noutpt,1020) uheadx(1:j3),ustr24(1:j2)
                    write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
                    qrderr = .true.
                    go to 999
                end if

                ustr = ufield(3)
                call chreal(nttyo,qrderr,ustr,var)

                if (qrderr) then
                    go to 999
                end if

                if (qokay) then
                    xhfgex(n,je,ne) = var
                    uhfgex(n,je,ne) = ufield(4)(1:8)
                end if

                ! Read the Delta V data for the reaction from a two-line
                ! header. The separator line completes the sub-sub-sub-block.
                uheadx = '----->'
                nfldtx = 5
                call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

                if (qrderr) then
                    go to 999
                end if

                write (noutpt,1016) uline1,uline2
                ustr24 = ufield(2)(1:24)
                uheadx = 'DelV0r'
                call locase(ustr24)
                call locase(uheadx)
                j2 = ilnobl(ustr24)
                j3 = ilnobl(uheadx)

                if (ustr24(1:j2) .ne. uheadx(1:j3)) then
                    write (noutpt,1020) uheadx(1:j3),ustr24(1:j2)
                    write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
                    qrderr = .true.
                    go to 999
                end if

                ustr = ufield(3)
                call chreal(nttyo,qrderr,ustr,var)

                if (qrderr) then
                    go to 999
                end if

                if (qokay) then
                    xvfgex(n,je,ne) = var
                    uvfgex(n,je,ne) = ufield(4)(1:8)
                end if
            end do

220 continue
        end do

230 continue
    end do

240 continue

    ! Add comment lines listing valid option strings.
    write (noutpt,1682)
1682 format('* Valid units strings (uxkgex(i,j,n)/uhfgex(i,j,n)/','uvfgex(i,j,n)) are:',t80,'*')

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
            write (noutpt,1684) ux80(1:80)
1684 format(a)

            ux80(2:79) = ' '
            j3 = 6
            k = 0
        end if
    end do

    if (k .gt. 0) then
        write (noutpt,1684) ux80(1:80)
    end if

    write (noutpt,2052)

575 continue

    ! Specified generic ion exchanger compositions.
    neti = 0
    qnonep = .false.

    ! Read a two-line header for the block.
    uheadx = 'Ion Exchanger Compositions'
    nfldtx = 2
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    write (noutpt,1016) uline1,uline2

    ! Loop on exchanger phases.
    nei = 0

    do nn = 1,netmax + 1
        ! Read a line. If the block has not been completely read,
        ! this contains the name of an exchanger phase, and a sub-block
        ! for that phase follows. Otherwise, this line is the first line
        ! of the next block.
        nfldtx = 0
        call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)
        ustr = ufield(1)
        uheadx = 'Exchanger phase'
        call locase(ustr)
        call locase(uheadx)
        j2 = ilnobl(ustr)
        j3 = ilnobl(uheadx)

        if (ustr(1:j2) .ne. uheadx(1:j3)) then
            ! Back up.
            backspace ninpts
            go to 340
        end if

        write (noutpt,1014) uline1
        ustr = ufield(2)
        ustrn = ustr(1:24)
        call locase(ustrn)

        if (ustrn(1:5).eq.'none ' .or. ustr(1:1).eq.' ') then
            ustr = 'None'
        end if

        qnonep = ustr(1:5).eq.'None '

        ! XX
        if (qnonep) then
            qgexef = .true.
        else
            neti = neti + 1
            nei = nei + 1

            if (nei .gt. netmax) then
                write (noutpt,1930) netmax,nei
                write (nttyo,1930) netmax,nei
1930 format(/' * Error - (EQ3NR/rd3ind) Have exceeded the',' maximum number of ',i3,/7x,'generic ion exchange phases',' while reading the data for concentrations',7x,'and',' compositions of such phases. Increase the dimensioning',/7x,'parameter netpar to at least ',i3,'.')

                go to 990
            end if

            ugexpi(nei) = ustr(1:24)

            ! XX
            j3 = ilnobl(ugexpi(nei))
            jgexti(nei) = 0

            ! Find the corresponding ne index.
            do ne = 1,net
                j2 = ilnobl(ugexp(ne))

                if (ugexp(nei)(1:j3) .eq. ugexp(ne)(1:j2)) then
                    go to 250
                end if
            end do

            write (noutpt,1932) ugexpi(nei)(1:j3)
            write (nttyo,1932) ugexpi(nei)(1:j3)
1932 format(/' * Error - (EQ3NR/rd3ind) Data are present on the',' input file for the',/7x,'generic ion exchange phase "',a,'", but this phase',/7x,"hasn't been previously defined."," Can't determine the requisite model,",/7x,"therefore can't",' finish reading the current data for this phase.')

            qrderr = .true.
            stop

250 continue
            j2 = ilnobl(ugexmo(ne))

            if (ugexmo(ne)(1:j2) .eq. 'Gapon' .or.      ugexmo(ne)(1:6) .eq. 'Gapon-' .or.      ugexmo(ne)(1:j2) .eq. 'Vanselow' .or.      ugexmo(ne)(1:9) .eq. 'Vanselow-') then
                ! Input composition is described in terms of equivalent
                ! fractions on the sites.
                qgexef = .true.
            else
                ! Input composition is described in terms of mole fractions
                ! on the sites.
                qgexef = .false.
            end if
        end if

        ! Read the separator line following the data line containing
        ! the name of an exchanger phase.
        nfldtx = 1
        call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

        if (qrderr) then
            go to 999
        end if

        write (noutpt,1014) uline1
        ustr24 = ufield(1)(1:24)
        uheadx = '--------'

        if (ustr24(1:8) .ne. uheadx(1:8)) then
            j2 = ilnobl(uline1)
            j2 = min(j2,70)
            write (noutpt,1070) uline1(1:j2)
            write (nttyo,1070) uline1(1:j2)
            qrderr = .true.
            go to 999
        end if

        ! Read the concentration of the exchanger phase (moles/kg.H2O)
        ! from a two-line header.
        uheadx = '->'
        nfldtx = 4
        call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

        if (qrderr) then
            go to 999
        end if

        write (noutpt,1016) uline1,uline2
        ustr24 = ufield(2)(1:24)
        uheadx = 'Moles/kg.H2O'
        call locase(ustr24)
        call locase(uheadx)
        j2 = ilnobl(ustr24)
        j3 = ilnobl(uheadx)

        if (ustr24(1:j2) .ne. uheadx(1:j3)) then
            write (noutpt,1020) uheadx(1:j3),ustr24(1:j2)
            write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
            qrderr = .true.
            go to 999
        end if

        ustr = ufield(3)
        call chreal(nttyo,qrderr,ustr,var)

        if (qrderr) then
            go to 999
        end if

        if (.not.qnonep) then
            cgexpi(nei) = var
        end if

        ! Loop on exchange sites.
        jei = 0

        do jj = 1,jetmax + 1
            ! Read a line. If the sub-block (for the current exchanger
            ! phase) has not been completely read, this contains the name
            ! of an exchange site, and a sub-sub-block for that site
            ! follows. Otherwise, this line is the first line of the next
            ! sub-block (for the next exchanger phase).
            nfldtx = 0
            call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)
            ustr = ufield(1)
            uheadx = '->'
            call locase(ustr)
            call locase(uheadx)
            j2 = ilnobl(ustr)
            j3 = ilnobl(uheadx)

            if (ustr(1:j2) .eq. uheadx(1:j3)) then
                write (noutpt,1014) uline1
                ustr24 = ufield(2)(1:24)
                uheadx = 'Exchange site'
                call locase(ustr24)
                call locase(uheadx)
                j2 = ilnobl(ustr24)
                j3 = ilnobl(uheadx)

                if (ustr24(1:j2) .ne. uheadx(1:j3)) then
                    write (noutpt,1020) uheadx(1:j3),ustr24(1:j2)
                    write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
                    qrderr = .true.
                    go to 999
                end if

                ustr = ufield(3)
                ustrn = ustr(1:24)
                call locase(ustrn)

                if (ustrn(1:5).eq.'none ' .or. ustr(1:1).eq.' ') then
                    ustr = 'None'
                end if

                qnones = ustr(1:5).eq.'None '

                qokay = .not.qnones .and. .not.qnonep

                if (qokay) then
                    ! Have found another exchange site sub-sub-block.
                    jgexti(ne) = jgexti(ne) + 1
                    jei = jei + 1

                    if (jei .gt. jetmax) then
                        j2 = ilnobl(ugexpi(nei))
                        write (noutpt,2000) jetmax,ugexpi(nei)(1:j2),jei
                        write (nttyo,2000) jetmax,ugexpi(nei)(1:j2),jei
2000 format(/' * Error - (EQ3NR/rd3ind) Have exceeded the',' maximum number of ',i3,/7x,'exchange sites on a',' generic ion exchange phase while reading',/7x,'the',' data for the concentration and composition of',/7x,a,'. Increase the dimensioning parameter jetpar',/7x,'to at least ',i3,'.')

                        go to 990
                    end if

                    ugexji(jei,nei) = ufield(3)(1:8)
                    ngexti(jei,nei) = 0
                end if
            else
                ! Back up.
                backspace ninpts

                ustr = ufield(1)
                uheadx = 'Exchanger phase'
                call locase(ustr)
                call locase(uheadx)
                j2 = ilnobl(ustr)
                j3 = ilnobl(uheadx)

                if (ustr(1:j2) .eq. uheadx(1:j3)) then
                    ! Have found another exchanger phase.
                    go to 330
                end if

                ! Have found the end of the current block.
                go to 340
            end if

            ! Read the separator line following the data line containing
            ! the name of an exchange site.
            nfldtx = 1
            call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

            if (qrderr) then
                go to 999
            end if

            write (noutpt,1014) uline1
            ustr24 = ufield(1)(1:24)
            uheadx = '--------'

            if (ustr24(1:8) .ne. uheadx(1:8)) then
                j2 = ilnobl(uline1)
                j2 = min(j2,70)
                write (noutpt,1070) uline1(1:j2)
                write (nttyo,1070) uline1(1:j2)
                qrderr = .true.
                go to 999
            end if

            ! Read a table header for the exchange site composition from
            ! from a two-line header.
            uheadx = '--->'
            nfldtx = 4
            call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

            if (qrderr) then
                go to 999
            end if

            write (noutpt,1016) uline1,uline2

            ustr24 = ufield(2)(1:24)
            uheadx = 'Exchange species'
            call locase(ustr24)
            call locase(uheadx)
            j2 = ilnobl(ustr24)
            j3 = ilnobl(uheadx)

            if (ustr24(1:j2) .ne. uheadx(1:j3)) then
                write (noutpt,1020) uheadx(1:j3),ustr24(1:j2)
                write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
                qrderr = .true.
                go to 999
            end if

            ! Check the measure of concentration: equivalent fraction
            ! versus mole fraction.
            ustr24 = ufield(3)(1:24)

            if (qnonep) then
                if (ustr24(1:10) .eq. 'Mole frac.') then
                    ustr24 = 'Eq. frac.'
                end if
            end if

            if (qgexef) then
                uheadx = 'Eq. frac.'
            else
                uheadx = 'Mole frac.'
            end if

            j2 = ilnobl(ustr24)
            j3 = ilnobl(uheadx)

            if (ustr24(1:j2) .ne. uheadx(1:j3)) then
                write (noutpt,2010) uheadx(1:j3),ustr24(1:j2)
                write (nttyo,2010) uheadx(1:j3),ustr24(1:j2)
2010 format(/' * Error - (EQ3NR/rd3ind) Was expecting to find',' a line containing',/7x,'"',a,'", instead found one',' containing "',a,'".')

                if (.not.qnonep) then
                    j4 = ilnobl(ugexmo(ne))
                    j5 = ilnobl(ugexpi(nei))
                    write (noutpt,2012) ugexmo(ne)(1:j4),ugexpi(nei)(1:j5)
                    write (nttyo,2012) ugexmo(ne)(1:j4),ugexpi(nei)(1:j5)
2012 format(7x,'The former is required for the ',a,' exchange model,',/7x,'which was specified for',' the generic ion exchange phase',/7x,a,'.')
                end if

                qrderr = .true.
                go to 999
            end if

            ! Loop on exchange species.
            iei = 0

            do ii = 1,ietmax + 1
                ! Read a line. If the sub-sub-block (for the current
                ! exchange site) has not been completely read, this contains
                ! the name of the exchange species and its mole fraction.
                ! Otherwise, this line is the first line of the next
                ! sub-sub-block (for the next site) or the first line of the
                ! next sub-block for the next exchanger phase).
                nfldtx = 0
                call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

                if (qrderr) then
                    go to 999
                end if

                write (noutpt,1014) uline1
                ustr = ufield(1)
                uheadx = '--->'
                call locase(ustr)
                call locase(uheadx)
                j2 = ilnobl(ustr)
                j3 = ilnobl(uheadx)

                if (ustr(1:j2) .eq. uheadx(1:j3)) then
                    ustr = ufield(2)
                    ustrn = ustr(1:24)
                    call locase(ustrn)

                    if (ustrn(1:5).eq.'none ' .or. ustr(1:1).eq.' ') then
                        ustr = 'None'
                    end if

                    qnonei = ustr(1:5).eq.'None '

                    qokay = .not.qnonei .and. .not.qnones .and. .not.qnonep

                    if (qokay) then
                        ! Have found another exchange species line.
                        ngexti(jei,nei) = ngexti(jei,nei) + 1
                        iei = iei + 1

                        if (iei .gt. ietmax) then
                            write (noutpt,2040) ietmax,ugexji(jei,nei),ugexpi(nei),iei
                            write (nttyo,2040) ietmax,ugexji(jei,nei),ugexpi(nei),iei
2040 format(/' * Error - (EQ3NR/rd3ind) Have exceeded the',' maximum number of ',i3,/7x,'species on an exchange',' site of a generic ion exchange phase',/7x,'while',' reading',/7x,'the data for the composition of',' site ',a,' of',/7x,'exchange phase ',a,'. Increase',' the dimensioning parameter ietpar to at least ',i3,'.')

                            go to 990
                        end if

                        ugexsi(iei,jei,nei) = ufield(2)(1:24)
                        ustr = ufield(3)
                        call chreal(nttyo,qrderr,ustr,var)

                        if (qrderr) then
                            go to 999
                        end if

                        ! Check the array name noted: egexsi versus xgexsi.
                        ustr24 = ' '
                        k2 = index(uline1,'egexsi')
                        k3 = index(uline1,'xgexsi')

                        if (k2 .gt. 0) then
                            ustr24 = 'egexsi'
                        end if

                        if (k3 .gt. 0) then
                            ustr24 = 'xgexsi'
                        end if

                        if (qgexef) then
                            uheadx = 'egexsi'
                        else
                            uheadx = 'xgexsi'
                        end if

                        j2 = ilnobl(ustr24)
                        j3 = ilnobl(uheadx)
                        j4 = ilnobl(ugexmo(ne))
                        j5 = ilnobl(ugexpi(nei))

                        if (ustr24(1:j2) .ne. uheadx(1:j3)) then
                            write (noutpt,2014) uheadx(1:j3),ustr24(1:j2),ugexmo(ne)(1:j4),ugexpi(nei)(1:j5)
                            write (nttyo,2014) uheadx(1:j3),ustr24(1:j2),ugexmo(ne)(1:j4),ugexpi(nei)(1:j5)
2014 format(/' * Warning - (EQ3NR/rd3ind) Was expecting',' to find a line containing',/7x,'"',a,'", instead',' found one containing "',a,'".',/7x,'The former',' matches the ',a,' exchange model, which was',/7x,'specified for the generic ion exchange phase ',a,'.')
                        end if

                        if (qgexef) then
                            egexsi(iei,jei,nei) = var
                        else
                            xgexsi(iei,jei,nei) = var
                        end if
                    end if
                else
                    uheadx = '--------'

                    if (ustr(1:8) .eq. uheadx(1:8)) then
                        ! Have found the end of the sub-sub-block for the
                        ! composition of the current exchange site
                        go to 320
                    else
                        ! Have unrecognized input.
                        j3 = ilnobl(ugexji(jei,nei))
                        j4 = ilnobl(ugexpi(nei))
                        write (noutpt,1360) ustr(1:j2),ugexji(jei,nei)(1:j3),ugexpi(nei)(1:j4)
                        write (nttyo,1360) ustr(1:j2),ugexji(jei,nei)(1:j3),ugexpi(nei)(1:j4)
1360 format(/' * Error - (EQ3NR/rd3ind) Have unrecognized',' input on the line',/7x,'beginning with "',a,'".',' This line should contain the',/7x,'name of an',' exchanger species and corresponding mole fraction',/7x,'for site ',a,' of exchanger phase ',a,', else it',/7x,'should be a separator line terminating a sequence',' of such lines.')
                    end if
                end if
            end do

320 continue
        end do

330 continue
    end do

340 continue

    ! Specified solid solution compositions.
    nxti = 0
    nxic = 0
    qnonep = .false.

    ! Read a two-line header for the block.
    uheadx = 'Solid Solution Compositions'
    nfldtx = 2
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    write (noutpt,1016) uline1,uline2

    ! Loop on solid solution phases.
    nxi = 0

    do nn = 1,nxtimx + 1
        ! Read a line. If the block has not been completely read,
        ! this contains the name of a solid solution, and a sub-block
        ! for that phase follows. Otherwise, this line is the first line
        ! of the next block.
        nfldtx = 0
        call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)
        ustr = ufield(1)
        uheadx = 'Solid Solution'
        call locase(ustr)
        call locase(uheadx)
        j2 = ilnobl(ustr)
        j3 = ilnobl(uheadx)

        if (ustr(1:j2) .ne. uheadx(1:j3)) then
            ! Back up.
            backspace ninpts
            go to 440
        end if

        write (noutpt,1014) uline1
        ustr = ufield(2)
        ustrn = ustr(1:24)
        call locase(ustrn)

        if (ustrn(1:5).eq.'none ' .or. ustr(1:1).eq.' ') then
            ustr = 'None'
        end if

        qnonep = ustr(1:5).eq.'None '

        if (.not.qnonep) then
            nxti = nxti + 1
            nxi = nxi + 1

            if (nxi .gt. nxtimx) then
                write (noutpt,1430) nxtimx,ustr(1:j2)
                write (nttyo,1430) nxtimx,ustr(1:j2)
1430 format(/' * Error - (EQ3NR/rd3ind) Have exceeded the',' maximum',/7x,'number of ',i5,' solid solutions for which',/7x,'compositions may be specified on the input file. This',/7x,'occurred while reading data for ',a,'.',/7x,'Increase the dimensioning parameter nxtipa.')

                go to 990
            end if

            usoli(nxi) = ustr(1:24)
            ncmpri(1,nxi) = nxic + 1
        end if

        ! Read the separator line following the data line containing
        ! the name of a solid solution.
        nfldtx = 1
        call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

        if (qrderr) then
            go to 999
        end if

        write (noutpt,1014) uline1
        ustr24 = ufield(1)(1:24)
        uheadx = '--------'

        if (ustr24(1:8) .ne. uheadx(1:8)) then
            j2 = ilnobl(uline1)
            j2 = min(j2,70)
            write (noutpt,1070) uline1(1:j2)
            write (nttyo,1070) uline1(1:j2)
            qrderr = .true.
            go to 999
        end if

        ! Read a table header for the solid solution composition from
        ! from a two-line header.
        uheadx = '->'
        nfldtx = 4
        call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

        if (qrderr) then
            go to 999
        end if

        write (noutpt,1016) uline1,uline2
        ustr24 = ufield(2)(1:24)
        uheadx = 'Component'
        call locase(ustr24)
        call locase(uheadx)
        j2 = ilnobl(ustr24)
        j3 = ilnobl(uheadx)

        if (ustr24(1:j2) .ne. uheadx(1:j3)) then
            write (noutpt,1020) uheadx(1:j3),ustr24(1:j2)
            write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
            qrderr = .true.
            go to 999
        end if

        ! Loop on component species.
        do nnn = 1,nxicmx + 1
            ! Read a line. If the sub-block (for the current
            ! solid solution) has not been completely read, this contains
            ! the name of the component species and its mole fraction.
            ! Otherwise, this line is the first line of the next
            ! sub-block (for the next solid solution).
            nfldtx = 0
            call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)
            write (noutpt,1014) uline1
            ustr = ufield(1)
            uheadx = '->'
            call locase(ustr)
            call locase(uheadx)
            j2 = ilnobl(ustr)
            j3 = ilnobl(uheadx)

            if (ustr(1:j2) .eq. uheadx(1:j3)) then
                ustr = ufield(2)
                ustrn = ustr(1:24)
                call locase(ustrn)

                if (ustrn(1:5).eq.'none ' .or. ustr(1:1).eq.' ') then
                    ustr = 'None'
                end if

                qnonei = ustr(1:5).eq.'None '

                qokay = .not.qnonep .and. .not.qnonei

                if (qokay) then
                    ! Have found another component species line.
                    nxic = nxic + 1

                    if (nxic .gt. nxicmx) then
                        j2 = ilnobl(ustr)
                        j3 = ilnobl(usoli(nxi))
                        write (noutpt,1460) nxicmx,ustr(1:j2),usoli(nxi)(1:j3)
                        write (nttyo,1460) nxicmx,ustr(1:j2),usoli(nxi)(1:j3)
1460 format(/' * Error - (EQ3NR/rd3ind) Have exceeded the',' maximum number',/7x,'of ',i5,' solid solution',' species for which mole fractions',/7x,'are specified',' on the input file. This occurred while reading',/7x,'data for the species ',a,'(',a,').',/7x,'Increase the dimensioning parameter nxipar.')

                        go to 990
                    end if

                    umemi(nxic) = ufield(2)(1:24)
                    ustr = ufield(3)
                    call chreal(nttyo,qrderr,ustr,var)

                    if (qrderr) then
                        go to 999
                    end if

                    xbari(nxic) = var
                end if
            else
                uheadx = '--------'

                if (ustr(1:8) .eq. uheadx(1:8)) then
                    ! Have found the end of the sub-block for the
                    ! composition of the current solid solution.
                    go to 420
                else
                    ! Have unrecognized input.
                    j3 = ilnobl(usoli(nxi))
                    write (noutpt,1470) ustr(1:j2),usoli(nxi)(1:j3)
                    write (nttyo,1470) ustr(1:j2),usoli(nxi)(1:j3)
1470 format(/' * Error - (EQ3NR/rd3ind) Have unrecognized',' input on the line',/7x,'beginning with "',a,'".',' This line should contain the',/7x,'name of a',' component species and corresponding mole fraction',/7x,'for solid solution ',a,', else it should be a',' separator line terminating a sequence of such lines.')
                end if
            end if
        end do

420 continue
        if (.not.qnonep) then
            ncmpri(2,nxi) = nxic
        end if
    end do

440 continue

    ! Nxmod options.
    uheadx = 'Alter/Suppress options'
    nfldtx = 2
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    write (noutpt,1016) uline1,uline2

    ! Read the first part of a table header from a one-line header.
    uheadx = 'Species'
    nfldtx = 3
    call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)

    if (qrderr) then
        go to 999
    end if

    write (noutpt,1014) uline1

    ! Read the second part of the table header from a two-line header.
    uheadx = '(uxmod(n))'
    nfldtx = 3
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    write (noutpt,1016) uline1,uline2

    n = 0

    ! Loop on the number of alter/suppress options.
    do nn = 1,nxmdmx + 1
        ! Read a line. This contains an alter/suppress option, else it is
        ! a separator line marking the end of the current block.
        nfldtx = 0
        call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

        if (qrderr) then
            go to 999
        end if

        write (noutpt,1014) uline1

        ustr = ufield(1)

        if (ustr(1:8) .eq. '--------') then
            go to 570
        end if

        call locase(ustr)
        ustrn = ustr(1:24)
        call locase(ustrn)

        if (ustrn(1:5).eq.'none ' .or. ustr(1:1).eq.' ') then
            ustr = 'None'
        end if

        if (ustr(1:5) .eq. 'None ') then
            go to 560
        end if

        n = n + 1

        if (n .gt. nxmdmx) then
            write (noutpt,1500) nxmdmx
            write (nttyo,1500) nxmdmx
1500 format(/' * Error - (EQ3NR/rd3ind) Have too many nxmod',/7x,'alter/suppress options. The code is only dimensioned',/7x,'for ',i3,' such options. Reduce the number of such',/7x,'options or increase the dimensioning parameter nxmdpa.')

            go to 990
        end if

        uxmod(n) = ufield(1)(1:48)

        ustr = ufield(2)
        call locase(ustr)

        do kxmd = -1,2
            uheadx = ukxm(kxmd)
            call locase(uheadx)

            if (ustr(1:16) .eq. uheadx(1:16)) then
                kxmod(n) = kxmd
                go to 152
            end if
        end do

        j2 = ilnobl(ustr)
        write (nttyo,1510) ustr(1:j2)
1510 format(/" * Error - (EQ3NR/rd3ind) Don't recognize the",' alter/suppress option',/7x,'string "',a,'". This should',' be one of the strings',/7x,'defined in the ukxm array.',' The valid strings are:',/)

        do kxmd = -1,2
            j3 = ilnobl(ukxm(kxmd))
            write (nttyo,1514) ukxm(kxmd)(1:j3)
1514 format(9x,a)
        end do

        go to 990

152 continue
        ustr = ufield(3)
        call chreal(nttyo,qrderr,ustr,var)

        if (qrderr) then
            go to 999
        end if

        xlkmod(n) = var

560 continue
    end do

570 continue
    nxmod = n

    ! Add comment lines listing valid option strings.
    write (noutpt,3092)
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
            write (noutpt,3094) ux80(1:80)
3094 format(a)

            ux80(2:79) = ' '
            j3 = 6
            k = 0
        end if
    end do

    if (k .gt. 0) then
        write (noutpt,3094) ux80(1:80)
    end if

    write (noutpt,2052)

    ! Iopt Model Option Switches.
    ! Note: iopt(1) = iopt1, etc.
    uheadx = 'Iopt Model Option Switches ("( 0)" marks default choices)'
    nfldtx = 1
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    write (noutpt,1016) uline1,uline2

    ! Loop on option switches.
    do nn = 1,noptmx
        ! Read the option title string from a one-line header.
        nfldtx = 1
        call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

        if (qrderr) then
            go to 999
        end if

        uheadx = ufield(1)
        j2 = ilnobl(uheadx)

        ! Check for the end of the block.
        if (uheadx(1:4) .ne. 'iopt') then
            ! Back up.
            backspace ninpts
            go to 650
        end if

        write (noutpt,1014) uline1
        k1 = index(uheadx,'(')
        k2 = index(uheadx,')')
        k3 = index(uheadx,'- ')

        if (k1.le.0 .or. k2.le.0 .or. k2.lt.(k1 + 2) .or. k3.le.0 .or.  k3.le.k2) then
            write (noutpt,1600) uheadx(1:j2)
            write (nttyo,1600) uheadx(1:j2)
1600 format(/' * Error - (EQ3NR/rd3ind) The iopt option switch',' title line',/7x,'"',a,'"',/7x,'read from the'," input file isn't in the required format.")

            go to 990
        end if

        ! Get the index of the option.
        ustr = uheadx(k1 + 1:k2 - 1)
        call chrint(ivar,nttyo,qrderr,ustr)

        if (qrderr) then
            go to 999
        end if

        n = ivar

        if (n.lt.1 .or. n.gt.nptxpa) then
            write (noutpt,1610) uheadx(1:j2),n
            write (nttyo,1610) uheadx(1:j2),n
1610 format(/' * Error - (EQ3NR/rd3ind) The iopt option switch',' title line',/7x,'"',a,'"',/7x,'read from the',' input file references an option switch index',/7x,'of ',i3,', which is out of range.')

            go to 990
        end if

        ! Check the full title string.
        ustr = uopttx(n)
        j3 = ilnobl(ustr)

        if (index(uheadx(1:j2),ustr(1:j3)) .le. 0) then
            write (noutpt,1612) uheadx(1:j2),ustr(1:j3)
            write (nttyo,1612) uheadx(1:j2),ustr(1:j3)
1612 format(/' * Error - (EQ3NR/rd3ind) The iopt option switch',' title string',/7x,'"',a,'"',/7x,"read from the input file doesn't contain the",' matching defined string',/7x,'"',a,'".')

            go to 990
        end if

        iopt(n) = 0
        nmark = 0

        ! Loop on option choices.
        do jj = 1,jptxpa + 1
            ! Read a line. This contains an option choice line, else it
            ! is a separator line marking the end of the option choice
            ! lines for the current option.
            nfldtx = 1
            call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

            if (qrderr) then
                go to 999
            end if

            write (noutpt,1014) uline1

            uheadx = ufield(1)
            j2 = ilnobl(uheadx)

            if (uheadx(1:8) .eq. '--------') then
                go to 620
            end if

            if (jj .gt. jptxpa) then
                j3 = ilnobl(uopttx(n))
                write (noutpt,1630) uopttx(n)(1:j3),jptxpa
                write (nttyo,1630) uopttx(n)(1:j3),jptxpa
1630 format(/' * Error - (EQ3NR/rd3ind) Have too many option',/7x,'choice lines for the iopt option switch whose title',' string is',/7x,'"',a,'".',/7x,' The code is only',' dimensioned for ',i3,' such lines. Reduce the',/7x,'number of such lines or increase the dimensioning',' parameter jptxpa.')

                go to 990
            end if

            k1 = index(uheadx,'[')
            k2 = index(uheadx,']')
            k3 = index(uheadx,'(')
            k4 = index(uheadx,')')

            if (k1.le.0 .or. k2.le.0 .or. k2.lt.(k1 + 2) .or.    k3.le.0 .or. k4.le.0 .or. k4.lt.(k3 + 2)) then
                write (noutpt,1640) uheadx(1:j2)
                write (nttyo,1640) uheadx(1:j2)
1640 format(/' * Error - (EQ3NR/rd3ind) The following iopt',/7x,'option switch choice line read from the input file',/7x,"isn't in the required format:",/7x,'"',a,'"')

                go to 990
            end if

            ! Get the index of the option choice.
            ustr = uheadx(k3 + 1:k4 - 1)
            call chrint(ivar,nttyo,qrderr,ustr)

            if (qrderr) then
                go to 999
            end if

            ival = ivar

            do j = 1,jptxpa
                ii = ioptox(j,n)

                if (ival .eq. ii) then
                    go to 600
                end if
            end do

            write (noutpt,1650) uheadx(1:j2)
            write (nttyo,1650) uheadx(1:j2)
1650 format(/' * Error - (EQ3NR/rd3ind) The iopt option switch',/7x,'choice line',/7x,'"',a,'"',/7x,'read from the input file references an out-of-range',' option',/7x,'choice index.')

            go to 990

600 continue

            ! Check the full option choice string.
            ustr = uoptox(j,n)
            j3 = ilnobl(ustr)
            j4 = ilnobl(uopttx(n))

            if (index(uheadx(1:j2),ustr(1:j3)) .le. 0) then
                write (noutpt,1660) uheadx(1:j2),ustr(1:j3),uopttx(n)(1:j4)
                write (nttyo,1660) uheadx(1:j2),ustr(1:j3),uopttx(n)(1:j4)
1660 format(/' * Error - (EQ3NR/rd3ind) The iopt option switch',' choice line',/7x,'"',a,'"',/7x," read from the input file doesn't contain the",' matching defined string',/7x,'"',a,'".',/7x,'This line belongs to the option switch whose title',' string is',/7x,'"',a,'".')

                go to 990
            end if

            k1 = index(uheadx,'[')
            k2 = index(uheadx,']')
            ustr24 = uheadx(k1 + 1:k2 -1)
            call lejust(ustr24)
            j3 = ilnobl(ustr24)
            qmark = .false.

            if (j3 .gt. 0) then
                if (index(ustr24,'*') .ge. 1) then
                    qmark = .true.
                else if (index(ustr24,'x') .ge. 1) then
                    qmark = .true.
                else if (index(ustr24,'X') .ge. 1) then
                    qmark = .true.
                else
                    j4 = ilnobl(uopttx(n))
                    write (noutpt,1670) ustr24(1:j3),uheadx(1:j2),uopttx(n)(1:j4)
                    write (nttyo,1670) ustr24(1:j3),uheadx(1:j2),uopttx(n)(1:j4)
1670 format(/" * Error - (EQ3NR/rd3ind) Don't recognize the",' string "',a,'"',/7x,'that appears on the iopt',' option switch choice line',/7x,'"',a,'"',/7x,'read from the input file. An option choice should',' be chosen by',/7x,'placing a "*", "x", or "X" in the',' checkbox ("[ ]"). This',/7x,'choice line belongs to',' the option switch whose title string is',/7x,'"',a,'".')

                    go to 990
                end if
            end if

            if (qmark) then
                nmark = nmark + 1
                iopt(n) = ival
                jlast = j
            end if
        end do

620 continue

        if (nmark .le. 0) then
            j2 = ilnobl(uopttx(n))
            write (noutpt,1680) uopttx(n)(1:j2)
            write (nttyo,1680) uopttx(n)(1:j2)
1680 format(/' * Warning - (EQ3NR/rd3ind) No option choice was',' checked on the input file',/7x,'for the iopt option',' switch whose title string is',/7x,'"',a,'".')

            do j = 1,jptxpa
                ival = ioptox(j,n)

                if (ival .eq. 0) then
                    go to 630
                end if
            end do

            write (noutpt,1690)
            write (nttyo,1690)
1690 format(/7x,'A default value of 0 has been applied, but no',/7x,'matching option choice string is defined.')

            go to 640

630 continue
            j3 = ilnobl(uoptox(j,n))
            write (noutpt,1692) uoptox(j,n)(1:j3)
            write (nttyo,1692) uoptox(j,n)(1:j3)
1692 format(/7x,'A default value of 0 has been applied. The',' matching string is',/7x,'"',a,'".')

640 continue
        end if

        if (nmark .gt. 1) then
            j2 = ilnobl(uopttx(n))
            j = jlast
            j3 = ilnobl(uoptox(j,n))
            write (noutpt,1694) uopttx(n)(1:j2),uoptox(j,n)(1:j3)
            write (nttyo,1694) uopttx(n)(1:j2),uoptox(j,n)(1:j3)
1694 format(/' * Warning - (EQ3NR/rd3ind) More than one option',' choice was checked',/7x,'on the input file for the iopt',' option switch whose title string is',/7x,'"',a,'".',/7x,'The last choice checked will be used. The',' matching string is',/7x,'"',a,'".')
        end if
    end do

650 continue

    ! Iopg Print Option Switches.
    ! Note: iopg(1) = iopg1, etc.
    uheadx = 'Iopg Activity Coefficient Option Switches ("( 0)" marks default choices)'
    nfldtx = 1
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    write (noutpt,1016) uline1,uline2

    ! Loop on option switches.
    do nn = 1,nopgmx
        ! Read the option title string from a one-line header.
        nfldtx = 1
        call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

        if (qrderr) then
            go to 999
        end if

        uheadx = ufield(1)
        j2 = ilnobl(uheadx)

        ! Check for the end of the block.
        if (uheadx(1:4) .ne. 'iopg') then
            ! Back up.
            backspace ninpts
            go to 750
        end if

        write (noutpt,1014) uline1
        k1 = index(uheadx,'(')
        k2 = index(uheadx,')')
        k3 = index(uheadx,'- ')

        if (k1.le.0 .or. k2.le.0 .or. k2.lt.(k1 + 2) .or. k3.le.0 .or.  k3.le.k2) then
            write (noutpt,1700) uheadx(1:j2)
            write (nttyo,1700) uheadx(1:j2)
1700 format(/' * Error - (EQ3NR/rd3ind) The iopg option switch',' title line',/7x,'"',a,'"',/7x,'read from the'," input file isn't in the required format.")

            go to 990
        end if

        ! Get the index of the option.
        ustr = uheadx(k1 + 1:k2 - 1)
        call chrint(ivar,nttyo,qrderr,ustr)

        if (qrderr) then
            go to 999
        end if

        n = ivar

        if (n.lt.1 .or. n.gt.npgxpa) then
            write (noutpt,1712) uheadx(1:j2),n
            write (nttyo,1712) uheadx(1:j2),n
1712 format(/' * Error - (EQ3NR/rd3ind) The iopg option switch',' title line',/7x,'"',a,'"',/7x,'read from the',' input file references an option switch index',/7x,'of ',i3,', which is out of range.')

            go to 990
        end if

        ! Check the full title string.
        ustr = uopgtx(n)
        j3 = ilnobl(ustr)

        if (index(uheadx(1:j2),ustr(1:j3)) .le. 0) then
            write (noutpt,1720) uheadx(1:j2),ustr(1:j3)
            write (nttyo,1720) uheadx(1:j2),ustr(1:j3)
1720 format(/' * Error - (EQ3NR/rd3ind) The iopg option switch',' title string',/7x,'"',a,'"',/7x,"read from the input file doesn't contain the",' matching defined string',/7x,'"',a,'".')

            go to 990
        end if

        iopg(n) = 0
        nmark = 0

        ! Loop on option choices.
        do jj = 1,jpgxpa + 1
            ! Read a line. This contains an option choice line, else it
            ! is a separator line marking the end of the option choice
            ! lines for the current option.
            nfldtx = 1
            call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

            if (qrderr) then
                go to 999
            end if

            write (noutpt,1014) uline1

            uheadx = ufield(1)
            j2 = ilnobl(uheadx)

            if (uheadx(1:8) .eq. '--------') then
                go to 720
            end if

            if (jj .gt. jpgxpa) then
                j3 = ilnobl(uopgtx(n))
                write (noutpt,1730) uopgtx(n)(1:j3),jpgxpa
                write (nttyo,1730) uopgtx(n)(1:j3),jpgxpa
1730 format(/' * Error - (EQ3NR/rd3ind) Have too many option',/7x,'choice lines for the iopg option switch whose title',' string is',/7x,'"',a,'".',/7x,' The code is only',' dimensioned for ',i3,' such lines. Reduce the',/7x,'number of such lines or increase the dimensioning',' parameter jpgxpa.')

                go to 990
            end if

            k1 = index(uheadx,'[')
            k2 = index(uheadx,']')
            k3 = index(uheadx,'(')
            k4 = index(uheadx,')')

            if (k1.le.0 .or. k2.le.0 .or. k2.lt.(k1 + 2) .or.    k3.le.0 .or. k4.le.0 .or. k4.lt.(k3 + 2)) then
                write (noutpt,1740) uheadx(1:j2)
                write (nttyo,1740) uheadx(1:j2)
1740 format(/' * Error - (EQ3NR/rd3ind) The following iopg',/7x,'option switch choice line read from the input file',/7x,"isn't in the required format:",/7x,'"',a,'"')

                go to 990
            end if

            ! Get the index of the option choice.
            ustr = uheadx(k3 + 1:k4 - 1)
            call chrint(ivar,nttyo,qrderr,ustr)

            if (qrderr) then
                go to 999
            end if

            ival = ivar

            do j = 1,jpgxpa
                ii = iopgox(j,n)

                if (ival .eq. ii) then
                    go to 700
                end if
            end do

            write (noutpt,1750) uheadx(1:j2)
            write (nttyo,1750) uheadx(1:j2)
1750 format(/' * Error - (EQ3NR/rd3ind) The iopg option switch',/7x,'choice line',/7x,'"',a,'"',/7x,'read from the input file references an out-of-range',' option',/7x,'choice index.')

            go to 990

700 continue

            ! Check the full option choice string.
            ustr = uopgox(j,n)
            j3 = ilnobl(ustr)
            j4 = ilnobl(uopgtx(n))

            if (index(uheadx(1:j2),ustr(1:j3)) .le. 0) then
                write (noutpt,1760) uheadx(1:j2),ustr(1:j3),uopgtx(n)(1:j4)
                write (nttyo,1760) uheadx(1:j2),ustr(1:j3),uopgtx(n)(1:j4)
1760 format(/' * Error - (EQ3NR/rd3ind) The iopg option switch',' choice line',/7x,'"',a,'"',/7x," read from the input file doesn't contain the",' matching defined string',/7x,'"',a,'".',/7x,'This line belongs to the option switch whose title',' string is',/7x,'"',a,'".')

                go to 990
            end if

            k1 = index(uheadx,'[')
            k2 = index(uheadx,']')
            ustr24 = uheadx(k1 + 1:k2 -1)
            call lejust(ustr24)
            j3 = ilnobl(ustr24)
            qmark = .false.

            if (j3 .gt. 0) then
                if (index(ustr24,'*') .ge. 1) then
                    qmark = .true.
                else if (index(ustr24,'x') .ge. 1) then
                    qmark = .true.
                else if (index(ustr24,'X') .ge. 1) then
                    qmark = .true.
                else
                    j4 = ilnobl(uopgtx(n))
                    write (noutpt,1770) ustr24(1:j3),uheadx(1:j2),uopgtx(n)(1:j4)
                    write (nttyo,1770) ustr24(1:j3),uheadx(1:j2),uopgtx(n)(1:j4)
1770 format(/" * Error - (EQ3NR/rd3ind) Don't recognize the",' string "',a,'"',/7x,'that appears on the iopg',' option switch choice line',/7x,'"',a,'"',/7x,'read from the input file. An option choice should',' be chosen by',/7x,'placing a "*", "x", or "X" in the',' checkbox ("[ ]"). This',/7x,'choice line belongs to',' the option switch whose title string is',/7x,'"',a,'".')

                    go to 990
                end if
            end if

            if (qmark) then
                nmark = nmark + 1
                iopg(n) = ival
                jlast = j
            end if
        end do

720 continue

        if (nmark .le. 0) then
            j2 = ilnobl(uopgtx(n))
            write (noutpt,1780) uopgtx(n)(1:j2)
            write (nttyo,1780) uopgtx(n)(1:j2)
1780 format(/' * Warning - (EQ3NR/rd3ind) No option choice was',' checked on the input file',/7x,'for the iopg option',' switch whose title string is',/7x,'"',a,'".')

            do j = 1,jpgxpa
                ival = iopgox(j,n)

                if (ival .eq. 0) then
                    go to 730
                end if
            end do

            write (noutpt,1790)
            write (nttyo,1790)
1790 format(/7x,'A default value of 0 has been applied, but no',/7x,'matching option choice string is defined.')

            go to 740

730 continue
            j3 = ilnobl(uopgox(j,n))
            write (noutpt,1792) uopgox(j,n)(1:j3)
            write (nttyo,1792) uopgox(j,n)(1:j3)
1792 format(/7x,'A default value of 0 has been applied. The',' matching string is',/7x,'"',a,'".')

740 continue
        end if

        if (nmark .gt. 1) then
            j2 = ilnobl(uopgtx(n))
            j = jlast
            j3 = ilnobl(uopgox(j,n))
            write (noutpt,1794) uopgtx(n)(1:j2),uopgox(j,n)(1:j3)
            write (nttyo,1794) uopgtx(n)(1:j2),uopgox(j,n)(1:j3)
1794 format(/' * Warning - (EQ3NR/rd3ind) More than one option',' choice was checked',/7x,'on the input file for the iopg',' option switch whose title string is',/7x,'"',a,'".',/7x,'The last choice checked will be used. The',' matching string is',/7x,'"',a,'".')
        end if
    end do

750 continue

    ! Iopr Print Option Switches.
    ! Note: iopr(1) = iopt1, etc.
    uheadx = 'Iopr Print Option Switches ("( 0)" marks default choices)'
    nfldtx = 1
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    write (noutpt,1016) uline1,uline2

    ! Loop on option switches.
    do nn = 1,noprmx
        ! Read the option title string from a one-line header.
        nfldtx = 1
        call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

        if (qrderr) then
            go to 999
        end if

        uheadx = ufield(1)
        j2 = ilnobl(uheadx)

        ! Check for the end of the block.
        if (uheadx(1:4) .ne. 'iopr') then
            ! Back up.
            backspace ninpts
            go to 850
        end if

        write (noutpt,1014) uline1
        k1 = index(uheadx,'(')
        k2 = index(uheadx,')')
        k3 = index(uheadx,'- ')

        if (k1.le.0 .or. k2.le.0 .or. k2.lt.(k1 + 2) .or. k3.le.0 .or.  k3.le.k2) then
            write (noutpt,1800) uheadx(1:j2)
            write (nttyo,1800) uheadx(1:j2)
1800 format(/' * Error - (EQ3NR/rd3ind) The iopr option switch',' title line',/7x,'"',a,'"',/7x,'read from the'," input file isn't in the required format.")

            go to 990
        end if

        ! Get the index of the option.
        ustr = uheadx(k1 + 1:k2 - 1)
        call chrint(ivar,nttyo,qrderr,ustr)

        if (qrderr) then
            go to 999
        end if

        n = ivar

        if (n.lt.1 .or. n.gt.nprxpa) then
            write (noutpt,1810) uheadx(1:j2),n
            write (nttyo,1810) uheadx(1:j2),n
1810 format(/' * Error - (EQ3NR/rd3ind) The iopr option switch',' title line',/7x,'"',a,'"',/7x,'read from the',' input file references an option switch index',/7x,'of ',i3,', which is out of range.')

            go to 990
        end if

        ! Check the full title string.
        ustr = uoprtx(n)
        j3 = ilnobl(ustr)

        if (index(uheadx(1:j2),ustr(1:j3)) .le. 0) then
            write (noutpt,1812) uheadx(1:j2),ustr(1:j3)
            write (nttyo,1812) uheadx(1:j2),ustr(1:j3)
1812 format(/' * Error - (EQ3NR/rd3ind) The iopr option switch',' title string',/7x,'"',a,'"',/7x,"read from the input file doesn't contain the",' matching defined string',/7x,'"',a,'".')

            go to 990
        end if

        iopr(n) = 0
        nmark = 0

        ! Loop on option choices.
        do jj = 1,jprxpa + 1
            ! Read a line. This contains an option choice line, else it
            ! is a separator line marking the end of the option choice
            ! lines for the current option.
            nfldtx = 1
            call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

            if (qrderr) then
                go to 999
            end if

            write (noutpt,1014) uline1

            uheadx = ufield(1)
            j2 = ilnobl(uheadx)

            if (uheadx(1:8) .eq. '--------') then
                go to 820
            end if

            if (jj .gt. jprxpa) then
                j3 = ilnobl(uoprtx(n))
                write (noutpt,1830) uoprtx(n)(1:j3),jprxpa
                write (nttyo,1830) uoprtx(n)(1:j3),jprxpa
1830 format(/' * Error - (EQ3NR/rd3ind) Have too many option',/7x,'choice lines for the iopr option switch whose title',' string is',/7x,'"',a,'".',/7x,' The code is only',' dimensioned for ',i3,' such lines. Reduce the',/7x,'number of such lines or increase the dimensioning',' parameter jprxpa.')

                go to 990
            end if

            k1 = index(uheadx,'[')
            k2 = index(uheadx,']')
            k3 = index(uheadx,'(')
            k4 = index(uheadx,')')

            if (k1.le.0 .or. k2.le.0 .or. k2.lt.(k1 + 2) .or.    k3.le.0 .or. k4.le.0 .or. k4.lt.(k3 + 2)) then
                write (noutpt,1840) uheadx(1:j2)
                write (nttyo,1840) uheadx(1:j2)
1840 format(/' * Error - (EQ3NR/rd3ind) The following iopr',/7x,'option switch choice line read from the input file',/7x,"isn't in the required format:",/7x,'"',a,'"')

                go to 990
            end if

            ! Get the index of the option choice.
            ustr = uheadx(k3 + 1:k4 - 1)
            call chrint(ivar,nttyo,qrderr,ustr)

            if (qrderr) then
                go to 999
            end if

            ival = ivar

            do j = 1,jprxpa
                ii = ioprox(j,n)

                if (ival .eq. ii) then
                    go to 800
                end if
            end do

            write (noutpt,1850) uheadx(1:j2)
            write (nttyo,1850) uheadx(1:j2)
1850 format(/' * Error - (EQ3NR/rd3ind) The iopr option switch',/7x,'choice line',/7x,'"',a,'"',/7x,'read from the input file references an out-of-range',' option',/7x,'choice index.')

            go to 990

800 continue

            ! Check the full option choice string.
            ustr = uoprox(j,n)
            j3 = ilnobl(ustr)
            j4 = ilnobl(uoprtx(n))

            if (index(uheadx(1:j2),ustr(1:j3)) .le. 0) then
                write (noutpt,1860) uheadx(1:j2),ustr(1:j3),uoprtx(n)(1:j4)
                write (nttyo,1860) uheadx(1:j2),ustr(1:j3),uoprtx(n)(1:j4)
1860 format(/' * Error - (EQ3NR/rd3ind) The iopr option switch',' choice line',/7x,'"',a,'"',/7x," read from the input file doesn't contain the",' matching defined string',/7x,'"',a,'".',/7x,'This line belongs to the option switch whose title',' string is',/7x,'"',a,'".')

                go to 990
            end if

            k1 = index(uheadx,'[')
            k2 = index(uheadx,']')
            ustr24 = uheadx(k1 + 1:k2 -1)
            call lejust(ustr24)
            j3 = ilnobl(ustr24)
            qmark = .false.

            if (j3 .gt. 0) then
                if (index(ustr24,'*') .ge. 1) then
                    qmark = .true.
                else if (index(ustr24,'x') .ge. 1) then
                    qmark = .true.
                else if (index(ustr24,'X') .ge. 1) then
                    qmark = .true.
                else
                    j4 = ilnobl(uoprtx(n))
                    write (noutpt,1870) ustr24(1:j3),uheadx(1:j2),uoprtx(n)(1:j4)
                    write (nttyo,1870) ustr24(1:j3),uheadx(1:j2),uoprtx(n)(1:j4)
1870 format(/" * Error - (EQ3NR/rd3ind) Don't recognize the",' string "',a,'"',/7x,'that appears on the iopr',' option switch choice line',/7x,'"',a,'"',/7x,'read from the input file. An option choice should',' be chosen by',/7x,'placing a "*", "x", or "X" in the',' checkbox ("[ ]"). This',/7x,'choice line belongs to',' the option switch whose title string is',/7x,'"',a,'".')

                    go to 990
                end if
            end if

            if (qmark) then
                nmark = nmark + 1
                iopr(n) = ival
                jlast = j
            end if
        end do

820 continue

        if (nmark .le. 0) then
            j2 = ilnobl(uoprtx(n))
            write (noutpt,1880) uoprtx(n)(1:j2)
            write (nttyo,1880) uoprtx(n)(1:j2)
1880 format(/' * Warning - (EQ3NR/rd3ind) No option choice was',' checked on the input file',/7x,'for the iopr option',' switch whose title string is',/7x,'"',a,'".')

            do j = 1,jprxpa
                ival = ioprox(j,n)

                if (ival .eq. 0) then
                    go to 830
                end if
            end do

            write (noutpt,1890)
            write (nttyo,1890)
1890 format(/7x,'A default value of 0 has been applied, but no',/7x,'matching option choice string is defined.')

            go to 840

830 continue
            j3 = ilnobl(uoprox(j,n))
            write (noutpt,1892) uoprox(j,n)(1:j3)
            write (nttyo,1892) uoprox(j,n)(1:j3)
1892 format(/7x,'A default value of 0 has been applied. The',' matching string is',/7x,'"',a,'".')

840 continue
        end if

        if (nmark .gt. 1) then
            j2 = ilnobl(uoprtx(n))
            j = jlast
            j3 = ilnobl(uoprox(j,n))
            write (noutpt,1894) uoprtx(n)(1:j2),uoprox(j,n)(1:j3)
            write (nttyo,1894) uoprtx(n)(1:j2),uoprox(j,n)(1:j3)
1894 format(/' * Warning - (EQ3NR/rd3ind) More than one option',' choice was checked',/7x,'on the input file for the iopr',' option switch whose title string is',/7x,'"',a,'".',/7x,'The last choice checked will be used. The',' matching string is',/7x,'"',a,'".')
        end if
    end do

850 continue

    ! Iodb debugging print option switches.
    ! Note: iodb(1) = iodb1, etc.
    uheadx = 'Iodb Debugging Print Option Switches ("( 0)" marks default choices)'
    nfldtx = 1
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    write (noutpt,1016) uline1,uline2

    ! Loop on option switches.
    do nn = 1,nodbmx
        ! Read the option title string from a one-line header.
        nfldtx = 1
        call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

        if (qrderr) then
            go to 999
        end if

        uheadx = ufield(1)
        j2 = ilnobl(uheadx)

        ! Check for the end of the block.
        if (uheadx(1:4) .ne. 'iodb') then
            ! Back up.
            backspace ninpts
            go to 950
        end if

        write (noutpt,1014) uline1
        k1 = index(uheadx,'(')
        k2 = index(uheadx,')')
        k3 = index(uheadx,'- ')

        if (k1.le.0 .or. k2.le.0 .or. k2.lt.(k1 + 2) .or. k3.le.0 .or.  k3.le.k2) then
            write (noutpt,1900) uheadx(1:j2)
            write (nttyo,1900) uheadx(1:j2)
1900 format(/' * Error - (EQ3NR/rd3ind) The iodb option switch',' title line',/7x,'"',a,'"',/7x,'read from the'," input file isn't in the required format.")

            go to 990
        end if

        ! Get the index of the option.
        ustr = uheadx(k1 + 1:k2 - 1)
        call chrint(ivar,nttyo,qrderr,ustr)

        if (qrderr) then
            go to 999
        end if

        n = ivar

        if (n.lt.1 .or. n.gt.ndbxpa) then
            write (noutpt,1910) uheadx(1:j2),n
            write (nttyo,1910) uheadx(1:j2),n
1910 format(/' * Error - (EQ3NR/rd3ind) The iodb option switch',' title line',/7x,'"',a,'"',/7x,'read from the',' input file references an option switch index',/7x,'of ',i3,', which is out of range.')

            go to 990
        end if

        ! Check the full title string.
        ustr = uodbtx(n)
        j3 = ilnobl(ustr)

        if (index(uheadx(1:j2),ustr(1:j3)) .le. 0) then
            write (noutpt,1920) uheadx(1:j2),ustr(1:j3)
            write (nttyo,1920) uheadx(1:j2),ustr(1:j3)
1920 format(/' * Error - (EQ3NR/rd3ind) The iodb option switch',' title string',/7x,'"',a,'"',/7x,"read from the input file doesn't contain the",' matching defined string',/7x,'"',a,'".')

            go to 990
        end if

        iodb(n) = 0
        nmark = 0

        ! Loop on option choices.
        do jj = 1,jdbxpa + 1
            ! Read a line. This contains an option choice line, else it
            ! is a separator line marking the end of the option choice
            ! lines for the current option.
            nfldtx = 1
            call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

            if (qrderr) then
                go to 999
            end if

            write (noutpt,1014) uline1

            uheadx = ufield(1)
            j2 = ilnobl(uheadx)

            if (uheadx(1:8) .eq. '--------') then
                go to 920
            end if

            if (jj .gt. jdbxpa) then
                j3 = ilnobl(uodbtx(n))
                write (noutpt,1922) uodbtx(n)(1:j3),jdbxpa
                write (nttyo,1922) uodbtx(n)(1:j3),jdbxpa
1922 format(/' * Error - (EQ3NR/rd3ind) Have too many option',/7x,'choice lines for the iodb option switch whose title',' string is',/7x,'"',a,'".',/7x,' The code is only',' dimensioned for ',i3,' such lines. Reduce the',/7x,'number of such lines or increase the dimensioning',' parameter jdbxpa.')

                go to 990
            end if

            k1 = index(uheadx,'[')
            k2 = index(uheadx,']')
            k3 = index(uheadx,'(')
            k4 = index(uheadx,')')

            if (k1.le.0 .or. k2.le.0 .or. k2.lt.(k1 + 2) .or.    k3.le.0 .or. k4.le.0 .or. k4.lt.(k3 + 2)) then
                write (noutpt,1940) uheadx(1:j2)
                write (nttyo,1940) uheadx(1:j2)
1940 format(/' * Error - (EQ3NR/rd3ind) The following iodb',/7x,'option switch choice line read from the input file',/7x,"isn't in the required format:",/7x,'"',a,'"')

                go to 990
            end if

            ! Get the index of the option choice.
            ustr = uheadx(k3 + 1:k4 - 1)
            call chrint(ivar,nttyo,qrderr,ustr)

            if (qrderr) then
                go to 999
            end if

            ival = ivar

            do j = 1,jdbxpa
                ii = iodbox(j,n)

                if (ival .eq. ii) then
                    go to 900
                end if
            end do

            write (noutpt,1950) uheadx(1:j2)
            write (nttyo,1950) uheadx(1:j2)
1950 format(/' * Error - (EQ3NR/rd3ind) The iodb option switch',/7x,'choice line',/7x,'"',a,'"',/7x,'read from the input file references an out-of-range',' option',/7x,'choice index.')

            go to 990

900 continue

            ! Check the full option choice string.
            ustr = uodbox(j,n)
            j3 = ilnobl(ustr)
            j4 = ilnobl(uodbtx(n))

            if (index(uheadx(1:j2),ustr(1:j3)) .le. 0) then
                write (noutpt,1960) uheadx(1:j2),ustr(1:j3),uodbtx(n)(1:j4)
                write (nttyo,1960) uheadx(1:j2),ustr(1:j3),uodbtx(n)(1:j4)
1960 format(/' * Error - (EQ3NR/rd3ind) The iodb option switch',' choice line',/7x,'"',a,'"',/7x," read from the input file doesn't contain the",' matching defined string',/7x,'"',a,'".',/7x,'This line belongs to the option switch whose title',' string is',/7x,'"',a,'".')

                go to 990
            end if

            k1 = index(uheadx,'[')
            k2 = index(uheadx,']')
            ustr24 = uheadx(k1 + 1:k2 -1)
            call lejust(ustr24)
            j3 = ilnobl(ustr24)
            qmark = .false.

            if (j3 .gt. 0) then
                if (index(ustr24,'*') .ge. 1) then
                    qmark = .true.
                else if (index(ustr24,'x') .ge. 1) then
                    qmark = .true.
                else if (index(ustr24,'X') .ge. 1) then
                    qmark = .true.
                else
                    j4 = ilnobl(uodbtx(n))
                    write (noutpt,1970) ustr24(1:j3),uheadx(1:j2),uodbtx(n)(1:j4)
                    write (nttyo,1970) ustr24(1:j3),uheadx(1:j2),uodbtx(n)(1:j4)
1970 format(/" * Error - (EQ3NR/rd3ind) Don't recognize the",' string "',a,'"',/7x,'that appears on the iodb',' option switch choice line',/7x,'"',a,'"',/7x,'read from the input file. An option choice should',' be chosen by',/7x,'placing a "*", "x", or "X" in the',' checkbox ("[ ]"). This',/7x,'choice line belongs to',' the option switch whose title string is',/7x,'"',a,'".')

                    go to 990
                end if
            end if

            if (qmark) then
                nmark = nmark + 1
                iodb(n) = ival
                jlast = j
            end if
        end do

920 continue

        if (nmark .le. 0) then
            j2 = ilnobl(uodbtx(n))
            write (noutpt,1980) uodbtx(n)(1:j2)
            write (nttyo,1980) uodbtx(n)(1:j2)
1980 format(/' * Warning - (EQ3NR/rd3ind) No option choice was',' checked on the input file',/7x,'for the iodb option',' switch whose title string is',/7x,'"',a,'".')

            do j = 1,jdbxpa
                ival = iodbox(j,n)

                if (ival .eq. 0) then
                    go to 930
                end if
            end do

            write (noutpt,1990)
            write (nttyo,1990)
1990 format(/7x,'A default value of 0 has been applied, but no',/7x,'matching option choice string is defined.')

            go to 940

930 continue
            j3 = ilnobl(uodbox(j,n))
            write (noutpt,1992) uodbox(j,n)(1:j3)
            write (nttyo,1992) uodbox(j,n)(1:j3)
1992 format(/7x,'A default value of 0 has been applied. The',' matching string is',/7x,'"',a,'".')

940 continue
        end if

        if (nmark .gt. 1) then
            j2 = ilnobl(uodbtx(n))
            j = jlast
            j3 = ilnobl(uodbox(j,n))
            write (noutpt,1994) uodbtx(n)(1:j2),uodbox(j,n)(1:j3)
            write (nttyo,1994) uodbtx(n)(1:j2),uodbox(j,n)(1:j3)
1994 format(/' * Warning - (EQ3NR/rd3ind) More than one option',' choice was checked',/7x,'on the input file for the iodb',' option switch whose title string is',/7x,'"',a,'".',/7x,'The last choice checked will be used. The',' matching string is',/7x,'"',a,'".')
        end if
    end do

950 continue

    ! Numerical parameters.
    ! Read the block title from a two-line header.
    uheadx = 'Numerical parameters'
    nfldtx = 1
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    write (noutpt,1016) uline1,uline2

    ! Read the beta convergence tolerance (tolbt) from a one-line
    ! header.
    uheadx = 'Beta convergence tolerance'
    nfldtx = 3
    call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)

    if (qrderr) then
        go to 999
    end if

    write (noutpt,1014) uline1
    ustr = ufield(2)
    call chreal(nttyo,qrderr,ustr,var)

    if (qrderr) then
        go to 999
    end if

    tolbt = var

    ! Read the del convergence tolerance (toldl) from a one-line header.
    uheadx = 'Del convergence tolerance'
    nfldtx = 3
    call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)

    if (qrderr) then
        go to 999
    end if

    write (noutpt,1014) uline1
    ustr = ufield(2)
    call chreal(nttyo,qrderr,ustr,var)

    if (qrderr) then
        go to 999
    end if

    toldl = var

    ! Read the maximum number of Newton-Raphson iterations (itermx)
    ! from a two-line header (the separator line is the end of the
    ! current block).
    uheadx = 'Max. Number of N-R Iterations'
    nfldtx = 3
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    write (noutpt,1016) uline1,uline2
    ustr = ufield(2)
    call chrint(ivar,nttyo,qrderr,ustr)

    if (qrderr) then
        go to 999
    end if

    itermx = ivar

    ! Ordinary basis switches.
    ! Read the block title from a two-line header.
    uheadx = 'Ordinary Basis Switches'
    nfldtx = 2
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    write (noutpt,1016) uline1,uline2

    nobsw = 0

    do n = 1,nbtmax + 1
        ! Read a line. If the block has not been completely read,
        ! this contains the name of a species to "Replace", and a
        ! sub-block for that species follows. Otherwise, this line is
        ! the first line of the next block.
        nfldtx = 0
        call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

        if (qrderr) then
            go to 999
        end if

        ustr = ufield(1)
        uheadx = 'Replace'
        call locase(ustr)
        call locase(uheadx)
        j2 = ilnobl(ustr)
        j3 = ilnobl(uheadx)

        if (ustr(1:j2) .ne. uheadx(1:j3)) then
            ! Back up.
            backspace ninpts
            go to 525
        end if

        write (noutpt,1014) uline1
        ustr = ufield(2)
        ustrn = ustr(1:24)
        call locase(ustrn)

        if (ustrn(1:5).eq.'none ' .or. ustr(1:1).eq.' ') then
            ustr = 'None'
        end if

        qnone = ustr(1:5).eq.'None '

        if (.not.qnone) then
            nobsw = nobsw + 1
            uobsw(1,nobsw) = ufield(2)(1:48)
        end if

        ! Read the name of the "with" species from a two-line header.
        uheadx = 'with'
        nfldtx = 3
        call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

        if (qrderr) then
            go to 999
        end if

        write (noutpt,1016) uline1,uline2

        if (.not.qnone) then
            uobsw(2,nobsw) = ufield(2)(1:48)
        end if
    end do

525 continue
    nobswt = nobsw

    ! Saturation flag tolerance (tolspf).
    tolspf = 0.

    ! Read the data from a two-line header.
    uheadx = 'Sat. flag tolerance'
    nfldtx = 3
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    write (noutpt,1016) uline1,uline2
    ustr = ufield(2)
    call chreal(nttyo,qrderr,ustr,var)

    if (qrderr) then
        go to 999
    end if

    tolspf = var

    ! Scale factor for the mass of aqueous solution to write
    ! on the pickup file.
    scamas = 0.

    ! Read the data from a two-line header.
    uheadx = 'Aq. Phase Scale Factor'
    nfldtx = 3
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    write (noutpt,1016) uline1,uline2
    ustr = ufield(2)
    call chreal(nttyo,qrderr,ustr,var)

    if (qrderr) then
        go to 999
    end if

    scamas = var

    ! Read a two-line header marking the end of the current problem
    ! input.
    uheadx = 'End of problem'
    nfldtx = 1
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    write (noutpt,1016) uline1,uline2

    write (noutpt,3050) nprob
    write (nttyo,3050) nprob
3050 format(/'   Done reading problem ',i3,'.',/)

    go to 999

990 continue
    qrderr = .true.

999 continue
end subroutine rd3ind
