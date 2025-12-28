program xcon3
    !! This is the main program of the XCON3 code. Configuration
    !! identification, the copyright statement, legal disclaimers,
    !! and similar statements are contained in XCON3/aaaxc3.f, the
    !! lead-off subroutine in the XCON3 source code. A short description
    !! of this program is also contained in that subroutine.
    implicit none

    include 'xcon3/x3op7.h'

    include 'eqlib/eqlj8.h'
    include 'eqlib/eqlo8.h'

    integer :: ietpar
    integer :: iktpar
    integer :: jetpar
    integer :: nbtpar
    integer :: netpar
    integer :: nodbpa
    integer :: nopgpa
    integer :: noprpa
    integer :: noptpa
    integer :: nsqpar
    integer :: ntitpa
    integer :: nxicpa
    integer :: nxmdpa
    integer :: nxtipa
    integer :: nxtpar

    parameter(ietpar = 10,iktpar = 20,jetpar = 3,nbtpar = 250)
    parameter(netpar = 12,ntitpa = 100,nxmdpa = 40,nxtpar = 50)
    parameter(nodbpa = 20,nopgpa = 20,noprpa = 20,noptpa = 20)

    parameter(nsqpar = nbtpar,nxtipa = nxtpar,nxicpa = iktpar*nxtpar)

    integer :: ietmax
    integer :: iktmax
    integer :: jetmax
    integer :: nbtmax
    integer :: netmax
    integer :: nodbmx
    integer :: nopgmx
    integer :: noprmx
    integer :: noptmx
    integer :: nsqmax
    integer :: ntitmx
    integer :: nxicmx
    integer :: nxmdmx
    integer :: nxtimx
    integer :: nxtmax

    integer :: newin
    integer :: ninpt
    integer :: ninpts
    integer :: noutpt
    integer :: nttyo

    integer :: iodb(nodbpa)
    integer :: iopg(nopgpa)
    integer :: iopr(noprpa)
    integer :: iopt(noptpa)
    integer :: jflagb(nsqpar)
    integer :: jflgi(nbtpar)
    integer :: jgext(netpar)
    integer :: jgexti(netpar)
    integer :: jxmod(nxmdpa)
    integer :: kxmod(nxmdpa)
    integer :: ncmpri(2,nxtipa)
    integer :: ncompb(nxtpar)
    integer :: ngexrt(jetpar,netpar)
    integer :: ngexti(jetpar,netpar)

    integer :: itermx
    integer :: jpres3
    integer :: net
    integer :: neti
    integer :: nprob
    integer :: ntitl
    integer :: nxcon
    integer :: nxmod
    integer :: nsq
    integer :: nxtb
    integer :: nbti
    integer :: nobswt
    integer :: nsbswt
    integer :: nxti

    integer :: i
    integer :: icount
    integer :: iebal3
    integer :: irdxc3
    integer :: itdsf3
    integer :: j
    integer :: jfli
    integer :: j2
    integer :: ik
    integer :: iktb
    integer :: k
    integer :: kl
    integer :: n
    integer :: nb
    integer :: nerr
    integer :: nmax
    integer :: nn
    integer :: nobsw
    integer :: nr1
    integer :: nr2
    integer :: ns
    integer :: ntest
    integer :: nx
    integer :: nxi
    integer :: nxic

    integer :: ilnobl

    logical :: qend
    logical :: qex
    logical :: qgexsh
    logical :: qrderr
    logical :: qxcon
    logical :: q8bchk
    logical :: q8beta

    character(len=80) :: utitl(ntitpa)
    character(len=56) :: ugexr(ietpar,jetpar,netpar)
    character(len=48) :: uxmod(nxmdpa)
    character(len=48) :: ucospi(nbtpar)
    character(len=48) :: uobsw(2,nbtpar)
    character(len=48) :: usbsw(2,nbtpar)
    character(len=48) :: uspeci(nbtpar)
    character(len=24) :: ugexmo(netpar)
    character(len=24) :: ugexp(netpar)
    character(len=24) :: ugexpi(netpar)
    character(len=24) :: ugexsi(ietpar,jetpar,netpar)
    character(len=24) :: usoli(nxtipa)
    character(len=24) :: umemi(nxicpa)
    character(len=24) :: ubasis(nsqpar)
    character(len=24) :: umemb(iktpar,nxtpar)
    character(len=24) :: uphas1(nsqpar)
    character(len=24) :: uphas2(nsqpar)
    character(len=24) :: usolb(nxtpar)
    character(len=24) :: uspecb(nsqpar)
    character(len=24) :: uxmd24(nxmdpa)
    character(len=8) :: ugexj(jetpar,netpar)
    character(len=8) :: ugexji(jetpar,netpar)
    character(len=8) :: uhfgex(ietpar,jetpar,netpar)
    character(len=8) :: uvfgex(ietpar,jetpar,netpar)
    character(len=8) :: uxkgex(ietpar,jetpar,netpar)

    character(len=80) :: uline
    character(len=24) :: uacion
    character(len=24) :: ublk24
    character(len=24) :: uebal
    character(len=24) :: uredox
    character(len=24) :: ustr
    character(len=24) :: ustr1
    character(len=24) :: ux24
    character(len=8) :: uplatc
    character(len=8) :: uplatm
    character(len=8) :: ustelu
    character(len=8) :: ustxc3
    character(len=8) :: uveelu
    character(len=8) :: uvexc3
    character(len=8) :: uv
    character(len=3) :: uoldv
    character(len=3) :: uoldvd
    character(len=3) :: unewv
    character(len=1) :: uoldf
    character(len=1) :: unewf
    character(len=1) :: ux1

    real(kind=8) :: apresh(5,2)
    real(kind=8) :: cspb(nsqpar)
    real(kind=8) :: xbarb(iktpar,nxtpar)
    real(kind=8) :: covali(nbtpar)
    real(kind=8) :: xbari(nxicpa)
    real(kind=8) :: cgexj(jetpar,netpar)
    real(kind=8) :: cgexpi(netpar)
    real(kind=8) :: egexsi(ietpar,jetpar,netpar)
    real(kind=8) :: mwtges(netpar)
    real(kind=8) :: tgexp(netpar)
    real(kind=8) :: xgexsi(ietpar,jetpar,netpar)
    real(kind=8) :: xhfgex(ietpar,jetpar,netpar)
    real(kind=8) :: xlkgex(ietpar,jetpar,netpar)
    real(kind=8) :: xvfgex(ietpar,jetpar,netpar)
    real(kind=8) :: xlkmod(nxmdpa)
    real(kind=8) :: zgexj(jetpar,netpar)

    real(kind=8) :: ehi
    real(kind=8) :: epstst
    real(kind=8) :: fep
    real(kind=8) :: fo2lgi
    real(kind=8) :: pei
    real(kind=8) :: presh
    real(kind=8) :: press
    real(kind=8) :: rho
    real(kind=8) :: scamas
    real(kind=8) :: tempc
    real(kind=8) :: tdspkg
    real(kind=8) :: tdspl
    real(kind=8) :: tolbt
    real(kind=8) :: toldl
    real(kind=8) :: tolsat
    real(kind=8) :: tolspf

    real(kind=8) :: verold
    real(kind=8) :: vernew

    ! BEGIN_MACHINE_DEPENDENT_CODE
    !   On some systems, a BLOCK DATA subroutine must be declared in an
    !   EXTERNAL statement to assure proper loading. On some other
    !   systems, this is not necessary, but neither it is not harmful.
    !   On yet some other systems, the EXTERNAL statement below may
    !   cause a problem. If so, try commenting it out. If you still
    !   have trouble, consult your local system documentation or
    !   experiment to find out how to correctly handle a BLOCK DATA
    !   subroutine on your system. The EXTERNAL statement below should
    !   not cause a problem if you are using a compiler which is fully
    !   compliant with the Fortran 90 standard. However, there is
    !   no guarantee that it will be adequate to assure correct loading
    !   of the BLOCK DATA subroutine.
    external bkdxc3

    ! END_MACHINE_DEPENDENT_CODE
    data noutpt /0/

    data ublk24 /'                        '/

    data apresh(1,1) / 1.013200000E+00/,apresh(2,1) / 0./,apresh(3,1) / 0./,apresh(4,1) / 0./,apresh(5,1) / 0./,apresh(1,2) /-4.345000000E-01/,apresh(2,2) / 7.632333333E-03/,apresh(3,2) / 5.514000000E-05/,apresh(4,2) /-1.263733333E-06/,apresh(5,2) / 1.396800000E-08/

    data epstst /1.e-12/

    ! BEGIN_MACHINE_DEPENDENT_CODE
    !   Define the console device number.
    !   The following works on UNIX, PC, and VAX machines.
    data nttyo  /6/

    !   BEGIN_MAC_DEPENDENT_CODE
    !     data nttyo  /9/
    !   END_MAC_DEPENDENT_CODE
    ! END_MACHINE_DEPENDENT_CODE
    data ninpt  /9/,ninpts /10/,newin /11/,nxcon /12/

    ! Get configuration identification data.
    call aaaxc3(ustxc3,uvexc3)
    call aaaelu(ustelu,uveelu)
    call platfd(uplatc,uplatm)

    ! Set dimensioning variables.
    ietmax = ietpar
    iktmax = iktpar
    jetmax = jetpar
    nbtmax = nbtpar
    netmax = netpar
    nodbmx = nodbpa
    nopgmx = nopgpa
    noprmx = noprpa
    noptmx = noptpa
    nsqmax = nsqpar
    ntitmx = ntitpa
    nxicmx = nxicpa
    nxmdmx = nxmdpa
    nxtimx = nxtipa
    nxtmax = nxtpar

    ! Does the INPUT file (the old input file) exist?
    inquire(file="input",exist=qex)

    if (.not.qex) then
        write (nttyo,1000)
1000 format(' * Error- (XCON3/xcon3) The INPUT file (the old input',' file).',/7x,"doesn't exist. Check the name that was",' specified.')

        go to 999
    end if

    ! Open the old input file.
    open(ninpt,file='input',status='old',err=700)
    go to 710

700 continue
    write (nttyo,1010)
1010 format(" * Error - (XCON3/xcon3) Can't open the INPUT file (the",' old input file).')

    go to 999

710 continue

    ! Determine the old input file format.
    !   W = compact
    !   D = menu-style
    k = 0
    kl = 0
    ntest = 10
    uoldf = 'W'

    do n = 1,100
        read (ninpt,1060) ux1
1060 format(a1)

        if (ux1(1:1) .ne. 'c') then
            kl = kl + 1
        end if

        if (ux1(1:1) .eq. '|') then
            k = k + 1
        end if

        if (kl .ge. ntest) then
            if (k .ge. ntest) then
                uoldf = 'D'
            end if

            go to 810
        end if
    end do

810 continue

    rewind(ninpt)

    ! Does the INPUTS file (the stripped input file) exist?
    ! If so, kill it. This file will be used to contain a copy of the
    ! old input file, stripped of comment lines.
    inquire(file="inputs",exist=qex)

    if (qex) then
        open(ninpts,file='inputs',status='old',err=720)
        close(ninpts,status='delete')
    end if

    ! Open the INPUTS file.
    open(ninpts,file='inputs',status='new',err=720)
    go to 730

720 continue
    write (nttyo,1020)
1020 format(" * Error - (XCON3/xcon3) Can't open the INPUTS file",/7x,'(the stripped old input file).')

    close(ninpt)
    go to 999

    ! Strip comment lines from the input file. Also strip any blank
    ! lines from a "D" format file.
730 continue
    if (uoldf(1:1) .eq. 'W') then
        do n = 1,10000
            read (ninpt,1030,end=750) uline
1030 format(a80)

            if (uline(1:1) .ne. '*') then
                write (ninpts,1030) uline
            end if
        end do

        write (nttyo,1040)
1040 format(" * Error - (XCON3/xcon3) The old input file is too  ' long.")

        close(ninpt)
        close(ninpts,status='delete')
        go to 999

750 continue
    else if (uoldf(1:1) .eq. 'D') then
        do n = 1,10000
            read (ninpt,1030,end=757) uline
            j2 = ilnobl(uline)

            if (j2 .gt. 0) then
                if (uline(1:1).ne.'*' .and. uline(1:1).ne.'c') then
                    write (ninpts,1030) uline
                end if
            end if
        end do

        write (nttyo,1040)
        close(ninpt)
        close(ninpts,status='delete')
        go to 999

757 continue
    end if

    ! Rewind the stripped input file.
    rewind(ninpts)

    ! Does the NEWIN file (to be the new input file) exist?
    ! If so, kill it.
    inquire(file="newin",exist=qex)

    if (qex) then
        open(newin,file='newin',status='old',err=760)
        close(newin,status='delete')
    end if

    ! Open the new input file.
    open(newin,file='newin',status='new',err=760)
    go to 770

760 continue
    write (nttyo,1050)
1050 format(" * Error - (XCON3/xcon3) Can't open the NEWIN file",' (to be the',/7x,'new input file.)')

    close(ninpt)
    close(ninpts,status='delete')
    go to 999

770 continue

    ! Input file formats:
    !   W = compact
    !   D = menu-style
    !    uoldf  = format of the old file
    !    unewf  = format of the new file
    unewf = ' '

    ! Input file version levels.
    !   '6.0' = version 6.0 (including version 6.1)
    !   '7.0' = version 7.0 (including versions 7.0x and 7.1)
    !   '7.2' = version 7.2 (including versions 7.2a and 7.2b)
    !   '8.0' = version 8.0
    !     uoldv  = version level of the old file
    !     uoldvd = default version level of the old file
    !     unewv  = version level of the new file
    uoldv = '   '
    uoldvd = '   '
    unewv = '   '

    ! Set ultimate default values.
    uoldvd = '7.0'
    unewv = '7.2'
    unewf = 'W'

    ! Does an options file exist?
    inquire(file="ixcon",exist=qxcon)

    if (qxcon) then
        open(nxcon,file='ixcon',status='old',err=780)
        go to 790

780 continue
        write (nttyo,1055)
1055 format(/" * Error - (XCON3/xcon3) Can't open the IXCON options",' file.')

        close(ninpt)
        close(ninpts,status='delete')
        close(newin)
        go to 999

        ! Read the IXCON options file.
790 continue
        call rddixc(nxcon,uoldvd,unewf,unewv)
    else
        write (nttyo,1057)
1057 format(/' * Error - (XCON3/xcon3) The IXCON options file'," doesn't exist.")

        close(ninpt)
        close(ninpts,status='delete')
        close(newin)
        go to 999
    end if

    ! Determine the old input file version level.
    !   '6.0' = version 6.0 (including version 6.1)
    !   '7.0' = version 7.0 (including versions 7.0x and 7.1)
    !   '7.2' = version 7.2 (including versions 7.2a and 7.2b)
    !   '8.0' = version 8.0
    i = 0

    do n = 1,ntitmx + 4
        read (ninpts,1030) uline

        if (uoldf(1:1) .eq. 'W') then
            if (uline(1:8) .eq. 'endit.  ') then
                go to 825
            end if
        else if (uoldf(1:1) .eq. 'D') then
            if (uline(1:8) .eq. '|-------') then
                i = i + 1
            end if

            if (i .ge. 3) then
                go to 825
            end if
        end if

        j = index(uline,'Version level=')

        if (j .le. 0) then
            j = index(uline,'version level=')
        end if

        if (j .gt. 0) then
            uv = uline(j + 14:j + 21)
            go to 840
        end if
    end do

825 continue

    i = 0
    rewind(ninpts)

    do n = 1,ntitmx
        read (ninpts,1030) uline

        if (uoldf(1:1) .eq. 'W') then
            if (uline(1:8) .eq. 'endit.  ') then
                go to 835
            end if
        else if (uoldf(1:1) .eq. 'D') then
            if (uline(1:8) .eq. '|-------') then
                i = i + 1
            end if

            if (i .ge. 3) then
                go to 835
            end if
        end if

        j = index(uline,'Version number=')

        if (j .le. 0) then
            j = index(uline,'version number=')
        end if

        if (j .gt. 0) then
            uv = uline(j + 15:j + 22)
            go to 840
        end if
    end do

835 continue

    write (nttyo,1070)
1070 format(/' * Warning - (XCON3/xcon3) The title on the old input',/7x,"file doesn't contain a version level tag. If this code is",/7x,'unable to determine the correct version level by other',/7x,'means, add a version tag by putting "Version level= X.X",',/7x,'where "X.X" is the version level, in the title, preferably',/7x,'on line 3. If there is more than problem on a single input',/7x,'file, the tag need be put only in the title of the first',' problem.')

    go to 880

840 continue
    call lejust(uv)

    if (uv(1:3) .eq. '6.0') then
        uoldv = '6.0'
    else if (uv(1:1) .eq. '6') then
        uoldv = '6.0'
    else if (uv(1:3) .eq. '7.0') then
        uoldv = '7.0'
    else if (uv(1:3) .eq. '7.1') then
        uoldv = '7.0'
    else if (uv(1:3) .eq. '7.2') then
        uoldv = '7.2'
    else if (uv(1:1) .eq. '7') then
        uoldv = '7.0'
    else if (uv(1:3) .eq. '8.0') then
        uoldv = '8.0'
    else if (uv(1:1) .eq. '8') then
        uoldv = '8.0'
    else
        j2 = ilnobl(uv)
        write (nttyo,1080) uv(1:j2)
1080 format(/" * Warning - (XCON3/xcon3) Can't determine the",/7x,'version level of the old input file. The version level',/7x,'is specified by a tag in the title of the first problem',/7x,'on the file as "',a,'". This is not a valid version',/7x,'level descriptor. Valid descriptors include "6.0", "7.0",',/7x,'"7.2", and "8.0". If the code is unable to determine the',/7x,'correct version level, replace the invalid descriptor.')
    end if

880 continue
    if (uoldv(1:3) .eq. '   ') then
        rewind(ninpts)

        do n = 1,1000
            read (ninpts,1030,end=838) uline
            j = index(uline,'uacion= ')

            if (j .gt. 0) then
                uoldv = '6.0'
                go to 838
            end if
        end do

838 continue
    end if

    ! If necessary, use the default for the version level of the old
    ! input file.
    if (uoldv(1:3) .eq. '   ') then
        uoldv = uoldvd
        write (nttyo,1090) uoldvd
1090 format(/' * Warning - (XCON3/xcon3) Taking the default',/7x,'value of "',a3,'" for the version level of the old',/7x,'input file.')
    end if

    ! Process the old input file, working from the stripped copy.
    rewind(ninpts)
    q8bchk = .false.
    q8beta = .false.
    nprob = 0
100 continue
    nprob = nprob + 1

    ! Zero or null various variables.
    call initcb(utitl,ntitmx)

    call initiz(iopt,noptmx)
    call initiz(iopg,nopgmx)
    call initiz(iopr,noprmx)
    call initiz(iodb,nodbmx)

    call initcb(uxmod,nxmdmx)
    call initcb(uxmd24,nxmdmx)
    call initiz(jxmod,nxmdmx)
    call initiz(kxmod,nxmdmx)
    call initaz(xlkmod,nxmdmx)

    call initcb(uspecb,nsqmax)
    call initcb(ubasis,nsqmax)
    call initiz(jflagb,nsqmax)
    call initaz(cspb,nsqmax)
    call initcb(uphas1,nsqmax)
    call initcb(uphas2,nsqmax)

    call initcb(usolb,nxtmax)

    nmax = iktmax*nxtmax
    call initcb(umemb,nmax)
    call initaz(xbarb,nmax)

    call initcb(uspeci,nbtmax)
    call initiz(jflgi,nbtmax)
    call initaz(covali,nbtmax)
    call initcb(ucospi,nbtmax)

    nmax = 2*nbtmax
    call initcb(usbsw,nmax)
    call initcb(uobsw,nmax)

    call initcb(usoli,nxtimx)

    call initcb(umemi,nxicmx)
    call initaz(xbari,nxicmx)

    nmax = 2*nxtimx
    call initiz(ncmpri,nmax)

    ! The uacion variable only appears on version level '6.0' input
    ! files. Provide a blank default value.
    uacion = ' '

    ! Read the current problem on the stripped input file.
    if (uoldf(1:1) .eq. 'W') then
        if (uoldv(1:3) .eq. '6.0') then
            call rd3w6(cspb,fep,iktmax,iodb,iopg,iopr,iopt,itermx,jflagb,jxmod,kxmod,ncompb,ninpts,nodbmx,nopgmx,noprmx,noptmx,nsq,nsqmax,ntitl,ntitmx,nttyo,nxmdmx,nxmod,nxtb,nxtmax,qend,qrderr,rho,tempc,tdspkg,tdspl,tolbt,toldl,tolsat,uacion,ubasis,uebal,umemb,uphas1,uphas2,uredox,usolb,uspecb,utitl,xbarb,uxmd24,xlkmod)

            if (qrderr) then
                go to 990
            end if
        else if (uoldv(1:3) .eq. '7.0') then
            call rd3w7(cspb,fep,iktmax,iodb,iopg,iopr,iopt,itermx,jflagb,jxmod,kxmod,ncompb,ninpts,nodbmx,nopgmx,noprmx,noptmx,nsq,nsqmax,ntitl,ntitmx,nttyo,nxmdmx,nxmod,nxtb,nxtmax,qend,qrderr,rho,tempc,tdspkg,tdspl,tolbt,toldl,tolsat,ubasis,uebal,umemb,uphas1,uphas2,uredox,usolb,uspecb,utitl,xbarb,uxmd24,xlkmod)

            if (qrderr) then
                go to 990
            end if
        else if (uoldv(1:3) .eq. '7.2') then
            ! Note: version level 7.2 is identical to version level 7.0
            ! for this format.
            call rd3w7(cspb,fep,iktmax,iodb,iopg,iopr,iopt,itermx,jflagb,jxmod,kxmod,ncompb,ninpts,nodbmx,nopgmx,noprmx,noptmx,nsq,nsqmax,ntitl,ntitmx,nttyo,nxmdmx,nxmod,nxtb,nxtmax,qend,qrderr,rho,tempc,tdspkg,tdspl,tolbt,toldl,tolsat,ubasis,uebal,umemb,uphas1,uphas2,uredox,usolb,uspecb,utitl,xbarb,uxmd24,xlkmod)

            if (qrderr) then
                go to 990
            end if
        else if (uoldv(1:3) .eq. '8.0') then
            call rd3w8(cgexj,cgexpi,covali,ehi,egexsi,fo2lgi,iebal3,ietmax,iodb,iopg,iopr,iopt,irdxc3,itdsf3,itermx,jgext,jetmax,jflgi,jgexti,jpres3,kxmod,mwtges,nbti,nbtmax,ncmpri,net,neti,netmax,ngexti,ninpts,ngexrt,nobswt,nodbmx,nopgmx,noprmx,noptmx,noutpt,nprob,nsbswt,ntitl,ntitmx,nttyo,nxicmx,nxmdmx,nxmod,nxti,nxtimx,pei,press,qend,qgexsh,qrderr,rho,scamas,tdspkg,tdspl,tempc,tgexp,tolbt,toldl,tolspf,ucospi,uebal,ugexj,ugexji,ugexmo,ugexp,ugexpi,ugexr,ugexsi,umemi,uobsw,uredox,usbsw,usoli,uspeci,utitl,uhfgex,uvfgex,uxkgex,uxmod,xbari,xgexsi,xhfgex,xlkgex,xvfgex,xlkmod,zgexj)

            if (qrderr) then
                go to 990
            end if
        else
            write (nttyo,1200) uoldv
1200 format(/' * Error - (XCON3/xcon3) Coding to implement',/7x,'reading an input file in "W" format has not been',/7x,'implemented for version level "',a3,'."')

            go to 990
        end if
    else if (uoldf(1:1) .eq. 'D') then
        if (uoldv(1:3).eq.'8.0' .and. .not.q8bchk) then
            ! Distinguish 8.0 from 8.0 beta. 8.0 will include the
            ! "Advisory:" and the "Option: on further procesing" strings,
            ! which 8.0 beta will not. Search for these only between
            ! the "Default redox constraint (irdxc3):" line and the
            ! "Alter/Suppress Options  | (nxmod)" line.
            icount = 0
120 continue
            read(ninpts,1030,end=140) uline
            j2 = index(uline(2:80),'|') - 1
            i = index(uline(2:j2),'Default redox constraint (irdxc3):')

            if (i .le. 0) then
                go to 120
            end if

130 continue
            read(ninpts,1030,end=140) uline
            j2 = index(uline(2:80),'|') - 1
            i = index(uline(2:j2),'Advisory:')

            if (i .gt. 0) then
                icount = icount + 1
            end if

            i = index(uline(2:j2),'Option: on further processing')

            if (i .gt. 0) then
                icount = icount + 1
            end if

            i = index(uline(2:j2),'Alter/Suppress Options  | (nxmod)')

            if (i .gt. 0) then
                go to 140
            end if

            if (icount .lt. 2) then
                go to 130
            end if

140 continue
            q8beta = icount.lt.2
            q8bchk = .true.
            rewind(ninpts)
        end if

        if (uoldv(1:3) .eq. '6.0') then
            write (nttyo,1205) uoldv
1205 format(/' * Error - (XCON3/xcon3) There is no "D" format',/7x,'for version level "6.0", hence',"can't read an input",/7x,'file in this format for this version level.')

            go to 990
        else if (uoldv(1:3) .eq. '7.0') then
            call rd3d7(cspb,fep,iktmax,iodb,iopg,iopr,iopt,itermx,jflagb,jxmod,kxmod,ncompb,ninpts,nodbmx,nopgmx,noprmx,noptmx,nsq,nsqmax,ntitl,ntitmx,nttyo,nxmdmx,nxmod,nxtb,nxtmax,qend,qrderr,rho,tempc,tdspkg,tdspl,tolbt,toldl,tolsat,ubasis,uebal,umemb,uphas1,uphas2,uredox,usolb,uspecb,utitl,xbarb,uxmd24,xlkmod)

            if (qrderr) then
                go to 990
            end if
        else if (uoldv(1:3) .eq. '7.2') then
            ! Note: XCON3/rd3d7.f can read a "D" format input file at
            ! version level '7.2' as well as at version level '7.0'.
            call rd3d7(cspb,fep,iktmax,iodb,iopg,iopr,iopt,itermx,jflagb,jxmod,kxmod,ncompb,ninpts,nodbmx,nopgmx,noprmx,noptmx,nsq,nsqmax,ntitl,ntitmx,nttyo,nxmdmx,nxmod,nxtb,nxtmax,qend,qrderr,rho,tempc,tdspkg,tdspl,tolbt,toldl,tolsat,ubasis,uebal,umemb,uphas1,uphas2,uredox,usolb,uspecb,utitl,xbarb,uxmd24,xlkmod)

            if (qrderr) then
                go to 990
            end if
        else if (uoldv(1:3).eq.'8.0' .and. q8beta) then
            call rd3d8b(cgexj,cgexpi,covali,ehi,egexsi,fo2lgi,iebal3,ietmax,iodb,iopg,iopr,iopt,irdxc3,itdsf3,itermx,jgext,jetmax,jflgi,jgexti,jpres3,kxmod,mwtges,nbti,nbtmax,ncmpri,net,neti,netmax,ngexti,ninpts,ngexrt,nobswt,nodbmx,nopgmx,noprmx,noptmx,noutpt,nprob,nsbswt,ntitl,ntitmx,nttyo,nxicmx,nxmdmx,nxmod,nxti,nxtimx,pei,press,qend,qrderr,rho,scamas,tdspkg,tdspl,tempc,tgexp,tolbt,toldl,tolspf,ucospi,uebal,ugexj,ugexji,ugexmo,ugexp,ugexpi,ugexr,ugexsi,umemi,uobsw,uredox,usbsw,usoli,uspeci,utitl,uhfgex,uvfgex,uxkgex,uxmod,xbari,xgexsi,xhfgex,xlkgex,xvfgex,xlkmod,zgexj)

            if (qrderr) then
                go to 990
            end if
        else if (uoldv(1:3) .eq. '8.0') then
            call rd3d8(cgexj,cgexpi,covali,ehi,egexsi,fo2lgi,iebal3,ietmax,iodb,iopg,iopr,iopt,irdxc3,itdsf3,itermx,jgext,jetmax,jflgi,jgexti,jpres3,kxmod,mwtges,nbti,nbtmax,ncmpri,net,neti,netmax,ngexti,ninpts,ngexrt,nobswt,nodbmx,nopgmx,noprmx,noptmx,noutpt,nprob,nsbswt,ntitl,ntitmx,nttyo,nxicmx,nxmdmx,nxmod,nxti,nxtimx,pei,press,qend,qgexsh,qrderr,rho,scamas,tdspkg,tdspl,tempc,tgexp,tolbt,toldl,tolspf,ucospi,uebal,ugexj,ugexji,ugexmo,ugexp,ugexpi,ugexr,ugexsi,umemi,uobsw,uredox,usbsw,usoli,uspeci,utitl,uhfgex,uvfgex,uxkgex,uxmod,xbari,xgexsi,xhfgex,xlkgex,xvfgex,xlkmod,zgexj)

            if (qrderr) then
                go to 990
            end if
        else
            write (nttyo,1210) uoldv
1210 format(/' * Error - (XCON3/xcon3) Coding to implement',/7x,'reading an input file in "D" format has not been',/7x,'implemented for version level "',a3,'."')

            go to 990
        end if
    else
        write (nttyo,1110) uoldf
1110 format(/' * Error - (XCON3/xcon3) Have unknown format',/7x,'specifier "',a1,'" for the old input file.')

        go to 990
    end if

    if (qend) then
        close(ninpt)
        close(ninpts,status='delete')
        close(newin)
        close(nxcon)
        go to 999
    end if

    ! Make changes in the content of the input file problem.
    ! Look for the string "Version level= " or Version number= " in the
    ! title. This identifies the version level of the input file.
    ! If present, change it. If not, insert a new line containing
    ! this at the beginning of the title.
    do n = 1,ntitl
        j = index(utitl(n),'Version level=')

        if (j .eq. 0) then
            j = index(utitl(n),'version level=')
        end if

        if (j .gt. 0) then
            utitl(n)(j + 14:j + 14) = ' '
            utitl(n)(j + 15:j + 17) = unewv
            utitl(n)(j + 18:80) = ' '
            go to 220
        end if

        j = index(utitl(n),'Version number=')

        if (j .eq. 0) then
            j = index(utitl(n),'version number=')
        end if

        if (j .gt. 0) then
            utitl(n)(j + 8:j + 14) = 'level= '
            utitl(n)(j + 15:j + 17) = unewv
            utitl(n)(j + 18:80) = ' '
            go to 220
        end if
    end do

    if ((ntitl + 1) .gt. ntitmx) then
        write (nttyo,1115) nprob,ntitmx
1115 format(/" * Error - (XCON3/xcon3) Can't add a version",/7x,'level marker to the first input file title of,'  /7x,'problem ',i2,' because this title already has the',/7x,'maximum length of ',i3,' lines.')

        go to 990
    end if

    do n = ntitl,1,-1
        utitl(n + 1) = utitl(n)
    end do

    ntitl = ntitl + 1

    utitl(1)(1:40) = '                                        '
    utitl(1)(41:80) = '                                        '
    utitl(1)(1:15)  = 'Version level= '
    utitl(1)(16:18) = unewv
    utitl(1)(21:27) = '(XCON3)'

220 continue

    ! Get version levels in floating point format.
    ! Calling sequence substitutions:
    !   uoldv for ustr
    !   verold for var
    call chreal(nttyo,qrderr,uoldv,verold)

    ! Calling sequence substitutions:
    !   unewv for ustr
    !   vernew for var
    call chreal(nttyo,qrderr,unewv,vernew)

    ! Patch up possible incompatibilities between "D" and "W" formats.
    ! The relevant variables here do not appear on EQ6 input files.
    if (vernew .lt. 8.0) then
        if (unewf(1:1) .eq. 'D') then
            if (uebal(1:5) .eq. 'none ') then
                uebal = ' '
            end if

            if (uebal(1:5) .eq. 'None ') then
                uebal = ' '
            end if

            if (uredox(1:5) .eq. 'none ') then
                uredox = ' '
            end if

            if (uredox(1:5) .eq. 'None ') then
                uredox = ' '
            end if
        else if (unewf(1:1) .eq. 'W') then
            if (uebal(1:5) .eq. 'None ') then
                uebal = 'none'
            end if

            if (uebal(1:1) .eq. ' ') then
                uebal = 'none'
            end if

            if (uredox(1:5) .eq. 'None ') then
                uredox = 'none'
            end if

            if (uredox(1:1) .eq. ' ') then
                uredox = 'none'
            end if
        end if
    else
        if (uebal(1:5) .eq. 'none ') then
            uebal = 'None'
        end if

        if (uebal(1:1) .eq. ' ') then
            uebal = 'None'
        end if

        if (uredox(1:5) .eq. 'none ') then
            uredox = 'None'
        end if

        if (uredox(1:1) .eq. ' ') then
            uredox = 'None'
        end if
    end if

    if (verold.lt.8.0 .and. vernew.ge.8.0) then
        ! Make translations from pre-Version 8.0 structures to Version 8.0
        ! and up structures. All pre-Version 8.0 basis switching is mapped
        ! to ordinary basis switching in Version 8.0 and up (none is
        ! mapped to special basis switching).
        nerr = 0

        jpres3 = 0

        if (tdspkg .gt. 0.) then
            itdsf3 = 0
        else if (tdspl .gt. 0.) then
            itdsf3 = 1
        else
            itdsf3 = 0
        end if

        j2 = ilnobl(uebal)

        if (j2 .eq. 0) then
            iebal3 = 0
        else if (uebal(1:j2) .eq. 'None') then
            iebal3 = 0
        else if (uebal(1:j2) .eq. 'none') then
            iebal3 = 0
        else if (uebal(1:j2) .ne. 'pick1.') then
            iebal3 = 1
        else
            write (nttyo,1608)
1608 format(/" * Error - (XCON3/xcon3) Can't translate this",/7x,'pre-version level 8.0 input file to version level',/7x,'8.0 or higher because it specifies the option that',/7x,'the code is to pick a species for electrical balancing',/7x,'(uebal= "pick1."). This option does not exist in',/7x,'version level 8.0 or higher. Change the uebal input',/7x,'on the old input file to blank or the name of a species',/7x,'to balance on.')

            nerr = nerr + 1
        end if

        scamas = 1.
        net = 0

        iopt(11) = iopt(2)
        iopt(2) = 0
        iopt(17) = iopt(3)
        iopt(3) = 0
        irdxc3 = iopt(1)
        iopt(1) = 0

        if (iopt(4) .ge. 2) then
            iopt(4) = 1
        end if

        iopr(10) = iopr(9)
        iopr(9) = iopr(6)
        iopr(6) = iopr(5)
        iopr(5) = 0
        iopr(3) = iopr(8)
        iopr(8) = 0

        iodb(3) = iodb(2)
        iodb(2) = 0

        nbti = nsq
        nsbswt = 0
        nobswt = 0

        if (irdxc3 .eq. -2) then
            pei = fep
        else if (irdxc3 .eq. -1) then
            ehi = fep
        else if (irdxc3 .eq. 0) then
            fo2lgi = fep
        end if

        j2 = ilnobl(uredox)

        if (j2 .le. 0) then
            uredox = 'None'
        end if

        j2 = ilnobl(uebal)

        if (j2 .le. 0) then
            uebal = 'None'
        end if

        do n = 1,nxmod
            uxmod(n)(1:24) = uxmd24(n)
        end do

        do ns = 1,nsq
            ! Basis species and basis switching.
            nb = ns
            uspeci(nb) = uspecb(ns)
            jflgi(nb) = jflagb(ns)
            covali(nb) = cspb(ns)

            if (uphas2(ns)(1:3) .eq. '   ') then
                ucospi(nb)(1:24) = uphas1(ns)
                ucospi(nb)(25:48) = uphas2(ns)
            else
                ucospi(nb)(1:24) = uphas2(ns)
                ucospi(nb)(25:48) = uphas1(ns)
            end if

            if (ubasis(ns)(1:3) .ne. '   ') then
                nobswt = nobswt + 1
                uobsw(1,nobswt) = uspecb(ns)
                uobsw(2,nobswt) = ubasis(ns)
            end if
        end do

        do nb = 1,nbti
            jfli = jflgi(nb)

            if (jfli.eq.19 .or. jfli.eq.20 .or. jfli.eq.21) then
                jflgi(nb) = 25
            end if
        end do

        do nb = 1,nbti
            jfli = jflgi(nb)

            if (jfli .eq. 31) then
                jflgi(nb) = 20
            end if

            if (jfli .eq. 32) then
                jflgi(nb) = 21
            end if
        end do

        do nb = 1,nbti
            jfli = jflgi(nb)

            if (jfli.ge.4 .and. jfli.le.8) then
                j2 = ilnobl(uspeci(nb)(1:24))
                write (nttyo,1610) jfli,uspeci(nb)(1:j2)
1610 format(/" * Error - (XCON3/xcon3) Can't translate this",/7x,'pre-version level 8.0 input file to version level',/7x,'8.0 or higher because it employs a jflag value of',/7x,i3,' for ',a,'. This jflag option for specifying',/7x,'a "free" concentration" is obsolete at the higher',/7x,'version level.')

                nerr = nerr + 1
            end if
        end do

        qgexsh = .false.

        nxti = nxtb
        nxic = 0

        do nx = 1,nxtb
            nxi = nx
            usoli(nxi) = usolb(nx)
            iktb = ncompb(nx)
            ncmpri(1,nxi) = nxic + 1
            ncmpri(2,nxi) = nxic + iktb

            do ik = 1,iktb
                nxic = nxic + 1
                umemi(nxic) = umemb(ik,nx)
                xbari(nxic) = xbarb(ik,nx)
            end do
        end do

        tolspf = tolsat

        if (nerr .gt. 0) then
            stop
        end if
    end if

    if ((abs(verold - 8.0).le.epstst .and. q8beta) .and. vernew.ge.8.0) then
        ! Make additions if converting from Version 8.0 beta to
        ! Version 8.0 and up.
        qgexsh = .false.
    end if

    if (verold.ge.8.0 .and. vernew.lt.8.0) then
        ! Make translations from Version 8.0 and up structures to
        ! pre-Version 8.0 structures, and vice versa. Special basis
        ! switching in Version 8.0 and up can't be mapped to any
        ! pre-Version 8.0 structure.
        nerr = 0

        if (tempc .le. 100.) then
            presh = apresh(1,1)
        else
            presh = 0.

            do nn = 1,5
                n = 6 - nn
                presh = apresh(n,2) + tempc*presh
            end do
        end if

        if (jpres3 .le. 0) then
            press = 0.
        else if (jpres3 .eq. 1) then
            write (nttyo,1692)
1692 format(/' * Note - (XCON3/xcon3) This version level 8.0',/7x,'or higher input file specifies the option to compute',/7x,'the pressure from the temperature according to the',/7x,'1.013-bar/steam-saturation curve (jpres3 = 1). This',/7x,'option does not directly translate to pre-Version 8.0.',/7x,'level The pressure at that level implicitly lies on the',/7x,'data file reference pressure curve, which may or may',/7x,'not coincide with the 1.013-bar/steam-saturation curve.',/7x,'To avoid getting this message, change the option on the',/7x,'old input file to specify computing the pressure',/7x,'according to the data file reference pressure curve',/7x,'(jpres3 = 0).')

            press = 0.
        else if (jpres3 .eq. 2) then
            if (press.gt.0 .and. abs(press - presh).gt.1.e-4) then
                write (nttyo,1700) press,presh
1700 format(/" * Error - (XCON3/xcon3) Can't translate this",/7x,'version level 8.0 or higher input file to a lower',/7x,'version level because it specifies a pressure of',/7x,f9.4,' bars, which differs from the pre-Version 8',/7x,'standard grid pressure of ',f9.4,' bars. To convert',/7x,'this input file, set press = 0., which defaults to',/7x,'the standard grid pressure.')

                nerr = nerr + 1
            end if
        else
            write (nttyo,1702) jpres3
1702 format(/" * Error - (XCON3/xcon3) Can't translate this",/7x,'version level 8.0 or higher input file to a lower',/7x,'version level because it specifies an unknown',/7x,'pressure option (jpres3= ',i2,'.')

            nerr = nerr + 1
        end if

        if (itdsf3 .le. 0) then
            continue
        else if (itdsf3 .ge. 1) then
            continue
        end if

        if (iebal3 .le. 0) then
            uebal = ' '
        else
            continue
        end if

        if (abs(scamas - 1.0) .gt. 1.e-6) then
            write (nttyo,1710) scamas
1710 format(/" * Error - (XCON3/xcon3) Can't translate this",/7x,'version level 8.0 or higher input file to a lower',/7x,'version level because it specifies a  scale factor',/7x,'of ',1pe11.4,' for the mass of aqueous solution to',/7x,'write on the PICKUP file. Pre-version 8 input files',/7x,'lack this parameter, effectively assuming a value of',' 1.0.',/7x,'To convert this input file, set scamas = 1.0.')

            nerr = nerr + 1
        end if

        if (net .gt. 0) then
            write (nttyo,1720) net
1720 format(/" * Error - (XCON3/xcon3) Can't translate this",/7x,'version level 8.0 or higher input file to a lower',/7x,'version level because it defines ',i2,' generic',/7x,'ion exchangers.Pre-version 8 input files lack the',/7x,'generic ion exchange capability.')

            nerr = nerr + 1
        end if

        iopt(2) = iopt(11)
        iopt(11) = 0
        iopt(3) = iopt(17)
        iopt(17) = 0
        iopt(1) = irdxc3

        if (iopt(4).ge.1 .and. nxti.gt.0) then
            iopt(4) = 2
        end if

        if (iopt(19) .gt. 0) then
            write (nttyo,1730)
1730 format(/" * Warning - (XCON3/xcon3) Can't translate the",/7x,'"Advanced EQ3NR PICKUP File Options" (iopt(19)) option',/7x,'specified on this version level 8.0 or higher input',/7x,'file to a lower version level because this option does',/7x,'not exist at the lower version level.')
        end if

        iopr(5) = iopr(6)
        iopr(6) = iopr(9)
        iopr(9) = iopr(10)
        iopr(10) = 0
        iopr(8) = iopr(3)
        iopr(3) = 0

        iodb(2) = iodb(3)
        iodb(3) = 0

        if (iopt(1) .eq. -2) then
            fep = pei
        else if (iopt(1) .eq. -1) then
            fep = ehi
        else if (iopt(1) .eq. 0) then
            fep = fo2lgi
        end if

        ! XX
        do nb = 1,nbti
            if (jflgi(nb) .eq. 22) then
                ! Translate the pmH option to a pH option.
                ! This is done at the version 8 level. Additional
                ! translation may follow.
                jflgi(nb) = 20
                iopg(2) = 1
            end if
        end do

        do nb = 1,nbti
            if (jflgi(nb) .eq. 23) then
                ! The pmX option is not translatable.
                j2 = ilnobl(uspeci(nb))
                write (nttyo,1732) uspeci(nb)(1:j2)
1732 format(/" * Error - (XCON3/xcon3) Can't translate this",/7x,'version level 8.0 or higher input file to a lower',/7x,'version level because it contains a pmX input for',/7x,a,'. The pmX option does not exist at a lower',/7x,'version level.')

                nerr = nerr + 1
            end if
        end do

        ! XX
        nsq = nbti

        do nb = 1,nbti
            ! Basis species.
            ns = nb
            uspecb(ns) = uspeci(nb)(1:24)
            jflagb(ns) = jflgi(nb)
            cspb(ns) = covali(nb)

            if (ucospi(nb)(25:27) .eq. '   ') then
                uphas1(ns) = ucospi(nb)(1:24)
                uphas2(ns) = ' '
            else if (ucospi(nb)(25:33) .eq. 'Aqueous solution ') then
                uphas1(ns) = ucospi(nb)(1:24)
                uphas2(ns) = ' '
            else if (ucospi(nb)(25:28) .eq. 'Gas ') then
                uphas1(ns) = ucospi(nb)(1:24)
                uphas2(ns) = ' '
            else
                uphas2(ns) = ucospi(nb)(1:24)
                uphas1(ns) = ucospi(nb)(25:48)
            end if
        end do

        do ns = 1,nsq
            if (jflagb(ns) .eq. 19) then
                jflagb(ns) = 16
                cspb(ns) = -cspb(ns)
            end if
        end do

        if (unewf(1:1) .eq. 'W') then
            do ns = 1,nsq
                if (jflagb(ns) .eq. 20) then
                    jflagb(ns) = 16
                    cspb(ns) = -cspb(ns)
                else if (jflagb(ns) .eq. 21) then
                    uphas1(ns) = 'Cl-'
                    jflagb(ns) = 17
                    cspb(ns) = -cspb(ns)
                end if
            end do
        else
            do ns = 1,nsq
                if (jflagb(ns) .eq. 20) then
                    jflagb(ns) = 31
                else if (jflagb(ns) .eq. 21) then
                    jflagb(ns) = 32
                end if
            end do
        end if

        do ns = 1,nsq
            if (jflagb(ns) .eq. 25) then
                jflagb(ns) = 19

                if (uphas2(ns)(1:24).ne.ublk24(1:24) .and.      uphas2(ns)(1:24).ne.uphas1(ns)(1:24)) then
                    jflagb(ns) = 20
                end if

                j = index(uphas1(ns),'(g)')

                if (j .gt. 0) then
                    jflagb(ns) = 21
                end if
            end if
        end do

        if (nsbswt .gt. 0) then
            write (nttyo,1740) nsbswt
1740 format(/" * Error - (XCON3/xcon3) Can't translate this",/7x,'version level 8.0 or higher input file to a lower',/7x,'version level because it contains ',i3,' directives',/7x,'for special basis switching.')

            nerr = nerr + 1
        end if

        do nobsw = 1,nobswt
            ! Basis switching.
            do nb = 1,nbti
                if (uobsw(1,nobsw)(1:48) .eq. uspeci(nb)(1:48)) then
                    ns = nb
                    ubasis(ns) = uobsw(2,nobsw)
                    go to 244
                end if
            end do

244 continue
        end do

        do n = 1,nxmod
            uxmd24(n) = uxmod(n)(1:24)

            ! Set the jxmod value. This is a flag for the type of species:
            !   0 = aqueous species
            !   1 = pure mineral
            !   2 = gas
            !   3 = solid solution
            ! The jxmod array is not used in Version 8.0 and up. It is
            ! not possible to design a perfect logic for determining the
            ! correct value. Here there may be a problem in determining
            ! whether a species is a pure mineral or a solid solution.
            ! Also, in some cases, the Version 8.0 and up option may not
            ! map back to the pre-Version 8.0 format.
            if (uxmod(n)(25:40) .eq. 'Aqueous solution') then
                jxmod(n) = 0
                go to 250
            end if

            if (uxmod(n)(25:28) .eq. 'Gas ') then
                jxmod(n) = 2
                go to 250
            end if

            if (uxmod(n)(25:48).ne.uxmod(n)(1:24) .and.      uxmod(n)(25:27).ne.'   ') then
                nerr = nerr + 1
                j2 = ilnobl(uxmod(n))
                write (nttyo,1117) uxmod(n)(1:j2)
1117 format(/" * Error - (XCON3/xcon3) Can't translate a",' suppress/alter option',/7x,'for "',a,'"',/7x,'to pre-Version 8.0 format.')

                go to 250
            end if

            i = index(uxmod(n),'(aq)')

            if (i .gt. 0) then
                jxmod(n) = 0
                go to 250
            end if

            ux24 = uxmod(n)(1:24)
            j2 = ilnobl(ux24)

            if (j2 .gt. 0) then
                if (ux24(j2:j2).eq.'+' .or. ux24(j2:j2).eq.'-') then
                    jxmod(n) = 0
                    go to 250
                end if
            end if

            i = index(uxmod(n),'(g)')

            if (i .le. 0) then
                i = index(uxmod(n),'Hydrogen')
            end if

            if (i .le. 0) then
                i = index(uxmod(n),'Oxygen')
            end if

            if (i .le. 0) then
                i = index(uxmod(n),'Nitrogen')
            end if

            if (i .le. 0) then
                i = index(uxmod(n),'Steam')
            end if

            if (i .gt. 0) then
                jxmod(n) = 2
                go to 250
            end if

            i = index(uxmod(n),'(ss)')

            if (i .le. 0) then
                i = index(uxmod(n),'-ss')
            end if

            if (i .le. 0) then
                i = index(uxmod(n),'Carbonate-Calcite')
            end if

            if (i .le. 0) then
                i = index(uxmod(n),'Biotite')
            end if

            if (i .le. 0) then
                i = index(uxmod(n),'Olivine')
            end if

            if (i .le. 0) then
                i = index(uxmod(n),'Plagioclase')
            end if

            if (i .le. 0) then
                i = index(uxmod(n),'Orthopyroxene')
            end if

            if (i .le. 0) then
                i = index(uxmod(n),'Smectite-di')
            end if

            if (i .le. 0) then
                i = index(uxmod(n),'Saponite-tri')
            end if

            if (i .gt. 0) then
                jxmod(n) = 3
                go to 250
            end if

            i = index(uxmod(n),'ite ')

            if (i .le. 0) then
                i = index(uxmod(n),'ime ')
            end if

            if (i .le. 0) then
                i = index(uxmod(n),'ine ')
            end if

            if (i .le. 0) then
                i = index(uxmod(n),'Quartz')
            end if

            if (i .gt. 0) then
                jxmod(n) = 1
                go to 250
            end if

            j2 = ilnobl(uxmod(n))
            write (nttyo,1119) uxmod(n)(1:j2)
1119 format(/" * Warning - (XCON3/xcon3) Can't unambiguously",/7x,'determine what kind of species in an alter/suppress',/7x,'option is "',a,'".',/7x,'Setting jxmod to 1. Check to see that this is correct',/7x,'(0 = aqueous species, 1 = pure mineral,',' 3 = gas species',/7x,'3 = solid solution).')

            jxmod(n) = 1
250 continue
        end do

        nxtb = nxti

        do nxi = 1,nxti
            nx = nxi
            usolb(nx) = usoli(nxi)
            nr1 = ncmpri(1,nxi)
            nr2 = ncmpri(2,nxi)
            ncompb(nx) =  nr2 - nr1 + 1
            ik = 0

            do nxic = nr1,nr2
                ik = ik + 1
                umemb(ik,nx) = umemi(nxic)
                xbarb(ik,nx) = xbari(nxic)
            end do
        end do

        tolsat = tolspf

        if (nerr .gt. 0) then
            stop
        end if
    end if

    if (vernew.lt.8.0 .and. unewf(1:1).eq."D") then
        ! Map any inputs equivalent to pH to pH, and any inputs equivalent
        ! to pHCl to pHCl, but only if the new input file is to be in "D"
        ! format. Use the hidden jflag values 31 (pH) and 32 (pHCl), which
        ! do not exist in the case of "W" format.
        do ns = 1,nsq
            ustr = uspecb(ns)

            if (ustr(1:3).eq.'h+ ' .or. ustr(1:3).eq.'H+ ') then
                if (jflagb(ns) .eq. 16) then
                    jflagb(ns) = 31
                    cspb(ns) = -cspb(ns)
                else if (jflagb(ns) .eq. 17) then
                    ustr1 = uphas1(ns)

                    if (ustr1(1:4).eq.'cl- ' .or. ustr1(1:4).eq.'Cl- ')          then
                        uphas1(ns) = ' '
                        jflagb(ns) = 32
                        cspb(ns) = -cspb(ns)
                    end if
                end if
            end if
        end do
    end if

    if (vernew .ge. 8.0) then
        ! Map any inputs equivalent to pH to pH, and any inputs equivalent
        ! to pHCl to pHCl.
        do nb = 1,nbti
            ustr = uspeci(nb)

            if (ustr(1:3) .eq. 'H+ ') then
                jfli = jflgi(nb)

                if (jfli .eq. 16) then
                    jflgi(nb) = 20
                    covali(nb) = -covali(nb)
                else if (jfli .eq. 19) then
                    jflgi(nb) = 20
                else if (jfli .eq. 17) then
                    ustr1 = ucospi(nb)

                    if (ustr1(1:4).eq.'Cl- ') then
                        ucospi(nb) = ' '
                        jflgi(nb) = 21
                        covali(nb) = -covali(nb)
                    end if
                end if
            end if
        end do
    end if

    ! Write the current problem on the new input file.
    if (unewf(1:1) .eq. 'W') then
        if (unewv(1:3) .eq. '6.0') then
            call wr3w6(cspb,fep,iktmax,iodb,iopg,iopr,iopt,itermx,jflagb,jxmod,kxmod,ncompb,newin,nodbmx,nopgmx,noprmx,noptmx,nsq,nsqmax,ntitl,ntitmx,nxmdmx,nxmod,nxtb,nxtmax,rho,tempc,tdspkg,tdspl,tolbt,toldl,tolsat,uacion,ubasis,uebal,umemb,uphas1,uphas2,uredox,usolb,uspecb,utitl,xbarb,uxmd24,xlkmod)
        else if (unewv(1:3) .eq. '7.0') then
            call wr3w7(cspb,fep,iktmax,iodb,iopg,iopr,iopt,itermx,jflagb,jxmod,kxmod,ncompb,newin,nodbmx,nopgmx,noprmx,noptmx,nsq,nsqmax,ntitl,ntitmx,nxmdmx,nxmod,nxtb,nxtmax,rho,tempc,tdspkg,tdspl,tolbt,toldl,tolsat,ubasis,uebal,umemb,uphas1,uphas2,uredox,usolb,uspecb,utitl,xbarb,uxmd24,xlkmod)
        else if (unewv(1:3) .eq. '7.2') then
            ! Note: version level 7.2 is identical to version level 7.0
            ! for this format.
            call wr3w7(cspb,fep,iktmax,iodb,iopg,iopr,iopt,itermx,jflagb,jxmod,kxmod,ncompb,newin,nodbmx,nopgmx,noprmx,noptmx,nsq,nsqmax,ntitl,ntitmx,nxmdmx,nxmod,nxtb,nxtmax,rho,tempc,tdspkg,tdspl,tolbt,toldl,tolsat,ubasis,uebal,umemb,uphas1,uphas2,uredox,usolb,uspecb,utitl,xbarb,uxmd24,xlkmod)
        else if (unewv(1:3) .eq. '8.0') then
            call wr3w8(cgexj,cgexpi,covali,ehi,egexsi,fo2lgi,iebal3,ietmax,iodb,iopg,iopr,iopt,irdxc3,itdsf3,itermx,jgext,jetmax,jflgi,jgexti,jpres3,kxmod,mwtges,nbti,nbtmax,ncmpri,net,neti,netmax,newin,ngexti,ngexrt,nobswt,nodbmx,nopgmx,noprmx,noptmx,nsbswt,ntitl,ntitmx,nttyo,nxicmx,nxmdmx,nxmod,nxti,nxtimx,pei,press,qgexsh,rho,scamas,tdspkg,tdspl,tempc,tgexp,tolbt,toldl,tolspf,ucospi,uebal,ugexj,ugexji,ugexmo,ugexp,ugexpi,ugexr,ugexsi,umemi,uobsw,uredox,usbsw,usoli,uspeci,utitl,uhfgex,uvfgex,uxkgex,uxmod,xbari,xgexsi,xhfgex,xlkgex,xvfgex,xlkmod,zgexj)
        else
            write (nttyo,1220) unewv
1220 format(/' * Error - (XCON3/xcon3) Coding to implement',/7x,'writing an input file in "W" format has not been',/7x,'implemented for version level "',a3,'."')

            go to 990
        end if
    else if (unewf(1:1) .eq. 'D') then
        if (unewv(1:3) .eq. '6.0') then
            write (nttyo,1230) unewv
1230 format(/' * Error - (XCON3/xcon3) There is no "D" format',/7x,'for version level "6.0", hence',"can't write an input",/7x,'file in this format for this version level.')

            go to 990
        else if (unewv(1:3) .eq. '7.0') then
            call wr3d7(cspb,fep,iktmax,iodb,iopg,iopr,iopt,itermx,jflagb,jxmod,kxmod,ncompb,newin,nodbmx,nopgmx,noprmx,noptmx,nsq,nsqmax,ntitl,ntitmx,nttyo,nxmdmx,nxmod,nxtb,nxtmax,rho,tempc,tdspkg,tdspl,tolbt,toldl,tolsat,ubasis,uebal,umemb,uphas1,uphas2,uredox,usolb,uspecb,utitl,xbarb,uxmd24,xlkmod)
        else if (unewv(1:3) .eq. '7.2') then
            call wr3d72(cspb,fep,iktmax,iodb,iopg,iopr,iopt,itermx,jflagb,jxmod,kxmod,ncompb,newin,nodbmx,nopgmx,noprmx,noptmx,nsq,nsqmax,ntitl,ntitmx,nttyo,nxmdmx,nxmod,nxtb,nxtmax,rho,tempc,tdspkg,tdspl,tolbt,toldl,tolsat,ubasis,uebal,umemb,uphas1,uphas2,uredox,usolb,uspecb,utitl,xbarb,uxmd24,xlkmod)
        else if (unewv(1:3) .eq. '8.0') then
            call wr3d8(cgexj,cgexpi,covali,ehi,egexsi,fo2lgi,iebal3,ietmax,iodb,iopg,iopr,iopt,irdxc3,itdsf3,itermx,jgext,jetmax,jflgi,jgexti,jpres3,kxmod,mwtges,nbti,nbtmax,ncmpri,net,neti,netmax,newin,ngexti,ngexrt,nobswt,nodbmx,nopgmx,noprmx,noptmx,nsbswt,ntitl,ntitmx,nttyo,nxicmx,nxmdmx,nxmod,nxti,nxtimx,pei,press,qgexsh,rho,scamas,tdspkg,tdspl,tempc,tgexp,tolbt,toldl,tolspf,ucospi,uebal,ugexj,ugexji,ugexmo,ugexp,ugexpi,ugexr,ugexsi,umemi,uobsw,uredox,usbsw,usoli,uspeci,utitl,uhfgex,uvfgex,uxkgex,uxmod,xbari,xgexsi,xhfgex,xlkgex,xvfgex,xlkmod,zgexj)
        else
            write (nttyo,1240) unewv
1240 format(/' * Error - (XCON3/xcon3) Coding to implement',/7x,'writing an input file in "D" format has not been',/7x,'implemented for version level "',a3,'."')

            go to 990
        end if
    else
        write (nttyo,1130) unewf
1130 format(/' * Error - (XCON3/xcon3) Have unknown format',/7x,'specifier "',a1,'" for the new input file.')

        go to 990
    end if

    ! Go back to do another problem on the same input file.
    go to 100

990 continue
    write (nttyo,1140) uoldv
1140 format(/' * Error - (XCON3/xcon3) Have encountered a read',/7x,'error while reading the old input file. First check the',/7x,'version level. The version level was taken to be "',a3,'".',/7x,'If the version level is marked incorrectly in the primary',/7x,'title of the first problem, correct this, as the default',/7x,'value read from the IXCON file or set in the code itself',/7x,'does not override this. If the problem is not due to an',/7x,'error in the assumed version level, there is probably a',/7x,'format error on the old input file. Check the stripped',/7x,'copy of the old input file (INPUTS).')

    ! Close all files. Don't delete the inputs file. Delete the NEWIN
    ! file, which is not valid.
    close(ninpt)
    close(ninpts)
    close(newin,status='delete')
    close(nxcon)

999 continue
end program xcon3
