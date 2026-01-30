subroutine rd6d8(awmaxi,awmini,cbsri,cdac,cesri,cgexj,csigma,dlaplo,dlaprn,dleplo,dleprn,dlhplo,dlhprn,dloplo,dloprn,dltpll,dltplo,dltprl,dltprn,dlxdmp,dlxmx0,dlxpll,dlxplo,dlxprl,dlxprn,eact,egersi,ehmaxi,ehmini,electr,fkrc,iact,ibsrti,iesrti,ietmax,iktmax,imchmx,imech,iodb,iopg,iopr,iopt,igerti,itermx,ixrti,jcode,jgerti,jetmax,jflgi,jgext,jpress,jreac,jtemp,kbt,kct,kdim,kmax,kmt,kprs,ksplmx,ksppmx,kstpmx,kxmod,kxt,hact,modr,moffg,morr,mprphi,mprspi,mtbaqi,mtbi,mwtges,nbti,nbtmax,nbt1mx,nctmax,ndact,ndctmx,nert,nertmx,net,netmax,nffg,nffgmx,ngexrt,ninpts,nobswt,nodbmx,nopgmx,noprmx,noptmx,nordmx,noutpt,nprob,nprpmx,nprpti,nprsmx,nprsti,nptkmx,nrct,nrctmx,nrk,nsbswt,nsk,nsrt,nsrtmx,ntitl1,ntitl2,ntitmx,ntrymx,nttkmx,nttyo,nxmdmx,nxmod,nxopex,nxopmx,nxopt,nxpemx,nxrt,nxrtmx,o2maxi,o2mini,phmaxi,phmini,pressb,pressi,ptk,qend,qgexsh,qrderr,rkb,rxbari,sfcar,ssfcar,tempcb,tempci,tgexp,timmxi,tistti,tolbt,toldl,tolsat,tolxsf,trkb,ttk,ubmtbi,ubsri,ucxri,udac,uesri,uffg,ugerji,ugermo,ugersi,ugexj,ugexmo,ugexp,ugexr,uhfgex,uobsw,uprphi,uprspi,ureac,usbsw,utitl1,utitl2,uvfgex,uxcat,uxkgex,uxmod,uxopex,uxopt,uzveci,vreac,xgersi,xhfgex,ximaxi,xistti,xlkffg,xlkgex,xlkmod,xvfgex,zgexj,zvclgi)
    !! This subroutine reads the EQ6 input file in menu-style ("D")
    !! format for version 8.0.
    !! This subroutine is a near-clone of EQ6/rd6ind.f. However, the
    !! present subroutine embodies only a pure read function (it does
    !! only minimal checking of what is read to ensure that what
    !! follows is readable). EQ6/rd6ind.f differs in that it also
    !! writes an instant echo of what is read to the EQ6 output file.
    !! The calling sequence of this subroutine is identical to that of
    !! EQ6/rd6ind.f, EQ6/rd6inw.f, and XCON6/rd6w8.f.
    !! This subroutine is called by:
    !!   XCON6/xcon6.f
    !! Principal input:
    !!   ninpts = unit number of the stripped input file
    !!   noutpt = unit number of the output file
    !!   nttyo  = unit number of the screen file
    !! Principal output:
    !!   qrderr = flag denoting a read error
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

    integer :: ninpts
    integer :: noutpt
    integer :: nttyo

    integer :: iact(imchmx,2,nrctmx)
    integer :: ibsrti(nsrtmx)
    integer :: igerti(jetmax,nertmx)
    integer :: iesrti(nsrtmx)
    integer :: imech(2,nrctmx)
    integer :: iodb(nodbmx)
    integer :: iopg(nopgmx)
    integer :: iopr(noprmx)
    integer :: iopt(noptmx)
    integer :: ixrti(nxrtmx)
    integer :: jcode(nrctmx)
    integer :: jgerti(nertmx)
    integer :: jflgi(nbtmax)
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
    integer :: nert
    integer :: net
    integer :: nffg
    integer :: nobswt
    integer :: nprob
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

    logical :: qend
    logical :: qgexsh
    logical :: qrderr

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

    include 'eqlib/eqlo8.h'
    include 'eqlib/eqlk8.h'

    ! Local parameter declarations.
    !   nfldpa = maximum number of fields per line
    !   nlchpa = character length of a line
    integer :: nfldpa
    integer :: nlchpa

    parameter (nfldpa = 8,nlchpa = 80)

    ! Local variable declarations.
    integer :: i
    integer :: icount
    integer :: ii
    integer :: iki
    integer :: imh
    integer :: ival
    integer :: ivar
    integer :: j
    integer :: jcox
    integer :: jd
    integer :: je
    integer :: jei
    integer :: jj
    integer :: jlast
    integer :: jrex
    integer :: jnrk
    integer :: j2
    integer :: j3
    integer :: j4
    integer :: krow
    integer :: k1
    integer :: k2
    integer :: k3
    integer :: k4
    integer :: kxmd
    integer :: n
    integer :: nbi
    integer :: ndt
    integer :: ne
    integer :: nei
    integer :: ner
    integer :: nfi
    integer :: nfldmx
    integer :: nfldt
    integer :: nfldtx
    integer :: nlchmx
    integer :: nmark
    integer :: nn
    integer :: nnn
    integer :: nobsw
    integer :: npi
    integer :: nsbsw
    integer :: nsr
    integer :: nsi
    integer :: nxi
    integer :: nxic
    integer :: nrc
    integer :: nxr
    integer :: nci

    integer :: ilnobl

    logical :: qgexbf
    logical :: qgexbs
    logical :: qgexrd
    logical :: qmark
    logical :: qnone
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
    character(len=48) :: ux48
    character(len=24) :: ustr24
    character(len=24) :: ustrn
    character(len=8) :: ux8
    character(len=1) :: ux1

    real(kind=8) :: var

    ! The following is nonsense to avoid compiler "unused variable"
    ! warnings. Here noutpt and nprob are not actually used. They are
    ! included in the calling sequence to allow it to match that of
    ! EQ6/rd6inw.f.
    noutpt = nttyo
    i = nprob

    nfldmx = nfldpa
    nlchmx = nlchpa

    qrderr = .false.

    ! Check some dimensioning parameters.
    qstop = .false.

    if (ndbxpa .ne. nodbmx) then
        write (nttyo,3000) ndbxpa,nodbmx
3000 format(/' * Error - (XCON6/rd6d8) The dimensioning parameter',' for the',/7x,'number of iodb debugging print option switches',' with string definitions',/7x,'(ndbxpa) has a value of ',i3,', but the dimensioning',/7x,'parameter for the number',' of such switches (nodbpa) has a',/7x,'value of ',i3,'.')

        qstop = .true.
    end if

    if (npgxpa .ne. nopgmx) then
        write (nttyo,3010) npgxpa,nopgmx
3010 format(/' * Error - (XCON6/rd6d8) The dimensioning parameter',' for the',/7x,'number of iopg activity coefficient option',' switches with string definitions',/7x,'(npgxpa) has a value',' of ' ,i3,', but the dimensioning',/7x,'parameter for the',' number of such switches (nopgpa) has a',/7x,'value of ',i3,'.')

        qstop = .true.
    end if

    if (nprxpa .ne. noprmx) then
        write (nttyo,3020) nprxpa,noprmx
3020 format(/' * Error - (XCON6/rd6d8) The dimensioning parameter',' for the',/7x,'number of iopr print option switches',' with string definitions',/7x,'(nprxpa) has a value of ',i3,', but the dimensioning',/7x,'parameter for the number',' of such switches (noprpa) has a',/7x,'value of ',i3,'.')

        qstop = .true.
    end if

    if (nptxpa .ne. noptmx) then
        write (nttyo,3030) nptxpa,noptmx
3030 format(/' * Error - (XCON6/rd6d8) The dimensioning parameter',' for the',/7x,'number of iopt model option switches',' with string definitions',/7x,'(nptxpa) has a value of ',i3,', but the dimensioning',/7x,'parameter for the number',' of such switches (noptpa) has a',/7x,'value of ',i3,'.')

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
        write (nttyo,1010) uline1(1:j2)
1010 format(/' * Error - (XCON6/rd6d8) The first line of a "D"',' format input file',/7x,'should be a separator line;',' therefore, it should begin with',/7x,'"|-------". The first',' line begins instead with',/7x,'"',a,'".')

        go to 990
    end if

    go to 105

100 continue
    qend = .true.
    go to 999

105 continue

    ! Read the block title ("Main Title") from a two-line header.
    uheadx = 'Main Title'
    nfldtx = 2
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    ! Now read the main title itself.
    n = 0

    do nn = 1,ntitmx + 1
        read (ninpts,1000,err=990) uline1
        call parsln(nfldmx,nfldt,nlchmx,ufield,uline1,ulscr)
        ustr = ufield(1)

        ! A separator line terminates the this block. It is not part
        ! of the title itself.
        if (ustr(1:8) .eq. '--------') then
            go to 120
        end if

        n = n + 1

        if (n .gt. ntitmx) then
            write (nttyo,1015) ntitmx
1015 format(/' * Error - (XCON6/rd6d8) Have too many lines in',/7x,'the main title. The code is only dimensioned for ',i4,/7x,'lines. Reduce the size of the main title or increase',/7x,'the dimensioning parameter ntitpa.')

            go to 990
        end if

        utitl1(n) = ufield(1)
    end do

120 continue
    ntitl1 = n

    ! Temperature.
    jtemp = -1
    icount = 0
    tempcb = 0.
    ttk(1) = 0.
    ttk(2) = 0.

    ! Read a one-line header.
    uheadx = 'Temperature option (jtemp):'
    nfldtx = 1
    call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)

    if (qrderr) then
        go to 999
    end if

    ! Read the first option (constant temperature) from a one-line
    ! header.
    nfldtx = 1
    call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

    if (qrderr) then
        go to 999
    end if

    ustr = ufield(1)
    ustr24 = ustr(5:29)
    uheadx = '( 0) Constant temperatur'
    call locase(ustr24)
    call locase(uheadx)
    j2 = ilnobl(ustr24)
    j3 = ilnobl(uheadx)

    if (ustr24(1:j2) .ne. uheadx(1:j3)) then
        write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
1020 format(/' * Error - (XCON6/rd6d8) Was expecting to find a',' line beginning with',/7x,'"',a,'", instead found one',' beginning with',/7x,'"',a,'".')

        qrderr = .true.
        go to 999
    end if

    if (ustr(1:1).ne.'[' .or. ustr(3:3).ne.']') then
        j2 = ilnobl(ustr)
        j2 = min(j2,70)
        write (nttyo,1022) ustr(1:j2)
1022 format(/' * Error - (XCON6/rd6d8) Was expecting to find an',/7x,'option check box "[ ]" in the line beginning with',/7x,'"',a,'".')

        qrderr = .true.
        go to 999
    end if

    ux1 = ustr(2:2)

    if (ux1.eq.'*' .or. ux1.eq.'x' .or. ux1.eq.'X') then
        jtemp = 0
        icount = icount + 1
    end if

    ! Read the associated temperature ("tempcb", C) from a one-line
    ! header.
    uheadx = 'Value (C)'
    nfldtx = 3
    call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)

    if (qrderr) then
        go to 999
    end if

    if (jtemp .eq. 0) then
        ustr = ufield(2)
        call chreal(nttyo,qrderr,ustr,var)

        if (qrderr) then
            go to 999
        end if

        tempcb = var
    end if

    ! Read the second option (linear tracking in Xi) from a one-line
    ! header.
    nfldtx = 1
    call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

    if (qrderr) then
        go to 999
    end if

    ustr = ufield(1)
    ustr24 = ustr(5:29)
    uheadx = '( 1) Linear tracking in '
    call locase(ustr24)
    call locase(uheadx)
    j2 = ilnobl(ustr24)
    j3 = ilnobl(uheadx)

    if (ustr24(1:j2) .ne. uheadx(1:j3)) then
        write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
        qrderr = .true.
        go to 999
    end if

    if (ustr(1:1).ne.'[' .or. ustr(3:3).ne.']') then
        j2 = ilnobl(ustr)
        j2 = min(j2,70)
        write (nttyo,1022) ustr(1:j2)
        qrderr = .true.
        go to 999
    end if

    ux1 = ustr(2:2)

    if (ux1.eq.'*' .or. ux1.eq.'x' .or. ux1.eq.'X') then
        jtemp = 1
        icount = icount + 1
    end if

    ! Read the associated base temperature ("tempcb", C) from a one-line
    ! header.
    uheadx = 'Base Value (C)'
    nfldtx = 3
    call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)

    if (qrderr) then
        go to 999
    end if

    if (jtemp .eq. 1) then
        ustr = ufield(2)
        call chreal(nttyo,qrderr,ustr,var)

        if (qrderr) then
            go to 999
        end if

        tempcb = var
    end if

    ! Read the associated derivative ("ttk(1)", dT/dXi, C/mol) from a
    ! one-line header.
    uheadx = 'Derivative'
    nfldtx = 3
    call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)

    if (qrderr) then
        go to 999
    end if

    if (jtemp .eq. 1) then
        ustr = ufield(2)
        call chreal(nttyo,qrderr,ustr,var)

        if (qrderr) then
            go to 999
        end if

        ttk(1) = var
    end if

    ! Read the third option (linear tracking in time) from a one-line
    ! header.
    nfldtx = 1
    call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

    if (qrderr) then
        go to 999
    end if

    ustr = ufield(1)
    ustr24 = ustr(5:29)
    uheadx = '( 2) Linear tracking in '
    call locase(ustr24)
    call locase(uheadx)
    j2 = ilnobl(ustr24)
    j3 = ilnobl(uheadx)

    if (ustr24(1:j2) .ne. uheadx(1:j3)) then
        write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
        qrderr = .true.
        go to 999
    end if

    if (ustr(1:1).ne.'[' .or. ustr(3:3).ne.']') then
        j2 = ilnobl(ustr)
        j2 = min(j2,70)
        write (nttyo,1022) ustr(1:j2)
        qrderr = .true.
        go to 999
    end if

    ux1 = ustr(2:2)

    if (ux1.eq.'*' .or. ux1.eq.'x' .or. ux1.eq.'X') then
        jtemp = 2
        icount = icount + 1
    end if

    ! Read the associated base temperature ("tempcb", C) from a one-line
    ! header.
    uheadx = 'Base Value (C)'
    nfldtx = 3
    call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)

    if (qrderr) then
        go to 999
    end if

    if (jtemp .eq. 2) then
        ustr = ufield(2)
        call chreal(nttyo,qrderr,ustr,var)

        if (qrderr) then
            go to 999
        end if

        tempcb = var
    end if

    ! Read the associated derivative ("ttk(1)", dT/dt, C/sec) from a
    ! one-line header.
    uheadx = 'Derivative'
    nfldtx = 3
    call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)

    if (qrderr) then
        go to 999
    end if

    if (jtemp .eq. 2) then
        ustr = ufield(2)
        call chreal(nttyo,qrderr,ustr,var)

        if (qrderr) then
            go to 999
        end if

        ttk(1) = var
    end if

    ! Read the fourth option (fluid mixing tracking) option from a
    ! one-line header.
    nfldtx = 1
    call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

    if (qrderr) then
        go to 999
    end if

    ustr = ufield(1)
    ustr24 = ustr(5:29)
    uheadx = '( 3) Fluid mixing tracki'
    call locase(ustr24)
    call locase(uheadx)
    j2 = ilnobl(ustr24)
    j3 = ilnobl(uheadx)

    if (ustr24(1:j2) .ne. uheadx(1:j3)) then
        write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
        go to 999
    end if

    if (ustr(1:1).ne.'[' .or. ustr(3:3).ne.']') then
        j2 = ilnobl(ustr)
        j2 = min(j2,70)
        write (nttyo,1022) ustr(1:j2)
        qrderr = .true.
        go to 999
    end if

    ux1 = ustr(2:2)

    if (ux1.eq.'*' .or. ux1.eq.'x' .or. ux1.eq.'X') then
        jtemp = 3
        icount = icount + 1
    end if

    ! Read the temperature of fluid 1 ("tempcb", C) from a one-line
    ! header.
    uheadx = 'T of fluid 1 (C)'
    nfldtx = 3
    call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)

    if (qrderr) then
        go to 999
    end if

    if (jtemp .eq. 3) then
        ustr = ufield(2)
        call chreal(nttyo,qrderr,ustr,var)

        if (qrderr) then
            go to 999
        end if

        tempcb = var
    end if

    ! Read the temperature of fluid 2 ("ttk(2)", C) from a one-line
    ! header.
    uheadx = 'T of fluid 2 (C)'
    nfldtx = 3
    call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)

    if (qrderr) then
        go to 999
    end if

    if (jtemp .eq. 3) then
        ustr = ufield(2)
        call chreal(nttyo,qrderr,ustr,var)

        if (qrderr) then
            go to 999
        end if

        ttk(2) = var
    end if

    ! Read the mass ratio factor ("ttk(1)") from a two-line header.
    ! This this the mass ratio of fluid2/fluid1 at Xi = 1.
    uheadx = 'Mass ratio factor'
    nfldtx = 3
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    if (jtemp .eq. 3) then
        ustr = ufield(2)
        call chreal(nttyo,qrderr,ustr,var)

        if (qrderr) then
            go to 999
        end if

        ttk(1) = var
    end if

    if (icount.eq.0 .or. jtemp.lt.0) then
        write (nttyo,1023) tempcb
1023 format(/' * Warning - (XCON6/rd6d8) No option was selected for',/7x,'the temperature. The temperature will be set at a fixed',' value',/7x,'of ',f7.2,'C.')

        jtemp = 0
        ttk(1) = 0.
        ttk(2) = 0.
    else if (icount .gt. 1) then
        write (nttyo,1024) tempcb
1024 format(/' * Warning - (XCON6/rd6d8) Multiple options were',' selected for',/7x,'the temperature. The temperature will be',' set at a fixed value',/7x,'of ',f7.2,'C.')

        jtemp = 0
        ttk(1) = 0.
        ttk(2) = 0.
    end if

    ! Pressure.
    jpress = -1
    icount = 0
    pressb = 0.
    ptk(1) = 0.
    ptk(2) = 0.

    ! Read a one-line header.
    uheadx = 'Pressure option (jpress):'
    nfldtx = 1
    call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)

    if (qrderr) then
        go to 999
    end if

    ! Read the first option (follow the data file reference pressure
    ! curve) from a one-line header.
    nfldtx = 1
    call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

    if (qrderr) then
        go to 999
    end if

    ustr = ufield(1)
    ustr24 = ustr(5:29)
    uheadx = '( 0) Follow the data fil'
    call locase(ustr24)
    call locase(uheadx)
    j2 = ilnobl(ustr24)
    j3 = ilnobl(uheadx)

    if (ustr24(1:j2) .ne. uheadx(1:j3)) then
        write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
        qrderr = .true.
        go to 999
    end if

    if (ustr(1:1).ne.'[' .or. ustr(3:3).ne.']') then
        j2 = ilnobl(ustr)
        j2 = min(j2,70)
        write (nttyo,1022) ustr(1:j2)
        qrderr = .true.
        go to 999
    end if

    ux1 = ustr(2:2)

    if (ux1.eq.'*' .or. ux1.eq.'x' .or. ux1.eq.'X') then
        jpress = 0
        pressb = 0.
        icount = icount + 1
    end if

    ! Read the second option (follow the 1.013-bar/steam-saturation
    ! curve) from a one-line header.
    nfldtx = 1
    call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

    if (qrderr) then
        go to 999
    end if

    ustr = ufield(1)
    ustr24 = ustr(5:29)
    uheadx = '( 1) Follow the 1.013-ba'
    call locase(ustr24)
    call locase(uheadx)
    j2 = ilnobl(ustr24)
    j3 = ilnobl(uheadx)

    if (ustr24(1:j2) .ne. uheadx(1:j3)) then
        write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
        qrderr = .true.
        go to 999
    end if

    if (ustr(1:1).ne.'[' .or. ustr(3:3).ne.']') then
        j2 = ilnobl(ustr)
        j2 = min(j2,70)
        write (nttyo,1022) ustr(1:j2)
        qrderr = .true.
        go to 999
    end if

    ux1 = ustr(2:2)

    if (ux1.eq.'*' .or. ux1.eq.'x' .or. ux1.eq.'X') then
        jpress = 1
        pressb = 0.
        icount = icount + 1
    end if

    ! Read the third option (specified constant pressure) from a
    ! one-line header.
    nfldtx = 1
    call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

    if (qrderr) then
        go to 999
    end if

    ustr = ufield(1)
    ustr24 = ustr(5:29)
    uheadx = '( 2) Constant pressure: '
    call locase(ustr24)
    call locase(uheadx)
    j2 = ilnobl(ustr24)
    j3 = ilnobl(uheadx)

    if (ustr24(1:j2) .ne. uheadx(1:j3)) then
        write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
        qrderr = .true.
        go to 999
    end if

    if (ustr(1:1).ne.'[' .or. ustr(3:3).ne.']') then
        j2 = ilnobl(ustr)
        j2 = min(j2,70)
        write (nttyo,1022) ustr(1:j2)
        qrderr = .true.
        go to 999
    end if

    ux1 = ustr(2:2)

    if (ux1.eq.'*' .or. ux1.eq.'x' .or. ux1.eq.'X') then
        jpress = 2
        icount = icount + 1
    end if

    ! Read the associated pressure ("pressb", bars) from a one-line
    ! header.
    uheadx = 'Value (bars)'
    nfldtx = 3
    call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)

    if (qrderr) then
        go to 999
    end if

    if (jpress .eq. 2) then
        ustr = ufield(2)
        call chreal(nttyo,qrderr,ustr,var)

        if (qrderr) then
            go to 999
        end if

        pressb = var
    end if

    ! Read the fourth option (linear tracking in Xi) from a one-line
    ! header.
    nfldtx = 1
    call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

    if (qrderr) then
        go to 999
    end if

    ustr = ufield(1)
    ustr24 = ustr(5:29)
    uheadx = '( 3) Linear tracking in '
    call locase(ustr24)
    call locase(uheadx)
    j2 = ilnobl(ustr24)
    j3 = ilnobl(uheadx)

    if (ustr24(1:j2) .ne. uheadx(1:j3)) then
        write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
        go to 999
    end if

    if (ustr(1:1).ne.'[' .or. ustr(3:3).ne.']') then
        j2 = ilnobl(ustr)
        j2 = min(j2,70)
        write (nttyo,1022) ustr(1:j2)
        qrderr = .true.
        go to 999
    end if

    ux1 = ustr(2:2)

    if (ux1.eq.'*' .or. ux1.eq.'x' .or. ux1.eq.'X') then
        jpress = 3
        icount = icount + 1
    end if

    ! Read the base pressure ("pressb", bars) from a one-line header.
    uheadx = 'Base Value (bars)'
    nfldtx = 3
    call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)

    if (qrderr) then
        go to 999
    end if

    if (jpress .eq. 3) then
        ustr = ufield(2)
        call chreal(nttyo,qrderr,ustr,var)

        if (qrderr) then
            go to 999
        end if

        pressb = var
    end if

    ! Read the associated derivative ("ptk(1)", dP/dXi, bars/mol) from
    ! a one-line header.
    uheadx = 'Derivative'
    nfldtx = 3
    call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)

    if (qrderr) then
        go to 999
    end if

    if (jpress.eq. 3) then
        ustr = ufield(2)
        call chreal(nttyo,qrderr,ustr,var)

        if (qrderr) then
            go to 999
        end if

        ptk(1) = var
    end if

    ! Read the fifth option (linear tracking in time) from a one-line
    ! header.
    nfldtx = 1
    call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

    if (qrderr) then
        go to 999
    end if

    ustr = ufield(1)
    ustr24 = ustr(5:29)
    uheadx = '( 4) Linear tracking in '
    call locase(ustr24)
    call locase(uheadx)
    j2 = ilnobl(ustr24)
    j3 = ilnobl(uheadx)

    if (ustr24(1:j2) .ne. uheadx(1:j3)) then
        write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
        go to 999
    end if

    if (ustr(1:1).ne.'[' .or. ustr(3:3).ne.']') then
        j2 = ilnobl(ustr)
        j2 = min(j2,70)
        write (nttyo,1022) ustr(1:j2)
        qrderr = .true.
        go to 999
    end if

    ux1 = ustr(2:2)

    if (ux1.eq.'*' .or. ux1.eq.'x' .or. ux1.eq.'X') then
        jpress = 4
        icount = icount + 1
    end if

    ! Read the base pressure ("pressb", bars) from a one-line header.
    uheadx = 'Base Value (bars)'
    nfldtx = 3
    call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)

    if (qrderr) then
        go to 999
    end if

    if (jpress .eq. 4) then
        ustr = ufield(2)
        call chreal(nttyo,qrderr,ustr,var)

        if (qrderr) then
            go to 999
        end if

        pressb = var
    end if

    ! Read the associated derivative ("ptk(1)", dP/dt, bars/sec) from
    ! a two-line header.
    uheadx = 'Derivative'
    nfldtx = 3
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    if (jpress .eq. 4) then
        ustr = ufield(2)
        call chreal(nttyo,qrderr,ustr,var)

        if (qrderr) then
            go to 999
        end if

        ptk(1) = var
    end if

    if (icount.eq.0 .or. jpress.lt.0) then
        write (nttyo,1030)
1030 format(/' * Warning - (XCON6/rd6d8) No option was selected for',' the pressure.',/7x,'The pressure will be set to be in',' accord with the',/7x,'data file reference pressure curve.')

        jpress = 0
        pressb = 0.
        ptk(1) = 0.
        ptk(2) = 0.
    else if (icount .gt. 1) then
        write (nttyo,1026)
1026 format(/' * Warning - (XCON6/rd6d8) Multiple options were',' selected for',/7x,'the pressure. The pressure will be set',' to be in accord with the',/7x,'data file reference pressure',' curve.')

        jpress = 0
        pressb = 0.
        ptk(1) = 0.
        ptk(2) = 0.
    end if

    ! Reactants.
    nrct = 0
    nsrt = 0
    nxrt = 0
    nert = 0

    ! Read a two-line header for the block.
    uheadx = 'Reactants (Irreversible Reactions)'
    nfldtx = 2
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    ! Loop on Reactants.
    nrc = 0
    nsr = 0
    nxr = 0
    ner = 0

    do nn = 1,nrctmx + 1
        ! Read a line. If the block has not been completely read, this
        ! contains the name of a reactant. Otherwise, this line is the
        ! first line of the block following the reactants super-block.
        nfldtx = 0
        call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

        if (qrderr) then
            go to 999
        end if

        ustr = ufield(1)
        uheadx = 'Reactant'
        call locase(ustr)
        call locase(uheadx)
        j2 = ilnobl(ustr)
        j3 = ilnobl(uheadx)

        if (ustr(1:j2) .ne. uheadx(1:j3)) then
            ! Back up.
            backspace ninpts

            ! Done reading the reactants super-block.
            go to 240
        else
            ! Back up.
            backspace ninpts

            ! Re-read the reactant line with the test on the number
            ! of fields.
            nfldtx = 3
            call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

            if (qrderr) then
                go to 999
            end if
        end if

        ustr = ufield(2)
        ustrn = ustr(1:24)
        call locase(ustrn)

        if (ustrn(1:5).eq.'none ' .or. ustr(1:1).eq.' ') then
            ustr = 'None'
        end if

        if (ustr(1:5) .eq. 'None ') then
            nrc = 0
        else
            nrc = nrc + 1
        end if

        if (nrc .gt. nrctmx) then
            write (nttyo,1520) nrctmx
1520 format(/' * Error - (XCON6/rd6d8) Have too many reactants',/7x,'The code is only dimensioned for ',i4,' reactants.',/7x,'Reduce the number of reactants or increase the',/7x,'dimensioning parameter nrctpa.')

            go to 990
        end if

        if (nrc .gt. 0) then
            ureac(nrc) = ustr(1:24)
        end if

        ! Read the separator line following the line containing the name
        ! of a reactant.
        nfldtx = 1
        call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

        if (qrderr) then
            go to 999
        end if

        ustr24 = ufield(1)(1:24)
        uheadx = '--------'

        if (ustr24(1:8) .ne. uheadx(1:8)) then
            j2 = ilnobl(uline1)
            j2 = min(j2,70)
            write (nttyo,1070) uline1(1:j2)
            qrderr = .true.
            go to 999
        end if

        ! Read the type of reactant from a two-line header.
        uheadx = '->'
        nfldtx = 4
        call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

        if (qrderr) then
            go to 999
        end if

        ustr24 = ufield(2)(1:24)
        uheadx = 'Type'
        call locase(ustr24)
        call locase(uheadx)
        j2 = ilnobl(ustr24)
        j3 = ilnobl(uheadx)

        if (ustr24(1:j2) .ne. uheadx(1:j3)) then
            write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
            qrderr = .true.
            go to 999
        end if

        ustr = ufield(3)
        call locase(ustr)

        if (nrc .eq. 0) then
            go to 150
        end if

        do jcox = 0,5
            uheadx = urcjco(jcox)
            call locase(uheadx)

            if (ustr(1:16) .eq. uheadx(1:16)) then
                jcode(nrc) = jcox
                go to 150
            end if
        end do

        j2 = ilnobl(ustr)
        write (nttyo,1060) ustr(1:j2)
1060 format(/" * Error - (XCON6/rd6d8) Don't recognize the",' jcode option string',/7x,' "',a,'". This should',' be one of the strings',/7x,'defined in the urcjco array.',' The valid strings are:',/)

        do jcox = 0,5
            j3 = ilnobl(urcjco(jcox))
            write (nttyo,1062) urcjco(jcox)(1:j3)
1062 format(9x,a)
        end do

        go to 990

        ! Read the reactant status from a two-line header.
150 continue
        uheadx = '->'
        nfldtx = 4
        call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

        if (qrderr) then
            go to 999
        end if

        ustr24 = ufield(2)(1:24)
        uheadx = 'Status'
        call locase(ustr24)
        call locase(uheadx)
        j2 = ilnobl(ustr24)
        j3 = ilnobl(uheadx)

        if (ustr24(1:j2) .ne. uheadx(1:j3)) then
            write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
            qrderr = .true.
            go to 999
        end if

        if (nrc .eq. 0) then
            go to 160
        end if

        ustr = ufield(3)
        call locase(ustr)

        do jrex = -1,3
            uheadx = urcjre(jrex)
            call locase(uheadx)

            if (ustr(1:16) .eq. uheadx(1:16)) then
                jreac(nrc) = jrex
                go to 160
            end if
        end do

        j2 = ilnobl(ustr)
        write (nttyo,1065) ustr(1:j2)
1065 format(/" * Error - (XCON6/rd6d8) Don't recognize the",' jreac option string',/7x,' "',a,'". This should',' be one of the strings',/7x,'defined in the urcjre array.',' The valid strings are:',/)

        do jrex = -1,3
            j3 = ilnobl(urcjre(jrex))
            write (nttyo,1067) urcjre(jrex)(1:j3)
1067 format(9x,a)
        end do

        go to 990

        ! Read the amount of reactant remaining (moles)
        ! from a two-line header.
160 continue
        uheadx = '->'
        nfldtx = 4
        call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

        if (qrderr) then
            go to 999
        end if

        ustr24 = ufield(2)(1:24)
        uheadx = 'Amount remaining (moles)'
        call locase(ustr24)
        call locase(uheadx)
        j2 = ilnobl(ustr24)
        j3 = ilnobl(uheadx)

        if (ustr24(1:j2) .ne. uheadx(1:j3)) then
            write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
            qrderr = .true.
            go to 999
        end if

        if (nrc .gt. 0) then
            ustr = ufield(3)
            call chreal(nttyo,qrderr,ustr,var)

            if (qrderr) then
                go to 999
            end if

            morr(nrc) = var
        end if

        ! Read the amount of destroyed remaining (moles)
        ! from a two-line header.
        uheadx = '->'
        nfldtx = 4
        call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

        if (qrderr) then
            go to 999
        end if

        ustr24 = ufield(2)(1:24)
        uheadx = 'Amount destroyed (moles)'
        call locase(ustr24)
        call locase(uheadx)
        j2 = ilnobl(ustr24)
        j3 = ilnobl(uheadx)

        if (ustr24(1:j2) .ne. uheadx(1:j3)) then
            write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
            qrderr = .true.
            go to 999
        end if

        if (nrc .gt. 0) then
            ustr = ufield(3)
            call chreal(nttyo,qrderr,ustr,var)

            if (qrderr) then
                go to 999
            end if

            modr(nrc) = var
        end if

        ! Skip jcode (reactant type) tests if no reactant is present.
        if (nrc .eq. 0) then
            go to 235
        end if

        if (jcode(nrc) .eq. 1) then
            ! Have a solid solution.
            nxr = nxr + 1

            if (nxr .gt. nxrtmx) then
                write (nttyo,1550) nxrtmx
1550 format(/' * Error - (XCON6/rd6d8) Have too many solid',/7x,'solution reactants. The code is only dimensioned',/7x,'for ',i4,' such reactants. Reduce the number of',/7x,'such reactants or increase the dimensioning',' parameter nxrtpa.')

                go to 990
            end if

            ! Composition sub-block header.
            uheadx = '->'
            nfldtx = 2
            call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

            if (qrderr) then
                go to 999
            end if

            ustr24 = ufield(2)(1:24)
            uheadx = 'Composition'
            call locase(ustr24)
            call locase(uheadx)
            j2 = ilnobl(ustr24)
            j3 = ilnobl(uheadx)

            if (ustr24(1:j2) .ne. uheadx(1:j3)) then
                write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
                qrderr = .true.
                go to 999
            end if

            ! Composition table header.
            uheadx = '--->'
            nfldtx = 4
            call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

            if (qrderr) then
                go to 999
            end if

            ustr24 = ufield(2)(1:24)
            uheadx = 'Component'
            call locase(ustr24)
            call locase(uheadx)
            j2 = ilnobl(ustr24)
            j3 = ilnobl(uheadx)

            if (ustr24(1:j2) .ne. uheadx(1:j3)) then
                write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
                qrderr = .true.
                go to 999
            end if

            ! Loop on end-member components.
            iki = 0

            do jj = 1,iktmax + 1
                ! Read a line. If the sub-block (for the components) has not
                ! been completely read, this contains the name of a component
                ! and its mole fraction. Otherwise, this line is the separator
                ! line before the next sub-block (surface area).
                nfldtx = 0
                call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

                if (qrderr) then
                    go to 999
                end if

                ustr = ufield(1)

                if (ustr(1:8) .eq. '--------') then
                    go to 210
                end if

                ! Have found another component.
                iki = iki + 1

                if (iki .gt. iktmax) then
                    j2 = ilnobl(ureac(nrc))
                    write (nttyo,1590) ureac(nrc)(1:j2),iktmax
1590 format(/' * Error - (XCON6/rd6d8) Have too many',' end-members',/7x,'in the solid solution reactant',a,'.',/7x,'The code is only dimensioned for ',i4,' end-members per',/7x,'solid solution. Reduce',' the number of end-members or',/7x,'increase the dimensioning parameter iktpar.')

                    go to 990
                end if

                ustr = ufield(3)
                call chreal(nttyo,qrderr,ustr,var)

                if (qrderr) then
                    go to 999
                end if

                rxbari(iki,nxr) = var
                ucxri(iki,nxr) = ufield(2)(1:24)
            end do

210 continue
            ixrti(nxr) = iki
        end if

        if (jcode(nrc) .eq. 2) then
            ! Have a special reactant. Read its molar volume, composition,
            ! and reaction.
            nsr = nsr + 1

            if (nsr .gt. nsrtmx) then
                write (nttyo,1600) nsrtmx
1600 format(/' * Error - (XCON6/rd6d8) Have too many special',/7x,'reactants. The code is only dimensioned for ',i4,/7x,'such reactants. Reduce the number of such reactants',/7x,'or increase the dimensioning parameter nsrtpa.')

                go to 990
            end if

            ! Read the molar volume (cm3/mol) from a two-line header.
            uheadx = '->'
            nfldtx = 4
            call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

            if (qrderr) then
                go to 999
            end if

            ustr24 = ufield(2)(1:24)
            uheadx = 'Molar volume (cm3/mol)'
            call locase(ustr24)
            call locase(uheadx)
            j2 = ilnobl(ustr24)
            j3 = ilnobl(uheadx)

            if (ustr24(1:j2) .ne. uheadx(1:j3)) then
                write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
                qrderr = .true.
                go to 999
            end if

            ustr = ufield(3)
            call chreal(nttyo,qrderr,ustr,var)

            if (qrderr) then
                go to 999
            end if

            vreac(nrc) = var

            ! Composition sub-block header.
            uheadx = '->'
            nfldtx = 2
            call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

            if (qrderr) then
                go to 999
            end if

            ustr24 = ufield(2)(1:24)
            uheadx = 'Composition'
            call locase(ustr24)
            call locase(uheadx)
            j2 = ilnobl(ustr24)
            j3 = ilnobl(uheadx)

            if (ustr24(1:j2) .ne. uheadx(1:j3)) then
                write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
                qrderr = .true.
                go to 999
            end if

            ! Composition table header.
            uheadx = '--->'
            nfldtx = 4
            call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

            if (qrderr) then
                go to 999
            end if

            ustr24 = ufield(2)(1:24)
            uheadx = 'Element'
            call locase(ustr24)
            call locase(uheadx)
            j2 = ilnobl(ustr24)
            j3 = ilnobl(uheadx)

            if (ustr24(1:j2) .ne. uheadx(1:j3)) then
                write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
                qrderr = .true.
                go to 999
            end if

            ! Loop on elements.
            nci = 0

            do jj = 1,nctmax + 1
                ! Read a line. If the sub-block (for the elements) has not
                ! been completely read, this contains the name of an element
                ! and  the stoichiometric number. Otherwise, this line is the
                ! separator line before the next sub-block (reaction).
                nfldtx = 0
                call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

                if (qrderr) then
                    go to 999
                end if

                ustr = ufield(1)

                if (ustr(1:8) .eq. '--------') then
                    go to 220
                end if

                ! Have found another element.
                nci = nci + 1

                if (nci .gt. nctmax) then
                    j2 = ilnobl(ureac(nrc))
                    write (nttyo,1660) ureac(nrc)(1:j2),nctmax
1660 format(/' * Error - (XCON6/rd6d8) Have too many','  chemical',/7x,'elements in the special reactant ',a,'.',/7x,'The code is only dimensioned for ',i4,' elements.',/7x,'Reduce the number of elements or',' increase the',/7x,'dimensioning parameter nctpar.')

                    go to 990
                end if

                ustr = ufield(3)
                call chreal(nttyo,qrderr,ustr,var)

                if (qrderr) then
                    go to 999
                end if

                cesri(nci,nsr) = var
                uesri(nci,nsr) = ufield(2)(1:8)
            end do

220 continue
            iesrti(nsr) = nci

            ! Reaction header.
            uheadx = '->'
            nfldtx = 2
            call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

            if (qrderr) then
                go to 999
            end if

            ustr24 = ufield(2)(1:24)
            uheadx = 'Reaction'
            call locase(ustr24)
            call locase(uheadx)
            j2 = ilnobl(ustr24)
            j3 = ilnobl(uheadx)

            if (ustr24(1:j2) .ne. uheadx(1:j3)) then
                write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
                qrderr = .true.
                go to 999
            end if

            ! Species in the reaction.
            uheadx = '--->'
            nfldtx = 4
            call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

            if (qrderr) then
                go to 999
            end if

            ustr24 = ufield(2)(1:24)
            uheadx = 'Species'
            call locase(ustr24)
            call locase(uheadx)
            j2 = ilnobl(ustr24)
            j3 = ilnobl(uheadx)

            if (ustr24(1:j2) .ne. uheadx(1:j3)) then
                write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
                qrderr = .true.
                go to 999
            end if

            ! Loop on species.
            nbi = 0

            do n = 1,nbt1mx+ 1
                ! Read a line. If the sub-block (for the species) has not
                ! been completely read, this contains the name of a species
                ! and  the reaction coefficient. Otherwise, this line is the
                ! separator line before the next sub-block (reaction).
                nfldtx = 0
                call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

                if (qrderr) then
                    go to 999
                end if

                ustr = ufield(1)

                if (ustr(1:8) .eq. '--------') then
                    go to 230
                end if

                ! Have found another species.
                nbi = nbi + 1

                if (nbi .gt. nbt1mx) then
                    j2 = ilnobl(ureac(nrc))
                    write (nttyo,1700) ureac(nrc)(1:j2),nbtmax
1700 format(/' * Error - (XCON6/rd6d8) Have too many basis',' basis species in the',/7x,'in the reaction for the',' special reactant ',a,'.',/7x,'The code is only',' dimensioned for ',i4,' basis species. Increase',/7x,'the dimensioning parameter nbtpar.')

                    go to 990
                end if

                ustr = ufield(3)
                call chreal(nttyo,qrderr,ustr,var)

                if (qrderr) then
                    go to 999
                end if

                cbsri(n,nsr) = var
                ubsri(n,nsr) = ufield(2)(1:24)
            end do

230 continue
            ibsrti(nsr) = nbi
        end if

        if (jcode(nrc) .eq. 5) then
            ! Have a generic ion exchanger.
            ner = ner + 1

            if (ner .gt. nertmx) then
                write (nttyo,1710) nertmx
1710 format(/' * Error - (XCON6/rd6d8) Have too many generic',' ion exchanger',/7x,'reactants. The code is only',' dimensioned for ',i4,' such reactants.',/7x,'Reduce',' the number of such reactants or increase the',/7x,'dimensioning parameter nertpa.')

                go to 990
            end if

            ! Read the exchange model name for the exchanger reactant
            ! from a two-line header.
            uheadx = '->'
            nfldtx = 4
            call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

            if (qrderr) then
                go to 999
            end if

            ustr24 = ufield(2)(1:24)
            uheadx = 'Exchange model'
            call locase(ustr24)
            call locase(uheadx)
            j2 = ilnobl(ustr24)
            j3 = ilnobl(uheadx)

            if (ustr24(1:j2) .ne. uheadx(1:j3)) then
                write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
                qrderr = .true.
                go to 999
            end if

            ugermo(ner) = ufield(3)(1:24)

            ! Composition sub-block header.
            uheadx = '->'
            nfldtx = 2
            call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

            if (qrderr) then
                go to 999
            end if

            ustr24 = ufield(2)(1:24)
            uheadx = 'Composition'
            call locase(ustr24)
            call locase(uheadx)
            j2 = ilnobl(ustr24)
            j3 = ilnobl(uheadx)

            if (ustr24(1:j2) .ne. uheadx(1:j3)) then
                write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
                qrderr = .true.
                go to 999
            end if

            ! Loop on exchange sites.
            jei = 0

            do jj = 1,jetmax + 1
                ! Read a line. If the sub-block (for the composition of the
                ! current exchanger reactant) has not been completely read,
                ! this contains the name of an exchange site, and a sub-sub-
                ! block for that site follows. Otherwise, this line is the
                ! first line following the composition sub-block for the
                ! current exchanger reactant.
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
                    ustr24 = ufield(2)(1:24)
                    uheadx = 'Exchange site'
                    call locase(ustr24)
                    call locase(uheadx)
                    j2 = ilnobl(ustr24)
                    j3 = ilnobl(uheadx)

                    if (ustr24(1:j2) .ne. uheadx(1:j3)) then
                        write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
                        qrderr = .true.
                        go to 999
                    end if

                    ! Have found another exchange site sub-sub-block.
                    jei = jei + 1

                    if (jei .gt. jetmax) then
                        write (nttyo,1716) ureac(nrc)(1:j2),jetmax
1716 format(/' * Error - (XCON6/rd6d8) Have too many',' exchange sites',/7x,'in the generic ion exchanger',' reactant ',a,'.',/7x,'The code is only dimensioned',' for ',i4,' exchange sites per',/7x,'generic ion',' exchanger. Reduce the number of exchange sites or',/7x,'increase the dimensioning parameter jetpar.')

                        go to 990
                    end if

                    ugerji(jei,ner) = ufield(3)(1:8)
                else
                    ! Should have found the first line after the sub-block for
                    ! the composition of the current exchanger reactant.
                    ! Back up.
                    backspace ninpts

                    go to 234
                end if

                ! Read the separator line following the line containing
                ! the name of an exchange site.
                nfldtx = 1
                call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

                if (qrderr) then
                    go to 999
                end if

                ustr24 = ufield(1)(1:24)
                uheadx = '--------'

                if (ustr24(1:8) .ne. uheadx(1:8)) then
                    j2 = ilnobl(uline1)
                    j2 = min(j2,70)
                    write (nttyo,1070) uline1(1:j2)
                    qrderr = .true.
                    go to 999
                end if

                ! Read the sub-sub-sub-block for the composition of the
                ! current site. This consists of a table header, followed
                ! by a table.
                ! Composition table header.
                uheadx = '----->'
                nfldtx = 4
                call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

                if (qrderr) then
                    go to 999
                end if

                ustr24 = ufield(2)(1:24)
                uheadx = 'Component'
                call locase(ustr24)
                call locase(uheadx)
                j2 = ilnobl(ustr24)
                j3 = ilnobl(uheadx)

                if (ustr24(1:j2) .ne. uheadx(1:j3)) then
                    write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
                    qrderr = .true.
                    go to 999
                end if

                ! Loop on components for the current site.
                nei = 0

                do nnn = 1,netmax + 1
                    ! Read a line. If the sub-sub-sub-block (for the components
                    ! for the current site) has not been completely read, this
                    ! contains the name of a component and its mole fraction on
                    ! the current site. Otherwise, this line is the separator
                    ! line before the next sub-sub-block (for the next site
                    ! of the exchanger reactant) or the sub-block (surface
                    ! area) following the composition sub-block.
                    nfldtx = 0
                    call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

                    if (qrderr) then
                        go to 999
                    end if

                    ustr = ufield(1)

                    if (ustr(1:8) .eq. '--------') then
                        go to 232
                    end if

                    ! Have found another component.
                    nei = nei + 1

                    if (nei .gt. netmax) then
                        j2 = ilnobl(ureac(nrc))
                        j3 = ilnobl(ugerji(jei,ner))
                        write (nttyo,1724) ureac(nrc)(1:j2),ugerji(jei,ner)(1:j3),netmax
1724 format(/' * Error - (XCON6/rd6d8) Have too many',' species on exchange',/7x,'site ',a,' of the generic',' ion exchanger reactant',/7x,a,'. The code is only',' dimensioned for species per',/7x,'exchange site.',' Reduce the number of end-members or increase',/7x,'the dimensioning parameter netpar.')

                        go to 990
                    end if

                    ustr = ufield(3)
                    call chreal(nttyo,qrderr,ustr,var)

                    if (qrderr) then
                        go to 999
                    end if

                    ugersi(nei,jei,ner) = ufield(2)(1:24)
                    egersi(nei,jei,ner) = var
                    xgersi(nei,jei,ner) = var
                end do

232 continue
                igerti(jei,ner) = nei
            end do

234 continue
            jgerti(ner) = jei
        end if

        ! Surface area option.
        icount = 0
        nsk(nrc) = 0.
        sfcar(nrc) = 0.
        ssfcar(nrc) = 0.

        ! Read a one-line header.
        uheadx = '->'
        nfldtx = 2
        call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)

        if (qrderr) then
            go to 999
        end if

        ustr24 = ufield(2)(1:24)
        uheadx = 'Surface area option (nsk'
        call locase(ustr24)
        call locase(uheadx)
        j2 = ilnobl(ustr24)
        j3 = ilnobl(uheadx)

        if (ustr24(1:j2) .ne. uheadx(1:j3)) then
            write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
            qrderr = .true.
            go to 999
        end if

        ! Read the first option (constant surface area) from a one-line
        ! header.
        nfldtx = 2
        call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

        if (qrderr) then
            go to 999
        end if

        ustr = ufield(2)
        ustr24 = ustr(5:29)
        uheadx = '( 0) Constant surface ar'
        call locase(ustr24)
        call locase(uheadx)
        j2 = ilnobl(ustr24)
        j3 = ilnobl(uheadx)

        if (ustr24(1:j2) .ne. uheadx(1:j3)) then
            write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
            qrderr = .true.
            go to 999
        end if

        if (ustr(1:1).ne.'[' .or. ustr(3:3).ne.']') then
            j2 = ilnobl(ustr)
            j2 = min(j2,70)
            write (nttyo,1022) ustr(1:j2)
            qrderr = .true.
            go to 999
        end if

        ux1 = ustr(2:2)

        if (ux1.eq.'*' .or. ux1.eq.'x' .or. ux1.eq.'X') then
            nsk(nrc) = 0
            icount = icount + 1
        end if

        ! Read the associated constant surface area value (sfcar(n),
        ! cm2) from a one-line header.
        uheadx = '->'
        nfldtx = 4
        call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)

        if (qrderr) then
            go to 999
        end if

        ustr24 = ufield(2)(1:24)
        uheadx = 'Value (cm2)'
        call locase(ustr24)
        call locase(uheadx)
        j2 = ilnobl(ustr24)
        j3 = ilnobl(uheadx)

        if (ustr24(1:j2) .ne. uheadx(1:j3)) then
            write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
            qrderr = .true.
            go to 999
        end if

        ustr = ufield(3)
        call chreal(nttyo,qrderr,ustr,var)

        if (qrderr) then
            go to 999
        end if

        if (nsk(nrc) .le. 0) then
            sfcar(nrc) = var
        end if

        ! Read the second option (constant specific surface area)
        ! from a one-line header.
        nfldtx = 2
        call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

        if (qrderr) then
            go to 999
        end if

        ustr = ufield(2)
        ustr24 = ustr(5:29)
        uheadx = '( 1) Constant specific s'
        call locase(ustr24)
        call locase(uheadx)
        j2 = ilnobl(ustr24)
        j3 = ilnobl(uheadx)

        if (ustr24(1:j2) .ne. uheadx(1:j3)) then
            write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
            go to 999
        end if

        if (ustr(1:1).ne.'[' .or. ustr(3:3).ne.']') then
            j2 = ilnobl(ustr)
            j2 = min(j2,70)
            write (nttyo,1022) ustr(1:j2)
            qrderr = .true.
            go to 999
        end if

        ux1 = ustr(2:2)

        if (ux1.eq.'*' .or. ux1.eq.'x' .or. ux1.eq.'X') then
            nsk(nrc) = 1
            icount = icount + 1
        end if

        ! Read the associated constant specific surface area value
        ! (ssfcar(n), cm2/g) from a one-line header.
        uheadx = '->'
        nfldtx = 4
        call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)

        if (qrderr) then
            go to 999
        end if

        ustr24 = ufield(2)(1:24)
        uheadx = 'Value (cm2/g)'
        call locase(ustr24)
        call locase(uheadx)
        j2 = ilnobl(ustr24)
        j3 = ilnobl(uheadx)

        if (ustr24(1:j2) .ne. uheadx(1:j3)) then
            write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
            qrderr = .true.
            go to 999
        end if

        ustr = ufield(3)
        call chreal(nttyo,qrderr,ustr,var)

        if (qrderr) then
            go to 999
        end if

        if (nsk(nrc) .eq. 1) then
            ssfcar(nrc) = var
        end if

        ! Read the third option (n**2/3 growth law: current surface
        ! area) from a two-line header.
        nfldtx = 2
        call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

        if (qrderr) then
            go to 999
        end if

        ustr = ufield(2)
        ustr24 = ustr(5:29)
        uheadx = '( 2) n**2/3 growth law-'
        call locase(ustr24)
        call locase(uheadx)
        j2 = ilnobl(ustr24)
        j3 = ilnobl(uheadx)

        if (ustr24(1:j2) .ne. uheadx(1:j3)) then
            write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
            qrderr = .true.
            go to 999
        end if

        if (ustr(1:1).ne.'[' .or. ustr(3:3).ne.']') then
            j2 = ilnobl(ustr)
            j2 = min(j2,70)
            write (nttyo,1022) ustr(1:j2)
            qrderr = .true.
            go to 999
        end if

        ux1 = ustr(2:2)

        if (ux1.eq.'*' .or. ux1.eq.'x' .or. ux1.eq.'X') then
            nsk(nrc) = 2
            icount = icount + 1
        end if

        ! Read the associated current surface area value (sfcar(n),
        ! cm2) from a two-line header.
        uheadx = '->'
        nfldtx = 4
        call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

        if (qrderr) then
            go to 999
        end if

        ustr24 = ufield(2)(1:24)
        uheadx = 'Value (cm2)'
        call locase(ustr24)
        call locase(uheadx)
        j2 = ilnobl(ustr24)
        j3 = ilnobl(uheadx)

        if (ustr24(1:j2) .ne. uheadx(1:j3)) then
            write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
            qrderr = .true.
            go to 999
        end if

        ustr = ufield(3)
        call chreal(nttyo,qrderr,ustr,var)

        if (qrderr) then
            go to 999
        end if

        if (nsk(nrc) .eq. 2) then
            sfcar(nrc) = var
        end if

        if (icount.eq.0) then
            j2 = ilnobl(ureac(nrc))
            write (nttyo,1033) ureac(nrc)(1:j2)
1033 format(/' * Warning - (XCON6/rd6d8) No option was selected',' for the',/7x,'surface area of the reactant ',a,'. The',/7x,'surface area and specific surface area will each be set',/7x,'to zero.')

            sfcar(nrc) = 0.0
            ssfcar(nrc) = 0.0
        else if (icount .gt. 1) then
            j2 = ilnobl(ureac(nrc))
            write (nttyo,1035) ureac(nrc)(1:j2),sfcar(nrc)
1035 format(/' * Warning - (XCON6/rd6d8) Multiple options were',' selected for',/7x,'the surface area of the reactant ',a,'. The',/7x,'surface area will be set to a fixed value of ',1pe10.3,' cm3.')

            ssfcar(nrc) = 0.0
        end if

        ! Read the ratio of active to total surface area (f) from a
        ! two-line header.
        uheadx = '->'
        nfldtx = 4
        call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

        if (qrderr) then
            go to 999
        end if

        ustr24 = ufield(2)(1:24)
        uheadx = 'Surface area factor'
        call locase(ustr24)
        call locase(uheadx)
        j2 = ilnobl(ustr24)
        j3 = ilnobl(uheadx)

        if (ustr24(1:j2) .ne. uheadx(1:j3)) then
            write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
            qrderr = .true.
            go to 999
        end if

        ustr = ufield(3)
        call chreal(nttyo,qrderr,ustr,var)

        if (qrderr) then
            go to 999
        end if

        fkrc(nrc)= var

        ! Loop over forward and backward directions. The forward
        ! direction corresponds to the disappearance (e.g., dissolution,
        ! dissociation) of the associated species. The backward direction
        ! corresponds to its formation (e.g., precipitation, association).
        do jd = 1,2
            if (jd .eq. 1) then
                ux8 = 'forward'
            else
                ux8 = 'backward'
            end if

            ! Read the defining rate law string from a two-line header.
            uheadx = '->'
            nfldtx = 4
            call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

            if (qrderr) then
                go to 999
            end if

            ustr24 = ufield(2)(1:24)

            if (jd .eq. 1) then
                uheadx = 'Forward rate law'
            else
                uheadx = 'Backward rate law'
            end if

            call locase(ustr24)
            call locase(uheadx)
            j2 = ilnobl(ustr24)
            j3 = ilnobl(uheadx)

            if (ustr24(1:j2) .ne. uheadx(1:j3)) then
                write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
                qrderr = .true.
                go to 999
            end if

            ustr = ufield(3)
            call locase(ustr)

            do jnrk = -1,3
                uheadx = urcnrk(jnrk,jd)
                call locase(uheadx)

                if (ustr(1:16) .eq. uheadx(1:16)) then
                    nrk(jd,nrc) = jnrk
                    go to 170
                end if
            end do

165 continue
            j2 = ilnobl(ustr)
            write (nttyo,1075) ustr(1:j2)
1075 format(/" * Error - (XCON6/rd6d8) Don't recognize the",' rate law option string',/7x,' "',a,'". This should',' be one of the strings',/7x,'defined in the urcnrk array.',' The valid strings are:',/)

            do jnrk = -1,3
                j3 = ilnobl(urcnrk(jnrk,jd))
                write (nttyo,1077) urcnrk(jnrk,jd)(1:j3)
1077 format(9x,a)
            end do

            go to 990

170 continue

            if (jd .eq. 1) then
                ! Trap illegal option for the forward direction.
                if (nrk(jd,nrc) .eq. 0) then
                    go to 165
                end if
            end if

            ! Trap options requiring no direct input (e.g., use of the
            ! rate law for the opposite direction, precipitation according
            ! to instantaneous partial equilibrium).
            if (nrk(jd,nrc) .le. 0) then
                go to 195
            end if

            ! Continue with options requiring direct input.
            if (nrk(jd,nrc) .eq. 1) then
                ! Arbitrary kinetics (relative rates, indifferent to time).
                imech(jd,nrc) = 3

                if (imech(jd,nrc) .gt. imchmx) then
                    j2 = ilnobl(ux8)
                    j3 = ilnobl(ureac(nrc))
                    write (nttyo,1960) ux8(1:j2),ureac(nrc)(1:j3),imchmx
1960 format(/' * Error - (XCON6/rd6d8) Have too many rate',/7x,'constants or corresponding mechanisms in the ',a,/7x,'rate law for reactant ',a,'. The code is only',/7x,'dimensioned for ',i2,' rate constants per rate law.',/7x,'Reduce the number of rate constants or increase the',/7x,'dimensioning parameter imchpa.')

                    go to 990
                end if

                do i = 1,3
                    uheadx = '--->'
                    nfldtx = 4
                    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

                    if (qrderr) then
                        go to 999
                    end if

                    ustr24 = ufield(2)(1:24)
                    uheadx = urcrel(i,1)
                    call locase(ustr24)
                    call locase(uheadx)
                    j2 = ilnobl(ustr24)
                    j3 = ilnobl(uheadx)

                    if (ustr24(1:j2) .ne. uheadx(1:j3)) then
                        write (nttyo,1080) ustr24(1:j2),uheadx(1:j3)
1080 format(/" * Error - (XCON6/rd6d8) Have an incorrect",' relative rate law header',/7x,'"',a,'".',' Was expecting the string "',a,'".')

                        go to 990
                    end if

                    ustr = ufield(3)
                    call chreal(nttyo,qrderr,ustr,var)

                    if (qrderr) then
                        go to 999
                    end if

                    rkb(i,jd,nrc) = var
                end do
            else if (nrk(jd,nrc) .eq. 2) then
                ! Transition state theory (TST) rate law. Up to imchmx
                ! parallel mechanisms are allowed.
                ! Loop on mechanisms.
                imh = 0

                do ii = 1,imchmx + 1
                    ! Read a line. If the sub-block for a mechanism has not
                    ! been completely read, this contains the name of a
                    ! mechanism, and a sub-sub-block follows. Otherwise, if
                    ! jd = 1, this line is the first line of the next sub-block
                    ! (backward rate law), or, if jd = 2, the first line of the
                    ! block for the next reactant, if any, else the first line
                    ! of the block following the reactants super-block.
                    uheadx = '--->'
                    nfldtx = 0
                    call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

                    if (qrderr) then
                        go to 999
                    end if

                    ustr = ufield(2)
                    uheadx = 'Mechanism'
                    call locase(ustr)
                    call locase(uheadx)

                    if (ustr(1:9) .ne. uheadx(1:9)) then
                        ! Done reading mechanisms.
                        if (imh .le. 0) then
                            j2 = ilnobl(ux8)
                            j3 = ilnobl(ureac(nrc))
                            write (nttyo,1958) ux8(1:j2),ureac(nrc)(1:j3)
1958 format(/' * Error - (XCON6/rd6d8) Have no',' "mechanisms" specified in the ',a,/7x,'rate law for reactant ',a,'.')

                            go to 990
                        end if

                        ! Set the number of mechanisms for the last
                        ! direction and reactant processed.
                        imech(jd,nrc) = imh

                        ! Back up.
                        backspace ninpts

                        if (jd .eq. 1) then
                            ! Go read the backward rate law data.
                            go to 195
                        else
                            ! Go read the data for the next reactant, if any.
                            go to 235
                        end if
                    end if

                    ! Read the separator line following the mechanism header.
                    nfldtx = 1
                    call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

                    if (qrderr) then
                        go to 999
                    end if

                    ustr24 = ufield(1)(1:24)
                    uheadx = '--------'

                    if (ustr24(1:8) .ne. uheadx(1:8)) then
                        j2 = ilnobl(uline1)
                        j2 = min(j2,70)
                        write (nttyo,1070) uline1(1:j2)
                        qrderr = .true.
                        go to 999
                    end if

                    imh = imh + 1

                    if (imh .gt. imchmx) then
                        j2 = ilnobl(ux8)
                        j3 = ilnobl(ureac(nrc))
                        write (nttyo,1960) ux8(1:j2),ureac(nrc)(1:j3),imchmx
                        go to 990
                    end if

                    ! Read the sigma (stoichiometric correction) factor from
                    ! a two-line header. In some treatments, m appears instead
                    ! (sigma = 1/m). Here sigma is represented by the "csigma"
                    ! variable.
                    uheadx = '----->'
                    nfldtx = 4
                    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

                    if (qrderr) then
                        go to 999
                    end if

                    ustr24 = ufield(2)(1:24)

                    if (jd .eq. 1) then
                        uheadx = 'sigma(i,+,n)'
                    else if (jd .eq. 2) then
                        uheadx = 'sigma(i,-,n)'
                    end if

                    call locase(ustr24)
                    call locase(uheadx)
                    j2 = ilnobl(ustr24)
                    j3 = ilnobl(uheadx)

                    if (ustr24(1:j2) .ne. uheadx(1:j3)) then
                        write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
                        qrderr = .true.
                        go to 999
                    end if

                    ustr = ufield(3)
                    call chreal(nttyo,qrderr,ustr,var)

                    if (qrderr) then
                        go to 999
                    end if

                    csigma(imh,jd,nrc) = var

                    ! Read the rate constant (k, "rkb", mol/cm2/sec) at the
                    ! base or reference temperature ("trkb", see below) from a
                    ! two-line header.
                    uheadx = '----->'
                    nfldtx = 4
                    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

                    if (qrderr) then
                        go to 999
                    end if

                    ustr24 = ufield(2)(1:24)

                    if (jd .eq. 1) then
                        uheadx = 'k(i,+,n) (mol/cm2/sec)'
                    else if (jd .eq. 2) then
                        uheadx = 'k(i,-,n) (mol/cm2/sec)'
                    end if

                    call locase(ustr24)
                    call locase(uheadx)
                    j2 = ilnobl(ustr24)
                    j3 = ilnobl(uheadx)

                    if (ustr24(1:j2) .ne. uheadx(1:j3)) then
                        write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
                        qrderr = .true.
                        go to 999
                    end if

                    ustr = ufield(3)
                    call chreal(nttyo,qrderr,ustr,var)

                    if (qrderr) then
                        go to 999
                    end if

                    rkb(imh,jd,nrc) = var

                    ! Read the base (or reference) temperature ("trkb", C)
                    ! from a two-line header.
                    uheadx = '----->'
                    nfldtx = 4
                    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

                    if (qrderr) then
                        go to 999
                    end if

                    ustr24 = ufield(2)(1:24)
                    uheadx = 'Ref. Temperature (C)'
                    call locase(ustr24)
                    call locase(uheadx)
                    j2 = ilnobl(ustr24)
                    j3 = ilnobl(uheadx)

                    if (ustr24(1:j2) .ne. uheadx(1:j3)) then
                        write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
                        qrderr = .true.
                        go to 999
                    end if

                    ustr = ufield(3)
                    call chreal(nttyo,qrderr,ustr,var)

                    if (qrderr) then
                        go to 999
                    end if

                    trkb(imh,jd,nrc) = var

                    ! Temperature dependence option.
                    icount = 0
                    iact(imh,jd,nrc) = 0

                    ! Read a one-line header.
                    uheadx = '----->'
                    nfldtx = 2
                    call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

                    if (qrderr) then
                        go to 999
                    end if

                    ustr24 = ufield(2)(1:24)
                    uheadx = 'Temperature dependence o'
                    call locase(ustr24)
                    call locase(uheadx)
                    j2 = ilnobl(ustr24)
                    j3 = ilnobl(uheadx)

                    if (ustr24(1:j2) .ne. uheadx(1:j3)) then
                        write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
                        qrderr = .true.
                        go to 999
                    end if

                    ! Read the first option (no temperature dependence) from
                    ! a one-line header.
                    uheadx = '----->'
                    nfldtx = 2
                    call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

                    if (qrderr) then
                        go to 999
                    end if

                    ustr = ufield(2)
                    ustr24 = ustr(5:29)
                    uheadx = '( 0) No temperature depe'
                    call locase(ustr24)
                    call locase(uheadx)
                    j2 = ilnobl(ustr24)
                    j3 = ilnobl(uheadx)

                    if (ustr24(1:j2) .ne. uheadx(1:j3)) then
                        write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
                        qrderr = .true.
                        go to 999
                    end if

                    if (ustr(1:1).ne.'[' .or. ustr(3:3).ne.']') then
                        j2 = ilnobl(ustr)
                        j2 = min(j2,70)
                        write (nttyo,1022) ustr(1:j2)
                        qrderr = .true.
                        go to 999
                    end if

                    ux1 = ustr(2:2)

                    if (ux1.eq.'*' .or. ux1.eq.'x' .or. ux1.eq.'X') then
                        iact(imh,jd,nrc) = 0
                        icount = icount + 1
                    end if

                    ! Read the second option (constant activation energy)
                    ! from a one-line header.
                    uheadx = '----->'
                    nfldtx = 2
                    call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

                    if (qrderr) then
                        go to 999
                    end if

                    ustr = ufield(2)
                    ustr24 = ustr(5:29)
                    uheadx = '( 1) Constant activation'
                    call locase(ustr24)
                    call locase(uheadx)
                    j2 = ilnobl(ustr24)
                    j3 = ilnobl(uheadx)

                    if (ustr24(1:j2) .ne. uheadx(1:j3)) then
                        write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
                        qrderr = .true.
                        go to 999
                    end if

                    if (ustr(1:1).ne.'[' .or. ustr(3:3).ne.']') then
                        j2 = ilnobl(ustr)
                        j2 = min(j2,70)
                        write (nttyo,1022) ustr(1:j2)
                        qrderr = .true.
                        go to 999
                    end if

                    ux1 = ustr(2:2)

                    if (ux1.eq.'*' .or. ux1.eq.'x' .or. ux1.eq.'X') then
                        iact(imh,jd,nrc) = 1
                        icount = icount + 1
                    end if

                    ! Read the associated activation energy from a one-line
                    ! header.
                    uheadx = '----->'
                    nfldtx = 4
                    call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

                    if (qrderr) then
                        go to 999
                    end if

                    ustr24 = ufield(2)(1:24)
                    uheadx = 'Value (kcal/mol)'
                    call locase(ustr24)
                    call locase(uheadx)
                    j2 = ilnobl(ustr24)
                    j3 = ilnobl(uheadx)

                    if (ustr24(1:j2) .ne. uheadx(1:j3)) then
                        write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
                        qrderr = .true.
                        go to 999
                    end if

                    ustr = ufield(3)
                    call chreal(nttyo,qrderr,ustr,var)

                    if (qrderr) then
                        go to 999
                    end if

                    eact(imh,jd,nrc) = var

                    ! Read the third option (constant activation enthalpy)
                    ! from a one-line header.
                    uheadx = '----->'
                    nfldtx = 2
                    call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

                    if (qrderr) then
                        go to 999
                    end if

                    ustr = ufield(2)
                    ustr24 = ustr(5:29)
                    uheadx = '( 2) Constant activation'
                    call locase(ustr24)
                    call locase(uheadx)
                    j2 = ilnobl(ustr24)
                    j3 = ilnobl(uheadx)

                    if (ustr24(1:j2) .ne. uheadx(1:j3)) then
                        write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
                        qrderr = .true.
                        go to 999
                    end if

                    if (ustr(1:1).ne.'[' .or. ustr(3:3).ne.']') then
                        j2 = ilnobl(ustr)
                        j2 = min(j2,70)
                        write (nttyo,1022) ustr(1:j2)
                        qrderr = .true.
                        go to 999
                    end if

                    ux1 = ustr(2:2)

                    if (ux1.eq.'*' .or. ux1.eq.'x' .or. ux1.eq.'X') then
                        iact(imh,jd,nrc) = 2
                        icount = icount + 1
                    end if

                    ! Read the associated activation enthalpy from a two-line
                    ! header.
                    uheadx = '----->'
                    nfldtx = 4
                    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

                    if (qrderr) then
                        go to 999
                    end if

                    ustr24 = ufield(2)(1:24)
                    uheadx = 'Value (kcal/mol)'
                    call locase(ustr24)
                    call locase(uheadx)
                    j2 = ilnobl(ustr24)
                    j3 = ilnobl(uheadx)

                    if (ustr24(1:j2) .ne. uheadx(1:j3)) then
                        write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
                        qrderr = .true.
                        go to 999
                    end if

                    ustr = ufield(3)
                    call chreal(nttyo,qrderr,ustr,var)

                    if (qrderr) then
                        go to 999
                    end if

                    hact(imh,jd,nrc) = var

                    if (icount .ne. 1) then
                        j2 = ilnobl(ux8)
                        j3 = ilnobl(ureac(nrc))

                        if (icount .le. 0) then
                            write (nttyo,1118) ux8(1:j2),imh,ureac(nrc)(1:j3)
1118 format(/' * Warning - (XCON6/rd6d8) No option was',' selected for',/7x,'treating the temperature',' dependence of the ',a,' rate constant',/7x,'for',' term ',i2,' of the rate equation for the reactant',' ',a,'.',/7x,'The temperature dependence will be',' set to zero.')
                        else
                            write (nttyo,1120) ux8(1:j2),imh,ureac(nrc)(1:j3)
1120 format(/' * Warning - (XCON6/rd6d8) Multiple options',' were selected for',/7x,'treating the temperature',' dependence of the ',a,' rate constant',/7x,'for',' term ',i2,' of the rate equation for the reactant',' ',a,'.',/7x,'The temperature dependence will be',' set to zero.')
                        end if

                        iact(imh,jd,nrc) = 0
                        eact(imh,jd,nrc) = 0.
                        hact(imh,jd,nrc) = 0.
                    end if

                    ! Read the species and associated stoichiometric numbers,
                    ! if any, appearing in the kinetic activity product.
                    ! Read the kinetic activity product header from a
                    ! two-line header.
                    uheadx = '----->'
                    nfldtx = 2
                    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

                    if (qrderr) then
                        go to 999
                    end if

                    ustr24 = ufield(2)(1:24)
                    uheadx = 'Kinetic activity product'
                    call locase(ustr24)
                    call locase(uheadx)
                    j2 = ilnobl(ustr24)
                    j3 = ilnobl(uheadx)

                    if (ustr24(1:j2) .ne. uheadx(1:j3)) then
                        write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
                        qrderr = .true.
                        go to 999
                    end if

                    ! Read the first line of the table header for the
                    ! kinetic activity product.
                    uheadx = '------->'
                    nfldtx = 3
                    call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

                    if (qrderr) then
                        go to 999
                    end if

                    ustr24 = ufield(2)(1:24)
                    uheadx = 'Species'
                    call locase(ustr24)
                    call locase(uheadx)
                    j2 = ilnobl(ustr24)
                    j3 = ilnobl(uheadx)

                    if (ustr24(1:j2) .ne. uheadx(1:j3)) then
                        write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
                        qrderr = .true.
                        go to 999
                    end if

                    ! Read the second line of the table header for the
                    ! kinetic activity product.
                    uheadx = '------->'
                    nfldtx = 3
                    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

                    if (qrderr) then
                        go to 999
                    end if

                    ustr24 = ufield(2)(1:24)

                    if (jd .eq. 1) then
                        uheadx = '(udac(j,i,1,n))'
                    else
                        uheadx = '(udac(j,i,2,n))'
                    end if

                    call locase(ustr24)
                    call locase(uheadx)
                    j2 = ilnobl(ustr24)
                    j3 = ilnobl(uheadx)

                    if (ustr24(1:j2) .ne. uheadx(1:j3)) then
                        write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
                        qrderr = .true.
                        go to 999
                    end if

                    ! Loop on species appearing in the kinetic activity
                    ! product.
                    ndt = 0

                    do jj = 1,ndctmx + 1
                        ! Read a line. If the sub-sub-block (for the current
                        ! mechanism) has not been completely read, this contains
                        ! the name of a species. Otherwise, this line is the
                        ! first line of the next sub-block (backward rate law).
                        uheadx = '------->'
                        nfldtx = 0
                        call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

                        if (qrderr) then
                            go to 999
                        end if

                        ustr24 = ufield(2)(1:24)
                        uheadx = 'None'
                        call locase(ustr24)
                        call locase(uheadx)
                        j2 = ilnobl(ustr24)
                        j3 = ilnobl(uheadx)

                        if (ustr24(1:j2) .eq. uheadx(1:j3)) then
                            ! Have no species appearing in the kinetic activity
                            ! product.
                            udac(1,imh,jd,nrc) = ufield(2)(1:24)
                            cdac(1,imh,jd,nrc) = 0.
                            uheadx = '--------'

                            ! Read a separator line.
                            nfldtx = 1
                            call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

                            if (qrderr) then
                                go to 999
                            end if

                            go to 175
                        end if

                        ndt = ndt + 1

                        if (ndt .gt. ndctmx) then
                            j2 = ilnobl(ux8)
                            j3 = ilnobl(ureac(nrc))
                            write (nttyo,1130) i,ux8(1:j2),ureac(nrc)(1:j3),ndctmx
1130 format(/' * Error - (XCON6/rd6d8) Have too many',/7x,'species in the activity product in term ',i2,/7x,'of the ',a,' direction rate law for reactant',/7x,a,'. The code is only dimensioned for ',i3,/7x,'such species. Reduce the number of such species',/7x,'or increase the dimensioning parameter ndctpa.')

                            go to 990
                        end if

                        udac(ndt,imh,jd,nrc) = ufield(2)(1:24)
                        ustr = ufield(3)
                        call chreal(nttyo,qrderr,ustr,var)

                        if (qrderr) then
                            go to 999
                        end if

                        cdac(ndt,imh,jd,nrc) = var

                        ! Read a line. If the sub-sub-sub-block (for the current
                        ! activity product) has been completely read, this line
                        ! is a separator line. Otherwise, it contains data for
                        ! another species in the current kinetic activity product.
                        nfldtx = 0
                        call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

                        if (qrderr) then
                            go to 999
                        end if

                        ustr = ufield(1)
                        uheadx = '--------'

                        if (ustr(1:8) .eq. uheadx(1:8)) then
                            go to 175
                        end if
                    end do

175 continue
                    ndact(imh,jd,nrc) = ndt
                end do
            else if (nrk(jd,nrc) .eq. 3) then
                ! Linear rate law.
                ! Loop on mechanisms. However, only one "mechanism" is
                ! currently allowed.
                imh = 0

                do ii = 1,imchmx + 1
                    ! Read a line. If the sub-block for a mechanism has not
                    ! been completely read, this contains the name of a
                    ! mechanism, and a sub-sub-block follows. Otherwise, if
                    ! jd = 1, this line is the first line of the next sub-block
                    ! (backward rate law), or, if jd = 2, the first line of the
                    ! block for the next reactant, if any, else the first line
                    ! of the block following the reactants super-block.
                    uheadx = '--->'
                    nfldtx = 0
                    call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

                    if (qrderr) then
                        go to 999
                    end if

                    ustr = ufield(2)
                    uheadx = 'Mechanism'
                    call locase(ustr)
                    call locase(uheadx)

                    if (ustr(1:9) .ne. uheadx(1:9)) then
                        ! Done reading mechanisms.
                        if (imh .le. 0) then
                            j2 = ilnobl(ux8)
                            j3 = ilnobl(ureac(nrc))
                            write (nttyo,1958) ux8(1:j2),ureac(nrc)(1:j3)
                            go to 990
                        end if

                        ! Set the number of mechanisms for the last
                        ! direction and reactant processed.
                        imech(jd,nrc) = imh

                        ! Back up.
                        backspace ninpts

                        if (jd .eq. 1) then
                            ! Go read the backward rate law data.
                            go to 195
                        else
                            ! Go read the data for the next reactant, if any.
                            go to 235
                        end if
                    end if

                    ! Read the separator line following the mechanism header.
                    nfldtx = 1
                    call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

                    if (qrderr) then
                        go to 999
                    end if

                    ustr24 = ufield(1)(1:24)
                    uheadx = '--------'

                    if (ustr24(1:8) .ne. uheadx(1:8)) then
                        j2 = ilnobl(uline1)
                        j2 = min(j2,70)
                        write (nttyo,1070) uline1(1:j2)
                        qrderr = .true.
                        go to 999
                    end if

                    imh = imh + 1

                    if (imh .gt. imchmx) then
                        j2 = ilnobl(ux8)
                        j3 = ilnobl(ureac(nrc))
                        write (nttyo,1960) ux8(1:j2),ureac(nrc)(1:j3),imchmx
                        go to 990
                    end if

                    if (imh .gt. 1) then
                        j2 = ilnobl(ux8)
                        j3 = ilnobl(ureac(nrc))
                        write (nttyo,1964) ux8(1:j2),ureac(nrc)(1:j3)
1964 format(/' * Error - (XCON6/rd6d8) Have more than the',' allowed one "mechanism"',/7x,'in the ',a,' rate','law for reactant ',a,'.')

                        go to 990
                    end if

                    ! Read the rate constant (k, "rkb", mol/cm2/sec) at the
                    ! base or reference temperature ("trkb", see below) from a
                    ! two-line header.
                    uheadx = '----->'
                    nfldtx = 4
                    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

                    if (qrderr) then
                        go to 999
                    end if

                    ustr24 = ufield(2)(1:24)

                    if (jd .eq. 1) then
                        uheadx = 'k(i,+,n) (mol/cm2/sec)'
                    else if (jd .eq. 2) then
                        uheadx = 'k(i,-,n) (mol/cm2/sec)'
                    end if

                    call locase(ustr24)
                    call locase(uheadx)
                    j2 = ilnobl(ustr24)
                    j3 = ilnobl(uheadx)

                    if (ustr24(1:j2) .ne. uheadx(1:j3)) then
                        write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
                        qrderr = .true.
                        go to 999
                    end if

                    ustr = ufield(3)
                    call chreal(nttyo,qrderr,ustr,var)

                    if (qrderr) then
                        go to 999
                    end if

                    rkb(imh,jd,nrc) = var

                    ! Read the base (or reference) temperature ("trkb", C)
                    ! from a two-line header.
                    uheadx = '----->'
                    nfldtx = 4
                    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

                    if (qrderr) then
                        go to 999
                    end if

                    ustr24 = ufield(2)(1:24)
                    uheadx = 'Ref. Temperature (C)'
                    call locase(ustr24)
                    call locase(uheadx)
                    j2 = ilnobl(ustr24)
                    j3 = ilnobl(uheadx)

                    if (ustr24(1:j2) .ne. uheadx(1:j3)) then
                        write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
                        qrderr = .true.
                        go to 999
                    end if

                    ustr = ufield(3)
                    call chreal(nttyo,qrderr,ustr,var)

                    if (qrderr) then
                        go to 999
                    end if

                    trkb(imh,jd,nrc) = var

                    ! Temperature dependence option.
                    icount = 0
                    iact(imh,jd,nrc) = 0

                    ! Read a one-line header.
                    uheadx = '----->'
                    nfldtx = 2
                    call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

                    if (qrderr) then
                        go to 999
                    end if

                    ustr24 = ufield(2)(1:24)
                    uheadx = 'Temperature dependence o'
                    call locase(ustr24)
                    call locase(uheadx)
                    j2 = ilnobl(ustr24)
                    j3 = ilnobl(uheadx)

                    if (ustr24(1:j2) .ne. uheadx(1:j3)) then
                        write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
                        qrderr = .true.
                        go to 999
                    end if

                    ! Read the first option (no temperature dependence) from
                    ! a one-line header.
                    uheadx = '----->'
                    nfldtx = 2
                    call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

                    if (qrderr) then
                        go to 999
                    end if

                    ustr = ufield(2)
                    ustr24 = ustr(5:29)
                    uheadx = '( 0) No temperature depe'
                    call locase(ustr24)
                    call locase(uheadx)
                    j2 = ilnobl(ustr24)
                    j3 = ilnobl(uheadx)

                    if (ustr24(1:j2) .ne. uheadx(1:j3)) then
                        write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
                        qrderr = .true.
                        go to 999
                    end if

                    if (ustr(1:1).ne.'[' .or. ustr(3:3).ne.']') then
                        j2 = ilnobl(ustr)
                        j2 = min(j2,70)
                        write (nttyo,1022) ustr(1:j2)
                        qrderr = .true.
                        go to 999
                    end if

                    ux1 = ustr(2:2)

                    if (ux1.eq.'*' .or. ux1.eq.'x' .or. ux1.eq.'X') then
                        iact(imh,jd,nrc) = 0
                        icount = icount + 1
                    end if

                    ! Read the second option (constant activation energy)
                    ! from a one-line header.
                    uheadx = '----->'
                    nfldtx = 2
                    call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

                    if (qrderr) then
                        go to 999
                    end if

                    ustr = ufield(2)
                    ustr24 = ustr(5:29)
                    uheadx = '( 1) Constant activation'
                    call locase(ustr24)
                    call locase(uheadx)
                    j2 = ilnobl(ustr24)
                    j3 = ilnobl(uheadx)

                    if (ustr24(1:j2) .ne. uheadx(1:j3)) then
                        write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
                        qrderr = .true.
                        go to 999
                    end if

                    if (ustr(1:1).ne.'[' .or. ustr(3:3).ne.']') then
                        j2 = ilnobl(ustr)
                        j2 = min(j2,70)
                        write (nttyo,1022) ustr(1:j2)
                        qrderr = .true.
                        go to 999
                    end if

                    ux1 = ustr(2:2)

                    if (ux1.eq.'*' .or. ux1.eq.'x' .or. ux1.eq.'X') then
                        iact(imh,jd,nrc) = 1
                        icount = icount + 1
                    end if

                    ! Read the associated activation energy from a one-line
                    ! header.
                    uheadx = '----->'
                    nfldtx = 4
                    call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

                    if (qrderr) then
                        go to 999
                    end if

                    ustr24 = ufield(2)(1:24)
                    uheadx = 'Value (kcal/mol)'
                    call locase(ustr24)
                    call locase(uheadx)
                    j2 = ilnobl(ustr24)
                    j3 = ilnobl(uheadx)

                    if (ustr24(1:j2) .ne. uheadx(1:j3)) then
                        write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
                        qrderr = .true.
                        go to 999
                    end if

                    ustr = ufield(3)
                    call chreal(nttyo,qrderr,ustr,var)

                    if (qrderr) then
                        go to 999
                    end if

                    eact(imh,jd,nrc) = var

                    ! Read the third option (constant activation enthalpy)
                    ! from a one-line header.
                    uheadx = '----->'
                    nfldtx = 2
                    call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

                    if (qrderr) then
                        go to 999
                    end if

                    ustr = ufield(2)
                    ustr24 = ustr(5:29)
                    uheadx = '( 2) Constant activation'
                    call locase(ustr24)
                    call locase(uheadx)
                    j2 = ilnobl(ustr24)
                    j3 = ilnobl(uheadx)

                    if (ustr24(1:j2) .ne. uheadx(1:j3)) then
                        write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
                        qrderr = .true.
                        go to 999
                    end if

                    if (ustr(1:1).ne.'[' .or. ustr(3:3).ne.']') then
                        j2 = ilnobl(ustr)
                        j2 = min(j2,70)
                        write (nttyo,1022) ustr(1:j2)
                        qrderr = .true.
                        go to 999
                    end if

                    ux1 = ustr(2:2)

                    if (ux1.eq.'*' .or. ux1.eq.'x' .or. ux1.eq.'X') then
                        iact(imh,jd,nrc) = 2
                        icount = icount + 1
                    end if

                    ! Read the associated activation enthalpy from a two-line
                    ! header.
                    uheadx = '----->'
                    nfldtx = 4
                    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

                    if (qrderr) then
                        go to 999
                    end if

                    ustr24 = ufield(2)(1:24)
                    uheadx = 'Value (kcal/mol)'
                    call locase(ustr24)
                    call locase(uheadx)
                    j2 = ilnobl(ustr24)
                    j3 = ilnobl(uheadx)

                    if (ustr24(1:j2) .ne. uheadx(1:j3)) then
                        write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
                        qrderr = .true.
                        go to 999
                    end if

                    ustr = ufield(3)
                    call chreal(nttyo,qrderr,ustr,var)

                    if (qrderr) then
                        go to 999
                    end if

                    hact(imh,jd,nrc) = var

                    if (icount .ne. 1) then
                        j3 = ilnobl(ux8)
                        j2 = ilnobl(ureac(nrc))

                        if (icount .le. 0) then
                            write (nttyo,1118) ux8(1:j3),imh,ureac(nrc)(1:j2)
                        else
                            write (nttyo,1120) ux8(1:j3),imh,ureac(nrc)(1:j2)
                        end if

                        iact(imh,jd,nrc) = 0
                        eact(imh,jd,nrc) = 0.
                        hact(imh,jd,nrc) = 0.
                    end if
                end do
            else
                ! Error, unknown rate law code.
                j3 = ilnobl(ux8)
                j2 = ilnobl(ureac(nrc))
                write (nttyo,1952) ux8(1:j3),nrk(jd,nrc),ureac(nrc)(1:j2)
1952 format(/' * Error - (XCON6/rd6d8) Programming error trap:',' The ',a,/7x,'rate law code has an unrecognized value',' of ',i2,' for reactant',/7x,a,'.')

                go to 990
            end if

195 continue
        end do

        ! End of the loop on the forward and backward directions of the
        ! irreversible reaction corresponding to a reactant.
235 continue
    end do

    ! End of the loop on reactants.
240 continue
    nrct = nrc
    nsrt = nsr
    nxrt = nxr
    nert = ner

    ! Starting, minimum, and maximum values of key run parameters.
    ! Read superblock header.
    ! Read the data from a two-line header.
    uheadx = 'Starting, minimum, and maximum values of key run parameters.'
    nfldtx = 1
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    ! Starting value of Xi.
    xistti = 0.

    ! Read the data from a two-line header.
    uheadx = 'Starting Xi value'
    nfldtx = 3
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    ustr = ufield(2)
    call chreal(nttyo,qrderr,ustr,var)

    if (qrderr) then
        go to 999
    end if

    xistti = var

    ! Maximum value of Xi.
    ximaxi = 0.

    ! Read the data from a two-line header.
    uheadx = 'Maximum Xi value'
    nfldtx = 3
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    ustr = ufield(2)
    call chreal(nttyo,qrderr,ustr,var)

    if (qrderr) then
        go to 999
    end if

    ximaxi = var

    ! Starting value of time.
    tistti = 0.

    ! Read the data from a two-line header.
    uheadx = 'Starting time (seconds)'
    nfldtx = 3
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    ustr = ufield(2)
    call chreal(nttyo,qrderr,ustr,var)

    if (qrderr) then
        go to 999
    end if

    tistti = var

    ! Maximum value of time.
    timmxi = 0.

    ! Read the data from a two-line header.
    uheadx = 'Maximum time (seconds)'
    nfldtx = 3
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    ustr = ufield(2)
    call chreal(nttyo,qrderr,ustr,var)

    if (qrderr) then
        go to 999
    end if

    timmxi = var

    ! Minimum value of pH.
    phmini = 0.

    ! Read the data from a two-line header.
    uheadx = 'Minimum value of pH'
    nfldtx = 3
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    ustr = ufield(2)
    call chreal(nttyo,qrderr,ustr,var)

    if (qrderr) then
        go to 999
    end if

    phmini = var

    ! Maximum value of pH.
    phmaxi = 0.

    ! Read the data from a two-line header.
    uheadx = 'Maximum value of pH'
    nfldtx = 3
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    ustr = ufield(2)
    call chreal(nttyo,qrderr,ustr,var)

    if (qrderr) then
        go to 999
    end if

    phmaxi = var

    ! Minimum value of Eh (v).
    ehmini = 0.

    ! Read the data from a two-line header.
    uheadx = 'Minimum value of Eh (v)'
    nfldtx = 3
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    ustr = ufield(2)
    call chreal(nttyo,qrderr,ustr,var)

    if (qrderr) then
        go to 999
    end if

    ehmini = var

    ! Maximum value of Eh (v).
    ehmaxi = 0.

    ! Read the data from a two-line header.
    uheadx = 'Maximum value of Eh (v)'
    nfldtx = 3
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    ustr = ufield(2)
    call chreal(nttyo,qrderr,ustr,var)

    if (qrderr) then
        go to 999
    end if

    ehmaxi = var

    ! Minimum value of log fO2.
    o2mini = 0.

    ! Read the data from a two-line header.
    uheadx = 'Minimum value of log fO2'
    nfldtx = 3
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    ustr = ufield(2)
    call chreal(nttyo,qrderr,ustr,var)

    if (qrderr) then
        go to 999
    end if

    o2mini = var

    ! Maximum value of log fO2.
    o2maxi = 0.

    ! Read the data from a two-line header.
    uheadx = 'Maximum value of log fO2'
    nfldtx = 3
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    ustr = ufield(2)
    call chreal(nttyo,qrderr,ustr,var)

    if (qrderr) then
        go to 999
    end if

    o2maxi = var

    ! Minimum value of aw.
    awmini = 0.

    ! Read the data from a two-line header.
    uheadx = 'Minimum value of aw'
    nfldtx = 3
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    ustr = ufield(2)
    call chreal(nttyo,qrderr,ustr,var)

    if (qrderr) then
        go to 999
    end if

    awmini = var

    ! Maximum value of aw.
    awmaxi = 0.

    ! Read the data from a two-line header.
    uheadx = 'Maximum value of aw'
    nfldtx = 3
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    ustr = ufield(2)
    call chreal(nttyo,qrderr,ustr,var)

    if (qrderr) then
        go to 999
    end if

    awmaxi = var

    ! Maximum number of steps.
    kstpmx = 0.

    ! Read the data from a two-line header.
    uheadx = 'Maximum number of steps'
    nfldtx = 3
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    ustr = ufield(2)
    call chrint(ivar,nttyo,qrderr,ustr)

    if (qrderr) then
        go to 999
    end if

    kstpmx = ivar

    ! Print interval parameters. Read superblock header.
    ! Read the data from a two-line header.
    uheadx = 'Print interval parameters.'
    nfldtx = 1
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    ! Xi print interval.
    dlxprn = 0.

    ! Read the data from a two-line header.
    uheadx = 'Xi print interval'
    nfldtx = 3
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    ustr = ufield(2)
    call chreal(nttyo,qrderr,ustr,var)

    if (qrderr) then
        go to 999
    end if

    dlxprn = var

    ! Log Xi print interval.
    dlxprl = 0.

    ! Read the data from a two-line header.
    uheadx = 'Log Xi print interval'
    nfldtx = 3
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    ustr = ufield(2)
    call chreal(nttyo,qrderr,ustr,var)

    if (qrderr) then
        go to 999
    end if

    dlxprl = var

    ! Time print interval.
    dltprn = 0.

    ! Read the data from a two-line header.
    uheadx = 'Time print interval'
    nfldtx = 3
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    ustr = ufield(2)
    call chreal(nttyo,qrderr,ustr,var)

    if (qrderr) then
        go to 999
    end if

    dltprn = var

    ! Log time print interval.
    dltprl = 0.

    ! Read the data from a two-line header.
    uheadx = 'Log time print interval'
    nfldtx = 3
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    ustr = ufield(2)
    call chreal(nttyo,qrderr,ustr,var)

    if (qrderr) then
        go to 999
    end if

    dltprl = var

    ! pH print interval.
    dlhprn = 0.

    ! Read the data from a two-line header.
    uheadx = 'pH print interval'
    nfldtx = 3
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    ustr = ufield(2)
    call chreal(nttyo,qrderr,ustr,var)

    if (qrderr) then
        go to 999
    end if

    dlhprn = var

    ! Eh (v) print interval.
    dleprn = 0.

    ! Read the data from a two-line header.
    uheadx = 'Eh (v) print interval'
    nfldtx = 3
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    ustr = ufield(2)
    call chreal(nttyo,qrderr,ustr,var)

    if (qrderr) then
        go to 999
    end if

    dleprn = var

    ! Log fO2 print interval.
    dloprn = 0.

    ! Read the data from a two-line header.
    uheadx = 'Log fO2 print interval'
    nfldtx = 3
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    ustr = ufield(2)
    call chreal(nttyo,qrderr,ustr,var)

    if (qrderr) then
        go to 999
    end if

    dloprn = var

    ! Activity of water print interval.
    dlaprn = 0.

    ! Read the data from a two-line header.
    uheadx = 'aw print interval'
    nfldtx = 3
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    ustr = ufield(2)
    call chreal(nttyo,qrderr,ustr,var)

    if (qrderr) then
        go to 999
    end if

    dlaprn = var

    ! Steps print interval.
    ksppmx = 0.

    ! Read the data from a two-line header.
    uheadx = 'Steps print interval'
    nfldtx = 3
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    ustr = ufield(2)
    call chrint(ivar,nttyo,qrderr,ustr)

    if (qrderr) then
        go to 999
    end if

    ksppmx = ivar

    ! Plot interval parameters. Read superblock header.
    ! Read the data from a two-line header.
    uheadx = 'Plot interval parameters.'
    nfldtx = 1
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    ! Xi plot interval.
    dlxplo = 0.

    ! Read the data from a two-line header.
    uheadx = 'Xi plot interval'
    nfldtx = 3
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    ustr = ufield(2)
    call chreal(nttyo,qrderr,ustr,var)

    if (qrderr) then
        go to 999
    end if

    dlxplo = var

    ! Log Xi plot interval.
    dlxpll = 0.

    ! Read the data from a two-line header.
    uheadx = 'Log Xi plot interval'
    nfldtx = 3
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    ustr = ufield(2)
    call chreal(nttyo,qrderr,ustr,var)

    if (qrderr) then
        go to 999
    end if

    dlxpll = var

    ! Time plot interval.
    dltplo = 0.

    ! Read the data from a two-line header.
    uheadx = 'Time plot interval'
    nfldtx = 3
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    ustr = ufield(2)
    call chreal(nttyo,qrderr,ustr,var)

    if (qrderr) then
        go to 999
    end if

    dltplo = var

    ! Log time plot interval.
    dltpll = 0.

    ! Read the data from a two-line header.
    uheadx = 'Log time plot interval'
    nfldtx = 3
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    ustr = ufield(2)
    call chreal(nttyo,qrderr,ustr,var)

    if (qrderr) then
        go to 999
    end if

    dltpll = var

    ! pH plot interval.
    dlhplo = 0.

    ! Read the data from a two-line header.
    uheadx = 'pH plot interval'
    nfldtx = 3
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    ustr = ufield(2)
    call chreal(nttyo,qrderr,ustr,var)

    if (qrderr) then
        go to 999
    end if

    dlhplo = var

    ! Eh (v) plot interval.
    dleplo = 0.

    ! Read the data from a two-line header.
    uheadx = 'Eh (v) plot interval'
    nfldtx = 3
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    ustr = ufield(2)
    call chreal(nttyo,qrderr,ustr,var)

    if (qrderr) then
        go to 999
    end if

    dleplo = var

    ! Log fO2 plot interval.
    dloplo = 0.

    ! Read the data from a two-line header.
    uheadx = 'Log fO2 plot interval'
    nfldtx = 3
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    ustr = ufield(2)
    call chreal(nttyo,qrderr,ustr,var)

    if (qrderr) then
        go to 999
    end if

    dloplo = var

    ! Activity of water plot interval.
    dlaplo = 0.

    ! Read the data from a two-line header.
    uheadx = 'aw plot interval'
    nfldtx = 3
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    ustr = ufield(2)
    call chreal(nttyo,qrderr,ustr,var)

    if (qrderr) then
        go to 999
    end if

    dlaplo = var

    ! Steps plot interval.
    ksplmx = 0.

    ! Read the data from a two-line header.
    uheadx = 'Steps plot interval'
    nfldtx = 3
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    ustr = ufield(2)
    call chrint(ivar,nttyo,qrderr,ustr)

    if (qrderr) then
        go to 999
    end if

    ksplmx = ivar

    ! Iopt Model Option Switches.
    ! Note: iopt(1) = iopt1, etc.
    uheadx = 'Iopt Model Option Switches ("( 0)" marks default choices)'
    nfldtx = 1
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

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

        k1 = index(uheadx,'(')
        k2 = index(uheadx,')')
        k3 = index(uheadx,'- ')

        if (k1.le.0 .or. k2.le.0 .or. k2.lt.(k1 + 2) .or. k3.le.0 .or.  k3.le.k2) then
            write (nttyo,1607) uheadx(1:j2)
1607 format(/' * Error - (XCON6/rd6d8) The iopt option switch',' title line',/7x,'"',a,'"',/7x,'read from the'," input file isn't in the required format.")

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
            write (nttyo,1610) uheadx(1:j2),n
1610 format(/' * Error - (XCON6/rd6d8) The iopt option switch',' title line',/7x,'"',a,'"',/7x,'read from the',' input file references an option switch index',/7x,'of ',i3,', which is out of range.')

            go to 990
        end if

        ! Check the full title string.
        ustr = uopttx(n)
        j3 = ilnobl(ustr)

        if (index(uheadx(1:j2),ustr(1:j3)) .le. 0) then
            write (nttyo,1612) uheadx(1:j2),ustr(1:j3)
1612 format(/' * Error - (XCON6/rd6d8) The iopt option switch',' title string',/7x,'"',a,'"',/7x,"read from the input file doesn't contain the",' matching defined string',/7x,'"',a,'".')

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

            uheadx = ufield(1)
            j2 = ilnobl(uheadx)

            if (uheadx(1:8) .eq. '--------') then
                go to 620
            end if

            if (jj .gt. jptxpa) then
                j3 = ilnobl(uopttx(n))
                write (nttyo,1630) uopttx(n)(1:j3),jptxpa
1630 format(/' * Error - (XCON6/rd6d8) Have too many option',/7x,'choice lines for the iopt option switch whose title',' string is',/7x,'"',a,'".',/7x,' The code is only',' dimensioned for ',i3,' such lines. Reduce the',/7x,'number of such lines or increase the dimensioning',' parameter jptxpa.')

                go to 990
            end if

            k1 = index(uheadx,'[')
            k2 = index(uheadx,']')
            k3 = index(uheadx,'(')
            k4 = index(uheadx,')')

            if (k1.le.0 .or. k2.le.0 .or. k2.lt.(k1 + 2) .or.    k3.le.0 .or. k4.le.0 .or. k4.lt.(k3 + 2)) then
                write (nttyo,1640) uheadx(1:j2)
1640 format(/' * Error - (XCON6/rd6d8) The following iopt',/7x,'option switch choice line read from the input file',/7x,"isn't in the required format:",/7x,'"',a,'"')

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

            write (nttyo,1650) uheadx(1:j2)
1650 format(/' * Error - (XCON6/rd6d8) The iopt option switch',/7x,'choice line',/7x,'"',a,'"',/7x,'read from the input file references an out-of-range',' option',/7x,'choice index.')

            go to 990

600 continue

            ! Check the full option choice string.
            ustr = uoptox(j,n)
            j3 = ilnobl(ustr)

            if (index(uheadx(1:j2),ustr(1:j3)) .le. 0) then
                j4 = ilnobl(uopttx(n))
                write (nttyo,1662) uheadx(1:j2),ustr(1:j3),uopttx(n)(1:j4)
1662 format(/' * Error - (XCON6/rd6d8) The iopt option switch',' choice line',/7x,'"',a,'"',/7x," read from the input file doesn't contain the",' matching defined string',/7x,'"',a,'".',/7x,'This line belongs to the option switch whose title',' string is',/7x,'"',a,'".')

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
                    write (nttyo,1670) ustr24(1:j3),uheadx(1:j2),uopttx(n)(1:j4)
1670 format(/" * Error - (XCON6/rd6d8) Don't recognize the",' string "',a,'"',/7x,'that appears on the iopt',' option switch choice line',/7x,'"',a,'"',/7x,'read from the input file. An option choice should',' be chosen by',/7x,'placing a "*", "x", or "X" in the',' checkbox ("[ ]"). This',/7x,'choice line belongs to',' the option switch whose title string is',/7x,'"',a,'".')

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
            write (nttyo,1680) uopttx(n)(1:j2)
1680 format(/' * Warning - (XCON6/rd6d8) No option choice was',' checked on the input file',/7x,'for the iopt option',' switch whose title string is',/7x,'"',a,'".')

            do j = 1,jptxpa
                ival = ioptox(j,n)

                if (ival .eq. 0) then
                    go to 630
                end if
            end do

            write (nttyo,1690)
1690 format(/7x,'A default value of 0 has been applied, but no',/7x,'matching option choice string is defined.')

            go to 640

630 continue
            j3 = ilnobl(uoptox(j,n))
            write (nttyo,1692) uoptox(j,n)(1:j3)
1692 format(/7x,'A default value of 0 has been applied. The',' matching string is',/7x,'"',a,'".')

640 continue
        end if

        if (nmark .gt. 1) then
            j2 = ilnobl(uopttx(n))
            j = jlast
            j3 = ilnobl(uoptox(j,n))
            write (nttyo,1694) uopttx(n)(1:j2),uoptox(j,n)(1:j3)
1694 format(/' * Warning - (XCON6/rd6d8) More than one option',' choice was checked',/7x,'on the input file for the iopt',' option switch whose title string is',/7x,'"',a,'".',/7x,'The last choice checked will be used. The',' matching string is',/7x,'"',a,'".')
        end if
    end do

650 continue

    ! Iopr Print Option Switches.
    ! Note: iopr(1) = iopr1, etc.
    uheadx = 'Iopr Print Option Switches ("( 0)" marks default choices)'
    nfldtx = 1
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

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

        k1 = index(uheadx,'(')
        k2 = index(uheadx,')')
        k3 = index(uheadx,'- ')

        if (k1.le.0 .or. k2.le.0 .or. k2.lt.(k1 + 2) .or. k3.le.0 .or.  k3.le.k2) then
            write (nttyo,1800) uheadx(1:j2)
1800 format(/' * Error - (XCON6/rd6d8) The iopr option switch',' title line',/7x,'"',a,'"',/7x,'read from the'," input file isn't in the required format.")

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
            write (nttyo,1810) uheadx(1:j2),n
1810 format(/' * Error - (XCON6/rd6d8) The iopr option switch',' title line',/7x,'"',a,'"',/7x,'read from the',' input file references an option switch index',/7x,'of ',i3,', which is out of range.')

            go to 990
        end if

        ! Check the full title string.
        ustr = uoprtx(n)
        j3 = ilnobl(ustr)

        if (index(uheadx(1:j2),ustr(1:j3)) .le. 0) then
            write (nttyo,1812) uheadx(1:j2),ustr(1:j3)
1812 format(/' * Error - (XCON6/rd6d8) The iopr option switch',' title string',/7x,'"',a,'"',/7x,"read from the input file doesn't contain the",' matching defined string',/7x,'"',a,'".')

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

            uheadx = ufield(1)
            j2 = ilnobl(uheadx)

            if (uheadx(1:8) .eq. '--------') then
                go to 820
            end if

            if (jj .gt. jprxpa) then
                j3 = ilnobl(uoprtx(n))
                write (nttyo,1830) uoprtx(n)(1:j3),jprxpa
1830 format(/' * Error - (XCON6/rd6d8) Have too many option',/7x,'choice lines for the iopr option switch whose title',' string is',/7x,'"',a,'".',/7x,' The code is only',' dimensioned for ',i3,' such lines. Reduce the',/7x,'number of such lines or increase the dimensioning',' parameter jprxpa.')

                go to 990
            end if

            k1 = index(uheadx,'[')
            k2 = index(uheadx,']')
            k3 = index(uheadx,'(')
            k4 = index(uheadx,')')

            if (k1.le.0 .or. k2.le.0 .or. k2.lt.(k1 + 2) .or.    k3.le.0 .or. k4.le.0 .or. k4.lt.(k3 + 2)) then
                write (nttyo,1840) uheadx(1:j2)
1840 format(/' * Error - (XCON6/rd6d8) The following iopr',/7x,'option switch choice line read from the input file',/7x,"isn't in the required format:",/7x,'"',a,'"')

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

            write (nttyo,1850) uheadx(1:j2)
1850 format(/' * Error - (XCON6/rd6d8) The iopr option switch',/7x,'choice line',/7x,'"',a,'"',/7x,'read from the input file references an out-of-range',' option',/7x,'choice index.')

            go to 990

800 continue

            ! Check the full option choice string.
            ustr = uoprox(j,n)
            j3 = ilnobl(ustr)

            if (index(uheadx(1:j2),ustr(1:j3)) .le. 0) then
                j4 = ilnobl(uoprtx(n))
                write (nttyo,1860) uheadx(1:j2),ustr(1:j3),uoprtx(n)(1:j4)
1860 format(/' * Error - (XCON6/rd6d8) The iopr option switch',' choice line',/7x,'"',a,'"',/7x," read from the input file doesn't contain the",' matching defined string',/7x,'"',a,'".',/7x,'This line belongs to the option switch whose title',' string is',/7x,'"',a,'".')

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
                    write (nttyo,1870) ustr24(1:j3),uheadx(1:j2),uoprtx(n)(1:j4)
1870 format(/" * Error - (XCON6/rd6d8) Don't recognize the",' string "',a,'"',/7x,'that appears on the iopr',' option switch choice line',/7x,'"',a,'"',/7x,'read from the input file. An option choice should',' be chosen by',/7x,'placing a "*", "x", or "X" in the',' checkbox ("[ ]"). This',/7x,'choice line belongs to',' the option switch whose title string is',/7x,'"',a,'".')

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
            write (nttyo,1880) uoprtx(n)(1:j2)
1880 format(/' * Warning - (XCON6/rd6d8) No option choice was',' checked on the input file',/7x,'for the iopr option',' switch whose title string is',/7x,'"',a,'".')

            do j = 1,jprxpa
                ival = ioprox(j,n)

                if (ival .eq. 0) then
                    go to 830
                end if
            end do

            write (nttyo,1890)
1890 format(/7x,'A default value of 0 has been applied, but no',/7x,'matching option choice string is defined.')

            go to 840

830 continue
            j3 = ilnobl(uoprox(j,n))
            write (nttyo,1892) uoprox(j,n)(1:j3)
1892 format(/7x,'A default value of 0 has been applied. The',' matching string is',/7x,'"',a,'".')

840 continue
        end if

        if (nmark .gt. 1) then
            j2 = ilnobl(uoprtx(n))
            j = jlast
            j3 = ilnobl(uoprox(j,n))
            write (nttyo,1894) uoprtx(n)(1:j2),uoprox(j,n)(1:j3)
1894 format(/' * Warning - (XCON6/rd6d8) More than one option',' choice was checked',/7x,'on the input file for the iopr',' option switch whose title string is',/7x,'"',a,'".',/7x,'The last choice checked will be used. The',' matching string is',/7x,'"',a,'".')
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

    ! Loop on option switches.
    do nn = 1,nodbmx
        ! Read the option title string from a one-line header.
        nfldtx = 0
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

        k1 = index(uheadx,'(')
        k2 = index(uheadx,')')
        k3 = index(uheadx,'- ')

        if (k1.le.0 .or. k2.le.0 .or. k2.lt.(k1 + 2) .or. k3.le.0 .or.  k3.le.k2) then
            write (nttyo,1900) uheadx(1:j2)
1900 format(/' * Error - (XCON6/rd6d8) The iodb option switch',' title line',/7x,'"',a,'"',/7x,'read from the'," input file isn't in the required format.")

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
            write (nttyo,1910) uheadx(1:j2),n
1910 format(/' * Error - (XCON6/rd6d8) The iodb option switch',' title line',/7x,'"',a,'"',/7x,'read from the',' input file references an option switch index',/7x,'of ',i3,', which is out of range.')

            go to 990
        end if

        ! Check the full title string.
        ustr = uodbtx(n)
        j3 = ilnobl(ustr)

        if (index(uheadx(1:j2),ustr(1:j3)) .le. 0) then
            write (nttyo,1920) uheadx(1:j2),ustr(1:j3)
1920 format(/' * Error - (XCON6/rd6d8) The iodb option switch',' title string',/7x,'"',a,'"',/7x,"read from the input file doesn't contain the",' matching defined string',/7x,'"',a,'".')

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

            uheadx = ufield(1)
            j2 = ilnobl(uheadx)

            if (uheadx(1:8) .eq. '--------') then
                go to 920
            end if

            if (jj .gt. jdbxpa) then
                j3 = ilnobl(uodbtx(n))
                write (nttyo,1922) uodbtx(n)(1:j3),jdbxpa
1922 format(/' * Error - (XCON6/rd6d8) Have too many option',/7x,'choice lines for the iodb option switch whose title',' string is',/7x,'"',a,'".',/7x,' The code is only',' dimensioned for ',i3,' such lines. Reduce the',/7x,'number of such lines or increase the dimensioning',' parameter jdbxpa.')

                go to 990
            end if

            k1 = index(uheadx,'[')
            k2 = index(uheadx,']')
            k3 = index(uheadx,'(')
            k4 = index(uheadx,')')

            if (k1.le.0 .or. k2.le.0 .or. k2.lt.(k1 + 2) .or.    k3.le.0 .or. k4.le.0 .or. k4.lt.(k3 + 2)) then
                write (nttyo,1940) uheadx(1:j2)
1940 format(/' * Error - (XCON6/rd6d8) The following iodb',/7x,'option switch choice line read from the input file',/7x,"isn't in the required format:",/7x,'"',a,'"')

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

            write (nttyo,1950) uheadx(1:j2)
1950 format(/' * Error - (XCON6/rd6d8) The iodb option switch',/7x,'choice line',/7x,'"',a,'"',/7x,'read from the input file references an out-of-range',' option',/7x,'choice index.')

            go to 990

900 continue

            ! Check the full option choice string.
            ustr = uodbox(j,n)
            j3 = ilnobl(ustr)

            if (index(uheadx(1:j2),ustr(1:j3)) .le. 0) then
                j4 = ilnobl(uodbtx(n))
                write (nttyo,1962) uheadx(1:j2),ustr(1:j3),uodbtx(n)(1:j4)
1962 format(/' * Error - (XCON6/rd6d8) The iodb option switch',' choice line',/7x,'"',a,'"',/7x," read from the input file doesn't contain the",' matching defined string',/7x,'"',a,'".',/7x,'This line belongs to the option switch whose title',' string is',/7x,'"',a,'".')

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
                    write (nttyo,1970) ustr24(1:j3),uheadx(1:j2),uodbtx(n)(1:j4)
1970 format(/" * Error - (XCON6/rd6d8) Don't recognize the",' string "',a,'"',/7x,'that appears on the iodb',' option switch choice line',/7x,'"',a,'"',/7x,'read from the input file. An option choice should',' be chosen by',/7x,'placing a "*", "x", or "X" in the',' checkbox ("[ ]"). This',/7x,'choice line belongs to',' the option switch whose title string is',/7x,'"',a,'".')

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
            write (nttyo,1980) uodbtx(n)(1:j2)
1980 format(/' * Warning - (XCON6/rd6d8) No option choice was',' checked on the input file',/7x,'for the iodb option',' switch whose title string is',/7x,'"',a,'".')

            do j = 1,jdbxpa
                ival = iodbox(j,n)

                if (ival .eq. 0) then
                    go to 930
                end if
            end do

            write (nttyo,1990)
1990 format(/7x,'A default value of 0 has been applied, but no',/7x,'matching option choice string is defined.')

            go to 940

930 continue
            j3 = ilnobl(uodbox(j,n))
            write (nttyo,1992) uodbox(j,n)(1:j3)
1992 format(/7x,'A default value of 0 has been applied. The',' matching string is',/7x,'"',a,'".')

940 continue
        end if

        if (nmark .gt. 1) then
            j2 = ilnobl(uodbtx(n))
            j = jlast
            j3 = ilnobl(uodbox(j,n))
            write (nttyo,1994) uodbtx(n)(1:j2),uodbox(j,n)(1:j3)
1994 format(/' * Warning - (XCON6/rd6d8) More than one option',' choice was checked',/7x,'on the input file for the iodb',' option switch whose title string is',/7x,'"',a,'".',/7x,'The last choice checked will be used. The',' matching string is',/7x,'"',a,'".')
        end if
    end do

950 continue

    ! Mineral sub-set selection suppression options.
    ! Read the block title from a two-line header.
    uheadx = 'Mineral Sub-Set Selection Suppression Options'
    nfldtx = 2
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    ! Read the option title from a two-line header.
    uheadx = 'Option'
    nfldtx = 0
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    ustr24 = ufield(1)(1:24)
    call locase(ustr24)
    call locase(uheadx)
    j2 = ilnobl(ustr24)
    j3 = ilnobl(uheadx)

    if (ustr24(1:j2) .ne. uheadx(1:j3)) then
        write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
        qrderr = .true.
        go to 999
    end if

    nxi = 0

    ! Read the first line.
    nfldtx = 0
    call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

    if (qrderr) then
        go to 999
    end if

    ! The label below marks a return point for processing subsequent
    ! lines in the current data block.
330 continue
    ustr = ufield(1)

    ! There are no data remaining in the current block if a
    ! separator line has been encountered.
    if (ustr(1:8) .eq. '--------') then
        go to 380
    end if

    ustrn = ustr(1:24)
    call locase(ustrn)

    if (ustrn(1:5).eq.'none ' .or. ustr(1:1).eq.' ') then
        ustr = 'None'
    end if

    qnone = ustr(1:5).eq.'None '

    if (.not.qnone)  then
        nxi = nxi + 1

        if (nxi .gt. nxopmx) then
            write (nttyo,1320) nxopmx
1320 format(/' * Error - (XCON6/rd6d8) Have too many mineral',/7x,'subset-selection suppression options. The code is',/7x,'only dimensioned for ',i3,' such options. Reduce the',/7x,'number of options or increase the dimensioning',/7x,'parameter nxoppa.')

            go to 990
        end if

        ustrn = ustr(1:24)
        call locase(ustrn)

        do n = 1,4
            uheadx = uxopti(n)
            call locase(uheadx)

            if (ustrn(1:16) .eq. uheadx(1:16)) then
                uxopt(nxi) = uxopti(n)(1:8)
                go to 350
            end if
        end do

        j2 = ilnobl(ustr)
        write (nttyo,3060) ustr(1:j2)
3060 format(/" * Error - (XCON6/rd6d8) Don't recognize the",' uxopt option string',/7x,'"',a,'". This should',' be one of the strings',/7x,'defined in the uxopti array.',' The valid strings are:',/)

        do n = 1,4
            j3 = ilnobl(uxopti(n))
            write (nttyo,3062) uxopti(n)(1:j3)
3062 format(9x,a)
        end do

        go to 990

350 continue
        uxcat(nxi) = ufield(2)(1:24)
    end if

    ! Read the next line. Go back to process it.
    nfldtx = 0
    call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

    if (qrderr) then
        go to 999
    end if

    go to 330

380 continue
    nxopt = nxi

    ! Exceptions to the mineral sub-set selection suppression options.
    ! Read the block title from a two-line header.
    uheadx = 'Exceptions to the Mineral Sub-Set Selection Suppression'
    nfldtx = 2
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    ! Read the mineral title from a two-line header.
    uheadx = 'Mineral'
    nfldtx = 0
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    ustr24 = ufield(1)(1:24)
    call locase(ustr24)
    call locase(uheadx)
    j2 = ilnobl(ustr24)
    j3 = ilnobl(uheadx)

    if (ustr24(1:j2) .ne. uheadx(1:j3)) then
        write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
        qrderr = .true.
        go to 999
    end if

    nxic = 0

    ! Read the first line.
    nfldtx = 0
    call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

    if (qrderr) then
        go to 999
    end if

    ! The label below marks a return point for processing subsequent
    ! lines in the current data block.
410 continue
    ustr = ufield(1)

    ! There are no data remaining in the current block if a
    ! separator line has been encountered.
    if (ustr(1:8) .eq. '--------') then
        go to 420
    end if

    ustrn = ustr(1:24)
    call locase(ustrn)

    if (ustrn(1:5).eq.'none ' .or. ustr(1:1).eq.' ') then
        ustr = 'None'
    end if

    qnone = ustr(1:5).eq.'None '

    if (.not.qnone)  then
        nxic = nxic + 1

        if (nxic .gt. nxpemx) then
            write (nttyo,1370) nxpemx
1370 format(/' * Error - (XCON6/rd6d8) Have too many',/7x,'exceptions specified to the mineral subset-selection',/7x,'suppression options. The code is only dimensioned',/7x,'for ',i3,'exceptions. Reduce the number of exceptions',/7x,'or increase the dimensioning parameter nxpepa.')

            go to 990
        end if

        uxopex(nxic) = ufield(1)(1:24)
    end if

    ! Read the next line. Go back to process it.
    nfldtx = 0
    call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

    if (qrderr) then
        go to 999
    end if

    go to 410

420 continue
    nxopex = nxic

    ! Fixed fugacity options.
    ! Read the block title from a two-line header.
    uheadx = 'Fixed Fugacity Options'
    nfldtx = 2
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    ! Read the second part of the block title from two lines.
    uheadx = 'Gas'
    nfldtx = 0
    call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

    if (qrderr) then
        go to 999
    end if

    ustr24 = ufield(1)(1:24)
    call locase(ustr24)
    call locase(uheadx)
    j2 = ilnobl(ustr24)
    j3 = ilnobl(uheadx)

    if (ustr24(1:j2) .ne. uheadx(1:j3)) then
        write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
        qrderr = .true.
        go to 999
    end if

    uheadx = '(uffg(n))'
    nfldtx = 0
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    ustr24 = ufield(1)(1:24)
    call locase(ustr24)
    call locase(uheadx)
    j2 = ilnobl(ustr24)
    j3 = ilnobl(uheadx)

    if (ustr24(1:j2) .ne. uheadx(1:j3)) then
        write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
        qrderr = .true.
        go to 999
    end if

    nfi = 0

    ! Read the first line.
    nfldtx = 0
    call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

    if (qrderr) then
        go to 999
    end if

    ! The label below marks a return point for processing subsequent
    ! lines in the current data block.
430 continue
    ustr = ufield(1)

    ! There are no data remaining in the current block if a
    ! separator line has been encountered.
    if (ustr(1:8) .eq. '--------') then
        go to 440
    end if

    ustrn = ustr(1:24)
    call locase(ustrn)

    if (ustrn(1:5).eq.'none ' .or. ustr(1:1).eq.' ') then
        ustr = 'None'
    end if

    qnone = ustr(1:5).eq.'None '

    if (.not.qnone) then
        nfi = nfi + 1

        if (nfi .gt. nffgmx) then
            write (nttyo,1420) nffgmx
1420 format(/' * Error - (XCON6/rd6d8) Have too many gases whose',' fugacities are to be fixed.',/7x,'The code is only',' dimensioned for ',i4,' such gases. Reduce the number',/7x,'of gases or increase the dimensioning parameter nffgpa.')

            go to 990
        end if

        uffg(nfi) = ufield(1)(1:24)
        ustr = ufield(2)
        call chreal(nttyo,qrderr,ustr,var)

        if (qrderr) then
            go to 999
        end if

        moffg(nfi) = var
        ustr = ufield(3)
        call chreal(nttyo,qrderr,ustr,var)

        if (qrderr) then
            go to 999
        end if

        xlkffg(nfi) = var
    end if

    ! Read the next line. Go back to process it.
    nfldtx = 0
    call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

    if (qrderr) then
        go to 999
    end if

    go to 430

440 continue
    nffg = nfi

    ! Numerical parameters.
    ! Read the block title from a two-line header.
    uheadx = 'Numerical parameters'
    nfldtx = 1
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    ! Read the maximum finite-difference order (nordmx)
    ! from a one-line header.
    uheadx = 'Max. finite-difference order'
    nfldtx = 3
    call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)

    if (qrderr) then
        go to 999
    end if

    ustr = ufield(2)
    call chrint(ivar,nttyo,qrderr,ustr)

    if (qrderr) then
        go to 999
    end if

    nordmx = ivar

    ! Read the beta convergence tolerance (tolbt) from a one-line
    ! header.
    uheadx = 'Beta convergence tolerance'
    nfldtx = 3
    call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)

    if (qrderr) then
        go to 999
    end if

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

    ustr = ufield(2)
    call chreal(nttyo,qrderr,ustr,var)

    if (qrderr) then
        go to 999
    end if

    toldl = var

    ! Read the maximum number of Newton-Raphson iterations (itermx)
    ! from a one-line header.
    uheadx = 'Max. No. of N-R iterations'
    nfldtx = 3
    call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)

    if (qrderr) then
        go to 999
    end if

    ustr = ufield(2)
    call chrint(ivar,nttyo,qrderr,ustr)

    if (qrderr) then
        go to 999
    end if

    itermx = ivar

    ! Read the search/find convergence tolerance (tolxsf) from
    ! a one-line header.
    uheadx = 'Search/find convergence tolerance'
    nfldtx = 3
    call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)

    if (qrderr) then
        go to 999
    end if

    ustr = ufield(2)
    call chreal(nttyo,qrderr,ustr,var)

    if (qrderr) then
        go to 999
    end if

    tolxsf = var

    ! Read the saturation tolerance (tolsat) from a one-line header.
    uheadx = 'Saturation tolerance'
    nfldtx = 3
    call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)

    if (qrderr) then
        go to 999
    end if

    ustr = ufield(2)
    call chreal(nttyo,qrderr,ustr,var)

    if (qrderr) then
        go to 999
    end if

    tolsat = var

    ! Read the maximum number of phase assemblage tries (ntrymx) from
    ! a one-line header.
    uheadx = 'Max. No. of Phase Assemblage Tries'
    nfldtx = 3
    call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)

    if (qrderr) then
        go to 999
    end if

    ustr = ufield(2)
    call chrint(ivar,nttyo,qrderr,ustr)

    if (qrderr) then
        go to 999
    end if

    ntrymx = ivar

    ! Read the zero order step size (in Xi) (dlxmx0) from
    ! a one-line header.
    uheadx = 'Zero order step size (in Xi)'
    nfldtx = 3
    call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)

    if (qrderr) then
        go to 999
    end if

    ustr = ufield(2)
    call chreal(nttyo,qrderr,ustr,var)

    if (qrderr) then
        go to 999
    end if

    dlxmx0 = var

    ! Read the maximum interval in Xi between PRS transfers
    ! from a two-line header (the separator line is the end of the
    ! current block).
    uheadx = 'Max. interval in Xi between PRS transfers'
    nfldtx = 3
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    ustr = ufield(2)
    call chreal(nttyo,qrderr,ustr,var)

    if (qrderr) then
        go to 999
    end if

    dlxdmp = var

    ! Secondary title.
    ! Read the block title ("Secondary Title") from a two-line header.
    uheadx = 'Secondary Title'
    nfldtx = 2
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    ! Now read the secondary title itself.
    n = 0

    do nn = 1,ntitmx + 1
        read (ninpts,1000,err=990) uline1
        call parsln(nfldmx,nfldt,nlchmx,ufield,uline1,ulscr)
        ustr = ufield(1)

        ! A separator line terminates the this block. It is not part
        ! of the title itself.
        if (ustr(1:8) .eq. '--------') then
            go to 520
        end if

        n = n + 1

        if (n .gt. ntitmx) then
            write (nttyo,5015) ntitmx
5015 format(/' * Error - (XCON6/rd6d8) Have too many lines in',/7x,'the secondary title. The code is only dimensioned for ',i4,/7x,'lines. Reduce the size of the secondary title or ','increase',/7x,'the dimensioning parameter ntitpa.')

            go to 990
        end if

        utitl2(n) = ufield(1)
    end do

520 continue
    ntitl2 = n

    ! Special basis switches.
    ! Read the block title from a two-line header.
    uheadx = 'Special Basis Switches'
    nfldtx = 2
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

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
            go to 525
        end if

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

        if (.not.qnone) then
            usbsw(2,nsbsw) = ufield(2)(1:48)
        end if
    end do

525 continue
    nsbswt = nsbsw

    ! Original temperature.
    tempci = 0.

    ! Read the data from a two-line header.
    uheadx = 'Original temperature (C)'
    nfldtx = 3
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    ustr = ufield(2)
    call chreal(nttyo,qrderr,ustr,var)

    if (qrderr) then
        go to 999
    end if

    tempci = var

    ! Original pressure.
    pressi = 0.

    ! Read the data from a two-line header.
    uheadx = 'Original pressure (bars)'
    nfldtx = 3
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    ustr = ufield(2)
    call chreal(nttyo,qrderr,ustr,var)

    if (qrderr) then
        go to 999
    end if

    pressi = var

    ! Ion exchanger creation.
    net = 0

    ! Read a two-line header for the block.
    uheadx = 'Create Ion Exchangers'
    nfldtx = 2
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

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

    ustr = ufield(1)
    call lejust(ustr)
    uheadx = 'Advisory:'
    call locase(ustr)
    call locase(uheadx)

    if (ustr(1:9) .ne. uheadx(1:9)) then
        qrderr = .true.
        j3 = ilnobl(uline1(1:70))
        write (nttyo,1210) uheadx(1:9),uline1(1:j3)
1210 format(/' * Error - (XCON6/rd6d8) Was expecting to find the',/7x,'string "',a,'" at the start of the line beginning with',/7x,'"',a,'".')

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
1220 format(/' * Warning - (XCON6/rd6d8) Could not interpret the',' advisory on whether',/7x,'or not one or more exchanger',' blocks follow. The problem is with',/7x,'the line',' beginning with',/7x,'"',a,'".',/7x,'The presence of such',' blocks will be determined directly.')
        end if
    end if

    ! Read the first line of the qgexsh option (on processing, show
    ! at least one exchanger creation block on a "D" format input or
    ! pickup file).
    uheadx = 'Option: on further processing (writing a pickup file or running XCON6 on the'
    nfldtx = 1
    call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)

    if (qrderr) then
        go to 999
    end if

    ! Read the second line of the qgexsh option.
    uheadx = 'present file), force the inclusion of at least one such block (qgexsh):'
    nfldtx = 1
    call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)

    if (qrderr) then
        go to 999
    end if

    ! Read the third and last line of the qgexsh option, plus the
    ! following separator line.
    nfldtx = 1
    call rdd2l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    k1 = index(uline1,'[')
    k2 = index(uline1,']')
    j3 = ilnobl(uline1)

    if (k1.le.0 .or. k2.le.0 .or. k2.lt.(k1 + 2)) then
        write (nttyo,1230) uline1(1:j3)
1230 format(/' * Error - (XCON6/rd6d8) Was expecting to find an',/7x,'option check box "[ ]" in the line beginning with',/7x,'"',a,'"',/7x,'The expected line contains the box for',' the qgexsh flag.')

        go to 990
    end if

    ! Check the full option choice string.
    ustr = '(.true.)'

    if (index(uline1(1:j3),ustr(1:8)) .le. 0) then
        write (nttyo,1240) uline1(1:j2),ustr(1:8)
1240 format(/' * Error - (XCON6/rd6d8) The qgexsh flag line',/7x,'"',a,'"',/7x,"read from the input file doesn't contain",' the expected string',/7x,'"',a,'".')

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
1250 format(/" * Error - (XCON6/rd6d8) Don't recognize the",' string "',a,'"',/7x,'that appears on the qgexsh',' option switch choice line',/7x,'"',a,'"',/7x,'read from the input file. An option choice should',' be chosen by',/7x,'placing a "*", "x", or "X" in the',' checkbox ("[ ]").')

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
1260 format(/' * Note - (XCON6/rd6d8) The advisory on the',' presence of',/7x,'ion exchanger blocks was negative,',' but one or more blocks is',/7x,'actually present.')
        end if
    else
        ! An exchanger block is not present.
        qgexrd = .false.

        if (qgexbs .and. qgexbf) then
            write (nttyo,1270)
1270 format(/' * Note - (XCON6/rd6d8) The advisory on the',' presence of',/7x,'ion exchanger blocks was positive,',' but no blocks are actually',/7x,'present.')
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
            go to 570
        end if

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
                write (nttyo,1620) netmax,ne
1620 format(/' * Error - (XCON6/rd6d8) Have exceeded the',' maximum number of ',i3,/7x,'generic ion exchange phases',' while reading the data to create',/7x,'such phases.',' Increase the dimensioning parameter netpar',/7x,'to at',' least ',i3,'.')

                go to 990
            end if

            ugexp(ne) = ustr(1:24)
            jgext(ne) = 0
        end if

        ! Read the separator line following the line containing
        ! the name of an exchanger phase.
        nfldtx = 1
        call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

        if (qrderr) then
            go to 999
        end if

        ustr24 = ufield(1)(1:24)
        uheadx = '--------'

        if (ustr24(1:8) .ne. uheadx(1:8)) then
            j2 = ilnobl(uline1)
            j2 = min(j2,70)
            write (nttyo,1070) uline1(1:j2)
1070 format(/' * Error - (XCON6/rd6d8) Read the following line',' where a separator line',/7x,'("|------- ... -------|")',' was expected:',/7x,'"',a,'".')

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

        ustr24 = ufield(2)(1:24)
        uheadx = 'Mol. Wt. (Z)'
        call locase(ustr24)
        call locase(uheadx)
        j2 = ilnobl(ustr24)
        j3 = ilnobl(uheadx)

        if (ustr24(1:j2) .ne. uheadx(1:j3)) then
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

        ustr24 = ufield(2)(1:24)
        uheadx = 'Exchange model'
        call locase(ustr24)
        call locase(uheadx)
        j2 = ilnobl(ustr24)
        j3 = ilnobl(uheadx)

        if (ustr24(1:j2) .ne. uheadx(1:j3)) then
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
                ustr24 = ufield(2)(1:24)
                uheadx = 'Exchange site'
                call locase(ustr24)
                call locase(uheadx)
                j2 = ilnobl(ustr24)
                j3 = ilnobl(uheadx)

                if (ustr24(1:j2) .ne. uheadx(1:j3)) then
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
                        write (nttyo,1290) jetmax,ugexp(ne)(1:j2),je
1290 format(/' * Error - (XCON6/rd6d8) Have exceeded the',' maximum number of ',i3,/7x,'exchange sites on a',' generic ion exchange phase while reading',/7x,'the',' data to create ',a,'. Increase the',/7x,'dimensioning parameter jetpar to at least ',i3,'.')

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
                    go to 560
                end if

                ! Have found the end of the current block.
                go to 570
            end if

            ! Read the separator line following the line containing
            ! the name of an exchange site.
            nfldtx = 1
            call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

            if (qrderr) then
                go to 999
            end if

            ustr24 = ufield(1)(1:24)
            uheadx = '--------'

            if (ustr24(1:8) .ne. uheadx(1:8)) then
                j2 = ilnobl(uline1)
                j2 = min(j2,70)
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

            ustr24 = ufield(2)(1:24)
            uheadx = 'Stoich. number'
            call locase(ustr24)
            call locase(uheadx)
            j2 = ilnobl(ustr24)
            j3 = ilnobl(uheadx)

            if (ustr24(1:j2) .ne. uheadx(1:j3)) then
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

            ustr24 = ufield(2)(1:24)
            uheadx = 'Electr. charge'
            call locase(ustr24)
            call locase(uheadx)
            j2 = ilnobl(ustr24)
            j3 = ilnobl(uheadx)

            if (ustr24(1:j2) .ne. uheadx(1:j3)) then
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
                    ustr24 = ufield(2)(1:24)
                    uheadx = 'Reaction'
                    call locase(ustr24)
                    call locase(uheadx)
                    j2 = ilnobl(ustr24)
                    j3 = ilnobl(uheadx)

                    if (ustr24(1:j2) .ne. uheadx(1:j3)) then
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
                            j3 = ilnobl(ugexj(je,ne))
                            write (nttyo,1820) netmax,ugexj(je,ne)(1:j3),ugexp(ne)(1:j2),n
1820 format(/' * Error - (XCON6/rd6d8) Have exceeded the',' maximum number of ',i3,/7x,'reactions for a site',' belonging to a generic ion exchange',/7x,'phase',' while reading the data for site ',a,' of exchange',' phase',/7x,a,'. Increase the dimensioning',' parameter',/7x,'ietpar to at least ',i3,'.')

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
                            go to 550
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
                        go to 560
                    end if

                    ! Have found the end of the current block.
                    go to 570
                end if

                ! Read the separator line following the line containing
                ! the string containing an exchange reaction in condensed
                ! format.
                nfldtx = 1
                call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

                if (qrderr) then
                    go to 999
                end if

                ustr24 = ufield(1)(1:24)
                uheadx = '--------'

                if (ustr24(1:8) .ne. uheadx(1:8)) then
                    j2 = ilnobl(uline1)
                    j2 = min(j2,70)
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

                ustr24 = ufield(2)(1:24)
                uheadx = 'Parameter'
                call locase(ustr24)
                call locase(uheadx)
                j2 = ilnobl(ustr24)
                j3 = ilnobl(uheadx)

                if (ustr24(1:j2) .ne. uheadx(1:j3)) then
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

                ustr24 = ufield(2)(1:24)
                uheadx = 'K func.'
                call locase(ustr24)
                call locase(uheadx)
                j2 = ilnobl(ustr24)
                j3 = ilnobl(uheadx)

                if (ustr24(1:j2) .ne. uheadx(1:j3)) then
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

                ustr24 = ufield(2)(1:24)
                uheadx = 'DelH0r'
                call locase(ustr24)
                call locase(uheadx)
                j2 = ilnobl(ustr24)
                j3 = ilnobl(uheadx)

                if (ustr24(1:j2) .ne. uheadx(1:j3)) then
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

                ustr24 = ufield(2)(1:24)
                uheadx = 'DelV0r'
                call locase(ustr24)
                call locase(uheadx)
                j2 = ilnobl(ustr24)
                j3 = ilnobl(uheadx)

                if (ustr24(1:j2) .ne. uheadx(1:j3)) then
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

550 continue
        end do

560 continue
    end do

570 continue

575 continue

    ! Nxmod options.
    uheadx = 'Alter/Suppress options'
    nfldtx = 2
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    ! Read the first part of a table header from a one-line header.
    uheadx = 'Species'
    nfldtx = 3
    call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)

    if (qrderr) then
        go to 999
    end if

    ! Read the second part of the table header from a two-line header.
    uheadx = '(uxmod(n))'
    nfldtx = 3
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

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

        ustr = ufield(1)

        if (ustr(1:8) .eq. '--------') then
            go to 590
        end if

        call locase(ustr)
        ustrn = ustr(1:24)
        call locase(ustrn)

        if (ustrn(1:5).eq.'none ' .or. ustr(1:1).eq.' ') then
            ustr = 'None'
        end if

        if (ustr(1:5) .eq. 'None ') then
            go to 580
        end if

        n = n + 1

        if (n .gt. nxmdmx) then
            write (nttyo,1500) nxmdmx
1500 format(/' * Error - (XCON6/rd6d8) Have too many nxmod',/7x,'alter/suppress options. The code is only dimensioned',/7x,'for ',i3,' such options. Reduce the number of such',/7x,'options or increase the dimensioning parameter nxmdpa.')

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
1510 format(/" * Error - (XCON6/rd6d8) Don't recognize the",' alter/suppress option',/7x,'string "',a,'". This should',' be one of the strings',/7x,'defined in the ukxm array.',' The valid strings are:',/)

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

580 continue
    end do

590 continue
    nxmod = n

    ! Iopg Print Option Switches.
    ! Note: iopg(1) = iopg1, etc.
    uheadx = 'Iopg Activity Coefficient Option Switches ("( 0)" marksdefault choices)'
    nfldtx = 1
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

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

        k1 = index(uheadx,'(')
        k2 = index(uheadx,')')
        k3 = index(uheadx,'- ')

        if (k1.le.0 .or. k2.le.0 .or. k2.lt.(k1 + 2) .or. k3.le.0 .or.  k3.le.k2) then
            write (nttyo,1705) uheadx(1:j2)
1705 format(/' * Error - (XCON6/rd6d8) The iopg option switch',' title line',/7x,'"',a,'"',/7x,'read from the'," input file isn't in the required format.")

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
            write (nttyo,1712) uheadx(1:j2),n
1712 format(/' * Error - (XCON6/rd6d8) The iopg option switch',' title line',/7x,'"',a,'"',/7x,'read from the',' input file references an option switch index',/7x,'of ',i3,', which is out of range.')

            go to 990
        end if

        ! Check the full title string.
        ustr = uopgtx(n)
        j3 = ilnobl(ustr)

        if (index(uheadx(1:j2),ustr(1:j3)) .le. 0) then
            write (nttyo,1720) uheadx(1:j2),ustr(1:j3)
1720 format(/' * Error - (XCON6/rd6d8) The iopg option switch',' title string',/7x,'"',a,'"',/7x,"read from the input file doesn't contain the",' matching defined string',/7x,'"',a,'".')

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

            uheadx = ufield(1)
            j2 = ilnobl(uheadx)

            if (uheadx(1:8) .eq. '--------') then
                go to 720
            end if

            if (jj .gt. jpgxpa) then
                j3 = ilnobl(uopgtx(n))
                write (nttyo,1730) uopgtx(n)(1:j3),jpgxpa
1730 format(/' * Error - (XCON6/rd6d8) Have too many option',/7x,'choice lines for the iopg option switch whose title',' string is',/7x,'"',a,'".',/7x,' The code is only',' dimensioned for ',i3,' such lines. Reduce the',/7x,'number of such lines or increase the dimensioning',' parameter jpgxpa.')

                go to 990
            end if

            k1 = index(uheadx,'[')
            k2 = index(uheadx,']')
            k3 = index(uheadx,'(')
            k4 = index(uheadx,')')

            if (k1.le.0 .or. k2.le.0 .or. k2.lt.(k1 + 2) .or.    k3.le.0 .or. k4.le.0 .or. k4.lt.(k3 + 2)) then
                write (nttyo,1740) uheadx(1:j2)
1740 format(/' * Error - (XCON6/rd6d8) The following iopg',/7x,'option switch choice line read from the input file',/7x,"isn't in the required format:",/7x,'"',a,'"')

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

            write (nttyo,1750) uheadx(1:j2)
1750 format(/' * Error - (XCON6/rd6d8) The iopg option switch',/7x,'choice line',/7x,'"',a,'"',/7x,'read from the input file references an out-of-range',' option',/7x,'choice index.')

            go to 990

700 continue

            ! Check the full option choice string.
            ustr = uopgox(j,n)
            j3 = ilnobl(ustr)

            if (index(uheadx(1:j2),ustr(1:j3)) .le. 0) then
                j4 = ilnobl(uopgtx(n))
                write (nttyo,1760) uheadx(1:j2),ustr(1:j3),uopgtx(n)(1:j4)
1760 format(/' * Error - (XCON6/rd6d8) The iopg option switch',' choice line',/7x,'"',a,'"',/7x," read from the input file doesn't contain the",' matching defined string',/7x,'"',a,'".',/7x,'This line belongs to the option switch whose title',' string is',/7x,'"',a,'".')

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
                    write (nttyo,1770) ustr24(1:j3),uheadx(1:j2),uopgtx(n)(1:j4)
1770 format(/" * Error - (XCON6/rd6d8) Don't recognize the",' string "',a,'"',/7x,'that appears on the iopg',' option switch choice line',/7x,'"',a,'"',/7x,'read from the input file. An option choice should',' be chosen by',/7x,'placing a "*", "x", or "X" in the',' checkbox ("[ ]"). This',/7x,'choice line belongs to',' the option switch whose title string is',/7x,'"',a,'".')

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
            write (nttyo,1780) uopgtx(n)(1:j2)
1780 format(/' * Warning - (XCON6/rd6d8) No option choice was',' checked on the input file',/7x,'for the iopg option',' switch whose title string is',/7x,'"',a,'".')

            do j = 1,jpgxpa
                ival = iopgox(j,n)

                if (ival .eq. 0) then
                    go to 730
                end if
            end do

            write (nttyo,1790)
1790 format(/7x,'A default value of 0 has been applied, but no',/7x,'matching option choice string is defined.')

            go to 740

730 continue
            j3 = ilnobl(uopgox(j,n))
            write (nttyo,1792) uopgox(j,n)(1:j3)
1792 format(/7x,'A default value of 0 has been applied. The',' matching string is',/7x,'"',a,'".')

740 continue
        end if

        if (nmark .gt. 1) then
            j2 = ilnobl(uopgtx(n))
            j = jlast
            j3 = ilnobl(uopgox(j,n))
            write (nttyo,1794) uopgtx(n)(1:j2),uopgox(j,n)(1:j3)
1794 format(/' * Warning - (XCON6/rd6d8) More than one option',' choice was checked',/7x,'on the input file for the iopg',' option switch whose title string is',/7x,'"',a,'".',/7x,'The last choice checked will be used. The',' matching string is',/7x,'"',a,'".')
        end if
    end do

750 continue

    ! Matrix Index Limits.
    ! Read the block title from a two-line header.
    uheadx = 'Matrix Index Limits'
    nfldtx = 1
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    ! Read the number of chemical elements (kct) from a one-line header.
    uheadx = 'No. of chem. elements'
    nfldtx = 3
    call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)

    if (qrderr) then
        go to 999
    end if

    ustr = ufield(2)
    call chrint(ivar,nttyo,qrderr,ustr)

    if (qrderr) then
        go to 999
    end if

    kct = ivar

    if (kct .gt. nctmax) then
        write (nttyo,3240) nctmax
3240 format(/' * Error - (XCON6/rd6d8) Have too many chemical',/7x,'elements present. The code is only dimensioned',/7x,'for ',i3,' elements. Reduce the number of elements',/7x,'or increase the dimensioning parameter nctpar.')

        go to 990
    end if

    ! Read the number of basis species (kbt) from a one-line header.
    uheadx = 'No. of basis species'
    nfldtx = 3
    call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)

    if (qrderr) then
        go to 999
    end if

    ustr = ufield(2)
    call chrint(ivar,nttyo,qrderr,ustr)

    if (qrderr) then
        go to 999
    end if

    kbt = ivar

    if (kbt .gt. nbtmax) then
        write (nttyo,3250) nbtmax
3250 format(/' * Error - (XCON6/rd6d8) Have too many basis',/7x,'species present. The code is only dimensioned',/7x,'for ',i3,' basis species. Reduce the number of elements',/7x,'or increase the dimensioning parameter nctpar.')

        go to 990
    end if

    ! Read the matrix index of the last pure min. (kmt) from a
    ! one-line header.
    uheadx = 'Index of last pure min.'
    nfldtx = 3
    call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)

    if (qrderr) then
        go to 999
    end if

    ustr = ufield(2)
    call chrint(ivar,nttyo,qrderr,ustr)

    if (qrderr) then
        go to 999
    end if

    kmt = ivar

    ! Read the matrix index of the last solid solution component
    ! from a one-line header.
    uheadx = 'Index of last sol-sol.'
    nfldtx = 3
    call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)

    if (qrderr) then
        go to 999
    end if

    ustr = ufield(2)
    call chrint(ivar,nttyo,qrderr,ustr)

    if (qrderr) then
        go to 999
    end if

    kxt = ivar

    ! Read the matrix size (kdim) from a one-line header.
    uheadx = 'Matrix size'
    nfldtx = 3
    call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)

    if (qrderr) then
        go to 999
    end if

    ustr = ufield(2)
    call chrint(ivar,nttyo,qrderr,ustr)

    if (qrderr) then
        go to 999
    end if

    kdim = ivar

    if (kdim .gt. kmax) then
        write (nttyo,3260) kmax
3260 format(/' * Error - (XCON6/rd6d8) Have too many matrix',/7x,'variables. The code is only dimensioned for ',i3,/7x,'matrix variables. Reduce the number of such variables',/7x,'or increase the dimensioning parameter kpar.')

        go to 990
    end if

    ! Read the PRS data flag (kprs) from a two-line header
    ! (the separator line is the end of the current block).
    uheadx = 'PRS data flag'
    nfldtx = 3
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    ustr = ufield(2)
    call chrint(ivar,nttyo,qrderr,ustr)

    if (qrderr) then
        go to 999
    end if

    kprs = ivar

    ! Species for which mass balances are defined.
    ! Read the first part of the block title from a one-line header.
    uheadx = 'Mass Balance Species (Matrix Row Variables)'
    nfldtx = 3
    call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)

    if (qrderr) then
        go to 999
    end if

    ! Read the second part of the block title from a two-line header.
    uheadx = '(ubmtbi(n))'
    nfldtx = 3
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    nbi = 0

    ! Read the first line.
    nfldtx = 0
    call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

    if (qrderr) then
        go to 999
    end if

    ! The label below marks a return point for processing subsequent
    ! lines in the current data block.
910 continue
    ustr = ufield(1)

    ! There are no data remaining in the current block if a
    ! separator line has been encountered.
    if (ustr(1:8) .eq. '--------') then
        go to 915
    end if

    ! Check for no species for which mass balances are defined.
    ustrn = ustr(1:24)
    call locase(ustrn)

    if (ustrn(1:5).eq.'none ' .or. ustr(1:1).eq.' ') then
        ustr = 'None'
    end if

    if (ustr(1:5) .eq. 'None ') then
        nbi = 0
        nfldtx = 0
        call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

        if (qrderr) then
            go to 999
        end if

        ustr = ufield(1)

        if (ustr(1:8) .eq. '--------') then
            go to 915
        end if
    end if

    nfldtx = 3

    if (nfldt .ne. nfldtx) then
        j2 = ilnobl(uline1)
        j2 = min(j2,70)
        write (nttyo,5010) nfldt,nfldtx,uline1(1:j2)
5010 format(/' * Error - (XCON6/rd6d8) Found ',i2,' fields',' where ',i2,/7x,'were expected on the line which begins with',/7x,'"',a,'".')

        go to 990
    end if

    ux48 = ufield(1)(1:48)

    nbi = nbi + 1

    if (nbi .gt. nbtmax) then
        j2 = ilnobl(ux48)
        write (nttyo,5020) nbtmax,ux48(1:j2)
5020 format(/' * Error - (XCON6/rd6d8) The number of mass balance',' species read',/7x,'from the input file exceeded the',' maximum of ',i3,' while',/7x,'trying to read data for',' the species ',a,'.',/7x,'Increase the dimensioning',' parameter nbtpar.')

        go to 990
    end if

    ubmtbi(nbi) = ux48

    ustr = ufield(2)
    ustrn = ustr(1:24)
    call locase(ustrn)

    if (ustrn(1:16) .eq. 'moles           ') then
        jflgi(nbi) = 0
    else if (ustrn(1:16) .eq. 'make non-basis  ') then
        jflgi(nbi) = 30
    else
        j2 = ilnobl(ustr)
        write (nttyo,5030) ustr(1:j2)
5030 format(/" * Error - (XCON6/rd6d8) Can't identify the",/7x,'following Mass-Balance-Species Units/constraints string- ',a,'.',/7x,'This must be one of "Moles", or "Make non-basis".')

        go to 990
    end if

    ! Read the next line. Go back to process it.
    nfldtx = 0
    call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

    if (qrderr) then
        go to 999
    end if

    go to 910

915 continue

    nbti = nbi

    ! Mass Balance Totals (moles).
    ! Read the first part of the block title from a two-line header.
    uheadx = 'Mass Balance Totals (moles)'
    nfldtx = 1
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    ! Read the next part of the block title from a one-line header.
    uheadx = 'Basis species (info. only)'
    nfldtx = 3
    call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

    if (qrderr) then
        go to 999
    end if

    ! Read the second part of the block title from a two-line header.
    uheadx = '(ubmtbi(n))'
    nfldtx = 3
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    ! Read the lines.
    do nbi = 1,nbti
        nfldtx = 0
        call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

        if (qrderr) then
            go to 999
        end if

        ustr = ufield(1)
        ux48 = ubmtbi(nbi)
        j2 = ilnobl(ustr)
        j3 = ilnobl(ux48)

        if (ustr(1:32) .ne. ux48(1:32)) then
            write (nttyo,5040) ux48(1:j3),ustr(1:j2)
5040 format(/' * Error - (XCON6/rd6d8) Found that the mass',' balance species name',/7x,'"',a,'".',/7x,"doesn't match",' the mass balance totals species name tag',/7x,'"',a,'"',/7x,'(which has a maximum length of only 32 characters).')

            qrderr = .true.
            go to 999
        end if

        ustr = ufield(2)
        call chreal(nttyo,qrderr,ustr,var)

        if (qrderr) then
            go to 999
        end if

        mtbi(nbi) = var
        ustr = ufield(3)
        call chreal(nttyo,qrderr,ustr,var)

        if (qrderr) then
            go to 999
        end if

        mtbaqi(nbi) = var
    end do

    ! Read the electrical imbalance from a two-line header.
    uheadx = 'Electrical imbalance'
    nfldtx = 3
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    ustr24 = ufield(1)(1:24)
    call locase(ustr24)
    call locase(uheadx)
    j2 = ilnobl(ustr24)
    j3 = ilnobl(uheadx)

    if (ustr24(1:j2) .ne. uheadx(1:j3)) then
        write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
        qrderr = .true.
        go to 999
    end if

    ustr = ufield(2)
    call chreal(nttyo,qrderr,ustr,var)

    if (qrderr) then
        go to 999
    end if

    electr = var

    ! Ordinary basis switches.
    ! Read the block title from a two-line header.
    uheadx = 'Ordinary Basis Switches'
    nfldtx = 2
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

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
            go to 925
        end if

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

        if (.not.qnone) then
            uobsw(2,nobsw) = ufield(2)(1:48)
        end if
    end do

925 continue
    nobswt = nobsw

    ! Matrix column variables and values.
    ! Read the first part of the block title from a two-line header.
    uheadx = 'Matrix Column Variables and Values'
    nfldtx = 1
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    ! Read the second part of the block title from a two-line header.
    uheadx = 'Basis species (uzveci(n))'
    nfldtx = 3
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    krow = 0

    ! Read the first line.
    nfldtx = 0
    call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

    if (qrderr) then
        go to 999
    end if

    ! The label below marks a return point for processing subsequent
    ! lines in the current data block.
960 continue
    ustr = ufield(1)

    ! There are no data remaining in the current block if a
    ! separator line has been encountered.
    if (ustr(1:8) .eq. '--------') then
        go to 965
    end if

    nfldtx = 3

    if (nfldt .ne. nfldtx) then
        j2 = ilnobl(uline1)
        j2 = min(j2,70)
        write (nttyo,5010) nfldt,nfldtx,uline1(1:j2)
        go to 990
    end if

    ux48 = ufield(1)(1:48)

    krow = krow + 1

    if (krow .gt. kmax) then
        j2 = ilnobl(ux48)
        write (nttyo,5050) kmax,ux48(1:j2)
5050 format(/' * Error - (XCON6/rd6d8) The number of matrix',' column variables',/7x,'and values read from the input file',' exceeded the maximum',/7x,'of ',i3,' while trying to read',' data for the basis species',/7x,a,'. Increase the',' dimensioning parameter kmax.')

        go to 990
    end if

    uzveci(krow) = ux48
    ustr = ufield(2)
    call chreal(nttyo,qrderr,ustr,var)

    if (qrderr) then
        go to 999
    end if

    zvclgi(krow) = var

    ! Read the next line. Go back to process it.
    nfldtx = 0
    call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

    if (qrderr) then
        go to 999
    end if

    go to 960

965 continue

    kdim = krow

    ! Phases and species in the PRS.
    ! Read a two-line header for the block.
    uheadx = 'Phases and Species in the PRS'
    nfldtx = 1
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    ! Loop on phases.
    nsi = 0
    npi = 0

    do nn = 1,nprpmx + 1
        ! Read a line. If the block has not been completely read,
        ! this contains the name of a phase, and a sub-block
        ! for that phase follows. Otherwise, this line is the first line
        ! of the next block.
        nfldtx = 0
        call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

        if (qrderr) then
            go to 999
        end if

        ustr = ufield(1)
        uheadx = 'Phase'
        call locase(ustr)
        call locase(uheadx)
        j2 = ilnobl(ustr)
        j3 = ilnobl(uheadx)

        if (ustr(1:j2) .ne. uheadx(1:j3)) then
            ! Back up.
            backspace ninpts
            go to 340
        end if

        ustr = ufield(2)
        ustrn = ustr(1:24)
        call locase(ustrn)

        if (ustrn(1:5).eq.'none ' .or. ustr(1:1).eq.' ') then
            ustr = 'None'
        end if

        qnonep = ustr(1:5).eq.'None '

        if (.not.qnonep) then
            npi = npi + 1

            if (npi .gt. nprpmx) then
                write (nttyo,3640) nprpmx
3640 format(/' * Error - (XCON6/rd6d8) Have too many phases',/7x,'in the physically removed system (PRS) on the',/7x,'input file. The code is only dimensioned for ',i3,/7x,'such phases. Reduce the number of such phases',/7x,'or increase the dimensioning parameter nprppa.')

                go to 990
            end if

            uprphi(npi )= ustr(1:24)
        end if

        ! Read the separator line following the line containing
        ! the phase name.
        nfldtx = 1
        call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

        if (qrderr) then
            go to 999
        end if

        ustr24 = ufield(1)(1:24)
        uheadx = '--------'

        if (ustr24(1:8) .ne. uheadx(1:8)) then
            j2 = ilnobl(uline1)
            j2 = min(j2,70)
            write (nttyo,1070) uline1(1:j2)
            qrderr = .true.
            go to 999
        end if

        ! Read the number of moles from a two-line header.
        uheadx = '->'
        nfldtx = 4
        call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

        if (qrderr) then
            go to 999
        end if

        ustr24 = ufield(2)(1:24)
        uheadx = 'No. of Moles'
        call locase(ustr24)
        call locase(uheadx)
        j2 = ilnobl(ustr24)
        j3 = ilnobl(uheadx)

        if (ustr24(1:j2) .ne. uheadx(1:j3)) then
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
            mprphi(npi) = var
        end if

        ! Read the first title line.
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
            ustr24 = ufield(2)(1:24)
            uheadx = 'Species'
            call locase(ustr24)
            call locase(uheadx)
            j2 = ilnobl(ustr24)
            j3 = ilnobl(uheadx)

            if (ustr24(1:j2) .ne. uheadx(1:j3)) then
                write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
                qrderr = .true.
                go to 999
            end if
        end if

        ! Read the second title line.
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
            ustr24 = ufield(2)(1:24)
            uheadx = '(uprspi(i,n))'
            call locase(ustr24)
            call locase(uheadx)
            j2 = ilnobl(ustr24)
            j3 = ilnobl(uheadx)

            if (ustr24(1:j2) .ne. uheadx(1:j3)) then
                write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
                qrderr = .true.
                go to 999
            end if
        end if

        ! Read the separator line following the second title line.
        nfldtx = 1
        call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

        if (qrderr) then
            go to 999
        end if

        ustr24 = ufield(1)(1:24)
        uheadx = '--------'

        if (ustr24(1:8) .ne. uheadx(1:8)) then
            j2 = ilnobl(uline1)
            j2 = min(j2,70)
            write (nttyo,1070) uline1(1:j2)
            qrderr = .true.
            go to 999
        end if

        ! Loop on component species.
        iki = 0

        do jj = 1,iktmax + 1
            ! Read a line.
            nfldtx = 0
            call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)

            if (qrderr) then
                go to 999
            end if

            ustr = ufield(1)
            uheadx = '--------'

            if (ustr(1:8) .ne. uheadx(1:8)) then
                ustr24 = ufield(1)(1:24)
                uheadx = '--->'
                call locase(ustr24)
                call locase(uheadx)
                j2 = ilnobl(ustr24)
                j3 = ilnobl(uheadx)

                if (ustr24(1:j2) .ne. uheadx(1:j3)) then
                    write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
                    qrderr = .true.
                    go to 999
                end if

                ustr = ufield(2)
                ustrn = ustr(1:24)
                call locase(ustrn)

                if (ustrn(1:5).eq.'none ' .or. ustr(1:1).eq.' ') then
                    ustr = 'None'
                end if

                qnones = ustr(1:5).eq.'None '

                qokay = .not.qnones .and. .not.qnonep

                if (qokay) then
                    ! Have found another species.
                    iki = iki + 1

                    if (iki .gt. iktmax) then
                        j2 = ilnobl(uprphi(npi))
                        write (nttyo,3680) uprphi(npi)(1:j2),iktmax
3680 format(/' * Error - (XCON6/rd6d8) Have too many',' end-members',/7x,'in the PRS solid solution ',a,'.',/7x,'The code is only dimensioned for ',i4,' end-members per',/7x,'solid solution. Reduce',' the number of end-members or',/7x,'increase the dimensioning parameter iktpar.')

                        go to 990
                    end if

                    nsi = nsi + 1

                    if (nsi .gt. nprsmx) then
                        write (nttyo,2000) nprsmx,nsi
2000 format(/' * Error - (XCON6/rd6d8) Have exceeded the',' maximum number of ',i3,/7x,'species while reading',' the data for phases and species in the PRS',/7x,' Increase the dimensioning parameter',' nprsmx to at least ',i3,'.')

                        go to 990
                    end if

                    uprspi(nsi) = ufield(2)(1:48)
                    ustr = ufield(3)
                    call chreal(nttyo,qrderr,ustr,var)

                    if (qrderr) then
                        go to 999
                    end if

                    mprspi(nsi) = var
                end if
            else
                ! Have found end of species.
                go to 335
            end if
        end do

335 continue
    end do

340 continue

    nprpti = npi
    nprsti = nsi

    ! Read a two-line header marking the end of the current problem
    ! input.
    uheadx = 'End of problem'
    nfldtx = 1
    call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)

    if (qrderr) then
        go to 999
    end if

    go to 999

990 continue
    qrderr = .true.

999 continue
end subroutine rd6d8