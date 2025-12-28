subroutine pcraq(aamatr,apr,atwt,avgrid,cdrs,cdrsi,cess,cessi,cof,dhfe,dhfs,dvfe,dvfs,eps100,gmmatr,ipch,ipchmx,ipcv,ipcvmx,ipivot,itgenf,mtotr,nacdpr,narxmx,narxt,nat,natmax,nbt,nbtmx1,nbtmx2,nch,nco,nct,nctmax,ndata1,ndat0s,ndat1f,ndbmax,ndbptg,ndbptl,nentei,nentri,nerr,nmodwr,noutpt,nsb,nslist,ntprmx,ntprt,nttyo,nwarn,qelect,q500fl,tempc,tempcs,tmpcmx,uaqsp,udbfmt,udbval,udrsi,uelem,uessi,uspec,xdbval,xhfe,xhfs,xlke,xlks,xvfe,xvfs,xvec,yvec,zaqsp,zchar)
    !! This subroutine reads data on aqueous species from the stripped
    !! DATA0 file, processes this data, and writes the results on the
    !! DATA1 and DATA1F files. The counter nerr is incremented by
    !! one for each error encountered. The counter nwarn is incremented
    !! similarly for each warning.
    !! This subroutine is called by:
    !!   EQPT/eqpt.f
    !! Principal input:
    !!   narxt  = array of numbers of coefficients in the temperature
    !!              ranges
    !!   nbt    = the number of basis species
    !!   nct    = the number of chemical elements
    !!   ndata1 = unit number of the DATA1 file
    !!   ndat0s = unit number of the stripped DATA0 file
    !!   ndat1f = unit number of the DATA1F file
    !!   ndbptg = the number of distinct points on the temperature grid
    !!   ndbptl = the maximum number of points on the temperature grid
    !!              per line
    !!   ntprt  = the number of temperature ranges on the standard
    !!              temperature grid
    !!   qelect = flag denoting the use of "e-" instead of "O2(g)"
    !!              in writing chemical reactions
    !!   tempc  = array of temperatures (on the temperature grid)
    !!   tempcs = array of scaled temperatures (on the temperature grid)
    !!   tmpcmx = the max norm of the tempeatures on the temperature grid
    !!   udbfmt = the format for reading a line of data on the
    !!              temperature grid
    !! Principal output:
    !!   atwt   = array of atomic weights
    !!   cdrs   = array of reaction coefficients
    !!   cess   = array of elemental composition coefficients
    !!   nat    = the number of aqueous species
    !!   nerr   = cumulative error counter
    !!   nwarn  = cumulative warning counter
    !!   uaqsp  = array of aqueous species names
    !!   uelem  = array of names of chemical elements
    !!   uspec  = array of names of species
    !!   xlke   = array of log K values for the "Eh" reaction (on the
    !!              temperature grid)
    !!   xlks   = array of log K values for reactions
    !!   zaqsp  = array of charge numbers of aqueous species
    !!   zchar  = array of electrical charge numbers of species
    !! Workspace:
    !!   aamatr = matrix used to calculate the polynomial coefficients
    !!   apr    = array of polynomial coefficients (for all ranges)
    !!   avgrid = array containing the data on the temperature grid
    !!   cof    = array of fitted polynomial coefficients (for a
    !!              single temperature range)
    !!   gmmatr = a copy of the amatr matrix
    !!   ipivot = the pivot vector, used in solving matrix equations
    !!   nacdpr = array containging the number of actual data points
    !!              by range on the "log K" temperature grid; excludes
    !!   udbval = string array for reading data on the temperature grid
    !!   xdbval = holding space array for data on the temperature grid
    !!   xvec   = array of scaled temperatures corresponding to the
    !!              data in the yvec array
    !!   yvec   = array of data to be fitted (for a single temperature
    !!              range)
    implicit none

    ! Calling sequence variable declarations.
    integer :: ipchmx
    integer :: ipcvmx
    integer :: narxmx
    integer :: natmax
    integer :: nbtmx1
    integer :: nbtmx2
    integer :: nctmax
    integer :: ndbmax
    integer :: ntprmx

    integer :: ndata1
    integer :: ndat0s
    integer :: ndat1f
    integer :: noutpt
    integer :: nslist
    integer :: nttyo

    integer :: ipivot(narxmx)
    integer :: nacdpr(ntprmx)
    integer :: narxt(ntprmx)
    integer :: nentei(nctmax)
    integer :: nentri(nbtmx1)

    integer :: ipch
    integer :: ipcv
    integer :: itgenf
    integer :: nat
    integer :: nbt
    integer :: nch
    integer :: nco
    integer :: nct
    integer :: ndbptg
    integer :: ndbptl
    integer :: nerr
    integer :: nmodwr
    integer :: nsb
    integer :: ntprt
    integer :: nwarn

    logical :: qelect
    logical :: q500fl

    character(len=24) :: uaqsp(natmax)
    character(len=24) :: udrsi(nbtmx1)
    character(len=24) :: uspec(nbtmx1)
    character(len=16) :: udbval(ndbmax)
    character(len=8) :: uelem(nctmax)
    character(len=8) :: uessi(nctmax)

    character(len=16) :: udbfmt

    real(kind=8) :: atwt(nctmax)
    real(kind=8) :: cdrs(nbtmx2,nbtmx1)
    real(kind=8) :: cdrsi(nbtmx1)
    real(kind=8) :: cess(nctmax,nbtmx1)
    real(kind=8) :: cessi(nctmax)
    real(kind=8) :: dhfe(narxmx,ntprmx,ipchmx)
    real(kind=8) :: dhfs(narxmx,ntprmx,ipchmx,nbtmx1)
    real(kind=8) :: dvfe(narxmx,ntprmx,ipcvmx)
    real(kind=8) :: dvfs(narxmx,ntprmx,ipcvmx,nbtmx1)
    real(kind=8) :: xhfe(narxmx,ntprmx)
    real(kind=8) :: xhfs(narxmx,ntprmx,nbtmx1)
    real(kind=8) :: xlke(narxmx,ntprmx)
    real(kind=8) :: xlks(narxmx,ntprmx,nbtmx1)
    real(kind=8) :: xvfe(narxmx,ntprmx)
    real(kind=8) :: xvfs(narxmx,ntprmx,nbtmx1)
    real(kind=8) :: zchar(nbtmx1)
    real(kind=8) :: zaqsp(natmax)

    real(kind=8) :: apr(narxmx,ntprmx)
    real(kind=8) :: avgrid(narxmx,ntprmx)
    real(kind=8) :: tempc(narxmx,ntprmx)
    real(kind=8) :: tempcs(narxmx,ntprmx)
    real(kind=8) :: tmpcmx(ntprmx)
    real(kind=8) :: aamatr(narxmx,narxmx)
    real(kind=8) :: gmmatr(narxmx,narxmx)
    real(kind=8) :: cof(narxmx)
    real(kind=8) :: xvec(narxmx)
    real(kind=8) :: yvec(narxmx)
    real(kind=8) :: xdbval(ndbmax)
    real(kind=8) :: mtotr(nctmax)

    real(kind=8) :: eps100

    ! Local variable declarations.
    integer :: i
    integer :: ii
    integer :: ipc
    integer :: j
    integer :: j2
    integer :: j3
    integer :: j4
    integer :: j5
    integer :: k
    integer :: n
    integer :: nbt1
    integer :: nbw
    integer :: nc
    integer :: ncount
    integer :: ncts
    integer :: ndrsts
    integer :: nmodx
    integer :: nn
    integer :: nnx
    integer :: ns
    integer :: nse
    integer :: nsm1
    integer :: nspnx
    integer :: nt
    integer :: ntpr
    integer :: nxm

    integer :: ilnobl

    logical :: qblkes
    logical :: qblkrs
    logical :: qend
    logical :: qerr
    logical :: qnofes
    logical :: qnofrs
    logical :: qzeres
    logical :: qzerrs
    logical :: q500nd

    character(len=24) :: uspn(2)

    character(len=80) :: ulbufa
    character(len=80) :: ulbufb
    character(len=80) :: uline
    character(len=80) :: ux80
    character(len=72) :: uterm
    character(len=72) :: utermc
    character(len=56) :: ustrgr
    character(len=24) :: uaq24
    character(len=24) :: uend24
    character(len=24) :: ublk24
    character(len=24) :: unone
    character(len=24) :: usblkf
    character(len=24) :: ux24
    character(len=16) :: ux16
    character(len=8) :: uaqu
    character(len=8) :: uaux
    character(len=8) :: ux8
    character(len=8) :: ux8a
    character(len=8) :: ux8b
    character(len=8) :: ux8c
    character(len=8) :: uendit
    character(len=8) :: usol

    real(kind=8) :: mwtss
    real(kind=8) :: zx

    data ublk24 /'                        '/
    data uaq24 /'aqueous                 '/
    data unone /'none                    '/
    data uendit /'endit.  '/
    data uaux /'auxiliar'/
    data uaqu /'aqueous '/
    data usol /'solids  '/

    uterm(1:48) = '+-----------------------------------------------'
    uterm(49:72) = '------------------------'
    utermc = uterm
    utermc(1:1) = '*'

    ! The variable nn must be set to zero here to provide a proper
    ! diagnostic if an early end-of-file or read format error occurs.
    nn = 0

    q500nd = q500fl

    ! Set the expected index of H2O.
    nbw = 1

    ! Skip 2 lines to the 'basis species' delimiter.
    read (ndat0s,1000,end=990,err=995) uline
    read (ndat0s,1000,end=990,err=995) uline
1000 format(a)

    ! Write label.
    usblkf = uaq24
    j3 = ilnobl(uaq24)
    write (ndata1) uaq24
    write (ndat1f,1010) uaq24(1:j3)
1010 format(a)

    write (ndat1f,1010) utermc(1:72)
    write (noutpt,1020) uaq24(1:j3)
    write (nttyo,1020) uaq24(1:j3)
    write (nslist,1020) uaq24(1:j3)
1020 format(//1x,a,/)

    ! Initialize some variables.
    qelect = .false.
    uend24(1:8) = uendit(1:8)
    uend24(9:24) = ublk24(9:24)
    uspn(1) = unone
    uspn(2) = ublk24
    nspnx = 1
    nnx = -1
    nn = 1
    nmodx = 0
    nbt1 = nbt + 1
    nat = 0

    ! Note: the main loop returns here.
100 continue
    ns = min(nn,nbt1)
    ncts = 0
    ndrsts = 0

    do nc = 1,nct
        cess(nc,ns) = 0.
        cessi(nc) = 0.
    end do

    do nse = 1,nbt1
        cdrs(nse,ns) = 0.
        cdrsi(nse) = 0.
    end do

    cdrs(nbt + 2,ns) = 0.

    call initcv(uessi,nct,' ')
    call initcv(udrsi,nbt1,' ')

    ! Read the first line of the block.
110 continue

    read (ndat0s,1000,end=990,err=995) uline
    ux24 = uline(1:24)
    ux8 = ux24(1:8)

    if (ux8(1:8).eq.uaux(1:8) .or. ux8(1:8).eq.uaqu(1:8)) then
        ! Skip the terminator line, read the next line.
        read (ndat0s,1000,end=990,err=995) uline
        go to 110
    else if (ux8(1:8) .eq. usol(1:8)) then
        ! Skip the terminator line.
        read (ndat0s,1000,end=990,err=995) uline

        ! Finish and exit here.
        if (uspn(1)(1:24) .eq. unone(1:24)) then
            nnx = 0
            nn = 0
            write (nttyo,1140) nnx,unone(1:4)
            write (nslist,1140) nnx,unone(1:4)
1140 format(1x,i5,2x,a)

            write (noutpt,1150) nn,unone(1:4)
1150 format(1x,i5,1x,a)
        else
            nxm = nspnx - 1

            if (nxm .ge. 1) then
                nnx = nnx + 2

                if (nxm .eq. 1) then
                    j3 = ilnobl(uspn(1))
                    write (nttyo,1140) nnx,uspn(1)(1:j3)
                    write (nslist,1140) nnx,uspn(1)(1:j3)
                else
                    j3 = ilnobl(uspn(2))
                    write (nttyo,1160) nnx,uspn(1),uspn(2)(1:j3)
                    write (nslist,1160) nnx,uspn(1),uspn(2)(1:j3)
1160 format(1x,i5,2x,a24,2x,a)
                end if
            else if (nmodx .ne. 1) then
                j3 = ilnobl(uspn(2))
                write (nttyo,1160) nnx,uspn(1),uspn(2)(1:j3)
            end if
        end if

        write (ndata1) uend24,ublk24,ublk24
        j3 = ilnobl(uend24)
        write (ndat1f,1010) uend24(1:j3)
        write (ndat1f,1010) utermc(1:72)

        if (nerr .gt. 0) then
            stop
        end if

        go to 999
    end if

    ! Have a species block.
    nat = nat + 1

    j2 = ilnobl(ux24)

    if (j2 .le. 0) then
        write (ux8a,'(i5)') nat
        call lejust(ux8a)
        j3 = ilnobl(ux8a)
        j4 = ilnobl(usblkf)
        write (noutpt,1240) ux8a(1:j3),usblkf(1:j4)
        write (nttyo,1240) ux8a(1:j3),usblkf(1:j4)
1240 format(/' * Error - (EQPT/pcraq) Have encountered a blank',' species name',/7x,'for species block ',a,' of the ',a,' superblock.')

        if (nat .gt. 1) then
            ux24 = uaqsp(nat - 1)
            j5 = ilnobl(ux24)

            if (j5 .gt. 0) then
                write (noutpt,1250) ux24(1:j5)
                write (nttyo,1250) ux24(1:j5)
1250 format(7x,'This block follows the one for ',a,'.')
            end if
        end if

        ux24 = '<blank>'
        nerr = nerr + 1
    end if

    uaqsp(nat) = ux24
    uspec(ns) = ux24

    ! Save the species name for screen and SLIST output.
    uspn(nspnx) = uspec(ns)

    ! Write the species name on the OUTPUT file.
    j2 = ilnobl(uspec(ns))
    write (noutpt,1150) nn,uspec(ns)(1:j2)

    ! Write species names two per line on the SLIST file.
    nspnx = nspnx + 1

    if (nspnx .gt. 2) then
        nspnx = 1
        nnx = nnx + 2
        nmodx = mod(nnx,nmodwr)

        if (nmodwr .eq. 1) then
            nmodx = 1
        end if

        j3 = ilnobl(uspn(2))

        if (nmodx .eq. 1) then
            write (nttyo,1160) nnx,uspn(1),uspn(2)(1:j3)
        end if

        write (nslist,1160) nnx,uspn(1),uspn(2)(1:j3)
    end if

    ! Skip the obsolete 'sp.type =' and 'revised =' lines.
    ! Read the electrical charge.
130 continue
    read (ndat0s,1000,end=990,err=995) uline
    ii = index(uline,'charge')

    if (ii .eq. 0) then
        go to 130
    end if

    ux80 = uline(ii + 6:80)
    ii = index(ux80,'=')
    ux16 = ux80(ii + 1:80)
    call lejust(ux16)
    read (ux16,'(f5.1)',err=995) zx

    ! 1280 format(14x,f5.1)
    zaqsp(nat) = zx
    zchar(ns) = zx

    ! Read the number of chemical elements composing the species.
    read (ndat0s,1000,end=990,err=995) uline
    ux80 = uline
    call lejust(ux80)
    ii = index(ux80,'element(s)')

    if (ii .gt. 1) then
        ux16 = ux80(1:ii - 1)
        read (ux16,'(i5)',err=995) ncts

        ! 1330   format(4x,i2)
    else
        ii = index(ux80,' ')

        if (ii .gt. 1) then
            ux16 = ux80(1:ii - 1)
            read (ux16,'(i5)',err=995) ncts
        else
            ncts = 0
        end if
    end if

    if (ncts .gt. nct) then
        j2 = ilnobl(uspec(ns))
        write (ux8a,'(i5)') ncts
        call lejust(ux8a)
        j3 = ilnobl(ux8a)
        write (ux8c,'(i5)') nct
        call lejust(ux8c)
        j5 = ilnobl(ux8c)
        write (noutpt,1350) uspec(ns)(1:j2),ux8a(1:j3),ux8c(1:j5)
        write (nttyo,1350) uspec(ns)(1:j2),ux8a(1:j3),ux8c(1:j5)
1350 format(/' * Error - (EQPT/pcraq) Species "',a,'" is',' composed of ',a,/7x,'chemical elements, but there',' are only ',a,' elements on the data file.')

        nerr = nerr + 1
    end if

    if (ncts .le. 0) then
        if (uspec(ns)(1:3).ne.'e- ' .and. uspec(ns)(1:3).ne.'E- ') then
            j2 = ilnobl(uspec(ns))
            write (noutpt,1360) uspec(ns)(1:j2)
            write (nttyo,1360) uspec(ns)(1:j2)
1360 format(/' * Error - (EQPT/pcraq) Species "',a,'" is not',' composed of',/7x,'any chemical elements.')

            nerr = nerr + 1
        end if
    end if

    ! Read the names of the elements and the corresponding
    ! compositional coefficients.
    n = 0

    do i = 1,ncts/3
        read (ndat0s,1000,end=990,err=995) uline
        ulbufa = uline

        do k = 1,3
            call lejust(ulbufa)
            ii = index(ulbufa,' ')

            if (ii .gt. 1) then
                ux16 = ulbufa(1:ii - 1)
                read (ux16,'(f8.4)',err=995) cessi(n + k)
            else
                cessi(n + k) = 0.
            end if

            ulbufb = ulbufa(ii + 1:80)
            call lejust(ulbufb)
            uessi(n + k) = ulbufb(1:8)

            if (k .lt. 3) then
                ulbufa = ulbufb(9:80)
            end if
        end do

        !        read (uline,1370,err=995) (cessi(n + k),uessi(n + k), k = 1,3)
        ! 1370   format((4x,3(f8.4,1x,a8,5x)))
        n = n + 3
    end do

    j = mod(ncts,3)

    if (j .gt. 0) then
        read (ndat0s,1000,end=990,err=995) uline
        ulbufa = uline

        do k = 1,j
            call lejust(ulbufa)
            ii = index(ulbufa,' ')

            if (ii .gt. 1) then
                ux16 = ulbufa(1:ii - 1)
                read (ux16,'(f8.4)',err=995) cessi(n + k)
            else
                cessi(n + k) = 0.
            end if

            ulbufb = ulbufa(ii + 1:80)
            call lejust(ulbufb)
            uessi(n + k) = ulbufb(1:8)

            if (k .lt. j) then
                ulbufa = ulbufb(9:80)
            end if
        end do

        ! read (uline,1370,err=995) (cessi(n + k),uessi(n + k), k = 1,j)
        n = n + j
    end if

    ! Check for blank element names, duplicate element names, and
    ! zero-valued stoichiometric coefficients.
    call elesck(cessi,nbtmx1,nctmax,ncts,nentei,nerr,noutpt,ns,nttyo,qblkes,qzeres,uessi,usblkf,uspec)

    if (ns .eq. nbw) then
        ! Make sure that H2O is the first aqueous species.
        if (uspec(nbw)(1:4) .ne. 'H2O ') then
            j2 = ilnobl(uspec(nbw))
            write (noutpt,1420) uspec(nbw)(1:j2)
            write (nttyo,1420) uspec(nbw)(1:j2)
1420 format(/' * Error - (EQPT/pcraq) The first strict basis',' species on the data file',/7x,'is ',a,'. The first',' strict basis species must be H2O.',/7x,'This is an',' idiosyncrasy of EQ3/6.')

            nerr = nerr + 1
        end if
    end if

    if (ns .le. nsb) then
        ! A strict basis species may be composed of only one chemical
        ! element other than O or H.
        ncount = 0

        do i = 1,ncts
            if (cessi(i) .gt. 0.) then
                if (uessi(i)(1:2).ne.'O ' .and. uessi(i)(1:2).ne.'H ' .and.      uessi(i)(1:7).ne.'<blank>') then
                    ncount = ncount + 1
                end if
            end if
        end do

        if (ncount .gt. 1) then
            j2 = ilnobl(uspec(ns))
            write (noutpt,1430) uspec(ns)(1:j2)
            write (nttyo,1430) uspec(ns)(1:j2)
1430 format(/' * Error - (EQPT/pcraq) The strict basis species ',a,' is composed',/7x,'of more than one chemical element',' other than O and H.')

            nerr = nerr + 1
        end if
    end if

    if (ns .le. nct) then
        ! Make sure that there is a 1:1 mapping to the chemical elements
        ! for the first nct members of the strict basis set.
        ! H2O (solvent water) is a special case. H2O always maps to O.
        ! Since the first chemical element is required to be O, the first
        ! aqueous species is required to  be H2O. Here it is sufficient
        ! to make sure that that indeed is the case.
        if (ns .eq. 1) then
            if (uspec(ns)(1:4) .ne. 'H2O ') then
                j2 = ilnobl(uspec(ns))
                write (noutpt,1440) uspec(ns)(1:j2)
                write (nttyo,1440) uspec(ns)(1:j2)
1440 format(/' * Error - (EQPT/pcraq) The first aqueous species',' is ',a,'.',/7x,'The first aqueous species must be H2O.')

                nerr = nerr + 1
            end if
        end if

        ! The hydrogen ion is also a special case. H+ always maps to H.
        if (uspec(ns)(1:3) .eq. 'H+ ') then
            nc = ns

            if (uelem(nc)(1:2) .ne. 'H ') then
                j2 = ilnobl(uspec(ns))
                write (noutpt,1450)
                write (nttyo,1450)
1450 format(/' * Error - (EQPT/pcraq) The species H+ is',' required to match',/7x,'the chemical element H. If the',' species is the n-th member of',/7x,'the strict',' basis set, the corresponding chemical element must be',/7x,'the n-th chemical element (there is one exception',' to this rule,',/7x,"the redox species, which doesn't",' correspond to any chemical element).')

                nerr = nerr + 1
            end if
        end if

        ! Any strict basis species, other than the redox species,
        ! must be composed of the chemical element to which it
        ! formally corresponds.
        nc = ns

        do i = 1,ncts
            if (uessi(i)(1:8) .eq. uelem(nc)(1:8)) then
                if (cessi(i) .gt. 0.) then
                    go to 120
                end if
            end if
        end do

        j2 = ilnobl(uspec(ns))
        j3 = ilnobl(uelem(nc))
        write (noutpt,1460) uspec(ns)(1:j2),uelem(nc)(1:j3)
        write (nttyo,1460) uspec(ns)(1:j2),uelem(nc)(1:j3)
1460 format(/' * Error - (EQPT/pcraq) The strict basis',' species ',a,' is not',/7x,'composed of the chemical',' element ',a,' to which it must correspond.',/7x,'If the species is the n-th member of the strict',' basis set,',/7x,'the, corresponding chemical',' element must be the n-th chemical element.',/7x,'There is one exception to this rule, the',' redox species. It must be',/7x,'the last strict'," basis species, if present. It doesn't correspond",/7x,'to any chemical element.')

        nerr = nerr + 1
120 continue
    end if

    ! Check out the redox species.
    if (ns .eq. nsb) then
        if (uspec(ns)(1:6).ne.'O2(g) ' .and.    uspec(ns)(1:3).ne.'e- ' .and.    uspec(ns)(1:3).ne.'E- ') then
            j2 = ilnobl(uspec(ns))
            write (noutpt,1470) uspec(ns)(1:j2)
            write (nttyo,1470) uspec(ns)(1:j2)
1470 format(/' * Error - (EQPT/pcraq) The redox species is ',a,'. This must',/7x,'be either O2(g) or e-.')

            nerr = nerr + 1
        end if
    end if

    if (ns .eq. nsb) then
        if (nerr .gt. 0) then
            stop
        end if
    end if

    if (ns .eq. nbt) then
        if (nerr .gt. 0) then
            stop
        end if
    end if

    if (ns .gt. nsb) then
        ! Read number of species in the associated reaction.
        read (ndat0s,1000,end=990,err=995) uline
        ux80 = uline
        call lejust(ux80)
        ii = index(ux80,'species')

        if (ii .gt. 1) then
            ux16 = ux80(1:ii - 1)
            read (ux16,'(i5)',err=995) ndrsts

            ! 1500   format(4x,i2)
        else
            ii = index(ux80,' ')

            if (ii .gt. 1) then
                ux16 = ux80(1:ii - 1)
                read (ux16,'(i5)',err=995) ndrsts
            else
                ndrsts = 0
            end if
        end if

        if (ndrsts .gt. nbt1) then
            j2 = ilnobl(uspec(ns))
            write (ux8a,'(i5)') ndrsts
            call lejust(ux8a)
            j3 = ilnobl(ux8a)
            write (ux8b,'(i5)') nbt
            call lejust(ux8b)
            j4 = ilnobl(ux8b)
            write (ux8c,'(i5)') nbt1
            call lejust(ux8c)
            j5 = ilnobl(ux8c)
            write (noutpt,1520) uspec(ns)(1:j2),ux8a(1:j3),ux8b(1:j4),ux8c(1:j5)
            write (nttyo,1520) uspec(ns)(1:j2),ux8a(1:j3),ux8b(1:j4),ux8c(1:j5)
1520 format(/' * Error - (EQPT/pcraq) The reaction for the',/7x,'destruction of species "',a,'" includes ',a,/7x,'species, but there are only ',a,' basis species on the',/7x,'data file, so only ',a,' species may appear in the ','reaction.')

            nerr = nerr + 1
        end if

        if (ndrsts .lt. 2) then
            j2 = ilnobl(uspec(ns))
            j4 = ilnobl(usblkf)
            write (noutpt,1530) uspec(ns)(1:j2),usblkf(1:j4)
            write (nttyo,1530) uspec(ns)(1:j2),usblkf(1:j4)
1530 format(/' * Error - (EQPT/pcraq) The species ',a,' appearing',/7x,'on the data file in the ',a,' superblock',' has fewer than',/7x,'two species in its associated',' reaction. This is not a',/7x,'valid reaction.')

            nerr = nerr + 1
        end if

        ! Read the names of the species in the reaction and the
        ! corresponding reaction coefficients.
        n = 0

        do i = 1,ndrsts/2
            read (ndat0s,1000,end=990,err=995) uline
            ulbufa = uline

            do k = 1,2
                call lejust(ulbufa)
                ii = index(ulbufa,' ')

                if (ii .gt. 1) then
                    ux16 = ulbufa(1:ii - 1)
                    read (ux16,'(f10.4)',err=995) cdrsi(n + k)
                else
                    cdrsi(n + k) = 0.
                end if

                ulbufb = ulbufa(ii + 1:80)
                call lejust(ulbufb)
                udrsi(n + k) = ulbufb(1:24)

                if (k .lt. 2) then
                    ulbufa = ulbufb(25:80)
                end if
            end do

            !          read (uline,1540,err=995) (cdrsi(n + k),udrsi(n + k), k = 1,2)
            ! 1540     format((2(1x,f10.4,2x,a24)))
            n = n + 2
        end do

        j = mod(ndrsts,2)

        if (j .gt. 0) then
            read (ndat0s,1000,end=990,err=995) uline
            ulbufa = uline

            do k = 1,j
                call lejust(ulbufa)
                ii = index(ulbufa,' ')

                if (ii .gt. 1) then
                    ux16 = ulbufa(1:ii - 1)
                    read (ux16,'(f10.4)',err=995) cdrsi(n + k)
                else
                    cdrsi(n + k) = 0.
                end if

                ulbufb = ulbufa(ii + 1:80)
                call lejust(ulbufb)
                udrsi(n + k) = ulbufb(1:24)

                if (k .lt. j) then
                    ulbufa = ulbufb(25:80)
                end if
            end do

            ! read (uline,1540,err=995) (cdrsi(n + k),udrsi(n + k), k = 1,j)
            n = n + j
        end if

        ! Check for blank species names, duplicate species names, and
        ! zero-valued reaction coefficients. Check that the reaction
        ! coefficient of the species with which the reaction is
        ! associated has a negative value.
        call rxnsck(nbtmx1,cdrsi,nct,ndrsts,nentri,nerr,noutpt,ns,nsb,nttyo,qblkrs,qzerrs,udrsi,usblkf,uspec)

        ! Read the log K grid for the current species.
        ! Return the data in the xdbval holding array.
        call rdgrid(ndat0s,ndbmax,ndbptg,ndbptl,qend,qerr,q500nd,udbval,xdbval)

        if (qend) then
            go to 990
        end if

        if (qerr) then
            go to 995
        end if

        ! Load the data into the xlks array.
        ! Calling sequence substitutions:
        !   ns for ipc
        !   nbtmx1 for ipcmax
        !   xlks for zdbval
        call ldbar3(ns,nbtmx1,nacdpr,narxmx,narxt,ndbmax,ntprmx,ntprt,xdbval,xlks)

        if (itgenf .ge. 0) then
            ! Test the grid ranges for sparseness of actual data.
            ustrgr = 'the log K for ' // uspec(ns)
            call tegrid(itgenf,nacdpr,narxt,nerr,noutpt,ntprmx,ntprt,nttyo,nwarn,ustrgr)
        end if
    end if

    ! In the case of a strict basis species (a species which has no
    ! real associated reaction), xhfs is the partial molar enthalpy
    ! of formation and xvfs is the partial molar volume. In the case
    ! of all other species, there is a real associated reaction and
    ! xhfs is then the partial molar enthalpy of reaction and xvfs
    ! is the partial molar volume of reaction.
    if (ipch .ge. 0) then
        ! Read the enthalpy function grid for the current species.
        ! Return the data in the xdbval holding array.
        call rdgrid(ndat0s,ndbmax,ndbptg,ndbptl,qend,qerr,q500nd,udbval,xdbval)

        if (qend) then
            go to 990
        end if

        if (qerr) then
            go to 995
        end if

        ! Load the data into the xhfs array.
        ! Calling sequence substitutions:
        !   ns for ipc
        !   nbtmx1 for ipcmax
        !   xhfs for zdbval
        call ldbar3(ns,nbtmx1,nacdpr,narxmx,narxt,ndbmax,ntprmx,ntprt,xdbval,xhfs)

        if (itgenf .ge. 0) then
            ! Test the grid ranges for sparseness of actual data.
            ustrgr = 'the enthalpy function for ' // uspec(ns)
            call tegrid(itgenf,nacdpr,narxt,nerr,noutpt,ntprmx,ntprt,nttyo,nwarn,ustrgr)
        end if

        do ipc = 1,ipch
            ! Read the enthalpy function derivative grid for the current
            ! species. Return the data in the xdbval holding array.
            call rdgrid(ndat0s,ndbmax,ndbptg,ndbptl,qend,qerr,q500nd,udbval,xdbval)

            if (qend) then
                go to 990
            end if

            if (qerr) then
                go to 995
            end if

            ! Load the data into the dhfs array.
            ! Calling sequence substitutions:
            !   ipchmx for ipcmax
            !   dhfs for zdbval
            call ldbar4(ipc,ipchmx,nacdpr,narxmx,narxt,nbtmx1,ndbmax,ns,ntprmx,ntprt,xdbval,dhfs)
        end do
    end if

    if (ipcv .ge. 0) then
        ! Read the volume function grid for the current species.
        ! Return the data in the xdbval holding array.
        call rdgrid(ndat0s,ndbmax,ndbptg,ndbptl,qend,qerr,q500nd,udbval,xdbval)

        if (qend) then
            go to 990
        end if

        if (qerr) then
            go to 995
        end if

        ! Load the data into the xvfs array.
        ! Calling sequence substitutions:
        !   ns for ipc
        !   nbtmx1 for ipcmax
        !   xvfs for zdbval
        call ldbar3(ns,nbtmx1,nacdpr,narxmx,narxt,ndbmax,ntprmx,ntprt,xdbval,xvfs)

        if (itgenf .ge. 0) then
            ! Test the grid ranges for sparseness of actual data.
            ustrgr = 'the volume function for ' // uspec(ns)
            call tegrid(itgenf,nacdpr,narxt,nerr,noutpt,ntprmx,ntprt,nttyo,nwarn,ustrgr)
        end if

        do ipc = 1,ipcv
            ! Read the volume function derivative grid for the current
            ! species. Return the data in the xdbval holding array.
            call rdgrid(ndat0s,ndbmax,ndbptg,ndbptl,qend,qerr,q500nd,udbval,xdbval)

            if (qend) then
                go to 990
            end if

            if (qerr) then
                go to 995
            end if

            ! Load the data into the dvfs array.
            ! Calling sequence substitutions:
            !   ipcvmx for ipcmax
            !   dvfs for zdbval
            call ldbar4(ipc,ipcvmx,nacdpr,narxmx,narxt,nbtmx1,ndbmax,ns,ntprmx,ntprt,xdbval,dvfs)
        end do
    end if

    ! Check the redox species in the strict basis set. Set the
    ! appropriate flags.
    if (ns .eq. nsb) then
        qelect = uspec(ns)(1:3).eq.'e- '
    end if

    ! If the strict basis redox species is e-, change it to O2(g).
    if (ns .eq. nsb) then
        if (qelect) then
            uspec(ns) = 'O2(g)'
            ncts = 1
            uessi(1) = 'O'
            cessi(1) = 2
            zchar(ns) = 0.
        end if
    end if

    ! If the reaction is written in terms of e-, rewrite it in terms
    ! of O2(g).
    if (qelect) then
        if (ns .gt. nsb) then
            call etoo2(cdrsi,dhfe,dhfs,dvfe,dvfs,ipch,ipchmx,ipcv,ipcvmx,narxmx,narxt,nbtmx1,ndrsts,ns,ntprmx,ntprt,udrsi,xhfe,xhfs,xlke,xlks,xvfe,xvfs)
        end if
    end if

    ! Process the species data.
    mwtss = 0.
    qnofes = .false.

    do n = 1,ncts
        ! Search for element name in the uelem array.
        ux8 = uessi(n)(1:8)

        if (ux8(1:7) .ne. '<blank>') then
            do nc = 1,nct
                if (ux8(1:8) .eq. uelem(nc)(1:8)) then
                    cess(nc,ns) = cess(nc,ns) + cessi(n)
                    mwtss = mwtss + atwt(nc)*cessi(n)
                    go to 220
                end if
            end do

            ! Error, not found.
            if (ns .gt. nsb) then
                j2 = ilnobl(uspec(ns))
                j3 = ilnobl(ux8)
                write (noutpt,1760) ux8(1:j3),uspec(ns)(1:j2)
                write (nttyo,1760) ux8(1:j3),uspec(ns)(1:j2)
1760 format(/' * Error - (EQPT/pcraq) Unrecognized chemical',' element "',a,'" is listed',/7x,'in the composition of',' the species ',a,'.')

                nerr = nerr + 1
                qnofes = .true.
            end if
        end if

220 continue
    end do

    write (ndata1) uspec(ns),ncts,ndrsts
    write (ndat1f,1800) uspec(ns),ncts,ndrsts
1800 format(a24,2x,2i5)

    write (ndata1) mwtss,zchar(ns)
    write (ndat1f,1810) mwtss,zchar(ns)
1810 format(5x,f10.3,f5.0)

    if (ncts .gt. 0) then
        write (ndata1) (cessi(n),uessi(n), n = 1,ncts)
        n = 1
230 continue

        if (n .eq. ncts) then
            j4 = ilnobl(uessi(n))
            write (ndat1f,1820) cessi(n),uessi(n)(1:j4)
1820 format(1x,f10.4,2x,a)

            n = n + 1
        else
            j4 = ilnobl(uessi(n + 1))
            write (ndat1f,1830) cessi(n),uessi(n),cessi(n + 1),uessi(n + 1)(1:j4)
1830 format(1x,f10.4,2x,a8,1x,f10.4,2x,a)

            n = n + 2
        end if

        if (n .le. ncts) then
            go to 230
        end if
    end if

    if (ns .le. nsb) then
        go to 270
    end if

    write (ndata1) (cdrsi(n),udrsi(n), n = 1,ndrsts)
    n = 1

240 continue

    if (n .eq. ndrsts) then
        j4 = ilnobl(udrsi(n))
        write (ndat1f,1860) cdrsi(n),udrsi(n)(1:j4)
1860 format(1x,f10.4,2x,a)

        n = n + 1
    else
        j4 = ilnobl(udrsi(n + 1))
        write (ndat1f,1870) cdrsi(n),udrsi(n),cdrsi(n + 1),udrsi(n + 1)(1:j4)
1870 format(1x,f10.4,2x,a24,1x,f10.4,2x,a)

        n = n + 2
    end if

    if (n .le. ndrsts) then
        go to 240
    end if

    qnofrs = .false.
    cdrs(nbt1,ns) = cdrsi(1)

    do n = 2,ndrsts
        ! Search for species name in the uspec array.
        ux24 = udrsi(n)

        if (ux24(1:7) .ne. '<blank>') then
            do nse = 1,nbt
                if (ux24 .eq. uspec(nse)) then
                    cdrs(nse,ns) = cdrs(nse,ns) + cdrsi(n)
                    go to 250
                end if
            end do

            ! Error, not found.
            j2 = ilnobl(uspec(ns))
            j3 = ilnobl(ux24)
            write (noutpt,1880) uspec(ns)(1:j2),ux24(1:j3)
            write (nttyo,1880) uspec(ns)(1:j2),ux24(1:j3)
1880 format(/' * Error - (EQPT/pcraq) The reaction which destroys',/7x,'non-basis species ',a,' is written in terms of an',/7x,'unrecognized basis species "',a,'".')

            nerr = nerr + 1
            qnofrs = .true.
250 continue
        end if
    end do

    ! Test the reaction for mass and charge balance.
    ! Skip if there are already obvious problems.
    if (ndrsts .ge. 2) then
        if (.not.qblkes .and. .not.qnofes) then
            if (.not.qblkrs .and. .not.qnofrs) then
                call rxnchk(cdrs,cess,mtotr,nbt,nbtmx1,nbtmx2,nco,nct,nctmax,nerr,noutpt,ns,nsb,nttyo,uelem,uspec,zchar)
            end if
        end if
    end if

    ! Fit interpolating polynomials to the log K grid.
    do ntpr = 1,ntprt
        do n = 1,narxt(ntpr)
            avgrid(n,ntpr) = xlks(n,ntpr,ns)
        end do
    end do

    call intrp(aamatr,apr,avgrid,cof,eps100,gmmatr,ipivot,narxmx,narxt,noutpt,ntprmx,ntprt,nttyo,tempc,tempcs,tmpcmx,xvec,yvec)

    ux24 = 'Log K'
    j2 = ilnobl(ux24)
    write (ndata1) ux24
    write (ndat1f,1010) ux24(1:j2)

    do ntpr = 1,ntprt
        nt = narxt(ntpr)
        write (ndata1) (apr(i,ntpr), i = 1,nt)
        write (ndat1f,1890) (apr(i,ntpr), i = 1,nt)
1890 format( 5(1pe16.9) )
    end do

270 continue

    if (ipch .ge. 0) then
        ! Fit interpolating polynomials to the enthalpy function grid.
        do ntpr = 1,ntprt
            do n = 1,narxt(ntpr)
                avgrid(n,ntpr) = xhfs(n,ntpr,ns)
            end do
        end do

        call intrp(aamatr,apr,avgrid,cof,eps100,gmmatr,ipivot,narxmx,narxt,noutpt,ntprmx,ntprt,nttyo,tempc,tempcs,tmpcmx,xvec,yvec)

        ux24 = 'Enthalpy'
        j2 = ilnobl(ux24)
        write (ndata1) ux24
        write (ndat1f,1010) ux24(1:j2)

        do ntpr = 1,ntprt
            nt = narxt(ntpr)
            write (ndata1) (apr(i,ntpr), i = 1,nt)
            write (ndat1f,1890) (apr(i,ntpr), i = 1,nt)
        end do

        ! Fit interpolating polynomials to the grids for its
        ! pressure derivatives.
        do ipc = 1,ipch
            do ntpr = 1,ntprt
                do n = 1,narxt(ntpr)
                    avgrid(n,ntpr) = dhfs(n,ntpr,ipc,ns)
                end do
            end do

            call intrp(aamatr,apr,avgrid,cof,eps100,gmmatr,ipivot,narxmx,narxt,noutpt,ntprmx,ntprt,nttyo,tempc,tempcs,tmpcmx,xvec,yvec)

            ux24 = 'dhfs( )'
            write (ux24(6:6),'(i1)') ipc
            j2 = ilnobl(ux24)
            write (ndata1) ux24
            write (ndat1f,1010) ux24(1:j2)

            do ntpr = 1,ntprt
                nt = narxt(ntpr)
                write (ndata1) (apr(i,ntpr), i = 1,nt)
                write (ndat1f,1890) (apr(i,ntpr), i = 1,nt)
            end do
        end do
    end if

    if (ipcv .ge. 0) then
        ! Fit interpolating polynomials to the volume function grid.
        do ntpr = 1,ntprt
            do n = 1,narxt(ntpr)
                avgrid(n,ntpr) = xvfs(n,ntpr,ns)
            end do
        end do

        call intrp(aamatr,apr,avgrid,cof,eps100,gmmatr,ipivot,narxmx,narxt,noutpt,ntprmx,ntprt,nttyo,tempc,tempcs,tmpcmx,xvec,yvec)

        ux24 = 'Volume'
        j2 = ilnobl(ux24)
        write (ndata1) ux24
        write (ndat1f,1010) ux24(1:j2)

        do ntpr = 1,ntprt
            nt = narxt(ntpr)
            write (ndata1) (apr(i,ntpr), i = 1,nt)
            write (ndat1f,1890) (apr(i,ntpr), i = 1,nt)
        end do

        ! Fit interpolating polynomials to the grids for its
        ! pressure derivatives.
        do ipc = 1,ipcv
            do ntpr = 1,ntprt
                do n = 1,narxt(ntpr)
                    avgrid(n,ntpr) = dvfs(n,ntpr,ipc,ns)
                end do
            end do

            call intrp(aamatr,apr,avgrid,cof,eps100,gmmatr,ipivot,narxmx,narxt,noutpt,ntprmx,ntprt,nttyo,tempc,tempcs,tmpcmx,xvec,yvec)

            ux24 = 'dvfs( )'
            write (ux24(6:6),'(i1)') ipc
            j2 = ilnobl(ux24)
            write (ndata1) ux24
            write (ndat1f,1010) ux24(1:j2)

            do ntpr = 1,ntprt
                nt = narxt(ntpr)
                write (ndata1) (apr(i,ntpr), i = 1,nt)
                write (ndat1f,1890) (apr(i,ntpr), i = 1,nt)
            end do
        end do
    end if

    write (ndat1f,1010) utermc(1:72)

    ! Look for another species block or terminator line.
    read (ndat0s,1000,end=990,err=995) uline
    nn = nn + 1
    go to 100

990 continue
    write (noutpt,2000)
    write (nttyo,2000)
2000 format(/' * Error - (EQPT/pcraq) Unexpectedly encountered',/7x,'end-of-file while reading the DATA0 file.')

    write (noutpt,2010) nn
    write (nttyo,2010) nn
2010 format(7x,'The value of the species block counter is ',i3,'.')

    if (ns .gt. 0) then
        j2 = ilnobl(uspec(ns))

        if (j2 .gt. 0) then
            write (noutpt,2020) uspec(ns)(1:j2)
            write (nttyo,2020) uspec(ns)(1:j2)
2020 format(7x,'The last species name read was "',a,'".')
        else
            nsm1 = ns - 1

            if ((nsm1) .gt. 0) then
                j2 = ilnobl(uspec(nsm1))

                if (j2 .gt. 0) then
                    write (noutpt,2020) uspec(nsm1)(1:j2)
                    write (nttyo,2020) uspec(nsm1)(1:j2)
                end if
            end if
        end if
    end if

    j2 = ilnobl(uline)

    if (j2 .gt. 0) then
        j2 = min(j2,70)
        write (noutpt,2030) uline(1:j2)
        write (nttyo,2030) uline(1:j2)
2030 format(7x,'The last line read was the following:',/7x,'"',a,'"')
    end if

    stop

995 continue
    write (noutpt,2040)
    write (nttyo,2040)
2040 format(/' * Error - (EQPT/pcraq) Encountered a read format',/7x,'error while reading the DATA0 file.')

    write (noutpt,2010) nn
    write (nttyo,2010) nn

    if (ns .gt. 0) then
        j2 = ilnobl(uspec(ns))

        if (j2 .gt. 0) then
            write (noutpt,2020) uspec(ns)(1:j2)
            write (nttyo,2020) uspec(ns)(1:j2)
        else
            nsm1 = ns - 1

            if ((nsm1) .gt. 0) then
                j2 = ilnobl(uspec(nsm1))

                if (j2 .gt. 0) then
                    write (noutpt,2020) uspec(nsm1)(1:j2)
                    write (nttyo,2020) uspec(nsm1)(1:j2)
                end if
            end if
        end if
    end if

    j2 = ilnobl(uline)

    if (j2 .gt. 0) then
        j2 = min(j2,67)
        write (noutpt,2030) uline(1:j2)
        write (nttyo,2030) uline(1:j2)
    end if

    stop

999 continue
end subroutine pcraq