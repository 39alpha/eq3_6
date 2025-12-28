subroutine rd3w8(cgexj,cgexpi,covali,ehi,egexsi,fo2lgi,iebal3,ietmax,iodb,iopg,iopr,iopt,irdxc3,itdsf3,itermx,jgext,jetmax,jflgi,jgexti,jpres3,kxmod,mwtges,nbti,nbtmax,ncmpri,net,neti,netmax,ngexti,ninpts,ngexrt,nobswt,nodbmx,nopgmx,noprmx,noptmx,noutpt,nprob,nsbswt,ntitl,ntitmx,nttyo,nxicmx,nxmdmx,nxmod,nxti,nxtimx,pei,press,qend,qgexsh,qrderr,rho,scamas,tdspkg,tdspl,tempc,tgexp,tolbt,toldl,tolspf,ucospi,uebal,ugexj,ugexji,ugexmo,ugexp,ugexpi,ugexr,ugexsi,umemi,uobsw,uredox,usbsw,usoli,uspeci,utitl,uhfgex,uvfgex,uxkgex,uxmod,xbari,xgexsi,xhfgex,xlkgex,xvfgex,xlkmod,zgexj)
    !! This suboutine reads the EQ3NR INPUT file in compact ("W") format
    !! for version 8.0.
    !! This suboutine is a near-clone of EQ3NR/rd3inw.f. However, the
    !! present suboutine embodies only a pure read function (it does
    !! only mimimal checking of what is read, to ensure that what
    !! follows is readable). EQ3NR/rd3inw.f differs in that it also
    !! writes an instant echo of what is read to the EQ3NR output file.
    !! The calling sequence of this suboutine is identical to that of
    !! EQ3NR/rd3inw.f, EQ3NR/rd3ind.f, and XCON3/rd3d8.f.
    !! This suboutine is called by:
    !!   XCON3/xcon3.f
    !! Principal input:
    !!   ninpts = unit number of the stripped INPUT file
    !!   noutpt = unit number of the OUTPUT file
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

    ! Local variable declarations.
    integer :: i
    integer :: iei
    integer :: je
    integer :: jei
    integer :: j2
    integer :: j3
    integer :: n
    integer :: nbi
    integer :: ne
    integer :: nei
    integer :: nn
    integer :: nr1
    integer :: nr2
    integer :: nssoti
    integer :: nxi
    integer :: nxic

    integer :: ilnobl

    logical :: qgexef

    character(len=80) :: uline
    character(len=48) :: ux48
    character(len=8) :: uendit
    character(len=8) :: ux8

    data uendit /'endit.  '/

    ! The following is nonsense to avoid compiler "unused variable"
    ! warnings. Here noutpt and nprob are not actually used. They are
    ! included in the calling sequence to allow it to match that of
    ! EQ3NR/rd3inw.f.
    noutpt = nttyo
    n = nprob
    nprob = n

    ! qend   = .true if the end of the INPUT file has been encountered
    ! qrderr = .true if the current problem can't be read because of
    !            a read format error or a dimensional overflow
    qend = .false.
    qrderr = .false.

    ! Title.
    n = 0
    read (ninpts,1000,end=100,err=990) uline
1000 format(a80)

    go to 110

100 continue
    qend = .true.
    go to 999

110 continue
    if (uline(1:8) .eq. uendit(1:8)) then
        go to 120
    end if

    n = 1
    utitl(n) = uline

    do nn = 2,ntitmx + 1
        read (ninpts,1000,err=990) uline

        if (uline(1:8) .eq. uendit(1:8)) then
            go to 120
        end if

        n = n + 1

        if (n .gt. ntitmx) then
            write (nttyo,1015) ntitmx
1015 format(/' * Error - (XCON3/rd3w8) Have too many lines in',/7x,'the title. The code is only dimensioned for ',i4,/7x,'lines. Reduce the size of the title or increase the',/7x,'dimensioning parameter ntitpa.')

            go to 990
        end if

        utitl(n) = uline
    end do

120 continue
    ntitl = n

    ! Special basis switches.
    read (ninpts,1110,err=990) nsbswt
1110 format(12x,i3)

    do n = 1,nsbswt
        read (ninpts,1130,err=990) usbsw(1,n)
1130 format(9x,a48)

        read (ninpts,1150,err=990) usbsw(2,n)
1150 format(15x,a48)
    end do

    ! Temperature.
    read (ninpts,1170,err=990) tempc
1170 format(12x,e12.5)

    ! Pressure.
    read (ninpts,1190,err=990) jpres3,press
1190 format(12x,i3,/12x,e12.5)

    ! Density.
    read (ninpts,1220,err=990) rho
1220 format(12x,e12.5)

    ! Total dissolved solutes.
    read (ninpts,1224,err=990) itdsf3,tdspkg,tdspl
1224 format(12x,i3,/12x,e12.5,12x,e12.5)

    ! Species to adjust for electrical balance.
    read (ninpts,1500,err=990) iebal3,uebal
1500 format(12x,i3,/12x,a24)

    ! Default redox constraint.
    read (ninpts,1232,err=990) irdxc3
1232 format(12x,i3)

    read (ninpts,1240,err=990) fo2lgi,ehi
1240 format(12x,e12.5,12x,e12.5)

    read (ninpts,1250,err=990) pei,uredox
1250 format(12x,e12.5,12x,a24)

    ! Aqueous basis species and associated constraints.
    nbi = 0

    do nn = 1,nbtmax + 1
        read (ninpts,2500,err=990) ux8,ux48
2500 format(a8,1x,a48)

        j2 = ilnobl(ux48)

        if (ux8(1:8) .eq. uendit(1:8)) then
            go to 400
        end if

        nbi = nbi + 1

        if (nbi .gt. nbtmax) then
            write (nttyo,2520) nbtmax,ux48(1:j2)
2520 format(/' * Error - (XCON3/rd3w8) The number of basis',' species read',/7x,'from the INPUT file exceeded the',' maximum of ',i3,' while',/7x,'trying to read data for',' the species ',a,'.',/7x,'Increase the dimensioning',' parameter nbtpar.')

            go to 990
        end if

        uspeci(nbi) = ux48

        read (ninpts,2540,err=990) jflgi(nbi),covali(nbi)
2540 format(10x,i2,12x,e12.5)

        if (jflgi(nbi).eq.17 .or. jflgi(nbi).eq.18  .or. jflgi(nbi).eq.25) then
            read (ninpts,2570,err=990) ucospi(nbi)
2570 format(10x,a48)
        end if
    end do

400 continue
    nbti = nbi

    ! Ion exchanger creation.
    read (ninpts,1600,err=990) net
1600 format(12x,i3)

    read (ninpts,1610) ux8
1610 format(12x,a8)

    call lejust(ux8)
    call locase(ux8)
    qgexsh = ux8(1:2).eq.'.t' .or. ux8(1:1).eq. 't'

    if (net .gt. netmax) then
        write (nttyo,1620) netmax,net
1620 format(/' * Error - (XCON3/rd3w8) Have exceeded the maximum',' number of ',i3,/7x,'generic ion exchange phases while',' reading the data to create',/7x,'such phases. Increase',' the dimensioning parameter netpar',/7x,'to at least ',i3,'.')

        go to 990
    end if

    do ne = 1,net
        read (ninpts,1630,err=990) ugexp(ne)
1630 format(12x,a24)

        j2 = ilnobl(ugexp(ne))

        read (ninpts,1650,err=990) mwtges(ne)
1650 format(12x,e12.5)

        read (ninpts,1630,err=990) ugexmo(ne)

        read (ninpts,1650,err=990) tgexp(ne)

        read (ninpts,1690,err=990) jgext(ne)
1690 format(12x,i3)

        if (jgext(ne) .gt. jetmax) then
            write (nttyo,1710) jetmax,ugexp(ne)(1:j2),jgext(ne)
1710 format(/' * Error - (XCON3/rd3w8) Have exceeded the maximum',' number of ',i3,/7x,'exchange sites on a generic ion',' exchange phase while reading',/7x,'the data to create',a,'. Increase the',/7x,'dimensioning parameter jetpar to',' at least ',i3,'.')

            go to 990
        end if

        do je = 1,jgext(ne)
            read (ninpts,1730,err=990) ugexj(je,ne)
1730 format(12x,a8)

            j3 = ilnobl(ugexj(je,ne))

            read (ninpts,1750,err=990) cgexj(je,ne),zgexj(je,ne)
1750 format(12x,e12.5,12x,e12.5)

            read (ninpts,1690,err=990) ngexrt(je,ne)

            if (ngexrt(je,ne) .gt. ietmax) then
                write (nttyo,1820) netmax,ugexj(je,ne)(1:j3),ugexp(ne)(1:j2),ngexrt(je,ne)
1820 format(/' * Error - (XCON3/rd3w8) Have exceeded the',' maximum number of ',i3,/7x,'reactions for a site',' belonging to a generic ion exchange',/7x,'phase while',' reading the data for site ',a,' of exchange phase',/7x,a,'. Increase the dimensioning parameter',/7x,'ietpar to at least ',i3,'.')

                go to 990
            end if

            do n = 1,ngexrt(je,ne)
                read (ninpts,1830,err=990) ugexr(n,je,ne)
1830 format(12x,a56)

                read (ninpts,1850,err=990) xlkgex(n,je,ne),uxkgex(n,je,ne)
1850 format(12x,e12.5,12x,a8)

                read (ninpts,1850,err=990) xhfgex(n,je,ne),uhfgex(n,je,ne)
                read (ninpts,1850,err=990) xvfgex(n,je,ne),uvfgex(n,je,ne)
            end do
        end do
    end do

    ! Specified generic ion exchanger compositions.
    read (ninpts,1910,err=990) neti
1910 format(12x,i3)

    if (neti .le. 0) then
        go to 250
    end if

    if (neti .gt. netmax) then
        write (nttyo,1930) netmax,neti
1930 format(/' * Error - (XCON3/rd3w8) Have exceeded the maximum',' number of ',i3,/7x,'generic ion exchange phases while',' reading the data for concentrations',7x,'and compositions',' of such phases. Increase the dimensioning',/7x,'parameter',' netpar to at least ',i3,'.')

        go to 990
    end if

    do nei = 1,neti
        read (ninpts,1940,err=990) ugexpi(nei)
1940 format(12x,a24)

        j2 = ilnobl(ugexpi(nei))

        read (ninpts,1960,err=990) cgexpi(nei)
1960 format(12x,e12.5)

        read (ninpts,1910,err=990) jgexti(nei)

        if (jgexti(nei) .gt. jetmax) then
            write (nttyo,2000) jetmax,ugexpi(nei)(1:j2),jgexti(nei)
2000 format(/' * Error - (XCON3/rd3w8) Have exceeded the maximum',' number of ',i3,/7x,'exchange sites on a generic ion',' exchange phase while reading',/7x,'the data for the',' concentration and composition of',/7x,a,'. Increase the',' dimensioning parameter jetpar',/7x,'to at least ',i3,'.')

            go to 990
        end if

        ! Find the corresponding ne index.
        do ne = 1,net
            j2 = ilnobl(ugexp(ne))
            j3 = ilnobl(ugexpi(nei))

            if (ugexp(nei)(1:j3) .eq. ugexp(ne)(1:j2)) then
                go to 240
            end if
        end do

        write (nttyo,2002) ugexpi(nei)(1:j3)
2002 format(/' * Error - (XCON3/rd3w8) Data are present on the',' input file for the',/7x,'generic ion exchange phase "',a,', but this phase',/7x,"hasn't been previously defined."," Can't determine the requisite model,",/7x,"therefore can't",' finish reading the current data for this phase.')

        stop

240 continue
        j2 = ilnobl(ugexmo(ne))

        if (ugexmo(ne)(1:j2) .eq. 'Gapon' .or.    ugexmo(ne)(1:6) .eq. 'Gapon-' .or.    ugexmo(ne)(1:j2) .eq. 'Vanselow' .or.    ugexmo(ne)(1:9) .eq. 'Vanselow-') then
            ! Input composition is described in terms of equivalent
            ! fractions on the sites.
            qgexef = .true.
        else
            ! Input composition is described in terms of mole fractions
            ! on the sites.
            qgexef = .false.
        end if

        do jei = 1,jgexti(nei)
            read (ninpts,2010,err=990) ugexji(jei,nei)
2010 format(12x,a8)

            j3 = ilnobl(ugexji(jei,nei))

            read (ninpts,1910,err=990) ngexti(jei,nei)

            if (ngexti(jei,nei) .gt. ietmax) then
                write (nttyo,2040) ietmax,ugexji(jei,nei),ugexpi(nei),ngexti(jei,nei)
2040 format(/' * Error - (XCON3/rd3w8) Have exceeded the',' maximum number of ',i3,/7x,'species on an exchange',' site of a generic ion exchange phase',/7x,'while reading',/7x,'the data for the composition of site ',a,' of',/7x,'exchange phase ',a,'. Increase the dimensioning',' parameter ietpar to at least ',i3,'.')

                go to 990
            end if

            if (qgexef) then
                ! The data are in equivalent fractions.
                do iei = 1,ngexti(jei,nei)
                    read (ninpts,2050,err=990) ugexsi(iei,jei,nei),egexsi(iei,jei,nei)
2050 format(12x,a24,10x,e12.5)
                end do
            else
                ! The data are in mole fractions.
                do iei = 1,ngexti(jei,nei)
                    read (ninpts,2050,err=990) ugexsi(iei,jei,nei),xgexsi(iei,jei,nei)
                end do
            end if
        end do
    end do

250 continue

    ! Specified solid solution compositions.
    read (ninpts,2102,err=990) nxti
2102 format(12x,i3)

    if (nxti .gt. nxtimx) then
        write (nttyo,2107) nxtimx,nxti
2107 format(/' * Error - (EQ3NR/rd3w8) Have exceeded the maximum',/7x,'number of ',i5,' solid solutions for which',' compositions',/7x,'may be specified on the INPUT file.',' Increase the dimensioning',/7x,'parameter nxtipa to at',' least ',i3,'.')

        go to 990
    end if

    nxic = 0

    do nxi = 1,nxti
        read (ninpts,2120,err=990) usoli(nxi)
2120 format(12x,a24)

        j2 = ilnobl(usoli(nxi))

        read (ninpts,2140,err=990) nssoti
2140 format(12x,i3)

        ncmpri(1,nxi) = nxic + 1
        ncmpri(2,nxi) = nxic + nssoti
        nr1 = ncmpri(1,nxi)
        nr2 = ncmpri(2,nxi)

        if (nr2 .gt. nxicmx) then
            write (nttyo,2150) nxicmx,usoli(nxi)(1:j2),nr2
2150 format(/' * Error - (EQ3NR/rd3w8) Have exceeded the',' maximum number',/7x,'of ',i5,' solid solution species',' for which mole fractions',/7x,'are specified on the',' INPUT file. This occurred while reading ',/7x,'data',' for the solid solution ',a,'. Increase',/7x,'the',' dimensioning parameter nxicpa to at least ',i3,'.')

            go to 990
        end if

        do nxic = nr1,nr2
            read (ninpts,2160,err=990) umemi(nxic),xbari(nxic)
2160 format(12x,a24,12x,e12.5)
        end do
    end do

    ! Number of nxmod options.
    read (ninpts,2340,err=990) nxmod
2340 format(12x,i3)

    if (nxmod .gt. nxmdmx) then
        write (nttyo,2360) nxmdmx
2360 format(/' * Error - (XCON3/rd3w8) Have too many nxmod',/7x,'alter/suppress options. The code is only dimensioned',/7x,'for ',i3,' such options. Reduce the number of such',/7x,'options or increase the dimensioning parameter nxmdpa.')

        go to 990
    end if

    ! Nxmod options.
    do n = 1,nxmod
        read (ninpts,2370,err=990) uxmod(n),kxmod(n),xlkmod(n)
2370 format(12x,a48,/12x,i2,22x,e12.5)
    end do

    ! Iopt model option switches.
    ! Note: iopt(1) = iopt1, etc.
    read (ninpts,1380,err=990) (iopt(i), i = 1,20)
1380 format(12x,10i5)

    ! Iopg activity coefficient option switches.
    ! Note: iopg(1) = iopg1, etc.
    read (ninpts,1410,err=990) (iopg(i), i = 1,20)
1410 format(12x,10i5)

    ! Iopr print option switches.
    ! Note: iopr(1) = iopr1, etc.
    read (ninpts,1430,err=990) (iopr(i), i = 1,20)
1430 format(12x,10i5)

    ! Iodb debugging print option switches.
    ! Note: iodb(1) = iodb1, etc.
    read (ninpts,1470,err=990) (iodb(i), i = 1,20)
1470 format(12x,10i5)

    ! Numerical parameters.
    read (ninpts,1280,err=990) tolbt,toldl
1280 format(12x,e12.5,12x,e12.5)

    read (ninpts,1300,err=990) itermx
1300 format(12x,i3)

    ! Ordinary basis switches.
    read (ninpts,1110,err=990) nobswt

    do n = 1,nobswt
        read (ninpts,1130,err=990) uobsw(1,n)
        read (ninpts,1150,err=990) uobsw(2,n)
    end do

    ! Saturation flag tolerance.
    read (ninpts,1294,err=990) tolspf
1294 format(12x,e12.5)

    ! Scale factor for the mass of aqueous solution to write
    ! on the PICKUP file.
    read (ninpts,1320,err=990) scamas
1320 format(12x,e12.5)

    ! End of input for one problem.
    go to 999

990 continue
    qrderr = .true.

999 continue
end subroutine rd3w8