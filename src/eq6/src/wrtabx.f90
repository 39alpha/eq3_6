subroutine wrtabx(actlg,afrc1,aft1,alk,cteaq,dvoso,dwoso,eh,fo2lg,iindx1,iktmax,iopt,ipndx1,kmax,km1,kmt,kstep,kx1,kxt,loph,ncmpr,modr,mopht,narn1,mosp,nct,nctmax,noptmx,nptmax,nrct,nrctmx,nstmax,ntabx,ntidmx,ntitl2,ntitld,ntitmx,nxtmax,pe,ph,ppmwe,prcinf,press,prminf,qbye,qmod,qriinf,tempc,time1,uelem,uphase,uplatm,ureac,uspec,usteq6,utitl2,utitld,uveeq6,vodrt,vosoct,wodrt,woh2o,wosoct,xbar,xi1)
    !! This subroutine writes the scratch tab file tabx. The length of
    !! any line should not exceed 129 characters.
    !! This subroutine is called by:
    !!   EQ6/path.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: iktmax
    integer :: kmax
    integer :: nctmax
    integer :: noptmx
    integer :: nptmax
    integer :: nrctmx
    integer :: nstmax
    integer :: ntidmx
    integer :: ntitmx
    integer :: nxtmax

    integer :: ntabx

    integer :: iindx1(kmax)
    integer :: iopt(noptmx)
    integer :: ipndx1(kmax)
    integer :: ncmpr(2,nptmax)

    integer :: km1
    integer :: kmt
    integer :: kstep
    integer :: kx1
    integer :: kxt
    integer :: narn1
    integer :: nct
    integer :: nrct
    integer :: ntitl2
    integer :: ntitld

    logical :: qbye
    logical :: qmod
    logical :: qriinf

    character(len=80) :: utitl2(ntitmx)
    character(len=80) :: utitld(ntidmx)
    character(len=48) :: uspec(nstmax)
    character(len=24) :: uphase(nptmax)
    character(len=24) :: ureac(nrctmx)
    character(len=8) :: uelem(nctmax)
    character(len=8) :: uplatm
    character(len=8) :: usteq6
    character(len=8) :: uveeq6

    real(kind=8) :: actlg(nstmax)
    real(kind=8) :: afrc1(nrctmx)
    real(kind=8) :: cteaq(nctmax)
    real(kind=8) :: loph(nptmax)
    real(kind=8) :: modr(nrctmx)
    real(kind=8) :: mopht(nptmax)
    real(kind=8) :: mosp(nstmax)
    real(kind=8) :: ppmwe(nctmax)
    real(kind=8) :: xbar(nstmax)

    real(kind=8) :: aft1
    real(kind=8) :: alk
    real(kind=8) :: dvoso
    real(kind=8) :: dwoso
    real(kind=8) :: eh
    real(kind=8) :: fo2lg
    real(kind=8) :: wodrt
    real(kind=8) :: wosoct
    real(kind=8) :: pe
    real(kind=8) :: ph
    real(kind=8) :: prcinf
    real(kind=8) :: press
    real(kind=8) :: prminf
    real(kind=8) :: tempc
    real(kind=8) :: time1
    real(kind=8) :: vodrt
    real(kind=8) :: vosoct
    real(kind=8) :: woh2o
    real(kind=8) :: xi1

    ! Local variable declarations with global dimensioning.
    character(len=16), dimension(:,:), allocatable :: uprss
    character(len=8), dimension(:,:), allocatable :: uprmn
    character(len=8), dimension(:), allocatable :: uelac

    real(kind=8), dimension(:), allocatable :: lcteaq
    real(kind=8), dimension(:), allocatable :: ppmaq
    real(kind=8), dimension(:), allocatable :: xfrac
    real(kind=8), dimension(:), allocatable :: zvclgx

    ! Saved values of local array sizes.
    integer :: isv_kmax
    integer :: isv_nctmax
    integer :: isv_nxcmax

    SAVE isv_kmax,isv_nctmax,isv_nxcmax

    SAVE lcteaq,ppmaq,uelac,uprmn,uprss,xfrac,zvclgx

    ! Local variable declarations.
    integer :: nxcmax

    integer :: i
    integer :: j2
    integer :: j3
    integer :: j4
    integer :: k
    integer :: kcol
    integer :: ktot
    integer :: ktotlm
    integer :: n
    integer :: nc
    integer :: nctpr
    integer :: nlim
    integer :: nlim1
    integer :: nlim2
    integer :: np
    integer :: np1
    integer :: nrc
    integer :: ns

    integer :: ilnobl

    logical :: qphasc
    logical :: qstart
    logical :: qtitl

    character(len=1) :: ua
    character(len=1) :: ub
    character(len=1) :: uc
    character(len=1) :: ud
    character(len=1) :: ue
    character(len=1) :: uf
    character(len=1) :: ug
    character(len=1) :: uh
    character(len=1) :: uj
    character(len=1) :: uk
    character(len=1) :: ul
    character(len=1) :: um
    character(len=1) :: un
    character(len=1) :: uo
    character(len=1) :: up
    character(len=1) :: uq
    character(len=1) :: ur
    character(len=1) :: us
    character(len=1) :: ut
    character(len=1) :: uu
    character(len=1) :: uv
    character(len=1) :: ux
    character(len=1) :: uy

    real(kind=8) :: cx
    real(kind=8) :: dx1
    real(kind=8) :: dx2
    real(kind=8) :: dx3
    real(kind=8) :: lalk
    real(kind=8) :: mx
    real(kind=8) :: tdays
    real(kind=8) :: tlogd
    real(kind=8) :: wkgh2o
    real(kind=8) :: xilog

    real(kind=8) :: tlg

    data ua,ub,uc,ud,ue,uf,ug,uh,uj,uk,ul,um,un,uo,up,uq,ur,us,ut,uu,uv,ux,uy/'a','b','c','d','e','f','g','h','j','k','l','m','n','o','p','q','r','s','t','u','v','x','y'/

    nxcmax = nxtmax*iktmax

    ! Allocate or reallocate local work arrays as needed.
    if (.not.ALLOCATED(uprss)) then
        ! Local work arrays are not allocated. Zero the saved
        ! array size variables. Note that only one array is tested
        ! to see if it is allocated. It is assumed that all local
        ! work arrays are either allocated or not.
        isv_kmax = 0
        isv_nctmax = 0
        isv_nxcmax = 0
    else
        ! Local work arrays are allocated. Check to see if any of the
        ! array size variables have changed. If so, deallocate
        ! the corresponding local work arrays and zero the corresponding
        ! saved size variables.
        if (kmax .ne. isv_kmax) then
            DEALLOCATE(uprmn)
            DEALLOCATE(uprss)
            DEALLOCATE(zvclgx)
            isv_kmax = 0
        end if

        if (nctmax .ne. isv_nctmax) then
            DEALLOCATE(uelac)
            DEALLOCATE(lcteaq,ppmaq)
            isv_nctmax = 0
        end if

        if (nxcmax .ne. isv_nxcmax) then
            DEALLOCATE(xfrac)
            isv_nxcmax = 0
        end if
    end if

    ! At this point, the saved array size values are zero if the
    ! corresponding arrays need to be allocated.
    if (isv_kmax .eq. 0) then
        ALLOCATE(uprmn(3,kmax))
        ALLOCATE(uprss(2,kmax))
        ALLOCATE(zvclgx(kmax))
        isv_kmax = kmax
    end if

    if (isv_nctmax .eq. 0) then
        ALLOCATE(uelac(nctmax))
        ALLOCATE(lcteaq(nctmax),ppmaq(nctmax))
        isv_nctmax = nctmax
    end if

    if (isv_nxcmax .eq. 0) then
        ALLOCATE(xfrac(nxcmax))
        isv_nxcmax = nxcmax
    end if

    ! Zero the contents of the local work arrays.
    do k = 1,kmax
        uprmn(1,k) = ' '
        uprmn(2,k) = ' '
        uprmn(3,k) = ' '
        uprss(1,k) = ' '
        uprss(2,k) = ' '
    end do

    do k = 1,kmax
        zvclgx(k) = -99999.
    end do

    do n = 1,nctmax
        uelac(n) = ' '
    end do

    do n = 1,nctmax
        lcteaq(n) = -99999.
        ppmaq(n) = 0.
    end do

    do n = 1,nxcmax
        xfrac(n) = 0.
    end do

    ! Write output tables of some selected parameters. These are
    ! identified by the letters a, b, c, etc. They are written on
    ! the file tabx in scrambled order (i.e., the lines of any given
    ! tables are interspersed with those of the other tables). Each
    ! line is marked in column one with the letter of the table to
    ! which it belongs. When the run is complete, EQLIBU/dscram.f
    ! writes the lines in descrambled form on file tab. File tabx
    ! is preserved. The maximum line length for the tabx file is 129
    ! characters. The maximum line length for the tab file is 128
    ! characters.
    !  Write table a. Header data.
    !           Run title
    !           Code version identification
    !           Data file title
    ! Note: the original primary title (utitl1) has been rolled over
    ! to the current secondary title (utitl2).
    qstart = kstep .eq. 0

    if (qstart) then
        do i = 1,3
            write (ntabx,1000) ua
1000 format(a1,1x)
        end do

        do n = 1,ntitl2
            j2 = ilnobl(utitl2(n))
            write (ntabx,1005) ua,utitl2(n)(1:j2)
1005 format(a1,1x,a)
        end do

        write (ntabx,1000) ua
        j2 = ilnobl(uveeq6)
        j3 = ilnobl(usteq6)
        j4 = ilnobl(uplatm)
        write (ntabx,1010) ua,uveeq6(1:j2),usteq6(1:j3),uplatm(1:j4)
1010 format(a1,1x,'Running EQ3/6-V',a,'-EQ6-EXE-',a,'-',a)

        write (ntabx,1015) ua
1015 format(a1)

        write (ntabx,1000) ua

        do n = 1,ntitld
            j2 = ilnobl(utitld(n))
            write (ntabx,1005) ua,utitld(n)(1:j2)
        end do
    end if

    qtitl = qstart

    if (iopt(18) .ge. 1) then
        qtitl = .false.
    end if

    ! Write table b. Major run parameters.
    !           Overall reaction progress
    !           Log of overall reaction progress
    !           Time, days
    !           Log days
    !           Temperature, C
    !           Pressure, bars
    !           pH
    !           log fO2
    !           Eh, volts
    !           pe
    !           Mass of solvent, kg
    !           Total affinity, kcal
    if (qtitl) then
        do i = 1,3
            write (ntabx,1000) ub
        end do

        write (ntabx,1020) ub
1020 format(a1,'       xi     log xi    time, d  log days','    tempc  ','   press       ph     log fo2      eh','        pe      kg h2o   tot aff')

        write (ntabx,1000) ub
    end if

    ! XX  Will eventually have to revise the tab file (-99999. vs. -999.).
    ! XX   xilog = -99999.
    xilog = -999.

    if (xi1 .gt. 0.) then
        xilog = tlg(xi1)
    end if

    wkgh2o = 0.001*woh2o

    if (qriinf) then
        tdays = prcinf

        ! XX     tlogd = +99999.
        tlogd = +999.
    else
        tdays = time1/86400.

        if (time1 .gt. prminf) then
            tlogd = tlg(tdays)
        else
            ! XX       tlogd = -99999.
            tlogd = -999.
        end if
    end if

    ! XX
    if (eh .le. -99999.) then
        eh = -999.
    end if

    if (pe .le. -99999.) then
        pe = -999.
    end if

    if (aft1 .ge. 99999.) then
        aft1 = 999.
    end if

    ! XX
    write (ntabx,1025) ub,xi1,xilog,tdays,tlogd,tempc,press,ph,fo2lg,eh,pe,wkgh2o,aft1
1025 format(a1,1x,1pe10.3,0pf10.4,1pe10.3,9(0pf10.4))

    ! Write table c. Miscellaneous aqueous solution composition
    ! parameters.
    !           Log of activity of water
    !           Alkalinity, eq/kg.H2O (not defined for
    !             temperatures greater than 50 C)
    if (qtitl) then
        do i = 1,3
            write (ntabx,1000) uc
        end do

        write (ntabx,1030) uc
1030 format(a1,1x,'   log xi   time, d   log days  log alk','  log tot   log tot   log tot   log a h2o')

        write (ntabx,1035) uc
1035 format(a1,41x,'  co3--     so4--     s--')

        write (ntabx,1000) uc
    end if

    ! XX   lalk = -99999.
    lalk = -999.

    if (alk .gt. 0.) then
        lalk = tlg(alk)
    end if

    dx1 = 0.
    dx2 = 0.
    dx3 = 0.
    write (ntabx,1040) uc,xilog,tdays,tlogd,lalk,dx1,dx2,dx3,actlg(narn1)
1040 format(a1,1x,0pf10.4,1pe10.3,10(0pf10.4))

    ! Write tables d, e, f, g, and h. Aqueous solution composition
    ! in terms of molalities of chemical elements. Write
    ! corresponding tables j, k, l, m, and n (aqueous solution
    ! composition in terms of ppm of chemical elements).
    nctpr = 0

    do nc = 1,nct
        if (uelem(nc)(1:2).ne.'H ' .and. uelem(nc)(1:2) .ne. 'O ') then
            nctpr = nctpr + 1
            uelac(nctpr) = uelem(nc)
            cx = cteaq(nc)
            lcteaq(nctpr) = tlg(cx)

            ! XX  Will eventually have to revise the tab file (-99999. vs. -999.).
            if (lcteaq(nctpr) .lt. -999.) then
                lcteaq(nctpr) = -999.
            end if

            ppmaq(nctpr) = ppmwe(nc)
        end if
    end do

    if (nctpr .le. 0) then
        go to 100
    end if

    nlim2 = min(9,nctpr)

    if (qtitl) then
        do i = 1,3
            write (ntabx,1000) ud
            write (ntabx,1000) uj
        end do

        write (ntabx,1045) ud
1045 format(a1,20x,'log molality of dissolved elements')

        write (ntabx,1050) uj
1050 format(a1,17x,'ppm (mg/kg.sol) of dissolved elements')

        write (ntabx,1000) ud
        write (ntabx,1000) uj
        write (ntabx,1055) ud,(uelac(n), n = 1,nlim2)
        write (ntabx,1055) uj,(uelac(n), n = 1,nlim2)
1055 format(a1,1x,'   log xi   time, d   log days ',10(4x,a3,3x))

        write (ntabx,1000) ud
        write (ntabx,1000) uj
    end if

    write (ntabx,1040) ud,xilog,tdays,tlogd,(lcteaq(n), n = 1,nlim2)
    write (ntabx,1065) uj,xilog,tdays,tlogd,(ppmaq(n), n = 1,nlim2)
1065 format(a1,1x,0pf10.4,1pe10.3,0pf10.4,9(1x,g8.3,1x))

    if (nctpr .le. nlim2) then
        go to 100
    end if

    nlim1 = nlim2 + 1
    nlim2 = min(18,nctpr)

    if (qtitl) then
        do i = 1,3
            write (ntabx,1000) ue
            write (ntabx,1000) uk
        end do

        write (ntabx,1055) ue,(uelac(n), n = nlim1,nlim2)
        write (ntabx,1055) uk,(uelac(n), n = nlim1,nlim2)
        write (ntabx,1000) ue
        write (ntabx,1000) uk
    end if

    write (ntabx,1040) ue,xilog,tdays,tlogd,(lcteaq(n), n = nlim1,nlim2)
    write (ntabx,1065) uk,xilog,tdays,tlogd,(ppmaq(n), n = nlim1,nlim2)

    if (nctpr .le. nlim2) then
        go to  100
    end if

    nlim1 = nlim2 + 1
    nlim2 = min(27,nctpr)

    if (qtitl) then
        do i = 1,3
            write (ntabx,1000) uf
            write (ntabx,1000) ul
        end do

        write (ntabx,1055) uf,(uelac(n), n = nlim1,nlim2)
        write (ntabx,1055) ul,(uelac(n), n = nlim1,nlim2)
        write (ntabx,1000) uf
        write (ntabx,1000) ul
    end if

    write (ntabx,1040) uf,xilog,tdays,tlogd,(lcteaq(n), n = nlim1,nlim2)
    write (ntabx,1065) ul,xilog,tdays,tlogd,(ppmaq(n), n = nlim1,nlim2)

    if (nctpr .le. nlim2) then
        go to  100
    end if

    nlim1 = nlim2 + 1
    nlim2 = min(36,nctpr)

    if (qtitl) then
        do i = 1,3
            write (ntabx,1000) ug
            write (ntabx,1000) um
        end do

        write (ntabx,1055) ug,(uelac(n), n = nlim1,nlim2)
        write (ntabx,1055) um,(uelac(n), n = nlim1,nlim2)
        write (ntabx,1000) ug
        write (ntabx,1000) um
    end if

    write (ntabx,1040) ug,xilog,tdays,tlogd,(lcteaq(n), n = nlim1,nlim2)
    write (ntabx,1065) um,xilog,tdays,tlogd,(ppmaq(n), n = nlim1,nlim2)

    if (nctpr .le. nlim2) then
        go to  100
    end if

    nlim1 = nlim2 + 1
    nlim2 = min(45,nctpr)

    if (qtitl) then
        do i = 1,3
            write (ntabx,1000) uh
            write (ntabx,1000) un
        end do

        write (ntabx,1055) uh,(uelac(n), n = nlim1,nlim2)
        write (ntabx,1055) un,(uelac(n), n = nlim1,nlim2)
        write (ntabx,1000) uh
        write (ntabx,1000) un
    end if

    write (ntabx,1040) uh,xilog,tdays,tlogd,(lcteaq(n), n = nlim1,nlim2)
    write (ntabx,1065) un,xilog,tdays,tlogd,(ppmaq(n), n = nlim1,nlim2)

100 continue

    ! Write table o. Product solid solution compositions.
    !              Endmember compositions of product phases.
    ! Write header on first pass through.
    if (qtitl) then
        do i = 1,3
            write (ntabx,1000) uo
        end do

        write (ntabx,2100) uo
2100 format(a1,20x,'solid solution product compositions')

        write (ntabx,1000) uo
    end if

    ! First check if solid solution phases are present in matrix.
    ktot = kxt - kx1 + 1

    if (ktot .le. 0) then
        go to 400
    end if

    ! Go through solid solutions in matrix.
    n = 0

    do kcol = kx1,kxt
        ns = iindx1(kcol)
        np = ipndx1(kcol)
        n = n + 1
        uprss(1,n) = uspec(ns)(1:10)
        uprss(2,n) = uspec(ns)(11:20)
        xfrac(n) = xbar(ns)
    end do

    ! Write out solid solution name and composition.
    nlim = min(6,n)
    write (ntabx,2006) uo,(uprss(1,i),i = 1,nlim)
    write (ntabx,2011) uo,(uprss(2,i),i = 1,nlim)
2006 format(a1,1x,'   log xi ',8x,6(1x,a10,1x))

2011 format(a1,19x,6(1x,a10,1x))

    write (ntabx,2007) uo,xilog,(xfrac(i),i = 1,nlim)
2007 format(a1,4x,f9.4,6x,6(1x,f10.8,1x))

    if (n .gt. nlim) then
        nlim = min(12,n)
        write (ntabx,2008) uo,(uprss(1,i),i = 7,nlim)
        write (ntabx,2008) uo,(uprss(2,i),i = 7,nlim)
        write (ntabx,2009) uo,(xfrac(i),i = 7,nlim)
2008 format(/,a1,19x,5(1x,a10,1x),1x,a10)

2009 format(a1,19x,6(1x,f10.8,1x))
    end if

    if (n .gt. nlim) then
        nlim = min(18,n)
        write (ntabx,2008) uo,(uprss(1,i),i = 13,nlim)
        write (ntabx,2008) uo,(uprss(2,i),i = 13,nlim)
        write (ntabx,2009) uo,(xfrac(i),i = 13,nlim)
    end if

    if (n .gt. nlim) then
        nlim = min(24,n)
        write (ntabx,2008) uo,(uprss(1,i),i = 19,nlim)
        write (ntabx,2008) uo,(uprss(2,i),i = 19,nlim)
        write (ntabx,2009) uo,(xfrac(i),i = 19,nlim)
    end if

    if (n .gt. 24) then
        write (ntabx,2051) uo
    end if

2051 format(a1,1x,'number of solid solution product phases > 24')

400 continue

    ! Write tables p, q, r, and s. Product minerals.
    !           Log of mass (moles) for closed system
    !           Log of cumulative mass (moles) for flow-through
    !              system
    if (qtitl) then
        do i = 1,3
            write (ntabx,1000) up
        end do

        if (iopt(1) .eq. 2) then
            write (ntabx,2000) up
2000 format(a1,20x,'log of moles of product minerals (cumulative)')
        else
            write (ntabx,2002) up
2002 format(a1,20x,'log of moles of product minerals')
        end if

        write (ntabx,1000) up
    end if

    ktot = kxt - km1 + 1

    if (ktot .le. 0) then
        go to 200
    end if

    ! ktotlm = no. of minerals that can be printed in the tables.
    ktotlm = 36

    if (kmt .ge. km1) then
        do kcol = km1,kmt
            n = kcol - km1 + 1

            if (n .gt. ktotlm) then
                go to 125
            end if

            np = ipndx1(kcol)
            uprmn(1,n) = uphase(np)(1:8)
            uprmn(2,n) = uphase(np)(9:16)
            uprmn(3,n) = uphase(np)(17:24)

            if (iopt(1) .eq. 2) then
                zvclgx(n) = tlg(mopht(np))
            else
                zvclgx(n) = loph(np)
            end if

            ! XX  Will eventually have to revise the tab file (-99999. vs. -999.).
            if (zvclgx(n) .lt. -999.) then
                zvclgx(n) = -999.
            end if
        end do
    else
        ! If all products present are solid solutions, reset counter to 0.
        n = 0
    end if

    if (kxt .ge. kx1) then
        np1 = 0

        do kcol = kx1,kxt
            ns = iindx1(kcol)
            np = ipndx1(kcol)

            ! The values of kx1 - kxt correspond to all the endmembers of
            ! solid solutions in the PRS. The following line checks to
            ! see if the next solid solution phase has been reached.
            if (np .ne. np1) then
                np1 = np
                n = n + 1

                if (n .gt. ktotlm) then
                    go to 125
                end if

                uprmn(1,n) = uphase(np)(1:8)
                uprmn(2,n) = uphase(np)(9:16)
                uprmn(3,n) = uphase(np)(17:24)

                if (iopt(1) .eq. 2) then
                    zvclgx(n) = tlg(mopht(np))
                else
                    zvclgx(n) = loph(np)
                end if

                if (zvclgx(n) .lt. -999.) then
                    zvclgx(n) = -999.
                end if
            end if
        end do
    end if

125 continue
    ktot = n

    nlim2 = min(9,ktot)

    qphasc = qmod .or. qbye .or. qtitl

    if (qphasc) then
        do i = 1,3
            write (ntabx,1000) up
        end do

        write (ntabx,2005) up,(uprmn(1,n), n = 1,nlim2)
2005 format(a1,1x,'   log xi   time, d   log days ',9(1x,a8,1x))

        write (ntabx,2010) up,(uprmn(2,n), n = 1,nlim2)
2010 format(a1,32x,9(1x,a8,1x),1x,a8)

        write (ntabx,2010) up,(uprmn(3,n), n = 1,nlim2)
        write (ntabx,1000) up
    end if

    write (ntabx,1040) up,xilog,tdays,tlogd,(zvclgx(n), n = 1,nlim2)

    if (ktot .le. nlim2) then
        go to 200
    end if

    nlim1 = nlim2 + 1
    nlim2 = min(18,ktot)

    if (qphasc) then
        do i = 1,3
            write (ntabx,1000) uq
        end do

        write (ntabx,2005) uq,(uprmn(1,n), n = nlim1,nlim2)
        write (ntabx,2010) uq,(uprmn(2,n), n = nlim1,nlim2)
        write (ntabx,2010) uq,(uprmn(3,n), n = nlim1,nlim2)
        write (ntabx,1000) uq
    end if

    write (ntabx,1040) uq,xilog,tdays,tlogd,(zvclgx(n), n = nlim1,nlim2)

    if (ktot .le. nlim2) then
        go to 200
    end if

    nlim1 = nlim2 + 1
    nlim2 = min(27,ktot)

    if (qphasc) then
        do i = 1,3
            write (ntabx,1000) ur
        end do

        write (ntabx,2005) ur,(uprmn(1,n), n = nlim1,nlim2)
        write (ntabx,2010) ur,(uprmn(2,n), n = nlim1,nlim2)
        write (ntabx,2010) ur,(uprmn(3,n), n = nlim1,nlim2)
        write (ntabx,1000) ur
    end if

    write (ntabx,1040) ur,xilog,tdays,tlogd,(zvclgx(n), n = nlim1,nlim2)

    if (ktot .le. nlim2) then
        go to 200
    end if

    nlim1 = nlim2 + 1
    nlim2 = min(36,ktot)

    if (qphasc) then
        do i = 1,3
            write (ntabx,1000) us
        end do

        write (ntabx,2005) us,(uprmn(1,n), n = nlim1,nlim2)
        write (ntabx,2010) us,(uprmn(2,n), n = nlim1,nlim2)
        write (ntabx,2010) us,(uprmn(3,n), n = nlim1,nlim2)
        write (ntabx,1000) us
    end if

    write (ntabx,1040) us,xilog,tdays,tlogd,(zvclgx(n), n = nlim1,nlim2)

200 continue

    ! Write tables t and u. Reactants.
    !           Log of mass (moles) of reactant destroyed.
    if (nrct .le. 0) then
        go to 215
    end if

    do nrc = 1,nrct
        mx = abs(modr(nrc))

        ! XX     zvclgx(nrc) = -99999.
        zvclgx(nrc) = -999.

        if (mx .gt. 0.) then
            zvclgx(nrc) = tlg(mx)
        end if
    end do

    nlim2 = min(9,nrct)

    if (qtitl) then
        do i = 1,3
            write (ntabx,1000) ut
        end do

        write (ntabx,2020) ut
2020 format(a1,20x,'log of destroyed moles of reactants')

        write (ntabx,1000) ut
        write (ntabx,2005) ut,(ureac(nrc)(1:8), nrc = 1,nlim2)
        write (ntabx,2010) ut,(ureac(nrc)(9:16), nrc = 1,nlim2)
        write (ntabx,2010) ut,(ureac(nrc)(17:24), nrc = 1,nlim2)
        write (ntabx,1000) ut
    end if

    write (ntabx,1040) ut,xilog,tdays,tlogd,(zvclgx(nrc), nrc = 1,nlim2)

    if (nrct .le. nlim2) then
        go to 215
    end if

    nlim1 = nlim2 + 1
    nlim2 = min(18,nrct)

    if (qtitl) then
        do i = 1,3
            write (ntabx,1000) uu
        end do

        write (ntabx,2005) uu,(ureac(nrc)(1:8), nrc = nlim1,nlim2)
        write (ntabx,2010) uu,(ureac(nrc)(9:16), nrc = nlim1,nlim2)
        write (ntabx,2010) uu,(ureac(nrc)(17:24), nrc = nlim1,nlim2)
        write (ntabx,1000) uu
    end if

    write (ntabx,1040) uu,xilog,tdays,tlogd,(zvclgx(nrc), nrc = nlim1,nlim2)

215 continue

    ! Write table v. Overall masses and volumes.
    !           Total mass of solids destroyed, grams
    !           Total mass of solids created, grams
    !           Net mass of solids created, grams
    !           Total volume of solids destroyed, cc
    !           Total volume of solids created, cc
    !           Net volume of solids created, cc
    if (qtitl) then
        do i = 1,3
            write (ntabx,1000) uv
        end do

        write (ntabx,2030) uv
2030 format(a1,1x,'      xi     log xi   time, d   log days   g des','     g cre     g net     cc des    cc cre    cc net')

        write (ntabx,1000) uv
    end if

    write (ntabx,2035) uv,xi1,xilog,tdays,tlogd,wodrt,wosoct,dwoso,vodrt,vosoct,dvoso
2035 format(a1,1x,1pe10.3,0pf10.4,1pe10.3,0pf10.4,8(1pe10.3))

    ! Write tables x and y. Affinities of irreversible reactions.
    !           Affinities of individual reactants, kcal.
    if (nrct .le. 0) then
        go to 300
    end if

    nlim2 = min(9,nrct)

    if (qtitl) then
        do i = 1,3
            write (ntabx,1000) ux
        end do

        write (ntabx,2040) ux
2040 format(a1,20x,'affinities of irreversible reactions')

        write (ntabx,1000) ux
        write (ntabx,2005) ux,(ureac(nrc)(1:8), nrc = 1,nlim2)
        write (ntabx,2010) ux,(ureac(nrc)(9:16), nrc = 1,nlim2)
        write (ntabx,2010) ux,(ureac(nrc)(17:24), nrc = 1,nlim2)
        write (ntabx,1000) ux
    end if

    write (ntabx,1040) ux,xilog,tdays,tlogd,(afrc1(nrc), nrc = 1,nlim2)

    if (nrct .le. nlim2) then
        go to 300
    end if

    nlim1 = nlim2 + 1
    nlim2 = min(18,nrct)

    if (qtitl) then
        do i = 1,3
            write (ntabx,1000) uy
        end do

        write (ntabx,2005) uy,(ureac(nrc)(1:8), nrc = nlim1,nlim2)
        write (ntabx,2010) uy,(ureac(nrc)(9:16), nrc = nlim1,nlim2)
        write (ntabx,2010) uy,(ureac(nrc)(17:24), nrc = nlim1,nlim2)
        write (ntabx,1000) uy
    end if

    write (ntabx,1040) uy,xilog,tdays,tlogd,(afrc1(nrc), nrc = nlim1,nlim2)

300 continue
end subroutine wrtabx