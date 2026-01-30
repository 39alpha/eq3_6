subroutine wrtabc(acflg,actlg,actw,afrc1,aft1,alk,conclg,cteaq,ctb,dvoso,dwoso,eh,fje,fo2lg,fugac,fxi,iktmax,iopt,jflag,jsflag,kmax,kstep,kx1,kxt,mrmlra,modr,mosp,mospt,moph,mopht,narn1,narn2,nbasp,nbaspd,nbt,nbtmax,ncmpr,nct,nctmax,nelect,ngrn1,ngrn2,ngtmax,nhydr,nhydx,nllnmx,no2gaq,noptmx,noutpt,npt,nptmax,nrct,nrctmx,nstmax,ntabx,ntidmx,ntitl2,ntitld,ntitmx,nttyo,nxrn1,nxrn2,nxtmax,pe,ph,phmes,ppmwb,ppmwe,prcinf,press,prminf,qrho,qriinf,rho,rhowc,sidrph,sigmam,tdsgks,tdsglw,tempc,time1,uelem,ulinex,uphase,uplatm,ureac,uspec,usteq6,utitl2,utitld,uveeq6,vodrt,vosoct,wkgh2o,wodrt,wosoct,xbar,xbarlg,xi1)
    use iso_fortran_env, only: dp => real64
    !! This subroutine writes to TABX (the scrambled TAB file) using
    !! a csv (comma separated value) format. A .csv file can be opened
    !! by spreadsheets and other software that support data analysis
    !! and plotting.
    !! This subroutine is called by:
    !!   EQ6/path.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: iktmax
    integer :: kmax
    integer :: nbtmax
    integer :: nctmax
    integer :: ngtmax
    integer :: nllnmx
    integer :: noptmx
    integer :: nptmax
    integer :: nrctmx
    integer :: nstmax
    integer :: ntidmx
    integer :: ntitmx
    integer :: nxtmax

    integer :: noutpt
    integer :: ntabx
    integer :: nttyo

    integer :: iopt(noptmx)
    integer :: jflag(nstmax)
    integer :: jsflag(nstmax)
    integer :: nbasp(nbtmax)
    integer :: nbaspd(nbtmax)
    integer :: ncmpr(2,nptmax)

    integer :: kstep
    integer :: kx1
    integer :: kxt
    integer :: narn1
    integer :: narn2
    integer :: nbt
    integer :: nct
    integer :: nelect
    integer :: ngrn1
    integer :: ngrn2
    integer :: nhydr
    integer :: nhydx
    integer :: no2gaq
    integer :: npt
    integer :: nrct
    integer :: ntitl2
    integer :: ntitld
    integer :: nxrn1
    integer :: nxrn2

    logical :: qrho
    logical :: qriinf

    character(len=nllnmx) :: ulinex
    character(len=80) :: utitl2(ntitmx)
    character(len=80) :: utitld(ntidmx)
    character(len=48) :: uspec(nstmax)
    character(len=24) :: uphase(nptmax)
    character(len=24) :: ureac(nrctmx)
    character(len=8) :: uelem(nctmax)
    character(len=8) :: uplatm
    character(len=8) :: usteq6
    character(len=8) :: uveeq6

    real(kind=8) :: acflg(nstmax)
    real(kind=8) :: actlg(nstmax)
    real(kind=8) :: afrc1(nrctmx)
    real(kind=8) :: conclg(nstmax)
    real(kind=8) :: ctb(nbtmax)
    real(kind=8) :: cteaq(nctmax)
    real(kind=8) :: fugac(ngtmax)
    real(kind=8) :: modr(nrctmx)
    real(kind=8) :: moph(nptmax)
    real(kind=8) :: mopht(nptmax)
    real(kind=8) :: mosp(nstmax)
    real(kind=8) :: mospt(nstmax)
    real(kind=8) :: ppmwb(nbtmax)
    real(kind=8) :: ppmwe(nctmax)
    real(kind=8) :: sidrph(nptmax)
    real(kind=8) :: xbar(nstmax)
    real(kind=8) :: xbarlg(nstmax)

    real(kind=8) :: actw
    real(kind=8) :: aft1
    real(kind=8) :: alk
    real(kind=8) :: dvoso
    real(kind=8) :: dwoso
    real(kind=8) :: eh
    real(kind=8) :: fje
    real(kind=8) :: fo2lg
    real(kind=8) :: fxi
    real(kind=8) :: mrmlra
    real(kind=8) :: pe
    real(kind=8) :: ph
    real(kind=8) :: phmes
    real(kind=8) :: prcinf
    real(kind=8) :: press
    real(kind=8) :: prminf
    real(kind=8) :: rho
    real(kind=8) :: rhowc
    real(kind=8) :: sigmam
    real(kind=8) :: tdsgks
    real(kind=8) :: tdsglw
    real(kind=8) :: tempc
    real(kind=8) :: time1
    real(kind=8) :: vodrt
    real(kind=8) :: vosoct
    real(kind=8) :: wkgh2o
    real(kind=8) :: wodrt
    real(kind=8) :: wosoct
    real(kind=8) :: xi1

    ! Local variable declarations with variable global dimensioning.
    integer :: ngfmax
    integer :: npmmax
    integer :: npsmax
    integer :: nxcmax

    SAVE ngfmax,npmmax,npsmax,nxcmax

    integer, dimension(:), allocatable :: ngasfx

    SAVE ngasfx

    character(len=24), dimension(:,:), allocatable :: ussphx
    character(len=24), dimension(:), allocatable :: ubaspx
    character(len=24), dimension(:), allocatable :: ubasqx
    character(len=24), dimension(:), allocatable :: usidrx
    character(len=24), dimension(:), allocatable :: uphasx
    character(len=8), dimension(:), allocatable :: uelacx

    SAVE ubaspx,ubasqx,uelacx,uphasx,usidrx,ussphx

    real(kind=8), dimension(:), allocatable :: acflbx
    real(kind=8), dimension(:), allocatable :: actlbx
    real(kind=8), dimension(:), allocatable :: afrcx
    real(kind=8), dimension(:), allocatable :: conlbx
    real(kind=8), dimension(:), allocatable :: ctbaqx
    real(kind=8), dimension(:), allocatable :: cteaqx
    real(kind=8), dimension(:), allocatable :: molrbx
    real(kind=8), dimension(:), allocatable :: molrex
    real(kind=8), dimension(:), allocatable :: mophx
    real(kind=8), dimension(:), allocatable :: ppmvbx
    real(kind=8), dimension(:), allocatable :: ppmvex
    real(kind=8), dimension(:), allocatable :: ppmwbx
    real(kind=8), dimension(:), allocatable :: ppmwex
    real(kind=8), dimension(:), allocatable :: sidrx
    real(kind=8), dimension(:), allocatable :: xfrac

    SAVE acflbx,actlbx,afrcx,conlbx,ctbaqx,cteaqx,molrbx,molrex,mophx,ppmvbx,ppmvex,ppmwbx,ppmwex,sidrx,xfrac

    ! Local variable declarations required to be saved between calls.
    integer :: ngft
    integer :: npmt
    integer :: npst
    integer :: nxct
    integer :: nbtpr
    integer :: nctpr
    integer :: nqtpr

    SAVE ngft,npmt,npst,nxct
    SAVE nbtpr,nctpr,nqtpr

    integer :: ilstmx

    SAVE ilstmx

    logical :: qgftag
    logical :: qpstag

    SAVE qgftag,qpstag

    ! Local variable declarations.
    integer :: i
    integer :: jj
    integer :: j1
    integer :: j2
    integer :: j3
    integer :: j4
    integer :: jc
    integer :: jlead
    integer :: jline
    integer :: jq
    integer :: jsum
    integer :: jt1
    integer :: jt2
    integer :: jt3
    integer :: jthv1
    integer :: jts
    integer :: n
    integer :: nb
    integer :: nc
    integer :: ng
    integer :: nlim
    integer :: np
    integer :: nr1
    integer :: nr2
    integer :: ns
    integer :: ns1
    integer :: ns2

    integer :: ilnobl

    logical :: qnewph
    logical :: qnewss
    logical :: qpsext
    logical :: qwrthd

    character(len=48) :: ux48
    character(len=16) :: ux16
    character(len=8) :: ux8a
    character(len=8) :: ux8b
    character(len=8) :: ultag1
    character(len=8) :: ultag2
    character(len=8) :: ultag3
    character(len=8) :: ultags
    character(len=8) :: ultgv1
    character(len=8) :: unotav

    real(kind=8) :: sidrcv
    real(kind=8) :: tdays

    data unotav/'N/A'/

    data sidrcv/-10./

    ! Allocate or reallocate local work arrays as needed.
    if (.not.ALLOCATED(uphasx)) then
        ! Local work arrays are not allocated. Note that only one array
        ! is tested to see if all are allocated or not.
        ilstmx = (nllnmx/25) - 3

        npsmax = min(ilstmx,nptmax)

        ngfmax = min(ilstmx,ngtmax)

        nxcmax = nxtmax*iktmax
        nxcmax = min(ilstmx,nxcmax)

        npmmax = 8*kmax
        npmmax = min(ilstmx,npmmax)

        npmt = 0
        npst = 0
        ngft = 0
        nxct = 0

        qgftag = .false.
        qpstag = .false.

        nctpr = 0
        nqtpr = 0
        nbtpr = 0

        ALLOCATE(uphasx(npmmax))
        ALLOCATE(mophx(npmmax))

        ALLOCATE(uelacx(nctmax))
        ALLOCATE(cteaqx(nctmax),ppmwex(nctmax))

        ALLOCATE(ubasqx(nbtmax))
        ALLOCATE(ctbaqx(nbtmax),ppmwbx(nbtmax))

        if (qrho) then
            ALLOCATE(molrex(nctmax),ppmvex(nctmax))
            ALLOCATE(molrbx(nbtmax),ppmvbx(nbtmax))
        end if

        ALLOCATE(ubaspx(nbtmax))
        ALLOCATE(acflbx(nbtmax),actlbx(nbtmax),conlbx(nbtmax))

        ALLOCATE(afrcx(nrctmx))

        ALLOCATE(usidrx(npsmax))
        ALLOCATE(sidrx(npsmax))

        ALLOCATE(ngasfx(ngfmax))

        ALLOCATE(ussphx(2,nxcmax))
        ALLOCATE(xfrac(nxcmax))

        ! Zero the contents of the local work arrays.
        do n = 1,nctmax
            uelacx(n) = ' '
            cteaqx(n) = 0.
            ppmwex(n) = 0.
        end do

        do n = 1,nbtmax
            ubasqx(n) = ' '
            ctbaqx(n) = 0.
            ppmwbx(n) = 0.
        end do

        if (qrho) then
            do n = 1,nctmax
                molrex(n) = 0.
                ppmvex(n) = 0.
            end do

            do n = 1,nbtmax
                molrbx(n) = 0.
                ppmvbx(n) = 0.
            end do
        end if

        do n = 1,nbtmax
            ubaspx(n) = ' '
            acflbx(n) = 0.
            actlbx(n) = 0.
            conlbx(n) = 0.
        end do

        do n = 1,nrctmx
            afrcx(n) = 0.
        end do

        do n = 1,npmmax
            uphasx(n) = ' '
            mophx(n) = 0.
        end do

        do n = 1,npsmax
            usidrx(n) = ' '
            sidrx(n) = 0.
        end do

        do n = 1,ngfmax
            ngasfx(n) = 0
        end do

        do n = 1,nxcmax
            ussphx(1,n) = ' '
            ussphx(2,n) = ' '
            xfrac(n) = 0.
        end do

        ! Construct the local element list.
        nctpr = 0

        do nc = 1,nct
            if (uelem(nc)(1:2).ne.'H ' .and. uelem(nc)(1:2) .ne. 'O ')      then
                nctpr = nctpr + 1
                uelacx(nctpr) = uelem(nc)
            end if
        end do

        ! Construct the local solute basis species list.
        nqtpr = 0

        do nb = 1,nbt
            ns1 = nbaspd(nb)
            ns2 = nbasp(nb)

            if (ns1.ge.narn1 .and. ns1.le.narn2) then
                if (jsflag(ns2).lt.2 .and. jflag(ns2).ne.30) then
                    if (ns1.ne.narn1 .and. ns1.ne.nelect .and.          ns1.ne.no2gaq) then
                        nqtpr = nqtpr + 1
                        ubasqx(nqtpr) = uspec(ns1)(1:24)
                    end if
                end if
            end if
        end do

        ! Construct the local complete basis species list.
        nbtpr = 0

        do nb = 1,nbt
            ns1 = nbaspd(nb)
            ns2 = nbasp(nb)

            if (ns1.ge.narn1 .and. ns1.le.narn2) then
                if (jsflag(ns2).lt.2 .and. jflag(ns2).ne.30) then
                    nbtpr = nbtpr + 1
                    ubaspx(nbtpr) = uspec(ns1)(1:24)
                end if
            end if
        end do

        ! Construct the local gas species index map list.
        ngft = 0

        do ns = ngrn1,ngrn2
            if (jsflag(ns) .lt. 2) then
                if (ngft .ge. ngfmax) then
                    ! Prepare to write a "More omitted" label
                    qgftag = .true.
                    go to 100
                end if

                ngft = ngft + 1
                ngasfx(ngft) = ns
            end if
        end do

100 continue
    end if

    ! Write output tables of some selected parameters. These are
    ! identified by short strings, commonly letters or letters and
    ! numbers. They are written on the TABX file in scrambled order
    ! (i.e., the lines of any given tables are interspersed with those
    ! of the other tables). Each line is marked at the beginning with
    ! an id string marking the table to which it belongs. The lines
    ! will later be descrambled so that all lines belonging to a given
    ! table are contiguous.
    ! Write Table A:
    !   Run title (first 15 lines max)
    !   Code version identification
    !   Data file title (first 10 lines max)
    ! Note: the original primary title (utitl1) has already been rolled
    ! over to the current secondary title (utitl2).
    qwrthd = .false.

    if (kstep.eq.0 .and. iopt(18).eq.0) then
        ultag1 = 'A       '
        jt1 = len_trim(ultag1)
        ultags = ultag1
        jts = jt1

        ! Write the table header.
        write (ntabx,1000) ultag1(1:jt1)
1000 format(a,',Table A,Title information (filtered),')

        qwrthd = .true.

        write (ntabx,'(a,",")') ultag1(1:jt1)
        j2 = ilnobl(uveeq6)
        j3 = ilnobl(usteq6)
        j4 = ilnobl(uplatm)
        write (ntabx,1020) ultag1(1:jt1),uveeq6(1:j2),usteq6(1:j3),uplatm(1:j4)
1020 format(a,',Running EQ6 from EQ3/6-V',a,' (',a,') for ',a,',')

        write (ntabx,'(a,",---,")') ultag1(1:jt1)
        write (ntabx,'(a,",")') ultag1(1:jt1)

        ! Write the table data.
        nlim = min(ntitl2,15)

        do n = 1,nlim
            jj = min(80,nllnmx)
            ulinex = utitl2(n)(1:jj)
            j1 = len_trim(ulinex)

            ! Check for double quotes and commas.
            jq = index(ulinex,'"')
            jc = index(ulinex,',')
            jsum = jq + jc

            if (jsum .gt. 0) then
                ! Replace any double quotes in the title text by
                ! single quotes.
                if (jq .gt. 0) then
110 continue
                    ulinex(jq:jq) = "'"
                    jq = index(ulinex,'"')

                    if (jq .gt. 0) then
                        go to 110
                    end if
                end if

                ! Replace any commas in the title text by semicolons.
                jc = index(ulinex,',')

                if (jc .gt. 0) then
120 continue
                    ulinex(jc:jc) = ";"
                    jc = index(ulinex,',')

                    if (jc .gt. 0) then
                        go to 120
                    end if
                end if
            end if

            write (ntabx,'(a,",",a,",")') ultag1(1:jt1),ulinex(1:j1)
        end do

        write (ntabx,'(a,",---,")') ultag1(1:jt1)
        write (ntabx,'(a,",")') ultag1(1:jt1)

        nlim = min(ntitld,10)

        do n = 1,nlim
            jj = min(80,nllnmx)
            ulinex = utitld(n)(1:jj)
            j1 = len_trim(ulinex)

            ! Check for double quotes and commas.
            jq = index(ulinex,'"')
            jc = index(ulinex,',')
            jsum = jq + jc

            if (jsum .gt. 0) then
                ! Replace any double quotes in the title text by
                ! single quotes.
                if (jq .gt. 0) then
130 continue
                    ulinex(jq:jq) = "'"
                    jq = index(ulinex,'"')

                    if (jq .gt. 0) then
                        go to 130
                    end if
                end if

                ! Replace any commas in the title text by semicolons.
                jc = index(ulinex,',')

                if (jc .gt. 0) then
140 continue
                    ulinex(jc:jc) = ";"
                    jc = index(ulinex,',')

                    if (jc .gt. 0) then
                        go to 140
                    end if
                end if
            end if

            write (ntabx,'(a,",",a,",")') ultag1(1:jt1),ulinex(1:j1)
        end do

        write (ntabx,'(a,",---,")') ultag1(1:jt1)
    end if

    ! Write Table B1:
    !   Xi
    !   t(days)
    !   Temperature(C)
    !   Pressure(bars)
    !   pH
    !   pmH
    !   log fO2
    !   Eh(v)
    !   pe
    !   Activity of water
    ! Prepare the data.
    if (.not.qriinf) then
        tdays = time1/86400.
    end if

    if (eh .le. -99999.) then
        eh = -99999.
    end if

    if (pe .le. -99999.) then
        pe = -99999.
    end if

    if (aft1 .ge. 99999.) then
        aft1 = 99999.
    end if

    ultag1 = 'B1      '
    jt1 = len_trim(ultag1)

    if (qwrthd) then
        ! Write an end marker for the previous table.
        write (ntabx,'(a,",EndTable:,",a,",")')  ultag1(1:jt1),ultags(1:jts)
        ultags = ultag1
        jts = jt1

        ! Write two empty lines before the table header.
        do i = 1,2
            write (ntabx,'(a,",")') ultag1(1:jt1)
        end do

        ! Write the table header.
        write (ntabx,1040) ultag1(1:jt1)
1040 format(a,',Table B1,Miscellaneous parameters I,')

        ! Write the column header.
        ulinex = 'Xi,t(days),Temp(C),Press(bars),pH,pmH,'
        j1 = len_trim(ulinex)
        ux48 = 'log fO2,Eh(v),pe,aw'
        j2 = len_trim(ux48)
        jsum = j1 + j2
        ulinex(j1 + 1:jsum) = ux48(1:j2)
        j1 = jsum

        write (ntabx,'(a,",",a,",")') ultag1(1:jt1),ulinex(1:j1)
    end if

    ! Write the table data.
    ulinex = ''
    j1 = 0

    write (ux16,'(1pg11.4)') xi1
    call lejust(ux16)
    j2 = len_trim(ux16)
    jsum = j1 + j2
    ulinex(j1 + 1:jsum) = ux16(1:j2)
    jsum = jsum + 1
    ulinex(jsum:jsum) = ','
    j1 = jsum

    write (ux16,'(1pg11.4)') tdays
    call lejust(ux16)
    j2 = len_trim(ux16)
    jsum = j1 + j2
    ulinex(j1 + 1:jsum) = ux16(1:j2)
    jsum = jsum + 1
    ulinex(jsum:jsum) = ','
    j1 = jsum

    write (ux16,'(f11.4)') tempc
    call lejust(ux16)
    j2 = len_trim(ux16)
    jsum = j1 + j2
    ulinex(j1 + 1:jsum) = ux16(1:j2)
    jsum = jsum + 1
    ulinex(jsum:jsum) = ','
    j1 = jsum

    write (ux16,'(f11.4)') press
    call lejust(ux16)
    j2 = len_trim(ux16)
    jsum = j1 + j2
    ulinex(j1 + 1:jsum) = ux16(1:j2)
    jsum = jsum + 1
    ulinex(jsum:jsum) = ','
    j1 = jsum

    write (ux16,'(f11.4)') ph
    call lejust(ux16)
    j2 = len_trim(ux16)
    jsum = j1 + j2
    ulinex(j1 + 1:jsum) = ux16(1:j2)
    jsum = jsum + 1
    ulinex(jsum:jsum) = ','
    j1 = jsum

    write (ux16,'(f11.4)') phmes
    call lejust(ux16)
    j2 = len_trim(ux16)
    jsum = j1 + j2
    ulinex(j1 + 1:jsum) = ux16(1:j2)
    jsum = jsum + 1
    ulinex(jsum:jsum) = ','
    j1 = jsum

    write (ux16,'(f11.4)') fo2lg
    call lejust(ux16)
    j2 = len_trim(ux16)
    jsum = j1 + j2
    ulinex(j1 + 1:jsum) = ux16(1:j2)
    jsum = jsum + 1
    ulinex(jsum:jsum) = ','
    j1 = jsum

    write (ux16,'(f11.4)') eh
    call lejust(ux16)
    j2 = len_trim(ux16)
    jsum = j1 + j2
    ulinex(j1 + 1:jsum) = ux16(1:j2)
    jsum = jsum + 1
    ulinex(jsum:jsum) = ','
    j1 = jsum

    write (ux16,'(f11.4)') pe
    call lejust(ux16)
    j2 = len_trim(ux16)
    jsum = j1 + j2
    ulinex(j1 + 1:jsum) = ux16(1:j2)
    jsum = jsum + 1
    ulinex(jsum:jsum) = ','
    j1 = jsum

    write (ux16,'(f11.4)') actw
    call lejust(ux16)
    j2 = len_trim(ux16)
    jsum = j1 + j2
    ulinex(j1 + 1:jsum) = ux16(1:j2)
    j1 = jsum

    write (ntabx,'(a,",",a,",")') ultag1(1:jt1),ulinex(1:j1)

    ! Write Table B2:
    !   Xi
    !   t(days)
    !   Mass of solvent H2O, kg
    !   Alkalinity(eq/kg H2O)
    !   Sigma m
    !   Ionic strength
    !   Ionic asymmetry
    !   TDS, g/kg.sol
    !   TDS, g/L
    !   Aqueous solution density, g/L
    !   Molarity/molality ratio
    ! Prepare the data. None to prepare.
    ultag1 = 'B2      '
    jt1 = len_trim(ultag1)

    if (qwrthd) then
        ! Write an end marker for the previous table.
        write (ntabx,'(a,",EndTable:,",a,",")')  ultag1(1:jt1),ultags(1:jts)
        ultags = ultag1
        jts = jt1

        ! Write two empty lines before the table header.
        do i = 1,2
            write (ntabx,'(a,",")') ultag1(1:jt1)
        end do

        ! Write the table header.
        write (ntabx,1060) ultag1(1:jt1)
1060 format(a,',Table B2,Miscellaneous parameters II,')

        ! Write the column header.
        ulinex = 'Xi,t(days),H2O(kg),Alk(eq/kg.H2O),Signam(m),I(m),'
        j1 = len_trim(ulinex)
        ux48 = 'J(m),TDS(g/kg.sol),TDS(g/L),density(g/L),M/m'
        j2 = len_trim(ux48)
        jsum = j1 + j2
        ulinex(j1 + 1:jsum) = ux48(1:j2)
        j1 = jsum

        write (ntabx,'(a,",",a,",")') ultag1(1:jt1),ulinex(1:j1)
    end if

    ! Write the table data.
    ulinex = ''
    j1 = 0

    write (ux16,'(1pg11.4)') xi1
    call lejust(ux16)
    j2 = len_trim(ux16)
    jsum = j1 + j2
    ulinex(j1 + 1:jsum) = ux16(1:j2)
    jsum = jsum + 1
    ulinex(jsum:jsum) = ','
    j1 = jsum

    write (ux16,'(1pg11.4)') tdays
    call lejust(ux16)
    j2 = len_trim(ux16)
    jsum = j1 + j2
    ulinex(j1 + 1:jsum) = ux16(1:j2)
    jsum = jsum + 1
    ulinex(jsum:jsum) = ','
    j1 = jsum

    write (ux16,'(1pg11.4)') wkgh2o
    call lejust(ux16)
    j2 = len_trim(ux16)
    jsum = j1 + j2
    ulinex(j1 + 1:jsum) = ux16(1:j2)
    jsum = jsum + 1
    ulinex(jsum:jsum) = ','
    j1 = jsum

    write (ux16,'(1pg11.4)') alk
    call lejust(ux16)
    j2 = len_trim(ux16)
    jsum = j1 + j2
    ulinex(j1 + 1:jsum) = ux16(1:j2)
    jsum = jsum + 1
    ulinex(jsum:jsum) = ','
    j1 = jsum

    write (ux16,'(1pg11.4)') sigmam
    call lejust(ux16)
    j2 = len_trim(ux16)
    jsum = j1 + j2
    ulinex(j1 + 1:jsum) = ux16(1:j2)
    jsum = jsum + 1
    ulinex(jsum:jsum) = ','
    j1 = jsum

    write (ux16,'(1pg11.4)') fxi
    call lejust(ux16)
    j2 = len_trim(ux16)
    jsum = j1 + j2
    ulinex(j1 + 1:jsum) = ux16(1:j2)
    jsum = jsum + 1
    ulinex(jsum:jsum) = ','
    j1 = jsum

    write (ux16,'(1pg11.4)') fje
    call lejust(ux16)
    j2 = len_trim(ux16)
    jsum = j1 + j2
    ulinex(j1 + 1:jsum) = ux16(1:j2)
    jsum = jsum + 1
    ulinex(jsum:jsum) = ','
    j1 = jsum

    write (ux16,'(f11.4)') tdsgks
    call lejust(ux16)
    j2 = len_trim(ux16)
    jsum = j1 + j2
    ulinex(j1 + 1:jsum) = ux16(1:j2)
    jsum = jsum + 1
    ulinex(jsum:jsum) = ','
    j1 = jsum

    write (ux16,'(f11.4)') tdsglw
    call lejust(ux16)
    j2 = len_trim(ux16)
    jsum = j1 + j2
    ulinex(j1 + 1:jsum) = ux16(1:j2)
    jsum = jsum + 1
    ulinex(jsum:jsum) = ','
    j1 = jsum

    write (ux16,'(f11.4)') rhowc
    call lejust(ux16)
    j2 = len_trim(ux16)
    jsum = j1 + j2
    ulinex(j1 + 1:jsum) = ux16(1:j2)
    jsum = jsum + 1
    ulinex(jsum:jsum) = ','
    j1 = jsum

    write (ux16,'(f11.4)') mrmlra
    call lejust(ux16)
    j2 = len_trim(ux16)
    jsum = j1 + j2
    ulinex(j1 + 1:jsum) = ux16(1:j2)
    j1 = jsum

    write (ntabx,'(a,",",a,",")') ultag1(1:jt1),ulinex(1:j1)

    ! Write Table C1:
    !   Xi
    !   t(days)
    !   Molalities of dissolved elements
    ! and Table C2:
    !   Xi
    !   t(days)
    !   ppm (mg/kg.sol) of dissolved elements
    ! Prepare the data.
    nctpr = 0

    do nc = 1,nct
        if (uelem(nc)(1:2).ne.'H ' .and. uelem(nc)(1:2) .ne. 'O ') then
            nctpr = nctpr + 1
            uelacx(nctpr) = uelem(nc)
            cteaqx(nctpr) = cteaq(nc)
            ppmwex(nctpr) = ppmwe(nc)
        end if
    end do

    ultag1 = 'C1      '
    ultag2 = 'C2      '
    jt1 = len_trim(ultag1)
    jt2 = len_trim(ultag2)

    if (qwrthd) then
        ! Write an end marker for the table previous to C1.
        write (ntabx,'(a,",EndTable:,",a,",")')  ultag1(1:jt1),ultags(1:jts)
        ultags = ultag1
        jts = jt1

        ! Write an end marker for the table previous to C2.
        write (ntabx,'(a,",EndTable:,",a,",")')  ultag2(1:jt2),ultags(1:jts)
        ultags = ultag2
        jts = jt2

        ! Write two empty lines before the table header.
        do i = 1,2
            write (ntabx,'(a,",")') ultag1(1:jt1)
            write (ntabx,'(a,",")') ultag2(1:jt2)
        end do

        ! Write the table header.
        write (ntabx,1100) ultag1(1:jt1)
1100 format(a,',Table C1,Dissolved elements(molality),')

        write (ntabx,1110) ultag2(1:jt2)
1110 format(a,',Table C2,Dissolved elements(ppm: mg/kg.sol),')

        ! Write the column headers.
        ulinex = 'Xi,t(days)'
        j1 = len_trim(ulinex)

        do n = 1,nctpr
            j1 = j1 + 1
            ulinex(j1:j1) = ","
            j2 = len_trim(uelacx(n))
            jsum = j1 + j2
            ulinex(j1 + 1:jsum) = uelacx(n)(1:j2)
            j1 = jsum
        end do

        jline = j1
        write (ntabx,'(a,",",a,",")') ultag1(1:jt1),ulinex(1:j1)
        write (ntabx,'(a,",",a,",")') ultag2(1:jt2),ulinex(1:j1)
    end if

    ! Write the table data. First, the common data (Xi, t( days)).
    ulinex = ''
    j1 = 0

    write (ux16,'(1pg11.4)') xi1
    call lejust(ux16)
    j2 = len_trim(ux16)
    jsum = j1 + j2
    ulinex(j1 + 1:jsum) = ux16(1:j2)
    jsum = jsum + 1
    ulinex(jsum:jsum) = ','
    j1 = jsum

    write (ux16,'(1pg11.4)') tdays
    call lejust(ux16)
    j2 = len_trim(ux16)
    jsum = j1 + j2
    ulinex(j1 + 1:jsum) = ux16(1:j2)
    j1 = jsum
    jlead = j1

    ! Write the data specific to Table C1.
    do n = 1,nctpr
        j1 = j1 + 1
        ulinex(j1:j1) = ","
        write (ux16,'(1pg11.4)') cteaqx(n)
        call lejust(ux16)
        j2 = len_trim(ux16)
        jsum = j1 + j2
        ulinex(j1 + 1:jsum) = ux16(1:j2)
        j1 = jsum
    end do

    write (ntabx,'(a,",",a,",")') ultag1(1:jt1),ulinex(1:j1)

    ! Write the data specific to Table C2.
    ulinex(jlead + 1:j1) = ' '
    j1 = jlead

    do n = 1,nctpr
        j1 = j1 + 1
        ulinex(j1:j1) = ","
        write (ux16,'(1pg11.4)') ppmwex(n)
        call lejust(ux16)
        j2 = len_trim(ux16)
        jsum = j1 + j2
        ulinex(j1 + 1:jsum) = ux16(1:j2)
        j1 = jsum
    end do

    write (ntabx,'(a,",",a,",")') ultag2(1:jt2),ulinex(1:j1)

    if (qrho) then
        ! Write Table C3:
        !   Xi
        !   t(days)
        !   Molarities of dissolved elements
        ! and Table C4:
        !   Xi
        !   t(days)
        !   ppm (mg/L) of dissolved elements
        ! Prepare the data.
        nctpr = 0

        do nc = 1,nct
            if (uelem(nc)(1:2).ne.'H ' .and. uelem(nc)(1:2) .ne. 'O ')      then
                nctpr = nctpr + 1
                uelacx(nctpr) = uelem(nc)
                molrex(nctpr) = cteaq(nc)*mrmlra
                ppmvex(nctpr) = ppmwe(nc)*rho
            end if
        end do

        ultag1 = 'C3       '(1:8)
        ultag2 = 'C4       '(1:8)
        jt1 = len_trim(ultag1)
        jt2 = len_trim(ultag2)

        if (qwrthd) then
            ! Write an end marker for the table previous to C3.
            write (ntabx,'(a,",EndTable:,",a,",")')    ultag1(1:jt1),ultags(1:jts)
            ultags = ultag1
            jts = jt1

            ! Write an end marker for the table previous to C4.
            write (ntabx,'(a,",EndTable:,",a,",")')    ultag2(1:jt2),ultags(1:jts)
            ultags = ultag2
            jts = jt2

            ! Write two empty lines before the table header.
            do i = 1,2
                write (ntabx,'(a,",")') ultag1(1:jt1)
                write (ntabx,'(a,",")') ultag2(1:jt2)
            end do

            ! Write the table header.
            write (ntabx,1120) ultag1(1:jt1)
1120 format(a,',Table C3,Dissolved elements(molarity),')

            write (ntabx,1130) ultag2(1:jt2)
1130 format(a,',Table C4,Dissolved elements(ppm: mg/L),')

            ! Write the column headers.
            ulinex = 'Xi,t(days)'
            j1 = len_trim(ulinex)

            do n = 1,nctpr
                j1 = j1 + 1
                ulinex(j1:j1) = ","
                j2 = len_trim(uelacx(n))
                jsum = j1 + j2
                ulinex(j1 + 1:jsum) = uelacx(n)(1:j2)
                j1 = jsum
            end do

            jline = j1
            write (ntabx,'(a,",",a,",")') ultag1(1:jt1),ulinex(1:j1)
            write (ntabx,'(a,",",a,",")') ultag2(1:jt2),ulinex(1:j1)
        end if

        ! Write the table data. First, the common data (Xi, t( days)).
        ulinex = ''
        j1 = 0

        write (ux16,'(1pg11.4)') xi1
        call lejust(ux16)
        j2 = len_trim(ux16)
        jsum = j1 + j2
        ulinex(j1 + 1:jsum) = ux16(1:j2)
        jsum = jsum + 1
        ulinex(jsum:jsum) = ','
        j1 = jsum

        write (ux16,'(1pg11.4)') tdays
        call lejust(ux16)
        j2 = len_trim(ux16)
        jsum = j1 + j2
        ulinex(j1 + 1:jsum) = ux16(1:j2)
        j1 = jsum
        jlead = j1

        ! Write the data specific to Table C3.
        do n = 1,nctpr
            j1 = j1 + 1
            ulinex(j1:j1) = ","
            write (ux16,'(1pg11.4)') molrex(n)
            call lejust(ux16)
            j2 = len_trim(ux16)
            jsum = j1 + j2
            ulinex(j1 + 1:jsum) = ux16(1:j2)
            j1 = jsum
        end do

        write (ntabx,'(a,",",a,",")') ultag1(1:jt1),ulinex(1:j1)

        ! Write the data specific to Table C4.
        ulinex(jlead + 1:j1) = ' '
        j1 = jlead

        do n = 1,nctpr
            j1 = j1 + 1
            ulinex(j1:j1) = ","
            write (ux16,'(1pg11.4)') ppmvex(n)
            call lejust(ux16)
            j2 = len_trim(ux16)
            jsum = j1 + j2
            ulinex(j1 + 1:jsum) = ux16(1:j2)
            j1 = jsum
        end do

        write (ntabx,'(a,",",a,",")') ultag2(1:jt2),ulinex(1:j1)
    end if

    ! Write Table D1:
    !   Xi
    !   t(days)
    !   Total molalities of solute basis species
    ! and Table D2:
    !   Xi
    !   t(days)
    !   Total ppm (mg/kg.sol) of solute basis species
    ! Prepare the data.
    n = 0

    do nb = 1,nbt
        ns1 = nbaspd(nb)
        ns2 = nbasp(nb)

        if (ns1.ge.narn1 .and. ns1.le.narn2) then
            if (jsflag(ns2).lt.2 .and. jflag(ns2).ne.30) then
                if (ns1.ne.narn1 .and. ns1.ne.nelect .and.        ns1.ne.no2gaq) then
                    n = n + 1
                    ctbaqx(n) = ctb(nb)
                    ppmwbx(n) = ppmwb(nb)
                end if
            end if
        end if
    end do

    if (n .ne. nqtpr) then
        ! Programming error trap.
        write (ux8a,'(i8)') n
        write (ux8b,'(i8)') nqtpr
        call lejust(ux8a)
        call lejust(ux8b)
        j1 = len_trim(ux8a)
        j2 = len_trim(ux8b)
        write (noutpt,1140) ux8a(1:j1),ux8b(1:j2)
        write (nttyo,1140) ux8a(1:j1),ux8b(1:j2)
1140 format(/' * Error - (EQ6/wrtabc) Programming error trap: ','The solute',/7x,'basis species counts for Tables D1 and D2',' do not match',/7x,'(',a,' versus ',a,').')

        stop
    end if

    ultag1 = 'D1      '
    ultag2 = 'D2      '
    jt1 = len_trim(ultag1)
    jt2 = len_trim(ultag2)

    if (qwrthd) then
        ! Write an end marker for the table previous to D1.
        write (ntabx,'(a,",EndTable:,",a,",")')  ultag1(1:jt1),ultags(1:jts)
        ultags = ultag1
        jts = jt1

        ! Write an end marker for the table previous to D2.
        write (ntabx,'(a,",EndTable:,",a,",")')  ultag2(1:jt2),ultags(1:jts)
        ultags = ultag2
        jts = jt2

        ! Write two empty lines before the table header.
        do i = 1,2
            write (ntabx,'(a,",")') ultag1(1:jt1)
            write (ntabx,'(a,",")') ultag2(1:jt2)
        end do

        ! Write the table header.
        write (ntabx,1220) ultag1(1:jt1)
1220 format(a,',Table D1,Solute basis species(total molality),')

        write (ntabx,1230) ultag2(1:jt2)
1230 format(a,',Table D2,Solute basis species','(total ppm: mg/kg.sol),')

        ! Write the column headers.
        ulinex = 'Xi,t(days)'
        j1 = len_trim(ulinex)

        do n = 1,nqtpr
            j1 = j1 + 1
            ulinex(j1:j1) = ","
            j2 = len_trim(ubasqx(n))
            jsum = j1 + j2
            ulinex(j1 + 1:jsum) = ubasqx(n)(1:j2)
            j1 = jsum
        end do

        write (ntabx,'(a,",",a,",")') ultag1(1:jt1),ulinex(1:j1)
        write (ntabx,'(a,",",a,",")') ultag2(1:jt2),ulinex(1:j1)
    end if

    ! Write the table data. First, the common data (Xi, t( days)).
    ulinex = ''
    j1 = 0

    write (ux16,'(1pg11.4)') xi1
    call lejust(ux16)
    j2 = len_trim(ux16)
    jsum = j1 + j2
    ulinex(j1 + 1:jsum) = ux16(1:j2)
    jsum = jsum + 1
    ulinex(jsum:jsum) = ','
    j1 = jsum

    write (ux16,'(1pg11.4)') tdays
    call lejust(ux16)
    j2 = len_trim(ux16)
    jsum = j1 + j2
    ulinex(j1 + 1:jsum) = ux16(1:j2)
    j1 = jsum
    jlead = j1

    ! Write the data specific to Table D1.
    do n = 1,nqtpr
        j1 = j1 + 1
        ulinex(j1:j1) = ","
        write (ux16,'(1pg11.4)') ctbaqx(n)
        call lejust(ux16)
        j2 = len_trim(ux16)
        jsum = j1 + j2
        ulinex(j1 + 1:jsum) = ux16(1:j2)
        j1 = jsum
    end do

    write (ntabx,'(a,",",a,",")') ultag1(1:jt1),ulinex(1:j1)

    ! Write the data specific to Table D2.
    ulinex(jlead + 1:j1) = ' '
    j1 = jlead

    do n = 1,nqtpr
        j1 = j1 + 1
        ulinex(j1:j1) = ","
        write (ux16,'(f11.4)') ppmwbx(n)
        call lejust(ux16)
        j2 = len_trim(ux16)
        jsum = j1 + j2
        ulinex(j1 + 1:jsum) = ux16(1:j2)
        j1 = jsum
    end do

    write (ntabx,'(a,",",a,",")') ultag2(1:jt2),ulinex(1:j1)

    if (qrho) then
        ! Write Table D3:
        !   Xi
        !   t(days)
        !   Total molarities of solute basis species
        ! and Table D4:
        !   Xi
        !   t(days)
        !   Total ppm (mg/L) solute file basis species
        ! Prepare the data.
        n = 0

        do nb = 1,nbt
            ns1 = nbaspd(nb)
            ns2 = nbasp(nb)

            if (ns1.ge.narn1 .and. ns1.le.narn2) then
                if (jsflag(ns2).lt.2 .and. jflag(ns2).ne.30) then
                    if (ns1.ne.narn1 .and. ns1.ne.nelect .and.          ns1.ne.no2gaq) then
                        n = n + 1
                        molrbx(n) = ctb(nb)*mrmlra
                        ppmvbx(n) = ppmwb(nb)*rho
                    end if
                end if
            end if
        end do

        if (n .ne. nqtpr) then
            ! Programming error trap.
            write (ux8a,'(i8)') n
            write (ux8b,'(i8)') nqtpr
            call lejust(ux8a)
            call lejust(ux8b)
            j1 = len_trim(ux8a)
            j2 = len_trim(ux8b)
            write (noutpt,1235) ux8a(1:j1),ux8b(1:j2)
            write (nttyo,1235) ux8a(1:j1),ux8b(1:j2)
1235 format(/' * Error - (EQ6/wrtabc) Programming error trap: ','The solute',/7x,'basis species counts for Tables D3 and D4',' do not match',/7x,'(',a,' versus ',a,').')

            stop
        end if

        ultag1 = 'D3      '
        ultag2 = 'D4      '
        jt1 = len_trim(ultag1)
        jt2 = len_trim(ultag2)

        if (qwrthd) then
            ! Write an end marker for the table previous to D3.
            write (ntabx,'(a,",EndTable:,",a,",")')    ultag1(1:jt1),ultags(1:jts)
            ultags = ultag1
            jts = jt1

            ! Write an end marker for the table previous to D4.
            write (ntabx,'(a,",EndTable:,",a,",")')    ultag2(1:jt2),ultags(1:jts)
            ultags = ultag2
            jts = jt2

            ! Write two empty lines before the table header.
            do i = 1,2
                write (ntabx,'(a,",")') ultag1(1:jt1)
                write (ntabx,'(a,",")') ultag2(1:jt2)
            end do

            ! Write the table header.
            write (ntabx,1240) ultag1(1:jt1)
1240 format(a,',Table D3,Solute basis species(total molarity),')

            write (ntabx,1250) ultag2(1:jt2)
1250 format(a,',Table D4,Solute basis species(total ppm: mg/L),')

            ! Write the column headers.
            ulinex = 'Xi,t(days)'
            j1 = len_trim(ulinex)

            do n = 1,nqtpr
                j1 = j1 + 1
                ulinex(j1:j1) = ","
                j2 = len_trim(ubasqx(n))
                jsum = j1 + j2
                ulinex(j1 + 1:jsum) = ubasqx(n)(1:j2)
                j1 = jsum
            end do

            write (ntabx,'(a,",",a,",")') ultag1(1:jt1),ulinex(1:j1)
            write (ntabx,'(a,",",a,",")') ultag2(1:jt2),ulinex(1:j1)
        end if

        ! Write the table data. First, the common data (Xi, t( days)).
        ulinex = ''
        j1 = 0

        write (ux16,'(1pg11.4)') xi1
        call lejust(ux16)
        j2 = len_trim(ux16)
        jsum = j1 + j2
        ulinex(j1 + 1:jsum) = ux16(1:j2)
        jsum = jsum + 1
        ulinex(jsum:jsum) = ','
        j1 = jsum

        write (ux16,'(1pg11.4)') tdays
        call lejust(ux16)
        j2 = len_trim(ux16)
        jsum = j1 + j2
        ulinex(j1 + 1:jsum) = ux16(1:j2)
        j1 = jsum
        jlead = j1

        ! Write the data specific to Table D3.
        do n = 1,nqtpr
            j1 = j1 + 1
            ulinex(j1:j1) = ","
            write (ux16,'(1pg11.4)') molrbx(n)
            call lejust(ux16)
            j2 = len_trim(ux16)
            jsum = j1 + j2
            ulinex(j1 + 1:jsum) = ux16(1:j2)
            j1 = jsum
        end do

        write (ntabx,'(a,",",a,",")') ultag1(1:jt1),ulinex(1:j1)

        ! Write the data specific to Table D4.
        ulinex(jlead + 1:j1) = ' '
        j1 = jlead

        do n = 1,nqtpr
            j1 = j1 + 1
            ulinex(j1:j1) = ","
            write (ux16,'(f11.4)') ppmvbx(n)
            call lejust(ux16)
            j2 = len_trim(ux16)
            jsum = j1 + j2
            ulinex(j1 + 1:jsum) = ux16(1:j2)
            j1 = jsum
        end do

        write (ntabx,'(a,",",a,",")') ultag2(1:jt2),ulinex(1:j1)
    end if

    ! Write Table E1:
    !   Xi
    !   t(days)
    !   Log true molalities of basis species (log mole fraction for H2O,
    !   (log fugacity for aqueous O2(g))
    ! Table E2:
    !   Xi
    !   t(days)
    !   Log activities of basis species (log fugacity for aqueous O2(g))
    ! and Table E3:
    !   Xi
    !   t(days)
    !   Log true activity coefficients of basis species
    ! Prepare the data.
    n = 0

    do nb = 1,nbt
        ns1 = nbaspd(nb)
        ns2 = nbasp(nb)

        if (ns1.ge.narn1 .and. ns1.le.narn2) then
            if (jsflag(ns2).lt.2 .and. jflag(ns2).ne.30) then
                n = n + 1

                if (ns1 .eq. narn1) then
                    conlbx(n) = xbarlg(ns1)
                    actlbx(n) = actlg(ns1)
                    acflbx(n) = acflg(ns1)
                else if (ns1 .eq. no2gaq) then
                    conlbx(n) = fo2lg
                    actlbx(n) = fo2lg
                    acflbx(n) = acflg(ns1)
                else
                    conlbx(n) = conclg(ns1)
                    actlbx(n) = actlg(ns1)
                    acflbx(n) = acflg(ns1)
                end if
            end if
        end if
    end do

    if (n .ne. nbtpr) then
        ! Programming error trap.
        write (ux8a,'(i8)') n
        write (ux8b,'(i8)') nbtpr
        call lejust(ux8a)
        call lejust(ux8b)
        j1 = len_trim(ux8a)
        j2 = len_trim(ux8b)
        write (noutpt,1260) ux8a(1:j1),ux8b(1:j2)
        write (nttyo,1260) ux8a(1:j1),ux8b(1:j2)
1260 format(/' * Error - (EQ6/wrtabc) Programming error trap: ','The local complete',/7x,'basis species counts for Tables E1,',' E2, and E3 do not match',/7x,'(',a,' versus ',a,').')

        stop
    end if

    ultag1 = 'E1      '
    ultag2 = 'E2      '
    ultag3 = 'E3      '
    jt1 = len_trim(ultag1)
    jt2 = len_trim(ultag2)
    jt3 = len_trim(ultag3)

    if (qwrthd) then
        ! Write an end marker for the table previous to E1.
        write (ntabx,'(a,",EndTable:,",a,",")')  ultag1(1:jt1),ultags(1:jts)
        ultags = ultag1
        jts = jt1

        ! Write an end marker for the table previous to E2.
        write (ntabx,'(a,",EndTable:,",a,",")')  ultag2(1:jt2),ultags(1:jts)
        ultags = ultag2
        jts = jt2

        ! Write an end marker for the table previous to E3.
        write (ntabx,'(a,",EndTable:,",a,",")')  ultag3(1:jt3),ultags(1:jts)
        ultags = ultag3
        jts = jt3

        ! Write two empty lines before the table header.
        do i = 1,2
            write (ntabx,'(a,",")') ultag1(1:jt1)
            write (ntabx,'(a,",")') ultag2(1:jt2)
            write (ntabx,'(a,",")') ultag3(1:jt3)
        end do

        ! Write the table header.
        write (ntabx,1270) ultag1(1:jt1)
1270 format(a,',Table E1,Basis species(log true molality;',' log mole fraction for H2O, log fugacity for O2(g)),')

        write (ntabx,1280) ultag2(1:jt2)
1280 format(a,',Table E2,Basis species(log activity;',' log fugacity for O2(g)),')

        write (ntabx,1290) ultag3(1:jt3)
1290 format(a,',Table E3,Basis species(log activity coefficient),')

        ! Write the column headers.
        ulinex = 'Xi,t(days)'
        j1 = len_trim(ulinex)

        do n = 1,nbtpr
            j1 = j1 + 1
            ulinex(j1:j1) = ","
            j2 = len_trim(ubaspx(n))
            jsum = j1 + j2
            ulinex(j1 + 1:jsum) = ubaspx(n)(1:j2)
            j1 = jsum
        end do

        write (ntabx,'(a,",",a,",")') ultag1(1:jt1),ulinex(1:j1)
        write (ntabx,'(a,",",a,",")') ultag2(1:jt2),ulinex(1:j1)
        write (ntabx,'(a,",",a,",")') ultag3(1:jt3),ulinex(1:j1)
    end if

    ! Write the table data. First, the common data (Xi, t( days)).
    ulinex = ''
    j1 = 0

    write (ux16,'(1pg11.4)') xi1
    call lejust(ux16)
    j2 = len_trim(ux16)
    jsum = j1 + j2
    ulinex(j1 + 1:jsum) = ux16(1:j2)
    jsum = jsum + 1
    ulinex(jsum:jsum) = ','
    j1 = jsum

    write (ux16,'(1pg11.4)') tdays
    call lejust(ux16)
    j2 = len_trim(ux16)
    jsum = j1 + j2
    ulinex(j1 + 1:jsum) = ux16(1:j2)
    j1 = jsum
    jlead = j1

    ! Write the data specific to Table E1.
    do n = 1,nbtpr
        j1 = j1 + 1
        ulinex(j1:j1) = ","
        write (ux16,'(1pg11.4)') conlbx(n)
        call lejust(ux16)
        j2 = len_trim(ux16)
        jsum = j1 + j2
        ulinex(j1 + 1:jsum) = ux16(1:j2)
        j1 = jsum
    end do

    write (ntabx,'(a,",",a,",")') ultag1(1:jt1),ulinex(1:j1)

    ! Write the data specific to Table E2.
    ulinex(jlead + 1:j1) = ' '
    j1 = jlead

    do n = 1,nbtpr
        j1 = j1 + 1
        ulinex(j1:j1) = ","
        write (ux16,'(f11.4)') actlbx(n)
        call lejust(ux16)
        j2 = len_trim(ux16)
        jsum = j1 + j2
        ulinex(j1 + 1:jsum) = ux16(1:j2)
        j1 = jsum
    end do

    write (ntabx,'(a,",",a,",")') ultag2(1:jt2),ulinex(1:j1)

    ! Write the data specific to Table E3.
    ulinex(jlead + 1:j1) = ' '
    j1 = jlead

    do n = 1,nbtpr
        j1 = j1 + 1
        ulinex(j1:j1) = ","
        write (ux16,'(f11.4)') acflbx(n)
        call lejust(ux16)
        j2 = len_trim(ux16)
        jsum = j1 + j2
        ulinex(j1 + 1:jsum) = ux16(1:j2)
        j1 = jsum
    end do

    write (ntabx,'(a,",",a,",")') ultag3(1:jt3),ulinex(1:j1)

    ! Write Table J:
    !   Xi
    !   t(days)
    !   Moles of reactants destroyed/created
    ! Prepare the data. None to prepare.
    ultag1 = 'J       '
    jt1 = len_trim(ultag1)

    if (qwrthd) then
        ! Write an end marker for the previous table.
        write (ntabx,'(a,",EndTable:,",a,",")')  ultag1(1:jt1),ultags(1:jts)
        ultags = ultag1
        jts = jt1

        ! Write two empty lines before the table header.
        do i = 1,2
            write (ntabx,'(a,",")') ultag1(1:jt1)
        end do

        ! Write the table header.
        write (ntabx,1340) ultag1(1:jt1)
1340 format(a,',Table J,Moles of reactants destroyed/created,')

        ! Write the column header.
        ulinex = 'Xi,t(days)'
        j1 = len_trim(ulinex)

        do n = 1,nrct
            j1 = j1 + 1
            ulinex(j1:j1) = ","
            j2 = len_trim(ureac(n))
            jsum = j1 + j2
            ulinex(j1 + 1:jsum) = ureac(n)(1:j2)
            j1 = jsum
        end do

        write (ntabx,'(a,",",a,",")') ultag1(1:jt1),ulinex(1:j1)
    end if

    ! Write the table data.
    ulinex = ''
    j1 = 0

    write (ux16,'(1pg11.4)') xi1
    call lejust(ux16)
    j2 = len_trim(ux16)
    jsum = j1 + j2
    ulinex(j1 + 1:jsum) = ux16(1:j2)
    jsum = jsum + 1
    ulinex(jsum:jsum) = ','
    j1 = jsum

    write (ux16,'(1pg11.4)') tdays
    call lejust(ux16)
    j2 = len_trim(ux16)
    jsum = j1 + j2
    ulinex(j1 + 1:jsum) = ux16(1:j2)
    j1 = jsum

    do n = 1,nrct
        j1 = j1 + 1
        ulinex(j1:j1) = ","
        write (ux16,'(1pg11.4)') modr(n)
        call lejust(ux16)
        j2 = len_trim(ux16)
        jsum = j1 + j2
        ulinex(j1 + 1:jsum) = ux16(1:j2)
        j1 = jsum
    end do

    write (ntabx,'(a,",",a,",")') ultag1(1:jt1),ulinex(1:j1)

    ! Write Table K:
    !   Xi
    !   t(days)
    !   Total Affinity, kcal
    !   Affinities of reactants, kcal
    ! Prepare the data.
    do n = 1,nrct
        afrcx(n) = min(99999.0_dp,afrc1(n))
    end do

    ultag1 = 'K       '
    jt1 = len_trim(ultag1)

    if (qwrthd) then
        ! Write an end marker for the previous table.
        write (ntabx,'(a,",EndTable:,",a,",")')  ultag1(1:jt1),ultags(1:jts)
        ultags = ultag1
        jts = jt1

        ! Write two empty lines before the table header.
        do i = 1,2
            write (ntabx,'(a,",")') ultag1(1:jt1)
        end do

        ! Write the table header.
        write (ntabx,1360) ultag1(1:jt1)
1360 format(a,',Table K,Affinities of reactants (kcal),')

        ! Write the column header.
        ulinex = 'Xi,t(days)'
        j1 = len_trim(ulinex)
        ux48 = ',Total'
        j2 = len_trim(ux48)
        jsum = j1 + j2
        ulinex(j1 + 1:jsum) = ux48(1:j2)
        j1 = jsum

        do n = 1,nrct
            j1 = j1 + 1
            ulinex(j1:j1) = ","
            j2 = len_trim(ureac(n))
            jsum = j1 + j2
            ulinex(j1 + 1:jsum) = ureac(n)(1:j2)
            j1 = jsum
        end do

        write (ntabx,'(a,",",a,",")') ultag1(1:jt1),ulinex(1:j1)
    end if

    ! Write the table data.
    ulinex = ''
    j1 = 0

    write (ux16,'(1pg11.4)') xi1
    call lejust(ux16)
    j2 = len_trim(ux16)
    jsum = j1 + j2
    ulinex(j1 + 1:jsum) = ux16(1:j2)
    jsum = jsum + 1
    ulinex(jsum:jsum) = ','
    j1 = jsum

    write (ux16,'(1pg11.4)') tdays
    call lejust(ux16)
    j2 = len_trim(ux16)
    jsum = j1 + j2
    ulinex(j1 + 1:jsum) = ux16(1:j2)
    jsum = jsum + 1
    ulinex(jsum:jsum) = ','
    j1 = jsum

    write (ux16,'(f11.4)') aft1
    call lejust(ux16)
    j2 = len_trim(ux16)
    jsum = j1 + j2
    ulinex(j1 + 1:jsum) = ux16(1:j2)
    j1 = jsum

    do n = 1,nrct
        j1 = j1 + 1
        ulinex(j1:j1) = ","
        write (ux16,'(f11.4)') afrcx(n)
        call lejust(ux16)
        j2 = len_trim(ux16)
        jsum = j1 + j2
        ulinex(j1 + 1:jsum) = ux16(1:j2)
        j1 = jsum
    end do

    write (ntabx,'(a,",",a,",")') ultag1(1:jt1),ulinex(1:j1)

    ! Write Table P:
    !   Xi
    !   t(days)
    !   Moles of product phases
    ! Note: in these tables, moles of product phase refers to what
    ! is present in the Equilibrium System for Titration and Closed
    ! system models, and to cumulative totals for the Fluid-Centered
    ! Flow-Through system.
    ! Prepare the data.
    do n = 1,npmt
        mophx(n) = 0.
    end do

    qnewph = .false.

    if (iopt(1) .ne. 2) then
        ! Closed or titration system.
        do np = 2,npt
            ! **       if (moph(np).gt.0. .and. uphase(np)(1:5).ne."fix_f") then
            if (moph(np) .gt. 0.) then
                do n = 1,npmt
                    if (uphase(np) .eq. uphasx(n)) then
                        ! Have found a product mineral in the existing list.
                        mophx(n) = moph(np)
                        go to 310
                    end if
                end do

                ! Did not find a current product mineral in the existing
                ! list. Add it to the list.
                qnewph = .true.
                npmt = npmt + 1
                n = npmt
                uphasx(n) = uphase(np)
                mophx(n) = moph(np)
            end if

310 continue
        end do
    else
        ! Fluid-centered flow-through system.
        do np = 2,npt
            ! **       if (moph(np).gt.0. .and. uphase(np)(1:5).ne."fix_f") then
            if (moph(np) .gt. 0.) then
                do n = 1,npmt
                    if (uphase(np) .eq. uphasx(n)) then
                        ! Have found a product mineral in the existing list.
                        mophx(n) = mopht(np)
                        go to 320
                    end if
                end do

                ! Did not find a current product mineral in the existing
                ! list. Add it to the list.
                qnewph = .true.
                npmt = npmt + 1
                n = npmt
                uphasx(n) = uphase(np)
                mophx(n) = mopht(np)
            end if

320 continue
        end do
    end if

    ultag1 = 'P       '
    jt1 = len_trim(ultag1)

    if (qwrthd) then
        ! Write an end marker for the previous table.
        write (ntabx,'(a,",EndTable:,",a,",")')  ultag1(1:jt1),ultags(1:jts)
        ultags = ultag1
        jts = jt1

        ! Write two empty lines before the table header.
        do i = 1,2
            write (ntabx,'(a,",")') ultag1(1:jt1)
        end do

        ! Write the table header.
        if (iopt(1) .ne. 2) then
            write (ntabx,1370) ultag1(1:jt1)
1370 format(a,',Table P,Moles of product minerals,')
        else
            write (ntabx,1380) ultag1(1:jt1)
1380 format(a,',Table P,Moles of product minerals (cumulative),')
        end if
    end if

    if (qnewph .or. kstep.le.0) then
        ! Write the column header. This is an example of a variable
        ! header.
        ulinex = 'Xi,t(days)'
        j1 = len_trim(ulinex)

        do n = 1,npmt
            j1 = j1 + 1
            ulinex(j1:j1) = ","
            j2 = len_trim(uphasx(n))
            jsum = j1 + j2
            ulinex(j1 + 1:jsum) = uphasx(n)(1:j2)
            j1 = jsum
        end do

        ultgv1(1:6) = ultag1(1:6)
        ultgv1(7:8) = 'vh'
        jthv1 = len_trim(ultgv1)
        write (ntabx,'(a,",",a,",")') ultgv1(1:jthv1),ulinex(1:j1)
    end if

    ! Write the table data.
    ulinex = ''
    j1 = 0

    write (ux16,'(1pg11.4)') xi1
    call lejust(ux16)
    j2 = len_trim(ux16)
    jsum = j1 + j2
    ulinex(j1 + 1:jsum) = ux16(1:j2)
    jsum = jsum + 1
    ulinex(jsum:jsum) = ','
    j1 = jsum

    write (ux16,'(1pg11.4)') tdays
    call lejust(ux16)
    j2 = len_trim(ux16)
    jsum = j1 + j2
    ulinex(j1 + 1:jsum) = ux16(1:j2)
    j1 = jsum

    do n = 1,npmt
        j1 = j1 + 1
        ulinex(j1:j1) = ","
        write (ux16,'(1pg11.4)') mophx(n)
        call lejust(ux16)
        j2 = len_trim(ux16)
        jsum = j1 + j2
        ulinex(j1 + 1:jsum) = ux16(1:j2)
        j1 = jsum
    end do

    write (ntabx,'(a,",",a,",")') ultag1(1:jt1),ulinex(1:j1)

    ! Write Table Q:
    !   Xi
    !   t(days)
    !   Saturation indices of potential product phases
    ! Note: in these tables, the saturation index (SI) is represented
    ! by the function log(Q/K). Values are contained in the sidrph
    ! array. A potential product phase is considered to be one which
    ! at some point in the run has a saturation index above a certain
    ! threshhold (sidrcv).
    ! Prepare the data.
    do n = 1,npst
        sidrx(n) = 0.
    end do

    qpsext = .false.

    do np = 2,npt
        if (uphase(np)(1:5).ne."fix_f") then
            do n = 1,npst
                if (uphase(np) .eq. usidrx(n)) then
                    ! Have found a potential product phase in the existing list.
                    sidrx(n) = sidrph(np)
                    go to 350
                end if
            end do

            if (sidrph(np).gt.sidrcv) then
                ! Did not find a potential product phase in the existing
                ! list. Add a phase, write a "More omitted" tag, or skip if
                ! a tag has been previously written.
                if (qpstag) then
                    ! The list has previously maxed out and a "More omitted"
                    ! tag has been written. Do nothing.
                    continue
                else if (npst .eq. npsmax) then
                    ! The list has previously maxed out, but a tag has
                    ! not been written. Prepare to write a tag.
                    qpsext = .true.
                    qpstag = .true.
                else
                    ! Add a new potential product phase to the list.
                    qpsext = .true.
                    npst = npst + 1
                    n = npst
                    usidrx(n) = uphase(np)
                    sidrx(n) = sidrph(np)
                end if
            end if
        end if

350 continue
    end do

    ultag1 = 'Q       '
    jt1 = len_trim(ultag1)

    if (qwrthd) then
        ! Write an end marker for the previous table.
        write (ntabx,'(a,",EndTable:,",a,",")')  ultag1(1:jt1),ultags(1:jts)
        ultags = ultag1
        jts = jt1

        ! Write two empty lines before the table header.
        do i = 1,2
            write (ntabx,'(a,",")') ultag1(1:jt1)
        end do

        ! Write the table header.
        write (ntabx,1392) ultag1(1:jt1)
1392 format(a,',Table Q,Saturation indices of potential',' product phases,')
    end if

    if (qpsext .or. kstep.le.0) then
        ! Write the column header. This is an example of a variable
        ! header.
        ulinex = 'Xi,t(days)'
        j1 = len_trim(ulinex)

        do n = 1,npst
            j1 = j1 + 1
            ulinex(j1:j1) = ","
            j2 = len_trim(usidrx(n))
            jsum = j1 + j2
            ulinex(j1 + 1:jsum) = usidrx(n)(1:j2)
            j1 = jsum
        end do

        if (qpstag) then
            ! Have more phases than can fit.
            j1 = j1 + 1
            ulinex(j1:j1) = ","
            ux16 = 'More omitted'
            j2 = len_trim(ux16)
            jsum = j1 + j2
            ulinex(j1 + 1:jsum) = ux16(1:j2)
            j1 = jsum
        end if

        ultgv1(1:6) = ultag1(1:6)
        ultgv1(7:8) = 'vh'
        jthv1 = len_trim(ultgv1)
        write (ntabx,'(a,",",a,",")') ultgv1(1:jthv1),ulinex(1:j1)
    end if

    ! Write the table data.
    ulinex = ''
    j1 = 0

    write (ux16,'(1pg11.4)') xi1
    call lejust(ux16)
    j2 = len_trim(ux16)
    jsum = j1 + j2
    ulinex(j1 + 1:jsum) = ux16(1:j2)
    jsum = jsum + 1
    ulinex(jsum:jsum) = ','
    j1 = jsum

    write (ux16,'(1pg11.4)') tdays
    call lejust(ux16)
    j2 = len_trim(ux16)
    jsum = j1 + j2
    ulinex(j1 + 1:jsum) = ux16(1:j2)
    j1 = jsum

    do n = 1,npst
        j1 = j1 + 1
        ulinex(j1:j1) = ","
        write (ux16,'(f11.4)') sidrx(n)
        call lejust(ux16)
        j2 = len_trim(ux16)
        jsum = j1 + j2
        ulinex(j1 + 1:jsum) = ux16(1:j2)
        j1 = jsum
    end do

    write (ntabx,'(a,",",a,",")') ultag1(1:jt1),ulinex(1:j1)

    ! Write Table T:
    !   Xi
    !   t(days)
    !   Fugacities of gases, bars
    ! Prepare the data. Nothing to prepare.
    ultag1 = 'T       '
    jt1 = len_trim(ultag1)

    if (qwrthd) then
        ! Write an end marker for the previous table.
        write (ntabx,'(a,",EndTable:,",a,",")')  ultag1(1:jt1),ultags(1:jts)
        ultags = ultag1
        jts = jt1

        ! Write two empty lines before the table header.
        do i = 1,2
            write (ntabx,'(a,",")') ultag1(1:jt1)
        end do

        ! Write the table header.
        write (ntabx,1400) ultag1(1:jt1)
1400 format(a,',Table T,Fugacities (bars),')

        ! Write the column header.
        ulinex = 'Xi,t(days)'
        j1 = len_trim(ulinex)

        do n = 1,ngft
            j1 = j1 + 1
            ulinex(j1:j1) = ","
            ns = ngasfx(n)
            j2 = len_trim(uspec(ns)(1:24))
            jsum = j1 + j2
            ulinex(j1 + 1:jsum) = uspec(ns)(1:j2)
            j1 = jsum
        end do

        if (qgftag) then
            ! Have more gas species than can fit.
            j1 = j1 + 1
            ulinex(j1:j1) = ","
            ux16 = 'More omitted'
            j2 = len_trim(ux16)
            jsum = j1 + j2
            ulinex(j1 + 1:jsum) = ux16(1:j2)
            j1 = jsum
        end if

        write (ntabx,'(a,",",a,",")') ultag1(1:jt1),ulinex(1:j1)
    end if

    ! Write the table data.
    ulinex = ''
    j1 = 0

    write (ux16,'(1pg11.4)') xi1
    call lejust(ux16)
    j2 = len_trim(ux16)
    jsum = j1 + j2
    ulinex(j1 + 1:jsum) = ux16(1:j2)
    jsum = jsum + 1
    ulinex(jsum:jsum) = ','
    j1 = jsum

    write (ux16,'(1pg11.4)') tdays
    call lejust(ux16)
    j2 = len_trim(ux16)
    jsum = j1 + j2
    ulinex(j1 + 1:jsum) = ux16(1:j2)
    j1 = jsum

    do n = 1,ngft
        j1 = j1 + 1
        ulinex(j1:j1) = ","
        ns = ngasfx(n)
        ng = ns - ngrn1 + 1
        write (ux16,'(1pg11.4)') fugac(ng)
        call lejust(ux16)
        j2 = len_trim(ux16)
        jsum = j1 + j2
        ulinex(j1 + 1:jsum) = ux16(1:j2)
        j1 = jsum
    end do

    write (ntabx,'(a,",",a,",")') ultag1(1:jt1),ulinex(1:j1)

    ! Write Table W:
    !   Xi
    !   t(days)
    !   Total mass of solids destroyed, grams
    !   Total mass of solids created, grams
    !   Net mass of solids created, grams
    !   Total volume of solids destroyed, cc
    !   Total volume of solids created, cc
    !   Net volume of solids created, cc
    ultag1 = 'W       '
    jt1 = len_trim(ultag1)

    if (qwrthd) then
        ! Write an end marker for the previous table.
        write (ntabx,'(a,",EndTable:,",a,",")')  ultag1(1:jt1),ultags(1:jts)
        ultags = ultag1
        jts = jt1

        ! Write two empty lines before the table header.
        do i = 1,2
            write (ntabx,'(a,",")') ultag1(1:jt1)
        end do

        ! Write the table header.
        write (ntabx,1410) ultag1(1:jt1)
1410 format(a,',Table W,Overall mass and volume changes,')

        ! Write the column header.
        ulinex = 'Xi,t(days),g destroyed,g created,g net,'
        j1 = len_trim(ulinex)
        ux48 = 'cc destroyed,cc created,cc net'
        j2 = len_trim(ux48)
        jsum = j1 + j2
        ulinex(j1 + 1:jsum) = ux48(1:j2)
        j1 = jsum

        write (ntabx,'(a,",",a,",")') ultag1(1:jt1),ulinex(1:j1)
    end if

    ! Write the table data.
    ulinex = ''
    j1 = 0

    write (ux16,'(1pg11.4)') xi1
    call lejust(ux16)
    j2 = len_trim(ux16)
    jsum = j1 + j2
    ulinex(j1 + 1:jsum) = ux16(1:j2)
    jsum = jsum + 1
    ulinex(jsum:jsum) = ','
    j1 = jsum

    write (ux16,'(1pg11.4)') tdays
    j2 = len_trim(ux16)
    jsum = j1 + j2
    ulinex(j1 + 1:jsum) = ux16(1:j2)
    jsum = jsum + 1
    ulinex(jsum:jsum) = ','
    j1 = jsum

    write (ux16,'(1pg11.4)') wodrt
    call lejust(ux16)
    j2 = len_trim(ux16)
    jsum = j1 + j2
    ulinex(j1 + 1:jsum) = ux16(1:j2)
    jsum = jsum + 1
    ulinex(jsum:jsum) = ','
    j1 = jsum

    write (ux16,'(1pg11.4)') wosoct
    call lejust(ux16)
    j2 = len_trim(ux16)
    jsum = j1 + j2
    ulinex(j1 + 1:jsum) = ux16(1:j2)
    jsum = jsum + 1
    ulinex(jsum:jsum) = ','
    j1 = jsum

    write (ux16,'(1pg11.4)') dwoso
    call lejust(ux16)
    j2 = len_trim(ux16)
    jsum = j1 + j2
    ulinex(j1 + 1:jsum) = ux16(1:j2)
    jsum = jsum + 1
    ulinex(jsum:jsum) = ','
    j1 = jsum

    write (ux16,'(1pg11.4)') vodrt
    call lejust(ux16)
    j2 = len_trim(ux16)
    jsum = j1 + j2
    ulinex(j1 + 1:jsum) = ux16(1:j2)
    jsum = jsum + 1
    ulinex(jsum:jsum) = ','
    j1 = jsum

    write (ux16,'(1pg11.4)') vosoct
    call lejust(ux16)
    j2 = len_trim(ux16)
    jsum = j1 + j2
    ulinex(j1 + 1:jsum) = ux16(1:j2)
    jsum = jsum + 1
    ulinex(jsum:jsum) = ','
    j1 = jsum

    write (ux16,'(1pg11.4)') dvoso
    call lejust(ux16)
    j2 = len_trim(ux16)
    jsum = j1 + j2
    ulinex(j1 + 1:jsum) = ux16(1:j2)
    j1 = jsum

    write (ntabx,'(a,",",a,",")') ultag1(1:jt1),ulinex(1:j1)

    if (iopt(4) .ge. 1) then
        ! Write Table X:
        !   Xi
        !   t(days)
        !   Mole fractions of solid solutions.
        ! Note: cumulative averages are presented for the case of the
        ! Fluid-Centered Flow-Through system.
        ! Prepare the data.
        do n = 1,nxct
            xfrac(n) = 0.
        end do

        qnewss = .false.

        if (iopt(1) .ne. 2) then
            ! Closed or titration system.
            do np = nxrn1,nxrn2
                if (moph(np) .gt. 0.) then
                    do n = 1,nxct
                        if (uphase(np) .eq. ussphx(1,n)) then
                            ! Have found a solid solution in the existing list.
                            nr1 = ncmpr(1,np)
                            nr2 = ncmpr(2,np)

                            do ns = nr1,nr2
                                if (uspec(ns) .eq. ussphx(2,n)) then
                                    xfrac(n) = xbar(ns)
                                    go to 410
                                end if
                            end do
                        end if
                    end do

                    ! Did not find a current solid solution in the existing
                    ! list. Add it to the list.
                    qnewss = .true.
                    nr1 = ncmpr(1,np)
                    nr2 = ncmpr(2,np)

                    do ns = nr1,nr2
                        nxct = nxct + 1
                        n = nxct
                        ussphx(1,n) = uphase(np)
                        ussphx(2,n) = uspec(ns)(1:24)
                        xfrac(n) = xbar(ns)
                    end do
                end if

410 continue
            end do
        else
            ! Fluid-centered flow-through system.
            do np = nxrn1,nxrn2
                if (mopht(np) .gt. 0.) then
                    do n = 1,nxct
                        if (uphase(np) .eq. ussphx(1,n)) then
                            ! Have found a solid solution in the existing list.
                            nr1 = ncmpr(1,np)
                            nr2 = ncmpr(2,np)

                            do ns = nr1,nr2
                                if (uspec(ns) .eq. ussphx(2,n)) then
                                    xfrac(n) = mospt(ns)/mopht(np)
                                    go to 420
                                end if
                            end do
                        end if
                    end do

                    ! Did not find a current solid solution in the existing
                    ! list. Add it to the list.
                    qnewss = .true.
                    nr1 = ncmpr(1,np)
                    nr2 = ncmpr(2,np)

                    do ns = nr1,nr2
                        nxct = nxct + 1
                        n = nxct
                        ussphx(1,n) = uphase(np)
                        ussphx(2,n) = uspec(ns)(1:24)
                        xfrac(n) = mospt(ns)/mopht(np)
                    end do
                end if

420 continue
            end do
        end if

        ultag1 = 'X       '
        jt1 = len_trim(ultag1)

        if (qwrthd) then
            ! Write an end marker for the previous table.
            write (ntabx,'(a,",EndTable:,",a,",")')    ultag1(1:jt1),ultags(1:jts)
            ultags = ultag1
            jts = jt1

            ! Write two empty lines before the table header.
            do i = 1,2
                write (ntabx,'(a,",")') ultag1(1:jt1)
            end do

            ! Write the table header.
            write (ntabx,1440) ultag1(1:jt1)
1440 format(a,',Table X,Solid solution mole fractions,')
        end if

        if (qnewss .or. kstep.le.0) then
            ! Write the column header.
            ulinex = 'Xi,t(days)'
            j1 = len_trim(ulinex)

            do n = 1,nxct
                j1 = j1 + 1
                ulinex(j1:j1) = ","
                j2 = len_trim(ussphx(1,n))
                jsum = j1 + j2
                ulinex(j1 + 1:jsum) = ussphx(1,n)(1:j2)
                j1 = jsum
            end do

            write (ntabx,'(a,",",a,",")') ultag1(1:jt1),ulinex(1:j1)
        end if

        ! Write the table data.
        ulinex = ''
        j1 = 0

        write (ux16,'(1pg11.4)') xi1
        call lejust(ux16)
        j2 = len_trim(ux16)
        jsum = j1 + j2
        ulinex(j1 + 1:jsum) = ux16(1:j2)
        jsum = jsum + 1
        ulinex(jsum:jsum) = ','
        j1 = jsum

        write (ux16,'(1pg11.4)') tdays
        call lejust(ux16)
        j2 = len_trim(ux16)
        jsum = j1 + j2
        ulinex(j1 + 1:jsum) = ux16(1:j2)
        j1 = jsum

        do n = 1,nxct
            j1 = j1 + 1
            ulinex(j1:j1) = ","
            write (ux16,'(1pg11.4)') xfrac(n)
            call lejust(ux16)
            j2 = len_trim(ux16)
            jsum = j1 + j2
            ulinex(j1 + 1:jsum) = ux16(1:j2)
            j1 = jsum
        end do

        write (ntabx,'(a,",",a,",")') ultag1(1:jt1),ulinex(1:j1)
    end if

    ! Write Table Z:
    !   End of data marker
    if (qwrthd) then
        ! Write an end marker for the previous table.
        write (ntabx,'("Z,EndTable:,",a,",")') ultags(1:jts)

        write (ntabx,'("Z,")')
        write (ntabx,'("Z,")')

        ! Write the table header.
        write (ntabx,'("Z,Table Z,Endfile,")')

        write (ntabx,'("Z,Endtable:,Z,")')
    end if
end subroutine wrtabc
