subroutine rd6w7(cdac,cesrb,cplim,csigma,dlzmx1,dlzmx2,dlzidp,dzpllg,dzplot,dzprlg,dzprnt,eact,electr,fk,iact,ifile,iktbt,iktmax,imchmx,imech,iodb,iopg,iopr,iopt,ioscan,itermx,jcode,jreac,jtemp,jxmod,kct,kdim,kmax,kmt,kprs,ksq,ksplmx,ksppmx,kstpmx,kxmod,kxt,hact,modr,moffg,morr,mprs,mteaqb,mteb,nctmax,ndact,ndctmx,nesrbt,nffg,nffgmx,ninpts,nmodl1,nmodl2,nodbmx,nopgmx,noprmx,noptmx,nordlm,npslmx,nprmn,nprmx,nprsmx,nrct,nrctmx,nrk,nsk,nsrtmx,nsscmx,nsslmx,ntitl1,ntitl2,ntitmx,ntrymx,nttkmx,nttyo,nxmdmx,nxmod,nxopex,nxopmx,nxopt,nxpemx,nxrtmx,qend,qrderr,rk0,rxbarb,sscrew,sk,tempci,tempc0,timemx,tolbt,toldl,tolsat,tolsst,tolx,trk0,tstrt,ttk,udac,uelemb,uendb,uesrb,uffg,undms,unrms,uprs,ureac,utitl1,utitl2,uxct16,uxmd24,uxopex,uxopt,vreac,xlkffg,xlkmod,zimax,zistrt,zkfac,zklogl,zklogu,zvclgi)
    !! This subroutine reads the EQ6 input file in compact ("W") format
    !! for version 7 (versions 7.0-7.1).
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

    integer :: ifile
    integer :: ioscan
    integer :: itermx
    integer :: jtemp
    integer :: kct
    integer :: kdim
    integer :: ksq
    integer :: kmt
    integer :: kprs
    integer :: ksplmx
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

    real(kind=8) :: cplim
    real(kind=8) :: dlzmx1
    real(kind=8) :: dlzmx2
    real(kind=8) :: dlzidp
    real(kind=8) :: dzpllg
    real(kind=8) :: dzplot
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

    ! Local variable declarations.
    integer :: i
    integer :: iktb
    integer :: j
    integer :: kc
    integer :: kcol
    integer :: n
    integer :: ncb
    integer :: nrc
    integer :: nsrt
    integer :: nxrt

    character(len=80) :: uline
    character(len=24) :: ux
    character(len=24) :: uxx
    character(len=8) :: uendit

    real(kind=8) :: xx

    data uendit /'endit.  '/

    qrderr = .false.

    ! Main title.
    ! Note: if there are exactly ntitpa lines in the title, the title
    ! need not be terminated by an 'endit.'. The 'endit.' if present
    ! is here not considered to be part of the title. The secondary
    ! title is treated in the same manner.
    qend = .false.
    read (ninpts,1100,end=100,err=990) uline
1100 format(a80)

    go to 105

100 continue
    qend = .true.
    go to 999

105 continue
    n = 1

    if (uline(1:8) .eq. uendit(1:8)) then
        go to 120
    end if

    utitl1(1) = uline

    do 110 n = 2,ntitmx
        read (ninpts,1100,err=990) uline

        if (uline(1:8) .eq. uendit(1:8)) then
            go to 120
        end if

        utitl1(n) = uline
110 continue

        n = n + 1

        read (ninpts,1100,err=990) uline
        call locase(uline)
        j = index(uline,'nmodl1=')

        if (j .gt. 0) then
            backspace(ninpts)
        else
            write (nttyo,1015) ntitmx
1015 format(/' * Error - (XCON6/rd6w7) Have too many lines in the',/7x,'main title. The code is only dimensioned for ',i4,/7x,'lines. Reduce the size of the title or increase the',/7x,'dimensioning parameter ntitpa.')

            go to 990
        end if

120 continue
        ntitl1 = n - 1

        ! Nmodl option switches.
        read (ninpts,1110,err=990) nmodl1,nmodl2
1110 format(12x,i2,22x,i2)

        ! Temperature parameters.
        ! Note: ttk(1) = tk1, etc.
        read (ninpts,1130,err=990) tempc0,jtemp,(ttk(i), i = 1,3)
1130 format(12x,e12.5,12x,i2,/3(12x,e12.5))

        ! Zi and time parameters.
        read (ninpts,1150,err=990) zistrt,zimax,tstrt,timemx,kstpmx,cplim
1150 format(2(12x,e12.5),/2(12x,e12.5),/12x,i12,12x,e12.5)

        ! Print interval parameters.
        read (ninpts,1170,err=990) dzprnt,dzprlg,ksppmx
1170 format (2(12x,e12.5),12x,i5)

        ! Plot interval parameters.
        read (ninpts,1190,err=990) dzplot,dzpllg,ksplmx
1190 format (2(12x,e12.5),12x,i5)

        ! Ifile.
        read (ninpts,1210,err=990) ifile
1210 format(12x,i2,22x,i2)

        ! Iopt option switches.
        ! Note: iopt(1) = iopt1, etc.
        read (ninpts,1240,err=990) (iopt(i), i = 1,20)
1240 format(12x,10i5)

        ! Iopr option switches.
        ! Note: iopr(1) = iopr1, etc.
        read (ninpts,1260,err=990) (iopr(i), i = 1,20)
1260 format(12x,10i5)

        ! Iodb option switches.
        ! Note: iodb(1) = iodb1, etc.
        read (ninpts,1280,err=990) (iodb(i), i = 1,20)
1280 format(12x,10i5)

        ! Number of nxopt options.
        read (ninpts,1300,err=990) nxopt
1300 format(12x,i2)

        if (nxopt .gt. nxopmx) then
            write (nttyo,1315) nxopmx
1315 format(/' * Error - (XCON6/rd6w7) Have too many mineral',/7x,'subset-selection suppression options. The code is',/7x,'only dimensioned for ',i3,' such options. Reduce the',/7x,'number of options or increase the dimensioning',/7x,'parameter nxoppa.')

            go to 990
        end if

        ! Nxopt options.
        if (nxopt .gt. 0) then
            do 140 n = 1,nxopt
                read (ninpts,1320,err=990) uxopt(n),uxct16(n)
1320 format(12x,a6,1x,a8)

140 continue

                read (ninpts,1340,err=990) nxopex
1340 format(12x,i2)

                if (nxopex .gt. nxpemx) then
                    write (nttyo,1345) nxpemx
1345 format(/' * Error - (XCON6/rd6w7) Have too many',/7x,'exceptions specified to the mineral subset-selection',/7x,'suppression options. The code is only dimensioned',/7x,'for ',i3,'exceptions. Reduce the number of exceptions',/7x,'or increase the dimensioning parameter nxoppa.')

                    go to 990
                end if

                if (nxopex .gt. 0) then
                    do 150 n = 1,nxopex
                        read (ninpts,1360,err=990) uxopex(n)
1360 format(12x,a24)

150 continue
                    end if
                end if

                ! Number of nffg options.
                read (ninpts,1380,err=990) nffg
1380 format(12x,i2)

                if (nffg .gt. nffgmx) then
                    write (nttyo,1385) nffgmx
1385 format(/' * Error - (XCON6/rd6w7) Have too many gases whose',/7x,'fugacities are to be fixed. The code is only dimensioned',/7x,'for ',i4,' such gases. Reduce the number of gases or',/7x,'increase the dimensioning parameter nffgpa.')

                    go to 990
                end if

                ! Nffg options.
                if (nffg .gt. 0) then
                    do 160 n = 1,nffg
                        read (ninpts,1400,err=990) uffg(n),moffg(n),xlkffg(n)
1400 format(12x,a12,12x,e12.5,12x,e12.5)

160 continue
                    end if

                    ! Number of reactants.
                    read (ninpts,1420,err=990) nrct
1420 format(12x,i2)

                    if (nrct .gt. nrctmx) then
                        write (nttyo,1425) nrctmx
1425 format(/' * Error - (XCON6/rd6w7) Have too many reactants',/7x,'The code is only dimensioned for ',i4,' reactants.',/7x,'Reduce the number of reactants or increase the',/7x,'dimensioning parameter nrctpa.')

                        go to 990
                    end if

                    ! Reactants.
                    nsrt = 0
                    nxrt = 0

                    if (nrct .gt. 0) then
                        do 390 nrc = 1,nrct
                            ! Name, flags, and masses.
                            read (ninpts,1440,err=990) ureac(nrc),jcode(nrc),jreac(nrc),morr(nrc),modr(nrc)
1440 format(12x,a24,/12x,i2,22x,i2,/2(12x,e12.5))

                            if (jcode(nrc) .eq. 1) then
                                ! Solid solution compositions.
                                nxrt = nxrt + 1

                                if (nxrt .gt. nxrtmx) then
                                    write (nttyo,1445) nxrtmx
1445 format(/' * Error - (XCON6/rd6w7) Have too many solid',/7x,'solution reactants. The code is only dimensioned',/7x,'for ',i4,' such reactants. Reduce the number of',/7x,'such reactants or increase the dimensioning',' parameter nxrtpa.')

                                    go to 990
                                end if

                                iktb = 0

                                do 170 i = 1,iktmax + 1
                                    read (ninpts,1460,err=990) ux,xx
1460 format(3x,a24,3x,e12.5)

                                    if (ux(1:8) .eq. uendit(1:8)) then
                                        go to 180
                                    end if

                                    iktb = iktb + 1

                                    if (iktb .gt. iktmax) then
                                        write (nttyo,1446) ureac(nrc),iktmax
1446 format(/' * Error - (XCON6/rd6w7) Have too many',' end-members',/7x,'in the solid solution reactant',' "',a24,'".',/7x,'The code is only dimensioned for ',i4,' end-members per',/7x,'solid solution. Reduce',' the number of end-members or',/7x,'increase the dimensioning parameter iktpar.')

                                        go to 990
                                    end if

                                    uendb(iktb,nxrt) = ux
                                    rxbarb(iktb,nxrt) = xx
170 continue

180 continue
                                    iktbt(nxrt) = iktb
                                else if (jcode(nrc) .eq. 2) then
                                    ! Special reactant compositions.
                                    nsrt = nsrt + 1

                                    if (nsrt .gt. nsrtmx) then
                                        write (nttyo,1447) nsrtmx
1447 format(/' * Error - (XCON6/rd6w7) Have too many special',/7x,'reactants. The code is only dimensioned for ',i4,/7x,'such reactants. Reduce the number of such reactants',/7x,'or increase the dimensioning parameter nsrtpa.')

                                        go to 990
                                    end if

                                    read (ninpts,1480,err=990) vreac(nrc)
1480 format(12x,e12.5)

                                    ncb = 0

                                    do 190 n = 1,nctmax + 1
                                        read (ninpts,1500,err=990) ux,xx
1500 format(3x,a8,1x,e25.15)

                                        if (ux(1:8) .eq. uendit(1:8)) then
                                            go to 200
                                        end if

                                        ncb = ncb + 1

                                        if (ncb .gt. nctmax) then
                                            write (nttyo,1510) ureac(nrc),nctmax
1510 format(/' * Error - (XCON6/rd6w7) Have too many','  chemical',/7x,'elements in the special  reactant "',a24,'".',/7x,'The code is only dimensioned for ',i4,' elements.',/7x,'Reduce the number of elements or increase the',/7x,'dimensioning parameter nctpar.')

                                            go to 990
                                        end if

                                        uesrb(ncb,nsrt) = ux(1:8)
                                        cesrb(ncb,nsrt) = xx
190 continue

200 continue
                                        nesrbt(nsrt) = ncb
                                    end if

                                    ! Surface area parameters.
                                    read (ninpts,1520,err=990) nsk(nrc),sk(nrc),fk(nrc)
1520 format(12x,i2,22x,e12.5,12x,e12.5)

                                    ! Rate law codes.
                                    read (ninpts,1540,err=990) nrk(1,nrc),nrk(2,nrc)
1540 format(12x,i2,22x,i2)

                                    ! Rate law parameters, forward direction.
                                    if (nrk(1,nrc) .eq. 1) then
                                        imech(1,nrc) = 3

                                        if (imech(1,nrc) .gt. imchmx) then
                                            write (nttyo,1550) ureac(nrc),imchmx
1550 format(/' * Error - (XCON6/rd6w7) Have too many rate',/7x,'constants in the forward direction rate law for',/7x,'reactant "',a24,'". The code is only',/7x,'dimensioned for ',i2,' rate constants per rate law.',/7x,'Reduce the number of rate constants or increase the',/7x,'dimensioning parameter imchpa.')

                                            go to 990
                                        end if

                                        read (ninpts,1560,err=990) (rk0(j,1,nrc), j = 1,3)
1560 format (3(12x,e12.5))
                                    else if (nrk(1,nrc) .eq. 2) then
                                        read (ninpts,1580,err=990) imech(1,nrc)
1580 format(12x,i2)

                                        if (imech(1,nrc) .gt. imchmx) then
                                            write (nttyo,1550) ureac(nrc),imchmx
                                            go to 990
                                        end if

                                        do 220 i = 1,imech(1,nrc)
                                            read (ninpts,1600,err=990) rk0(i,1,nrc),trk0(i,1,nrc),iact(i,1,nrc)
1600 format(12x,e12.5,12x,e12.5,12x,i2)

                                            read (ninpts,1620,err=990) eact(i,1,nrc),hact(i,1,nrc)
1620 format(12x,e12.5,12x,e12.5)

                                            read (ninpts,1650,err=990) ndact(i,1,nrc),csigma(i,1,nrc)
1650 format(12x,i2,22x,e12.5)

                                            if (ndact(i,1,nrc) .gt. ndctmx) then
                                                write (nttyo,1655) ureac(nrc),i,ndctmx
1655 format(/' * Error - (XCON6/rd6w7) Have too many',/7x,'species in the activity product in term ',i2,/7x,'of the forward direction rate law for reactant',/7x,'"',a24,'". The code is only dimensioned for ',i3,/7x,'such species. Reduce the number of such species',/7x,'or increase the dimensioning parameter ndctpa.')

                                                go to 990
                                            end if

                                            if (ndact(i,1,nrc) .gt. 0) then
                                                do 210 n = 1,ndact(i,1,nrc)
                                                    read (ninpts,1670,err=990) udac(n,i,1,nrc),cdac(n,i,1,nrc)
1670 format(12x,a8,16x,e12.5)

210 continue
                                                end if

220 continue
                                            else if (nrk(1,nrc) .eq. 3) then
                                                imech(1,nrc) = 1

                                                if (imech(1,nrc) .gt. imchmx) then
                                                    write (nttyo,1550) ureac(nrc),imchmx
                                                    go to 990
                                                end if

                                                i = 1
                                                read (ninpts,1600,err=990) rk0(i,1,nrc),trk0(i,1,nrc),iact(i,1,nrc)
                                                read (ninpts,1620,err=990) eact(i,1,nrc),hact(i,1,nrc)
                                            else if (nrk(1,nrc) .eq. 4) then
                                                read (ninpts,1580,err=990) imech(1,nrc)

                                                if (imech(1,nrc) .gt. imchmx) then
                                                    write (nttyo,1550) ureac(nrc),imchmx
                                                    go to 990
                                                end if

                                                do 240 i = 1,imech(1,nrc)
                                                    read (ninpts,1600,err=990) rk0(i,1,nrc),trk0(i,1,nrc),iact(i,1,nrc)
                                                    read (ninpts,1620,err=990) eact(i,1,nrc),hact(i,1,nrc)
                                                    read (ninpts,1650,err=990) ndact(i,1,nrc)

                                                    if (ndact(i,1,nrc) .gt. ndctmx) then
                                                        write (nttyo,1655) ureac(nrc),i,ndctmx
                                                        go to 990
                                                    end if

                                                    if (ndact(i,1,nrc) .gt. 0) then
                                                        do 230 n = 1,ndact(i,1,nrc)
                                                            read (ninpts,1670,err=990) udac(n,i,1,nrc),cdac(n,i,1,nrc)
230 continue
                                                        end if

240 continue
                                                    end if

                                                    ! Rate law parameters, backward direction.
                                                    if (nrk(2,nrc) .eq. 1) then
                                                        imech(2,nrc) = 3

                                                        if (imech(2,nrc) .gt. imchmx) then
                                                            write (nttyo,1673) ureac(nrc),imchmx
1673 format(/' * Error - (XCON6/rd6w7) Have too many rate',/7x,'constants in the backward direction rate law for',/7x,'reactant "',a24,'". The code is only',/7x,'dimensioned for ',i2,' rate constants per rate law.',/7x,'Reduce the number of rate constants or increase the',/7x,'dimensioning parameter imchpa.')

                                                            go to 990
                                                        end if

                                                        read (ninpts,1560,err=990) (rk0(j,2,nrc), j = 1,3)
                                                    else if (nrk(2,nrc) .eq. 2) then
                                                        read (ninpts,1580,err=990) imech(2,nrc)

                                                        if (imech(2,nrc) .gt. imchmx) then
                                                            write (nttyo,1673) ureac(nrc),imchmx
                                                            go to 990
                                                        end if

                                                        do 260 i = 1,imech(2,nrc)
                                                            read (ninpts,1600,err=990) rk0(i,2,nrc),trk0(i,2,nrc),iact(i,2,nrc)
                                                            read (ninpts,1620,err=990) eact(i,2,nrc),hact(i,2,nrc)
                                                            read (ninpts,1650,err=990) ndact(i,2,nrc),csigma(i,2,nrc)

                                                            if (ndact(i,2,nrc) .gt. ndctmx) then
                                                                write (nttyo,1675) ureac(nrc),i,ndctmx
1675 format(/' * Error - (XCON6/rd6w7) Have too many',/7x,'species in the activity product in term ',i2,/7x,'of the backward direction rate law for reactant',/7x,'"',a24,'". The code is only dimensioned for ',i3,/7x,'such species. Reduce the number of such species',/7x,'or increase the dimensioning parameter ndctpa.')

                                                                go to 990
                                                            end if

                                                            if (ndact(i,2,nrc) .gt. 0) then
                                                                do 250 n = 1,ndact(i,2,nrc)
                                                                    read (ninpts,1670,err=990) udac(n,i,2,nrc),cdac(n,i,2,nrc)
250 continue
                                                                end if

260 continue
                                                            else if (nrk(2,nrc) .eq. 3) then
                                                                imech(2,nrc) = 1

                                                                if (imech(2,nrc) .gt. imchmx) then
                                                                    write (nttyo,1673) ureac(nrc),imchmx
                                                                    go to 990
                                                                end if

                                                                i = 1
                                                                read (ninpts,1600,err=990) rk0(i,2,nrc),trk0(i,2,nrc),iact(i,2,nrc)
                                                                read (ninpts,1620,err=990) eact(i,2,nrc),hact(i,2,nrc)
                                                            else if (nrk(2,nrc) .eq. 4) then
                                                                read (ninpts,1580,err=990) imech(2,nrc)

                                                                if (imech(2,nrc) .gt. imchmx) then
                                                                    write (nttyo,1673) ureac(nrc),imchmx
                                                                    go to 990
                                                                end if

                                                                do 280 i = 1,imech(2,nrc)
                                                                    read (ninpts,1600,err=990) rk0(i,2,nrc),trk0(i,2,nrc),iact(i,2,nrc)
                                                                    read (ninpts,1620,err=990) eact(i,2,nrc),hact(i,2,nrc)
                                                                    read (ninpts,1650,err=990) ndact(i,2,nrc)

                                                                    if (ndact(i,2,nrc) .gt. ndctmx) then
                                                                        write (nttyo,1675) ureac(nrc),i,ndctmx
                                                                        go to 990
                                                                    end if

                                                                    if (ndact(i,2,nrc) .gt. 0) then
                                                                        do 270 n = 1,ndact(i,2,nrc)
                                                                            read (ninpts,1670,err=990) udac(n,i,2,nrc),cdac(n,i,2,nrc)
270 continue
                                                                        end if

280 continue
                                                                    end if

390 continue
                                                                end if

                                                                ! Dump interval.
                                                                read (ninpts,2000,err=990) dlzidp
2000 format(3(12x,e12.5))

                                                                ! Tolerances.
                                                                read (ninpts,2020,err=990) tolbt,toldl,tolx,tolsat,tolsst
2020 format(3(12x,e12.5))

                                                                ! Setscrew parameters.
                                                                ! Note: sscrew(1) = screw1, etc.
                                                                read (ninpts,2040,err=990) (sscrew(i), i = 1,6)
2040 format(3(12x,e12.5))

                                                                ! Z parameters.
                                                                read (ninpts,2060,err=990) zklogu,zklogl,zkfac
2060 format(3(12x,e12.5))

                                                                ! Step size limits and maximum order.
                                                                read (ninpts,2080,err=990) dlzmx1,dlzmx2,nordlm
2080 format (2(12x,e12.5),12x,i2)

                                                                ! Maximum number of iterations and maximum number of phase
                                                                ! assemblage tries.
                                                                read (ninpts,2100,err=990) itermx,ntrymx
2100 format(12x,i2,22x,i2)

                                                                ! Slide and scan control parameters.
                                                                read (ninpts,2120,err=990) npslmx,nsslmx,ioscan
2120 format(12x,i2,22x,i2,22x,i2)

                                                                ! Process the bottom half of the current output file.
                                                                ! Secondary title.
                                                                do 400 n = 1,ntitmx
                                                                    read (ninpts,3000,err=990) uline
3000 format(a80)

                                                                    if (uline(1:8) .eq. uendit(1:8)) then
                                                                        go to 410
                                                                    end if

                                                                    utitl2(n) = uline
400 continue

                                                                    n = n + 1

                                                                    read (ninpts,1100,err=990) uline
                                                                    call locase(uline)
                                                                    j = index(uline,'tempci=')

                                                                    if (j .gt. 0) then
                                                                        backspace(ninpts)
                                                                    else
                                                                        write (nttyo,3002) ntitmx
3002 format(/' * Error - (XCON6/rd6w7) Have too many lines in the',/7x,'secondary title. The code is only dimensioned for ',i4,/7x,'lines. Reduce the size of the title or increase the',/7x,'dimensioning parameter ntitpa.')

                                                                        go to 990
                                                                    end if

410 continue
                                                                    ntitl2 = n - 1

                                                                    ! Original temperature.
                                                                    read (ninpts,3010,err=990) tempci
3010 format(12x,e12.5)

                                                                    ! Number of nxmod options.
                                                                    read (ninpts,3030,err=990) nxmod
3030 format(12x,i2)

                                                                    if (nxmod .gt. nxmdmx) then
                                                                        write (nttyo,3035) nxmdmx
3035 format(/' * Error - (XCON6/rd6w7) Have too many nxmod',/7x,'alter/suppress options. The code is only dimensioned',/7x,'for ',i3,' such options. Reduce the number of such',' options',/7x,'or increase the dimensioning parameter',' nxmdpa.')

                                                                        go to 990
                                                                    end if

                                                                    ! Nxmod options.
                                                                    if (nxmod .gt. 0) then
                                                                        do 420 n = 1,nxmod
                                                                            read (ninpts,3050,err=990) uxmd24(n),jxmod(n),kxmod(n),xlkmod(n)
3050 format(12x,a24,/12x,i2,22x,i2,22x,e12.5)

420 continue
                                                                        end if

                                                                        ! Iopg options.
                                                                        ! Note: iopg(1) = iopg1, etc.
                                                                        read (ninpts,3070,err=990) (iopg(i), i = 1,10)
3070 format (12x,i2,22x,i2,22x,i2)

                                                                        ! Index limits.
                                                                        read (ninpts,3090,err=990) kct,ksq,kmt,kxt,kdim,kprs
3090 format(3(12x,i2,10x),/3(12x,i2,10x))

                                                                        if (kct .gt. nctmax) then
                                                                            write (nttyo,3095) nctmax
3095 format(/' * Error - (XCON6/rd6w7) Have too many chemical',/7x,'elements present. The code is only dimensioned',/7x,'for ',i3,' elements. Reduce the number of elements',/7x,'or increase the dimensioning parameter nctpar.')

                                                                            go to 990
                                                                        end if

                                                                        if (kdim .gt. kmax) then
                                                                            write (nttyo,3100) kmax
3100 format(/' * Error - (XCON6/rd6w7) Have too many master',/7x,'variables. The code is only dimensioned for ',i3,/7x,'master variables. Reduce the number of such variables',/7x,'or increase the dimensioning parameter kpar.')

                                                                            go to 990
                                                                        end if

                                                                        ! Balance totals.
                                                                        do 430 kc = 1,kct
                                                                            read (ninpts,3120,err=990) uelemb(kc),mteb(kc),mteaqb(kc)
3120 format(3x,a8,7x,e25.15,3x,e25.15)

430 continue

                                                                            read (ninpts,3140,err=990) ux,electr
3140 format(3x,a8,7x,e25.15)

                                                                            ! Basis variable data.
                                                                            do 440 kcol = 1,kdim
                                                                                read (ninpts,3160,err=990) unrms(kcol),undms(kcol),zvclgi(kcol)
3160 format(3x,a18,1x,a18,1x,e25.15)

440 continue

                                                                                ! Physically removed system data.
                                                                                nprmn = 0
                                                                                nprmx = 0

                                                                                if (kprs .gt. 0) then
                                                                                    do 450 n = 1,nprsmx
                                                                                        read (ninpts,3180,err=990) ux,xx
3180 format(3x,a18,1x,e25.15)

                                                                                        if (ux(1:8) .eq. uendit(1:8)) then
                                                                                            go to 460
                                                                                        end if

                                                                                        nprmn = nprmn + 1

                                                                                        if (nprmn .gt. nprsmx) then
                                                                                            write (nttyo,3190) nprsmx
3190 format(/' * Error - (XCON6/rd6w7) Have too many mineral',/7x,'species in the physically removed system. The code is',/7x,'only dimensioned for ',i3,' such species. Reduce the',/7x,'number of such species or increase the dimensioning',/7x,'parameter nprspa.')

                                                                                            go to 990
                                                                                        end if

                                                                                        uprs(nprmn)(1:24) = ux
                                                                                        uprs(nprmn)(25:48) = ux
                                                                                        mprs(nprmn) = xx
450 continue

460 continue

                                                                                        nprmx = nprmn

                                                                                        do 480 n = 1,nprsmx
                                                                                            read (ninpts,3200,err=990) ux
3200 format(3x,a24)

                                                                                            if (ux(1:8) .eq. uendit(1:8)) then
                                                                                                go to 490
                                                                                            end if

                                                                                            do 470 i = 1,iktmax + 1
                                                                                                read (ninpts,3220,err=990) uxx,xx
3220 format(3x,a18,1x,e25.15)

                                                                                                if (uxx(1:8) .eq. uendit(1:8)) then
                                                                                                    go to 480
                                                                                                end if

                                                                                                nprmx = nprmx + 1

                                                                                                if (nprmn .gt. nprsmx) then
                                                                                                    write (nttyo,3190) nprsmx
                                                                                                    go to 990
                                                                                                end if

                                                                                                uprs(nprmx)(1:24) = ux
                                                                                                uprs(nprmx)(25:48) = uxx
                                                                                                mprs(nprmx) = xx
470 continue

480 continue

490 continue
                                                                                            end if

                                                                                            go to 999

990 continue
                                                                                            qrderr = .true.

999 continue
                                                                                        end subroutine rd6w7