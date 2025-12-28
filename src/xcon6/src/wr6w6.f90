subroutine wr6w6(cdac,cesrb,cplim,csigma,dlzmx1,dlzmx2,dlzidp,dzpllg,dzplot,dzprlg,dzprnt,electr,fk,ifile,iktbt,iktmax,imchmx,imech,iodb,iopg,iopr,iopt,ioscan,itermx,jcode,jreac,jtemp,jxmod,kct,kdim,kmax,kmt,kprs,ksq,ksplmx,ksppmx,kstpmx,kxmod,kxt,modr,moffg,morr,mprs,mteaqb,mteb,nctmax,ndact,ndctmx,nesrbt,newin,nffg,nffgmx,nmodl1,nmodl2,nodbmx,nopgmx,noprmx,noptmx,nordlm,npslmx,nprmn,nprmx,nprsmx,nrct,nrctmx,nrk,nsk,nsrtmx,nsscmx,nsslmx,ntitl1,ntitl2,ntitmx,ntrymx,nttkmx,nxmdmx,nxmod,nxopex,nxopmx,nxopt,nxpemx,nxrtmx,rk0,rxbarb,sscrew,sk,tempci,tempc0,timemx,tolbt,toldl,tolsat,tolsst,tolx,tstrt,ttk,uacion,ucode,udac,uelemb,uendb,ueqlrn,ueqlst,uesrb,uffg,undms,unrms,uprs,ureac,urelno,ustage,utitl1,utitl2,uxct16,uxmd24,uxopex,uxopt,vreac,xlkffg,xlkmod,zimax,zistrt,zkfac,zklogl,zklogu,zvclgi)
    !! This subroutine writes the EQ6 input file in compact ("W") format
    !! for version 6 (versions 6.0-6.1).
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
    integer :: newin
    integer :: nffg
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
    integer :: nxmod
    integer :: nxopex
    integer :: nxopt

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
    character(len=24) :: uacion
    character(len=16) :: uxct16(nxopmx)
    character(len=8) :: uelemb(nctmax)
    character(len=8) :: uesrb(nctmax,nsrtmx)
    character(len=8) :: uxopt(nxopmx)
    character(len=8) :: ucode
    character(len=8) :: urelno
    character(len=8) :: ustage
    character(len=8) :: ueqlrn
    character(len=8) :: ueqlst

    real(kind=8) :: cdac(ndctmx,imchmx,2,nrctmx)
    real(kind=8) :: cesrb(nctmax,nsrtmx)
    real(kind=8) :: csigma(imchmx,2,nrctmx)
    real(kind=8) :: fk(nrctmx)
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
    integer :: npr
    integer :: nrc
    integer :: nsrt
    integer :: nxrt

    character(len=24) :: uxx
    character(len=8) :: uendit

    data uendit /'endit.  '/

    ! Write the new input file.
    ! Title.
    do 110 n = 1,ntitl1
        write (newin,1100) utitl1(n)
1100 format(a80)

110 continue

        if (ntitl1 .lt. ntitmx) then
            write (newin,1110)
        end if

1110 format('endit.')

        ! Nmodl option switches.
        write (newin,1120) nmodl1,nmodl2
1120 format(4x,'nmodl1= ',i2,14x,'nmodl2= ',i2)

        ! Temperature parameters.
        ! Note: ttk(1) = tk1, etc.
        write (newin,1140) tempc0,jtemp,(ttk(i), i = 1,3)
1140 format(4x,'tempc0= ',1pe12.5,5x,'jtemp= ',i2,/7x,'tk1= ',1pe12.5,7x,'tk2= ',1pe12.5,7x,'tk3= ',1pe12.5)

        ! Zi and time parameters.
        write (newin,1160) zistrt,zimax,tstrt,timemx,kstpmx,cplim
1160 format(4x,'zistrt= ',1pe12.5,5x,'zimax= ',1pe12.5,/5x,'tstrt= ',1pe12.5,4x,'timemx= ',1pe12.5,/4x,'kstpmx= ',i12,5x,'cplim= ',1pe12.5)

        ! Print interval parameters.
        write (newin,1180) dzprnt,dzprlg,ksppmx
1180 format(4x,'dzprnt= ',1pe12.5,4x,'dzprlg= ',1pe12.5,4x,'ksppmx= ',i5)

        ! Plot interval parameters.
        write (newin,1200) dzplot,dzpllg,ksplmx
1200 format (4x,'dzplot= ',1pe12.5,4x,'dzpllg= ',1pe12.5,4x,'ksplmx= ',i5)

        ! Ifile.
        write (newin,1220) ifile
1220 format(5x,'ifile= ',i2)

        ! Comment header for iopt, iopr, and iodb option switches.
        write (newin,1230)
1230 format('*',15x,'1    2    3    4    5    6    7    8    9   10')

        ! Iopt option switches.
        ! Note: iopt(1) = iopt1, etc.
        write (newin,1250) (iopt(i), i = 1,20)
1250 format(2x,'iopt1-10= ',10i5,/1x,'iopt11-20= ',10i5)

        ! Iopr option switches.
        ! Note: iopr(1) = iopY1, etc.
        write (newin,1270) (iopr(i), i = 1,20)
1270 format(2x,'iopr1-10= ',10i5,/1x,'iopr11-20= ',10i5)

        ! Iodb option switches.
        ! Note: iodb(1) = iodb1, etc.
        write (newin,1290) (iodb(i), i = 1,20)
1290 format(2x,'iodb1-10= ',10i5,/1x,'iodb11-20= ',10i5)

        ! Number of nxopt options.
        write (newin,1310) nxopt
1310 format(5x,'nxopt= ',i2)

        ! Nxopt options.
        if (nxopt .gt. 0) then
            do 140 n = 1,nxopt
                write (newin,1330) uxopt(n),uxct16(n)
1330 format(4x,'option= ',a6,1x,a8)

140 continue

                write (newin,1350) nxopex
1350 format(4x,'nxopex= ',i2)

                if (nxopex .gt. 0) then
                    do 150 n = 1,nxopex
                        write (newin,1370) uxopex(n)
1370 format(1x,'exception= ',a8)

150 continue
                    end if
                end if

                ! Number of nffg options.
                write (newin,1390) nffg
1390 format(6x,'nffg= ',i2)

                ! Nffg options.
                if (nffg .gt. 0) then
                    do 160 n = 1,nffg
                        write (newin,1410) uffg(n),moffg(n),xlkffg(n)
1410 format(3x,'species= ',a12,5x,'moffg= ',1pe12.5,4x,'xlkffg= ',1pe12.5)

160 continue
                    end if

                    ! Number of reactants.
                    write (newin,1430) nrct
1430 format(6x,'nrct= ',i2)

                    write (newin,1435)
1435 format('*-------------------------------------------------','----------------------------')

                    ! Reactants.
                    nsrt = 0
                    nxrt = 0

                    if (nrct .gt. 0) then
                        do 390 nrc = 1,nrct
                            ! Name, flags, and masses.
                            write (newin,1450) ureac(nrc),jcode(nrc),jreac(nrc),morr(nrc),modr(nrc)
1450 format(2x,'reactant= ',a24,/5x,'jcode= ',i2,15x,'jreac= ',i2,/6x,'morr= ',1pe12.5,6x,'modr= ',1pe12.5)

                            if (jcode(nrc) .eq. 1) then
                                ! Solid solution compositions.
                                nxrt = nxrt + 1

                                do 170 iktb = 1,iktbt(nxrt)
                                    write (newin,1470) uendb(iktb,nxrt),rxbarb(iktb,nxrt)
1470 format(3x,a18,3x,1pe12.5)

170 continue

                                    write (newin,1475) uendit
1475 format(3x,a8)
                                else if (jcode(nrc) .eq. 2) then
                                    ! Special reactant compositions.
                                    nsrt = nsrt + 1
                                    write (newin,1490) vreac(nrc)
1490 format(5x,'vreac= ',1pe12.5)

                                    do 190 ncb = 1,nesrbt(nsrt)
                                        write (newin,1510) uesrb(ncb,nsrt),cesrb(ncb,nsrt)
1510 format(3x,a8,1x,1pe25.15)

190 continue

                                        write (newin,1475) uendit
                                    end if

                                    ! Surface area parameters.
                                    write (newin,1530) nsk(nrc),sk(nrc),fk(nrc)
1530 format(7x,'nsk= ',i2,18x,'sk= ',1pe12.5,8x,'fk= ',1pe12.5)

                                    write (newin,1550) nrk(1,nrc),nrk(2,nrc)
1550 format(7x,'nrk= ',i2,16x,'nrpk= ',i2)

                                    ! Rate law parameters, forward direction.
                                    if (nrk(1,nrc) .eq. 1) then
                                        i = imech(1,nrc)
                                        write (newin,1570) (rk0(j,1,nrc), j = 1,i)
1570 format(7x,'rk1= ',1pe12.5,7x,'rk2= ',1pe12.5,7x,'rk3= ',1pe12.5)
                                    else if (nrk(1,nrc) .eq. 2) then
                                        write (newin,1590) imech(1,nrc)
1590 format(5x,'imech= ',i2)

                                        do 220 i = 1,imech(1,nrc)
                                            write (newin,1610) rk0(i,1,nrc),ndact(i,1,nrc),csigma(i,1,nrc)
1610 format(7x,'rk0= ',1pe12.5,5x,'ndact= ',i2,14x,'csigma= ',1pe12.5)

                                            if (ndact(i,1,nrc) .gt. 0) then
                                                do 210 n = 1,ndact(i,1,nrc)
                                                    write (newin,1680) udac(n,i,1,nrc),cdac(n,i,1,nrc)
1680 format(6x,'udac= ',a8,10x,'cdac= ',1pe12.5)

210 continue
                                                end if

220 continue
                                            else if (nrk(1,nrc) .eq. 3) then
                                                write (newin,1570) (rk0(j,1,nrc), j = 1,3)
                                            else if (nrk(1,nrc) .eq. 4) then
                                                write (newin,1590) imech(1,nrc)

                                                do 240 i = 1,imech(1,nrc)
                                                    write (newin,1690) rk0(i,1,nrc),ndact(i,1,nrc)
1690 format(7x,'rk0= ',1pe12.5,5x,'ndact= ',i2)

1800 format(5x,'ndact= ',i2)

                                                    if (ndact(i,1,nrc) .gt. 0) then
                                                        do 230 n = 1,ndact(i,1,nrc)
                                                            write (newin,1680) udac(n,i,1,nrc),cdac(n,i,1,nrc)
230 continue
                                                        end if

240 continue
                                                    else if (nrk(1,nrc) .eq. 5) then
                                                        write (newin,1810) (rk0(j,1,nrc), j = 1,2)
1810 format(7x,'rk1= ',1pe12.5,7x,'rk2= ',1pe12.5)
                                                    end if

                                                    ! Rate law parameters, backward direction.
                                                    if (nrk(2,nrc) .eq. 1) then
                                                        i = imech(2,nrc)
                                                        write (newin,1830) (rk0(j,2,nrc), j = 1,i)
1830 format(6x,'rpk1= ',1pe12.5,6x,'rpk2= ',1pe12.5,6x,'rpk3= ',1pe12.5)
                                                    else if (nrk(2,nrc) .eq. 2) then
                                                        write (newin,1840) imech(2,nrc)
1840 format(4x,'ipmech= ',i2)

                                                        do 260 i = 1,imech(2,nrc)
                                                            write (newin,1850) rk0(i,2,nrc),ndact(i,2,nrc),csigma(i,2,nrc)
1850 format(7x,'rpk= ',1pe12.5,4x,'npdact= ',i2,13x,'cpsigma= ',1pe12.5)

                                                            if (ndact(i,2,nrc) .gt. 0) then
                                                                do 250 n = 1,ndact(i,2,nrc)
                                                                    write (newin,1880) udac(n,i,2,nrc),cdac(n,i,2,nrc)
1880 format(5x,'updac= ',a8,9x,'cpdac= ',1pe12.5)

250 continue
                                                                end if

260 continue
                                                            else if (nrk(2,nrc) .eq. 3) then
                                                                write (newin,1830) rk0(1,2,nrc)
1920 format(6x,'rpk1= ',1pe12.5)
                                                            else if (nrk(2,nrc) .eq. 4) then
                                                                write (newin,1840) imech(2,nrc)

                                                                do 280 i = 1,imech(2,nrc)
                                                                    write (newin,1930) rk0(i,2,nrc),ndact(i,2,nrc)
1930 format(7x,'rpk= ',1pe12.5,4x,'npdact= ',i2)

1940 format(4x,'npdact= ',i2)

                                                                    if (ndact(i,2,nrc) .gt. 0) then
                                                                        do 270 n = 1,ndact(i,2,nrc)
                                                                            write (newin,1880) udac(n,i,2,nrc),cdac(n,i,2,nrc)
270 continue
                                                                        end if

280 continue
                                                                    else if (nrk(2,nrc) .eq. 5) then
                                                                        write (newin,1810) (rk0(j,2,nrc), j = 1,2)
                                                                    end if

                                                                    write (newin,1435)
390 continue
                                                                end if

                                                                ! Dump interval.
                                                                write (newin,2010) dlzidp
2010 format(4x,'dlzidp= ',1pe12.5)

                                                                ! Tolerances.
                                                                write (newin,2030) tolbt,toldl,tolx,tolsat,tolsst
2030 format(5x,'tolbt= ',1pe12.5,5x,'toldl= ',1pe12.5,6x,'tolx= ',1pe12.5,/4x,'tolsat= ',1pe12.5,4x,'tolsst= ',1pe12.5)

                                                                ! Setscrew parameters.
                                                                ! Note: sscrew(1) = screw1, etc.
                                                                write (newin,2050) (sscrew(i), i = 1,6)
2050 format(4x,'screw1= ',1pe12.5,4x,'screw2= ',1pe12.5,4x,'screw3= ',1pe12.5,/4x,'screw4= ',1pe12.5,4x,'screw5= ',1pe12.5,4x,'screw6= ',1pe12.5)

                                                                ! Z parameters.
                                                                write (newin,2070) zklogu,zklogl,zkfac
2070 format(4x,'zklogu= ',f12.3,4x,'zklogl= ',f12.3,5x,'zkfac= ',f12.3)

                                                                ! Step size limits and maximum order.
                                                                write (newin,2090) dlzmx1,dlzmx2,nordlm
2090 format(4x,'dlzmx1= ',1pe12.5,4x,'dlzmx2= ',1pe12.5,4x,'nordlm= ',i2)

                                                                ! Maximum number of iterations and maximum number of phase
                                                                ! assemblage tries.
                                                                write (newin,2110) itermx,ntrymx
2110 format(4x,'itermx= ',i2,14x,'ntrymx= ',i2)

                                                                ! Slide and scan control parameters.
                                                                write (newin,2130) npslmx,nsslmx,ioscan
2130 format(4x,'npslmx= ',i2,14x,'nsslmx= ',i2,14x,'ioscan= ',i2)

                                                                write (newin,1435)

                                                                ! Process the bottom half of the current output file.
                                                                ! Write recaptured data in comments on the originating code, if
                                                                ! known.
                                                                if (ucode(1:3).eq.'EQ6' .or. ucode(1:3).eq.'eq6') then
                                                                    write (newin,2180) urelno,ustage,ueqlrn,ueqlst
2180 format('* pickup file written by eq6.',a4,a6,/'*  supported by eqlib.',a4,a6)
                                                                end if

                                                                if (ucode(1:5).eq.'EQ3NR' .or. ucode(1:5).eq.'eq3nr') then
                                                                    write (newin,2190) urelno,ustage,ueqlrn,ueqlst
2190 format('* pickup file written by eq3nr.',a4,a6,/'*  supported by eqlib.',a4,a6)
                                                                end if

                                                                ! Title.
                                                                do 400 n = 1,ntitl2
                                                                    write (newin,3000) utitl2(n)
3000 format(a80)

400 continue

                                                                    if (ntitl2 .lt. ntitmx) then
                                                                        write (newin,1110)
                                                                    end if

                                                                    ! Ion that defines the equivalent stoichiometric ionic strength.
                                                                    write (newin,3015) uacion
3015 format(4x,'uacion= ',a24)

                                                                    ! Original temperature.
                                                                    write (newin,3020) tempci
3020 format(4x,'tempci= ',1pe12.5)

                                                                    ! Number of nxmod options.
                                                                    write (newin,3040) nxmod
3040 format(5x,'nxmod= ',i2)

                                                                    ! Nxmod options.
                                                                    if (nxmod .gt. 0) then
                                                                        do 420 n = 1,nxmod
                                                                            write (newin,3060) uxmd24(n),jxmod(n),kxmod(n),xlkmod(n)
3060 format(3x,'species= ',a24,/6x,'type= ',i2,14x,'option= ',i2,14x,'xlkmod= ',1pe12.5)

420 continue
                                                                        end if

                                                                        ! Iopg options.
                                                                        ! Note: iopg(1) = iopg1, etc.
                                                                        write(newin,3080) (iopg(i), i = 1,10)
3080 format(5x,'iopg1= ',i2,15x,'iopg2= ',i2,15x,'iopg3= ',i2,/5x,'iopg4= ',i2,15x,'iopg5= ',i2,15x,'iopg6= ',i2,/5x,'iopg7= ',i2,15x,'iopg8= ',i2,15x,'iopg9= ',i2,/4x,'iopg10= ',i2)

                                                                        ! Index limits.
                                                                        write (newin,3100) kct,ksq,kmt,kxt,kdim,kprs
3100 format(7x,'kct= ',i2,17x,'ksq= ',i2,17x,'kmt= ',i2,/7x,'kxt= ',i2,16x,'kdim= ',i2,16x,'kprs= ',i2)

                                                                        ! Balance totals.
                                                                        write (newin,3110)
3110 format('*  Component',13x,'Moles Total',17x,'Moles Aqueous')

                                                                        do 430 kc = 1,kct
                                                                            write (newin,3130) uelemb(kc),mteb(kc),mteaqb(kc)
3130 format(3x,a8,7x,1pe25.15,3x,1pe25.15)

430 continue

                                                                            write (newin,3150) electr
3150 format(3x,'electr  ',7x,1pe25.15)

                                                                            ! Basis variable data.
                                                                            do 440 kcol = 1,kdim
                                                                                write (newin,3170) unrms(kcol),undms(kcol),zvclgi(kcol)
3170 format(3x,a18,1x,a18,1x,1pe25.15)

440 continue

                                                                                ! Physically removed system data.
                                                                                if (kprs .gt. 0) then
                                                                                    do 450 npr = 1,nprmn
                                                                                        write (newin,3190) uprs(npr)(1:24),mprs(npr)
3190 format(3x,a18,1x,1pe25.15)

450 continue

                                                                                        write (newin,3200) uendit
3200 format(3x,a8)

                                                                                        uxx  = '        '

                                                                                        do 480 npr = nprmn + 1,nprmx
                                                                                            if (uprs(npr)(1:8) .eq. '        ') then
                                                                                                go to 490
                                                                                            end if

                                                                                            if (uprs(npr)(1:24) .eq. uxx(1:24)) then
                                                                                                write (newin,3200) uendit
                                                                                                write (newin,3210) uprs(npr)(1:24)
3210 format(3x,a24)
                                                                                            end if

                                                                                            write (newin,3230) uprs(npr)(25:48),mprs(npr)
3230 format(3x,a18,1x,1pe25.15)

480 continue

490 continue

                                                                                            write (newin,3200) uendit
                                                                                            write (newin,3200) uendit
                                                                                        end if

999 continue
                                                                                    end subroutine wr6w6