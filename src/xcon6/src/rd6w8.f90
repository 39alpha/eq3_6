subroutine rd6w8(awmaxi,awmini,cbsri,cdac,cesri,cgexj,csigma,dlaplo,dlaprn,dleplo,dleprn,dlhplo,dlhprn,dloplo,dloprn,dltpll,dltplo,dltprl,dltprn,dlxdmp,dlxmx0,dlxpll,dlxplo,dlxprl,dlxprn,eact,egersi,ehmaxi,ehmini,electr,fkrc,iact,ibsrti,iesrti,ietmax,iktmax,imchmx,imech,iodb,iopg,iopr,iopt,igerti,itermx,ixrti,jcode,jgerti,jetmax,jflgi,jgext,jpress,jreac,jtemp,kbt,kct,kdim,kmax,kmt,kprs,ksplmx,ksppmx,kstpmx,kxmod,kxt,hact,modr,moffg,morr,mprphi,mprspi,mtbaqi,mtbi,mwtges,nbti,nbtmax,nbt1mx,nctmax,ndact,ndctmx,nert,nertmx,net,netmax,nffg,nffgmx,ngexrt,ninpts,nobswt,nodbmx,nopgmx,noprmx,noptmx,nordmx,noutpt,nprob,nprpmx,nprpti,nprsmx,nprsti,nptkmx,nrct,nrctmx,nrk,nsbswt,nsk,nsrt,nsrtmx,ntitl1,ntitl2,ntitmx,ntrymx,nttkmx,nttyo,nxmdmx,nxmod,nxopex,nxopmx,nxopt,nxpemx,nxrt,nxrtmx,o2maxi,o2mini,phmaxi,phmini,pressb,pressi,ptk,qend,qgexsh,qrderr,rkb,rxbari,sfcar,ssfcar,tempcb,tempci,tgexp,timmxi,tistti,tolbt,toldl,tolsat,tolxsf,trkb,ttk,ubmtbi,ubsri,ucxri,udac,uesri,uffg,ugerji,ugermo,ugersi,ugexj,ugexmo,ugexp,ugexr,uhfgex,uobsw,uprphi,uprspi,ureac,usbsw,utitl1,utitl2,uvfgex,uxcat,uxkgex,uxmod,uxopex,uxopt,uzveci,vreac,xgersi,xhfgex,ximaxi,xistti,xlkffg,xlkgex,xlkmod,xvfgex,zgexj,zvclgi)
    !! This subroutine reads the EQ6 input file in compact ("W") format
    !! for version 8.0.
    !! This subroutine is a near-clone of EQ6/rd6inw.f. However, the
    !! present subroutine embodies only a pure read function (it does
    !! only minimal checking of what is read to ensure that what
    !! follows is readable). EQ6/rd6inw.f differs in that it also
    !! writes an instant echo of what is read to the EQ6 output file.
    !! The calling sequence of this subroutine is identical to that of
    !! EQ6/rd6inw.f, EQ6/rd6ind.f, and XCON6/rd6d8.f.
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
    integer :: nffg
    integer :: nert
    integer :: net
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

    ! Local variable declarations.
    integer :: i
    integer :: iki
    integer :: ikti
    integer :: je
    integer :: jei
    integer :: jeti
    integer :: jj
    integer :: j2
    integer :: j3
    integer :: kcol
    integer :: krow
    integer :: n
    integer :: nbi
    integer :: nci
    integer :: ne
    integer :: nei
    integer :: ner
    integer :: neti
    integer :: nn
    integer :: npi
    integer :: nrc

    integer :: ilnobl

    character(len=80) :: uline
    character(len=24) :: uxf
    character(len=24) :: uxg
    character(len=24) :: uxs
    character(len=8) :: uendit
    character(len=8) :: uxe
    character(len=8) :: ux8

    real(kind=8) :: mx
    real(kind=8) :: xx

    data uendit /'endit.  '/

    ! The following is nonsense to avoid compiler "unused variable"
    ! warnings. Here noutpt and nprob are not actually used. They are
    ! included in the calling sequence to allow it to match that of
    ! EQ6/rd6inw.f.
    noutpt = nttyo
    i = nprob

    ! qend   = .true if the end of the input file has been encountered
    ! qrderr = .true if the current problem can't be read because of
    !            a read format error or a dimensional overflow
    qend = .false.
    qrderr = .false.

    ! Main title.
    ! Note: if there are exactly ntitpa lines in the title, the title
    ! need not be terminated by an 'endit.'. The 'endit.' if present
    ! is here not considered to be part of the title.
    ntitl1 = 0
    read (ninpts,1000,end=100,err=990) uline
1000 format(a80)

    go to 110

100 continue
    qend = .true.
    go to 999

110 continue

    if (uline(1:8) .ne. uendit(1:8)) then
        ntitl1 = 1
        utitl1(1) = uline

        do n = 2,ntitmx
            read (ninpts,1000,err=990) uline

            if (uline(1:8) .eq. uendit(1:8)) then
                go to 120
            end if

            ntitl1 = n
            utitl1(n) = uline
        end do

120 continue
    end if

    ! Temperature parameters.
    ! Note: ttk(1) = ttk1, etc.
    read (ninpts,1100,err=990) jtemp,tempcb,ttk(1),ttk(2)
1100 format(12x,i2,/12x,e12.5,/2(12x,e12.5))

    ! Pressure parameters.
    ! Note: ptk(1) = ptk1, etc.
    read (ninpts,1120,err=990) jpress,pressb,ptk(1),ptk(2)
1120 format(12x,i2,/12x,e12.5,/2(12x,e12.5))

    ! Number of reactants.
    read (ninpts,1480,err=990) nrct
1480 format(12x,i2)

    if (nrct .gt. nrctmx) then
        write (nttyo,1520) nrctmx
1520 format(/' * Error - (XCON6/rd6w8) Have too many reactants',/7x,'The code is only dimensioned for ',i4,' reactants.',/7x,'Reduce the number of reactants or increase the',/7x,'dimensioning parameter nrctpa.')

        go to 990
    end if

    ! Reactants.
    nsrt = 0
    nxrt = 0
    nert = 0

    do nrc = 1,nrct
        ! Name, flags, and masses.
        read (ninpts,1530,err=990) ureac(nrc),jcode(nrc),jreac(nrc),morr(nrc),modr(nrc)
1530 format(12x,a24,/12x,i2,22x,i2,/2(12x,e12.5))

        j2 = ilnobl(ureac(nrc))

        if (jcode(nrc) .eq. 1) then
            ! Have a solid solution reactant.
            nxrt = nxrt + 1

            if (nxrt .gt. nxrtmx) then
                write (nttyo,1550) nxrtmx
1550 format(/' * Error - (XCON6/rd6w8) Have too many solid',/7x,'solution reactants. The code is only dimensioned',/7x,'for ',i4,' such reactants. Reduce the number of',/7x,'such reactants or increase the dimensioning',' parameter nxrtpa.')

                go to 990
            end if

            ikti = 0

            do i = 1,iktmax + 1
                read (ninpts,1560,err=990) uxs,xx
1560 format(3x,a24,3x,e12.5)

                if (uxs(1:8) .eq. uendit(1:8)) then
                    go to 150
                end if

                ikti = ikti + 1

                if (ikti .gt. iktmax) then
                    write (nttyo,1590) ureac(nrc)(1:j2),iktmax
1590 format(/' * Error - (XCON6/rd6w8) Have too many',' end-members',/7x,'in the solid solution reactant ',a,'.',/7x,'The code is only dimensioned for ',i4,' end-members per',/7x,'solid solution. Reduce',' the number of end-members or',/7x,'increase the dimensioning parameter iktpar.')

                    go to 990
                end if

                ucxri(ikti,nxrt) = uxs
                rxbari(ikti,nxrt) = xx
            end do

150 continue
            ixrti(nxrt) = ikti
        else if (jcode(nrc) .eq. 2) then
            ! Have a special reactant.
            nsrt = nsrt + 1

            if (nsrt .gt. nsrtmx) then
                write (nttyo,1600) nsrtmx
1600 format(/' * Error - (XCON6/rd6w8) Have too many special',/7x,'reactants. The code is only dimensioned for ',i4,/7x,'such reactants. Reduce the number of such reactants',/7x,'or increase the dimensioning parameter nsrtpa.')

                go to 990
            end if

            read (ninpts,1610,err=990) vreac(nrc)
1610 format(12x,e12.5)

            nci = 0

            do n = 1,nctmax + 1
                read (ninpts,1630,err=990) uxe,xx
1630 format(3x,a8,3x,e22.15)

                if (uxe(1:8) .eq. uendit(1:8)) then
                    go to 160
                end if

                nci = nci + 1

                if (nci .gt. nctmax) then
                    write (nttyo,1660) ureac(nrc)(1:j2),nctmax
1660 format(/' * Error - (XCON6/rd6w8) Have too many',' chemical',/7x,'elements in the special reactant ',a,'.',/7x,'The code is only dimensioned for ',i4,' elements.',/7x,'Reduce the number of elements or',' increase the',/7x,'dimensioning parameter nctpar.')

                    go to 990
                end if

                uesri(nci,nsrt) = uxe
                cesri(nci,nsrt) = xx
            end do

160 continue
            iesrti(nsrt) = nci

            nbi = 0

            do n = 1,nbt1mx + 1
                read (ninpts,1670,err=990) uxs,xx
1670 format(3x,a24,3x,e22.15)

                if (uxs(1:8) .eq. uendit(1:8)) then
                    go to 170
                end if

                nbi = nbi + 1

                if (nbi .gt. nbt1mx) then
                    write (nttyo,1700) ureac(nrc)(1:j2),nbtmax
1700 format(/' * Error - (XCON6/rd6w8) Have too many basis',' basis species in the',/7x,'in the reaction for the',' special reactant ',a,'.',/7x,'The code is only',' dimensioned for ',i4,' basis species. Increase',/7x,'the dimensioning parameter nbtpar.')

                    go to 990
                end if

                ubsri(nbi,nsrt) = uxs
                cbsri(nbi,nsrt) = xx
            end do

170 continue
            ibsrti(nsrt) = nbi
        else if (jcode(nrc) .eq. 5) then
            ! Have a generic ion exchanger reactant.
            nert = nert + 1
            ner = nert

            if (nert .gt. nertmx) then
                write (nttyo,1710) nertmx
1710 format(/' * Error - (XCON6/rd6w8) Have too many generic',' ion exchanger',/7x,'reactants. The code is only',' dimensioned for ',i4,' such reactants.',/7x,'Reduce',' the number of such reactants or increase the',/7x,'dimensioning parameter nertpa.')

                go to 990
            end if

            read (ninpts,1712,err=990) ugermo(ner)
1712 format(12x,a24)

            jeti = 0

            do jj = 1,jetmax + 1
                read (ninpts,1716,err=990) uxs
1716 format(3x,a24)

                if (uxs(1:8) .eq. uendit(1:8)) then
                    go to 154
                end if

                jeti = jeti + 1
                jei = jeti

                if (jeti .gt. jetmax) then
                    write (nttyo,1718) ureac(nrc)(1:j2),jetmax
1718 format(/' * Error - (XCON6/rd6w8) Have too many',' exchange sites',/7x,'in the generic ion exchanger',' reactant ',a,'.',/7x,'The code is only dimensioned',' for ',i4,' exchange sites per',/7x,'generic ion',' exchanger. Reduce the number of exchange sites or',/7x,'increase the dimensioning parameter jetpar.')

                    go to 990
                end if

                ugerji(jei,ner) = uxs

                neti = 0

                do nn = 1,netmax + 1
                    read (ninpts,1720,err=990) uxs,xx
1720 format(6x,a24,3x,e12.5)

                    if (uxs(1:8) .eq. uendit(1:8)) then
                        go to 152
                    end if

                    neti = neti + 1
                    nei = neti

                    if (neti .gt. netmax) then
                        j3 = ilnobl(ugerji(jei,nert))
                        write (nttyo,1724) ureac(nrc)(1:j2),ugerji(jei,ner)(1:j3),netmax
1724 format(/' * Error - (XCON6/rd6w8) Have too many',' species on exchange',/7x,'site ',a,' of the generic',' ion exchanger reactant',/7x,a,'. The code is only',' dimensioned for species per',/7x,'exchange site.',' Reduce the number of end-members or increase',/7x,'the dimensioning parameter netpar.')

                        go to 990
                    end if

                    ugersi(nei,jei,ner) = uxs
                    egersi(nei,jei,ner) = xx
                    xgersi(nei,jei,ner) = xx
                end do

152 continue
                igerti(jei,ner) = neti
            end do

154 continue
            jgerti(ner) = jeti
        end if

        ! Surface area parameters.
        read (ninpts,1740,err=990) nsk(nrc),sfcar(nrc),ssfcar(nrc),fkrc(nrc)
1740 format(12x,i2,22x,e12.5,12x,e12.5,/12x,e12.5)

        ! Rate law codes: nrk(1,nrc) = forward rate law code,
        ! nrk(2,nrc) = backward rate law code.
        read (ninpts,1760,err=990) nrk(1,nrc),nrk(2,nrc)
1760 format(12x,i2,22x,i2)

        ! Read forward (dissolution, dissociation) rate data.
        if (nrk(1,nrc) .eq. -1) then
            ! Use the backward rate law.
            continue
        else if (nrk(1,nrc) .eq. 0) then
            ! Illegal value.
            write (nttyo,1775) ureac(nrc)(1:j2)
1775 format(/' * Error - (XCON6/rd6w8) The forward rate law code',/7x,'has an illegal value of 0 for ',a,'.')

            go to 990
        else if (nrk(1,nrc) .eq. 1) then
            ! Arbitrary kinetics (relative rates, indifferent to time).
            imech(1,nrc) = 3

            if (imech(1,nrc) .gt. imchmx) then
                write (nttyo,1780) ureac(nrc)(1:j2),imchmx
1780 format(/' * Error - (XCON6/rd6w8) Have too many rate',/7x,'constants or corresponding mechanisms in the forward',/7x,'rate law for reactant ',a,'. The code is only',/7x,'dimensioned for ',i2,' rate constants per rate law.',/7x,'Reduce the number of rate constants or increase the',/7x,'dimensioning parameter imchpa.')

                go to 990
            end if

            read (ninpts,1790,err=990) (rkb(i,1,nrc), i = 1,3)
1790 format(3(12x,e12.5))
        else if (nrk(1,nrc) .eq. 2) then
            ! Transition state rate law. Up to imchmx parallel mechanisms
            ! are allowed.
            read (ninpts,1810,err=990) imech(1,nrc)
1810 format(12x,i2)

            if (imech(1,nrc) .gt. imchmx) then
                write (nttyo,1780) ureac(nrc)(1:j2),imchmx
                go to 990
            end if

            do i = 1,imech(1,nrc)
                read (ninpts,1830,err=990) rkb(i,1,nrc),trkb(i,1,nrc),iact(i,1,nrc)
1830 format(12x,e12.5,12x,e12.5,12x,i2)

                read (ninpts,1850,err=990) eact(i,1,nrc),hact(i,1,nrc)
1850 format(12x,e12.5,12x,e12.5)

                read (ninpts,1870,err=990) ndact(i,1,nrc),csigma(i,1,nrc)
1870 format(12x,i2,22x,e12.5)

                if (ndact(i,1,nrc) .gt. ndctmx) then
                    write (nttyo,1890) i,ureac(nrc)(1:j2),ndctmx
1890 format(/' * Error - (XCON6/rd6w8) Have too many',/7x,'species in the activity product in term ',i2,/7x,'of the forward direction rate law for reactant',/7x,a,'. The code is only dimensioned for ',i3,/7x,'such species. Reduce the number of such species',/7x,'or increase the dimensioning parameter ndctpa.')

                    go to 990
                end if

                do n = 1,ndact(i,1,nrc)
                    read (ninpts,1910,err=990) udac(n,i,1,nrc),cdac(n,i,1,nrc)
1910 format(12x,a24,12x,e12.5)
                end do
            end do
        else if (nrk(1,nrc) .eq. 3) then
            ! Linear rate law.
            imech(1,nrc) = 1

            if (imech(1,nrc) .gt. imchmx) then
                write (nttyo,1780) ureac(nrc)(1:j2),imchmx
                go to 990
            end if

            i = 1
            read (ninpts,1830,err=990) rkb(i,1,nrc),trkb(i,1,nrc),iact(i,1,nrc)
            read (ninpts,1850,err=990) eact(i,1,nrc),hact(i,1,nrc)
        else
            write (nttyo,1950) nrk(1,nrc),ureac(nrc)(1:j2)
1950 format(/' * Error - (XCON6/rd6w8) The forward rate law code',/7x,'has an unrecognized value of ',i2,' for ',a,'.')

            go to 990
        end if

        ! Read backward (precipitation, formation) rate data.
        if (nrk(2,nrc) .eq. -1) then
            ! Use the forward rate law.
            continue
        else if (nrk(2,nrc) .eq. 0) then
            ! Use instantaneous partial equilibrium.
            continue
        else if (nrk(2,nrc) .eq. 1) then
            ! Arbitrary kinetics.
            imech(2,nrc) = 3

            if (imech(2,nrc) .gt. imchmx) then
                write (nttyo,1960) ureac(nrc)(1:j2),imchmx
1960 format(/' * Error - (XCON6/rd6w8) Have too many rate',/7x,'constants or corresponding mechanisms in the backward',/7x,'rate law for reactant ',a,'. The code is only',/7x,'dimensioned for ',i2,' rate constants per rate law.',/7x,'Reduce the number of rate constants or increase the',/7x,'dimensioning parameter imchpa.')

                go to 990
            end if

            read (ninpts,1790,err=990) (rkb(i,2,nrc), i = 1,3)
        else if (nrk(2,nrc) .eq. 2) then
            ! Transition state rate law.
            read (ninpts,1810,err=990) imech(2,nrc)

            if (imech(2,nrc) .gt. imchmx) then
                write (nttyo,1960) ureac(nrc)(1:j2),imchmx
                go to 990
            end if

            do i = 1,imech(2,nrc)
                read (ninpts,1830,err=990) rkb(i,2,nrc),trkb(i,2,nrc),iact(i,2,nrc)

                read (ninpts,1850,err=990) eact(i,2,nrc),hact(i,2,nrc)

                read (ninpts,1870,err=990) ndact(i,2,nrc),csigma(i,2,nrc)

                if (ndact(i,2,nrc) .gt. ndctmx) then
                    write (nttyo,1970) i,ureac(nrc)(1:j2),ndctmx
1970 format(/' * Error - (XCON6/rd6w8) Have too many',/7x,'species in the activity product in term ',i2,/7x,'of the backward direction rate law for reactant',/7x,a,'. The code is only dimensioned for ',i3,/7x,'such species. Reduce the number of such species',/7x,'or increase the dimensioning parameter ndctpa.')

                    go to 990
                end if

                do n = 1,ndact(i,2,nrc)
                    read (ninpts,1910,err=990) udac(n,i,2,nrc),cdac(n,i,2,nrc)
                end do
            end do
        else if (nrk(2,nrc) .eq. 3) then
            ! Linear rate law.
            imech(2,nrc) = 1

            if (imech(2,nrc) .gt. imchmx) then
                write (nttyo,1960) ureac(nrc)(1:j2),imchmx
                go to 990
            end if

            i = 1
            read (ninpts,1830,err=990) rkb(i,2,nrc),trkb(i,2,nrc),iact(i,2,nrc)
            read (ninpts,1850,err=990) eact(i,2,nrc),hact(i,2,nrc)
        else
            write (nttyo,1980) nrk(2,nrc),ureac(nrc)(1:j2)
1980 format(/' * Error - (XCON6/rd6w8) The backward rate law code',/7x,'has an unrecognized value of ',i2,' for ',a,'.')

            go to 990
        end if
    end do

    ! Starting, minimum, and maximum values of key run parameters.
    read (ninpts,1150,err=990) xistti,ximaxi,tistti,timmxi
1150 format(2(12x,e12.5),/2(12x,e12.5))

    read (ninpts,1160,err=990) phmini,phmaxi,ehmini,ehmaxi,o2mini,o2maxi,awmini,awmaxi,kstpmx
1160 format(2(12x,e12.5),/2(12x,e12.5),/2(12x,e12.5),/2(12x,e12.5),/12x,i12)

    ! Print interval parameters.
    read (ninpts,1170,err=990) dlxprn,dlxprl,dltprn,dltprl
1170 format(2(12x,e12.5),12x,/2(12x,e12.5))

    read (ninpts,1180,err=990) dlhprn,dleprn,dloprn,dlaprn,ksppmx
1180 format(2(12x,e12.5),12x,/2(12x,e12.5),/12x,i12)

    ! Plot interval parameters.
    read (ninpts,1190,err=990) dlxplo,dlxpll,dltplo,dltpll
1190 format(2(12x,e12.5),/2(12x,e12.5))

    read (ninpts,1195,err=990) dlhplo,dleplo,dloplo,dlaplo,ksplmx
1195 format(2(12x,e12.5),/2(12x,e12.5),/12x,i12)

    ! Iopt model option switches.
    ! Note: iopt(1) = iopt1, etc.
    read (ninpts,1210,err=990) (iopt(i), i = 1,20)
1210 format(12x,10i5)

    ! Iopr print option switches.
    ! Note: iopr(1) = iopr1, etc.
    read (ninpts,1230,err=990) (iopr(i), i = 1,20)
1230 format(12x,10i5)

    ! Iodb debugging print option switches.
    ! Note: iodb(1) = iodb1, etc.
    read (ninpts,1250,err=990) (iodb(i), i = 1,20)
1250 format(12x,10i5)

    ! Number of nxopt options.
    read (ninpts,1300,err=990) nxopt
1300 format(12x,i2)

    if (nxopt .gt. nxopmx) then
        write (nttyo,1320) nxopmx
1320 format(/' * Error - (XCON6/rd6w8) Have too many mineral',/7x,'subset-selection suppression options. The code is',/7x,'only dimensioned for ',i3,' such options. Reduce the',/7x,'number of options or increase the dimensioning',/7x,'parameter nxoppa.')

        go to 990
    end if

    ! Nxopt options.
    do n = 1,nxopt
        read (ninpts,1330,err=990) uxopt(n),uxcat(n)
1330 format(12x,a6,1x,a24)
    end do

    if (nxopt .gt. 0) then
        read (ninpts,1350,err=990) nxopex
1350 format(12x,i2)

        if (nxopex .gt. nxpemx) then
            write (nttyo,1370) nxpemx
1370 format(/' * Error - (XCON6/rd6w8) Have too many',/7x,'exceptions specified to the mineral subset-selection',/7x,'suppression options. The code is only dimensioned',/7x,'for ',i3,'exceptions. Reduce the number of exceptions',/7x,'or increase the dimensioning parameter nxpepa.')

            go to 990
        end if

        do n = 1,nxopex
            read (ninpts,1380,err=990) uxopex(n)
1380 format(12x,a24)
        end do
    end if

    ! Number of nffg options.
    read (ninpts,1400,err=990) nffg
1400 format(12x,i2)

    if (nffg .gt. nffgmx) then
        write (nttyo,1420) nffgmx
1420 format(/' * Error - (XCON6/rd6w8) Have too many gases whose',/7x,'fugacities are to be fixed. The code is only dimensioned',/7x,'for ',i4,' such gases. Reduce the number of gases or',/7x,'increase the dimensioning parameter nffgpa.')

        go to 990
    end if

    ! Nffg options.
    do n = 1,nffg
        read (ninpts,1430,err=990) uffg(n),moffg(n),xlkffg(n)
1430 format(12x,a24,/12x,e12.5,12x,e12.5)
    end do

    ! Maximum finite-difference order.
    read (ninpts,2010,err=990) nordmx
2010 format(12x,i3)

    ! Newton-Raphson convergence tolerances.
    read (ninpts,2020,err=990) tolbt,toldl
2020 format(2(12x,e12.5))

    ! Maximum number of Newton-Raphson iterations.
    read (ninpts,2010,err=990) itermx

    ! Search/find tolerance.
    read (ninpts,2022,err=990) tolxsf
2022 format(12x,e12.5)

    ! Saturation tolerance.
    read (ninpts,2024,err=990) tolsat
2024 format(12x,e12.5)

    ! Maximum number of phase assemblage tries.
    read (ninpts,2102,err=990) ntrymx
2102 format(12x,i3)

    ! Zero-order step size (in Xi).
    read (ninpts,2080,err=990) dlxmx0
2080 format(12x,e12.5)

    ! PRS transfer interval in Xi.
    read (ninpts,2000,err=990) dlxdmp
2000 format(12x,e12.5)

    ! Process the bottom half of the current input file.
    ! Secondary title.
    ntitl2 = 0
    read (ninpts,1000,err=990) uline

    if (uline(1:8) .ne. uendit(1:8)) then
        ntitl2 = 1
        utitl2(1) = uline

        do n = 2,ntitmx
            read (ninpts,1000,err=990) uline

            if (uline(1:8) .eq. uendit(1:8)) then
                go to 220
            end if

            utitl2(n) = uline
            ntitl2 = n
        end do

220 continue
    end if

    ! Special basis switches.
    read (ninpts,2620,err=990) nsbswt
2620 format(12x,i3)

    do n = 1,nsbswt
        read (ninpts,2640,err=990) usbsw(1,n)
2640 format(9x,a48)

        read (ninpts,2660,err=990) usbsw(2,n)
2660 format(15x,a48)
    end do

    ! Original temperature.
    read (ninpts,2680,err=990) tempci
2680 format(12x,e12.5)

    ! Original pressure.
    read (ninpts,2700,err=990) pressi
2700 format(12x,e12.5)

    ! Ion exchanger creation. This section reads the directives for
    ! creating ion exchange phases and species and their associated
    ! intrinsic (including thermodynamic) properties. The data is
    ! essentially that which might be read from a supporting data file.
    read (ninpts,2795) ux8
2795 format(12x,a8)

    call lejust(ux8)
    call locase(ux8)
    qgexsh = ux8(1:2).eq.'.t' .or. ux8(1:1).eq. 't'

    read (ninpts,2800,err=990) net
2800 format(12x,i3)

    if (net .gt. netmax) then
        write (nttyo,2820) netmax,net
2820 format(/' * Error - (XCON6/rd6w8) Have exceeded the maximum',' number of ',i3,/7x,'generic ion exchange phases while',' reading the data to create',/7x,'such phases. Increase',' the dimensioning parameter netpar',/7x,'to at least ',i3,'.')

        go to 990
    end if

    do ne = 1,net
        read (ninpts,2830,err=990) ugexp(ne)
2830 format(12x,a24)

        j2 = ilnobl(ugexp(ne))

        read (ninpts,2850,err=990) mwtges(ne)
2850 format(12x,e12.5)

        read (ninpts,2830,err=990) ugexmo(ne)
        j3 = ilnobl(ugexmo(ne))

        read (ninpts,2850,err=990) tgexp(ne)

        read (ninpts,2890,err=990) jgext(ne)
2890 format(12x,i3)

        if (jgext(ne) .gt. jetmax) then
            write (nttyo,2910) jetmax,ugexp(ne)(1:j2),jgext(ne)
2910 format(/' * Error - (XCON6/rd6w8) Have exceeded the maximum',' number of ',i3,/7x,'exchange sites on a generic ion',' exchange phase while reading',/7x,'the data to create',a,'. Increase the',/7x,'dimensioning parameter jetpar to',' at least ',i3,'.')

            go to 990
        end if

        do je = 1,jgext(ne)
            read (ninpts,2930,err=990) ugexj(je,ne)
2930 format(12x,a8)

            j3 = ilnobl(ugexj(je,ne))

            read (ninpts,2950,err=990) cgexj(je,ne),zgexj(je,ne)
2950 format(12x,e12.5,12x,e12.5)

            read (ninpts,2890,err=990) ngexrt(je,ne)

            if (ngexrt(je,ne) .gt. ietmax) then
                write (nttyo,3020) netmax,ugexj(je,ne)(1:j3),ugexp(ne)(1:j2),ngexrt(je,ne)
3020 format(/' * Error - (XCON6/rd6w8) Have exceeded the',' maximum number of ',i3,/7x,'reactions for a site',' belonging to a generic ion exchange',/7x,'phase while',' reading the data for site ',a,' of exchange phase'      /7x,a,'. Increase the dimensioning parameter',/7x,'ietpar to at least ',i3,'.')

                go to 990
            end if

            do n = 1,ngexrt(je,ne)
                read (ninpts,3030,err=990) ugexr(n,je,ne)
3030 format(12x,a56)

                read (ninpts,3050,err=990) xlkgex(n,je,ne),uxkgex(n,je,ne)
3050 format(12x,e12.5,12x,a8)

                read (ninpts,3050,err=990) xhfgex(n,je,ne),uhfgex(n,je,ne)
                read (ninpts,3050,err=990) xvfgex(n,je,ne),uvfgex(n,je,ne)
            end do
        end do
    end do

    ! Number of nxmod options.
    read (ninpts,3120,err=990) nxmod
3120 format(12x,i2)

    if (nxmod .gt. nxmdmx) then
        write (nttyo,3140) nxmdmx
3140 format(/' * Error - (XCON6/rd6w8) Have too many nxmod',/7x,'alter/suppress options. The code is only dimensioned',/7x,'for ',i2,' such options. Reduce the number of such',/7x,'options or increase the dimensioning parameter nxmdpa.')

        go to 990
    end if

    ! Nxmod options.
    do n = 1,nxmod
        read (ninpts,3150,err=990) uxmod(n),kxmod(n),xlkmod(n)
3150 format(12x,a48,/12x,i2,22x,e12.5)
    end do

    ! Iopg options.
    ! Note: iopg(1) = iopg1, etc.
    read (ninpts,3170,err=990) (iopg(i), i = 1,20)
3170 format(12x,10i5)

    ! Index limits.
    read (ninpts,3220,err=990) kct,kbt,kmt,kxt,kdim,kprs
3220 format(3(12x,i2,10x),/3(12x,i2,10x))

    nbti = kbt

    if (kct .gt. nctmax) then
        write (nttyo,3240) nctmax
3240 format(/' * Error - (XCON6/rd6w8) Have too many chemical',/7x,'elements present. The code is only dimensioned',/7x,'for ',i3,' elements. Reduce the number of elements',/7x,'or increase the dimensioning parameter nctpar.')

        go to 990
    end if

    if (kbt .gt. nbtmax) then
        write (nttyo,3250) nbtmax
3250 format(/' * Error - (XCON6/rd6w8) Have too many basis',/7x,'species present. The code is only dimensioned',/7x,'for ',i3,' basis species. Reduce the number of elements',/7x,'or increase the dimensioning parameter nctpar.')

        go to 990
    end if

    if (kdim .gt. kmax) then
        write (nttyo,3260) kmax
3260 format(/' * Error - (XCON6/rd6w8) Have too many matrix',/7x,'variables. The code is only dimensioned for ',i3,/7x,'matrix variables. Reduce the number of such variables',/7x,'or increase the dimensioning parameter kpar.')

        go to 990
    end if

    ! Species for which mass balances are defined.
    !   ubmtbi = names of the data file basis species for which mass
    !              balances are defined
    !   jflgi  = jflag input for basis species
    !      0 = Retain as an active basis species
    !     30 = Convert to a dependent species; fold the mass balance
    !            total for this species into the mass balance totals
    !            of basis species which remain active
    do nbi = 1,nbti
        read (ninpts,3340,err=990) ubmtbi(nbi),jflgi(nbi)
3340 format(3x,a48,3x,i2)
    end do

    ! Mass balance totals.
    !   mtbi   = total number of moles of basis species in the
    !              equilibrium system (ES)
    !   mtbaqi = total number of moles of basis species in the
    !              aqueous solution
    do nbi = 1,nbti
        read (ninpts,3410,err=990) mtbi(nbi),mtbaqi(nbi)
3410 format(6x,e22.15,6x,e22.15)
    end do

    read (ninpts,3430,err=990) electr
3430 format(34x,e22.15)

    ! Ordinary basis switches.
    read (ninpts,2620,err=990) nobswt

    do n = 1,nobswt
        read (ninpts,2640,err=990) uobsw(1,n)

        read (ninpts,2660,err=990) uobsw(2,n)
    end do

    ! Matrix column variables.
    !   uzveci = names of species or properties associated with the
    !              matrix column variables
    do krow = 1,kdim
        read (ninpts,3480,err=990) uzveci(krow)
3480 format(3x,a48)
    end do

    ! Matrix variable values.
    !   zvclg1 = values of matrix variables
    do kcol = 1,kdim
        read (ninpts,3520,err=990) zvclgi(kcol)
3520 format(3x,e22.15)
    end do

    ! Initialize any non-zero values for the physically removed
    ! system. Here mprphi and mprspi are arrays for the number of
    ! moles of phases and species, respectively, in the PRS.
    nprpti = 0
    nprsti = 0

    if (kprs .gt. 0) then
        do npi = 1,nprpmx + 1
            read (ninpts,3610,err=990) uxf,mx
3610 format(1x,a24,6x,e22.15)

            if (uxf(1:8) .eq. uendit(1:8)) then
                go to 310
            end if

            nprpti = nprpti + 1

            if (nprpti .gt. nprpmx) then
                write (nttyo,3640) nprpmx
3640 format(/' * Error - (XCON6/rd6w8) Have too many phases',/7x,'in the physically removed system (PRS) on the',/7x,'input file. The code is only dimensioned for ',i3,/7x,'such phases. Reduce the number of such phases',/7x,'or increase the dimensioning parameter nprppa.')

                go to 990
            end if

            uprphi(nprpti)(1:24) = uxf
            mprphi(nprpti) = mx

            do iki = 1,iktmax + 1
                read (ninpts,3650,err=990) uxg,mx
3650 format(3x,a24,6x,e22.15)

                if (uxg(1:8) .eq. uendit(1:8)) then
                    go to 300
                end if

                if (iki .gt. iktmax) then
                    j2 = ilnobl(uxf)
                    write (nttyo,3680) uxf(1:j2),iktmax
3680 format(/' * Error - (XCON6/rd6w8) Have too many',' end-members',/7x,'in the PRS solid solution ',a,'.',/7x,'The code is only dimensioned for ',i4,' end-members per',/7x,'solid solution. Reduce',' the number of end-members or',/7x,'increase the dimensioning parameter iktpar.')

                    go to 990
                end if

                nprsti = nprsti + 1

                if (nprsti .gt. nprsmx) then
                    write (nttyo,3690) nprsmx
3690 format(/' * Error - (XCON6/rd6w8) Have too many species',/7x,'in the physically removed system (PRS) on the',/7x,'input file. The code is only dimensioned for ',i3,/7x,'such species. Reduce the number of such species',/7x,'or increase the dimensioning parameter nprspa.')

                    go to 990
                end if

                uprspi(nprsti)(1:24) = uxg
                uprspi(nprsti)(25:48) = uxf
                mprspi(nprsti) = mx
            end do

300 continue
        end do

310 continue
    end if

    go to 999

990 continue
    qrderr = .true.

999 continue
end subroutine rd6w8