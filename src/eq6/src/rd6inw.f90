subroutine rd6inw(awmaxi,awmini,cbsri,cdac,cesri,cgexj,csigma,dlaplo,dlaprn,dleplo,dleprn,dlhplo,dlhprn,dloplo,dloprn,dltpll,dltplo,dltprl,dltprn,dlxdmp,dlxmx0,dlxpll,dlxplo,dlxprl,dlxprn,eact,egersi,ehmaxi,ehmini,electr,fkrc,iact,ibsrti,iesrti,ietmax,iktmax,imchmx,imech,iodb,iopg,iopr,iopt,igerti,itermx,ixrti,jcode,jgerti,jetmax,jflgi,jgext,jpress,jreac,jtemp,kbt,kct,kdim,kmax,kmt,kprs,ksplmx,ksppmx,kstpmx,kxmod,kxt,hact,modr,moffg,morr,mprphi,mprspi,mtbaqi,mtbi,mwtges,nbti,nbtmax,nbt1mx,nctmax,ndact,ndctmx,nert,nertmx,net,netmax,nffg,nffgmx,ngexrt,ninpts,nobswt,nodbmx,nopgmx,noprmx,noptmx,nordmx,noutpt,nprob,nprpmx,nprpti,nprsmx,nprsti,nptkmx,nrct,nrctmx,nrk,nsbswt,nsk,nsrt,nsrtmx,ntitl1,ntitl2,ntitmx,ntrymx,nttkmx,nttyo,nxmdmx,nxmod,nxopex,nxopmx,nxopt,nxpemx,nxrt,nxrtmx,o2maxi,o2mini,phmaxi,phmini,pressb,pressi,ptk,qend,qgexsh,qrderr,rkb,rxbari,sfcar,ssfcar,tempcb,tempci,tgexp,timmxi,tistti,tolbt,toldl,tolsat,tolxsf,trkb,ttk,ubmtbi,ubsri,ucxri,udac,uesri,uffg,ugerji,ugermo,ugersi,ugexj,ugexmo,ugexp,ugexr,uhfgex,uobsw,uprphi,uprspi,ureac,usbsw,utitl1,utitl2,uvfgex,uxcat,uxkgex,uxmod,uxopex,uxopt,uzveci,vreac,xgersi,xhfgex,ximaxi,xistti,xlkffg,xlkgex,xlkmod,xvfgex,zgexj,zvclgi)
    !! This subroutine reads the EQ6 input file in compact ("W") format.
    !! This subroutine is a near-clone of XCON6/rd6w8.f. The present
    !! subroutine differs from the latter, in that in addition to the
    !! pure read function, it writes an instant echo of what is read to
    !! the output file, and sandwiches this between prefatory and
    !! ending messages. To maintain close consistency with XCON6/rd6w8.f,
    !! this subroutine assigns no default values and only performs such
    !! checking of what is read to ensure that what follows is readable.
    !! The calling sequence of this subroutine is identical to that of
    !! EQ6/rd6ind.f, XCON6/rd6w8.f, and XCON6/rd6d8.f.
    !! This subroutine is called by:
    !!   EQ6/eq6.f
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
    integer :: jflgi(nbtmax)
    integer :: jgerti(nertmx)
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
    integer :: j4
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

    ! qend   = .true if the end of the input file has been encountered
    ! qrderr = .true if the current problem can't be read because of
    !            a read format error or a dimensional overflow
    qend = .false.
    qrderr = .false.

    ! Main title.
    !   utitl1 = main title
    !   ntitl1 = number of lines in the main title
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
    write (noutpt,1010) nprob
    write (nttyo,1010) nprob
1010 format(//' Reading problem ',i3,' from the input file ...',/)

    j2 = ilnobl(uline)
    j2 = min(j2,79)
    write (noutpt,1020) uline
1020 format(1x,a)

    if (uline(1:8) .ne. uendit(1:8)) then
        ntitl1 = 1
        utitl1(1) = uline
    end if

    if (ntitl1 .gt. 0) then
        do n = 2,ntitmx
            read (ninpts,1000,err=990) uline
            j2 = ilnobl(uline)
            j2 = min(j2,79)
            write (noutpt,1020) uline(1:j2)

            if (uline(1:8) .eq. uendit(1:8)) then
                go to 120
            end if

            ntitl1 = n
            utitl1(n) = uline
        end do

120 continue

        ! Write the first 5 lines of the input problem title to the
        ! screen file.
        write (nttyo,1040)
1040 format(3x,'The input problem title is (first 5 lines',' maximum):',/)

        i = min(ntitl1,5)

        do n = 1,i
            j2 = ilnobl(utitl1(n))
            j2 = min(j2,74)
            write (nttyo,1050) utitl1(n)(1:j2)
1050 format(5x,a)
        end do
    end if

    write (nttyo,1060)
1060 format(/3x,'Continuing to read the problem input ...')

    ! Temperature parameters.
    ! Note: ttk(1) = ttk1, etc.
    !   jtemp  = temperature tracking code:
    !     0 = Constant temperature:
    !           tempc = tempcb
    !     1 = Linear tracking in Xi:
    !           tempc = tempcb + ttk(1)*xi1
    !     2 = Linear tracking in time (iopt(2) must be 1):
    !           tempc = tempcb + ttk(1)*time1
    !     3 = Fluid mixing tracking:
    !           tempc= (tempcb*ttk(1) + xi*ttk(2))/(xi + ttk(1)), where:
    !             ttk(1) = ratio of the mass of the starting aqueous
    !               to that of the aqueous solution being treated as a
    !               reactant; usually the mass of each will be close to
    !               1 kilogram, hence the value of ttk(1) will then be
    !               about 1.0.
    !             ttk(2) = temperature of the fluid being treated as a
    !               reactant ("Fluid 2")
    !   tempcb = the base temperature, C
    !   ttk    = temperature tracking coefficients
    read (ninpts,1100,err=990) jtemp,tempcb,ttk(1),ttk(2)
1100 format(12x,i2,/12x,e12.5,/2(12x,e12.5))

    write (noutpt,1110) jtemp,tempcb,ttk(1),ttk(2)
1110 format(6x,'jtemp= ',i2,/5x,'tempcb= ',1pe12.5,/7x,'ttk1= ',e12.5,6x,'ttk2= ',e12.5)

    ! Pressure parameters.
    ! Note: ptk(1) = ptk1, etc.
    !   jpress = pressure tracking code:
    !     0 = Follow the data file's reference pressure curve:
    !           press = func(data file, T)
    !     1 = Follow the 1.013-bar/steam-saturation curve:
    !           press = built-in func(T)
    !     2 = Constant pressure:
    !           press = pressb
    !     3 = Linear tracking in Xi:
    !           press = pressb + ptk(1)*xi1
    !     4 = Linear tracking in time (iopt(2) must be 1):
    !           press = pressb + ptk(1)*time1
    !   pressb = the base pressure, bars
    !   ptk    = pressure tracking coefficients
    read (ninpts,1120,err=990) jpress,pressb,ptk(1),ptk(2)
1120 format(12x,i2,/12x,e12.5,/2(12x,e12.5))

    write (noutpt,1130) jpress,pressb,ptk(1),ptk(2)
1130 format(4x,'jpress= ',i2,/5x,'pressb= ',1pe12.5,/7x,'ptk1= ',e12.5,6x,'ptk2= ',e12.5)

    ! Number of reactants.
    !   nrct   = number of reactants
    read (ninpts,1480,err=990) nrct
1480 format(12x,i2)

    write (noutpt,1490) nrct
1490 format(7x,'nrct= ',i2)

    write (noutpt,1500)
1500 format(' *-------------------------------------------------','----------------------------')

    if (nrct .gt. nrctmx) then
        write (noutpt,1520) nrctmx
        write (nttyo,1520) nrctmx
1520 format(/' * Error - (EQ6/rd6inw) Have too many reactants',/7x,'The code is only dimensioned for ',i4,' reactants.',/7x,'Reduce the number of reactants or increase the',/7x,'dimensioning parameter nrctpa.')

        go to 990
    end if

    ! Reactants.
    !   nxrt  = number of solid solution reactants
    !   nsrt  = number of special reactants
    !   nert  = number of ion exchanger reactants
    !   ureac = name of reactant
    !   morr  = moles of reactant remaining when Xi = xistti
    !             (the starting number of moles for the process
    !             corresponds to Xi = 0)
    !   modr  = moles of reactant irreversibly destroyed when
    !             Xi = xistti (usually zero at Xi = 0)
    !   jcode = reactant type flag:
    !             0 = mineral
    !             1 = solid solution
    !             2 = special reactant
    !             3 = aqueous species
    !             4 = gas
    !   jreac = reactant status flag:
    !              0 = set to react (this is the only value a user
    !                    should set; the code will write other values
    !                    on pickup files)
    !             -1 = saturated, but the remaining reactant mass
    !                    continues to react irreversibly; there is
    !                    usually also a secondary product mass of the
    !                    same species, so that the net rate of reaction
    !                    is zero and the solution does not saturate
    !              1 = exhausted
    !              2 = saturated; the status of any remaining reactant
    !                    mass is changed to that of a product phase;
    !                    this only happens if iopt(1) = 0
    nxrt = 0
    nsrt = 0
    nert = 0

    do nrc = 1,nrct
        ! Name, flags, and masses.
        read (ninpts,1530,err=990) ureac(nrc),jcode(nrc),jreac(nrc),morr(nrc),modr(nrc)
1530 format(12x,a24,/12x,i2,22x,i2,/2(12x,e12.5))

        j2 = ilnobl(ureac(nrc))
        write (noutpt,1540) ureac(nrc)(1:j2),jcode(nrc),jreac(nrc),morr(nrc),modr(nrc)
1540 format(3x,'reactant= ',a,/6x,'jcode= ',i2,15x,'jreac= ',i2,/7x,'morr= ',1pe12.5,6x,'modr= ',e12.5)

        if (jcode(nrc) .eq. 1) then
            ! Have a solid solution reactant.
            !   ucxri  = name of end-member of a solid solution reactant
            !   rxbari = mole fraction of end-member of a solid solution
            !              reactant
            nxrt = nxrt + 1

            if (nxrt .gt. nxrtmx) then
                write (noutpt,1550) nxrtmx
                write (nttyo,1550) nxrtmx
1550 format(/' * Error - (EQ6/rd6inw) Have too many solid',/7x,'solution reactants. The code is only dimensioned',/7x,'for ',i4,' such reactants. Reduce the number of',/7x,'such reactants or increase the dimensioning',' parameter nxrtpa.')

                go to 990
            end if

            ikti = 0

            do i = 1,iktmax + 1
                read (ninpts,1560,err=990) uxs,xx
1560 format(3x,a24,3x,e12.5)

                if (uxs(1:8) .eq. uendit(1:8)) then
                    j3 = ilnobl(uxs)
                    write (noutpt,1570) uxs(1:j3)
1570 format(4x,a)

                    go to 150
                end if

                write (noutpt,1580) uxs,xx
1580 format(4x,a24,3x,1pe12.5)

                ikti = ikti + 1

                if (ikti .gt. iktmax) then
                    write (noutpt,1590) ureac(nrc)(1:j2),iktmax
                    write (nttyo,1590) ureac(nrc)(1:j2),iktmax
1590 format(/' * Error - (EQ6/rd6inw) Have too many',' end-members',/7x,'in the solid solution reactant ',a,'.',/7x,'The code is only dimensioned for ',i4,' end-members per',/7x,'solid solution. Reduce',' the number of end-members or',/7x,'increase the dimensioning parameter iktpar.')

                    go to 990
                end if

                ucxri(ikti,nxrt) = uxs
                rxbari(ikti,nxrt) = xx
            end do

150 continue
            ixrti(nxrt) = ikti
        else if (jcode(nrc) .eq. 2) then
            ! Have a special reactant.
            !   vreac  = volume, cc/mol
            !   uesri  = name of element
            !   cesri  =  moles chemical element per mole reactant
            !   ubsri  = name of basis species in the chemical reaction
            !   cbsri =  corresponding reaction coefficient
            nsrt = nsrt + 1

            if (nsrt .gt. nsrtmx) then
                write (noutpt,1600) nsrtmx
                write (nttyo,1600) nsrtmx
1600 format(/' * Error - (EQ6/rd6inw) Have too many special',/7x,'reactants. The code is only dimensioned for ',i4,/7x,'such reactants. Reduce the number of such reactants',/7x,'or increase the dimensioning parameter nsrtpa.')

                go to 990
            end if

            read (ninpts,1610,err=990) vreac(nrc)
1610 format(12x,e12.5)

            write (noutpt,1620) vreac(nrc)
1620 format(6x,'vreac= ',1pe12.5)

            write (noutpt,1622)
1622 format (' * Elemental composition')

            nci = 0

            do n = 1,nctmax + 1
                read (ninpts,1630,err=990) uxe,xx
1630 format(3x,a8,3x,e22.15)

                if (uxe(1:8) .eq. uendit(1:8)) then
                    j3 = ilnobl(uxe)
                    write (noutpt,1640) uxe(1:j3)
1640 format(4x,a)

                    go to 160
                end if

                write (noutpt,1650) uxe,xx
1650 format(4x,a8,3x,1pe22.15)

                nci = nci + 1

                if (nci .gt. nctmax) then
                    write (noutpt,1660) ureac(nrc)(1:j2),nctmax
                    write (nttyo,1660) ureac(nrc)(1:j2),nctmax
1660 format(/' * Error - (EQ6/rd6inw) Have too many',' chemical',/7x,'elements in the special reactant ',a,'.',/7x,'The code is only dimensioned for ',i4,' elements.',/7x,'Reduce the number of elements or',' increase the',/7x,'dimensioning parameter nctpar.')

                    go to 990
                end if

                uesri(nci,nsrt) = uxe
                cesri(nci,nsrt) = xx
            end do

160 continue
            iesrti(nsrt) = nci

            write (noutpt,1662)
1662 format (' * Reaction')

            ! Note: the special reactant itself should appear first in
            ! the associated chemical reaction. Ordinarily, its reaction
            ! coefficient should be -1.0.
            nbi = 0

            do n = 1,nbt1mx + 1
                read (ninpts,1670,err=990) uxs,xx
1670 format(3x,a24,3x,e22.15)

                if (uxs(1:8) .eq. uendit(1:8)) then
                    j3 = ilnobl(uxs)
                    write (noutpt,1680) uxs(1:j3)
1680 format(4x,a)

                    go to 170
                end if

                write (noutpt,1690) uxs,xx
1690 format(4x,a24,3x,1pe22.15)

                nbi = nbi + 1

                if (nbi .gt. nbt1mx) then
                    write (noutpt,1700) ureac(nrc)(1:j2),nbtmax
                    write (nttyo,1700) ureac(nrc)(1:j2),nbtmax
1700 format(/' * Error - (EQ6/rd6inw) Have too many basis',' basis species in the',/7x,'in the reaction for the',' special reactant ',a,'.',/7x,'The code is only',' dimensioned for ',i4,' basis species. Increase',/7x,'the dimensioning parameter nbtpar.')

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
                write (noutpt,1710) nertmx
                write (nttyo,1710) nertmx
1710 format(/' * Error - (EQ6/rd6inw) Have too many generic',' ion exchanger',/7x,'reactants. The code is only',' dimensioned for ',i4,' such reactants.',/7x,'Reduce',' the number of such reactants or increase the',/7x,'dimensioning parameter nertpa.')

                go to 990
            end if

            read (ninpts,1712,err=990) ugermo(ner)
1712 format(12x,a24)

            j3 = ilnobl(ugermo(ner))
            write (noutpt,1714) ugermo(ner)(1:j3)
1714 format(5x,'ugermo= ',a)

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
                    write (noutpt,1718) ureac(nrc)(1:j2),jetmax
                    write (nttyo,1718) ureac(nrc)(1:j2),jetmax
1718 format(/' * Error - (EQ6/rd6inw) Have too many',' exchange sites',/7x,'in the generic ion exchanger',' reactant ',a,'.',/7x,'The code is only dimensioned',' for ',i4,' exchange sites per',/7x,'generic ion',' exchanger. Reduce the number of exchange sites or',/7x,'increase the dimensioning parameter jetpar.')

                    go to 990
                end if

                ugerji(jei,ner) = uxs(1:8)

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
                        write (noutpt,1724) ureac(nrc)(1:j2),ugerji(jei,ner)(1:j3),netmax
                        write (nttyo,1724) ureac(nrc)(1:j2),ugerji(jei,ner)(1:j3),netmax
1724 format(/' * Error - (EQ6/rd6inw) Have too many',' species on exchange',/7x,'site ',a,' of the generic',' ion exchanger reactant',/7x,a,'. The code is only',' dimensioned for species per',/7x,'exchange site.',' Reduce the number of end-members or increase',/7x,'the dimensioning parameter netpar.')

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
        !   nsk    = surface area flag
        !     0 = fixed surface area
        !     1 = fixed specific surface area
        !     2 = n**2/3 growth law- current surface area
        !   sfcar  = surface area, cm2
        !   ssfcar = specific surface area, cm2/mol
        !   fkrc   = ratio of effective surface area to total surface
        !              area of a reactant; this parameter can also be
        !              used as a generalized correction factor for the
        !              rate constant
        read (ninpts,1740,err=990) nsk(nrc),sfcar(nrc),ssfcar(nrc),fkrc(nrc)
1740 format(12x,i2,22x,e12.5,12x,e12.5,/12x,e12.5)

        write (noutpt,1750) nsk(nrc),sfcar(nrc),fkrc(nrc)
1750 format(8x,'nsk= ',i2,15x,'sfcar= ',1pe12.5,4x,'ssfcar= ',e12.5,/7x,'fkrc= ',e12.5)

        ! Rate law codes: nrk(1,nrc) = forward rate law code,
        ! nrk(2,nrc) = backward rate law code.
        !   Forward rate law codes:
        !     -1 = Use the backward rate law form (legal only if
        !            nrk(2,nrc) = 2)
        !      0 = Illegal value
        !      1 = Relative rate
        !      2 = Transition state theory net rate
        !      3 = Linear rate
        !   For the case nrk(1,nrc) = 2, the reactant must be
        !   either a pure mineral (jcode = 0) or a solid solution
        !   (jcode = 1).
        !   Backward rate law codes:
        !     -1 = Use the forward rate law form (legal only if
        !            nrk(1,nrc) = 2)
        !      0 = No rate law specified; the reaction may be controlled
        !            by partial equilibrium
        !      1 = Relative rate
        !      2 = Transition state theory net rate
        !      3 = Linear rate
        !   For the case nrk(2,nrc) = 2, the reactant must be
        !   either a pure mineral (jcode = 0) or a solid solution
        !   (jcode = 1)
        read (ninpts,1760,err=990) nrk(1,nrc),nrk(2,nrc)
1760 format(12x,i2,22x,i2)

        write (noutpt,1770) nrk(1,nrc),nrk(2,nrc)
1770 format(7x,'nrk1= ',i2,16x,'nrk2= ',i2)

        ! Read forward (dissolution, dissociation) rate data.
        if (nrk(1,nrc) .eq. -1) then
            ! Use the backward rate law.
            continue
        else if (nrk(1,nrc) .eq. 0) then
            ! Illegal value.
            write (noutpt,1775) ureac(nrc)(1:j2)
            write (nttyo,1775) ureac(nrc)(1:j2)
1775 format(/' * Error - (EQ6/rd6inw) The forward rate law code',/7x,'has an illegal value of 0 for ',a,'.')

            go to 990
        else if (nrk(1,nrc) .eq. 1) then
            ! Arbitrary kinetics (relative rates, indifferent to time).
            !   relative rate = rkb(1,1,nrc) + rkb(2,1,nrc)*xi1
            !                     + (rkb(3,1,nrc)/2.)*xi1**2
            ! where
            !   rkb    = relative rate constants
            !   xi1    = the reaction progress variable
            imech(1,nrc) = 3

            if (imech(1,nrc) .gt. imchmx) then
                write (noutpt,1780) ureac(nrc)(1:j2),imchmx
                write (nttyo,1780) ureac(nrc)(1:j2),imchmx
1780 format(/' * Error - (EQ6/rd6inw) Have too many rate',/7x,'constants or corresponding mechanisms in the forward',/7x,'rate law for reactant ',a,'. The code is only',/7x,'dimensioned for ',i2,' rate constants per rate law.',/7x,'Reduce the number of rate constants or increase the',/7x,'dimensioning parameter imchpa.')

                go to 990
            end if

            read (ninpts,1790,err=990) (rkb(i,1,nrc), i = 1,3)
1790 format(3(12x,e12.5))

            write (noutpt,1800) (rkb(i,1,nrc), i = 1,3)
1800 format(7x,'rkb1= ',1pe12.5,6x,'rkb2= ',e12.5,6x,'rkb3= ',e12.5)
        else if (nrk(1,nrc) .eq. 2) then
            ! Transition state rate law. Up to imchmx parallel mechanisms
            ! are allowed.
            !   imech  = number of parallel transition state mechanisms
            !   rkb    = rate constants, one for each rate law term
            !   ndact  = number of aqueous species that appear in each
            !              rate law term
            !   udac   = the name of such an aqueous species
            !   cdac   = the power the activity of such a species is raised
            !              to in the rate law term
            !   csigma = ratio of affinity of a macroscopic reaction to the
            !              affinity of the microcopic reaction for
            !              destruction of an activated complex
            read (ninpts,1810,err=990) imech(1,nrc)
1810 format(12x,i2)

            write (noutpt,1820) imech(1,nrc)
1820 format(6x,'imech= ',i2)

            if (imech(1,nrc) .gt. imchmx) then
                write (noutpt,1780) ureac(nrc)(1:j2),imchmx
                write (nttyo,1780) ureac(nrc)(1:j2),imchmx
                go to 990
            end if

            do i = 1,imech(1,nrc)
                read (ninpts,1830,err=990) rkb(i,1,nrc),trkb(i,1,nrc),iact(i,1,nrc)
1830 format(12x,e12.5,12x,e12.5,12x,i2)

                write (noutpt,1840) rkb(i,1,nrc),trkb(i,1,nrc),iact(i,1,nrc)
1840 format(8x,'rkb= ',1pe12.5,6x,'trkb= ',e12.5,6x,'iact= ',i2)

                read (ninpts,1850,err=990) eact(i,1,nrc),hact(i,1,nrc)
1850 format(12x,e12.5,12x,e12.5)

                write (noutpt,1860) eact(i,1,nrc),hact(i,1,nrc)
1860 format(7x,'eact= ',1pe12.5,6x,'hact= ',e12.5)

                read (ninpts,1870,err=990) ndact(i,1,nrc),csigma(i,1,nrc)
1870 format(12x,i2,22x,e12.5)

                write (noutpt,1880) ndact(i,1,nrc),csigma(i,1,nrc)
1880 format(6x,'ndact= ',i2,14x,'csigma= ',1pe12.5)

                if (ndact(i,1,nrc) .gt. ndctmx) then
                    write (noutpt,1890) i,ureac(nrc)(1:j2),ndctmx
                    write (nttyo,1890) i,ureac(nrc)(1:j2),ndctmx
1890 format(/' * Error - (EQ6/rd6inw) Have too many',/7x,'species in the activity product in term ',i2,/7x,'of the forward direction rate law for reactant',/7x,a,'. The code is only dimensioned for ',i3,/7x,'such species. Reduce the number of such species',/7x,'or increase the dimensioning parameter ndctpa.')

                    go to 990
                end if

                do n = 1,ndact(i,1,nrc)
                    read (ninpts,1910,err=990) udac(n,i,1,nrc),cdac(n,i,1,nrc)
1910 format(12x,a24,12x,e12.5)

                    write (noutpt,1920) udac(n,i,1,nrc),cdac(n,i,1,nrc)
1920 format(7x,'udac= ',a24,6x,'cdac= ',1pe12.5)
                end do
            end do
        else if (nrk(1,nrc) .eq. 3) then
            ! Linear rate law.
            !   Rate = fkrc(nrc)*sfcar(nrc)*rkb(1,1,nrc)
            imech(1,nrc) = 1

            if (imech(1,nrc) .gt. imchmx) then
                write (noutpt,1780) ureac(nrc)(1:j2),imchmx
                write (nttyo,1780) ureac(nrc)(1:j2),imchmx
                go to 990
            end if

            i = 1
            read (ninpts,1830,err=990) rkb(i,1,nrc),trkb(i,1,nrc),iact(i,1,nrc)
            write (noutpt,1840) rkb(i,1,nrc),trkb(i,1,nrc),iact(i,1,nrc)
            read (ninpts,1850,err=990) eact(i,1,nrc),hact(i,1,nrc)
            write (noutpt,1860) eact(i,1,nrc),hact(i,1,nrc)
        else
            write (noutpt,1950) nrk(1,nrc),ureac(nrc)(1:j2)
            write (nttyo,1950) nrk(1,nrc),ureac(nrc)(1:j2)
1950 format(/' * Error - (EQ6/rd6inw) The forward rate law code',/7x,'has an unrecognized value of ',i2,' for ',a,'.')

            go to 990
        end if

        ! Read backward (precipitation, formation) rate data.
        !   If nrk(2,nrc) = 0, then the backward may be governed by
        !   instantaneous partial equilibrium.
        !   If nrk(2,nrc) = -1, then the forward rate law will be used.
        !   This is legal only for the cases of nrk(1,nrc) = 2
        !   (transition state theory net rate).
        !   If nrk(2,nrc) <= 0, don't enter any additional backward
        !   rate law data.
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
                write (noutpt,1960) ureac(nrc)(1:j2),imchmx
                write (nttyo,1960) ureac(nrc)(1:j2),imchmx
1960 format(/' * Error - (EQ6/rd6inw) Have too many rate',/7x,'constants or corresponding mechanisms in the backward',/7x,'rate law for reactant ',a,'. The code is only',/7x,'dimensioned for ',i2,' rate constants per rate law.',/7x,'Reduce the number of rate constants or increase the',/7x,'dimensioning parameter imchpa.')

                go to 990
            end if

            read (ninpts,1790,err=990) (rkb(i,2,nrc), i = 1,3)
            write (noutpt,1800) (rkb(i,2,nrc), i = 1,3)
        else if (nrk(2,nrc) .eq. 2) then
            ! Transition state rate law.
            read (ninpts,1810,err=990) imech(2,nrc)
            write (noutpt,1820) imech(2,nrc)

            if (imech(2,nrc) .gt. imchmx) then
                write (noutpt,1960) ureac(nrc)(1:j2),imchmx
                write (nttyo,1960) ureac(nrc)(1:j2),imchmx
                go to 990
            end if

            do i = 1,imech(2,nrc)
                read (ninpts,1830,err=990) rkb(i,2,nrc),trkb(i,2,nrc),iact(i,2,nrc)
                write (noutpt,1840) rkb(i,2,nrc),trkb(i,2,nrc),iact(i,2,nrc)

                read (ninpts,1850,err=990) eact(i,2,nrc),hact(i,2,nrc)
                write (noutpt,1860) eact(i,2,nrc),hact(i,2,nrc)

                read (ninpts,1870,err=990) ndact(i,2,nrc),csigma(i,2,nrc)
                write (noutpt,1880) ndact(i,2,nrc),csigma(i,2,nrc)

                if (ndact(i,2,nrc) .gt. ndctmx) then
                    write (noutpt,1970) i,ureac(nrc)(1:j2),ndctmx
                    write (nttyo,1970) i,ureac(nrc)(1:j2),ndctmx
1970 format(/' * Error - (EQ6/rd6inw) Have too many',/7x,'species in the activity product in term ',i2,/7x,'of the backward direction rate law for reactant',/7x,a,'. The code is only dimensioned for ',i3,/7x,'such species. Reduce the number of such species',/7x,'or increase the dimensioning parameter ndctpa.')

                    go to 990
                end if

                do n = 1,ndact(i,2,nrc)
                    read (ninpts,1910,err=990) udac(n,i,2,nrc),cdac(n,i,2,nrc)
                    write (noutpt,1920) udac(n,i,2,nrc),cdac(n,i,2,nrc)
                end do
            end do
        else if (nrk(2,nrc) .eq. 3) then
            ! Linear rate law.
            imech(2,nrc) = 1

            if (imech(2,nrc) .gt. imchmx) then
                write (noutpt,1960) ureac(nrc)(1:j2),imchmx
                write (nttyo,1960) ureac(nrc)(1:j2),imchmx
                go to 990
            end if

            i = 1
            read (ninpts,1830,err=990) rkb(i,2,nrc),trkb(i,2,nrc),iact(i,2,nrc)
            write (noutpt,1840) rkb(i,2,nrc),trkb(i,2,nrc),iact(i,2,nrc)
            read (ninpts,1850,err=990) eact(i,2,nrc),hact(i,2,nrc)
            write (noutpt,1860) eact(i,2,nrc),hact(i,2,nrc)
        else
            write (noutpt,1980) nrk(2,nrc),ureac(nrc)(1:j2)
            write (nttyo,1980) nrk(2,nrc),ureac(nrc)(1:j2)
1980 format(/' * Error - (EQ6/rd6inw) The backward rate law code',/7x,'has an unrecognized value of ',i2,' for ',a,'.')

            go to 990
        end if

        write (noutpt,1500)
    end do

    ! Starting, minimum, and maximum values of key run parameters.
    !   xistti  = starting value of reaction progress for this run
    !   ximaxi  = maximum value of Xi for this run
    !   tistti  = time at Xi = xistti (sec)
    !   timmxi  = maximum time (sec) for this run
    read (ninpts,1150,err=990) xistti,ximaxi,tistti,timmxi
1150 format(2(12x,e12.5),/2(12x,e12.5))

    write (noutpt,1155) xistti,ximaxi,tistti,timmxi
1155 format(5x,'xistti= ',1pe12.5,4x,'ximaxi= ',e12.5,/ 5x,'tistti= ',1pe12.5,4x,'timmxi= ',e12.5)

    ! phmini  = minimum pH for this run
    ! phmaxi  = maximum pH for this run
    ! ehmini  = minimum Eh (v) for this run
    ! ehmaxi  = maximum Eh (v) for this run
    ! o2mini  = minimum log fO2 for this run
    ! o2maxi  = maximum log fO2 for this run
    ! awmini  = minimum aw for this run
    ! awmaxi  = maximum aw for this run
    ! kstpmx  = maximum number of steps for this run
    read (ninpts,1160,err=990) phmini,phmaxi,ehmini,ehmaxi,o2mini,o2maxi,awmini,awmaxi,kstpmx
1160 format(2(12x,e12.5),/2(12x,e12.5),/2(12x,e12.5),/2(12x,e12.5),/12x,i12)

    write (noutpt,1165) phmini,phmaxi,ehmini,ehmaxi,o2mini,o2maxi,awmini,awmaxi,kstpmx
1165 format(5x,'phmini= ',1pe12.5,4x,'phmaxi= ',e12.5,/ 5x,'ehmini= ',1pe12.5,4x,'ehmaxi= ',e12.5,/ 5x,'o2mini= ',1pe12.5,4x,'o2maxi= ',e12.5,/ 5x,'awmini= ',1pe12.5,4x,'awmaxi= ',e12.5,/ 5x,'kstpmx= ',i12)

    ! Print interval parameters.
    !   dlxprn = print increment in terms of Xi
    !   dlxprl = log Xi print interval
    !   dltprn = print increment in terms of time (seconds)
    !   dltprl = log time print interval
    read (ninpts,1170,err=990) dlxprn,dlxprl,dltprn,dltprl
1170 format(2(12x,e12.5),12x,/2(12x,e12.5))

    write (noutpt,1175) dlxprn,dlxprl,dltprn,dltprl
1175 format(5x,'dlxprn= ',1pe12.5,4x,'dlxprl= ',e12.5,/5x,'dltprn= ',e12.5,4x,'dltprl= ',e12.5)

    ! dlhprn = print increment in terms of pH
    ! dleprn = print increment in terms of Eh (v)
    ! dloprn = print increment in terms of log fO2
    ! dlaprn = print increment in terms of aw
    ! ksppmx = limit on number of steps from the last print point
    !            at which a print will be made
    read (ninpts,1180,err=990) dlhprn,dleprn,dloprn,dlaprn,ksppmx
1180 format(2(12x,e12.5),12x,/2(12x,e12.5),/12x,i12)

    write (noutpt,1185) dlhprn,dleprn,dloprn,dlaprn,ksppmx
1185 format(5x,'dlhprn= ',1pe12.5,4x,'dleprn= ',e12.5,/5x,'dloprn= ',e12.5,4x,'dlaprn= ',e12.5,/5x,'ksppmx= ',i12)

    ! Plot interval parameters.
    !   dlxplo = plot increment in terms of Xi
    !   dlxpll = log Xi plot interval
    !   dltplo = plot increment in terms of time (seconds)
    !   dltpll = log time plot interval
    read (ninpts,1190,err=990) dlxplo,dlxpll,dltplo,dltpll
1190 format(2(12x,e12.5),/2(12x,e12.5))

    write (noutpt,1192) dlxplo,dlxpll,dltplo,dltpll
1192 format(5x,'dlxplo= ',1pe12.5,4x,'dlxpll= ',e12.5,/5x,'dltplo= ',e12.5,4x,'dltpll= ',e12.5)

    ! dlhplo = plot increment in terms of pH
    ! dleplo = plot increment in terms of Eh (v)
    ! dloplo = plot increment in terms of log fO2
    ! dlaplo = plot increment in terms of aw
    ! ksplmx = limit on number of steps from the last plot point
    !            at which a plot will be made
    read (ninpts,1195,err=990) dlhplo,dleplo,dloplo,dlaplo,ksplmx
1195 format(2(12x,e12.5),/2(12x,e12.5),/12x,i12)

    write (noutpt,1197) dlhplo,dleplo,dloplo,dlaplo,ksplmx
1197 format(5x,'dlhplo= ',1pe12.5,4x,'dleplo= ',e12.5,/5x,'dloplo= ',e12.5,4x,'dlaplo= ',e12.5,/5x,'ksplmx= ',i12)

    ! Iopt model option switches.
    ! Note: iopt(1) = iopt1, etc.
    !   iopt(1) = Switch to choose the physical model:
    !      0 = Closed system
    !      1 = Titration system
    !      2 = Fluid-centered flow-through open system
    !   iopt(2) = Kinetic mode:
    !      0 = Reaction progress mode (arbitrary kinetics)
    !      1 = Reaction progress/time mode (true kinetics)
    !   iopt(3) = Phase boundary searches:
    !      0 =  The step size is constrained by the locations of
    !             the predicted phase boundaries (iopt(3) is reset
    !             to 0 if other options require this)
    !      1 =  The location of phase boundaries is estimated from
    !             Taylor's series and printed, but the step size is
    !             unconstrained
    !      2  = The locations of phase boundaries are ignored
    !   iopt(4) = Solid solution products:
    !      0 = Solid solutions are ignored
    !      1 = Solid solutions are permitted
    !   iopt(5) = Clear the ES solids read from the input file:
    !      0 = Don't do it
    !      1 = Do it
    !   iopt(6) = Clear the ES solids at the initial value of reaction
    !               progress; this is done after calculating the state
    !               of the system at this point:
    !      0 = Don't do it
    !      1 = Do it
    !   iopt(7) = Clear the ES solids at the end of the run:
    !      0 = Don't do it
    !      1 = Do it, unless the run terminates early due to
    !            calculational difficulties
    !   iopt(8) = Not used
    !   iopt(9) = Clear the PRS solids read from the input file:
    !      0 = Don't do it
    !      1 = Do it
    !   iopt(10) = Clear the PRS solids at the end of the run:
    !      0 = Don't do it
    !      1 = Do it, unless the run terminates early due to
    !            calculational difficulties
    !   iopt(11) = Auto basis switching, in pre-N-R optimization:
    !      0 = Off
    !      1 = On
    !   iopt(12) = Auto basis switching, after N-R iteration:
    !      0 = Off
    !      1 = On
    !   iopt(13) = Switch to choose calculational mode:
    !      0 = Normal path tracing
    !      1 = Economy mode (if permissible)
    !      2 = Super economy mode (if permissible)
    !   iopt(14) = ODE integrator corrector mode:
    !      0 = Allow stiff and simple correctors
    !      1 = Allow only the simple corrector
    !      2 = Allow only the stiff corrector
    !      3 = Allow no correctors
    !   iopt(15) = Force suppression of all redox reactions:
    !      0 = Don't do it
    !      1 = Do it (this is potentially dangerous)
    !   iopt(16) = Backup file options:
    !     -1 = Don't write backup files
    !      0 = Write backup files
    !      1 = Write a sequential backup file
    !   iopt(17) = Pickup file options:
    !     -1 = Don't write a pickup file
    !      0 = Write a pickup file
    !   iopt(18) = Tab file options:
    !     -1 = Don't write the tab file
    !      0 = Write the tab file
    !      1 = The tab file output is appended to the tabx file from a
    !             previous run; if more than one problem is stacked on
    !             the input file, this option applies to only the first
    !             problem, for any subsequent any subsequent problems,
    !             iopt(18) is reset to 0
    !   iopt(19) = Used only by EQ3NR
    !   iopt(20) = Advanced EQ6 pickup file options:
    !      0 = Write a normal EQ6 pickup file
    !      1 = Write an EQ6 input file with "Fluid 2" set up as a
    !            special reactant for mixing with another fluid or
    !            reaction with a more generalized equilibrium system
    write (noutpt,1200)
1200 format(' *',15x,'1    2    3    4    5    6    7    8    9   10')

    read (ninpts,1210,err=990) (iopt(i), i = 1,20)
1210 format(12x,10i5)

    write (noutpt,1220) (iopt(i), i = 1,20)
1220 format(3x,'iopt1-10= ',10i5,/2x,'iopt11-20= ',10i5)

    ! Iopr print option switches.
    ! Note: iopr(1) = iopr1, etc.
    !   iopr(1) = List the names of all species read from the
    !               supporting data file:
    !      0 = Don't list
    !      1 = List
    !   iopr(2) = List all reactions:
    !      0 = Don't list
    !      1 = List all reactions (this can be quite lengthy)
    !      2 = Also print the log K values
    !      3 = Also print the coefficients of the interpolating
    !            polynomials
    !   iopr(3) = List the hard core diameters of the aqueous species:
    !      0 = Don't list
    !      1 = List
    !   iopr(4) = Print a table at each print point of the
    !               concentrations, activities, and activity
    !               coefficients of the aqueous species:
    !     -3  = Omit species with molalities < 1.e-8
    !     -2 =  Omit species with molalities < 1.e-12
    !     -1 =  Omit species with molalities < 1.e-20
    !      0 =  Omit species with molalities < 1.e-100
    !      1 =  Include all species
    !   iopr(5) = Print a table at each print point of the cation/H+
    !               activity ratios, anion-H+ activity products, and
    !               neutral species activities:
    !      0 = Don't print
    !      1 = Print cation/H+ activity ratios only
    !      2 = Print cation/H+ activity ratios and anion-H+ activity
    !            products only
    !      3 = Print cation/H+ activity ratios, anion-H+ activity
    !            products, and neutral species activities
    !   iopr(6) = At each print point, print a table of the percentage
    !                contributions for each aqueous mass balance total:
    !     -1 = Don't print any tables
    !      0 = Print tables including 99% of all contributing species
    !      1 = Print tables including all contributing species
    !   iopr(7) = Print a table at each print point of the saturation
    !               indices and affinities of the various non-aqueous
    !               phases:
    !     -1 = Don't print
    !      0 = Print for those phases not undersaturated by
    !            more than 10 kcal
    !      1 = Print for all phases
    !   iopr(8) = Print a table at each print point of the fugacities
    !               of the gas species:
    !     -1 = Don't print
    !      0 = Print
    !      1 = Print
    !   iopr(9) = Print a table at each print of the mean molal
    !               activity coefficients:
    !     -1 = Don't print
    !      0 = Don't print
    !      1 = Print
    !   iopr(10) = Print a tabulation at the start of running the
    !                current problem of the Pitzer interaction
    !                coefficients:
    !      0 = Don't print
    !      1 = Print a summary of the names of the species present and
    !            the number of Pitzer interaction coefficients
    !      2 = Print a summary of the names of the species present and
    !            the number of Pitzer interaction coefficients
    !   iopr(11) - iopr(16) = Not used
    !   iopr(17) = Pickup file format:
    !      0 = Use the same format ("D" or "W") as the input file
    !      1 = Use "W" format
    !      2 = Use "D" format
    !   iopr(18) - iopr(20) = Not used
    read (ninpts,1230,err=990) (iopr(i), i = 1,20)
1230 format(12x,10i5)

    write (noutpt,1240) (iopr(i), i = 1,20)
1240 format(3x,'iopr1-10= ',10i5,/2x,'iopr11-20= ',10i5)

    ! Iodb debugging print option switches.
    ! Note: iodb(1) = iodb1, etc.
    !   iodb(1) = General diagnostic messages:
    !      0 = Don't print
    !      1 = Print Level 1 diagnostic messages
    !      2 = Print Level 1 and Level 2 diagnostic messages
    !   iodb(2) = Kinetics-related diagnostic messages:
    !      0 = Don't print
    !      1 = Print Level 1 kinetics diagnostic messages
    !      2 = Print Level 1 and Level 2 kinetics diagnostic messages
    !   iodb(3) = Pre-Newton-Raphson optimization:
    !      0 = Don't print
    !      1 = Print summary information
    !      2 = Print detailed information
    !      3 = Print more detailed information
    !      4 = Also print changes to activity coefficients
    !   iodb(4) = Information describing Newton-Raphson iteration:
    !      0 = Don't print
    !      1 = Print summary information
    !      2 = Print detailed information including the residual and
    !            correction vectors
    !      3 = Also print the Jacobian matrix
    !      4 = Also print changes to activity coefficients
    !   iodb(5) = Step size/order selection:
    !      0 = Don't print
    !      1 = Print the chosen scale factor
    !      2 = Print the orders under consideration and their respective
    !            step size scaling factors
    !   iodb(6) = Hypothetical affinity calculations (for solid
    !               solutions, etc.):
    !      0 = Don't print
    !      1 = Print summary information
    !      2 = Print detailed information
    !   iodb(7) = Search iterations (for phase boundaries, etc.):
    !      0 = Don't print
    !      1 = Print summary information
    !   iodb(8) = Information describing ODE corrector iteration:
    !      0 = Don't print
    !      1 = Print summary information
    !      2 = Print detailed information
    !   iodb(9) - iodb(20) = Not used
    read (ninpts,1250,err=990) (iodb(i), i = 1,20)
1250 format(12x,10i5)

    write (noutpt,1260) (iodb(i), i = 1,20)
1260 format(3x,'iodb1-10= ',10i5,/2x,'iodb11-20= ',10i5)

    ! Number of nxopt options.
    !   nxopt = number of mineral subset-selection suppression options
    read (ninpts,1300,err=990) nxopt
1300 format(12x,i2)

    write (noutpt,1310) nxopt
1310 format(6x,'nxopt= ',i2)

    if (nxopt .gt. nxopmx) then
        write (noutpt,1320) nxopmx
        write (nttyo,1320) nxopmx
1320 format(/' * Error - (EQ6/rd6inw) Have too many mineral',/7x,'subset-selection suppression options. The code is',/7x,'only dimensioned for ',i3,' such options. Reduce the',/7x,'number of options or increase the dimensioning',/7x,'parameter nxoppa.')

        go to 990
    end if

    ! Nxopt options.
    !   uxopt  = 'all' or 'alwith'
    !   uxcat  = name of a chemical element
    !   nxopex = number of exceptions
    !   uxopex = name of the minerals that are the exceptions
    do n = 1,nxopt
        read (ninpts,1330,err=990) uxopt(n),uxcat(n)
1330 format(12x,a6,1x,a24)

        j2 = ilnobl(uxcat(n))
        write (noutpt,1340) uxopt(n),uxcat(n)(1:j2)
1340 format(5x,'option= ',a6,1x,a)
    end do

    if (nxopt .gt. 0) then
        read (ninpts,1350,err=990) nxopex
1350 format(12x,i2)

        write (noutpt,1360) nxopex
1360 format(5x,'nxopex= ',i2)

        if (nxopex .gt. nxpemx) then
            write (noutpt,1370) nxpemx
            write (nttyo,1370) nxpemx
1370 format(/' * Error - (EQ6/rd6inw) Have too many',/7x,'exceptions specified to the mineral subset-selection',/7x,'suppression options. The code is only dimensioned',/7x,'for ',i3,'exceptions. Reduce the number of exceptions',/7x,'or increase the dimensioning parameter nxpepa.')

            go to 990
        end if

        do n = 1,nxopex
            read (ninpts,1380,err=990) uxopex(n)
1380 format(12x,a24)

            j2 = ilnobl(uxopex(n))
            write (noutpt,1390) uxopex(n)(1:j2)
1390 format(2x,'exception= ',a)
        end do
    end if

    ! Number of nffg options.
    !   nffg   = number of gases with fixed fugacities
    read (ninpts,1400,err=990) nffg
1400 format(12x,i2)

    write (noutpt,1410) nffg
1410 format(7x,'nffg= ',i2)

    if (nffg .gt. nffgmx) then
        write (noutpt,1420) nffgmx
        write (nttyo,1420) nffgmx
1420 format(/' * Error - (EQ6/rd6inw) Have too many gases whose',/7x,'fugacities are to be fixed. The code is only dimensioned',/7x,'for ',i4,' such gases. Reduce the number of gases or',/7x,'increase the dimensioning parameter nffgpa.')

        go to 990
    end if

    ! Nffg options.
    !   uffg   = name of species
    !   moffg  = moles of gas species
    !   xlkffg = log fugacity
    do n = 1,nffg
        read (ninpts,1430,err=990) uffg(n),moffg(n),xlkffg(n)
1430 format(12x,a24,/12x,e12.5,12x,e12.5)

        write (noutpt,1440) uffg(n),moffg(n),xlkffg(n)
1440 format(4x,'species= ',a24,/6x,'moffg= ',1pe12.5,4x,'xlkffg= ',e12.5)
    end do

    ! Under normal circumstances, users should enter zeros and take
    ! the code default values for dlxdmp, tolbt, toldl, tolxsf, tolsat,
    ! dlxmx0, itermx, and ntrymx.
    ! Maximum finite-difference order.
    read (ninpts,2010,err=990) nordmx
2010 format(12x,i3)

    ! Newton-Raphson convergence tolerances.
    !   tolbt = convergence tolerance on residual magnitude
    !             (in Newton-Raphson iteration)
    !   toldl = convergence tolerance on correction magnitude
    !             (in Newton-Raphson iteration)
    read (ninpts,2020,err=990) tolbt,toldl
2020 format(2(12x,e12.5))

    write (noutpt,2030) tolbt,toldl
2030 format(6x,'tolbt= ',1pe12.5,5x,'toldl= ',e12.5)

    ! Maximum number of Newton-Raphson iterations.
    !   itermx = limit on the number of Newton-Raphson iterations.
    !              Recommended values lie in the range 50 to 200.
    read (ninpts,2100,err=990) itermx
2100 format(12x,i3)

    write (noutpt,2110) itermx
2110 format(5x,'itermx= ',i3)

    ! Search/find tolerance.
    !   tolxsf = search/find tolerance
    read (ninpts,2022,err=990) tolxsf
2022 format(12x,e12.5)

    write (noutpt,2032) tolxsf
2032 format(5x,'tolxsf= ', 1pe12.5)

    ! Saturation tolerance.
    !   tolsat = supersaturation tolerance below which no attempt
    !              is made to precipitate a phase. A good value is
    !              0.0001 kcal. Don't make it too small.
    read (ninpts,2024,err=990) tolsat
2024 format(12x,e12.5)

    write (noutpt,2034) tolsat
2034 format(5x,'tolxsf= ', 1pe12.5)

    ! Maximum number of phase assemblage tries.
    !   ntrymx = limit on the number of attempted phase assemblages
    !              at a point of reaction progress. Recommended
    !              values lie in the range of 50 to 100.
    read (ninpts,2102,err=990) ntrymx
2102 format(12x,i3)

    write (noutpt,2112) ntrymx
2112 format(5x,'ntrymx= ',i3)

    ! Zero-order step size (in Xi).
    !   dlxmx0 = the normal step size for order zero; this is not the
    !              minimum step size (dlxmin). Recommended values lie
    !              in the range 1.e-10 to 1.e-6.
    read (ninpts,2080,err=990) dlxmx0
2080 format(12x,e12.5)

    write (noutpt,2090) dlxmx0
2090 format(5x,'dlxmx0= ',1pe12.5)

    ! PRS transfer interval in Xi.
    !   dlxdmp = fixed reaction progress interval for the transfer of
    !              phases from the equilibrium system (ES) to the
    !              physically removed system (PRS); this only works
    !              in conjunction with the fluid-centered flow-through
    !              open system model (iopt(1) = 2)
    read (ninpts,2000,err=990) dlxdmp
2000 format(12x,e12.5)

    write (noutpt,2005) dlxdmp
2005 format(5x,'dlxdmp= ',1pe12.5)

    write (noutpt,1500)

    ! Process the bottom half of the current input file.
    ! Secondary title.
    ! utitl2 = secondary title
    ! ntitl2 = number of lines in the secondary title.
    ntitl2 = 0
    read (ninpts,1000,err=990) uline

    j2 = ilnobl(uline)
    j2 = min(j2,79)
    write (noutpt,1020) uline

    if (uline(1:8) .ne. uendit(1:8)) then
        ntitl2 = 1
        utitl2(1) = uline

        do n = 2,ntitmx
            read (ninpts,1000,err=990) uline
            j2 = ilnobl(uline)
            j2 = min(j2,79)
            write (noutpt,1020) uline(1:j2)

            if (uline(1:8) .eq. uendit(1:8)) then
                go to 220
            end if

            utitl2(n) = uline
            ntitl2 = n
        end do

220 continue
    end if

    ! Special basis switches.
    !   nsbswt = the number of special basis switches
    !   usbsw(1,n) = the species to be switched from the strict
    !                  basis to the auxiliary basis set in the
    !                  in the n-th switch
    !   usbsw(2,n) = the species to be switched from the auxiliary
    !                  basis set to the strict basis set in the
    !                  n-th switch
    write (noutpt,2610)
2610 format(' *   Special basis switches')

    read (ninpts,2620,err=990) nsbswt
2620 format(12x,i3)

    write (noutpt,2630) nsbswt
2630 format(5x,'nsbswt= ',i3)

    do n = 1,nsbswt
        read (ninpts,2640,err=990) usbsw(1,n)
2640 format(9x,a48)

        j2 = ilnobl(usbsw(1,n))
        write (noutpt,2650) usbsw(1,n)(1:j2)
2650 format(1x,'species= ',a)

        read (ninpts,2660,err=990) usbsw(2,n)
2660 format(15x,a48)

        j2 = ilnobl(usbsw(2,n))
        write (noutpt,2670) usbsw(2,n)(1:j2)
2670 format(3x,'switch with= ',a)
    end do

    ! Original temperature.
    !   tempci = temperature, C, at the end of the previous
    !              EQ3NR or EQ6 run
    read (ninpts,2680,err=990) tempci
2680 format(12x,e12.5)

    write (noutpt,2690) tempci
2690 format(5x,'tempci= ',1pe12.5)

    ! Original pressure.
    !   pressi = pressure, bars, at the end of the previous
    !              EQ3NR or EQ6 run
    read (ninpts,2700,err=990) pressi
2700 format(12x,e12.5)

    write (noutpt,2710) pressi
2710 format(5x,'pressi= ',1pe12.5)

    ! Ion exchanger creation. This section reads the directives for
    ! creating ion exchange phases and species and their associated
    ! intrinsic (including thermodynamic) properties. The data is
    ! essentially that which might be read from a supporting data file.
    !   qgexsh = logical flag. If .true., force the writing of at
    !              least one exchanger block on a "D" format input
    !              file (when running XCON6) or pickup file.
    !   net    = the number of ion exchangers
    !   ugexp  = array of ion exchange phase names
    !   mwtges = array of molecular weights of substrates of ion
    !              exchange phases
    !   ugexmo = array of strings identifying models for composing
    !              exact ion exchange species and corresponding
    !              reactions; examples include:
    !                'Vanselow' = Vanselow model
    !                'Gapon'    = Gapon model
    !                'Site-mixing' = the general site-mixing model
    !   tgexp  = array of reference temperatures (C)  for the
    !              thermodynamic data for the reactions specified
    !              for the exchanger phases
    !   jgext  = array of numbers of exchange sites in ion exchange
    !              phases
    !   ugexj  = array of names specified for the sites on the various
    !             exchangers
    !   cgexj  = array of numbers of moles of exchange sites per mole
    !              of substrate in ion exchange phases
    !   zgexj  = array of electrical charges of exchange sites
    !              in ion exchange phases. The total charge number
    !              for a site is the product of this charge and
    !              and the number of moles of the site per mole
    !              of exchanger. The sign of the charge on a site
    !              is generally the opposite of that of the ions
    !              which exchange on it.
    !   ngexrt = array of numbers of ion exchange reactions specified
    !              for a given site and exchanger
    !   ugexr  = array of strings containing compact representations
    !              of the exchange reactions; e.g., 'Na+ = Ca++' for a
    !              reaction in which Na+ on the exchanger is replaced
    !              by Ca++. One may make specifications such as
    !              'Na+ = ' in which case the ion goes into solution
    !              leaving a bare substrate. All reactions are
    !              normalized to the exchange (or loss) of one
    !              equivalent. The exact form of the reaction is
    !              otherwise dependent on the mixing law specifed in
    !              the element of the ugexmo array for the current
    !              exchanger.
    !   xlkgex = array of equilibrium constants or related data for
    !              the reaction specified in the above string
    !   uxkgex = array of strings denoting the kind of data in the
    !              corresponding entry of the xlkgex array:
    !                ' '       = log K per equivalent
    !                'LogK/eq' = log K per equivalent
    !                'kcal/eq' = DeltaG0r, kcal, per equivalent
    !                'kJ/eq'   = DeltaG0r, kJ, per equivalent
    !   xhfgex = array of enthalpy of reaction data for the reaction
    !              specified in the above string
    !   uhfgex = array of strings denoting the kind of data in the
    !              corresponding entry of the xhfgex array:
    !                ' '       = DeltaG0r, kcal, per equivalent
    !                'kcal/eq' = DeltaG0r, kcal, per equivalent
    !                'kJ/eq'   = DeltaG0r, kJ, per equivalent
    !   xvfgex = array of volume of reaction data for the reaction
    !              specified in the above string
    !   uvfgex = array of strings denoting the kind of data in the
    !              corresponding entry of the xhfgex array:
    !                ' '      = cm3 per equivalent
    !                'cm3/eq' = cm3 per equivalent
    write (noutpt,2790)
2790 format(' * Ion exchanger creation')

    read (ninpts,2795) ux8
2795 format(12x,a8)

    call lejust(ux8)
    call locase(ux8)
    qgexsh = ux8(1:2).eq.'.t' .or. ux8(1:1).eq. 't'
    write (noutpt,2797) qgexsh
2797 format(5x,'qgexsh= ',l8)

    read (ninpts,2800,err=990) net
2800 format(12x,i3)

    write (noutpt,2810) net
2810 format(8x,'net= ',i3)

    if (net .gt. netmax) then
        write (noutpt,2820) netmax,net
        write (nttyo,2820) netmax,net
2820 format(/' * Error - (EQ6/rd6inw) Have exceeded the maximum',' number of ',i3,/7x,'generic ion exchange phases while',' reading the data to create',/7x,'such phases. Increase',' the dimensioning parameter netpar',/7x,'to at least ',i3,'.')

        go to 990
    end if

    do ne = 1,net
        read (ninpts,2830,err=990) ugexp(ne)
2830 format(12x,a24)

        j2 = ilnobl(ugexp(ne))
        write (noutpt,2840) ugexp(ne)(1:j2)
2840 format(6x,'ugexp= ',a)

        read (ninpts,2850,err=990) mwtges(ne)
2850 format(12x,e12.5)

        write (noutpt,2860) mwtges(ne)
2860 format(5x,'mwtges= ',1pe12.5)

        read (ninpts,2830,err=990) ugexmo(ne)
        j3 = ilnobl(ugexmo(ne))
        write (noutpt,2870) ugexmo(ne)(1:j3)
2870 format(5x,'ugexmo= ',a)

        read (ninpts,2850,err=990) tgexp(ne)
        write (noutpt,2880) tgexp(ne)
2880 format(6x,'tgexp= ',1pe12.5)

        read (ninpts,2890,err=990) jgext(ne)
2890 format(12x,i3)

        write (noutpt,2900) jgext(ne)
2900 format(6x,'jgext= ',i3)

        if (jgext(ne) .gt. jetmax) then
            write (noutpt,2910) jetmax,ugexp(ne)(1:j2),jgext(ne)
            write (nttyo,2910) jetmax,ugexp(ne)(1:j2),jgext(ne)
2910 format(/' * Error - (EQ6/rd6inw) Have exceeded the maximum',' number of ',i3,/7x,'exchange sites on a generic ion',' exchange phase while reading',/7x,'the data to create',a,'. Increase the',/7x,'dimensioning parameter jetpar to',' at least ',i3,'.')

            go to 990
        end if

        do je = 1,jgext(ne)
            read (ninpts,2930,err=990) ugexj(je,ne)
2930 format(12x,a8)

            j3 = ilnobl(ugexj(je,ne))
            write (noutpt,2940) ugexj(je,ne)(1:j3)
2940 format(6x,'ugexj= ',a)

            read (ninpts,2950,err=990) cgexj(je,ne),zgexj(je,ne)
2950 format(12x,e12.5,12x,e12.5)

            write (noutpt,2960) cgexj(je,ne),zgexj(je,ne)
2960 format(6x,'cgexj= ',1pe12.5,5x,'zgexj= ',e12.5)

            read (ninpts,2890,err=990) ngexrt(je,ne)
            write (noutpt,3010) ngexrt(je,ne)
3010 format(5x,'ngexrt= ',i3)

            if (ngexrt(je,ne) .gt. ietmax) then
                write (noutpt,3020) netmax,ugexj(je,ne)(1:j3),ugexp(ne)(1:j2),ngexrt(je,ne)
                write (nttyo,3020) netmax,ugexj(je,ne)(1:j3),ugexp(ne)(1:j2),ngexrt(je,ne)
3020 format(/' * Error - (EQ6/rd6inw) Have exceeded the',' maximum number of ',i3,/7x,'reactions for a site',' belonging to a generic ion exchange',/7x,'phase while',' reading the data for site ',a,' of exchange phase'      /7x,a,'. Increase the dimensioning parameter',/7x,'ietpar to at least ',i3,'.')

                go to 990
            end if

            do n = 1,ngexrt(je,ne)
                read (ninpts,3030,err=990) ugexr(n,je,ne)
3030 format(12x,a56)

                j2 = ilnobl(ugexr(n,je,ne))
                write (noutpt,3040) ugexr(n,je,ne)(1:j2)
3040 format(6x,'ugexr= ',a)

                read (ninpts,3050,err=990) xlkgex(n,je,ne),uxkgex(n,je,ne)
3050 format(12x,e12.5,12x,a8)

                j4 = ilnobl(uxkgex(n,je,ne))
                write (noutpt,3060) xlkgex(n,je,ne),uxkgex(n,je,ne)(1:j4)
3060 format(5x,'xlkgex= ',1pe12.5,5x,'units= ',a)

                read (ninpts,3050,err=990) xhfgex(n,je,ne),uhfgex(n,je,ne)
                j4 = ilnobl(uhfgex(n,je,ne))
                write (noutpt,3070) xhfgex(n,je,ne),uhfgex(n,je,ne)(1:j4)
3070 format(5x,'xhfgex= ',1pe12.5,5x,'units= ',a)

                read (ninpts,3050,err=990) xvfgex(n,je,ne),uvfgex(n,je,ne)
                j4 = ilnobl(uvfgex(n,je,ne))
                write (noutpt,3080) xvfgex(n,je,ne),uvfgex(n,je,ne)(1:j4)
3080 format(5x,'xvfgex= ',1pe12.5,5x,'units= ',a)
            end do
        end do
    end do

    ! Number of nxmod options.
    !   nxmod = the number of suppressed/altered species/reactions
    read (ninpts,3120,err=990) nxmod
3120 format(12x,i2)

    write (noutpt,3130) nxmod
3130 format(6x,'nxmod= ',i2)

    if (nxmod .gt. nxmdmx) then
        write (noutpt,3140) nxmdmx
        write (nttyo,3140) nxmdmx
3140 format(/' * Error - (EQ6/rd6inw) Have too many nxmod',/7x,'alter/suppress options. The code is only dimensioned',/7x,'for ',i2,' such options. Reduce the number of such',/7x,'options or increase the dimensioning parameter nxmdpa.')

        go to 990
    end if

    ! Nxmod options.
    !   uxmod = the name (all 48 letters) of the species. If the
    !             phase part is not given, the option is applied to
    !             every species for which the species part of its name
    !             generates a match.
    !   kxmod  = alter/suppress code
    !     -1 = the species is suppressed
    !      0 = the log K is replaced by xlkmod
    !      1 = the log K is augmented by xlkmod
    !      2 = same as kxmod=1, but xlkmod is input in units of kcal
    !            per mole of the associated species
    !   xlkmod = log K value alteration function as defined above
    do n = 1,nxmod
        read (ninpts,3150,err=990) uxmod(n),kxmod(n),xlkmod(n)
3150 format(12x,a48,/12x,i2,22x,e12.5)

        j2 = ilnobl(uxmod(n))
        write (noutpt,3160) uxmod(n)(1:j2),kxmod(n),xlkmod(n)
3160 format(4x,'species= ',a,/5x,'option= ',i2,14x,'xlkmod= ',1pe12.5)
    end do

    ! Iopg options.
    ! Note: iopg(1) = iopg1, etc.
    !   iopg(1) = Model for aqueous species activity coefficients:
    !     -1 = The  Davies equation
    !      0 = The B-dot equation
    !      1 = Pitzer's  equations
    !      2 = HC + DH equations
    !   iopg(2) = Rescaling of aqueous ionic activity coefficients for
    !               consistency with a desired pH scale:
    !     -1 = "Internal" pH scale; no scaling (e.g., the "Davies"
    !            scale if iopg(1) = -1, the "B-dot" scale if
    !            iopg(1) = 0, or the Pitzer scale if iopg(1) = 1)
    !      0 = The NBS pH scale (log gamma(Cl-) is defined by the
    !            Bates-Guggenheim equation)
    !      1 = The Mesmer pH scale (log gamma(H+) = 0)
    !   iopg(3) = iopg(10) = Not used
    read (ninpts,3170,err=990) (iopg(i), i = 1,20)
3170 format(12x,10i5)

    write (noutpt,1200)
    write (noutpt,3180) (iopg(i), i = 1,20)
3180 format(3x,'iopg1-10= ',10i5,/2x,'iopg11-20= ',10i5)

    ! Index limits.
    !   kct    = the number of elements in the matrix
    !   kbt    = the number of basis species in the matrix
    !   kmt    = position of last pure mineral in the matrix
    !   kxt    = position of last solid solution in the matrix
    !   kdim   = size of the matrix
    !   kprs   = flag to read input for initializing the physically
    !              removed system (PRS)
    read (ninpts,3220,err=990) kct,kbt,kmt,kxt,kdim,kprs
3220 format(3(12x,i2,10x),/3(12x,i2,10x))

    write (noutpt,3230) kct,kbt,kmt,kxt,kdim,kprs
3230 format(8x,'kct= ',i2,17x,'kbt= ',i2,17x,'kmt= ',i2,/ 8x,'kxt= ',i2,16x,'kdim= ',i2,16x,'kprs= ',i2)

    nbti = kbt

    if (kct .gt. nctmax) then
        write (noutpt,3240) nctmax
        write (nttyo,3240) nctmax
3240 format(/' * Error - (EQ6/rd6inw) Have too many chemical',/7x,'elements present. The code is only dimensioned',/7x,'for ',i3,' elements. Reduce the number of elements',/7x,'or increase the dimensioning parameter nctpar.')

        go to 990
    end if

    if (kbt .gt. nbtmax) then
        write (noutpt,3250) nbtmax
        write (nttyo,3250) nbtmax
3250 format(/' * Error - (EQ6/rd6inw) Have too many basis',/7x,'species present. The code is only dimensioned',/7x,'for ',i3,' basis species. Reduce the number of elements',/7x,'or increase the dimensioning parameter nctpar.')

        go to 990
    end if

    if (kdim .gt. kmax) then
        write (noutpt,3260) kmax
        write (nttyo,3260) kmax
3260 format(/' * Error - (EQ6/rd6inw) Have too many matrix',/7x,'variables. The code is only dimensioned for ',i3,/7x,'matrix variables. Reduce the number of such variables',/7x,'or increase the dimensioning parameter kpar.')

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
    write (noutpt,3330)
3330 format(' * Data file basis species and jflag values')

    do nbi = 1,nbti
        read (ninpts,3340,err=990) ubmtbi(nbi),jflgi(nbi)
3340 format(3x,a48,3x,i2)

        write (noutpt,3350) ubmtbi(nbi),jflgi(nbi)
3350 format(4x,a48,3x,i2)
    end do

    ! Mass balance totals.
    !   mtbi   = total number of moles of basis species in the
    !              equilibrium system (ES)
    !   mtbaqi = total number of moles of basis species in the
    !              aqueous solution
    !   electr = electrical imbalance
    write (noutpt,3400)
3400 format(' *',7x,'Mass balance totals     Aqueous mass balance totals')

    do nbi = 1,nbti
        read (ninpts,3410,err=990) mtbi(nbi),mtbaqi(nbi)
3410 format(6x,e22.15,6x,e22.15)

        write (noutpt,3420) mtbi(nbi),mtbaqi(nbi)
3420 format(7x,1pe22.15,6x,e22.15)
    end do

    read (ninpts,3430,err=990) electr
3430 format(34x,e22.15)

    write (noutpt,3440) electr
3440 format(10x,'Electrical imbalance= ',3x,1pe22.15)

    ! Ordinary basis switches.
    !   nobswt = the number of ordinary basis switches
    !   uobsw(1,n) = the species to be switched from the basis
    !                  set in the n-th switch
    !   uobsw(2,n) = the species to be switched into the basis
    !                  set in the n-th switch
    write (noutpt,3450)
3450 format(' *   Ordinary basis switches')

    read (ninpts,2620,err=990) nobswt
    write (noutpt,3460) nobswt
3460 format(5x,'nobswt= ',i3)

    do n = 1,nobswt
        read (ninpts,2640,err=990) uobsw(1,n)
        j2 = ilnobl(uobsw(1,n))
        write (noutpt,2650) uobsw(1,n)(1:j2)

        read (ninpts,2660,err=990) uobsw(2,n)
        j2 = ilnobl(uobsw(2,n))
        write (noutpt,2670) uobsw(2,n)(1:j2)
    end do

    ! Matrix column variables.
    !   uzveci = names of species or properties associated with the
    !              matrix column variables
    write (noutpt,3470)
3470 format(' * Matrix species or entities')

    do krow = 1,kdim
        read (ninpts,3480,err=990) uzveci(krow)
3480 format(3x,a48)

        j2 = ilnobl(uzveci(krow))
        write (noutpt,3490) uzveci(krow)(1:j2)
3490 format(4x,a)
    end do

    ! Matrix variable values.
    !   zvclg1 = values of matrix variables
    write (noutpt,3510)
3510 format(' *   Values of matrix variables')

    do kcol = 1,kdim
        read (ninpts,3520,err=990) zvclgi(kcol)
3520 format(3x,e22.15)

        write (noutpt,3530) zvclgi(kcol)
3530 format(4x,1pe22.15)
    end do

    ! Initialize any non-zero values for the physically removed
    ! system. Here mprphi and mprspi are arrays for the number of
    ! moles of phases and species, respectively, in the PRS.
    nprpti = 0
    nprsti = 0

    if (kprs .gt. 0) then
        write (noutpt,3600)
3600 format(' *  Phases and species in the PRS')

        j3 = ilnobl(uendit)

        do npi = 1,nprpmx + 1
            read (ninpts,3610,err=990) uxf,mx
3610 format(1x,a24,6x,e22.15)

            if (uxf(1:8) .eq. uendit(1:8)) then
                write (noutpt,3620) uendit(1:j3)
3620 format(2x,a)

                go to 310
            else
                write (noutpt,3630) uxf,mx
3630 format(2x,a24,6x,1pe22.15)
            end if

            nprpti = nprpti + 1

            if (nprpti .gt. nprpmx) then
                write (noutpt,3640) nprpmx
                write (nttyo,3640) nprpmx
3640 format(/' * Error - (EQ6/rd6inw) Have too many phases',/7x,'in the physically removed system (PRS) on the',/7x,'input file. The code is only dimensioned for ',i3,/7x,'such phases. Reduce the number of such phases',/7x,'or increase the dimensioning parameter nprppa.')

                go to 990
            end if

            uprphi(nprpti)(1:24) = uxf
            mprphi(nprpti) = mx

            do iki = 1,iktmax + 1
                read (ninpts,3650,err=990) uxg,mx
3650 format(3x,a24,6x,e22.15)

                if (uxg(1:8) .eq. uendit(1:8)) then
                    write (noutpt,3660) uendit(1:j3)
3660 format(4x,a)

                    go to 300
                else
                    write (noutpt,3670) uxg,mx
3670 format(4x,a24,6x,1pe22.15)
                end if

                if (iki .gt. iktmax) then
                    j2 = ilnobl(uxf)
                    write (noutpt,3680) uxf(1:j2),iktmax
                    write (nttyo,3680) uxf(1:j2),iktmax
3680 format(/' * Error - (EQ6/rd6inw) Have too many',' end-members',/7x,'in the PRS solid solution ',a,'.',/7x,'The code is only dimensioned for ',i4,' end-members per',/7x,'solid solution. Reduce',' the number of end-members or',/7x,'increase the dimensioning parameter iktpar.')

                    go to 990
                end if

                nprsti = nprsti + 1

                if (nprsti .gt. nprsmx) then
                    write (noutpt,3690) nprsmx
                    write (nttyo,3690) nprsmx
3690 format(/' * Error - (EQ6/rd6inw) Have too many species',/7x,'in the physically removed system (PRS) on the',/7x,'input file. The code is only dimensioned for ',i3,/7x,'such species. Reduce the number of such species',/7x,'or increase the dimensioning parameter nprspa.')

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

    write (noutpt,4000) nprob
    write (nttyo,4000) nprob
4000 format(/'   Done reading problem ',i3,'.',/)

    go to 999

990 continue
    qrderr = .true.

999 continue
end subroutine rd6inw