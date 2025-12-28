subroutine stpk36(awmaxi,awmini,cbsri,cbsr1,cdac,cesri,cesr1,csigma,dlaplo,dlaprn,dleplo,dleprn,dlhplo,dlhprn,dloplo,dloprn,dltpll,dltplo,dltprl,dltprn,dlxdmp,dlxmx0,dlxpll,dlxplo,dlxprl,dlxprn,eact,ehmaxi,ehmini,fkrc,hact,iact,ibsrti,ibsrt1,iesrti,iesrt1,iktmax,imchmx,imech,iopt,itermx,ixrti,jcode,jpress,jreac,jtemp,ksplmx,ksppmx,kstpmx,modr,moffg,morr,mprphi,mprspi,nbt1mx,nctmax,ndact,ndctmx,nffg,nffgmx,noptmx,nordmx,noutpt,nprob,nprpmx,nprpti,nprsmx,nprsti,nptkmx,nrct,nrctmx,nrk,nsk,nsrt,nsrtmx,ntitl1,ntitmx,ntrymx,nttkmx,nttyo,nxopex,nxopmx,nxopt,nxpemx,nxrt,nxrtmx,o2maxi,o2mini,phmaxi,phmini,ptk,press,pressb,rkb,rxbari,sfcar,ssfcar,tempc,tempcb,tempc1,timmxi,tistti,tolsat,tolxsf,trkb,ttk,ubsri,ubsr1,ucxri,udac,uesri,uesr1,uffg,uprphi,uprspi,ureac,ureac1,utitl1,uxcat,uxopex,uxopt,vreac,ximaxi,xistti,xlkffg)
    !! This subroutine sets up certain variables and arrays for writing
    !! the top half of an EQ6 input file. There are currently three
    !! options for the form and content of this part of this file.
    !! (see below). EQ3NR/setpk3.f performs a somewhat similiar function
    !! for data to be written on the normal EQ3NR pickup file, or the
    !! bottom half of an EQ6 input file.
    !! This subroutine should only be called if iopt(19) is greater than
    !! zero.
    !! This subroutine is called by:
    !!   EQ3NR/eq3nr.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: iktmax
    integer :: imchmx
    integer :: nbt1mx
    integer :: nctmax
    integer :: ndctmx
    integer :: nffgmx
    integer :: noptmx
    integer :: nordmx
    integer :: nprpmx
    integer :: nprsmx
    integer :: nptkmx
    integer :: nsrtmx
    integer :: ntitmx
    integer :: nttkmx
    integer :: nrctmx
    integer :: nxopmx
    integer :: nxpemx
    integer :: nxrtmx

    integer :: noutpt
    integer :: nttyo

    integer :: iact(imchmx,2,nrctmx)
    integer :: ibsrti(nsrtmx)
    integer :: iesrti(nsrtmx)
    integer :: imech(2,nrctmx)
    integer :: iopt(noptmx)
    integer :: ixrti(nxrtmx)
    integer :: jcode(nrctmx)
    integer :: jreac(nrctmx)
    integer :: ndact(imchmx,2,nrctmx)
    integer :: nrk(2,nrctmx)
    integer :: nsk(nrctmx)

    integer :: ibsrt1
    integer :: iesrt1
    integer :: itermx
    integer :: jpress
    integer :: jtemp
    integer :: ksplmx
    integer :: ksppmx
    integer :: kstpmx
    integer :: nffg
    integer :: nprob
    integer :: nprpti
    integer :: nprsti
    integer :: nrct
    integer :: nsrt
    integer :: ntitl1
    integer :: ntrymx
    integer :: nxopex
    integer :: nxopt
    integer :: nxrt

    character(len=80) :: utitl1(ntitmx)
    character(len=48) :: uprspi(nprsmx)
    character(len=24) :: ubsri(nbt1mx,nsrtmx)
    character(len=24) :: ubsr1(nbt1mx)
    character(len=24) :: ucxri(iktmax,nxrtmx)
    character(len=24) :: udac(ndctmx,imchmx,2,nrctmx)
    character(len=24) :: uffg(nffgmx)
    character(len=24) :: uprphi(nprpmx)
    character(len=24) :: ureac(nrctmx)
    character(len=24) :: uxcat(nxopmx)
    character(len=24) :: uxopex(nxpemx)
    character(len=8) :: uesri(nctmax,nsrtmx)
    character(len=8) :: uesr1(nctmax)
    character(len=8) :: uxopt(nxopmx)

    character(len=24) :: ureac1

    real(kind=8) :: cbsri(nbt1mx,nsrtmx)
    real(kind=8) :: cbsr1(nbt1mx)
    real(kind=8) :: cdac(ndctmx,imchmx,2,nrctmx)
    real(kind=8) :: cesri(nctmax,nsrtmx)
    real(kind=8) :: cesr1(nctmax)
    real(kind=8) :: csigma(imchmx,2,nrctmx)
    real(kind=8) :: eact(imchmx,2,nrctmx)
    real(kind=8) :: fkrc(nrctmx)
    real(kind=8) :: hact(imchmx,2,nrctmx)
    real(kind=8) :: modr(nrctmx)
    real(kind=8) :: moffg(nffgmx)
    real(kind=8) :: morr(nrctmx)
    real(kind=8) :: mprphi(nprpmx)
    real(kind=8) :: mprspi(nprsmx)
    real(kind=8) :: ptk(nptkmx)
    real(kind=8) :: rkb(imchmx,2,nrctmx)
    real(kind=8) :: rxbari(iktmax,nxrtmx)
    real(kind=8) :: sfcar(nrctmx)
    real(kind=8) :: ssfcar(nrctmx)
    real(kind=8) :: trkb(imchmx,2,nrctmx)
    real(kind=8) :: ttk(nttkmx)
    real(kind=8) :: vreac(nrctmx)
    real(kind=8) :: xlkffg(nffgmx)

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
    real(kind=8) :: o2maxi
    real(kind=8) :: o2mini
    real(kind=8) :: phmaxi
    real(kind=8) :: phmini
    real(kind=8) :: press
    real(kind=8) :: pressb
    real(kind=8) :: tempc
    real(kind=8) :: tempcb
    real(kind=8) :: tempc1
    real(kind=8) :: timmxi
    real(kind=8) :: tistti
    real(kind=8) :: tolsat
    real(kind=8) :: tolxsf
    real(kind=8) :: ximaxi
    real(kind=8) :: xistti

    ! Local variable declarations.
    integer :: j2
    integer :: n
    integer :: nmax

    integer :: ilnobl

    character(len=8) :: ux8

    if (iopt(19) .le. 0) then
        go to 999
    end if

    ! Set up the initial part of the new main title.
    call initcb(utitl1,ntitmx)
    ntitl1 = 25
    utitl1(1) = 'EQ6 input file name= sample.6i'
    utitl1(2) = 'Description= "Sample"'
    utitl1(3) = 'Version level= 8.0'
    utitl1(4) = 'Revised mm/dd/yy    Revisor= Username'
    utitl1(6) = '  This is a sample EQ6 input file written as an EQ3NR pickup file.'
    utitl1(7) = 'The EQ3NR input file used to generate this output is identified in'
    utitl1(8) = 'the second title given below.'
    utitl1(10) = '  You are expected to modify this EQ6 input file to meet your own needs.'

    ! Set various EQ6 input file variables.
    ksplmx = 0
    ksppmx = 0
    kstpmx = 500
    ntrymx = 0
    itermx = 0

    xistti = 0.
    ximaxi = 1.0

    tistti = 0.
    timmxi = 1.e+38

    awmaxi = 1.e+38
    awmini = -1.e+38
    ehmaxi = 1.e+38
    ehmini = -1.e+38
    o2maxi = 1.e+38
    o2mini = -1.e+38
    phmaxi = 1.e+38
    phmini = -1.e+38

    dlxdmp = 0.
    dlxmx0 = 0.

    dlxprn = 0.
    dlxprl = 0.5
    dlxplo = 0.
    dlxpll = 0.

    dltprn = 1.e+38
    dltprl = 1.e+38
    dltplo = 1.e+38
    dltpll = 1.e+38

    dlaprn = 1.e+38
    dleprn = 1.e+38
    dlhprn = 1.e+38
    dloprn = 1.e+38
    dlaplo = 1.e+38
    dleplo = 1.e+38
    dlhplo = 1.e+38
    dloplo = 1.e+38

    tolsat = 0.
    tolxsf = 0.
    nordmx = 6

    jtemp = 0
    call initaz(ttk,nttkmx)
    tempcb = tempc

    if (iopt(19) .eq. 3) then
        if (abs(tempc - tempc1) .gt. 1.e-6) then
            ! Set up for non-isothermal fluid-mixing.
            jtemp = 3
            ttk(2) = tempc1
            ttk(1) = 1.0
        end if
    end if

    jpress = 0
    call initaz(ptk,nptkmx)
    pressb = press

    nrct = 0
    call initiz(jcode,nrctmx)
    call initiz(jreac,nrctmx)
    call initiz(nsk,nrctmx)

    nsrt = 0
    call initiz(ibsrti,nsrtmx)
    call initiz(iesrti,nsrtmx)

    nmax = nbt1mx*nsrtmx
    call initaz(cbsri,nmax)
    call initcb(ubsri,nmax)

    nmax = nctmax*nsrtmx
    call initaz(cesri,nmax)
    call initcb(uesri,nmax)

    call initiz(ixrti,nxrtmx)

    nxrt = 0
    nmax = iktmax*nxrtmx
    call initaz(rxbari,nmax)
    call initcb(ucxri,nmax)

    call initcb(ureac,nrctmx)
    call initaz(fkrc,nrctmx)
    call initaz(sfcar,nrctmx)
    call initaz(ssfcar,nrctmx)

    nmax = 2*nrctmx
    call initiz(imech,nmax)
    call initiz(nrk,nmax)

    nmax = imchmx*2*nrctmx
    call initaz(csigma,nmax)
    call initaz(rkb,nmax)
    call initaz(trkb,nmax)
    call initaz(eact,nmax)
    call initaz(hact,nmax)
    call initiz(iact,nmax)
    call initiz(ndact,nmax)

    nmax = ndctmx*imchmx*2*nrctmx
    call initaz(cdac,nmax)
    call initcb(udac,nmax)

    nffg = 0
    call initcb(uffg,nffgmx)
    call initaz(moffg,nffgmx)
    call initaz(xlkffg,nffgmx)

    nxopt = 0
    call initcb(uxopt,nxopmx)
    call initcb(uxcat,nxopmx)
    nxopex = 0
    call initcb(uxopex,nxpemx)

    nprpti = 0
    call initaz(mprphi,nprpmx)
    call initcb(uprphi,nprpmx)

    nprsti = 0
    call initaz(mprspi,nprsmx)
    call initcb(uprspi,nprsmx)

    if (iopt(19) .eq. 1) then
        ! Write an EQ6 input file with a single dissolving reactant,
        ! Quartz. Its dissolution is controlled by a specified
        ! relative rate of 1.0. The basic purpose of this option is
        ! to provide a template for a file on which the reactants
        ! react according to specified relative rates. It is intended
        ! that the user then modify this template to describe the
        ! problem actually desired.
        utitl1(12) = '  This sample file has Quartz dissolving at a relative rate of 1.0.'
        utitl1(13) = 'It may or may not actually run. For example, Quartz may not appear'
        utitl1(14) = 'on the supporting data file.'

        nrct = 1
        ureac(1) = 'Quartz'
        jcode(1) = 0
        jreac(1) = 0
        morr(1) = 1.0
        modr(1) = 0.
        nrk(1,1) = 1
        rkb(1,1,1) = 1.0
    else if (iopt(19) .eq. 2) then
        ! Write an EQ6 input file with a single dissolving reactant,
        ! Albite. Its dissolution is controlled by a specified
        ! TST rate law. The basic purpose of this option is to
        ! provide a template for a file on which the reactants
        ! react according to specified TST rate laws. It is intended
        ! that the user then modify this template to describe the
        ! problem actually desired.
        utitl1(12) = '  This sample file has Albite dissolving according to a TST rate law.'
        utitl1(13) = 'It may or may not actually run. For example, Albite may not appear'
        utitl1(14) = 'on the supporting data file.'
        utitl1(16) = '  The kinetic data shown here are taken from Knauss and Wolery (1986). You'
        utitl1(17) = 'may or may not wish to use these particular data.In any case, you are'
        utitl1(18) = 'entirely responsible for the data that you do use.'
        utitl1(20) = '                            References'
        utitl1(22) = 'Knauss, K.G., and Wolery, T.J., 1986, Dependence of albite dissolution'
        utitl1(23) = '  kinetics on pH and time at 25C and 70C, Geochimica et Cosmochimica Acta,'

        utitl1(24) = '  v. 50, p. 2481-2497.'
        nrct = 1
        ureac(1) = 'Albite'
        jcode(1) = 0
        jreac(1) = 0
        morr(1) = 1.0
        modr(1) = 0.
        nsk(1) = 0
        sfcar(1) = 1.0e+4
        fkrc(1) = 1.0

        nrk(1,1) = 1
        nrk(2,1) = -1
        imech(1,1) = 3

        ! Acid range.
        rkb(1,1,1) = 1.00e-15
        trkb(1,1,1) = 25.0
        iact(1,1,1) = 2
        hact(1,1,1) = 28.47
        csigma(1,1,1) = 1.0
        ndact(1,1,1) = 1
        udac(1,1,1,1) = 'H+'
        cdac(1,1,1,1) = 1.0

        ! Neutral range.
        rkb(2,1,1) = 3.98e-17
        trkb(2,1,1) = 25.0
        iact(2,1,1) = 2
        hact(2,1,1) = 12.89
        csigma(2,1,1) = 1.0
        ndact(2,1,1) = 0

        ! Alkaline range.
        rkb(3,1,1) = 5.01e-21
        trkb(3,1,1) = 25.0
        iact(3,1,1) = 2
        hact(3,1,1) = 7.68
        csigma(3,1,1) = 1.0
        ndact(3,1,1) = 1
        udac(1,3,1,1) = 'H+'
        cdac(1,3,1,1) = -0.5
    else if (iopt(19) .eq. 3) then
        ! Write an EQ6 input file for a fluid mixing calculation in
        ! which "Fluid 2" is mixed into "Fluid 1". Note that "Fluid 2"
        ! appears first on the EQ3NR input file, beflore "Fluid 1".
        ! It appears as a special reactant ont he top half of the EQ6
        ! input file. "Fluid 1" is the starting aqueous solution in the
        ! equilibrium system, and appears in the bottom half. If this
        ! option is chosen for the first and only problem on the EQ3NR
        ! input file, the pickup file produced is set up for mixing of
        ! this fluid with itself. This would normally be useful only
        ! as a template.
        ! Set iopt(17) for the first fluid ("Fluid 2") to -1, so no
        ! pickup file is produced (iopt(19) can be set to 0). For the
        ! second fluid ("Fluid 1"), set iopt(17) equal to 0 and
        ! iopt(19) equal to 3.
        utitl1(12) = '  This sample file has "Fluid 2" mixing into "Fluid 1". "Fluid 2" was the'
        utitl1(13) = 'first aqueous solution on the EQ3NR input file, "Fluid 1" was the second.'
        utitl1(14) = 'The option switch iopt(19) was set to 3 for the latter to produce the'
        utitl1(15) = 'present EQ6 input file.'

        if (jtemp .eq. 0) then
            utitl1(16) = ' '
            utitl1(17) = '  Mixing will be isothermal as specified below.'
        else if (jtemp .eq. 3) then
            utitl1(16) = ' '
            utitl1(17) = '  Mixing will be non-isothermal as specified below.'
        end if

        if (nprob .eq. 1) then
            utitl1(18) = ' '
            utitl1(19) = '  Warning: This EQ6 input file is set up to mix a fluid with itself.'
        end if

        nrct = 1
        nsrt = 1
        ureac(1) = ureac1
        jcode(1) = 2
        jreac(1) = 0
        morr(1) = 1.0
        modr(1) = 0.
        nrk(1,1) = 1
        rkb(1,1,1) = 1.0

        vreac(1) = 1.0

        iesrti(1) = iesrt1

        do n = 1,iesrt1
            uesri(n,1) = uesr1(n)
            cesri(n,1) = cesr1(n)
        end do

        ibsrti(1) = ibsrt1

        do n = 1,ibsrt1
            ubsri(n,1) = ubsr1(n)
            cbsri(n,1) = cbsr1(n)
        end do
    else
        ! Unknown option switch value.
        write (ux8,"(i5)") iopt(19)
        call lejust(ux8)
        j2 = ilnobl(ux8)
        write (noutpt,1000) ux8(1:j2)
        write (nttyo,1000) ux8(1:j2)
1000 format(/' * Error - (EQ3NR/stpk36) Programming error trap:',' The "Advanced EQ3NR',/7x,'PICKUP File Options" switch',' iopt(19) has an unknown value of ',a,'.')

        stop
    end if

999 continue
end subroutine stpk36
