subroutine rd3inw(cgexj,cgexpi,covali,ehi,egexsi,fo2lgi,iebal3,ietmax,iodb,iopg,iopr,iopt,irdxc3,itdsf3,itermx,jgext,jetmax,jflgi,jgexti,jpres3,kxmod,mwtges,nbti,nbtmax,ncmpri,net,neti,netmax,ngexti,ninpts,ngexrt,nobswt,nodbmx,nopgmx,noprmx,noptmx,noutpt,nprob,nsbswt,ntitl,ntitmx,nttyo,nxicmx,nxmdmx,nxmod,nxti,nxtimx,pei,press,qend,qgexsh,qrderr,rho,scamas,tdspkg,tdspl,tempc,tgexp,tolbt,toldl,tolspf,ucospi,uebal,ugexj,ugexji,ugexmo,ugexp,ugexpi,ugexr,ugexsi,umemi,uobsw,uredox,usbsw,usoli,uspeci,utitl,uhfgex,uvfgex,uxkgex,uxmod,xbari,xgexsi,xhfgex,xlkgex,xvfgex,xlkmod,zgexj)
    !! This subroutine reads the EQ3NR input file in compact ("W")
    !! format.
    !! This subroutine is a near-clone of XCON3/rd3w8.f. The present
    !! subroutine differs from the latter, in that in addition to the
    !! pure read function, it writes an instant echo of what is read to
    !! the output file. To maintain close consistency with XCON3/rd3w8.f,
    !! this subroutine assigns no default values and only performs such
    !! checking of what is read to ensure that what follows is readable.
    !! The calling sequence of this subroutine is identical to that of
    !! EQ3NR/rd3ind.f, XCON3/rd3w8.f, and XCON3/rd3d8.f.
    !! This subroutine is called by:
    !!   EQ3NR/eq3nr.f
    !! Principal input:
    !!   ninpts = unit number of the stripped input file
    !!   noutpt = unit number of the output file
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
    integer :: j4
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

    ! qend   = .true if the end of the input file has been encountered
    ! qrderr = .true if the current problem can't be read because of
    !            a read format error or a dimensional overflow
    qend = .false.
    qrderr = .false.

    ! Title.
    !   utitl1 = Title
    !   ntitl1 = number of lines in the title
    ! Note: if there are exactly ntitpa lines in the title, the title
    ! need not be terminated by an 'endit.'. The 'endit.' if present
    ! is here not considered to be part of the title.
    n = 0
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

    if (uline(1:8) .eq. uendit(1:8)) then
        write (noutpt,1030)
        write (nttyo,1030)
1030 format(3x,'This input file has no title.')

        go to 120
    else
        n = 1
        utitl(n) = uline
    end if

    do nn = 2,ntitmx + 1
        read (ninpts,1000,err=990) uline
        j2 = ilnobl(uline)
        j2 = min(j2,79)
        write (noutpt,1020) uline(1:j2)

        if (uline(1:8) .eq. uendit(1:8)) then
            go to 120
        end if

        n = n + 1

        if (n .gt. ntitmx) then
            write (nttyo,1015) ntitmx
1015 format(/' * Error - (EQ3NR/rd3inw) Have too many lines in',/7x,'the title. The code is only dimensioned for ',i4,/7x,'lines. Reduce the size of the title or increase the',/7x,'dimensioning parameter ntitpa.')

            go to 990
        end if

        utitl(n) = uline
    end do

120 continue
    ntitl = n

    if (ntitl .gt. 0) then
        ! Write the first 5 lines of the input problem title to the
        ! screen file.
        write (nttyo,1040)
1040 format(3x,'The input problem title is (first 5 lines',' maximum):',/)

        i = min(ntitl,5)

        do n = 1,i
            j2 = ilnobl(utitl(n))
            j2 = min(j2,74)
            write (nttyo,1050) utitl(n)(1:j2)
1050 format(5x,a)
        end do
    end if

    write (nttyo,1060)
1060 format(/3x,'Continuing to read the problem input ...')

    ! Special basis switches.
    !   nsbswt = the number of special basis switches
    !   usbsw(1,n) = the species to be switched from the strict
    !                  basis to the auxiliary basis set in the
    !                  in the n-th switch
    !   usbsw(2,n) = the species to be switched from the auxiliary
    !                  basis set to the strict basis set in the
    !                  n-th switch
    write (noutpt,1100)
1100 format(' * Special basis switches')

    read (ninpts,1110,err=990) nsbswt
1110 format(12x,i3)

    write (noutpt,1120) nsbswt
1120 format(5x,'nsbswt= ',i3)

    do n = 1,nsbswt
        read (ninpts,1130,err=990) usbsw(1,n)
1130 format(9x,a48)

        j2 = ilnobl(usbsw(1,n))
        write (noutpt,1140) usbsw(1,n)(1:j2)
1140 format(1x,'species= ',a)

        read (ninpts,1150,err=990) usbsw(2,n)
1150 format(15x,a48)

        j2 = ilnobl(usbsw(2,n))
        write (noutpt,1160) usbsw(2,n)(1:j2)
1160 format(3x,'switch with= ',a)
    end do

    ! Temperature.
    !   tempc  = temperature, C
    write (noutpt,1165)
1165 format(' * General')

    read (ninpts,1170,err=990) tempc
1170 format(12x,e12.5)

    write (noutpt,1180) tempc
1180 format(6x,'tempc= ',1pe12.5)

    ! Pressure.
    !   jpres3 = pressure code:
    !     0 = Compute the pressure from the temperauture, using the
    !         data file's reference pressure curve:
    !           press = func(data file, T)
    !     1 = Compute the pressure from the temperature, using the
    !         1.013-bar/steam-saturation curve:
    !           press = built-in func(T)
    !     2 = Constant pressure:
    !           press = pressb
    !   press  = pressure, bars
    read (ninpts,1190,err=990) jpres3,press
1190 format(12x,i3,/12x,e12.5)

    write (noutpt,1200) jpres3,press
1200 format(5x,'jpres3= ',i3,/6x,'press= ',1pe12.5)

    ! Density.
    !   rho    = density of aqueous solution in g/ml
    read (ninpts,1220,err=990) rho
1220 format(12x,e12.5)

    write (noutpt,1222) rho
1222 format(8x,'rho= ',1pe12.5)

    ! Total dissolved solutes.
    !   itdsf3 = TDS code:
    !     0 = Use the input value in mg/kg.sol (tdspkg)
    !     1 = Use the input value in mg/L (tdspl)
    !   tdspkg = total dissolved solutes, mg/kg.sol
    !   tdspl  = total dissolved solutes, mg/L
    read (ninpts,1224,err=990) itdsf3,tdspkg,tdspl
1224 format(12x,i3,/12x,e12.5,12x,e12.5)

    write (noutpt,1227) itdsf3,tdspkg,tdspl
1227 format(5x,'itdsf3= ',i3,/5x,'tdspkg= ',1pe12.5,5x,'tdspl= ',e12.5)

    ! Species to adjust for electrical balance.
    !   iebal3 = electrical balancing code:
    !     0 = No balancing
    !     1 = Balance on the species whose name is entered as uredox
    !   uebal = name of a basis species whose concentration is to be
    !             adjusted to achieve electrical balance
    read (ninpts,1500,err=990) iebal3,uebal
1500 format(12x,i3,/12x,a24)

    j2 = ilnobl(uebal)
    write (noutpt,1510) iebal3,uebal(1:j2)
1510 format(5x,'iebal3= ',i3/6x,'uebal= ',a)

    ! Default redox constraint.
    !   irdxc3 = default redox option switch:
    !     -3 = Read a constraint on aqueous "O2(g)" in the regular
    !            basis species constraint block
    !     -2 = The pe is fixed at the value entered as "pei"
    !     -1 = The Eh (V) is fixed at the value entered as "ehi"
    !      0 = The log fO2 is fixed at the value entered as "fo2lgi"
    !      1 = The log Fo2, Eh, pe, etc., is determined by a redox
    !            couple, one member of which is defined by the name
    !            of an auxiliary basis species entered as "uredox",
    !            the member other of which is defined by its appearance
    !            in the reaction for the former. Usually Fe++ is a
    !            strict basis species, and Fe++ is in the auxiliary
    !            basis set, so entering uredox = "Fe+++" specifies
    !            the Fe+++/Fe++ couple. Concentration or other data
    !            must be entered below for both members of the couple.
    !            Couples involving non-aqueous phases can be defined
    !            indirectly, using the jflgi = 25 option described
    !            below. For example, to use the Hematite/Magnetite
    !            couple, enter uredox = "Fe+++". In the basis species
    !            constraint block, for Fe++, enter jflgi = 25 and
    !            ucospi = "Magnetite", and for Fe+++, enter jflgi = 25
    !            and ucospi = "Hematite".
    !   fo2lgi = log fO2
    !   ehi    = Eh, oxidation-reduction parameter
    !   pei    = pe, electron activity
    !   uredox = name (all 24 letters) of an auxiliary basis species
    !              that is part of a redox couple that will be used
    !              to  distribute other redox couples
    !   Note: heterogeneous equilibria can be used also by setting
    !     irdxc3 = 0, fep = 0., and jflgi(O2(g)) = 25
    read (ninpts,1232,err=990) irdxc3
1232 format(12x,i3)

    write (noutpt,1234) irdxc3
1234 format(5x,'irdxc3= ',i3)

    read (ninpts,1240,err=990) fo2lgi,ehi
1240 format(12x,e12.5,12x,e12.5)

    write (noutpt,1242) fo2lgi,ehi
1242 format(5x,'fo2lgi= ',1pe12.5,7x,'ehi= ',e12.5)

    read (ninpts,1250,err=990) pei,uredox
1250 format(12x,e12.5,12x,a24)

    j2 = ilnobl(uredox)
    write (noutpt,1252) pei,uredox(1:j2)
1252 format(8x,'pei= ',1pe12.5,4x,'uredox= ',a)

    ! Aqueous basis species and associated constraints.
    !   uspeci = name of an aqueous species in the basis set. If
    !              the species is not defined as a basis species on
    !              the input file, it will be treated as a basis species
    !              in the current run. The name is constructed in
    !              species/phase format (the first 24 characters
    !              comprise the species part,the second 24 characters,
    !              the name of the corresponding phase). Only the
    !              species part is input; the code provides the
    !              phase part.
    !   jflgi  = flag variable defining the type of constraint placed
    !              on a basis species. It determines the meaning of the
    !              "covali" parameter.
    !   covali = a floating point datum (concentration, activity, etc.)
    !              whose meaning is specified by the jflgi parameter.
    !   ucospi = the name of a species/phase which must be designated
    !              to complete a specification for the jflgi = 17, 18,
    !              or 25 options. For jflgi = 17 or 18, the species is
    !              an aqueous ion. For jflgi = 25, the species
    !              specified must belong to a phase other than the
    !              aqueous solution. The name is constructed in
    !              species/phase format. For the jflgi = 17 or 18
    !              options, the phase part may be left blank, as the
    !              phase by implication is 'aqueous solution.'
    !              For the jflgi = 25 option, the phase part may also
    !              be left blank. The phase will then be taken as that
    !              which first produces a name match on the species
    !              part. If the phase is intended to be a solid or
    !              liquid solution, failure to specify the phase part
    !              will result in matching the pure component phase
    !              instead.
    !   Summary of the jflgi conventions:
    !     -1 = Not present in the system. This value is assigned
    !            as a default to strict basis species.
    !      0 = Total concentration, molal ("covali")
    !            Note: enter jflgi = 0 for water and O2(g). The
    !            "covali" parameter for O2(g) is always log fugacity.
    !            A line for O2(g) is not entered here unless
    !            irdxc3 = -3.
    !      1 = Total concentration, molar ("covali")
    !      2 = Total concentration, mg/L ("covali")
    !      3 = Total concentration, mg/kg.sol ("covali")
    !      7 = Total alkalinity, eq/kg.H2O ("covali")
    !      8 = Total alkalinity, eq/L ("covali")
    !      9 = Total alkalinity, eq/kg.sol ("covali")
    !     10 = Total alkalinity, mg/L CaCO3 ("covali")
    !     11 = Total alkalinity, mg/L HCO3- ("covali")
    !     16 = Log activity (molal scale) ("covali")
    !     17 = The function abs(zj) log ai +/- abs(zi) log aj,
    !            where i and j are the ions specified by "uspeci"
    !            and "ucospi". The sign of the second term is positive
    !            for cation-anion and anion-cation combinations,
    !            negative for like-charge combinations. To input the
    !            pHCl defined as pH + pCl, enter H+ as "uspeci", -pHCl
    !            for "covali", and Cl- for uphas1.
    !     18 = Log mean activity of the neutral electrolyte composed
    !            of the ions specified by "uspeci" and "ucospi".
    !     19 = pX, or -log activity. This is very similar to the
    !            jflgi= 16 option, except that the "covali" input
    !            has the opposite sign. This can be used to enter
    !            pH, but see jflag= 20 below.
    !     20 = pH. This is very similar to the jflgi= 19 (pX) option.
    !            It can only be applied to H+, however.
    !     21 = pHCl. This is very similar to the jflgi= 17 (log
    !            activity combination) option. It can only be applied
    !            to H+ or Cl-, however. The other ion is then the
    !            counterion. Note that the counterion is not specified
    !            in the "ucospi" input.
    !     25 = Heterogeneous equilibrium. Here "usp" contains the name
    !            of the relevant species (first 24 characters) and
    !            phase (last 24 characters). If the relevant species
    !            is a gas species, its fugacity is specified by the
    !            "covali" input. This is a dangerous option.
    !     27 = Homogeneous equilibrium - the species must be in the
    !            auxiliary basis when the calculation commences;
    !            it and its own dependent species are treated as
    !            depenent species of the corresponding species in the
    !            strict basis, but are not limited by the total
    !            concentration, if one is specified, for that
    !            corresponding species. No "covali" input is required.
    !            This is a dangerous option.
    !            For example, suppose you have a sulfate analysis
    !            of a water and you expect that sulfide species
    !            may be significant and possibly in equilibrium with
    !            sulfate. If you think the sulfate analysis
    !            measured both sulfate and sulfide, enter the
    !            analytical value as the total concentration of
    !            sulfate and specify jflgi(HS-) = 30. If you think
    !            the analysis measured only sulfate proper, then
    !            enter jflgi(HS-) = 27 instead. If the redox state of
    !            the water is very reducing you may get an infinite
    !            or near-infinite amount of calculated sulfide.
    !     16 = Log activity (molal scale) ("covali")
    !     30 = Make non-basis. The species is treated as above,
    !            but it and its own dependent species are constrained
    !            by the total concentration entered for the
    !            corresponding strict basis species. No "covali" input
    !            is required.
    write (noutpt,2490)
2490 format(' * Aqueous basis species')

    nbi = 0

    do nn = 1,nbtmax + 1
        read (ninpts,2500,err=990) ux8,ux48
2500 format(a8,1x,a48)

        j2 = ilnobl(ux48)

        if (ux8(1:8) .eq. uendit(1:8)) then
            j3 = ilnobl(uendit)
            write (noutpt,2510) uendit(1:j3)
2510 format(1x,a)

            go to 400
        end if

        nbi = nbi + 1

        if (nbi .gt. nbtmax) then
            write (noutpt,2520) nbtmax,ux48(1:j2)
            write (nttyo,2520) nbtmax,ux48(1:j2)
2520 format(/' * Error - (EQ3NR/rd3inw) The number of basis',' species read',/7x,'from the input file exceeded the',' maximum of ',i3,' while',/7x,'trying to read data for',' the species ',a,'.',/7x,'Increase the dimensioning',' parameter nbtpar.')

            go to 990
        end if

        write (noutpt,2530) ux48(1:j2)
2530 format(1x,'species= ',a)

        uspeci(nbi) = ux48

        read (ninpts,2540,err=990) jflgi(nbi),covali(nbi)
2540 format(10x,i2,12x,e12.5)

        write (noutpt,2550) jflgi(nbi),covali(nbi)
2550 format(4x,'jflgi= ',i2,4x,'covali= ',1pe12.5)

        if (jflgi(nbi).eq.17 .or. jflgi(nbi).eq.18  .or. jflgi(nbi).eq.25) then
            read (ninpts,2570,err=990) ucospi(nbi)
2570 format(10x,a48)

            j3 = ilnobl(ucospi(nbi))
            write (noutpt,2580) ucospi(nbi)(1:j3)
2580 format(3x,'ucospi= ',a)
        end if
    end do

400 continue
    nbti = nbi

    ! Ion exchanger creation. This section reads the directives for
    ! creating ion exchange phases and species and their associated
    ! intrinsic (including thermodynamic) properties. The data is
    ! essentially that which might be read from a supporting data file.
    ! Scenario-specific data (e.g., a specified amount of exchanger or
    ! a specific composition for an exchanger phase) are not included
    ! in this section.
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
    !                ' '       = DeltaH0r, kcal, per equivalent
    !                'kcal/eq' = DeltaH0r, kcal, per equivalent
    !                'kJ/eq'   = DeltaH0r, kJ, per equivalent
    !   xvfgex = array of volume of reaction data for the reaction
    !              specified in the above string
    !   uvfgex = array of strings denoting the kind of data in the
    !              corresponding entry of the xhfgex array:
    !                ' '      = DeltaV0r, cm3 per equivalent
    !                'cm3/eq' = DeltaV0r, cm3 per equivalent
    write (noutpt,1590)
1590 format(' * Ion exchangers')

    read (ninpts,1595) ux8
1595 format(12x,a8)

    call lejust(ux8)
    call locase(ux8)
    qgexsh = ux8(1:2).eq.'.t' .or. ux8(1:1).eq. 't'
    write (noutpt,1597)
1597 format(5x,'qgexsh= ',l8)

    read (ninpts,1600,err=990) net
1600 format(12x,i3)

    write (noutpt,1610) net
1610 format(8x,'net= ',i3)

    if (net .gt. netmax) then
        write (noutpt,1620) netmax,net
        write (nttyo,1620) netmax,net
1620 format(/' * Error - (EQ3NR/rd3inw) Have exceeded the maximum',' number of ',i3,/7x,'generic ion exchange phases while',' reading the data to create',/7x,'such phases. Increase',' the dimensioning parameter netpar',/7x,'to at least ',i3,'.')

        go to 990
    end if

    do ne = 1,net
        read (ninpts,1630,err=990) ugexp(ne)
1630 format(12x,a24)

        j2 = ilnobl(ugexp(ne))
        write (noutpt,1640) ugexp(ne)(1:j2)
1640 format(6x,'ugexp= ',a)

        read (ninpts,1650,err=990) mwtges(ne)
1650 format(12x,e12.5)

        write (noutpt,1660) mwtges(ne)
1660 format(5x,'mwtges= ',1pe12.5)

        read (ninpts,1630,err=990) ugexmo(ne)
        j3 = ilnobl(ugexmo(ne))
        write (noutpt,1670) ugexmo(ne)(1:j3)
1670 format(5x,'ugexmo= ',a)

        read (ninpts,1650,err=990) tgexp(ne)
        write (noutpt,1680) tgexp(ne)
1680 format(6x,'tgexp= ',1pe12.5)

        read (ninpts,1690,err=990) jgext(ne)
1690 format(12x,i3)

        write (noutpt,1700) jgext(ne)
1700 format(6x,'jgext= ',i3)

        if (jgext(ne) .gt. jetmax) then
            write (noutpt,1710) jetmax,ugexp(ne)(1:j2),jgext(ne)
            write (nttyo,1710) jetmax,ugexp(ne)(1:j2),jgext(ne)
1710 format(/' * Error - (EQ3NR/rd3inw) Have exceeded the maximum',' number of ',i3,/7x,'exchange sites on a generic ion',' exchange phase while reading',/7x,'the data to create',a,'. Increase the',/7x,'dimensioning parameter jetpar to',' at least ',i3,'.')

            go to 990
        end if

        do je = 1,jgext(ne)
            read (ninpts,1730,err=990) ugexj(je,ne)
1730 format(12x,a8)

            j3 = ilnobl(ugexj(je,ne))
            write (noutpt,1740) ugexj(je,ne)(1:j3)
1740 format(6x,'ugexj= ',a)

            read (ninpts,1750,err=990) cgexj(je,ne),zgexj(je,ne)
1750 format(12x,e12.5,12x,e12.5)

            write (noutpt,1760) cgexj(je,ne),zgexj(je,ne)
1760 format(6x,'cgexj= ',1pe12.5,5x,'zgexj= ',e12.5)

            read (ninpts,1690,err=990) ngexrt(je,ne)
            write (noutpt,1810) ngexrt(je,ne)
1810 format(5x,'ngexrt= ',i3)

            if (ngexrt(je,ne) .gt. ietmax) then
                write (noutpt,1820) netmax,ugexj(je,ne)(1:j3),ugexp(ne)(1:j2),ngexrt(je,ne)
                write (nttyo,1820) netmax,ugexj(je,ne)(1:j3),ugexp(ne)(1:j2),ngexrt(je,ne)
1820 format(/' * Error - (EQ3NR/rd3inw) Have exceeded the',' maximum number of ',i3,/7x,'reactions for a site',' belonging to a generic ion exchange',/7x,'phase while',' reading the data for site ',a,' of exchange phase',/7x,a,'. Increase the dimensioning parameter',/7x,'ietpar to at least ',i3,'.')

                go to 990
            end if

            do n = 1,ngexrt(je,ne)
                read (ninpts,1830,err=990) ugexr(n,je,ne)
1830 format(12x,a56)

                j2 = ilnobl(ugexr(n,je,ne))
                write (noutpt,1840) ugexr(n,je,ne)(1:j2)
1840 format(6x,'ugexr= ',a)

                read (ninpts,1850,err=990) xlkgex(n,je,ne),uxkgex(n,je,ne)
1850 format(12x,e12.5,12x,a8)

                j4 = ilnobl(uxkgex(n,je,ne))
                write (noutpt,1860) xlkgex(n,je,ne),uxkgex(n,je,ne)(1:j4)
1860 format(5x,'xlkgex= ',1pe12.5,5x,'units= ',a)

                read (ninpts,1850,err=990) xhfgex(n,je,ne),uhfgex(n,je,ne)
                j4 = ilnobl(uhfgex(n,je,ne))
                write (noutpt,1870) xhfgex(n,je,ne),uhfgex(n,je,ne)(1:j4)
1870 format(5x,'xhfgex= ',1pe12.5,5x,'units= ',a)

                read (ninpts,1850,err=990) xvfgex(n,je,ne),uvfgex(n,je,ne)
                j4 = ilnobl(uvfgex(n,je,ne))
                write (noutpt,1880) xvfgex(n,je,ne),uvfgex(n,je,ne)(1:j4)
1880 format(5x,'xvfgex= ',1pe12.5,5x,'units= ',a)
            end do
        end do
    end do

    ! Specified generic ion exchanger compositions.
    !   neti   = the number of generic exchangers for which compositions
    !              are entered
    !   ugexpi = array of names of generic ion exchange phases for
    !              which compositions are entered
    !   cgexpi = concentration (moles/kg H2O) of generic ion exchange
    !              phases for which compositions are entered
    !   jgexti  = array of numbers of sites for the various generic
    !              exchanger phases
    !   ugexji = array of names of exchange sites of generic ion
    !              exchange phases
    !   ngexti  = array of numbers of exchange species for sites of
    !               the various generic exchanger phasess
    !   ugexsi = array of names of species for sites of generic ion
    !              exchange phases
    !   egexsi = array of equivalent fractions of ions on exchange sites
    !              in an ion exchange phase
    !   xgexsi = array of moles fractions of ions on exchange sites
    !              in an ion exchange phase
    write (noutpt,1900)
1900 format(' * Ion exchanger compositions')

    read (ninpts,1910,err=990) neti
1910 format(12x,i3)

    write (noutpt,1920) neti
1920 format(7x,'neti= ',i3)

    if (neti .le. 0) then
        go to 250
    end if

    if (neti .gt. netmax) then
        write (noutpt,1930) netmax,neti
        write (nttyo,1930) netmax,neti
1930 format(/' * Error - (EQ3NR/rd3inw) Have exceeded the maximum',' number of ',i3,/7x,'generic ion exchange phases while',' reading the data for concentrations',7x,'and compositions',' of such phases. Increase the dimensioning',/7x,'parameter',' netpar to at least ',i3,'.')

        go to 990
    end if

    do nei = 1,neti
        read (ninpts,1940,err=990) ugexpi(nei)
1940 format(12x,a24)

        j2 = ilnobl(ugexpi(nei))
        write (noutpt,1950) ugexpi(nei)(1:j2)
1950 format(5x,'ugexpi= ',a)

        read (ninpts,1960,err=990) cgexpi(nei)
1960 format(12x,e12.5)

        write (noutpt,1970) cgexpi(nei)
1970 format(5x,'cgexpi= ',1pe12.5)

        read (ninpts,1910,err=990) jgexti(nei)
        write (noutpt,1990) jgexti(nei)
1990 format(5x,'jgexti= ',i3)

        if (jgexti(nei) .gt. jetmax) then
            write (noutpt,2000) jetmax,ugexpi(nei)(1:j2),jgexti(nei)
            write (nttyo,2000) jetmax,ugexpi(nei)(1:j2),jgexti(nei)
2000 format(/' * Error - (EQ3NR/rd3inw) Have exceeded the maximum',' number of ',i3,/7x,'exchange sites on a generic ion',' exchange phase while reading',/7x,'the data for the',' concentration and composition of',/7x,a,'. Increase the',' dimensioning parameter jetpar',/7x,'to at least ',i3,'.')

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

        write (noutpt,2002) ugexpi(nei)(1:j3)
        write (nttyo,2002) ugexpi(nei)(1:j3)
2002 format(/' * Error - (EQ3NR/rd3inw) Data are present on the',' input file for the',/7x,'generic ion exchange phase "',a,', but this phase',/7x,"hasn't been previously defined."," Can't determine the requisite model,",/7x,"therefore can't",' finish reading the current data for this phase.')

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
            write (noutpt,2020) ugexji(jei,nei)(1:j3)
2020 format(5x,'ugexji= ',a)

            read (ninpts,1910,err=990) ngexti(jei,nei)
            write (noutpt,2030) ngexti(jei,nei)
2030 format(5x,'ngexti= ',i3)

            if (ngexti(jei,nei) .gt. ietmax) then
                write (noutpt,2040) ietmax,ugexji(jei,nei),ugexpi(nei),ngexti(jei,nei)
                write (nttyo,2040) ietmax,ugexji(jei,nei),ugexpi(nei),ngexti(jei,nei)
2040 format(/' * Error - (EQ3NR/rd3inw) Have exceeded the',' maximum number of ',i3,/7x,'species on an exchange',' site of a generic ion exchange phase',/7x,'while reading',/7x,'the data for the composition of site ',a,' of',/7x,'exchange phase ',a,'. Increase the dimensioning',' parameter ietpar to at least ',i3,'.')

                go to 990
            end if

            if (qgexef) then
                ! The data are in equivalent fractions.
                do iei = 1,ngexti(jei,nei)
                    read (ninpts,2050,err=990) ugexsi(iei,jei,nei),egexsi(iei,jei,nei)
2050 format(12x,a24,10x,e12.5)

                    write (noutpt,2052) ugexsi(iei,jei,nei),egexsi(iei,jei,nei)
2052 format(5x,'ugexsi= ',a24,2x,'egexsi= ',1pe12.5)
                end do
            else
                ! The data are in mole fractions.
                do iei = 1,ngexti(jei,nei)
                    read (ninpts,2050,err=990) ugexsi(iei,jei,nei),xgexsi(iei,jei,nei)
                    write (noutpt,2060) ugexsi(iei,jei,nei),xgexsi(iei,jei,nei)
2060 format(5x,'ugexsi= ',a24,2x,'xgexsi= ',1pe12.5)
                end do
            end if
        end do
    end do

250 continue

    ! Specified solid solution compositions.
    !   nxti   = the number of solid solutions for which compositions
    !              are entered
    !   usoli  = array of names of solid solutions for which
    !              compositions are entered
    !   ncmpri = range pointer array for solid solutions for which
    !              compositions are entered; defines the range in
    !              the umemi and xbari arrays which corresponds to the
    !              end-member components of a given solid solution
    !   nssoti = number of end-member components in a given solid
    !              solution
    !   umemi  = array of names of end-member components of solid
    !              solutions for which compositions are entered
    !   xbari  = array of moles fractions of end-member components
    !              of solid solutions for which compositions are entered
    write (noutpt,2100)
2100 format(' * Solid solution compositions')

    read (ninpts,2102,err=990) nxti
2102 format(12x,i3)

    write (noutpt,2104) nxti
2104 format(7x,'nxti= ',i3)

    if (nxti .gt. nxtimx) then
        write (noutpt,2107) nxtimx,nxti
        write (nttyo,2107) nxtimx,nxti
2107 format(/' * Error - (EQ3NR/rd3inw) Have exceeded the maximum',/7x,'number of ',i5,' solid solutions for which',' compositions',/7x,'may be specified on the input file.',' Increase the dimensioning',/7x,'parameter nxtipa to at',' least ',i3,'.')

        go to 990
    end if

    nxic = 0

    do nxi = 1,nxti
        read (ninpts,2120,err=990) usoli(nxi)
2120 format(12x,a24)

        j2 = ilnobl(usoli(nxi))
        write (noutpt,2130) usoli(nxi)(1:j2)
2130 format(6x,'usoli= ',a24)

        read (ninpts,2140,err=990) nssoti
2140 format(12x,i3)

        write (noutpt,2142) nssoti
2142 format(5x,'nssoti= ',i3)

        ncmpri(1,nxi) = nxic + 1
        ncmpri(2,nxi) = nxic + nssoti
        nr1 = ncmpri(1,nxi)
        nr2 = ncmpri(2,nxi)

        if (nr2 .gt. nxicmx) then
            write (noutpt,2150) nxicmx,usoli(nxi)(1:j2),nr2
            write (nttyo,2150) nxicmx,usoli(nxi)(1:j2),nr2
2150 format(/' * Error - (EQ3NR/rd3inw) Have exceeded the',' maximum number',/7x,'of ',i5,' solid solution species',' for which mole fractions',/7x,'are specified on the',' input file. This occurred while reading ',/7x,'data',' for the solid solution ',a,'. Increase',/7x,'the',' dimensioning parameter nxicpa to at least ',i3,'.')

            go to 990
        end if

        do nxic = nr1,nr2
            read (ninpts,2160,err=990) umemi(nxic),xbari(nxic)
2160 format(12x,a24,12x,e12.5)

            write (noutpt,2170) umemi(nxic),xbari(nxic)
2170 format(6x,'umemi= ',a24,5x,'xbari= ',1pe12.5)
        end do
    end do

    ! Number of nxmod options.
    !   nxmod = the number of suppressed/altered species/reactions
    write (noutpt,2330)
2330 format(' * Alter/suppress options')

    read (ninpts,2340,err=990) nxmod
2340 format(12x,i3)

    write (noutpt,2350) nxmod
2350 format(6x,'nxmod= ',i3)

    if (nxmod .gt. nxmdmx) then
        write (noutpt,2360) nxmdmx
        write (nttyo,2360) nxmdmx
2360 format(/' * Error - (EQ3NR/rd3inw) Have too many nxmod',/7x,'alter/suppress options. The code is only dimensioned',/7x,'for ',i3,' such options. Reduce the number of such',/7x,'options or increase the dimensioning parameter nxmdpa.')

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
        read (ninpts,2370,err=990) uxmod(n),kxmod(n),xlkmod(n)
2370 format(12x,a48,/12x,i2,22x,e12.5)

        j2 = ilnobl(uxmod(n))
        write (noutpt,2380) uxmod(n)(1:j2),kxmod(n),xlkmod(n)
2380 format(4x,'species= ',a,/5x,'option= ',i2,14x,'xlkmod= ',1pe12.5)
    end do

    ! Iopt model option switches.
    ! Note: iopt(1) = iopt1, etc.
    !   iopt(1) - iopt(3) = Used only by EQ6
    !   iopt(4) = Solid solution products:
    !      0 = Solid solutions are ignored
    !      1 = Solid solutions are permitted
    !      2 = Solid solutions are permitted and compositions to
    !            use in computing saturation indices are read
    !   iopt(5) - iopt(7) = Used only by EQ6
    !   iopt(8) = Not used
    !   iopt(9) - iopt(10) = Used only by EQ6
    !   iopt(11) = Auto basis switching, in pre-N-R optimization:
    !      0 = Off
    !      1 = On
    !   iopt(12) - iopt(13) Used only by EQ6
    !   iopt(14) = Not used
    !   iopt(15) - iopt(16) = Used only by EQ6
    !   iopt(17) = pickup file options:
    !     -1 = Don't write a pickup file
    !      0 = Write a pickup file at the end of the run
    !   iopt(18) = Used only by EQ6
    !   iopt(19) = Advanced EQ3NR pickup file options:
    !      0 = Write a normal EQ3NR pickup file
    !      1 = Write an EQ6 input file with Quartz dissolving, using
    !            a relative rate law
    !      2 = Write an EQ6 input file with Albite dissolving, using
    !            a TST rate law
    !      3 = Write an EQ6 input file to mix two fluids appearing
    !            on a single EQ3NR input file; the first fluid on that
    !            file is set up as a special reactant called "Fluid 2"
    !            that will be added to the second fluid on that
    !            file, which in effect becomes "Fluid 1". Set iopt(19)
    !            to 3 only for the second fluid.
    !   iopt(20) = Used only by EQ6
    write (noutpt,1365)
1365 format(' * Iopt, iopg, iopr, and iodb options')

    write (noutpt,1370)
1370 format(' *',15x,'1    2    3    4    5    6    7    8    9   10')

    read (ninpts,1380,err=990) (iopt(i), i = 1,20)
1380 format(12x,10i5)

    write (noutpt,1390) (iopt(i), i = 1,20)
1390 format(3x,'iopt1-10= ',10i5,/2x,'iopt11-20= ',10i5)

    ! Iopg activity coefficient option switches.
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
    !   iopg(3) = iopg(20) = Not used
    read (ninpts,1410,err=990) (iopg(i), i = 1,20)
1410 format(12x,10i5)

    write (noutpt,1420) (iopg(i), i = 1,20)
1420 format(3x,'iopg1-10= ',10i5,/2x,'iopg11-20= ',10i5)

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
    !   iopr(17) = pickup file format:
    !      0 = Use the same format ("D" or "W") as the input file
    !      1 = Use "W" format
    !      2 = Use "D" format
    !   iopr(18) - iopr(20) = Not used
    read (ninpts,1430,err=990) (iopr(i), i = 1,20)
1430 format(12x,10i5)

    write (noutpt,1440) (iopr(i), i = 1,20)
1440 format(3x,'iopr1-10= ',10i5,/2x,'iopr11-20= ',10i5)

    ! Iodb debugging print option switches.
    ! Note: iodb(1) = iodb1, etc.
    !   iodb(1) = General diagnostic messages:
    !      0 = Don't print
    !      1 = Print Level 1 diagnostic messages
    !      2 = Print Level 1 and Level 2 diagnostic messages
    !   iodb(2)  = Used only by EQ6
    !   iodb(3) = Pre-Newton-Raphson optimization:
    !      0 = Don't print
    !      1 = Print summary information
    !      2 = Print detailed information
    !      3 = Print more detailed information
    !      4 = Also print changes to activity coefficients
    !   iodb(4) = Newton-Raphson iteration:
    !      0 = Don't print
    !      1 = Print summary information
    !      2 = Print detailed information including the residual and
    !            correction vectors
    !      3 = Also print the Jacobian matrix
    !      4 = Also print changes to activity coefficients
    !   iodb(5)  = Used only by EQ6
    !   iodb(6) = Hypothetical affinity calculations (for solid
    !               solutions, etc.):
    !      0 = Don't print
    !      1 = Print summary information
    !      2 = Print detailed information
    !   iodb(7)  = Used only by EQ6
    !   iodb(8) -iodb(20) = Not used
    read (ninpts,1470,err=990) (iodb(i), i = 1,20)
1470 format(12x,10i5)

    write (noutpt,1480) (iodb(i), i = 1,20)
1480 format(3x,'iodb1-10= ',10i5,/2x,'iodb11-20= ',10i5)

    ! Numerical parameters.
    !   tolbt = convergence tolerance on residual magnitude
    !             (in Newton-Raphson iteration)
    !   toldl = convergence tolerance on correction magnitude
    !             (in Newton-Raphson iteration)
    !   itermx = maximum number of Newton-Raphson iterations
    write (noutpt,1270)
1270 format(' * Numerical parameters')

    read (ninpts,1280,err=990) tolbt,toldl
1280 format(12x,e12.5,12x,e12.5)

    write (noutpt,1290) tolbt,toldl
1290 format(6x,'tolbt= ',1pe12.5,5x,'toldl= ',e12.5)

    read (ninpts,1300,err=990) itermx
1300 format(12x,i3)

    write (noutpt,1310) itermx
1310 format(5x,'itermx= ',i3)

    ! Ordinary basis switches.
    !   nobswt = the number of ordinary basis switches
    !   uobsw(1,n) = the species to be switched from the basis
    !                  set in the n-th switch
    !   uobsw(2,n) = the species to be switched into the basis
    !                  set in the n-th switch
    write (noutpt,2700)
2700 format(' * Ordinary basis switches')

    read (ninpts,1110,err=990) nobswt
    write (noutpt,2710) nobswt
2710 format(5x,'nobswt= ',i3)

    do n = 1,nobswt
        read (ninpts,1130,err=990) uobsw(1,n)
        j2 = ilnobl(uobsw(1,n))
        write (noutpt,1140) uobsw(1,n)(1:j2)
        read (ninpts,1150,err=990) uobsw(2,n)
        j2 = ilnobl(uobsw(2,n))
        write (noutpt,1160) uobsw(2,n)(1:j2)
    end do

    ! Saturation flag tolerance.
    !   tolspf = saturation flag tolerance for minerals and
    !              various other phases, in kcal. This parameter
    !              controls the printing of the "SATD" and "SSATD"
    !              flags in the saturation state tables. It does not
    !              have any other effect in EQ3NR.
    write (noutpt,1292)
1292 format(' * Saturation flag tolerance')

    read (ninpts,1294,err=990) tolspf
1294 format(12x,e12.5)

    write (noutpt,1296) tolspf
1296 format(5x,'tolspf= ',1pe12.5)

    ! Scale factor for the mass of aqueous solution to write
    ! on the pickup file.
    !   scamas = scale factor (the default is 1.0, which corresponds
    !              to a mass of aqueous solution containing 1 kg
    !              of solvent water)
    write (noutpt,1315)
1315 format(' * Aqueous phase scale factor')

    read (ninpts,1320,err=990) scamas
1320 format(12x,e12.5)

    write (noutpt,1330) scamas
1330 format(5x,'scamas= ',1pe12.5)

    write (noutpt,3000) nprob
    write (nttyo,3000) nprob
3000 format(/'   Done reading problem ',i3,'.',/)

    go to 999

990 continue
    qrderr = .true.

999 continue
end subroutine rd3inw