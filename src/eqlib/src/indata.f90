subroutine indata(aadh,aadhh,aadhv,aaphi,abdh,abdhh,abdhv,abdot,abdoth,abdotv,adadhh,adadhv,adbdhh,adbdhv,adbdth,adbdtv,adhfe,adhfsa,advfe,advfsa,amua,aprehw,apresg,apxa,aslma,atwta,axhfe,axhfsa,axlke,axlksa,axvfe,axvfsa,azeroa,bpxa,cco2,cdrsa,cessa,eps100,iapxa_asv,iapxta,iaqsla,ibpxa_asv,ibpxta,ielam,igas,ikta_asv,insgfa,ipbt_asv,ipch,ipch_asv,ipcv,ipcv_asv,ixrn1a,ixrn2a,jpdblo,jpfc_asv,jptffl,jsola,mwtspa,nad1,nalpaa,napa_asv,napta,narn1a,narn2a,narx_asv,narxt,nata,nata_asv,nbaspa,nbta,nbta_asv,nbtafd,nbta1_asv,ncmpra,ncta,ncta_asv,ndrsa,ndrsa_asv,ndrsra,nessa,nessa_asv,nessra,ngrn1a,ngrn2a,ngta,ngta_asv,nlrn1a,nlrn2a,nlta,nlta_asv,nmrn1a,nmrn2a,nmta,nmta_asv,nmuta,nmuta_asv,nmuxa,noutpt,npta,npta_asv,nslta,nslta_asv,nslxa,nsta,nsta_asv,ntid_asv,ntitld,ntpr_asv,ntprt,nttyo,nxrn1a,nxrn2a,nxta,nxta_asv,qclnsa,palpaa,tdamax,tdamin,tempcu,ubasp,udakey,udatfi,uelema,uspeca,uphasa,uptypa,utitld,vosp0a,zchara)
    !! This subroutine reads the supporting data file DATA1, starting
    !! with its title. This subroutine should be called after subroutine
    !! indath.f has been called to read the header section (which
    !! mostly contains array size allocation data required to allocate
    !! sufficient array space to store the data read by the present
    !! subroutine.
    !! This subroutine is called by:
    !!   EQ3NR/eq3nr.f
    !!   EQ6/eq6.f
    !! Principal input:
    !!   noutpt = unit number of the output file
    !!   nttyo  = unit number of the screen file
    !! Principal output:
    !!   aadh   = array of polynomial coefficients for computing
    !!              the Debye-Huckel A(gamma,10) parameter
    !!   aadhh  = array of polynomial coefficients for computing
    !!              the Debye-Huckel A(H) parameter
    !!   aadhv  = array of polynomial coefficients for computing
    !!              the Debye-Huckel A(V) parameter
    !!   aaphi  = array of polynomial coefficients for computing
    !!              the Debye-Huckel A(phi) parameter
    !!   abdh   = array of polynomial coefficients for computing
    !!              the Debye-Huckel B parameter
    !!   abdhh  = array of polynomial coefficients for computing
    !!              the Debye-Huckel B(H) parameter
    !!   abdhv  = array of polynomial coefficients for computing
    !!              the Debye-Huckel B(V) parameter
    !!   abdot  = array of polynomial coefficients for computing
    !!              the Helgeson (1969) B-dot parameter
    !!   abdoth = array of polynomial coefficients for computing
    !!              the B-dot(H) parameter
    !!   abdotv = array of polynomial coefficients for computing
    !!              the B-dot(V) parameter
    !!   adadhh = array of polynomial coefficients for computing
    !!              the the pressure derivatives of the Debye-Huckel
    !!              A(H) parameter
    !!   adadhv = array of polynomial coefficients for computing
    !!              the the pressure derivatives of the Debye-Huckel
    !!              A(V) parameter
    !!   adbdhh = array of polynomial coefficients for computing
    !!              the the pressure derivatives of the Debye-Huckel
    !!              B(H) parameter
    !!   adbdhv = array of polynomial coefficients for computing
    !!              the the pressure derivatives of the Debye-Huckel
    !!              B(V) parameter
    !!   adbdth = array of polynomial coefficients for computing
    !!              the the pressure derivatives of the B-dot(H)
    !!              parameter
    !!   adbdtv = array of polynomial coefficients for computing
    !!              the the pressure derivatives of the B-dot(V)
    !!              parameter
    !!   adhfsa = array of polynomial coefficients for computing
    !!              the pressure derivatives of the enthalpy of reaction
    !!              for reactions read from the data file
    !!   adhfe  = array of polynomial coefficients for computing
    !!              the pressure derivatives of the enthalpy of
    !!              reaction for the "Eh" reaction
    !!   advfe  = array of polynomial coefficients for computing
    !!              the pressure derivatives of the volume of
    !!              reaction for the "Eh" reaction
    !!   advfsa = array of polynomial coefficients for computing
    !!              the pressure derivatives of the volume of reaction
    !!              for reactions read from the data file
    !!   amua   = coefficients for calculating Pitzer mu parameters
    !!              as a function of temperature, read from the
    !!              data file
    !!   aprehw = array of polynomial coefficients for computing
    !!              the recommended pressure envelope half-width (bars)
    !!   apresg = array of polynomial coefficients for computing
    !!              the pressure (bars) on the 1-atm:steam saturation
    !!              curve
    !!   apxa   = array of interaction coefficients for computing
    !!               activity coefficients in solid solutions
    !!   aslm   = coefficients for calculating Pitzer S-lambda(n)
    !!              parameters as a function of temperature, read
    !!              read from the data file
    !!   atwta  = array of atomic weights of the elements read from
    !!              the data file
    !!   axhfe  = array of polynomial coefficients for computing
    !!              the enthalpy of reaction for the "Eh" reaction
    !!   axhfsa = array of polynomial coefficients for computing
    !!              the enthalpy of reaction for reactions read
    !!              from the data file
    !!   axlke  = array of polynomial coefficients for computing
    !!              the log K for the "Eh" reaction
    !!   axlksa = array of polynomial coefficients for computing
    !!              the log K functions of reactions read from the
    !!              data file
    !!   axvfe  = array of polynomial coefficients for computing
    !!              the volume of reaction for the "Eh" reaction
    !!   axvfsa = array of polynomial coefficients for computing
    !!              the volume of reaction for reactions read
    !!              from the data file
    !!   azeroa = array of hard core diameters of the species read
    !!              from the data file
    !!   bpxa   = array of site-mixing parameters for computing activity
    !!              coefficients in solid solutions
    !!   cco2   = array of coefficients for computing the activity
    !!              coefficient of CO2(aq) from the Drummond (1981)
    !!              equation
    !!   cdrsa  = array of reaction coefficients for the reactions
    !!              read from the data file
    !!   cessa  = array of elemental composition coefficients for the
    !!              species read from the data file
    !!   iapxta = array of the numbers of non-zero interaction
    !!              coefficients for computing activity coefficients in
    !!              solid solutions
    !!   iaqsla = index of the aqueous solution phase, as read from
    !!              the data file
    !!   ibpxta = array of the numbers of non-zero site-mixing
    !!              paramters for computing activity coefficients in
    !!              solid solutions
    !!   ielam  = flag indicating no use or use of E-lambda terms
    !!   igas   = index of the gas phase
    !!   insgfa = array of flags for treating the activity coefficients
    !!              of species that happen to be electrically neutral,
    !!              as read from the data file
    !!   ixrn1a = index of the first non-aqueous solution phase,
    !!              as read from the data file
    !!   ixrn2a = index of the last non-aqueous solution phase,
    !!              as read from the data file
    !!   jpdblo = Pitzer data block organization flag
    !!   jptffl = Pitzer parameter temperature function flag
    !!   jsola  = array of flags for the activity coefficient models
    !!              to use for solid solutions, as read from the
    !!              data file
    !!   mwtspa = array of the molecular weights of the species read
    !!              from the data file
    !!   nalpaa = pointer array giving the index of the set of alpha
    !!              coeffcients for a given set of S-lambda
    !!              coefficients, as read from the data file
    !!   narn1a = start of the range of aqueous species, as read from
    !!              the data file
    !!   narn1a = end of the range of aqueous species, as read from
    !!              the data file
    !!   nbaspa = array of indices of species in basis set, as read
    !!              from the data file
    !!   ncmpra = array giving the start and the end of the range of the
    !!              species belonging to a given phase, as read from the
    !!              data file
    !!   ndrsa  = array of indices of the species corresponding to
    !!              the reaction coefficients in the cdrsa array
    !!   ndrsra = array giving the start and end of the range in the
    !!              cdrsa and ndrsa arrays of the species in the reaction
    !!              which belongs to a given species
    !!   nessa  = array of indices of the chemical elements
    !!              corresponding to the composition coefficients in
    !!              the cessa array
    !!   nessra = array giving the start and end of the range in the
    !!              cessa and nessa arrays of the chemical elements
    !!              comprising a given species
    !!   ngrn1a = index of the first species in the gas range, as read
    !!              from the data file
    !!   ngrn2a = index of the last species in the gas range, as read
    !!              from the data file
    !!   nlrn1a = index of the first species in the pure liquid range,
    !!              as read from the data file
    !!   nlrn2a = index of the last species in the pure liquid range,
    !!              as read from the data file
    !!   nmrn1a = index of the first species in the pure mineral range,
    !!              as read from the data file
    !!   nmrn2a = index of the last species in the pure mineral range,
    !!              as read from the data file
    !!   nmuta  = the number of species triplets for which mu
    !!              coefficients are defined, as read from the data file
    !!   nmuxa  = indices of species in triplets for which mu
    !!              coefficients are defined, as read from the data file
    !!   nslta  = number of species pairs for which S-lambda
    !!              coefficients are defined, as read from the data file
    !!   nslxa  = indices of species in pairs for which S-lambda
    !!              coefficients are defined, as read from the data file
    !!   nxrn1a = index of the first species in the solid solution range,
    !!              as read from the data file
    !!   nxrn2a = index of the last species in the solid solution range,
    !!              as read from the data file
    !!   palpaa = array of sets of alpha coefficients read from the
    !!              data file
    !!   qclnsa = array of flags indicating if a species was created
    !!              by cloning
    !!   tdamax = the nominal upper temperature limit of the data file
    !!   tdamin = the nominal lower temperature limit of the data file
    !!   ubasp  = array of names of the species in the basis set
    !!   uelema = array of names of the chemical elements read from
    !!              the data file
    !!   uphasa = array of names of phases read from the data file
    !!   uptypa = array giving the type of a given phase (e.g.,
    !!              liquid, solid, gas) read from the data file
    !!   uspeca = array of names of species read from the data file
    !!   utitld = the title read from the data file
    !!   vosp0a = array of standard state molar volumes of the species
    !!              read from the data file
    !!   zchara = array of electrical charge numbers of the species
    !!              read from the data file
    implicit none

    ! Calling sequence variable declarations.
    integer :: iapxa_asv
    integer :: ibpxa_asv
    integer :: ikta_asv
    integer :: ipbt_asv
    integer :: ipch_asv
    integer :: ipcv_asv
    integer :: jpfc_asv
    integer :: napa_asv
    integer :: narx_asv
    integer :: nata_asv
    integer :: nbta_asv
    integer :: nbta1_asv
    integer :: ncta_asv
    integer :: ndrsa_asv
    integer :: nessa_asv
    integer :: ngta_asv
    integer :: nlta_asv
    integer :: nmta_asv
    integer :: nmuta_asv
    integer :: npta_asv
    integer :: nslta_asv
    integer :: nsta_asv
    integer :: ntid_asv
    integer :: ntpr_asv
    integer :: nxta_asv

    integer :: nad1
    integer :: noutpt
    integer :: nttyo

    integer :: iapxta(nxta_asv)
    integer :: ibpxta(nxta_asv)
    integer :: insgfa(nata_asv)
    integer :: jsola(nxta_asv)
    integer :: nalpaa(nslta_asv)
    integer :: narxt(ntpr_asv)
    integer :: nbaspa(nbta_asv)
    integer :: ncmpra(2,npta_asv)
    integer :: ndrsa(ndrsa_asv)
    integer :: ndrsra(2,nsta_asv)
    integer :: nessa(nessa_asv)
    integer :: nessra(2,nsta_asv)
    integer :: nmuxa(3,nmuta_asv)
    integer :: nslxa(2,nslta_asv)

    integer :: iaqsla
    integer :: ielam
    integer :: igas
    integer :: ipch
    integer :: ipcv
    integer :: ixrn1a
    integer :: ixrn2a
    integer :: jpdblo
    integer :: jptffl
    integer :: napta
    integer :: narn1a
    integer :: narn2a
    integer :: nata
    integer :: nbta
    integer :: nbtafd
    integer :: ncta
    integer :: ngrn1a
    integer :: ngrn2a
    integer :: ngta
    integer :: nlrn1a
    integer :: nlrn2a
    integer :: nlta
    integer :: nmrn1a
    integer :: nmrn2a
    integer :: nmta
    integer :: nmuta
    integer :: npta
    integer :: nslta
    integer :: nsta
    integer :: ntitld
    integer :: ntprt
    integer :: nxrn1a
    integer :: nxrn2a
    integer :: nxta

    logical :: qclnsa(nsta_asv)

    character(len=8) :: uelema(ncta_asv)
    character(len=24) :: uphasa(npta_asv)
    character(len=24) :: uptypa(npta_asv)
    character(len=48) :: ubasp(nbta_asv)
    character(len=48) :: uspeca(nsta_asv)
    character(len=80) :: utitld(ntid_asv)

    character(len=8) :: udakey
    character(len=8) :: udatfi

    real(kind=8) :: tempcu(ntpr_asv)

    real(kind=8) :: aadh(narx_asv,ntpr_asv)
    real(kind=8) :: aadhh(narx_asv,ntpr_asv)
    real(kind=8) :: aadhv(narx_asv,ntpr_asv)
    real(kind=8) :: aaphi(narx_asv,ntpr_asv)
    real(kind=8) :: abdh(narx_asv,ntpr_asv)
    real(kind=8) :: abdhh(narx_asv,ntpr_asv)
    real(kind=8) :: abdhv(narx_asv,ntpr_asv)
    real(kind=8) :: abdot(narx_asv,ntpr_asv)
    real(kind=8) :: abdoth(narx_asv,ntpr_asv)
    real(kind=8) :: abdotv(narx_asv,ntpr_asv)
    real(kind=8) :: adadhh(narx_asv,ntpr_asv,ipch_asv)
    real(kind=8) :: adadhv(narx_asv,ntpr_asv,ipcv_asv)
    real(kind=8) :: adbdhh(narx_asv,ntpr_asv,ipch_asv)
    real(kind=8) :: adbdhv(narx_asv,ntpr_asv,ipcv_asv)
    real(kind=8) :: adbdth(narx_asv,ntpr_asv,ipch_asv)
    real(kind=8) :: adbdtv(narx_asv,ntpr_asv,ipcv_asv)
    real(kind=8) :: adhfe(narx_asv,ntpr_asv,ipch_asv)
    real(kind=8) :: advfe(narx_asv,ntpr_asv,ipcv_asv)
    real(kind=8) :: adhfsa(narx_asv,ntpr_asv,ipch_asv,nsta_asv)
    real(kind=8) :: advfsa(narx_asv,ntpr_asv,ipcv_asv,nsta_asv)
    real(kind=8) :: aprehw(narx_asv,ntpr_asv)
    real(kind=8) :: apresg(narx_asv,ntpr_asv)
    real(kind=8) :: apxa(iapxa_asv,nxta_asv)
    real(kind=8) :: atwta(ncta_asv)
    real(kind=8) :: axhfe(narx_asv,ntpr_asv)
    real(kind=8) :: axhfsa(narx_asv,ntpr_asv,nsta_asv)
    real(kind=8) :: axlke(narx_asv,ntpr_asv)
    real(kind=8) :: axlksa(narx_asv,ntpr_asv,nsta_asv)
    real(kind=8) :: axvfe(narx_asv,ntpr_asv)
    real(kind=8) :: axvfsa(narx_asv,ntpr_asv,nsta_asv)

    real(kind=8) :: amua(jpfc_asv,nmuta_asv)
    real(kind=8) :: aslma(jpfc_asv,0:ipbt_asv,nslta_asv)
    real(kind=8) :: azeroa(nata_asv)
    real(kind=8) :: bpxa(ibpxa_asv,nxta_asv)
    real(kind=8) :: cco2(5)
    real(kind=8) :: cdrsa(ndrsa_asv)
    real(kind=8) :: cessa(nessa_asv)
    real(kind=8) :: mwtspa(nsta_asv)
    real(kind=8) :: palpaa(ipbt_asv,napa_asv)
    real(kind=8) :: vosp0a(nsta_asv)
    real(kind=8) :: zchara(nsta_asv)

    real(kind=8) :: eps100

    real(kind=8) :: tdamax
    real(kind=8) :: tdamin

    ! Local variables with global dimensioning.
    ! None of these variables need to be SAVEd.
    character(len=24), dimension(:), allocatable :: udrsv
    character(len=8), dimension(:), allocatable :: uessv

    real(kind=8), dimension(:), allocatable :: cdrsv
    real(kind=8), dimension(:), allocatable :: cessv

    ! Local variable declarations.
    integer :: ipc
    integer :: j2
    integer :: j3
    integer :: n
    integer :: nc
    integer :: ndrsn
    integer :: nerr
    integer :: nessn
    integer :: nmax
    integer :: nn
    integer :: np
    integer :: ns
    integer :: nsc
    integer :: ntpr
    integer :: nwatra

    integer :: ilnobl

    logical :: qx

    character(len=80) :: ux80
    character(len=56) :: ux56
    character(len=24) :: uphasv
    character(len=24) :: uwater
    character(len=24) :: uaqsln
    character(len=24) :: uptsld
    character(len=24) :: uptliq
    character(len=24) :: uptgas
    character(len=24) :: usblkf
    character(len=24) :: ux24
    character(len=8) :: uendit
    character(len=8) :: usedh
    character(len=8) :: upitz
    character(len=8) :: ustr
    character(len=8) :: ustr2
    character(len=8) :: ustr3
    character(len=8) :: uterm

    data uwater /'H2O                     '/
    data uaqsln /'Aqueous solution        '/
    data uptsld /'Solid                   '/
    data uptliq /'Liquid                  '/
    data uptgas /'Gas                     '/

    data uendit /'endit.  '/,uterm  /'+-------'/,usedh  /'SEDH    '/,upitz  /'Pitzer  '/

    ! Allocate scratch arrays.
    ALLOCATE(cessv(ncta_asv))
    ALLOCATE(uessv(ncta_asv))

    ALLOCATE(cdrsv(nbta1_asv))
    ALLOCATE(udrsv(nbta1_asv))

    ! Zero some arrays.
    nmax = nbta_asv
    call initiz(nbaspa,nmax)

    nmax = ndrsa_asv
    call initiz(ndrsa,nmax)
    call initaz(cdrsa,nmax)

    nmax = nsta_asv
    call initaz(vosp0a,nmax)

    nmax = 2*nsta_asv
    call initiz(ndrsra,nmax)

    nmax = narx_asv*ntpr_asv*nsta_asv
    call initaz(axhfsa,nmax)
    call initaz(axlksa,nmax)
    call initaz(axvfsa,nmax)

    ! SEDH arrays.
    nmax = nata_asv
    call initaz(azeroa,nmax)
    call initiz(insgfa,nmax)

    ! Pitzer arrays.
    nmax = ipbt_asv*napa_asv
    call initaz(palpaa,nmax)
    nmax = jpfc_asv*(ipbt_asv + 1)*nslta_asv
    call initaz(aslma,nmax)
    nmax = jpfc_asv*nmuta_asv
    call initaz(amua,nmax)

    ! Don't rewind the DATA1 file (nad1). It should be positioned
    ! properly by the preceding call to subroutine indath.f.
    write (noutpt,1000)
    write (nttyo,1000)
1000 format(/' Reading the rest of the DATA1 file ...')

    ! Initialize the error counter.
    nerr = 0

    ! Set the number of elements on the data file. This equal to
    ! the dimensioned limit. Note that nbta, the number of basis
    ! species, is determined in this subroutine as it is not usually
    ! equal to the corresponding dimensioned limit nbta_asv.
    ! In fact, nbta as returned by the subroutine may be subsequently
    ! increased due to "promotions" implied by the contents of the
    ! input file.
    ncta = ncta_asv

    ! Read the data file title.
    ntitld = 0

    do n = 1,ntid_asv
        read (nad1) utitld(n)
        ntitld = ntitld + 1

        if (utitld(n)(1:8) .eq. uterm(1:8)) then
            go to 140
        end if
    end do

140 continue

    ! Write the first 5 lines of the data file title to the screen
    ! file.
    write (nttyo,1060)
1060 format(/3x,'The data file title is (first 5 lines maximum):',/)

    nn = min(ntitld,5)

    do n = 1,nn
        if (utitld(n)(1:8) .eq. uterm(1:8)) then
            go to 150
        end if

        j2 = ilnobl(utitld(n))
        j2 = min(j2,74)
        write (nttyo,1070) utitld(n)(1:j2)
1070 format(5x,a)
    end do

150 continue

    ! Write the complete data file title, less terminator, if any,
    ! to the output file.
    write (noutpt,1080)
1080 format(/3x,'The data file title is:',/)

    do n = 1,ntitld
        if (utitld(n)(1:8) .eq. uterm(1:8)) then
            go to 160
        end if

        j2 = ilnobl(utitld(n))
        j2 = min(j2,74)
        write (noutpt,1070) utitld(n)(1:j2)
    end do

160 continue

    ! Search the data file title for the usual identifying three-letter
    ! keystring (e.g., com, hmw, etc.).
    udatfi = 'current'

    do n = 1,ntitld
        ux80 = utitld(n)
        call lejust(ux80)

        if (ux80(1:6) .eq. 'data0.') then
            udatfi = ux80(7:9)
            go to 170
        end if
    end do

170 continue

    write (noutpt,1090)
    write (nttyo,1090)
1090 format(/3x,'Continuing to read the DATA1 file ...')

    ! Read the actual number of temperature ranges in the temperature
    ! grid.
    read (nad1) ux56
    read (nad1) ntprt

    ! Read the actual number of coefficients in each range.
    do ntpr = 1,ntprt
        read (nad1) ux56
        read (nad1) narxt(ntpr)
    end do

    ! Read a flag indicating the degree of representation of enthalpy
    ! functions.
    read (nad1) ux56
    read (nad1) ipch

    ! Read a flag indicating the degree of representation of volume
    ! functions.
    read (nad1) ux56
    read (nad1) ipcv

    if (udakey(1:8) .eq. upitz(1:8)) then
        ! Read the Pitzer data block organization flag.
        read (nad1) ux56
        read (nad1) jpdblo

        ! Read the Pitzer parameter temperature function flag.
        read (nad1) ux56
        read (nad1) jptffl
    end if

    ! Read chemical elements data.
    read (nad1) ux24

    do nc = 1,ncta
        read (nad1) uelema(nc),atwta(nc)
    end do

    ! Read nominal temperature limits (Celsius) for the data file.
    read (nad1) ux56
    read (nad1) tdamin,tdamax

    ! Read the maximum temperatures (tempcu) of the temperature ranges.
    ! These are used to find the temperature range flag (ntpr) for any
    ! specified temperature (see EQLIB/gntpr.f).
    read (nad1) ux56

    do n = 1,ntprt
        read (nad1) tempcu(n)
    end do

    ! Read interpolating polynomial coefficients for computing the
    ! standard grid pressure (bars) as a function of temperature
    ! (Celsius). Often this grid corresponds to 1.013 bar up to 100C
    ! and the steam-liquid water saturation pressure at higher
    ! temperatures.
    ! Calling sequence substitutions:
    !   apresg for arr
    call indatc(apresg,nad1,narx_asv,narxt,ntpr_asv,ntprt,ux24)

    if (ipcv .ge. 0) then
        ! Read interpolating polynomial coefficients for computing the
        ! recommended pressure envelope half width (bars).
        ! Calling sequence substitutions:
        !   aprehw for arr
        call indatc(aprehw,nad1,narx_asv,narxt,ntpr_asv,ntprt,ux24)
    end if

    if (udakey(1:8) .eq. usedh(1:8)) then
        ! The aqeuous species activity coefficient formalism is
        ! consistent with the B-dot equation, the Davies equation,
        ! or similar equations.
        ! Read interpolating polynomial coefficients for the Debye-Huckel
        ! A(gamma,10) parameter and related parameters.
        ! Calling sequence substitutions:
        !   aadh for arr
        call indatc(aadh,nad1,narx_asv,narxt,ntpr_asv,ntprt,ux24)

        if (ipch .ge. 0) then
            ! Calling sequence substitutions:
            !   aadhh for arr
            call indatc(aadhh,nad1,narx_asv,narxt,ntpr_asv,ntprt,ux24)

            do ipc = 1,ipch
                ! Calling sequence substitutions:
                !   adadhh for arr
                !   ipch_asv for ipcx_asv
                call indatd(adadhh,ipc,ipch_asv,nad1,narxt,narx_asv,ntprt,ntpr_asv,ux24)
            end do
        end if

        if (ipcv .ge. 0) then
            ! Calling sequence substitutions:
            !   aadhv for arr
            call indatc(aadhv,nad1,narx_asv,narxt,ntpr_asv,ntprt,ux24)

            do ipc = 1,ipcv
                ! Calling sequence substitutions:
                !   adadhv for arr
                !   ipcv_asv for ipcx_asv
                call indatd(adadhv,ipc,ipcv_asv,nad1,narxt,narx_asv,ntprt,ntpr_asv,ux24)
            end do
        end if

        ! Read interpolating polynomial coefficients for the Debye-Huckel
        ! B(gamma) and related parameters.
        ! Calling sequence substitutions:
        !   abdh for arr
        call indatc(abdh,nad1,narx_asv,narxt,ntpr_asv,ntprt,ux24)

        if (ipch .ge. 0) then
            ! Calling sequence substitutions:
            !   abdhh for arr
            call indatc(abdhh,nad1,narx_asv,narxt,ntpr_asv,ntprt,ux24)

            do ipc = 1,ipch
                ! Calling sequence substitutions:
                !   adbdhh for arr
                !   ipch_asv for ipcx_asv
                call indatd(adbdhh,ipc,ipch_asv,nad1,narxt,narx_asv,ntprt,ntpr_asv,ux24)
            end do
        end if

        if (ipcv .ge. 0) then
            ! Calling sequence substitutions:
            !   abdhv for arr
            call indatc(abdhv,nad1,narx_asv,narxt,ntpr_asv,ntprt,ux24)

            do ipc = 1,ipcv
                ! Calling sequence substitutions:
                !   adbdhv for arr
                !   ipcv_asv for ipcx_asv
                call indatd(adbdhv,ipc,ipcv_asv,nad1,narxt,narx_asv,ntprt,ntpr_asv,ux24)
            end do
        end if

        ! Read interpolating polynomial coefficients for the B-dot
        ! parameter and related parameters.
        ! Calling sequence substitutions:
        !   abdot for arr
        call indatc(abdot,nad1,narx_asv,narxt,ntpr_asv,ntprt,ux24)

        if (ipch .ge. 0) then
            ! Calling sequence substitutions:
            !   abdoth for arr
            call indatc(abdoth,nad1,narx_asv,narxt,ntpr_asv,ntprt,ux24)

            do ipc = 1,ipch
                ! Calling sequence substitutions:
                !   adbdth for arr
                !   ipch_asv for ipcx_asv
                call indatd(adbdth,ipc,ipch_asv,nad1,narxt,narx_asv,ntprt,ntpr_asv,ux24)
            end do
        end if

        if (ipcv .ge. 0) then
            ! Calling sequence substitutions:
            !   abdotv for arr
            call indatc(abdotv,nad1,narx_asv,narxt,ntpr_asv,ntprt,ux24)

            do ipc = 1,ipcv
                ! Calling sequence substitutions:
                !   adbdtv for arr
                !   ipcv_asv for ipcx_asv
                call indatd(adbdtv,ipc,ipcv_asv,nad1,narxt,narx_asv,ntprt,ntpr_asv,ux24)
            end do
        end if

        ! Read coefficients for the Drummond (1981) equation.
        ! This equation represents the activity coeficient
        ! of CO2(aq) as a function of temperature and ionic
        ! strength in sodium chloride solutions.
        read (nad1) ux24
        read (nad1) (cco2(n), n = 1,5)
    end if

    if (udakey(1:8) .eq. upitz(1:8)) then
        ! The aqeuous species activity coefficient formalism is
        ! consistent with Pitzer's equations. Read interpolating
        ! polynomial coefficients for Debye-Huckel A(phi) and
        ! related parameters.
        ! Calling sequence substitutions:
        !   aaphi for arr
        call indatc(aaphi,nad1,narx_asv,narxt,ntpr_asv,ntprt,ux24)

        if (ipch .ge. 0) then
            ! Calling sequence substitutions:
            !   aadhh for arr
            call indatc(aadhh,nad1,narx_asv,narxt,ntpr_asv,ntprt,ux24)

            do ipc = 1,ipch
                ! Calling sequence substitutions:
                !   adadhh for arr
                !   ipch_asv for ipcx_asv
                call indatd(adadhh,ipc,ipch_asv,nad1,narxt,narx_asv,ntprt,ntpr_asv,ux24)
            end do
        end if

        if (ipcv .ge. 0) then
            ! Calling sequence substitutions:
            !   aadhv for arr
            call indatc(aadhv,nad1,narx_asv,narxt,ntpr_asv,ntprt,ux24)

            do ipc = 1,ipcv
                ! Calling sequence substitutions:
                !   adadhv for arr
                !   ipcv_asv for ipcx_asv
                call indatd(adadhv,ipc,ipcv_asv,nad1,narxt,narx_asv,ntprt,ntpr_asv,ux24)
            end do
        end if
    end if

    ! Read interpolating polynomial coefficients for the log K of
    ! the "Eh" reaction:
    !   2 H2O(l) = 2 O2(g) + 4 H+ + 4 e-
    ! Calling sequence substitutions:
    !   axlke for arr
    call indatc(axlke,nad1,narx_asv,narxt,ntpr_asv,ntprt,ux24)

    if (ipch .ge. 0) then
        ! Calling sequence substitutions:
        !   axhfe for arr
        call indatc(axhfe,nad1,narx_asv,narxt,ntpr_asv,ntprt,ux24)

        do ipc = 1,ipch
            ! Calling sequence substitutions:
            !   adhfe for arr
            !   ipch_asv for ipcx_asv
            call indatd(adhfe,ipc,ipch_asv,nad1,narxt,narx_asv,ntprt,ntpr_asv,ux24)
        end do
    end if

    if (ipcv .ge. 0) then
        ! Calling sequence substitutions:
        !   axvfe for arr
        call indatc(axvfe,nad1,narx_asv,narxt,ntpr_asv,ntprt,ux24)

        do ipc = 1,ipcv
            ! Calling sequence substitutions:
            !   advfe for arr
            !   ipcv_asv for ipcx_asv
            call indatd(advfe,ipc,ipcv_asv,nad1,narxt,narx_asv,ntprt,ntpr_asv,ux24)
        end do
    end if

    ! Read the phases and species and assocated data from the
    ! data file. The following variables are used as counters:
    !   np     = phase index
    !   ns     = species index
    ! Initialize the qclnsa array to .false.
    qx = .false.

    call initlv(qclnsa,nsta_asv,qx)

    ! Initialize the following variables, all of which are used
    ! below as counters.
    !   nessn  = counter for the cessa and nessa arrays
    !   ndrsn  = counter for the cdrsa and ndrsa arrays
    !   ns     = species index
    !   nbta   = number of basis species
    !   nata   = number of aqueous species
    !   ngta   = number of gas species
    !   nmta   = number of pure mineral species
    !   nlta   = number of pure liquid species
    !   nxta   = number of solid solution phases
    nessn = 0
    ndrsn = 0
    ns = 0
    nbta = 0
    nata = 0
    ngta = 0
    nmta = 0
    nlta = 0
    nxta = 0

    ! Process the aqueous species superblock. Set up the aqueous
    ! solution as the first phase.
    np = 1
    uphasv(1:24) = uaqsln(1:24)
    iaqsla = np
    uphasa(np)(1:24) = uphasv(1:24)
    uptypa(np)(1:24) = uaqsln(1:24)
    usblkf(1:24) = uaqsln(1:24)
    ncmpra(1,1) = 1
    narn1a = 1

    ! Read the superblock header.
    read (nad1) ustr,ustr2,ustr3

    ! Read the blocks in the current superblock.
    call indats(adhfsa,advfsa,axhfsa,axlksa,axvfsa,cdrsa,cdrsv,cessa,cessv,ipch,ipch_asv,ipcv,ipcv_asv,mwtspa,nad1,narxt,narx_asv,nata,nata_asv,nbta,nbta_asv,nbta1_asv,nbtafd,ncmpra,ncta,ncta_asv,ndrsa,ndrsa_asv,ndrsn,ndrsra,nerr,nessa,nessa_asv,nessn,nessra,ngta,ngta_asv,nlta,nlta_asv,nmta,nmta_asv,noutpt,np,npta_asv,ns,nsta_asv,ntprt,ntpr_asv,nttyo,uaqsln,ubasp,udrsv,uelema,uendit,uessv,uphasa,uphasv,uptgas,uptliq,uptsld,uptypa,usblkf,uspeca,vosp0a,zchara)

    ncmpra(2,np) = ns
    narn2a = ns

    ! Check to see that water is the first species in the phase
    ! "aqueous solution". This is required so that all solute species
    ! can be cleanly referenced in simple loop; that is, in a loop of
    ! the form "do ns = narn1a + 1,narn2a".
    nwatra = 0

    if (uspeca(narn1a)(1:24) .eq. uwater(1:24)) then
        nwatra = 1
    end if

    if (nwatra .eq. 0) then
        j3 = ilnobl(uwater)
        write (noutpt,1150) uwater(1:j3),uaqsln(1:j2)
        write (nttyo,1150) uwater(1:j3),uaqsln(1:j2)
1150 format(/' * Error - (EQLIB/indata) The species ',/7x,a,' (',a,') must be listed on the',/7x,'data file as the first species in its phase.')

        stop
    end if

    ! Process the pure minerals superblock. Each pure mineral is both
    ! a species and a phase. The phase name variable uphasv will be set
    ! to the name of each mineral as its block is read.
    usblkf(1:24) = uptsld(1:24)
    nmrn1a = ns + 1

    ! Read the superblock header.
    read (nad1) ustr,ustr2,ustr3

    ! Read the blocks in the current superblock.
    call indats(adhfsa,advfsa,axhfsa,axlksa,axvfsa,cdrsa,cdrsv,cessa,cessv,ipch,ipch_asv,ipcv,ipcv_asv,mwtspa,nad1,narxt,narx_asv,nata,nata_asv,nbta,nbta_asv,nbta1_asv,nbtafd,ncmpra,ncta,ncta_asv,ndrsa,ndrsa_asv,ndrsn,ndrsra,nerr,nessa,nessa_asv,nessn,nessra,ngta,ngta_asv,nlta,nlta_asv,nmta,nmta_asv,noutpt,np,npta_asv,ns,nsta_asv,ntprt,ntpr_asv,nttyo,uaqsln,ubasp,udrsv,uelema,uendit,uessv,uphasa,uphasv,uptgas,uptliq,uptsld,uptypa,usblkf,uspeca,vosp0a,zchara)

    nmrn2a = ns

    if (nmrn2a .le. narn2a) then
        nmrn1a = 0
    end if

    ! There is presently no superblock corresponding to pure liquids.
    ! Set up water as a pure liquid.
    usblkf(1:24) = uptliq(1:24)
    nlrn1a = 0
    nlrn2a = 0
    nlta = 0

    if (nwatra .gt. 0) then
        np = np + 1
        ns = ns + 1
        uphasv(1:24) = uspeca(nwatra)(1:24)
        uphasa(np)(1:24) = uphasv(1:24)
        uptypa(np)(1:24) = uptliq(1:24)
        ncmpra(1,np) = ns
        ncmpra(2,np) = ns
        nsc = 1
        call clones(axlksa,cdrsa,cessa,mwtspa,narx_asv,ndrsa,ndrsa_asv,ndrsn,ndrsra,nessa,nessa_asv,nessn,nessra,np,npta_asv,ns,nsc,nsta_asv,ntpr_asv,uphasa,uspeca,zchara)
        qclnsa(ns) = .true.
        nlrn1a = ns
        nlrn2a = ns
        nlta = 1
    end if

    ! Process the gas species superblock. Set up the gas phase
    ! as the next phase.
    uphasv(1:24) = uptgas(1:24)
    usblkf(1:24) = uptgas(1:24)
    np = np + 1
    igas = np
    uphasa(np)(1:24) = uphasv(1:24)
    uptypa(np)(1:24) = uptgas(1:24)
    ncmpra(1,np) = ns + 1
    ngrn1a = ns + 1

    ! Read the superblock header.
    read (nad1) ustr,ustr2,ustr3

    ! Read the blocks in the current superblock.
    call indats(adhfsa,advfsa,axhfsa,axlksa,axvfsa,cdrsa,cdrsv,cessa,cessv,ipch,ipch_asv,ipcv,ipcv_asv,mwtspa,nad1,narxt,narx_asv,nata,nata_asv,nbta,nbta_asv,nbta1_asv,nbtafd,ncmpra,ncta,ncta_asv,ndrsa,ndrsa_asv,ndrsn,ndrsra,nerr,nessa,nessa_asv,nessn,nessra,ngta,ngta_asv,nlta,nlta_asv,nmta,nmta_asv,noutpt,np,npta_asv,ns,nsta_asv,ntprt,ntpr_asv,nttyo,uaqsln,ubasp,udrsv,uelema,uendit,uessv,uphasa,uphasv,uptgas,uptliq,uptsld,uptypa,usblkf,uspeca,vosp0a,zchara)

    ngrn2a = ns

    if (ngrn2a .le. nlrn2a) then
        ngrn1a = 0
    end if

    ncmpra(2,np) = ns

    ! Process the solid solution superblock.
    read (nad1) ustr,ustr2,ustr3
    nxrn1a = ns + 1
    ixrn1a = np + 1

    call indatp(apxa,axlksa,bpxa,cdrsa,cessa,iapxa_asv,iapxta,ibpxa_asv,ibpxta,ikta_asv,jsola,mwtspa,nad1,narx_asv,ncmpra,ndrsa,ndrsa_asv,ndrsn,ndrsra,nerr,nessa,nessa_asv,nessn,nessra,nmrn1a,nmrn2a,noutpt,np,npta_asv,ns,nsta_asv,ntpr_asv,nttyo,nxta,nxta_asv,qclnsa,uendit,uspeca,uphasa,uptsld,uptypa,zchara)

    npta = np
    nsta = ns
    nxrn2a = ns
    ixrn2a = np

    if (nxrn2a .le. nlrn2a) then
        nxrn1a = 0
        nxrn2a = 0
        ixrn1a = 0
        ixrn2a = 0
    end if

    ! Set up the nbaspa array and expand the elements of the ubasp array
    ! to the full 48 characters. Convert basis species indices in the
    ! ndrsa array to species indices.
    call dfbasp(nbaspa,nbta,nbta_asv,ndrsa,ndrsa_asv,ndrsra,nerr,noutpt,nsta,nsta_asv,nttyo,ubasp,uspeca)

    if (udakey(1:8) .eq. usedh(1:8)) then
        ! The aqeuous species activity coefficient formalism is
        ! consistent with the B-dot equation, the Davies equation,
        ! or similar equations. Read the appropriate data for
        ! this formalism.
        call inbdot(azeroa,insgfa,nad1,narn1a,narn2a,nata,nata_asv,nerr,noutpt,nsta_asv,nttyo,uspeca)
    end if

    if (udakey(1:8) .eq. upitz(1:8)) then
        ! The aqeuous species activity coefficient formalism is
        ! consistent with Pitzer's equations. Read the appropriate
        ! data for this formalism.
        call inupt(amua,aslma,ielam,ipbt_asv,jpdblo,jpfc_asv,nad1,nalpaa,napa_asv,napta,narn1a,narn2a,nerr,nmuta,nmuta_asv,nmuxa,noutpt,nslta,nslta_asv,nslxa,nsta_asv,nttyo,palpaa,uspeca,zchara)
    end if

    ! The following is a bit of nonsense so compiler warnings will
    ! not be generated saying that ux80, ustr, ustr2, and ustr3 are
    ! not used (they are all used to read unused data).
    ux80(1:56) = ux56
    ustr3 = ux80(1:8)
    ustr = ustr2
    ustr2 = ustr3
    ustr3 = ustr

    ! Deallocate the scratch arrays used by calls to EQLIB/indats.f.
    DEALLOCATE(cessv)
    DEALLOCATE(uessv)

    DEALLOCATE(cdrsv)
    DEALLOCATE(udrsv)

    if (nerr .gt. 0) then
        stop
    end if

    write (noutpt,1200)
    write (nttyo,1200)
1200 format(/'   Done reading the DATA1 file.')
end subroutine indata