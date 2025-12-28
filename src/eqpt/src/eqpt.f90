program eqpt
    !! This is the main program of the EQPT code. Configuration
    !! identification, the copyright statement, legal disclaimers,
    !! and similar statements are contained in EQPT/aaaeqt.f, the
    !! lead-off subroutine in the EQPT source code. A short description
    !! of this program is also contained in that subroutine.
    implicit none

    ! File path parameters
    integer :: numargs
    character(len=1024) :: temppath
    character(len=:), allocatable :: data0path
    integer :: pathindices(2)
    character(len=:), allocatable :: basename
    character(len=:), allocatable :: ofile
    character(len=:), allocatable :: d1file
    character(len=:), allocatable :: d1ffile
    character(len=:), allocatable :: sfile
    character(len=:), allocatable :: d0sfile

    ! Dimensioning parameters:
    !   iapxpa = maximum number of interaction coefficients for
    !              computing activity coefficients in solid solutions
    !   ibpxpa = maximum number of site-mixing parameters for
    !              computing activity coefficients in solid solutions
    integer :: iapxpa
    integer :: ibpxpa
    integer :: iktpar

    parameter(iapxpa = 20,ibpxpa = 20)

    ! Array allocation size variables used in EQPT.
    !   ipch_asv = The order for pressure corrections to enthalpy
    !                functions
    !   ipcv_asv = The order for pressure corrections to volume
    !                functions; the maximum order for pressure
    !                corrections to log K and other Gibbs-energy-based
    !                functions is one greater than this
    !   nap_asv  = The maximum number of distinct sets of Pitzer alpha
    !                parameters
    !   narx_asv = The maximum number of coefficients per temperature
    !                range
    !   nat_asv  = The number of aqueous species on the data file
    !   nazt_asv = the number of aqeuous species on the data file for
    !                which hard core diameters are specified
    !   nbt_asv  = The number of basis species on the data file
    !   nct_asv  = The number of chemical elements on the data file
    !   ngt_asv  = The number of gas species on the data file
    !   nlt_asv  = The number of pure liquids on the data file
    !   nmt_asv  = The number of pure minerals on the data file
    !   npx2_asv = The number of pairs of species not of the same
    !                charge sign for which Pitzer parameters are
    !                defined; typically, one of the pair is a cation
    !                and the other is an anion, but one or both
    !                species may also be electrically neutral
    !   npx3_asv = The number of triplets of species corresponding to
    !                aqueous electrolyte mixtures for which Pitzer
    !                parameters are defined; generally, no more than
    !                two of these may have an electrical charge number
    !                that is postive, negative, or zero
    !   ntid_asv = The number of lines in the data file title
    !   ntpr_asv = The number of temperature ranges
    !   ipbt_asv = The greatest index for a Pitzer beta(n) coefficient
    !                 (the first such index is zero); also the number
    !                 number of Pitzer alpha coefficients in a set
    !   jpfc_asv = The number of coefficients in a function for
    !                representing a Pitzer coefficient
    !   ndb_asv  = maximum number of distinct points on the temperature
    !                grid; ndb_asv = ntpr_asv*(narx_asv - 1) + 1
    integer :: ipbt_asv
    integer :: ipch_asv
    integer :: ipcv_asv
    integer :: jpfc_asv
    integer :: narx_asv
    integer :: nat_asv
    integer :: nazt_asv
    integer :: nbt_asv
    integer :: nct_asv
    integer :: ngt_asv
    integer :: nlt_asv
    integer :: nmt_asv
    integer :: npx2_asv
    integer :: npx3_asv
    integer :: ntid_asv
    integer :: ntpr_asv

    integer :: ndb_asv

    ! Array allocation size variables which are determined in EQPT
    ! but only written on the DATA1 file for use by EQ3NR and EQ6.
    !   ikt_asv  = The maximum number of end-member component species
    !                of any solid solution on the data filen
    !   nlat_asv = The number of members in the set of Pitzer lambda
    !                coefficients
    !   nmut_asv = The number of members in the set of Pitzer mu
    !                coefficients
    !   npt_asv  = The number of phases of all types on the data file
    !   nst_asv  = The number of species of all types on the data file
    !   nxt_asv  = The number of solid-solution phases on the data file
    integer :: ikt_asv
    integer :: nap_asv
    integer :: nlat_asv
    integer :: nmut_asv
    integer :: npt_asv
    integer :: nst_asv
    integer :: nxt_asv

    ! Global variable declarations.
    ! Dimensioning variables.
    integer :: iapxmx
    integer :: ibpxmx
    integer :: iktmax
    integer :: ipbtmx
    integer :: ipchmx
    integer :: ipcvmx
    integer :: jpfcmx
    integer :: narxmx
    integer :: natmax
    integer :: naztmx
    integer :: nbtmax
    integer :: nctmax
    integer :: ngtmax
    integer :: nltmax
    integer :: nmtmax
    integer :: npx2mx
    integer :: npx3mx
    integer :: ntidmx
    integer :: ntprmx
    integer :: nxtmax
    integer :: nbtmx1
    integer :: nbtmx2
    integer :: ndbmax

    ! File unit numbers.
    integer :: ndata0
    integer :: ndata1
    integer :: ndat0s
    integer :: ndat1f
    integer :: noutpt
    integer :: nslist
    integer :: nttyo

    ! Other variables.
    integer, dimension(:), allocatable :: insgf
    integer, dimension(:), allocatable :: ipivot
    integer, dimension(:), allocatable :: issot
    integer, dimension(:), allocatable :: nacdpr
    integer, dimension(:), allocatable :: narxt
    integer, dimension(:), allocatable :: nentei
    integer, dimension(:), allocatable :: nentri

    integer :: ier
    integer :: ipch
    integer :: ipcv
    integer :: irang
    integer :: nazt
    integer :: nat
    integer :: nbt
    integer :: nct
    integer :: ncvaz
    integer :: ndbptg
    integer :: ndbptl
    integer :: nerr
    integer :: ngt
    integer :: nlt
    integer :: nmt
    integer :: nmodwr
    integer :: nsb
    integer :: nthdt
    integer :: ntitld
    integer :: ntprt
    integer :: nwarn
    integer :: nxt

    logical, dimension(:), allocatable :: qpdaz

    logical :: qelect

    character(len=80), dimension(:), allocatable :: utitld
    character(len=24), dimension(:), allocatable :: uaqsp
    character(len=24), dimension(:), allocatable :: uazp
    character(len=24), dimension(:), allocatable :: udrsi
    character(len=24), dimension(:), allocatable :: ugassp
    character(len=24), dimension(:), allocatable :: uliqsp
    character(len=24), dimension(:), allocatable :: uminsp
    character(len=24), dimension(:), allocatable :: uspec
    character(len=24), dimension(:), allocatable :: ussoph
    character(len=24), dimension(:,:), allocatable :: ussosp
    character(len=16), dimension(:), allocatable :: udbval
    character(len=8), dimension(:), allocatable :: uelem
    character(len=8), dimension(:), allocatable :: uessi

    character(len=16) :: udbfmt
    character(len=8) :: uplatc
    character(len=8) :: uplatm
    character(len=8) :: uveeqt
    character(len=8) :: usteqt
    character(len=8) :: uveelu
    character(len=8) :: ustelu
    character(len=8) :: uakey
    character(len=8) :: uethfl

    real(kind=8), dimension(:), allocatable :: atwt
    real(kind=8), dimension(:), allocatable :: azero
    real(kind=8), dimension(:), allocatable :: cdrsi
    real(kind=8), dimension(:), allocatable :: cessi
    real(kind=8), dimension(:), allocatable :: mtotr
    real(kind=8), dimension(:), allocatable :: zaqsp
    real(kind=8), dimension(:), allocatable :: zchar
    real(kind=8), dimension(:,:), allocatable :: cdrs
    real(kind=8), dimension(:,:), allocatable :: cess

    real(kind=8) :: cco2(5)

    real(kind=8), dimension(:,:,:), allocatable :: xhfs
    real(kind=8), dimension(:,:,:), allocatable :: xlks
    real(kind=8), dimension(:,:,:), allocatable :: xvfs
    real(kind=8), dimension(:,:,:,:), allocatable :: dhfs
    real(kind=8), dimension(:,:,:,:), allocatable :: dvfs

    real(kind=8), dimension(:,:), allocatable :: adh
    real(kind=8), dimension(:,:), allocatable :: adhh
    real(kind=8), dimension(:,:), allocatable :: adhv
    real(kind=8), dimension(:,:), allocatable :: aphi
    real(kind=8), dimension(:,:), allocatable :: bdh
    real(kind=8), dimension(:,:), allocatable :: bdhh
    real(kind=8), dimension(:,:), allocatable :: bdhv
    real(kind=8), dimension(:,:), allocatable :: bdot
    real(kind=8), dimension(:,:), allocatable :: bdoth
    real(kind=8), dimension(:,:), allocatable :: bdotv
    real(kind=8), dimension(:,:), allocatable :: prehw
    real(kind=8), dimension(:,:), allocatable :: presg
    real(kind=8), dimension(:,:), allocatable :: xhfe
    real(kind=8), dimension(:,:), allocatable :: xlke
    real(kind=8), dimension(:,:), allocatable :: xvfe
    real(kind=8), dimension(:,:,:), allocatable :: dadhh
    real(kind=8), dimension(:,:,:), allocatable :: dadhv
    real(kind=8), dimension(:,:,:), allocatable :: dbdhh
    real(kind=8), dimension(:,:,:), allocatable :: dbdhv
    real(kind=8), dimension(:,:,:), allocatable :: dbdth
    real(kind=8), dimension(:,:,:), allocatable :: dbdtv
    real(kind=8), dimension(:,:,:), allocatable :: dhfe
    real(kind=8), dimension(:,:,:), allocatable :: dvfe

    real(kind=8), dimension(:), allocatable :: tmpcmx
    real(kind=8), dimension(:), allocatable :: xdbval
    real(kind=8), dimension(:,:), allocatable :: apr
    real(kind=8), dimension(:,:), allocatable :: avgrid
    real(kind=8), dimension(:,:), allocatable :: tempc
    real(kind=8), dimension(:,:), allocatable :: tempcs

    real(kind=8), dimension(:), allocatable :: tvec
    real(kind=8), dimension(:), allocatable :: tvecs
    real(kind=8), dimension(:), allocatable :: cof
    real(kind=8), dimension(:), allocatable :: xvec
    real(kind=8), dimension(:), allocatable :: yvec
    real(kind=8), dimension(:,:), allocatable :: aamatr
    real(kind=8), dimension(:,:), allocatable :: gmmatr

    real(kind=8) :: pcvaz

    real(kind=8) :: eps100
    real(kind=8) :: smp100
    real(kind=8) :: tdamax
    real(kind=8) :: tdamin
    real(kind=8) :: tvecmx

    ! Pitzer interaction parameter arrays.
    integer :: naapr
    integer :: ncapr
    integer :: nccpr
    integer :: nnapr
    integer :: nncpr
    integer :: nnnpr
    integer :: nn2pr
    integer :: naactr
    integer :: na2ctr
    integer :: nccatr
    integer :: nc2atr
    integer :: nncatr
    integer :: nn2ntr
    integer :: nn3tr

    integer, dimension(:), allocatable :: in2pr
    integer, dimension(:), allocatable :: in3tr
    integer, dimension(:,:), allocatable :: iaapr
    integer, dimension(:,:), allocatable :: icapr
    integer, dimension(:,:), allocatable :: iccpr
    integer, dimension(:,:), allocatable :: inapr
    integer, dimension(:,:), allocatable :: incpr
    integer, dimension(:,:), allocatable :: innpr
    integer, dimension(:,:), allocatable :: iaactr
    integer, dimension(:,:), allocatable :: ia2ctr
    integer, dimension(:,:), allocatable :: iccatr
    integer, dimension(:,:), allocatable :: ic2atr
    integer, dimension(:,:), allocatable :: incatr
    integer, dimension(:,:), allocatable :: in2ntr

    integer :: jassan
    integer :: jassca
    integer :: jassne
    integer :: npx2t
    integer :: npx3t

    integer :: ncvaa
    integer :: ncvca
    integer :: ncvcc
    integer :: ncvna
    integer :: ncvnc
    integer :: ncvnn
    integer :: ncvn2
    integer :: ncvaac
    integer :: ncvcca
    integer :: ncvn2n
    integer :: ncvnca

    integer :: npxca
    integer :: npxth
    integer :: npxni
    integer :: npxnn
    integer :: npxn2
    integer :: npxpsi
    integer :: npxzet
    integer :: npxn2n

    logical, dimension(:), allocatable :: qpdaa
    logical, dimension(:), allocatable :: qpdca
    logical, dimension(:), allocatable :: qpdcc
    logical, dimension(:), allocatable :: qpdna
    logical, dimension(:), allocatable :: qpdnc
    logical, dimension(:), allocatable :: qpdnn
    logical, dimension(:), allocatable :: qpdn2

    logical, dimension(:), allocatable :: qpdaac
    logical, dimension(:), allocatable :: qpdcca
    logical, dimension(:), allocatable :: qpdnca
    logical, dimension(:), allocatable :: qpdn2n

    character(len=24), dimension(:,:), allocatable :: upair
    character(len=24), dimension(:,:), allocatable :: uthdtr
    character(len=24), dimension(:,:), allocatable :: utripl

    real(kind=8), dimension(:,:), allocatable :: alpha
    real(kind=8), dimension(:,:), allocatable :: alphca
    real(kind=8), dimension(:,:), allocatable :: acphi
    real(kind=8), dimension(:,:), allocatable :: amua2c
    real(kind=8), dimension(:,:), allocatable :: amuaac
    real(kind=8), dimension(:,:), allocatable :: amucca
    real(kind=8), dimension(:,:), allocatable :: amuc2a
    real(kind=8), dimension(:,:), allocatable :: amunca
    real(kind=8), dimension(:,:), allocatable :: amun2n
    real(kind=8), dimension(:,:), allocatable :: amun3
    real(kind=8), dimension(:,:,:), allocatable :: abeta
    real(kind=8), dimension(:,:,:), allocatable :: alamaa
    real(kind=8), dimension(:,:,:), allocatable :: alamca
    real(kind=8), dimension(:,:,:), allocatable :: alamcc
    real(kind=8), dimension(:,:,:), allocatable :: alamna
    real(kind=8), dimension(:,:,:), allocatable :: alamnc
    real(kind=8), dimension(:,:,:), allocatable :: alamnn
    real(kind=8), dimension(:,:,:), allocatable :: alamn2

    real(kind=8), dimension(:,:), allocatable :: atheta
    real(kind=8), dimension(:,:), allocatable :: apsi

    real(kind=8) :: pcvaa
    real(kind=8) :: pcvca
    real(kind=8) :: pcvcc
    real(kind=8) :: pcvna
    real(kind=8) :: pcvnc
    real(kind=8) :: pcvnn
    real(kind=8) :: pcvn2

    real(kind=8) :: pcvaac
    real(kind=8) :: pcvcca
    real(kind=8) :: pcvn2n
    real(kind=8) :: pcvnca

    ! XX
    real(kind=8) :: apx(iapxpa)
    real(kind=8) :: bpx(ibpxpa)

    ! Local variable declarations.
    integer :: i
    integer :: iexec0
    integer :: itgenf
    integer :: j
    integer :: jexec0
    integer :: jpdblo
    integer :: jptffl
    integer :: j2
    integer :: j3
    integer :: j4
    integer :: k
    integer :: n
    integer :: natm1
    integer :: nch
    integer :: nco
    integer :: nmax
    integer :: npx2
    integer :: npx2r1
    integer :: npx2r2
    integer :: nrecl
    integer :: nthd
    integer :: ntpr

    integer :: ilnobl

    logical :: q500fl

    character(len=16) :: uacfst
    character(len=11) :: utime0
    character(len=11) :: utime1
    character(len=8) :: ux8
    character(len=8) :: ux8a
    character(len=8) :: ux8b
    character(len=9) :: udate0
    character(len=9) :: udate1
    character(len=80) :: uline

    real(kind=8) :: tcpu
    real(kind=8) :: texec0
    real(kind=8) :: trun
    real(kind=8) :: tuser

    data nrecl /0/

    ! BEGIN_MACHINE_DEPENDENT_CODE
    !   On some systems, a BLOCK DATA subroutine must be declared in
    !   an EXTERNAL statement to assure proper loading. On some other
    !   systems, this is not necessary, but neither it is not harmful.
    !   On yet some other systems, the EXTERNAL statement below may
    !   cause a problem. If so, try commenting it out. If you still
    !   have trouble, consult your local system documentation or
    !   experiment to find out how to correctly handle a BLOCK DATA
    !   subroutine on your system. The EXTERNAL statement below should
    !   not cause a problem if you are using a compiler which is fully
    !   compliant with the Fortran 90 standard. However, there is
    !   no guarantee that it will be adequate to assure correct loading
    !   of the BLOCK DATA subroutine.
    external bkdeqp

    ! END_MACHINE_DEPENDENT_CODE
    ! BEGIN_MACHINE_DEPENDENT_CODE
    !   Define the console output device number. A value of 6 applies
    !   in most cases.
    data nttyo  /6/

    ! END_MACHINE_DEPENDENT_CODE
    ! Get time and date at the start of execution. This information
    ! will be used for time and date stamping.
    call initim(iexec0,jexec0,texec0,noutpt,nttyo,udate0,utime0)

    ! Disable underflow trapping, if any.
    call undflw

    ! Set file unit numbers to zero.
    noutpt = 0
    ndata1 = 0
    ndat1f = 0
    nslist = 0
    ndata0 = 0
    ndat0s = 0

    ! Open all output files.
    numargs = COMMAND_ARGUMENT_COUNT()

    if (numargs.eq.1) then
        call GET_COMMAND_ARGUMENT(1,temppath)
        data0path = TRIM(temppath)
    else
        write (0, *) 'usage: eqpt <data0>'
        stop
    end if

    call getbasename(data0path,pathindices)
    basename = data0path(pathindices(1):pathindices(2))

    call openin(noutpt,nttyo,data0path,'formatted',ndata0)

    ofile = basename // '.po'
    d1file = basename // '.d1'
    d1ffile = basename // '.d1f'
    sfile = basename // '.s'
    d0sfile = basename // '.d0s'

    call openou(noutpt,nttyo,ofile,'formatted',nrecl,noutpt)
    call openou(noutpt,nttyo,d1file,'unformatted',nrecl,ndata1)
    call openou(noutpt,nttyo,d1ffile,'formatted',nrecl,ndat1f)
    call openou(noutpt,nttyo,sfile,'formatted',nrecl,nslist)
    call openou(noutpt,nttyo,d0sfile,'formatted',nrecl,ndat0s)

    ! Make a copy of the DATA0 file, stripped of comments.
    call stripl(ndata0,ndat0s)
    close(ndata0)
    ndata0 = 0
    rewind ndat0s

    ! Get configuration identification data.
    call aaaeqt(usteqt,uveeqt)
    call aaaelu(ustelu,uveelu)
    call platfd(uplatc,uplatm)

    ! Write configuration identification data, the copyright statement,
    ! and any remaining statements or disclaimers.
    i = index(uveeqt,' ') - 1
    j = index(uveeqt,'.') - 1
    k = index(uplatc,' ') - 1

    if (i .le. 0) then
        i = 8
    end if

    if (j .le. 0) then
        j = 8
    end if

    if (k .le. 0) then
        k = 8
    end if

    write (nttyo,1000) uveeqt(1:i),uveeqt(1:j),uveeqt(1:i),uplatc(1:k)
    write (noutpt,1000) uveeqt(1:i),uveeqt(1:j),uveeqt(1:i),uplatc(1:k)
1000 format(/' EQ3/6, Version ',a,' (EQ3/6-V',a,'-REL-V',a,'-',a,')')

    i = j
    j = index(usteqt,' ') - 1
    k = index(uplatm,' ') - 1

    if (j .le. 0) then
        j = 8
    end if

    if (k .le. 0) then
        k = 8
    end if

    write (nttyo,1010) uveeqt(1:i),usteqt(1:j),uplatm(1:k)
    write (noutpt,1010) uveeqt(1:i),usteqt(1:j),uplatm(1:k)
1010 format(' EQPT Data File Preprocessor Code (EQ3/6-V',a,'-EQPT-EXE-',a,'-',a,')')

    write (noutpt,1020)
    write (nttyo,1020)
1020 format(' Supported by the following EQ3/6 libraries:')

    i = index(uveelu,'.') - 1
    j = index(ustelu,' ') - 1

    if (i .le. 0) then
        i = 8
    end if

    if (j .le. 0) then
        j = 8
    end if

    write (nttyo,1030) uveelu(1:i),ustelu(1:j),uplatm(1:k)
    write (noutpt,1030) uveelu(1:i),ustelu(1:j),uplatm(1:k)
1030 format('   EQLIBU (EQ3/6-V',a,'-EQLIBU-LIB-',a,'-',a,')',/)

    write (nttyo,1040)
    write (noutpt,1040)
1040 format(' Copyright (c) 1987, 1990-1993, 1995, 1997, 2002 The',' Regents of the',/' University of California, Lawrence',' Livermore National Laboratory.',/' All rights reserved.',/)

    ! Write additional statements and disclaimers.
    call prcndi(noutpt,nttyo)

    ! Write the time and date on the output.
    j2 = ilnobl(udate0)
    write(noutpt,1070) utime0,udate0(1:j2)
    write(nttyo,1070) utime0,udate0(1:j2)
1070 format(' Run',2x,a8,2x,a,//)

    ! Get the platform's real*8 floating-point parameters.
    call flpars(eps100,irang,noutpt,nttyo,smp100)

    ! Initialize array dimension variables.
    iapxmx = iapxpa
    ibpxmx = ibpxpa

    ! Write header on the slist file.
    write (nslist,1090)
1090 format(' EQPT Species List (SLIST) File:',//)

    ! Initialize the cumulative error counter.
    nerr = 0

    ! Initialize the cumulative warning counter.
    nwarn = 0

    ! Check the first line of the data file to ensure that the
    ! required header is present.
    call hdrchk(ndat0s,noutpt,nttyo)

    ! Get the aqueous species activity coefficient model type
    ! (Extended Debye-Huckel, Pitzer, etc.) associated with this
    ! data file. This subroutine scans the data file and counts
    ! keywords to make the determination.
    call gakey(ndat0s,noutpt,nttyo,uakey)

    ! Scan the data file to determine the necessary array dimensions
    ! and any corresponding structural information that will be
    ! required before "reading" the data file.
    ! First scan the data file title for embedded data.
    call ggridp(ipch_asv,ipcv_asv,itgenf,jpdblo,jpfc_asv,jptffl,narx_asv,ndat0s,ndb_asv,noutpt,ntid_asv,ntpr_asv,nttyo,q500fl,uakey)

    ! Set the number of parameters in a Pitzer alpha set.
    ipbt_asv = 1

    if (uakey(1:8) .eq. 'Pitzer  ') then
        ipbt_asv = 2
    end if

    ! Rescan the data file to get other array dimensions.
    ! Also get the actual numbers of basis species and chemical
    ! elements on the data file.
    call gnenb(ipbt_asv,ikt_asv,jpdblo,jpfc_asv,nap_asv,nat_asv,nazt_asv,nbt_asv,nct_asv,ndat0s,ngt_asv,nlt_asv,nmt_asv,noutpt,npt_asv,npx2_asv,npx3_asv,nsb,nst_asv,nttyo,nxt_asv,uakey)

    nbt = nbt_asv
    nct = nct_asv

    if (nsb .ne. (nct + 1)) then
        write (ux8a,'(i5)') nsb
        call lejust(ux8a)
        j2 = ilnobl(ux8a)
        write (ux8b,'(i5)') nct
        call lejust(ux8b)
        j3 = ilnobl(ux8b)
        write(noutpt,1100) ux8a(1:j2),ux8b(1:j3)
        write(nttyo,1100) ux8a(1:j2),ux8b(1:j3)
1100 format(/' * Error - (EQPT/eqpt) The number of strict basis',' species',/7x,'must be the number of chemical elements plus',' one. The number of',/7x,'strict basis species is ',a,',',' but the number of chemical elements',/7x,'is ',a,'.')

        stop
    end if

    ipbt_asv = 1

    if (uakey(1:8) .eq. 'Pitzer  ') then
        ipbt_asv = 2

        ! The following array allocation variables are only used by
        ! EQ3NR and EQ6.
        nlat_asv = npx2_asv + npx3_asv
        nmut_asv = npx2_asv + 2*npx3_asv
        nlat_asv = max(1,nlat_asv)
        nmut_asv = max(1,nmut_asv)
    end if

    iktmax = ikt_asv
    ipbtmx = ipbt_asv
    ipchmx = ipch_asv
    ipcvmx = ipcv_asv
    jpfcmx = jpfc_asv
    narxmx = narx_asv
    natmax = nat_asv
    naztmx = nazt_asv
    nbtmax = nbt_asv
    nctmax = nct_asv
    ndbmax = ndb_asv
    ngtmax = ngt_asv
    nltmax = nlt_asv
    nmtmax = nmt_asv
    npx2mx = npx2_asv
    npx3mx = npx3_asv
    ntidmx = ntid_asv
    ntprmx = ntpr_asv
    nxtmax = nxt_asv

    nbtmx1 = nbtmax + 1
    nbtmx2 = nbtmax + 2

    ! Allocate the necessary arrays.
    ALLOCATE(ipivot(narx_asv))
    ALLOCATE(issot(nxt_asv))
    ALLOCATE(nacdpr(ntpr_asv))
    ALLOCATE(narxt(ntpr_asv))
    ALLOCATE(nentei(nct_asv))
    ALLOCATE(nentri(nbt_asv + 1))

    ALLOCATE(qpdaz(nat_asv))

    ALLOCATE(uaqsp(nat_asv))
    ALLOCATE(zaqsp(nat_asv))

    ALLOCATE(ugassp(ngt_asv))
    ALLOCATE(uliqsp(nlt_asv))
    ALLOCATE(uminsp(nmt_asv))
    ALLOCATE(ussoph(nxt_asv))

    ALLOCATE(utitld(ntid_asv))

    ALLOCATE(uspec(nbt_asv + 1))
    ALLOCATE(udrsi(nbt_asv + 1))

    ALLOCATE(ussosp(ikt_asv,nxt_asv))

    ALLOCATE(uelem(nct_asv))
    ALLOCATE(uessi(nct_asv))

    ALLOCATE(atwt(nct_asv))
    ALLOCATE(cdrsi(nbt_asv + 1))
    ALLOCATE(cessi(nct_asv))
    ALLOCATE(zchar(nbt_asv + 1))
    ALLOCATE(mtotr(nct_asv))

    ALLOCATE(cdrs(nbt_asv + 2,nbt_asv + 1))
    ALLOCATE(cess(nct_asv,nbt_asv + 1))

    ALLOCATE(dhfs(narx_asv,ntpr_asv,ipch_asv,nbt_asv + 1))
    ALLOCATE(dvfs(narx_asv,ntpr_asv,ipcv_asv,nbt_asv + 1))
    ALLOCATE(xhfs(narx_asv,ntpr_asv,nbt_asv + 1))
    ALLOCATE(xlks(narx_asv,ntpr_asv,nbt_asv + 1))
    ALLOCATE(xvfs(narx_asv,ntpr_asv,nbt_asv + 1))

    ALLOCATE(adh(narx_asv,ntpr_asv))
    ALLOCATE(adhh(narx_asv,ntpr_asv))
    ALLOCATE(adhv(narx_asv,ntpr_asv))
    ALLOCATE(aphi(narx_asv,ntpr_asv))
    ALLOCATE(bdh(narx_asv,ntpr_asv))
    ALLOCATE(bdhh(narx_asv,ntpr_asv))
    ALLOCATE(bdhv(narx_asv,ntpr_asv))
    ALLOCATE(bdot(narx_asv,ntpr_asv))
    ALLOCATE(bdoth(narx_asv,ntpr_asv))
    ALLOCATE(bdotv(narx_asv,ntpr_asv))
    ALLOCATE(prehw(narx_asv,ntpr_asv))
    ALLOCATE(presg(narx_asv,ntpr_asv))
    ALLOCATE(xhfe(narx_asv,ntpr_asv))
    ALLOCATE(xlke(narx_asv,ntpr_asv))
    ALLOCATE(xvfe(narx_asv,ntpr_asv))

    ALLOCATE(dadhh(narx_asv,ntpr_asv,ipch_asv))
    ALLOCATE(dbdhh(narx_asv,ntpr_asv,ipch_asv))
    ALLOCATE(dbdth(narx_asv,ntpr_asv,ipch_asv))
    ALLOCATE(dhfe(narx_asv,ntpr_asv,ipch_asv))

    ALLOCATE(dadhv(narx_asv,ntpr_asv,ipcv_asv))
    ALLOCATE(dbdhv(narx_asv,ntpr_asv,ipcv_asv))
    ALLOCATE(dbdtv(narx_asv,ntpr_asv,ipcv_asv))
    ALLOCATE(dvfe(narx_asv,ntpr_asv,ipcv_asv))

    ALLOCATE(apr(narx_asv,ntpr_asv))
    ALLOCATE(avgrid(narx_asv,ntpr_asv))
    ALLOCATE(tempc(narx_asv,ntpr_asv))
    ALLOCATE(tempcs(narx_asv,ntpr_asv))
    ALLOCATE(tmpcmx(ntpr_asv))
    ALLOCATE(udbval(ndb_asv))
    ALLOCATE(xdbval(ndb_asv))

    ALLOCATE(tvec(narx_asv))
    ALLOCATE(tvecs(narx_asv))
    ALLOCATE(cof(narx_asv))
    ALLOCATE(xvec(narx_asv))
    ALLOCATE(yvec(narx_asv))

    ALLOCATE(aamatr(narx_asv,narx_asv))
    ALLOCATE(gmmatr(narx_asv,narx_asv))

    if (uakey(1:8) .eq. 'SEDH    ') then
        ALLOCATE(uazp(nazt_asv))
        ALLOCATE(azero(nazt_asv))
        ALLOCATE(insgf(nazt_asv))
    end if

    if (uakey(1:8) .eq. 'Pitzer  ') then
        ALLOCATE(alpha(ipbt_asv,npx2_asv))
        ALLOCATE(abeta(jpfc_asv,0:ipbt_asv,npx2_asv))
        ALLOCATE(acphi(jpfc_asv,npx2_asv))

        ALLOCATE(upair(2,npx2_asv))
        ALLOCATE(uthdtr(3,npx3_asv))
        ALLOCATE(utripl(3,npx3_asv))

        ALLOCATE(atheta(jpfc_asv,npx3_asv))
        ALLOCATE(apsi(jpfc_asv,npx3_asv))
    end if

    ! Write some tallies on the OUTPUT file.
    write (noutpt,1120) nct,nbt
    write (nslist,1120) nct,nbt
1120 format(11x,'Number of elements = ',i5,/11x,'Number of basis species = ',i5,/)

    ! Write the header on the DATA1 and DATA1f files. This begins with
    ! the string 'data1' at the top of the files, followed by the
    ! keystring for the type of aqueous species activity coefficient
    ! model, and the array dimensions necessary to read the rest of
    ! the DATA1 file.
    call wrhdr(ikt_asv,ipch_asv,ipcv_asv,jpfc_asv,nap_asv,narx_asv,nat_asv,nbt_asv,nct_asv,ndata1,ndat1f,ngt_asv,nlat_asv,nlt_asv,nmt_asv,nmut_asv,npt_asv,nst_asv,ntid_asv,ntpr_asv,nxt_asv,uakey)

    ! Write the data file title on the various output files. Determine
    ! certain embedded options apart from those that set dimensioning
    ! parameters.
    ntitld = ntid_asv
    call rdwttl(ipch,ipcv,jpdblo,jpfcmx,jptffl,narxt,ndata1,ndat0s,ndat1f,noutpt,nslist,ntitld,ntidmx,ntprmx,ntprt,nttyo,uakey,utitld)

    ! Read the miscellaneous parameters (write them after the chemical
    ! elements block). Write the nominal temperature limits and the
    ! upper limits of the temperature ranges here, however.
    call rdpar(adh,adhh,adhv,aphi,bdh,bdhh,bdhv,bdot,bdoth,bdotv,cco2,dadhh,dadhv,dbdhh,dbdhv,dbdth,dbdtv,dhfe,dvfe,ipch,ipchmx,ipcv,ipcvmx,itgenf,nacdpr,narxmx,narxt,ndat0s,ndbmax,ndbptg,ndbptl,nerr,noutpt,ntprmx,ntprt,nttyo,nwarn,prehw,presg,q500fl,tdamax,tdamin,tempc,uakey,udbfmt,udbval,xdbval,xhfe,xlke,xvfe)

    ! Read and write the element data.
    call rdwele(atwt,nch,nco,nct,nctmax,ndata1,ndat0s,ndat1f,nerr,noutpt,nslist,nttyo,uelem)

    ! Check the element names for uniqueness.
    call neleck(nct,nctmax,nerr,noutpt,nttyo,uelem)

    if (nerr .gt. 0) then
        write (ux8,'(i5)') nerr
        call lejust(ux8)
        j2 = ilnobl(ux8)
        write (noutpt,1160) ux8(1:j2)
        write (nttyo,1160) ux8(1:j2)
1160 format(//' ',a,' errors(s) were encountered.',/)

        stop
    end if

    ! Scale the grid temperature values for subsequent interpolation.
    do ntpr = 1,ntprt
        nmax = narxt(ntpr)

        do n = 1,nmax
            tvec(n) = tempc(n,ntpr)
        end do

        ! Calling sequence substitutions:
        !   tvec for avx
        !   tvecmx for avxmax
        !   tvecs for avxs
        call scalx1(tvec,tvecmx,tvecs,ier,nmax)

        do n = 1,nmax
            tempcs(n,ntpr) = tvecs(n)
        end do

        tmpcmx(ntpr) = tvecmx
    end do

    ! Set the interval parameter (nmodwr) for writing species and
    ! phase names to nttyo. This is scaled to the number of chemical
    ! elements.
    if (nct .gt. 50) then
        nmodwr = 20
    else if (nct .gt. 30) then
        nmodwr = 10
    else if (nct .gt. 20) then
        nmodwr = 5
    else
        nmodwr = 1
    end if

    ! Interpolate and write the miscellaneous parameters.
    call wrpar(aamatr,adh,adhh,adhv,aphi,apr,avgrid,bdh,bdhh,bdhv,bdot,bdoth,bdotv,cco2,cof,dadhh,dadhv,dbdhh,dbdhv,dbdth,dbdtv,dhfe,dvfe,eps100,gmmatr,ipch,ipchmx,ipcv,ipcvmx,ipivot,narxmx,narxt,ndata1,ndat1f,noutpt,ntprmx,ntprt,nttyo,presg,prehw,tdamax,tdamin,tempc,tempcs,tmpcmx,uakey,xhfe,xlke,xvfe,xvec,yvec)

    ! Read and write aqueous species.
    call pcraq(aamatr,apr,atwt,avgrid,cdrs,cdrsi,cess,cessi,cof,dhfe,dhfs,dvfe,dvfs,eps100,gmmatr,ipch,ipchmx,ipcv,ipcvmx,ipivot,itgenf,mtotr,nacdpr,narxmx,narxt,nat,natmax,nbt,nbtmx1,nbtmx2,nch,nco,nct,nctmax,ndata1,ndat0s,ndat1f,ndbmax,ndbptg,ndbptl,nentei,nentri,nerr,nmodwr,noutpt,nsb,nslist,ntprmx,ntprt,nttyo,nwarn,qelect,q500fl,tempc,tempcs,tmpcmx,uaqsp,udbfmt,udbval,udrsi,uelem,uessi,uspec,xdbval,xhfe,xhfs,xlke,xlks,xvfe,xvfs,xvec,yvec,zaqsp,zchar)

    ! Read and write minerals, liquids, gases.
    call pcrsg(aamatr,apr,atwt,avgrid,cdrs,cdrsi,cess,cessi,cof,dhfe,dhfs,dvfe,dvfs,eps100,gmmatr,ipch,ipchmx,ipcv,ipcvmx,ipivot,itgenf,mtotr,nacdpr,narxmx,narxt,nat,natmax,nbt,nbtmx1,nbtmx2,nch,nco,nct,nctmax,ndata1,ndat0s,ndat1f,ndbmax,ndbptg,ndbptl,nentei,nentri,nerr,ngt,ngtmax,nlt,nltmax,nmodwr,nmt,nmtmax,noutpt,nsb,nslist,ntprmx,ntprt,nttyo,nwarn,qelect,q500fl,tempc,tempcs,tmpcmx,udbfmt,udbval,udrsi,uelem,uessi,ugassp,uliqsp,uminsp,uspec,xdbval,xhfe,xhfs,xlke,xlks,xvfe,xvfs,xvec,yvec,zchar)

    ! Read and write solid solutions.
    call pcrss(apx,bpx,iapxmx,ibpxmx,iktmax,issot,nbtmx1,ndata1,ndat0s,ndat1f,nerr,nmodwr,nmt,nmtmax,noutpt,nslist,nttyo,nxt,nxtmax,uminsp,ussosp,ussoph)

    ! Make sure that all aqueous species names are unique.
    call naqsck(nat,natmax,nerr,noutpt,nttyo,uaqsp)

    ! Make sure that all pure mineral names are unique.
    call nminck(nerr,nmt,nmtmax,noutpt,nttyo,uminsp)

    ! Make sure that all gas species names are unique.
    call ngasck(nerr,ngt,ngtmax,noutpt,nttyo,ugassp)

    ! Make sure that all solid solution names are unique.
    call nssock(nerr,noutpt,nttyo,nxt,nxtmax,ussoph)

    ! Make sure that all solid solution end-member names are unique
    ! within each solid solution.
    call nxspck(iktmax,issot,nerr,noutpt,nttyo,nxt,nxtmax,ussoph,ussosp)

    ! Validate the end-member names for the solid solutions. Each
    ! end-member name must match that of a pure mineral.
    call vxspck(iktmax,issot,nerr,nmt,nmtmax,noutpt,nttyo,nxt,nxtmax,uminsp,ussoph,ussosp)

    if (nerr .gt. 0) then
        write (ux8,'(i5)') nerr
        call lejust(ux8)
        j2 = ilnobl(ux8)
        write (noutpt,1160) ux8(1:j2)
        write (nttyo,1160) ux8(1:j2)
        stop
    end if

    ! Rewind and position to read aqueous species activity coefficient
    ! parameters.
    rewind(ndat0s)

    if (uakey(1:8) .eq. 'SEDH    ') then
        uacfst = 'bdot parameters '
    else if (uakey(1:8) .eq. 'Pitzer  ') then
        if (jpdblo .eq. -1) then
            uacfst = 'single-salt para'
        else
            uacfst = 'ca combinations:'
        end if
    else
        j2 = ilnobl(uakey)
        write(noutpt,1200) uakey(1:j2)
        write(nttyo,1200) uakey(1:j2)
1200 format(/' * Error - (EQPT/eqpt) Unrecognized data file',' type = "',a,'".',/7x,'Allowed types are "SEDH" (Simple',' Extended Debye-Huckel) and',/7x,'"Pitzer".')

        stop
    end if

270 continue
    read(ndat0s,1210,end=280,err=280) uline
1210 format(a)

    if (uline(1:16) .ne. uacfst(1:16)) then
        go to 270
    end if

    go to 290

280 continue
    j3 = ilnobl(uacfst)
    write(noutpt,1220) uacfst(1:j3)
    write(nttyo,1220) uacfst(1:j3)
1220 format(/' * Error - (EQPT/eqpt) End-of-file hit or other read',/7x,'error occurred while searching for the block header',/7x,'beginning with "',a,'".')

    stop

290 continue
    backspace(ndat0s)

    ! Zero arrays for activity coefficient parameters.
    if (uakey(1:8) .eq. 'SEDH    ') then
        ! Zero the azero array.
        call initaz(azero,naztmx)

        ! Zero the insgf array.
        call initiz(insgf,naztmx)
    end if

    if (uakey(1:8) .eq. 'Pitzer  ') then
        ! Zero the alpha array.
        nmax = ipbtmx*npx2mx
        call initaz(alpha,nmax)

        ! Zero the abeta array.
        nmax = jpfcmx*(ipbtmx + 1)*npx2mx
        call initaz(abeta,nmax)

        ! Zero the acphi array.
        nmax = jpfcmx*npx2mx
        call initaz(acphi,nmax)

        ! Zero the atheta and apsi arrays.
        nmax = jpfcmx*npx3mx
        call initaz(atheta,nmax)
        call initaz(apsi,nmax)
    end if

    ! Read the SEDH or Pitzer section.
    if (uakey(1:8) .eq. 'SEDH    ') then
        ! Read hard core diameters and related parameters used
        ! in the B-dot equation.
        call rdazp(azero,insgf,nazt,naztmx,ndat0s,nerr,noutpt,nttyo,uazp)
    end if

    if (uakey(1:8) .eq. 'Pitzer  ') then
        ! Read Pitzer parameters.
        if (jpdblo .eq. -1) then
            ! Read the data according to the old Pitzer data block
            ! organization. There is one superblock for data
            ! pertaining to pure aqueous electrolytes, another
            ! for mixtures of such with a common ion. Parameters
            ! involving neutral aqueous species are handled (awkwardly,
            ! but with appropriate checks) within this structure as
            ! noted below.
            ! Read the beta(n) and Cphi parameters.
            !   - normally for ca (cation-anion) combinations.
            !   - lamda(0) values may be read in place of beta(0) values
            !     for nn, nn', nc, and na combinations.
            !   - mu values may be read in place of Cphi values for
            !     nnn combinations.
            call rdpz2(abeta,alpha,acphi,ipbtmx,jpfcmx,nat,natmax,ndat0s,nerr,noutpt,npx2mx,npx2t,nttyo,nwarn,uaqsp,upair,zaqsp)

            ! Read the S-theta and psi parameters.
            !   - normally for cc'a and aa'c combinations.
            !   - zeta values may be read in place of psi values for
            !     nca combinations.
            !   - mu values may be read in place of psi values for
            !     nnn' and n'n'n combinations.
            call rdpz3(apsi,atheta,jpfcmx,nat,natmax,ndat0s,nerr,noutpt,npx3mx,npx3t,nthdt,nttyo,nwarn,uaqsp,uethfl,uthdtr,utripl,zaqsp)
        else
            ! Read the data according to the new Pitzer data block
            ! organization. There is one superblock for each allowed
            ! species pair or triplet type. Note that for the "new"
            ! organization, the E-theta flag is hardwired to 'on'.
            ! That is because basically all modern implementations
            ! of Pitzer's equations use the theoretically-based
            ! higher-order electrostatic term formalism (which is
            ! turned on by uethfl = 'on').
            uethfl = 'on'
            npx2t = 0
            npx3t = 0

            ! Read ca (cation-anion) data.
            call rdpca(abeta,alpha,acphi,ipbtmx,jpfcmx,nat,natmax,ndat0s,nerr,noutpt,npxca,npx2mx,npx2t,nttyo,nwarn,uaqsp,upair,zaqsp)

            ! Read cc' (cation-different cation) and aa' (anion-
            ! different anion) data.
            call rdpth(abeta,ipbtmx,jpfcmx,nat,natmax,ndat0s,nerr,noutpt,npxth,npx2mx,npx2t,nttyo,nwarn,uaqsp,upair,zaqsp)

            ! Read nc (neutral-cation) and na (neutral-anion) data.
            call rdpni(abeta,ipbtmx,jpfcmx,nat,natmax,ndat0s,nerr,noutpt,npxni,npx2mx,npx2t,nttyo,nwarn,uaqsp,upair,zaqsp)

            ! Read nn (neutral-same neutral) data.
            call rdpn2(abeta,acphi,ipbtmx,jpfcmx,nat,natmax,ndat0s,nerr,noutpt,npxn2,npx2mx,npx2t,nttyo,nwarn,uaqsp,upair,zaqsp)

            ! Read nn'(neutral-different neutral) data.
            call rdpnn(abeta,ipbtmx,jpfcmx,nat,natmax,ndat0s,nerr,noutpt,npxnn,npx2mx,npx2t,nttyo,nwarn,uaqsp,upair,zaqsp)

            ! Read cc'a (cation-different cation-anion) and aa'c
            ! (anion-different anion-cation) data.
            call rdppsi(apsi,jpfcmx,nat,natmax,ndat0s,nerr,noutpt,npxpsi,npx3mx,npx3t,nttyo,nwarn,uaqsp,utripl,zaqsp)

            ! Read nca (neutral-cation--anion) data.
            call rdpzet(apsi,jpfcmx,nat,natmax,ndat0s,nerr,noutpt,npxzet,npx3mx,npx3t,nttyo,nwarn,uaqsp,utripl,zaqsp)

            ! Read nnn' (neutral-neutral-different neutral) data.
            call rdpn2n(apsi,jpfcmx,nat,natmax,ndat0s,nerr,noutpt,npxn2n,npx3mx,npx3t,nttyo,nwarn,uaqsp,utripl,zaqsp)

            ! Copy the theta data into arrays used with the old
            ! Pitzer data block organization.
            nthdt = npxth
            nthd = 0
            npx2r1 = npxca + 1
            npx2r2 = npxca + npxth

            do npx2 = npx2r1,npx2r2
                nthd = nthd + 1
                uthdtr(1,nthd) = upair(1,npx2)
                uthdtr(2,nthd) = upair(2,npx2)
                uthdtr(3,nthd) = '<dummy>'

                do j = 1,jpfcmx
                    atheta(j,nthd) = abeta(j,0,npx2)
                end do
            end do
        end if
    end if

    if (nerr .gt. 0) then
        write (ux8,'(i5)') nerr
        call lejust(ux8)
        j2 = ilnobl(ux8)
        write (noutpt,1160) ux8(1:j2)
        write (nttyo,1160) ux8(1:j2)
        stop
    end if

    ! Test and write the SEDH data section of DATA1.
    if (uakey(1:8) .eq. 'SEDH    ') then
        ! Validate the names of aqueous species used to specify azero
        ! ('bdot') data. Each such name should correspond to an aqueous
        ! species for which there is a species block on the data file.
        ! Write a note for any exceptions.
        call vazpck(nat,natmax,nazt,naztmx,noutpt,nttyo,uaqsp,uazp)

        ! Test and process all azero and insgf ('bdot') data.
        call tpraz(nat,natmax,nazt,naztmx,ncvaz,nerr,noutpt,nttyo,pcvaz,qpdaz,uaqsp,uazp)

        ! Write the azero and insgf data summary.
        ux8 = ' '
        write (ux8,'(f6.2)') pcvaz
        call lejust(ux8)
        j2 = ilnobl(ux8)
        write (noutpt,1230)
1230 format(//2x,'Aqueous species hard core diameter coverage:')

        natm1 = nat - 1
        write (noutpt,1240) ncvaz,natm1,ux8(1:j2)
1240 format(/4x,i5,' aqueous species have hard core diameters',' specified on the data file',/4x,i5,' aqueous solute',' species are present on this file',/4x,'Coverage is ',a,' per cent')

        if (nerr .gt. 0) then
            write (ux8,'(i5)') nerr
            call lejust(ux8)
            j2 = ilnobl(ux8)
            write (noutpt,1160) ux8(1:j2)
            write (nttyo,1160) ux8(1:j2)
            stop
        end if

        ! Write hard core diameters and related parameters used
        ! for example in the B-dot equation.
        call wrazp(azero,insgf,nazt,naztmx,ndata1,ndat1f,noutpt,nttyo,uazp)
    end if

    ! Test and write the Pitzer data section of DATA1.
    if (uakey(1:8) .eq. 'Pitzer  ') then
        ! Count the numbers of aqueous solution cations, anions, and
        ! neutrals, excluding any fictive redox species. These will
        ! be used to compute the species pairs and triplets relevant
        ! to Pitzer parameters.
        call coasst(jassan,jassca,jassne,nat,natmax,uaqsp,zaqsp)

        ! Compute the number of each of the relevant types of these
        ! species pairs and triplets (e.g., ca, cc', aa', nc, na, nn,
        ! nn', cca, aac, nnn, nnn', and nca).
        call cpcomb(jassan,jassca,jassne,naapr,nat,ncapr,nccpr,nnapr,nncpr,nnnpr,nn2pr,naactr,na2ctr,nncatr,nccatr,nc2atr,nn2ntr,nn3tr)

        ! Allocate the associated arrays.
        ALLOCATE(alphca(ipbt_asv,ncapr))

        ALLOCATE(alamaa(jpfc_asv,0:ipbt_asv,naapr))
        ALLOCATE(alamca(jpfc_asv,0:ipbt_asv,ncapr))
        ALLOCATE(alamcc(jpfc_asv,0:ipbt_asv,nccpr))
        ALLOCATE(alamna(jpfc_asv,0:ipbt_asv,nnapr))
        ALLOCATE(alamnc(jpfc_asv,0:ipbt_asv,nncpr))
        ALLOCATE(alamnn(jpfc_asv,0:ipbt_asv,nnnpr))
        ALLOCATE(alamn2(jpfc_asv,0:ipbt_asv,nn2pr))

        ALLOCATE(amua2c(jpfc_asv,na2ctr))
        ALLOCATE(amuc2a(jpfc_asv,nc2atr))
        ALLOCATE(amun3(jpfc_asv,nn3tr))

        ALLOCATE(amuaac(jpfc_asv,naactr))
        ALLOCATE(amucca(jpfc_asv,nccatr))
        ALLOCATE(amunca(jpfc_asv,nncatr))
        ALLOCATE(amun2n(jpfc_asv,nn2ntr))

        ALLOCATE(qpdaa(naapr))
        ALLOCATE(qpdca(ncapr))
        ALLOCATE(qpdcc(nccpr))
        ALLOCATE(qpdna(nnapr))
        ALLOCATE(qpdnc(nncpr))
        ALLOCATE(qpdnn(nnnpr))
        ALLOCATE(qpdn2(nn2pr))

        ALLOCATE(qpdaac(naactr))
        ALLOCATE(qpdcca(nccatr))
        ALLOCATE(qpdnca(nncatr))
        ALLOCATE(qpdn2n(nn2ntr))

        ALLOCATE(in2pr(nn2pr))
        ALLOCATE(in3tr(nn3tr))

        ALLOCATE(iaapr(2,naapr))
        ALLOCATE(icapr(2,ncapr))
        ALLOCATE(iccpr(2,nccpr))
        ALLOCATE(inapr(2,nnapr))
        ALLOCATE(incpr(2,nncpr))
        ALLOCATE(innpr(2,nnnpr))

        ALLOCATE(in2ntr(2,nn2ntr))
        ALLOCATE(ia2ctr(2,na2ctr))
        ALLOCATE(ic2atr(2,nc2atr))
        ALLOCATE(iaactr(3,naactr))
        ALLOCATE(iccatr(3,nccatr))
        ALLOCATE(incatr(3,nncatr))

        ! Construct index arrays for those pairs and triplets.
        call bldspc(iaapr,icapr,iccpr,inapr,incpr,innpr,in2pr,iaactr,ia2ctr,iccatr,ic2atr,incatr,in2ntr,in3tr,jassan,jassca,jassne,nat,natmax,naapr,ncapr,nccpr,nnapr,nncpr,nnnpr,nn2pr,naactr,na2ctr,nncatr,nccatr,nc2atr,nn2ntr,nn3tr,uaqsp,zaqsp)

        ! Note on subroutine and array naming for species pairs and
        ! triplets associated with Pitzer coefficients:
        !   ca   is represented by "ca"
        !   cc'a is represented by "cca"
        !   cca  is represented by "c2a"
        !   aa'c is represented by "a2c"
        !   aac  is represented by "aac"
        !   cc' is represented by "cc"
        ! Test and process all ca (cation-anion) parameters. This includes
        ! the cca (repeated cation-anion) and aac (cation-repeated anion)
        ! cases. The basic transformations are:
        !   beta(n)(ca) -> lambda(n)(ca)   (n = 0,2)
        !   Cphi(ca)    -> mu(cca) and mu(aac)
        call tprca(abeta,acphi,alamca,alpha,alphca,amua2c,amuc2a,icapr,ipbtmx,jpfcmx,natmax,na2ctr,ncapr,ncvca,nc2atr,nerr,noutpt,npx2mx,npx2t,nttyo,nwarn,pcvca,qpdca,uaqsp,upair,zaqsp)

        ! Test and process the cc' (cation-different cation) parameters.
        ! cases. The basic transformation is:
        !   theta(cc') -> lambda(cc')
        call tprcc(alamcc,atheta,iccpr,ipbtmx,jpfcmx,natmax,nccpr,ncvcc,nerr,noutpt,npx3mx,nthdt,nttyo,nwarn,pcvcc,qpdcc,uaqsp,uthdtr)

        ! Test and process the aa' (anion-different anion) parameters.
        ! cases. The basic transformation is:
        !   theta(aa') -> lambda(aa')
        call tpraa(alamaa,atheta,iaapr,ipbtmx,jpfcmx,naapr,natmax,ncvaa,nerr,noutpt,npx3mx,nthdt,nttyo,nwarn,pcvaa,qpdaa,uaqsp,uthdtr)

        ! Test and process the nn (repeated-neutral) parameters. This
        ! includes the nnn (doubly repeated neutral) cases. Note on
        ! subroutine and array naming:
        !   nn  is represented by "n2"
        !   nn' is represented by "nn"
        !   lambda(nn) -> lambda(nn)
        !   mu(nnn)    -> mu(nnn)
        call tprn2(abeta,acphi,alamn2,amun3,in2pr,ipbtmx,jpfcmx,natmax,ncvn2,nerr,nn2pr,nn3tr,noutpt,npx2mx,npx2t,nttyo,nwarn,pcvn2,qpdn2,uaqsp,upair)

        ! Test and process the nn' (neutral-different neutral) parameters.
        ! See the above note on subroutine and array naming.
        !   lambda(nn') -> lambda(nn')
        call tprnn(abeta,acphi,alamnn,innpr,in2pr,ipbtmx,jpfcmx,natmax,ncvnn,nerr,nnnpr,nn2pr,noutpt,npx2mx,npx2t,nttyo,nwarn,pcvnn,qpdnn,qpdn2,uaqsp,upair)

        ! Test and process the nc (neutral-cation) parameters.
        !   lambda(nc) -> lambda(nc)
        call tprnc(abeta,alamnc,incpr,ipbtmx,jpfcmx,natmax,ncvnc,nerr,nncpr,noutpt,npx2mx,npx2t,nttyo,nwarn,pcvnc,qpdnc,uaqsp,upair)

        ! Test and process the na (neutral-anion) parameters.
        !   lambda(na) -> lambda(na)
        call tprna(abeta,alamna,inapr,ipbtmx,jpfcmx,natmax,ncvna,nerr,nnapr,noutpt,npx2mx,npx2t,nttyo,nwarn,pcvna,qpdna,uaqsp,upair)

        ! Test and process the cc'a (cation-different cation-anion)
        ! parameters. Note that mu(cc'a) depends not only on psi(cc'a),
        ! but also on mu(cca) and mu(c'c'a), hence indirectly on Cphi(ca)
        ! and Cphi(c'a). Note on subroutine and array naming:
        !   cc'a is represented by "cca"
        !   cca  is represented by "c2a"
        !   psi(cc'a) -> mu(cc'a)
        call tprcca(amucca,amuc2a,apsi,icapr,iccatr,ipbtmx,jpfcmx,natmax,ncapr,nccatr,ncvcca,nc2atr,nerr,noutpt,npx3mx,npx3t,nttyo,nwarn,pcvcca,qpdca,qpdcca,uaqsp,utripl,zaqsp)

        ! Test and process the aa'c (anion-different anion-cation)
        ! parameters. Note that mu(aa'c) depends not only on psi(aa'c),
        ! but also on mu(aac) and mu(a'a'c), hence indirectly on Cphi(ca)
        ! and Cphi(ca'). Note on subroutine and array naming:
        !   aa'c is represented by "aac"
        !   aac  is represented by "a2c"
        !   psi(aa'c) -> mu(aa'c)
        call tpraac(amuaac,amua2c,apsi,iaactr,icapr,ipbtmx,jpfcmx,naactr,natmax,na2ctr,ncapr,ncvaac,nerr,noutpt,npx3mx,npx3t,nttyo,nwarn,pcvaac,qpdaac,qpdca,uaqsp,utripl,zaqsp)

        ! Test and process the nnn' (repeated neutral-different neutral)
        ! parameters. These are not processed with the nn' parameters
        ! because there are two mu coefficients (nnn' and n'n'n) for
        ! each lambda coefficient (nn'). Note on subroutine and array
        ! naming:
        !   nnn' is represented by "n2n"
        !   mu(nnn') -> mu(nnn')
        call tprn2n(amun2n,apsi,innpr,in2pr,in2ntr,ipbtmx,jpfcmx,natmax,ncvn2n,nerr,nnnpr,nn2pr,nn2ntr,noutpt,npx3mx,npx3t,nttyo,nwarn,pcvn2n,qpdnn,qpdn2,qpdn2n,uaqsp,utripl)

        ! Test and process the nca (neutral-cation-anion) parameters.
        !   mu(nca) -> mu(nca)
        call tprnca(amunca,apsi,inapr,incatr,incpr,ipbtmx,jpfcmx,natmax,ncvnca,nerr,nnapr,nncatr,nncpr,noutpt,npx3mx,npx3t,nttyo,nwarn,pcvnca,qpdna,qpdnca,qpdnc,uaqsp,utripl)

        ! What about the following combinations?
        !   cc     (repeated cation)
        !   aa     (repeated neutral)
        !   ccc    (doubly repeated cation)
        !   aaa    (doubly repeated cation)
        !   ccc'   (repeated cation and a distinct cation)
        !   aaa'   (repeated anion and a distinct anion)
        !   cc'c'' (three distinct cations)
        !   aa'a'' (three distinct anions)
        !   nn'n'' (three distinct neutrals)
        !   ncc'   (neutral and two distinct cations)
        !   naa'   (neutral and two distinct anions)
        !   nn'c   (two distinct neutrals and a cation)
        !   nn'a   (two distinct neutrals and an anion)
        ! Some (e.g., cc,  aa) are defined to be zero by convention.
        ! Others (e.g., cc'c'' and nn'n'') correspond to quaternary
        ! or higher mixtures. In any case, they are not used in the
        ! normal Pitzer treatment of electrolyte solutions.
        ! Write the cation-anion (ca) data summary.
        ux8 = ' '
        write (ux8,'(f6.2)') pcvca
        call lejust(ux8)
        j2 = ilnobl(ux8)
        write (noutpt,1300)
1300 format(//2x,'Cation-anion (ca) pair coverage:')

        write (noutpt,1310) ncvca,ncapr,ux8(1:j2)
1310 format(/4x,i5,' pairs have Pitzer parameters specified on',' the data file',/4x,i5,' pairs can be constructed from',' the species present on this file',/4x,'Coverage is ',a,' per cent')

        ! Write the cation-distinct cation (cc') data summary.
        ux8 = ' '
        write (ux8,'(f6.2)') pcvcc
        call lejust(ux8)
        j2 = ilnobl(ux8)
        write (noutpt,1320)
1320 format(//2x,"Cation-distinct cation (cc') pair coverage:")

        write (noutpt,1310) ncvcc,nccpr,ux8(1:j2)

        ! Write the anion-distinct anion (aa') data summary.
        ux8 = ' '
        write (ux8,'(f6.2)') pcvaa
        call lejust(ux8)
        j2 = ilnobl(ux8)
        write (noutpt,1330)
1330 format(//2x,"Anion-distinct anion (aa') pair coverage:")

        write (noutpt,1310) ncvaa,naapr,ux8(1:j2)

        ! Write the repeated-neutral (nn) data summary.
        ux8 = ' '
        write (ux8,'(f6.2)') pcvn2
        call lejust(ux8)
        j2 = ilnobl(ux8)
        write (noutpt,1340)
1340 format(//2x,'Repeated-neutral (nn) pair coverage:')

        write (noutpt,1310) ncvn2,nn2pr,ux8(1:j2)

        ! Write the neutral-distinct neutral (nn') data summary.
        ux8 = ' '
        write (ux8,'(f6.2)') pcvnn
        call lejust(ux8)
        j2 = ilnobl(ux8)
        write (noutpt,1350)
1350 format(//2x,"Neutral-distinct neutral (nn') pair coverage:")

        write (noutpt,1310) ncvnn,nnnpr,ux8(1:j2)

        ! Write the neutral-cation (nc) data summary.
        ux8 = ' '
        write (ux8,'(f6.2)') pcvnc
        call lejust(ux8)
        j2 = ilnobl(ux8)
        write (noutpt,1360)
1360 format(//2x,'Neutral-cation (nc) pair coverage:')

        write (noutpt,1310) ncvnc,nncpr,ux8(1:j2)

        ! Write the neutral-anion (na) data summary.
        ux8 = ' '
        write (ux8,'(f6.2)') pcvna
        call lejust(ux8)
        j2 = ilnobl(ux8)
        write (noutpt,1370)
1370 format(//2x,'Neutral-anion (na) pair coverage:')

        write (noutpt,1310) ncvna,nnapr,ux8(1:j2)

        ! Write the cation-distinct cation-anion (cc'a) data summary.
        ux8 = ' '
        write (ux8,'(f6.2)') pcvcca
        call lejust(ux8)
        j2 = ilnobl(ux8)
        write (noutpt,1380)
1380 format(//2x,"Cation-distinct cation-anion (cc'a) triplet",' coverage:')

        write (noutpt,1390) ncvcca,nccatr,ux8(1:j2)
1390 format(/4x,i5,' triplets have Pitzer parameters specified on',' the data file',/4x,i5,' triplets can be constructed from',' the species present on this file',/4x,'Coverage is ',a,' per cent')

        ! Write the anion-distinct anion-cation (aa'c) data summary.
        ux8 = ' '
        write (ux8,'(f6.2)') pcvaac
        call lejust(ux8)
        j2 = ilnobl(ux8)
        write (noutpt,1400)
1400 format(//2x,"Anion-distinct anion-cation (aa'c) triplet",' coverage:')

        write (noutpt,1390) ncvaac,naactr,ux8(1:j2)

        ! Write the repeated neutral-distinct neutral (nnn') data
        ! summary.
        ux8 = ' '
        write (ux8,'(f6.2)') pcvn2n
        call lejust(ux8)
        j2 = ilnobl(ux8)
        write (noutpt,1410)
1410 format(//2x,"Repeated neutral-distinct neutral (nnn')",' triplet coverage:')

        write (noutpt,1390) ncvn2n,nn2ntr,ux8(1:j2)

        ! Write the neutral-cation-anion data summary.
        ux8 = ' '
        write (ux8,'(f6.2)') pcvnca
        call lejust(ux8)
        j2 = ilnobl(ux8)
        write (noutpt,1420)
1420 format(//2x,'Neutral-cation-anion (nca) triplet coverage:')

        write (noutpt,1390) ncvnca,nncatr,ux8(1:j2)

        if (nerr .gt. 0) then
            write (ux8,'(i5)') nerr
            call lejust(ux8)
            j2 = ilnobl(ux8)
            write (noutpt,1160) ux8(1:j2)
            write (nttyo,1160) ux8(1:j2)
            stop
        end if

        ! Write Pitzer parameters.
        call wrpz23(alphca,alamaa,alamca,alamcc,alamna,alamnc,alamnn,alamn2,amuaac,amua2c,amucca,amuc2a,amunca,amun2n,amun3,ipbtmx,iaapr,icapr,iccpr,inapr,incpr,innpr,in2pr,iaactr,ia2ctr,iccatr,ic2atr,incatr,in2ntr,in3tr,jpdblo,jpfcmx,natmax,naapr,ncapr,nccpr,nnapr,nncpr,nnnpr,nn2pr,naactr,na2ctr,nccatr,nc2atr,nncatr,nn2ntr,nn3tr,ndata1,ndat1f,noutpt,nttyo,uaqsp,uethfl)
    end if

    j2 = ilnobl(uakey)
    j3 = ilnobl(utitld(1))
    write (nttyo,1460) uakey(1:j2),utitld(1)(1:j3)
    write (noutpt,1460) uakey(1:j2),utitld(1)(1:j3)
1460 format(//' Completed processing the ',a,' data file ',a,'.')

    write (noutpt,1470)
    write (nttyo,1470)
1470 format(//' No errors were encountered.')

    if (nwarn .le. 0) then
        write (noutpt,1480)
        write (nttyo,1480)
1480 format(/' No warnings were encountered.',//)
    else
        write (ux8,'(i5)') nwarn
        call lejust(ux8)
        j4 = ilnobl(ux8)
        write (noutpt,1490) ux8(1:j4)
        write (nttyo,1490) ux8(1:j4)
1490 format(/' ',a,' warning(s) were encountered.',//)
    end if

    ! Get end time and date. Also get the run time. Optionally, on
    ! a Unix platform, get the user and cpu times.
    call runtim(iexec0,jexec0,texec0,noutpt,nttyo,trun,tuser,tcpu,udate1,utime1)

    j2 = ilnobl(udate0)
    j3 = ilnobl(udate1)
    write (noutpt,1500) utime0,udate0(1:j2),utime1,udate1(1:j3)
    write (nttyo,1500) utime0,udate0(1:j2),utime1,udate1(1:j3)
1500 format(10x,'Start time = ',a8,2x,a,/12x,'End time = ',a8,2x,a)

    ! Print the run, user, and cpu times.
    write (noutpt,1510) trun
    write (nttyo,1510)  trun
1510 format(/10x,' Run time = ',g10.3,' seconds')

    if (tuser .gt. 0.) then
        write (noutpt,1520) tuser
        write (nttyo,1520)  tuser
1520 format(10x,'User time = ',g10.3,' seconds')
    end if

    if (tcpu .gt. 0.) then
        write (noutpt,1530) tcpu
        write (nttyo,1530)  tcpu
1530 format(10x,' Cpu time = ',g10.3,' seconds')
    end if

    write (noutpt,1540)
    write (nttyo,1540)
1540 format(/' Normal exit')

    ! Clear the IEEE flag for floating-point underflow, if such a
    ! flag is present, to avoid getting an unnecessary system
    ! warning message. Underflow is a normal condition in EQ3/6.
    ! Make porting changes in the EQLIBU subroutine that is called
    ! in this section. Do not make the porting changes here.
    call cliefu()

    ! Close files.
    close(ndat0s,status='delete')

    close(noutpt)
    close(nslist)
    close(ndata1)
    close(ndat1f)
end program eqpt