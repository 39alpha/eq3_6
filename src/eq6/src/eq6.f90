program eq6
    !! This is the main program of the EQ6 code. Configuration
    !! identification, the copyright statement, legal disclaimers,
    !! and similar statements are contained in EQ6/aaaeq6.f, the
    !! lead-off subroutine in the EQ3NR source code. A short description
    !! of this program is also contained in that subroutine.
    !! Modules.
    !! The module mod6pt contains data required to evaluate Pitzer's
    !! equations.
    use mod6pt

    ! The module mod6xf contains most of the standard-state
    ! thermodynamic data.
    use mod6xf

    implicit none

    include 'eqlib/eqlpar.h'

    include 'eqlib/eqldv.h'
    include 'eqlib/eqlge.h'
    include 'eqlib/eql1s.h'

    include 'eqlib/eqlo8.h'

    ! File path parameters
    integer :: numargs
    character(len=1024) :: temppath
    character(len=:), allocatable :: data1path
    character(len=:), allocatable :: sixipath
    integer :: pathindices(2)
    character(len=:), allocatable :: basename
    character(len=:), allocatable :: ofile
    character(len=:), allocatable :: pfile
    character(len=:), allocatable :: bafile
    character(len=:), allocatable :: bbfile
    character(len=:), allocatable :: ifile
    character(len=:), allocatable :: tfile
    character(len=:), allocatable :: txfile
    character(len=:), allocatable :: tsfile

    integer, parameter :: nllnpa = 7200

    ! Array allocation size variables used in EQ6.
    !   General size variables:
    !     ntid_asv  = the number of lines in the data file title
    !     ntit_asv  = the number of lines in any title appearing
    !                   on an input or pickup file
    !   Size variables associated with the description of thermodynamic
    !   data:
    !     ipch_asv  = the order for pressure corrections to enthalpy
    !                   functions
    !     ipcv_asv  = the order for pressure corrections to volume
    !                   functions; the maximum order for pressure
    !                   corrections to log K and other Gibbs-energy-
    !                   based functions is one greater than this
    !     ipbt_asv  = the maximum number of Pitzer alpha parameters
    !                   for any species pair
    !     jpfc_asv  = the number of coefficients in the Pitzer parameter
    !                   temperature function
    !     narx_asv  = maximum number of coefficients per temperature
    !                   range
    !     ntpr_asv  = number of temperature ranges
    !   Size variables required to read other data from the
    !   supporting data file:
    !   Size variables required to read other data from the
    !   supporting data file:
    !     ikta_asv  = the maximum number of end-member components
    !                   of any solid solution on the data file
    !     nata_asv  = the number of aqueous species on the data file
    !     nbta_asv  = the number of basis species on the data file
    !     ncta_asv  = the number of chemical elements on the data file
    !     ngta_asv  = the number of gas species on the data file
    !     nlta_asv  = the number of pure liquids on the data file
    !     nmta_asv  = the number of pure minerals on the data file
    !     napa_asv  = number of distinct sets of Pitzer alpha
    !                   parameters on the data file
    !     nmuta_asv = the number of triplets of ions on the data file
    !                   for which distinct Pitzer mu coefficients
    !                   are defined
    !     nslta_asv = the number of pairs of ions on the data file
    !                   for which distinct Pitzer lambda coefficients
    !                   are defined
    !   Size variables required to hold the compressed set of such
    !   data:
    !     ikt_asv   = the maximum number of end-member components
    !                   of any solid solution in the compressed set
    !     nat_asv   = the number of aqueous species in the compressed
    !                   set
    !     nbt_asv   = the number of basis species in the compressed set
    !     nct_asv   = the number of chemical elements in the compressed
    !                   set
    !     ngt_asv   = the number of gas species in the compressed set
    !     nlt_asv   = the number of pure liquids in the compressed set
    !     nmt_asv   = the number of pure minerals in the compressed set
    !     nap_asv   = the number of distinct sets of Pitzer alpha
    !                   parameters in the compressed set
    !     nmut_asv  = the number of triplets of ions in the compressed
    !                   set for which distinct Pitzer mu coefficients
    !                   are defined
    !     nslt_asv  = the number of pairs of ions in the compressed set
    !                   for which distinct Pitzer lambda coefficients
    !                   are defined
    integer :: ntid_asv

    integer :: ipch_asv
    integer :: ipcv_asv
    integer :: ipbt_asv
    integer :: jpfc_asv
    integer :: narx_asv
    integer :: ntpr_asv

    integer :: nazm_asv
    integer :: nazp_asv

    integer :: ikta_asv
    integer :: napa_asv
    integer :: nata_asv
    integer :: nbta_asv
    integer :: ncta_asv
    integer :: ndrsa_asv
    integer :: nessa_asv
    integer :: ngta_asv
    integer :: nlta_asv
    integer :: nmta_asv
    integer :: nmuta_asv
    integer :: nslta_asv
    integer :: npta_asv
    integer :: nsta_asv
    integer :: nxta_asv

    integer :: iapxa_asv
    integer :: ibpxa_asv

    integer :: nbta1_asv

    integer :: ntf1a_asv
    integer :: ntf2a_asv

    integer :: ikt_asv
    integer :: nap_asv
    integer :: nat_asv
    integer :: nbt_asv
    integer :: nct_asv
    integer :: ndrs_asv
    integer :: ness_asv
    integer :: ngt_asv
    integer :: nlt_asv
    integer :: nmt_asv
    integer :: nmut_asv
    integer :: nslt_asv
    integer :: npt_asv
    integer :: nst_asv
    integer :: nxt_asv

    integer :: iapx_asv
    integer :: ibpx_asv

    integer :: ntit_asv
    integer :: nxmd_asv

    integer :: nbt1_asv
    integer :: nmx_asv
    integer :: nsts_asv
    integer :: nsx_asv

    integer :: ntf1_asv
    integer :: ntf2_asv

    integer :: k_asv
    integer :: npet_asv
    integer :: nset_asv

    integer :: imch_asv
    integer :: ndct_asv
    integer :: nert_asv
    integer :: nffg_asv
    integer :: nprp_asv
    integer :: nprs_asv
    integer :: nptk_asv
    integer :: nrct_asv
    integer :: nsrt_asv
    integer :: nttk_asv
    integer :: nxop_asv
    integer :: nxpe_asv
    integer :: nxrt_asv

    integer :: nad1
    integer :: nbkupa
    integer :: nbkupb
    integer :: newin
    integer :: ninpt
    integer :: ninpts
    integer :: noutpt
    integer :: ntab
    integer :: ntabs
    integer :: ntabx
    integer :: nttyo

    integer :: iodb(nodb_par)
    integer :: iopg(nopg_par)
    integer :: iopr(nopr_par)
    integer :: iopt(nopt_par)

    integer :: igerti(jet_par,nert_par)
    integer :: jgerti(nert_par)

    integer, dimension(:), allocatable :: iapxt
    integer, dimension(:), allocatable :: ibpxt
    integer, dimension(:), allocatable :: ibsrti
    integer, dimension(:), allocatable :: iesrti
    integer, dimension(:), allocatable :: iffg
    integer, dimension(:), allocatable :: iindx1
    integer, dimension(:), allocatable :: insgf
    integer, dimension(:), allocatable :: ipndx1
    integer, dimension(:), allocatable :: ixrti
    integer, dimension(:), allocatable :: jcode
    integer, dimension(:), allocatable :: jffg
    integer, dimension(:), allocatable :: jflag
    integer, dimension(:), allocatable :: jflagd
    integer, dimension(:), allocatable :: jflgi
    integer, dimension(:), allocatable :: jpflag
    integer, dimension(:), allocatable :: jreac
    integer, dimension(:), allocatable :: jsflag
    integer, dimension(:), allocatable :: jsitex
    integer, dimension(:), allocatable :: jsol
    integer, dimension(:), allocatable :: kxmod

    integer, dimension(:), allocatable :: narxt
    integer, dimension(:), allocatable :: nbasp
    integer, dimension(:), allocatable :: nbaspd
    integer, dimension(:), allocatable :: nbaspi
    integer, dimension(:), allocatable :: nbaspx
    integer, dimension(:), allocatable :: nbmap
    integer, dimension(:), allocatable :: ncmap
    integer, dimension(:), allocatable :: ndecsp
    integer, dimension(:), allocatable :: ndrs
    integer, dimension(:), allocatable :: ndrsd
    integer, dimension(:), allocatable :: ndrsx
    integer, dimension(:), allocatable :: ness
    integer, dimension(:), allocatable :: npchk
    integer, dimension(:), allocatable :: nphasx
    integer, dimension(:), allocatable :: nrndex
    integer, dimension(:), allocatable :: nsk
    integer, dimension(:), allocatable :: nsmap
    integer, dimension(:), allocatable :: nsts
    integer, dimension(:), allocatable :: ntf1
    integer, dimension(:), allocatable :: ntf2
    integer, dimension(:), allocatable :: nxridx

    integer, dimension(:,:), allocatable :: imech
    integer, dimension(:,:), allocatable :: ncmpr
    integer, dimension(:,:), allocatable :: ndrsr
    integer, dimension(:,:), allocatable :: ndrsrd
    integer, dimension(:,:), allocatable :: ndrsrx
    integer, dimension(:,:), allocatable :: nessr
    integer, dimension(:,:), allocatable :: nrk
    integer, dimension(:,:), allocatable :: nstsr

    integer, dimension(:,:,:), allocatable :: iact
    integer, dimension(:,:,:), allocatable :: ndact

    integer, dimension(:,:,:,:), allocatable :: ndac

    integer :: iern1
    integer :: iern2
    integer :: ifrn1
    integer :: ifrn2
    integer :: ilrn1
    integer :: ilrn2
    integer :: imrn1
    integer :: imrn2
    integer :: ixrn1
    integer :: ixrn2

    integer :: ifcphi1
    integer :: ifcphi2
    integer :: ifnnn
    integer :: ifn2n
    integer :: ifpsi1
    integer :: ifpsi2
    integer :: ifzeta
    integer :: ilcphi1
    integer :: ilcphi2
    integer :: ilnnn
    integer :: iln2n
    integer :: ilpsi1
    integer :: ilpsi2
    integer :: ilzeta

    integer :: nat
    integer :: nbt
    integer :: nct
    integer :: net
    integer :: ngt
    integer :: nlt
    integer :: nmt
    integer :: npt
    integer :: nst
    integer :: nxt

    integer :: narn1
    integer :: narn2
    integer :: nern1
    integer :: nern2
    integer :: nfrn1
    integer :: nfrn2
    integer :: ngrn1
    integer :: ngrn2
    integer :: nlrn1
    integer :: nlrn2
    integer :: nmrn1
    integer :: nmrn2
    integer :: nxrn1
    integer :: nxrn2

    integer :: kprs
    integer :: nbti
    integer :: nprpti
    integer :: nprsti

    integer :: i
    integer :: ie
    integer :: ier
    integer :: iaqsla
    integer :: iaqsln
    integer :: ielam
    integer :: iexec0
    integer :: igas
    integer :: ilevel
    integer :: ipch
    integer :: ipcv
    integer :: irang
    integer :: itermx
    integer :: ixrn1a
    integer :: ixrn2a
    integer :: izmax
    integer :: j
    integer :: je
    integer :: jexec0
    integer :: jlen
    integer :: jpdblo
    integer :: jpress
    integer :: jptffl
    integer :: jtemp
    integer :: j2
    integer :: j3
    integer :: k
    integer :: kbt
    integer :: kcol
    integer :: kct
    integer :: kdim
    integer :: kelect
    integer :: khydr
    integer :: khydx
    integer :: klim
    integer :: km1
    integer :: kmt
    integer :: kobswt
    integer :: ko2gaq
    integer :: krdxsp
    integer :: ksplmx
    integer :: ksppmx
    integer :: kstpmx
    integer :: kwater
    integer :: kx1
    integer :: kxt

    integer :: n
    integer :: napta
    integer :: narn1a
    integer :: narn2a
    integer :: nart
    integer :: nata
    integer :: nb
    integer :: nbta
    integer :: nbtafd
    integer :: nbtd
    integer :: nbw
    integer :: nbwa
    integer :: nb1
    integer :: nb2
    integer :: nchloa
    integer :: nchlor
    integer :: ncta
    integer :: ne
    integer :: neleca
    integer :: nelect
    integer :: ner
    integer :: nerr
    integer :: nert
    integer :: nffg
    integer :: ngrn1a
    integer :: ngrn2a
    integer :: ngrt
    integer :: ngta
    integer :: nhydr
    integer :: nhydra
    integer :: nhydx
    integer :: nhydxa
    integer :: nllnmx
    integer :: nlrn1a
    integer :: nlrn2a
    integer :: nlta
    integer :: nmax
    integer :: nmrn1a
    integer :: nmrn2a
    integer :: nmrt
    integer :: nmta
    integer :: nmuta
    integer :: nobswt
    integer :: no2gaa
    integer :: no2gaq
    integer :: np
    integer :: npi
    integer :: npslmx
    integer :: npta
    integer :: nprob
    integer :: nrc
    integer :: nrct
    integer :: nrdxsa
    integer :: nrdxsp
    integer :: nrecl
    integer :: nrr1
    integer :: nrr2
    integer :: nr1
    integer :: nr2
    integer :: nsbsw
    integer :: nsbswt
    integer :: ns
    integer :: nse
    integer :: nsi
    integer :: nslta
    integer :: nsrt
    integer :: nss
    integer :: nsslmx
    integer :: nsta
    integer :: ns1
    integer :: ns2
    integer :: nt
    integer :: ntf1t
    integer :: ntf1ta
    integer :: ntf2t
    integer :: ntf2ta
    integer :: ntitld
    integer :: ntitl1
    integer :: ntitl2
    integer :: ntpr
    integer :: ntprh
    integer :: ntprt
    integer :: ntrymx
    integer :: nxmod
    integer :: nxopex
    integer :: nxopt
    integer :: nxrn1a
    integer :: nxrn2a
    integer :: nxrt
    integer :: nxta

    integer :: ilnobl
    integer :: nbasis

    logical :: qbassw
    logical :: qbswok
    logical :: qchlor
    logical :: qcnpre
    logical :: qcntmp
    logical :: qcwrpj
    logical :: qdwipp
    logical :: qecon
    logical :: qelim
    logical :: qend
    logical :: qex
    logical :: qgexsh
    logical :: qhawep
    logical :: qloffg
    logical :: qop
    logical :: qoptmz
    logical :: qpit75
    logical :: qrderr
    logical :: qredox
    logical :: qscon
    logical :: qtatxt

    character(len=24) :: ugersi(iet_par,jet_par,nert_par)
    character(len=24) :: ugermo(nert_par)
    character(len=8) :: ugerji(jet_par,nert_par)

    character(len=80), dimension(:), allocatable :: utitl1
    character(len=80), dimension(:), allocatable :: utitl2
    character(len=48), dimension(:), allocatable :: ubasp
    character(len=48), dimension(:), allocatable :: ubmtbi
    character(len=48), dimension(:), allocatable :: uprspi
    character(len=48), dimension(:), allocatable :: uspec
    character(len=48), dimension(:), allocatable :: uxmod
    character(len=48), dimension(:), allocatable :: uzveci
    character(len=48), dimension(:), allocatable :: uzvec1
    character(len=24), dimension(:), allocatable :: uffg
    character(len=24), dimension(:), allocatable :: uphase
    character(len=24), dimension(:), allocatable :: uprphi
    character(len=24), dimension(:), allocatable :: uptype
    character(len=24), dimension(:), allocatable :: ureac
    character(len=24), dimension(:), allocatable :: uxcat
    character(len=24), dimension(:), allocatable :: uxopex
    character(len=8), dimension(:), allocatable :: uelem
    character(len=8), dimension(:), allocatable :: uxopt

    character(len=48), dimension(:,:), allocatable :: uobsw
    character(len=48), dimension(:,:), allocatable :: usbsw
    character(len=24), dimension(:,:), allocatable :: ubsri
    character(len=24), dimension(:,:), allocatable :: ucxri
    character(len=8), dimension(:,:), allocatable :: uesri

    character(len=24), dimension(:,:,:,:), allocatable :: udac

    character(len=nllnpa) :: ulinex
    character(len=80) :: ux80
    character(len=56) :: uspn56
    character(len=32) :: uactop
    character(len=24) :: unamsp
    character(len=24) :: unamph
    character(len=24) :: ux
    character(len=24) :: uaqsln
    character(len=24) :: ublk24
    character(len=11) :: utime0
    character(len=11) :: utime1
    character(len=9) :: udate0
    character(len=9) :: udate1
    character(len=8) :: udatfi
    character(len=8) :: ufixf
    character(len=8) :: udakey
    character(len=8) :: uinfor
    character(len=8) :: uplatc
    character(len=8) :: uplatm
    character(len=8) :: ustelg
    character(len=8) :: ustelu
    character(len=8) :: usteql
    character(len=8) :: usteq6
    character(len=8) :: uveelg
    character(len=8) :: uveelu
    character(len=8) :: uveeql
    character(len=8) :: uveeq6
    character(len=8) :: ux8

    real(kind=8) :: sscrew(nssc_par)
    real(kind=8) :: cco2(5)

    real(kind=8) :: egers(iet_par,jet_par,nert_par)
    real(kind=8) :: egersi(iet_par,jet_par,nert_par)
    real(kind=8) :: mrgers(iet_par,jet_par,nert_par)
    real(kind=8) :: xgers(iet_par,jet_par,nert_par)
    real(kind=8) :: xgersi(iet_par,jet_par,nert_par)

    real(kind=8) :: cegexs(iet_par,jet_par,net_par)
    real(kind=8) :: cpgexs(iet_par,jet_par,net_par)
    real(kind=8) :: egexjf(jet_par,net_par)
    real(kind=8) :: mrgexs(iet_par,jet_par,net_par)

    real(kind=8), dimension(:), allocatable :: atwt
    real(kind=8), dimension(:), allocatable :: azero
    real(kind=8), dimension(:), allocatable :: cdrs
    real(kind=8), dimension(:), allocatable :: cdrsd
    real(kind=8), dimension(:), allocatable :: cdrsx
    real(kind=8), dimension(:), allocatable :: cess
    real(kind=8), dimension(:), allocatable :: cscale
    real(kind=8), dimension(:), allocatable :: csts
    real(kind=8), dimension(:), allocatable :: elecsr
    real(kind=8), dimension(:), allocatable :: fkrc
    real(kind=8), dimension(:), allocatable :: loph
    real(kind=8), dimension(:), allocatable :: losp
    real(kind=8), dimension(:), allocatable :: modr
    real(kind=8), dimension(:), allocatable :: moffg
    real(kind=8), dimension(:), allocatable :: moph
    real(kind=8), dimension(:), allocatable :: morr
    real(kind=8), dimension(:), allocatable :: mosp
    real(kind=8), dimension(:), allocatable :: mprph
    real(kind=8), dimension(:), allocatable :: mprphi
    real(kind=8), dimension(:), allocatable :: mprsp
    real(kind=8), dimension(:), allocatable :: mprspi
    real(kind=8), dimension(:), allocatable :: mtb
    real(kind=8), dimension(:), allocatable :: mtbaq
    real(kind=8), dimension(:), allocatable :: mtbaqi
    real(kind=8), dimension(:), allocatable :: mtbi
    real(kind=8), dimension(:), allocatable :: mte
    real(kind=8), dimension(:), allocatable :: mteaq
    real(kind=8), dimension(:), allocatable :: mwtrc
    real(kind=8), dimension(:), allocatable :: mwtsp
    real(kind=8), dimension(:), allocatable :: ptk
    real(kind=8), dimension(:), allocatable :: sfcar
    real(kind=8), dimension(:), allocatable :: ssfcar
    real(kind=8), dimension(:), allocatable :: tempcu
    real(kind=8), dimension(:), allocatable :: tf1
    real(kind=8), dimension(:), allocatable :: tf2
    real(kind=8), dimension(:), allocatable :: ttk
    real(kind=8), dimension(:), allocatable :: xlkffg
    real(kind=8), dimension(:), allocatable :: xlkmod
    real(kind=8), dimension(:), allocatable :: vosp0
    real(kind=8), dimension(:), allocatable :: vreac
    real(kind=8), dimension(:), allocatable :: zchar
    real(kind=8), dimension(:), allocatable :: zchcu6
    real(kind=8), dimension(:), allocatable :: zchsq2
    real(kind=8), dimension(:), allocatable :: zvclgi
    real(kind=8), dimension(:), allocatable :: zvclg1
    real(kind=8), dimension(:), allocatable :: zvec1

    real(kind=8), dimension(:,:), allocatable :: cbsr
    real(kind=8), dimension(:,:), allocatable :: cesr
    real(kind=8), dimension(:,:), allocatable :: rxbar
    real(kind=8), dimension(:,:), allocatable :: cbsri
    real(kind=8), dimension(:,:), allocatable :: cesri
    real(kind=8), dimension(:,:), allocatable :: rxbari
    real(kind=8), dimension(:,:,:), allocatable :: csigma
    real(kind=8), dimension(:,:,:), allocatable :: eact
    real(kind=8), dimension(:,:,:), allocatable :: hact
    real(kind=8), dimension(:,:,:), allocatable :: rk
    real(kind=8), dimension(:,:,:), allocatable :: rkb
    real(kind=8), dimension(:,:,:), allocatable :: trkb
    real(kind=8), dimension(:,:,:,:), allocatable :: cdac

    real(kind=8), dimension(:), allocatable :: dadhh
    real(kind=8), dimension(:), allocatable :: dadhv
    real(kind=8), dimension(:), allocatable :: dbdhh
    real(kind=8), dimension(:), allocatable :: dbdhv
    real(kind=8), dimension(:), allocatable :: dbdth
    real(kind=8), dimension(:), allocatable :: dbdtv

    real(kind=8), dimension(:,:), allocatable :: apx
    real(kind=8), dimension(:,:), allocatable :: bpx
    real(kind=8), dimension(:,:), allocatable :: wfac

    real(kind=8), dimension(:,:), allocatable :: aadh
    real(kind=8), dimension(:,:), allocatable :: aadhh
    real(kind=8), dimension(:,:), allocatable :: aadhv
    real(kind=8), dimension(:,:), allocatable :: aaphi
    real(kind=8), dimension(:,:), allocatable :: abdh
    real(kind=8), dimension(:,:), allocatable :: abdhh
    real(kind=8), dimension(:,:), allocatable :: abdhv
    real(kind=8), dimension(:,:), allocatable :: abdot
    real(kind=8), dimension(:,:), allocatable :: abdoth
    real(kind=8), dimension(:,:), allocatable :: abdotv

    real(kind=8), dimension(:,:,:), allocatable :: adadhh
    real(kind=8), dimension(:,:,:), allocatable :: adadhv
    real(kind=8), dimension(:,:,:), allocatable :: adbdhh
    real(kind=8), dimension(:,:,:), allocatable :: adbdhv
    real(kind=8), dimension(:,:,:), allocatable :: adbdth
    real(kind=8), dimension(:,:,:), allocatable :: adbdtv

    real(kind=8) :: adh
    real(kind=8) :: adhh
    real(kind=8) :: adhv
    real(kind=8) :: aphi
    real(kind=8) :: bdh
    real(kind=8) :: bdhh
    real(kind=8) :: bdhv
    real(kind=8) :: bdot
    real(kind=8) :: bdoth
    real(kind=8) :: bdotv

    real(kind=8) :: afcnst
    real(kind=8) :: aftarg
    real(kind=8) :: al10
    real(kind=8) :: avcnst
    real(kind=8) :: awmax
    real(kind=8) :: awmaxi
    real(kind=8) :: awmin
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
    real(kind=8) :: dlxmax
    real(kind=8) :: dlxdmp
    real(kind=8) :: dlxmin
    real(kind=8) :: dlxmx0
    real(kind=8) :: dlxpll
    real(kind=8) :: dlxplo
    real(kind=8) :: dlxprl
    real(kind=8) :: dlxprn
    real(kind=8) :: ehfac
    real(kind=8) :: ehmax
    real(kind=8) :: ehmaxi
    real(kind=8) :: ehmin
    real(kind=8) :: ehmini
    real(kind=8) :: electr
    real(kind=8) :: eps100
    real(kind=8) :: farad
    real(kind=8) :: o2max
    real(kind=8) :: o2maxi
    real(kind=8) :: o2min
    real(kind=8) :: o2mini
    real(kind=8) :: phmax
    real(kind=8) :: phmaxi
    real(kind=8) :: phmin
    real(kind=8) :: phmini
    real(kind=8) :: prcinf
    real(kind=8) :: press
    real(kind=8) :: pressb
    real(kind=8) :: pressd
    real(kind=8) :: pressi
    real(kind=8) :: rconst
    real(kind=8) :: rcnstv
    real(kind=8) :: rtcnst
    real(kind=8) :: smp100
    real(kind=8) :: tcpu
    real(kind=8) :: tdamax
    real(kind=8) :: tdamin
    real(kind=8) :: tempc
    real(kind=8) :: tempcb
    real(kind=8) :: tempcd
    real(kind=8) :: tempci
    real(kind=8) :: tempk
    real(kind=8) :: texec0
    real(kind=8) :: timemx
    real(kind=8) :: timmxi
    real(kind=8) :: time1
    real(kind=8) :: tistrt
    real(kind=8) :: tistti
    real(kind=8) :: tolaft
    real(kind=8) :: tolbt
    real(kind=8) :: toldl
    real(kind=8) :: tolsat
    real(kind=8) :: tolsst
    real(kind=8) :: tolxsf
    real(kind=8) :: tolxst
    real(kind=8) :: tolxsu
    real(kind=8) :: trun
    real(kind=8) :: tuser
    real(kind=8) :: ximax
    real(kind=8) :: ximaxi
    real(kind=8) :: xistrt
    real(kind=8) :: xistti
    real(kind=8) :: xi1
    real(kind=8) :: x10
    real(kind=8) :: zkfac
    real(kind=8) :: zklgmn
    real(kind=8) :: zklogl
    real(kind=8) :: zklogu

    real(kind=8) :: av
    real(kind=8) :: azch
    real(kind=8) :: azchmx
    real(kind=8) :: dp
    real(kind=8) :: dt
    real(kind=8) :: ex
    real(kind=8) :: gx
    real(kind=8) :: xx
    real(kind=8) :: zx

    ! XX   real*8 vpgstp
    !      Variable declarations: Data file original contents. These
    !      variables and arrays remain unchanged after being filled by
    !      reading the data file. These data comprise the 'a' set.
    integer, dimension(:), allocatable :: iapxta
    integer, dimension(:), allocatable :: ibpxta
    integer, dimension(:), allocatable :: insgfa
    integer, dimension(:), allocatable :: jsola
    integer, dimension(:), allocatable :: nbaspa
    integer, dimension(:), allocatable :: ndrsa
    integer, dimension(:), allocatable :: nessa

    integer, dimension(:,:), allocatable :: ncmpra
    integer, dimension(:,:), allocatable :: ndrsra
    integer, dimension(:,:), allocatable :: nessra

    real(kind=8), dimension(:,:,:), allocatable :: axhfsa
    real(kind=8), dimension(:,:,:), allocatable :: axlksa
    real(kind=8), dimension(:,:,:), allocatable :: axvfsa

    real(kind=8), dimension(:,:,:,:), allocatable :: adhfsa
    real(kind=8), dimension(:,:,:,:), allocatable :: advfsa

    real(kind=8), dimension(:), allocatable :: atwta
    real(kind=8), dimension(:), allocatable :: azeroa
    real(kind=8), dimension(:), allocatable :: cdrsa
    real(kind=8), dimension(:), allocatable :: cessa
    real(kind=8), dimension(:), allocatable :: mwtspa
    real(kind=8), dimension(:), allocatable :: vosp0a
    real(kind=8), dimension(:), allocatable :: zchara

    real(kind=8), dimension(:,:), allocatable :: apxa
    real(kind=8), dimension(:,:), allocatable :: bpxa

    character(len=80), dimension(:), allocatable :: utitld
    character(len=48), dimension(:), allocatable :: uspeca
    character(len=24), dimension(:), allocatable :: uphasa
    character(len=24), dimension(:), allocatable :: uptypa
    character(len=8), dimension(:), allocatable :: uelema

    ! Pitzer's equations subset.
    integer, dimension(:), allocatable :: nalpaa
    integer, dimension(:,:), allocatable :: nmuxa
    integer, dimension(:,:), allocatable :: nslxa

    real(kind=8), dimension(:,:), allocatable :: amua
    real(kind=8), dimension(:,:), allocatable :: palpaa
    real(kind=8), dimension(:,:,:), allocatable :: aslma

    ! Variable declarations: Variables closely associated with the
    ! 'a' set.
    logical, dimension(:), allocatable :: qclnsa

    ! Variable declarations: Pseudo-data file original contents.
    integer, dimension(:), allocatable :: ntf1a
    integer, dimension(:), allocatable :: ntf2a

    real(kind=8), dimension(:), allocatable :: tf1a
    real(kind=8), dimension(:), allocatable :: tf2a

    ! BEGIN_MACHINE_DEPENDENT_CODE
    !   On some systems, a BLOCK DATA subroutine must be declared in an
    !   EXTERNAL statement to assure proper loading. On some other
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
    external bkdeq6

    ! END_MACHINE_DEPENDENT_CODE
    ! Special data base of coefficients for the HCO3-CO3-OH total
    ! alkalinity.
    include 'eqlib/eqltf1.h'

    ! Special data base of coefficients for the extended total
    ! alkalinity.
    include 'eqlib/eqltf2.h'

    ! Other data statements.
    data uaqsln /'Aqueous solution        '/
    data ublk24 /'                        '/
    data ufixf  /'fix_f   '/

    ! Set practical infinity.
    data prcinf/1.e+38/

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
    nad1 = 0
    nbkupa = 0
    nbkupb = 0
    ninpt = 0
    ninpts = 0
    newin = 0
    ntab = 0
    ntabs = 0
    ntabx = 0
    nrecl = 0

    ! Open all files execpt pickup, tab, tabs, and tabx.
    numargs = COMMAND_ARGUMENT_COUNT()

    if (numargs.eq.2) then
        call GET_COMMAND_ARGUMENT(1,temppath)
        data1path = TRIM(temppath)
        temppath(:)=''
        call GET_COMMAND_ARGUMENT(2,temppath)
        sixipath = TRIM(temppath)
    else
        write (0, *) 'usage: eq6 <data1> <6i>'
        stop
    end if

    call getbasename(sixipath, pathindices)
    basename = sixipath(pathindices(1):pathindices(2))

    ofile = basename // '.6o'
    pfile = basename // '.6p'
    bafile = basename // '.6ba'
    bbfile = basename // '.6bb'
    ifile = basename // '.6ib'
    tfile = basename // '.6t'
    txfile = basename // '.6tx'
    tsfile = basename // '.6ts'

    call openin(noutpt,nttyo,data1path,'unformatted',nad1)
    call openin(noutpt,nttyo,sixipath,'formatted',ninpt)

    call openou(noutpt,nttyo,ofile,'formatted',nrecl,noutpt)
    call openou(noutpt,nttyo,ifile,'formatted',nrecl,ninpts)

    ! Make a copy of the input file, stripped of comments.
    call stripl(ninpt,ninpts)
    close (ninpt)
    ninpt = 0
    rewind ninpts

    ! Get configuration identification data.
    call aaaeq6(usteq6,uveeq6)
    call aaaeql(usteql,uveeql)
    call aaaelg(ustelg,uveelg)
    call aaaelu(ustelu,uveelu)
    call platfd(uplatc,uplatm)

    ! Write configuration identification data, the copyright statement,
    ! and any remaining statements or disclaimers.
    i = index(uveeq6,' ') - 1
    j = index(uveeq6,'.') - 1
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

    write (nttyo,1000) uveeq6(1:i),uveeq6(1:j),uveeq6(1:i),uplatc(1:k)
    write (noutpt,1000) uveeq6(1:i),uveeq6(1:j),uveeq6(1:i),uplatc(1:k)
1000 format(/' EQ3/6, Version ',a,' (EQ3/6-V',a,'-REL-V',a,'-',a,')')

    i = j
    j = index(usteq6,' ') - 1
    k = index(uplatm,' ') - 1

    if (j .le. 0) then
        j = 8
    end if

    if (k .le. 0) then
        k = 8
    end if

    write (nttyo,1010) uveeq6(1:i),usteq6(1:j),uplatm(1:k)
    write (noutpt,1010) uveeq6(1:i),usteq6(1:j),uplatm(1:k)
1010 format(' EQ6 Reaction-Path Code (EQ/36-V',a,'-EQ6-EXE-',a,'-',a,')')

    write (noutpt,1020)
    write (nttyo,1020)
1020 format(' Supported by the following EQ3/6 libraries:')

    i = index(uveeql,'.') - 1
    j = index(usteql,' ') - 1

    if (i .le. 0) then
        i = 8
    end if

    if (j .le. 0) then
        j = 8
    end if

    write (nttyo,1030) uveeql(1:i),usteql(1:j),uplatm(1:k)
    write (noutpt,1030) uveeql(1:i),usteql(1:j),uplatm(1:k)
1030 format('   EQLIB (EQ3/6-V',a,'-EQLIB-LIB-',a,'-',a,')')

    i = index(uveelg,'.') - 1
    j = index(ustelg,' ') - 1

    if (i .le. 0) then
        i = 8
    end if

    if (j .le. 0) then
        j = 8
    end if

    write (nttyo,1040) uveelg(1:i),ustelg(1:j),uplatm(1:k)
    write (noutpt,1040) uveelg(1:i),ustelg(1:j),uplatm(1:k)
1040 format('   EQLIBG (EQ3/6-V',a,'-EQLIBG-LIB-',a,'-',a,')')

    i = index(uveelu,'.') - 1
    j = index(ustelu,' ') - 1

    if (i .le. 0) then
        i = 8
    end if

    if (j .le. 0) then
        j = 8
    end if

    write (nttyo,1050) uveelu(1:i),ustelu(1:j),uplatm(1:k)
    write (noutpt,1050) uveelu(1:i),ustelu(1:j),uplatm(1:k)
1050 format('   EQLIBU (EQ3/6-V',a,'-EQLIBU-LIB-',a,'-',a,')',/)

    write (nttyo,1060)
    write (noutpt,1060)
1060 format(' Copyright (c) 1987, 1990-1993, 1995, 1997, 2002 The',' Regents of the',/' University of California, Lawrence',' Livermore National Laboratory.',/' All rights reserved.',/)

    ! Write additional statements and disclaimers.
    call prcndi(noutpt,nttyo)

    ! Write the time and date on the output.
    j2 = ilnobl(udate0)
    write(noutpt,1080) utime0,udate0(1:j2)
    write(nttyo,1080) utime0,udate0(1:j2)
1080 format(' Run',2(2x,a8,2x,a),//)

    ! Get the platform's real*8 floating-point parameters.
    call flpars(eps100,irang,noutpt,nttyo,smp100)

    ! Initialize array dimension variables.
    ! The following are common with EQ3NR.
    ietmax = iet_par
    jetmax = jet_par
    ketmax = ket_par
    netmax = net_par

    jsomax = jso_par
    nodbmx = nodb_par
    nopgmx = nopg_par
    noprmx = nopr_par
    noptmx = nopt_par
    nvetmx = nvet_par

    nxmd_asv = nxmd_par
    nxmdmx = nxmd_asv

    ntit_asv = ntit_par
    ntitmx = ntit_asv

    iapxa_asv = iapxa_par
    ibpxa_asv = ibpxa_par

    ! The following are common with EQ3NR, but are only needed by that
    ! code to write a full EQ6 input file as its pickup file.
    nrct_asv = nrct_par
    nert_asv = nert_par
    nsrt_asv = nsrt_par
    nxrt_asv = nxrt_par
    imch_asv = imch_par
    ndct_asv = ndct_par
    nffg_asv = nffg_par
    nprp_asv = nprp_par
    nprs_asv = nprs_par
    nptk_asv = nptk_par
    nttk_asv = nttk_par
    nxop_asv = nxop_par
    nxpe_asv = nxpe_par

    nrctmx = nrct_asv
    nertmx = nert_asv
    nsrtmx = nsrt_asv
    nxrtmx = nxrt_asv
    imchmx = imch_asv
    ndctmx = ndct_asv
    nffgmx = nffg_asv
    nprpmx = nprp_asv
    nprsmx = nprs_asv
    nptkmx = nptk_asv
    nttkmx = nttk_asv
    nxopmx = nxop_asv
    nxpemx = nxpe_asv

    ! The following are unique to EQ6.
    nsscmx = nssc_par
    nllnmx = nllnpa

    ! Initialize the following floating-point constants:
    !   rconst = the gas constant: 1.98726 cal/mole-K
    !   rcnstv = the gas constant: 83.14510 bar-cm3/mol-K
    !   vpgstp  = the volume of a perfect gas at STP: 22413.6 cm3
    !   farad  = the Faraday constant: 23062.3 cal/equiv-volt
    !   al10   = ln 10; =~ 2.3026
    rconst = 1.98726
    rcnstv = 83.14510

    ! XX   vpgstp = 22413.6
    farad = 23062.3
    x10 = 10.
    al10 = log(x10)

    ! Read the header section of the DATA1 file. This section consists
    ! of a record containing the string 'data1' (to ensure that the file
    ! is indeed a DATA1 file), a record containing the keystring for the
    ! activity coefficient model for aqueous species to which this data
    ! file corresponds, and the array allocation size variables required
    ! to allocate sufficient array space to store the rest of the data
    ! on this data file.
    call indath(ikta_asv,ipbt_asv,ipch_asv,ipcv_asv,jpfc_asv,nad1,napa_asv,narx_asv,nata_asv,nbta_asv,ncta_asv,ngta_asv,nlta_asv,nmta_asv,npta_asv,nmuta_asv,noutpt,nslta_asv,nsta_asv,ntid_asv,ntpr_asv,nttyo,nxta_asv,udakey)

    ! The value of nbta_asv at this point matches the number of
    ! formally declared basis species on the data file. That is
    ! the number of species in the "strict" and "auxiliary" basis
    ! set sections of the aqueous species superblock. Save this
    ! as nbtafd. Other species can be added to the basis set in
    ! this version of EQ3/6. Thus, nbta_asv at this point is not
    ! the "final answer."
    nbtafd = nbta_asv

    ! Increase the values of some allocation size variables to allow
    ! for some species, phases, etc., to created by this software.
    ! Examples: generic ion exchanger phases and species.
    nbta_asv = nbta_asv + jet_par*net_par + 10
    npta_asv = npta_asv + net_par + 10
    nsta_asv = nsta_asv + iet_par*jet_par*net_par + 10

    ! Set the values of some secondary allocation size variables.
    nbta1_asv = nbta_asv + 1
    ndrsa_asv = 7*nsta_asv
    nessa_asv = 5*nsta_asv

    ! Allocate arrays to store the data read from the rest of the DATA1
    ! file.
    ALLOCATE(iapxta(nxta_asv))
    ALLOCATE(ibpxta(nxta_asv))
    ALLOCATE(insgfa(nata_asv))
    ALLOCATE(jsola(nxta_asv))
    ALLOCATE(narxt(ntpr_asv))
    ALLOCATE(nbaspa(nbta_asv))
    ALLOCATE(ndrsa(ndrsa_asv))
    ALLOCATE(nessa(nessa_asv))

    ALLOCATE(ncmpra(2,npta_asv))
    ALLOCATE(ndrsra(2,nsta_asv))
    ALLOCATE(nessra(2,nsta_asv))

    ALLOCATE(qclnsa(nsta_asv))

    ALLOCATE(utitld(ntid_asv))
    ALLOCATE(uspeca(nsta_asv))
    ALLOCATE(uphasa(npta_asv))
    ALLOCATE(uptypa(npta_asv))
    ALLOCATE(uelema(ncta_asv))
    ALLOCATE(ubasp(nbta_asv))

    ALLOCATE(tempcu(ntpr_asv))

    ALLOCATE(aprehw(narx_asv,ntpr_asv))
    ALLOCATE(apresg(narx_asv,ntpr_asv))

    ALLOCATE(aadh(narx_asv,ntpr_asv))
    ALLOCATE(aadhh(narx_asv,ntpr_asv))
    ALLOCATE(aadhv(narx_asv,ntpr_asv))
    ALLOCATE(aaphi(narx_asv,ntpr_asv))
    ALLOCATE(abdh(narx_asv,ntpr_asv))
    ALLOCATE(abdhh(narx_asv,ntpr_asv))
    ALLOCATE(abdhv(narx_asv,ntpr_asv))
    ALLOCATE(abdot(narx_asv,ntpr_asv))
    ALLOCATE(abdoth(narx_asv,ntpr_asv))
    ALLOCATE(abdotv(narx_asv,ntpr_asv))

    ALLOCATE(axhfe(narx_asv,ntpr_asv))
    ALLOCATE(axlke(narx_asv,ntpr_asv))
    ALLOCATE(axvfe(narx_asv,ntpr_asv))

    ALLOCATE(adadhh(narx_asv,ntpr_asv,ipch_asv))
    ALLOCATE(adadhv(narx_asv,ntpr_asv,ipcv_asv))
    ALLOCATE(adbdhh(narx_asv,ntpr_asv,ipch_asv))
    ALLOCATE(adbdhv(narx_asv,ntpr_asv,ipcv_asv))
    ALLOCATE(adbdth(narx_asv,ntpr_asv,ipch_asv))
    ALLOCATE(adbdtv(narx_asv,ntpr_asv,ipcv_asv))
    ALLOCATE(adhfe(narx_asv,ntpr_asv,ipch_asv))
    ALLOCATE(advfe(narx_asv,ntpr_asv,ipcv_asv))

    ALLOCATE(axhfsa(narx_asv,ntpr_asv,nsta_asv))
    ALLOCATE(axlksa(narx_asv,ntpr_asv,nsta_asv))
    ALLOCATE(axvfsa(narx_asv,ntpr_asv,nsta_asv))

    ALLOCATE(adhfsa(narx_asv,ntpr_asv,ipch_asv,nsta_asv))
    ALLOCATE(advfsa(narx_asv,ntpr_asv,ipcv_asv,nsta_asv))

    ALLOCATE(atwta(ncta_asv))
    ALLOCATE(cdrsa(ndrsa_asv))
    ALLOCATE(cessa(nessa_asv))
    ALLOCATE(mwtspa(nsta_asv))
    ALLOCATE(vosp0a(nsta_asv))
    ALLOCATE(zchara(nsta_asv))

    ALLOCATE(apxa(iapxa_asv,nxta_asv))
    ALLOCATE(bpxa(ibpxa_asv,nxta_asv))

    ! SEDH data.
    ALLOCATE(azeroa(nata_asv))

    ! Pitzer data.
    ALLOCATE(nalpaa(nslta_asv))
    ALLOCATE(nmuxa(3,nmuta_asv))
    ALLOCATE(nslxa(2,nslta_asv))

    ALLOCATE(amua(jpfc_asv,nmuta_asv))

    ALLOCATE(aslma(jpfc_asv,0:ipbt_asv,nslta_asv))
    ALLOCATE(palpaa(ipbt_asv,napa_asv))

    ! Read the remainder of the DATA1 file. The image of the data file
    ! is stored in arrays and variables that typically end in 'a.'
    ! A subset of this is the reaction data. The variable nbta and
    ! the arrays nbaspa, cdrsa, ndrsa, ndrsra, and axlksa comprise
    ! what is called the 'a' set of reaction data. These are not
    ! modified by any code manipulations.
    call indata(aadh,aadhh,aadhv,aaphi,abdh,abdhh,abdhv,abdot,abdoth,abdotv,adadhh,adadhv,adbdhh,adbdhv,adbdth,adbdtv,adhfe,adhfsa,advfe,advfsa,amua,aprehw,apresg,apxa,aslma,atwta,axhfe,axhfsa,axlke,axlksa,axvfe,axvfsa,azeroa,bpxa,cco2,cdrsa,cessa,eps100,iapxa_asv,iapxta,iaqsla,ibpxa_asv,ibpxta,ielam,igas,ikta_asv,insgfa,ipbt_asv,ipch,ipch_asv,ipcv,ipcv_asv,ixrn1a,ixrn2a,jpdblo,jpfc_asv,jptffl,jsola,mwtspa,nad1,nalpaa,napa_asv,napta,narn1a,narn2a,narx_asv,narxt,nata,nata_asv,nbaspa,nbta,nbta_asv,nbtafd,nbta1_asv,ncmpra,ncta,ncta_asv,ndrsa,ndrsa_asv,ndrsra,nessa,nessa_asv,nessra,ngrn1a,ngrn2a,ngta,ngta_asv,nlrn1a,nlrn2a,nlta,nlta_asv,nmrn1a,nmrn2a,nmta,nmta_asv,nmuta,nmuta_asv,nmuxa,noutpt,npta,npta_asv,nslta,nslta_asv,nslxa,nsta,nsta_asv,ntid_asv,ntitld,ntpr_asv,ntprt,nttyo,nxrn1a,nxrn2a,nxta,nxta_asv,qclnsa,palpaa,tdamax,tdamin,tempcu,ubasp,udakey,udatfi,uelema,uspeca,uphasa,uptypa,utitld,vosp0a,zchara)

    ! Allocate arrays to store the data extracted from the pseudo-data
    ! files.
    ntf1a_asv = ntf1a_par
    ALLOCATE(ntf1a(ntf1a_asv))
    ALLOCATE(tf1a(ntf1a_asv))

    ntf2a_asv = ntf2a_par
    ALLOCATE(ntf2a(ntf2a_asv))
    ALLOCATE(tf2a(ntf2a_asv))

    ntf1mx = ntf1a_asv
    ntf2mx = ntf2a_asv
    ntfxmx = ntf1mx

    ! Get the coefficients for calculating the HCO3-CO3-OH total
    ! alkalinity. These are actually set in a pseudo-data file in the
    ! EQLIB INCLUDE file eqltf1.h. The subroutine called here assigns
    ! the coefficients to the species read from the data file.
    ! Calling sequence substitutions:
    !   ntf1a for ntfxa
    !   ntf1mx for ntfxmx
    !   ntf1ta for ntfxta
    !   tf1a for tfxa
    !   utf1xd for utfxxd
    call inttfx(narn1a,narn2a,noutpt,nsta_asv,ntf1a,ntf1mx,ntf1ta,nttyo,tf1a,uspeca,utf1xd)

    ! Get the coefficients for calculating the extended total
    ! alkalinity. These are actually set in a pseudo-data file in the
    ! EQLIB INCLUDE file eqltf2.h. The subroutine called here assigns
    ! the coefficients to the species read from the data file.
    ! Calling sequence substitutions:
    !   ntf2a for ntfxa
    !   ntf2mx for ntfxmx
    !   ntf2ta for ntfxta
    !   tf2a for tfxa
    !   utf2xd for utfxxd
    call inttfx(narn1a,narn2a,noutpt,nsta_asv,ntf2a,ntf2mx,ntf2ta,nttyo,tf2a,uspeca,utf2xd)

    ! Get the indices of H+, OH-, Cl-, fictive aqueous O2(g), and
    ! fictive aqueous e-. This will be repeated after data compression.
    ! Calling sequence substitutions:
    !   narn1a for narn1
    !   narn2a for narn2
    !   nchloa for nchlor
    !   neleca for nelect
    !   nhydra for nhydr
    !   nhydxa for nhydx
    !   no2gaa for no2gaq
    !   nsta_asv for nstmax
    !   uspeca for uspec
    call gspion(narn1a,narn2a,nchloa,neleca,nhydra,nhydxa,noutpt,no2gaa,nsta_asv,nttyo,uspeca)

    ! Calling sequence substitutions:
    !   nbaspa for nbasp
    !   nbta_asv for nbtmax
    !   nbta for nbt
    !   ncta for nct
    !   ndrsra for ndrsr
    !   nrdxsa for nrdxsp
    !   nsta_asv for nstmax
    !   uspeca for uspec
    call grdxsp(nbaspa,nbta,nbta_asv,ncta,ndrsra,noutpt,nrdxsa,nsta_asv,nttyo,uspeca)

    ! Get the basis index of water.
    ! Calling sequence substitutions:
    !   nbaspa for nbasp
    !   nbta for nbt
    !   nbta_asv for nbtmax
    !   narn1a for ns
    nbwa = nbasis(nbaspa,nbta,nbta_asv,narn1a)

    ! Setup up the remaining array allocation variables required to
    ! allocate variables needed to store the dta to be read from the
    ! input file.
    ! The following value of k_asv is special to EQ6.
    ! A smaller value is adequate for EQ3NR.
    k_asv = 2*nbta_asv + 5

    ! Allocate arrays to store the data to be read from the input
    ! file.
    ALLOCATE(jflgi(nbta_asv))
    ALLOCATE(nbaspi(nbta_asv))

    ALLOCATE(mtbaqi(nbta_asv))
    ALLOCATE(mtbi(nbta_asv))
    ALLOCATE(zvclgi(k_asv))

    ! Setup up the remaining required array allocation variables.
    ! At the present time, just set the regular variables which
    ! correspond to specific variables required for the 'a' set
    ! of data (e.g., nbt_asv to nbta_asv) to the values of the
    ! corresponding 'a' set variables.
    ! Primary allocation variables.
    nat_asv = nata_asv
    nbt_asv = nbta_asv
    nct_asv = ncta_asv
    ngt_asv = ngta_asv
    nlt_asv = nlta_asv
    nmt_asv = nmta_asv
    npt_asv = npta_asv
    nst_asv = nsta_asv
    nxt_asv = nxta_asv

    nap_asv = napa_asv
    nslt_asv = nslta_asv
    nmut_asv = nmuta_asv

    ikt_asv = ikta_asv
    iapx_asv = iapxa_asv
    ibpx_asv = ibpxa_asv

    ntf1_asv = ntf1a_asv
    ntf2_asv = ntf2a_asv

    ! Secondary allocation variables.
    nbt1_asv = nbt_asv + 1
    ndrs_asv = 7*nst_asv
    ness_asv = 5*nst_asv
    nmx_asv = 3*nmut_asv
    nsx_asv = 2*nslt_asv

    npet_asv = nbt_asv + 5
    nset_asv = 20*npet_asv

    ! Note: nazp_asv and nazm_asv could be better calculated after
    ! compression of the set of aqueous species.
    azchmx = 0.

    do ns = narn1a,narn2a
        azch = abs(zchara(ns))
        azchmx = max(azch,azchmx)
    end do

    nazp_asv = nint(azchmx)
    nazm_asv = -nazp_asv

    ! Initialize array dimension variables corresponding to allocation
    ! size variables associated with run-time dimensioning.
    iapxmx = iapx_asv
    ibpxmx = ibpx_asv
    iktmax = ikt_asv
    ipchmx = ipch_asv
    ipcvmx = ipcv_asv
    ipbtmx = ipbt_asv
    jpfcmx = jpfc_asv

    kmax = k_asv

    napmax = nap_asv
    narxmx = narx_asv
    natmax = nat_asv
    nazpmx = nazp_asv
    nazmmx = nazm_asv
    nbtmax = nbt_asv
    nbt1mx = nbt1_asv
    nctmax = nct_asv
    ndrsmx = ndrs_asv
    nessmx = ness_asv
    ngtmax = ngt_asv
    nltmax = nlt_asv
    nmtmax = nmt_asv
    nmutmx = nmut_asv
    npetmx = npet_asv
    nptmax = npt_asv
    nsetmx = nset_asv
    nsltmx = nslt_asv
    nstmax = nst_asv
    ntidmx = ntid_asv
    ntprmx = ntpr_asv
    nxtmax = nxt_asv

    nmxmax = nmx_asv
    nsxmax = nsx_asv

    ! Allocate the rest of the necessary arrays.
    ALLOCATE(iindx1(k_asv))
    ALLOCATE(ipndx1(k_asv))

    ALLOCATE(insgf(nat_asv))

    ALLOCATE(nbasp(nbt_asv))
    ALLOCATE(nbaspd(nbt_asv))
    ALLOCATE(nbaspx(nbt_asv))
    ALLOCATE(nbmap(nbt_asv))
    ALLOCATE(ndecsp(nbt_asv))

    ALLOCATE(ncmap(nct_asv))

    ALLOCATE(ndrs(ndrs_asv))
    ALLOCATE(ndrsd(ndrs_asv))
    ALLOCATE(ndrsx(ndrs_asv))

    ALLOCATE(ness(ness_asv))

    ALLOCATE(jpflag(npt_asv))
    ALLOCATE(ncmpr(2,npt_asv))

    ALLOCATE(jflag(nst_asv))
    ALLOCATE(jflagd(nst_asv))
    ALLOCATE(jsflag(nst_asv))
    ALLOCATE(jsitex(nst_asv))
    ALLOCATE(nphasx(nst_asv))
    ALLOCATE(nsmap(nst_asv))

    ALLOCATE(ntf1(ntf1_asv))
    ALLOCATE(ntf2(ntf2_asv))

    ALLOCATE(kxmod(nxmd_asv))

    ALLOCATE(iapxt(nxt_asv))
    ALLOCATE(ibpxt(nxt_asv))
    ALLOCATE(jsol(nxt_asv))

    ALLOCATE(ndrsr(2,nst_asv))
    ALLOCATE(ndrsrd(2,nst_asv))
    ALLOCATE(ndrsrx(2,nst_asv))
    ALLOCATE(nessr(2,nst_asv))

    ALLOCATE(utitl1(ntit_asv))
    ALLOCATE(utitl2(ntit_asv))

    ALLOCATE(ubmtbi(nbt_asv))

    ALLOCATE(uobsw(2,nbt_asv))
    ALLOCATE(usbsw(2,nbt_asv))

    ALLOCATE(uspec(nst_asv))
    ALLOCATE(uxmod(nxmd_asv))
    ALLOCATE(uzveci(k_asv))
    ALLOCATE(uzvec1(k_asv))

    ALLOCATE(uphase(npt_asv))
    ALLOCATE(uptype(npt_asv))

    ALLOCATE(uelem(nct_asv))

    ALLOCATE(apx(iapx_asv,nxt_asv))
    ALLOCATE(bpx(ibpx_asv,nxt_asv))
    ALLOCATE(wfac(ikt_asv,nxt_asv))

    ALLOCATE(loph(npt_asv))
    ALLOCATE(moph(npt_asv))

    ALLOCATE(azero(nat_asv))
    ALLOCATE(atwt(nct_asv))

    ALLOCATE(zvclg1(k_asv))
    ALLOCATE(zvec1(k_asv))

    ALLOCATE(mtb(nbt_asv))
    ALLOCATE(mtbaq(nbt_asv))

    ALLOCATE(mte(nct_asv))
    ALLOCATE(mteaq(nct_asv))

    ALLOCATE(cdrs(ndrs_asv))
    ALLOCATE(cdrsd(ndrs_asv))
    ALLOCATE(cdrsx(ndrs_asv))
    ALLOCATE(losp(nst_asv))
    ALLOCATE(mosp(nst_asv))
    ALLOCATE(mwtsp(nst_asv))
    ALLOCATE(vosp0(nst_asv))
    ALLOCATE(zchar(nst_asv))
    ALLOCATE(zchcu6(nst_asv))
    ALLOCATE(zchsq2(nst_asv))

    ALLOCATE(cess(ness_asv))
    ALLOCATE(tf1(ntf1_asv))
    ALLOCATE(tf2(ntf2_asv))
    ALLOCATE(xlkmod(nxmd_asv))

    ALLOCATE(dadhh(ipch_asv))
    ALLOCATE(dadhv(ipcv_asv))
    ALLOCATE(dbdhh(ipch_asv))
    ALLOCATE(dbdhv(ipcv_asv))
    ALLOCATE(dbdth(ipch_asv))
    ALLOCATE(dbdtv(ipcv_asv))

    ALLOCATE(dhfe(ipch_asv))
    ALLOCATE(dvfe(ipcv_asv))

    ! Allocate standard-state data arrays for module mod6xf.
    ! The coefficient arrays here correspond to the 'a' set
    ! arrays read from the data file. The present coefficient
    ! arrays (and the other arrays) correspond to a compressed
    ! set of species.
    ALLOCATE(xhfs(nst_asv))
    ALLOCATE(xhfsd(nst_asv))
    ALLOCATE(xlks(nst_asv))
    ALLOCATE(xlksd(nst_asv))
    ALLOCATE(xvfs(nst_asv))
    ALLOCATE(xvfsd(nst_asv))

    ALLOCATE(dhfs(ipch_asv,nst_asv))
    ALLOCATE(dhfsd(ipch_asv,nst_asv))
    ALLOCATE(dvfs(ipcv_asv,nst_asv))
    ALLOCATE(dvfsd(ipcv_asv,nst_asv))

    ALLOCATE(axhfs(narx_asv,ntpr_asv,nst_asv))
    ALLOCATE(axhfsd(narx_asv,ntpr_asv,nst_asv))
    ALLOCATE(axhfsx(narx_asv,ntpr_asv,nst_asv))
    ALLOCATE(axlks(narx_asv,ntpr_asv,nst_asv))
    ALLOCATE(axlksd(narx_asv,ntpr_asv,nst_asv))
    ALLOCATE(axlksx(narx_asv,ntpr_asv,nst_asv))
    ALLOCATE(axvfs(narx_asv,ntpr_asv,nst_asv))
    ALLOCATE(axvfsd(narx_asv,ntpr_asv,nst_asv))
    ALLOCATE(axvfsx(narx_asv,ntpr_asv,nst_asv))

    ALLOCATE(adhfs(narx_asv,ntpr_asv,ipch_asv,nst_asv))
    ALLOCATE(adhfsd(narx_asv,ntpr_asv,ipch_asv,nst_asv))
    ALLOCATE(adhfsx(narx_asv,ntpr_asv,ipch_asv,nst_asv))
    ALLOCATE(advfs(narx_asv,ntpr_asv,ipcv_asv,nst_asv))
    ALLOCATE(advfsd(narx_asv,ntpr_asv,ipcv_asv,nst_asv))
    ALLOCATE(advfsx(narx_asv,ntpr_asv,ipcv_asv,nst_asv))

    ! Allocate Pitzer data arrays for module mod6pt.
    ALLOCATE(nalpha(nslt_asv))

    ALLOCATE(nmux(3,nmut_asv))
    ALLOCATE(nmxi(2,nat_asv))
    ALLOCATE(nmxx(3,nmx_asv))
    ALLOCATE(nslx(2,nslt_asv))
    ALLOCATE(nsxi(2,nat_asv))
    ALLOCATE(nsxx(2,nsx_asv))

    ALLOCATE(amu(jpfc_asv,nmut_asv))
    ALLOCATE(pmu(nmut_asv))

    ALLOCATE(aslm(jpfc_asv,0:ipbtmx,nslt_asv))
    ALLOCATE(pslm(nslt_asv))
    ALLOCATE(dpslm(2,nslt_asv))
    ALLOCATE(pslamn(0:ipbt_asv,nslt_asv))

    ALLOCATE(gpit(ipbt_asv,nap_asv))
    ALLOCATE(dgpit(2,ipbt_asv,nap_asv))
    ALLOCATE(palpha(ipbt_asv,nap_asv))

    ALLOCATE(elam(nazp_asv,nazp_asv))
    ALLOCATE(delam(2,nazp_asv,nazp_asv))
    ALLOCATE(pelm(nazp_asv,nazp_asv))
    ALLOCATE(dpelm(2,nazp_asv,nazp_asv))
    ALLOCATE(selm(nazm_asv:nazp_asv))
    ALLOCATE(dselm(2,nazm_asv:nazp_asv))

    ! Allocate still more arrays.
    ALLOCATE(iffg(nffg_asv))
    ALLOCATE(jffg(nffg_asv))
    ALLOCATE(npchk(npt_asv))

    ALLOCATE(cscale(nst_asv))

    ALLOCATE(jcode(nrct_asv))
    ALLOCATE(jreac(nrct_asv))
    ALLOCATE(nsk(nrct_asv))
    ALLOCATE(nrndex(nrct_asv))
    ALLOCATE(nxridx(nrct_asv))

    ALLOCATE(imech(2,nrct_asv))
    ALLOCATE(nrk(2,nrct_asv))

    ALLOCATE(iact(imch_asv,2,nrct_asv))
    ALLOCATE(ndact(imch_asv,2,nrct_asv))

    ALLOCATE(ibsrti(nsrt_asv))
    ALLOCATE(iesrti(nsrt_asv))

    ALLOCATE(ixrti(nxrt_asv))
    ALLOCATE(ndac(ndct_asv,imch_asv,2,nrct_asv))

    ALLOCATE(uprspi(nprs_asv))

    ALLOCATE(uprphi(nprp_asv))
    ALLOCATE(ubsri(nbt1_asv,nsrt_asv))
    ALLOCATE(ucxri(ikt_asv,nxrt_asv))
    ALLOCATE(udac(ndct_asv,imch_asv,2,nrct_asv))
    ALLOCATE(uffg(nffg_asv))
    ALLOCATE(ureac(nrct_asv))
    ALLOCATE(uxcat(nxop_asv))
    ALLOCATE(uxopt(nxop_asv))
    ALLOCATE(uxopex(nxpe_asv))
    ALLOCATE(uesri(nct_asv,nsrt_asv))

    ALLOCATE(cbsr(nbt1_asv,nsrt_asv))
    ALLOCATE(cbsri(nbt1_asv,nsrt_asv))
    ALLOCATE(cdac(ndct_asv,imch_asv,2,nrct_asv))
    ALLOCATE(cesr(nct_asv,nsrt_asv))
    ALLOCATE(cesri(nct_asv,nsrt_asv))
    ALLOCATE(csigma(imch_asv,2,nrct_asv))
    ALLOCATE(eact(imch_asv,2,nrct_asv))
    ALLOCATE(elecsr(nsrt_asv))
    ALLOCATE(hact(imch_asv,2,nrct_asv))
    ALLOCATE(rk(imch_asv,2,nrct_asv))
    ALLOCATE(rkb(imch_asv,2,nrct_asv))
    ALLOCATE(trkb(imch_asv,2,nrct_asv))
    ALLOCATE(fkrc(nrct_asv))
    ALLOCATE(modr(nrct_asv))
    ALLOCATE(morr(nrct_asv))
    ALLOCATE(mwtrc(nrct_asv))
    ALLOCATE(sfcar(nrct_asv))
    ALLOCATE(ssfcar(nrct_asv))
    ALLOCATE(vreac(nrct_asv))
    ALLOCATE(moffg(nffg_asv))
    ALLOCATE(xlkffg(nffg_asv))
    ALLOCATE(mprph(npt_asv))
    ALLOCATE(mprphi(nprp_asv))
    ALLOCATE(mprsp(nst_asv))
    ALLOCATE(mprspi(nprs_asv))
    ALLOCATE(ptk(nptk_asv))
    ALLOCATE(ttk(nttk_asv))
    ALLOCATE(rxbar(ikt_asv,nxrt_asv))
    ALLOCATE(rxbari(ikt_asv,nxrt_asv))

    ! Allocate arrays to handle mass balance coefficients.
    nsts_asv = 10*nst_asv
    nstsmx = nsts_asv

    ALLOCATE(nsts(nsts_asv))
    ALLOCATE(nstsr(2,nst_asv))
    ALLOCATE(csts(nsts_asv))

    ! Initialize the problem counter.
    nprob = 0

    ! The label below is a return point for processing an additional
    ! problem specifed on the input file.
20 continue
    nprob = nprob + 1

    ! Set special species indices.
    nchlor = nchloa
    nhydr = nhydra
    nhydx = nhydxa
    no2gaq = no2gaa
    nelect = neleca
    nrdxsp = nrdxsa
    nbw = nbwa

    ! Initialize certain arrays to zeros, blanks, etc.
    call initiz(iopt,noptmx)
    call initiz(iopg,nopgmx)
    call initiz(iopr,noprmx)
    call initiz(iodb,nodbmx)

    call initcb(uzveci,kmax)
    call initcb(uzvec1,kmax)
    call initaz(zvclg1,kmax)
    call initaz(zvec1,kmax)

    ! Note: cspi and uspi are not used in EQ6.
    call initiz(jflgi,nbtmax)

    call initcb(utitl1,ntitmx)
    call initcb(utitl2,ntitmx)

    call initcb(ubmtbi,nbtmax)

    nmax = 2*nbtmax
    call initcb(usbsw,nmax)
    call initcb(uobsw,nmax)

    nmax = nbt_asv
    call initiz(nbasp,nmax)
    call initiz(nbaspd,nmax)
    call initiz(nbaspx,nmax)

    nmax = nst_asv
    call initiz(jflag,nmax)
    call initiz(jflagd,nmax)
    call initiz(jsflag,nmax)
    call initcb(uspec,nmax)
    call initaz(vosp0,nmax)
    call initaz(xlks,nmax)
    call initaz(zchar,nmax)
    call initaz(zchsq2,nmax)
    call initaz(zchcu6,nmax)

    ! xxxxxxxxxxxx
    call initaz(cdrs,ndrsmx)
    call initaz(cdrsx,ndrsmx)
    call initiz(ndrs,ndrsmx)
    call initiz(ndrsx,ndrsmx)

    ! xxxxxxxxxxxx
    nmax = 2*nst_asv
    call initiz(ndrsr,nmax)
    call initiz(ndrsrx,nmax)

    nmax = narx_asv*ntpr_asv*nst_asv
    call initaz(axlks,nmax)
    call initaz(axhfs,nmax)
    call initaz(axvfs,nmax)

    ! xxxxxxxxxxxx
    call initaz(axlksx,nmax)
    call initaz(axhfsx,nmax)
    call initaz(axvfsx,nmax)

    nmax = narxmx*ntprmx*ipchmx*nstmax
    call initaz(adhfs,nmax)
    call initaz(adhfsx,nmax)

    nmax = narxmx*ntprmx*ipcvmx*nstmax
    call initaz(advfs,nmax)
    call initaz(advfsx,nmax)

    ! xxxxxxxxxxxx
    call initiz(iffg,nffgmx)
    call initiz(jffg,nffgmx)

    call initcb(uffg,nffgmx)
    call initaz(moffg,nffgmx)
    call initaz(xlkffg,nffgmx)

    call initiz(jcode,nrctmx)
    call initiz(jreac,nrctmx)
    call initiz(nsk,nrctmx)

    call initcb(ureac,nrctmx)
    call initaz(fkrc,nrctmx)
    call initaz(modr,nrctmx)
    call initaz(morr,nrctmx)
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

    call initcb(ugexp,netmax)

    do ne = 1,netmax
        jgext(ne) = 0
        kern1(ne) = 0
        kern2(ne) = 0
    end do

    nmax = jetmax*netmax
    call initiz(jern1,nmax)
    call initiz(jern2,nmax)
    call initiz(ngext,nmax)
    call initiz(ngexrt,nmax)
    call initaz(egexjf,nmax)

    nmax = ietmax*jetmax*netmax
    call initaz(cegexs,nmax)
    call initaz(cpgexs,nmax)
    call initaz(mrgexs,nmax)
    call initiz(ngexro,nmax)
    call initiz(ngexso,nmax)
    call initiz(ngexsa,nmax)

    nmax = ietmax*jetmax*nertmx
    call initaz(egersi,nmax)
    call initaz(xgersi,nmax)

    nmax = ketmax*netmax
    call initiz(kgexsa,nmax)

    nmax = iktmax*nxrtmx
    call initaz(rxbari,nmax)

    nmax = iapxmx*nxtmax
    call initaz(apx,nmax)

    nmax = ibpxmx*nxtmax
    call initaz(bpx,nmax)

    ! Determine the problem input format.
    read (ninpts,1090,end=105,err=107) ux8
1090 format(a8)

    backspace ninpts
    uinfor = 'W'

    if (ux8(1:8) .eq. '|-------') then
        uinfor = 'D'
    end if

    ! Read the problem input.
    if (uinfor(1:1) .eq. 'W') then
        ! Compact (W) format.
        call rd6inw(awmaxi,awmini,cbsri,cdac,cesri,cgexj,csigma,dlaplo,dlaprn,dleplo,dleprn,dlhplo,dlhprn,dloplo,dloprn,dltpll,dltplo,dltprl,dltprn,dlxdmp,dlxmx0,dlxpll,dlxplo,dlxprl,dlxprn,eact,egersi,ehmaxi,ehmini,electr,fkrc,iact,ibsrti,iesrti,ietmax,iktmax,imchmx,imech,iodb,iopg,iopr,iopt,igerti,itermx,ixrti,jcode,jgerti,jetmax,jflgi,jgext,jpress,jreac,jtemp,kbt,kct,kdim,kmax,kmt,kprs,ksplmx,ksppmx,kstpmx,kxmod,kxt,hact,modr,moffg,morr,mprphi,mprspi,mtbaqi,mtbi,mwtges,nbti,nbtmax,nbt1mx,nctmax,ndact,ndctmx,nert,nertmx,net,netmax,nffg,nffgmx,ngexrt,ninpts,nobswt,nodbmx,nopgmx,noprmx,noptmx,nordmx,noutpt,nprob,nprpmx,nprpti,nprsmx,nprsti,nptkmx,nrct,nrctmx,nrk,nsbswt,nsk,nsrt,nsrtmx,ntitl1,ntitl2,ntitmx,ntrymx,nttkmx,nttyo,nxmdmx,nxmod,nxopex,nxopmx,nxopt,nxpemx,nxrt,nxrtmx,o2maxi,o2mini,phmaxi,phmini,pressb,pressi,ptk,qend,qgexsh,qrderr,rkb,rxbari,sfcar,ssfcar,tempcb,tempci,tgexp,timmxi,tistti,tolbt,toldl,tolsat,tolxsf,trkb,ttk,ubmtbi,ubsri,ucxri,udac,uesri,uffg,ugerji,ugermo,ugersi,ugexj,ugexmo,ugexp,ugexr,uhfgex,uobsw,uprphi,uprspi,ureac,usbsw,utitl1,utitl2,uvfgex,uxcat,uxkgex,uxmod,uxopex,uxopt,uzveci,vreac,xgersi,xhfgex,ximaxi,xistti,xlkffg,xlkgex,xlkmod,xvfgex,zgexj,zvclgi)
    else
        ! Menu-style (D) format.
        call rd6ind(awmaxi,awmini,cbsri,cdac,cesri,cgexj,csigma,dlaplo,dlaprn,dleplo,dleprn,dlhplo,dlhprn,dloplo,dloprn,dltpll,dltplo,dltprl,dltprn,dlxdmp,dlxmx0,dlxpll,dlxplo,dlxprl,dlxprn,eact,egersi,ehmaxi,ehmini,electr,fkrc,iact,ibsrti,iesrti,ietmax,iktmax,imchmx,imech,iodb,iopg,iopr,iopt,igerti,itermx,ixrti,jcode,jgerti,jetmax,jflgi,jgext,jpress,jreac,jtemp,kbt,kct,kdim,kmax,kmt,kprs,ksplmx,ksppmx,kstpmx,kxmod,kxt,hact,modr,moffg,morr,mprphi,mprspi,mtbaqi,mtbi,mwtges,nbti,nbtmax,nbt1mx,nctmax,ndact,ndctmx,nert,nertmx,net,netmax,nffg,nffgmx,ngexrt,ninpts,nobswt,nodbmx,nopgmx,noprmx,noptmx,nordmx,noutpt,nprob,nprpmx,nprpti,nprsmx,nprsti,nptkmx,nrct,nrctmx,nrk,nsbswt,nsk,nsrt,nsrtmx,ntitl1,ntitl2,ntitmx,ntrymx,nttkmx,nttyo,nxmdmx,nxmod,nxopex,nxopmx,nxopt,nxpemx,nxrt,nxrtmx,o2maxi,o2mini,phmaxi,phmini,pressb,pressi,ptk,qend,qgexsh,qrderr,rkb,rxbari,sfcar,ssfcar,tempcb,tempci,tgexp,timmxi,tistti,tolbt,toldl,tolsat,tolxsf,trkb,ttk,ubmtbi,ubsri,ucxri,udac,uesri,uffg,ugerji,ugermo,ugersi,ugexj,ugexmo,ugexp,ugexr,uhfgex,uobsw,uprphi,uprspi,ureac,usbsw,utitl1,utitl2,uvfgex,uxcat,uxkgex,uxmod,uxopex,uxopt,uzveci,vreac,xgersi,xhfgex,ximaxi,xistti,xlkffg,xlkgex,xlkmod,xvfgex,zgexj,zvclgi)
    end if

    go to 109
105 continue
    qend = .true.
    go to 109
107 continue
    qrderr = .true.
109 continue

    if (qrderr) then
        write (noutpt,1100)
        write (nttyo,1100)
1100 format(/' * Error - (EQ6/eq6) Encountered an error while',' reading the input file.',/7x,'The line associated with',' the problem should be the last one',/7x,'echoed on the',' output file or the first line immediately thereafter.')

        stop
    end if

    if (qend) then
        write (noutpt,1110)
        write (nttyo,1110)
1110 format(/' No further input found.',/)

        ! Make porting changes in the EQLIBU routines that are called in
        ! this section. Do not make the porting changes here.
        ! Get end time and date. Also get the run, user, and cpu times.
        call runtim(iexec0,jexec0,texec0,noutpt,nttyo,trun,tuser,tcpu,udate1,utime1)

        j2 = ilnobl(udate0)
        j3 = ilnobl(udate1)
        write (noutpt,1120) utime0,udate0(1:j2),utime1,udate1(1:j3)
        write (nttyo,1120) utime0,udate0(1:j2),utime1,udate1(1:j3)
1120 format(10x,'Start time = ',a8,2x,a,/12x,'End time = ',a8,2x,a)

        ! Print the run, user, and cpu times.
        write (noutpt,1140) trun
        write (nttyo,1140)  trun
1140 format(/10x,' Run time = ',g10.3,' seconds')

        if (tuser .gt. 0.) then
            write (noutpt,1150) tuser
            write (nttyo,1150)  tuser
1150 format(10x,'User time = ',g10.3,' seconds')
        end if

        if (tcpu .gt. 0.) then
            write (noutpt,1160) tcpu
            write (nttyo,1160)  tcpu
1160 format(10x,' Cpu time = ',g10.3,' seconds')
        end if

        write (noutpt,1170)
        write (nttyo,1170)
1170 format(/' Normal exit')

        ! BEGIN_MACHINE_DEPENDENT_CODE
        !   Clear the IEEE flag for floating-point underflow, if such a
        !   flag is present, to avoid getting an unnecessary system
        !   warning message. Underflow is a normal condition in EQ3/6.
        !   Make porting changes in the EQLIBU subroutine that is called
        !   in this section. Do not make the porting changes here.
        call cliefu()

        ! END_MACHINE_DEPENDENT_CODE
        ! Close and delete the stripped input file and the tabs file.
        close (ninpts,status='delete')

        if (ntabs .gt. 0) then
            close (ntabs,status='delete')
        end if

        stop
    end if

    ! Scan the input file title for option strings.
    ! Scan the input file title for the USEOLDPITZERMU string.
    ! If found and if iopg(1) = 1 (use Pitzer's equations), set up
    ! set to calculate activity coefficients using recalculated
    ! observable third-order parameters (Cphi, psi, zeta) instead
    ! of their conventional mu equivalents.
    qhawep = .true.

    if (iopg(1) .eq. 1) then
        do n = 1,ntitl1
            ux80 = utitl1(n)
            call locase(ux80)
            i = index(ux80,'useoldpitzermu')

            if (i .gt. 0) then
                qhawep = .false.
                go to 106
            end if
        end do

        do n = 1,ntitl2
            ux80 = utitl2(n)
            call locase(ux80)
            i = index(ux80,'useoldpitzermu')

            if (i .gt. 0) then
                qhawep = .false.

                if (ntitl1 .lt. ntitmx) then
                    ntitl1 = ntitl1 + 1
                    utitl1(ntitl1) = utitl2(n)
                else
                    write (noutpt,1180)
                    write (nttyo,1180)
1180 format(/' * Note - (EQ6/eq6) Found the USEOLDPITZERMU',' option string',/7x,'in the secondary input file title',' but cannot add it',/7x,'to the main title. The option',' will not carry over to any',/7x,'PICKUP file that may',' be written at the end of the present run.')
                end if

                go to 106
            end if
        end do

106 continue
    end if

    if (.not.qhawep) then
        write (noutpt,1190)
        write (nttyo,1190)
1190 format(/' * Note - (EQ6/eq6) Found the USEOLDPITZERMU',' option string',/7x,'in one of the input file titles.',' Will evaluate the Pitzer',/7x,'equations using the',' lambda-mu format. Implied psi coefficients',/7x,'not',' on the supporting data file will not be effectively',' treated',/7x,'as having zero values. Rather, the',' corresponding mu coefficients',/7x,'will be treated as',' having zero values. In effect, an implied',/7x,'psi',' coefficient value will include contributions from',' the',/7x,'two cognate Cphi coefficients that appear',' in the psi-mu',/7x,'breakdown equation.')
    end if

    ! Scan the input file title for the USEOLDPITZER75 string.
    ! If found and if iopg(1) = 1 (use Pitzer's equations), use
    ! the older approximation given by Pitzer (1975) for calculating
    ! the higher-order electrostatic terms.
    qpit75 = .false.

    if (iopg(1) .eq. 1) then
        do n = 1,ntitl1
            ux80 = utitl1(n)
            call locase(ux80)
            i = index(ux80,'useoldpitzer75')

            if (i .gt. 0) then
                qpit75 = .true.
                go to 170
            end if
        end do

        do n = 1,ntitl2
            ux80 = utitl2(n)
            call locase(ux80)
            i = index(ux80,'useoldpitzer75')

            if (i .gt. 0) then
                qpit75 = .true.

                if (ntitl1 .lt. ntitmx) then
                    ntitl1 = ntitl1 + 1
                    utitl1(ntitl1) = utitl2(n)
                else
                    write (noutpt,1192)
                    write (nttyo,1192)
1192 format(/' * Warning - (EQ6/eq6) Found the USEOLDPITZER75',' option string',/7x,'in the secondary input file title',' but cannot add it',/7x,'to the main title. The option',' will not carry over to any',/7x,'PICKUP file that may',' be written at the end of the present run.')
                end if

                go to 170
            end if
        end do

170 continue
    end if

    if (qpit75) then
        write (noutpt,1193)
        write (nttyo,1193)
1193 format(/' * Warning - (EQ6/eq6) Found the USEOLDPITZER75',' option string',/7x,'in one of the input file titles.',' Will evaluate the Pitzer',/7x,'equations using the',' old Pitzer (1975) approximation for',/7x,'higher order',' electrostatic terms, not the later Harvie (1981)',/7x,'approximation that is now used nearly universally',' to evaluate',/7x,'the Pitzer equations.')
    end if

    ! Scan the input file title for the WRITEPITZERJTABLES string.
    ! If found and if iopg(1) = 1 (use Pitzer's equations), calculate
    ! and output tables of the J(x) and J'(x) functions. This is a
    ! one-time only option. It is not carried forward on the PICKUP
    ! file.
    qcwrpj = .false.

    if (iopg(1) .eq. 1) then
        do n = 1,ntitl1
            ux80 = utitl1(n)
            call locase(ux80)
            i = index(ux80,'writepitzerjtables')

            if (i .gt. 0) then
                qcwrpj = .true.
                go to 180
            end if
        end do

        do n = 1,ntitl2
            ux80 = utitl2(n)
            call locase(ux80)
            i = index(ux80,'writepitzerjtables')

            if (i .gt. 0) then
                qcwrpj = .true.
                go to 180
            end if
        end do

180 continue
    end if

    if (qcwrpj) then
        write (noutpt,1194)
        write (nttyo,1194)
1194 format(/' * Note - (EQ6/eq6) Found the WRITEPITZERJTABLES',' option string',/7x,'in one of the input file titles.',' Will calculate the Pitzer',/7x,"J(x) and J'(x) functions",' for higher order electrostatic terms and',/7x,'write output',' tables for both the Pitzer (1975) and Harvie (1981)',/7x,'approximations.')
    end if

    ! Scan the input file title for the TURNOFFOPTIMIZER string.
    ! If found, turn off the pre-Newton-Raphson optimizer.
    qoptmz = .true.

    do n = 1,ntitl1
        ux80 = utitl1(n)
        call locase(ux80)
        i = index(ux80,'turnoffoptimizer')

        if (i .gt. 0) then
            qoptmz = .false.
            go to 110
        end if
    end do

    do n = 1,ntitl2
        ux80 = utitl2(n)
        call locase(ux80)
        i = index(ux80,'turnoffoptimizer')

        if (i .gt. 0) then
            qoptmz = .false.

            if (ntitl1 .lt. ntitmx) then
                ntitl1 = ntitl1 + 1
                utitl1(ntitl1) = utitl2(n)
            else
                write (noutpt,1195)
                write (nttyo,1195)
1195 format(/' * Warning - (EQ6/eq6) Found the TURNOFFOPTIMIZER',' option string',/7x,'in the secondary input file title',' but cannot add it',/7x,'to the main title. The option',' will not carry over to any',/7x,'PICKUP file that may',' be written at the end of the present run.')
            end if

            go to 110
        end if
    end do

110 continue

    if (.not.qoptmz) then
        write (noutpt,1197)
        write (nttyo,1197)
1197 format(/' * Note - (EQ6/eq6) Found the TURNOFFOPTIMIZER',' option string',/7x,'in one of the input file titles.',' Will turn off the',/7x,'pre-Newton-Raphson optimizer.',' This may help if Newton-Raphson',/7x,'iteration fails',' to converge due to insufficient optimizer',/7x,'performance. This is most likely to be a problem',' when the',/7x,'phase assemblage changes.')
    end if

    ! Scan the input file title for the TABFILEASTXT string. If found,
    ! found, make the TAB file an ordinary text file (the old style
    ! TAB file) instead of as a comma separated value (.csv) file.
    qtatxt = .false.

    do n = 1,ntitl1
        ux80 = utitl1(n)
        call locase(ux80)
        i = index(ux80,'tabfileastxt')

        if (i .gt. 0) then
            qtatxt = .true.
            go to 210
        end if
    end do

    do n = 1,ntitl2
        ux80 = utitl2(n)
        call locase(ux80)
        i = index(ux80,'tabfileastxt')

        if (i .gt. 0) then
            qtatxt = .true.

            if (ntitl1 .lt. ntitmx) then
                ntitl1 = ntitl1 + 1
                utitl1(ntitl1) = utitl2(n)
            else
                write (noutpt,1196)
                write (nttyo,1196)
1196 format(/' * Warning - (EQ6/eq6) Found the TABFILEASTXT',' option string',/7x,'in the secondary input file title',' but cannot add it',/7x,'to the main title. The option',' will not carry over to any',/7x,'PICKUP file that may',' be written at the end of the present run.')
            end if

            go to 210
        end if
    end do

210 continue

    if (qtatxt) then
        write (noutpt,1198)
        write (nttyo,1198)
1198 format(/' * Note - (EQ6/eq6) Found the TABFILEASTXT',' option string',/7x,'in one of the input file titles.',' Will write the TAB',/7x,'file as an ordinary text file',' instead of the default csv',/7x,'(comma-separated-variable)',' file. The data content may not',/7x,'be the same for both',' files.')
    end if

    ! Check consistency between the activity coefficient option and the
    ! data1 file.
    call cdakey(iopg,nopgmx,noutpt,nttyo,udakey,udatfi)

    ! Get the name of the option for calculating the activity
    ! coefficients of aqueous species.
    call nactop(iopg,nopgmx,noutpt,nttyo,uactop)

    if (iopt(16) .ge. 0) then
        inquire(file=bafile,opened=qop)

        if (.not.qop) then
            call openou(noutpt,nttyo,bafile,'formatted',nrecl,nbkupa)
        end if

        if (iopt(16) .eq. 0) then
            inquire(file=bbfile,opened=qop)

            if (.not.qop) then
                call openou(noutpt,nttyo,bbfile,'formatted',nrecl,nbkupb)
            end if
        end if
    end if

    if (iopt(17) .ge. 0) then
        inquire(file=pfile,opened=qop)

        if (.not.qop) then
            call openou(noutpt,nttyo,pfile,'formatted',nrecl,newin)
        end if
    end if

    if (nprob.gt.1 .and. iopt(18).ge.1) then
        iopt(18) = 0
    end if

    if (iopt(18) .ge. 1) then
        inquire(file=txfile,exist=qex)

        if (qex) then
            call openin(noutpt,nttyo,txfile,'formatted',ntabx)

            do i = 1,10000
                read (ntabx,1200,end=200) ux8
1200 format(a8)
            end do

200 continue
        else
            write (noutpt,1210)
            write (nttyo,1210)
1210 format(/' * Error - (EQ6/eq6) Have iopt(18)= 1, but there',/7x,'is no existing tabx file to which to append.')

            stop
        end if
    end if

    if (iopt(18) .ge. 0) then
        ! Set up the TAB, TABS, and TABX files.
        ! For the .csv form, the record length (nrecl) must match the
        ! maximum character lengths of the table tag (8 characters)
        ! and the ulinex variable (nllnmx characters), plus two
        ! (for one comma trailing the table tag and another trailing
        ! The ulinex variable).
        nrecl = nllnmx + 10

        if (qtatxt) then
            nrecl = 128
        end if

        inquire(file=tfile,opened=qop)

        if (.not.qop) then
            call openou(noutpt,nttyo,tfile,'formatted',nrecl,ntab)
        end if

        inquire(file=tsfile,opened=qop)

        if (.not.qop) then
            call openou(noutpt,nttyo,tsfile,'formatted',nrecl,ntabs)
        else
            rewind ntabs
        end if

        nrecl = nllnmx + 10

        if (qtatxt) then
            nrecl = 129
        end if

        inquire(file=txfile,opened=qop)

        if (.not.qop) then
            call openou(noutpt,nttyo,txfile,'formatted',nrecl,ntabx)
        end if

        if (iopt(18) .eq. 0) then
            rewind ntabx
        end if
    end if

    ! qecon = .true. denotes economy mode.
    ! qscon = .true. denotes super economy mode.
    qecon = .false.
    qscon = .false.

    if (iopt(1).lt.2 .and. iopt(2).le.0) then
        if (iopt(13) .ge. 1) then
            qecon = .true.
        end if

        if (iopt(13) .ge. 2) then
            qscon = .true.
        end if
    end if

    ! Set default values of run control parameters read from the
    ! data file. Put any out-of-range values into acceptable range.
    tolxst = 0.
    tolxsu = 0.
    call dfaltz(dlaplo,dlaprn,dleplo,dleprn,dlhplo,dlhprn,dloplo,dloprn,dltplo,dltpll,dltprl,dltprn,dlxdmp,dlxmx0,dlxplo,dlxpll,dlxprl,dlxprn,iopt,itermx,ksplmx,ksppmx,kstpmx,net,noptmx,nordmx,nrct,noutpt,ntrymx,nttyo,prcinf,qecon,qscon,timmxi,tistti,tolbt,toldl,tolsat,tolxsf,tolxst,tolxsu,ximaxi,xistti)
    nrd1mx = nordmx + 1

    if (iopt(2) .gt. 0) then
        ! Disallow print and plot intervals in log time space if the
        ! starting time is less than zero.
        if (tistti.lt.0. .or. timmxi.lt.0.) then
            dltprl = prcinf
            dltpll = prcinf
        end if
    end if

    ! Set other run control parameters.
    call setrcp(aftarg,dlxmax,dlxmin,dlxmx0,npslmx,nsscmx,nsslmx,prcinf,sscrew,tolaft,tolsat,tolsst,zkfac,zklgmn,zklogl,zklogu)

    ! Set some parameters using the data read from the input file.
    nart = 0
    ngrt = 0
    nmrt = 0

    do nrc = 1,nrct
        if (jcode(nrc) .eq. 0) then
            nmrt = nmrt + 1
        else if (jcode(nrc) .eq. 3) then
            nart = nart + 1
        else if (jcode(nrc) .eq. 4) then
            ngrt = ngrt + 1
        end if
    end do

    km1 = kbt + 1
    kx1 = kmt + 1

    ! Set the starting and maximum values of reaction progress and
    ! time to the input values. Note that the values of the input
    ! variables will be changed prior to writing data on a backup
    ! or pickup file.
    xistrt = xistti
    ximax = ximaxi
    tistrt = tistti
    timemx = timmxi

    xi1 = xistrt
    time1 = tistrt

    ! Set the values of other run control parameters. These parameters
    ! are also associated with mechanisms for terminating the current
    ! reaction path; for example, if pH decreases to phmin, or increases
    ! to phmax.
    phmin = phmini
    phmax = phmaxi
    ehmin = ehmini
    ehmax = ehmaxi
    o2min = o2mini
    o2max = o2maxi
    awmin = awmini
    awmax = awmaxi

    ! Check the input reactant data for various kinds of errors
    ! and inconsistencies.
    call chkinz(ier,imchmx,imech,iopt,jcode,kmax,kxt,nelect,noptmx,noutpt,no2gaq,nrct,nrctmx,nrk,nstmax,nttyo,rkb,ureac,uspeca,uzveci,zvclgi)

    if (ier .gt. 0) then
        stop
    end if

    ! Set some defaults.
    do nrc = 1,nrct
        if (fkrc(nrc) .le. 0.) then
            fkrc(nrc) = 1.0
        end if

        if (sfcar(nrc) .lt. 0.) then
            sfcar(nrc) = 0.
        end if

        if (ssfcar(nrc) .lt. 0.) then
            ssfcar(nrc) = 0.
        end if

        do j = 1,2
            if (nrk(j,nrc) .eq. 2) then
                if (imech(j,nrc) .le. 0) then
                    imech(j,nrc) = 1
                end if

                do i = 1,imech(j,nrc)
                    if (csigma(i,j,nrc) .le. 0.) then
                        csigma(i,j,nrc) = 1.
                    end if
                end do
            end if
        end do
    end do

    do n = 1,nsbswt
        if (usbsw(1,n)(1:24) .eq. ublk24(1:24)) then
            usbsw(1,n)(1:24) = uaqsln
        end if

        if (usbsw(2,n)(1:24) .eq. ublk24(1:24)) then
            usbsw(2,n)(1:24) = uaqsln
        end if
    end do

    do n = 1,nobswt
        if (uobsw(1,n)(1:24) .eq. ublk24(1:24)) then
            uobsw(1,n)(1:24) = uaqsln
        end if

        if (uobsw(2,n)(1:24) .eq. ublk24(1:24)) then
            uobsw(2,n)(1:24) = uaqsln
        end if
    end do

    ! Copy the 'a' set of reaction data into the 'd' set. The 'a'
    ! set remains a faithful image of the data file. The 'd' set
    ! at this stage is subject to modification. The number of basis
    ! species may be increased in response to input file directives.
    ! Thus, nbtd may become greater than nbta, and the basis species
    ! pointer array nbaspd lengthened with respect to nbaspa. The
    ! 'd' set at this stage may also be modified by directives on the
    ! input file to execute special basis switching. This redefines
    ! which of the basis species are the strict basis species, and
    ! according causes reactions and attendant thermodynamic data to
    ! be rewritten. Modifications of the 'd' set at this stage
    ! support options affecting problem definition. After the data
    ! have been compressed to eliminate unnecessary species, the 'd'
    ! set will be redfined. This second stage form is equivalent to
    ! a compressed version of the first stage form with additional
    ! input file directed modifications.
    ! Calling sequence substitutions:
    !   adhfsa for adhfs
    !   advfsa for advfs
    !   axhfsa for axhfs
    !   axlksa for axlks
    !   axvfsa for axvfs
    !   cdrsa for cdrs
    !   nbaspa for nbasp
    !   ndrsa for ndrs
    !   ndrsra for ndrsr
    call cdrssd(adhfsa,adhfsd,advfsa,advfsd,axhfsa,axhfsd,axlksa,axlksd,axvfsa,axvfsd,cdrsa,cdrsd,ipch,ipchmx,ipcv,ipcvmx,narxmx,nbaspa,nbaspd,nbtmax,ndrsa,ndrsd,ndrsmx,ndrsra,ndrsrd,nstmax,ntprmx)
    nbtd = nbta

    ! Get the basis index of water (nbw).
    ! Calling sequence substitutions:
    !   nbaspd for nbasp
    !   nbtd for nbt
    !   narn1a for ns
    nbw = nbasis(nbaspd,nbtd,nbtmax,narn1a)

    ! Execute special basis switching to redefine which basis species
    ! are the strict basis species in the 'd' set of reaction data
    ! (nbtd/nbaspd/cdrsd/ndrsd/ndrsrd/axlksd). This type of basis
    ! switching is executed to assist in problem definition. Ordinary
    ! basis switching, executed later in this code, is done only
    ! for numerical purposes.
    if (nsbswt .gt. 0) then
        do nsbsw = 1,nsbswt
            ! Interpret the directive for the nsbsw-th switch.
            call intsbs(nb1,nb2,nbaspd,nbtd,nbtmax,noutpt,ns1,ns2,nsbsw,nstmax,nttyo,usbsw,uspeca)

            ! Execute the switch.
            call swtchb(adhfsd,adhfsx,advfsd,advfsx,axhfsd,axhfsx,axlksd,axlksx,axvfsd,axvfsx,cdrsd,cdrsx,ipch,ipchmx,ipcv,ipcvmx,narxmx,nbaspd,nbtmax,nbw,nb1,nb2,ndrsd,ndrsmx,ndrsx,ndrsrd,ndrsrx,noutpt,ns1,ns2,nsta,nstmax,ntprmx,nttyo,uspeca)
        end do
    end if

    ! Initialize some arrays.
    ! Note: csp and nsp are not used in EQ6.
    call initiz(ndecsp,nbtmax)

    ! Note: mtb, mtbaq, mte, and mteaq are not used in EQ3NR.
    call initaz(mtb,nbtmax)
    call initaz(mtbaq,nbtmax)
    call initaz(mte,nctmax)
    call initaz(mteaq,nctmax)

    ! Note: conc and conclg are initialized in EQ6/path.f.
    call initaz(moph,nptmax)
    call initaz(mosp,nstmax)
    av = -99999.
    call initav(loph,nptmax,av)
    call initav(losp,nstmax,av)

    ! Note: mprsp is not used in EQ3NR.
    call initaz(mprsp,nstmax)

    ! Note: acflg, acflgo, act, actlg, xbar, and xbarlg, are initialized
    ! in EQ6/path.f.
    ! Note: npchk and mprph are not used in EQ3NR.
    call initiz(npchk,nptmax)
    call initaz(mprph,nptmax)

    ! Note: zvec0 and zvclg0 are initialized in EQ6/path.f.
    call initiz(iindx1,kmax)
    call initiz(ipndx1,kmax)

    call initaz(zvec1,kmax)

    av = -99999.
    call initav(zvclg1,kmax,av)

    ! Note: the following are not used in EQ3NR.
    nmax = imchmx*2*nrctmx
    call initaz(rk,nmax)

    nmax = ndctmx*imchmx*2*nrctmx
    call initiz(ndac,nmax)

    nmax = iktmax*nxrtmx
    call initaz(rxbar,nmax)

    nmax = nctmax*nsrtmx
    call initaz(cesr,nmax)

    nmax = nbt1mx*nsrtmx
    call initaz(cbsr,nmax)

    nmax = ietmax*jetmax*nertmx
    call initaz(egers,nmax)
    call initaz(xgers,nmax)
    call initaz(mrgers,nmax)

    ! Interpret the data file basis species listed on the input file.
    ! "Data file" basis species to be created according to instructions
    ! read from the input file (e.g., for generic ion exchangers)
    ! will be ignored at this point, as creation occurs later in
    ! this code.
    call intbs6(jflag,jflgi,kmax,narn1a,narn2a,nbaspd,nbtd,nbti,nbtmax,ndrsrd,ndecsp,noutpt,nsta,nstmax,nttyo,uspeca,ubmtbi)

    ! Set jflag to -1 for species that can not appear in the system.
    call jflaux(jflag,nbaspd,nbtd,nbtmax,ndrsd,ndrsmx,ndrsrd,nstmax)

    ! Test the temperature-tracking input.
    if (jtemp .eq. 0) then
        continue
    else if (jtemp .eq. 1) then
        if (ttk(1) .eq. 0.) then
            write (noutpt,2110)
            write (nttyo,2110)
2110 format(/' * Note (EQ6/eq6) The specified temperature',/7x,'derivative for computing the temperature as a linear',/7x,'function of Xi is zero. The temperature is therefore',/7x,'constant. Will reset the temperature-tracking flag',/7x,'jtemp from 1 to 0.')

            jtemp = 0
        end if
    else if (jtemp .eq. 2) then
        if (ttk(1) .eq. 0.) then
            write (noutpt,2120)
            write (nttyo,2120)
2120 format(/' * Note (EQ6/eq6) The specified temperature',/7x,'derivative for computing the temperature as a linear',/7x,'function of time is zero. The temperature is therefore',/7x,'constant. Will reset the temperature-tracking flag',/7x,'jtemp from 2 to 0.')

            jtemp = 0
        end if
    else if (jtemp .eq. 3) then
        if (ttk(1) .eq. 0.) then
            write (noutpt,2130)
            write (nttyo,2130)
2130 format(/' * Note (EQ6/eq6) The specified mass factor',/7x,'ratio (ttk(1) for computing the temperature from',/7x,'a fluid mixing ratio (jtemp = 3) was zero. This is',/7x,'not a valid value. Will change this to 1.0. This',/7x,'parameter affects the relation between the fluid',/7x,'mixing ratio and Xi. A value of 1.0 corresponds to',/7x,'a 1:1 mixing ratio.')
        end if
    else
        write (noutpt,2140) jtemp
        write (nttyo,2140) jtemp
2140 format(/' * Error - (EQ6/eq6) The temperature tracking',/7x,'flag (jtemp) has an unknown value of ',i2,'.')

        stop
    end if

    ! Test the pressure-tracking input.
    if (jpress .eq. 0) then
        continue
    else if (jpress .eq. 1) then
        continue
    else if (jpress .eq. 2) then
        continue
    else if (jpress .eq. 3) then
        if (ptk(1) .eq. 0.) then
            write (noutpt,2150)
            write (nttyo,2150)
2150 format(/' * Note (EQ6/eq6) The specified pressure',/7x,'derivative for computing the pressure as a linear',/7x,'function of Xi is zero. The pressure is therefore',/7x,'constant. Will reset the pressure-tracking flag',/7x,'jpress from 3 to 2.')

            jpress = 0
        end if
    else if (jpress .eq. 4) then
        if (ptk(1) .eq. 0.) then
            write (noutpt,2160)
            write (nttyo,2160)
2160 format(/' * Note (EQ6/eq6) The specified pressure',/7x,'derivative for computing the pressure as a linear',/7x,'function of time is zero. The pressure is therefore',/7x,'constant. Will reset the pressure-tracking flag',/7x,'jpress from 4 to 2.')

            jpress = 0
        end if
    else
        write (noutpt,2170) jpress
        write (nttyo,2170) jpress
2170 format(/' * Error - (EQ6/eq6) The pressure tracking',/7x,'flag (jpress) has an unknown value of ',i2,'.')

        stop
    end if

    ! Compute the initial temperature.
    call gtemp(afcnst,al10,iopt,jtemp,noptmx,noutpt,nttkmx,nttyo,rconst,rtcnst,tempc,tempcb,tempk,time1,ttk,xi1)

    ! Determine the corresponding temperature range flag.
    call gntpr(ntpr,ntprmx,ntprt,tempc,tempcu)

    ! Check for constant temperature.
    qcntmp = jtemp .eq. 0

    ! Check for a temperature jump.
    dt = tempc - tempci

    if (abs(dt) .gt. 1.e-4) then
        write (noutpt,1300) tempc,tempci
        write (nttyo,1300) tempc,tempci
1300 format(/' * Note - (EQ6/eq6) The temperature is jumping from',/7x,'the previous run:',/,/'     Initial temperature      = ',f10.4,' C',/'     Previous run temperature = ',f10.4,' C')
    end if

    ! Compute the data file reference grid pressure at the initial
    ! temperature.
    ! Calling sequence substitutions:
    !   apresg for arr
    !   presg for prop
    call evdat2(apresg,narxmx,narxt,ntpr,ntprmx,presg,tempc)

    ! Compute the 1.013-bar/steam-saturation curve pressure at the
    ! initial temperature.
    if (tempc .le. 100.) then
        ntprh = 1
    else
        ntprh = 2
    end if

    ! Calling sequence substitutions:
    !   apresh for arr
    !   5 for narxmx
    !   narxth for narxt
    !   ntprh for ntpr
    !   2 for ntprmx
    !   presh for prop
    call evdat2(apresh,5,narxth,ntprh,2,presh,tempc)

    ! Compute the initial pressure.
    call gpress(iopt,jpress,noptmx,noutpt,nptkmx,nttyo,presg,presh,press,pressb,time1,ptk,xi1)

    ! Check for constant pressure.
    qcnpre = jpress.eq.2 .or. (jpress.eq.0 .and. qcntmp) .or. (jpress.eq.1 .and. qcntmp)

    ! Check for a pressure jump.
    dp = press - pressi

    if (abs(dp) .gt. 1.e-4) then
        write (noutpt,1310) press,pressi
        write (nttyo,1310) press,pressi
1310 format(/' * Note - (EQ6/eq6) The pressure is jumping from',/7x,'the previous run:',/,/'     Initial pressure      = ',1pg12.5,' bars',/'     Previous run pressure = ',1pg12.5,' bars')
    end if

    ! Set species status flags.
    call flgset(axlksd,iopt,jflag,jpflag,jsflag,kxmod,narn1a,narn2a,narxmx,nbaspd,nbtd,nbtmax,ncmpra,ncta,ndrsd,ndrsmx,ndrsrd,noptmx,noutpt,npta,nptmax,nrdxsp,nsta,nstmax,ntpr,ntprmx,nttyo,nxmdmx,nxmod,uphasa,uptypa,uspeca,uxmod)

    if (iopt(15) .ge. 1) then
        ! Execute the option to suppress all redox reactions.
        call suprdx(jflag,jsflag,narn1a,narn2a,ndrsd,ndrsmx,ndrsrd,nrdxsp,nsta,nstmax,uspeca)
    end if

    ! Execute the pure mineral subset selection suppression options.
    call mincsp(cdrsd,jpflag,jsflag,nbaspd,nbtd,nbtmax,ncmpra,ndrsd,ndrsmx,ndrsrd,nmrn1a,nmrn2a,noutpt,npta,nptmax,nstmax,nttyo,nxopex,nxopmx,nxopt,nxpemx,uspeca,uxcat,uxopex,uxopt)

    ! Check the jpflag array to make sure it is consistent with the
    ! jsflag array.
    call flgchk(jpflag,jsflag,ncmpra,npta,nptmax,nstmax,qclnsa)

    ! Examine each active auxiliary basis species. Print a warning if
    ! any other species in the corresponding dissociation reaction is
    ! not present in the model.
    call bspchk(jsflag,nbaspd,nbtd,nbtmax,ndrsd,ndrsmx,ndrsrd,noutpt,nrdxsp,nstmax,nttyo,uspeca)

    ! Do data array compression. Write working data arrays that
    ! don't include phases and species that aren't necessary
    ! for the current problem.
    call cmpdat(adhfs,adhfsd,advfs,advfsd,amu,amua,apx,apxa,aslm,aslma,atwt,atwta,axhfs,axhfsd,axlks,axlksd,axvfs,axvfsd,azero,azeroa,bpx,bpxa,cdrs,cdrsd,cess,cessa,iapxmx,iapxt,iapxta,iaqsla,iaqsln,ibpxmx,ibpxt,ibpxta,ilrn1,ilrn2,imrn1,imrn2,ipbtmx,ipch,ipchmx,ipcv,ipcvmx,insgf,insgfa,iopg,ixrn1,ixrn1a,ixrn2,ixrn2a,jflag,jpfcmx,jpflag,jsflag,jsitex,jsol,jsola,mwtsp,mwtspa,nalpaa,nalpha,napmax,napt,napta,narn1,narn1a,narn2,narn2a,narxmx,nat,natmax,nbasp,nbaspd,nbmap,nbt,nbtd,nbti,nbtmax,nchlor,ncmap,ncmpr,ncmpra,nct,ncta,nctmax,ndecsp,ndrs,ndrsd,ndrsmx,ndrsr,ndrsrd,ness,nessa,nessmx,nessr,nessra,ngrn1,ngrn1a,ngrn2,ngrn2a,ngt,nlrn1,nlrn1a,nlrn2,nlrn2a,nlt,nmrn1,nmrn1a,nmrn2,nmrn2a,nmt,nmut,nmuta,nmutmx,nmux,nmuxa,nopgmx,nslt,nsltmx,nslta,nslx,nslxa,nphasx,npt,npta,nptmax,nsmap,nsta,nst,nstmax,ntf1,ntf1a,ntf1mx,ntf1t,ntf1ta,ntf2,ntf2a,ntf2mx,ntf2t,ntf2ta,ntprmx,nxrn1,nxrn1a,nxrn2,nxrn2a,nxt,nxtmax,palpaa,palpha,qchlor,tf1,tf1a,tf2,tf2a,uelem,uelema,uphasa,uphase,uptypa,uptype,uspec,uspeca,vosp0,vosp0a,zchar,zchara)

    ! Get the indices of special species after compression.
    ! Get the indices of H+, OH-, Cl-, fictive aqueous O2(g), and
    ! fictive aqueous e-.
    call gspion(narn1,narn2,nchlor,nelect,nhydr,nhydx,noutpt,no2gaq,nstmax,nttyo,uspec)

    ! Get the index of the redox basis species.
    ! Calling sequence substitutions:
    !   nbaspd for nbasp
    !   nbtd for nbt
    !   ndrsrd for ndrsr
    call grdxsp(nbaspd,nbtd,nbtmax,nct,ndrsrd,noutpt,nrdxsp,nstmax,nttyo,uspec)

    ! Get the basis index of water (nbw).
    ! Calling sequence substitutions:
    !   nbaspd for nbasp
    !   nbtd for nbt
    !   narn1 for ns
    nbw = nbasis(nbaspd,nbtd,nbtmax,narn1)

    ! Examine each active auxiliary basis species. Change the
    ! corresponding log K polynomial coefficients so that log K is
    ! fixed at a value of -99999 if any other species in the
    ! corresponding dissociation reaction is not in the model.
    call bsplkp(axlks,narxmx,nbasp,nbt,nbtmax,ndrs,ndrsmx,ndrsr,nstmax,ntprmx)

    ! Interpret input file directives to create generic ion-exchange
    ! phases and species.
    call intexi(al10,axhfs,axlks,axvfs,cegexs,cess,cdrs,cgexj,cpgexs,egexjf,iern1,iern2,ietmax,jern1,jern2,jetmax,jflag,jgext,jpflag,jsflag,jsitex,kern1,kern2,ketmax,kgexsa,mwtges,mwtsp,narn1,narn2,narxmx,narxt,nbasp,nbt,nbtmax,ncmpr,ndrs,ndrsmx,ndrsr,ness,nessmx,nessr,netmax,net,nern1,nern2,ngexro,ngexrt,ngexsa,ngexso,ngext,noutpt,nphasx,npt,nptmax,nst,nstmax,ntprt,ntprmx,nttyo,nvetmx,rconst,tgexp,ugexj,ugexmo,ugexmv,ugexp,ugexr,ugexs,ugexsr,uhfgex,uphase,uspec,uvfgex,uxkgex,xhfgex,xlkgex,xvfgex,zchar,zgexj)

    if (net .gt. 0) then
        ! Set jflag values for the exchanger species.
        do ne = 1,net
            do je = 1,jgext(ne)
                nrr1 = jern1(je,ne)
                nrr2 = jern2(je,ne)

                do ns = nrr1,nrr2
                    ! Calling sequence substitutions:
                    !   nbaspd for nbasp
                    nb = nbasis(nbasp,nbt,nbtmax,ns)

                    if (nb .eq. 0) then
                        jflag(ns) = 30
                    end if
                end do
            end do
        end do
    end if

    if (net .gt. 0) then
        ! Echo a table for the generic ion exchangers, describing the
        ! setup of species, reactions, and corresponding thermodynamic
        ! data.
        call echgex(axlks,cdrs,cgexj,iern1,iern2,jern1,jern2,jetmax,jgext,jpflag,jsflag,narxmx,narxt,ndrs,ndrsmx,ndrsr,netmax,noutpt,nptmax,ntprmx,ntprt,nstmax,press,tempc,ugexj,ugexmo,uphase,uspec,xlks)
    end if

    ! Alter any log K values as directed by the input file
    ! (actually, it is the set of interpolating polynomial
    ! coefficients which is altered).
    if (nxmod .gt. 0) then
        call alters(afcnst,apresg,axlks,cdrs,kxmod,narxmx,narxt,ndrs,ndrsmx,ndrsr,noutpt,npt,nptmax,nst,nstmax,ntpr,ntprmx,nttyo,nxmdmx,nxmod,tempc,uphase,uspec,uxmod,xlkmod)
    end if

    ! Create fictive minerals for the fixed fugacity option.
    nfrn1 = nst + 1
    nfrn2 = nst
    ifrn1 = npt + 1
    ifrn2 = npt

    if (nffg .gt. 0) then
        call nlkffg(axlks,cess,cdrs,iffg,ifrn1,ifrn2,jffg,jpflag,jsflag,mwtsp,narxmx,ncmpr,ndrs,ndrsmx,ndrsr,ness,nessmx,nessr,nffg,nffgmx,nfrn1,nfrn2,ngrn1,ngrn2,noutpt,nphasx,npt,nptmax,nst,nstmax,ntpr,ntprmx,nttyo,qcntmp,uffg,ufixf,uphase,uspec,vosp0,xlkffg)
    end if

    ! Make up (z**2)/2 and (z**3)/6 for later use.
    do ns = 1,nst
        zx = 0.5*zchar(ns)*zchar(ns)
        zchsq2(ns) = zx
        zchcu6(ns) = (zx*zchar(ns))/3.
    end do

    ! Get the max norm of the electrical charges of the aqueous
    ! species (izmax).
    call zsrt(izmax,narn1,narn2,nstmax,zchar)

    ! Determine if the model to be calculated has a redox aspect.
    ! Note: the iopt(15) option may be been excecuted above.
    call tstrdx(cdrs,iodb,iopt,jflag,jsflag,narn1,narn2,ndrs,ndrsmx,ndrsr,nodbmx,noptmx,noutpt,nrdxsp,nstmax,qredox,uspec)

    ! If there is no redox aspect, make sure that all active redox
    ! reactions are suppressed.
    if (.not.qredox .and. iopt(15).le.0) then
        ! Calling sequence substitutions:
        !   narn1 for narn1a
        !   narn2 for narn2a
        !   ndrs for ndrsd
        !   ndrsr for ndrsrd
        !   nst for nsta
        !   uspec for uspeca
        call suprdx(jflag,jsflag,narn1,narn2,ndrs,ndrsmx,ndrsr,nrdxsp,nst,nstmax,uspec)
    end if

    ! Interpret the mass balance totals. Construct the mtb and mtbaq
    ! arrays.
    call intmtb(mtb,mtbaq,mtbaqi,mtbi,nbasp,nbt,nbti,nbtmax,noutpt,nstmax,nttyo,ubmtbi,uspec)

    ! Interpret matrix variables. Construct the iindx1 and zvclg1
    ! arrays. Note that new fictive fugacity-fixing phases will not
    ! yet appear in the matrix after this has been done.
    call intmat(iaqsln,iindx1,ipndx1,kbt,kdim,kelect,khydr,khydx,kmax,km1,kmt,ko2gaq,kwater,kx1,kxt,narn1,narn2,nbasp,nbt,nbti,nbtmax,ncmpr,nelect,nern1,nern2,nhydr,nhydx,nobswt,noutpt,no2gaq,nphasx,npt,nptmax,nstmax,nttyo,qloffg,ubmtbi,ufixf,uobsw,uphase,uspec,uzveci,uzvec1,zvclgi,zvclg1,zvec1)

    krdxsp = ko2gaq
    nrdxsp = no2gaq

    if (krdxsp .eq. 0) then
        krdxsp = kelect
        nrdxsp = nelect
    end if

    if (nsrt .gt. 0) then
        ! Construct reactions for special reactants in cases for
        ! which a reaction was not read from the input file.
        call makrsr(cbsri,cesri,cess,eps100,ibsrti,iesrti,jcode,nbt1mx,nct,nctmax,ness,nessmx,nessr,noutpt,nrct,nrctmx,nsrtmx,nstmax,nttyo,ubsri,uelem,uesri,ureac,uspec,zchar)
    end if

    ! Interpret reactant names and associated properties read from the
    ! input file.
    call intrct(cbsr,cbsri,cesr,cesri,egers,egersi,ibsrti,iern1,iern2,iesrti,ietmax,igerti,iktmax,imrn1,imrn2,ixrn1,ixrn2,ixrti,jcode,jetmax,jgerti,jgext,narn1,narn2,nbaspd,nbt,nbtmax,nbt1mx,ncmpr,nct,nctmax,nertmx,netmax,ngexsa,ngext,ngrn1,ngrn2,noutpt,nptmax,nrct,nrctmx,nrndex,nsrtmx,nstmax,nttyo,nxridx,nxrtmx,rxbar,rxbari,ubsri,ucxri,uelem,uesri,ugerji,ugermo,ugersi,ugexj,ugexmo,uphase,ureac,uspec,xgers,xgersi)

    if (nert .gt. 0) then
        ! Calculate the mole ratios (mrgers) for species on exchange
        ! sites of generic ion exchanger reactants. The mole ratio is the
        ! moles of exchanger species (specific to a site) per mole of
        ! exchanger substrate. These mole ratios are used to compute mass
        ! balance increments. They are calculated from the corresponding
        ! equivalent fractions (egers). Note that these mole ratios and
        ! equivalent fractions are model-independent, unlike the
        ! corresponding mole fractions.
        do nrc = 1,nrct
            if (jcode(nrc) .eq. 5) then
                ner = nxridx(nrc)
                np = nrndex(nrc)
                ne = np - iern1 + 1

                do je = 1,jgext(ne)
                    ex = zgexj(je,ne)*cgexj(je,ne)

                    do ie = 1,ngext(je,ne)
                        nss = ngexsa(ie,je,ne)

                        if (nss .gt. 0) then
                            mrgers(ie,je,ner) = -ex*egers(ie,je,ner)/zchar(nss)
                        end if
                    end do
                end do
            end if
        end do

        ! Calculate the mole fractions (xgers) for species on exchange
        ! sites of generic ion exchanger reactants.
        do nrc = 1,nrct
            if (jcode(nrc) .eq. 5) then
                ner = nxridx(nrc)
                np = nrndex(nrc)
                ne = np - iern1 + 1

                do je = 1,jgext(ne)
                    xx = 0.

                    do ie = 1,ngext(je,ne)
                        xx = xx  + mrgers(ie,je,ner)
                    end do

                    do ie = 1,ngext(je,ne)
                        xgers(ie,je,ner) = mrgers(ie,je,ner)/xx
                    end do
                end do
            end if
        end do
    end if

    if (nsrt .gt. 0) then
        ! Make sure that reactions for special reactants are charge
        ! balanced.
        call chzrsr(cbsr,elecsr,eps100,jcode,nbasp,nbt,nbtmax,nbt1mx,noutpt,nrct,nrctmx,nrndex,nsrtmx,nstmax,nttyo,ureac,uspec,zchar)
    end if

    ! Find any basis species in the matrix with jflag = 30 and
    ! eliminate them from the matrix. Adjust as necessary the total
    ! masses of those basis species which remain in the matrix.
    call combmb(cdrs,iindx1,ipndx1,jflag,kbt,kdim,km1,kmax,kmt,kx1,kxt,mtb,mtbaq,ndrsmx,nbasp,nbt,nbtmax,ndrs,ndrsr,noutpt,nstmax,nttyo,uspec,uzvec1,zvclg1,zvec1)

    ! Initialize the error counter.
    nerr = 0

    ! Make sure that no mineral species in the equilibrium system (ES)
    ! has been suppressed.
    do kcol = km1,kxt
        np = ipndx1(kcol)
        ns = iindx1(kcol)

        if (jsflag(ns) .gt. 0) then
            ! Calling sequence substitutions:
            !   uspec(ns) for unam48
            call fmspnm(jlen,uspec(ns),uspn56)
            write (noutpt,1350) uspn56(1:jlen)
            write (nttyo,1350) uspn56(1:jlen)
1350 format(/' * Error - (EQ6/eq6) The species ',a," can't be",/7x,'suppressed because it is present in the equilibrium',' system.')

            nerr = nerr + 1
        else
            jpflag(np) = -1
            jsflag(ns) = -1
        end if
    end do

    ! Check to see that no species in the strict basis has a jflag
    ! value of 30 at this point.
    do nb = 1,nbt
        ns = nbasp(nb)
        nt = ndrsr(2,ns) - ndrsr(1,ns) + 1

        if (nt .lt. 2) then
            if (jflag(ns) .eq. 30) then
                ux = uspec(ns)(1:24)
                j2 = ilnobl(ux)
                write (ux8,'(i5)') jflag(ns)
                call lejust(ux8)
                j3 = ilnobl(ux8)
                write (noutpt,1360) uspec(ns)(1:j2),ux8(1:j3)
                write (nttyo,1360) uspec(ns)(1:j2),ux8(1:j3)
1360 format(/' * Error - (EQ6/eq6) The species ',a," can't have",/7x,'a jflag value of ',a,', because this species is in',' the strict basis set.')

                nerr = nerr + 1
            end if
        end if
    end do

    if (nerr .gt. 0) then
        stop
    end if

    ! Redefine the 'd' set of reactions and attendant data to
    ! match the ordinary set as it presently exists. The ordinary
    ! set will be futher subjected to eliminations from the active
    ! basis set and ordinary basis switching.
    call cdrssd(adhfs,adhfsd,advfs,advfsd,axhfs,axhfsd,axlks,axlksd,axvfs,axvfsd,cdrs,cdrsd,ipch,ipchmx,ipcv,ipcvmx,narxmx,nbasp,nbaspd,nbtmax,ndrs,ndrsd,ndrsmx,ndrsr,ndrsrd,nstmax,ntprmx)
    nbtd = nbt

    ! Compute affinity scaling factors. Note that these are
    ! defined in terms of the reactions before any eliminations
    ! or ordinary basis switches are made.
    call gafscl(cdrsd,cscale,ndrsmx,ndrsrd,nst,nstmax)

    if (iopg(1) .eq. 1) then
        ! Build the S-lambda index arrays nsxi and nsxx.
        call bdslx(narn1,narn2,natmax,noutpt,nslt,nsltmx,nslx,nsxi,nsxx,nsxmax,nttyo)

        ! Build the mu index arrays nmxi and nmxx.
        call bdmlx(narn1,narn2,natmax,nmut,nmutmx,nmux,nmxi,nmxmax,nmxx,noutpt,nttyo)

        if (iopr(10) .gt. 0) then
            ! Write tables concerning Pitzer interaction coefficients.
            call ptztab(iopr,narn1,narn2,natmax,nmutmx,nmux,nmxi,nmxmax,nmxx,noprmx,noutpt,nsltmx,nslx,nstmax,nsxi,nsxmax,nsxx,uspec)
        end if

        ! Print warnings for species lacking Pitzer coefficients.
        call ptzchk(narn1,narn2,natmax,nmxi,noutpt,nstmax,nsxi,nttyo,uspec)

        !        Transform conventional mu data to corresponding C, psi,
        !        and zeta data (data originally defined in mu form is
        !        not affected).
        ! xxxxxxxxxxx
        if (qhawep) then
            call rc3ocf(amu,jpfcmx,ifcphi1,ifcphi2,ifnnn,ifn2n,ifpsi1,ifpsi2,ifzeta,ilcphi1,ilcphi2,ilnnn,iln2n,ilpsi1,ilpsi2,ilzeta,iodb,nmux,nmut,nmutmx,nodbmx,noutpt,nstmax,nttyo,uspec,zchar)
        end if
    end if

    ! Compute the thermodynamic parameters that are functions of
    ! temperature.
    call evdata(aadh,aadhh,aadhv,aaphi,abdh,abdhh,abdhv,abdot,abdoth,abdotv,adadhh,adadhv,adbdhh,adbdhv,adbdth,adbdtv,adh,adhfe,adhh,adhv,adhfs,adhfsd,advfe,advfs,advfsd,afcnst,al10,amu,aslm,aphi,aprehw,apresg,apresh,apx,avcnst,axhfe,axhfs,axhfsd,axlke,axlks,axlksd,axvfe,axvfs,axvfsd,bdh,bdhh,bdhv,bdot,bdoth,bdotv,dadhh,dadhv,dbdhh,dbdhv,dbdth,dbdtv,dhfe,dhfs,dhfsd,dvfe,dvfs,dvfsd,ehfac,farad,iapxmx,iktmax,iopg,iopt,ipbtmx,ipch,ipchmx,ipcv,ipcvmx,ixrn1,ixrn2,jpfcmx,jptffl,jsol,narxmx,narxt,narxth,ncmpr,nmut,nmutmx,nopgmx,noptmx,noutpt,nptmax,nslt,nsltmx,nst,nstmax,ntpr,ntprmx,nttyo,nxt,nxtmax,pmu,prehw,presg,presh,press,pslamn,rconst,rcnstv,rtcnst,tempc,tempk,uphase,uspec,wfac,xhfe,xhfs,xhfsd,xlke,xlks,xlksd,xvfe,xvfs,xvfsd)

    ! Compute the rate parameters that are functions of temperature.
    call evratc(eact,hact,iact,imchmx,imech,nrct,nrctmx,nrk,rk,rkb,rtcnst,tempk,trkb)

    tempcd = tempc

    ! Check for thermodynamic pressure corrections.
    dp = press - presg

    if (abs(dp) .ge. 1.e-4) then
        if (ipcv .lt. 0) then
            ! There are no data to support needed pressure corrections.
            write (noutpt,1370) press,presg,dp
            write (nttyo,1370) press,presg,dp
1370 format(/' * Warning - (EQ6/eq6) The supporting data file',/7x,'contains no data to support making thermodynamic',/7x,'pressure corrections. No such corrections will be made.',/7x,'The current pressure is ',1pg12.5,' bars, the standard',/7x,'grid pressure is ',g12.5,' bars, and the pressure',/7x,'difference is ',g12.5,' bars.')
        else
            ! Make pressure corrections to the thermodynamic data.
            call pcorrm(adh,adhh,adhv,al10,aphi,avcnst,bdh,bdhh,bdhv,bdot,bdoth,bdotv,dadhh,dadhv,dbdhh,dbdhv,dbdth,dbdtv,dhfe,dvfe,iopg,ipch,ipchmx,ipcv,ipcvmx,nopgmx,presg,press,rcnstv,tempk,xhfe,xlke,xvfe)

            call pcorrx(avcnst,dhfs,dvfs,ipch,ipchmx,ipcv,ipcvmx,nbasp,nbt,nbtmax,ndrsr,nst,nstmax,presg,press,xhfs,xlks,xvfs)

            ! Calling sequence substitutions:
            !   dhfsd for dhfs
            !   dvfsd for dvfs
            !   nbaspd for nbasp
            !   nbtd for nbt
            !   ndrsrd for ndrsr
            !   xhfsd for xhfs
            !   xlksd for xlks
            !   xvfsd for xvfs
            call pcorrx(avcnst,dhfsd,dvfsd,ipch,ipchmx,ipcv,ipcvmx,nbaspd,nbtd,nbtmax,ndrsrd,nst,nstmax,presg,press,xhfsd,xlksd,xvfsd)
        end if
    end if

    pressd = press

    write (noutpt,1400)
1400 format(//21x,'--- Inactive Species ---',/)

    k = 0
    klim = 20

    do ns = 1,nst
        if (jsflag(ns) .eq. 1) then
            if (ns .ne. no2gaq) then
                k = k + 1

                if (k .le. klim) then
                    ! Calling sequence substitutions:
                    !   uspec(ns) for unam48
                    call fmspnm(jlen,uspec(ns),uspn56)
                    write (noutpt,1410) uspn56(1:jlen)
1410 format(4x,a)
                end if
            end if
        end if
    end do

    if (k .gt. klim) then
        i = k - klim
        write (noutpt,1420) i
1420 format(/6x,'plus ',i4,' others',/)
    else if (k .le. 0) then
        write (noutpt,1430)
1430 format(4x,'None',/)
    else
        write (noutpt,1440)
1440 format(1x)
    end if

    ! Call EQLIB/elim.f to rewrite the reaction equations (cdrs/ndrs/
    ! ndrsr arrays) so that auxiliary basis variables with jflag = 30
    ! are eliminated from the active basis set.
    qelim = .false.

    do nb = 1,nbt
        nse = nbasp(nb)
        nt = ndrsr(2,nse) - ndrsr(1,nse) + 1

        if (nt.ge.2 .and. jflag(nse).eq.30) then
            if (jsflag(nse) .gt. 0) then
                mtb(nb) = 0.
                mtbaq(nb) = 0.
            else
                qelim = .true.
                call elim(adhfs,adhfsx,advfs,advfsx,axhfs,axhfsx,axlks,axlksx,axvfs,axvfsx,cdrs,cdrsx,eps100,ipch,ipchmx,ipcv,ipcvmx,jsflag,narxmx,ndrs,ndrsmx,ndrsx,ndrsr,ndrsrx,nse,nst,nstmax,ntprmx,noutpt,nttyo,uspec)
            end if
        end if
    end do

    ! Compute mass balance coefficients using the active part of the
    ! original basis set.
    call gcsts(cdrs,csts,jflag,nbaspd,nbt,nbtmax,ndrs,ndrsmx,ndrsr,noutpt,nsts,nstsmx,nstsr,nst,nstmax,nttyo,uspec)

    ! Copy the existing nbasp array into the nbaspi array. The
    ! main purpose of this is to recall what the basis set was
    ! when mass balance relationships were defined. In particular,
    ! this array holds the basis set as it was prior to the switching
    ! out of bare exchangerspecies. This information will be used
    ! in describing mass balance totals. In particular, it will be
    ! needed to write the pickup file.
    do nb = 1,nbt
        nbaspi(nb) = nbasp(nb)
    end do

    if (net .gt. 0) then
        ! Switch bare site species of generic ion exchange phases with
        ! certain models such as Gapon and Vanselow out of the basis set.
        ! Then suppress these bare site species so that they appear in the
        ! model only in association with corresponding mass balances.
        ! Changes are made here to the ion exchanger sections of both the
        ! ordinary and the 'd' sets of reactions and associated data.
        call chsgex(adhfs,adhfsd,adhfsx,advfs,advfsd,advfsx,axhfs,axhfsd,axhfsx,axlks,axlksd,axlksx,axvfs,axvfsd,axvfsx,cdrs,cdrsd,cdrsx,eps100,iern1,ipch,ipchmx,ipcv,ipcvmx,jern1,jetmax,jflag,jgext,jsflag,narn1,narxmx,narxt,nbasp,nbaspd,nbaspx,nbt,nbtmax,nbw,ndrs,ndrsd,ndrsmx,ndrsr,ndrsrd,ndrsrx,ndrsx,nern1,nern2,net,netmax,ngext,noutpt,nphasx,nst,nstmax,ntprmx,ntprt,nttyo,qbassw,qbswok,ugexmo,uspec)

        do kcol = 1,kbt
            nb = iindx1(kcol)
            ns = nbasp(nb)

            if (ns.ge.nern1 .and. ns.le.nern2) then
                uzvec1(kcol) = uspec(ns)
            end if
        end do
    end if

    if (nsrt .gt. 0) then
        ! Check the reactions for the special reactants. Rewrite these
        ! reactions to eliminate any basis species for which jflag = 30.
        call ckfrsr(cbsr,csts,jcode,jflag,nbaspd,nbtd,nbtmax,nbt1mx,noutpt,nrct,nrctmx,nrndex,nsrtmx,nstmax,nsts,nstsmx,nstsr,nttyo,ureac,uspec)
    end if

    if (iopt(5) .gt. 0) then
        ! Clear the ES solids read from the input file. Fictive
        ! fugacity-fixing minerals are not cleared.
        write (noutpt,1470)
        write (nttyo,1470)
1470 format(/' * Note (EQ6/eq6) Clearing equilibrium system',' (ES) solids',/7x,'read from the input file.')

        call clress(csts,iindx1,ipndx1,jpflag,jsflag,kdim,kmax,km1,kmt,kx1,kxt,loph,losp,moph,mosp,mtb,mtbaq,nbt,nbtmax,nptmax,nstmax,nsts,nstsmx,nstsr,ufixf,uzvec1,zvec1,zvclg1)
    end if

    ! Save the present configuration of the jflag array. This will
    ! be used to support calculation of saturation indices and
    ! affinities for the 'd' set of reactions.
    call copyia(jflag,jflagd,nbtmax)

    ! Execute any ordinary basis switching directives from the input
    ! file.
    kobswt = 0

    if (nobswt .gt. 0) then
        ! Copy the existing nbasp set into the nbaspx array.
        call copyia(nbasp,nbaspx,nbt)

        ! Interpret the switches.
        call intbsw(nbasp,nbaspx,nbt,nbtmax,nobswt,noutpt,nst,nstmax,nttyo,uobsw,uspec)

        ! Execute the switches.
        do nb = 1,nbt
            ns1 = nbaspx(nb)
            ns2 = nbasp(nb)

            if (ns1 .ne. ns2) then
                call switch(adhfs,adhfsx,advfs,advfsx,axhfs,axhfsx,axlks,axlksx,axvfs,axvfsx,cdrs,cdrsx,eps100,ipch,ipchmx,ipcv,ipcvmx,jflag,jsflag,narn1,narxmx,nbasp,nbaspd,nbaspx,nb,nbt,nbtmax,nbw,ndrs,ndrsmx,ndrsx,ndrsr,ndrsrx,noutpt,ns2,nst,nstmax,ntprmx,nttyo,qbassw,qbswok,uspec)
                kobswt = kobswt + 1
            end if
        end do
    end if

    ! Put new fixed fugacity minerals into the matrix if the
    ! corresponding added masses are positive. Keep old ones
    ! in if the corresponding net masses are positive. If necessary,
    ! recalculate the mass balance totals. Note: this call is
    ! necessary even if nffg = 0, because it completes the
    ! elimination of any fictive fugacity-fixing phases left over
    ! from a previous run.
    call setffg(csts,iindx1,iffg,ipndx1,jpflag,jsflag,kbt,kdim,kmax,km1,kmt,kx1,kxt,losp,moffg,mtb,mtbaq,nbaspd,nbt,nbtmax,ncmpr,nffg,nffgmx,noutpt,nphasx,npt,nptmax,nstmax,nsts,nstsmx,nstsr,nttyo,qloffg,uffg,uspec,uzvec1,zvclg1,zvec1)

    ! Recompute equilibrium constants, etc., for the ordinary set
    ! if the reactions have been modified by elimination of
    ! auxiliary basis species or by basis switching. The presence
    ! of generic ion exchangers (net > 0) implies intrinsic or
    ! hidden basis switching.
    if (qelim .or. kobswt.gt.0 .or. net.gt.0) then
        call evdatr(adhfs,advfs,axhfs,axlks,axvfs,dhfs,dvfs,ipch,ipchmx,ipcv,ipcvmx,narxmx,narxt,nst,nstmax,ntpr,ntprmx,tempc,xhfs,xlks,xvfs)

        if (ipcv .ge. 0) then
            call pcorrx(avcnst,dhfs,dvfs,ipch,ipchmx,ipcv,ipcvmx,nbasp,nbt,nbtmax,ndrsr,nst,nstmax,presg,press,xhfs,xlks,xvfs)
        end if
    end if

    kobswt = 0

    ! Recompute equilibrium constants, etc., for the 'd' set
    ! if generic ion exchangers are present. If so, the subset
    ! of the 'd' set reactions dealing with such exchangers
    ! was modified above by the call to EQLIB/chsgex.f.
    if (net .gt. 0) then
        ! Calling sequence substitutions:
        !   adhfsd for adhfs
        !   advfsd for advfs
        !   axhfsd for axhfs
        !   axlksd for axlks
        !   axvfsd for avhfs
        !   dhfsd for dhfs
        !   dvfsd for dvfs
        !   xhfsd for xhfs
        !   xlksd for xlks
        !   xvfsd for xvfs
        call evdatr(adhfsd,advfsd,axhfsd,axlksd,axvfsd,dhfsd,dvfsd,ipch,ipchmx,ipcv,ipcvmx,narxmx,narxt,nst,nstmax,ntpr,ntprmx,tempc,xhfsd,xlksd,xvfsd)

        ! Calling sequence substitutions:
        !   dhfsd for dhfs
        !   dvfsd for dvfs
        !   nbaspd for nbasp
        !   nbtd for nbt
        !   ndrsrd for ndrsr
        !   xhfsd for xhfs
        !   xlksd for xlks
        !   xvfsd for xvfs
        if (ipcv .ge. 0) then
            call pcorrx(avcnst,dhfsd,dvfsd,ipch,ipchmx,ipcv,ipcvmx,nbaspd,nbtd,nbtmax,ndrsrd,nst,nstmax,presg,press,xhfsd,xlksd,xvfsd)
        end if
    end if

    if (iopr(2) .ge. 1) then
        ! Calling sequence substitutions:
        !   noutpt for nf
        ilevel = iopr(2)
        call echolk(axlks,cdrs,ilevel,jsflag,narxmx,ndrs,ndrsmx,ndrsr,noutpt,nst,ntprmx,nstmax,press,tempc,uspec,xlks)
    end if

    if (nrct .le. 0) then
        go to 350
    end if

    ! Process the input reactant data.
    ! Get local indices of reactants, check their names, and set up
    ! their molecular weights and molar volumes.
    call rsetup(atwt,cbsr,cesr,iern1,ietmax,iindx1,iktmax,jcode,jern1,jern2,jetmax,jgext,kbt,kmax,mwtges,mwtrc,mwtsp,nbaspd,nbt,nbtmax,nbt1mx,ncmpr,nct,nctmax,nertmx,netmax,ngext,noutpt,nptmax,nrct,nrctmx,nrndex,nsrtmx,nstmax,nsts,nstsmx,nstsr,nttyo,nxridx,nxrtmx,rxbar,ureac,uspec,vosp0,vreac,xgers)

    ! Calculate surface area parameters.
    do nrc = 1,nrct
        gx = morr(nrc)*mwtrc(nrc)

        if (nsk(nrc).eq.0 .or. nsk(nrc).eq.2) then
            ! Calculate the surface area (cm2) from the specific surface
            ! area (cm2/g).
            ssfcar(nrc) = 0.

            if (gx .gt. 0.) then
                ssfcar(nrc) = sfcar(nrc)/gx
            end if
        else if (nsk(nrc) .eq. 1) then
            ! Calculate the specific surface area (cm2/g) from the
            ! surface area (cm2).
            sfcar(nrc) = ssfcar(nrc)*gx
        else
            j2 = ilnobl(ureac(nrc))
            write (noutpt,3110) nsk(nrc),ureac(nrc)(1:j2)
            write (nttyo,3110) nsk(nrc),ureac(nrc)(1:j2)
3110 format(/' * Error - (EQ6/eq6) The surface area code nsk',/7x,'has an unrecognized value of ',i3,' for reactant',/7x,a,'.')

            nerr = nerr + 1
        end if
    end do

    ! Set backward/precipitation kinetics flag array, npchk.
    do nrc = 1,nrct
        if (nrk(2,nrc) .ne. 0) then
            if (jcode(nrc).eq.0 .or. jcode(nrc).eq.1) then
                do np = 1,npt
                    if (ureac(nrc) .eq. uphase(np)) then
                        npchk(np) = 1
                        go to 300
                    end if
                end do
            end if
        end if

300 continue
    end do

    ! Check names for activity terms in rate equations.
    do nrc = 1,nrct
        if (nrk(1,nrc) .eq. 2) then
            do i = 1,imech(1,nrc)
                if (ndact(i,1,nrc) .gt. 0) then
                    do n = 1,ndact(i,1,nrc)
                        unamsp = udac(n,i,1,nrc)

                        do ns = 1,nst
                            if (unamsp(1:24) .eq. uspec(ns)(1:24)) then
                                go to 310
                            end if
                        end do

                        j3 = ilnobl(unamsp)
                        j2 = ilnobl(ureac(nrc))
                        write (noutpt,1500) unamsp(1:j3),ureac(nrc)(1:j2)
                        write (nttyo,1500) unamsp(1:j3),ureac(nrc)(1:j2)
1500 format(/" * Error - (EQ6/eq6) Can't match the species",/7x,'name "',a,'" appearing in the rate law for',/7x,'with any species present in the model system.')

                        nerr = nerr + 1
310 continue
                        ndac(n,i,1,nrc) = ns
                    end do
                end if
            end do
        end if

        if (nrk(2,nrc) .eq. 2) then
            do i = 1,imech(2,nrc)
                if (ndact(i,2,nrc) .gt. 0) then
                    do n = 1,ndact(i,2,nrc)
                        unamsp = udac(n,i,2,nrc)

                        do ns = 1,nst
                            if (unamsp(1:24) .eq. uspec(ns)(1:24)) then
                                go to 320
                            end if
                        end do

                        j3 = ilnobl(unamsp)
                        j2 = ilnobl(ureac(nrc))
                        write (noutpt,1500) unamsp(1:j3),ureac(nrc)(1:j2)
                        write (noutpt,1500) unamsp(1:j3),ureac(nrc)(1:j2)
                        nerr = nerr + 1
320 continue
                        ndac(n,i,2,nrc) = ns
                    end do
                end if
            end do
        end if
    end do

    if (nerr .gt. 0) then
        stop
    end if

350 continue

    ! Initialize non-aqueous phases in the physically removed system
    ! (PRS).
    if (nprpti .gt. 0) then
        nerr = 0
        npi = 0
        unamph = ' '

        do nsi = 1,nprsti
            unamsp(1:24) = uprspi(nsi)(1:24)

            if (unamph(1:24) .ne. uprspi(nsi)(25:48)) then
                ! Have found another phase.
                npi = npi + 1
                unamph(1:24) = uprphi(npi)(1:24)

                do np = 1,npt
                    if (unamph(1:24) .eq. uphase(np)(1:24)) then
                        go to 370
                    end if
                end do

                j2 = ilnobl(unamph)
                write (noutpt,1530) unamph(1:j2)
                write (nttyo,1530) unamph(1:j2)
1530 format(/' * Error - (EQ6/eq6) The phase ',a,' is listed',/7x,'as present in the physically removed system (PRS)',' on the input file,',/7x,"but it isn't in the data file",' set (after compression).',/7x,'You may be using the',' wrong data file.')

                nerr = nerr + 1
                go to 380

370 continue
                mprph(np) = mprphi(npi)
                nr1 = ncmpr(1,np)
                nr2 = ncmpr(2,np)

                do ns = nr1,nr2
                    if (unamsp(1:24) .eq. uspec(ns)(1:24)) then
                        mprsp(ns) = mprspi(nsi)
                        go to 380
                    end if
                end do

                j2 = ilnobl(unamph)
                j3 = ilnobl(unamsp)
                write (noutpt,1540) unamsp(1:j3),unamph(1:j2)
                write (nttyo,1540) unamsp(1:j3),unamph(1:j2)
1540 format(/' * Error - (EQ6/eq6) The species ',a,' of phase',/7x,a,' is listed as present in the physically removed',/7x,'system (PRS) on the input file, but this species'," isn't in the",/7x,'data file set (after compression).',' You may be using the wrong',/7x,'data file.')

                nerr = nerr + 1
            else
                ! Have found another species belonging to the same phase.
                do ns = nr1,nr2
                    if (unamsp(1:24) .eq. uspec(ns)(1:24)) then
                        mprsp(ns) = mprspi(nsi)
                        go to 380
                    end if
                end do

                j2 = ilnobl(unamph)
                j3 = ilnobl(unamsp)
                write (noutpt,1540) unamsp(1:j3),unamph(1:j2)
                write (nttyo,1540) unamsp(1:j3),unamph(1:j2)
                nerr = nerr + 1
            end if

380 continue
        end do

        if (nerr .gt. 0) then
            stop
        end if
    end if

    if (iopt(9) .gt. 0) then
        ! Clear the PRS solids read from the input file.
        write (noutpt,1570)
        write (nttyo,1570)
1570 format(/' * Note (EQ6/eq6) Clearing physically removed',' system (PRS) solids',/7x,'read from the input file.')

        call initaz(mprsp,nstmax)
        call initaz(mprph,nptmax)
    end if

    ! Write an echo of the input problem, including defaults, on
    ! the output file.
    call echoz(axlks,awmaxi,awmini,azero,cbsr,cdac,cdrs,cesr,csigma,dlaplo,dlaprn,dleplo,dleprn,dlhplo,dlhprn,dloplo,dloprn,dltpll,dltplo,dltprl,dltprn,dlxdmp,dlxmax,dlxmx0,dlxpll,dlxplo,dlxprl,dlxprn,eact,ehmaxi,ehmini,iktmax,imchmx,imech,iodb,iopg,iopr,iopt,itermx,jcode,jpress,jtemp,jsflag,ksplmx,ksppmx,kstpmx,kxmod,mwtrc,narn1,narn2,narxmx,nat,nata,natmax,nbaspd,nbt,nbta,nbtd,nbtmax,nbt1mx,ncmpr,nct,ncta,nctmax,ndact,ndctmx,ndrs,ndrsmx,ndrsr,nffgmx,ngt,ngta,ngtmax,nlt,nlta,nltmax,nmt,nmta,nmtmax,nodbmx,nopgmx,noprmx,noptmx,nordmx,noutpt,npslmx,npt,npta,nptkmx,nptmax,nrct,nrctmx,nrk,nrndex,nsk,nsrt,nsrtmx,nsscmx,nsslmx,nst,nsta,nstmax,ntprmx,ntrymx,nttkmx,nttyo,nxmdmx,nxmod,nxopmx,nxopex,nxopt,nxpemx,nxridx,nxrt,nxrtmx,nxt,nxta,nxtmax,o2maxi,o2mini,phmaxi,phmini,press,pressb,ptk,qredox,rkb,rxbar,sscrew,tempc,tempcb,tempk,timmxi,tistti,tolbt,toldl,tolsat,tolsst,tolxsf,tolxst,tolxsu,trkb,ttk,uactop,udac,uelem,uffg,ureac,uspec,uxcat,uxmod,uxopex,uxopt,vreac,xistti,ximaxi,xlkmod,xlks,zkfac,zklgmn,zklogl,zklogu)

    ! Roll over the main title to the secondary title.
    ntitl2 = ntitl1

    do n = 1,ntitl1
        utitl2(n) = utitl1(n)
    end do

    if (utitl1(1)(7:37) .ne. 'This main title is a carry-over') then
        if (ntitl1 .lt. ntitmx) then
            do n = ntitl1,1,-1
                utitl1(n + 1) = utitl1(n)
            end do

            utitl1(1) =    'Note: This main title is a carry-over from a previous run.'
        end if
    end if

    if (qcwrpj) then
        ! Write tables of the Pitzer J(x) and J'(x) functions.
        call cwrpjt(noutpt)
    end if

    ! Trace the reaction path.
    call path(aadh,aadhh,aadhv,aaphi,abdh,abdhh,abdhv,abdot,abdoth,abdotv,adadhh,adadhv,adbdhh,adbdhv,adbdth,adbdtv,adh,adhh,adhv,afcnst,aftarg,al10,aphi,apx,atwt,avcnst,awmax,awmaxi,awmin,awmini,azero,bdh,bdhh,bdhv,bdot,bdoth,bdotv,bpx,cbsr,cbsri,cco2,cdac,cdrs,cdrsd,cdrsx,cegexs,cesr,cesri,cess,cpgexs,cscale,csigma,csts,dadhh,dadhv,dbdhh,dbdhv,dbdth,dbdtv,dlaplo,dlaprn,dleplo,dleprn,dlhplo,dlhprn,dloplo,dloprn,dltplo,dltpll,dltprl,dltprn,dlxdmp,dlxmax,dlxmin,dlxmx0,dlxplo,dlxpll,dlxprl,dlxprn,eact,egers,egersi,egexjf,ehfac,ehmax,ehmaxi,ehmin,ehmini,elecsr,electr,eps100,farad,fkrc,hact,iact,iapxt,iaqsln,ibpxt,ibsrti,ielam,iern1,iern2,iesrti,ifcphi1,ifcphi2,ifnnn,ifn2n,ifpsi1,ifpsi2,ifrn1,ifrn2,ifzeta,igas,igerti,iindx1,ilcphi1,ilcphi2,ilnnn,iln2n,ilpsi1,ilpsi2,ilrn1,ilrn2,ilzeta,imech,imrn1,imrn2,insgf,iodb,iopg,iopr,iopt,ipch,ipndx1,ipcv,irang,itermx,ixrn1,ixrn2,ixrti,izmax,jcode,jffg,jflag,jflagd,jflgi,jgerti,jpflag,jpress,jptffl,jreac,jsflag,jsitex,jsol,jtemp,kbt,kct,kdim,kelect,khydr,khydx,km1,kmt,ko2gaq,kprs,krdxsp,ksplmx,ksppmx,kstpmx,kwater,kxmod,kx1,kxt,loph,losp,modr,moffg,moph,morr,mosp,mprph,mprphi,mprsp,mprspi,mrgers,mrgexs,mtb,mtbi,mtbaq,mtbaqi,mte,mteaq,mwtrc,mwtsp,narn1,narn2,narxt,nat,nbasp,nbaspd,nbaspi,nbaspx,nbkupa,nbkupb,nbt,nbtd,nbti,nbw,nchlor,ncmpr,nct,ndac,ndact,nelect,nern1,nern2,ness,nessr,net,ndrs,ndrsd,ndrsx,ndrsr,ndrsrd,ndrsrx,nert,newin,nffg,nfrn1,nfrn2,ngrn1,ngrn2,ngt,nhydr,nhydx,nllnmx,nlrn1,nlrn2,nlt,nmrn1,nmrn2,nmrt,nmt,nobswt,noutpt,no2gaq,npchk,nphasx,nprob,nprpti,nprsti,npslmx,npt,nrct,nrdxsp,nrk,nrndex,nsbswt,nsk,nsrt,nsslmx,nst,nsts,nstsr,ntabx,ntf1,ntf1t,ntf2,ntf2t,ntitl1,ntitl2,ntitld,ntpr,ntprt,ntrymx,nttyo,nxmod,nxopex,nxopt,nxridx,nxrn1,nxrn2,nxrt,nxt,o2max,o2maxi,o2min,o2mini,phmax,phmaxi,phmin,phmini,prcinf,press,pressb,pressd,pressi,ptk,qcnpre,qcntmp,qdwipp,qecon,qgexsh,qhawep,qoptmz,qpit75,qredox,qscon,qtatxt,rconst,rcnstv,rk,rkb,rtcnst,rxbar,rxbari,sfcar,smp100,sscrew,ssfcar,tempc,tempcb,tempcd,tempci,tempcu,tempk,tf1,tf2,timemx,time1,timmxi,tistrt,tistti,trkb,ttk,tolaft,tolbt,toldl,tolsat,tolsst,tolxsf,tolxst,tolxsu,uaqsln,ubmtbi,ubsri,ucxri,udac,uelem,uesri,uffg,ufixf,ugerji,ugermo,ugersi,uinfor,ulinex,uobsw,uphase,uplatm,uprphi,uprspi,ureac,usbsw,uspec,usteq6,utitl1,utitl2,utitld,uveeq6,uxcat,uxmod,uxopex,uxopt,uzvec1,uzveci,vosp0,vreac,wfac,xgers,xgersi,ximax,ximaxi,xistrt,xistti,xi1,xlkffg,xlkmod,zchar,zchcu6,zchsq2,zklgmn,zklogl,zklogu,zvclgi,zvclg1,zvec1)

    ! Descramble the tabx file onto the scratch tab file, then copy
    ! the scratch tab file onto the tab file.
    if (iopt(18) .ge. 0) then
        if (qtatxt) then
            ! The TAB file is an ordinary text file.
            ! Calling sequence substitutions:
            !   ntabx for nf1
            !   ntabs for nf2
            call dscrax(ntabx,ntabs,nllnmx,ulinex)

            ! Calling sequence substitutions:
            !   ntabs for nf1
            !   ntab for nf2
            call fcopyx(ntabs,ntab,nllnmx,ulinex)
        else
            ! The TAB file is a .csv file.
            call dscramc(ntabx,ntabs,nllnmx,ulinex)

            ! Calling sequence substitutions:
            !   ntabs for nf1
            !   ntab for nf2
            call fcopyx(ntabs,ntab,nllnmx,ulinex)
        end if
    end if

    go to 20
end program eq6