program eq3nr
    !! This is the main program of the EQ3NR code. Configuration
    !! identification, the copyright statement, legal disclaimers, and
    !! similar statements are contained in EQ3NR/aaaeq3.f, the lead-off
    !! subroutine in the EQ3NR source code. A short description of this
    !! program is also contained in that subroutine.
    implicit none

    include 'eqlib/eqldef.h'
    include 'eqlib/eqlpar.h'
    include 'eqlib/eqldv.h'
    include 'eqlib/eqlge.h'
    include 'eqlib/eql1s.h'
    include 'eqlib/eqlwd.h'

    include 'eqlib/eqlj8.h'
    include 'eqlib/eqlk8.h'
    include 'eqlib/eqlo8.h'

    ! File path parameters
    integer :: numargs
    character(len=1024) :: temppath
    character(len=:), allocatable :: data1path
    character(len=:), allocatable :: threeipath
    integer :: pathindices(2)
    character(len=:), allocatable :: basename
    character(len=:), allocatable :: ofile
    character(len=:), allocatable :: ifile
    character(len=:), allocatable :: pfile

    ! Array allocation size variables used in EQ3NR.
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
    integer :: nxic_asv

    integer :: ntf1_asv
    integer :: ntf2_asv
    integer :: ntfx_asv

    integer :: k_asv

    integer :: nxti_asv

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

    integer :: newin
    integer :: ninpt
    integer :: ninpts
    integer :: noutpt
    integer :: nttyo

    integer :: iodb(nodb_par)
    integer :: iopg(nopg_par)
    integer :: iopr(nopr_par)
    integer :: iopt(nopt_par)

    integer :: jgexti(net_par)
    integer :: ngexpi(net_par)
    integer :: ngexti(jet_par,net_par)

    integer, dimension(:), allocatable :: iapxt
    integer, dimension(:), allocatable :: ibpxt
    integer, dimension(:), allocatable :: ibswx
    integer, dimension(:), allocatable :: iction
    integer, dimension(:), allocatable :: igstak
    integer, dimension(:), allocatable :: iindx1
    integer, dimension(:), allocatable :: insgf
    integer, dimension(:), allocatable :: ipivot
    integer, dimension(:), allocatable :: ipndx1
    integer, dimension(:), allocatable :: istack
    integer, dimension(:), allocatable :: ixbasp
    integer, dimension(:), allocatable :: jcsort
    integer, dimension(:), allocatable :: jflag
    integer, dimension(:), allocatable :: jflagd
    integer, dimension(:), allocatable :: jgsort
    integer, dimension(:), allocatable :: jgstak
    integer, dimension(:), allocatable :: jjndex
    integer, dimension(:), allocatable :: jjsort
    integer, dimension(:), allocatable :: jpflag
    integer, dimension(:), allocatable :: jsflag
    integer, dimension(:), allocatable :: jsitex
    integer, dimension(:), allocatable :: jsol
    integer, dimension(:), allocatable :: jssort
    integer, dimension(:), allocatable :: jstack
    integer, dimension(:), allocatable :: kction
    integer, dimension(:), allocatable :: kkndex
    integer, dimension(:), allocatable :: kxmod

    integer, dimension(:), allocatable :: jflgi

    integer, dimension(:), allocatable :: narxt
    integer, dimension(:), allocatable :: nbasp
    integer, dimension(:), allocatable :: nbaspd
    integer, dimension(:), allocatable :: nbaspx
    integer, dimension(:), allocatable :: nbmap
    integer, dimension(:), allocatable :: ncmap
    integer, dimension(:), allocatable :: ncosp
    integer, dimension(:), allocatable :: ndecsp
    integer, dimension(:), allocatable :: ndrs
    integer, dimension(:), allocatable :: ndrsd
    integer, dimension(:), allocatable :: ndrsx
    integer, dimension(:), allocatable :: ness
    integer, dimension(:), allocatable :: nfac
    integer, dimension(:), allocatable :: nphasx
    integer, dimension(:), allocatable :: npnxp
    integer, dimension(:), allocatable :: nsmap
    integer, dimension(:), allocatable :: nsts
    integer, dimension(:), allocatable :: ntfx
    integer, dimension(:), allocatable :: ntf1
    integer, dimension(:), allocatable :: ntf2

    integer, dimension(:), allocatable :: nbaspi

    integer, dimension(:,:), allocatable :: ncmpr
    integer, dimension(:,:), allocatable :: ndrsr
    integer, dimension(:,:), allocatable :: ndrsrd
    integer, dimension(:,:), allocatable :: ndrsrx
    integer, dimension(:,:), allocatable :: nessr
    integer, dimension(:,:), allocatable :: nstsr

    integer, dimension(:,:), allocatable :: ncmpri

    integer, dimension(:), allocatable :: nalpha

    integer, dimension(:,:), allocatable :: nmux
    integer, dimension(:,:), allocatable :: nmxi
    integer, dimension(:,:), allocatable :: nmxx
    integer, dimension(:,:), allocatable :: nslx
    integer, dimension(:,:), allocatable :: nsxi
    integer, dimension(:,:), allocatable :: nsxx

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

    integer :: i
    integer :: iaqsla
    integer :: iaqsln
    integer :: ibetmx
    integer :: icarb
    integer :: idelmx
    integer :: iebal
    integer :: iebal3
    integer :: ielam
    integer :: ier
    integer :: iern1
    integer :: iern2
    integer :: iexec0
    integer :: ifrn1
    integer :: ifrn2
    integer :: igas
    integer :: ilevel
    integer :: ilrn1
    integer :: ilrn2
    integer :: imrn1
    integer :: imrn2
    integer :: io2gaq
    integer :: ipch
    integer :: ipcv
    integer :: irang
    integer :: irdxc3
    integer :: itdsf3
    integer :: iter
    integer :: itermx
    integer :: itmxsv
    integer :: ixrn1
    integer :: ixrn2
    integer :: ixrn1a
    integer :: ixrn2a
    integer :: izmax
    integer :: j
    integer :: je
    integer :: jexec0
    integer :: jfl
    integer :: jfleba
    integer :: jlen
    integer :: jpdblo
    integer :: jpres3
    integer :: jptffl
    integer :: j1
    integer :: j2
    integer :: j3
    integer :: j4
    integer :: k
    integer :: ka1
    integer :: kat
    integer :: kbt
    integer :: kcarb
    integer :: kcol
    integer :: kct
    integer :: kdim
    integer :: kebal
    integer :: kelect
    integer :: ker
    integer :: ke1
    integer :: ket
    integer :: khydr
    integer :: km1
    integer :: kmt
    integer :: kobswt
    integer :: ko2gaq
    integer :: kprs
    integer :: kwater
    integer :: kx1
    integer :: kxt

    integer :: n
    integer :: nad1
    integer :: napt
    integer :: napta
    integer :: narn1
    integer :: narn2
    integer :: narn1a
    integer :: narn2a
    integer :: nat
    integer :: nata
    integer :: nb
    integer :: nbi
    integer :: nbt
    integer :: nbta
    integer :: nbtafd
    integer :: nbtd
    integer :: nbti
    integer :: nbw
    integer :: nbwa
    integer :: nb1
    integer :: nb2
    integer :: nc
    integer :: ncarb
    integer :: nchloa
    integer :: nchlor
    integer :: nct
    integer :: ncta
    integer :: ne
    integer :: neleca
    integer :: nelect
    integer :: nern1
    integer :: nern2
    integer :: nerr
    integer :: net
    integer :: neti
    integer :: ngrn1
    integer :: ngrn1a
    integer :: ngrn2a
    integer :: ngrn2
    integer :: ngt
    integer :: ngta
    integer :: nhydr
    integer :: nhydra
    integer :: nhydx
    integer :: nhydxa
    integer :: nlrn1
    integer :: nlrn2
    integer :: nlrn1a
    integer :: nlrn2a
    integer :: nlt
    integer :: nlta
    integer :: nmax
    integer :: nmrn1
    integer :: nmrn2
    integer :: nmrn1a
    integer :: nmrn2a
    integer :: nmt
    integer :: nmta
    integer :: nmut
    integer :: nmuta
    integer :: nobswt
    integer :: no2gaa
    integer :: no2gai
    integer :: no2gaq
    integer :: np
    integer :: npt
    integer :: npta
    integer :: nprob
    integer :: nrecl
    integer :: nrdxsa
    integer :: nrdxsp
    integer :: nredox
    integer :: nrr1
    integer :: nrr2
    integer :: nr1
    integer :: nr2
    integer :: ns
    integer :: nsbsw
    integer :: nsbswt
    integer :: nse
    integer :: nslt
    integer :: nslta
    integer :: nss
    integer :: nst
    integer :: nsta
    integer :: ns1
    integer :: ns2
    integer :: nt
    integer :: ntfxt
    integer :: ntf1t
    integer :: ntf1ta
    integer :: ntf2t
    integer :: ntf2ta
    integer :: ntitl
    integer :: ntitl2
    integer :: ntitld
    integer :: ntpr
    integer :: ntprt
    integer :: nxi
    integer :: nxic
    integer :: nxmod
    integer :: nxrn1
    integer :: nxrn2
    integer :: nxrn1a
    integer :: nxrn2a
    integer :: nxt
    integer :: nxta
    integer :: nxti

    integer :: ilnobl
    integer :: nbasis

    logical, dimension(:), allocatable :: qxknph

    logical :: qbassw
    logical :: qchlor
    logical :: qbswok
    logical :: qcwrpj
    logical :: qdwipp
    logical :: qelim
    logical :: qend
    logical :: qhawep
    logical :: qop
    logical :: qpit75
    logical :: qredox
    logical :: qrderr
    logical :: qrho
    logical :: q6mode

    character(len=80), dimension(:), allocatable :: utitl
    character(len=80), dimension(:), allocatable :: utitl2
    character(len=48), dimension(:), allocatable :: ubasp
    character(len=48), dimension(:), allocatable :: ubmtbi
    character(len=48), dimension(:), allocatable :: ucospi
    character(len=48), dimension(:), allocatable :: uspec
    character(len=48), dimension(:), allocatable :: uspeci
    character(len=48), dimension(:), allocatable :: uxmod
    character(len=48), dimension(:), allocatable :: uzveci
    character(len=48), dimension(:), allocatable :: uzvec1
    character(len=48), dimension(:,:), allocatable :: uobsw
    character(len=48), dimension(:,:), allocatable :: usbsw
    character(len=24), dimension(:), allocatable :: umemi
    character(len=24), dimension(:), allocatable :: usoli
    character(len=24), dimension(:), allocatable :: uphase
    character(len=24), dimension(:), allocatable :: uptype
    character(len=8), dimension(:), allocatable :: uelem
    character(len=8), dimension(:), allocatable :: uldel
    character(len=8), dimension(:), allocatable :: ulbeta

    character(len=32) :: ujflls(0:njf_par)
    character(len=32) :: uxtype(jso_par)
    character(len=24) :: ugexpi(net_par)
    character(len=24) :: ugexsi(iet_par,jet_par,net_par)
    character(len=8) :: ugexji(jet_par,net_par)

    character(len=56) :: uspn56
    character(len=48) :: ubacmx
    character(len=48) :: ubbig
    character(len=48) :: ubetmx
    character(len=48) :: ubgamx
    character(len=48) :: ubneg
    character(len=32) :: uactop
    character(len=24) :: uaqsln
    character(len=24) :: ublk24
    character(len=24) :: uebal
    character(len=24) :: uredox
    character(len=24) :: ux24
    character(len=16) :: ux16a
    character(len=16) :: ux16b
    character(len=11) :: utime0
    character(len=11) :: utime1
    character(len=9) :: udate0
    character(len=9) :: udate1
    character(len=8) :: udatfi
    character(len=8) :: udakey
    character(len=8) :: uinfor
    character(len=8) :: upkfor
    character(len=8) :: uplatc
    character(len=8) :: uplatm
    character(len=8) :: ustelg
    character(len=8) :: ustelu
    character(len=8) :: usteql
    character(len=8) :: usteq3
    character(len=8) :: uveelg
    character(len=8) :: uveelu
    character(len=8) :: uveeql
    character(len=8) :: uveeq3

    character(len=8) :: uv
    character(len=8) :: ux8
    character(len=8) :: ux8a
    character(len=8) :: ux8b
    character(len=80) :: ux80

    real(kind=8), dimension(:), allocatable :: covali

    real(kind=8), dimension(:), allocatable :: acflg
    real(kind=8), dimension(:), allocatable :: acflgo
    real(kind=8), dimension(:), allocatable :: act
    real(kind=8), dimension(:), allocatable :: actlg
    real(kind=8), dimension(:), allocatable :: affpd
    real(kind=8), dimension(:), allocatable :: affsd
    real(kind=8), dimension(:), allocatable :: ahrc
    real(kind=8), dimension(:), allocatable :: alpha
    real(kind=8), dimension(:), allocatable :: amtb
    real(kind=8), dimension(:), allocatable :: atwt
    real(kind=8), dimension(:), allocatable :: azero
    real(kind=8), dimension(:), allocatable :: a3bars
    real(kind=8), dimension(:), allocatable :: beta
    real(kind=8), dimension(:), allocatable :: betao
    real(kind=8), dimension(:), allocatable :: bfac
    real(kind=8), dimension(:), allocatable :: cdrs
    real(kind=8), dimension(:), allocatable :: cdrsd
    real(kind=8), dimension(:), allocatable :: cdrsx
    real(kind=8), dimension(:), allocatable :: cdrtw
    real(kind=8), dimension(:), allocatable :: cdrw
    real(kind=8), dimension(:), allocatable :: cess
    real(kind=8), dimension(:), allocatable :: cjbasp
    real(kind=8), dimension(:), allocatable :: cnufac
    real(kind=8), dimension(:), allocatable :: conc
    real(kind=8), dimension(:), allocatable :: conclg
    real(kind=8), dimension(:), allocatable :: coval
    real(kind=8), dimension(:), allocatable :: csts
    real(kind=8), dimension(:), allocatable :: ctb
    real(kind=8), dimension(:), allocatable :: cteaq
    real(kind=8), dimension(:), allocatable :: delvco
    real(kind=8), dimension(:), allocatable :: delvec
    real(kind=8), dimension(:), allocatable :: dlogxw
    real(kind=8), dimension(:), allocatable :: efac
    real(kind=8), dimension(:), allocatable :: ehrc
    real(kind=8), dimension(:), allocatable :: fo2lrc
    real(kind=8), dimension(:), allocatable :: fsort
    real(kind=8), dimension(:), allocatable :: fugac
    real(kind=8), dimension(:), allocatable :: fugalg
    real(kind=8), dimension(:), allocatable :: loph
    real(kind=8), dimension(:), allocatable :: losp
    real(kind=8), dimension(:), allocatable :: lsort
    real(kind=8), dimension(:), allocatable :: moph
    real(kind=8), dimension(:), allocatable :: mosp
    real(kind=8), dimension(:), allocatable :: mtb
    real(kind=8), dimension(:), allocatable :: mtbaq
    real(kind=8), dimension(:), allocatable :: mtbaqi
    real(kind=8), dimension(:), allocatable :: mtbi
    real(kind=8), dimension(:), allocatable :: mte
    real(kind=8), dimension(:), allocatable :: mteaq
    real(kind=8), dimension(:), allocatable :: mwtsp
    real(kind=8), dimension(:), allocatable :: perc
    real(kind=8), dimension(:), allocatable :: ppmwe
    real(kind=8), dimension(:), allocatable :: rhsvec
    real(kind=8), dimension(:), allocatable :: sidrph
    real(kind=8), dimension(:), allocatable :: sidrsp
    real(kind=8), dimension(:), allocatable :: tempcu
    real(kind=8), dimension(:), allocatable :: tfx
    real(kind=8), dimension(:), allocatable :: tf1
    real(kind=8), dimension(:), allocatable :: tf2
    real(kind=8), dimension(:), allocatable :: vosp0
    real(kind=8), dimension(:), allocatable :: weight
    real(kind=8), dimension(:), allocatable :: xbar
    real(kind=8), dimension(:), allocatable :: xbari
    real(kind=8), dimension(:), allocatable :: xbarlg
    real(kind=8), dimension(:), allocatable :: xlkmod
    real(kind=8), dimension(:), allocatable :: zchar
    real(kind=8), dimension(:), allocatable :: zchcu6
    real(kind=8), dimension(:), allocatable :: zchsq2
    real(kind=8), dimension(:), allocatable :: zvclgi
    real(kind=8), dimension(:), allocatable :: zvclg1
    real(kind=8), dimension(:), allocatable :: zvec1

    real(kind=8), dimension(:,:), allocatable :: amu
    real(kind=8), dimension(:,:,:), allocatable :: aslm

    real(kind=8), dimension(:), allocatable :: pmu
    real(kind=8), dimension(:), allocatable :: pslm
    real(kind=8), dimension(:,:), allocatable :: dpslm
    real(kind=8), dimension(:,:), allocatable :: gpit
    real(kind=8), dimension(:,:), allocatable :: palpha
    real(kind=8), dimension(:,:), allocatable :: pslamn
    real(kind=8), dimension(:,:,:), allocatable :: dgpit

    real(kind=8), dimension(:,:), allocatable :: elam
    real(kind=8), dimension(:,:), allocatable :: pelm
    real(kind=8), dimension(:,:,:), allocatable :: delam
    real(kind=8), dimension(:,:,:), allocatable :: dpelm

    real(kind=8), dimension(:), allocatable :: selm
    real(kind=8), dimension(:,:), allocatable :: dselm

    real(kind=8), dimension(:,:), allocatable :: aprehw
    real(kind=8), dimension(:,:), allocatable :: apresg

    real(kind=8), dimension(:), allocatable :: dadhh
    real(kind=8), dimension(:), allocatable :: dadhv
    real(kind=8), dimension(:), allocatable :: dbdhh
    real(kind=8), dimension(:), allocatable :: dbdhv
    real(kind=8), dimension(:), allocatable :: dbdth
    real(kind=8), dimension(:), allocatable :: dbdtv

    real(kind=8), dimension(:), allocatable :: xhfs
    real(kind=8), dimension(:), allocatable :: xhfsd
    real(kind=8), dimension(:), allocatable :: xlks
    real(kind=8), dimension(:), allocatable :: xlksd
    real(kind=8), dimension(:), allocatable :: xvfs
    real(kind=8), dimension(:), allocatable :: xvfsd
    real(kind=8), dimension(:), allocatable :: dhfe
    real(kind=8), dimension(:), allocatable :: dvfe
    real(kind=8), dimension(:,:), allocatable :: axhfe
    real(kind=8), dimension(:,:), allocatable :: axlke
    real(kind=8), dimension(:,:), allocatable :: axvfe
    real(kind=8), dimension(:,:,:), allocatable :: adhfe
    real(kind=8), dimension(:,:,:), allocatable :: advfe

    real(kind=8), dimension(:,:), allocatable :: dhfs
    real(kind=8), dimension(:,:), allocatable :: dhfsd
    real(kind=8), dimension(:,:), allocatable :: dvfs
    real(kind=8), dimension(:,:), allocatable :: dvfsd
    real(kind=8), dimension(:,:,:), allocatable :: axhfs
    real(kind=8), dimension(:,:,:), allocatable :: axhfsd
    real(kind=8), dimension(:,:,:), allocatable :: axhfsx
    real(kind=8), dimension(:,:,:), allocatable :: axlks
    real(kind=8), dimension(:,:,:), allocatable :: axlksd
    real(kind=8), dimension(:,:,:), allocatable :: axlksx
    real(kind=8), dimension(:,:,:), allocatable :: axvfs
    real(kind=8), dimension(:,:,:), allocatable :: axvfsd
    real(kind=8), dimension(:,:,:), allocatable :: axvfsx
    real(kind=8), dimension(:,:,:,:), allocatable :: adhfs
    real(kind=8), dimension(:,:,:,:), allocatable :: adhfsd
    real(kind=8), dimension(:,:,:,:), allocatable :: adhfsx
    real(kind=8), dimension(:,:,:,:), allocatable :: advfs
    real(kind=8), dimension(:,:,:,:), allocatable :: advfsd
    real(kind=8), dimension(:,:,:,:), allocatable :: advfsx

    real(kind=8), dimension(:,:), allocatable :: apx
    real(kind=8), dimension(:,:), allocatable :: bpx
    real(kind=8), dimension(:,:), allocatable :: wfac
    real(kind=8), dimension(:,:), allocatable :: aamatr
    real(kind=8), dimension(:,:), allocatable :: gmmatr

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

    real(kind=8) :: cco2(5)
    real(kind=8) :: cegexs(iet_par,jet_par,net_par)
    real(kind=8) :: cgexpi(net_par)
    real(kind=8) :: cgexp(net_par)
    real(kind=8) :: cpgexs(iet_par,jet_par,net_par)
    real(kind=8) :: egexjc(jet_par,net_par)
    real(kind=8) :: egexjf(jet_par,net_par)
    real(kind=8) :: egexpa(net_par)
    real(kind=8) :: egexpc(net_par)
    real(kind=8) :: egexs(iet_par,jet_par,net_par)
    real(kind=8) :: egexsi(iet_par,jet_par,net_par)
    real(kind=8) :: egexw(ket_par,net_par)
    real(kind=8) :: mrgexs(iet_par,jet_par,net_par)
    real(kind=8) :: xgexsi(iet_par,jet_par,net_par)
    real(kind=8) :: xgexw(ket_par,net_par)

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

    real(kind=8) :: prehw
    real(kind=8) :: presg
    real(kind=8) :: xhfe
    real(kind=8) :: xlke
    real(kind=8) :: xvfe

    real(kind=8) :: abar
    real(kind=8) :: actwlc
    real(kind=8) :: afcnst
    real(kind=8) :: alki
    real(kind=8) :: al10
    real(kind=8) :: avcnst
    real(kind=8) :: a3bar
    real(kind=8) :: azch
    real(kind=8) :: azchmx
    real(kind=8) :: bacfmx
    real(kind=8) :: bbig
    real(kind=8) :: betamx
    real(kind=8) :: bfje
    real(kind=8) :: bgamx
    real(kind=8) :: bneg
    real(kind=8) :: bsigmm
    real(kind=8) :: bfxi
    real(kind=8) :: delmax
    real(kind=8) :: eh
    real(kind=8) :: ehfac
    real(kind=8) :: ehi
    real(kind=8) :: electr
    real(kind=8) :: eps100
    real(kind=8) :: farad
    real(kind=8) :: fje
    real(kind=8) :: fjeo
    real(kind=8) :: fo2
    real(kind=8) :: fo2lg
    real(kind=8) :: fo2lgi
    real(kind=8) :: fxi
    real(kind=8) :: fxio
    real(kind=8) :: omega
    real(kind=8) :: omeglg
    real(kind=8) :: pe
    real(kind=8) :: pei
    real(kind=8) :: presmx
    real(kind=8) :: press
    real(kind=8) :: pressi
    real(kind=8) :: rconst
    real(kind=8) :: rcnstv
    real(kind=8) :: rho
    real(kind=8) :: rtcnst
    real(kind=8) :: scamas
    real(kind=8) :: screwd
    real(kind=8) :: screwn
    real(kind=8) :: sigmam
    real(kind=8) :: sigmmo
    real(kind=8) :: sigzi
    real(kind=8) :: smp100
    real(kind=8) :: tcpu
    real(kind=8) :: tdamax
    real(kind=8) :: tdamin
    real(kind=8) :: tdspkg
    real(kind=8) :: tdspl
    real(kind=8) :: tempc
    real(kind=8) :: tempci
    real(kind=8) :: tempk
    real(kind=8) :: texec0
    real(kind=8) :: tolbt
    real(kind=8) :: toldl
    real(kind=8) :: tolspf
    real(kind=8) :: trun
    real(kind=8) :: tuser
    real(kind=8) :: vosol
    real(kind=8) :: wfh2o
    real(kind=8) :: wftds
    real(kind=8) :: wkgwi
    real(kind=8) :: woh2o
    real(kind=8) :: wosol
    real(kind=8) :: wotds
    real(kind=8) :: xbarw
    real(kind=8) :: xbarwc
    real(kind=8) :: xbrwlc
    real(kind=8) :: xbrwlg
    real(kind=8) :: x10
    real(kind=8) :: zx

    ! XX   real(8) vpgstp
    real(kind=8) :: mlmrra
    real(kind=8) :: mrmlra
    real(kind=8) :: rhoc
    real(kind=8) :: rhowc
    real(kind=8) :: tdsgks
    real(kind=8) :: tdsglw
    real(kind=8) :: tdspkc
    real(kind=8) :: tdsplc

    real(kind=8) :: axx
    real(kind=8) :: av
    real(kind=8) :: dp
    real(kind=8) :: pxl
    real(kind=8) :: pxu
    real(kind=8) :: xx

    real(kind=8) :: texp
    real(kind=8) :: tlg

    ! Variable declarations: Data file original contents. These
    ! variables and arrays remain unchanged after being filled by
    ! reading the data file. These data comprise the 'a' set.
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

    ! Variable declarations: Data needed to write a complete EQ6 input
    ! file as the EQ3NR pickup file.
    integer :: igerti(jet_par,nert_par)
    integer :: jgerti(nert_par)

    integer, dimension(:), allocatable :: ibsrti
    integer, dimension(:), allocatable :: iesrti
    integer, dimension(:), allocatable :: ixrti
    integer, dimension(:), allocatable :: jcode
    integer, dimension(:), allocatable :: jreac
    integer, dimension(:), allocatable :: nsk

    integer, dimension(:,:), allocatable :: imech
    integer, dimension(:,:), allocatable :: nrk

    integer, dimension(:,:,:), allocatable :: iact
    integer, dimension(:,:,:), allocatable :: ndact

    integer :: jpress
    integer :: jtemp
    integer :: ksplmx
    integer :: ksppmx
    integer :: kstpmx
    integer :: nert
    integer :: nffg
    integer :: nprpti
    integer :: nprsti
    integer :: nrct
    integer :: nsrt
    integer :: ntitl1
    integer :: ntrymx
    integer :: nxopex
    integer :: nxopt
    integer :: nxrt

    logical :: qgexsh

    character(len=24) :: ugersi(iet_par,jet_par,nert_par)
    character(len=24) :: ugermo(nert_par)
    character(len=8) :: ugerji(jet_par,nert_par)

    character(len=80), dimension(:), allocatable :: utitl1
    character(len=48), dimension(:), allocatable :: uprspi
    character(len=24), dimension(:), allocatable :: uffg
    character(len=24), dimension(:), allocatable :: uprphi
    character(len=24), dimension(:), allocatable :: ureac
    character(len=24), dimension(:), allocatable :: uxcat
    character(len=24), dimension(:), allocatable :: uxopex
    character(len=8), dimension(:), allocatable :: uxopt

    character(len=24), dimension(:,:), allocatable :: ubsri
    character(len=24), dimension(:,:), allocatable :: ucxri
    character(len=8), dimension(:,:), allocatable :: uesri

    character(len=24), dimension(:,:,:,:), allocatable :: udac

    real(kind=8), dimension(:), allocatable :: fkrc
    real(kind=8), dimension(:), allocatable :: modr
    real(kind=8), dimension(:), allocatable :: moffg
    real(kind=8), dimension(:), allocatable :: morr
    real(kind=8), dimension(:), allocatable :: mprphi
    real(kind=8), dimension(:), allocatable :: mprspi
    real(kind=8), dimension(:), allocatable :: ptk
    real(kind=8), dimension(:), allocatable :: sfcar
    real(kind=8), dimension(:), allocatable :: ssfcar
    real(kind=8), dimension(:), allocatable :: ttk
    real(kind=8), dimension(:), allocatable :: vreac
    real(kind=8), dimension(:), allocatable :: xlkffg
    real(kind=8), dimension(:,:), allocatable :: cbsri
    real(kind=8), dimension(:,:), allocatable :: cesri
    real(kind=8), dimension(:,:), allocatable :: rxbari
    real(kind=8), dimension(:,:,:), allocatable :: csigma
    real(kind=8), dimension(:,:,:), allocatable :: eact
    real(kind=8), dimension(:,:,:), allocatable :: hact
    real(kind=8), dimension(:,:,:), allocatable :: rkb
    real(kind=8), dimension(:,:,:), allocatable :: trkb
    real(kind=8), dimension(:,:,:,:), allocatable :: cdac

    real(kind=8) :: egersi(iet_par,jet_par,nert_par)
    real(kind=8) :: xgersi(iet_par,jet_par,nert_par)

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
    real(kind=8) :: pressb
    real(kind=8) :: tempcb
    real(kind=8) :: tempc1
    real(kind=8) :: timmxi
    real(kind=8) :: tistti
    real(kind=8) :: tolsat
    real(kind=8) :: tolxsf
    real(kind=8) :: ximaxi
    real(kind=8) :: xistti

    ! Variable declarations: Local data needed to assist in writing a
    ! complete EQ6 input file as the EQ3NR pickup file. These variables
    ! do not appear on that file itself.
    integer :: ibsrt1
    integer :: iesrt1

    character(len=24), dimension(:), allocatable :: ubsr1
    character(len=8), dimension(:), allocatable :: uesr1

    character(len=24) :: ureac1

    real(kind=8), dimension(:), allocatable :: cbsr1
    real(kind=8), dimension(:), allocatable :: cesr1

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
    external bkdeq3

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

    data ujflls(0)  /'Total molality                  '/
    data ujflls(1)  /'Total molarity                  '/
    data ujflls(2)  /'Total mg/L                      '/
    data ujflls(3)  /'Total mg/kg.sol                 '/
    data ujflls(4)  /'ERROR                           '/
    data ujflls(5)  /'ERROR                           '/
    data ujflls(6)  /'ERROR                           '/
    data ujflls(7)  /'Total alkalinity, eq/kg.H2O     '/
    data ujflls(8)  /'Total alkalinity, eq/L          '/
    data ujflls(9)  /'Total alkalinity, eq/kg.sol     '/
    data ujflls(10) /'Total alkalinity, mg/L CaCO3    '/
    data ujflls(11) /'Total alkalinity, mg/L HCO3-    '/
    data ujflls(12) /'ERROR                           '/
    data ujflls(13) /'ERROR                           '/
    data ujflls(14) /'ERROR                           '/
    data ujflls(15) /'ERROR                           '/
    data ujflls(16) /'Log activity                    '/
    data ujflls(17) /'|zj| log ai +/- |zi| log aj     '/
    data ujflls(18) /'Log a(+/-,ij)                   '/
    data ujflls(19) /'pX                              '/
    data ujflls(20) /'pH                              '/
    data ujflls(21) /'pHCl                            '/
    data ujflls(22) /'pmH                             '/
    data ujflls(23) /'pmX                             '/
    data ujflls(24) /'ERROR                           '/
    data ujflls(25) /'Heterogenous equilibrium        '/
    data ujflls(26) /'ERROR                           '/
    data ujflls(27) /'Homogenous equilibrium          '/
    data ujflls(28) /'ERROR                           '/
    data ujflls(29) /'ERROR                           '/
    data ujflls(30) /'Make non-basis                  '/

    data uxtype(1)  /'Ideal solution                  '/
    data uxtype(2)  /'Binary, third-order Maclaurin   '/
    data uxtype(3)  /'Binary, parabolic Maclaurin     '/
    data uxtype(4)  /'Binary, cubic Maclaurin (P,T)   '/
    data uxtype(5)  /'Binary, Guggenheim  (T)         '/
    data uxtype(6)  /'Ternary, regular                '/

    data nrecl /0/

    ! Under-relaxation control parameters:
    !   screwd = bound on max norm of applied part of the delvec vector
    !   screwn = factor bounding the increase in the max norm of the
    !            beta vector
    data screwd/2.0/,screwn/0.50/

    ! Fixed fugacity phase range delimiters (not currently used
    ! by EQ3NR):
    data ifrn1,ifrn2 /0,0/

    ! Maximum allowed pressure (bars):
    data presmx /10000./

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
    ninpt = 0
    ninpts = 0
    newin = 0

    ! Open all files execpt pickup.
    numargs = COMMAND_ARGUMENT_COUNT()

    if (numargs.eq.2) then
        call GET_COMMAND_ARGUMENT(1,temppath)
        data1path = TRIM(temppath)
        temppath(:)='\0'
        call GET_COMMAND_ARGUMENT(2,temppath)
        threeipath = TRIM(temppath)
    else
        write (0, *) 'usage: eq3nr <data1> <3i>'
        stop
    end if

    call getbasename(threeipath, pathindices)
    basename = threeipath(pathindices(1):pathindices(2))

    ofile = basename // '.3o'
    ifile = basename // '.3ib'
    pfile = basename // '.3p'

    call openin(noutpt,nttyo,data1path,'unformatted',nad1)
    call openin(noutpt,nttyo,threeipath,'formatted',ninpt)

    call openou(noutpt,nttyo,ofile,'formatted',nrecl,noutpt)
    call openou(noutpt,nttyo,ifile,'formatted',nrecl,ninpts)

    ! Make a copy of the input file, stripped of comments.
    call stripl(ninpt,ninpts)
    close (ninpt)
    ninpt = 0
    rewind ninpts

    ! Get configuration identification data.
    call aaaeq3(usteq3,uveeq3)
    call aaaeql(usteql,uveeql)
    call aaaelg(ustelg,uveelg)
    call aaaelu(ustelu,uveelu)
    call platfd(uplatc,uplatm)

    ! Write configuration identification data, the copyright statement,
    ! and any remaining statements or disclaimers.
    i = index(uveeq3,' ') - 1
    j = index(uveeq3,'.') - 1
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

    write (nttyo,1000) uveeq3(1:i),uveeq3(1:j),uveeq3(1:i),uplatc(1:k)
    write (noutpt,1000) uveeq3(1:i),uveeq3(1:j),uveeq3(1:i),uplatc(1:k)
1000 format(/' EQ3/6, Version ',a,' (EQ3/6-V',a,'-REL-V',a,'-',a,')')

    i = j
    j = index(usteq3,' ') - 1
    k = index(uplatm,' ') - 1

    if (j .le. 0) then
        j = 8
    end if

    if (k .le. 0) then
        k = 8
    end if

    write (nttyo,1001) uveeq3(1:i),usteq3(1:j),uplatm(1:k)
    write (noutpt,1001) uveeq3(1:i),usteq3(1:j),uplatm(1:k)
1001 format(' EQ3NR Speciation-Solubility Code (EQ3/6-V',a,'-EQ3NR-EXE-',a,'-',a,')')

    write (noutpt,1002)
    write (nttyo,1002)
1002 format(' Supported by the following EQ3/6 libraries:')

    i = index(uveeql,'.') - 1
    j = index(usteql,' ') - 1

    if (i .le. 0) then
        i = 8
    end if

    if (j .le. 0) then
        j = 8
    end if

    write (nttyo,1003) uveeql(1:i),usteql(1:j),uplatm(1:k)
    write (noutpt,1003) uveeql(1:i),usteql(1:j),uplatm(1:k)
1003 format('   EQLIB (EQ3/6-V',a,'-EQLIB-LIB-',a,'-',a,')')

    i = index(uveelg,'.') - 1
    j = index(ustelg,' ') - 1

    if (i .le. 0) then
        i = 8
    end if

    if (j .le. 0) then
        j = 8
    end if

    write (nttyo,1004) uveelg(1:i),ustelg(1:j),uplatm(1:k)
    write (noutpt,1004) uveelg(1:i),ustelg(1:j),uplatm(1:k)
1004 format('   EQLIBG (EQ3/6-V',a,'-EQLIBG-LIB-',a,'-',a,')')

    i = index(uveelu,'.') - 1
    j = index(ustelu,' ') - 1

    if (i .le. 0) then
        i = 8
    end if

    if (j .le. 0) then
        j = 8
    end if

    write (nttyo,1005) uveelu(1:i),ustelu(1:j),uplatm(1:k)
    write (noutpt,1005) uveelu(1:i),ustelu(1:j),uplatm(1:k)
1005 format('   EQLIBU (EQ3/6-V',a,'-EQLIBU-LIB-',a,'-',a,')',/)

    write (nttyo,1010)
    write (noutpt,1010)
1010 format(' Copyright (c) 1987, 1990-1993, 1995, 1997, 2002 The',' Regents of the',/' University of California, Lawrence',' Livermore National Laboratory.',/' All rights reserved.',/)

    ! Write additional statements and disclaimers.
    call prcndi(noutpt,nttyo)

    ! Write the time and date on the output.
    j2 = ilnobl(udate0)
    write (noutpt,1080) utime0,udate0(1:j2)
    write (nttyo,1080) utime0,udate0(1:j2)
1080 format(' Run',2x,a8,2x,a,/)

    ! Get the platform's real(8) floating-point parameters.
    call flpars(eps100,irang,noutpt,nttyo,smp100)

    ! Initialize array dimension variables corresponding to fixed
    ! parameters.
    ! The following are common with EQ6.
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

    ! The following are common with EQ6, but are only needed to write
    ! a full EQ6 input file as the EQ3NR pickup file.
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

    ! Special over-rides for EQ3NR.
    nrct_asv = 2
    nsrt_asv = 1
    nxrt_asv = 1
    nprp_asv = 1
    nprs_asv = 1
    nxop_asv = 1
    nxpe_asv = 1

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

    ! The following are unique to EQ3NR.
    nxti_asv = nxti_par

    njfmax = njf_par

    nxtimx = nxti_asv

    ! Initialize the following floating-point constants:
    !   rconst = the gas constant: 1.98726 cal/mol-K
    !   rcnstv = the gas constant: 83.14510 bar-cm3/mol-K
    !   vpgstp = the molar volume of a perfect gas at STP:
    !              22413.6 cm3/mol
    !   farad  = the Faraday constant: 23062.3 cal/equiv-volt
    !   al10   = ln 10; =~ 2.3026
    rconst = 1.98726
    rcnstv = 83.14510

    ! XX   vpgstp = 22413.6
    farad = 23062.3
    x10 = 10.
    al10 = log(x10)

    ! Set the EQ6 calculational mode flag (.false. in EQ3NR).
    q6mode = .false.

    ! Read the header section of the data1 file. This section consists
    ! of a record containing the string 'data1' (to ensure that the file
    ! is indeed a data1 file), a record containing the keystring for the
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

    ! Bump up the values of some allocation size variables to allow
    ! for some species, phases, etc., to created by this software.
    ! Examples: generic ion exchanger phases and species.
    nbta_asv = nbta_asv + 10
    npta_asv = npta_asv + net_par + 10
    nsta_asv = nsta_asv + iet_par*jet_par*net_par + 10

    ! Set the values of some secondary allocation size variables.
    nbta1_asv = nbta_asv + 1
    ndrsa_asv = 7*nsta_asv
    nessa_asv = 5*nsta_asv

    ! Allocate arrays to store the data read from the rest of the data1
    ! file.
    ALLOCATE(iapxta(nxta_asv))
    ALLOCATE(ibpxta(nxta_asv))
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

    ALLOCATE(axhfe(narx_asv,ntpr_asv))
    ALLOCATE(axlke(narx_asv,ntpr_asv))
    ALLOCATE(axvfe(narx_asv,ntpr_asv))

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
    ALLOCATE(insgfa(nata_asv))

    ! Pitzer data.
    ALLOCATE(nalpaa(nslta_asv))
    ALLOCATE(nmuxa(3,nmuta_asv))
    ALLOCATE(nslxa(2,nslta_asv))

    ALLOCATE(amua(jpfc_asv,nmuta_asv))

    ALLOCATE(aslma(jpfc_asv,0:ipbt_asv,nslta_asv))
    ALLOCATE(palpaa(ipbt_asv,napa_asv))

    ! Read the remainder of the data1 file. The image of the data file
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

    ! Get the indices of H+, OH-, Cl-, aqueous O2(g), and aqueous e-.
    ! This will be repeated after data compression.
    ! Calling sequence substitutions:
    !   narn1a for narn1
    !   narn2a for narn2
    !   nchloa for nchlor
    !   nhydra for nhydr
    !   neleca for nelect
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

    ! Get some constants that depend on the molecular weight of
    ! water.
    omega = 1000./mwtspa(narn1a)
    omeglg = log10(omega)

    ! Setup up the remaining array allocation variables required to
    ! allocate variables needed to store the dta to be read from the
    ! input file.
    ! The following value of k_asv is special to EQ3NR.
    ! A larger value is necessary for EQ6.
    k_asv = nbta_asv + 5

    ! The following is unique to EQ3NR.
    nxic_asv = ikta_asv*nxti_asv
    nxicmx = nxic_asv

    ! Allocate arrays to store the data to be read from the input
    ! file.
    ALLOCATE(jflgi(nbta_asv))
    ALLOCATE(nbaspi(nbta_asv))

    ALLOCATE(npnxp(nxti_asv))

    ALLOCATE(ncmpri(2,nxti_asv))

    ALLOCATE(uspeci(nbta_asv))
    ALLOCATE(umemi(nxic_asv))
    ALLOCATE(usoli(nxti_asv))

    ALLOCATE(covali(nbta_asv))

    ALLOCATE(mtbaqi(nbta_asv))
    ALLOCATE(mtbi(nbta_asv))
    ALLOCATE(xbari(nxic_asv))
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
    ntfx_asv = max(ntf1_asv,ntf2_asv)
    nmx_asv = 3*nmut_asv
    nsx_asv = 2*nslt_asv

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
    nptmax = npt_asv
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
    ALLOCATE(ipivot(k_asv))

    ALLOCATE(insgf(nat_asv))

    ALLOCATE(ibswx(nbt_asv))
    ALLOCATE(ixbasp(nbt_asv))
    ALLOCATE(nbasp(nbt_asv))
    ALLOCATE(nbaspd(nbt_asv))
    ALLOCATE(nbaspx(nbt_asv))
    ALLOCATE(nbmap(nbt_asv))
    ALLOCATE(ncosp(nbt_asv))
    ALLOCATE(ndecsp(nbt_asv))
    ALLOCATE(nfac(nbt_asv))

    ALLOCATE(ncmap(nct_asv))

    ALLOCATE(iction(nbt_asv))
    ALLOCATE(kction(nbt_asv))
    ALLOCATE(jjndex(nbt_asv))
    ALLOCATE(kkndex(nbt_asv))

    ALLOCATE(ndrs(ndrs_asv))
    ALLOCATE(ndrsd(ndrs_asv))
    ALLOCATE(ndrsx(ndrs_asv))

    ALLOCATE(ness(ness_asv))

    ALLOCATE(igstak(ngt_asv))
    ALLOCATE(jgstak(ngt_asv))
    ALLOCATE(jgsort(ngt_asv))

    ALLOCATE(jpflag(npt_asv))
    ALLOCATE(ncmpr(2,npt_asv))

    ALLOCATE(istack(nst_asv))
    ALLOCATE(jstack(nst_asv))
    ALLOCATE(jcsort(nst_asv))
    ALLOCATE(jjsort(nst_asv))
    ALLOCATE(jssort(nst_asv))
    ALLOCATE(jflag(nst_asv))
    ALLOCATE(jflagd(nst_asv))
    ALLOCATE(jsflag(nst_asv))
    ALLOCATE(jsitex(nst_asv))
    ALLOCATE(nphasx(nst_asv))
    ALLOCATE(nsmap(nst_asv))

    ALLOCATE(ntfx(ntfx_asv))
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

    ALLOCATE(qxknph(npt_asv))

    ALLOCATE(utitl(ntit_asv))
    ALLOCATE(utitl2(ntit_asv))

    ALLOCATE(ubmtbi(nbt_asv))
    ALLOCATE(ucospi(nbt_asv))

    ALLOCATE(uobsw(2,nbt_asv))
    ALLOCATE(usbsw(2,nbt_asv))

    ALLOCATE(uspec(nst_asv))
    ALLOCATE(uxmod(nxmd_asv))
    ALLOCATE(uzveci(k_asv))
    ALLOCATE(uzvec1(k_asv))

    ALLOCATE(uphase(npt_asv))
    ALLOCATE(uptype(npt_asv))

    ALLOCATE(uelem(nct_asv))

    ALLOCATE(uldel(k_asv))
    ALLOCATE(ulbeta(k_asv))

    ALLOCATE(apx(iapx_asv,nxt_asv))
    ALLOCATE(bpx(ibpx_asv,nxt_asv))
    ALLOCATE(wfac(ikt_asv,nxt_asv))

    ALLOCATE(aamatr(k_asv,k_asv))
    ALLOCATE(gmmatr(k_asv,k_asv))

    ALLOCATE(acflg(nst_asv))
    ALLOCATE(acflgo(nst_asv))
    ALLOCATE(act(nst_asv))
    ALLOCATE(actlg(nst_asv))
    ALLOCATE(affsd(nst_asv))

    ALLOCATE(affpd(npt_asv))
    ALLOCATE(loph(npt_asv))
    ALLOCATE(moph(npt_asv))
    ALLOCATE(sidrph(npt_asv))

    ALLOCATE(azero(nat_asv))
    ALLOCATE(a3bars(nat_asv))

    ALLOCATE(ahrc(nbt_asv))
    ALLOCATE(amtb(nbt_asv))
    ALLOCATE(cjbasp(nbt_asv))

    ALLOCATE(atwt(nct_asv))
    ALLOCATE(cteaq(nct_asv))
    ALLOCATE(mte(nct_asv))
    ALLOCATE(mteaq(nct_asv))
    ALLOCATE(ppmwe(nct_asv))

    ALLOCATE(alpha(k_asv))
    ALLOCATE(beta(k_asv))
    ALLOCATE(betao(k_asv))
    ALLOCATE(delvco(k_asv))
    ALLOCATE(delvec(k_asv))
    ALLOCATE(rhsvec(k_asv))
    ALLOCATE(zvclg1(k_asv))
    ALLOCATE(zvec1(k_asv))

    ALLOCATE(bfac(nbt_asv))
    ALLOCATE(coval(nbt_asv))
    ALLOCATE(ctb(nbt_asv))
    ALLOCATE(dlogxw(nbt_asv))
    ALLOCATE(efac(nbt_asv))
    ALLOCATE(ehrc(nbt_asv))
    ALLOCATE(fo2lrc(nbt_asv))
    ALLOCATE(mtb(nbt_asv))
    ALLOCATE(mtbaq(nbt_asv))
    ALLOCATE(perc(nbt_asv))

    ALLOCATE(cdrs(ndrs_asv))
    ALLOCATE(cdrsd(ndrs_asv))
    ALLOCATE(cdrsx(ndrs_asv))
    ALLOCATE(cdrtw(nst_asv))
    ALLOCATE(cdrw(nst_asv))
    ALLOCATE(cnufac(nst_asv))
    ALLOCATE(conc(nst_asv))
    ALLOCATE(conclg(nst_asv))
    ALLOCATE(losp(nst_asv))
    ALLOCATE(lsort(nst_asv))
    ALLOCATE(mosp(nst_asv))
    ALLOCATE(mwtsp(nst_asv))
    ALLOCATE(sidrsp(nst_asv))
    ALLOCATE(vosp0(nst_asv))
    ALLOCATE(weight(nst_asv))
    ALLOCATE(xbar(nst_asv))
    ALLOCATE(xbarlg(nst_asv))
    ALLOCATE(zchar(nst_asv))
    ALLOCATE(zchcu6(nst_asv))
    ALLOCATE(zchsq2(nst_asv))

    ALLOCATE(cess(ness_asv))
    ALLOCATE(tfx(ntfx_asv))
    ALLOCATE(tf1(ntf1_asv))
    ALLOCATE(tf2(ntf2_asv))
    ALLOCATE(xlkmod(nxmd_asv))

    ALLOCATE(fsort(ngt_asv))
    ALLOCATE(fugac(ngt_asv))
    ALLOCATE(fugalg(ngt_asv))

    ALLOCATE(xhfs(nst_asv))
    ALLOCATE(xhfsd(nst_asv))
    ALLOCATE(xlks(nst_asv))
    ALLOCATE(xlksd(nst_asv))
    ALLOCATE(xvfs(nst_asv))
    ALLOCATE(xvfsd(nst_asv))

    ALLOCATE(dhfe(ipch_asv))
    ALLOCATE(dvfe(ipcv_asv))

    ALLOCATE(dadhh(ipch_asv))
    ALLOCATE(dadhv(ipcv_asv))
    ALLOCATE(dbdhh(ipch_asv))
    ALLOCATE(dbdhv(ipcv_asv))
    ALLOCATE(dbdth(ipch_asv))
    ALLOCATE(dbdtv(ipcv_asv))

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

    ! Pitzer data.
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
    ALLOCATE(pslamn(0:ipbtmx,nslt_asv))

    ALLOCATE(gpit(ipbtmx,nap_asv))
    ALLOCATE(dgpit(2,ipbtmx,nap_asv))
    ALLOCATE(palpha(ipbtmx,nap_asv))

    ALLOCATE(elam(nazp_asv,nazp_asv))
    ALLOCATE(delam(2,nazp_asv,nazp_asv))
    ALLOCATE(pelm(nazp_asv,nazp_asv))
    ALLOCATE(dpelm(2,nazp_asv,nazp_asv))
    ALLOCATE(selm(nazm_asv:nazp_asv))
    ALLOCATE(dselm(2,nazm_asv:nazp_asv))

    ! Allocate arrays needed to write a complete EQ6 input file as
    ! the EQ3NR pickup file.
    ALLOCATE(jcode(nrct_asv))
    ALLOCATE(jreac(nrct_asv))
    ALLOCATE(nsk(nrct_asv))

    ALLOCATE(imech(2,nrct_asv))
    ALLOCATE(nrk(2,nrct_asv))

    ALLOCATE(iact(imch_asv,2,nrct_asv))
    ALLOCATE(ndact(imch_asv,2,nrct_asv))

    ALLOCATE(ibsrti(nsrt_asv))
    ALLOCATE(iesrti(nsrt_asv))

    ALLOCATE(ixrti(nxrt_asv))

    ALLOCATE(utitl1(ntit_asv))

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

    ALLOCATE(cbsri(nbt1_asv,nsrt_asv))
    ALLOCATE(cdac(ndct_asv,imch_asv,2,nrct_asv))
    ALLOCATE(cesri(nct_asv,nsrt_asv))
    ALLOCATE(csigma(imch_asv,2,nrct_asv))
    ALLOCATE(eact(imch_asv,2,nrct_asv))
    ALLOCATE(hact(imch_asv,2,nrct_asv))
    ALLOCATE(rkb(imch_asv,2,nrct_asv))
    ALLOCATE(trkb(imch_asv,2,nrct_asv))
    ALLOCATE(fkrc(nrct_asv))
    ALLOCATE(modr(nrct_asv))
    ALLOCATE(morr(nrct_asv))
    ALLOCATE(sfcar(nrct_asv))
    ALLOCATE(ssfcar(nrct_asv))
    ALLOCATE(vreac(nrct_asv))
    ALLOCATE(moffg(nffg_asv))
    ALLOCATE(xlkffg(nffg_asv))
    ALLOCATE(mprphi(nprp_asv))
    ALLOCATE(mprspi(nprs_asv))
    ALLOCATE(ptk(nptk_asv))
    ALLOCATE(ttk(nttk_asv))
    ALLOCATE(rxbari(ikt_asv,nxrt_asv))

    ! Allocate some arrays needed to assist in writing a complete
    ! EQ6 input file as the EQ3NR pickup file. These arrays do not
    ! appear on that file itself.
    ALLOCATE(uesr1(nct_asv))
    ALLOCATE(cesr1(nct_asv))
    ALLOCATE(cbsr1(nbt1_asv))
    ALLOCATE(ubsr1(nbt1_asv))

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
    qbassw = .false.

    call initiz(iopt,noptmx)
    call initiz(iopg,nopgmx)
    call initiz(iopr,noprmx)
    call initiz(iodb,nodbmx)

    call initiz(jflgi,nbtmax)
    call initaz(covali,nbtmax)

    call initcb(uspeci,nbtmax)
    call initcb(ucospi,nbtmax)

    nmax = 2*nbtmax
    call initcb(uobsw,nmax)
    call initcb(usbsw,nmax)

    call initiz(ncosp,nbtmax)
    call initiz(ndecsp,nbtmax)
    call initaz(coval,nbtmax)

    call initiz(nbasp,nbtmax)

    call initiz(jflag,nstmax)
    call initiz(jflagd,nstmax)
    call initcb(uspec,nstmax)
    call initaz(xlks,nstmax)
    call initaz(zchar,nstmax)
    call initaz(zchsq2,nstmax)
    call initaz(zchcu6,nstmax)

    ! xxxxxxxxxxx
    call initaz(cdrs,ndrsmx)
    call initaz(cdrsx,ndrsmx)
    call initiz(ndrs,ndrsmx)
    call initiz(ndrsx,ndrsmx)

    ! xxxxxxxxxxx
    nmax = 2*nstmax
    call initiz(ndrsr,nmax)
    call initiz(ndrsrx,nmax)

    nmax = narxmx*ntprmx*nstmax
    call initaz(axlks,nmax)
    call initaz(axhfs,nmax)
    call initaz(axvfs,nmax)

    ! xxxxxxxxxxx
    call initaz(axlksx,nmax)
    call initaz(axhfsx,nmax)
    call initaz(axvfsx,nmax)

    nmax = narxmx*ntprmx*ipchmx*nstmax
    call initaz(adhfs,nmax)
    call initaz(adhfsx,nmax)

    nmax = narxmx*ntprmx*ipcvmx*nstmax
    call initaz(advfs,nmax)
    call initaz(advfsx,nmax)

    ! xxxxxxxxxxx
    call initcb(ugexp,netmax)
    call initaz(cgexpi,netmax)
    call initiz(jgext,netmax)

    nmax = jetmax*netmax
    call initiz(ngexrt,nmax)

    nmax = ietmax*jetmax*netmax
    call initaz(xgexsi,nmax)

    nmax = iapxmx*nxtmax
    call initaz(apx,nmax)

    nmax = ibpxmx*nxtmax
    call initaz(bpx,nmax)

    ! Determine the problem input format.
    read (ninpts,1030,end=105,err=107) ux8
1030 format(a8)

    backspace ninpts
    uinfor = 'W'

    if (ux8(1:8) .eq. '|-------') then
        uinfor = 'D'
    end if

    ! Read the problem input.
    if (uinfor(1:1) .eq. 'W') then
        ! Compact (W) format.
        call rd3inw(cgexj,cgexpi,covali,ehi,egexsi,fo2lgi,iebal3,ietmax,iodb,iopg,iopr,iopt,irdxc3,itdsf3,itermx,jgext,jetmax,jflgi,jgexti,jpres3,kxmod,mwtges,nbti,nbtmax,ncmpri,net,neti,netmax,ngexti,ninpts,ngexrt,nobswt,nodbmx,nopgmx,noprmx,noptmx,noutpt,nprob,nsbswt,ntitl,ntitmx,nttyo,nxicmx,nxmdmx,nxmod,nxti,nxtimx,pei,press,qend,qgexsh,qrderr,rho,scamas,tdspkg,tdspl,tempc,tgexp,tolbt,toldl,tolspf,ucospi,uebal,ugexj,ugexji,ugexmo,ugexp,ugexpi,ugexr,ugexsi,umemi,uobsw,uredox,usbsw,usoli,uspeci,utitl,uhfgex,uvfgex,uxkgex,uxmod,xbari,xgexsi,xhfgex,xlkgex,xvfgex,xlkmod,zgexj)
    else
        ! Menu-style (D) format.
        call rd3ind(cgexj,cgexpi,covali,ehi,egexsi,fo2lgi,iebal3,ietmax,iodb,iopg,iopr,iopt,irdxc3,itdsf3,itermx,jgext,jetmax,jflgi,jgexti,jpres3,kxmod,mwtges,nbti,nbtmax,ncmpri,net,neti,netmax,ngexti,ninpts,ngexrt,nobswt,nodbmx,nopgmx,noprmx,noptmx,noutpt,nprob,nsbswt,ntitl,ntitmx,nttyo,nxicmx,nxmdmx,nxmod,nxti,nxtimx,pei,press,qend,qgexsh,qrderr,rho,scamas,tdspkg,tdspl,tempc,tgexp,tolbt,toldl,tolspf,ucospi,uebal,ugexj,ugexji,ugexmo,ugexp,ugexpi,ugexr,ugexsi,umemi,uobsw,uredox,usbsw,usoli,uspeci,utitl,uhfgex,uvfgex,uxkgex,uxmod,xbari,xgexsi,xhfgex,xlkgex,xvfgex,xlkmod,zgexj)
    end if

    go to 109
105 continue
    qend = .true.
    go to 109
107 continue
    qrderr = .true.
109 continue

    if (qrderr) then
        write (noutpt,1203)
        write (nttyo,1203)
1203 format(/' * Error - (EQ3NR/eq3nr) Encountered an error while',' reading the input file.',/7x,'The line associated with',' the problem should be the last one',/7x,'echoed on the',' output file or the first line immediately thereafter.')

        stop
    end if

    if (qend) then
        write (noutpt,1205)
        write (nttyo,1205)
1205 format(/' No further input found.',/)

        ! Get end time and date. Also get the run, user, and cpu times.
        call runtim(iexec0,jexec0,texec0,noutpt,nttyo,trun,tuser,tcpu,udate1,utime1)

        j2 = ilnobl(udate0)
        j3 = ilnobl(udate1)
        write (noutpt,1208) utime0,udate0(1:j2),utime1,udate1(1:j3)
        write (nttyo,1208) utime0,udate0(1:j2),utime1,udate1(1:j3)
1208 format(/10x,'Start time = ',a8,2x,a,/12x,'End time = ',a8,2x,a)

        ! Print the run, user, and cpu times.
        write (noutpt,1020) trun
        write (nttyo,1020)  trun
1020 format(/10x,' Run time = ',g10.3,' seconds')

        if (tuser .gt. 0.) then
            write (noutpt,1022) tuser
            write (nttyo,1022)  tuser
1022 format(10x,'User time = ',g10.3,' seconds')
        end if

        if (tcpu .gt. 0.) then
            write (noutpt,1024) tcpu
            write (nttyo,1024)  tcpu
1024 format(10x,' Cpu time = ',g10.3,' seconds')
        end if

        write (noutpt,1090)
        write (nttyo,1090)
1090 format(/' Normal exit')

        ! Clear the IEEE flag for floating-point underflow, if such a
        ! flag is present, to avoid getting an unnecessary system
        ! warning message. Underflow is a normal condition in EQ3/6.
        ! Make porting changes in the EQLIBU subroutine that is called
        ! in this section. Do not make the porting changes here.
        call cliefu()

        ! Close and delete the stripped input file.
        close (ninpts,status='delete')

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
        do n = 1,ntitl
            ux80 = utitl(n)
            call locase(ux80)
            i = index(ux80,'useoldpitzermu')

            if (i .gt. 0) then
                qhawep = .false.
                go to 106
            end if
        end do

106 continue
    end if

    if (.not.qhawep) then
        write (noutpt,1190)
        write (nttyo,1190)
1190 format(/' * Note - (EQ3NR/eq3nr) Found the USEOLDPITZERMU',' option string',/7x,'in one of the input file titles.',' Will evaluate the Pitzer',/7x,'equations using the',' lambda-mu format. Implied psi coefficients',/7x,'not',' on the supporting data file will not be effectively',' treated',/7x,'as having zero values. Rather, the',' corresponding mu coefficients',/7x,'will be treated as',' having zero values. In effect, an implied',/7x,'psi',' coefficient value will include contributions from',' the',/7x,'two cognate Cphi coefficients that appear',' in the psi-mu',/7x,'breakdown equation.')
    end if

    ! Scan the input file title for the USEOLDPITZER75 string.
    ! If found and if iopg(1) = 1 (use Pitzer's equations), use
    ! the older approximation given by Pitzer (1975) for calculating
    ! the higher-order electrostatic terms.
    qpit75 = .false.

    if (iopg(1) .eq. 1) then
        do n = 1,ntitl
            ux80 = utitl(n)
            call locase(ux80)
            i = index(ux80,'useoldpitzer75')

            if (i .gt. 0) then
                qpit75 = .true.
                go to 170
            end if
        end do

170 continue
    end if

    if (qpit75) then
        write (noutpt,1193)
        write (nttyo,1193)
1193 format(/' * Warning - (EQ3NR/eq3nr) Found the USEOLDPITZER75',' option string',/7x,'in one of the input file titles.',' Will evaluate the Pitzer',/7x,'equations using the',' old Pitzer (1975) approximation for',/7x,'higher order',' electrostatic terms, not the later Harvie (1981)',/7x,'approximation that is now used nearly universally',' to evaluate',/7x,'the Pitzer equations.')
    end if

    ! Scan the input file title for the WRITEPITZERJTABLES string.
    ! If found and if iopg(1) = 1 (use Pitzer's equations), calculate
    ! and output tables of the J(x) and J'(x) functions. This is a
    ! one-time only option. It is not carried forward on the PICKUP
    ! file.
    qcwrpj = .false.

    if (iopg(1) .eq. 1) then
        do n = 1,ntitl
            ux80 = utitl(n)
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
1194 format(/' * Note - (EQ3NR/eq3nr) Found the WRITEPITZERJTABLES',' option string',/7x,'in the input file title.',' Will calculate the Pitzer',/7x,"J(x) and J'(x) functions",' for higher order electrostatic terms and',/7x,'write output',' tables for both the Pitzer (1975) and Harvie (1981)',/7x,'approximations.')
    end if

    ! If the phase part of certain species names read from the
    ! input file is blank, make this part 'Aqueous solution'.
    do n = 1,nsbswt
        if (usbsw(1,n)(25:48) .eq. ublk24(1:24)) then
            usbsw(1,n)(25:48) = uaqsln(1:24)
        end if

        if (usbsw(2,n)(25:48) .eq. ublk24(1:24)) then
            usbsw(2,n)(25:48) = uaqsln(1:24)
        end if
    end do

    do nbi = 1,nbti
        if (uspeci(nbi)(25:48) .eq. ublk24(1:24)) then
            uspeci(nbi)(25:48) = uaqsln(1:24)
        end if
    end do

    do n = 1,nobswt
        if (uobsw(1,n)(25:48) .eq. ublk24(1:24)) then
            uobsw(1,n)(25:48) = uaqsln(1:24)
        end if

        if (uobsw(2,n)(25:48) .eq. ublk24(1:24)) then
            uobsw(2,n)(25:48) = uaqsln(1:24)
        end if
    end do

    ! Set the redox variables.
    qredox = .true.

    if (irdxc3 .ne. -3) then
        do nbi = 1,nbti
            if (uspeci(nbi)(1:6).eq.'O2(g) ' .and.      uspeci(nbi)(25:48).eq.uaqsln(1:24)) then
                write (noutpt,1100)
                write (nttyo,1100)
1100 format(/' * Note - (EQ3NR/eq3nr) The input line for O2(g)',/7x,'will be ignored because irdxc3 is not equal to -3.')

                jflgi(nbi) = 0
                covali(nbi) = 0.
                go to 102
            end if
        end do

102 continue
    end if

    fo2lg = -99999.
    fo2 = 0.

    if (irdxc3 .eq. 0) then
        fo2lg = fo2lgi
        fo2 = texp(fo2lg)
    else if (irdxc3 .eq. -1) then
        eh = ehi
    else if (irdxc3 .eq. -2) then
        pe = pei
    end if

    ! Normalize the mole fractions of any solid solution components.
    nerr = 0

    do nxi = 1,nxti
        nr1 = ncmpri(1,nxi)
        nr2 = ncmpri(2,nxi)

        xx = 0.

        do nxic = nr1,nr2
            xx = xx + xbari(nxic)
        end do

        if (xx .le. 0.) then
            j2 = ilnobl(usoli(nxi))
            write (noutpt,1110) usoli(nxi)(1:j2)
            write (nttyo,1110) usoli(nxi)(1:j2)
1110 format(/' * Error - (EQ3NR/eq3nr) The solid solution',a,/7x,'was given on the input file with a null composition.')

            nerr = nerr + 1
            go to 104
        end if

        do nxic = nr1,nr2
            xbari(nxic) = xbari(nxic)/xx
        end do

104 continue
    end do

    if (nerr .gt. 0) then
        go to 20
    end if

    ! Check consistency between the activity coefficient option and the
    ! data1 file.
    call cdakey(iopg,nopgmx,noutpt,nttyo,udakey,udatfi)

    ! Get the name of the option for calculating the activity
    ! coefficients of aqueous species.
    call nactop(iopg,nopgmx,noutpt,nttyo,uactop)

    ! Open the pickup file if this is to be used.
    if (iopt(17) .ne. -1) then
        inquire(file=pfile,opened=qop)

        if (.not.qop) then
            call openou(noutpt,nttyo,pfile,'formatted',nrecl,newin)
        end if
    end if

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
    ! set will be redefined. This second stage form is equivalent to
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

    ! Force jflag to default to 27 instead of 30 for the species
    ! O2(aq) and H2(aq).
    call jfloha(jflgi,nbtd,nbti,nbtmax,noutpt,nttyo,nxmdmx,nxmod,ubasp,uspeci,uxmod)

    ! Interpret the basis species information listed on the input file.
    call intbs3(covali,ier,jflag,jflgi,narn1a,narn2a,nbaspd,nbtd,nbti,nbtmax,ndrsrd,ndecsp,noutpt,nrdxsp,nsta,nstmax,nttyo,uspeca,uspeci)

    if (ier .gt. 0) then
        go to 20
    end if

    ! Set jflag to -1 for species that can't appear in the system.
    call jflaux(jflag,nbaspd,nbtd,nbtmax,ndrsd,ndrsmx,ndrsrd,nstmax)

    ! Set up temperature variables.
    tempk = tempc + 273.15
    call gntpr(ntpr,ntprmx,ntprt,tempc,tempcu)

    ! Set species status flags.
    call flgset(axlksd,iopt,jflag,jpflag,jsflag,kxmod,narn1a,narn2a,narxmx,nbaspd,nbtd,nbtmax,ncmpra,ncta,ndrsd,ndrsmx,ndrsrd,noptmx,noutpt,npta,nptmax,nrdxsp,nsta,nstmax,ntpr,ntprmx,nttyo,nxmdmx,nxmod,uphasa,uptypa,uspeca,uxmod)

    ! Check the jpflag array to make sure it is consistent with the
    ! jsflag array.
    call flgchk(jpflag,jsflag,ncmpra,npta,nptmax,nstmax,qclnsa)

    ! Look at each active auxiliary basis species. Write a warning if
    ! any other species in the corresponding dissociation reaction is
    ! not present in the model.
    call bspchk(jsflag,nbaspd,nbtd,nbtmax,ndrsd,ndrsmx,ndrsrd,noutpt,nrdxsp,nstmax,nttyo,uspeca)

    ! Do data array compression. Write working data arrays that
    ! do not include phases and species that are not necessary
    ! for the current problem.
    call cmpdat(adhfs,adhfsd,advfs,advfsd,amu,amua,apx,apxa,aslm,aslma,atwt,atwta,axhfs,axhfsd,axlks,axlksd,axvfs,axvfsd,azero,azeroa,bpx,bpxa,cdrs,cdrsd,cess,cessa,iapxmx,iapxt,iapxta,iaqsla,iaqsln,ibpxmx,ibpxt,ibpxta,ilrn1,ilrn2,imrn1,imrn2,ipbtmx,ipch,ipchmx,ipcv,ipcvmx,insgf,insgfa,iopg,ixrn1,ixrn1a,ixrn2,ixrn2a,jflag,jpfcmx,jpflag,jsflag,jsitex,jsol,jsola,mwtsp,mwtspa,nalpaa,nalpha,napmax,napt,napta,narn1,narn1a,narn2,narn2a,narxmx,nat,natmax,nbasp,nbaspd,nbmap,nbt,nbtd,nbti,nbtmax,nchlor,ncmap,ncmpr,ncmpra,nct,ncta,nctmax,ndecsp,ndrs,ndrsd,ndrsmx,ndrsr,ndrsrd,ness,nessa,nessmx,nessr,nessra,ngrn1,ngrn1a,ngrn2,ngrn2a,ngt,nlrn1,nlrn1a,nlrn2,nlrn2a,nlt,nmrn1,nmrn1a,nmrn2,nmrn2a,nmt,nmut,nmuta,nmutmx,nmux,nmuxa,nopgmx,nslt,nsltmx,nslta,nslx,nslxa,nphasx,npt,npta,nptmax,nsmap,nsta,nst,nstmax,ntf1,ntf1a,ntf1mx,ntf1t,ntf1ta,ntf2,ntf2a,ntf2mx,ntf2t,ntf2ta,ntprmx,nxrn1,nxrn1a,nxrn2,nxrn2a,nxt,nxtmax,palpaa,palpha,qchlor,tf1,tf1a,tf2,tf2a,uelem,uelema,uphasa,uphase,uptypa,uptype,uspec,uspeca,vosp0,vosp0a,zchar,zchara)

    do ne = 1,netmax
        cgexp(ne) = 0.
    end do

    do ns = 1,nstmax
        acflg(ns) = 0.
        acflgo(ns) = 0.
        act(ns) = 0.
        conc(ns) = 0.
        mosp(ns) = 0.
        xbar(ns) = 0.
        zchsq2(ns) = 0.
        zchcu6(ns) = 0.
    end do

    av = -99999.
    call initav(conclg,nstmax,av)
    call initav(losp,nstmax,av)
    call initav(actlg,nstmax,av)
    call initav(xbarlg,nstmax,av)

    do np = 1,nptmax
        moph(np) = 0.
    end do

    av = -99999.
    call initav(loph,nptmax,av)

    do kcol = 1,kmax
        zvec1(kcol) = 0.
    end do

    av = -99999.
    call initav(zvclg1,kmax,av)

    uv = 'conc'
    call initcv(ulbeta,kmax,uv)
    call initcv(uldel,kmax,uv)

    do ne = 1,netmax
        kern1(ne) = 0
        kern2(ne) = 0
    end do

    do ne = 1,netmax
        egexpa(ne) = 0.
        egexpc(ne) = 0.
    end do

    nmax = jetmax*netmax
    call initiz(jern1,nmax)
    call initiz(jern2,nmax)
    call initiz(ngext,nmax)
    call initaz(egexjc,nmax)
    call initaz(egexjf,nmax)
    call initaz(mgext,nmax)

    nmax = ietmax*jetmax*netmax
    call initaz(cegexs,nmax)
    call initaz(cpgexs,nmax)
    call initaz(mrgexs,nmax)
    call initiz(ngexro,nmax)
    call initiz(ngexso,nmax)
    call initiz(ngexsa,nmax)
    call initaz(egexs,nmax)

    nmax = ietmax*jetmax*nertmx
    call initaz(egersi,nmax)
    call initaz(xgersi,nmax)

    nmax = ketmax*netmax
    call initiz(kgexsa,nmax)
    call initaz(egexw,nmax)

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

    ! Look at each active auxiliary basis species. Change the
    ! corresponding log K polynomial coefficients so that log K is
    ! fixed at a value of -9999999. if any other species in the
    ! corresponding dissociation reaction is not in the model.
    call bsplkp(axlks,narxmx,nbasp,nbt,nbtmax,ndrs,ndrsmx,ndrsr,nstmax,ntprmx)

    ! Initialize the mole fractions and activities of pure solids and
    ! liquids.
    do ns = nlrn1,nlrn2
        act(ns) = 1.
        actlg(ns) = 0.
        xbar(ns) = 1.
        xbarlg(ns) = 0.
    end do

    do ns = nmrn1,nmrn2
        act(ns) = 1.
        actlg(ns) = 0.
        xbar(ns) = 1.
        xbarlg(ns) = 0.
    end do

    ! Load coefficients for interpreting the input alkalinity, if any.
    ntfxt = ntf1t
    ntfxmx = ntf1mx

    do n = 1,ntf1t
        ntfx(n) = ntf1(n)
        tfx(n) = tf1(n)
    end do

    ! Interpret input file directives to create generic ion-exchange
    ! phases and species.
    call intexi(al10,axhfs,axlks,axvfs,cegexs,cess,cdrs,cgexj,cpgexs,egexjf,iern1,iern2,ietmax,jern1,jern2,jetmax,jflag,jgext,jpflag,jsflag,jsitex,kern1,kern2,ketmax,kgexsa,mwtges,mwtsp,narn1,narn2,narxmx,narxt,nbasp,nbt,nbtmax,ncmpr,ndrs,ndrsmx,ndrsr,ness,nessmx,nessr,netmax,net,nern1,nern2,ngexro,ngexrt,ngexsa,ngexso,ngext,noutpt,nphasx,npt,nptmax,nst,nstmax,ntprt,ntprmx,nttyo,nvetmx,rconst,tgexp,ugexj,ugexmo,ugexmv,ugexp,ugexr,ugexs,ugexsr,uhfgex,uphase,uspec,uvfgex,uxkgex,xhfgex,xlkgex,xvfgex,zchar,zgexj)

    if (ier .gt. 0) then
        go to 20
    end if

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

    ! Alter any log K values as directed by the input file.
    if (nxmod .gt. 0) then
        rtcnst = 0.001*rconst*tempk
        afcnst = al10*rtcnst
        call alters(afcnst,apresg,axlks,cdrs,kxmod,narxmx,narxt,ndrs,ndrsmx,ndrsr,noutpt,npt,nptmax,nst,nstmax,ntpr,ntprmx,nttyo,nxmdmx,nxmod,tempc,uphase,uspec,uxmod,xlkmod)
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

    ! Interpret input file data which defines the indices of species
    ! required to evaluate such constraints as a mean activity or a
    ! heterogeneous equilibrium placed on a basis species.
    call intnsp(coval,covali,ier,jflag,narn1,narn2,nbasp,nbt,nbti,nbtmax,nchlor,ncosp,ndecsp,nhydr,noutpt,nst,nstmax,nttyo,ucospi,uspec)

    if (ier .gt. 0) then
        go to 20
    end if

    if (neti .gt. 0) then
        ! Interpret input file data for concentrations and compositions
        ! of generic ion exchange phases.
        call intge3(cgexp,cgexpi,ier,iern1,iern2,ietmax,jern1,jern2,jetmax,jgext,jgexti,net,neti,netmax,ngexpi,ngexti,noutpt,nptmax,nstmax,nttyo,ugexj,ugexji,ugexp,ugexpi,ugexsi,uphase,uspec,xbar,xbarlg,xgexsi)

        if (ier .gt. 0) then
            go to 20
        end if
    end if

    if (net .gt. 0) then
        ! Set default values for the concentrations (mol/kg.H2O) of
        ! generic exchange phases. These values represent trace amounts.
        ! Currently, the default concentration of an exchanger phase
        ! is 1.e-12 molal.
        do ne = 1,net
            if (cgexp(ne) .le. 0.) then
                cgexp(ne) = 1.e-12
            end if
        end do

        ! Map the concentrations of generic exchange phases into the
        ! coval array. Set correpsonding values in the mtb, moph, and
        ! loph arrays.
        do nb = 1,nbt
            ns = nbasp(nb)

            if (ns.ge.nern1 .and. ns.le.nern2) then
                np = nphasx(ns)
                ne = np - iern1 + 1
                coval(nb) = cgexp(ne)
                mtb(nb) = cgexp(ne)
                moph(np) = cgexp(ne)
                loph(np) = tlg(moph(np))
            end if
        end do
    end if

    ! Interpret input file data for compositions of solid solutions.
    ! file.
    if (iopt(4).ge.1 .and. nxti.gt.0) then
        call intinx(ier,ixrn1,ixrn2,ncmpr,ncmpri,noutpt,npnxp,nptmax,nstmax,nttyo,nxicmx,nxti,nxtimx,umemi,uphase,usoli,uspec,xbar,xbari,xbarlg)

        if (ier .gt. 0) then
            go to 20
        end if
    end if

    ! Check to see that no species in the strict basis has a jflag
    ! value of 30 at this point.
    nerr = 0

    do nb = 1,nbt
        ns = nbasp(nb)
        nt = ndrsr(2,ns) - ndrsr(1,ns) + 1

        if (nt .lt. 2) then
            if (jflag(ns) .eq. 30) then
                ! Calling sequence substitutions:
                !   uspec(ns) for unam48
                call fmspnm(jlen,uspec(ns),uspn56)
                write (noutpt,1240) uspn56(1:jlen)
                write (nttyo,1240) uspn56(1:jlen)
1240 format(/' * Error - (EQ3NR/eq3nr) The strict basis species'      /7x,a,' has a jflag value of 30. This implies',/7x,'that this species is to be treated as a dependent',' species. However,',/7x,"a strict basis species can't",' be treated in this manner.')

                nerr = nerr + 1
            end if
        end if
    end do

    if (nerr .gt. 0) then
        go to 20
    end if

    ! Redefine the 'd' set of reactions and attendant data to
    ! match the ordinary set as it presently exists. The ordinary
    ! set will be futher subjected to eliminations from the active
    ! basis set and ordinary basis switching.
    call cdrssd(adhfs,adhfsd,advfs,advfsd,axhfs,axhfsd,axlks,axlksd,axvfs,axvfsd,cdrs,cdrsd,ipch,ipchmx,ipcv,ipcvmx,narxmx,nbasp,nbaspd,nbtmax,ndrs,ndrsd,ndrsmx,ndrsr,ndrsrd,nstmax,ntprmx)
    nbtd = nbt

    if (iopg(1) .eq. 1) then
        ! Build the S-lambda index arrays nsxi and nsxx.
        call bdslx(narn1,narn2,natmax,noutpt,nslt,nsltmx,nslx,nsxi,nsxx,nsxmax,nttyo)

        ! Build the mu index arrays nmxi and nmxx.
        call bdmlx(narn1,narn2,natmax,nmut,nmutmx,nmux,nmxi,nmxmax,nmxx,noutpt,nttyo)

        if (iopr(10) .gt. 0) then
            ! Write tables concerning Pitzer coefficients.
            call ptztab(iopr,narn1,narn2,natmax,nmutmx,nmux,nmxi,nmxmax,nmxx,noprmx,noutpt,nsltmx,nslx,nstmax,nsxi,nsxmax,nsxx,uspec)
        end if

        ! Write warnings for species lacking Pitzer coefficients.
        call ptzchk(narn1,narn2,natmax,nmxi,noutpt,nstmax,nsxi,nttyo,uspec)

        ! Transform conventional mu data to corresponding C, psi,
        ! and zeta data (data originally defined in mu form is
        ! not affected).
        if (qhawep) then
            call rc3ocf(amu,jpfcmx,ifcphi1,ifcphi2,ifnnn,ifn2n,ifpsi1,ifpsi2,ifzeta,ilcphi1,ilcphi2,ilnnn,iln2n,ilpsi1,ilpsi2,ilzeta,iodb,nmux,nmut,nmutmx,nodbmx,noutpt,nstmax,nttyo,uspec,zchar)
        end if
    end if

    ! Compute the thermodynamic parameters that are functions of
    ! temperature and pressure.
    call evdata(aadh,aadhh,aadhv,aaphi,abdh,abdhh,abdhv,abdot,abdoth,abdotv,adadhh,adadhv,adbdhh,adbdhv,adbdth,adbdtv,adh,adhfe,adhh,adhv,adhfs,adhfsd,advfe,advfs,advfsd,afcnst,al10,amu,aslm,aphi,aprehw,apresg,apresh,apx,avcnst,axhfe,axhfs,axhfsd,axlke,axlks,axlksd,axvfe,axvfs,axvfsd,bdh,bdhh,bdhv,bdot,bdoth,bdotv,dadhh,dadhv,dbdhh,dbdhv,dbdth,dbdtv,dhfe,dhfs,dhfsd,dvfe,dvfs,dvfsd,ehfac,farad,iapxmx,iktmax,iopg,iopt,ipbtmx,ipch,ipchmx,ipcv,ipcvmx,ixrn1,ixrn2,jpfcmx,jptffl,jsol,narxmx,narxt,narxth,ncmpr,nmut,nmutmx,nopgmx,noptmx,noutpt,nptmax,nslt,nsltmx,nst,nstmax,ntpr,ntprmx,nttyo,nxt,nxtmax,pmu,prehw,presg,presh,press,pslamn,rconst,rcnstv,rtcnst,tempc,tempk,uphase,uspec,wfac,xhfe,xhfs,xhfsd,xlke,xlks,xlksd,xvfe,xvfs,xvfsd)

    if (jpres3 .eq. 0) then
        ! The pressure corresponds to the data file reference presssure
        ! curve.
        press = presg
    else if (jpres3 .eq. 1) then
        ! The pressure corresponds to the 1.013-bar/steam-saturation
        ! presssure curve.
        press = presh
    else if (jpres3 .eq. 2) then
        ! Constant pressure.
        continue
    else
        ! Error.
        write (noutpt,1242) jpres3
        write (nttyo,1242) jpres3
1242 format(/' * Error - (EQ3NR/eq3nr) The pressure option flag',' (jpres3)',/7x,'has an unknown value of ',i2,'.')
    end if

    ! Stop if the pressure isn't greater than 0.
    if (press .le. 0) then
        write (noutpt,1244) press
        write (nttyo,1244) press
1244 format(/' * Error - (EQ3NR/eq3nr) The pressure must be',' greater than zero.',/7x,'The current pressure is ',1pg12.5,' bars.')

        stop
    end if

    ! Stop if press is greater than presmx.
    if (press .gt. presmx) then
        write (noutpt,1246) press,presmx
        write (nttyo,1246) press,presmx
1246 format(/' * Error - (EQ3NR/eq3nr) The calculated pressure is',' ',1pg12.5,' bars,',/7x,"greater than the code's built-in",' maximum value of ',1pg12.5,' bars.',/7x,'Other limits',' associated with the supporting data file may also apply.')

        stop
    end if

    ! Note if press is less than the 1.013-bar/steam-saturation curve
    ! value.
    if ((press - presh) .lt. -1.e-4) then
        write (noutpt,1248) press,presh
        write (nttyo,1248) press,presh
1248 format(/' * Note - (EQ3NR/eq3nr) The calculated pressure is',' ',1pg12.5,' bars,',/7x,'less than the 1.013-bar',' steam-saturation curve value',/7x,'of ',g12.5,' bars.')
    end if

    ! If the pressure is nearly identical to the data file reference
    ! pressure curve value, set it equal to that value.
    if (abs(press - presg) .le. 1.e-4) then
        press = presg
    end if

    ! Check for thermodynamic pressure corrections.
    dp = press - presg

    if (abs(dp) .ge. 1.e-4) then
        if (ipcv .lt. 0) then
            ! There are no data to support needed pressure corrections.
            write (noutpt,1250) press,presg,dp
            write (nttyo,1250) press,presg,dp
1250 format(/' * Warning - (EQ3NR/eq3nr) The supporting data file',/7x,'contains no data to support making thermodynamic',/7x,'pressure corrections. No such corrections will be made.',/7x,'The current pressure is ',1pg12.5,' bars, the standard',/7x,'grid pressure is ',g12.5,' bars, and the pressure',/7x,'difference is ',g12.5,' bars.')
        else
            if (press .le. 0.) then
                ! The pressure is zero or negative.
                write (noutpt,1260) press
                write (nttyo,1260) press
1260 format(/' * Error - (EQ3NR/eq3nr) The pressure must be',/7x,'greater than zero. The current pressure is ',1pg12.5,' bars.')

                stop
            end if

            if (abs(dp) .gt. prehw) then
                ! The pressure is outside the recommended envelope.
                pxu = presg + prehw
                pxl = presg - prehw

                if (pxl .le. 0.) then
                    pxl = 0.
                end if

                write (noutpt,1270) press,pxl,pxu,tempc
                write (nttyo,1270) press,pxl,pxu,tempc
1270 format(/' * Warning - (EQ3NR/eq3nr) The current pressure',/7x,'of ',1pg12.5,' bars is outside the recommended',/7x,'pressure envelope of ',g12.5,' to ',g12.5,' bars',/7x,'at ',0pf6.2,' C.')
            end if

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

    ! Check default redox constraints.
    if (irdxc3 .eq. -3) then
        continue
    else if (irdxc3 .eq. -2) then
        continue
    else if (irdxc3 .eq. -1) then
        continue
    else if (irdxc3 .eq. 0) then
        continue
    else if (irdxc3 .eq. 1) then
        continue
    else
        ! Error.
        write (noutpt,1280) irdxc3
        write (nttyo,1280) irdxc3
1280 format(/' * Error - (EQ3NR/eq3nr) The default redox option',/7x,'flag (irdxc3) has an unknown value of ',i2,'.')

        go to 20
    end if

    ! Make sure that aqueous O2(g) was included in the list of active
    ! basis species on the data file.
    no2gai = 0

    do nbi = 1,nbti
        nb = ndecsp(nbi)

        if (nb .gt. 0) then
            ns = nbaspd(nb)

            if (ns .eq. no2gaq) then
                no2gai = nbi
                go to 110
            end if
        end if
    end do

110 continue

    if (irdxc3 .ne.-3 .and. jflag(no2gaq).ne.0) then
        write (noutpt,1300)
1300 format(/' * Note - (EQ3NR/eq3nr) Have conflicting redox',' options:')

        write (noutpt,1310) irdxc3,jflag(no2gaq)
1310 format(9x,'irdxc3 = ',i2,' overrides jflag(O2(g)) = ',i2,/)

        jflag(no2gaq) = 0

        ! Calling sequence substitutions:
        !   no2gaq for ns
        io2gaq = nbasis(nbasp,nbt,nbtmax,no2gaq)
        coval(io2gaq) = -99999.
        ncosp(io2gaq) = 0

        if (no2gai .gt. 0) then
            jflgi(no2gai) = 0
            covali(no2gai) = -99999.
            ucospi(no2gai)(1:48) = ' '
        end if
    end if

    j2 = ilnobl(uredox)

    if (irdxc3.ne.1 .and. uredox(1:5).ne.'None ' .and. j2.gt.0) then
        write (noutpt,1300)
        write (noutpt,1320) irdxc3,uredox(1:j2)
1320 format(9x,'irdxc3 = ',i2,' overrides uredox = ',a,/)

        uredox(1:24) = 'None'
    end if

    ! If irdxc3 .ge. 1, find the index of the species whose name is
    ! uredox.
    nredox = 0

    if (irdxc3 .eq. 1) then
        j2 = ilnobl(uredox)

        if (j2.le.0 .or. uspec(ns)(1:4).eq.'None') then
            write (noutpt,1330)
            write (nttyo,1330)
1330 format(/' * Error - (EQ3NR/eq3nr) The default redox state is',/7x,'controlled by a couple, but the auxiliary basis species',/7x,"which defines the couple (uredox) isn't specified on",/7x,'the input file.')

            go to 20
        else
            do nb = 1,nbt
                ns = nbasp(nb)

                if (uspec(ns)(1:24) .eq. uredox(1:24)) then
                    nredox  = ns
                    go to 120
                end if
            end do

            write (noutpt,1340) uredox(1:j2)
            write (nttyo,1340) uredox(1:j2)
1340 format(/' * Error - (EQ3NR/eq3nr) The default redox state is',/7x,'controlled by a couple, but the species ',a,', whose',/7x,"associated reaction defines this couple, isn't",/7x,'in the basis set, as is required.')

            go to 20
        end if
    end if

120 continue

    ! Check the total dissolved solutes option.
    if (itdsf3 .eq. 0) then
        continue
    else if (itdsf3 .eq. 1) then
        continue
    else
        ! Error.
        write (noutpt,1342) itdsf3
        write (nttyo,1342) itdsf3
1342 format(/' * Error - (EQ3NR/eq3nr) The total dissolved solutes',/7x,'option flag (itdsf3) has an unknown value of ',i2,'.')

        go to 20
    end if

    ! Check the electrical balancing option.
    if (iebal3 .eq. 0) then
        continue
    else if (iebal3 .eq. 1) then
        continue
    else
        ! Error.
        write (noutpt,1344) iebal3
        write (nttyo,1344) iebal3
1344 format(/' * Error - (EQ3NR/eq3nr) The electrical balancing',/7x,'option flag (iebal3) has an unknown value of ',i2,'.')

        go to 20
    end if

    ! Set default values.
    qrho = rho .gt. 0.
    call dfaltx(itermx,rho,scamas,tdspkg,tdspl,tolbt,toldl,tolspf)

    if (net .gt. 0) then
        ! Echo a table for the generic ion exchangers, describing the
        ! setup of species, reactions, and corresponding thermodynamic
        ! data.
        call echgex(axlks,cdrs,cgexj,iern1,iern2,jern1,jern2,jetmax,jgext,jpflag,jsflag,narxmx,narxt,ndrs,ndrsmx,ndrsr,netmax,noutpt,nptmax,ntprmx,ntprt,nstmax,press,tempc,ugexj,ugexmo,uphase,uspec,xlks)
    end if

    ! Print an echo of the processed input.
    call echox(azero,cdrs,covali,eh,fo2lg,iebal3,iodb,iopg,iopr,iopt,irdxc3,itdsf3,itermx,ixrn1,ixrn2,jflag,jflgi,jpres3,narn1,narn2,nat,nata,natmax,nbt,nbta,nbti,nbtmax,nct,ncta,nctmax,ncmpr,ncosp,ndecsp,ndrs,ndrsmx,ndrsr,nelect,ngt,ngta,ngtmax,nhydr,nhydx,njfmax,nlt,nlta,nltmax,nmt,nmta,nmtmax,nodbmx,nopgmx,noprmx,noptmx,noutpt,no2gaq,npt,npta,nptmax,nredox,nst,nsta,nstmax,nxt,nxta,nxti,nxtmax,pe,presg,press,rho,scamas,tdspkg,tdspl,tempc,tolbt,toldl,tolspf,uactop,ucospi,uebal,ujflls,uphase,uredox,uspec,uspeci,xbar)

    if (qcwrpj) then
        ! Write tables of the Pitzer J(x) and J'(x) functions.
        call cwrpjt(noutpt)
    end if

    ! Convert input concentration data not on the molal scale to
    ! that scale. Convert pe, if input, to Eh.
    call setup(coval,eh,ehfac,ier,irdxc3,itdsf3,jflag,mwtsp,narn1,nbaspd,nbt,nbtmax,noutpt,nstmax,nttyo,pe,rho,tdspkg,tdspl,uspec)

    if (ier .gt. 0) then
        go to 20
    end if

    ! Echo input data after modification.
    write (noutpt,1400)
1400 format(/21x,'--- Modified Input Constraints ---',//5x,'Species',20x,'coval   jflag   Type of Input',/)

    do nb = 1,nbt
        ns1 = nbasp(nb)

        if (jsflag(ns1) .le. 0) then
            jfl = jflag(ns1)
            j2 = ilnobl(ujflls(jfl))

            if (jfl.eq.17 .or. jfl.eq.18) then
                nse = ncosp(nb)
                write (noutpt,1410) uspec(ns1),coval(nb),jfl,ujflls(jfl)(1:j2)
1410 format(2x,a24,2x,1pe12.5,2x,i2,2x,a)

                j3 = ilnobl(uspec(nse)(1:24))
                write (noutpt,1420) uspec(nse)(1:j3)
1420 format(39x,'Counterion= ',a)
            else if (jfl .eq. 25) then
                nse = ncosp(nb)
                write (noutpt,1422) uspec(ns1),jfl,ujflls(jfl)(1:j2)
1422 format(2x,a24,16x,i2,2x,a)

                j3 = ilnobl(uspec(nse)(1:24))
                ux24 = uspec(nse)(25:48)
                j4 = ilnobl(ux24)
                write (noutpt,1425) uspec(nse)(1:j3),ux24(1:j4)
1425 format(46x,'Species= ',a,/48x,'Phase= ',a)

                ! Calling sequence substitutions:
                !   noutpt for nf
                !   nse for ns
                call prreac(cdrs,ndrs,ndrsmx,ndrsr,noutpt,nse,nstmax,uspec)
                write (noutpt,1426)
1426 format(1x)
            else if (jfl.eq.27 .or. jfl.eq.30) then
                write (noutpt,1422) uspec(ns1),jfl,ujflls(jfl)(1:j2)
            else if (ns1 .eq. no2gaq) then
                write (noutpt,1427) uspec(ns1),coval(nb),jfl
1427 format(2x,a24,2x,1pe12.5,2x,i2,2x,'Log fO2')
            else
                if (jfl .ge. 0) then
                    write (noutpt,1410) uspec(ns1),coval(nb),jfl,ujflls(jfl)(1:j2)
                else
                    write (noutpt,1430) uspec(ns1)
1430 format(2x,a24,16x,'Not present in the model')
                end if
            end if
        end if
    end do

    write (noutpt,1440)
1440 format(1x)

    ! Find the index of the species to be adjusted for electrical
    ! balance, if any.
    call intieb(iebal,iebal3,ier,jsflag,nbasp,nbt,nbtmax,noutpt,nstmax,nttyo,uebal,uspec,zchar)

    if (ier .gt. 0) then
        go to 20
    end if

    if (uebal(1:5) .ne. 'None ') then
        j2 = ilnobl(uebal)
        write (noutpt,1450) uebal(1:j2)
1450 format(/' Electrical balance will be achieved by adjusting',/'   the concentration of ',a,'.',/)
    end if

    ! Write a list of inactive species.
    write (noutpt,1470)
1470 format(/21x,'--- Inactive Species ---',/)

    k = 0

    do ns = 1,nst
        if (jsflag(ns) .eq. 1) then
            if (ns .ne. no2gaq) then
                k = k + 1
                j2 = ilnobl(uspec(ns)(1:24))
                write (noutpt,1480) uspec(ns)(1:j2)
1480 format(4x,a)
            end if
        end if
    end do

    if (k .le. 0) then
        write (noutpt,1490)
    end if

1490 format(4x,'None',/)

    ! Call EQLIB/elim.f to rewrite the reaction equations (cdrs/ndrs/
    ! ndrsr arrays) so that auxiliary basis variables with jflag = 30
    ! are eliminated from the active basis set.
    qelim = .false.

    do nb = 1,nbt
        nse = nbasp(nb)
        nt = ndrsr(2,nse) - ndrsr(1,nse) + 1

        if (nt.ge.2 .and. jflag(nse).eq.30) then
            if (jsflag(nse) .gt. 0) then
                coval(nb) = 0.
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

    ! Set up the ixbasp and cjbasp arrays. The former is a flag
    ! array, each member of which denotes whether the
    ! thermodynamic activity of the corresponding basis species
    ! is defined in terms of molality (= 0) or mole fraction (= 1).
    ! The cjbasp array contains any site stoichiometric factors
    ! associated with the operational basis species.
    call gibasp(cgexj,cjbasp,iern1,ixbasp,jern1,jern2,jetmax,jgext,narn1,narn2,nbasp,nbt,nbtmax,nern1,nern2,netmax,nphasx,nstmax)

    ! Call lambda to compute activity coefficients for components of
    ! solid solutions for which compositions were specified on the
    ! input file.
    if (iopt(4).ge.1 .and. nxti.gt.0) then
        do np = ixrn1,ixrn2
            ! Calling sequence substitutions:
            !   acflg for acflgc
            call lambda(acflg,afcnst,bpx,ibpxmx,ibpxt,iktmax,ixrn1,ixrn2,jsol,ncmpr,noutpt,np,nptmax,nstmax,nttyo,nxtmax,wfac,xbar,xbarlg,uphase,uspec)
        end do
    end if

    ! Find which of HCO3- and CO3-- is in the basis set. If both are
    ! in the basis set (potentially a bad situation), take the first.
    ncarb = 0

    do nb = 1,nbt
        ns = nbasp(nb)

        if (uspec(ns)(1:6).eq.'HCO3- ' .or.    uspec(ns)(1:6).eq.'CO3-- ') then
            icarb = nb
            ncarb = ns
            go to 210
        end if
    end do

210 continue

    ! Check the input file for inconsistencies and errors.
    ! Here ier accumulates the number of errors caught.
    call chkinx(cdrs,coval,ier,irdxc3,jflag,jsflag,narn1,narn2,nbasp,nbt,nbtmax,ncosp,ndrs,ndrsmx,ndrsr,nelect,nhydr,nhydx,noutpt,no2gaq,nredox,nstmax,nttyo,tempc,ucospi,uredox,uspec,zchar)

    if (ier .ge. 1) then
        write (noutpt,1510) ier
        write (nttyo,1510) ier
1510 format(/' * Error - (EQ3NR/eq3nr) The input has failed ',i4,' checks.')

        go to 20
    end if

    ! Save the input value of the alkalinity, if any.
    alki = 0.

    if (ncarb .gt. 0) then
        jfl = jflag(ncarb)

        if (jfl .eq. 7) then
            alki = coval(icarb)
        end if
    end if

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
        ilevel = iopr(2)

        ! Calling sequence substitutions:
        !   noutpt for nf
        call echolk(axlks,cdrs,ilevel,jsflag,narxmx,ndrs,ndrsmx,ndrsr,noutpt,nst,ntprmx,nstmax,presg,tempc,uspec,xlks)
    end if

    ! Initialize the jcsort, jjsort, jssort, and jgsort arrays by
    ! setting each element equal to its index.
    call initii(jcsort,nst)
    call initii(jjsort,nst)
    call initii(jssort,nst)
    call initii(jgsort,ngt)

    ! Set up the cdrw array. This provides a fast way to get the
    ! reaction coefficient of H2O (Aqueous solution) in any reaction.
    call gcdrw(cdrs,cdrw,narn1,ndrs,ndrsmx,ndrsr,nst,nstmax)

    ! The cdrtw array is not actually used by EQ3NR, only EQ6.
    ! It is included here for the sake of consistency.
    call gcdrtw(cdrs,cdrtw,narn1,narn2,ndrs,ndrsmx,ndrsr,nelect,no2gaq,nst,nstmax)

    write (noutpt,1520)
1520 format(/' - - BEGIN ITERATIVE CALCULATIONS  - - - - - - - - - -',' - - - - - - - - - - - -',/)

    ! Set up the matrix structure (key index array iindx1) for hybrid
    ! Newton-Raphson iteration and compute starting values of the matrix
    ! variables.
    call arrset(aamatr,abar,acflg,acflgo,act,actlg,adh,adhh,adhv,adhfs,adhfsx,advfs,advfsx,afcnst,alpha,al10,amtb,aphi,avcnst,azero,a3bar,a3bars,axhfs,axhfsx,axlks,axlksx,axvfs,axvfsx,bacfmx,bbig,bdh,bdhh,bdhv,bdot,bdoth,bdotv,beta,betamx,bfac,bgamx,bneg,bpx,bsigmm,bfje,bfxi,cco2,cdrs,cdrsx,cdrtw,cdrw,cegexs,cgexj,cjbasp,cnufac,conc,conclg,coval,cpgexs,csts,delam,delvec,dgpit,dhfs,dlogxw,dpelm,dpslm,dselm,dvfs,efac,egexjc,egexjf,egexs,eh,ehfac,elam,eps100,fje,fjeo,fo2,fo2lg,fsort,fugac,fugalg,fxi,fxio,gmmatr,gpit,ibetmx,ibpxt,ibswx,iction,iebal,ielam,iern1,iern2,ifcphi1,ifcphi2,ifnnn,ifn2n,ifpsi1,ifpsi2,ifrn1,ifrn2,ifzeta,igas,igstak,iindx1,ilcphi1,ilcphi2,ilnnn,iln2n,ilpsi1,ilpsi2,ilrn1,ilrn2,ilzeta,imrn1,imrn2,insgf,iopg,iodb,iopt,ipch,ipcv,ipivot,ipndx1,irdxc3,istack,ixbasp,ixrn1,ixrn2,izmax,jcsort,jern1,jern2,jflag,jgext,jgsort,jgstak,jjndex,jjsort,jpflag,jsflag,jsitex,jsol,jssort,jstack,ka1,kat,kbt,kct,kction,kdim,kebal,kelect,ker,ke1,ket,khydr,kkndex,km1,ko2gaq,kwater,kx1,kxt,loph,losp,lsort,mgext,moph,mosp,mrgexs,mtb,nalpha,napt,narn1,narn2,narxt,nbasp,nbaspd,nbaspx,nbt,nbtd,nbti,nbw,nchlor,ncmpr,ncosp,nct,ndecsp,ndrs,ndrsx,ndrsr,ndrsrd,ndrsrx,nelect,nern1,nern2,net,nfac,ngexsa,ngext,ngrn1,ngrn2,ngt,nhydr,nhydx,nmut,nmux,nmxi,nmxx,noutpt,no2gaq,nphasx,npt,nredox,nslt,nslx,nst,nsts,nstsr,nsxi,nsxx,ntfx,ntfxt,ntpr,nttyo,omega,omeglg,palpha,pe,pelm,pmu,presg,press,pslamn,pslm,qbassw,qchlor,qhawep,qpit75,qredox,q6mode,rhsvec,selm,sigmam,sigmmo,smp100,tempc,tempk,tfx,ubacmx,ubbig,ubgamx,ubneg,ubetmx,ucospi,ugexj,ugexmo,ujflls,uphase,uspec,uzvec1,weight,wfac,xbar,xbarlg,xbarw,xbarwc,xbrwlc,xbrwlg,xhfs,xlke,xlks,xvfs,zchar,zchsq2,zchcu6,zgexj,zvclg1,zvec1)

    if (ker .ge. 2) then
        go to 20
    end if

    ! Find the matrix position of HCO3- or CO3--, if any.
    kcarb = 0

    if (ncarb .gt. 0) then
        do kcol = 1,kbt
            nb = iindx1(kcol)
            ns = nbasp(nb)

            if (ns .eq. ncarb) then
                kcarb = kcol
                go to 220
            end if
        end do

220 continue
    end if

    ! Calculate the speciation-solubility model.
    write (noutpt,1600)
    write (nttyo,1600)
1600 format(/,' Starting hybrid Newton-Raphson iteration.',/)

    call newton(aamatr,abar,acflg,acflgo,act,actlg,actwlc,adh,adhh,adhv,afcnst,alpha,al10,amtb,aphi,azero,a3bar,a3bars,bacfmx,bbig,beta,betamx,betao,bdh,bdhh,bdhv,bdot,bdoth,bdotv,bfje,bfxi,bgamx,bneg,bpx,bsigmm,cco2,cdrs,cdrtw,cdrw,cegexs,cgexj,cjbasp,cnufac,conc,conclg,coval,cpgexs,csts,delam,delmax,delvco,delvec,dgpit,dlogxw,dpelm,dpslm,dselm,egexjc,egexjf,egexs,eh,ehfac,elam,eps100,fje,fjeo,fo2,fo2lg,fsort,fugac,fugalg,fxi,fxio,gmmatr,gpit,ibpxt,idelmx,iebal,ielam,ier,iern1,iern2,ifcphi1,ifcphi2,ifnnn,ifn2n,ifpsi1,ifpsi2,ifrn1,ifrn2,ifzeta,igas,igstak,iindx1,ilcphi1,ilcphi2,ilnnn,iln2n,ilpsi1,ilpsi2,ilrn1,ilrn2,ilzeta,imrn1,imrn2,insgf,iodb,iopg,ipivot,ipndx1,irdxc3,istack,iter,itermx,ixbasp,ixrn1,ixrn2,izmax,jcsort,jern1,jern2,jflag,jgext,jgsort,jgstak,jjsort,jpflag,jsflag,jsitex,jsol,jssort,jstack,kbt,kction,kdim,kelect,khydr,km1,kmt,ko2gaq,kwater,kx1,kxt,loph,losp,lsort,mgext,moph,mosp,mrgexs,mtb,nalpha,napt,narn1,narn2,nbasp,nbt,nbw,nchlor,ncmpr,ncosp,ndrs,ndrsr,nelect,nern1,nern2,net,ngexsa,ngext,ngrn1,ngrn2,ngt,nhydr,nmut,nmux,nmxi,nmxx,noutpt,no2gaq,nphasx,npt,nredox,nslt,nslx,nst,nsts,nstsr,nsxi,nsxx,ntfx,ntfxt,nttyo,omega,omeglg,palpha,pelm,pmu,press,pslamn,pslm,qhawep,qpit75,qredox,q6mode,rhsvec,screwd,screwn,selm,sigmam,sigmmo,tempk,tfx,tolbt,toldl,ubacmx,ubbig,ubetmx,ubgamx,ubneg,ugexj,ugexmo,ulbeta,uldel,uphase,uspec,uzvec1,weight,wfac,xbar,xbarlg,xbarw,xbarwc,xbrwlc,xbrwlg,xlke,xlks,zchar,zchsq2,zchcu6,zgexj,zvclg1,zvec1)

    if (ier .gt. 0) then
        write (noutpt,1610)
        write (nttyo,1610)
1610 format(/' * Error - (EQ3NR/eq3nr) Hybrid Newton-Raphson',' iteration failed')

        ux8a = ' '
        write (ux8a,1620) iter
1620 format(i7)

        call lejust(ux8a)
        j2 = ilnobl(ux8a)

        if (ier .eq. 1) then
            write (noutpt,1630) ux8a(1:j2)
            write (nttyo,1630) ux8a(1:j2)
1630 format(7x,'after ',a,' iterations because a zero matrix',' was encountered. This is',/7x,'probably due to a',' programming error.')
        else if (ier .eq. 2) then
            write (noutpt,1640) ux8a(1:j2)
            write (nttyo,1640) ux8a(1:j2)
1640 format(7x,'after ',a,' iterations because a non-zero,',' computationally singular',/7x,'matrix was encountered.')
        else if (ier .eq. 3) then
            write (noutpt,1650) ux8a(1:j2)
            write (nttyo,1650) ux8a(1:j2)
1650 format(7x,'after ',a,' iterations because the code',' detected that',/7x,'iteration was diverging.')
        else if (ier .eq. 4) then
            write (noutpt,1660) ux8a(1:j2)
            write (nttyo,1660) ux8a(1:j2)
1660 format(7x,'after ',a,' iterations because the maximum',' number of iterations',/7x,'was done.')
        else
            ux8b = ' '
            write (ux8b,1620) ier
            call lejust(ux8b)
            j3 = ilnobl(ux8b)
            write (noutpt,1670) ux8a(1:j2),ux8b(1:j3)
            write (nttyo,1670) ux8a(1:j2),ux8b(1:j3)
1670 format(7x,'after ',a,' iterations because an unknown event',' occurred. The ier',/7x,'error code has the unknown value',' ',a,'. This condition is a',/7x,'programming error.')
        end if

        call ndiagx(delmax,delvec,eps100,idelmx,iebal,iindx1,irdxc3,jflag,kcarb,kebal,khydr,kmax,ko2gaq,nbasp,nbtmax,nhydr,noutpt,nstmax,nttyo,screwd,uspec)
        go to 20
    end if

    write (noutpt,1690) iter
    write (nttyo,1690) iter
1690 format('   Done. Hybrid Newton-Raphson iteration converged in ',i3,' iterations.',/)

    ! Save the jflag value of the species being adjusted for electrical
    ! balance.
    jfleba = 0

    if (iebal .gt. 0) then
        ns = nbaspd(iebal)
        jfleba = jflag(ns)
    end if

    ! Reset jflag values which are not 0 or 30 to one of these values.
    ! Species for which jflag will be set to 0 now have associated
    ! mass balance totals, those for which jflag will be set to 30
    ! do not.
    if (nredox .ge. 1) then
        jflag(nredox) = 0
    end if

    do nb = 1,nbt
        ns = nbasp(nb)
        jflgi(nb) = jflag(ns)

        if (jflgi(nb) .eq. 27) then
            jflgi(nb) = 30
        end if

        if (jflag(ns) .ne. 30) then
            jflag(ns) = 0
        end if
    end do

    ! Get the weights (actually, the masses) of total dissolved
    ! solutes, the solvent, and the solution, and also some related
    ! quantities. Also get the solution density (rhowc) from the
    ! WIPP brine density model if the temperature lies in the range
    ! 20-30 C. This model is a fit to NaCl solution data having
    ! the form:
    !   density (g/L) = a + b x TDS (g/L)
    ! where a = 1000.96 and b = 0.639963. Also get the molarity/
    ! molality conversion ratio (mrmlra, molarity = mrmlra x molality)
    ! and its inverse.
    call gwdenp(adwipp,bdwipp,jcsort,mlmrra,mosp,mrmlra,mwtsp,narn1,narn2,nstmax,qdwipp,rhoc,rhowc,tdsgks,tdsglw,tdspkc,tdsplc,tempc,vosol,wfh2o,wftds,wkgwi,woh2o,wosol,wotds)

    ! Check potential difference between input and calculated
    ! densities.
    if (rho.gt.0. .and. rhoc.gt.0.) then
        axx = abs((rhoc - rho)/rho)

        if (axx .gt. 0.01) then
            write (ux16a,'(1pg12.5)') rhoc
            call lejust(ux16a)
            j1 = ilnobl(ux16a)
            write (ux16b,'(1pg12.5)') rho
            call lejust(ux16b)
            j2 = ilnobl(ux16b)
            write (noutpt,1720) ux16a(1:j1),ux16b(1:j2)
            write (nttyo,1720) ux16a(1:j1),ux16b(1:j2)
1720 format(/' * Warning - (EQ3NR/eq3nr) The calculated density',' of ',a,' g/mL',/7x,'differs from the input file/default',' value of ',a,' g/mL',/7x,'by more than 1%.',' The calculated value will be used',/7x,'in subsequent',' calculations.')
        end if
    end if

    ! Check potential difference between input and calculated
    ! TDS values.
    if (itdsf3 .le. 0) then
        if (tdspkg.gt.0. .and. tdspkc.gt.0.) then
            axx = abs((tdspkg - tdspkc)/tdspkc)

            if (axx .gt. 0.01) then
                write (ux16a,'(1pg12.5)') tdspkc
                call lejust(ux16a)
                j1 = ilnobl(ux16a)
                write (ux16b,'(1pg12.5)') tdspkg
                call lejust(ux16b)
                j2 = ilnobl(ux16b)
                write (noutpt,1730) ux16a(1:j1),ux16b(1:j2)
                write (nttyo,1730) ux16a(1:j1),ux16b(1:j2)
1730 format(/' * Warning - (EQ3NR/eq3nr) The calculated TDS',' of ',a,' mg/kg.sol',/7x,'differs from the input',' file/default value of ',a,' mg/kg.sol',/7x,'by more',' than 1%. The calculated value will be used in',/7x,'subsequent calculations.')
            end if
        else if (tdspkg.le.0. .and. tdspkc.gt.0.) then
            write (ux16a,'(1pg12.5)') tdspkc
            call lejust(ux16a)
            j1 = ilnobl(ux16a)
            write (ux16b,'(1pg12.5)') tdspkg
            call lejust(ux16b)
            j2 = ilnobl(ux16b)
            write (noutpt,1732) ux16a(1:j1),ux16b(1:j2)
            write (nttyo,1732) ux16a(1:j1),ux16b(1:j2)
1732 format(/' * Warning - (EQ3NR/eq3nr) The calculated TDS',' of ',a,' mg/kg.sol',/7x,'differs from the input',' file/default value of ',a,' mg/kg.sol.',/7x,'The calculated value will be used in subsequent',' calculations.')
        end if
    else
        if (tdspl.gt.0. .and. tdsplc.gt.0.) then
            axx = abs((tdspl - tdsplc)/tdsplc)

            if (axx .gt. 0.01) then
                write (ux16a,'(1pg12.5)') tdsplc
                call lejust(ux16a)
                j1 = ilnobl(ux16a)
                write (ux16b,'(1pg12.5)') tdspl
                call lejust(ux16b)
                j2 = ilnobl(ux16b)
                write (noutpt,1740) ux16a(1:j1),ux16b(1:j2)
                write (nttyo,1740) ux16a(1:j1),ux16b(1:j2)
1740 format(/' * Warning - (EQ3NR/eq3nr) The calculated TDS',' of ',a,' mg/L',/7x,'differs from the input',' file value of ',a,' mg/L',/7x,'by more',' than 1%. The calculated value will be used in',/7x,'subsequent calculations.')
            end if
        else if (tdspl.le.0. .and. tdsplc.gt.0.) then
            write (ux16a,'(1pg12.5)') tdsplc
            call lejust(ux16a)
            j1 = ilnobl(ux16a)
            write (ux16b,'(1pg12.5)') tdspl
            call lejust(ux16b)
            j2 = ilnobl(ux16b)
            write (noutpt,1742) ux16a(1:j1),ux16b(1:j2)
            write (nttyo,1742) ux16a(1:j1),ux16b(1:j2)
1742 format(/' * Warning - (EQ3NR/eq3nr) The calculated TDS',' of ',a,' mg/L',/7x,'differs from the input',' file/default value of ',a,' mg/L.',/7x,'The calculated value will be used in subsequent',' calculations.')
        end if
    end if

    ! If there is a density calculated from the WIPP brine density
    ! model, use it as the density in subsequent calculations. If not
    ! and a density value greater than zero was entered on the
    ! input file, use the input value.
    if (qdwipp) then
        qrho = .true.
        rho = rhoc
        tdspkg = tdspkc
        tdspl = tdsplc
    else if (qrho) then
        mrmlra = 0.001*wfh2o*rho
        mlmrra = 1./mrmlra
    end if

    if (.not.qrho) then
        rho = 0.0
    end if

    ! Compute mass balance totals for chemical elements.
    do nc = 1,nctmax
        mte(nc) = 0.
        mteaq(nc) = 0.
        cteaq(nc) = 0.
    end do

    do nss = narn1,narn2
        ns = jcsort(nss)

        if (jsflag(ns) .le. 0) then
            nr1 = nessr(1,ns)
            nr2 = nessr(2,ns)

            do n = nr1,nr2
                nc = ness(n)
                mteaq(nc) = mteaq(nc) + cess(n)*mosp(ns)
            end do
        end if
    end do

    do nc = 1,nct
        mte(nc) = mteaq(nc)
        cteaq(nc) = wkgwi*mteaq(nc)
        ppmwe(nc) = 1000.*cteaq(nc)*atwt(nc)*wfh2o
    end do

    do nss = nern1,nern2
        ns = jcsort(nss)

        if (jsflag(ns) .le. 0) then
            nr1 = nessr(1,ns)
            nr2 = nessr(2,ns)

            do n = nr1,nr2
                nc = ness(n)

                if (nc .gt. 0) then
                    mte(nc) = mte(nc) + cess(n)*mosp(ns)
                end if
            end do
        end if
    end do

    ! Compute mass balance totals for active basis species.
    do nb = 1,nbtmax
        mtb(nb) = 0.
        mtbaq(nb) = 0.
        ctb(nb) = 0.
    end do

    do nss = narn1,narn2
        ns = jcsort(nss)

        if (jsflag(ns) .le. 0) then
            nr1 = nstsr(1,ns)
            nr2 = nstsr(2,ns)

            do n = nr1,nr2
                nb = nsts(n)
                mtbaq(nb) = mtbaq(nb) + csts(n)*mosp(ns)
            end do
        end if
    end do

    do nb = 1,nbt
        mtb(nb) = mtbaq(nb)
        ctb(nb) = wkgwi*mtbaq(nb)
    end do

    do nss = nern1,nern2
        ns = jcsort(nss)

        if (jsflag(ns) .le. 0) then
            nr1 = nstsr(1,ns)
            nr2 = nstsr(2,ns)

            do n = nr1,nr2
                nb = nsts(n)
                mtb(nb) = mtb(nb) + csts(n)*mosp(ns)
            end do
        end if
    end do

    ! Compute apparent "whole-phase" equivalent fractions and mole
    ! fractions of the exchange ions present in generic ion exchanger
    ! phases. Cations and anions are treated separately in these
    ! calculations.
    call gegexw(cegexs,egexpc,egexpa,egexw,iern1,iern2,ietmax,jern1,jetmax,jgext,kern1,kern2,ketmax,kgexsa,moph,mosp,netmax,ngexsa,ngext,noutpt,nptmax,nstmax,nttyo,xgexw,zchar)

    ! Initialize the flag array which marks if compositions are known
    ! which maximize the affinities of the phases. Here assume that
    ! this is the case for all phases except solid solutions. The
    ! compositions for these phases will be determined by EQLIB/hpsat.f,
    ! which will be called below by scripx.f. Note that input solid
    ! solution compositions do not in general maximize the affinities.
    ! Thus, they don't count here.
    do np = 1, npt
        qxknph(np) = .true.
    end do

    do np = ixrn1,ixrn2
        qxknph(np) = .false.
    end do

    ! Print the computed description of the aqueous solution.
    call scripx(abar,acflg,act,actlg,adh,afcnst,affpd,affsd,ahrc,alki,apx,atwt,a3bar,a3bars,bpx,cdrsd,cegexs,cess,conc,conclg,coval,csts,ctb,cteaq,egexjc,egexjf,egexpa,egexpc,egexs,egexw,ehfac,ehrc,eps100,eh,farad,fje,fo2,fo2lg,fo2lrc,fugac,fugalg,fxi,iapxmx,ibpxmx,iebal,iern1,iern2,ietmax,igas,iktmax,iopg,iopr,iopt,ilrn1,ilrn2,imrn1,imrn2,ixrn1,ixrn2,jcsort,jern1,jern2,jetmax,jflag,jflagd,jflgi,jfleba,jgext,jgsort,jpflag,jsflag,jsol,jsomax,kern1,kern2,ketmax,kgexsa,mlmrra,mrmlra,moph,mosp,mte,mteaq,mwtsp,narn1,narn2,natmax,nbasp,nbaspd,nbt,nbtmax,nchlor,ncmpr,nct,nctmax,ndrsd,ndrsmx,ndrsrd,nelect,ness,nessmx,nessr,net,neti,netmax,ngexpi,ngexsa,ngext,ngrn1,ngrn2,ngt,ngtmax,nhydr,nhydx,nopgmx,noprmx,noptmx,noutpt,no2gaq,npnxp,npt,nptmax,nrdxsp,nst,nstmax,nsts,nstsmx,nstsr,ntf1,ntf1mx,ntf1t,ntf2,ntf2mx,ntf2t,nttyo,nxrn1,nxrn2,nxt,nxti,nxtimx,nxtmax,omega,pe,perc,ppmwe,qrho,qxknph,rho,rhoc,rhowc,sidrph,sidrsp,sigmam,sigzi,tdsglw,tdspkc,tdspkg,tdspl,tdsplc,tempc,tf1,tf2,tolspf,uelem,ugexj,ugexmo,uphase,uspec,uxtype,vosol,wfac,wfh2o,wftds,wkgwi,woh2o,wosol,wotds,xbar,xbarlg,xbarw,xbrwlg,xgexw,xlke,xlksd,zchar,zchcu6,zchsq2)

    if (nprob .le. 1) then
        ! Calculate data to describe the aqueous solution in the first
        ! problem on the input file as a special reactant. This defines
        ! the so-called 'Fluid 2' to be described on the EQ3NR pickup
        ! file under the iopt(19) = 3 option.
        ureac1 = 'Fluid 2'
        tempc1 = tempc

        do nc = 1,nct
            n = nc
            uesr1(n) = uelem(nc)
            cesr1(n) = mteaq(nc)
        end do

        iesrt1 = nct

        n = 1
        ubsr1(n) = ureac1
        cbsr1(n) = -1.0

        do nb = 1,nbt
            ns = nbaspd(nb)

            if (jflag(ns) .lt. 30) then
                n = n + 1
                ubsr1(n) = uspec(ns)
                cbsr1(n) = mtbaq(nb)
            end if
        end do

        ibsrt1 = n
    end if

    if (iopt(17) .ge. 0) then
        ! Write the pickup file.
        ! First set up some needed data.
        call setpk3(electr,iindx1,jflag,jflgi,kbt,kdim,kmax,kmt,kprs,kwater,kxt,mtb,mtbi,mtbaq,mtbaqi,narn1,narn2,nbasp,nbaspd,nbaspi,nbti,nbtmax,ndrsrd,nern1,nern2,nobswt,nstmax,ntitl,ntitl2,ntitmx,omeglg,press,pressi,scamas,sigzi,tempc,tempci,ubmtbi,uobsw,uspec,utitl,utitl2,uzveci,uzvec1,zvclgi,zvclg1)

        if (iopr(17) .le. 0) then
            upkfor = uinfor
        else if (iopr(17) .eq. 1) then
            upkfor = 'W'
        else if (iopr(17) .ge. 2) then
            upkfor = 'D'
        end if

        if (iopt(19) .le. 0) then
            ! Write a normal EQ3NR pickup file. This corresponds to the
            ! bottom half of a full EQ6 input file. To use this pickup
            ! file in an EQ6 input file, you must concatenate it with the
            ! top half. Usually this top half is extracted from an
            ! existing EQ6 input file and edited as necessary.
            if (upkfor(1:1) .eq. 'W') then
                ! Compact (W) format.
                call wr3pkw(electr,cgexj,ietmax,iopg,jetmax,jflgi,jgext,kbt,kct,kdim,kmax,kmt,kxmod,kxt,mtbi,mtbaqi,mwtges,nbti,nbtmax,net,netmax,newin,ngexrt,nobswt,nopgmx,nsbswt,ntitl2,ntitmx,nxmdmx,nxmod,qgexsh,pressi,tempci,tgexp,ubmtbi,ugexj,ugexmo,ugexp,ugexr,uhfgex,uobsw,usbsw,utitl2,uvfgex,uxkgex,uxmod,uzveci,xhfgex,xlkgex,xlkmod,xvfgex,zgexj,zvclgi)
            else
                ! Menu-style (D) format.
                call wr3pkd(electr,cgexj,ietmax,iopg,jetmax,jflgi,jgext,kbt,kct,kdim,kmax,kmt,kxmod,kxt,mtbi,mtbaqi,mwtges,nbti,nbtmax,net,netmax,newin,ngexrt,nobswt,nopgmx,nsbswt,ntitl2,ntitmx,nxmdmx,nxmod,qgexsh,pressi,tempci,tgexp,ubmtbi,ugexj,ugexmo,ugexp,ugexr,uhfgex,uobsw,usbsw,utitl2,uvfgex,uxkgex,uxmod,uzveci,xhfgex,xlkgex,xlkmod,xvfgex,zgexj,zvclgi)
            end if
        else
            ! Write an advanced EQ3NR pickup file. This corresponds to
            ! a full EQ6 input file. There are currently three options
            ! for the form and content of the top half of this file.
            ! These are determined by the iopt(19) option switch. See
            ! comments in EQ3NR/stpk36.f.
            itmxsv = itermx

            call stpk36(awmaxi,awmini,cbsri,cbsr1,cdac,cesri,cesr1,csigma,dlaplo,dlaprn,dleplo,dleprn,dlhplo,dlhprn,dloplo,dloprn,dltpll,dltplo,dltprl,dltprn,dlxdmp,dlxmx0,dlxpll,dlxplo,dlxprl,dlxprn,eact,ehmaxi,ehmini,fkrc,hact,iact,ibsrti,ibsrt1,iesrti,iesrt1,iktmax,imchmx,imech,iopt,itermx,ixrti,jcode,jpress,jreac,jtemp,ksplmx,ksppmx,kstpmx,modr,moffg,morr,mprphi,mprspi,nbt1mx,nctmax,ndact,ndctmx,nffg,nffgmx,noptmx,nordmx,noutpt,nprob,nprpmx,nprpti,nprsmx,nprsti,nptkmx,nrct,nrctmx,nrk,nsk,nsrt,nsrtmx,ntitl1,ntitmx,ntrymx,nttkmx,nttyo,nxopex,nxopmx,nxopt,nxpemx,nxrt,nxrtmx,o2maxi,o2mini,phmaxi,phmini,ptk,press,pressb,rkb,rxbari,sfcar,ssfcar,tempc,tempcb,tempc1,timmxi,tistti,tolsat,tolxsf,trkb,ttk,ubsri,ubsr1,ucxri,udac,uesri,uesr1,uffg,uprphi,uprspi,ureac,ureac1,utitl1,uxcat,uxopex,uxopt,vreac,ximaxi,xistti,xlkffg)

            itermx = itmxsv

            if (upkfor(1:1) .eq. 'W') then
                ! Compact (W) format.
                call wr6pkw(awmaxi,awmini,cbsri,cdac,cesri,cgexj,csigma,dlaplo,dlaprn,dleplo,dleprn,dlhplo,dlhprn,dloplo,dloprn,dltpll,dltplo,dltprl,dltprn,dlxdmp,dlxmx0,dlxpll,dlxplo,dlxprl,dlxprn,eact,egersi,ehmaxi,ehmini,electr,fkrc,iact,ibsrti,iesrti,ietmax,igerti,iktmax,imchmx,imech,iodb,iopg,iopr,iopt,itermx,ixrti,jcode,jetmax,jflgi,jgerti,jgext,jpress,jreac,jtemp,kbt,kct,kdim,kmax,kmt,kprs,ksplmx,ksppmx,kstpmx,kxmod,kxt,hact,modr,moffg,morr,mprphi,mprspi,mtbaqi,mtbi,mwtges,nbti,nbtmax,nbt1mx,nctmax,ndact,ndctmx,nert,nertmx,net,netmax,newin,nffg,nffgmx,ngexrt,nobswt,nodbmx,nopgmx,noprmx,noptmx,nordmx,nprpmx,nprpti,nprsmx,nprsti,nptkmx,nrct,nrctmx,nrk,nsbswt,nsk,nsrt,nsrtmx,ntitl1,ntitl2,ntitmx,ntrymx,nttkmx,nxmdmx,nxmod,nxopex,nxopmx,nxopt,nxpemx,nxrt,nxrtmx,o2maxi,o2mini,phmaxi,phmini,pressb,pressi,ptk,qgexsh,rkb,rxbari,sfcar,ssfcar,tempcb,tempci,tgexp,timmxi,tistti,tolbt,toldl,tolsat,tolxsf,trkb,ttk,ubmtbi,ubsri,ucxri,udac,uesri,uffg,ugerji,ugermo,ugersi,ugexj,ugexmo,ugexp,ugexr,uhfgex,uobsw,uprphi,uprspi,ureac,usbsw,utitl1,utitl2,uvfgex,uxcat,uxkgex,uxmod,uxopex,uxopt,uzveci,vreac,xgersi,xhfgex,ximaxi,xistti,xlkffg,xlkgex,xlkmod,xvfgex,zgexj,zvclgi)
            else
                ! Menu-style (D) format.
                call wr6pkd(awmaxi,awmini,cbsri,cdac,cesri,cgexj,csigma,dlaplo,dlaprn,dleplo,dleprn,dlhplo,dlhprn,dloplo,dloprn,dltpll,dltplo,dltprl,dltprn,dlxdmp,dlxmx0,dlxpll,dlxplo,dlxprl,dlxprn,eact,egersi,ehmaxi,ehmini,electr,fkrc,iact,ibsrti,iesrti,ietmax,igerti,iktmax,imchmx,imech,iodb,iopg,iopr,iopt,itermx,ixrti,jcode,jetmax,jflgi,jgerti,jgext,jpress,jreac,jtemp,kbt,kct,kdim,kmax,kmt,kprs,ksplmx,ksppmx,kstpmx,kxmod,kxt,hact,modr,moffg,morr,mprphi,mprspi,mtbaqi,mtbi,mwtges,nbti,nbtmax,nbt1mx,nctmax,ndact,ndctmx,nert,nertmx,net,netmax,newin,nffg,nffgmx,ngexrt,nobswt,nodbmx,nopgmx,noprmx,noptmx,nordmx,nprpmx,nprpti,nprsmx,nprsti,nptkmx,nrct,nrctmx,nrk,nsbswt,nsk,nsrt,nsrtmx,ntitl1,ntitl2,ntitmx,ntrymx,nttkmx,nxmdmx,nxmod,nxopex,nxopmx,nxopt,nxpemx,nxrt,nxrtmx,o2maxi,o2mini,phmaxi,phmini,pressb,pressi,ptk,qgexsh,rkb,rxbari,sfcar,ssfcar,tempcb,tempci,tgexp,timmxi,tistti,tolbt,toldl,tolsat,tolxsf,trkb,ttk,ubmtbi,ubsri,ucxri,udac,uesri,uffg,ugerji,ugermo,ugersi,ugexj,ugexmo,ugexp,ugexr,uhfgex,uobsw,uprphi,uprspi,ureac,usbsw,utitl1,utitl2,uvfgex,uxcat,uxkgex,uxmod,uxopex,uxopt,uzveci,vreac,xgersi,xhfgex,ximaxi,xistti,xlkffg,xlkgex,xlkmod,xvfgex,zgexj,zvclgi)
            end if
        end if

        write (noutpt,1800)
        write (nttyo,1800)
1800 format(/' The pickup file has been written.')
    else
        write (noutpt,1810)
        write (nttyo,1810)
1810 format(/' No pickup file was written.')
    end if

    ! Go look for another problem on the input file.
    go to 20
end program eq3nr
