! eqlpar.h
!     Explanations of the dimensioning parameters and equivalent
!     variables may be found in the file eqldef.h.
!     Some parameters are also defined and set in the include files
!     eqlj8.h and eqlo8.h. These are involved in the menu-style ("D")
!     input format for Version 8.
!     Dimensioning parameters and equivalent variables used in both
!     EQ3NR and EQ6 modes (fixed/variable mass of solvent water).
!       Primary parameters:
integer :: iet_par
integer :: jet_par
integer :: jso_par
integer :: net_par
integer :: nodb_par
integer :: nopg_par
integer :: nopr_par
integer :: nopt_par
integer :: ntit_par
integer :: nxmd_par

parameter(iet_par = 10,jet_par = 4,jso_par = 12,net_par = 12,nodb_par = 20,nopg_par = 20,nopr_par = 20,nopt_par = 20,ntit_par = 200,nxmd_par = 100)

integer :: iapxa_par
integer :: ibpxa_par

parameter(iapxa_par = 30,ibpxa_par = 10)

! Secondary parameters:
integer :: ket_par

parameter(ket_par = iet_par*net_par)

! Dimensioning parameters and equivalent variables used in only
! EQ3NR mode (fixed mass of solvent water).
!   Primary parameters:
integer :: njf_par
integer :: nxti_par

parameter(njf_par = 30,nxti_par = 50)

! Dimensioning parameters and equivalent variables used in only
! EQ6 mode (variable mass of solvent water).
!   Primary parameters:
integer :: imch_par
integer :: ndct_par
integer :: nert_par
integer :: nffg_par
integer :: nprp_par
integer :: nptk_par
integer :: nrct_par
integer :: nsrt_par
integer :: nttk_par
integer :: nxop_par
integer :: nxpe_par
integer :: nxrt_par
integer :: nssc_par

parameter(imch_par = 4,ndct_par = 4,nert_par = 40,nffg_par = 20,nprp_par = 100,nptk_par = 3,nrct_par = 40,nsrt_par = 40,nttk_par = 3,nxop_par = 40,nxpe_par = 100,nxrt_par = 40,nssc_par = 20)

! Secondary parameters:
integer :: nprs_par

parameter(nprs_par = 3*nprp_par)

