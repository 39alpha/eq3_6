! eqldv.h
!     Dimensioning variables.
!     Those used in both EQ3NR and EQ6 modes (fixed/variable mass of
!     solvent water).
integer :: iapxmx
integer :: ibpxmx
integer :: ietmax
integer :: iktmax
integer :: ipbtmx
integer :: ipchmx
integer :: ipcvmx
integer :: jetmax
integer :: jpfcmx
integer :: jsomax
integer :: ketmax
integer :: kmax
integer :: napmax
integer :: narxmx
integer :: natmax
integer :: nazmmx
integer :: nazpmx
integer :: nbtmax
integer :: nbt1mx
integer :: nctmax
integer :: netmax
integer :: ndrsmx
integer :: nessmx
integer :: ngtmax
integer :: nltmax
integer :: nmtmax
integer :: nmutmx
integer :: nmxmax
integer :: nodbmx
integer :: nopgmx
integer :: noprmx
integer :: noptmx
integer :: nptmax
integer :: nsltmx
integer :: nstmax
integer :: nstsmx
integer :: nsxmax
integer :: ntidmx
integer :: ntitmx
integer :: ntf1mx
integer :: ntf2mx
integer :: ntfxmx
integer :: ntprmx
integer :: nvetmx
integer :: nxmdmx
integer :: nxtmax

common /eqldv/ iapxmx,ibpxmx,ietmax,iktmax,ipbtmx,ipchmx,ipcvmx,jetmax,jpfcmx,jsomax,ketmax,kmax,napmax,narxmx,natmax,nazmmx,nazpmx,nbtmax,nbt1mx,nctmax,netmax,ndrsmx,nessmx,ngtmax,nltmax,nmtmax,nmutmx,nmxmax,nodbmx,nopgmx,noprmx,noptmx,nptmax,nsltmx,nstmax,nstsmx,nsxmax,ntidmx,ntitmx,ntf1mx,ntf2mx,ntfxmx,ntprmx,nvetmx,nxmdmx,nxtmax

! Those used in only EQ3NR mode (fixed mass of solvent water).
integer :: njfmax
integer :: nxicmx
integer :: nxtimx

common /eqldv1/ njfmax,nxicmx,nxtimx

! Those used in only EQ6 mode (variable mass of solvent water).
integer :: imchmx
integer :: ndctmx
integer :: nertmx
integer :: nffgmx
integer :: nordmx
integer :: npetmx
integer :: nprpmx
integer :: nprsmx
integer :: nptkmx
integer :: nrctmx
integer :: nrd1mx
integer :: nsetmx
integer :: nttkmx
integer :: nsrtmx
integer :: nxopmx
integer :: nxpemx
integer :: nxrtmx
integer :: nsscmx

common /eqldv2/ imchmx,ndctmx,nertmx,nffgmx,nordmx,npetmx,nprpmx,nprsmx,nptkmx,nrctmx,nrd1mx,nsetmx,nttkmx,nsrtmx,nxopmx,nxpemx,nxrtmx,nsscmx

