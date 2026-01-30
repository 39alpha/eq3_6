subroutine wrhdr(ikt_asv,ipch_asv,ipcv_asv,jpfc_asv,nap_asv,narx_asv,nat_asv,nbt_asv,nct_asv,ndata1,ndat1f,ngt_asv,nlat_asv,nlt_asv,nmt_asv,nmut_asv,npt_asv,nst_asv,ntid_asv,ntpr_asv,nxt_asv,uakey)
    !! This suboutine writes a header on the DATA1 and DATA1F files.
    !! The header includes the string "data1', followed by the keystring
    !! for the type of aqueous species activity coefficient model, and
    !! the dimensioning parameters required to read the rest of the data
    !! on the DATA1 file. The dimensioning parameters are used by EQ3NR
    !! and EQ6.
    !! This suboutine is called by:
    !!   EQPT/eqpt.f
    !! Principal input:
    !!   ndata1    = the unit number of the DATA1 file
    !!   ndat1f    = the unit number of the DATA1F file
    !!   uakey     = the keystring identifying the type of model for the
    !!                 aqueous species activity coefficients
    !!   ikt_asv  = the maximum number of end-member component species
    !!                in any solid solution on the data file
    !!   jpfc_asv = the number of terms (coefficients) in the
    !!                temperature function used to represent Pitzer
    !!                interaction parameters
    !!   nap_asv  = maximum number of distinct sets of Pitzer alpha
    !!                parameters
    !!   nat_asv  = the number of aqueous species on the data file
    !!   nbt_asv  = the number of basis species on the data file
    !!   nct_asv  = the number of chemical elements on the data file
    !!   ngt_asv  = the number of gas species on the data file
    !!   nlt_asv  = the number of pure liquid species on the data file
    !!   nmt_asv  = the number of gas species on the data file
    !!   npt_asv  = The number of phases of all types on the data file
    !!   nlat_asv = the number of members in the set of Pitzer lambda
    !!                coefficients
    !!   nmut_asv = the number of members in the set of Pitzer mu
    !!                coefficients
    !!   nst_asv  = the number of species of all types on the data file
    !!   ntid_asv = the number of lines in the data file title
    !!   nxt_asv  = the number of solid-solution phases on the data file
    !! Principal output:
    !!   None
    implicit none

    ! Calling sequence variable declarations.
    integer :: ndata1
    integer :: ndat1f

    integer :: ikt_asv
    integer :: ipch_asv
    integer :: ipcv_asv
    integer :: jpfc_asv
    integer :: nap_asv
    integer :: narx_asv
    integer :: nat_asv
    integer :: nbt_asv
    integer :: nct_asv
    integer :: ngt_asv
    integer :: nlat_asv
    integer :: nlt_asv
    integer :: nmt_asv
    integer :: nmut_asv
    integer :: npt_asv
    integer :: nst_asv
    integer :: ntid_asv
    integer :: ntpr_asv
    integer :: nxt_asv

    character(len=8) :: uakey

    ! Local variable declarations.
    integer :: j2
    integer :: j3

    integer :: ilnobl

    character(len=72) :: uterm
    character(len=72) :: utermc
    character(len=56) :: ux56
    character(len=8) :: udat1

    uterm(1:48) = '+-----------------------------------------------'
    uterm(49:72) = '------------------------'
    utermc = uterm
    utermc(1:1) = '*'

    ! Write 'data1' at the top of the DATA1 and DATA1F files.
    udat1 = 'data1'
    j2 = ilnobl(udat1)

    write (ndata1) udat1
    write (ndat1f,1000) udat1(1:j2)
1000 format(a)

    ! Next write the keystring for the type of model for the activity
    ! coefficients of the aqueous species.
    j3 = ilnobl(uakey)

    write (ndata1) uakey
    write (ndat1f,1000) uakey(1:j3)
    write (ndat1f,1000) utermc(1:72)

    ! Write the minimum dimensioning parameters required by the current
    ! data file. This will be used by EQ3NR and EQ6 to allocate array
    ! sizes for arrays which will contain the data read from this data
    ! file.
    ! Write the standard grid parameters on the output files.
    ux56 = 'Dim: chemical elements'
    write (ndata1) ux56
    j2 = ilnobl(ux56)
    write (ndat1f,1000) ux56(1:j2)
    write (ndata1) nct_asv
    write (ndat1f,1010) nct_asv
1010 format(i5)

    ux56 = 'Dim: basis species'
    write (ndata1) ux56
    j2 = ilnobl(ux56)
    write (ndat1f,1000) ux56(1:j2)
    write (ndata1) nbt_asv
    write (ndat1f,1010) nbt_asv

    ux56 = 'Dim: total species of all types'
    write (ndata1) ux56
    j2 = ilnobl(ux56)
    write (ndat1f,1000) ux56(1:j2)
    write (ndata1) nst_asv
    write (ndat1f,1010) nst_asv

    ux56 = 'Dim: total phases'
    write (ndata1) ux56
    j2 = ilnobl(ux56)
    write (ndat1f,1000) ux56(1:j2)
    write (ndata1) npt_asv
    write (ndat1f,1010) npt_asv

    ux56 = 'Dim: aqueous species'
    write (ndata1) ux56
    j2 = ilnobl(ux56)
    write (ndat1f,1000) ux56(1:j2)
    write (ndata1) nat_asv
    write (ndat1f,1010) nat_asv

    ux56 = 'Dim: pure minerals'
    write (ndata1) ux56
    j2 = ilnobl(ux56)
    write (ndat1f,1000) ux56(1:j2)
    write (ndata1) nmt_asv
    write (ndat1f,1010) nmt_asv

    ux56 = 'Dim: pure liquids'
    write (ndata1) ux56
    j2 = ilnobl(ux56)
    write (ndat1f,1000) ux56(1:j2)
    write (ndata1) nlt_asv
    write (ndat1f,1010) nlt_asv

    ux56 = 'Dim: gas species'
    write (ndata1) ux56
    j2 = ilnobl(ux56)
    write (ndat1f,1000) ux56(1:j2)
    write (ndata1) ngt_asv
    write (ndat1f,1010) ngt_asv

    ux56 = 'Dim: solid solutions'
    write (ndata1) ux56
    j2 = ilnobl(ux56)
    write (ndat1f,1000) ux56(1:j2)
    write (ndata1) nxt_asv
    write (ndat1f,1010) nxt_asv

    ux56 = 'Dim: max. components in a solid solution'
    write (ndata1) ux56
    j2 = ilnobl(ux56)
    write (ndat1f,1000) ux56(1:j2)
    write (ndata1) ikt_asv
    write (ndat1f,1010) ikt_asv

    ux56 = 'Dim: lines in the data file title'
    write (ndata1) ux56
    j2 = ilnobl(ux56)
    write (ndat1f,1000) ux56(1:j2)
    write (ndata1) ntid_asv
    write (ndat1f,1010) ntid_asv

    ux56 = 'Dim: ranges in the logK temperature grid'
    write (ndata1) ux56
    j2 = ilnobl(ux56)
    write (ndat1f,1000) ux56(1:j2)
    write (ndata1) ntpr_asv
    write (ndat1f,1010) ntpr_asv

    ux56 = 'Dim: max. points in a logK temperature grid range'
    write (ndata1) ux56
    j2 = ilnobl(ux56)
    write (ndat1f,1000) ux56(1:j2)
    write (ndata1) narx_asv
    write (ndat1f,1010) narx_asv

    ux56 = 'Dim: dH/dP order'
    write (ndata1) ux56
    j2 = ilnobl(ux56)
    write (ndat1f,1000) ux56(1:j2)
    write (ndata1) ipch_asv
    write (ndat1f,1010) ipch_asv

    ux56 = 'Dim: dV/dP order'
    write (ndata1) ux56
    j2 = ilnobl(ux56)
    write (ndat1f,1000) ux56(1:j2)
    write (ndata1) ipcv_asv
    write (ndat1f,1010) ipcv_asv

    if (uakey(1:8) .eq. 'Pitzer  ') then
        ux56 = 'Dim: distinct sets of Pitzer alpha parameters'
        write (ndata1) ux56
        j2 = ilnobl(ux56)
        write (ndat1f,1000) ux56(1:j2)
        write (ndata1) nap_asv
        write (ndat1f,1010) nap_asv

        ux56 = 'Dim: Pitzer lambda parameters'
        write (ndata1) ux56
        j2 = ilnobl(ux56)
        write (ndat1f,1000) ux56(1:j2)
        write (ndata1) nlat_asv
        write (ndat1f,1010) nlat_asv

        ux56 = 'Dim: Pitzer mu parameters'
        write (ndata1) ux56
        j2 = ilnobl(ux56)
        write (ndat1f,1000) ux56(1:j2)
        write (ndata1) nmut_asv
        write (ndat1f,1010) nmut_asv

        ux56 = 'Dim: coefficients in the Pitzer parameter temp. function'
        write (ndata1) ux56
        j2 = ilnobl(ux56)
        write (ndat1f,1000) ux56(1:j2)
        write (ndata1) jpfc_asv
        write (ndat1f,1010) jpfc_asv
    end if

    write (ndat1f,1000) utermc(1:72)
end subroutine wrhdr
