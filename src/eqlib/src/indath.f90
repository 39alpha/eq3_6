subroutine indath(ikta_asv,ipbt_asv,ipch_asv,ipcv_asv,jpfc_asv,nad1,napa_asv,narx_asv,nata_asv,nbta_asv,ncta_asv,ngta_asv,nlta_asv,nmta_asv,npta_asv,nmuta_asv,noutpt,nslta_asv,nsta_asv,ntid_asv,ntpr_asv,nttyo,nxta_asv,udakey)
    !! This subroutine reads the header section of the data1 file. This
    !! section consists of a record containing the string 'data1' (to
    !! ensure that the file is indeed a data1 file), a record containing
    !! the keystring for the activity coefficient model for aqueous
    !! species to which this data file corresponds, and the array
    !! allocation size variables required to allocate sufficient array
    !! space to store the rest of the data on this data file.
    !! This subroutine is called by:
    !!   EQ3NR/eq3nr.f
    !!   EQ6/eq6.f
    !! Principal input:
    !!   nad1   = unit number of the data1 file
    !!   noutpt = unit number of the output file
    !!   nttyo  = unit number of the screen file
    !! Principal output:
    !!   nbta_asv  = the number of basis species on the data file
    !!   ncta_asv  = the number of chemical elements on the data file
    !!   ntid_asv  = the number of lines in the data file title
    !!   ikta_asv  = the maximum number of end-member component species
    !!                 in any solid solution on the data file
    !!   nata_asv  = the number of aqueous species on the data file
    !!   ngta_asv  = the number of gas species on the data file
    !!   nlta_asv  = the number of pure liquids on the data file
    !!   nmta_asv  = the number of pure minerals on the data file
    !!   npta_asv  = The number of phases of all types on the data file
    !!   nsta_asv  = the number of species of all types on the data file
    !!   nxta_asv  = the number of solid-solution phases on the data
    !!                 file
    !!   napa_asv  = the number of distinct sets of Pitzer alpha values
    !!   nmuta_asv = the number of triplets of ions on the data file
    !!                 for which distinct Pitzer mu coefficients
    !!                 are defined
    !!   nslta_asv = the number of pairs of ions on the data file
    !!                 for which distinct Pitzer S-lambda coefficients
    !!                 are defined
    !!   ipch_asv = the maximum order for pressure corrections to
    !!                enthalpy functions
    !!   ipcv_asv = the maximum order for pressure corrections to
    !!                volume functions; the maximum order for pressure
    !!                corrections to log K and other Gibbs-energy-based
    !!                functions is one greater than this
    !!   ipbt_asv  = the maximum number of Pitzer alpha parameters
    !!                 for any species pair
    !!   jpfc_asv  = the number of coefficients in the Pitzer parameter
    !!                temperature function
    !!   narx_asv = maximum number of coefficients per temperature range
    !!   ntpr_asv = number of temperature ranges
    implicit none

    ! Calling sequence variable declarations.
    integer :: nad1
    integer :: noutpt
    integer :: nttyo

    integer :: ikta_asv
    integer :: ipbt_asv
    integer :: napa_asv
    integer :: nata_asv
    integer :: nbta_asv
    integer :: ncta_asv
    integer :: ngta_asv
    integer :: nlta_asv
    integer :: nmta_asv
    integer :: npta_asv
    integer :: nmuta_asv
    integer :: nslta_asv
    integer :: nsta_asv
    integer :: ntid_asv
    integer :: nxta_asv

    integer :: ipch_asv
    integer :: ipcv_asv
    integer :: jpfc_asv
    integer :: narx_asv
    integer :: ntpr_asv

    character(len=8) :: udakey

    ! Local variable declarations.
    character(len=56) :: ux56
    character(len=8) :: ux8

    write (noutpt,1000)
    write (nttyo,1000)
1000 format(/' Reading the data1 file header section ...')

    ! Rewind the data file. The unit number is nad1.
    rewind nad1

    ! Read the file header string.
    read (nad1) ux8

    if (ux8(1:5) .ne. 'data1') then
        write (noutpt,1010) ux8
        write (nttyo,1010) ux8
1010 format(/' * Error - (EQLIB/indath) Have wrong file header',' "',a5,'"',/7x,'on the data1 file. The first five characters',' must be "data1".')

        stop
    end if

    ! Read key string indicating the type of option supported by the
    ! data file for treating activity coefficients of aqueous species.
    ! This will be used to test the model specified on the input file
    ! for consistency with the data file. That test is made in
    ! EQLIB/cdakey.f.
    read (nad1) udakey

    ! Read the array dimension required for chemical elements.
    read (nad1) ux56
    read (nad1) ncta_asv

    ! Read the array dimension required for basis species.
    read (nad1) ux56
    read (nad1) nbta_asv

    ! Read the array dimension required for total species.
    read (nad1) ux56
    read (nad1) nsta_asv

    ! Read the array dimension required for total phases.
    read (nad1) ux56
    read (nad1) npta_asv

    ! Read the array dimension required for aqueous species.
    read (nad1) ux56
    read (nad1) nata_asv

    ! Read the array dimension required for pure minerals.
    read (nad1) ux56
    read (nad1) nmta_asv

    ! Read the array dimension required for pure liquids.
    read (nad1) ux56
    read (nad1) nlta_asv

    ! Read the array dimension required for gas species.
    read (nad1) ux56
    read (nad1) ngta_asv

    ! Read the array dimension required for solid-solution phases.
    read (nad1) ux56
    read (nad1) nxta_asv

    ! Read the array dimension required for solid-solution components
    ! (maximum number per solid solution).
    read (nad1) ux56
    read (nad1) ikta_asv

    ! Read the array dimension required for the data file title.
    read (nad1) ux56
    read (nad1) ntid_asv

    ! Read the array dimension required for the temperature ranges in
    ! the logK temperature grid.
    read (nad1) ux56
    read (nad1) ntpr_asv

    ! Read the array dimension required for the points in a temperature
    ! range on the logK temperature grid.
    read (nad1) ux56
    read (nad1) narx_asv

    ! Read the array dimension required for the dH/dP order represented
    ! on the data file.
    read (nad1) ux56
    read (nad1) ipch_asv

    ! Read the array dimension required for the dV/dP order represented
    ! on the data file.
    read (nad1) ux56
    read (nad1) ipcv_asv

    ! Read the array dimensions for parameters needed for activity
    ! coefficient models.
    if (udakey(1:8) .eq. 'Pitzer  ') then
        ! Read the array dimension required for the sets of Pitzer
        ! alpha coefficients.
        read (nad1) ux56
        read (nad1) napa_asv

        ! Read the array dimension required for the set of Pitzer
        ! S-lambda coefficients.
        read (nad1) ux56
        read (nad1) nslta_asv

        ! Read the array dimension required for the set of Pitzer mu
        ! coefficients.
        read (nad1) ux56
        read (nad1) nmuta_asv

        ! Read the number of coefficients for the Pitzer parameter
        ! temperature function.
        read (nad1) ux56
        read (nad1) jpfc_asv
    else
        ! Minimum values.
        napa_asv = 1
        nslta_asv = 1
        nmuta_asv = 1
        jpfc_asv = 1
    end if

    ! Set the number of parameters in a Pitzer alpha set.
    if (udakey(1:8) .eq. 'Pitzer  ') then
        ipbt_asv = 2
    else
        ipbt_asv = 1
    end if

    ! The following is a bit of nonsense so compiler warninings will
    ! not be generated that ux56 is not used.
    ux8 = ux56(1:8)
end subroutine indath