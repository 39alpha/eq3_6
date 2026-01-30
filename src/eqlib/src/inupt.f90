subroutine inupt(amua,aslma,ielam,ipbt_asv,jpdblo,jpfc_asv,nad1,nalpaa,napa_asv,napta,narn1a,narn2a,nerr,nmuta,nmuta_asv,nmuxa,noutpt,nslta,nslta_asv,nslxa,nsta_asv,nttyo,palpaa,uspeca,zchara)
    !! This subroutine reads the parameters for Pitzer's equations from
    !! the data file.
    !! This subroutine is called by:
    !!   EQLIB/indata.f
    !! Principal input:
    !!   nad1      = unit number of the data1 file
    !!   ipbt_asv  = dimensioning variable, the number of Pitzer
    !!                 alpha coefficients for any cation-anion pairnts
    !!   jpdblo    = the Pitzer data block organization flag:
    !!                 -1 = classical (pre-version 8)
    !!                  0 = new
    !!   jpfc_asv  = dimensioning variable, the number of coefficients
    !!                 in a Pitzer parameter temperature functionientsies
    !!   nmuta_asv = dimensioning variable, the maximum number of species
    !!                 triplets for which mu coefficients are defined
    !!   nslta_asv = dimensioning variable, the maximum number of species
    !!                 pairs for which S-lambda coefficients are defined
    !!   zchara    = array of electrical charge numbers of the species
    !!                 read from the data file
    !! Principal output:
    !!   amua   = coefficients for calculating Pitzer mu parameters
    !!              as a function of temperature
    !!   aslma  = coefficients for calculating Pitzer S-lamnda(n)
    !!              parameters as a function of temperature
    !!   nalpaa = pointer array giving the index of the set of alpha
    !!              coeffcients for a given set of S-lambda coefficients
    !!   nmuta  = number of species triplets for which mu
    !!              coefficients are defined
    !!   nmuxa  = indices of species in triplets for which mu
    !!              coefficients are defined
    !!   nslta  = number of species pairs for which S-lambda
    !!              coefficients are defined
    !!   nslxa  = indices of species in pairs for which S-lambda
    !!              coefficients are defined
    !!   palpaa = array of sets of alpha coefficients
    implicit none

    ! Calling sequence variable declarations.
    integer :: ipbt_asv
    integer :: jpfc_asv
    integer :: napa_asv
    integer :: nmuta_asv
    integer :: nslta_asv
    integer :: nsta_asv

    integer :: noutpt
    integer :: nttyo

    integer :: nalpaa(nslta_asv)
    integer :: nmuxa(3,nmuta_asv)
    integer :: nslxa(2,nslta_asv)

    integer :: ielam
    integer :: jpdblo
    integer :: nad1
    integer :: napta
    integer :: narn1a
    integer :: narn2a
    integer :: nerr
    integer :: nmuta
    integer :: nslta

    character(len=48) :: uspeca(nsta_asv)

    real(kind=8) :: amua(jpfc_asv,nmuta_asv)
    real(kind=8) :: aslma(jpfc_asv,0:ipbt_asv,nslta_asv)
    real(kind=8) :: palpaa(ipbt_asv,napa_asv)
    real(kind=8) :: zchara(nsta_asv)

    ! Local variable declarations.
    integer :: i
    integer :: j
    integer :: ja
    integer :: j2
    integer :: nblock
    integer :: nmu
    integer :: ns
    integer :: nsl

    integer :: ilnobl

    logical :: qx

    character(len=80) :: uline
    character(len=24) :: unam1
    character(len=24) :: unam2
    character(len=24) :: unam3
    character(len=16) :: ux16
    character(len=8) :: uelam
    character(len=8) :: uendit
    character(len=8) :: ux8

    real(kind=8) :: aax
    real(kind=8) :: z1
    real(kind=8) :: z2

    ! Local variable declarations with global dimensioning.
    ! These need not be SAVEd.
    real(kind=8), dimension(:), allocatable :: alphai

    data uendit /'endit.  '/

    ALLOCATE(alphai(ipbt_asv))

    if (jpdblo .eq. -1) then
        if (ipbt_asv .ne. 2) then
            write (ux8,'(i5)') ipbt_asv
            call lejust(ux8)
            j2 = ilnobl(ux8)
            write (noutpt,1000) ux8(1:j2)
            write (nttyo,1000) ux8(1:j2)
1000 format(/' * Error - (EQLIB/inupt) Have an illegal value',' of ',a,' for the array',/7x,'allocation size variable',' ipbt_asv (number of Pitzer alpha parameters',/7x,'for ',' any cation-anion or like pair. For the classical Pitzer',/7x,'data block organization present here, this variable',' must have',/7x,'a value of 2.')

            stop
        end if
    end if

    if (jpdblo .eq. -1) then
        ! In the "classical" Pitzer data block organization, the E-lambda
        ! (E-theta) flag may be "on" or "off", as specified on the data
        ! file. Read and decode this flag.
        uelam = ' '
        read (nad1) uline
        uelam(1:8) = uline(17:24)

        if (uelam(1:3) .eq. 'off') then
            ielam = -1
        else if (uelam(1:3) .eq. 'on ') then
            ielam = 0
        else
            j2 = ilnobl(uelam)
            write (noutpt,1010) uelam(1:j2)
            write (nttyo,1010) uelam(1:j2)
1010 format(/' * Error - (EQLIB/inupt) Have read an unrecognized',/7x,'value of "',a,'" for the E-lambda (E-theta) flag from',/7x,'the data file. Allowed values are "on" and "off".')

            stop
        end if
    else
        ! In the "new" Pitzer data block organization, the E-lambda
        ! flag is always "on".
        ielam = 0
    end if

    ! Read the species pairs for which S-lambda coefficients are
    ! defined. These data are comprised in two superblocks.
    ja = 0
    nsl = 0

    do nblock = 1,2
        ! Read a block.
        ! Read the names of the species in a pair.
15 continue
        read (nad1) unam1,unam2

        ! Test for end of superblock.
        if (unam1(1:6) .eq. uendit(1:6)) then
            go to 130
        end if

        ! Increment the S-lambda species pair counter.
        nsl = nsl + 1

        ! Calling sequence substitutions:
        !   narn1a for nrn1a
        !   narn2a for nrn2a
        !   unam1 for unam
        call srchn(narn1a,narn2a,ns,nsta_asv,unam1,uspeca)
        nslxa(1,nsl) = ns
        z1 = zchara(ns)

        ! Calling sequence substitutions:
        !   narn1a for nrn1a
        !   narn2a for nrn2a
        !   unam2 for unam
        call srchn(narn1a,narn2a,ns,nsta_asv,unam2,uspeca)
        nslxa(2,nsl) = ns
        z2 = zchara(ns)

        ! Read the temperature function coefficients for calculating
        ! Pitzer S-lambda(n) parameters. Read also the associated
        ! Pitzer alpha parameters.
        if (jpdblo .eq. -1) then
            ! Classical Pitzer data block organization.
            ! Read the lowest-order coefficients, which are here
            ! the 25C values of the corresponding interaction
            ! coefficients.
            read (nad1) (aslma(1,i,nsl), i = 0,2)

            ! Read the associated alpha coefficients.
            read (nad1) (alphai(i), i = 1,2)

            ! Read the higher-order coefficients, which here are
            ! first and second-order temperature derivatives.
            do i = 0,2
                read (nad1) aslma(2,i,nsl),aslma(3,i,nsl)
            end do
        else
            ! New Pitzer data block organization.
            ! Read the associated alpha coefficients.
            if ((z1*z2) .lt. 0.) then
                ! Read a larger set of data for cation-anion pairs.
                do i = 1,ipbt_asv
                    read (nad1) ux16,alphai(i)
                end do

                do i = 0,ipbt_asv
                    read (nad1) ux16

                    do j = 1,jpfc_asv
                        read (nad1) aslma(j,i,nsl)
                    end do
                end do
            else
                ! Read a more limited set of data for all other pair types.
                read (nad1) ux16

                do j = 1,jpfc_asv
                    read (nad1) aslma(j,0,nsl)
                end do
            end if
        end if

        ! Read the block terminator.
        read (nad1) uline(1:72)

        ! Test for alpha coefficient set already in the palpaa array.
        do j = 1,ja
            qx = .true.

            do i = 1,ipbt_asv
                aax = abs(palpaa(i,j) - alphai(i))

                if (aax .gt. 1.e-12) then
                    qx = .false.
                    go to 110
                end if
            end do

110 continue

            if (qx) then
                ! Found the index to an existing alpha pair. Store this
                ! index, then go process the data for another species pair.
                nalpaa(nsl) = j
                go to 15
            end if
        end do

        ! Not in the set; add it and store the index for this species
        ! pair.
        ja = ja + 1

        do i = 1,ipbt_asv
            palpaa(i,ja) = alphai(i)
        end do

        nalpaa(nsl) = ja

        ! Go process the data for another species pair.
        go to 15

130 continue
    end do

    nslta = nsl
    napta = ja

    DEALLOCATE(alphai)

    ! Read the species triplets for which mu coefficients are
    ! defined. These data are comprised in two superblocks.
    nmu = 0

    do nblock = 1,2
        ! Read a block.
        ! Read the names of the species in a triplet.
150 continue
        read (nad1) unam1,unam2,unam3

        ! Test for end of superblock.
        if (unam1(1:6) .eq. uendit(1:6)) then
            go to 200
        end if

        ! Increment the mu species triplet counter.
        nmu = nmu + 1

        ! Calling sequence substitutions:
        !   narn1a for nrn1a
        !   narn2a for nrn2a
        !   unam1 for unam
        call srchn(narn1a,narn2a,ns,nsta_asv,unam1,uspeca)
        nmuxa(1,nmu) = ns

        ! Calling sequence substitutions:
        !   narn1a for nrn1a
        !   narn2a for nrn2a
        !   unam2 for unam
        call srchn(narn1a,narn2a,ns,nsta_asv,unam2,uspeca)
        nmuxa(2,nmu) = ns

        ! Calling sequence substitutions:
        !   narn1a for nrn1a
        !   narn2a for nrn2a
        !   unam3 for unam
        call srchn(narn1a,narn2a,ns,nsta_asv,unam3,uspeca)
        nmuxa(3,nmu) = ns

        ! Read the coefficients for calculating Pitzer mu
        ! parameters as a function of temperature.
        if (jpdblo .eq. -1) then
            ! Classical Pitzer data block organization.
            read (nad1) amua(1,nmu),amua(2,nmu),amua(3,nmu)
        else
            ! New Pitzer data block organization.
            read (nad1) ux16

            do j = 1,jpfc_asv
                read (nad1) amua(j,nmu)
            end do
        end if

        ! Read the block terminator.
        read (nad1) uline(1:72)

        go to 150

200 continue
    end do

    nmuta = nmu
end subroutine inupt
