subroutine inbdot(azeroa,insgfa,nad1,narn1a,narn2a,nata,nata_asv,nerr,noutpt,nsta_asv,nttyo,uspeca)
    !! This subroutine reads from the data file the block of individual
    !! species parameters for the B-dot model of aqueous species
    !! activity coefficients. For each solute species, these parameters
    !! are its hard core diameter (azeroa) and its neutral species flag
    !! (insgfa). The latter is relevant only for electrically neutral
    !! species. If it is has a value of 0, the log activity coefficient
    !! is set to zero; if it has a value of -1, the log activity
    !! coefficient is represented by the Drummond (1981) polynomial.
    !! Water (the solvent) is not subject to control by the insgfa flag.
    !! This subroutine is called by:
    !!   EQLIB/indata.f
    !! Principal input:
    !!   noutpt = unit number of the output file
    !!   nttyo  = unit number of the screen file
    !!   uspeca = array of names of species read from the data file
    !! Principal output:
    !!   azeroa = array of hard core diameters read from the data file
    !!   insgfa = array of flags for treating the activity coefficients
    !!              of species that happen to be electrically neutral
    implicit none

    ! Calling sequence variable declarations.
    integer :: nata_asv
    integer :: nsta_asv

    integer :: insgfa(nata_asv)
    integer :: nad1
    integer :: narn1a
    integer :: narn2a
    integer :: nata
    integer :: nerr
    integer :: noutpt
    integer :: nttyo

    character(len=48) :: uspeca(nsta_asv)

    real(kind=8) :: azeroa(nata_asv)

    ! Local variable declarations.
    integer :: igx
    integer :: j2
    integer :: na
    integer :: naz
    integer :: ncount
    integer :: ns

    integer :: ilnobl

    real(kind=8) :: azdef
    real(kind=8) :: azx

    character(len=24) :: unam
    character(len=24) :: unam1
    character(len=24) :: unam2
    character(len=8) :: uendit

    data uendit /'endit.  '/

    ! The following variable is the default value for the hard core
    ! diameter of an aqueous solute species. Note that hard core
    ! diameters are given in units of 10**-8 cm.
    data azdef  /4.0/

    ! Initialize the azeroa and insgfa arrays to zero.
    do na = 1,nata_asv
        azeroa(na) = 0.
        insgfa(na) = 0
    end do

    ! Read the block contents.
    naz = 0

    ! Read the species name and azeroa and insgfa data.
110 continue
    read (nad1) unam,azx,igx

    ! Test for the end of the data block.
    if (unam(1:8) .eq. uendit(1:8)) then
        go to 120
    end if

    ! Find the species index, ns.
    ! Calling sequence substitutions:
    !   narn1a for nrn1a
    !   narn2a for nrn2a
    call srchn(narn1a,narn2a,ns,nsta_asv,unam,uspeca)

    ! If not found among the loaded species, skip.
    if (ns .le. 0) then
        j2 = ilnobl(unam)
        write (noutpt,1020) unam(1:j2)
        write (nttyo,1020) unam(1:j2)
1020 format(/' * Warning - (EQLIB/inbdot) Have B-dot model',' parameters listed',/7x,'for an aqueous species named ',a,', but there is',/7x,'no data block for such a species.')

        go to 110
    end if

    ! Compute the aqueous species index na.
    na = ns - narn1a + 1

    ! Test for duplicate input.
    if (azeroa(na) .gt. 0.) then
        write (noutpt,1030) unam(1:j2),azx,igx,azeroa(na),insgfa(na)
        write (nttyo,1030) unam(1:j2),azx,igx,azeroa(na),insgfa(na)
1030 format(/' * Warning - (EQLIB/inbdot) Have a duplicate entry',' on the data file for',/7x,'the B-dot model parameters of',' the species ',a,'.',/7x,'The duplicate entry values are:',//9x,'azero= ',f7.2,', insgf= ',i3,//7x,'The first entry values are:',//9x,'azero= ',f7.2,', insgf= ',i3,//7x,'The first entry values will be used.')

        go to 110
    end if

    naz = naz + 1

    if (naz .gt. nata) then
        write (noutpt,1040) nata
        write (nttyo,1040) nata
1040 format(/' * Error - (EQLIB/inbdot) There are more entries',/7x,'for species with  B-dot model parameters than there are',/7x,'aqueous species, which number ',i4,'.')

        nerr = nerr + 1
        go to 999
    end if

    azeroa(na) = azx
    insgfa(na) = igx
    go to 110

    ! Check for solute species with no entries.
120 continue
    do na = 1,nata
        if (azeroa(na) .le. 0.) then
            go to 140
        end if
    end do

    go to 999

    ! Assign the default value if needed.
140 continue
    write (noutpt,1050) azdef
1050 format(/' * Note - (EQLIB/inbdot) The following aqueous species',' have been assigned',/7x,'a default hard core diameter of',' ',f6.3,' x 10**-8 cm:',/)

    ncount = 0

    do ns = narn1a,narn2a
        na = ns - narn1a + 1

        if (azeroa(na) .le. 0.) then
            azeroa(na) = azdef

            if (ns .ne. narn1a) then
                if (ncount .le. 0) then
                    unam1 = uspeca(ns)(1:24)
                    ncount = 1
                else
                    unam2 = uspeca(ns)(1:24)
                    j2 = ilnobl(unam2)
                    write (noutpt,1060) unam1,unam2(1:j2)
1060 format(9x,a24,3x,a)

                    ncount = 0
                end if
            end if
        end if
    end do

    if (ncount .eq. 1) then
        j2 = ilnobl(unam1)
        write (noutpt,1070) unam1(1:j2)
1070 format(9x,a)
    end if

999 continue
end subroutine inbdot