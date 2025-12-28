subroutine tpraz(nat,natmax,nazt,naztmx,ncvaz,nerr,noutpt,nttyo,pcvaz,qpdaz,uaqsp,uazp)
    !! Test and process the hard core diameter (azero) and neutral
    !! activity coefficient flag (insgf) data (i.e., the 'bdot' data)
    !! read from the DATA0 file. Find and flag errors, such as duplication
    !! of data (e.g., two entries for the same aqueous species).
    !! Check the coverage of entered data against all aqueous solute
    !! species present on the data file.
    !! This subroutine is called by:
    !!   EQPT/eqpt.f
    !! Principal input:
    !!   nat    = the number of aqueous species
    !!   nazt   = the number of specified hard core diameters
    !!   uaqsp  = array of names of aqueous species
    !!   uazp   = array of aqueous species names used to specify
    !!              hard core diamters on the data file
    !! Principal output:
    !!   nerr   = error counter
    implicit none

    ! Calling sequence variable declarations.
    integer :: natmax
    integer :: naztmx

    integer :: noutpt
    integer :: nttyo

    integer :: nat
    integer :: nazt
    integer :: ncvaz
    integer :: nerr

    logical :: qpdaz(natmax)

    character(len=24) :: uaqsp(natmax)
    character(len=24) :: uazp(naztmx)

    real(kind=8) :: pcvaz

    ! Local variable declarations.
    integer :: jaz
    integer :: j2
    integer :: j3
    integer :: j4
    integer :: j5
    integer :: na
    integer :: naz
    integer :: ncount
    integer :: ndupl
    integer :: nlistl
    integer :: nn
    integer :: nodatc

    integer :: ilnobl

    character(len=24) :: unam
    character(len=8) :: ux8

    ! Limit on the list of aqueous species for which no data
    ! were found.
    data nlistl / 20 /

    ! Initialize the qpdaz array. A value of .false. for a species
    ! means that no azero or insgf data for this species were read
    ! from the data file. Solvent water (the first aqueous species)
    ! does not have such data, so set the flag to .true. for this
    ! species.
    qpdaz(1) = .true.

    do na = 2,nat
        qpdaz(na) = .false.
    end do

    ! Check the entered azero and insgf data for aqueous species.
    nodatc = 0

    do na = 2,nat
        ! Search for uaqsp(na) in the uazp array. That array contains
        ! the aqueous species names for which azero and insgf data
        ! were read from the data file.
        unam = uaqsp(na)

        do naz = 1,nazt
            if (unam(1:24) .eq. uazp(naz)(1:24)) then
                go to 100
            end if
        end do

        naz = 0
100 continue

        if (naz .gt. 0) then
            ! Have found an entry.
            qpdaz(na) = .true.

            ! Check for duplicate data sets.
            ndupl = 0

            do jaz = naz + 1,nazt
                if (unam(1:24) .eq. uazp(jaz)(1:24)) then
                    ndupl = ndupl + 1
                end if
            end do

            if (ndupl .gt. 0) then
                j2 = ilnobl(unam)

                if (ndupl .eq. 1) then
                    write (noutpt,1010) unam(1:j2)
                    write (nttyo,1010) unam(1:j2)
1010 format(/' * Error - (EQPT/tpraz) Have found a duplicate',' entry on the DATA0 file',/7x,'for the hard core',' diameter of ',a,'.')
                else
                    ux8 = ' '
                    write (ux8,'(i5)') ndupl
                    call lejust(ux8)
                    j5 = ilnobl(ux8)
                    write (noutpt,1020) ux8(1:j5),unam(1:j2)
                    write (nttyo,1020) ux8(1:j5),unam(1:j2)
1020 format(/' * Error - (EQPT/tpraz) Have found ',a,' duplicate entries on the DATA0 file',/7x,'for the hard core diameter of ',a,'.')
                end if

                nerr = nerr + ndupl
            end if
        else
            ! No data entry was found on the DATA0 file.
            ! Note that qpdaz(na) is left with a value of .false.
            nodatc = nodatc + 1
        end if
    end do

    if (nodatc .gt. 0) then
        write (noutpt,1040)
        write (nttyo,1040)
1040 format(/' * Warning - (EQPT/tpraz) Did not find a hard core',' diameter entry on the',/7x,'DATA0 file for any of the',' following aqueous species:',/)

        ncount = 0

        do na = 1,nat
            if (.not.qpdaz(na)) then
                ncount = ncount + 1
                unam = uaqsp(na)
                j4 = ilnobl(unam)
                write (noutpt,1050) unam(1:j4)
                write (nttyo,1050) unam(1:j4)
1050 format(9x,a)

                if (ncount .eq. nlistl) then
                    go to 200
                end if
            end if
        end do

200 continue

        nn = nodatc - ncount

        if (nn .gt. 0) then
            write (ux8,'(i5)') nn
            call lejust(ux8)
            j3 = ilnobl(ux8)
            write (noutpt,1060) ux8(1:j3)
            write (nttyo,1060) ux8(1:j3)
1060 format(/9x,'plus ',a,' others')
        end if

        write (noutpt,1070)
        write (nttyo,1070)
1070 format(1x)
    end if

    ncvaz = nat - nodatc - 1
    pcvaz = (100.*ncvaz)/float(nat - 1)
end subroutine tpraz