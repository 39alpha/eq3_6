subroutine rddixc(nxcon,uoldvd,unewf,unewv)
    !! This subroutine reads user input from the IXCON control file used
    !! by XCON3 and XCON6. The IXCON file should look as follows:
    !! 1 (Column 1) (This line is not part of the IXCON file.)
    !! IXCON
    !!   This file contains the user-controlled options for
    !!   the EQ3/6 input file converters XCON3 and XCON6.
    !!   The code will scan the old input file to determine
    !!   its format. It will also try to determine its
    !!   version level. If it can't do so, the version level
    !!   must be set using the default in this control file.
    !! OLD INPUT FILE VERSION LEVEL (DEFAULT ONLY)
    !!  | | 6.0 (versions 6.0-6.1)
    !!  |X| 7.0 (versions 7.0-7.1)
    !!  | | 7.2 (version 7.2)
    !!  | | 8.0 (version 8.0)
    !! NEW INPUT FILE FORMAT
    !!  | | W   (compact)
    !!  |X| D   (menu-style)
    !! NEW INPUT FILE VERSION LEVEL
    !!  | | 6.0 (versions 6.0-6.1)
    !!  | | 7.0 (versions 7.0-7.1)
    !!  |X| 7.2 (version 7.2)
    !!  | | 8.0 (version 8.0)
    !! END
    !! 1 (Column 1) (This line is not part of the IXCON file.)
    !! Note that the '|' characters defining the checkboxes must be
    !! in columns 2 and 4. Select a choice with a non-blank character.
    !! The first choice will be selected. Secondary choices will be
    !! ignored.
    !! This subroutine is called by:
    !!   XCON3/xcon3.f
    !!   XCON6/xcon6.f
    !! Input:
    !!   nxcon  = the unit number of the IXCON file
    !! Output:
    !!   uoldvd = the default version level of the old input file
    !!   unewf  = the desired format ("W" or "D") of the new input file
    !!   unewvd = the desired version level of the new input file
    implicit none

    ! Calling sequence variable declarations.
    integer :: nxcon

    character(len=3) :: uoldvd
    character(len=3) :: unewv
    character(len=1) :: unewf

    ! Local variable declarations.
    integer :: i
    integer :: ii
    integer :: j

    character(len=80) :: uline

    rewind(nxcon)

    do i = 1,1000
        read (nxcon,1000,end=999) uline
1000 format(a80)

        j = index(uline,'OLD INPUT FILE VERSION LEVEL')

        if (j .gt. 0) then
            do ii = 1,20
                read (nxcon,1000,end=999) uline

                if (uline(2:2).ne.'|' .or. uline(4:4).ne.'|') then
                    go to 120
                end if

                if (uline(3:3) .ne. ' ') then
                    uoldvd = uline(6:8)
                    go to 120
                end if
            end do

            go to 120
        end if
    end do

120 continue

    rewind(nxcon)

    do i = 1,1000
        read (nxcon,1000,end=999) uline
        j = index(uline,'NEW INPUT FILE FORMAT')

        if (j .gt. 0) then
            do ii = 1,20
                read (nxcon,1000,end=999) uline

                if (uline(2:2).ne.'|' .or. uline(4:4).ne.'|') then
                    go to 150
                end if

                if (uline(3:3) .ne. ' ') then
                    unewf = uline(6:6)
                    go to 150
                end if
            end do

            go to 150
        end if
    end do

150 continue

    rewind(nxcon)

    do i = 1,1000
        read (nxcon,1000,end=999) uline
        j = index(uline,'NEW INPUT FILE VERSION LEVEL')

        if (j .gt. 0) then
            do ii = 1,20
                read (nxcon,1000,end=999) uline

                if (uline(2:2).ne.'|' .or. uline(4:4).ne.'|') then
                    go to 180
                end if

                if (uline(3:3) .ne. ' ') then
                    unewv = uline(6:8)
                    go to 180
                end if
            end do

            go to 180
        end if
    end do

180 continue

999 continue
end subroutine rddixc