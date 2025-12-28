subroutine rd6wvz(ninpt,ucode,urelno,ustage,ueqlrn,ueqlst)
    !! This subroutine looks at the comments on the old EQ6 INPUT file
    !! in compact ("W") format, and tries to recover the code and
    !! version number data normally written on the bottom half of the
    !! PICKUP file to identify the code and version which wrote this
    !! part of that file.
    !! This subroutine is called by:
    !!   XCON6/xcon6.f
    implicit none

    ! Calling sequence variable declarations.
    integer :: ninpt

    character(len=8) :: ucode
    character(len=8) :: urelno
    character(len=8) :: ustage
    character(len=8) :: ueqlrn
    character(len=8) :: ueqlst

    ! Local variable declarations.
    integer :: i
    integer :: j

    character(len=80) :: uline

    ucode  = ' '
    urelno = ' '
    ustage = ' '
    ueqlrn = ' '
    ueqlst = ' '
    rewind(ninpt)

    do 100 i = 1,10000
        read (ninpt,'(a)',end=999) uline

        if (uline(1:1) .eq. '*') then
            j = index(uline,'pickup file written by')

            if (j .le. 0) then
                j = index(uline,'Pickup file written by')
            end if

            if (j .gt. 0) then
                go to 110
            end if
        end if

100 continue

        go to 999

110 continue
        j = index(uline,'eq3nr.')

        if (j .eq. 0) then
            j = index(uline,'EQ3NR.')
        end if

        if (j .gt. 0) then
            ucode = 'eq3nr'
            urelno= uline(j + 6:j + 9)
            ustage= uline(j + 10:j+ 13)
            go to 120
        end if

        j = index(uline,'EQ3NR, version ')

        if (j .gt. 0) then
            ucode = 'EQ3NR'
            urelno= uline(j + 15:j + 18)
            ustage= uline(j + 21:j+ 28)
            go to 120
        end if

        j = index(uline,'eq6.')

        if (j .eq. 0) then
            j = index(uline,'EQ6.')
        end if

        if (j .gt. 0) then
            ucode = 'eq6'
            urelno= uline(j + 4:j + 7)
            ustage= uline(j + 8:j + 11)
            go to 120
        end if

        j = index(uline,'EQ6, version ')

        if (j .gt. 0) then
            ucode = 'EQ6'
            urelno= uline(j + 13:j + 16)
            ustage= uline(j + 19:j + 26)
            go to 120
        end if

120 continue
        read (ninpt,'(a)',end=999) uline
        j = index(uline,'supported by eqlib.')

        if (j .eq. 0) then
            j = index(uline,'supported by EQLIB.')
        end if

        if (j .gt. 0) then
            ueqlrn = uline(j + 19:j + 22)
            ueqlst = uline(j + 23:j + 26)
            go to 999
        end if

        j = index(uline,'supported by EQLIB, version ')

        if (j .gt. 0) then
            ueqlrn = uline(j + 28:j + 31)
            ueqlst = uline(j + 34:j + 41)
            go to 999
        end if

        j = index(uline,'Supported by EQLIB, version ')

        if (j .gt. 0) then
            ueqlrn = uline(j + 28:j + 31)
            ueqlst = uline(j + 34:j + 41)
            go to 999
        end if

999 continue
    end subroutine rd6wvz