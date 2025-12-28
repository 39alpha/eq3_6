subroutine adgexj(je,noutpt,nttyo,ugexjd)
    !! This subroutine assigns a default name for the je-th site of a
    !! generic exchange phase. The general model is "S(n)", where n
    !! is the site number. However, if n is larger than 99999, the S
    !! and both parentheses are all dropped.
    !! This subroutine is called by:
    !!   EQLIB/intexi.f
    !!   EQ3NR/intge3.f
    !! Principal input:
    !!   je     = site index
    !! Principal output:
    !!   ugexjd = default name of the je-th site of a generic ion
    !!              exchange phase
    implicit none

    ! Calling sequence variable declarations.
    integer :: noutpt
    integer :: nttyo

    integer :: je

    character(len=8) :: ugexjd

    ! Local variable declarations.
    ! None
    ugexjd = ' '
    ugexjd(1:2) = 'S('

    if (je .le. 9) then
        write (ugexjd(3:3),'(i1)') je
        ugexjd(4:4) = ')'
    else if (je .le. 99) then
        write (ugexjd(3:4),'(i2)') je
        ugexjd(5:5) = ')'
    else if (je .le. 999) then
        write (ugexjd(3:5),'(i3)') je
        ugexjd(6:6) = ')'
    else if (je .le. 9999) then
        write (ugexjd(3:6),'(i4)') je
        ugexjd(7:7) = ')'
    else if (je .le. 99999) then
        write (ugexjd(3:7),'(i5)') je
        ugexjd(8:8) = ')'
    else if (je .le. 999999) then
        ugexjd = 'S '
        write (ugexjd(2:7),'(i6)') je
    else if (je .le. 9999999) then
        ugexjd = 'S '
        write (ugexjd(2:8),'(i7)') je
    else if (je .le. 99999999) then
        ugexjd = ' '
        write (ugexjd(1:8),'(i8)') je
    else
        write (noutpt,1000) je
        write (nttyo,1000) je
1000 format (/' * Error - (EQLIB/adgexj) Programming error trap:',' The site',/7x,'number ',i12,' for a generic ion exchange',' phase',/7x,'is too large to fit in the 8-character field',' reserved for it',/7x,'in the default site name.')

        stop
    end if
end subroutine adgexj