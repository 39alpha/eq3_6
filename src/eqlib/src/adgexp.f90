subroutine adgexp(ne,noutpt,nttyo,ugexpd)
    !! This subroutine assigns a default name to the ne-th generic ion
    !! exchange phase. The general model is "Exchanger(n)", where n
    !! is the index number ne.
    !! This subroutine is called by:
    !!   EQLIB/intexi.f
    !!   EQ3NR/intge3.f
    !! Principal input:
    !!   ne     = index number of the current generic ion exchange phase
    !! Principal output:
    !!   ugexpd = default name of the ne-th generic ion exchange phase
    implicit none

    ! Calling sequence variable declarations.
    integer :: noutpt
    integer :: nttyo

    integer :: ne

    character(len=24) :: ugexpd

    ! Local variable declarations.
    ! None
    ugexpd = ' '
    ugexpd(1:10) = 'Exchanger('

    if (ne .le. 9) then
        write (ugexpd(11:11),'(i1)') ne
        ugexpd(12:12) = ')'
    else if (ne .le. 99) then
        write (ugexpd(11:12),'(i2)') ne
        ugexpd(13:13) = ')'
    else if (ne .le. 999) then
        write (ugexpd(11:13),'(i3)') ne
        ugexpd(14:14) = ')'
    else if (ne .le. 9999) then
        write (ugexpd(11:14),'(i4)') ne
        ugexpd(15:15) = ')'
    else if (ne .le. 99999) then
        write (ugexpd(11:15),'(i5)') ne
        ugexpd(16:16) = ')'
    else if (ne .le. 999999) then
        write (ugexpd(11:16),'(i6)') ne
        ugexpd(17:17) = ')'
    else if (ne .le. 9999999) then
        write (ugexpd(11:17),'(i7)') ne
        ugexpd(18:18) = ')'
    else if (ne .le. 99999999) then
        write (ugexpd(11:18),'(i8)') ne
        ugexpd(19:19) = ')'
    else
        write (noutpt,1000) ne
        write (nttyo,1000) ne
1000 format (/' * Error - (EQLIB/adgexp) Programming error trap:',' The index',/7x,'number ',i12,' of a generic ion exchange',' phase',/7x,'is too large to fit in the 8-character field',' reserved for it',/7x,'in the default site name.')

        stop
    end if
end subroutine adgexp