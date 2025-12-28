subroutine realch(noutpt,nttyo,ustr,var)
    !! This subroutine converts writes the real*8 variable var into the
    !! string ustr, employing left justification and blank fill. If
    !! right justification is desired, one should just write the
    !! real*8 variable into the string.
    !! The length of the string variable is unknown. A buffer variable
    !! which is employed has a character length of ichpar. This limits
    !! the size of the substring available for the input real*8 number.
    !! This subroutine is called by:
    !!   Any
    !! Input:
    !!   noutpt = unit number of the output file
    !!   nttyo  = unit number of the screen file
    !!   var   = the input real*8 number
    !! Output:
    !!   ustr   = the string variable in which the integer is to
    !!              be written
    implicit none

    ! Calling sequence variable declarations.
    integer :: noutpt
    integer :: nttyo

    character(len=*) :: ustr

    real(kind=8) :: var

    ! Local parameter declarations.
    integer :: ichpar

    parameter (ichpar = 24)

    ! Local variable declarations.
    integer :: ilen
    integer :: j2

    integer :: ilnobl

    character(len=ichpar) :: ux

    ! Write the real*8 number into a character variable buffer.
    write (ux,1000,err = 100) var

    ! The format below should specify a real*8 field which matches
    ! the character length of the buffer variable (ichpar).
1000 format(g24.4)

    go to 110

100 continue
    write (noutpt,1010) var
    write (nttyo ,1010) var
1010 format(/" * Error - (EQLIBU/realch) Can't write the real*8",/7x,'number ",g24.4," into a character string due to formatted',/7x,'write error.')

    stop

    ! Left justify.
110 continue
    call lejust(ux)

    ! Get the length of the real*8 string.
    j2 = ilnobl(ux)

    ! Get the length of the string variable.
    ilen = len(ustr)

    if (j2 .gt. ilen) then
        write (noutpt,1020) var,j2,ilen
        write (nttyo ,1020) var,j2,ilen
1020 format(/" * Error - (EQLIBU/realch) Can't write the real*8",/7x,'number "',g24.4,'" into a string because the number',/7x,'requires',i3,' characters. The string variable',' has only ',i3,' characters.')

        stop
    end if

    ! Load the string.
    ustr = ux(1:j2)
end subroutine realch