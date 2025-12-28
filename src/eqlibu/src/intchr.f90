subroutine intchr(ivar,noutpt,nttyo,ustr)
    !! This subroutine converts writes the integer variable ivar into the
    !! character variable ustr, employing left justification and blank
    !! fill. If right justification is desired, one should just write
    !! the integer into the string.
    !! The length of the string variable is unknown. A buffer variable
    !! which is employed has a character length of ichpar. This limits
    !! the size of the substring available for the input integer.
    !! This subroutine is called by:
    !!   Any
    !! Input:
    !!   ivar   = the input integer
    !!   noutpt = unit number of the output file
    !!   nttyo  = unit number of the screen file
    !! Output:
    !!   ustr   = the string variable in which the integer is to
    !!              be written
    implicit none

    ! Calling sequence variable declarations.
    integer :: noutpt
    integer :: nttyo
    integer :: ivar

    character(len=*) :: ustr

    ! Local parameter declarations.
    integer :: ichpar

    parameter (ichpar = 16)

    ! Local variable declarations.
    integer :: j2
    integer :: ilen

    integer :: ilnobl

    character(len=ichpar) :: ux

    ! Write the integer into a character variable buffer.
    write (ux,1000,err = 100) ivar

    ! The format below should specify an integer field which matches
    ! the character length of the buffer variable (ichpar).
1000 format(i16)

    go to 110

100 continue
    write (noutpt,1010) ivar
    write (nttyo ,1010) ivar
1010 format(/" * Error - (EQLIBU/intchr) Can't write the integer",/7x,'"',i16,'" into a character string due to formatted',/7x,'write error.')

    stop

    ! Left justify.
110 continue
    call lejust(ux)

    ! Get the length of the integer string.
    j2 = ilnobl(ux)

    ! Get the length of the string variable.
    ilen = len(ustr)

    if (j2 .gt. ilen) then
        write (noutpt,1020) ivar,j2,ilen
        write (nttyo ,1020) ivar,j2,ilen
1020 format(/" * Error - (EQLIBU/intchr) Can't write the integer",/7x,'"',i16,'" into a string because the integer requires',i3,/7x,' characters, but the string variable has only ',i3,'.')

        stop
    end if

    ! Load the string.
    ustr = ux(1:j2)
end subroutine intchr