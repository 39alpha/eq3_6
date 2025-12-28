integer function ifnobl(ustr)
    !! This subroutine finds the position of the first non-blank
    !! character in the string ustr.
    !! This subroutine is called by:
    !!   Any
    !! Input:
    !!   ustr   = the input string variable
    !! Output:
    !!   ifnobl = the position of the first non-blank character
    implicit none

    ! Calling sequence variable declarations.
    character(len=*) :: ustr

    ! Local variable declarations.
    integer :: j
    integer :: nchars

    ! Get the length of the string variable.
    nchars = len(ustr)

    ! Find the first non-blank character.
    ifnobl = 0

    do j = 1,nchars
        if (ustr(j:j) .ne. ' ') then
            ifnobl = j
            go to 999
        end if
    end do

999 continue
end function ifnobl