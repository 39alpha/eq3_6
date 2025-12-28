integer function ilnobl(ustr)
    !! This subroutine finds the position of the last non-blank character
    !! in the string ustr.
    !! This subroutine is called by:
    !!   Any
    !! Input:
    !!   ustr   = the input string variable
    !! Output:
    !!   ilnobl = the position of the first non-blank character
    implicit none

    ! Calling sequence variable declarations.
    character(len=*) :: ustr

    ! Local variable declarations.
    integer :: j
    integer :: nchars

    ! Get the length of the string variable.
    nchars = len(ustr)

    ! Find the first non-blank character.
    ilnobl = 0

    do j = nchars,1,-1
        if (ustr(j:j) .ne. ' ') then
            ilnobl = j
            go to 999
        end if
    end do

999 continue
end function ilnobl