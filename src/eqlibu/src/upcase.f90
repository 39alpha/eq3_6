subroutine upcase(ustr)
    !! This subroutine converts a string from lower case to upper case.
    !! This subroutine is called by:
    !!   Any
    !! Input:
    !!   ustr   = the input string variable
    !! Output:
    !!   ustr   = the output string variable
    implicit none

    ! Calling sequence variable declarations.
    character(len=*) :: ustr

    ! Local variable declarations.
    integer :: idel
    integer :: j
    integer :: nchars

    character(len=1) :: u1

    idel = ichar('A') - ichar('a')

    if (idel .ne. 0) then
        nchars = len(ustr)

        do j = 1,nchars
            u1 = ustr(j:j)

            if (u1.ge.'a' .and. u1.le.'z') then
                ustr(j:j) = char(ichar(u1) + idel)
            end if
        end do
    end if
end subroutine upcase