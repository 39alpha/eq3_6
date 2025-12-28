subroutine rijust(ustr)
    !! This subroutine right-justifies the non-blank portion of the
    !! string ustr.
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
    integer :: j
    integer :: jj
    integer :: jbl
    integer :: j2
    integer :: nchars

    integer :: ilnobl

    ! Get the length of the string variable.
    nchars = len(ustr)

    ! Get the position of the last non-blank character and the number
    ! of blanks on the right-hand-side.
    j2 = ilnobl(ustr)
    jbl = nchars - j2

    if (jbl .gt. 0) then
        do jj = j2,1,-1
            j = jj + jbl
            ustr(j:j) = ustr(jj:jj)
        end do

        do j = 1,jbl
            ustr(j:j) = ' '
        end do
    end if
end subroutine rijust