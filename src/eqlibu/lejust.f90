subroutine lejust(ustr)
    !! This subroutine left-justifies the non-blank portion of the string
    !! ustr.
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
    integer :: j1
    integer :: nchars

    integer :: ifnobl

    ! Get the length of the string variable.
    nchars = len(ustr)

    ! Get the position of the first non-blank character and the number
    ! of blanks on the left-hand-side.
    j1 = ifnobl(ustr)
    jbl = j1 - 1

    if (jbl .gt. 0) then
        do jj = j1,nchars
            j = jj - jbl
            ustr(j:j) = ustr(jj:jj)
        end do

        do j = nchars - jbl + 1,nchars
            ustr(j:j) = ' '
        end do
    end if
end subroutine lejust