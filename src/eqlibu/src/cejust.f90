subroutine cejust(ustr)
    !! This subroutine centers the non-blank portion of the string ustr.
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
    integer :: jbll
    integer :: jblr
    integer :: jshl
    integer :: jshr
    integer :: j1
    integer :: j2
    integer :: nchars

    integer :: ifnobl
    integer :: ilnobl

    ! Get the length of the string variable.
    nchars = len(ustr)

    ! Get the position of the first non-blank character and the number
    ! of blanks on the left-hand-side.
    j1 = ifnobl(ustr)
    jbll = j1 - 1

    ! Get the position of the last non-blank character and the number
    ! of blanks on the right-hand-side.
    j2 = ilnobl(ustr)
    jblr = nchars - j2

    jbl = (jbll + jblr)/2
    jshl = jbll - jbl

    if (jshl .gt. 0) then
        ! Shift left.
        do jj = j1,nchars
            j = jj - jshl
            ustr(j:j) = ustr(jj:jj)
        end do

        do j = nchars - jshl,nchars
            ustr(j:j) = ' '
        end do
    end if

    jshr = - jshl

    if (jshr .gt. 1) then
        do jj = j2,1,-1
            j = jj + jshr
            ustr(j:j) = ustr(jj:jj)
        end do

        do j = 1,jshr
            ustr(j:j) = ' '
        end do
    end if
end subroutine cejust