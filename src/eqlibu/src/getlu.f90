subroutine getlu(nlu,nerr)
    !! This subroutine finds a currently unused unit number.
    !! This subroutine is called by:
    !!   EQLIBU/openin.f
    !!   EQLIBU/openou.f
    !! Input:
    !!   None
    !! Output:
    !!   nlu    = first currently unused unit number
    !!   nerr   = error flag:
    !!              = 0   Okay
    !!              = 1   Error
    implicit none

    ! Calling sequence variable declarations.
    integer :: nlu
    integer :: nerr

    ! Local variable declarations.
    integer :: iumax
    integer :: iumin

    logical :: qopen

    data iumax,iumin /40,0/

    nerr = 0

    ! Loop through all valid file numbers, beginning with the largest.
    do nlu = iumax,iumin,-1
        inquire(unit=nlu,opened=qopen)

        if (.not.qopen) then
            go to 999
        end if
    end do

    nerr = 1

999 continue
end subroutine getlu