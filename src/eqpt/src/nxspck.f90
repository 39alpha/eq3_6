subroutine nxspck(iktmax,issot,nerr,noutpt,nttyo,nxt,nxtmax,ussoph,ussosp)
    !! Check the names of solid solution end-members for uniqueness
    !! within each solid solution.
    !! This subroutine is called by:
    !!   EQPT/eqpt.f
    !! Principal input:
    !!   ussoph = array containing the names of solid solutions
    !!   ussosp = array containing the names of solid solution
    !!              end-members
    !! Principal output:
    !!   nerr   = incremented error counter
    implicit none

    ! Calling sequence variable declarations.
    integer :: noutpt
    integer :: nttyo

    integer :: iktmax
    integer :: nxtmax

    integer :: nxt
    integer :: nerr

    integer :: issot(nxtmax)

    character(len=24) :: ussoph(nxtmax)
    character(len=24) :: ussosp(iktmax,nxtmax)

    ! Local variable declarations.
    character(len=24) :: unam
    character(len=8) :: ux8

    integer :: i
    integer :: ikt
    integer :: ilist
    integer :: j
    integer :: j2
    integer :: j3
    integer :: j4
    integer :: k
    integer :: nx
    integer :: ncount

    integer :: ilnobl

    ! Loop over all solid solution phases.
    do nx = 1,nxt
        ikt = issot(nx)

        ! Note: ilist = the number of end-members that are duplicated
        ! for a given solid solution.
        ilist = 0

        ! Check each name for uniqueness.
        do i = 1,ikt
            unam = ussosp(i,nx)
            ncount = 1

            if (unam(1:7) .ne. '<blank>') then
                do j = i + 1,ikt
                    if (unam(1:24) .eq. ussosp(j,nx)(1:24)) then
                        ! Have found a duplicate entry moving down the list from
                        ! the entry now being tested for uniqueness. For the first
                        ! such duplicate entry only, make sure that there is not a
                        ! duplicate entry moving up the list from the entry now
                        ! being tested. If there is such a entry, then the
                        ! duplication in question has been noted previously and
                        ! should not be noted again.
                        if (ncount .eq. 1) then
                            do k = 1,i - 1
                                if (unam(1:24) .eq. ussosp(k,nx)(1:24)) then
                                    ! Have found a duplicate entry preceding the entry
                                    ! now being tested for uniqueness.
                                    go to 100
                                end if
                            end do
                        end if

                        ! Have found a duplication that has not been previously
                        ! noted.
                        ncount = ncount + 1
                    end if
                end do

100 continue
            end if

            if (ncount .gt. 1) then
                ! Have a duplicate entry or entries.
                ilist = ilist + 1

                if (ilist .eq. 1) then
                    ! Write error message header.
                    j2 = ilnobl(ussoph(nx))
                    write (noutpt,1000) ussoph(nx)(1:j2)
                    write (nttyo,1000) ussoph(nx)(1:j2)
1000 format(/' * Error - (EQPT/nxspck) The following',' end-members are listed more than',/7x,'once for solid',' solution ',a,':',/)
                end if

                write (ux8,'(i5)') ncount
                call lejust(ux8)
                j4 = ilnobl(ux8)
                j3 = ilnobl(unam)
                write (noutpt,1010) unam(1:j3),ux8(1:j4)
                write (nttyo,1010) unam(1:j3),ux8(1:j4)
1010 format(9x,a,' has ',a,' entries')
            end if
        end do

        nerr = nerr + ilist
    end do
end subroutine nxspck