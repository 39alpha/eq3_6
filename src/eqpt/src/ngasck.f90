subroutine ngasck(nerr,ngt,ngtmax,noutpt,nttyo,ugassp)
    !! Check the names of gas species for uniqueness.
    !! This subroutine is called by:
    !!   EQPT/eqpt.f
    !! Principal input:
    !!   ugassp = array containing the names of gas species
    !! Principal output:
    !!   nerr   = incremented error counter
    implicit none

    ! Calling sequence variable declarations.
    integer :: noutpt
    integer :: nttyo

    integer :: ngtmax

    integer :: ngt
    integer :: nerr

    character(len=24) :: ugassp(ngtmax)

    ! Local variable declarations.
    character(len=24) :: unam
    character(len=8) :: ux8

    integer :: i
    integer :: ilist
    integer :: j
    integer :: j2
    integer :: j3
    integer :: k
    integer :: ncount

    integer :: ilnobl

    ! Note: ilist = the number of species with duplicate species blocks
    ! present on the data file.
    ilist = 0

    ! Check each name for uniqueness.
    do i = 1,ngt
        unam = ugassp(i)
        ncount = 1

        if (unam(1:7) .ne. '<blank>') then
            do j = i + 1,ngt
                if (unam(1:24) .eq. ugassp(j)(1:24)) then
                    ! Have found a duplicate block moving down the list from
                    ! the block now being tested for uniqueness. For the first
                    ! such duplicate block only, make sure that there is not a
                    ! duplicate block moving up the list from the block now
                    ! being tested. If there is such a block, then the
                    ! duplication in question has been noted previously and
                    ! should not be noted again.
                    if (ncount .eq. 1) then
                        do k = 1,i - 1
                            if (unam(1:24) .eq. ugassp(k)(1:24)) then
                                ! Have found a duplicate block preceding the block
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
            ! Have a duplicate block or blocks.
            ilist = ilist + 1

            if (ilist .eq. 1) then
                ! Write error message header.
                write (noutpt,1000)
                write (nttyo,1000)
1000 format(/' * Error - (EQPT/ngasck) The following gas',' species have multiple',/7x,'species blocks:',/)
            end if

            write (ux8,'(i5)') ncount
            call lejust(ux8)
            j3 = ilnobl(ux8)
            j2 = ilnobl(unam)
            write (noutpt,1010) unam(1:j2),ux8(1:j3)
            write (nttyo,1010) unam(1:j2),ux8(1:j3)
1010 format(9x,a,' has ',a,' blocks')
        end if
    end do

    nerr = nerr + ilist
end subroutine ngasck