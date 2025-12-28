subroutine neleck(nct,nctmax,nerr,noutpt,nttyo,uelem)
    !! Check the names of chemical elements for uniqueness.
    !! This subroutine is called by:
    !!   EQPT/eqpt.f
    !! Principal input:
    !! Principal output:
    !!   nerr   = incremented error counter
    implicit none

    ! Calling sequence variable declarations.
    integer :: nctmax

    integer :: noutpt
    integer :: nttyo

    integer :: nct
    integer :: nerr

    character(len=8) :: uelem(nctmax)

    ! Local variable declarations.
    character(len=8) :: unam
    character(len=8) :: ux8

    integer :: i
    integer :: ilist
    integer :: j
    integer :: j2
    integer :: j3
    integer :: k
    integer :: ncount

    integer :: ilnobl

    ilist = 0

    ! Check each name for uniqueness.
    do i = 1,nct
        unam = uelem(i)
        ncount = 1

        if (unam(1:7) .ne. '<blank>') then
            do j = i + 1,nct
                if (unam(1:8) .eq. uelem(j)(1:8)) then
                    ! Have found a duplicate entry moving down the list from
                    ! the entry now being tested for uniqueness. For the first
                    ! such duplicate entry only, make sure that there is not a
                    ! duplicate entry moving up the list from the entry now
                    ! being tested. If there is such a entry, then the
                    ! duplication in question has been noted previously and
                    ! should not be noted again.
                    if (ncount .eq. 1) then
                        do k = 1,i - 1
                            if (unam(1:8) .eq. uelem(k)(1:8)) then
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
                write (noutpt,1000)
                write (nttyo,1000)
1000 format(/' * Error - (EQPT/neleck) The following chemical',' elements have multiple',/7x,'entries in the elements',' block:',/)
            end if

            write(ux8,'(i5)') ncount
            call lejust(ux8)
            j3 = ilnobl(ux8)
            j2 = ilnobl(unam)
            write (noutpt,1010) unam(1:j2),ux8(1:j3)
            write (nttyo,1010) unam(1:j2),ux8(1:j3)
1010 format(9x,a,' has ',a,' entries')
        end if
    end do

    nerr = nerr + ilist
end subroutine neleck