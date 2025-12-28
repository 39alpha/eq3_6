subroutine nrxnck(nbtmx1,ndrsts,nentri,nerr,qduprs,udrsi)
    !! Check the associated raection of a species to ensure that
    !! each species name appearing in the reaction is unique.
    !! This subroutine is called by:
    !!   EQPT/rxnsck.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: nbtmx1

    integer :: nentri(nbtmx1)

    integer :: ndrsts
    integer :: nerr

    logical :: qduprs

    character(len=24) :: udrsi(nbtmx1)

    ! Local variable declarations.
    character(len=24) :: unam

    integer :: i
    integer :: j
    integer :: k
    integer :: ncount

    qduprs = .false.

    ! Check each name for uniqueness.
    do i = 1,ndrsts
        unam = udrsi(i)
        ncount = 1

        if (unam(1:7) .ne. '<blank>') then
            do j = i + 1,ndrsts
                if (unam(1:24) .eq. udrsi(j)(1:24)) then
                    ! Have found a duplicate entry moving down the list from
                    ! the entry now being tested for uniqueness. For the first
                    ! such duplicate entry only, make sure that there is not a
                    ! duplicate entry moving up the list from the entry now
                    ! being tested. If there is such a entry, then the
                    ! duplication in question has been noted previously and
                    ! should not be noted again.
                    if (ncount .eq. 1) then
                        do k = 1,i - 1
                            if (unam(1:24) .eq. udrsi(k)(1:24)) then
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

        nentri(i) = ncount

        if (ncount .gt. 1) then
            qduprs = .true.
        end if
    end do
end subroutine nrxnck