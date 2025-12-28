subroutine nelcck(nctmax,ncts,nentei,nerr,qdupes,uessi)
    !! Check the elemental composition of a species to ensure that
    !! each chemical element name appearing in the composition is unique.
    !! This subroutine is called by:
    !!   EQPT/elesck.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: nctmax

    integer :: nentei(nctmax)

    integer :: ncts
    integer :: nerr

    logical :: qdupes

    character(len=8) :: uessi(nctmax)

    ! Local variable declarations.
    character(len=8) :: unam

    integer :: i
    integer :: j
    integer :: k
    integer :: ncount

    qdupes = .false.

    ! Check each name for uniqueness.
    do i = 1,ncts
        unam = uessi(i)
        ncount = 1

        if (unam(1:7) .ne. '<blank>') then
            do j = i + 1,ncts
                if (unam(1:8) .eq. uessi(j)(1:8)) then
                    ! Have found a duplicate entry moving down the list from
                    ! the entry now being tested for uniqueness. For the first
                    ! such duplicate entry only, make sure that there is not a
                    ! duplicate entry moving up the list from the entry now
                    ! being tested. If there is such a entry, then the
                    ! duplication in question has been noted previously and
                    ! should not be noted again.
                    if (ncount .eq. 1) then
                        do k = 1,i - 1
                            if (unam(1:8) .eq. uessi(k)(1:8)) then
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

        nentei(i) = ncount

        if (ncount .gt. 1) then
            qdupes = .true.
        end if
    end do
end subroutine nelcck