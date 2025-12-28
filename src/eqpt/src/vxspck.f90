subroutine vxspck(iktmax,issot,nerr,nmt,nmtmax,noutpt,nttyo,nxt,nxtmax,uminsp,ussoph,ussosp)
    !! Validate the names of solid solution end-members. Make sure that
    !! these names appear on the list of pure minerals.
    !! This subroutine is called by:
    !!   EQPT/eqpt.f
    !! Principal input:
    !!   nxt    = the number of solid solutions
    !!   issot  = array of numbers of end-members of solid solutions
    !!   ussoph = array of names of solid solutions
    !!   ussosp = array of names of end-members of solid solutions
    !!   uminsp = array of pure mineral names
    !! Principal output:
    !!   nerr   = error counter
    implicit none

    ! Calling sequence variable declarations.
    integer :: iktmax
    integer :: nmtmax
    integer :: nxtmax

    integer :: noutpt
    integer :: nttyo

    integer :: nerr
    integer :: nmt
    integer :: nxt

    integer :: issot(nxtmax)

    character(len=24) :: uminsp(nmtmax)
    character(len=24) :: ussoph(nxtmax)
    character(len=24) :: ussosp(iktmax,nxtmax)

    ! Local variable declarations.
    character(len=24) :: unam

    integer :: i
    integer :: ikt
    integer :: j2
    integer :: j3
    integer :: kerr
    integer :: nm
    integer :: nx

    integer :: ilnobl

    ! Loop over all solid solution phases.
    do nx = 1,nxt
        ikt = issot(nx)
        kerr = 0

        ! Check that each end-member appears on the list of pure minerals.
        do i = 1,ikt
            unam = ussosp(i,nx)

            if (unam(1:7) .ne. '<blank>') then
                ! Check that each end-member appears on the list of
                ! pure minerals.
                do nm = 1,nmt
                    if (unam(1:24) .eq. uminsp(nm)(1:24)) then
                        go to 100
                    end if
                end do

                ! Did not find this end-member on the list of pure minerals.
                kerr = kerr + 1
                nerr = nerr + 1

                if (kerr .eq. 1) then
                    ! Write a header for the solid solution.
                    j2 = ilnobl(ussoph(nx))
                    write (noutpt,1000) ussoph(nx)(1:j2)
                    write (nttyo,1000) ussoph(nx)(1:j2)
1000 format(/' * Error - (EQPT/vxspck) The following',' end-members of solid solution',/7x,a,' are not present',' on the data file as pure minerals:',/)
                end if

                j3 = ilnobl(unam)
                write (noutpt,1010) unam(1:j3)
                write (nttyo,1010) unam(1:j3)
1010 format(9x,a)

100 continue
            end if
        end do
    end do
end subroutine vxspck