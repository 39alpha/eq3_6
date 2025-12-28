subroutine hdrchk(ndat0s,noutpt,nttyo)
    !! This suboutine checks the first line of the DATA0 file to ensure
    !! that the mandatory header ("data0" beginning in column 1) is
    !! present.
    !! This suboutine is called by:
    !!   EQPT/eqpt.f
    !! Principal input:
    !!   ndat0s = unit number of the stripped DATA0 file
    !! Principal output:
    !!   None
    implicit none

    ! Calling sequence variable declarations.
    integer :: ndat0s
    integer :: noutpt
    integer :: nttyo

    ! Local variable declarations.
    integer :: j2

    integer :: ilnobl

    character(len=80) :: uline

    ! Check the DATA0 file header. This insures that the file that
    ! is supposed to be a DATA0 file really is one.
    read (ndat0s,1000,end=990,err=995) uline
1000 format(a)

    if (uline(1:5).ne.'data0' .and. uline(1:5).ne.'Data0' .and.  uline(1:5).ne.'DATA0') then
        j2 = ilnobl(uline)
        j2 = min(j2,50)
        write (noutpt,1010) uline(1:j2)
        write (nttyo,1010) uline(1:j2)
1010 format(/' * Error - (EQPT/hdrchk) The DATA0 must have "data0"',' beginning in',/7x,'column 1 of the first line. The first',' line of this file begins',/7x,'instead with: "',a,'".')

        stop
    end if

    rewind(ndat0s)
    go to 999

    ! Write a message for any read error.
990 continue
    write (noutpt,2000)
    write (nttyo,2000)
2000 format(/' * Error - (EQPT/hdrchk) The DATA0 file is an empty',' file.')

    stop

995 continue
    write (noutpt,2010)
    write (nttyo,2010)
2010 format(/' * Error - (EQPT/hdrchk) Encountered a read format',' error while',/7x,'checking for the mandatory header on the',' first line of the',/7x,'DATA0 file.')

    stop

999 continue
end subroutine hdrchk