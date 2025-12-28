subroutine openou(noutpt,nttyo,ufiln,uform,nrecl,ilu)
    !! This subroutine opens an output type file. If a file of the
    !! same name already exists, it is first destroyed. An unused
    !! logical unit number is obtained.
    !! This subroutine is called by:
    !!   EQ3NR/eq3nr.f
    !!   EQ6/eq6.f
    !! Input:
    !!   noutpt = the unit number of the output file
    !!   nttyo  = the unit number of the screen file
    !!   ufiln  = the name of the file to open
    !!   uform  = the file format, 'formatted' or 'unformatted'
    !!   nrecl  = the record length (number of characters per line)
    !! Output:
    !!   ilu     = the logical unit number of opened file
    implicit none

    ! Calling sequence variable declarations.
    integer :: noutpt
    integer :: nttyo

    integer :: ilu
    integer :: nrecl

    ! Note- the character length for the variables ufiln and uform
    ! can not be exactly specified.
    character(len=*) :: ufiln
    character(len=*) :: uform

    ! Local variable declarations.
    integer :: j2
    integer :: nerr

    integer :: ilnobl

    logical :: qex

    character(len=8) :: ustat

    j2 = ilnobl(ufiln)

    ! See if a file of the same name already exists. If so,
    ! destroy it. This makes the logical unit number available.
    ! If a file of the same name does not exist, get the next
    ! available logical unit number.
    inquire(file=ufiln,exist=qex)

    if (qex) then
        ustat = 'old'
        call getlu(ilu,nerr)

        if (nerr .ne. 0) then
            if (noutpt .gt. 0) then
                write (noutpt,1000) ustat,ufiln(1:j2)
            end if

            write (nttyo,1000) ustat,ufiln(1:j2)
1000 format(/' * Error - (EQLIBU/openou) No logical unit number',/7x,'is available to open the ',a3,' file "',a,'".')

            stop
        end if

        open(ilu,file=ufiln,status=ustat,err=10)
        close(ilu,status='delete',err=15)
    else
        ustat = 'new'
        call getlu(ilu,nerr)

        if (nerr .ne. 0) then
            if (noutpt .gt. 0) then
                write (noutpt,1000) ustat,ufiln(1:j2)
            end if

            write (nttyo,1000) ustat,ufiln(1:j2)
            stop
        end if
    end if

    ! Open the new file.
    if (nrecl .gt. 0) then
        ! Use the specified record length.
        open(ilu,file=ufiln,form=uform,status='new',recl=nrecl,err=10)
    else
        ! The record length is not specified. Open the file at the
        ! default record length.
        open(ilu,file=ufiln,form=uform,status='new',err=10)
    end if

    go to 999

10 continue
    if (noutpt .gt. 0) then
        write (noutpt,1010) ustat,ufiln(1:j2)
    end if

    write (nttyo ,1010) ustat,ufiln(1:j2)
1010 format(/" * Error - (EQLIBU/openou) Can't open the ",a3,' copy',/7x,'of the file "',a,'".')

    stop

15 continue
    if (noutpt .gt. 0) then
        write (noutpt,1020) ufiln(1:j2)
    end if

    write (nttyo,1020) ufiln(1:j2)
1020 format(/" * Error - (EQLIBU/openou) Can't delete the old copy",/7x,'of the file "',a,'".')

    stop

999 continue
end subroutine openou