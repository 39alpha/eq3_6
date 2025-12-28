subroutine openin(noutpt,nttyo,ufiln,uform,ilu)
    !! This subroutine opens an input type file. The file must already
    !! exist. An unused logical unit number is obtained.
    !! This subroutine is called by:
    !!   EQ3NR/eq3nr.f
    !!   EQ6/eq6.f
    !! Input:
    !!   noutpt = the unit number of the output file
    !!   nttyo  = the unit number of the screen file
    !!   ufiln  = the name of the file to open
    !!   uform  = the file format, 'formatted' or 'unformatted'
    !! Output:
    !!   ilu     = the logical unit number of opened file
    implicit none

    ! Calling sequence variable declarations.
    ! Note- the character length for the variables ufiln and uform
    ! can not be exactly specified.
    integer :: ilu
    integer :: noutpt
    integer :: nttyo

    character(len=*) :: ufiln
    character(len=*) :: uform

    ! Local variable declarations.
    integer :: j2
    integer :: nerr

    integer :: ilnobl

    logical :: qex

    character(len=8) :: uformo

    j2 = ilnobl(ufiln)

    ! Check to make sure the file exists.
    inquire(file=ufiln,exist=qex,formatted=uformo)

    if (.not.qex) then
        if (noutpt .gt. 0) then
            write (noutpt,1000) ufiln(1:j2)
        end if

        write (nttyo,1000) ufiln(1:j2)
1000 format(/' * Error - (EQLIBU/openin) The file "',a,'"',/7x,'does not exist.')

        stop
    end if

    ! Check the file format.
    if (uform.eq.'formatted' .and. uformo.eq.'no') then
        if (noutpt .gt. 0) then
            write (noutpt,1010) ufiln(1:j2)
        end if

        write (nttyo,1010) ufiln(1:j2)
1010 format(/' * Error - (EQLIBU/openin) The file "',a,'"',/7x,'should be formatted, but it is not.')

        stop
    end if

    if (uform.eq.'unformatted' .and. uformo.eq.'yes') then
        !        if (noutpt .gt. 0) write (noutpt,1020) ufiln(1:j2)
        !        write (nttyo,1020) ufiln(1:j2)
        ! 1020   format(/' * Error - (EQLIBU/openin) The file "',a,'"',
        !     $  /7x,'should be unformatted, but it is not.')
        !      stop
    end if

    ! Get the next available logical unit number.
    call getlu(ilu,nerr)

    if (nerr .ne. 0) then
        if (noutpt .gt. 0) then
            write (noutpt,1050) ufiln(1:j2)
        end if

        write (nttyo,1050) ufiln(1:j2)
1050 format(/' * Error - (EQLIBU/openin) No logical unit number',/7x,'is available for the file "',a,'".')

        stop
    end if

    open(ilu,file=ufiln,form=uform,status='old',err=10)
    go to 999

10 continue

    if (noutpt .gt. 0) then
        write (noutpt,1060) ufiln(1:j2)
    end if

    write (nttyo,1060) ufiln(1:j2)
1060 format(/" * Error - (EQLIBU/openin) Can't open the file",/7x,'"',a,'".')

    stop

999 continue
end subroutine openin
