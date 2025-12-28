subroutine rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)
    !! This subroutine reads one line containing an expected header
    !! or tag string from an EQ3/6 input in menu-style ("D")
    !! format. The line generally contains both the header and data.
    !! EQLIBU/rdd2lh.f is like the present subroutine, but also reads
    !! a second line, which is a separator containing a field of dashes
    !! between two separators. EQLIBU/rdd1l.f is also very much like
    !! the present subroutine. It differs in that there is no expected
    !! header.
    !! This subroutine is called by:
    !!   XCON6/rd6d7.f
    !! Input:
    !!   nfldmx = the maximum number of fields on line uline1
    !!   nfldtx = the expected number of fields on line uline1
    !!   ninpts = the unit number of the stripped input file
    !!   nlchmx = the maximum number of characters on line uline1
    !!   nttyo  = the unit number of the screen file
    !!   uheadx = the expected header (the expected contents of the
    !!              first of the two lines)
    !!   ulscr  = scratch character variable
    !! Output:
    !!   nfldt  = the number of fields on line uline1
    !!   qrderr = logical flag, true if an error was encountered
    !!   ufield = the array of fields on line uline1
    !!   uheadr = the header (contents of the first of the two lines)
    !!   uline1 = the line to be read
    implicit none

    ! Calling sequence variable declarations.
    integer :: nfldmx
    integer :: nlchmx

    integer :: nttyo

    integer :: nfldt
    integer :: nfldtx
    integer :: ninpts

    logical :: qrderr

    character(len=*) :: ufield(nfldmx)
    character(len=*) :: uline1
    character(len=*) :: ulscr
    character(len=*) :: uheadr
    character(len=*) :: uheadx

    ! Local variable declarations.
    integer :: j
    integer :: j2
    integer :: j3

    integer :: ilnobl

    qrderr = .false.

    ! Read and parse the first line.
    read (ninpts,1000,err=990) uline1
1000 format(a80)

    call parslj(nfldmx,nfldt,nlchmx,ufield,uline1,ulscr)

    ! Compare the number of fields found with that expected.
    if (nfldtx .gt. 0) then
        if (nfldt .ne. nfldtx) then
            j3 = ilnobl(uline1(1:70))
            write (nttyo,1010) nfldt,nfldtx,uline1(1:j3)
1010 format(/' * Warning - (EQLIBU/rdd1lh) Found ',i2,' fields',/7x,'where ',i2,' were expected on the line beginning with',/7x,'"',a,'".')
        end if
    end if

    ! Compare the header with that expected.
    uheadr = ufield(1)
    call locase(uheadr)
    call locase(uheadx)

    j2 = ilnobl(uheadx)
    j = index(uheadr,uheadx(1:j2))

    if (j .eq. 0) then
        qrderr = .true.
        j3 = ilnobl(uline1(1:70))
        j2 = ilnobl(uheadx(1:70))
        write (nttyo,1020) uheadx(1:j2),uline1(1:j3)
1020 format(/' * Error - (EQLIBU/rdd1lh) Was expecting to find the',/7x,'header beginning with',/7x,'"',a,'"',/7x,'on the line beginning with',/7x,'"',a,'".')

        go to 999
    end if

    go to 999

990 continue
    qrderr = .true.

999 continue
end subroutine rdd1lh