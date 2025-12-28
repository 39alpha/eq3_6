subroutine rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,nttyo,qrderr,ufield,uline1,ulscr)
    !! This subroutine reads single line containing an unknown header
    !! or tag string from an EQ3/6 input in menu-style ("D")
    !! format. The line generally contains both the header and data.
    !! EQLIBU/rdd2l.f is much like the present subroutine, but also reads
    !! a second line, which is a separator containing a field of dashes
    !! between two separators. EQLIBU/rdd1lh.f is also much like the
    !! present subroutine. It differs in that it has an expected header.
    !! This subroutine is called by:
    !!   XCON3/rd3d7.f
    !!   XCON6/rd6d7.f
    !! Input:
    !!   ninpts = the unit number of the stripped input file
    !!   nfldmx = the maximum number of fields on line uline1
    !!   nfldtx = the expected number of fields on line uline1
    !!   nlchmx = the maximum number of characters on line uline1
    !!   nttyo  = the unit number of the screen file
    !!   ulscr  = scratch character variable
    !! Output:
    !!   nfldt  = the number of fields on line uline1
    !!   qrderr = logical flag, true if an error was encountered
    !!   ufield = the array of fields on line uline1
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

    ! Local variable declarations.
    integer :: j2

    integer :: ilnobl

    qrderr = .false.

    ! Read and parse the first line.
    read (ninpts,1000,err=990) uline1
1000 format(a80)

    call parslj(nfldmx,nfldt,nlchmx,ufield,uline1,ulscr)

    ! Compare the number of fields found with that expected.
    if (nfldtx .gt. 0) then
        if (nfldt .ne. nfldtx) then
            j2 = ilnobl(uline1(1:70))
            write (nttyo,1010) nfldt,nfldtx,uline1(1:j2)
1010 format(/' * Warning - (EQLIBU/rdd1l) Found ',i2,' fields',/7x,'where ',i2,' were expected on the line beginning with',/7x,'"',a,'".')
        end if
    end if

    go to 999

990 continue
    qrderr = .true.

999 continue
end subroutine rdd1l