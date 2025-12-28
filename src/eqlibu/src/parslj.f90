subroutine parslj(nfldmx,nfldt,nlchmx,ufield,uline1,ulscr)
    !! This subroutine parses the input string uline1, based on the
    !! delimiter '|'. Character strings found between delimiters are
    !! copied into the output array ufield. The elements of this
    !! array are then left-justified.
    !! The scheme is exemplified as follows:
    !! |  field 1  |  field 2   |  field3  |
    !! where the first '|' is in column 1 or
    !!    field 1  |  field 2   |  field3  |
    !! The first field may or may not have a delimiter to its left.
    !! Each field must have a delimiter to its right. A field to
    !! the right of the last delimiter, if any, is ignored.
    !! This subroutine is called by:
    !!   EQLIBU/rdd1l.f
    !!   EQLIBU/rdd1lh.f
    !!   EQLIBU/rdd2l.f
    !!   EQLIBU/rdd2lh.f
    !! Input:
    !!   nfldmx = maximum number of fields
    !!   nlchmx = the dimension of the ufield array
    !!   uline1 = the line to be parsed
    !!   ulscr  = scratch character variable
    !! Output:
    !!   nfldt  = actual number of fields
    !!   ufield = array of fields
    !! Calling sequence variable declarations.
    !! Note: the character length of the variables uline1, ulstr, and
    !! ufield is actually nlchmx. The Sun Fortran compiler does not
    !! allow this to be directly specified.
    integer :: nfldmx
    integer :: nlchmx

    integer :: nfldt

    character(len=*) :: uline1
    character(len=*) :: ulscr
    character(len=*) :: ufield(nfldmx)

    ! Local variable declarations.
    integer :: i

    ! Do a simple parse.
    call parsln(nfldmx,nfldt,nlchmx,ufield,uline1,ulscr)

    ! Left-justify the elements of the ufield array.
    do i = 1,nfldt
        call lejust(ufield(i))
    end do
end subroutine parslj