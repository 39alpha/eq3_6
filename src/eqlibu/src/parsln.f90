subroutine parsln(nfldmx,nfldt,nlchmx,ufield,uline1,ulscr)
    !! This subroutine parses the input string uline1, based on the
    !! delimiter '|'. Character strings found between delimiters are
    !! copied into the output array ufield. The elements of this array
    !! are not justified. EQLIBU/parslj.f can be used to call this
    !! subroutine and follow up the parsing with left-justification.
    !! The scheme is exemplified as follows:
    !! |  field 1  |  field 2   |  field3  |
    !! where the first '|' is in column 1 or
    !!    field 1  |  field 2   |  field3  |
    !! The first field may or may not have a delimiter to its left.
    !! Each field must have a delimiter to its right. A field to
    !! the right of the last delimiter, if any, is ignored.
    !! This subroutine is called by:
    !!   EQLIBU/parslj.f
    !!   XCON3/rd3d7.f
    !!   XCON6/rd6d7.f
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
    integer :: icht
    integer :: i2
    integer :: nfld

    icht = 0
    nfldt = 0

    do i = 1,nfldmx
        ufield(i) = ' '
    end do

    ! Find the first delimiter. This is normally expected to appear in
    ! column 1. The first field then begins in column 2. If the first
    ! delimiter appears in a column beyond the first, the first
    ! field begins in column 1.
    i  = index(uline1,'|')

    if (i .eq. 1) then
        icht = 1
    else
        icht = 0
    end if

    do nfld = 1,nfldmx
        ulscr = uline1(icht + 1:nlchmx)
        i = index(ulscr,'|')

        if (i .eq. 0) then
            ! Have a last field with no delimiter to its right. Ignore it.
            go to 999
        else
            i2 = i - 1
        end if

        nfldt = nfldt + 1

        if (nfldt .le. nfldmx) then
            if (i2 .gt. 0) then
                ufield(nfld)(1:i2) = ulscr(1:i2)
            end if
        end if

        icht = icht + i2 + 1

        if (icht .ge. nlchmx) then
            go to 999
        end if
    end do

999 continue
end subroutine parsln