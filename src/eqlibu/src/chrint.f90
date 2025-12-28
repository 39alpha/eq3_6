subroutine chrint(ivar,nttyo,qrderr,ustr)
    !! This subroutine reads the integer ivar from the character string
    !! ustr. The string may contain non-blank characters other than the
    !! integer. If there is more than one integer in the string, only the
    !! first is read. If the string is completely blank, the integer is
    !! returned with a value of zero. If no integer is present and at
    !! least one blank is present, the integer is returned with a value
    !! of zero.
    !! The length of the string variable is unknown. A buffer variable
    !! which is employed has a character length of ichpar. This limits
    !! the size of the substring for the integer.
    !! This subroutine is called by:
    !!   Any
    !! Input:
    !!   nttyo  = unit number of the screen file
    !!   ustr   = the input string variable
    !! Output:
    !!   ivar   = the output integer
    !!   qrderr = logical flag, .true. if a read format error occurred
    implicit none

    ! Calling sequence variable declarations.
    integer :: nttyo
    integer :: ivar

    logical :: qrderr

    character(len=*) :: ustr

    ! Local parameter declarations.
    integer :: ichpar

    parameter (ichpar = 16)

    ! Local variable declarations.
    integer :: i
    integer :: ifirst
    integer :: ipos1
    integer :: istr
    integer :: ichars
    integer :: j2
    integer :: nchars

    integer :: ifnobl
    integer :: ilnobl

    character(len=ichpar) :: ux
    character(len=1) :: uc

    logical :: qsign

    qrderr = .false.

    ! Find the position of the first non-blank character in the
    ! buffer variable.
    ipos1 = ifnobl(ustr)

    if (ipos1 .le. 0) then
        ivar = 0
        go to 999
    end if

    ! Get the length of the string variable.
    nchars = len(ustr)

    ! Find the first integer substring. This must begin with one of the
    ! following characters: +-0123456789. The substring is considered
    ! found when a number is found. If the substring begins with a
    ! "+" or "-" sign, intervening blanks are permitted.
    ifirst = 0

    do i = ipos1,nchars
        uc(1:1) = ustr(i:i)

        if (uc.ge.'0' .and. uc.le.'9') then
            if (ifirst .eq. 0) then
                ifirst = i
            end if

            go to 100
        else if (uc .eq. '-') then
            ifirst = i
        else if (uc .eq. '+') then
            ifirst = i
        else if (uc .ne. ' ') then
            ifirst = 0
        end if
    end do

    ! No integer is present in the string. See if there is at
    ! least one blank present.
    i = index(ustr,' ')

    if (i .gt. 0) then
        ivar = 0
        go to 999
    end if

    j2 = ilnobl(ustr)
    write (nttyo ,1000) ustr(1:j2)
1000 format(/" * Error - (EQLIBU/chrint) Can't find an integer",/7x,'to read from the string "',a,'".')

    qrderr = .true.
    go to 999

    ! Copy the integer field into the buffer string.
    !   ichars = length of the substring containing the integer
    !   qsign  = logical flag denoting the reading of a substring
    !            containing an initial '-' or '+' sign followed
    !            by zero or more blanks
100 continue
    ichars = 0
    qsign = .false.

    do istr = ifirst,nchars
        uc(1:1) = ustr(istr:istr)

        if (uc.ge.'0' .and. uc.le.'9') then
            ! Have encountered a number character in the integer substring.
            ichars = ichars + 1
            ux(ichars:ichars) = uc
            qsign = .false.
        else if (uc.eq.'-' .or. uc.eq.'+') then
            ! Have encountered a '-' or '+' sign which begins the
            ! integer substring.
            ichars = ichars + 1
            ux(ichars:ichars) = uc
            qsign = .true.
        else if (qsign .and. uc.eq.' ') then
            ! Have encountered a blank between the '-' or '+' sign and
            ! the first non-blank character and the integer string.
            ! Allow it in the substring.
            ichars = ichars + 1
            ux(ichars:ichars) = uc
        else
            go to 110
        end if
    end do

    ! Load blank fill.
110 continue
    do i = ichars + 1,ichpar
        ux(i:i) = ' '
    end do

    ! Read the integer from the string buffer.
    read (ux,1010,err=995) ivar

    ! The format below should specify an integer field which matches
    ! the character length of the buffer variable (ichpar).
1010 format(i16)

    go to 999

995 continue
    j2 = ilnobl(ux)
    write (nttyo ,1020) ux(1:j2)
1020 format(/" * Error - (EQLIBU/chrint) Can't read an integer from",/7x,'the buffer string "',a,'" due to a formatted read error.')

    qrderr = .true.

999 continue
end subroutine chrint