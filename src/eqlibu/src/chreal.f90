subroutine chreal(nttyo,qrderr,ustr,var)
    !! This subroutine reads the real*8 var from the character string
    !! ustr. The string may contain non-blank characters other than the
    !! real*8 number. If there is more than one real*8 number in the
    !! string, only the first is read. If the string is completely blank,
    !! the real*8 number is returned with a value of zero. If no real*8
    !! number is present and at least one blank is present, the
    !! real*8 number is returned with a value of zero.
    !! The length of the string variable is unknown. A buffer variable
    !! which is employed has a character length of ichpar. This limits
    !! the size of the substring for the real*8 number.
    !! This subroutine is called by:
    !!   Any
    !! Input:
    !!   nttyo  = unit number of the screen file
    !!   ustr   = the input string variable
    !! Output:
    !!   var    = the output real*8 number
    !!   qrderr = logical flag, .true. if a read format error occurred
    implicit none

    ! Calling sequence variable declarations.
    integer :: nttyo

    logical :: qrderr

    character(len=*) :: ustr

    real(kind=8) :: var

    ! Local parameter declarations.
    integer :: ichpar

    parameter (ichpar = 24)

    ! Local variable declarations.
    integer :: i
    integer :: iexpf
    integer :: ifirst
    integer :: ilettr
    integer :: inum1
    integer :: ipoint
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
        var = 0.
        go to 999
    end if

    ! Get the length of the string variable.
    nchars = len(ustr)

    ! Find the first real*8 substring. This must begin with one of the
    ! following characters: .+-0123456789. The substring is considered
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
        else if (uc .eq. '.') then
            if (ifirst .eq. 0) then
                ifirst = i
            end if
        else if (uc .ne. ' ') then
            ifirst = 0
        end if
    end do

    ! No real*8 number is present in the string. See if there is at
    ! least one blank present.
    i = index(ustr,' ')

    if (i .gt. 0) then
        var = 0.
        go to 999
    end if

    j2 = ilnobl(ustr)
    write (nttyo,1000) ustr(1:j2)
1000 format(/" * Error - (EQLIBU/chreal) Can't find a real*8 number",/7x,'to read from the string "',a,'".')

    qrderr = .true.
    go to 999

    ! Copy the real*8 field into the buffer string.
    !   ipoint = original string position containing '.'
    !   iexpf  = original string position of the first character
    !              in the exponent field, normally 'e','E', 'd',
    !              or 'D', but possibly '-' or '+'
    !   ilettr = original string position containing 'e','E',
    !              'd', or 'D'
    !   inum1  = number of number characters in the non-exponent
    !              field
    !   ichars = length of the substring containing the real*8
    !              number
    !   qsign  = logical flag denoting the reading of a substring
    !            containing an initial '-' or '+' sign followed
    !            by zero or more blanks
100 continue
    ipoint = 0
    iexpf = 0
    ilettr = 0
    inum1 = 0
    ichars = 0
    qsign = .false.

    do istr = ifirst,nchars
        uc(1:1) = ustr(istr:istr)

        if (uc.ge.'0' .and. uc.le.'9') then
            ! Have encountered a number character in either the non-exponent
            ! or exponent field.
            ichars = ichars + 1
            ux(ichars:ichars) = uc

            if (iexpf .eq. 0) then
                inum1 = inum1 + 1
            end if

            qsign = .false.
        else if (uc.eq.'-' .or. uc.eq.'+') then
            if (iexpf .eq. 0) then
                if (inum1 .eq. 0) then
                    ! Have encountered a '-' or '+' sign which begins the
                    ! non-exponent field.
                    ichars = ichars + 1
                    ux(ichars:ichars) = uc
                    qsign = .true.
                else
                    ! Have encountered a '-' or '+' sign which starts
                    ! the exponent field.
                    ichars = ichars + 1
                    ux(ichars:ichars) = uc
                    iexpf = istr
                    qsign = .false.
                end if
            else
                if (istr .eq. (ilettr + 1)) then
                    ! Have a '-' or '+' sign following an 'e', 'E',
                    ! 'd', or 'D' (in the exponent field).
                    ichars = ichars + 1
                    ux(ichars:ichars) = uc
                    qsign = .false.
                else
                    ! Have encountered a second '-' or '+' sign in the
                    ! exponent field. Terminate the substring without it.
                    go to 110
                end if
            end if
        else if (uc.eq.' ') then
            if (qsign) then
                ! Have encountered a blank between the initial '-' or
                ! '+' sign and the first non-blank character in the
                ! substring. Allow it in the substring.
                ichars = ichars + 1
                ux(ichars:ichars) = uc
            else
                ! Have encountered an illegal blank. Terminate the
                ! substring without it.
                go to 110
            end if
        else if (uc .eq. '.') then
            if (iexpf .eq. 0) then
                if (ipoint .eq. 0) then
                    ichars = ichars + 1
                    ux(ichars:ichars) = uc
                    ipoint = istr
                    qsign = .false.
                else
                    ! Have encountered a second '.' in the non-exponent
                    ! field. Terminate the substring without it.
                    go to 110
                end if
            else
                ! Have encountered a '.' in the exponent field.
                ! Terminate the substring without it.
                go to 110
            end if
        else if (uc.eq.'e' .or. uc.eq.'E' .or. uc.eq.'d'    .or. uc.eq.'D') then
            if (iexpf .eq. 0) then
                ! Have encountered an exponent character.
                ichars = ichars + 1
                ux(ichars:ichars) = uc
                ilettr = istr
                iexpf = istr
                qsign = .false.
            else
                ! Have encountered a second 'e', 'E', 'd', or 'D' in
                ! the exponent field. Terminate the substring without it.
                go to 110
            end if
        else
            ! Have encountered an illegal character. Terminate the
            ! substring without it.
            go to 110
        end if
    end do

    ! Load blank fill.
110 continue
    do i = ichars + 1,ichpar
        ux(i:i) = ' '
    end do

    if (inum1 .le. 0) then
        j2 = ilnobl(ux)
        write (nttyo,1010) ux(1:j2)
1010 format(/" * Warning - (EQLIBU/chreal) The real*8 number in",/7x,'the buffer string "',a,'" has no numbers in the',/7x,'non-exponent part.')
    end if

    ! Read the real*8 from the string buffer.
    read (ux,1030,err=995) var

    ! The format below should specify a real*8 field which matches
    ! the character length of the buffer variable (ichpar).
1030 format(g24.0)

    go to 999

995 continue
    j2 = ilnobl(ux)
    write (nttyo,1040) ux(1:j2)
1040 format(/" * Error - (EQLIBU/chreal) Can't read a real*8 number",/7x,'from the buffer string "',a,'" due to a formatted read',' error.')

    qrderr = .true.

999 continue
end subroutine chreal