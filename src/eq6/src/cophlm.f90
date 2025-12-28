subroutine cophlm(actw,awmax,awmin,delxi,dlxmin,eh,ehmax,ehmin,fo2lg,iodb,nodbmx,noutpt,o2max,o2min,ph,phmax,phmin,qadjdx,qredox,tolxsu)
    !! This subroutine checks for oversteps with regard to specified
    !! minimum and maximum values for the pH, Eh, log fO2, and
    !! activity of water.
    !! This subroutine is called by:
    !!   EQ6/path.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: nodbmx

    integer :: noutpt

    integer :: iodb(nodbmx)

    logical :: qadjdx
    logical :: qredox

    real(kind=8) :: actw
    real(kind=8) :: awmax
    real(kind=8) :: awmin
    real(kind=8) :: delxi
    real(kind=8) :: dlxmin
    real(kind=8) :: eh
    real(kind=8) :: ehmax
    real(kind=8) :: ehmin
    real(kind=8) :: fo2lg
    real(kind=8) :: o2max
    real(kind=8) :: o2min
    real(kind=8) :: ph
    real(kind=8) :: phmax
    real(kind=8) :: phmin
    real(kind=8) :: tolxsu

    ! Local variable declarations.
    ! None
    qadjdx = .false.

    ! Is the pH less than the requested minimum value?
    if ((ph - phmin) .lt. -tolxsu) then
        if (iodb(1) .gt. 0) then
            write (noutpt,1000) ph,phmin
1000 format(/3x,'The pH is ',f7.4,', which is less than the',' minimum value',/5x,'of ',f7.4,'. This overstep exceeds',' the specified tolerance.',/)
        end if

        ! Determine whether or not to reduce delxi and go back and try
        ! again to better locate the event in question.
        call dadjdx(delxi,dlxmin,iodb,nodbmx,noutpt,qadjdx)

        if (qadjdx) then
            go to 999
        end if
    end if

    ! Is the pH greater than the requested maximum value?
    if ((ph - phmax) .gt. tolxsu) then
        if (iodb(1) .gt. 0) then
            write (noutpt,1010) ph,phmax
1010 format(/3x,'The pH is ',f7.4,', which is more than the',' maximum value',/5x,'of ',f7.4,'. This overstep exceeds',' the specified tolerance.',/)
        end if

        ! Determine whether or not to reduce delxi and go back and try
        ! again to better locate the event in question.
        call dadjdx(delxi,dlxmin,iodb,nodbmx,noutpt,qadjdx)

        if (qadjdx) then
            go to 999
        end if
    end if

    if (qredox) then
        ! Is the Eh less than the requested minimum value?
        if ((eh - ehmin) .lt. -tolxsu) then
            if (iodb(1) .gt. 0) then
                write (noutpt,1020) eh,ehmin
1020 format(/3x,'The Eh is ',f9.4,' v, which is less than the',' minimum value',/5x,'of ',f9.4,'. This overstep exceeds',' the specified tolerance.',/)
            end if

            ! Determine whether or not to reduce delxi and go back and try
            ! again to better locate the event in question.
            call dadjdx(delxi,dlxmin,iodb,nodbmx,noutpt,qadjdx)

            if (qadjdx) then
                go to 999
            end if
        end if

        ! Is the Eh greater than the requested maximum value?
        if ((eh - ehmax) .gt. tolxsu) then
            if (iodb(1) .gt. 0) then
                write (noutpt,1040) eh,ehmax
1040 format(/3x,'The Eh is ',f9.4,' v, which is more than the',' maximum value',/5x,'of ',f9.4,'. This overstep exceeds',' the specified tolerance.',/)
            end if

            ! Determine whether or not to reduce delxi and go back and try
            ! again to better locate the event in question.
            call dadjdx(delxi,dlxmin,iodb,nodbmx,noutpt,qadjdx)

            if (qadjdx) then
                go to 999
            end if
        end if

        ! Is the log fO2 less than the requested minimum value?
        if ((fo2lg - o2min) .lt. -tolxsu) then
            if (iodb(1) .gt. 0) then
                write (noutpt,1050) fo2lg,o2min
1050 format(/3x,'The log fO2 is ',f9.4,', which is less than',' the minimum value',/5x,'of ',f9.4,'. This overstep',' exceeds the specified tolerance.',/)
            end if

            ! Determine whether or not to reduce delxi and go back and try
            ! again to better locate the event in question.
            call dadjdx(delxi,dlxmin,iodb,nodbmx,noutpt,qadjdx)

            if (qadjdx) then
                go to 999
            end if
        end if

        ! Is the Eh greater than the requested maximum value?
        if ((fo2lg - o2max) .gt. tolxsu) then
            if (iodb(1) .gt. 0) then
                write (noutpt,1060) fo2lg,o2max
1060 format(/3x,'The log fO2 is ',f9.4,', which is more than',' the maximum value',/5x,'of ',f9.4,'. This overstep',' exceeds the specified tolerance.',/)
            end if

            ! Determine whether or not to reduce delxi and go back and try
            ! again to better locate the event in question.
            call dadjdx(delxi,dlxmin,iodb,nodbmx,noutpt,qadjdx)

            if (qadjdx) then
                go to 999
            end if
        end if
    end if

    ! Is the activity of water less than the requested minimum value?
    if ((actw - awmin) .lt. -tolxsu) then
        if (iodb(1) .gt. 0) then
            write (noutpt,1070) actw,awmin
1070 format(/3x,'The activity of water is ',f6.4,', which is',' less than the',/5x,'minimum value of ',f6.4,'. This',' overstep exceeds the',/5x,'specified tolerance.',/)
        end if

        ! Determine whether or not to reduce delxi and go back and try
        ! again to better locate the event in question.
        call dadjdx(delxi,dlxmin,iodb,nodbmx,noutpt,qadjdx)

        if (qadjdx) then
            go to 999
        end if
    end if

    ! Is the activity of water greater than the requested maximum
    ! value?
    if ((actw - awmax) .gt. tolxsu) then
        if (iodb(1) .gt. 0) then
            write (noutpt,1080) actw,awmax
1080 format(/3x,'The activity of water is ',f6.4,', which is',' greater than the',/5x,'maximum value of ',f6.4,'. This',' overstep exceeds the',/5x,'specified tolerance.',/)
        end if

        ! Determine whether or not to reduce delxi and go back and try
        ! again to better locate the event in question.
        call dadjdx(delxi,dlxmin,iodb,nodbmx,noutpt,qadjdx)

        if (qadjdx) then
            go to 999
        end if
    end if

999 continue
end subroutine cophlm