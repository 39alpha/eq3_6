subroutine cophpr(actw,aw0prn,aw1prn,delxi,dlxmin,eh,eh0prn,eh1prn,fo2lg,iodb,nodbmx,noutpt,o20prn,o21prn,ph,ph0prn,ph1prn,qadjdx,qredox,tolxsu)
    !! This subroutine checks for oversteps with regard to currently
    !! defined lesser and greater print point values for the pH, Eh,
    !! log fO2, and activity of water.
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
    real(kind=8) :: aw0prn
    real(kind=8) :: aw1prn
    real(kind=8) :: delxi
    real(kind=8) :: dlxmin
    real(kind=8) :: eh
    real(kind=8) :: eh0prn
    real(kind=8) :: eh1prn
    real(kind=8) :: fo2lg
    real(kind=8) :: o20prn
    real(kind=8) :: o21prn
    real(kind=8) :: ph
    real(kind=8) :: ph0prn
    real(kind=8) :: ph1prn
    real(kind=8) :: tolxsu

    ! Local variable declarations.
    ! None
    qadjdx = .false.

    ! Is the pH less than the current lesser print point value?
    if ((ph - ph0prn) .lt. -tolxsu) then
        if (iodb(1) .gt. 0) then
            write (noutpt,1000) ph,ph0prn
1000 format(/3x,'The pH is ',f7.4,', which is less than the',' current lesser print point',/5x,'value of ',f7.4,'. This',' overstep exceeds the specified tolerance.',/)
        end if

        ! Determine whether or not to reduce delxi and go back and try
        ! again to better locate the event in question.
        call dadjdx(delxi,dlxmin,iodb,nodbmx,noutpt,qadjdx)

        if (qadjdx) then
            go to 999
        end if
    end if

    ! Is the pH greater than the current greater print point value?
    if ((ph - ph1prn) .gt. tolxsu) then
        if (iodb(1) .gt. 0) then
            write (noutpt,1010) ph,ph1prn
1010 format(/3x,'The pH is ',f7.4,', which is greater than the',' current greater print point',/5x,'value of ',f7.4,'. This',' overstep exceeds the specified tolerance.',/)
        end if

        ! Determine whether or not to reduce delxi and go back and try
        ! again to better locate the event in question.
        call dadjdx(delxi,dlxmin,iodb,nodbmx,noutpt,qadjdx)

        if (qadjdx) then
            go to 999
        end if
    end if

    if (qredox) then
        ! Is the Eh less than the current lesser print point value?
        if ((eh - eh0prn) .lt. -tolxsu) then
            if (iodb(1) .gt. 0) then
                write (noutpt,1020) eh,eh0prn
1020 format(/3x,'The Eh is ',f7.4,' v, which is less than the',' current lesser print point',/5x,'value of ',f7.4,' v.',' This overstep exceeds the specified tolerance.',/)
            end if

            ! Determine whether or not to reduce delxi and go back and try
            ! again to better locate the event in question.
            call dadjdx(delxi,dlxmin,iodb,nodbmx,noutpt,qadjdx)

            if (qadjdx) then
                go to 999
            end if
        end if

        ! Is the Eh greater than the current greater print point value?
        if ((eh - eh1prn) .gt. tolxsu) then
            if (iodb(1) .gt. 0) then
                write (noutpt,1030) eh,eh1prn
1030 format(/3x,'The Eh is ',f7.4,' v, which is greater than',' the current greater print point',/5x,'value of ',f7.4,' v. This overstep exceeds the specified tolerance.',/)
            end if

            ! Determine whether or not to reduce delxi and go back and try
            ! again to better locate the event in question.
            call dadjdx(delxi,dlxmin,iodb,nodbmx,noutpt,qadjdx)

            if (qadjdx) then
                go to 999
            end if
        end if

        ! Is the log fO2 less than the current lesser print point value?
        if ((fo2lg - o20prn) .lt. -tolxsu) then
            if (iodb(1) .gt. 0) then
                write (noutpt,1040) fo2lg,o20prn
1040 format(/3x,'The log fO2 is ',f7.4,', which is less than',' the current lesser print point',/5x,'value of ',f7.4,'. This overstep exceeds the specified tolerance.',/)
            end if

            ! Determine whether or not to reduce delxi and go back and try
            ! again to better locate the event in question.
            call dadjdx(delxi,dlxmin,iodb,nodbmx,noutpt,qadjdx)

            if (qadjdx) then
                go to 999
            end if
        end if

        ! Is the log fO2 greater than the current greater print point
        ! value?
        if ((fo2lg - o21prn) .gt. tolxsu) then
            if (iodb(1) .gt. 0) then
                write (noutpt,1050) fo2lg,o21prn
1050 format(/3x,'The log fO2 is ',f7.4,', which is greater',' than the current greater print',/5x,'point value of ',f7.4,'. This overstep exceeds the specified tolerance.',/)
            end if

            ! Determine whether or not to reduce delxi and go back and try
            ! again to better locate the event in question.
            call dadjdx(delxi,dlxmin,iodb,nodbmx,noutpt,qadjdx)

            if (qadjdx) then
                go to 999
            end if
        end if
    end if

    ! Is the activity of water less than the current lesser print
    ! point value?
    if ((actw - aw0prn) .lt. -tolxsu) then
        if (iodb(1) .gt. 0) then
            write (noutpt,1060) actw,aw0prn
1060 format(/3x,'The activity of water is ',f7.4,', which is less',' than the current lesser',/5x,'print point value of ',f7.4,'. This overstep exceeds the specified',/5x,'tolerance.',/)
        end if

        ! Determine whether or not to reduce delxi and go back and try
        ! again to better locate the event in question.
        call dadjdx(delxi,dlxmin,iodb,nodbmx,noutpt,qadjdx)

        if (qadjdx) then
            go to 999
        end if
    end if

    ! Is the activity of water greater than the current greater print
    ! point value?
    if ((actw - aw1prn) .gt. tolxsu) then
        if (iodb(1) .gt. 0) then
            write (noutpt,1070) actw,aw1prn
1070 format(/3x,'The activity of water is ',f7.4,', which is',' greater than the current',/5x,'greater print point value',' of ',f7.4,'. This overstep exceeds the',/5x,'specified',' tolerance.',/)
        end if

        ! Determine whether or not to reduce delxi and go back and try
        ! again to better locate the event in question.
        call dadjdx(delxi,dlxmin,iodb,nodbmx,noutpt,qadjdx)

        if (qadjdx) then
            go to 999
        end if
    end if

999 continue
end subroutine cophpr