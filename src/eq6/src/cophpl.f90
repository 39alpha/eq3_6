subroutine cophpl(actw,aw0plo,aw1plo,delxi,dlxmin,eh,eh0plo,eh1plo,fo2lg,iodb,nodbmx,noutpt,o20plo,o21plo,ph,ph0plo,ph1plo,qadjdx,qredox,tolxsu)
    !! This subroutine checks for oversteps with regard to currently
    !! defined lesser and greater plot point values for the pH, Eh,
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
    real(kind=8) :: aw0plo
    real(kind=8) :: aw1plo
    real(kind=8) :: delxi
    real(kind=8) :: dlxmin
    real(kind=8) :: eh
    real(kind=8) :: eh0plo
    real(kind=8) :: eh1plo
    real(kind=8) :: fo2lg
    real(kind=8) :: o20plo
    real(kind=8) :: o21plo
    real(kind=8) :: ph
    real(kind=8) :: ph0plo
    real(kind=8) :: ph1plo
    real(kind=8) :: tolxsu

    ! Local variable declarations.
    ! None
    qadjdx = .false.

    ! Is the pH less than the current lesser plot point value?
    if ((ph - ph0plo) .lt. -tolxsu) then
        if (iodb(1) .gt. 0) then
            write (noutpt,1000) ph,ph0plo
1000 format(/3x,'The pH is ',f7.4,', which is less than the',' current lesser plot point',/5x,'value of ',f7.4,'. This',' overstep exceeds the specified tolerance.',/)
        end if

        ! Determine whether or not to reduce delxi and go back and try
        ! again to better locate the event in question.
        call dadjdx(delxi,dlxmin,iodb,nodbmx,noutpt,qadjdx)

        if (qadjdx) then
            go to 999
        end if
    end if

    ! Is the pH greater than the current greater plot point value?
    if ((ph - ph1plo) .gt. tolxsu) then
        if (iodb(1) .gt. 0) then
            write (noutpt,1010) ph,ph1plo
1010 format(/3x,'The pH is ',f7.4,', which is greater than the',' current greater plot point',/5x,'value of ',f7.4,'. This',' overstep exceeds the specified tolerance.',/)
        end if

        ! Determine whether or not to reduce delxi and go back and try
        ! again to better locate the event in question.
        call dadjdx(delxi,dlxmin,iodb,nodbmx,noutpt,qadjdx)

        if (qadjdx) then
            go to 999
        end if
    end if

    if (qredox) then
        ! Is the Eh less than the current lesser plot point value?
        if ((eh - eh0plo) .lt. -tolxsu) then
            if (iodb(1) .gt. 0) then
                write (noutpt,1020) eh,eh0plo
1020 format(/3x,'The Eh is ',f7.4,' v, which is less than the',' current lesser plot point',/5x,'value of ',f7.4,' v.',' This overstep exceeds the specified tolerance.',/)
            end if

            ! Determine whether or not to reduce delxi and go back and try
            ! again to better locate the event in question.
            call dadjdx(delxi,dlxmin,iodb,nodbmx,noutpt,qadjdx)

            if (qadjdx) then
                go to 999
            end if
        end if

        ! Is the Eh greater than the current greater plot point value?
        if ((eh - eh1plo) .gt. tolxsu) then
            if (iodb(1) .gt. 0) then
                write (noutpt,1030) eh,eh1plo
1030 format(/3x,'The Eh is ',f7.4,' v, which is greater than',' the current greater plot point',/5x,'value of ',f7.4,' v. This overstep exceeds the specified tolerance.',/)
            end if

            ! Determine whether or not to reduce delxi and go back and try
            ! again to better locate the event in question.
            call dadjdx(delxi,dlxmin,iodb,nodbmx,noutpt,qadjdx)

            if (qadjdx) then
                go to 999
            end if
        end if

        ! Is the log fO2 less than the current lesser plot point value?
        if ((fo2lg - o20plo) .lt. -tolxsu) then
            if (iodb(1) .gt. 0) then
                write (noutpt,1040) fo2lg,o20plo
1040 format(/3x,'The log fO2 is ',f7.4,', which is less than',' the current lesser plot point',/5x,'value of ',f7.4,'. This overstep exceeds the specified tolerance.',/)
            end if

            ! Determine whether or not to reduce delxi and go back and try
            ! again to better locate the event in question.
            call dadjdx(delxi,dlxmin,iodb,nodbmx,noutpt,qadjdx)

            if (qadjdx) then
                go to 999
            end if
        end if

        ! Is the log fO2 greater than the current greater plot point
        ! value?
        if ((fo2lg - o21plo) .gt. tolxsu) then
            if (iodb(1) .gt. 0) then
                write (noutpt,1050) fo2lg,o21plo
1050 format(/3x,'The log fO2 is ',f7.4,', which is greater',' than the current greater plot',/5x,'point value of ',f7.4,'. This overstep exceeds the specified tolerance.',/)
            end if

            ! Determine whether or not to reduce delxi and go back and try
            ! again to better locate the event in question.
            call dadjdx(delxi,dlxmin,iodb,nodbmx,noutpt,qadjdx)

            if (qadjdx) then
                go to 999
            end if
        end if
    end if

    ! Is the activity of water less than the current lesser plot
    ! point value?
    if ((actw - aw0plo) .lt. -tolxsu) then
        if (iodb(1) .gt. 0) then
            write (noutpt,1060) actw,aw0plo
1060 format(/3x,'The activity of water is ',f7.4,', which is less',' than the current lesser',/5x,'plot point value of ',f7.4,'. This overstep exceeds the specified',/5x,'tolerance.',/)
        end if

        ! Determine whether or not to reduce delxi and go back and try
        ! again to better locate the event in question.
        call dadjdx(delxi,dlxmin,iodb,nodbmx,noutpt,qadjdx)

        if (qadjdx) then
            go to 999
        end if
    end if

    ! Is the activity of water greater than the current greater plot
    ! point value?
    if ((actw - aw1plo) .gt. tolxsu) then
        if (iodb(1) .gt. 0) then
            write (noutpt,1070) actw,aw1plo
1070 format(/3x,'The activity of water is ',f7.4,', which is',' greater than the current',/5x,'greater plot point value',' of ',f7.4,'. This overstep exceeds the',/5x,'specified',' tolerance.',/)
        end if

        ! Determine whether or not to reduce delxi and go back and try
        ! again to better locate the event in question.
        call dadjdx(delxi,dlxmin,iodb,nodbmx,noutpt,qadjdx)

        if (qadjdx) then
            go to 999
        end if
    end if

999 continue
end subroutine cophpl