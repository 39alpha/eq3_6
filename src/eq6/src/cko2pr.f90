subroutine cko2pr(delxi,dlxmin,do20,dxo0pr,dxo1pr,dxval0,eps100,iodb,nodbmx,nord,noutpt,nrd1mx,nttyo,fo2lg0,fo2lg1,o20prn,o21prn,prcinf,qdump,tolxsu,xi0,xi1,xval0)
    !! This subroutine checks to see that the next log fO2-based print
    !! point is not exceeded. Because the log fO2 might be decreasing
    !! or increasing, two potential target points (o20prn and o21prn)
    !! must be addressed.
    !! This subroutine is called by:
    !!   EQ6/path.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: nodbmx
    integer :: nrd1mx

    integer :: noutpt
    integer :: nttyo

    integer :: iodb(nodbmx)

    integer :: nord

    logical :: qdump

    real(kind=8) :: do20(nrd1mx)
    real(kind=8) :: dxval0(nrd1mx)

    real(kind=8) :: delxi
    real(kind=8) :: dlxmin
    real(kind=8) :: dxo0pr
    real(kind=8) :: dxo1pr
    real(kind=8) :: eps100
    real(kind=8) :: prcinf
    real(kind=8) :: fo2lg0
    real(kind=8) :: o20prn
    real(kind=8) :: fo2lg1
    real(kind=8) :: o21prn
    real(kind=8) :: tolxsu
    real(kind=8) :: xi0
    real(kind=8) :: xi1
    real(kind=8) :: xval0
    real(kind=8) :: xtargv

    ! Local variable declarations.
    integer :: icount
    integer :: ier
    integer :: ilsign
    integer :: n

    character(len=48) :: usearch
    character(len=24) :: unam24

    real(kind=8) :: dxp
    real(kind=8) :: dxsv
    real(kind=8) :: tolsx

    real(kind=8) :: fctrl

    data usearch /'log fO2 reaches a print point value             '/
    data unam24  /'                       '/

    dxsv = delxi

    dxo0pr = prcinf
    dxo1pr = prcinf

    if (nord .le. 0) then
        go to 999
    end if

    ! The first search target is the currently defined log fO2-based
    ! print point associated with the lesser log fO2 value.
    ilsign = +1
    xtargv = o20prn
    tolsx = tolxsu

    icount = -1
100 continue
    xi1 = xi0 + delxi
    icount = icount + 1

    if (icount .gt. 10) then
        if (iodb(1).gt.0 .or. iodb(5).gt.0) then
            write (noutpt,1000)
1000 format(/3x,'A scan for where the log fO2 reaches the lesser',' of two print point values',/3x,'has failed. Dropping to',' the minimum step size.')
        end if

        delxi = dlxmin
        go to 110
    end if

    ! Estimate the log fO2 from a Taylor's series expansion.
    fo2lg1 = fo2lg0
    dxp = 1.

    do n = 1,nord
        dxp = dxp*delxi
        fo2lg1 = fo2lg1 + ( do20(n)/fctrl(n) )*dxp
    end do

    if (delxi .gt. dlxmin) then
        if (fo2lg0 .le. (1. + tolxsu)*o20prn) then
            ! The log fO2 at the base point is too close to the requested
            ! lesser print point value.
            delxi = dlxmin
            dxo0pr = delxi

            if (iodb(1).gt.0 .or. iodb(5).gt.0) then
                write (noutpt,1010) fo2lg0,o20prn,xi0,delxi
1010 format(/3x,'The base point log fO2 of ',1pe11.4,' is',' already very close',/3x,'to the requested lesser print',' point value of',e11.4,/3x,'at Xi= ',e11.4,'. Have set',' delxi equal to the minimum',/7x,'value of ',e11.4,'.')

                go to 100
            end if
        end if
    end if

    if (delxi .gt. dlxmin) then
        if ((fo2lg1 - o20prn) .lt. -tolsx) then
            xval0 = fo2lg0

            do n = 1,nord
                dxval0(n) = do20(n)
            end do

            call search(delxi,dlxmin,dxval0,eps100,ier,ilsign,iodb,nodbmx,nord,noutpt,nrd1mx,nttyo,tolsx,unam24,usearch,xtargv,xval0)

            dxo0pr = delxi

            if (ier .le. 0) then
                go to 100
            end if

            if (ier .ge. 2) then
                ! Note: if ier = 1, the returned "safe" value of delxi
                ! is used.
                delxi = dlxmin
            end if

            dxo0pr = delxi
            go to 110
        end if
    end if

    if (iodb(1).gt.0 .or. iodb(5).gt.0) then
        if (abs(fo2lg1 - o20prn) .le. tolsx) then
            write (noutpt,1020) o20prn,xi1,delxi
1020 format(/3x,"Taylor's series predict that the log fO2 is",' at the requested',/3x,'lesser print point value of ',1pe11.4,' at Xi= ',e11.4,',',/3x,'delxi= ',e11.4,'.')
        end if
    end if

110 continue

    ! The second search target is the currently defined log fO2-based
    ! print point associated with the greater log fO2 value.
    ilsign = -1
    xtargv = o21prn
    tolsx = tolxsu

    icount = -1
200 continue
    xi1 = xi0 + delxi
    icount = icount + 1

    if (icount .gt. 10) then
        if (iodb(1).gt.0 .or. iodb(5).gt.0) then
            write (noutpt,1100)
1100 format(/3x,'A scan for where the log fO2 reaches the greater',' of two print point values',/3x,'has failed. Dropping to',' the minimum step size.')
        end if

        delxi = dlxmin
        go to 210
    end if

    ! Estimate the log fO2 from a Taylor's series expansion.
    fo2lg1 = fo2lg0
    dxp = 1.

    do n = 1,nord
        dxp = dxp*delxi
        fo2lg1 = fo2lg1 + ( do20(n)/fctrl(n) )*dxp
    end do

    if (delxi .gt. dlxmin) then
        if (fo2lg0 .ge. (1. - tolxsu)*o21prn) then
            ! The log fO2 at the base point is too close to the requested
            ! greater print point value.
            delxi = dlxmin
            dxo1pr = delxi

            if (iodb(1).gt.0 .or. iodb(5).gt.0) then
                write (noutpt,1110) fo2lg0,o21prn,xi0,delxi
1110 format(/3x,'The base point log fO2 of ',1pe11.4,' is',' already very close',/3x,'to the requested greater print',' point value of',e11.4,/3x,'at Xi= ',e11.4,'. Have set',' delxi equal to the minimum',/7x,'value of ',e11.4,'.')

                go to 200
            end if
        end if
    end if

    if (delxi .gt. dlxmin) then
        if ((fo2lg1 - o21prn) .gt. tolsx) then
            xval0 = fo2lg0

            do n = 1,nord
                dxval0(n) = do20(n)
            end do

            call search(delxi,dlxmin,dxval0,eps100,ier,ilsign,iodb,nodbmx,nord,noutpt,nrd1mx,nttyo,tolsx,unam24,usearch,xtargv,xval0)

            dxo1pr = delxi

            if (ier .le. 0) then
                go to 200
            end if

            if (ier .ge. 2) then
                ! Note: if ier = 1, the returned "safe" value of delxi
                ! is used.
                delxi = dlxmin
            end if

            dxo1pr = delxi
            go to 210
        end if
    end if

    if (iodb(1).gt.0 .or. iodb(5).gt.0) then
        if (abs(fo2lg1 - o21prn) .le. tolsx) then
            write (noutpt,1120) o21prn,xi1,delxi
1120 format(/3x,"Taylor's series predict that the log fO2 is",' at the requested',/3x,'greater print point value of ',1pe11.4,' at Xi= ',e11.4,',',/3x,'delxi= ',e11.4,'.')
        end if
    end if

210 continue

    if (dxsv .gt. delxi) then
        qdump = .false.
    end if

999 continue
end subroutine cko2pr
