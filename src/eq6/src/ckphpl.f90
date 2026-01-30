subroutine ckphpl(delxi,dlxmin,dph0,dxh0pl,dxh1pl,dxval0,eps100,iodb,nodbmx,nord,noutpt,nrd1mx,nttyo,ph0,ph1,ph0plo,ph1plo,prcinf,qdump,tolxsu,xi0,xi1,xval0)
    !! This subroutine checks to see that the next pH-based plot
    !! point is not exceeded. Because the pH might be decreasing or
    !! increasing, two potential target points (ph0plo and ph1plo)
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

    real(kind=8) :: dph0(nrd1mx)
    real(kind=8) :: dxval0(nrd1mx)

    real(kind=8) :: delxi
    real(kind=8) :: dlxmin
    real(kind=8) :: dxh0pl
    real(kind=8) :: dxh1pl
    real(kind=8) :: eps100
    real(kind=8) :: prcinf
    real(kind=8) :: ph0
    real(kind=8) :: ph0plo
    real(kind=8) :: ph1
    real(kind=8) :: ph1plo
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

    data usearch /'pH reaches a plot point value                   '/
    data unam24  /'                       '/

    dxsv = delxi

    dxh0pl = prcinf
    dxh1pl = prcinf

    if (nord .le. 0) then
        go to 999
    end if

    ! The first search target is the currently defined pH-based plot
    ! point associated with the lesser pH value.
    ilsign = +1
    xtargv = ph0plo
    tolsx = tolxsu

    icount = -1
100 continue
    xi1 = xi0 + delxi
    icount = icount + 1

    if (icount .gt. 10) then
        if (iodb(1).gt.0 .or. iodb(5).gt.0) then
            write (noutpt,1000)
1000 format(/3x,'A scan for where the pH reaches the lesser of',' two plot point values',/3x,'has failed. Dropping to the',' minimum step size.')
        end if

        delxi = dlxmin
        go to 110
    end if

    ! Estimate the pH from a Taylor's series expansion.
    ph1 = ph0
    dxp = 1.

    do n = 1,nord
        dxp = dxp*delxi
        ph1 = ph1 + ( dph0(n)/fctrl(n) )*dxp
    end do

    if (delxi .gt. dlxmin) then
        if (ph0 .le. (1. + tolxsu)*ph0plo) then
            ! The pH at the base point is too close to the requested
            ! lesser plot point value.
            delxi = dlxmin
            dxh0pl = delxi

            if (iodb(1).gt.0 .or. iodb(5).gt.0) then
                write (noutpt,1010) ph0,ph0plo,xi0,delxi
1010 format(/3x,'The base point pH of ',1pe11.4,' is already',' very close',/3x,'to the requested lesser plot point',' value of',e11.4,/3x,'at Xi= ',e11.4,'. Have set delxi',' equal to the minimum',/7x,'value of ',e11.4,'.')

                go to 100
            end if
        end if
    end if

    if (delxi .gt. dlxmin) then
        if ((ph1 - ph0plo) .lt. -tolsx) then
            xval0 = ph0

            do n = 1,nord
                dxval0(n) = dph0(n)
            end do

            call search(delxi,dlxmin,dxval0,eps100,ier,ilsign,iodb,nodbmx,nord,noutpt,nrd1mx,nttyo,tolsx,unam24,usearch,xtargv,xval0)

            dxh0pl = delxi

            if (ier .le. 0) then
                go to 100
            end if

            if (ier .ge. 2) then
                ! Note: if ier = 1, the returned "safe" value of delxi
                ! is used.
                delxi = dlxmin
            end if

            dxh0pl = delxi
            go to 110
        end if
    end if

    if (iodb(1).gt.0 .or. iodb(5).gt.0) then
        if (abs(ph1 - ph0plo) .le. tolsx) then
            write (noutpt,1020) ph0plo,xi1,delxi
1020 format(/3x,"Taylor's series predict that the pH is",' at the requested',/3x,'lesser plot point value of ',1pe11.4,' at Xi= ',e11.4,',',/3x,'delxi= ',e11.4,'.')
        end if
    end if

110 continue

    ! The second search target is the currently defined pH-based plot
    ! point associated with the greater pH value.
    ilsign = -1
    xtargv = ph1plo
    tolsx = tolxsu

    icount = -1
200 continue
    xi1 = xi0 + delxi
    icount = icount + 1

    if (icount .gt. 10) then
        if (iodb(1).gt.0 .or. iodb(5).gt.0) then
            write (noutpt,1100)
1100 format(/3x,'A scan for where the pH reaches the greater of',' two plot point values',/3x,'has failed. Dropping to the',' minimum step size.')
        end if

        delxi = dlxmin
        go to 210
    end if

    ! Estimate the pH from a Taylor's series expansion.
    ph1 = ph0
    dxp = 1.

    do n = 1,nord
        dxp = dxp*delxi
        ph1 = ph1 + ( dph0(n)/fctrl(n) )*dxp
    end do

    if (delxi .gt. dlxmin) then
        if (ph0 .ge. (1. - tolxsu)*ph1plo) then
            ! The pH at the base point is too close to the requested
            ! greater plot point value.
            delxi = dlxmin
            dxh1pl = delxi

            if (iodb(1).gt.0 .or. iodb(5).gt.0) then
                write (noutpt,1110) ph0,ph1plo,xi0,delxi
1110 format(/3x,'The base point pH of ',1pe11.4,' is already',' very close',/3x,'to the requested greater plot point',' value of',e11.4,/3x,'at Xi= ',e11.4,'. Have set delxi',' equal to the minimum',/7x,'value of ',e11.4,'.')

                go to 200
            end if
        end if
    end if

    if (delxi .gt. dlxmin) then
        if ((ph1 - ph1plo) .gt. tolsx) then
            xval0 = ph0

            do n = 1,nord
                dxval0(n) = dph0(n)
            end do

            call search(delxi,dlxmin,dxval0,eps100,ier,ilsign,iodb,nodbmx,nord,noutpt,nrd1mx,nttyo,tolsx,unam24,usearch,xtargv,xval0)

            dxh1pl = delxi

            if (ier .le. 0) then
                go to 200
            end if

            if (ier .ge. 2) then
                ! Note: if ier = 1, the returned "safe" value of delxi
                ! is used.
                delxi = dlxmin
            end if

            dxh1pl = delxi
            go to 210
        end if
    end if

    if (iodb(1).gt.0 .or. iodb(5).gt.0) then
        if (abs(ph1 - ph1plo) .le. tolsx) then
            write (noutpt,1120) ph1plo,xi1,delxi
1120 format(/3x,"Taylor's series predict that the pH is",' at the requested',/3x,'greater plot point value of ',1pe11.4,' at Xi= ',e11.4,',',/3x,'delxi= ',e11.4,'.')
        end if
    end if

210 continue

    if (dxsv .gt. delxi) then
        qdump = .false.
    end if

999 continue
end subroutine ckphpl
