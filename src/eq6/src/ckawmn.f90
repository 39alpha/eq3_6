subroutine ckawmn(delxi,dlxmin,daw0,dxw0mx,dxval0,eps100,iodb,nodbmx,nord,noutpt,nrd1mx,nttyo,aw0,aw1,awmin,prcinf,qdump,tolxsu,xi0,xi1,xval0)
    !! This subroutine checks to see that the requested minimum value of
    !! activity of water (aw) is not exceeded.
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

    real(kind=8) :: daw0(nrd1mx)
    real(kind=8) :: dxval0(nrd1mx)

    real(kind=8) :: delxi
    real(kind=8) :: dlxmin
    real(kind=8) :: dxw0mx
    real(kind=8) :: eps100
    real(kind=8) :: prcinf
    real(kind=8) :: awmin
    real(kind=8) :: aw0
    real(kind=8) :: aw1
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

    data usearch /'aw reaches the requested minimum value          '/
    data unam24  /'                       '/

    dxsv = delxi

    dxw0mx = prcinf

    if (nord .le. 0) then
        go to 999
    end if

    ! The search target is the requested limit on the aw.
    ilsign = +1
    xtargv = awmin
    tolsx = tolxsu

    icount = -1
100 continue
    xi1 = xi0 + delxi
    icount = icount + 1

    if (icount .gt. 10) then
        if (iodb(1).gt.0 .or. iodb(5).gt.0) then
            write (noutpt,1000)
1000 format(/3x,'A scan for where the aw reaches the requested',' minimum value',/3x,'has failed. Dropping to the minimum',' step size.')
        end if

        delxi = dlxmin
        go to 990
    end if

    ! Estimate the aw from a Taylor's series expansion.
    aw1 = aw0
    dxp = 1.

    do n = 1,nord
        dxp = dxp*delxi
        aw1 = aw1 + ( daw0(n)/fctrl(n) )*dxp
    end do

    if (delxi .gt. dlxmin) then
        if (aw0 .le. (1. + tolxsu)*awmin) then
            ! The aw at the base point is too close to the requested
            ! minimum value.
            delxi = dlxmin
            dxw0mx = delxi

            if (iodb(1).gt.0 .or. iodb(5).gt.0) then
                write (noutpt,1010) aw0,awmin,xi0,delxi
1010 format(/3x,'The base point aw of ',1pe11.4,' is already',' very close',/3x,'to the requested value of ',e11.4,/3x,'at Xi= ',e11.4,'. Have set delxi',' equal to the minimum',/7x,'value of ',e11.4,'.')

                go to 100
            end if
        end if
    end if

    if (delxi .gt. dlxmin) then
        if ((aw1 - awmin) .lt. -tolsx) then
            xval0 = aw0

            do n = 1,nord
                dxval0(n) = daw0(n)
            end do

            call search(delxi,dlxmin,dxval0,eps100,ier,ilsign,iodb,nodbmx,nord,noutpt,nrd1mx,nttyo,tolsx,unam24,usearch,xtargv,xval0)

            dxw0mx = delxi

            if (ier .le. 0) then
                go to 100
            end if

            if (ier .ge. 2) then
                ! Note: if ier = 1, the returned "safe" value of delxi
                ! is used.
                delxi = dlxmin
            end if

            dxw0mx = delxi
            go to 990
        end if
    end if

    if (iodb(1).gt.0 .or. iodb(5).gt.0) then
        if (abs(aw1 - awmin) .le. tolsx) then
            write (noutpt,1020) awmin,xi1,delxi
1020 format(/3x,"Taylor's series predict that the aw is",' at the requested',/3x,'maximum value of ',1pe11.4,' at Xi= ',e11.4,',',/3x,'delxi= ',e11.4,'.')
        end if
    end if

990 continue
    if (dxsv .gt. delxi) then
        qdump = .false.
    end if

999 continue
end subroutine ckawmn