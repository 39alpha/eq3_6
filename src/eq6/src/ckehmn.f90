subroutine ckehmn(delxi,dlxmin,deh0,dxe0mx,dxval0,eps100,iodb,nodbmx,nord,noutpt,nrd1mx,nttyo,eh0,eh1,ehmin,prcinf,qdump,tolxsu,xi0,xi1,xval0)
    !! This subroutine checks to see that the requested minimum value of
    !! Eh is not exceeded.
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

    real(kind=8) :: deh0(nrd1mx)
    real(kind=8) :: dxval0(nrd1mx)

    real(kind=8) :: delxi
    real(kind=8) :: dlxmin
    real(kind=8) :: dxe0mx
    real(kind=8) :: eps100
    real(kind=8) :: prcinf
    real(kind=8) :: ehmin
    real(kind=8) :: eh0
    real(kind=8) :: eh1
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

    data usearch /'Eh reaches the requested minimum value          '/
    data unam24  /'                       '/

    dxsv = delxi

    dxe0mx = prcinf

    if (nord .le. 0) then
        go to 999
    end if

    ! The search target is the requested limit on the Eh.
    ilsign = +1
    xtargv = ehmin
    tolsx = tolxsu

    icount = -1
100 continue
    xi1 = xi0 + delxi
    icount = icount + 1

    if (icount .gt. 10) then
        if (iodb(1).gt.0 .or. iodb(5).gt.0) then
            write (noutpt,1000)
1000 format(/3x,'A scan for where the Eh reaches the requested',' minimum value',/3x,'has failed. Dropping to the minimum',' step size.')
        end if

        delxi = dlxmin
        go to 990
    end if

    ! Estimate the Eh from a Taylor's series expansion.
    eh1 = eh0
    dxp = 1.

    do n = 1,nord
        dxp = dxp*delxi
        eh1 = eh1 + ( deh0(n)/fctrl(n) )*dxp
    end do

    if (delxi .gt. dlxmin) then
        if (eh0 .le. (1. + tolxsu)*ehmin) then
            ! The Eh at the base point is too close to the requested
            ! minimum value.
            delxi = dlxmin
            dxe0mx = delxi

            if (iodb(1).gt.0 .or. iodb(5).gt.0) then
                write (noutpt,1010) eh0,ehmin,xi0,delxi
1010 format(/3x,'The base point Eh of ',1pe11.4,' v is already',' very close',/3x,'to the requested value of ',e11.4,/3x,'at Xi= ',e11.4,'. Have set delxi',' equal to the minimum',/7x,'value of ',e11.4,'.')

                go to 100
            end if
        end if
    end if

    if (delxi .gt. dlxmin) then
        if ((eh1 - ehmin) .lt. -tolsx) then
            xval0 = eh0

            do n = 1,nord
                dxval0(n) = deh0(n)
            end do

            call search(delxi,dlxmin,dxval0,eps100,ier,ilsign,iodb,nodbmx,nord,noutpt,nrd1mx,nttyo,tolsx,unam24,usearch,xtargv,xval0)

            dxe0mx = delxi

            if (ier .le. 0) then
                go to 100
            end if

            if (ier .ge. 2) then
                ! Note: if ier = 1, the returned "safe" value of delxi
                ! is used.
                delxi = dlxmin
            end if

            dxe0mx = delxi
            go to 990
        end if
    end if

    if (iodb(1).gt.0 .or. iodb(5).gt.0) then
        if (abs(eh1 - ehmin) .le. tolsx) then
            write (noutpt,1020) ehmin,xi1,delxi
1020 format(/3x,"Taylor's series predict that the Eh is",' at the requested',/3x,'minimum value of ',1pe11.4,' v at Xi= ',e11.4,',',/3x,'delxi= ',e11.4,'.')
        end if
    end if

990 continue
    if (dxsv .gt. delxi) then
        qdump = .false.
    end if

999 continue
end subroutine ckehmn