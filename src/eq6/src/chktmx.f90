subroutine chktmx(delxi,dlxmin,dlxtmx,drir0,dxval0,eps100,iodb,nodbmx,nord,noutpt,nrd1mx,nttyo,prcinf,qdump,qriinf,rirec0,timemx,time0,time1,tolxst,xi0,xi1,xval0)
    !! This subroutine checks to see that the requested maximum value of
    !! time is not exceeded.
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
    logical :: qriinf

    real(kind=8) :: drir0(nrd1mx)
    real(kind=8) :: dxval0(nrd1mx)

    real(kind=8) :: delxi
    real(kind=8) :: dlxmin
    real(kind=8) :: dlxtmx
    real(kind=8) :: eps100
    real(kind=8) :: prcinf
    real(kind=8) :: rirec0
    real(kind=8) :: timemx
    real(kind=8) :: time0
    real(kind=8) :: time1
    real(kind=8) :: tolxst
    real(kind=8) :: xi0
    real(kind=8) :: xi1
    real(kind=8) :: xval0
    real(kind=8) :: xtargv

    ! Local variable declarations.
    integer :: icount
    integer :: ier
    integer :: ilsign
    integer :: n
    integer :: nordp1

    character(len=48) :: usearch
    character(len=24) :: unam24

    real(kind=8) :: dxp
    real(kind=8) :: dxsv
    real(kind=8) :: tolsx

    real(kind=8) :: fctrl

    data usearch /'time reaches the requested maximum value        '/
    data unam24  /'                       '/

    dxsv = delxi

    nordp1 = nord + 1
    dlxtmx = prcinf

    if (timemx .le. 0.) then
        go to 999
    end if

    ! The search target is the requested limit on the time. The interval
    ! for convergence is ((1. - tolxst)*timemx,(1. + tolxst)*timemx)).
    ilsign = -1
    xtargv = timemx
    tolsx = tolxst*timemx

    icount = -1
100 continue
    xi1 = xi0 + delxi
    icount = icount + 1

    if (icount .gt. 10) then
        if (iodb(1).gt.0 .or. iodb(5).gt.0) then
            write (noutpt,1000)
1000 format(/3x,'A scan for where time reaches the requested',' maximum value',/3x,'has failed. Dropping to the minimum',' step size.')
        end if

        delxi = dlxmin
        go to 990
    end if

    ! Estimate the time from a Taylor's series expansion.
    if (qriinf) then
        time1 = prcinf
    else
        time1 = time0 + rirec0*delxi
        dxp = delxi

        do n = 1,nord
            dxp = dxp*delxi
            time1 = time1 + ( drir0(n)*dxp )/fctrl(n + 1)
        end do
    end if

    if (delxi .gt. dlxmin) then
        if (time0 .gt. (1. - tolxst)*timemx) then
            ! The time at the base point is too close to the requested
            ! maximum value.
            delxi = dlxmin
            dlxtmx = delxi

            if (iodb(1).gt.0 .or. iodb(5).gt.0) then
                write (noutpt,1010) time0,timemx,xi0,delxi
1010 format(/3x,'The base point time of ',1pe11.4,' seconds',' is already very close',/3x,'to the requested value of ',e11.4,' seconds',/3x,'at Xi= ',e11.4,'. Have set delxi',' equal to the minimum',/7x,'value of ',e11.4,'.')

                go to 100
            end if
        end if
    end if

    if (delxi .gt. dlxmin) then
        if ((time1 - timemx) .gt. tolsx) then
            xval0 = time0
            dxval0(1) = rirec0

            do n = 2,nordp1
                dxval0(n) = drir0(n - 1)
            end do

            ! Calling sequence substitutions:
            !   nordp1 for nord
            call search(delxi,dlxmin,dxval0,eps100,ier,ilsign,iodb,nodbmx,nordp1,noutpt,nrd1mx,nttyo,tolsx,unam24,usearch,xtargv,xval0)

            dlxtmx = delxi

            if (ier .le. 0) then
                go to 100
            end if

            if (ier .ge. 2) then
                ! Note: if ier = 1, the returned "safe" value of delxi
                ! is used.
                delxi = dlxmin
            end if

            dlxtmx = delxi
            go to 990
        end if
    end if

    if (iodb(1).gt.0 .or. iodb(5).gt.0) then
        if (abs(time1 - timemx) .le. tolsx) then
            write (noutpt,1020) timemx,xi1,delxi
1020 format(/3x,"Taylor's series predict that the model time is",' at the requested',/3x,'maximum value of ',1pe11.4,' seconds at Xi= ',e11.4,',',/3x,'delxi= ',e11.4,'.')
        end if
    end if

990 continue
    if (dxsv .gt. delxi) then
        qdump = .false.
    end if

999 continue
end subroutine chktmx