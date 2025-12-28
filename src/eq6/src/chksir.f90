subroutine chksir(delxi,dlxmin,drir0,dxval0,eps100,iodb,nodbmx,nord,noutpt,nrd1mx,nttyo,rirec0,rirecp,xi0,xi1,xval0)
    !! This subroutine checks the sign of the inverse rate. It finds the
    !! point of reaction progress at which the inverse rate becomes
    !! zero. Technically, this should never happen. The inverse rate
    !! should start with a positive value and increase toward infinity.
    !! The inverse rate is tracked using finite differences.
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

    real(kind=8) :: drir0(nrd1mx)
    real(kind=8) :: dxval0(nrd1mx)

    real(kind=8) :: delxi
    real(kind=8) :: dlxmin
    real(kind=8) :: eps100
    real(kind=8) :: rirec0
    real(kind=8) :: rirecp
    real(kind=8) :: xi0
    real(kind=8) :: xi1
    real(kind=8) :: xval0

    ! Local variable declarations.
    integer :: icount
    integer :: ier
    integer :: ilsign
    integer :: n

    character(len=48) :: usearch
    character(len=24) :: unam24

    real(kind=8) :: dxp
    real(kind=8) :: tolsx
    real(kind=8) :: xtargv

    real(kind=8) :: fctrl

    data usearch /'the inverse rate becomes zero                   '/
    data unam24 /'                        '/

    if (nord .le. 0) then
        go to 999
    end if

    ! Note that the search target is not where the inverse rate is zero,
    ! but eps100. The interval for convergence is (0,2*eps100).
    ilsign = 1
    xtargv = eps100
    tolsx = eps100

    icount = -1
100 continue
    xi1 = xi0 + delxi
    icount = icount + 1

    if (icount .gt. 10) then
        if (iodb(1).gt.0 .or. iodb(5).gt.0) then
            write (noutpt,1000)
1000 format(/3x,'A scan for the inverse rate becoming zero',' has failed.',/3x,'Dropping to the minimum step size.')
        end if

        delxi = dlxmin
        go to 999
    end if

    ! Estimate the inverse rate from a Taylor's series expansion.
    rirecp = rirec0
    dxp = 1.

    do n = 1,nord
        dxp = dxp*delxi
        rirecp = rirecp + ( drir0(n)/fctrl(n) )*dxp
    end do

    if (rirecp .ge. 0.) then
        go to 999
    end if

    if (delxi .gt. dlxmin) then
        if (rirec0 .le. (2.*eps100)) then
            ! The inverse rate is too small at the base point.
            delxi = dlxmin

            if (iodb(1).gt.0 .or. iodb(5).gt.0) then
                write (noutpt,1010) xi0,delxi
1010 format(/' --- The inverse rate is already very close to',/7x,'zero at Xi= ',1pe11.4,'. Have set delxi equal to the',/7x,'minimum value of ',e11.4,'.')

                go to 100
            end if
        end if
    end if

    if (delxi .gt. dlxmin) then
        xval0 = rirec0

        do n = 1,nord
            dxval0(n) = drir0(n)
        end do

        call search(delxi,dlxmin,dxval0,eps100,ier,ilsign,iodb,nodbmx,nord,noutpt,nrd1mx,nttyo,tolsx,unam24,usearch,xtargv,xval0)

        if (ier .le. 0) then
            go to 100
        end if

        if (ier .ge. 2) then
            ! Note: if ier = 1, the returned "safe" value of delxi
            ! is used.
            delxi = dlxmin
        end if

        go to 999
    end if

    if (iodb(1).gt.0 .or. iodb(5).gt.0) then
        if (abs(rirecp - xtargv) .le. tolsx) then
            write (noutpt,1020) xi1,delxi
1020 format(/3x,"Taylor's series predict that the inverse rate",/3x,'is zero at Xi= ',1pe11.4,', delxi= ',1pe11.4,'.')
        end if
    end if

999 continue
end subroutine chksir