subroutine chksrr(delxi,dlxmin,drer0,dxval0,eps100,iodb,jreac,nodbmx,noutpt,nord,nrct,nrctmx,nrd1mx,nttyo,rrelr0,rrelrp,tolsrr,ureac,xi0,xi1,xval0)
    !! This subroutine checks the signs of the relative rates. It finds
    !! the point of reaction progress at which any relative rate of an
    !! irreversible reaction becomes zero. The relative rates are
    !! tracked using finite differences.
    !! This subroutine is called by:
    !!   EQ6/path.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: nodbmx
    integer :: nrctmx
    integer :: nrd1mx

    integer :: noutpt
    integer :: nttyo

    integer :: iodb(nodbmx)
    integer :: jreac(nrctmx)

    integer :: nord
    integer :: nrct

    character(len=24) :: ureac(nrctmx)

    real(kind=8) :: drer0(nrd1mx,nrctmx)
    real(kind=8) :: dxval0(nrd1mx)
    real(kind=8) :: rrelr0(nrctmx)
    real(kind=8) :: rrelrp(nrctmx)

    real(kind=8) :: delxi
    real(kind=8) :: dlxmin
    real(kind=8) :: eps100
    real(kind=8) :: tolsrr
    real(kind=8) :: xi0
    real(kind=8) :: xi1
    real(kind=8) :: xval0

    ! Local variable declarations.
    integer :: icount
    integer :: ier
    integer :: ilsign
    integer :: j2
    integer :: krzero
    integer :: n
    integer :: nrc
    integer :: nrzero

    integer :: ilnobl

    character(len=48) :: usearch
    character(len=24) :: unam24

    real(kind=8) :: arrxp
    real(kind=8) :: atgsrr
    real(kind=8) :: dxp
    real(kind=8) :: rrx
    real(kind=8) :: rrx0
    real(kind=8) :: rrxp
    real(kind=8) :: tolsx
    real(kind=8) :: xtargv
    real(kind=8) :: xval

    real(kind=8) :: fctrl

    data usearch /'a relative rate changes sign                    '/

    if (nord .le. 0) then
        go to 999
    end if

    ! The search target is not where a relative rate is zero. The
    ! magnitude of the search target is atgsrr = 0.5*tolsrr. The search
    ! tolerance has the same value. Thus, in the case of a relative rate
    ! going from positive to negative values, the search target is
    ! -atgsrr and the interval for convergence is (-tolsrr,0.). In the
    ! case of a relative rate going from negative to positive values,
    ! the search target is atgsrr, and the interval for convergence is
    ! (0.,tolsrr).
    atgsrr = 0.5*tolsrr
    tolsx = atgsrr

    icount = -1
100 continue
    xi1 = xi0 + delxi
    icount = icount + 1

    if (icount .gt. 50) then
        if (iodb(1).gt.0 .or. iodb(5).gt.0) then
            write (noutpt,1000)
1000 format(/3x,'A scan for reactant relative rates changing sign',' has failed.',/3x,'Dropping to the minimum step size.')
        end if

        delxi = dlxmin
        go to 999
    end if

    ! Estimate the relative rates from Taylor's series expansions.
    do nrc = 1,nrct
        rrx = 0.

        if (jreac(nrc).eq.0 .or. jreac(nrc).eq.-1) then
            rrx = rrelr0(nrc)
            dxp = 1.

            do n = 1,nord
                dxp = dxp*delxi
                rrx = rrx + ( drer0(n,nrc)/fctrl(n) )*dxp
            end do
        end if

        rrelrp(nrc) = rrx
    end do

    ! Find any cross-overs. Note that there are two kinds, positive to
    ! negative, and negative to positive.
    xval = 0.
    nrzero = 0
    krzero = 0
    unam24 = ' '

    do nrc = 1,nrct
        if (jreac(nrc).eq.0 .or. jreac(nrc).eq.-1) then
            rrx0 = rrelr0(nrc)
            rrxp = rrelrp(nrc)

            if (rrx0 .gt. 0.) then
                if (rrxp .lt. -tolsrr) then
                    ! A case has been found in which the cross over tolerance
                    ! is exceeded for the positive to negative case.
                    krzero = krzero + 1
                    arrxp = -rrxp

                    if (arrxp .gt. xval) then
                        xval = arrxp
                        nrzero = nrc
                        unam24 = ureac(nrc)
                        ilsign = 1
                        xtargv = -atgsrr
                    end if
                end if
            else if (rrx0 .lt. 0.) then
                if (rrxp .gt. tolsrr) then
                    ! A case has been found in which the cross over tolerance
                    ! is exceeded for the negative to positive case.
                    krzero = krzero + 1
                    arrxp = rrxp

                    if (arrxp .gt. xval) then
                        xval = arrxp
                        nrzero = nrc
                        unam24 = ureac(nrc)
                        ilsign = -1
                        xtargv = atgsrr
                    end if
                end if
            end if
        end if
    end do

    if (krzero .gt. 0) then
        if (delxi .gt. dlxmin) then
            if (abs(xval - xtargv) .gt. tolsx) then
                xval0 = rrelr0(nrzero)

                do n = 1,nord
                    dxval0(n) = drer0(n,nrzero)
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
        end if
    end if

    if (iodb(1).gt.0 .or. iodb(5).gt.0) then
        if (krzero .gt. 0) then
            if (abs(xval - xtargv) .le. tolsx) then
                j2 = ilnobl(unam24)
                write (noutpt,1010) unam24(1:j2),xi1,delxi
1010 format(/" --- Taylor's series predict a zero relative",/7x,'rate for ',a,' at Xi= ',1pe11.4,',',/7x,'delxi= ',1pe11.4,' ---')
            end if
        end if
    end if

999 continue
end subroutine chksrr