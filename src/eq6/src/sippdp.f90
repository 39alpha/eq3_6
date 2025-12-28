subroutine sippdp(actw,aw0plo,aw0prn,aw1plo,aw1prn,dlaplo,dlaprn,dleplo,dleprn,dlhplo,dlhprn,dloplo,dloprn,dltpll,dltplo,dltprl,dltprn,dlxdmp,dlxmx0,dlxpll,dlxplo,dlxprl,dlxprn,eh,eh0plo,eh0prn,eh1plo,eh1prn,eps100,fo2lg,lprcin,o20plo,o20prn,o21plo,o21prn,ph,ph0plo,ph0prn,ph1plo,ph1prn,prcinf,qredox,tiplol,tiplot,tiprnl,tiprnt,tistsv,xidump,xiplol,xiplot,xiprnl,xiprnt,xistsv)
    !! This subroutine advances the print, plot, and dump points.
    !! This subroutine is called by:
    !!   EQ6/path.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    logical :: qredox

    real(kind=8) :: actw
    real(kind=8) :: aw0plo
    real(kind=8) :: aw0prn
    real(kind=8) :: aw1plo
    real(kind=8) :: aw1prn
    real(kind=8) :: dlaplo
    real(kind=8) :: dlaprn
    real(kind=8) :: dleplo
    real(kind=8) :: dleprn
    real(kind=8) :: dlhplo
    real(kind=8) :: dlhprn
    real(kind=8) :: dloplo
    real(kind=8) :: dloprn
    real(kind=8) :: dltpll
    real(kind=8) :: dltplo
    real(kind=8) :: dltprl
    real(kind=8) :: dltprn
    real(kind=8) :: dlxdmp
    real(kind=8) :: dlxmx0
    real(kind=8) :: dlxpll
    real(kind=8) :: dlxplo
    real(kind=8) :: dlxprl
    real(kind=8) :: dlxprn
    real(kind=8) :: eh
    real(kind=8) :: eh0plo
    real(kind=8) :: eh0prn
    real(kind=8) :: eh1plo
    real(kind=8) :: eh1prn
    real(kind=8) :: eps100
    real(kind=8) :: fo2lg
    real(kind=8) :: lprcin
    real(kind=8) :: o20plo
    real(kind=8) :: o20prn
    real(kind=8) :: o21plo
    real(kind=8) :: o21prn
    real(kind=8) :: ph
    real(kind=8) :: ph0plo
    real(kind=8) :: ph0prn
    real(kind=8) :: ph1plo
    real(kind=8) :: ph1prn
    real(kind=8) :: prcinf
    real(kind=8) :: tiplol
    real(kind=8) :: tiplot
    real(kind=8) :: tiprnl
    real(kind=8) :: tiprnt
    real(kind=8) :: tistsv
    real(kind=8) :: xidump
    real(kind=8) :: xiplol
    real(kind=8) :: xiplot
    real(kind=8) :: xiprnl
    real(kind=8) :: xiprnt
    real(kind=8) :: xistsv

    ! Local variable declarations.
    integer :: n

    real(kind=8) :: lloc
    real(kind=8) :: tsmall
    real(kind=8) :: t1day
    real(kind=8) :: xloc
    real(kind=8) :: yloc

    ! XX   real*8 t1sec,t1min,t1hr,t1year
    real(kind=8) :: texp
    real(kind=8) :: tlg

    xiprnt = prcinf
    xiprnl = prcinf
    xiplot = prcinf
    xiplol = prcinf

    xidump = prcinf

    tiprnt = prcinf
    tiprnl = prcinf
    tiplot = prcinf
    tiplol = prcinf

    ph0prn = -prcinf
    eh0prn = -prcinf
    o20prn = -prcinf
    aw0prn = -prcinf
    ph0plo = -prcinf
    eh0plo = -prcinf
    o20plo = -prcinf
    aw0plo = -prcinf

    ph1prn = prcinf
    eh1prn = prcinf
    o21prn = prcinf
    aw1prn = prcinf
    ph1plo = prcinf
    eh1plo = prcinf
    o21plo = prcinf
    aw1plo = prcinf

    ! XX   t1sec= 1.0
    ! XX   t1min= 60.
    ! XX   t1hr = 3600.
    t1day = 86400.

    ! XX   t1year = t1day*365.25
    !      Calculate the initial print point for the specified Xi interval.
    if (dlxprn .lt. prcinf) then
        n = int( xistsv/dlxprn ) + 1
        xiprnt = n*dlxprn

        if ((xiprnt + eps100) .le. xistsv) then
            xiprnt = xiprnt + dlxprn
        end if
    end if

    ! Calculate the initial print point for the specified log Xi
    ! interval.
    if (dlxprl .lt. lprcin) then
        if (xistsv .lt. dlxmx0) then
            xloc = tlg(dlxmx0)
            lloc = int(xloc)
            n = int( ( xloc - lloc )/dlxprl )
            yloc = lloc + n*dlxprl

            if ((yloc + eps100) .le. xloc) then
                yloc = yloc + dlxprl
            end if
        else
            xloc = tlg(xistsv)
            lloc = int(xloc)
            n = int( ( xloc - lloc )/dlxprl ) + 1
            yloc = lloc + n*dlxprl

            if ((yloc + eps100) .le. xloc) then
                yloc = yloc + dlxprl
            end if
        end if

        xiprnl = texp(yloc)
    end if

    ! Calculate the initial plot point for the specified Xi interval.
    if (dlxplo .lt. prcinf) then
        n = int( xistsv/dlxplo ) + 1
        xiplot = n*dlxplo

        if ((xiplot + eps100) .le. xistsv) then
            xiplot = xiplot + dlxplo
        end if
    end if

    ! Calculate the initial plot point for the specified log Xi
    ! interval.
    if (dlxpll .lt. lprcin) then
        if (xistsv .lt. dlxmx0) then
            xloc = tlg(dlxmx0)
            lloc = int(xloc)
            n = int( ( xloc - lloc )/dlxpll )
            yloc = lloc + n*dlxpll

            if ((yloc + eps100) .le. xloc) then
                yloc = yloc + dlxpll
            end if
        else
            xloc = tlg(xistsv)
            lloc = int(xloc)
            n = int( ( xloc - lloc )/dlxpll ) + 1
            yloc = lloc + n*dlxpll

            if ((yloc + eps100) .le. xloc) then
                yloc = yloc + dlxpll
            end if
        end if

        xiplol = texp(yloc)
    end if

    ! Calculate the initial dump point for the specified Xi
    ! interval.
    if (dlxdmp .lt. prcinf) then
        n = int( xistsv/dlxdmp ) + 1
        xidump = n*dlxdmp

        if ((xidump + eps100) .le. xistsv) then
            xidump = xidump + dlxdmp
        end if
    end if

    ! Calculate the initial print point for the specified time interval.
    if (dltprn .lt. prcinf) then
        n = int( tistsv/dltprn ) + 1
        tiprnt = n*dltprn

        if ((tiprnt + eps100) .le. tistsv) then
            tiprnt = tiprnt + dltprn
        end if
    end if

    ! Calculate the initial print point for the specified log time
    ! interval.
    if (dltprl .lt. lprcin) then
        tsmall = t1day

        if (tistsv .lt. tsmall) then
            xloc = tlg(tsmall)
            lloc = int(xloc)
            n = int( ( xloc - lloc )/dltprl )
            yloc = lloc + n*dltprl

            if ((yloc + eps100) .le. xloc) then
                yloc = yloc + dltprl
            end if
        else
            xloc = tlg(tistsv)
            lloc = int(xloc)
            n = int( ( xloc - lloc )/dltprl ) + 1
            yloc = lloc + n*dltprl

            if ((yloc + eps100) .le. xloc) then
                yloc = yloc + dltprl
            end if
        end if

        tiprnl = texp(yloc)
    end if

    ! Calculate the initial plot point for the specified time interval.
    if (dltplo .lt. prcinf) then
        n = int( tistsv/dltplo ) + 1
        tiplot = n*dltplo

        if ((tiplot + eps100) .le. tistsv) then
            tiplot = tiplot + dltplo
        end if
    end if

    ! Calculate the initial plot point for the specified log time
    ! interval.
    if (dltpll .lt. lprcin) then
        tsmall = t1day

        if (tistsv .lt. tsmall) then
            xloc = tlg(tsmall)
            lloc = int(xloc)
            n = int( ( xloc - lloc )/dltpll )
            yloc = lloc + n*dltpll

            if ((yloc + eps100) .le. xloc) then
                yloc = yloc + dltpll
            end if
        else
            xloc = tlg(tistsv)
            lloc = int(xloc)
            n = int( ( xloc - lloc )/dltpll ) + 1
            yloc = lloc + n*dltpll

            if ((yloc + eps100) .le. xloc) then
                yloc = yloc + dltpll
            end if
        end if

        tiplol = texp(yloc)
    end if

    ! Calculate the initial print points for the specified pH interval.
    ! Note there are two such points, one lower and one higher than the
    ! current value.
    if (dlhprn .lt. prcinf) then
        n = int( ph/dlhprn )
        ph0prn = n*dlhprn

        if ((ph0prn - eps100) .ge. ph) then
            ph0prn = ph0prn - dlhprn
        end if

        ph1prn = (n + 1)*dlhprn

        if ((ph1prn + eps100) .le. ph) then
            ph1prn = ph1prn + dlhprn
        end if
    end if

    ! Calculate the initial print points for the specified aw interval.
    ! Note there are two such points, one lower and one higher than the
    ! current value.
    if (dlaprn .lt. prcinf) then
        n = int( actw/dlaprn )
        aw0prn = n*dlaprn

        if ((aw0prn - eps100) .ge. actw) then
            aw0prn = aw0prn - dlaprn
        end if

        aw1prn = (n + 1)*dlaprn

        if ((aw1prn + eps100) .le. actw) then
            aw1prn = aw1prn + dlaprn
        end if
    end if

999 continue
end subroutine sippdp