subroutine adprpl(actw,aw0plo,aw0prn,aw1plo,aw1prn,dlaplo,dlaprn,dleplo,dleprn,dlhplo,dlhprn,dloplo,dloprn,dltpll,dltplo,dltprl,dltprn,dlxpll,dlxplo,dlxprl,dlxprn,dlxdmp,eh,eh0plo,eh0prn,eh1plo,eh1prn,fo2lg,fxprpl,o20plo,o20prn,o21plo,o21prn,ph,ph0plo,ph0prn,ph1plo,ph1prn,qplaw0,qplaw1,qpleh0,qpleh1,qplolt,qplolx,qplott,qplotx,qplo20,qplo21,qplph0,qplph1,qpraw0,qpraw1,qpreh0,qpreh1,qprnlx,qprnlt,qprntx,qprntt,qpro20,qpro21,qprph0,qprph1,qredox,time1,tiplol,tiplot,tiprnl,tiprnt,xidump,xiplol,xiplot,xiprnl,xiprnt,xi1)
    !! This subroutine adjusts print, plot, and PRS transfer points
    !! as needed for the next step along the reaction path.
    !! This subroutine is called by:
    !!   EQ6/path.f
    !! Principal input:
    !! Principal output:
    implicit none

    logical :: qplaw0
    logical :: qplaw1
    logical :: qpleh0
    logical :: qpleh1
    logical :: qplolt
    logical :: qplolx
    logical :: qplott
    logical :: qplotx
    logical :: qplo20
    logical :: qplo21
    logical :: qplph0
    logical :: qplph1
    logical :: qpraw0
    logical :: qpraw1
    logical :: qpreh0
    logical :: qpreh1
    logical :: qprnlx
    logical :: qprnlt
    logical :: qprntx
    logical :: qprntt
    logical :: qpro20
    logical :: qpro21
    logical :: qprph0
    logical :: qprph1
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
    real(kind=8) :: dlxpll
    real(kind=8) :: dlxplo
    real(kind=8) :: dlxprl
    real(kind=8) :: dlxprn
    real(kind=8) :: dlxdmp
    real(kind=8) :: eh
    real(kind=8) :: eh0plo
    real(kind=8) :: eh0prn
    real(kind=8) :: eh1plo
    real(kind=8) :: eh1prn
    real(kind=8) :: fo2lg
    real(kind=8) :: fxprpl
    real(kind=8) :: o20plo
    real(kind=8) :: o20prn
    real(kind=8) :: o21plo
    real(kind=8) :: o21prn
    real(kind=8) :: ph
    real(kind=8) :: ph0plo
    real(kind=8) :: ph0prn
    real(kind=8) :: ph1plo
    real(kind=8) :: ph1prn
    real(kind=8) :: time1
    real(kind=8) :: tiplol
    real(kind=8) :: tiplot
    real(kind=8) :: tiprnl
    real(kind=8) :: tiprnt
    real(kind=8) :: xidump
    real(kind=8) :: xiplol
    real(kind=8) :: xiplot
    real(kind=8) :: xiprnl
    real(kind=8) :: xiprnt
    real(kind=8) :: xi1

    ! Local variable declarations.
    real(kind=8) :: tolx
    real(kind=8) :: tx
    real(kind=8) :: xx

    real(kind=8) :: texp
    real(kind=8) :: tlg

    ! Reset print points as necessary.
    if (qprntx) then
        tolx = fxprpl*dlxprn
100 continue
        xiprnt = xiprnt + dlxprn

        if ((xi1 + tolx) .ge. xiprnt) then
            go to 100
        end if
    end if

    if (qprntt) then
        tolx = fxprpl*dlxprn
110 continue
        tiprnt = tiprnt + dltprn

        if ((time1 + tolx) .ge. tiprnt) then
            go to 110
        end if
    end if

    if (qprnlx) then
        tolx = fxprpl*dlxprl
120 continue
        xiprnl = texp(tlg(xiprnl) + dlxprl)
        xx = texp(tlg(xi1) + tolx)

        if (xx .ge. xiprnl) then
            go to 120
        end if
    end if

    if (qprnlt) then
        tolx = fxprpl*dltprl
130 continue
        tiprnl = texp(tlg(tiprnl) + dltprl)
        tx = texp(tlg(time1) + tolx)

        if (tx .ge. tiprnl) then
            go to 130
        end if
    end if

    if (qprph0) then
        tolx = fxprpl*dlhprn
140 continue
        ph0prn = ph0prn - dlhprn

        if ((ph - tolx) .le. ph0prn) then
            go to 140
        end if

        ph1prn = ph0prn
150 continue
        ph1prn = ph1prn + dlhprn

        if ((ph + tolx) .ge. ph1prn) then
            go to 150
        end if
    end if

    if (qprph1) then
        tolx = fxprpl*dlhprn
160 continue
        ph1prn = ph1prn + dlhprn

        if ((ph + tolx) .ge. ph1prn) then
            go to 160
        end if

        ph0prn = ph1prn
170 continue
        ph0prn = ph0prn - dlhprn

        if ((ph - tolx) .le. ph0prn) then
            go to 170
        end if
    end if

    if (qredox) then
        if (qpreh0) then
            tolx = fxprpl*dleprn
180 continue
            eh0prn = eh0prn - dleprn

            if ((eh - tolx) .le. eh0prn) then
                go to 180
            end if

            eh1prn = eh0prn
190 continue
            eh1prn = eh1prn + dleprn

            if ((eh + tolx) .ge. eh1prn) then
                go to 190
            end if
        end if

        if (qpreh1) then
            tolx = fxprpl*dleprn
200 continue
            eh1prn = eh1prn + dleprn

            if ((eh + tolx) .ge. eh1prn) then
                go to 200
            end if

            eh0prn = eh1prn
210 continue
            eh0prn = eh0prn - dleprn

            if ((eh - tolx) .le. eh0prn) then
                go to 210
            end if
        end if

        if (qpro20) then
            tolx = fxprpl*dloprn
220 continue
            o20prn = o20prn - dloprn

            if ((fo2lg - tolx) .le. o20prn) then
                go to 220
            end if

            o21prn = o20prn
230 continue
            o21prn = o21prn + dloprn

            if ((fo2lg + tolx) .ge. o21prn) then
                go to 230
            end if
        end if

        if (qpro21) then
            tolx = fxprpl*dloprn
240 continue
            o21prn = o21prn + dloprn

            if ((fo2lg + tolx) .ge. o21prn) then
                go to 240
            end if

            o20prn = o21prn
250 continue
            o20prn = o20prn - dloprn

            if ((fo2lg - tolx) .le. o20prn) then
                go to 250
            end if
        end if
    end if

    if (qpraw0) then
        tolx = fxprpl*dlaprn
260 continue
        aw0prn = aw0prn - dlaprn

        if ((actw - tolx) .le. aw0prn) then
            go to 260
        end if

        aw1prn = aw0prn
270 continue
        aw1prn = aw1prn + dlaprn

        if ((actw + tolx) .ge. aw1prn) then
            go to 270
        end if
    end if

    if (qpraw1) then
        tolx = fxprpl*dlaprn
280 continue
        aw1prn = aw1prn + dlaprn

        if ((actw + tolx) .ge. aw1prn) then
            go to 280
        end if

        aw0prn = aw1prn
290 continue
        aw0prn = aw0prn - dlaprn

        if ((actw - tolx) .le. aw0prn) then
            go to 290
        end if
    end if

    ! Reset plot points as necessary.
    if (qplotx) then
        tolx = fxprpl*dlxplo
400 continue
        xiplot = xiplot + dlxplo

        if ((xi1 + tolx)  .ge. xiplot) then
            go to 400
        end if
    end if

    if (qplott) then
        tolx = fxprpl*dltplo
410 continue
        tiplot = tiplot + dltplo

        if ((time1 + tolx)  .ge. tiplot) then
            go to 410
        end if
    end if

    if (qplolx) then
        tolx = fxprpl*dlxpll
420 continue
        xiplol = texp(tlg(xiplol) + dlxpll)
        xx = texp(tlg(xi1) + tolx)

        if (xx  .ge. xiplol) then
            go to 420
        end if
    end if

    if (qplolt) then
        tolx = fxprpl*dltpll
430 continue
        tiplol = texp(tlg(tiplol) + dltpll)
        tx = texp(tlg(time1) + tolx)

        if (tx .ge. tiplol) then
            go to 430
        end if
    end if

    if (qplph0) then
        tolx = fxprpl*dlhplo
440 continue
        ph0plo = ph0plo - dlhplo

        if ((ph - tolx) .le. ph0plo) then
            go to 440
        end if

        ph1plo = ph0plo
450 continue
        ph1plo = ph1plo + dlhplo

        if ((ph + tolx) .ge. ph1plo) then
            go to 450
        end if
    end if

    if (qplph1) then
        tolx = fxprpl*dlhplo
460 continue
        ph1plo = ph1plo + dlhplo

        if ((ph + tolx) .ge. ph1plo) then
            go to 460
        end if

        ph0plo = ph1plo
470 continue
        ph0plo = ph0plo - dlhplo

        if ((ph - tolx) .le. ph0plo) then
            go to 470
        end if
    end if

    if (qredox) then
        if (qpleh0) then
            tolx = fxprpl*dleplo
480 continue
            eh0plo = eh0plo - dleplo

            if ((eh - tolx) .le. eh0plo) then
                go to 480
            end if

            eh1plo = eh0plo
490 continue
            eh1plo = eh1plo + dleplo

            if ((eh + tolx) .ge. eh1plo) then
                go to 490
            end if
        end if

        if (qpleh1) then
            tolx = fxprpl*dleplo
500 continue
            eh1plo = eh1plo + dleplo

            if ((eh + tolx) .ge. eh1plo) then
                go to 500
            end if

            eh0plo = eh1plo
510 continue
            eh0plo = eh0plo - dleplo

            if ((eh - tolx) .le. eh0plo) then
                go to 510
            end if
        end if

        if (qplo20) then
            tolx = fxprpl*dloplo
520 continue
            o20plo = o20plo - dloplo

            if ((fo2lg - tolx) .le. o20plo) then
                go to 520
            end if

            o21plo = o20plo
530 continue
            o21plo = o21plo + dloplo

            if ((fo2lg + tolx) .ge. o21plo) then
                go to 530
            end if
        end if

        if (qplo21) then
            tolx = fxprpl*dloplo
540 continue
            o21plo = o21plo + dloplo

            if ((fo2lg + tolx) .ge. o21plo) then
                go to 540
            end if

            o20plo = o21plo
550 continue
            o20plo = o20plo - dloplo

            if ((fo2lg - tolx) .le. o20plo) then
                go to 550
            end if
        end if
    end if

    if (qplaw0) then
        tolx = fxprpl*dlaplo
560 continue
        aw0plo = aw0plo - dlaplo

        if ((actw - tolx) .le. aw0plo) then
            go to 560
        end if

        aw1plo = aw0plo
570 continue
        aw1plo = aw1plo + dlaplo

        if ((actw + tolx) .ge. aw1plo) then
            go to 570
        end if
    end if

    if (qplaw1) then
        tolx = fxprpl*dlaplo
580 continue
        aw1plo = aw1plo + dlaplo

        if ((actw + tolx) .ge. aw1plo) then
            go to 580
        end if

        aw0plo = aw1plo
590 continue
        aw0plo = aw0plo - dlaplo

        if ((actw - tolx) .le. aw0plo) then
            go to 590
        end if
    end if

    ! Reset the PRS transfer point as necessary.
    if (xidump .le. xi1) then
        tolx = fxprpl*dlxdmp
700 continue
        xidump = xidump + dlxdmp

        if ((xi1 + tolx)  .ge. xidump) then
            go to 700
        end if
    end if
end subroutine adprpl