subroutine chdxec(delxi,dlxmx0,dzvc0,iodb,kdim,kmax,km1,kxt,nodbmx,nord,noutpt,nrd1mx,qmin,qscon,scale,scalim,scnstd,scnsti,uzvec1,zvec0)
    !! This subroutine chooses the step size (delxi) for economy and
    !! super economy modes.
    !! This subroutine is called by:
    !!   EQ6/path.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: kmax
    integer :: nodbmx
    integer :: nrd1mx

    integer :: noutpt

    integer :: iodb(nodbmx)

    integer :: kdim
    integer :: km1
    integer :: kxt
    integer :: nord

    logical :: qmin
    logical :: qscon

    character(len=48) :: uzvec1(kmax)

    real(kind=8) :: dzvc0(nrd1mx,kmax)
    real(kind=8) :: zvec0(kmax)

    real(kind=8) :: delxi
    real(kind=8) :: dlxmx0
    real(kind=8) :: scale
    real(kind=8) :: scalim
    real(kind=8) :: scnstd
    real(kind=8) :: scnsti

    ! Local variable declarations.
    integer :: j2
    integer :: j3
    integer :: kcol
    integer :: kcolcn

    integer :: ilnobl

    character(len=24) :: ubas
    character(len=24) :: uzt
    character(len=24) :: uzn

    real(kind=8) :: aterm
    real(kind=8) :: d1
    real(kind=8) :: d2
    real(kind=8) :: d1sq
    real(kind=8) :: mxx0
    real(kind=8) :: rd1
    real(kind=8) :: rd2
    real(kind=8) :: ri1
    real(kind=8) :: ri2
    real(kind=8) :: rx
    real(kind=8) :: scfac
    real(kind=8) :: scnst
    real(kind=8) :: sx
    real(kind=8) :: termd
    real(kind=8) :: termi
    real(kind=8) :: thetad
    real(kind=8) :: thetai

    data ubas /'basis variable          '/

    if (qscon) then
        ! Choose the step size for super economy mode.
        nord = 0
        scale = dlxmx0/delxi
        go to 999
    end if

    ! Choose step size for ordinary economy mode.
    if (nord .le. 0) then
        scale = dlxmx0/delxi
        go to 999
    end if

    uzt = 'None'
    uzn = ' '
    kcolcn = 0
    scfac = scale

    do kcol = 1,kdim
        qmin = kcol.ge.km1 .and. kcol.le.kxt
        d1 = dzvc0(1,kcol)
        d2 = dzvc0(2,kcol)
        mxx0 = zvec0(kcol)

        if (nord.le.1 .or. d2.eq.0.) then
            ! First order bound.
            if (d1 .eq. 0.) then
                go to 100
            end if

            if (d1 .gt. 0.) then
                scnst = scnsti
            else
                scnst = scnstd
            end if

            sx = (mxx0*scnst) / (d1*delxi)

            if (iodb(5) .ge. 2) then
                if (sx .lt. scale) then
                    rx = sx*delxi
                    write (noutpt,1000) uzvec1(kcol),rx
1000 format(3x,'Element= ',a,/7x,'Root= ',1pe12.5)

                    write (noutpt,1010) dzvc0(1,kcol)
1010 format(7x,'Derivatives= ',4(1pe12.5,2x))
                end if
            end if
        else
            ! Second order bound. There are two quadratic equations to be
            ! solved, each of which has two roots. Only the smallest
            ! positive root is of interest. If one root of an equation is
            ! complex, the other will also be complex. If both roots are
            ! real, they could be positive and negative, both positive, or
            ! both negative.
            d1sq = d1*d1
            thetai = mxx0*scnsti
            thetad = mxx0*scnstd
            termi = 2.*d2*thetai
            termd = 2.*d2*thetad

            ! Here ri1, ri2, rd1, and rd2 are values of Xi which are roots
            ! of the equations.
            aterm = d1sq - termi

            if (aterm .ge. 0.) then
                ri1 = (d1 + sqrt(aterm))/d2

                if (ri1.le.0.) then
                    ri1 = 1.e+30
                end if

                ri2 = (d1 - sqrt(aterm))/d2

                if (ri2.le.0.) then
                    ri2=1.e+30
                end if
            else
                ri1 = 1.e+30
                ri2 = 1.e+30
            end if

            aterm = d1sq - termd

            if (aterm .ge. 0.) then
                rd1 = (d1 + sqrt(aterm))/d2

                if (rd1.le.0.) then
                    rd1 = 1.e+30
                end if

                rd2 = (d1 - sqrt(aterm))/d2

                if (rd2.le.0.) then
                    rd2 = 1.e+30
                end if
            else
                rd1 = 1.e+30
                rd2 = 1.e+30
            end if

            rx = min(ri1,ri2,rd1,rd2)

            if (rx .lt. 1.e+30) then
                sx = rx/delxi
            else
                sx = 1.e+30
            end if

            if (iodb(5) .ge. 2) then
                if (sx .lt. scale) then
                    j2 = ilnobl(uzvec1(kcol))
                    write (noutpt,1000) uzvec1(kcol)(1:j2),rx
                    write (noutpt,1010) dzvc0(1,kcol),dzvc0(2,kcol)
                end if
            end if
        end if

        if (sx .lt. scfac) then
            scfac = sx
            kcolcn = kcol
            uzt(1:24) = ubas(1:24)
        end if

100 continue
    end do

    if (kcolcn .gt. 0) then
        uzn(1:24) = uzvec1(kcolcn)(1:24)
    end if

    if (iodb(5) .ge. 1) then
        write (noutpt,1020) scfac
1020 format('   Economy mode: scale= ',f10.5)

        j2 = ilnobl(uzt)
        j3 = ilnobl(uzn)

        if (j3 .gt. 0) then
            write (noutpt,1030) scfac,uzt(1:j2),uzn(1:j3)
1030 format(5x,'Constrained by ',a,1x,a)
        else
            write (noutpt,1040) uzt(1:j2)
1040 format(5x,'Constrained by ',a)
        end if
    end if

    scale = min(scalim,scfac)

999 continue
end subroutine chdxec