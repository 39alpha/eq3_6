subroutine chkstc(actw,awmax,awmin,eh,ehmax,ehmin,fo2lg,iopt,jreac,kstep,kstpmx,noptmx,noutpt,nrct,nrctmx,nttyo,o2max,o2min,ph,phmax,phmin,prcinf,qaft1,qcnpre,qcntmp,qconst,qredox,qstop,qvhfxi,qvlsow,timemx,time1,tolxst,tolxsu,ximax,xi1)
    !! This subroutine checks for conditions which call for terminating
    !! the current reaction path calculation.
    !! This subroutine is called by:
    !!   EQ6/path.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: kstpmx
    integer :: noptmx
    integer :: nrctmx

    integer :: noutpt
    integer :: nttyo

    integer :: iopt(noptmx)
    integer :: jreac(nrctmx)

    integer :: kstep
    integer :: nrct

    logical :: qaft1
    logical :: qcnpre
    logical :: qcntmp
    logical :: qconst
    logical :: qredox
    logical :: qstop
    logical :: qvhfxi
    logical :: qvlsow

    real(kind=8) :: actw
    real(kind=8) :: awmax
    real(kind=8) :: awmin
    real(kind=8) :: eh
    real(kind=8) :: ehmax
    real(kind=8) :: ehmin
    real(kind=8) :: fo2lg
    real(kind=8) :: o2max
    real(kind=8) :: o2min
    real(kind=8) :: ph
    real(kind=8) :: phmax
    real(kind=8) :: phmin
    real(kind=8) :: prcinf
    real(kind=8) :: timemx
    real(kind=8) :: time1
    real(kind=8) :: tolxst
    real(kind=8) :: tolxsu
    real(kind=8) :: ximax
    real(kind=8) :: xi1

    ! Local variable declarations.
    integer :: nrc
    integer :: nrrct

    real(kind=8) :: tres

    qstop = .false.

    if (xi1 .ge. ximax) then
        write (noutpt,1000)
        write (nttyo,1000)
1000 format(/3x,'Have reached the maximum value of Xi.')

        qstop = .true.
    end if

    if (iopt(2) .ge. 1) then
        if (timemx .lt. prcinf) then
            tres = (time1 - timemx)/timemx

            if (abs(tres) .le. tolxst) then
                write (noutpt,1010)
                write (nttyo,1010)
1010 format(/3x,'Have reached the maximum value of time.')

                qstop = .true.
            end if
        end if
    end if

    if (phmin .gt. -prcinf) then
        if (abs(ph - phmin) .le. tolxsu) then
            write (noutpt,1020)
            write (nttyo,1020)
1020 format(/3x,'Have reached the minimum value of pH.')

            qstop = .true.
        else if (ph .lt. phmin) then
            write (noutpt,1030)
            write (nttyo,1030)
1030 format(/3x,'The pH is below the minimum value.')

            qstop = .true.
        end if
    end if

    if (phmax .lt. prcinf) then
        if (abs(ph - phmax) .le. tolxsu) then
            write (noutpt,1040)
            write (nttyo,1040)
1040 format(/3x,'Have reached the maximum value of pH.')

            qstop = .true.
        else if (ph .gt. phmax) then
            write (noutpt,1050)
            write (nttyo,1050)
1050 format(/3x,'The pH is above the maximum value.')

            qstop = .true.
        end if
    end if

    if (qredox) then
        if (ehmin .gt. -prcinf) then
            if (abs(eh - ehmin) .le. tolxsu) then
                write (noutpt,1060)
                write (nttyo,1060)
1060 format(/3x,'Have reached the minimum value of Eh (v).')

                qstop = .true.
            else if (eh .lt. ehmin) then
                write (noutpt,1070)
                write (nttyo,1070)
1070 format(/3x,'The Eh (v) is below the minimum value.')

                qstop = .true.
            end if
        end if

        if (ehmax .lt. prcinf) then
            if (abs(eh - ehmax) .le. tolxsu) then
                write (noutpt,1080)
                write (nttyo,1080)
1080 format(/3x,'Have reached the maximum value of Eh (v).')

                qstop = .true.
            else if (eh .gt. ehmax) then
                write (noutpt,1090)
                write (nttyo,1090)
1090 format(/3x,'The Eh (v) is above the maximum value.')

                qstop = .true.
            end if
        end if
    end if

    if (qredox) then
        if (o2min .gt. -prcinf) then
            if (abs(fo2lg - o2min) .le. tolxsu) then
                write (noutpt,1100)
                write (nttyo,1100)
1100 format(/3x,'Have reached the minimum value of log fO2.')

                qstop = .true.
            else if (fo2lg .lt. o2min) then
                write (noutpt,1110)
                write (nttyo,1110)
1110 format(/3x,'The log fO2 is below the minimum value.')

                qstop = .true.
            end if
        end if

        if (o2max .lt. prcinf) then
            if (abs(fo2lg - o2max) .le. tolxsu) then
                write (noutpt,1120)
                write (nttyo,1120)
1120 format(/3x,'Have reached the maximum value of log fO2.')

                qstop = .true.
            else if (fo2lg .gt. o2max) then
                write (noutpt,1130)
                write (nttyo,1130)
1130 format(/3x,'The log fO2 is above the maximum value.')

                qstop = .true.
            end if
        end if
    end if

    if (awmin .gt. -prcinf) then
        if (abs(actw - awmin) .le. tolxsu) then
            write (noutpt,1140)
            write (nttyo,1140)
1140 format(/3x,'Have reached the minimum value of the',' activity of water.')

            qstop = .true.
        else if (actw .lt. awmin) then
            write (noutpt,1150)
            write (nttyo,1150)
1150 format(/3x,'The activity of water is below the minimum',' value.')

            qstop = .true.
        end if
    end if

    if (awmax .lt. prcinf) then
        if (abs(actw - awmax) .le. tolxsu) then
            write (noutpt,1160)
            write (nttyo,1160)
1160 format(/3x,'Have reached the maximum value of the',' activity of water.')

            qstop = .true.
        else if (actw .gt. awmax) then
            write (noutpt,1170)
            write (nttyo,1170)
1170 format(/3x,'The activity of water is above the maximum',' value.')

            qstop = .true.
        end if
    end if

    if (kstep .ge. kstpmx) then
        write (noutpt,1200)
        write (nttyo,1200)
1200 format(/3x,'Have done the maximum number of steps.')

        qstop = .true.
    end if

    if (qconst) then
        write (noutpt,1400)
        write (nttyo,1400)
1400 format(/3x,'All irreversible reaction rates are now zero.')

        qstop=.true.
    end if

    if (qaft1) then
        write (noutpt,1410)
        write (nttyo,1410)
1410 format(/3x,'The total affinity is repeatedly nearly zero.')

        qstop=.true.
    end if

    ! Test to see if any reactants remain active. If none, terminate
    ! the run unless the temperature or pressure are changing.
    if (nrct .gt. 0) then
        if (qcntmp .and. qcnpre) then
            nrrct = 0

            do nrc = 1,nrct
                if (jreac(nrc) .eq. 0) then
                    nrrct = nrrct + 1
                end if
            end do

            if (nrrct .le. 0) then
                write (noutpt,1420)
                write (nttyo,1420)
1420 format(/3x,'Each reactant is now saturated or exhausted.')

                qstop = .true.
            end if
        end if
    end if

    if (qvlsow) then
        write (noutpt,1430)
        write (nttyo,1430)
1430 format(/3x,'Have very nearly exhausted solvent water.')

        qstop=.true.
    end if

    if (qvhfxi) then
        write (noutpt,1440)
        write (nttyo,1440)
1440 format(/3x,'Have extremely high ionic strength.')

        qstop=.true.
    end if
end subroutine chkstc