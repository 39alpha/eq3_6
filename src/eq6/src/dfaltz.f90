subroutine dfaltz(dlaplo,dlaprn,dleplo,dleprn,dlhplo,dlhprn,dloplo,dloprn,dltplo,dltpll,dltprl,dltprn,dlxdmp,dlxmx0,dlxplo,dlxpll,dlxprl,dlxprn,iopt,itermx,ksplmx,ksppmx,kstpmx,net,noptmx,nordmx,nrct,noutpt,ntrymx,nttyo,prcinf,qecon,qscon,timmxi,tistti,tolbt,toldl,tolsat,tolxsf,tolxst,tolxsu,ximaxi,xistti)
    !! This subroutine sets the defaults for various run parameters read
    !! from the input file. In some cases, it forces the parameters to
    !! take on certain values, or to fall in certain ranges.
    !! This subroutine is called by:
    !!   EQ6/eq6.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: noptmx

    integer :: noutpt
    integer :: nttyo

    integer :: iopt(noptmx)

    integer :: itermx
    integer :: ksplmx
    integer :: ksppmx
    integer :: kstpmx
    integer :: net
    integer :: nordmx
    integer :: nrct
    integer :: ntrymx

    logical :: qecon
    logical :: qscon

    real(kind=8) :: dlaplo
    real(kind=8) :: dlaprn
    real(kind=8) :: dleplo
    real(kind=8) :: dleprn
    real(kind=8) :: dlhplo
    real(kind=8) :: dlhprn
    real(kind=8) :: dloplo
    real(kind=8) :: dloprn
    real(kind=8) :: dltplo
    real(kind=8) :: dltpll
    real(kind=8) :: dltprl
    real(kind=8) :: dltprn
    real(kind=8) :: dlxdmp
    real(kind=8) :: dlxmx0
    real(kind=8) :: dlxplo
    real(kind=8) :: dlxpll
    real(kind=8) :: dlxprl
    real(kind=8) :: dlxprn
    real(kind=8) :: prcinf
    real(kind=8) :: timmxi
    real(kind=8) :: tistti
    real(kind=8) :: tolbt
    real(kind=8) :: toldl
    real(kind=8) :: tolsat
    real(kind=8) :: tolxsf
    real(kind=8) :: tolxst
    real(kind=8) :: tolxsu
    real(kind=8) :: ximaxi
    real(kind=8) :: xistti

    ! Local variable declarations.
    !   None
    ! Iopt(1): Physical system model option.
    if (iopt(1) .lt. 0) then
        write (noutpt,1000) iopt(1)
        write (nttyo,1000) iopt(1)
1000 format(/' * Note - (EQ6/dfaltz) The value of iopt(1) is ',i3,',',/7x,'which is out of range. Resetting iopt(1) to 0.')

        iopt(1) = 0
    end if

    if (iopt(1) .gt. 2) then
        write (noutpt,1010) iopt(1)
        write (nttyo,1010) iopt(1)
1010 format(/' * Note - (EQ6/dfaltz) The value of iopt(1) is ',i3,',',/7x,'which is out of range. Resetting iopt(1) to 2.')

        iopt(1) = 2
    end if

    ! Iopt(3): Phase boundary search option.
    if (iopt(3).gt.0 .and. iopt(1).ge.2) then
        write (noutpt,1020) iopt(3),iopt(1)
        write (nttyo,1020) iopt(3),iopt(1)
1020 format(/' * Note - (EQ6/dfaltz) The value of iopt(3) is ',i3,',',/7x,'which is not valid because iopt(1) is ',i3,'.',/7x,'Resetting iopt(3) = 0.')

        iopt(3) = 0
    end if

    if (iopt(3).gt.0 .and. iopt(2).ge.1) then
        write (noutpt,1030) iopt(3),iopt(2)
        write (nttyo,1030) iopt(3),iopt(2)
1030 format(/' * Note - (EQ6/dfaltz) The value of iopt(3) is ',i3,',',/7x,'which is not valid because iopt(2) is ',i3,'.',/7x,'Resetting iopt(3) = 0.')

        iopt(3) = 0
    end if

    ! Maximum value of Xi.
    if (ximaxi .lt. xistti) then
        write (noutpt,1040)
        write (nttyo,1040)
1040 format(/' * Note - (EQ6/dfaltz) The maximum value of reaction',/7x,'progress (ximaxi) is less than the starting value',' (xistti).',/7x,'The maximum value is being reset to',' infinity.')

        ximaxi = prcinf
    end if

    if (ximaxi .gt. prcinf) then
        ximaxi = prcinf
    end if

    ! Maximum value of time.
    if (timmxi .lt. tistti) then
        write (noutpt,1050)
        write (nttyo,1050)
1050 format(/' * Note - (EQ6/dfaltz) The maximum value of time',/7x,'(timmxi) is less than the starting value (tistti).',/7x,'The maximum value is being reset to infinity.')

        timmxi = prcinf
    end if

    if (timmxi .gt. prcinf) then
        timmxi = prcinf
    end if

    ! Maximum number of steps.
    if (kstpmx .lt. 0) then
        kstpmx = 0
    end if

    ! Newton-Raphson beta tolerance.
    if (iopt(2) .le. 0) then
        if (tolbt .le. 0.) then
            tolbt = 1.e-6
        end if

        if (tolbt .lt. 1.e-10) then
            tolbt = 1.e-10
        end if

        if (tolbt .gt. 1.e-4) then
            tolbt = 1.e-4
        end if
    else
        if (tolbt .le. 0.) then
            tolbt = 1.e-8
        end if

        if (tolbt .lt. 1.e-10) then
            tolbt = 1.e-10
        end if

        if (tolbt .gt. 1.e-6) then
            tolbt = 1.e-6
        end if
    end if

    ! Newton-Raphson del tolerance.
    if (iopt(2) .le. 0) then
        if (toldl .le. 0.) then
            toldl = 1.e-6
        end if

        if (toldl .lt. 1.e-10) then
            toldl = 1.e-10
        end if

        if (toldl .gt. 1.e-4) then
            toldl = 1.e-4
        end if
    else
        if (toldl .le. 0.) then
            toldl = 1.e-8
        end if

        if (toldl .lt. 1.e-10) then
            toldl = 1.e-10
        end if

        if (toldl .gt. 1.e-6) then
            toldl = 1.e-6
        end if
    end if

    ! Maximum finite-difference order.
    if (nordmx .eq. 0) then
        nordmx = 6
    end if

    if (nordmx .lt. 1) then
        nordmx = 1
    end if

    if (nordmx .gt. 8) then
        nordmx = 8
    end if

    ! Maximum number of Newton-Raphson iterations.
    if (itermx .le. 0) then
        itermx = 200
    end if

    ! Maximum number of tries to find the correct phase assemblage.
    if (ntrymx .le. 0) then
        ntrymx = 100
    end if

    ! Ordinary search/find tolerance.
    if (tolxsf .le. 0.) then
        tolxsf = 1.e-6
    end if

    ! Search tolerance on time (relative).
    if (tolxst .le. 0.) then
        tolxst = 1.e-8
    end if

    if (tolxst .lt. 1.e-12) then
        tolxst = 1.e-12
    end if

    if (tolxst .gt. 1.e-4) then
        tolxst = 1.e-4
    end if

    ! Search tolerance on pH, Eh, log fO2, and activity of water
    ! (linear).
    if (tolxsu .le. 0.) then
        tolxsu = 1.e-5
    end if

    if (tolxsu .lt. 1.e-8) then
        tolxsu = 1.e-8
    end if

    if (tolxsu .gt. 1.e-4) then
        tolxsu = 1.e-4
    end if

    ! Saturation tolerance.
    if (tolsat .le. 0.) then
        tolsat = 0.0005
    end if

    if (tolsat .lt. 0.00005) then
        tolsat = 0.00005
    end if

    if (tolsat .gt. 0.05) then
        tolsat = 0.05
    end if

    ! Limits on Newton-Raphson delta tolerance.
    if (tolxsf .lt. 1.e-10) then
        toldl = 1.e-10
    end if

    if (tolxsf .gt. 1.e-2) then
        toldl = 1.e-2
    end if

    ! Zero-order step size.
    if (dlxmx0 .le. 0.) then
        dlxmx0 = 1.e-8

        if (iopt(2) .ge. 1) then
            dlxmx0 = 1.e-9
        end if

        if (iopt(1) .ge. 2) then
            dlxmx0 = 1.e-9
        end if

        if (qecon) then
            dlxmx0 = 1.e-6
        end if

        if (nrct .le. 0) then
            dlxmx0 = 1.e-2
        end if

        if (qscon) then
            if (dlxprn.gt.0.) then
                dlxmx0 = dlxprn
            else if (ximaxi .gt. xistti) then
                dlxmx0 = 0.1*(ximaxi - xistti)
            end if
        end if
    end if

    if (dlxmx0 .lt. 1.e-12) then
        dlxmx0 = 1.e-12
    end if

    ! Print interval in Xi.
    if (dlxprn .le. 0.) then
        dlxprn = prcinf
    end if

    ! Print interval in log Xi.
    if (dlxprl .le. 0.) then
        dlxprl = 0.5
    end if

    ! Print interval in time (seconds).
    if (dltprn .le. 0.) then
        dltprn = prcinf
    end if

    ! Print interval in log time.
    if (dltprl .le. 0.) then
        dltprl = prcinf
    end if

    ! Print interval in pH.
    if (dlhprn .le. 0.) then
        dlhprn = prcinf
    end if

    ! Print interval in Eh (v).
    if (dleprn .le. 0.) then
        dleprn = prcinf
    end if

    ! Print interval in log fO2.
    if (dloprn .le. 0.) then
        dloprn = prcinf
    end if

    ! Print interval in activity of water.
    if (dlaprn .le. 0.) then
        dlaprn = prcinf
    end if

    ! Print interval in steps.
    if (ksppmx .le. 0) then
        ksppmx = 100
    end if

    ! Plot interval in Xi.
    if (dlxplo .le. 0.) then
        dlxplo = prcinf
    end if

    ! Plot interval in log Xi.
    if (dlxpll .le. 0.) then
        dlxpll = prcinf
    end if

    ! Plot interval in pH.
    if (dlhplo .le. 0.) then
        dlhplo = prcinf
    end if

    ! Plot interval in Eh (v).
    if (dleplo .le. 0.) then
        dleplo = prcinf
    end if

    ! Plot interval in log fO2.
    if (dloplo .le. 0.) then
        dloplo = prcinf
    end if

    ! Plot interval in activity of water.
    if (dlaplo .le. 0.) then
        dlaplo = prcinf
    end if

    ! Plot interval in time (seconds).
    if (dltplo .le. 0.) then
        dltplo = prcinf
    end if

    ! Plot interval in log time.
    if (dltpll .le. 0.) then
        dltpll = prcinf
    end if

    ! Plot interval in steps.
    if (ksplmx .le. 0.) then
        ksplmx = 10000
    end if

    ! PRS transfer interval in Xi.
    if (dlxdmp .le. 0.) then
        if (iopt(1) .le. 1) then
            ! Have the titration model or the closed system model.
            dlxdmp = prcinf
        else if (iopt(1) .eq. 2) then
            ! Have the fluid-centered, flow-through open system model.
            dlxdmp = prcinf

            if (iopt(4).ge.1 .or. net.gt.0) then
                if (ximaxi .gt. 0.) then
                    dlxdmp = ximaxi/25.
                else
                    dlxdmp = 1.e-6
                end if
            end if
        end if
    end if

    if (iopt(1) .eq. 2) then
        if (iopt(4).ge.1 .or. net.gt.0) then
            ! Have solid solutions or generic ion exchangers, or both.
            ! Make sure the PRS transfer interval is no greater then the
            ! linear plot interval.
            if (dlxdmp .gt. dlxplo) then
                dlxdmp = dlxplo
            end if
        end if
    end if
end subroutine dfaltz