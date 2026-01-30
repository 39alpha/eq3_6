subroutine ldlxrc(al10,delxi,dlxmin,dzvc0,iodb,iindx1,kbt,kdim,kelect,khydr,khydx,km1,kmax,ko2gaq,krdxsp,kwater,kxt,nbasp,nbtmax,nodbmx,nord,noutpt,nrd1mx,nstmax,nttyo,qrapch,uspec,zklogu,zvclg0,zvclg1,zvec0,zvec1)
    !! This subroutine limits delxi when a variable corresponding to a
    !! basis species is rapidly changing. Special forms of this
    !! constraint apply to the following species:
    !!   H2O
    !!   H+ (or OH- if that is in the basis set instead of H+)
    !!   O2(g,aq) (or any other auxiliary basis species such as e-,
    !!             HS-, or Fe+++, which might be used as the redox
    !!             defining basis species)
    !! In tracing a reaction path, it is undesirable, for example,
    !! to suddenly run out of solvent water, or to skip over a redox
    !! jump, without obtaining some level of resolution. This subroutine
    !! forces smaller steps to be taken to avoid such problems. Also,
    !! it insures that some level of descriptive detail is obtained
    !! when a signficant, abrupt change is occurring in the system.
    !! It is not the purpose of this subroutine to produce step size
    !! values that correspond exactly to the specified constraints.
    !! The methodology here is rather approximate. Changes in activity
    !! coefficients along the reaction path are not accounted for.
    !! In reducing delxi, the basic approach is based only on the
    !! value at the base point and the first derivative. The full power
    !! of the finite difference representations is therefore not
    !! brought to bear.
    !! This subroutine is called by:
    !!   EQ6/path.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: kmax
    integer :: nbtmax
    integer :: nodbmx
    integer :: nrd1mx
    integer :: nstmax

    integer :: noutpt
    integer :: nttyo

    integer :: iodb(nodbmx)
    integer :: iindx1(kmax)
    integer :: nbasp(nbtmax)

    integer :: kbt
    integer :: kdim
    integer :: kelect
    integer :: khydr
    integer :: khydx
    integer :: km1
    integer :: ko2gaq
    integer :: krdxsp
    integer :: kwater
    integer :: kxt
    integer :: nord

    logical :: qrapch

    character(len=48) :: uspec(nstmax)

    real(kind=8) :: dzvc0(nrd1mx,kmax)
    real(kind=8) :: zvclg0(kmax)
    real(kind=8) :: zvclg1(kmax)
    real(kind=8) :: zvec0(kmax)
    real(kind=8) :: zvec1(kmax)

    real(kind=8) :: al10
    real(kind=8) :: delxi
    real(kind=8) :: dlxmin
    real(kind=8) :: zklogu

    ! Local variable declarations.
    integer :: j2
    integer :: kcol
    integer :: nb
    integer :: ns
    integer :: nstep
    integer :: nstepl

    integer :: ilnobl

    logical :: qztayl

    character(len=24) :: unamsp

    real(kind=8) :: lobasp
    real(kind=8) :: lelect
    real(kind=8) :: lhydr
    real(kind=8) :: lo2gaq
    real(kind=8) :: lwater
    real(kind=8) :: lrdxsp

    real(kind=8) :: dlxic
    real(kind=8) :: dlxsv
    real(kind=8) :: zdif
    real(kind=8) :: zdzid
    real(kind=8) :: zdzspi
    real(kind=8) :: zdzwa
    real(kind=8) :: zdzwai
    real(kind=8) :: zx

    ! Maximum number of cycles cutting the step size in this subroutine.
    data nstepl /20/

    ! Change limit values.
    data lelect /8.0/,lhydr /0.5/,lo2gaq /8.0/,lobasp /8.0/,lrdxsp /8.0/,lwater /0.50/

    qrapch = .false.
    dlxsv = delxi
    unamsp = 'Error'
    nstep = -1

    ! Here is a return point if the value of delxi has been reduced by
    ! this subroutine. This provides an additional check, as the
    ! step-size reduction algorithms are only approximate.
100 continue
    nstep = nstep + 1

    if (nstep .ge. nstepl) then
        if (iodb(1).gt.0 .or. iodb(5).gt.0) then
            write (noutpt,1000) nstepl
            write (nttyo,1000) nstepl
1000 format(/' * Note - (EQ6/ldlxrc) The step size has been cut',' the maximum ',i2,' times',/7x,'to limit rapid change in',' one or more of the variables associated with the',/7x,'aqueous basis species. Cutting the step size to the',' minimum value.')
        end if

        delxi = dlxmin
        go to 990
    end if

    ! Make a Taylor's series expansion of the z vector, without applying
    ! change limits.
    qztayl = .false.
    call ztaylr(delxi,dzvc0,kdim,kmax,km1,kxt,nord,nrd1mx,qztayl,zklogu,zvclg0,zvclg1,zvec0,zvec1)

    if (kwater .gt. 0) then
        ! Number of moles of H2O. The limit is +/- lwater %.
        nb = iindx1(kwater)
        ns = nbasp(nb)
        dlxic = delxi
        zdif = zvec0(kwater) - zvec1(kwater)
        zx = zdif/zvec0(kwater)

        if (abs(zx) .gt. lwater) then
            ! The predicted change is too big. Reduce the step size.
            if (dzvc0(1,kwater) .ne. 0.) then
                ! Estimate the new step size from the finite-difference
                ! data.
                zdzwa = zvec0(kwater)/dzvc0(1,kwater)
                dlxic = lwater*abs(zdzwa)

                if (dlxic .ge. delxi) then
                    dlxic = 0.5*delxi
                end if
            else
                ! Just cut the step size.
                dlxic = 0.5*delxi
            end if

            dlxic = max(dlxic,dlxmin)

            if (dlxic .lt. delxi) then
                delxi = dlxic
                unamsp = uspec(ns)(1:24)
                go to 100
            end if
        end if
    end if

    zdzwai = dzvc0(1,kwater)/zvec0(kwater)

    if (kwater .gt. 0) then
        if (khydr .gt. 0) then
            ! The pH, if H+ is in the basis set. The limit is +/- lhydr
            ! pH unit.
            nb = iindx1(khydr)
            ns = nbasp(nb)
            dlxic = delxi
            zdif = zvclg0(khydr) - zvclg1(khydr)

            if (abs(zdif) .gt. lhydr) then
                ! The predicted change is too big. Reduce the step size.
                zdzspi = dzvc0(1,khydr)/zvec0(khydr)
                zdzid  = zdzspi - zdzwai

                if (zdzid .ne. 0.) then
                    ! Estimate the new step size from the finite-difference
                    ! data.
                    dlxic = al10*lhydr/abs(zdzid)

                    if (dlxic .ge. delxi) then
                        dlxic = 0.5*delxi
                    end if
                else
                    ! Just cut the step size.
                    dlxic = 0.5*delxi
                end if

                dlxic = max(dlxic,dlxmin)

                if (dlxic .lt. delxi) then
                    delxi = dlxic
                    unamsp = uspec(ns)(1:24)
                    go to 100
                end if
            end if
        else if (khydx .gt. 0) then
            ! The pH, if OH- is in the basis set instead of H+.
            ! The limit is +/- lhydr pH unit.
            nb = iindx1(khydx)
            ns = nbasp(nb)
            dlxic = delxi
            zdif = zvclg0(khydx) - zvclg1(khydx)

            if (abs(zdif) .gt. lhydr) then
                ! The predicted change is too big. Reduce the step size.
                zdzspi = dzvc0(1,khydx)/zvec0(khydx)
                zdzid  = zdzspi - zdzwai

                if (zdzid .ne. 0.) then
                    ! Estimate the new step size from the finite-difference
                    ! data.
                    dlxic = al10*lhydr/abs(zdzid)

                    if (dlxic .ge. delxi) then
                        dlxic = 0.5*delxi
                    end if
                else
                    ! Just cut the step size.
                    dlxic = 0.5*delxi
                end if

                dlxic = max(dlxic,dlxmin)

                if (dlxic .lt. delxi) then
                    delxi = dlxic
                    unamsp = uspec(ns)(1:24)
                    go to 100
                end if
            end if
        end if
    end if

    if (ko2gaq .gt. 0) then
        ! The log fO2. The limit is +/- lo2gaq log units.
        nb = iindx1(ko2gaq)
        ns = nbasp(nb)
        dlxic = delxi
        zdif = zvclg0(ko2gaq) - zvclg1(ko2gaq)

        if (abs(zdif) .gt. lo2gaq) then
            if (dzvc0(1,ko2gaq) .ne. 0.) then
                ! Estimate the new step size from the finite-difference
                ! data.
                dlxic = al10*lo2gaq*abs(zvec0(ko2gaq)/dzvc0(1,ko2gaq))

                if (dlxic .ge. delxi) then
                    dlxic = 0.5*delxi
                end if
            else
                ! Just cut the step size.
                dlxic = 0.5*delxi
            end if

            dlxic = max(dlxic,dlxmin)

            if (dlxic .lt. delxi) then
                delxi = dlxic
                unamsp = uspec(ns)(1:24)
                go to 100
            end if
        end if
    end if

    if (kelect .gt. 0) then
        ! The pe. The limit is +/- lelect log units.
        nb = iindx1(kelect)
        ns = nbasp(nb)
        dlxic = delxi
        zdif = zvclg0(kelect) - zvclg1(kelect)

        if (abs(zdif) .gt. lelect) then
            if (dzvc0(1,kelect) .ne. 0.) then
                ! Estimate the new step size from the finite-difference
                ! data.
                dlxic = al10*lelect*abs(zvec0(kelect)/dzvc0(1,kelect))

                if (dlxic .ge. delxi) then
                    dlxic = 0.5*delxi
                end if
            else
                ! Just cut the step size.
                dlxic = 0.5*delxi
            end if

            dlxic = max(dlxic,dlxmin)

            if (dlxic .lt. delxi) then
                delxi = dlxic
                unamsp = uspec(ns)(1:24)
                go to 100
            end if
        end if
    end if

    if (krdxsp.gt.0 .and. kwater.gt.0) then
        if (krdxsp.ne.ko2gaq .and. krdxsp.ne.kelect) then
            ! The log molality of an auxiliary basis species. The limit is
            ! +/- lrdxsp log units.
            nb = iindx1(krdxsp)
            ns = nbasp(nb)
            dlxic = delxi
            zdif = zvclg0(krdxsp) - zvclg1(krdxsp)

            if (abs(zdif) .gt. lrdxsp) then
                zdzspi = dzvc0(1,krdxsp)/zvec0(krdxsp)
                zdzid  = zdzspi - zdzwai

                if (zdzid .ne. 0.) then
                    ! Estimate the new step size from the finite-difference
                    ! data.
                    dlxic = al10*lrdxsp/abs(zdzid)

                    if (dlxic .ge. delxi) then
                        dlxic = 0.5*delxi
                    end if
                else
                    ! Just cut the step size.
                    dlxic = 0.5*delxi
                end if

                dlxic = max(dlxic,dlxmin)

                if (dlxic .lt. delxi) then
                    if (dlxic .gt. dlxmin) then
                        delxi = dlxic
                        unamsp = uspec(ns)(1:24)
                        go to 100
                    end if
                end if
            end if
        end if
    end if

    do kcol = 1,kbt
        if (kcol.ne.kwater .and. kcol.ne.khydr .and. kcol.ne.khydx    .and. kcol.ne.krdxsp) then
            ! The log molality of any other basis species. The limit is
            ! +/- lobasp log units.
            nb = iindx1(kcol)
            ns = nbasp(nb)
            dlxic = delxi
            zdif = zvclg0(kcol) - zvclg1(kcol)

            if (abs(zdif) .gt. lobasp) then
                zdif = zvclg0(kcol) - zvclg1(kcol)

                if (abs(zdif) .gt. lobasp) then
                    zdzspi = dzvc0(1,kcol)/zvec0(kcol)
                    zdzid  = zdzspi - zdzwai

                    if (zdzid .ne. 0.) then
                        ! Estimate the new step size from the finite-difference
                        ! data.
                        dlxic = al10*lobasp/abs(zdzid)

                        if (dlxic .ge. delxi) then
                            dlxic = 0.5*delxi
                        end if
                    else
                        ! Just cut the step size.
                        dlxic = 0.5*delxi
                    end if

                    dlxic = max(dlxic,dlxmin)

                    if (dlxic .lt. delxi) then
                        delxi = dlxic
                        unamsp = uspec(ns)(1:24)
                        go to 100
                    end if
                end if
            end if
        end if
    end do

990 continue
    qrapch = delxi .lt. dlxsv

    if (qrapch) then
        if (iodb(5) .gt. 0) then
            j2 = ilnobl(unamsp)
            write (noutpt,1010) dlxsv,delxi,unamsp(1:j2)
            write (nttyo,1010) dlxsv,delxi,unamsp(1:j2)
1010 format(/' * Note - (EQ6/ldlxrc) The step size has been cut',' from ',1pe12.5,/7x,'to ',e12.5,' to limit rapid change',' associated with ',a,'.')
        end if
    end if
end subroutine ldlxrc
