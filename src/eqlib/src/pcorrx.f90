subroutine pcorrx(avcnst,dhfs,dvfs,ipch,ipchmx,ipcv,ipcvmx,nbasp,nbt,nbtmax,ndrsr,nst,nstmax,presg,press,xhfs,xlks,xvfs)
    !! This subroutine makes pressure corrections for equilibrium
    !! constants and related thermodynamic functions. It normally
    !! corrects for pressures off the standard T-P grid (e.g., from
    !! "presg" to "press"). However, it can be used to correct from
    !! any pressure to another by substituting the value of the former
    !! for "presg". This is done in EQ6 to correct for changing pressure
    !! along an isothermal reaction path. In order to do this, the
    !! pressure derivatives of enthalpy and volume functions (dhfs,dvfs)
    !! are also corrected to the pressure of interest. The derivatives of
    !! highest order are necessarily treated as constants, so there
    !! is no correction in this case.
    !! EQLIB/pcorrm.f performs the same function for some miscellaneous
    !! thermodynamic functions, such as Debye-Huckel
    !! parameters.
    !! This subroutine is called by:
    !!   EQ3NR/arrset.f
    !!   EQ3NR/eq3nr.f
    !!   EQ6/absswb.f
    !!   EQ6/eq6.f
    !!   EQ6/tpadv.f
    !! Principal input:
    !!  avcnst = 2.303 RT, with R in units of bar-cm3/mol-K
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: ipchmx
    integer :: ipcvmx
    integer :: nbtmax
    integer :: nstmax

    integer :: ipch
    integer :: ipcv
    integer :: nbt
    integer :: nst

    integer :: nbasp(nbtmax)
    integer :: ndrsr(2,nstmax)

    real(kind=8) :: dhfs(ipchmx,nstmax)
    real(kind=8) :: dvfs(ipcvmx,nstmax)
    real(kind=8) :: xhfs(nstmax)
    real(kind=8) :: xlks(nstmax)
    real(kind=8) :: xvfs(nstmax)

    real(kind=8) :: avcnst
    real(kind=8) :: presg
    real(kind=8) :: press

    ! Local variable declarations.
    integer :: i
    integer :: ipc
    integer :: n
    integer :: nb
    integer :: ns
    integer :: nt

    real(kind=8) :: dp
    real(kind=8) :: xx

    real(kind=8) :: fctrl

    ! Calculate the pressure difference.
    dp = press - presg

    if (ipcv .eq. 0) then
        ! Correct the log K values to the current pressure, using
        ! a first order (constant volume of reaction) correction.
        ! This is expected to be the most common order for pressure
        ! corrections, hence this special block has been written
        ! for the sake of efficiency.
        do ns = 1,nst
            if (xlks(ns).gt.-9999999. .and. xlks(ns).lt.9999999.) then
                xx = -xvfs(ns)*dp/avcnst
                xlks(ns) = xlks(ns) + xx
            end if
        end do
    end if

    if (ipcv .gt. 0) then
        ! Correct the log K values to the current pressure, using
        ! a correction of order 2 or higher. The coding in this block
        ! could be used to handle the order 1 case. However, that
        ! would not be quite as efficient as the above block.
        do ns = 1,nst
            if (xlks(ns).gt.-9999999. .and. xlks(ns).lt.9999999.) then
                xx = -xvfs(ns)*dp

                do ipc = 1,ipcv
                    n = ipc + 1
                    xx = xx - ( dvfs(ipc,ns)*(dp**n) )/fctrl(n)
                end do

                xx = xx/avcnst
                xlks(ns) = xlks(ns) + xx
            end if
        end do
    end if

    if (ipcv .ge. 0) then
        ! Uncorrect the log K values for strict basis species,
        ! as these are fixed at zero and the volume function array
        ! contains the standard partial molar volume, not the
        ! standard partial molar volume of reaction.
        do nb = 1,nbt
            ns = nbasp(nb)
            nt = ndrsr(2,ns) - ndrsr(1,ns) + 1

            if (nt .lt. 2) then
                xlks(ns) = 0.
            end if
        end do
    end if

    if (ipch .gt. 0) then
        ! Correct the enthalpy function values to the current pressure.
        do ns = 1,nst
            xx = 0.

            do ipc = 1,ipch
                xx = xx  + ( dhfs(ipc,ns)*(dp**ipc) )/fctrl(ipc)
            end do

            xhfs(ns) = xhfs(ns) + xx

            do ipc = 1,ipch - 1
                ! Correct the pressure derivatives of the enthalpy functions
                ! to the current pressure.
                xx = 0.

                do i = ipc + 1,ipch
                    n = i - ipc
                    xx = xx  + ( dhfs(i,ns)*(dp**n) )/fctrl(n)
                end do

                dhfs(ipc,ns) = dhfs(ipc,ns) + xx
            end do
        end do
    end if

    if (ipcv .gt. 0) then
        ! Correct the volume function values to the current pressure.
        do ns = 1,nst
            xx = 0.

            do ipc = 1,ipcv
                xx = xx  + ( dvfs(ipc,ns)*(dp**ipc) )/fctrl(ipc)
            end do

            xvfs(ns) = xvfs(ns) + xx

            do ipc = 1,ipcv - 1
                ! Correct the pressure derivatives of the volume functions
                ! to the current pressure.
                xx = 0.

                do i = ipc + 1,ipcv
                    n = i - ipc
                    xx = xx  + ( dvfs(i,ns)*(dp**n) )/fctrl(n)
                end do

                dvfs(ipc,ns) = dvfs(ipc,ns) + xx
            end do
        end do
    end if
end subroutine pcorrx