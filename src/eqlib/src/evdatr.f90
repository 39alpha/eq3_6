subroutine evdatr(adhfs,advfs,axhfs,axlks,axvfs,dhfs,dvfs,ipch,ipchmx,ipcv,ipcvmx,narxmx,narxt,nst,nstmax,ntpr,ntprmx,tempc,xhfs,xlks,xvfs)
    !! This subroutine evaluates equilibrium constants and related
    !! thermodynamic functions as functions of temperature and
    !! standard grid pressure. To make pressure corrections off the
    !! standard T-P grid for these data, use EQLIB/pcorrx.f.
    !! This subroutine is called by:
    !!   EQLIB/absswa.f
    !!   EQLIB/evdata.f
    !!   EQ3NR/eq3nr.f
    !!   EQ6/absswb.f
    !!   EQ6/eq6.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: ipchmx
    integer :: ipcvmx
    integer :: narxmx
    integer :: nstmax
    integer :: ntprmx

    integer :: narxt(ntprmx)

    integer :: ipch
    integer :: ipcv
    integer :: nst
    integer :: ntpr

    real(kind=8) :: adhfs(narxmx,ntprmx,ipchmx,nstmax)
    real(kind=8) :: advfs(narxmx,ntprmx,ipcvmx,nstmax)
    real(kind=8) :: axhfs(narxmx,ntprmx,nstmax)
    real(kind=8) :: axlks(narxmx,ntprmx,nstmax)
    real(kind=8) :: axvfs(narxmx,ntprmx,nstmax)
    real(kind=8) :: dhfs(ipchmx,nstmax)
    real(kind=8) :: dvfs(ipcvmx,nstmax)
    real(kind=8) :: xhfs(nstmax)
    real(kind=8) :: xlks(nstmax)
    real(kind=8) :: xvfs(nstmax)

    real(kind=8) :: tempc

    ! Local variable declarations.
    integer :: ipc
    integer :: ns

    real(kind=8) :: prop

    ! Compute the log K values for all reactions.
    do ns = 1,nst
        ! Compute the log K values on the standard grid.
        if (axlks(1,ntpr,ns) .lt. 9999999.) then
            ! Calling sequence substitutions:
            !   axlks for arr
            !   ns for k
            !   nstmax for nmax
            call evdat3(axlks,ns,nstmax,narxmx,narxt,ntpr,ntprmx,prop,tempc)
            xlks(ns) = prop
        else
            xlks(ns) = 9999999.
        end if

        if (ipch .ge. 0) then
            ! Compute the enthalpy function values on the standard grid.
            if (axhfs(1,ntpr,ns) .lt. 9999999.) then
                ! Calling sequence substitutions:
                !   axhfs for arr
                !   ns for k
                !   nstmax for nmax
                call evdat3(axhfs,ns,nstmax,narxmx,narxt,ntpr,ntprmx,prop,tempc)
                xhfs(ns) = prop
            else
                xhfs(ns) = 9999999.
            end if

            do ipc = 1,ipch
                ! Compute the pressure derivatives of the enthalpy function
                ! values on the standard grid.
                ! Calling sequence substitutions:
                !   adhfs for arr
                !   ipchmx for ipcxmx
                !   ns for k
                !   nstmax for nmax
                call evdat4(adhfs,ipc,ipchmx,ns,nstmax,narxmx,narxt,ntpr,ntprmx,prop,tempc)
                dhfs(ipc,ns) = prop
            end do
        end if

        if (ipcv .ge. 0) then
            ! Compute the volume function values on the standard grid.
            if (axvfs(1,ntpr,ns) .lt. 9999999.) then
                ! Calling sequence substitutions:
                !   axvfs for arr
                !   ns for k
                !   nstmax for nmax
                call evdat3(axvfs,ns,nstmax,narxmx,narxt,ntpr,ntprmx,prop,tempc)
                xvfs(ns) = prop
            else
                xvfs(ns) = 9999999.
            end if

            do ipc = 1,ipcv
                ! Compute the pressure derivatives of the volume function
                ! values on the standard grid.
                ! Calling sequence substitutions:
                !   advfs for arr
                !   ipcvmx for ipcxmx
                !   ns for k
                !   nstmax for nmax
                call evdat4(advfs,ipc,ipcvmx,ns,nstmax,narxmx,narxt,ntpr,ntprmx,prop,tempc)
                dvfs(ipc,ns) = prop
            end do
        end if
    end do
end subroutine evdatr