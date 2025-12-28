subroutine autosw(adhfs,adhfsx,advfs,advfsx,axhfs,axhfsx,axlks,axlksx,axvfs,axvfsx,cdrs,cdrsx,eps100,ibswx,iindx1,ipch,ipchmx,ipcv,ipcvmx,jflag,jsflag,kbt,kmax,narn1,narxmx,nbasp,nbaspd,nbaspx,nbt,nbtmax,nbw,ndrs,ndrsmx,ndrsx,ndrsr,ndrsrx,noutpt,nst,nstmax,ntprmx,nttyo,qbassw,uspec,uzvec1)
    !! This subroutine executes automatic basis switching for the purpose
    !! of reducing mass balance residuals. EQLIB/fbassw.f finds
    !! candidates for basis switching, and EQLIB/gabswx.f resolves any
    !! conflicts.
    !! This subroutine is called by:
    !!   EQLIB/absswa.f
    !!   EQ6/absswb.f
    !! Principal input:
    !!   axlks  = array of polynomial coefficients for computing log K
    !!              values (altered by this subroutine)
    !!   ibswx  = array defining switches to be made
    !!   nbasp  = array defining the species in the basis set
    !!              (altered by this subroutine)
    !!   nbaspd = array defining the species in the data file basis set
    !! Principal output:
    !!   axlks  = array of polynomial coefficients for computing log K
    !!              values
    !!   cdrs   = array of reaction coefficients
    !!   nbasp  = array defining the species in the basis set
    !!   ndrs   = array of species indices corresponding to the reaction
    !!              coefficients in the cdrs array
    !!   ndrsr  = array defining the range in the cdrs and ndrs arrays
    !!              corresponding to the reaction for a given species
    implicit none

    ! Calling sequence variable declarations.
    integer :: ipchmx
    integer :: ipcvmx
    integer :: kmax
    integer :: narxmx
    integer :: nbtmax
    integer :: ndrsmx
    integer :: nstmax
    integer :: ntprmx

    integer :: noutpt
    integer :: nttyo

    integer :: ibswx(nbtmax)
    integer :: iindx1(kmax)
    integer :: jflag(nstmax)
    integer :: jsflag(nstmax)
    integer :: nbasp(nbtmax)
    integer :: nbaspd(nbtmax)
    integer :: nbaspx(nbtmax)
    integer :: ndrs(ndrsmx)
    integer :: ndrsx(ndrsmx)
    integer :: ndrsr(2,nstmax)
    integer :: ndrsrx(2,nstmax)

    integer :: ipch
    integer :: ipcv
    integer :: kbt
    integer :: narn1
    integer :: nbt
    integer :: nbw
    integer :: nst

    logical :: qbassw

    character(len=48) :: uspec(nstmax)
    character(len=48) :: uzvec1(kmax)

    real(kind=8) :: adhfs(narxmx,ntprmx,ipchmx,nstmax)
    real(kind=8) :: adhfsx(narxmx,ntprmx,ipchmx,nstmax)
    real(kind=8) :: advfs(narxmx,ntprmx,ipcvmx,nstmax)
    real(kind=8) :: advfsx(narxmx,ntprmx,ipcvmx,nstmax)
    real(kind=8) :: axhfs(narxmx,ntprmx,nstmax)
    real(kind=8) :: axhfsx(narxmx,ntprmx,nstmax)
    real(kind=8) :: axlks(narxmx,ntprmx,nstmax)
    real(kind=8) :: axlksx(narxmx,ntprmx,nstmax)
    real(kind=8) :: axvfs(narxmx,ntprmx,nstmax)
    real(kind=8) :: axvfsx(narxmx,ntprmx,nstmax)
    real(kind=8) :: cdrs(ndrsmx)
    real(kind=8) :: cdrsx(ndrsmx)

    real(kind=8) :: eps100

    ! Local variable declarations.
    integer :: nb
    integer :: nsd
    integer :: ns1
    integer :: ns2
    integer :: krow

    logical :: qbswok

    ! Save the nbasp array.
    call copyia(nbasp,nbaspx,nbt)

    do krow = 1,kbt
        nb = iindx1(krow)
        ns2 = ibswx(nb)

        if (ns2 .gt. 0) then
            ns1 = nbaspx(nb)
            nsd = nbaspd(nb)

            ! If the ns1-th species has not previously been switched with
            ! the nsd-th species, go make the currently requested switch.
            ! If it has been switched, first undo that switch. Take care
            ! not to undo it twice.
            if (nsd .ne. ns2) then
                if (ns1 .ne. nsd) then
                    ! First undo the existing switch.
                    nbasp(nb) = nsd

                    ! Calling sequence substitutions:
                    !   nsd for ns2
                    call switch(adhfs,adhfsx,advfs,advfsx,axhfs,axhfsx,axlks,axlksx,axvfs,axvfsx,cdrs,cdrsx,eps100,ipch,ipchmx,ipcv,ipcvmx,jflag,jsflag,narn1,narxmx,nbasp,nbaspd,nbaspx,nb,nbt,nbtmax,nbw,ndrs,ndrsmx,ndrsx,ndrsr,ndrsrx,noutpt,nsd,nst,nstmax,ntprmx,nttyo,qbassw,qbswok,uspec)
                end if
            end if

            nbasp(nb) = ns2

            call switch(adhfs,adhfsx,advfs,advfsx,axhfs,axhfsx,axlks,axlksx,axvfs,axvfsx,cdrs,cdrsx,eps100,ipch,ipchmx,ipcv,ipcvmx,jflag,jsflag,narn1,narxmx,nbasp,nbaspd,nbaspx,nb,nbt,nbtmax,nbw,ndrs,ndrsmx,ndrsx,ndrsr,ndrsrx,noutpt,ns2,nst,nstmax,ntprmx,nttyo,qbassw,qbswok,uspec)

            ! Update the names in the uzvec1 array.
            uzvec1(krow) = uspec(ns2)
        end if
    end do
end subroutine autosw