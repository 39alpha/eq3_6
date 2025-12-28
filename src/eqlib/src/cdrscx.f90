subroutine cdrscx(adhfs,adhfsx,advfs,advfsx,axhfs,axhfsx,axlks,axlksx,axvfs,axvfsx,cdrs,cdrsx,ipch,ipchmx,ipcv,ipcvmx,narxmx,ndrs,ndrsmx,ndrsx,ndrsr,ndrsrx,nstmax,ntprmx)
    !! This subroutine copies the scratch arrays for reactions and
    !! reaction properties into the standard arrays.
    !! This subroutine is called by:
    !!   EQLIB/elim.f
    !!   EQLIB/mdrgex.f
    !!   EQLIB/switch.f
    !!   EQLIB/swtchb.f
    !! Principal input:
    !!   adhfsx = scratch array of coefficients for computing
    !!              pressure derivatives of enthalpy functions as a
    !!              function of temperature
    !!   advfsx = scratch array of coefficients for computing
    !!              pressure derivatives of volume functions as a
    !!              function of temperature
    !!   axhfsx = scratch array of coefficients for computing
    !!              enthalpy functions as a function of temperature
    !!   axlksx = scratch array of coefficients for computing
    !!              equilibrium constants as a function of temperature
    !!   axvfsx = scratch array of coefficients for computing
    !!              volume functions as a function of temperature
    !!   cdrsx  = scratch array of reaction coefficients
    !!   narxmx = number of coefficient elements of axlks per species
    !!              per temperature range
    !!   ndrsx  = scratch array of indices of species corresponding to
    !!              reaction coefficients
    !!   ndrsrx = scratch pointer array for indices of species appearing
    !!              in reactions
    !!   nstmax = maximum number of species
    !!   ntprmx = maximum number of temperature ranges
    !! Principal output:
    !!   adhfs  = standard array of coefficients for computing
    !!              pressure derivatives of enthalpy functions as a
    !!              function of temperature
    !!   advfs  = standard array of coefficients for computing
    !!              pressure derivatives of volume functions as a
    !!              function of temperature
    !!   axhfs  = standard array of coefficients for computing
    !!              enthalpy functions as a function of temperature
    !!   axlks  = standard array of coefficients for computing
    !!              equilibrium constants as a function of temperature
    !!   axvfs  = standard array of coefficients for computing
    !!              volume functions as a function of temperature
    !!   cdrs   = standard array of reaction coefficients
    !!   ndrs   = standard  array of indices of species corresponding to
    !!              reaction coefficients
    !!   ndrsr  = standard pointer array for indices of species
    !!              appearing in reactions
    implicit none

    ! Calling sequence variable declarations.
    integer :: ipchmx
    integer :: ipcvmx
    integer :: narxmx
    integer :: ndrsmx
    integer :: nstmax
    integer :: ntprmx

    integer :: ndrs(ndrsmx)
    integer :: ndrsx(ndrsmx)
    integer :: ndrsr(2,nstmax)
    integer :: ndrsrx(2,nstmax)
    integer :: ipch
    integer :: ipcv

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

    ! Local variable declarations.
    integer :: jhmax
    integer :: jvmax
    integer :: nmax

    nmax = narxmx*ntprmx*nstmax
    jhmax = narxmx*ntprmx*ipchmx*nstmax
    jvmax = narxmx*ntprmx*ipcvmx*nstmax

    call copyaa(axlksx,axlks,nmax)

    if (ipch .ge. 0) then
        call copyaa(axhfsx,axhfs,nmax)

        if (ipch .ge. 1) then
            call copyaa(adhfsx,adhfs,jhmax)
        end if
    end if

    if (ipcv .ge. 0) then
        call copyaa(axvfsx,axvfs,nmax)

        if (ipcv .ge. 1) then
            call copyaa(advfsx,advfs,jvmax)
        end if
    end if

    nmax = 2*nstmax
    call copyia(ndrsrx,ndrsr,nmax)

    call copyaa(cdrsx,cdrs,ndrsmx)
    call copyia(ndrsx,ndrs,ndrsmx)
end subroutine cdrscx