subroutine cdrssd(adhfs,adhfsd,advfs,advfsd,axhfs,axhfsd,axlks,axlksd,axvfs,axvfsd,cdrs,cdrsd,ipch,ipchmx,ipcv,ipcvmx,narxmx,nbasp,nbaspd,nbtmax,ndrs,ndrsd,ndrsmx,ndrsr,ndrsrd,nstmax,ntprmx)
    !! This subroutine copies the reactions and reaction properties as
    !! they are currently written into the 'd' set.
    !! This subroutine is called by:
    !!   EQ3NR/eq3nr.f
    !!   EQ6/eq6.f
    !! Principal input:
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
    !!   narxmx = number of coefficient elements of axlks per species
    !!              per temperature range
    !!   nbtmax = maximum number of basis species
    !!   ndrs   = standard  array of indices of species corresponding to
    !!              reaction coefficients
    !!   ndrsr  = standard pointer array for indices of species appearing
    !!              in reactions
    !!   nstmax = maximum number of species
    !!   ntprmx = maximum number of temperature ranges
    !! Principal output:
    !!   adhfsd = 'data file' array of coefficients for computing
    !!              pressure derivatives of enthalpy functions as a
    !!              function of temperature
    !!   advfsd = 'data file' array of coefficients for computing
    !!              pressure derivatives of volume functions as a
    !!              function of temperature
    !!   axhfsd = 'data file' array of coefficients for computing
    !!              enthalpy functions as a function of temperature
    !!   axlksd = 'data file' array of coefficients for computing
    !!              equilibrium constants as a function of temperature
    !!   axvfsd = 'data file' array of coefficients for computing
    !!              volume functions as a function of temperature
    !!   cdrsd  = 'data file' array of reaction coefficients
    !!   ndrsd  = 'data file'  array of indices of species corresponding
    !!              to reaction coefficients
    !!   ndrsrd = 'data file' pointer array for indices of species
    !!              appearing in reactions
    implicit none

    ! Calling sequence variable declarations.
    integer :: ipchmx
    integer :: ipcvmx
    integer :: narxmx
    integer :: nbtmax
    integer :: ndrsmx
    integer :: nstmax
    integer :: ntprmx

    integer :: nbasp(nbtmax)
    integer :: nbaspd(nbtmax)
    integer :: ndrs(ndrsmx)
    integer :: ndrsd(ndrsmx)
    integer :: ndrsr(2,nstmax)
    integer :: ndrsrd(2,nstmax)
    integer :: ipch
    integer :: ipcv

    real(kind=8) :: adhfs(narxmx,ntprmx,ipchmx,nstmax)
    real(kind=8) :: adhfsd(narxmx,ntprmx,ipchmx,nstmax)
    real(kind=8) :: advfs(narxmx,ntprmx,ipcvmx,nstmax)
    real(kind=8) :: advfsd(narxmx,ntprmx,ipcvmx,nstmax)
    real(kind=8) :: axhfs(narxmx,ntprmx,nstmax)
    real(kind=8) :: axhfsd(narxmx,ntprmx,nstmax)
    real(kind=8) :: axlks(narxmx,ntprmx,nstmax)
    real(kind=8) :: axlksd(narxmx,ntprmx,nstmax)
    real(kind=8) :: axvfs(narxmx,ntprmx,nstmax)
    real(kind=8) :: axvfsd(narxmx,ntprmx,nstmax)
    real(kind=8) :: cdrs(ndrsmx)
    real(kind=8) :: cdrsd(ndrsmx)

    ! Local variable declarations.
    integer :: jhmax
    integer :: jvmax
    integer :: nmax

    nmax = narxmx*ntprmx*nstmax
    jhmax = narxmx*ntprmx*ipchmx*nstmax
    jvmax = narxmx*ntprmx*ipcvmx*nstmax

    call copyaa(axlks,axlksd,nmax)

    if (ipch .ge. 0) then
        call copyaa(axhfs,axhfsd,nmax)

        if (ipch .ge. 1) then
            call copyaa(adhfs,adhfsd,jhmax)
        end if
    end if

    if (ipcv .ge. 0) then
        call copyaa(axvfs,axvfsd,nmax)

        if (ipcv .ge. 1) then
            call copyaa(advfs,advfsd,jvmax)
        end if
    end if

    nmax = 2*nstmax
    call copyia(ndrsr,ndrsrd,nmax)

    call copyaa(cdrs,cdrsd,ndrsmx)
    call copyia(ndrs,ndrsd,ndrsmx)

    call copyia(nbasp,nbaspd,nbtmax)
end subroutine cdrssd