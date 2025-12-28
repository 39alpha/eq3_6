subroutine mdrgex(adhfs,adhfsd,adhfsx,advfs,advfsd,advfsx,axhfs,axhfsd,axhfsx,axlks,axlksd,axlksx,axvfs,axvfsd,axvfsx,cdrs,cdrsd,cdrsx,ipch,ipchmx,ipcv,ipcvmx,narxmx,nbasp,nbaspd,nbtmax,ndrs,ndrsd,ndrsx,ndrsmx,ndrsr,ndrsrd,ndrsrx,nern1,nern2,noutpt,nst,nstmax,ntprmx,nttyo)
    !! This subroutine folds the reactions and reaction properties for
    !! generic ion exchangers into the 'd' set. This is generally done
    !! after those reactions and properties have been manipulated.
    !! Note that the modified arrays must first be written into a set
    !! of scratch arrays. This is because the ion exchanger sections
    !! must generally go in the middle of the new arrays, and their
    !! lengths are generally such that the positioning of the start of
    !! the terminal "unmodified" sections is modified.
    !! This subroutine is called by:
    !!   EQLIB/chsgex.f
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

    integer :: noutpt
    integer :: nttyo

    integer :: nbasp(nbtmax)
    integer :: nbaspd(nbtmax)
    integer :: ndrs(ndrsmx)
    integer :: ndrsd(ndrsmx)
    integer :: ndrsx(ndrsmx)
    integer :: ndrsr(2,nstmax)
    integer :: ndrsrd(2,nstmax)
    integer :: ndrsrx(2,nstmax)

    integer :: ipch
    integer :: ipcv
    integer :: nern1
    integer :: nern2
    integer :: nst

    real(kind=8) :: adhfs(narxmx,ntprmx,ipchmx,nstmax)
    real(kind=8) :: adhfsd(narxmx,ntprmx,ipchmx,nstmax)
    real(kind=8) :: adhfsx(narxmx,ntprmx,ipchmx,nstmax)
    real(kind=8) :: advfs(narxmx,ntprmx,ipcvmx,nstmax)
    real(kind=8) :: advfsd(narxmx,ntprmx,ipcvmx,nstmax)
    real(kind=8) :: advfsx(narxmx,ntprmx,ipcvmx,nstmax)
    real(kind=8) :: axhfs(narxmx,ntprmx,nstmax)
    real(kind=8) :: axhfsd(narxmx,ntprmx,nstmax)
    real(kind=8) :: axhfsx(narxmx,ntprmx,nstmax)
    real(kind=8) :: axlks(narxmx,ntprmx,nstmax)
    real(kind=8) :: axlksd(narxmx,ntprmx,nstmax)
    real(kind=8) :: axlksx(narxmx,ntprmx,nstmax)
    real(kind=8) :: axvfs(narxmx,ntprmx,nstmax)
    real(kind=8) :: axvfsd(narxmx,ntprmx,nstmax)
    real(kind=8) :: axvfsx(narxmx,ntprmx,nstmax)
    real(kind=8) :: cdrs(ndrsmx)
    real(kind=8) :: cdrsd(ndrsmx)
    real(kind=8) :: cdrsx(ndrsmx)

    ! Local variable declarations.
    integer :: i
    integer :: ipc
    integer :: j
    integer :: jhmax
    integer :: jvmax
    integer :: j2
    integer :: j3
    integer :: n
    integer :: nmax
    integer :: nb
    integer :: nn
    integer :: nr1
    integer :: nr2
    integer :: ns
    integer :: nt

    integer :: ilnobl

    character(len=8) :: ux8a
    character(len=8) :: ux8b

    nmax = narxmx*ntprmx*(nern1 -1)
    jhmax = narxmx*ntprmx*ipchmx*(nern1 -1)
    jvmax = narxmx*ntprmx*ipcvmx*(nern1 -1)

    ! Coefficients for equilibrium constants.
    if (nern1 .gt. 1) then
        call copyaa(axlksd,axlksx,nmax)
    end if

    do ns = nern1,nern2
        do j = 1,ntprmx
            do i = 1,narxmx
                axlksx(i,j,ns) = axlks(i,j,ns)
            end do
        end do
    end do

    if (nern2 .lt. nst) then
        do ns = nern2 + 1,nst
            do j = 1,ntprmx
                do i = 1,narxmx
                    axlksx(i,j,ns) = axlksd(i,j,ns)
                end do
            end do
        end do
    end if

    ! Coefficients related to enthalpy.
    if (ipch .ge. 0) then
        ! Coefficients for enthalpy functions.
        if (nern1 .gt. 1) then
            call copyaa(axhfsd,axhfsx,nmax)
        end if

        do ns = nern1,nern2
            do j = 1,ntprmx
                do i = 1,narxmx
                    axhfsx(i,j,ns) = axhfs(i,j,ns)
                end do
            end do
        end do

        if (nern2 .lt. nst) then
            do ns = nern2 + 1,nst
                do j = 1,ntprmx
                    do i = 1,narxmx
                        axhfsx(i,j,ns) = axhfsd(i,j,ns)
                    end do
                end do
            end do
        end if

        if (ipch .ge. 1) then
            ! Coefficients for derivatives of enthalpy functions.
            if (nern1 .gt. 1) then
                call copyaa(adhfsd,adhfsx,jhmax)
            end if

            do ns = nern1,nern2
                do ipc = 1,ipchmx
                    do j = 1,ntprmx
                        do i = 1,narxmx
                            adhfsx(i,j,ipc,ns) = adhfs(i,j,ipc,ns)
                        end do
                    end do
                end do
            end do

            if (nern2 .lt. nst) then
                do ns = nern2 + 1,nst
                    do ipc = 1,ipchmx
                        do j = 1,ntprmx
                            do i = 1,narxmx
                                adhfsx(i,j,ipc,ns) = adhfsd(i,j,ipc,ns)
                            end do
                        end do
                    end do
                end do
            end if
        end if
    end if

    ! Coefficients related to volume.
    if (ipcv .ge. 0) then
        ! Coefficients for volume functions.
        if (nern1 .gt. 1) then
            call copyaa(axvfsd,axvfsx,nmax)
        end if

        do ns = nern1,nern2
            do j = 1,ntprmx
                do i = 1,narxmx
                    axvfsx(i,j,ns) = axvfs(i,j,ns)
                end do
            end do
        end do

        if (nern2 .lt. nst) then
            do ns = nern2 + 1,nst
                do j = 1,ntprmx
                    do i = 1,narxmx
                        axvfsx(i,j,ns) = axvfsd(i,j,ns)
                    end do
                end do
            end do
        end if

        if (ipcv .ge. 1) then
            ! Coefficients for derivatives of volume functions.
            if (nern1 .gt. 1) then
                call copyaa(advfsd,advfsx,jvmax)
            end if

            do ns = nern1,nern2
                do ipc = 1,ipcvmx
                    do j = 1,ntprmx
                        do i = 1,narxmx
                            advfsx(i,j,ipc,ns) = advfs(i,j,ipc,ns)
                        end do
                    end do
                end do
            end do

            if (nern2 .lt. nst) then
                do ns = nern2 + 1,nst
                    do ipc = 1,ipcvmx
                        do j = 1,ntprmx
                            do i = 1,narxmx
                                advfsx(i,j,ipc,ns) = advfsd(i,j,ipc,ns)
                            end do
                        end do
                    end do
                end do
            end if
        end if
    end if

    ! Index range pointers for reactions.
    if (nern1 .gt. 1) then
        nmax = 2*(nern1 - 1)
        call copyia(ndrsrd,ndrsrx,nmax)
    end if

    do ns = nern1,nern2
        ndrsrx(1,ns) = ndrsr(1,ns)
        ndrsrx(1,ns) = ndrsr(1,ns)
    end do

    if (nern2 .lt. nst) then
        do ns = nern2 + 1,nst
            ndrsrx(1,ns) = ndrsrd(1,ns)
            ndrsrx(1,ns) = ndrsrd(1,ns)
        end do
    end if

    ! Reaction coefficients and indices of associated species.
    nt = 0

    if (nern1 .gt. 1) then
        nt = ndrsrd(2,nern1 - 1)
    end if

    if (nern2 .ge. nern1) then
        nt = nt + ndrsr(2,nern2) - ndrsr(1,nern1) + 1
    end if

    if (nst .gt. nern2) then
        nt = nt + ndrsrd(2,nst) - ndrsrd(1,nern2 + 1) + 1
    end if

    if (nt .gt. ndrsmx) then
        write (ux8a,'(i5)') ndrsmx
        write (ux8b,'(i5)') nt
        j2 = ilnobl(ux8a)
        j3 = ilnobl(ux8b)
        write (noutpt,1000) ux8a(1:j2),ux8b(1:j3)
        write (nttyo,1000)  ux8a(1:j2),ux8b(1:j3)
1000 format(/' * Error - (EQLIB/mdrgex.f) Have insufficient array',' space to expand the',/7x,'reaction coefficient arrays to',' accommodate changes to the',/7x,'ion exchanger section'," of the 'd' set of reaction data. Increase",/7x,'the',' dimensioning parameter ndrspa from ',a,' to at least ',a,'.')

        stop
    end if

    if (nern1 .gt. 1) then
        nmax = ndrsrd(2,nern1 - 1)
        call copyaa(cdrsd,cdrsx,nmax)
        call copyia(ndrsd,ndrsx,nmax)
    end if

    if (nern1 .gt. 1) then
        nn = ndrsrd(2,nern1 - 1)
    else
        nn = 0
    end if

    do ns = nern1,nern2
        nr1 = ndrsr(1,ns)
        nr2 = ndrsr(2,ns)

        do n = nr1,nr2
            nn = nn + 1
            ndrsx(nn) = ndrs(n)
            cdrsx(nn) = cdrs(n)
        end do
    end do

    nn = ndrsrx(2,nern2)

    do ns = nern2 + 1,nst
        nr1 = ndrsrd(1,ns)
        nr2 = ndrsrd(2,ns)

        do n = nr1,nr2
            nn = nn + 1
            ndrsx(nn) = ndrsd(n)
            cdrsx(nn) = cdrsd(n)
        end do
    end do

    do nb = 1,nbtmax
        ns = nbasp(nb)

        if (ns.ge.nern1 .and. ns.le.nern2) then
            nbaspd(nb) = nbasp(nb)
        end if
    end do

    ! Copy the scratch arrays into the corresponding 'd' set arrays.
    ! Copy the new reactions in the 'x' set into the standard arrays
    ! in the 'd' set.
    ! Calling sequence substitutions:
    !   adhfsd for adhfs
    !   advfsd for advfs
    !   axhfsd for axhfs
    !   axlksd for axlks
    !   axvfsd for avhfs
    !   cdrsd for cdrs
    !   ndrsd for ndrs
    !   ndrsrd for ndrsr
    call cdrscx(adhfsd,adhfsx,advfsd,advfsx,axhfsd,axhfsx,axlksd,axlksx,axvfsd,axvfsx,cdrsd,cdrsx,ipch,ipchmx,ipcv,ipcvmx,narxmx,ndrsd,ndrsmx,ndrsx,ndrsrd,ndrsrx,nstmax,ntprmx)
end subroutine mdrgex