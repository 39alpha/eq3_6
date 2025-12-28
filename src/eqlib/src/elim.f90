subroutine elim(adhfs,adhfsx,advfs,advfsx,axhfs,axhfsx,axlks,axlksx,axvfs,axvfsx,cdrs,cdrsx,eps100,ipch,ipchmx,ipcv,ipcvmx,jsflag,narxmx,ndrs,ndrsmx,ndrsx,ndrsr,ndrsrx,nse,nst,nstmax,ntprmx,noutpt,nttyo,uspec)
    !! This subroutine rewrites reaction equations so that the auxiliary
    !! basis species with index nse and jflag = 30 is eliminated from
    !! all reactions except the one linking it with its corresponding
    !! strict basis variable. The polynomial coefficients for computing
    !! the equilibrium coefficients are recomputed accordingly.
    !! This subroutine is called by:
    !!   EQ3NR/eq3nr.f
    !!   EQ6/eq6.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: ipchmx
    integer :: ipcvmx
    integer :: narxmx
    integer :: ndrsmx
    integer :: nstmax
    integer :: ntprmx

    integer :: jsflag(nstmax)
    integer :: ndrs(ndrsmx)
    integer :: ndrsx(ndrsmx)
    integer :: ndrsr(2,nstmax)
    integer :: ndrsrx(2,nstmax)
    integer :: ipch
    integer :: ipcv
    integer :: nse
    integer :: nst
    integer :: noutpt
    integer :: nttyo

    character(len=48) :: uspec(nstmax)

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
    integer :: i
    integer :: ipc
    integer :: j
    integer :: jlen
    integer :: n
    integer :: nnx
    integer :: nre1
    integer :: nre2
    integer :: nr1
    integer :: nr2
    integer :: ns
    integer :: nss
    integer :: nt

    character(len=56) :: uspn56

    real(kind=8) :: cx
    real(kind=8) :: cxe
    real(kind=8) :: cxee
    real(kind=8) :: cxse
    real(kind=8) :: cxss
    real(kind=8) :: stofac

    real(kind=8) :: coefdr

    nt = ndrsr(2,nse) - ndrsr(1,nse) + 1

    if (nt .lt. 2) then
        ! Calling sequence substitutions:
        !   uspec(nse) for unam48
        call fmspnx(jlen,uspec(nse),uspn56)
        write (noutpt,1000) uspn56(1:jlen)
        write (nttyo,1000) uspn56(1:jlen)
1000 format(/' * Error - (EQLIB/elim) The species ',a,/7x,'is in the strict basis and therefore can not be',' eliminated',/7x,'from the working basis set.')

        stop
    end if

    ! Calling sequence substitutions:
    !   nse for ns
    cxee = coefdr(cdrs,ndrs,ndrsmx,ndrsr,nse,nse,nstmax)
    nre1 = ndrsr(1,nse)
    nre2 = ndrsr(2,nse)

    nnx = 0

    do ns = 1,nst
        nr1 = ndrsr(1,ns)
        nr2 = ndrsr(2,ns)
        ndrsrx(1,ns) = nnx + 1
        cxe = coefdr(cdrs,ndrs,ndrsmx,ndrsr,nse,ns,nstmax)

        if (ns.eq.nse .or. cxe.eq.0. .or. jsflag(ns).ge.2) then
            ! Have found a reaction that is not to be changed.
            do n = nr1,nr2
                nnx = nnx + 1

                if (nnx .gt. ndrsmx) then
                    ! Calling sequence substitutions:
                    !   uspec(nse) for unam48
                    call fmspnx(jlen,uspec(nse),uspn56)
                    write (noutpt,1005) ndrsmx,uspn56(1:jlen),ndrsmx
                    write (nttyo,1005) ndrsmx,uspn56(1:jlen),ndrsmx
1005 format(/' * Error - (EQLIB/elim) The maximum ',i7,' entries in the',/7x,'cdrs and ndrs arrays has been',' exceeded in trying to eliminate',/7x,'the species ',a,' from the working basis set.',/7x,'Increase the',' dimensioning  parameter ndrspa from its current',/7x,'value of ',i6,'.')

                    stop
                end if

                cdrsx(nnx) = cdrs(n)
                ndrsx(nnx) = ndrs(n)
            end do

            ndrsrx(2,ns) = nnx

            ! Log K coefficients.
            do j = 1,ntprmx
                do i = 1,narxmx
                    axlksx(i,j,ns) = axlks(i,j,ns)
                end do
            end do

            if (ipch .ge. 0) then
                ! Enthalpy function coefficients.
                do j = 1,ntprmx
                    do i = 1,narxmx
                        axhfsx(i,j,ns) = axhfs(i,j,ns)
                    end do
                end do

                do ipc = 1,ipch
                    do j = 1,ntprmx
                        do i = 1,narxmx
                            adhfsx(i,j,ipc,ns) = adhfs(i,j,ipc,ns)
                        end do
                    end do
                end do
            end if

            if (ipcv .ge. 0) then
                ! Volume function coefficients.
                do j = 1,ntprmx
                    do i = 1,narxmx
                        axvfsx(i,j,ns) = axvfs(i,j,ns)
                    end do
                end do

                do ipc = 1,ipcv
                    do j = 1,ntprmx
                        do i = 1,narxmx
                            advfsx(i,j,ipc,ns) = advfs(i,j,ipc,ns)
                        end do
                    end do
                end do
            end if
        else
            stofac = cxe/cxee

            ! Do species appearing in the existing reaction for the ns-th
            ! species.
            do n = nr1,nr2
                nss = ndrs(n)

                ! Calling sequence substitutions:
                !   nss for nse
                !   nse for ns
                cxse = coefdr(cdrs,ndrs,ndrsmx,ndrsr,nss,nse,nstmax)
                cx = cdrs(n) - stofac*cxse

                if (abs(cx) .le. eps100) then
                    cx = 0.
                end if

                if (cx .ne. 0.) then
                    nnx = nnx + 1

                    if (nnx .gt. ndrsmx) then
                        ! Calling sequence substitutions:
                        !   uspec(nse) for unam48
                        call fmspnx(jlen,uspec(nse),uspn56)
                        write (noutpt,1005) ndrsmx,uspn56(1:jlen),ndrsmx
                        write (nttyo,1005) ndrsmx,uspn56(1:jlen),ndrsmx
                        stop
                    end if

                    cdrsx(nnx) = cx
                    ndrsx(nnx) = nss
                end if
            end do

            ! Do species appearing in the existing reaction for the nse-th
            ! species but not in that for the ns-th species.
            do n = nre1,nre2
                nss = ndrs(n)

                ! Calling sequence substitutions:
                !   nss for nse
                !   nse for ns
                cxse = coefdr(cdrs,ndrs,ndrsmx,ndrsr,nss,nse,nstmax)

                ! Calling sequence substitutions:
                !   nss for nse
                cxss = coefdr(cdrs,ndrs,ndrsmx,ndrsr,nss,ns,nstmax)

                if (cxss .eq. 0.) then
                    cx = -stofac*cxse

                    if (abs(cx) .le. eps100) then
                        cx = 0.
                    end if

                    if (cx .ne. 0.) then
                        nnx = nnx + 1

                        if (nnx .gt. ndrsmx) then
                            ! Calling sequence substitutions:
                            !   uspec(nse) for unam48
                            call fmspnx(jlen,uspec(nse),uspn56)
                            write (noutpt,1005) ndrsmx,uspn56(1:jlen),ndrsmx
                            write (nttyo,1005) ndrsmx,uspn56(1:jlen),ndrsmx
                            stop
                        end if

                        cdrsx(nnx) = cx
                        ndrsx(nnx) = nss
                    end if
                end if
            end do

            ! If the new reaction has no entries, put in a null entry.
            if (nnx .lt. ndrsrx(1,ns)) then
                nnx = nnx + 1

                if (nnx .gt. ndrsmx) then
                    ! Calling sequence substitutions:
                    !   uspec(nse) for unam48
                    call fmspnx(jlen,uspec(nse),uspn56)
                    write (noutpt,1005) ndrsmx,uspn56(1:jlen),ndrsmx
                    write (nttyo,1005) ndrsmx,uspn56(1:jlen),ndrsmx
                    stop
                end if

                cdrsx(nnx) = 0.
                ndrsx(nnx) = 0
            end if

            ndrsrx(2,ns) = nnx

            ! Log K coefficients.
            do j = 1,ntprmx
                if (axlks(1,j,ns) .lt. 9999999.) then
                    do i = 1,narxmx
                        axlksx(i,j,ns) = axlks(i,j,ns) - stofac*axlks(i,j,nse)
                    end do
                end if
            end do

            if (ipch .ge. 0) then
                ! Enthalpy function coefficients.
                do j = 1,ntprmx
                    if (axhfs(1,j,ns) .lt. 9999999.) then
                        do i = 1,narxmx
                            axhfsx(i,j,ns) = axhfs(i,j,ns) - stofac*axhfs(i,j,nse)
                        end do
                    end if
                end do

                do ipc = 1,ipch
                    do j = 1,ntprmx
                        do i = 1,narxmx
                            adhfsx(i,j,ipc,ns) = adhfs(i,j,ipc,ns)            - stofac*adhfs(i,j,ipc,nse)
                        end do
                    end do
                end do
            end if

            if (ipcv .ge. 0) then
                ! Volume function coefficients.
                do j = 1,ntprmx
                    if (axvfs(1,j,ns) .lt. 9999999.) then
                        do i = 1,narxmx
                            axvfsx(i,j,ns) = axvfs(i,j,ns) - stofac*axvfs(i,j,nse)
                        end do
                    end if
                end do

                do ipc = 1,ipcv
                    do j = 1,ntprmx
                        do i = 1,narxmx
                            advfsx(i,j,ipc,ns) = advfs(i,j,ipc,ns)            - stofac*advfs(i,j,ipc,nse)
                        end do
                    end do
                end do
            end if
        end if
    end do

    ! Copy the new reactions into the standard arrays.
    call cdrscx(adhfs,adhfsx,advfs,advfsx,axhfs,axhfsx,axlks,axlksx,axvfs,axvfsx,cdrs,cdrsx,ipch,ipchmx,ipcv,ipcvmx,narxmx,ndrs,ndrsmx,ndrsx,ndrsr,ndrsrx,nstmax,ntprmx)
end subroutine elim