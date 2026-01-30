subroutine switch(adhfs,adhfsx,advfs,advfsx,axhfs,axhfsx,axlks,axlksx,axvfs,axvfsx,cdrs,cdrsx,eps100,ipch,ipchmx,ipcv,ipcvmx,jflag,jsflag,narn1,narxmx,nbasp,nbaspd,nbaspx,nb,nbt,nbtmax,nbw,ndrs,ndrsmx,ndrsx,ndrsr,ndrsrx,noutpt,ns2,nst,nstmax,ntprmx,nttyo,qbassw,qbswok,uspec)
    !! This subroutine executes an ordinary basis switch. The ns2-th
    !! species (not in the current active basis set) is exchanged into
    !! that set for the ns1-th species. All reactions are rewritten where
    !! necessary for consistency. Here nbaspx is the unmodified copy of
    !! the nbasp array. It contains a record of the basis set prior to
    !! basis switching.
    !! Basis switching is subject to the following rules:
    !!    1. Neither species may be suppressed.
    !!    2. The species switched in must have reaction that gives the
    !!       species to be switched out as a product.
    !!    3. The species to be switched into the active basis set must
    !!       not already occupy a position of its own in that set (e.g.,
    !!       be an active auxiliary basis species). To switch an active
    !!       auxiliary basis species into the strict basis set (e.g.,
    !!       make a special basis switch), use EQLIB/swtchb.f.
    !! The species indices are not interchanged by the switch. If both
    !! species are basis species, however, their basis indices are
    !! interchanged.
    !! This subroutine is called by:
    !!   EQLIB/autosw.f
    !!   EQ3NR/eq3nr.f
    !!   EQ6/eq6.f
    !!   EQ6/eqcalc.f
    !!   EQ6/eqphas.f
    !! Principal input:
    !! Principal output:
    !!   qbassw = logical flag:
    !!              = .true. if the basis set after the current switch
    !!                is not identical to the data file basis set
    !!   qbswok = logical flag:
    !!              = .true. if the current switch was completed okay
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
    integer :: narn1
    integer :: nb
    integer :: nbt
    integer :: nbw
    integer :: nst
    integer :: ns2

    logical :: qbassw

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
    integer :: jlen1
    integer :: jlen2
    integer :: jfl
    integer :: j2
    integer :: n
    integer :: nbb
    integer :: nb2
    integer :: nmax
    integer :: nrf1
    integer :: nrf2
    integer :: nrl1
    integer :: nrl2
    integer :: nr1
    integer :: nr2
    integer :: ns
    integer :: nse
    integer :: nsi
    integer :: ns1
    integer :: ntf
    integer :: nx

    integer :: ilnobl
    integer :: nbasis

    logical :: qbswok

    character(len=56) :: uspn56
    character(len=56) :: usp156
    character(len=56) :: usp256
    character(len=8) :: ux8

    real(kind=8) :: axx
    real(kind=8) :: cx
    real(kind=8) :: cxes
    real(kind=8) :: cxe1
    real(kind=8) :: cxe2
    real(kind=8) :: cx1s
    real(kind=8) :: cx11
    real(kind=8) :: cx12
    real(kind=8) :: stofac

    real(kind=8) :: coefdr

    ! Zero scratch arrays.
    nmax = narxmx*ntprmx*ipchmx*nstmax
    call initaz(adhfsx,nmax)

    nmax = narxmx*ntprmx*ipcvmx*nstmax
    call initaz(advfsx,nmax)

    nmax = narxmx*ntprmx*nstmax
    call initaz(axhfsx,nmax)
    call initaz(axlksx,nmax)
    call initaz(axvfsx,nmax)

    nmax = 2*nstmax
    call initiz(ndrsrx,nmax)

    nmax = ndrsmx
    call initaz(cdrsx,nmax)
    call initiz(ndrsx,nmax)

    qbswok = .false.

    ns1 = nbaspx(nb)

    ! Check to see if the switch is okay.
    call swtchk(cdrs,jflag,jsflag,nbaspx,nbt,nbtmax,ndrs,ndrsmx,ndrsr,noutpt,ns1,ns2,nstmax,nttyo,uspec)

    ! Calling sequence substitutions:
    !   jlen1 for jlen
    !   uspec(ns1) for unam48
    !   usp156 for uspn56
    call fmspnx(jlen1,uspec(ns1),usp156)

    ! Calling sequence substitutions:
    !   jlen2 for jlen
    !   uspec(ns2) for unam48
    !   usp256 for uspn56
    call fmspnx(jlen2,uspec(ns2),usp256)

    ! Check whether or not the species being switched in is not already
    ! in the active basis set (associated with some other mass balance).
    ! This is necessary to ensure that the switch to be made is an
    ! ordinary basis switch (the only type allowed in the present
    ! subroutine).
    ! Calling sequence substitutions:
    !   nbaspx for nbasp
    !   ns2 for ns
    nb2 = nbasis(nbaspx,nbt,nbtmax,ns2)

    if (nb2.gt.0 .and. jflag(ns2).ne.30) then
        write (noutpt,1000) usp156(1:jlen1),usp256(1:jlen2)
        write (nttyo,1000) usp156(1:jlen1),usp256(1:jlen2)
1000 format(/' * Error - (EQLIB/switch) Programming error trap:'," Can't replace",/7x,'the species ', a,' in the basis set with',' ',a,',',/7x,'doing an ordinary basis switch, because the',' former is already',/7x,'in the active basis set.')

        stop
    end if

    write (noutpt,1010) usp156(1:jlen1),usp256(1:jlen2)
    write (nttyo,1010) usp156(1:jlen1),usp256(1:jlen2)
1010 format(/' Making an ordinary basis switch: replacing ',a,/1x,'in the active basis set with ',a,'.',/)

    if (nb .eq. nbw) then
        ! Switching solvent water out of the basis set. Write a warning.
        write (noutpt,1020) usp156(1:jlen1)
        write (nttyo,1020) usp156(1:jlen1)
1020 format(/' * Warning - (EQLIB/switch) Making an ordinary',' basis switch',/7x,'which moves ',a,' out of the active',' basis set.')
    end if

    write (noutpt,1040)
1040 format(/'  Before the switch:',/)

    ! Print the existing reaction for the ns1-th species, if any.
    ! Calling sequence substitutions:
    !   noutpt for nf
    !   ns1 for ns
    call prreac(cdrs,ndrs,ndrsmx,ndrsr,noutpt,ns1,nstmax,uspec)

    ! Print the existing linking reaction.
    ! Calling sequence substitutions:
    !   noutpt for nf
    !   ns2 for ns
    call prreac(cdrs,ndrs,ndrsmx,ndrsr,noutpt,ns2,nstmax,uspec)

    write (noutpt,1050)
1050 format(1x)

    ! Calling sequence substitutions:
    !   ns1 for nse
    !   ns2 for ns
    cx12 = coefdr(cdrs,ndrs,ndrsmx,ndrsr,ns1,ns2,nstmax)
    nrl1 = ndrsr(1,ns2)
    nrl2 = ndrsr(2,ns2)

    nrf1 = ndrsr(1,ns1)
    nrf2 = ndrsr(2,ns1)
    cx11 = cdrs(nrf1)

    write (ux8,'(i5)') ndrsmx
    call lejust(ux8)
    j2 = ilnobl(ux8)

    nx = 0

    do ns = 1,nst
        nr1 = ndrsr(1,ns)
        nr2 = ndrsr(2,ns)
        ndrsrx(1,ns) = nx + 1

        if (ns .eq. ns1) then
            ! Invert the linking reaction. Put the coefficient for the
            ! ns1-th species first in the range.
            nx = nx + 1

            if (nx .gt. ndrsmx) then
                write (noutpt,1070) ux8(1:j2),usp256(1:jlen2),usp156(1:jlen1),usp156(1:jlen1)
                write (nttyo,1070) ux8(1:j2),usp256(1:jlen2),usp156(1:jlen1),usp156(1:jlen1)
1070 format(/' * Error - (EQLIB/switch) The maximum ',a,' entries in the cdrs and ndrs',/7x,'arrays has been',' exceeded in trying to switch ',a,/7x,'into the active',' basis set for ',a,' while writing',/7x,'the reaction',' for ',a,'. Increase the value of the',/7x,' dimensioning parameter ndrspa.')

                stop
            end if

            cdrsx(nx) = -cx12
            ndrsx(nx) = ns

            do n = nrl1,nrl2
                if (ndrs(n) .ne. ns) then
                    nx = nx + 1

                    if (nx .gt. ndrsmx) then
                        write (noutpt,1070) ux8(1:j2),usp256(1:jlen2),usp156(1:jlen1),usp156(1:jlen1)
                        write (nttyo,1070) ux8(1:j2),usp256(1:jlen2),usp156(1:jlen1),usp156(1:jlen1)
                        stop
                    end if

                    cdrsx(nx) = -cdrs(n)
                    ndrsx(nx) = ndrs(n)
                end if
            end do

            ! Log K coefficients.
            do j = 1,ntprmx
                do i = 1,narxmx
                    axlksx(i,j,ns) = -axlks(i,j,ns2)
                end do
            end do

            if (ipch .ge. 0) then
                ! Enthalpy function coefficients.
                do j = 1,ntprmx
                    do i = 1,narxmx
                        axhfsx(i,j,ns) = -axhfs(i,j,ns2)
                    end do
                end do

                do ipc = 1,ipch
                    do j = 1,ntprmx
                        do i = 1,narxmx
                            adhfsx(i,j,ipc,ns) = -adhfs(i,j,ipc,ns2)
                        end do
                    end do
                end do
            end if

            if (ipcv .ge. 0) then
                ! Volume function coefficients.
                do j = 1,ntprmx
                    do i = 1,narxmx
                        axvfsx(i,j,ns) = -axvfs(i,j,ns2)
                    end do
                end do

                do ipc = 1,ipcv
                    do j = 1,ntprmx
                        do i = 1,narxmx
                            advfsx(i,j,ipc,ns) = -advfs(i,j,ipc,ns2)
                        end do
                    end do
                end do
            end if
        else if (ns .eq. ns2) then
            ntf = nrf2 - nrf1 + 1

            if (ntf .lt. 2) then
                ! If the ns1-th species was a strict basis species, make
                ! the ns2-th species into one. Write a null reaction.
                nx = nx + 1

                if (nx .gt. ndrsmx) then
                    write (noutpt,1070) ux8(1:j2),usp256(1:jlen2),usp156(1:jlen1),usp256(1:jlen2)
                    write (nttyo,1070) ux8(1:j2),usp256(1:jlen2),usp156(1:jlen1),usp256(1:jlen2)
                    stop
                end if

                cdrsx(nx) = 0.
                ndrsx(nx) = 0

                ! Log K coefficients. Log K = 0 for a strict basis species.
                do j = 1,ntprmx
                    do i = 1,narxmx
                        axlksx(i,j,ns) = 0.
                    end do
                end do

                if (ipch .ge. 0) then
                    ! Enthalpy function coefficients. This function is the
                    ! partial molar enthalpy of formation for a strict
                    ! basis species.
                    do j = 1,ntprmx
                        do i = 1,narxmx
                            axx = axhfs(i,j,ns)

                            do n = nr1 + 1,nr2
                                nsi = ndrs(n)
                                cx = cdrs(n)
                                axx = axx - cx*axhfs(i,j,nsi)
                            end do

                            cx = cdrs(nr1)
                            axhfsx(i,j,ns) = axx/cx
                        end do
                    end do

                    do ipc = 1,ipch
                        do j = 1,ntprmx
                            do i = 1,narxmx
                                axx = adhfs(i,j,ipc,ns)

                                do n = nr1 + 1,nr2
                                    nsi = ndrs(n)
                                    cx = cdrs(n)
                                    axx = axx - cx*adhfs(i,j,ipc,nsi)
                                end do

                                cx = cdrs(nr1)
                                adhfsx(i,j,ipc,ns) = axx/cx
                            end do
                        end do
                    end do
                end if

                if (ipcv .ge. 0) then
                    ! Volume function coefficients. This function is the
                    ! partial molar volume for a strict basis species.
                    do j = 1,ntprmx
                        do i = 1,narxmx
                            axx = axvfs(i,j,ns)

                            do n = nr1 + 1,nr2
                                nsi = ndrs(n)
                                cx = cdrs(n)
                                axx = axx - cx*axvfs(i,j,nsi)
                            end do

                            cx = cdrs(nr1)
                            axvfsx(i,j,ns) = axx/cx
                        end do
                    end do

                    do ipc = 1,ipcv
                        do j = 1,ntprmx
                            do i = 1,narxmx
                                axx = advfs(i,j,ipc,ns)

                                do n = nr1 + 1,nr2
                                    nsi = ndrs(n)
                                    cx = cdrs(n)
                                    axx = axx - cx*advfs(i,j,ipc,nsi)
                                end do

                                cx = cdrs(nr1)
                                advfsx(i,j,ipc,ns) = axx/cx
                            end do
                        end do
                    end do
                end if
            else
                ! If the ns1-th species was an auxiliary basis species, make
                ! the ns2-th species into one. Write a reaction which is the
                ! original reaction for the breakdown of the the ns1-th
                ! species, but using the original reaction for the ns1-th
                ! species instead of the linking reaction to eliminate the
                ! ns1-th species.
                stofac = -cx12/cx11

                ! Go through the species of the existing reaction for the
                ! ns2-th species (this is the linking reaction).
                do n = nr1,nr2
                    nse = ndrs(n)

                    ! Calling sequence substitutions:
                    !   ns1 for ns
                    cxe1 = coefdr(cdrs,ndrs,ndrsmx,ndrsr,nse,ns1,nstmax)
                    cx = cdrs(n) + stofac*cxe1

                    if (abs(cx) .gt. eps100) then
                        nx = nx + 1

                        if (nx .gt. ndrsmx) then
                            write (noutpt,1070) ux8(1:j2),usp256(1:jlen2),usp156(1:jlen1),usp256(1:jlen2)
                            write (nttyo,1070) ux8(1:j2),usp256(1:jlen2),usp156(1:jlen1),usp256(1:jlen2)
                            stop
                        end if

                        cdrsx(nx) = cx
                        ndrsx(nx) = nse
                    end if
                end do

                ! Add to the new reaction any species in the original
                ! reaction for the breakdown of the ns1-th species.
                do n = nrf1,nrf2
                    nse = ndrs(n)
                    cxes = coefdr(cdrs,ndrs,ndrsmx,ndrsr,nse,ns,nstmax)

                    if (cxes .eq. 0.) then
                        cx = stofac*cdrs(n)
                        nx = nx + 1

                        if (nx .gt. ndrsmx) then
                            write (noutpt,1070) ux8(1:j2),usp256(1:jlen2),usp156(1:jlen1),usp256(1:jlen2)
                            write (nttyo,1070) ux8(1:j2),usp256(1:jlen2),usp156(1:jlen1),usp256(1:jlen2)
                            stop
                        end if

                        cdrsx(nx) = cx
                        ndrsx(nx) = nse
                    end if
                end do

                ! Log K coefficients.
                do j = 1,ntprmx
                    do i = 1,narxmx
                        axlksx(i,j,ns) = axlks(i,j,ns) + stofac*axlks(i,j,ns1)
                    end do
                end do

                if (ipch .ge. 0) then
                    ! Enthalpy function coefficients.
                    do j = 1,ntprmx
                        do i = 1,narxmx
                            axhfsx(i,j,ns) = axhfs(i,j,ns) + stofac*axhfs(i,j,ns1)
                        end do
                    end do

                    do ipc = 1,ipch
                        do j = 1,ntprmx
                            do i = 1,narxmx
                                adhfsx(i,j,ipc,ns) = adhfs(i,j,ipc,ns)              + stofac*adhfs(i,j,ipc,ns1)
                            end do
                        end do
                    end do
                end if

                if (ipcv .ge. 0) then
                    ! Volume function coefficients.
                    do j = 1,ntprmx
                        do i = 1,narxmx
                            axvfsx(i,j,ns) = axvfs(i,j,ns)            + stofac*axvfs(i,j,ns1)
                        end do
                    end do

                    do ipc = 1,ipcv
                        do j = 1,ntprmx
                            do i = 1,narxmx
                                advfsx(i,j,ipc,ns) = advfs(i,j,ipc,ns)              + stofac*advfs(i,j,ipc,ns1)
                            end do
                        end do
                    end do
                end if
            end if
        else
            ! Do reactions for all other species.
            ! Calling sequence substitutions:
            !   uspec(ns) for unam48
            call fmspnx(jlen,uspec(ns),uspn56)

            ! Calling sequence substitutions:
            !   ns1 for nse
            cx1s = coefdr(cdrs,ndrs,ndrsmx,ndrsr,ns1,ns,nstmax)

            if (cx1s .eq. 0.) then
                ! Have a reaction which is not connected. Copy as is.
                do n = nr1,nr2
                    nx = nx + 1

                    if (nx .gt. ndrsmx) then
                        write (noutpt,1070) ux8(1:j2),usp256(1:jlen2),usp156(1:jlen1),uspn56(1:jlen)
                        write (nttyo,1070) ux8(1:j2),usp256(1:jlen2),usp156(1:jlen1),uspn56(1:jlen)
                        stop
                    end if

                    cdrsx(nx) = cdrs(n)
                    ndrsx(nx) = ndrs(n)
                end do

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
                ! Have a reaction which is connected.
                stofac = -cx1s/cx12

                ! Go through the species of the existing reaction for the
                ! ns-th species.
                do n = nr1,nr2
                    nse = ndrs(n)

                    ! Calling sequence substitutions:
                    !   ns2 for ns
                    cxe2 = coefdr(cdrs,ndrs,ndrsmx,ndrsr,nse,ns2,nstmax)
                    cx = cdrs(n) + stofac*cxe2

                    if (abs(cx) .gt. eps100) then
                        nx = nx + 1

                        if (nx .gt. ndrsmx) then
                            write (noutpt,1070) ux8(1:j2),usp256(1:jlen2),usp156(1:jlen1),uspn56(1:jlen)
                            write (nttyo,1070) ux8(1:j2),usp256(1:jlen2),usp156(1:jlen1),uspn56(1:jlen)
                            stop
                        end if

                        cdrsx(nx) = cx
                        ndrsx(nx) = nse
                    end if
                end do

                ! Add to the new reaction any species in the linking reaction
                ! that were not in the original reaction for the breakdown of
                ! the ns-th species.
                do n = nrl1,nrl2
                    nse = ndrs(n)
                    cxes = coefdr(cdrs,ndrs,ndrsmx,ndrsr,nse,ns,nstmax)

                    if (cxes .eq. 0.) then
                        cx = stofac*cdrs(n)
                        nx = nx + 1

                        if (nx .gt. ndrsmx) then
                            write (noutpt,1070) ux8(1:j2),usp256(1:jlen2),usp156(1:jlen1),uspn56(1:jlen)
                            write (nttyo,1070) ux8(1:j2),usp256(1:jlen2),usp156(1:jlen1),uspn56(1:jlen)
                            stop
                        end if

                        cdrsx(nx) = cx
                        ndrsx(nx) = nse
                    end if
                end do

                ! Log K coefficients.
                do j = 1,ntprmx
                    if (axlks(1,j,ns) .lt. 9999999.) then
                        do i = 1,narxmx
                            axlksx(i,j,ns) = axlks(i,j,ns) + stofac*axlks(i,j,ns2)
                        end do
                    end if
                end do

                if (ipch .ge. 0) then
                    ! Enthalpy function coefficients.
                    do j = 1,ntprmx
                        if (axhfs(1,j,ns) .lt. 9999999.) then
                            do i = 1,narxmx
                                axhfsx(i,j,ns) = axhfs(i,j,ns)              + stofac*axhfs(i,j,ns2)
                            end do
                        end if
                    end do

                    do ipc = 1,ipch
                        do j = 1,ntprmx
                            do i = 1,narxmx
                                adhfsx(i,j,ipc,ns) = adhfs(i,j,ipc,ns)              + stofac*adhfs(i,j,ipc,ns2)
                            end do
                        end do
                    end do
                end if

                if (ipcv .ge. 0) then
                    ! Enthalpy function coefficients.
                    do j = 1,ntprmx
                        do i = 1,narxmx
                            if (axvfs(1,j,ns) .lt. 9999999.) then
                                axvfsx(i,j,ns) = axvfs(i,j,ns)              + stofac*axvfs(i,j,ns2)
                            end if
                        end do
                    end do

                    do ipc = 1,ipcv
                        do j = 1,ntprmx
                            do i = 1,narxmx
                                advfsx(i,j,ipc,ns) = advfs(i,j,ipc,ns)              + stofac*advfs(i,j,ipc,ns2)
                            end do
                        end do
                    end do
                end if
            end if
        end if

        ndrsrx(2,ns) = nx
    end do

    ! Copy the new reactions from the 'x' arrays into the standard
    ! arrays.
    call cdrscx(adhfs,adhfsx,advfs,advfsx,axhfs,axhfsx,axlks,axlksx,axvfs,axvfsx,cdrs,cdrsx,ipch,ipchmx,ipcv,ipcvmx,narxmx,ndrs,ndrsmx,ndrsx,ndrsr,ndrsrx,nstmax,ntprmx)

    write (noutpt,1090)
1090 format(/'  After the switch:',/)

    ! Print the new linking reaction.
    ! Calling sequence substitutions:
    !   noutpt for nf
    !   ns1 for ns
    call prreac(cdrs,ndrs,ndrsmx,ndrsr,noutpt,ns1,nstmax,uspec)

    ! Print the new reaction for the ns2-th species.
    ! Calling sequence substitutions:
    !   noutpt for nf
    !   ns2 for ns
    call prreac(cdrs,ndrs,ndrsmx,ndrsr,noutpt,ns2,nstmax,uspec)

    write (noutpt,1050)

    ! Exchange jflag values.
    jfl = jflag(ns1)
    jflag(ns1) = jflag(ns2)
    jflag(ns2) = jfl

    ! Clear the content in the nbasp and nbaspx arrays which constitutes
    ! the instruction to make the current switch. Change the nbaspx
    ! entry to clear the instruction for the current switch.
    nbaspx(nb) = ns2

    ! If necessary, reset the basis index of water (nbw).
    if (ns1 .eq. narn1) then
        nbw = 0
    else if (ns2 .eq. narn1) then
        nbw = nb
    end if

    ! Set the flag indicating that the switch completed successfully.
    qbswok = .true.

    ! Check the value of the qbassw flag. Set it to .true. if
    ! the basis set is not identical to the corresponding data file
    ! basis set.
    qbassw = .false.

    do nbb = 1,nbt
        ns1 = nbaspd(nbb)
        ns2 = nbasp(nbb)

        if (ns2 .ne. ns1) then
            go to 200
        end if
    end do

    qbassw = .true.
200 continue
end subroutine switch
