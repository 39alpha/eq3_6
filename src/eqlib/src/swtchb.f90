subroutine swtchb(adhfsd,adhfsx,advfsd,advfsx,axhfsd,axhfsx,axlksd,axlksx,axvfsd,axvfsx,cdrsd,cdrsx,ipch,ipchmx,ipcv,ipcvmx,narxmx,nbaspd,nbtmax,nbw,nb1,nb2,ndrsd,ndrsmx,ndrsx,ndrsrd,ndrsrx,noutpt,ns1,ns2,nsta,nstmax,ntprmx,nttyo,uspeca)
    !! This subroutine performs a special kind of basis switch. It
    !! interchanges the status of a strict basis species with an
    !! auxiliary basis species. The effect of this action is equivalent
    !! to modifying the data file so that a different basis species is
    !! the strict basis species corresponding to a chemical element.
    !! Here nb1 is the basis index of the species originally in the
    !! strict set, and nb2 is the basis index of the species originally
    !! in the auxiliary basis set. Here also ns1 and ns2 are the
    !! corresponding species indices. The species indices are not
    !! interchanged by the switch. The basis indices, as represented by
    !! the contents of the nbaspd array, are interchanged.
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
    integer :: nbtmax
    integer :: ndrsmx
    integer :: nstmax
    integer :: ntprmx

    integer :: noutpt
    integer :: nttyo

    integer :: nbaspd(nbtmax)
    integer :: ndrsd(ndrsmx)
    integer :: ndrsx(ndrsmx)
    integer :: ndrsrd(2,nstmax)
    integer :: ndrsrx(2,nstmax)
    integer :: ipch
    integer :: ipcv
    integer :: nbw
    integer :: nb1
    integer :: nb2
    integer :: ns1
    integer :: ns2
    integer :: nsta

    character(len=48) :: uspeca(nstmax)

    real(kind=8) :: adhfsd(narxmx,ntprmx,ipchmx,nstmax)
    real(kind=8) :: adhfsx(narxmx,ntprmx,ipchmx,nstmax)
    real(kind=8) :: advfsd(narxmx,ntprmx,ipcvmx,nstmax)
    real(kind=8) :: advfsx(narxmx,ntprmx,ipcvmx,nstmax)
    real(kind=8) :: axhfsd(narxmx,ntprmx,nstmax)
    real(kind=8) :: axhfsx(narxmx,ntprmx,nstmax)
    real(kind=8) :: axlksd(narxmx,ntprmx,nstmax)
    real(kind=8) :: axlksx(narxmx,ntprmx,nstmax)
    real(kind=8) :: axvfsd(narxmx,ntprmx,nstmax)
    real(kind=8) :: axvfsx(narxmx,ntprmx,nstmax)
    real(kind=8) :: cdrsd(ndrsmx)
    real(kind=8) :: cdrsx(ndrsmx)

    ! Local variable declarations.
    integer :: i
    integer :: ipc
    integer :: j
    integer :: jlen
    integer :: jlen1
    integer :: jlen2
    integer :: n
    integer :: nerr
    integer :: nrl1
    integer :: nrl2
    integer :: nr1
    integer :: nr2
    integer :: ns
    integer :: nsi
    integer :: nt1
    integer :: nt2
    integer :: nx

    character(len=56) :: uspn56
    character(len=56) :: usp156
    character(len=56) :: usp256

    real(kind=8) :: axx
    real(kind=8) :: cx
    real(kind=8) :: cx12

    real(kind=8) :: coefdr

    nerr = 0

    ! Calling sequence substitutions:
    !   jlen1 for jlen
    !   uspeca(ns1) for unam48
    !   usp156 for uspn56
    call fmspnx(jlen1,uspeca(ns1),usp156)

    ! Calling sequence substitutions:
    !   jlen2 for jlen
    !   uspeca(ns2) for unam48
    !   usp256 for uspn56
    call fmspnx(jlen2,uspeca(ns2),usp256)

    if (nb1 .eq. nb2) then
        write (noutpt,1004) usp156(1:jlen1)
        write (nttyo,1004) usp156(1:jlen1)
1004 format(/" * Error - (EQLIB/swtchb) Can't make a special basis",' switch',/7x,'replacing ',a,' with itself.')

        nerr = nerr + 1
    end if

    if (nb2 .lt. nb1) then
        write (noutpt,1006) usp156(1:jlen1),usp256(1:jlen2)
        write (nttyo,1006) usp156(1:jlen1),usp256(1:jlen2)
1006 format(/" * Error - (EQLIB/swtchb) Can't make a special",' switch',/7x,'replacing ',a,' with ',a,' because these',' species are',/7x,'not in the proper hierarchical order.',' The former must appear before',/7x,'the latter in the list',' of basis species at the time the special switch',/7x,'is executed. Check the ordering on the data file and any',/7x,'changes made by previously executed special basis',' switches.')

        nerr = nerr + 1
    end if

    ! Make sure that the first species is in the strict basis set.
    nt1 = ndrsrd(2,ns1) - ndrsrd(1,ns1) + 1

    if (nt1 .ge. 2) then
        write (noutpt,1010) usp156(1:jlen1),usp256(1:jlen2)
        write (nttyo,1010) usp156(1:jlen1),usp256(1:jlen2)
1010 format(/" * Error - (EQLIB/swtchb) Can't make a special",' switch',/7x,'replacing ',a,' with ',a,' because the',' former species',/7x,'is not currently in the strict',' basis set.',/7x,'Check the status of this species on',' the data file and any',/7x,'changes made by previously',' executed special basis switches.')

        nerr = nerr + 1
    end if

    ! Check the linking reaction.
    nt2 = ndrsrd(2,ns2) - ndrsrd(1,ns2) + 1

    if (nt2 .lt. 2) then
        ! The second species must not be in the strict basis, because
        ! there is then no possibility of a linking reaction.
        write (noutpt,1030) usp156(1:jlen1),usp256(1:jlen2)
        write (nttyo,1030) usp156(1:jlen1),usp256(1:jlen2)
1030 format(/" * Error - (EQLIB/swtchb) Can't make a special",' switch',/7x,'replacing ',a,' with ',a,' because the',' latter species',/7x,'is currently in the strict basis',' set, precluding.',/7x,'the possibility of it having a',' reaction linking it to the',/7x,'former. Check the status',' of the latter species on the',/7x,'data file and any',' changes made by previously executed',/7x,'special',/7x,'basis switches.')

        nerr = nerr + 1
    else
        ! Make sure that the first species appears as a product in
        ! the reaction belonging to the second species.
        ! Calling sequence substitutions:
        !   cdrsd for cdrs
        !   ndrsd for ndrs
        !   ndrsrd for ndrsr
        !   ns1 for nse
        !   ns2 for ns
        cx12 = coefdr(cdrsd,ndrsd,ndrsmx,ndrsrd,ns1,ns2,nstmax)

        if (cx12 .eq. 0.) then
            write (noutpt,1040) usp156(1:jlen1),usp256(1:jlen2)
            write (nttyo,1040) usp156(1:jlen1),usp256(1:jlen2)
1040 format(/" * Error - (EQLIB/swtchb) Can't make a special",' switch',/7x,'replacing ',a,' with ',a,' because the',' former species',/7x,'does not appear in the reaction',' for the former species. Check',/7x,'the reaction',' on the data file and any changes made by',/7x,'previously executed special basis switches.')

            ns = ns2

            ! Calling sequence substitutions:
            !   cdrsd for cdrs
            !   ndrsd for ndrs
            !   ndrsrd for ndrsr
            !   noutpt for nf
            !   uspeca for uspec
            call prreac(cdrsd,ndrsd,ndrsmx,ndrsrd,noutpt,ns,nstmax,uspeca)
            nerr = nerr + 1
        end if
    end if

    if (nerr .gt. 0) then
        stop
    end if

    write (noutpt,1050) usp156(1:jlen1),usp256(1:jlen2)
    write (nttyo,1050) usp156(1:jlen1),usp256(1:jlen2)
1050 format(/' Making a special basis switch: replacing ',a,' with ',a,'.'/)

    if (nb1 .eq. nbw) then
        write (noutpt,1052) usp156(1:jlen1)
        write (nttyo,1052) usp156(1:jlen1)
1052 format(/' * Warning - (EQLIB/swtchb) Making a special basis',' switch',/7x,'replacing ',a,' with another species.')
    end if

    if (nb2 .eq. nbw) then
        write (noutpt,1054) usp256(1:jlen2)
        write (nttyo,1054) usp256(1:jlen2)
1054 format(/' * Warning - (EQLIB/swtchb) Making a special basis',' switch',/7x,'replacing a species with ',a,'.')
    end if

    ! Print the existing linking reaction.
    ! Calling sequence substitutions:
    !   cdrsd for cdrs
    !   ndrsd for ndrs
    !   ndrsrd for ndrsr
    !   noutpt for nf
    !   ns2 for ns
    !   uspeca for uspec
    call prreac(cdrsd,ndrsd,ndrsmx,ndrsrd,noutpt,ns2,nstmax,uspeca)

    nrl1 = ndrsrd(1,ns2)
    nrl2 = ndrsrd(2,ns2)

    nx = 0

    do ns = 1,nsta
        nr1 = ndrsrd(1,ns)
        nr2 = ndrsrd(2,ns)
        ndrsrx(1,ns) = nx + 1

        ! Calling sequence substitutions:
        !   uspeca(ns) for unam48
        call fmspnx(jlen,uspeca(ns),uspn56)

        if (ns .eq. ns1) then
            ! Invert the linking reaction. Put the coefficient for the
            ! ns1-th species first in the range.
            nx = nx + 1

            if (nx .gt. ndrsmx) then
                write (noutpt,1060) ndrsmx,usp156(1:jlen1),usp256(1:jlen2),uspn56(1:jlen),ndrsmx
                write (nttyo,1060) ndrsmx,usp156(1:jlen1),usp256(1:jlen2),uspn56(1:jlen),ndrsmx
1060 format(/' * Error - (EQLIB/swtchb) The maximum ',i7,' entries in the',/7x,'cdrsd and ndrsd arrays has been',' exceeded in trying to',/7x,'make a special basis switch',' replacing ',a,' with',/7x,a,' while writing the reaction',' for ',a,'. Increase',/7x,'the dimensioning parameter',' ndrspa from its present value of ',i7,'.')

                stop
            end if

            cdrsx(nx) = -cx12
            ndrsx(nx) = ns

            do n = nrl1,nrl2
                if (ndrsd(n) .ne. ns) then
                    nx = nx + 1

                    if (nx .gt. ndrsmx) then
                        write (noutpt,1060) ndrsmx,usp156(1:jlen1),usp256(1:jlen2),uspn56(1:jlen),ndrsmx
                        write (nttyo,1060) ndrsmx,usp156(1:jlen1),usp256(1:jlen2),uspn56(1:jlen),ndrsmx
                        stop
                    end if

                    cdrsx(nx) = -cdrsd(n)
                    ndrsx(nx) = ndrsd(n)
                end if
            end do

            ndrsrx(2,ns) = nx

            ! Log K coefficients.
            do j = 1,ntprmx
                do i = 1,narxmx
                    axlksx(i,j,ns) = -axlksd(i,j,ns2)
                end do
            end do

            if (ipch .ge. 0) then
                ! Enthalpy function coefficients.
                do j = 1,ntprmx
                    do i = 1,narxmx
                        axhfsx(i,j,ns) = -axhfsd(i,j,ns2)
                    end do
                end do

                do ipc = 1,ipch
                    do j = 1,ntprmx
                        do i = 1,narxmx
                            adhfsx(i,j,ipc,ns) = -adhfsd(i,j,ipc,ns2)
                        end do
                    end do
                end do
            end if

            if (ipcv .ge. 0) then
                ! Volume function coefficients.
                do j = 1,ntprmx
                    do i = 1,narxmx
                        axvfsx(i,j,ns) = -axvfsd(i,j,ns2)
                    end do
                end do

                do ipc = 1,ipcv
                    do j = 1,ntprmx
                        do i = 1,narxmx
                            advfsx(i,j,ipc,ns) = -advfsd(i,j,ipc,ns2)
                        end do
                    end do
                end do
            end if
        else if (ns .eq. ns2) then
            ! Write a null reaction.
            nx = nx + 1

            if (nx .gt. ndrsmx) then
                write (noutpt,1060) ndrsmx,usp156(1:jlen1),usp256(1:jlen2),uspn56(1:jlen),ndrsmx
                write (nttyo,1060) ndrsmx,usp156(1:jlen1),usp256(1:jlen2),uspn56(1:jlen),ndrsmx
                stop
            end if

            cdrsx(nx) = 0.
            ndrsx(nx) = 0
            ndrsrx(2,ns) = nx

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
                        axx = axhfsd(i,j,ns)

                        do n = nr1 + 1,nr2
                            nsi = ndrsd(n)
                            cx = cdrsd(n)
                            axx = axx - cx*axhfsd(i,j,nsi)
                        end do

                        cx = cdrsd(nr1)
                        axhfsx(i,j,ns) = axx/cx
                    end do
                end do

                do ipc = 1,ipch
                    do j = 1,ntprmx
                        do i = 1,narxmx
                            axx = adhfsd(i,j,ipc,ns)

                            do n = nr1 + 1,nr2
                                nsi = ndrsd(n)
                                cx = cdrsd(n)
                                axx = axx - cx*adhfsd(i,j,ipc,nsi)
                            end do

                            cx = cdrsd(nr1)
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
                        axx = axvfsd(i,j,ns)

                        do n = nr1 + 1,nr2
                            nsi = ndrsd(n)
                            cx = cdrsd(n)
                            axx = axx - cx*axvfsd(i,j,nsi)
                        end do

                        cx = cdrsd(nr1)
                        axvfsx(i,j,ns) = axx/cx
                    end do
                end do

                do ipc = 1,ipcv
                    do j = 1,ntprmx
                        do i = 1,narxmx
                            axx = advfsd(i,j,ipc,ns)

                            do n = nr1 + 1,nr2
                                nsi = ndrsd(n)
                                cx = cdrsd(n)
                                axx = axx - cx*advfsd(i,j,ipc,nsi)
                            end do

                            cx = cdrsd(nr1)
                            advfsx(i,j,ipc,ns) = axx/cx
                        end do
                    end do
                end do
            end if
        else
            ! Copy as is reactions for all other species.
            do n = nr1,nr2
                nx = nx + 1

                if (nx .gt. ndrsmx) then
                    write (noutpt,1060) ndrsmx,usp156(1:jlen1),usp256(1:jlen2),uspn56(1:jlen),ndrsmx
                    write (nttyo,1060) ndrsmx,usp156(1:jlen1),usp256(1:jlen2),uspn56(1:jlen),ndrsmx
                    stop
                end if

                cdrsx(nx) = cdrsd(n)
                ndrsx(nx) = ndrsd(n)
            end do

            ndrsrx(2,ns) = nx

            ! Log K coefficients.
            do j = 1,ntprmx
                do i = 1,narxmx
                    axlksx(i,j,ns) = axlksd(i,j,ns)
                end do
            end do

            if (ipch .ge. 0) then
                ! Enthalpy function coefficients.
                do j = 1,ntprmx
                    do i = 1,narxmx
                        axhfsx(i,j,ns) = axhfsd(i,j,ns)
                    end do
                end do

                do ipc = 1,ipch
                    do j = 1,ntprmx
                        do i = 1,narxmx
                            adhfsx(i,j,ipc,ns) = adhfsd(i,j,ipc,ns)
                        end do
                    end do
                end do
            end if

            if (ipcv .ge. 0) then
                ! Volume function coefficients.
                do j = 1,ntprmx
                    do i = 1,narxmx
                        axvfsx(i,j,ns) = axvfsd(i,j,ns)
                    end do
                end do

                do ipc = 1,ipcv
                    do j = 1,ntprmx
                        do i = 1,narxmx
                            advfsx(i,j,ipc,ns) = advfsd(i,j,ipc,ns)
                        end do
                    end do
                end do
            end if
        end if
    end do

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

    ! Print the new linking reaction.
    ! Calling sequence substitutions:
    !   cdrsd for cdrs
    !   ndrsd for ndrs
    !   ndrsrd for ndrsr
    !   ns1 for ns
    !   noutpt for nf
    !   uspeca for uspec
    call prreac(cdrsd,ndrsd,ndrsmx,ndrsrd,noutpt,ns1,nstmax,uspeca)

    ! Interchange the basis indices.
    nbaspd(nb1) = ns2
    nbaspd(nb2) = ns1

    if (nb1 .eq. nbw) then
        nbw = nb2
    else if (nb2 .eq. nbw) then
        nbw = nb1
    end if

999 continue
end subroutine swtchb