subroutine chsgex(adhfs,adhfsd,adhfsx,advfs,advfsd,advfsx,axhfs,axhfsd,axhfsx,axlks,axlksd,axlksx,axvfs,axvfsd,axvfsx,cdrs,cdrsd,cdrsx,eps100,iern1,ipch,ipchmx,ipcv,ipcvmx,jern1,jetmax,jflag,jgext,jsflag,narn1,narxmx,narxt,nbasp,nbaspd,nbaspx,nbt,nbtmax,nbw,ndrs,ndrsd,ndrsmx,ndrsr,ndrsrd,ndrsrx,ndrsx,nern1,nern2,net,netmax,ngext,noutpt,nphasx,nst,nstmax,ntprmx,ntprt,nttyo,qbassw,qbswok,ugexmo,uspec)
    !! This subroutine changes the original setup of component species
    !! of generic ion exchanger phases for certain exchange models
    !! (e.g., Gapon, Vanselow). It switches bare site species out of
    !! the active basis set; e.g., it replaces __-Z with Na-Z. The
    !! thermodynamic properties and the jsflag status array element
    !! for each such base site species are changed to effectively take
    !! this species out of the model. It persists only as a reference
    !! for the mass balance of the corresponding exchanger substrate.
    !! This subroutine is called by:
    !!   EQ3NR/eq3nr.f
    !!   EQ6/eq6.f
    !! Principal input:
    !!   narxmx = number of coefficient elements of axlks per species
    !!              per temperature range
    !!   nbtmax = maximum number of basis species
    !!   nstmax = maximum number of species
    !!   ntprmx = maximum number of temperature ranges
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: ipchmx
    integer :: ipcvmx
    integer :: jetmax
    integer :: narxmx
    integer :: nbtmax
    integer :: ndrsmx
    integer :: netmax
    integer :: nstmax
    integer :: ntprmx

    integer :: noutpt
    integer :: nttyo

    integer :: jern1(jetmax,netmax)
    integer :: jflag(nstmax)
    integer :: jgext(netmax)
    integer :: jsflag(nstmax)
    integer :: narxt(ntprmx)
    integer :: nbasp(nbtmax)
    integer :: nbaspd(nbtmax)
    integer :: nbaspx(nbtmax)
    integer :: ndrs(ndrsmx)
    integer :: ndrsd(ndrsmx)
    integer :: ndrsr(2,nstmax)
    integer :: ndrsrd(2,nstmax)
    integer :: ndrsx(ndrsmx)
    integer :: ndrsrx(2,nstmax)
    integer :: ngext(jetmax,netmax)
    integer :: nphasx(nstmax)

    integer :: iern1
    integer :: ipch
    integer :: ipcv
    integer :: narn1
    integer :: nbt
    integer :: nbw
    integer :: nern1
    integer :: nern2
    integer :: net
    integer :: nst
    integer :: ntprt

    logical :: qbassw
    logical :: qbswok

    character(len=24) :: ugexmo(netmax)
    character(len=48) :: uspec(nstmax)

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

    real(kind=8) :: eps100

    ! Local variable declarations.
    integer :: i
    integer :: ie
    integer :: j
    integer :: je
    integer :: j2
    integer :: nb
    integer :: ne
    integer :: np
    integer :: ns
    integer :: ns1
    integer :: ns2

    integer :: ilnobl

    ! Copy the existing nbasp set into the nbaspx array.
    do nb = 1,nbt
        nbaspx(nb) = nbasp(nb)
    end do

    ! Switch bare site species of generic ion exchange phases with
    ! certain exchange models (e.g., Gapon, Vanselow) out of the
    ! basis set. E.g., use Na-Z instead of __-Z.
    do nb = 1,nbt
        ns1 = nbaspx(nb)

        if (ns1.ge.nern1 .and. ns1.le.nern2) then
            np = nphasx(ns1)
            ne = np - iern1 + 1
            j2 = ilnobl(ugexmo(ne))

            if (ugexmo(ne)(1:j2).eq.'Gapon' .or.      ugexmo(ne)(1:6).eq.'Gapon-' .or.      ugexmo(ne)(1:j2).eq.'Vanselow' .or.      ugexmo(ne)(1:9).eq.'Vanselow-') then
                ! Theoretically, the test below should be unnecessary,
                ! as the only exchanger species presently in the basis
                ! set should be bare site species.
                if (uspec(ns1)(1:3) .eq. '__ ') then
                    ! Set up to switch the bare exchanger species with the
                    ! next species in the same site.
                    ns2 = ns1 + 1
                    nbasp(nb) = ns2
                end if
            end if
        end if
    end do

    ! Make the switches.
    do nb = 1,nbt
        ns1 = nbaspx(nb)
        ns2 = nbasp(nb)

        if (ns1 .ne. ns2) then
            call switch(adhfs,adhfsx,advfs,advfsx,axhfs,axhfsx,axlks,axlksx,axvfs,axvfsx,cdrs,cdrsx,eps100,ipch,ipchmx,ipcv,ipcvmx,jflag,jsflag,narn1,narxmx,nbasp,nbaspd,nbaspx,nb,nbt,nbtmax,nbw,ndrs,ndrsmx,ndrsx,ndrsr,ndrsrx,noutpt,ns2,nst,nstmax,ntprmx,nttyo,qbassw,qbswok,uspec)
        end if
    end do

    ! Now for certain models (e.g., Gapon, Vanselow) adjust the
    ! thermodynamic data for the bare site species so that these
    ! species are effectively taken out of the model.
    do ne = 1,net
        j2 = ilnobl(ugexmo(ne))

        if (ugexmo(ne)(1:j2).eq.'Gapon' .or.    ugexmo(ne)(1:6).eq.'Gapon-' .or.    ugexmo(ne)(1:j2).eq.'Vanselow' .or.    ugexmo(ne)(1:9).eq.'Vanselow-') then
            do je = 1,jgext(ne)
                ns = jern1(je,ne) - 1

                do ie = 1,ngext(je,ne)
                    ns = ns + 1

                    if (uspec(ns)(1:3) .eq. '__ ') then
                        ! Null the thermodynamic properties of the
                        ! bare site species.
                        do j = 1,ntprt
                            axlks(1,j,ns) = +9999999.
                            axhfs(1,j,ns) = 0.
                            axvfs(1,j,ns) = 0.

                            do i = 2,narxt(j)
                                axlks(i,j,ns) = 0.
                                axhfs(i,j,ns) = 0.
                                axvfs(i,j,ns) = 0.
                            end do
                        end do

                        jsflag(ns) = 2
                    end if
                end do
            end do
        end if
    end do

    ! Make the 'd' set consistent with the ion exchanger part of the
    ! ordinary set as it now exists. Do not modify other parts of the
    ! 'd' set. This action will leave, say, Na-Z in the 'd' set in
    ! place of __-Z, which was the basis species on which the
    ! exchanger substrate mass balance was defined. The only memory
    ! of this definition is contained in the nbaspi array. Note
    ! that the 'd' set reactions for exchangers will now include
    ! the effects of eliminations of other basis species from the
    ! active basis set.
    ! The effects of eliminations on exchanger reactions could be
    ! overcome by saving the exchanger reactions and associated data
    ! prior to making the eliminations. However, new arrays would be
    ! required to do that. It doesn't seem worthwhile at the present
    ! time.
    call mdrgex(adhfs,adhfsd,adhfsx,advfs,advfsd,advfsx,axhfs,axhfsd,axhfsx,axlks,axlksd,axlksx,axvfs,axvfsd,axvfsx,cdrs,cdrsd,cdrsx,ipch,ipchmx,ipcv,ipcvmx,narxmx,nbasp,nbaspd,nbtmax,ndrs,ndrsd,ndrsx,ndrsmx,ndrsr,ndrsrd,ndrsrx,nern1,nern2,noutpt,nst,nstmax,ntprmx,nttyo)

999 continue
end subroutine chsgex