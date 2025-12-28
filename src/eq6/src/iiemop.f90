subroutine iiemop(iemop,iemos,iindx1,ipndx1,jsflag,kdim,kmax,ncmpe,ncmpr,noutpt,npet,npetmx,npt,nptmax,nset,nsetmx,nstmax,nttyo,uaqsln,uspec,uphase)
    !! This subroutine initializes the indexing used for tracking the
    !! numbers of moles of phases present in the Equilibrium System (ES).
    !! This indexing must be re-initialized whenever a change occurs
    !! in the ES phase assemblage. The iemop array contains the indices
    !! of the phases in the ES. The iemos array contains the indicies
    !! of the active species in these phases, except that only the
    !! species H2O(l) is included in the case of the aqueous solution
    !! phase. solution. The ncmpe array is a species range pointer
    !! array for the phases analogous to ncmpr.
    !! The iemop and emop arrays respectively contain the indices and
    !! numbers of moles of the phases currently in the ES. The iemos
    !! and emos array are the analogs for the species of these phases.
    !! However, only the species H2O(l) is represented for the aqueous
    !! solution.  The ncmpe array is a species range pointer array for
    !! the phases analogous to ncmpr.
    !! The fdpe0 and fdse0 arrays respectively contain the finite
    !! differences for the data in the emop and emos arrays. The demop
    !! and demos arrays contain the corresponding derivatives with
    !! respect to the reaction progress variable, xi.
    !! This subroutine is called by:
    !!   EQ6/dumpdp.f
    !!   EQ6/path.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: kmax
    integer :: npetmx
    integer :: nptmax
    integer :: nsetmx
    integer :: nstmax

    integer :: noutpt
    integer :: nttyo

    integer :: iemop(npetmx)
    integer :: iemos(nsetmx)
    integer :: iindx1(kmax)
    integer :: ipndx1(kmax)
    integer :: jsflag(nstmax)
    integer :: ncmpe(2,npetmx)
    integer :: ncmpr(2,nptmax)

    integer :: kdim
    integer :: npet
    integer :: npt
    integer :: nset

    character(len=48) :: uspec(nstmax)
    character(len=24) :: uphase(nptmax)

    character(len=24) :: uaqsln

    ! Local variable declarations.
    integer :: jlen
    integer :: j2
    integer :: kcol
    integer :: np
    integer :: npe
    integer :: nplast
    integer :: nr1
    integer :: nr2
    integer :: ns
    integer :: nse

    integer :: ilnobl

    character(len=56) :: uspn56
    character(len=8) :: ux8

    ! Loop over all phases present in the ES.
    npe = 0
    nse = 0
    nplast = 0

    do kcol = 1,kdim
        np = ipndx1(kcol)

        if (np .ne. nplast) then
            npe = npe + 1
            iemop(npe) = np
            nplast = np
            nr1 = ncmpr(1,np)
            nr2 = ncmpr(2,np)

            if (uphase(np)(1:24) .eq. uaqsln(1:24)) then
                nr2 = nr1
            end if

            ncmpe(1,npe) = nse + 1

            do ns = nr1,nr2
                if (jsflag(ns) .le. 0) then
                    nse = nse + 1

                    if (nse .gt. nsetmx) then
                        ! Calling sequence substitutions:
                        !   uspec(ns) for unam48
                        call fmspnm(jlen,uspec(ns),uspn56)
                        write (ux8,'(i5)') nsetmx
                        call lejust(ux8)
                        j2 = ilnobl(ux8)
                        write (noutpt,1000) ux8(1:j2),uspn56(1:jlen)
                        write (nttyo,1000) ux8(1:j2),uspn56(1:jlen)
1000 format(/' * Error - (EQ6/iiemop) Have too many species',' to track by finite differences.',/7x,'Exceeded the',' dimensioned limit of ',a,' upon trying to set up',/7x,'tracking for ',a,'. Increase the dimensioning',/7x,'parameter nsetpa.')
                    end if

                    iemos(nse) = ns
                end if
            end do

            ncmpe(2,npe) = nse
        end if
    end do

    npet = npe
    nset = nse
end subroutine iiemop