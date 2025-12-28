subroutine tstrdx(cdrs,iodb,iopt,jflag,jsflag,narn1,narn2,ndrs,ndrsmx,ndrsr,nodbmx,noptmx,noutpt,nrdxsp,nstmax,qredox,uspec)
    !! This subroutine determines if the chemical model to be
    !! computed has a redox aspect. This will be determined to be
    !! so if a species in the model has an associated reaction
    !! that is a redox reaction and the species is not in the
    !! active basis set.
    !! An auxiliary basis species (say Oxalate-) that is in the
    !! model but is included in the active basis set is effectively.
    !! treated as detached from a corresponding strict basis species
    !! (e.g., the concentration/activity of Oxalate- is not determined
    !! by the assumption of equilibrium for its associated reaction,
    !! which would link it to HCO3-). In effect, Oxalate- is treated
    !! as being composed of a pseudo-element, and its presence in the
    !! model does not require a redox variable.
    !! This subroutine is called by:
    !!   EQ6/eq6.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: ndrsmx
    integer :: nodbmx
    integer :: noptmx
    integer :: nstmax

    integer :: noutpt

    integer :: iodb(nodbmx)
    integer :: iopt(noptmx)
    integer :: jflag(nstmax)
    integer :: jsflag(nstmax)
    integer :: ndrs(ndrsmx)
    integer :: ndrsr(2,nstmax)

    integer :: narn1
    integer :: narn2
    integer :: nrdxsp

    logical :: qredox

    character(len=48) :: uspec(nstmax)

    real(kind=8) :: cdrs(ndrsmx)

    ! Local variable declarations.
    integer :: jlen
    integer :: ns

    character(len=56) :: uspn56

    real(kind=8) :: cx

    real(kind=8) :: coefdr

    qredox = .false.

    if (iopt(15) .le. 0) then
        do ns = narn1,narn2
            if (jflag(ns) .eq. 30) then
                ! The species is not in the active basis set. It is a
                ! "dependent" species whose concentration/activity is
                ! computed assuming its associated reaction is in a state
                ! of equilibrium.
                if (jsflag(ns) .lt. 2) then
                    ! The species not hard suppressed.
                    ! Calling sequence substitutions:
                    !   nrdxsp for nse
                    cx = coefdr(cdrs,ndrs,ndrsmx,ndrsr,nrdxsp,ns,nstmax)

                    if (cx .ne. 0.) then
                        qredox = .true.

                        if (iodb(1) .ge. 1) then
                            ! Calling sequence substitutions:
                            !   uspec(ns) for unam48
                            call fmspnm(jlen,uspec(ns),uspn56)
                            write (noutpt,1000) uspn56(1:jlen)
1000 format(/' * Note - (EQ6/tstrdx) The reaction',' associated with the species',/7x,a,' is a redox',' reaction.')

                            go to 999
                        end if
                    end if
                end if
            end if
        end do
    end if

999 continue
end subroutine tstrdx