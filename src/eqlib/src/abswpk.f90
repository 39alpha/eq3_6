subroutine abswpk(beta,cdrs,csts,efac,ibswx,iebal,iindx1,jcsort,jflag,jssort,kbt,kmax,mosp,narn1,narn2,nbasp,nbaspd,nbt,nbtmax,ndrs,ndrsmx,ndrsr,nelect,nhydr,no2gaq,nstmax,nsts,nstsmx,nstsr,qbswx,q6mode,weight)
    !! This subroutine determines the dominant species and the associated
    !! factors required for continued fraction corrections.
    !! This subroutine is called by:
    !!   EQLIB/absswa.f
    !! Principal input:
    !!   beta   = array of normalized Newton-Raphson residual functions
    !!   q6mode = flag denoting usage for EQ3NR or EQ6:
    !!              .false. = EQ3NR
    !!              .true.  = EQ6
    !! Principal output:
    !!   ibswx  = array of species indices of candidates for switching
    !!              into the active basis set
    !!   qbswx  = flag denoting that candidates were found for
    !!              basis switching
    implicit none

    ! Calling sequence variable declarations.
    integer :: kmax
    integer :: nbtmax
    integer :: ndrsmx
    integer :: nstmax
    integer :: nstsmx

    integer :: ibswx(nbtmax)
    integer :: iindx1(kmax)
    integer :: jcsort(nstmax)
    integer :: jflag(nstmax)
    integer :: jssort(nstmax)
    integer :: nbasp(nbtmax)
    integer :: nbaspd(nbtmax)
    integer :: ndrs(ndrsmx)
    integer :: ndrsr(2,nstmax)
    integer :: nsts(nstsmx)
    integer :: nstsr(2,nstmax)

    integer :: iebal
    integer :: kbt
    integer :: narn1
    integer :: narn2
    integer :: nbt
    integer :: nelect
    integer :: nhydr
    integer :: no2gaq

    logical :: qbswx
    logical :: q6mode

    real(kind=8) :: beta(kmax)
    real(kind=8) :: cdrs(ndrsmx)
    real(kind=8) :: csts(nstsmx)
    real(kind=8) :: efac(nbtmax)
    real(kind=8) :: mosp(nstmax)
    real(kind=8) :: weight(nstmax)

    real(kind=8) :: coefdr
    real(kind=8) :: coefst

    ! Local variable declarations.
    integer :: ix
    integer :: iy
    integer :: krow
    integer :: nb
    integer :: ns
    integer :: nse
    integer :: nsi
    integer :: nsj

    logical :: qskip

    real(kind=8) :: btest
    real(kind=8) :: cx
    real(kind=8) :: wsi

    ! Note: the following statements don't really do anything except
    ! cause the compiler not to complain that nbt and jssort
    ! are not used.
    ix = nbt
    iy = ix
    nbt = iy

    ix = jssort(1)
    iy = ix
    jssort(1) = iy

    ! Set flag indicating that candidates for basis switching exist.
    qbswx = .false.

    ! Clear the ibswx array.
    call initiz(ibswx,nbtmax)

    ! Loop over all mass balance relations.
    do krow = 1,kbt
        nb = iindx1(krow)
        nse = nbaspd(nb)
        nsj = nbasp(nb)

        if (nse.ne.narn1 .and. nse.ne.nhydr .and. nb.ne.iebal .and.    nse.ne.no2gaq .and.nse.ne.nelect) then
            ! If the current species is not H2O, H+, O2(g,aq), or e-,
            ! and is not being used for electrical balancing in EQ3NR,
            ! consider a switch.
            ! Pick a candidate for automatic basis switching if a
            ! continued fraction correction would be large. Store
            ! the index of the candidate, if any, in the array ibswx.
            btest = (beta(krow) + 1.)**efac(nb)

            if (btest .ge. 10.) then
                do ns = narn1,narn2
                    weight(ns) = coefst(csts,nsts,nstsmx,nstsr,nb,ns,nstmax)
                end do

                ! Screen out species in the mass balance that are not
                ! linked by current reaction with the current basis species.
                ! Do this by setting the corresonding weights to zero.
                do ns = narn1,narn2
                    if (jflag(ns) .eq. 30) then
                        ! Calling sequence substitutions:
                        !   nsj for nse
                        cx = coefdr(cdrs,ndrs,ndrsmx,ndrsr,nsj,ns,nstmax)

                        if (cx .eq. 0.) then
                            weight(ns) = 0.
                        end if
                    end if
                end do

                call fbassw(jcsort,jflag,mosp,narn1,narn2,nse,nsi,nsj,nstmax,weight,wsi)

                if (nsi .gt. 0) then
                    ibswx(nb) = nsi
                    qbswx = .true.
                end if
            end if
        end if
    end do
end subroutine abswpk