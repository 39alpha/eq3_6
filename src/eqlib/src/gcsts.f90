subroutine gcsts(cdrs,csts,jflag,nbaspd,nbt,nbtmax,ndrs,ndrsmx,ndrsr,noutpt,nsts,nstsmx,nstsr,nst,nstmax,nttyo,uspec)
    !! This subroutine computes the stoichiometric factors which relate
    !! each aqueous species to the nb-th member of the data file
    !! basis set. These stoichiometric factors permit calculation of
    !! sums which correspond to physically meaningful 'total' masses
    !! or concentrations of the basis species, except for three of
    !! these species. It is not possible to define physically
    !! meaningful mass balances for for water, hydrogen ion, and the
    !! aqueous species oxygen gas. The reactions input to this subroutine
    !! may be rewritten from the data file forms to reflect the
    !! elmination of one or more auxiliary basis species from the
    !! active basis set. They may not be rewritten to reflect basis
    !! switching, except for a switch which exchanges a strict basis
    !! species with an auxiliary basis species.
    !! This subroutine is called by:
    !!   EQ3NR/eq3nr.f
    !!   EQ6/eq6.f
    !! Principal input:
    !!   cdrs   = array of reaction coefficients
    !!   nbaspd = indices of the species in the 'd' basis set
    !!   ndrs   = array of indices of the species corresponding to the
    !!              coefficients in the cdrs array
    !!   ndrsr  = array giving the range in the cdrs and ndrs arrays
    !!              corresponding to the reaction for a given species
    !!   uspec  = array of species names
    !! Principal output:
    !!   csts   = array of stoichiometric coefficients appearing in
    !!              mass balance relations
    !!   nsts   = array of indices of the basis species corresponding
    !!              to the coefficients in the csts array
    !!   nstsr  = array giving the range in the csts and nsts arrays
    !!              corresponding to a given species
    implicit none

    ! Calling sequence variable declarations.
    integer :: nbtmax
    integer :: ndrsmx
    integer :: nstsmx
    integer :: nstmax

    integer :: noutpt
    integer :: nttyo

    integer :: jflag(nstmax)
    integer :: nbaspd(nbtmax)
    integer :: ndrs(ndrsmx)
    integer :: ndrsr(2,nstmax)
    integer :: nsts(nstsmx)
    integer :: nstsr(2,nstmax)
    integer :: nbt
    integer :: nst

    character(len=48) :: uspec(nstmax)

    real(kind=8) :: cdrs(ndrsmx)
    real(kind=8) :: csts(nstsmx)

    ! Local variable declarations.
    integer :: jlen
    integer :: jlene
    integer :: n
    integer :: nb
    integer :: nbb
    integer :: nerr
    integer :: nj
    integer :: nn
    integer :: nr1
    integer :: nr2
    integer :: ns
    integer :: nse
    integer :: nt
    integer :: nts

    integer :: nbasis

    logical :: qbasis

    character(len=56) :: uspe56
    character(len=56) :: uspn56

    real(kind=8) :: cxs

    nerr = 0
    n = 0

    do ns = 1,nst
        ! Set first element of the range pointer array.
        nstsr(1,ns) = n + 1

        nr1 = ndrsr(1,ns)
        nr2 = ndrsr(2,ns)
        nt = nr2 - nr1 + 1
        nts = nt - 1

        ! Check to see if the current species is a member of the active
        ! basis set. If so, set qbasis to .true. and nts to 1.
        ! Calling sequence substitutions:
        !   nbaspd for nbasp
        nb = nbasis(nbaspd,nbt,nbtmax,ns)

        qbasis = .false.

        if (nb.gt.0 .and. jflag(ns).ne.30) then
            qbasis = .true.
            nts = 1
        end if

        ! Check to see if the csts/nsts array size is sufficient.
        nj = n + nts

        if (nj .gt. nstsmx) then
            ! Calling sequence substitutions:
            !   uspec(ns) for unam48
            call fmspnx(jlen,uspec(ns),uspn56)
            write (noutpt,1000) nstsmx,uspn56(1:jlen)
            write (nttyo,1000) nstsmx,uspn56(1:jlen)
1000 format(/' * Error - (EQLIB/gcsts) The maximum ',i7,' entries have been exceeded',/7x,'computing the csts',' array of stoichiometric coefficients. The last',/7x,'species for which the coefficients were being computed',/7x,'was ',a,'. Increase the dimensioning',/7x,'parameter nstspa')

            stop
        end if

        if (qbasis) then
            ! Set mass balance coefficient for the ns-th species if it is
            ! in the active basis set.
            n = n + 1
            nsts(n) = nb
            csts(n) = 1.
        else
            ! Set mass balance coefficients for other species appearing in
            ! the reaction for the ns-th species if the ns-th species is
            ! not in the active basis set.
            if (nt .lt. 2) then
                ! Calling sequence substitutions:
                !   uspec(ns) for unam48
                call fmspnx(jlen,uspec(ns),uspn56)
                write (noutpt,1005) uspn56(1:jlen)
                write (nttyo,1005) uspn56(1:jlen)
1005 format(/' * Error - (EQLIB/gcsts) The species ',a,/7x,'has no associated reaction on the data file, but it',/7x,'is not a strict basis species.')

                nerr = nerr + 1

                ! Make a single null entry.
                n = n + 1
                nsts(n) = 0
                csts(n) = 0.
            else
                cxs = -cdrs(nr1)

                do nn = nr1 + 1,nr2
                    n = n + 1
                    nse = ndrs(nn)

                    ! Calling sequence substitutions:
                    !   nbaspd for nbasp
                    !   nse for ns
                    nbb = nbasis(nbaspd,nbt,nbtmax,nse)

                    if (nbb .le. 0) then
                        ! Calling sequence substitutions:
                        !   uspec(ns) for unam48
                        call fmspnx(jlen,uspec(ns),uspn56)

                        ! Calling sequence substitutions:
                        !   jlene for jlen
                        !   uspec(nse) for unam48
                        !   uspe56 for uspn56
                        call fmspnx(jlene,uspec(nse),uspe56)
                        write (noutpt,1010) uspe56(1:jlene),uspn56(1:jlen)
                        write (nttyo,1010) uspe56(1:jlene),uspn56(1:jlen)
1010 format(/' * Error - (EQLIB/gcsts) The species ',a,/7x,'appears in the data file reaction for ',a,',',/7x,'but it is not in the data file basis set.')

                        nerr = nerr + 1
                    end if

                    nsts(n) = nbb
                    csts(n) = cdrs(nn)/cxs
                end do
            end if
        end if

        ! Set the second element of the range pointer array.
        nstsr(2,ns) = n
    end do

    if (nerr .gt. 0) then
        write (noutpt,1015)
        write (nttyo,1015)
1015 format(/' * Error - (EQLIB/gcsts) One or more reactions',/7x,'on the data file are not consistent with the data file',/7x,'basis set.')

        stop
    end if
end subroutine gcsts