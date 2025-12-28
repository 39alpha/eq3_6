subroutine gegexw(cegexs,egexpc,egexpa,egexw,iern1,iern2,ietmax,jern1,jetmax,jgext,kern1,kern2,ketmax,kgexsa,moph,mosp,netmax,ngexsa,ngext,noutpt,nptmax,nstmax,nttyo,xgexw,zchar)
    !! This subroutine computes the apparent "whole-phase" equivalent
    !! fractions (egexw) and mole fractions (xgexw) of the exchange
    !! ions present in generic ion exchanger phases. The exchange
    !! fractions are calculated separately for cations and anions.
    !! For example, the exchange fraction of a cation is the number
    !! of equivalents of that cation divided by the sum of the
    !! equivalents of all cations in the same exchange phase. The
    !! mole fractions are calculated in similar fashion.
    !! Subroutine EQLIB/gegexs.f computes the equivalent fractions
    !! (egexs) and mole ratios (mrgexs) of exchanger species (moles
    !! of exchanger species per mole of exchanger phase). Those
    !! equivalent fractions are defined for individual sites.
    !! This subroutine is called by:
    !!   EQ3NR/eq3nr.f
    !!   EQ6/path.f
    !! Principal input:
    !!   moph   = array of numbers of moles of phases
    !!   mosp   = array of numbers of moles of species
    !! Principal output:
    !!   egexw  = array of apparent "whole-phase" equivalent fractions,
    !!              for exchange ions
    !!   egexpc = array of apparent cationic exchange capacities
    !!   egexpa = array of apparent anionic exchange capacities
    !!   xgexw  = array of apparent "whole-phase" mole fractions,
    !!              for exchange ions
    implicit none

    ! Calling sequence variable declarations.
    integer :: ietmax
    integer :: jetmax
    integer :: ketmax
    integer :: netmax
    integer :: nptmax
    integer :: nstmax

    integer :: noutpt
    integer :: nttyo

    integer :: jern1(jetmax,netmax)
    integer :: jgext(netmax)
    integer :: kern1(netmax)
    integer :: kern2(netmax)
    integer :: kgexsa(ketmax,netmax)
    integer :: ngext(jetmax,netmax)
    integer :: ngexsa(ietmax,jetmax,netmax)

    integer :: iern1
    integer :: iern2

    real(kind=8) :: cegexs(ietmax,jetmax,netmax)
    real(kind=8) :: egexpa(netmax)
    real(kind=8) :: egexpc(netmax)
    real(kind=8) :: egexw(ketmax,netmax)
    real(kind=8) :: moph(nptmax)
    real(kind=8) :: mosp(nstmax)
    real(kind=8) :: xgexw(ketmax,netmax)
    real(kind=8) :: zchar(nstmax)

    ! Local variable declarations.
    integer :: ie
    integer :: je
    integer :: ke
    integer :: ne
    integer :: np
    integer :: ns
    integer :: nss

    real(kind=8) :: ex
    real(kind=8) :: expac
    real(kind=8) :: expcc
    real(kind=8) :: mx
    real(kind=8) :: mxpac
    real(kind=8) :: mxpcc
    real(kind=8) :: zx

    ! Compute the apparent "whole-phase" equivalent fractions (egexw)
    ! and mole fractions (xgexw) for the exchange species in generic
    ! ion exchangers. Note that for a cationic exchange species
    ! (e.g., Na+), this equivalent fraction is defined as the total
    ! number of equivalents of the species on the exchanger (in all
    ! sites, regardless of charge) divided by the sum of the total
    ! number of equivalents of all such cationic species (expcc).
    ! For an anionic species, this equivalent fraction is defined
    ! analogously (expac being analogous to expcc). The formally
    ! declared exchange capacities for the sites and the phases are
    ! not used here.
    do np = iern1,iern2
        ne = np - iern1 + 1

        ! Loop on exchanging species (e.g., Na+). Get the total
        ! number of equivalents of each such species on the phase
        ! by summing over contributions from all sites.
        do ke = kern1(ne),kern2(ne)
            egexw(ke,ne) = 0.
            nss = kgexsa(ke,ne)

            if (nss .gt. 0) then
                ! Loop on sites.
                do je = 1,jgext(ne)
                    ns = jern1(je,ne) - 1

                    ! Loop on exchanger species (e.g., Na-Z).
                    do ie = 1,ngext(je,ne)
                        ns = ns + 1

                        if (nss .eq. ngexsa(ie,je,ne)) then
                            ! The current exchanger species (e.g., Na-Z) is
                            ! comprised of the current exchanging species
                            ! (e.g., Na+).
                            ex = mosp(ns)*cegexs(ie,je,ne)
                            egexw(ke,ne) = egexw(ke,ne) + ex
                        end if
                    end do
                end do
            end if
        end do

        ! Get the apparent cationic and anionic "whole-phase" exchange
        ! capacities (egexpc and egexpa, respectively).
        expcc = 0.
        expac= 0.

        do ke = kern1(ne),kern2(ne)
            nss = kgexsa(ke,ne)

            if (nss .gt. 0) then
                ex = egexw(ke,ne)

                if (ex .gt. 0.) then
                    expcc = expcc + ex
                else if (ex .lt. 0.) then
                    expac = expac + ex
                end if
            end if
        end do

        egexpc(ne) = expcc/moph(np)
        egexpa(ne) = expac/moph(np)

        ! Calculate the number of moles of exchange ions in the exchange
        ! phase from the corresponding number of equivalents. Get the
        ! total number of moles of exchange cations (mxpcc) and of
        ! exchange anions (mxpac).
        mxpcc = 0.
        mxpac= 0.

        do ke = kern1(ne),kern2(ne)
            nss = kgexsa(ke,ne)

            if (nss .gt. 0) then
                zx = zchar(nss)

                if (zx .ne. 0.) then
                    mx = egexw(ke,ne)/zx
                else
                    mx = 0.
                end if

                xgexw(ke,ne) = mx

                if (ex .gt. 0.) then
                    mxpcc = mxpcc + mx
                else if (ex .lt. 0.) then
                    mxpac = mxpac + mx
                end if
            end if
        end do

        ! Loop on exchanging species (e.g., Na+). This time calculate
        ! the apparent "whole-phase" equivalent fractions.
        do ke = kern1(ne),kern2(ne)
            nss = kgexsa(ke,ne)

            if (nss .gt. 0) then
                ex = egexw(ke,ne)

                if (ex .gt. 0.) then
                    egexw(ke,ne) = ex/expcc
                else if (ex .lt. 0.) then
                    egexw(ke,ne) = ex/expac
                end if
            end if
        end do

        ! Loop on exchanging species (e.g., Na+). This time calculate
        ! the apparent "whole-phase" mole fractions.
        do ke = kern1(ne),kern2(ne)
            nss = kgexsa(ke,ne)

            if (nss .gt. 0) then
                zx = zchar(nss)
                mx = xgexw(ke,ne)

                if (zx .gt. 0.) then
                    xgexw(ke,ne) = mx/mxpcc
                else if (zx .lt. 0.) then
                    xgexw(ke,ne) = mx/mxpac
                end if
            end if
        end do
    end do
end subroutine gegexw