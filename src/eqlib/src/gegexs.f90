subroutine gegexs(cegexs,cgexj,egexjc,egexjf,egexs,iern1,iern2,ietmax,jern1,jetmax,jgext,moph,mosp,mrgexs,netmax,ngexsa,ngext,noutpt,nptmax,nstmax,nttyo,zchar,zgexj)
    !! This subroutine computes the equivalent fractions (egexs) and
    !! mole ratios (mrgexs) of exchanger species belonging to generic
    !! ion exchange phases.
    !! The equivalent fractions for the exchanger species are the same
    !! as those for the corresponding exchange ions. The mole ratios
    !! are the number of moles of exchanger species per mole of
    !! the exchange phase. The number of moles of an exchanger species
    !! may or may not equal the number of moles of the corresponding
    !! exchange ion. Thus, these mole ratios do not necessarily give
    !! the number of moles of the exchange ions per mole of the exchange
    !! phase.
    !! The calculated equivalent fractions are based on formally
    !! declared exchange capacities (egexjf). Thus, they need not sum
    !! to unity, though they will do so for some models (e.g., Gapon,
    !! Vanselow). If an anion should exchange onto a negatively charged
    !! site, it will have a negative equivalent fraction.
    !! Subroutine EQLIB/gegexw.f computes the apparent "whole-phase"
    !! equivalent fractions (egexw) of cations and anions on generic ion
    !! exchangers.
    !! This subroutine is called by:
    !!   EQLIB/ncmpex.f
    !! Principal input:
    !!   moph   = array of numbers of moles of phases
    !!   mosp   = array of numbers of moles of species
    !! Principal output:
    !!   egexs  = array of equivalent fractions, by site, based on
    !!              formally declared exchange capacities
    !!   mrgexs = array of mole ratios
    implicit none

    ! Calling sequence variable declarations.
    integer :: ietmax
    integer :: jetmax
    integer :: netmax
    integer :: nptmax
    integer :: nstmax

    integer :: noutpt
    integer :: nttyo

    integer :: jern1(jetmax,netmax)
    integer :: jgext(netmax)
    integer :: ngexsa(ietmax,jetmax,netmax)
    integer :: ngext(jetmax,netmax)

    integer :: iern1
    integer :: iern2

    real(kind=8) :: cegexs(ietmax,jetmax,netmax)
    real(kind=8) :: cgexj(jetmax,netmax)
    real(kind=8) :: egexjc(jetmax,netmax)
    real(kind=8) :: egexjf(jetmax,netmax)
    real(kind=8) :: egexs(ietmax,jetmax,netmax)
    real(kind=8) :: mrgexs(ietmax,jetmax,netmax)
    real(kind=8) :: moph(nptmax)
    real(kind=8) :: mosp(nstmax)
    real(kind=8) :: zchar(nstmax)
    real(kind=8) :: zgexj(jetmax,netmax)

    ! Local variable declarations.
    integer :: ie
    integer :: je
    integer :: ne
    integer :: np
    integer :: ns
    integer :: nss

    real(kind=8) :: ex
    real(kind=8) :: exf
    real(kind=8) :: exsc

    ! Compute the equivalent fractions (egexs) for the sites of
    ! the generic ion exchangers. Note that these fractions are
    ! based on the formal declared exchange capacities (egexjf) of
    ! these sites. These equivalent fractions will only sum to unity
    ! for certain models, such as the Vanselow and Gapon models,
    ! which require the actual capacity (egexjc) to equal the
    ! corresponding formal capacity.
    do np = iern1,iern2
        ne = np - iern1 + 1

        do je = 1,jgext(ne)
            ! First calculate the number of equivalents for each
            ! component species from the corresponding number of moles.
            ! Also calculate the total number of equivalents (egexjc)
            ! actually on the current site of the current exchanger. This
            ! may or may not equal the formal exchange capacity of the site
            ! egexjf), depending on the exchange model. It will equal it in
            ! the case of simple exchange models which allow exchange of
            ! ions of only one charge sign on a given site. The Vanselow
            ! and Gapon models satisfy this condition.
            exsc = 0.
            ns = jern1(je,ne) - 1

            do ie = 1,ngext(je,ne)
                ns = ns + 1
                ex = mosp(ns)*cegexs(ie,je,ne)
                egexs(ie,je,ne) = ex
                exsc = exsc + ex
            end do

            ! Note that the sign of the capacity is the opposite of that
            ! of the charge on the site itself.
            egexjc(je,ne) = exsc/moph(np)

            if (egexjf(je,ne) .ne. 0.) then
                ! Calculate the equivalent fractions.
                exf = egexjf(je,ne)*moph(np)

                do ie = 1,ngext(je,ne)
                    egexs(ie,je,ne) = egexs(ie,je,ne)/exf
                end do
            else
                ! If the declared formal capacity is zero (not permitted
                ! for the Vanselow and Gapon models), do not calculate
                ! equivalent fractions. Set the corresponding elements of
                ! the egexs array to zero.
                do ie = 1,ngext(je,ne)
                    egexs(ie,je,ne) = 0.
                end do
            end if
        end do
    end do

    ! Calculate the mole ratios (mrgexs) for exchanger species of
    ! of generic ion exchanger reactants. The mole ratio is the
    ! moles of exchanger species (specific to a site) per mole of
    ! exchanger substrate.
    do np = iern1,iern2
        ne = np - iern1 + 1

        do je = 1,jgext(ne)
            ex = zgexj(je,ne)*cgexj(je,ne)

            do ie = 1,ngext(je,ne)
                nss = ngexsa(ie,je,ne)

                if (nss .gt. 0) then
                    mrgexs(ie,je,ne) = -ex*egexs(ie,je,ne)/zchar(nss)
                end if
            end do
        end do
    end do
end subroutine gegexs