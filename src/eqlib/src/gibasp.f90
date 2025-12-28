subroutine gibasp(cgexj,cjbasp,iern1,ixbasp,jern1,jern2,jetmax,jgext,narn1,narn2,nbasp,nbt,nbtmax,nern1,nern2,netmax,nphasx,nstmax)
    !! This subroutine sets up the ixbasp and cjbasp arrays. The former
    !! is a flag array, each member of which denotes whether the
    !! thermodynamic activity of the corresponding basis species is
    !! defined in terms of molality (= 0) or mole fraction (= 1). The
    !! cjbasp array contains any site stoichiometric factor associated
    !! with a given basis species. If there is no such factor, then
    !! an element of cjbasp may be assigned a value of unity.
    !! This subroutine is called by:
    !!   EQ3NR/eq3nr.f
    !!   EQ3NR/arrset.f
    !!   EQ6/eqphas.f
    !!   EQ6/eqshel.f
    !!   EQ6/optmzr.f
    !!   EQ6/path.f
    !! Principal input:
    !!   cgexj  = array of site stoichiometric factors for generic ion
    !!              exchanger species
    !!   jgext  = array giving the number of sites for a generic ion
    !!              exchanger species.
    !!   nbasp  = array containing the species indices of the basis
    !!              species
    !!   nbt    = the number of basis species
    !!   nphasx = array giving the phase to which a species belongs
    !! Principal output:
    !!   cjbasp = array of site stoichiometric factors corresponding to
    !!              sites occupied by basis species, if any
    !!   ixbasp = array of flag switches indicating whether the
    !!              thermodynamic activity of a basis species is based
    !!              on molality (= 0) or mole fraction (= 1)
    implicit none

    ! Calling sequence variable declarations.
    integer :: jetmax
    integer :: nbtmax
    integer :: netmax
    integer :: nstmax

    integer :: ixbasp(nbtmax)
    integer :: jern1(jetmax,netmax)
    integer :: jern2(jetmax,netmax)
    integer :: jgext(netmax)
    integer :: nbasp(nbtmax)
    integer :: nphasx(nstmax)

    integer :: iern1
    integer :: narn1
    integer :: narn2
    integer :: nbt
    integer :: nern1
    integer :: nern2

    real(kind=8) :: cjbasp(nbtmax)
    real(kind=8) :: cgexj(jetmax,netmax)

    ! Local variable declarations.
    integer :: je
    integer :: nb
    integer :: ne
    integer :: np
    integer :: nrr1
    integer :: nrr2
    integer :: ns
    integer :: nss

    do nb = 1,nbt
        ixbasp(nb) = 0
        cjbasp(nb) = 1.0

        ns = nbasp(nb)

        if (ns .eq. narn1) then
            ixbasp(nb) = 1
        else if (ns.ge.nern1 .and. ns.le.nern2) then
            np = nphasx(ns)
            ne = np - iern1 + 1
            ixbasp(nb) = 1

            do je = 1,jgext(ne)
                nrr1 = jern1(je,ne)
                nrr2 = jern2(je,ne)

                do nss = nrr1,nrr2
                    if (ns .eq. nss) then
                        cjbasp(nb) = cgexj(je,ne)
                    end if
                end do
            end do
        end if
    end do

999 continue
end subroutine gibasp