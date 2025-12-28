subroutine intmtb(mtb,mtbaq,mtbaqi,mtbi,nbasp,nbt,nbti,nbtmax,noutpt,nstmax,nttyo,ubmtbi,uspec)
    !! This subroutine interprets the mass balance totals read from the
    !! input file. It constructs the mtb and mtbaq arrays.
    !! This subroutine is called by:
    !!   EQ6/eq6.f
    !! Principal input:
    !!   mtbi   = total number of moles of basis species in the
    !!              equilibrium system
    !!   mtbaqi = total number of moles of basis species in the
    !!              aqueous phase
    !!   ubmtbi = array of names of the corresponding species
    !! Principal output:
    !!   mtb    = total number of moles of basis species in the
    !!              equilibrium system
    !!   mtbaq  = total number of moles of basis species in the
    !!              aqueous phase
    implicit none

    ! Calling sequence variable declarations.
    integer :: nbtmax
    integer :: nstmax

    integer :: noutpt
    integer :: nttyo

    integer :: nbasp(nbtmax)

    integer :: nbt
    integer :: nbti

    character(len=48) :: ubmtbi(nbtmax)
    character(len=48) :: uspec(nstmax)

    real(kind=8) :: mtb(nbtmax)
    real(kind=8) :: mtbaq(nbtmax)
    real(kind=8) :: mtbaqi(nbtmax)
    real(kind=8) :: mtbi(nbtmax)

    ! Local variable declarations.
    integer :: jlen2
    integer :: nb
    integer :: nbi
    integer :: nerr
    integer :: ns

    character(len=56) :: uspn56

    nerr = 0

    do nbi = 1,nbti
        do nb = 1,nbt
            ns = nbasp(nb)

            if (ubmtbi(nbi)(1:48) .eq. uspec(ns)(1:48)) then
                mtb(nb) = mtbi(nbi)
                mtbaq(nb) = mtbaqi(nbi)
                go to 100
            end if
        end do

        call fmspnm(jlen2,ubmtbi(nbi),uspn56)
        write (noutpt,1000) uspn56(1:jlen2)
        write (nttyo,1000) uspn56(1:jlen2)
1000 format(/' * Error - (EQ6/intmtb) A mass balance is defined',' on the input',/7x,'file for ',a,", but this species isn't",' in the',/7x,"currently active basis set. Either it isn't",' on the current data file',/7x,'or it has been suppressed',' as by an nxmod or iopt(15) option.')

        nerr = nerr + 1
100 continue
    end do

    if (nerr .gt. 0) then
        stop
    end if
end subroutine intmtb