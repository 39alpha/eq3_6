integer function nbasis(nbasp,nbt,nbtmax,ns)
    !! This subroutine finds the position of the ns-th species in the
    !! array nbasp, which defines the basis set. If the species
    !! is not in the basis set, nbasis is returned as zero.
    !! This subroutine is called by:
    !!   EQLIB/flgset.f
    !!   EQLIB/gcsts.f
    !!   EQLIB/switch.f
    !!   EQLIB/swtchk.f
    !!   EQ3NR/arrsim.f
    !!   EQ3NR/chkinx.f
    !!   EQ3NR/eq3nr.f
    !!   EQ3NR/intbs3.f
    !!   EQ6/combmb.f
    !!   EQ6/intbs6.f
    !! Principal input:
    !!   nbasp  = array of indices of species in the basis set
    !!   nbt    = the number of species in the basis set
    !!   ns     = the index of a species
    !! Principal output:
    !!   nbasis = the basis index of the ns-th species, if this species
    !!              is in the basis set
    implicit none

    ! Calling sequence variable declarations.
    integer :: nbtmax

    integer :: nbasp(nbtmax)
    integer :: nbt
    integer :: ns

    ! Local variable declarations.
    integer :: nb
    integer :: nss

    nbasis = 0

    do nb = 1,nbt
        nss = nbasp(nb)

        if (ns .eq. nss) then
            nbasis = nb
            go to 999
        end if
    end do

999 continue
end function nbasis