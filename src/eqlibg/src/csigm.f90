subroutine csigm(conc,jcsort,narn1,narn2,nstmax,sigmmc)
    !! This subroutine calculates the sum of the molalities of aqueous
    !! solute species (sigmmc):
    !!   Sigma(i) m(i)
    !! Note that a sorted summation is used.
    !! This subroutine is called by:
    !!   EQLIB/ngcadv.f
    !!   EQLIB/ncmpex.f
    !!   EQ3NR/arrset.f
    !!   EQ3NR/scripx.f
    !!   EQ6/scripz.f
    !! Principal input:
    !!   conc   = array of species concentrations
    !!   jcsort = array of species indices, in order of increasing
    !!              concentration, but with sorting restricted to within
    !!              phase ranges
    !!   narn1  = start of the range of aqueous species; this is
    !!              also the index of the solvent, water
    !!   narn2  = end of the range of aqueous species
    !! Principal output:
    !!   sigmmc = the sum of the molalities of the aqueous solute
    !!              species
    implicit none

    ! Calling sequence variable declarations.
    integer :: nstmax

    integer :: jcsort(nstmax)
    integer :: narn1
    integer :: narn2

    real(kind=8) :: conc(nstmax)
    real(kind=8) :: sigmmc

    ! Local variable declarations.
    integer :: ileft
    integer :: nss
    integer :: nval

    real(kind=8) :: cw
    real(kind=8) :: sx

    ! Logically, this could be represented by:
    !   sigmmc = 0.
    !   do ns = narn1 + 1,narn2
    !     sigmmc = sigmmc + conc(ns)
    !   enddo
    ! In doing a sorted summation, a complication arises in that
    ! water (ns = narn1) must be included in the range of the
    ! loop, because the jcsort array is organized to contain the
    ! species indices sorted according to increasing concentration
    ! within phase ranges. There is no available sorting among
    ! aqueous solute species only. This complication is dealt with
    ! by temporarily setting the concentration of water (conc(narn1))
    ! to zero. Note the use of a local variable (sx) within the loop.
    ! Note also that the loop is unrolled.
    sx = 0.
    cw = conc(narn1)
    conc(narn1) = 0.
    nval = narn2 - narn1 + 1
    ileft = (nval/8)*8 + narn1 - 1

    do nss = narn1,ileft,8
        sx = sx + conc(jcsort(nss))      + conc(jcsort(nss + 1))          + conc(jcsort(nss + 2))  + conc(jcsort(nss + 3))          + conc(jcsort(nss + 4))  + conc(jcsort(nss + 5))          + conc(jcsort(nss + 6))  + conc(jcsort(nss + 7))
    end do

    do nss = ileft + 1,narn2
        sx = sx + conc(jcsort(nss))
    end do

    sigmmc = sx
    conc(narn1) = cw
end subroutine csigm