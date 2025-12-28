subroutine cfxi(conc,fxic,jcsort,narn1,narn2,nstmax,zchsq2)
    !! This subroutine calculates the ionic strength. Note that a sorted
    !! summation is used.
    !! This subroutine is called by:
    !!   EQLIB/ngcadv.f
    !!   EQ3NR/arrset.f
    !!   EQ3NR/scripx.f
    !!   EQ6/scripz.f
    !! Principal input:
    !!   conc   = array of species concentrations
    !!   jcsort = array of species indices, in order of increasing
    !!              concentration, but with sorting restricted to within
    !!              phase ranges
    !!   narn1  = start of the range of aqueous species
    !!   narn2  = end of the range of aqueous species
    !!   zchsq2 = array of one-half the electrical charge squared
    !! Principal output:
    !!   fxic   = the ionic strength (the 2nd-order electrostatic
    !!              moment function I)
    implicit none

    ! Calling sequence variable declarations.
    integer :: nstmax

    integer :: jcsort(nstmax)
    integer :: narn1
    integer :: narn2

    real(kind=8) :: conc(nstmax)
    real(kind=8) :: zchsq2(nstmax)
    real(kind=8) :: fxic

    ! Local variable declarations.
    integer :: ileft
    integer :: ns
    integer :: nss
    integer :: ns0
    integer :: ns1
    integer :: ns2
    integer :: ns3
    integer :: ns4
    integer :: ns5
    integer :: ns6
    integer :: ns7
    integer :: nval

    real(kind=8) :: xx

    ! Logically, this could be represented by:
    !   fxic = 0.
    !   do ns = narn1 + 1,narn2
    !     fxic = fxic + conc(ns)*zchsq2(ns)
    !   enddo
    ! In doing a sorted summation, a complication arises in that
    ! water (ns = narn1) must be included in the range of the
    ! loop, because the jcsort array is organized to contain the
    ! species indices sorted according to increasing concentration
    ! within phase ranges. There is no available sorting among
    ! aqueous solute species only. Because the electrical charge
    ! of water is zero, the calculation can be made looping from
    ! narn1 to narn2. It is not necessary to temporarily set the
    ! concentration of water (conc(narn1)) to zero. Note the use
    ! of a local variable (xx) within the loop. Note also that the
    ! loop is unrolled.
    xx = 0.
    nval = narn2 - narn1 + 1
    ileft = (nval/8)*8 + narn1 - 1

    do nss = narn1,ileft,8
        ns0 = jcsort(nss)
        ns1 = jcsort(nss + 1)
        ns2 = jcsort(nss + 2)
        ns3 = jcsort(nss + 3)
        ns4 = jcsort(nss + 4)
        ns5 = jcsort(nss + 5)
        ns6 = jcsort(nss + 6)
        ns7 = jcsort(nss + 7)
        xx = xx + conc(ns0)*zchsq2(ns0) + conc(ns1)*zchsq2(ns1)          + conc(ns2)*zchsq2(ns2) + conc(ns3)*zchsq2(ns3)          + conc(ns4)*zchsq2(ns4) + conc(ns5)*zchsq2(ns5)          + conc(ns6)*zchsq2(ns6) + conc(ns7)*zchsq2(ns7)
    end do

    do nss = ileft + 1,narn2
        ns = jcsort(nss)
        xx = xx + conc(ns)*zchsq2(ns)
    end do

    fxic = xx
end subroutine cfxi