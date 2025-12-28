subroutine cfje(conc,fjec,jcsort,narn1,narn2,nstmax,zchcu6)
    !! This subroutine calculates the ionic asymmetry (the 3-rd order
    !! electrostatic moment function J). This is defined as:
    !!   J = 1/6 SUM(i) m(i)z(i)**3
    !! This is a higher order analogue of the ionic strength, I, which
    !! is the 2nd-order electrostatic moment function. For comparison,
    !! the ionic strength is defined as:
    !!   I = 1/2 SUM(i) m(i)z(i)**2
    !! This J is not to be confused with the functions J0(x), J1(x), and
    !! J2(x), which are also involved in higher-order electrostatic
    !! contributions to activity coefficients (see EQLIBG/ghj0.f).
    !! Note that a sorted summation is used.
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
    !!   zchcu6 = array of one-sixth the electrical charge cubed
    !! Principal output:
    !!   fjec   = the ionic asymmetry (the 3rd-order electrostatic
    !!              moment function J)
    implicit none

    ! Calling sequence variable declarations.
    integer :: nstmax

    integer :: jcsort(nstmax)
    integer :: narn1
    integer :: narn2

    real(kind=8) :: conc(nstmax)
    real(kind=8) :: zchcu6(nstmax)
    real(kind=8) :: fjec

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
    !   fjec = 0.
    !   do ns = narn1 + 1,narn2
    !     fjec = fjec + conc(ns)*zchcu6(ns)
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
        xx = xx + conc(ns0)*zchcu6(ns0) + conc(ns1)*zchcu6(ns1)          + conc(ns2)*zchcu6(ns2) + conc(ns3)*zchcu6(ns3)          + conc(ns4)*zchcu6(ns4) + conc(ns5)*zchcu6(ns5)          + conc(ns6)*zchcu6(ns6) + conc(ns7)*zchcu6(ns7)
    end do

    do nss = ileft + 1,narn2
        ns = jcsort(nss)
        xx = xx + conc(ns)*zchcu6(ns)
    end do

    fjec = xx
end subroutine cfje