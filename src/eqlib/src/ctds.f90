subroutine ctds(jcsort,mosp,mwtsp,narn1,narn2,nstmax,wotds)
    !! This subroutine calculates the total dissolved solute mass
    !! (wotds, g). Note that a sorted summation is used.
    !! This subroutine is called by:
    !!   EQLIB/gwdenp.f
    !! Principal input:
    !!   mosp   = array of numbers of moles of species
    !!   jcsort = array of species indices, in order of increasing
    !!              concentration/number of moles, but with sorting
    !!              restricted to within phase ranges
    !!   mwtsp  = array of molecular weights of species
    !!   narn1  = start of the range of aqueous species
    !!   narn2  = end of the range of aqueous species
    !! Principal output:
    !!   wotds  = the total dissolved solute mass (g)
    implicit none

    ! Calling sequence variable declarations.
    integer :: nstmax

    integer :: jcsort(nstmax)
    integer :: narn1
    integer :: narn2

    real(kind=8) :: mosp(nstmax)
    real(kind=8) :: mwtsp(nstmax)
    real(kind=8) :: wotds

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

    real(kind=8) :: mxw
    real(kind=8) :: wx

    ! Logically, this could be represented by:
    !   wotds = 0.
    !   do ns = narn1 + 1,narn2
    !     wotds = wotds + mosp(ns)*mwtsp(ns)
    !   enddo
    ! In doing a sorted summation, a complication arises in that
    ! water (ns = narn1) must be included in the range of the
    ! loop, because the jcsort array is organized to contain the
    ! species indices sorted according to increasing concentration
    ! within phase ranges. There is no available sorting among
    ! aqueous solute species only. Thus, the concentration of water
    ! (conc(nanrn1)) is temporarily set to zero. Note the use
    ! of a local variable (wx) within the loop. Note also that the
    ! loop is unrolled.
    wx = 0.
    mxw = mosp(narn1)
    mosp(narn1) = 0.
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
        wx = wx + mosp(ns0)*mwtsp(ns0) + mosp(ns1)*mwtsp(ns1)          + mosp(ns2)*mwtsp(ns2) + mosp(ns3)*mwtsp(ns3)          + mosp(ns4)*mwtsp(ns4) + mosp(ns5)*mwtsp(ns5)          + mosp(ns6)*mwtsp(ns6) + mosp(ns7)*mwtsp(ns7)
    end do

    do nss = ileft + 1,narn2
        ns = jcsort(nss)
        wx = wx + mosp(ns)*mwtsp(ns)
    end do

    wotds = wx
    mosp(narn1) = mxw
end subroutine ctds