subroutine cabar(abar,azero,conc,jcsort,fxi,narn1,narn2,natmax,nstmax,zchsq2)
    !! This subroutine calculates the average hard core diameter of
    !! aqueous ionic species (abar). This average is defined
    !! using "ionic strength weighting"; that is, the weighting factor
    !! for each species is 1/2 the square of its electrical charge.
    !! This thus excludes contributions from electrically neutral
    !! solute species. This "abar" is used in empirical ion size
    !! mixing models of the first-order Debye-Huckel contribution to
    !! aqueous species activity coefficients.
    !! This subroutine is called by:
    !!   EQLIBG/gcoeff.f
    !! Principal input:
    !!   azero  = array of hard core diameters
    !!   conc   = array of species concentrations
    !!   jcsort = array of species indices, in order of increasing
    !!              concentration, but with sorting restricted to within
    !!              phase ranges
    !!   fxi    = the ionic strength (the 2nd-order electrostatic
    !!              moment function I)
    !!   narn1  = start of the range of aqueous species
    !!   narn2  = end of the range of aqueous species
    !!   zchsq2 = array of z(i)**2/2 values
    !! Principal output:
    !!   abar  = average ionic hard core diameter
    implicit none

    ! Calling sequence variable declarations.
    integer :: natmax
    integer :: nstmax

    integer :: jcsort(nstmax)

    integer :: narn1
    integer :: narn2

    real(kind=8) :: azero(natmax)
    real(kind=8) :: conc(nstmax)
    real(kind=8) :: zchsq2(nstmax)
    real(kind=8) :: abar
    real(kind=8) :: fxi

    ! Local variable declarations.
    integer :: ileft
    integer :: nnn
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

    real(kind=8) :: sx

    abar = 0.

    if (fxi .le. 0.) then
        go to 999
    end if

    ! Logically, this could be represented by:
    !   sx = 0.
    !   do ns = narn1 + 1,narn2
    !     na = ns - narn1 + 1
    !     sx = sx + conc(ns)*zchsq2(ns)*azero(na)
    !   enddo
    !   abar = sx/fxi
    ! In doing a sorted summation, a complication arises in that
    ! water (ns = narn1) must be included in the range of the
    ! loop, because the jcsort array is organized to contain the
    ! species indices sorted according to increasing concentration
    ! within phase ranges. There is no available sorting among
    ! aqueous solute species only. Because the electrical charge
    ! of water is zero, the calculation can be made looping from
    ! narn1 to narn2. It is not necessary to temporarily set the
    ! concentration of water (conc(narn1)) to zero. Note the use
    ! of a local variable (sx) within the loop. Note also that the
    ! loop is unrolled.
    sx = 0.
    nnn = narn1 - 1
    nval = narn2 - nnn
    ileft = (nval/8)*8 + nnn

    do nss = narn1,ileft,8
        ns0 = jcsort(nss)
        ns1 = jcsort(nss + 1)
        ns2 = jcsort(nss + 2)
        ns3 = jcsort(nss + 3)
        ns4 = jcsort(nss + 4)
        ns5 = jcsort(nss + 5)
        ns6 = jcsort(nss + 6)
        ns7 = jcsort(nss + 7)
        sx = sx + conc(ns0)*zchsq2(ns0)*azero(ns0 - nnn)          + conc(ns1)*zchsq2(ns1)*azero(ns1 - nnn)          + conc(ns2)*zchsq2(ns2)*azero(ns2 - nnn)          + conc(ns3)*zchsq2(ns3)*azero(ns3 - nnn)          + conc(ns4)*zchsq2(ns4)*azero(ns4 - nnn)          + conc(ns5)*zchsq2(ns5)*azero(ns5 - nnn)          + conc(ns6)*zchsq2(ns6)*azero(ns6 - nnn)          + conc(ns7)*zchsq2(ns7)*azero(ns7 - nnn)
    end do

    do nss = ileft + 1,narn2
        ns = jcsort(nss)
        sx = sx + conc(ns)*zchsq2(ns)*azero(ns - nnn)
    end do

    abar = sx/fxi

999 continue
end subroutine cabar