subroutine ca3bar(azero,a3bar,a3bars,conc,jcsort,narn1,narn2,natmax,nstmax,sigmam)
    !! This subroutine calculates the characteristic average cubed
    !! distance of closest approach for each aqueous solute species
    !! (a3bars) and the average cubed distance of closest approach for
    !! all solute species (a3bar). These quantities appear in the
    !! first order term representing hard core repulsion in models
    !! of aqueous species activity coefficients. In the calculation
    !! of these quantities, the weighting factor corresponding to a
    !! species is its molality.
    !! This subroutine is called by:
    !!   EQLIBG/gcoeff.f
    !! Principal input:
    !!   azero  = array of hard core diameters
    !!   conc   = array of species concentrations
    !!   sigmam = sum of the molalities of the aqueous solute species
    !! Prinicpal output:
    !!   a3bars = characteristic average cubed distance of closest
    !!              approach array
    !!   a3bar  = average cubed distance of closest approach
    implicit none

    ! Calling sequence variable declarations.
    integer :: natmax
    integer :: nstmax

    integer :: jcsort(nstmax)
    integer :: narn1
    integer :: narn2

    real(kind=8) :: azero(natmax)
    real(kind=8) :: a3bars(natmax)
    real(kind=8) :: conc(nstmax)
    real(kind=8) :: a3bar
    real(kind=8) :: sigmam

    ! Local variable declarations.
    integer :: ileft
    integer :: na
    integer :: nnn
    integer :: ns
    integer :: nsj
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

    real(kind=8) :: aij
    real(kind=8) :: aij0
    real(kind=8) :: aij1
    real(kind=8) :: aij2
    real(kind=8) :: aij3
    real(kind=8) :: aij4
    real(kind=8) :: aij5
    real(kind=8) :: aij6
    real(kind=8) :: aij7
    real(kind=8) :: azi
    real(kind=8) :: cw
    real(kind=8) :: sx

    do na = 1,narn2 - narn1 + 1
        a3bars(na) = 0
    end do

    a3bar = 0

    if (sigmam .le. 0.) then
        go to 999
    end if

    ! Get the average value of the cube of the distance of closest
    ! approach for each aqueous solute species. Logically, this could be
    ! represented by:
    !   do ns = narn1 + 1,narn2
    !     na = ns - narn1 + 1
    !     sx = 0.
    !     do nsj = narn1 + 1,narn2
    !       naj = nsj - narn1 + 1
    !       aij = 0.5*(azero(na) + azero(naj))
    !       sx = sx + conc(nsj)*(aij**3)
    !     enddo
    !     a3bars(na) = sx/sigmam
    !   enddo
    ! In doing a sorted summation in the inner loop, a complication
    ! arises in that water (ns = narn1) must be included in the
    ! range of the loop, because the jcsort array is organized to
    ! contain the species indices sorted according to increasing
    ! concentration within phase ranges. There is no available
    ! sorting among aqueous solute species only. This complication
    ! is dealt with by temporarily setting the concentration of water
    ! (conc(narn1)) to zero. Note the use of a local variable (sx)
    ! within the loop. Note also that the loop is unrolled.
    cw = conc(narn1)
    conc(narn1) = 0.
    nnn = narn1 - 1
    nval = narn2 - nnn
    ileft = (nval/8)*8 + nnn

    do ns = narn1 + 1,narn2
        na = ns - narn1 + 1
        azi = azero(na)
        sx = 0.

        do nss = narn1,ileft,8
            ns0 = jcsort(nss)
            ns1 = jcsort(nss + 1)
            ns2 = jcsort(nss + 2)
            ns3 = jcsort(nss + 3)
            ns4 = jcsort(nss + 4)
            ns5 = jcsort(nss + 5)
            ns6 = jcsort(nss + 6)
            ns7 = jcsort(nss + 7)
            aij0 = 0.5*(azi + azero(ns0 - nnn))
            aij1 = 0.5*(azi + azero(ns1 - nnn))
            aij2 = 0.5*(azi + azero(ns2 - nnn))
            aij3 = 0.5*(azi + azero(ns3 - nnn))
            aij4 = 0.5*(azi + azero(ns4 - nnn))
            aij5 = 0.5*(azi + azero(ns5 - nnn))
            aij6 = 0.5*(azi + azero(ns6 - nnn))
            aij7 = 0.5*(azi + azero(ns7 - nnn))
            sx = sx + conc(ns0)*(aij0**3) + conc(ns1)*(aij1**3)            + conc(ns2)*(aij2**3) + conc(ns3)*(aij3**3)            + conc(ns4)*(aij4**3) + conc(ns5)*(aij5**3)            + conc(ns6)*(aij6**3) + conc(ns7)*(aij7**3)
        end do

        do nss = ileft + 1,narn2
            nsj = jcsort(nss)
            aij = 0.5*(azi + azero(nsj - nnn))
            sx = sx + conc(ns)*(aij**3)
        end do

        a3bars(na) = sx/sigmam
    end do

    conc(narn1) = cw

    ! Get the average value of the cube of the distance of closest
    ! approach for all aqeuous solute species. Logically, this could be
    ! represented by:
    !   sx = 0.
    !   do ns = narn1 + 1,narn2
    !     na = ns - narn1 + 1
    !     sx = sx + conc(ns)*a3bars(na)
    !   enddo
    !   a3bar = sx/sigmam
    ! Because a3bars(narn1) is defined to be zero, a sorted summation
    ! can be safely made by looping from narn1 to narn2. Note the use
    ! of a local variable (sx) within the loop. Note also that the loop
    ! is unrolled.
    do nss = narn1,ileft,8
        ns0 = jcsort(nss)
        ns1 = jcsort(nss + 1)
        ns2 = jcsort(nss + 2)
        ns3 = jcsort(nss + 3)
        ns4 = jcsort(nss + 4)
        ns5 = jcsort(nss + 5)
        ns6 = jcsort(nss + 6)
        ns7 = jcsort(nss + 7)
        sx = sx + conc(ns0)*a3bars(ns0 - nnn)          + conc(ns1)*a3bars(ns1 - nnn)          + conc(ns2)*a3bars(ns2 - nnn)          + conc(ns3)*a3bars(ns3 - nnn)          + conc(ns4)*a3bars(ns4 - nnn)          + conc(ns5)*a3bars(ns5 - nnn)          + conc(ns6)*a3bars(ns6 - nnn)          + conc(ns7)*a3bars(ns7 - nnn)
    end do

    do nss = ileft + 1,narn2
        ns = jcsort(nss)
        sx = sx + conc(ns)*a3bars(ns - nnn)
    end do

    a3bar = sx/sigmam

999 continue
end subroutine ca3bar