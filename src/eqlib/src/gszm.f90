subroutine gszm(conc,jcsort,narn1,narn2,nstmax,sigza,sigzc,sigzi,sigzm,zchar)
    !! This subroutine calculates the sums of equivalent concentrations
    !! and the charge imbalance. Note that a sorted summation is
    !! used.
    !! This subroutine is called by:
    !!   EQLIB/betas.f
    !!   EQ3NR/arrset.f
    !!   EQ3NR/scripx.f
    !!   EQ6/scripz.f
    !! Principal input:
    !!   conc   = array of species concentrations
    !!   jcsort = array of species indices, in order of increasing
    !!              concentration, but with sorting restricted to within
    !!              phase ranges
    !!   narn1  = start of the range of aqeuous species; also
    !!              the species index of solvent water
    !!   narn2  = end of the range of aqeuous species
    !!   zchar  = array of electrical charge numbers
    !! Principal output:
    !!   sigza  = the sum of equivalent concentrations of aqueous
    !!              anions, SUM(i) abs(z(i))m(i), for z(i) < 0
    !!   sigzc  = the sum of equivalent concentrations of aqueous
    !!              cations, SUM(i) z(i)m(i), for z(i) > 0
    !!   sigzi  = the calculated charge imbalance, sigzc - sigza
    !!   sigzm  = the sum of equivalent concentrations of aqueous
    !!              ions, sigzc + sigza
    implicit none

    ! Calling sequence variable declarations.
    integer :: nstmax

    integer :: jcsort(nstmax)
    integer :: narn1
    integer :: narn2

    real(kind=8) :: conc(nstmax)
    real(kind=8) :: zchar(nstmax)
    real(kind=8) :: sigza
    real(kind=8) :: sigzc
    real(kind=8) :: sigzi
    real(kind=8) :: sigzm

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

    real(kind=8) :: ec
    real(kind=8) :: sxa
    real(kind=8) :: sxc
    real(kind=8) :: sxi
    real(kind=8) :: sxm

    ! Logically, this could be represented by:
    !   sigzc = 0.
    !   sigza = 0.
    !   do ns = narn1,narn2
    !     if (zchar(ns) .gt. 0.) then
    !       sigzc = sigzc + conc(ns)*zchar(ns)
    !     elseif (zchar(ns) .lt. 0.) then
    !       sigza = sigza + abs(conc(ns)*zchar(ns))
    !     endif
    !   enddo
    !   sigzi = sigzc - sigza
    !   sigzm = sigzc + sigza
    sigzc = 0.
    sigza = 0.
    sigzm = 0.
    sigzi = 0.

    ! Note the use of local variables (sxc, sxa, sxi, and sxm) within
    ! the loop. Note also that the loop is unrolled.
    sxc = 0.
    sxa = 0.
    sxi = 0.
    sxm = 0.
    nnn = narn1 - 1
    nval = narn2 - nnn
    ileft = (nval/8)*8 + nnn

    do nss = narn1,ileft,8
        ns0 = jcsort(nss)

        if (zchar(ns0) .gt. 0.) then
            ec = zchar(ns0)*conc(ns0)
            sxc = sxc + ec
            sxi = sxi + ec
            sxm = sxm + ec
        else if (zchar(ns0) .lt. 0.) then
            ec = - zchar(ns0)*conc(ns0)
            sxa = sxa + ec
            sxi = sxi - ec
            sxm = sxm + ec
        end if

        ns1 = jcsort(nss + 1)

        if (zchar(ns1) .gt. 0.) then
            ec = zchar(ns1)*conc(ns1)
            sxc = sxc + ec
            sxi = sxi + ec
            sxm = sxm + ec
        else if (zchar(ns1) .lt. 0.) then
            ec = - zchar(ns1)*conc(ns1)
            sxa = sxa + ec
            sxi = sxi - ec
            sxm = sxm + ec
        end if

        ns2 = jcsort(nss + 2)

        if (zchar(ns2) .gt. 0.) then
            ec = zchar(ns2)*conc(ns2)
            sxc = sxc + ec
            sxi = sxi + ec
            sxm = sxm + ec
        else if (zchar(ns2) .lt. 0.) then
            ec = - zchar(ns2)*conc(ns2)
            sxa = sxa + ec
            sxi = sxi - ec
            sxm = sxm + ec
        end if

        ns3 = jcsort(nss + 3)

        if (zchar(ns3) .gt. 0.) then
            ec = zchar(ns3)*conc(ns3)
            sxc = sxc + ec
            sxi = sxi + ec
            sxm = sxm + ec
        else if (zchar(ns3) .lt. 0.) then
            ec = - zchar(ns3)*conc(ns3)
            sxa = sxa + ec
            sxi = sxi - ec
            sxm = sxm + ec
        end if

        ns4 = jcsort(nss + 4)

        if (zchar(ns4) .gt. 0.) then
            ec = zchar(ns4)*conc(ns4)
            sxc = sxc + ec
            sxi = sxi + ec
            sxm = sxm + ec
        else if (zchar(ns4) .lt. 0.) then
            ec = - zchar(ns4)*conc(ns4)
            sxa = sxa + ec
            sxi = sxi - ec
            sxm = sxm + ec
        end if

        ns5 = jcsort(nss + 5)

        if (zchar(ns5) .gt. 0.) then
            ec = zchar(ns5)*conc(ns5)
            sxc = sxc + ec
            sxi = sxi + ec
            sxm = sxm + ec
        else if (zchar(ns5) .lt. 0.) then
            ec = - zchar(ns5)*conc(ns5)
            sxa = sxa + ec
            sxi = sxi - ec
            sxm = sxm + ec
        end if

        ns6 = jcsort(nss + 6)

        if (zchar(ns6) .gt. 0.) then
            ec = zchar(ns6)*conc(ns6)
            sxc = sxc + ec
            sxi = sxi + ec
            sxm = sxm + ec
        else if (zchar(ns6) .lt. 0.) then
            ec = - zchar(ns6)*conc(ns6)
            sxa = sxa + ec
            sxi = sxi - ec
            sxm = sxm + ec
        end if

        ns7 = jcsort(nss + 7)

        if (zchar(ns7) .gt. 0.) then
            ec = zchar(ns7)*conc(ns7)
            sxc = sxc + ec
            sxi = sxi + ec
            sxm = sxm + ec
        else if (zchar(ns7) .lt. 0.) then
            ec = - zchar(ns7)*conc(ns7)
            sxa = sxa + ec
            sxi = sxi - ec
            sxm = sxm + ec
        end if
    end do

    do nss = ileft + 1,narn2
        ns = jcsort(nss)

        if (zchar(ns) .gt. 0.) then
            ec = zchar(ns)*conc(ns)
            sxc = sxc + ec
            sxi = sxi + ec
            sxm = sxm + ec
        else if (zchar(ns) .lt. 0.) then
            ec = - zchar(ns)*conc(ns)
            sxa = sxa + ec
            sxi = sxi - ec
            sxm = sxm + ec
        end if
    end do

    sigzc = sxc
    sigza = sxa
    sigzi = sxi
    sigzm = sxm
end subroutine gszm