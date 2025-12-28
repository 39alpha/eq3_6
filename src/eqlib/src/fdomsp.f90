subroutine fdomsp(jssort,mosp,nsi,nst,nstmax,weight,wsi)
    !! This subroutine finds the species that dominates a mass balance.
    !! The primary purpose of this subroutine is to support the
    !! pre-Newton-Raphson optimization algorithm. This subroutine is not
    !! intended to support optimizing the basis set.
    !! Note - to save time, it is assumed that the ratio of the
    !! largest value of weight to the smallest non-zero value is
    !! no greater than 100.
    !! This subroutine is called by:
    !!   EQLIB/cfracf.f
    !! Principal input:
    !!   jssort = array of species indices in order of increasing
    !!              mass
    !!   weight = array of stoichiometric weighting factors
    !!   mosp   = array of moles of species
    !!   nst    = number of species
    !! Principal output:
    !!   nsi    = index of the dominant species
    !!   wsi    = weight(nsi)
    implicit none

    ! Calling sequence variable declarations.
    integer :: nstmax

    integer :: jssort(nstmax)
    integer :: nsi
    integer :: nst

    real(kind=8) :: mosp(nstmax)
    real(kind=8) :: weight(nstmax)
    real(kind=8) :: wsi

    ! Local variable declarations.
    integer :: nsc
    integer :: nss
    integer :: ns

    real(kind=8) :: ap0
    real(kind=8) :: ap1
    real(kind=8) :: m0
    real(kind=8) :: m1
    real(kind=8) :: w0
    real(kind=8) :: w1
    real(kind=8) :: p1
    real(kind=8) :: rx

    m0 = 0.
    w0 = 1.
    ap0 = 0.
    nsi = 0

    do nsc = 1,nst
        nss = nst + 1 - nsc
        ns = jssort(nss)
        w1 = weight(ns)
        m1 = mosp(ns)

        if (m1 .le. 0.) then
            go to 20
        end if

        if (w1 .ne. 0.) then
            p1 = w1*m1
            ap1 = abs(p1)

            if (ap1 .gt. ap0) then
                m0 = m1
                w0 = w1
                ap0 = ap1
                nsi = ns
                go to 15
            end if
        end if

        rx = m0/m1

        if (rx .gt. 100.) then
            go to 20
        end if

15 continue
    end do

20 continue

    wsi = w0
end subroutine fdomsp