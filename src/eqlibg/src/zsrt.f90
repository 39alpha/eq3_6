subroutine zsrt(izmax,narn1,narn2,nstmax,zchar)
    !! This subroutine finds the largest absolute value of the electrical
    !! charge of any aqueous species (izmax). It is used as a limit in
    !! calculating higher-order electrostatic terms in Pitzer's
    !! equations.
    !! This subroutine is called by:
    !!   EQ3NR/eq3nr.f
    !!   EQ6/eq6.f
    !! Principal input:
    !!   zchar  = array of electrical charge numbers
    !!   narn1  = start of the range of aqueous species
    !!   narn2  = end of the range of aqueous species
    !! Principal output:
    !!   izmax  = max norm of the electrical charges of the
    !!            aqueous species
    implicit none

    ! Calling sequence variable declarations.
    integer :: nstmax

    integer :: izmax
    integer :: narn1
    integer :: narn2

    real(kind=8) :: zchar(nstmax)

    ! Local variable declarations.
    integer :: ileft
    integer :: jzmax
    integer :: jzmin
    integer :: ns
    integer :: nval

    real(kind=8) :: zmax
    real(kind=8) :: zmin
    real(kind=8) :: zx
    real(kind=8) :: zx0
    real(kind=8) :: zx1
    real(kind=8) :: zx2
    real(kind=8) :: zx3
    real(kind=8) :: zx4
    real(kind=8) :: zx5
    real(kind=8) :: zx6
    real(kind=8) :: zx7

    ! Find the extreme values of electrical charge among the aqueous
    ! species. Note that the loop is unrolled.
    zmax = 0.
    zmin = 0.
    nval = narn2 - narn1 + 1
    ileft = (nval/8)*8 + narn1 - 1

    do ns = narn1,ileft,8
        zx0 = zchar(ns)
        zx1 = zchar(ns + 1)
        zx2 = zchar(ns + 2)
        zx3 = zchar(ns + 3)
        zx4 = zchar(ns + 4)
        zx5 = zchar(ns + 5)
        zx6 = zchar(ns + 6)
        zx7 = zchar(ns + 7)
        zmin = min(zx0,zx1,zx2,zx3,zx4,zx5,zx6,zx7,zmin)
        zmax = max(zx0,zx1,zx2,zx3,zx4,zx5,zx6,zx7,zmax)
    end do

    do ns = ileft + 1,narn2
        zx = zchar(ns)
        zmin = min(zx,zmin)
        zmax = max(zx,zmax)
    end do

    ! Find the largest absolute value of electrical charge.
    jzmin = nint(zmin)
    jzmax = nint(zmax)
    izmax = max(-jzmin,jzmax)
end subroutine zsrt