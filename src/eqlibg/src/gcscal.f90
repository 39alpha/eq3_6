subroutine gcscal(acflgc,delacf,narn1,narn2,nref,nstmax,zchar)
    !! This subroutine rescales the activity coefficients of aqeuous
    !! ionic species to make them consistent with a given pH scale,
    !! which has been used to define the correction parameter "delacf".
    !! This subroutine is called by:
    !!   EQLIBG/gcoeff.f
    !! Input:
    !!   acflgc = calculated activity coefficient array  (unscaled)
    !!   delacf = log gamma of the reference ion (new scale) -
    !!              log gamma of the reference ion (old scale)
    !!   nref   = the species index of the reference ion
    !!   nstmax = the maximum number of species
    !!   zchar  = electrical charge array
    !! Output:
    !!   acflgc = calculated activity coefficient array (scaled)
    implicit none

    ! Calling sequence variable declarations.
    integer :: nstmax

    integer :: narn1
    integer :: narn2
    integer :: nref

    real(kind=8) :: acflgc(nstmax)
    real(kind=8) :: zchar(nstmax)
    real(kind=8) :: delacf

    ! Local variable declarations.
    integer :: ileft
    integer :: ns
    integer :: nval

    real(kind=8) :: zref

    zref = zchar(nref)

    ! Note that the loop can run from narn1 to narn2 without
    ! causing a problem in the case of water because water has no
    ! electrical charge. Note also that the loop is unrolled.
    nval = narn2 - narn1 + 1
    ileft = (nval/8)*8 + narn1 - 1

    do ns = narn1,ileft,8
        acflgc(ns) = acflgc(ns) + (zchar(ns)/zref)*delacf
        acflgc(ns + 1) = acflgc(ns + 1) + (zchar(ns + 1)/zref)*delacf
        acflgc(ns + 2) = acflgc(ns + 2) + (zchar(ns + 2)/zref)*delacf
        acflgc(ns + 3) = acflgc(ns + 3) + (zchar(ns + 3)/zref)*delacf
        acflgc(ns + 4) = acflgc(ns + 4) + (zchar(ns + 4)/zref)*delacf
        acflgc(ns + 5) = acflgc(ns + 5) + (zchar(ns + 5)/zref)*delacf
        acflgc(ns + 6) = acflgc(ns + 6) + (zchar(ns + 6)/zref)*delacf
        acflgc(ns + 7) = acflgc(ns + 7) + (zchar(ns + 7)/zref)*delacf
    end do

    do ns = ileft + 1,narn2
        acflgc(ns) = acflgc(ns) + (zchar(ns)/zref)*delacf
    end do
end subroutine gcscal