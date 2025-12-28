subroutine cfxist(ctb,fxistc,nbaspd,nbt,nbtmax,nstmax,zchsq2)
    !! This subroutine calculates the stoichiometric ionic strength
    !! (fxistc). Note that a negative value of total concentration
    !! (ctb) for H+ is treated as a positive value for OH- (and vice
    !! versa).
    !! This subroutine is called by:
    !!   EQ3NR/scripx.f
    !!   EQ6/scripz.f
    !! Principal input:
    !!   ctb    = array of total concentrations of the basis species in
    !!              the 'd' set
    !!   nbaspd = array containing the indices of the species in the 'd'
    !!              basis set
    !!   nbt    = the number of active basis species
    !!   zchsq2 =  array of one-half the charge squared
    !! Principal output:
    !!   fxistc = the stoichiometric ionic strength, defined in terms
    !!              of the total molalities of the basis species in the
    !!              'd' set
    implicit none

    ! Calling sequence variable declarations.
    integer :: nbtmax
    integer :: nstmax

    integer :: nbaspd(nbtmax)
    integer :: nbt

    real(kind=8) :: ctb(nbtmax)
    real(kind=8) :: zchsq2(nstmax)
    real(kind=8) :: fxistc

    ! Local variable declarations.
    integer :: nb
    integer :: ns

    ! Note: the relatively small likely value of nbt does not justify
    ! the use of an unrolled loop.
    fxistc = 0.

    do nb = 1,nbt
        ns = nbaspd(nb)
        fxistc = fxistc + zchsq2(ns)*abs(ctb(nb))
    end do
end subroutine cfxist