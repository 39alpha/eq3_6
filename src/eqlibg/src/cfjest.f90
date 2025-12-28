subroutine cfjest(ctb,fjestc,nbaspd,nbt,nbtmax,nstmax,zchcu6)
    !! This subroutine calculates the stoichiometric ionic asymmetry
    !! (fjestc). Note that a negative value of total concentration
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
    !!   zchcu6 =  array of one-sixth the charge cubed
    !! Principal output:
    !!   fjestc = the stoichiometric ionic asymmetry, defined in terms
    !!              of the total molalities of the basis species in the
    !!              'd' set
    implicit none

    ! Calling sequence variable declarations.
    integer :: nbtmax
    integer :: nstmax

    integer :: nbaspd(nbtmax)
    integer :: nbt

    real(kind=8) :: ctb(nbtmax)
    real(kind=8) :: zchcu6(nstmax)
    real(kind=8) :: fjestc

    ! Local variable declarations.
    integer :: nb
    integer :: ns

    ! Note: the relatively small likely value of nbt does not justify
    ! the use of an unrolled loop.
    fjestc = 0.

    do nb = 1,nbt
        ns = nbaspd(nb)
        fjestc = fjestc + zchcu6(ns)*abs(ctb(nb))
    end do
end subroutine cfjest