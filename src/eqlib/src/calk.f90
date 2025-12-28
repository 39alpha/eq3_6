subroutine calk(alkc,conc,nstmax,ntfx,ntfxmx,ntfxt,tfx)
    !! This subroutine calculates the alkalinity (eq/kg H2O). A sorted
    !! summation is not done here, because relatively few species
    !! contribute to alkalinity. Also, the structure of the titration
    !! factor arrays (tfx, ntfx) is not amenable to an efficient
    !! calculation of that sort.
    !! This subroutine is called by:
    !!   EQLIB/betas.f
    !!   EQLIB/prtalk.f
    !!   EQ3NR/scripx.f
    !!   EQ6/scripz.f
    !! Principal input:
    !!   conc   = array of species concentrations
    !!   ntfx   = array of species indices corresponding to the
    !!              alkalinity factors in the tfx array
    !!   ntfxt  = the number of species contributing to alkalinity
    !!   tfx    = array of alkalinity factors
    !! Principal output:
    !!   alkc   = the calculated alkalinity (eq/kg H20)
    implicit none

    ! Calling sequence variable declarations.
    integer :: nstmax
    integer :: ntfxmx

    integer :: ntfx(ntfxmx)
    integer :: ntfxt

    real(kind=8) :: conc(nstmax)
    real(kind=8) :: tfx(ntfxmx)
    real(kind=8) :: alkc

    ! Local variable declarations.
    integer :: ileft
    integer :: n

    real(kind=8) :: ax

    ! Note that the loop is unrolled.
    ax = 0.
    ileft = (ntfxt/8)*8

    do n = 1,ileft,8
        ax = ax + tfx(n)*conc(ntfx(n))          + tfx(n + 1)*conc(ntfx(n + 1))          + tfx(n + 2)*conc(ntfx(n + 2))          + tfx(n + 3)*conc(ntfx(n + 3))          + tfx(n + 4)*conc(ntfx(n + 4))          + tfx(n + 5)*conc(ntfx(n + 5))          + tfx(n + 6)*conc(ntfx(n + 6))          + tfx(n + 7)*conc(ntfx(n + 7))
    end do

    do n = ileft + 1,ntfxt
        ax = ax + tfx(n)*conc(ntfx(n))
    end do

    alkc = ax
end subroutine calk