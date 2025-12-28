subroutine lamgex(acflgc,cgexj,jern1,jern2,jetmax,jgext,net,netmax,nstmax,xbarlg)
    !! This subroutine computes the activity coefficients of generic ion
    !! exchanger species. The exchanger phase is taken to be ideal
    !! in the site-mixing sense. The activity (a) and activity coefficent
    !! (lambda) are related to the mole fraction (x) as follows:
    !!   log a = N log x + N log lambda
    !! where N is the site stoichiometric factor. This ensures that the
    !! activity in the ideal case is equal to x**N, a treatment that
    !! generally works better in hybrid Newton-Raphson iteration than
    !! a treatment based on:
    !!   log a = log x + log lambda
    !! whcih, in the ideal case, requires that:
    !!   log lambda = (N - 1) log x
    !! The treatment adopted here is also preferable in the the activity
    !! coefficient is reserved for treating actual non-ideality, rather
    !! than as a correction used to relate two different definitions of
    !! ideality.
    !! This subroutine is called by:
    !!   EQLIB/ngcadv.f
    !!   EQ3NR/arrset.f
    !!   EQ6/exivar.f
    !!   EQ6/raff.f
    !! Principal input:
    !!   cgexj  = array of site stoichiometry factors for site of generic
    !!              ion exchange phases
    !!   jern1  = array giving the start of the range in the species
    !!              list corresponding to species in the je-th site
    !!              of the ne-th exchanger.
    !!   jern2  = array giving the end of the range in the species
    !!              list corresponding to species in the je-th site
    !!              of the ne-th exchanger.
    !!   jgext  = array giving the number of exchange sites in each of
    !!              the ion exchange phases
    !!   xbarlg = array of mole fractions of species
    !!  Principal output:
    !!   acflgc = array of log activity coefficient values
    implicit none

    ! Calling sequence variable declarations.
    integer :: jetmax
    integer :: netmax
    integer :: nstmax

    integer :: jern1(jetmax,netmax)
    integer :: jern2(jetmax,netmax)
    integer :: jgext(netmax)

    integer :: net

    real(kind=8) :: acflgc(nstmax)
    real(kind=8) :: cgexj(jetmax,netmax)
    real(kind=8) :: xbarlg(nstmax)

    ! Local variable declarations.
    integer :: je
    integer :: ne
    integer :: nr1
    integer :: nr2
    integer :: ns

    ! Compute log lambda for generic exchanger species.
    do ne = 1,net
        do je = 1,jgext(ne)
            nr1 = jern1(je,ne)
            nr2 = jern2(je,ne)

            do ns = nr1,nr2
                acflgc(ns) = 0.
            end do
        end do
    end do
end subroutine lamgex