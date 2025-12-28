real(kind=8) function coefst(csts,nsts,nstsmx,nstsr,nb,ns,nstmax)
    !! This subroutine finds the stoichiometric coefficient of the nb-th
    !! basis species for the ns-th species. Note that 'nb' denotes
    !! position in the list of basis species, not the general list of
    !! species. The corresponding index on the latter list would be
    !! given by nbasp(nb).
    !! This subroutine is called by:
    !!   EQLIB/prtpct.f
    !!   EQ3NR/betas.f
    !!   EQ3NR/matrix.f
    !!   EQ3NR/scripx.f
    !!   EQ6/betaz.f
    !!   EQ6/matrxz.f
    !!   EQ6/pabssw.f
    !!   EQ6/setffg.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: nstsmx
    integer :: nstmax

    integer :: nsts(nstsmx)
    integer :: nstsr(2,nstmax)
    integer :: nb
    integer :: ns

    real(kind=8) :: csts(nstsmx)

    ! Local variable declarations.
    integer :: n
    integer :: n1
    integer :: n2

    coefst = 0.
    n1 = nstsr(1,ns)
    n2 = nstsr(2,ns)

    do n = n1,n2
        if (nb .eq. nsts(n)) then
            coefst = csts(n)
            go to 999
        end if
    end do

999 continue
end function coefst