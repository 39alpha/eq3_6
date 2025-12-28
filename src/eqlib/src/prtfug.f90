subroutine prtfug(jgsort,fugac,fugalg,jsflag,ngrn1,ngt,ngtmax,noutpt,nstmax,uspec)
    !! This subroutine prints a table giving the equilibrium fugacities
    !! of gas species.
    !! This subroutine is called by:
    !!   EQ3NR/scripx.f
    !!   EQ6/scripz.f
    !! Principal input:
    !!   fugac  = gas fugacity array
    !!   fugalg = gas log fugacity array
    !!   jsflag = species status flag array
    !!   ngrn1  = start of gas species range
    !!   ngt    = number of gas species
    !!   uspec  = species name array
    !! Principal output:
    !!   None
    implicit none

    ! Calling sequence variable declarations.
    integer :: ngtmax
    integer :: nstmax

    integer :: noutpt

    integer :: jgsort(ngtmax)
    integer :: jsflag(nstmax)
    integer :: ngrn1
    integer :: ngt

    character(len=48) :: uspec(nstmax)

    real(kind=8) :: fugac(ngtmax)
    real(kind=8) :: fugalg(ngtmax)

    ! Local variable declarations.
    integer :: ncount
    integer :: ng
    integer :: ngg
    integer :: ns

    write(noutpt,1000)
1000 format(/21x,'--- Fugacities ---',//5x,'Gas ',20x,'Log Fugacity',4x,'Fugacity',/)

    ncount = 0

    do ngg= ngt,1,-1
        ng = jgsort(ngg)
        ns = ngrn1 + ng - 1

        if (jsflag(ns) .lt. 2) then
            ncount = ncount + 1
            write (noutpt,1020) uspec(ns),fugalg(ng),fugac(ng)
1020 format(3x,a24,3x,f10.5,3x,1pe12.5)
        end if
    end do

    if (ncount .le. 0) then
        write (noutpt,1030)
    end if

1030 format(3x,'None')

    write (noutpt,1040)
1040 format(1x)
end subroutine prtfug