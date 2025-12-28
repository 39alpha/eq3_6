subroutine lambda(acflgc,afcnst,bpx,ibpxmx,ibpxt,iktmax,ixrn1,ixrn2,jsol,ncmpr,noutpt,np,nptmax,nstmax,nttyo,nxtmax,wfac,xbar,xbarlg,uphase,uspec)
    !! This subroutine computes activity coefficients for the components
    !! of the np-th phase, which is a solid solution.
    !! This subroutine is called by:
    !!   EQLIB/ngcadv.f
    !!   EQ3NR/eq3nr.f
    !!   EQ6/raff.f
    !! Principal input:
    !!   afcnst = affinity constant, 2.303RT/1000
    !!   bpx    = array of site-mixing parameters for solid solution
    !!              activity coefficients
    !!   ibpxt  = array of the numbers of non-zero site-mixing
    !!              paramters for computing activity coefficients in
    !!              solid solutions
    !!   jsol   = array of activity coefficient model flags
    !!   ncmpr  = species range pointer array for phases
    !!   noutpt = unit number of the output file
    !!   nttyo  = unity number of the tty output file
    !!   nx     = solid solution index
    !!   nxtmax = array dimension, maximum number of solid solutions
    !!   wfac   = array of non-ideal mixing parameters calculated from
    !!               the apx array
    !!   xbar   = array of mole fractions
    !!   xbarlg = array of log mole fraction values
    !!   wterm  = array of temperature-dependent coefficients for
    !!              solid solution activity coefficients
    !! Principal output:
    !!   acflgc  = array of activity coefficients
    implicit none

    ! Calling sequence variable declarations.
    integer :: ibpxmx
    integer :: iktmax
    integer :: nptmax
    integer :: nstmax
    integer :: nxtmax

    integer :: noutpt
    integer :: nttyo

    integer :: ibpxt(nxtmax)
    integer :: jsol(nxtmax)
    integer :: ncmpr(2,nptmax)
    integer :: ixrn1
    integer :: ixrn2
    integer :: np

    character(len=48) :: uspec(nstmax)
    character(len=24) :: uphase(nptmax)

    real(kind=8) :: acflgc(nstmax)
    real(kind=8) :: bpx(ibpxmx,nxtmax)
    real(kind=8) :: xbar(nstmax)
    real(kind=8) :: xbarlg(nstmax)
    real(kind=8) :: wfac(iktmax,nxtmax)
    real(kind=8) :: afcnst

    ! Local variable declarations.
    integer :: ix
    integer :: iy
    integer :: j2
    integer :: k
    integer :: nr1
    integer :: nr2
    integer :: ns
    integer :: nt
    integer :: nx

    integer :: ilnobl

    character(len=48) :: ux48
    character(len=48) :: uy48

    real(kind=8) :: stxm1
    real(kind=8) :: xx
    real(kind=8) :: xy

    ! Note: the following statements don't really do anything except
    ! cause the compiler not to complain that ixrn2, uspec, afcnst,
    ! xbar, and wfac are not used.
    ix = ixrn2
    iy = ix
    ixrn2 = iy

    ux48 = uspec(1)
    uy48 = ux48
    uspec(1) = uy48

    xx = afcnst
    xy = xx
    afcnst = xy

    xx = xbar(1)
    xy = xx
    xbar(1) = xy

    xx = wfac(1,1)
    xy = xx
    wfac(1,1) = xy

    nx = np - ixrn1 + 1

    ! XX   Temporarily make all solid solutions ideal (for version 8.0).
    ! XX   Non-ideal solid solutions will be reintroduced in some later
    ! XX   version.
    if (jsol(nx) .ne. 1) then
        jsol(nx) = 1
    end if

    ! XX
    k = jsol(nx)

    nr1 = ncmpr(1,np)
    nr2 = ncmpr(2,np)
    nt = nr2 - nr1 + 1

    if (nt .le. 0) then
        j2 = ilnobl(uphase(np))
        write (noutpt,1000) uphase(np)(1:j2)
        write (nttyo,1000) uphase(np)(1:j2)
1000 format(/' * Error - (EQLIBG/lambda) Programming error trap:',' Have no',/7x,'components present for solid solution ',a,'.')

        stop
    end if

    if (k .eq. 1) then
        if (ibpxt(nx) .eq. 0.) then
            ! Molecular mixing ideal solution.
            do ns = nr1,nr2
                acflgc(ns) = 0.
            end do
        else
            ! Site mixing ideal solution.
            stxm1 = bpx(1,nx) - 1.0

            do ns = nr1,nr2
                acflgc(ns) = stxm1*xbarlg(ns)
            end do
        end if
    else
        j2 = ilnobl(uphase(np))
        write (noutpt,1010) jsol(nx),nx,uphase(np)(1:j2)
        write (nttyo,1010) jsol(nx),nx,uphase(np)(1:j2)
1010 format(/' * Error - (EQLIBG/lambda) Programming error trap:',' Have an',/7x,'undefined value of ',i3,' for jsol(',i3,'),',' the activity'/7x,'coefficient model flag for solid',' solution ',a,'.')

        stop
    end if
end subroutine lambda