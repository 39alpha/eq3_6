subroutine wterm(apx,iapxmx,iktmax,ixrn1,ixrn2,jsol,ncmpr,noutpt,nptmax,nstmax,nttyo,nxt,nxtmax,press,tempk,uphase,uspec,wfac)
    !! This subroutine computes the wfac(i,nx) array, which contains
    !! the coefficients for the excess free energy function of
    !! solid solutions. If non-zero coefficients are lacking, the
    !! solid solution reduces to an ideal solution and the jsol value
    !! is changed to zero to reflect this.
    !! This subroutine is called by:
    !!   EQLIB/evdata.f
    !! Principal input:
    !!   apx    = array of non-ideal mixing parameters
    !!   jsol   = array of non-ideal mixing law flags
    !!   press  = pressure, bars
    !!   tempk  = temperature, K
    !!   uphase = array of phase names
    !! Principal input:
    !!   wfac   = array of non-ideal mixing parameters calculated from
    !!              the apx array
    implicit none

    ! Calling sequence variable declarations.
    integer :: iapxmx
    integer :: iktmax
    integer :: nptmax
    integer :: nstmax
    integer :: nxtmax

    integer :: jsol(nxtmax)
    integer :: ncmpr(2,nptmax)
    integer :: ixrn1
    integer :: ixrn2
    integer :: noutpt
    integer :: nttyo
    integer :: nxt

    character(len=48) :: uspec(nstmax)
    character(len=24) :: uphase(nptmax)

    real(kind=8) :: apx(iapxmx,nxtmax)
    real(kind=8) :: wfac(iktmax,nxtmax)
    real(kind=8) :: press
    real(kind=8) :: tempk

    ! Local variable declarations.
    integer :: i
    integer :: ix
    integer :: iy
    integer :: j2
    integer :: k
    integer :: np
    integer :: nx

    integer :: ilnobl

    character(len=48) :: ux48
    character(len=48) :: uy48

    real(kind=8) :: sum

    ! Note: the following statements don't really do anything except
    ! cause the compiler not to complain that uspec, ncmpr, and ixrn2
    ! are not used.
    ux48 = uspec(1)
    uy48 = ux48
    uspec(1) = uy48

    ix = ncmpr(1,1)
    iy = ix
    ncmpr(1,1) = iy

    ix = ixrn2
    iy = ix
    ixrn2 = iy

    np = ixrn1 - 1

    do nx = 1,nxt
        np = np + 1
        sum = 0.

        do i = 1,iapxmx
            sum = sum + apx(i,nx)*apx(i,nx)
        end do

        if (jsol(nx).eq.1 .and. sum.gt. 0.) then
            j2 = ilnobl(uphase(np))
            write (noutpt,1735) uphase(np)(1:j2)
            write (nttyo,1735) uphase(np)(1:j2)
1735 format(/' * Note - (EQLIBG/wterm) The phase ',a,' is',/7x,'listed as ideal. This is inconsistent with the',/7x,'presence of non-zero data in the apx array. These',/7x,'data will be ignored.')

            do i = 1,iapxmx
                apx(i,nx) = 0.
            end do
        end if

        if (jsol(nx).ne.1 .and. sum.eq.0.) then
            j2 = ilnobl(uphase(np))
            write (noutpt,1736) uphase(np)(1:j2)
            write (nttyo,1736) uphase(np)(1:j2)
1736 format(/' * Note - (EQLIBG/wterm) The phase ',a,' is',/7x,'listed as non-ideal, but no non-zero data are present'    /7x,'in the apx array. This phase will be treated as ideal.')

            jsol(nx) = 1
        end if

        ! Initialize wfac values to zero
        do i = 1,iktmax
            wfac(i,nx) = 0.
        end do
    end do

    np = ixrn1 - 1

    do nx = 1,nxt
        np = np + 1
        k = jsol(nx)

        if (k .eq. 1) then
            ! Ideal solution.
            continue
        else if (k.eq.2) then
            ! binary solution, third-order maclaurin expansion
            ! original pathi solid solution model.
            wfac(1,nx) = apx(1,nx)
            wfac(2,nx) = apx(2,nx)
            wfac(3,nx) = apx(3,nx)
            wfac(1,nx) = -wfac(2,nx)/2. - wfac(3,nx)/6.
        else if (k.eq.3) then
            ! binary solution, parabolic maclaurin expansion
            wfac(1,nx) = apx(1,nx)
        else if (k.eq.4) then
            ! binary solution, cubic maclaurin (p,t dependent)
            wfac(1,nx) = apx(1,nx) + apx(2,nx)*tempk + apx(3,nx)*press
            wfac(2,nx) = apx(4,nx) + apx(5,nx)*tempk + apx(6,nx)*press
        else if (k.eq.5) then
            ! binary solution, guggenheim polynomial (t dependent)
            wfac(1,nx) = apx(1,nx) + apx(2,nx)*tempk    + apx(3,nx)*tempk*tempk
            wfac(2,nx) = apx(4,nx) + apx(5,nx)*tempk    + apx(6,nx)*tempk*tempk
            wfac(3,nx) = apx(7,nx) + apx(8,nx)*tempk    + apx(9,nx)*tempk*tempk
        else if (k.eq.6) then
            ! ternary regular solution (see prigogine and defay, p. 257)
            wfac(1,nx) = apx(1,nx)
            wfac(2,nx) = apx(2,nx)
            wfac(3,nx) = apx(3,nx)
        else if (k.eq.7) then
            ! newton et al plagioclase model (gca vol 44 p. 933, 1980)
            ! 1 - albite; 2 - anorthite
            wfac(1,nx) = apx(1,nx)
            wfac(2,nx) = apx(2,nx)
        else
            j2 = ilnobl(uphase(np))
            write (noutpt,100) jsol(nx),uphase(np)(1:j2)
            write (nttyo,100) jsol(nx),uphase(np)(1:j2)
100 format(/' * Error - (EQLIBG/wterm) Have an Illegal jsol',/7x,'value of ',i2,' for solid solution ',a,'.')
        end if
    end do
end subroutine wterm
