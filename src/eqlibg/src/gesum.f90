subroutine gesum(conc,delam,elam,elsump,elsums,elsumw,fxi,narn1,narn2,nazpmx,nstmax,zchar)
    !! This subroutine calculates three second order sums involving the
    !! E-lambda function (elam) and its ionic strength derivatives
    !! (delam). These are used in Pitzer's equations. These
    !! sums are:
    !!   elsumw: SUM(ij) [E-lambda(ij) + I*E-lambda'(ij)]*m(i)*m(j)
    !!   elsums: SUM(ij) E-lambda'(ij)*m(i)*m(j)
    !!   elsump: SUM(ij) E-lambda''(ij)*m(i)*m(j) [Not used]
    !! Here ij refers to any species pair, but non-zero
    !! contributions only arise from cc' and aa' pairs. Note that:
    !!   elsumw = 2*
    !!         [SUM(c'>c) {E-theta(cc') + I*E-theta'(cc')}m(c)m(c')
    !!        + SUM(a'>a) {E-theta(aa') + I*E-theta'(aa')}m(a)m(a')]
    !!   elsums = 2*[SUM(c'>c) E-theta'(cc')m(c)m(c')
    !!             + SUM(a'>a) E-theta'(aa')m(a)m(a')]
    !!   elsump = 2*[SUM(c'>c) E-theta''(cc')m(c)m(c')
    !!             + SUM(a'>a) E-theta''(aa')m(a)m(a')] [Not used]
    !! This subroutine is called by:
    !!   EQLIBG/gcoeff.f
    !! Principal input:
    !!   conc   = array of concentrations
    !!   elam   = array of E-lambda functions
    !!   delam  = array of ionic strength derivatives of E-lambda
    !!              functions
    !!   fxi    = the ionic strength (the 2nd-order electrostatic
    !!              moment function I)
    !!   narn1  = start of species range for aqueous solution
    !!   narn2  = end of species range for aqueous solution
    !!   zchar  = array of charges
    !! Principal output:
    !!   elsumw = sum: SUM(ij) [ E-lambda(ij) + I*E-lambda'(ij) ]*mi*mj
    !!   elsums = sum: SUM(ij) E-lambda'(ij)*mi*mj
    !!   elsump = sum: SUM(ij) E-lambda''(ij)*mi*mj
    implicit none

    ! Calling sequence variable declarations.
    integer :: nazpmx
    integer :: nstmax

    integer :: narn1
    integer :: narn2

    real(kind=8) :: conc(nstmax)
    real(kind=8) :: delam(2,nazpmx,nazpmx)
    real(kind=8) :: elam(nazpmx,nazpmx)
    real(kind=8) :: zchar(nstmax)
    real(kind=8) :: elsump
    real(kind=8) :: elsums
    real(kind=8) :: elsumw
    real(kind=8) :: fxi

    ! Local variable declarations.
    integer :: ijz
    integer :: iz
    integer :: jz
    integer :: nsi
    integer :: nsj

    real(kind=8) :: cij
    real(kind=8) :: el
    real(kind=8) :: elp
    real(kind=8) :: elpp
    real(kind=8) :: zi
    real(kind=8) :: zj

    elsumw = 0.
    elsums = 0.
    elsump = 0.

    do nsi = narn1,narn2
        zi = zchar(nsi)
        iz = nint(abs(zi))

        if (iz .ne. 0) then
            do nsj = narn1,narn2
                if (nsi .ne. nsj) then
                    zj = zchar(nsj)
                    jz = nint(abs(zj))

                    if (jz .ne. 0) then
                        ijz = nint(zchar(nsi)*zchar(nsj))

                        if ((ijz) .gt. 0) then
                            el = elam(iz,jz)
                            elp = delam(1,iz,jz)
                            elpp = delam(2,iz,jz)
                            cij = conc(nsi)*conc(nsj)
                            elsumw = elsumw + (el + fxi*elp)*cij
                            elsums = elsums + elp*cij
                            elsump = elsump + elpp*cij
                        end if
                    end if
                end if
            end do
        end if
    end do
end subroutine gesum