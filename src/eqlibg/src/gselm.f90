subroutine gselm(conc,delam,dselm,elam,izmax,narn1,narn2,nazmmx,nazpmx,nstmax,selm,zchar)
    !! This subroutine computes the following first order sums used
    !! in Pitzer's equations:
    !!   selm(i):  SUM(j) E-lambda(ij)*m(j)
    !!   dselm(1,i): SUM(j) E-lambda'(ij)*m(j)
    !! This subroutine is called by:
    !!   EQLIBG/gcoeff.f
    !! Principal input:
    !!   conc   = array of concentrations
    !!   elam   = array of E-lambda functions
    !!   delam  = array of ionic strength derivatives of E-lambda
    !!              functions
    !!   narn1  = start of species range for aqueous solution
    !!   narn2  = end of species range for aqueous solution
    !!   zchar  = array of charges
    !! Principal output:
    !!   selm   = array of sums: 2 SUM(j) E-lambda(ij)*m(j)
    !!   dselm  = array of the corresponding derivatives with respect
    !!              to ionic strength
    implicit none

    ! Calling sequence variable declarations.
    integer :: nazmmx
    integer :: nazpmx
    integer :: nstmax

    integer :: izmax
    integer :: narn1
    integer :: narn2

    real(kind=8) :: delam(2,nazpmx,nazpmx)
    real(kind=8) :: dselm(2,nazmmx:nazpmx)
    real(kind=8) :: conc(nstmax)
    real(kind=8) :: elam(nazpmx,nazpmx)
    real(kind=8) :: selm(nazmmx:nazpmx)
    real(kind=8) :: zchar(nstmax)

    ! Local variable declarations.
    integer :: iajz
    integer :: iz
    integer :: j
    integer :: jz

    real(kind=8) :: sum1
    real(kind=8) :: sum1p
    real(kind=8) :: sum2
    real(kind=8) :: sum2p
    real(kind=8) :: zj

    ! Note: nazmmx = -nazpmx.
    do iz = nazmmx,nazpmx
        selm(iz) = 0.
        dselm(1,iz) = 0.
    end do

    ! Looping over iz from 1 to izmax, get the quantities for charges
    ! iz and -iz simultaneously.
    do iz = 1,izmax
        sum1 = 0.
        sum1p = 0.
        sum2 = 0.
        sum2p = 0.

        do j = narn1 + 1,narn2
            zj = zchar(j)
            jz = nint(zj)
            iajz = abs(jz)

            if (jz .gt. 0) then
                sum1 = sum1 + elam(iz,iajz)*conc(j)
                sum1p = sum1p + delam(1,iz,iajz)*conc(j)
            else if (jz .lt. 0) then
                sum2 = sum2 + elam(iz,iajz)*conc(j)
                sum2p = sum2p + delam(1,iz,iajz)*conc(j)
            end if
        end do

        selm(iz) = sum1
        dselm(1,iz) = sum1p
        selm(-iz) = sum2
        dselm(1,-iz) = sum2p
    end do
end subroutine gselm