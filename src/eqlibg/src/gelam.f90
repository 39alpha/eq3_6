subroutine gelam(aphi,delam,dpelm,elam,fxi,izmax,nazpmx,pelm,qpit75)
    !! This subroutine calculates the E-lambda function (elam(i,j)) and
    !! its first two ionic strength derivatives (delam(1,i,j) and
    !! delam(2,i,j)). Here i, j refers to the charge pair zi, zj or
    !! -zi, -zj, where zi and zj are both positive. The E-lambda
    !! function and its derivatives are used in Pitzer's equations
    !! to represent higher-order electrical interactions. Here izmax
    !! is the max norm of the electrical charges of the aqueous species
    !! and aphi is the Debye-Huckel A(phi) parameter. The array
    !! pelm(i,j) contains a set of "primitive" E-lambda functions
    !! corresponding to those in elam(i,j). The arrays dpelm(1,i,j)
    !! and dpelm(2,i,j) similarly parallel delam(1,i,j) and delam(2,i,j).
    !! This subroutine is called by:
    !!   EQLIBG/gcoeff.f
    !! Principal input:
    !!   aphi   = the Debye-Huckel A(phi) parameter
    !!   izmax  = the max norm of the electrical charge numbers of
    !!              the aqueous species
    !!   fxi    = the ionic strength (the 2nd-order electrostatic
    !!              moment function I)
    !! Principal output:
    !!   elam   = array of values of the E-lambda functions
    !!   delam  = array of values of the ionic strength derivatives
    !!              of the E-lambda functions
    !! Work space:
    !!   pelm   = array of values of primitive E-lambda functions
    !!   dpelm  = array of values of the ionic strength derivatives
    !!              of the primitive E-lambda functions
    implicit none

    ! Calling sequence variable declarations.
    integer :: nazpmx

    integer :: izmax

    logical :: qpit75

    real(kind=8) :: delam(2,nazpmx,nazpmx)
    real(kind=8) :: dpelm(2,nazpmx,nazpmx)
    real(kind=8) :: elam(nazpmx,nazpmx)
    real(kind=8) :: pelm(nazpmx,nazpmx)
    real(kind=8) :: aphi
    real(kind=8) :: fxi

    ! Local variable declarations.
    integer :: i
    integer :: ijz
    integer :: iz
    integer :: jz

    real(kind=8) :: el
    real(kind=8) :: elp
    real(kind=8) :: elpp
    real(kind=8) :: wi
    real(kind=8) :: wj

    ! First calculate "primitive" E-lambdas (pelm) and their
    ! derivatives (dpelm) for the various charge pairs. These depend
    ! only on the product of the electrical charges.
    do jz = 1,izmax
        do iz = jz,izmax
            ijz = iz*jz

            ! Get the primitive E-lambda and its derivatives for
            ! this current charge pair.
            call elmdd(aphi,el,elp,elpp,fxi,ijz,qpit75)
            pelm(iz,jz) = el
            dpelm(1,iz,jz) = elp
            dpelm(2,iz,jz) = elpp

            if (jz .ne. iz) then
                pelm(jz,iz) = el
                dpelm(1,jz,iz) = elp
                dpelm(2,jz,iz) = elpp
            end if
        end do
    end do

    ! Now calculate conventional E-lambdas (elam) and their derivatives
    ! (delam) for the various charge pairs.
    do jz = 1,izmax
        do iz = 1,izmax
            wj = 0.5*iz/jz
            wi = 0.5*jz/iz
            elam(iz,jz) = pelm(iz,jz) - wj*pelm(jz,jz) - wi*pelm(iz,iz)
            delam(1,iz,jz) = dpelm(1,iz,jz) - wj*dpelm(1,jz,jz)    - wi*dpelm(1,iz,iz)
            delam(2,iz,jz) = dpelm(2,iz,jz) - wj*dpelm(2,jz,jz)    - wi*dpelm(2,iz,iz)
        end do
    end do

    ! Set the diagonal elements to zero.
    do i = 1,izmax
        elam(i,i) = 0.
        delam(1,i,i) = 0.
        delam(2,i,i) = 0.
    end do
end subroutine gelam