subroutine gndacb(cdac,cdacb,cdrs,eps100,imech,imchmx,jflag,nbasp,nbt,nbtmax,ndac,ndacb,ndact,ndactb,ndctmx,ndrs,ndrsmx,ndrsr,nrct,nrctmx,nstmax)
    !! This subroutine computes the ndacb, ndactb, and cdacb arrays,
    !! which are used to support the higher-order stiff ODE corrector.
    !! Basically, these are the respective analogs of the ndac, ndact,
    !! and cdac arrays. The latter arrays partially describe the
    !! kinetic activity products of terms in rate expressions. The
    !! ndac array contains the species indices of the species whose
    !! activities appear in an activity product, the ndact array
    !! contains the total numbers of such species for the kinetic
    !! activity products, and the cdac array contains the respective
    !! exponents for the species whose indices are given in the
    !! ndac array. The species whose indices appear in the ndac
    !! array need not be basis species. The ndacb, ndactb, and
    !! cdacb arrays are transforms in which the species in ndacb
    !! are all members of the active basis set. These arrays must
    !! be adjusted if this set is changed by basis switching.
    !! This subroutine is called by:
    !!   EQ6/path.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: imchmx
    integer :: nbtmax
    integer :: ndctmx
    integer :: ndrsmx
    integer :: nrctmx
    integer :: nstmax

    integer :: nbt
    integer :: nrct

    integer :: imech(2,nrctmx)
    integer :: jflag(nstmax)
    integer :: nbasp(nbtmax)
    integer :: ndac(ndctmx,imchmx,2,nrctmx)
    integer :: ndacb(nbt,imchmx,2,nrct)
    integer :: ndact(imchmx,2,nrctmx)
    integer :: ndactb(imchmx,2,nrct)
    integer :: ndrs(ndrsmx)
    integer :: ndrsr(2,nstmax)

    real(kind=8) :: cdac(ndctmx,imchmx,2,nrctmx)
    real(kind=8) :: cdacb(nbt,imchmx,2,nrct)
    real(kind=8) :: cdrs(ndrsmx)
    real(kind=8) :: eps100

    ! Local variable declarations.
    integer :: i
    integer :: j
    integer :: k
    integer :: n
    integer :: nb
    integer :: nn
    integer :: nnn
    integer :: ns
    integer :: nse
    integer :: nr1
    integer :: nr2
    integer :: nsi

    real(kind=8) :: cxe

    do k = 1,nrct
        do j = 1,2
            do i = 1,imech(j,k)
                ! Process the current kinetic activity product.
                nn = 0

                do nb = 1,nbt
                    nse = nbasp(nb)

                    if (jflag(nse) .eq. 0) then
                        cxe = 0.

                        do n = 1,ndact(i,j,k)
                            ns = ndac(n,i,j,k)

                            if (ns .eq. nse) then
                                ! The current species in the kinetic activity product
                                ! is the basis species presently being examined.
                                cxe = cxe + cdac(n,i,j,k)
                            else if (jflag(ns) .eq. 30) then
                                ! The current species in the activity product
                                ! is a dependent species. Look for a dependence
                                ! on the basis species presently being examined.
                                nr1 = ndrsr(1,ns)
                                nr2 = ndrsr(2,ns)

                                do nnn = nr1 + 1,nr2
                                    nsi = ndrs(nnn)

                                    if (nsi .eq. nse) then
                                        cxe = cxe - cdac(n,i,j,k)*(cdrs(nnn)/cdrs(nr1))
                                    end if
                                end do
                            end if

                            if (abs(cxe) .ge. eps100) then
                                ! The basis species presently being examined
                                ! belongs in the current transformed kinetic
                                ! activity product.
                                nn = nn + 1
                                ndacb(nn,i,j,k) = nse
                                cdacb(nn,i,j,k) = cxe
                            end if
                        end do
                    end if
                end do

                ndactb(i,j,k) = nn

                do nn = ndactb(i,j,k) + 1,nbt
                    ndacb(nn,i,j,k) = 0
                    cdacb(nn,i,j,k) = 0.
                end do
            end do
        end do
    end do
end subroutine gndacb