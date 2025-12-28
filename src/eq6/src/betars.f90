subroutine betars(alphar,betar,btrfnc,btrmax,btrmxo,ibtrmx,nrct,nrct1,nrctmx,rirec1,rirecp,rrelr1,rrelrp)
    !! This subroutine calculates residual functions for the ODE
    !! integrator.
    !! This subroutine is called by:
    !!   EQ6/path.f
    !! Principal input:
    !! Principal output:
    !!   alphar = array of alpha[r] residuals
    !!   betar  = array of beta[r] residuals
    implicit none

    ! Calling sequence variable declarations.
    integer :: nrctmx

    integer :: ibtrmx
    integer :: nrct
    integer :: nrct1

    real(kind=8) :: alphar(nrct1)
    real(kind=8) :: betar(nrct1)
    real(kind=8) :: rrelr1(nrctmx)
    real(kind=8) :: rrelrp(nrctmx)

    real(kind=8) :: btrfnc
    real(kind=8) :: btrmax
    real(kind=8) :: btrmxo
    real(kind=8) :: rirec1
    real(kind=8) :: rirecp

    ! Local variable declarations.
    integer :: nrc

    integer :: iarmxn

    ! Calculate absolute (alphar) and normalized (betar) residuals.
    do nrc = 1,nrct
        alphar(nrc) = rrelr1(nrc) - rrelrp(nrc)
        betar(nrc) = 0.

        if (rrelrp(nrc) .ne. 0.) then
            betar(nrc) = alphar(nrc)/rrelrp(nrc)
        end if
    end do

    alphar(nrct1) = rirec1 - rirecp

    if (rirecp .ne. 0.) then
        betar(nrct1) = alphar(nrct1)/rirecp
    end if

    ! Find the max norm of the normalized residual vector (betar).
    ibtrmx = iarmxn(betar,nrct1)
    btrmax = 0.

    if (ibtrmx .gt. 0) then
        btrmax = abs(betar(ibtrmx))
    end if

    ! Calculate the betar improvement function (btrfnc).
    btrfnc = 0.

    if (btrmxo .gt. 0.) then
        btrfnc = (btrmxo -btrmax)/btrmxo
    end if

    btrmxo = btrmax
end subroutine betars