subroutine evratc(eact,hact,iact,imchmx,imech,nrct,nrctmx,nrk,rk,rkb,rtcnst,tempk,trkb)
    !! This subroutine evaluates rate constants as functions of
    !! temperature. There are two alternative treatments. One assumes
    !! a constant activation energy. The other assumes a constant
    !! activation enthalpy.
    !! This subroutine is called by:
    !!   EQ6/eq6.f
    !!   EQ6/tpadv.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: imchmx
    integer :: nrctmx

    integer :: iact(imchmx,2,nrctmx)
    integer :: imech(2,nrctmx)
    integer :: nrk(2,nrctmx)

    integer :: nrct

    real(kind=8) :: eact(imchmx,2,nrctmx)
    real(kind=8) :: hact(imchmx,2,nrctmx)
    real(kind=8) :: rk(imchmx,2,nrctmx)
    real(kind=8) :: rkb(imchmx,2,nrctmx)
    real(kind=8) :: trkb(imchmx,2,nrctmx)

    real(kind=8) :: rtcnst
    real(kind=8) :: tempk

    ! Local variable declarations.
    integer :: i
    integer :: nrc

    real(kind=8) :: tfx
    real(kind=8) :: tk1

    do nrc = 1,nrct
        ! Forward rate laws (forward equates to dissolution,
        ! dissociation, or other means of destruction).
        if (nrk(1,nrc) .ge. 1) then
            do i = 1,imech(1,nrc)
                if (iact(i,1,nrc) .eq. 0) then
                    ! No temperature dependence.
                    rk(i,1,nrc) = rkb(i,1,nrc)
                else if (iact(i,1,nrc) .eq. 1) then
                    ! Constant activation energy.
                    tk1 = trkb(i,1,nrc) + 273.15
                    tfx = (eact(i,1,nrc)/rtcnst) * ((tempk/tk1) - 1.)
                    rk(i,1,nrc) = rkb(i,1,nrc)*exp(tfx)
                else if (iact(i,1,nrc) .eq. 2) then
                    ! Constant activation enthalpy.
                    tk1 = trkb(i,1,nrc) + 273.15
                    tfx = (hact(i,1,nrc)/rtcnst) * ((tempk/tk1) - 1.)
                    rk(i,1,nrc) = rkb(i,1,nrc)*(tempk/tk1)*exp(tfx)
                end if
            end do
        end if

        ! Backward rate laws (backward equates to precipitation,
        ! association, or other means of formation).
        if (nrk(2,nrc) .ge. 1) then
            do i = 1,imech(2,nrc)
                if (iact(i,2,nrc) .eq. 0) then
                    ! No temperature dependence.
                    rk(i,2,nrc) = rkb(i,2,nrc)
                else if (iact(i,2,nrc) .eq. 1) then
                    ! Constant activation energy.
                    tk1 = trkb(i,2,nrc) + 273.15
                    tfx = (eact(i,2,nrc)/rtcnst) * ((tempk/tk1) - 1.)
                    rk(i,2,nrc) = rkb(i,2,nrc)*exp(tfx)
                else if (iact(i,2,nrc) .eq. 2) then
                    ! Constant activation enthalpy.
                    tk1 = trkb(i,2,nrc) + 273.15
                    tfx = (hact(i,2,nrc)/rtcnst) * ((tempk/tk1) - 1.)
                    rk(i,2,nrc) = rkb(i,2,nrc)*(tempk/tk1)*exp(tfx)
                end if
            end do
        end if
    end do
end subroutine evratc