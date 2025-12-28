subroutine csfar(afrc1,morr,morr0,mwtrc,noutpt,nrc,nrctmx,nsk,nttyo,prcinf,sfcar,sfcar0,ssfcar,ureac)
    !! This subroutine calculates the surface area required to calculate
    !! the rate for the nrc-th irreversible reaction.
    !! This subroutine is called by:
    !!   EQ6/rtcalc.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: nrctmx

    integer :: noutpt
    integer :: nttyo

    integer :: nsk(nrctmx)

    integer :: nrc

    character(len=24) :: ureac(nrctmx)

    real(kind=8) :: afrc1(nrctmx)
    real(kind=8) :: morr(nrctmx)
    real(kind=8) :: morr0(nrctmx)
    real(kind=8) :: mwtrc(nrctmx)
    real(kind=8) :: sfcar(nrctmx)
    real(kind=8) :: sfcar0(nrctmx)
    real(kind=8) :: ssfcar(nrctmx)

    real(kind=8) :: prcinf

    ! Local variable declarations.
    integer :: j2

    integer :: ilnobl

    real(kind=8) :: gx

    ! Compute surface area.
    gx = morr(nrc)*mwtrc(nrc)

    if (nsk(nrc) .eq. 0) then
        ! Constant surface area (cm2).
        ! Compute the specific surface area (cm2/g).
        ssfcar(nrc) = 0.

        if (gx .gt. 0.) then
            ssfcar(nrc) = sfcar(nrc)/gx
        end if
    else if (nsk(nrc) .eq. 1) then
        ! Constant specific surface area (cm2/g).
        ! Compute the surface area (cm2).
        if (morr(nrc) .gt. 0.) then
            sfcar(nrc) = ssfcar(nrc)*gx
        else if (afrc1(nrc) .le. 0.) then
            sfcar(nrc) = 1.e+5
            j2 = ilnobl(ureac(nrc))
            write (noutpt,1000) ureac(nrc)(1:j2),sfcar(nrc)
            write (nttyo,1000) ureac(nrc)(1:j2),sfcar(nrc)
1000 format(/' * Note - (EQ6/csfar) The surface area of',' reactant ',a,/7x,'has been temporarily set to ',1pe11.4,' cm2 to permit this phase',/7x,'to begin precipitating.',' This is necessary because the surface area',/7x,' model requires the phase to have a mass, which it',/7x,'presently lacks.')
        else
            sfcar(nrc) = 0.
        end if
    else if (nsk(nrc) .eq. 2) then
        ! Constant particle number surface area growth law.
        ! Compute the surface area (cm2).
        ! Compute the specific surface area (cm2/g).
        if (morr0(nrc) .gt. 0.) then
            sfcar(nrc) = sfcar0(nrc)*(morr(nrc)/morr0(nrc))**(2./3.)
            ssfcar(nrc) = 0.

            if (gx .gt. 0) then
                ssfcar(nrc) = sfcar(nrc)/gx
            end if
        else if (afrc1(nrc) .le. 0.) then
            sfcar(nrc) = 1.e+5
            ssfcar(nrc) = prcinf
            j2 = ilnobl(ureac(nrc))
            write (noutpt,1010) ureac(nrc)(1:j2),sfcar(nrc)
            write (nttyo,1010) ureac(nrc)(1:j2),sfcar(nrc)
1010 format(/' * Note - (EQ6/csfar) The surface area of',' reactant ',a,/7x,'has been temporarily set to ',1pe11.4,' cm2 to permit this phase',/7x,'to begin precipitating.',' This is necessary because the surface area',/7x,' model requires the phase to have a mass at the',/7x,'previous point, which it lacks.')
        else
            sfcar(nrc) = 0.
        end if
    end if
end subroutine csfar