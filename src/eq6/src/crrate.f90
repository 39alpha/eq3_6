subroutine crrate(act,afrc1,cdac,csigma,eps100,fkrc,idirec,imchmx,imech,iodb,jreac,morr,ndac,ndact,ndctmx,nodbmx,noutpt,nrc,nrctmx,nrk,nstmax,nttyo,rk,rreac1,rrelr1,rrxfi1,rtcnst,sfcar,udac,ureac)
    !! This subroutine calculates the rate (relative or absolute) from
    !! the specified rate law for the nrc-th irreversible reaction.
    !! This subroutine is called by:
    !!   EQ6/rtcalc.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: imchmx
    integer :: ndctmx
    integer :: nodbmx
    integer :: nrctmx
    integer :: nstmax

    integer :: noutpt
    integer :: nttyo

    integer :: nrc

    integer :: idirec(nrctmx)
    integer :: imech(2,nrctmx)
    integer :: iodb(nodbmx)
    integer :: jreac(nrctmx)
    integer :: ndac(ndctmx,imchmx,2,nrctmx)
    integer :: ndact(imchmx,2,nrctmx)
    integer :: nrk(2,nrctmx)

    character(len=24) :: udac(ndctmx,imchmx,2,nrctmx)
    character(len=24) :: ureac(nrctmx)

    real(kind=8) :: act(nstmax)
    real(kind=8) :: afrc1(nrctmx)
    real(kind=8) :: cdac(ndctmx,imchmx,2,nrctmx)
    real(kind=8) :: csigma(imchmx,2,nrctmx)
    real(kind=8) :: fkrc(nrctmx)
    real(kind=8) :: morr(nrctmx)
    real(kind=8) :: rk(imchmx,2,nrctmx)
    real(kind=8) :: rreac1(nrctmx)
    real(kind=8) :: rrelr1(nrctmx)
    real(kind=8) :: rrxfi1(imchmx,nrctmx)
    real(kind=8) :: sfcar(nrctmx)

    real(kind=8) :: eps100
    real(kind=8) :: rtcnst

    ! Local variable declarations.
    integer :: i
    integer :: j2
    integer :: j3
    integer :: n
    integer :: ns

    integer :: ilnobl

    logical :: qovstp

    character(len=8) :: ux8

    real(kind=8) :: affr
    real(kind=8) :: aprod
    real(kind=8) :: efx
    real(kind=8) :: fs
    real(kind=8) :: rx

    ! Calculate the net rate (relative or absolute) for each
    ! kinetically-governed reaction. Use either a forward rate constant
    ! plus affinity form or a backward rate constant plus affinity form.
    ! This subroutine does not use a form involving a forward rate
    ! constant plus backward rate constant.
    rreac1(nrc) = 0.
    rrelr1(nrc) = 0.

    do i = 1,imchmx
        rrxfi1(i,nrc) = 0.
    end do

    ! Check for exhausted or saturated reactant.
    if (jreac(nrc) .gt. 0) then
        go to 999
    end if

    ! Check the sign of the affinity. For minerals, the affinity
    ! here is the affinity to dissolve, so:
    !   Negative value -> precipitation
    !   Positive value -> dissolution
    ! Some reactants (special reactants, aqeuous species, and
    ! gases) currently do not have actual affinities associated
    ! with them. For these, the affinity is currently taken as
    ! +9999999.
    affr = afrc1(nrc)
    qovstp = affr.lt.0. .and. nrk(2,nrc).eq.0 .and. jreac(nrc).le.0

    if ( iodb(2) .ge. 2 ) then
        j2 = ilnobl(ureac(nrc))

        if (affr.gt.0. .or. qovstp) then
            write (noutpt,1010) ureac(nrc)(1:j2)
1010 format(/' --- Rate calculations for destruction of ',a,' ---',/)
        else if (affr .lt. 0.) then
            write (noutpt,1020) ureac(nrc)(1:j2)
1020 format(/' --- Rate calculations for formation of ',a,' ---',/)
        else
            write (noutpt,1030) sfcar(nrc)
1030 format(/5x,'Surface area= ',e12.5,' cm2.',/)
        end if
    end if

    fs = fkrc(nrc)*sfcar(nrc)

    if (idirec(nrc) .eq. 1) then
        ! Calculate the net rate using the forward (dissolution) form
        ! of rate law.
        if (nrk(1,nrc) .eq. 0) then
            j2 = ilnobl(ureac(nrc))
            write (noutpt,1100) ureac(nrc)(1:j2)
            write (nttyo,1100) ureac(nrc)(1:j2)
1100 format(/' * Error - (EQ6/crrate) Programming error trap:',/7x,'The forward rate law code (nrk(1,nrc)) is zero for',/7x,'reactant ',a,'.')

            stop
        else if (nrk(1,nrc) .eq. 1) then
            ! Specified relative rate (arbitrary kinetics).
            rrelr1(nrc) = rk(1,1,nrc)
            rrxfi1(1,nrc) = rrelr1(nrc)
        else if (nrk(1,nrc) .eq. 2) then
            ! Transition state theory (more properly, the TST-like form).
            ! The net rate is equal to an absolute forward rate minus an
            ! absolute backward rate. However, that is obvious only in the
            ! equivalent two rate constant form.
            rreac1(nrc) = 0.

            do i = 1,imech(1,nrc)
                efx = affr/(csigma(i,1,nrc)*rtcnst)
                aprod = 1.

                do n = 1,ndact(i,1,nrc)
                    ns = ndac(n,i,1,nrc)
                    aprod = aprod*act(ns)**cdac(n,i,1,nrc)
                end do

                rx = rk(i,1,nrc)*fs*aprod*(1.0 - exp(-efx))
                rreac1(nrc) = rreac1(nrc) + rx
                rrxfi1(i,nrc) = rx
            end do
        else if (nrk(1,nrc) .eq. 3) then
            ! Specified rate (no affinity dependence).
            rreac1(nrc) = fs*rk(1,1,nrc)
            rrxfi1(1,nrc) = rreac1(nrc)
        else
            ! Unrecognized rate law code.
            j2 = ilnobl(ureac(nrc))
            write (ux8,'(i5)') nrk(1,nrc)
            call lejust(ux8)
            j3 = ilnobl(ux8)
            write (noutpt,1110) ureac(nrc)(1:j2),ux8(1:j3)
            write (nttyo,1110) ureac(nrc)(1:j2),ux8(1:j3)
1110 format(/' * Error - (EQ6/crrate) The reaction for the',/7x,'destruction of reactant ',a,' has an unrecognized',/7x,'forward rate law code of ',a,'.')

            stop
        end if
    else
        ! Calculate the net rate using the backward (formation) form
        ! of rate law.
        ! Note: All rates are expressed in this code as the net of
        ! forward minus backward rates, so a postive formation rate
        ! is expressed as a negative number.
        if (nrk(2,nrc) .eq. 0) then
            ! Instantaneous equilibrium. Escape for case of a
            ! supersaturated reactant that is not in the matrix.
        else if (nrk(2,nrc) .eq. 1) then
            ! Specified relative rate (aribitrary kinetics).
            rrelr1(nrc) = -rk(1,2,nrc)
            rrxfi1(1,nrc) = rrelr1(nrc)
        else if (nrk(2,nrc) .eq. 2) then
            ! Transition state theory (more properly, the TST-like form).
            rreac1(nrc) = 0.

            do i = 1,imech(2,nrc)
                ! Note: if affr < 0, then efx < 0, affac < 0, and rreac1 < 0.
                efx = affr/(csigma(i,2,nrc)*rtcnst)
                aprod = 1.

                do n = 1,ndact(i,2,nrc)
                    ns = ndac(n,i,2,nrc)
                    aprod = aprod*act(ns)**cdac(n,i,2,nrc)
                end do

                rx = rk(i,2,nrc)*fs*aprod*(1.0 - exp(-efx))
                rreac1(nrc) = rreac1(nrc) + rx
                rrxfi1(i,nrc) = rx
            end do
        else if (nrk(2,nrc) .eq. 3) then
            ! Linear rate law.
            rreac1(nrc) = -fs*rk(1,2,nrc)
            rrxfi1(1,nrc) = rreac1(nrc)
        else
            j2 = ilnobl(ureac(nrc))
            write (ux8,'(i5)') nrk(2,nrc)
            call lejust(ux8)
            j3 = ilnobl(ux8)
            write (noutpt,250) ureac(nrc)(1:j2),ux8(1:j3)
            write (nttyo,250) ureac(nrc)(1:j2),ux8(1:j3)
250 format(/' * Error - (EQ6/crrate) The reaction for the',/7x,'destruction of reactant ',a,' has an unrecognized',/7x,'backward rate law code of ',a,'.')

            stop
        end if
    end if

999 continue
end subroutine crrate