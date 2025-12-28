subroutine rtcalc(act,afrc1,cdac,csigma,eps100,fkrc,idirec,imchmx,imech,iodb,iopt,jcode,jreac,morr,morr0,mwtrc,ndac,ndact,ndctmx,nodbmx,noptmx,nord,noutpt,nrk,nrct,nrctmx,nsk,nstmax,nttyo,prcinf,prminf,qriinf,rirec1,rk,rreac1,rrelr1,rtcnst,rrxfi1,sfcar,sfcar0,ssfcar,udac,ureac)
    !! This subroutine calculates the relative and absolute rates of
    !! the nrc-th irreversible reaction. This rate is computed from the
    !! specified rate expression.
    !! Under iopt(2) = 1, xi1 is now defined as the sum of the absolute
    !! values of the progress variables for all kinetically governed
    !! (nrk(1,nrc).gt.1 or nrk(2,nrc).gt.1) reactions. The inverse rate
    !! is the inverse of the sum of the absolute values of the rates of
    !! these reactions.
    !! This subroutine is called by:
    !!   EQ6/path.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: imchmx
    integer :: ndctmx
    integer :: nodbmx
    integer :: noptmx
    integer :: nrctmx
    integer :: nstmax

    integer :: idirec(nrctmx)
    integer :: imech(2,nrctmx)
    integer :: iodb(nodbmx)
    integer :: iopt(noptmx)
    integer :: jcode(nrctmx)
    integer :: jreac(nrctmx)
    integer :: ndac(ndctmx,imchmx,2,nrctmx)
    integer :: ndact(imchmx,2,nrctmx)
    integer :: nrk(2,nrctmx)
    integer :: nsk(nrctmx)

    integer :: noutpt
    integer :: nttyo

    integer :: nord
    integer :: nrct

    logical :: qriinf

    character(len=24) :: udac(ndctmx,imchmx,2,nrctmx)
    character(len=24) :: ureac(nrctmx)

    real(kind=8) :: act(nstmax)
    real(kind=8) :: afrc1(nrctmx)
    real(kind=8) :: cdac(ndctmx,imchmx,2,nrctmx)
    real(kind=8) :: csigma(imchmx,2,nrctmx)
    real(kind=8) :: fkrc(nrctmx)
    real(kind=8) :: morr(nrctmx)
    real(kind=8) :: morr0(nrctmx)
    real(kind=8) :: mwtrc(nrctmx)
    real(kind=8) :: rk(imchmx,2,nrctmx)
    real(kind=8) :: rreac1(nrctmx)
    real(kind=8) :: rrelr1(nrctmx)
    real(kind=8) :: rrxfi1(imchmx,nrctmx)
    real(kind=8) :: sfcar(nrctmx)
    real(kind=8) :: sfcar0(nrctmx)
    real(kind=8) :: ssfcar(nrctmx)

    real(kind=8) :: eps100
    real(kind=8) :: prcinf
    real(kind=8) :: prminf
    real(kind=8) :: rirec1
    real(kind=8) :: rtcnst

    ! Local variable declarations.
    integer :: id
    integer :: j2
    integer :: nrc

    integer :: ilnobl

    logical :: qovstp

    real(kind=8) :: affr
    real(kind=8) :: trate

    ! Loop over all kinetically-governed reactions. Calculate the
    ! relevant surface areas for surface-area-controlled reactions.
    do nrc = 1,nrct
        call csfar(afrc1,morr,morr0,mwtrc,noutpt,nrc,nrctmx,nsk,nttyo,prcinf,sfcar,sfcar0,ssfcar,ureac)
    end do

    ! Determine the rate law form (forward or backward) to use for
    ! each reactant.
    do nrc = 1,nrct
        affr = afrc1(nrc)
        qovstp = affr.lt.0. .and. nrk(2,nrc).eq.0 .and. jreac(nrc).le.0

        if (affr .ge. 0.) then
            ! The affinity favors the forward direction.
            idirec(nrc) = 1

            if (nrk(1,nrc) .eq. -1) then
                idirec(nrc) = 2
            end if
        else
            ! The affinity favors the backward direction.
            idirec(nrc) = 2

            if (nrk(2,nrc).eq.-1 .or. qovstp) then
                idirec(nrc) = 1
            end if
        end if
    end do

    ! Loop over all kinetically-governed reactions. Calculate the
    ! net rate (relative or absolute) for each one.
    do nrc = 1,nrct
        call crrate(act,afrc1,cdac,csigma,eps100,fkrc,idirec,imchmx,imech,iodb,jreac,morr,ndac,ndact,ndctmx,nodbmx,noutpt,nrc,nrctmx,nrk,nstmax,nttyo,rk,rreac1,rrelr1,rrxfi1,rtcnst,sfcar,udac,ureac)
    end do

    ! If in time mode, calculate the inverse rate. Calculate relative
    ! rates from the corresponding absolute rates. If any primary
    ! calculated rates were relative rates, calculate the corresponding
    ! absolute rates.
    if (iopt(2) .gt. 0) then
        ! Compute the total rate.
        trate = 0.

        do nrc = 1,nrct
            id = idirec(nrc)

            if (nrk(id,nrc) .gt. 1) then
                trate = trate + abs(rreac1(nrc))
            end if
        end do

        ! Compute absolute rates for reactants constrained by
        ! relative rates.
        do nrc = 1,nrct
            id = idirec(nrc)

            if (nrk(id,nrc) .eq. 1) then
                rreac1(nrc) = rrelr1(nrc)*trate
            end if
        end do

        ! Compute the inverse rate.
        qriinf = .false.

        if (trate .le. prminf) then
            qriinf = .true.
            nord = 0
            rirec1 = prcinf
        else
            rirec1 = 1./trate
        end if

        ! Compute relative rates for reactants constrained by
        ! absolute rates.
        do nrc = 1,nrct
            id = idirec(nrc)

            if (nrk(id,nrc) .gt. 1) then
                rrelr1(nrc) = rreac1(nrc)*rirec1
            end if
        end do
    end if

    ! Check if the rate and the affinity are compatible.
    do nrc = 1,nrct
        if (jreac(nrc).le.0 .and. jcode(nrc) .le. 1) then
            affr = afrc1(nrc)

            if (affr.gt.0. .and. rrelr1(nrc).lt.0.) then
                j2 = ilnobl(ureac(nrc))
                write (noutpt,1000) ureac(nrc)(1:j2)
                write (nttyo,1000) ureac(nrc)(1:j2)
1000 format(/' * Warning - (EQ6/rtcalc) The calculated rate',' for reactant',/7x,a,' is less than zero, but the',' affinity',/7x,'is greater than zero:')

                write (noutpt,1010) rreac1(nrc),rrelr1(nrc),affr
                write (nttyo,1010) rreac1(nrc),rrelr1(nrc),affr
1010 format(/9x,'Rate = ',1pe12.5,' mol/d',/9x,'Relative rate= ',e12.5,' mol/mol',/9x,'Affinity= ',e12.5,' kcal')

                write (noutpt,1020)
                write (nttyo,1020)
1020 format(/7x,'Check the signs of the rate constants.',/)
            end if

            qovstp = affr.lt.0. .and. nrk(2,nrc).eq.0 .and.    jreac(nrc).le.0

            if (affr.lt.0. .and. .not.qovstp .and.      rrelr1(nrc).gt.0.) then
                j2 = ilnobl(ureac(nrc))
                write (noutpt,1050) ureac(nrc)(1:j2)
                write (nttyo,1050) ureac(nrc)(1:j2)
1050 format(/' * Warning - (EQ6/rtcalc) The calculated rate',' for reactant',/7x,a,' is greater than zero, but the',' affinity',/7x,'is less than zero:')

                write (noutpt,1010) rreac1(nrc),rrelr1(nrc),affr
                write (nttyo,1010) rreac1(nrc),rrelr1(nrc),affr
                write (noutpt,1020)
                write (nttyo,1020)
            end if
        end if
    end do

    if (iodb(2) .gt. 1) then
        ! Write a summary of calculated rate information on the
        ! output file.
        do nrc = 1,nrct
            if (jreac(nrc) .le. 0) then
                j2 = ilnobl(ureac(nrc))
                write (noutpt,1070) ureac(nrc)(1:j2),rreac1(nrc),rrelr1(nrc),affr
1070 format(/5x,'Calculated Rate Summary',/5x,a,/7x,'Rate = ',1pe12.5,' mol/d',/7x,'Relative rate= ',e12.5,' mol/mol',/7x,'Affinity= ',e12.5,' kcal')
            end if
        end do
    end if
end subroutine rtcalc