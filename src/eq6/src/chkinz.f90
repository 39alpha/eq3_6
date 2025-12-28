subroutine chkinz(ier,imchmx,imech,iopt,jcode,kmax,kxt,nelect,noptmx,noutpt,no2gaq,nrct,nrctmx,nrk,nstmax,nttyo,rkb,ureac,uspeca,uzveci,zvclgi)
    !! This subroutine checks the code input for various kinds of
    !! errors and inconsistencies. Here ier accumulates the
    !! number of errors caught.
    !! This subroutine is called by:
    !!   EQ6/eq6.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: imchmx
    integer :: kmax
    integer :: noptmx
    integer :: nrctmx
    integer :: nstmax

    integer :: noutpt
    integer :: nttyo

    integer :: imech(2,nrctmx)
    integer :: iopt(noptmx)
    integer :: jcode(nrctmx)
    integer :: nrk(2,nrctmx)

    integer :: ier
    integer :: kxt
    integer :: nelect
    integer :: no2gaq
    integer :: nrct

    character(len=48) :: uspeca(nstmax)
    character(len=48) :: uzveci(kmax)
    character(len=24) :: ureac(nrctmx)

    real(kind=8) :: rkb(imchmx,2,nrctmx)
    real(kind=8) :: zvclgi(kmax)

    ! Local variable declarations.
    integer :: i
    integer :: jlen
    integer :: j2
    integer :: kcol
    integer :: nrc

    integer :: ilnobl

    character(len=56) :: uspn56

    ier = 0

    ! Rate law data.
    !   nrk(1,nrc)= forward rate law code:
    !     -1 = Use the backward rate law form (legal only if
    !          nrk(2,nrc) = 2)
    !      0 = Illegal value
    !      1 = Relative rate
    !      2 = Transition state theory net rate
    !      3 = Linear rate
    !   For the case nrk(1,nrc) = 2, the reactant must be
    !   either a pure mineral (jcode = 0) or a solid solution
    !   (jcode = 1).
    !   nrk(2,nrc)= backward rate law code:
    !     -1 = Use the forward rate law form (legal only if
    !          nrk(1,nrc) = 2)
    !      0 = No rate law specified; the reaction may be controlled
    !          by partial equilibrium.
    !      1 = Relative rate
    !      2 = Transition state theory net rate
    !      3 = Linear rate
    !   For the case nrk(2,nrc) = 2, the reactant must be
    !   either a pure mineral (jcode = 0) or a solid solution
    !   (jcode = 1)
    do nrc = 1,nrct
        j2 = ilnobl(ureac(nrc))

        if (iopt(2) .le. 0) then
            ! Check for rate law choices that can't be used without
            ! an explicit time frame.
            if (nrk(2,nrc) .ge. 2) then
                write (noutpt,1000) ureac(nrc)(1:j2),nrk(2,nrc)
                write (nttyo,1000) ureac(nrc)(1:j2),nrk(2,nrc)
1000 format(/' * Error - (EQ6/chkinz) The forward rate law code',/7x,'for ',a,' is ',i2,', which is not valid',/7x,'unless iopt(2) = 1 (the model has a time frame).')

                ier = ier + 1
            end if

            if (nrk(2,nrc) .ge. 2) then
                write (noutpt,1010) ureac(nrc)(1:j2),nrk(2,nrc)
                write (nttyo,1010) ureac(nrc)(1:j2),nrk(2,nrc)
1010 format(/' * Error - (EQ6/chkinz) The backward rate law',' code',/7x,'for ',a,' is ',i2,', which is not valid',/7x,'unless iopt(2) = 1 (the model has a time frame).')

                ier = ier + 1
            end if
        end if

        ! Check for invalid forward and backward rate law combinations.
        ! The following combinations are invalid:
        !   nrk(1,nrc) = 1 and nrk(2,nrc) = 1
        !   nrk(1,nrc) = 3 and nrk(2,nrc) = 3
        !   nrk(1,nrc) = -1 and nrk(2,nrc) is not 2
        !   nrk(2,nrc) = 0
        if (nrk(1,nrc).eq.1 .and. nrk(2,nrc).eq.1) then
            write (noutpt,1020) ureac(nrc)(1:j2)
            write (nttyo,1020) ureac(nrc)(1:j2)
1020 format(' * Error - (EQ6/chkinz) The forward and backward',/7x,'rate law codes for ',a,' are both 1, which is not a',/7x,'valid combination.')

            ier = ier + 1
        end if

        if (nrk(1,nrc).eq.3 .and. nrk(2,nrc).eq.3) then
            write (noutpt, 1030) ureac(nrc)(1:j2)
            write (nttyo, 1030) ureac(nrc)(1:j2)
1030 format(' * Error - (EQ6/chkinz) The forward and backward',/7x,'rate law codes for ',a,' are both 3, which is not a',/7x,'valid combination.')

            ier = ier + 1
        end if

        if (nrk(1,nrc) .eq. -1) then
            if (nrk(2,nrc) .ne. 2) then
                write (noutpt, 1040) ureac(nrc)(1:j2)
                write (nttyo, 1040) ureac(nrc)(1:j2)
1040 format(' * Error - (EQ6/chkinz) The forward rate law code',/7x,'for ',a,' is -1, which is not a valid value  unless',/7x,'the backward rate law code is 2.')

                ier = ier + 1
            end if
        end if

        if (nrk(1,nrc) .eq. 0) then
            write (noutpt,1050) ureac(nrc)(1:j2)
            write (nttyo,1050) ureac(nrc)(1:j2)
1050 format(' * Error - (EQ6/chkinz) The forward rate law code',/7x,'for ',a,' is 0, which is not a valid value.')

            ier = ier + 1
        end if

        ! Check the parameters for the forward rate laws.
        if (nrk(1,nrc) .eq. 2) then
            ! Transition state theory.
            do i = 1,imech(1,nrc)
                if (rkb(i,1,nrc) .lt. 0.) then
                    write (noutpt,1060) i,rkb(i,1,nrc),ureac(nrc)(1:j2),nrk(1,nrc)
                    write (nttyo,1060) i,rkb(i,1,nrc),ureac(nrc)(1:j2),nrk(1,nrc)
1060 format (/' * Error - (EQ6/chkinz) Forward rate constant',/7x,i2,' has a value of ',1pe11.4,' for reactant'        /7x,a,'. A negative rate constant is not valid for a',/7x,'forward rate code of ',i2,'.')

                    ier = ier + 1
                end if
            end do

            if (jcode(nrc) .gt. 1) then
                write (noutpt,1070) nrk(1,nrc),ureac(nrc)(1:j2)
                write (nttyo,1070) nrk(1,nrc),ureac(nrc)(1:j2)
1070 format(/' * Error - (EQ6/chkinz) Have specified a',/7x,'forward rate law code of ',i2,' for reactant',/7x,a,'. This rate law incorporates an affinity',/7x,'dependence and can only be used for a reactant which',/7x,'is a mineral or solid solution.')

                ier = ier + 1
            end if
        else if (nrk(1,nrc) .eq. 3) then
            ! Linear rate law
            if (rkb(1,1,nrc) .lt. 0.) then
                write (noutpt,1060) i,rkb(1,1,nrc),ureac(nrc)(1:j2),nrk(1,nrc)
                write (nttyo,1060) i,rkb(1,1,nrc),ureac(nrc)(1:j2),nrk(1,nrc)
                ier = ier + 1
            end if
        end if

        ! Check the parameters for the backward rate laws.
        if (nrk(2,nrc) .eq. -1) then
            ! Use the forward rate law for the backward case.
            if (nrk(1,nrc) .ne. 2) then
                write (noutpt,1080) ureac(nrc)(1:j2),nrk(1,nrc)
                write (nttyo,1080) ureac(nrc)(1:j2),nrk(1,nrc)
1080 format(' * Error - (EQ6/chkinz) The backward rate law',/7x,'code for ',a,' is  -1, which is invalid because',/7x,'the forward rate law code is ',i2,'.')

                ier = ier + 1
            end if
        else if (nrk(2,nrc) .eq. 2) then
            ! Transition state rate law.
            do i = 1,imech(2,nrc)
                if (rkb(i,2,nrc) .lt. 0.) then
                    write (noutpt,1090) i,rkb(i,2,nrc),ureac(nrc)(1:j2),nrk(2,nrc)
                    write (nttyo,1090) i,rkb(i,2,nrc),ureac(nrc)(1:j2),nrk(2,nrc)
1090 format (/' * Error - (EQ6/chkinz) Backward rate constant',/7x,i2,' has a value of ',1pe11.4,' for reactant'        /7x,a,'. A negative rate constant is not valid for a',/7x,'forward rate code of ',i2,'.')

                    write (noutpt,1060) ureac(nrc)(1:j2)
                    write (nttyo,1060) ureac(nrc)(1:j2)
                    ier = ier + 1
                end if
            end do

            if (jcode(nrc) .gt. 1) then
                write (noutpt,1100) nrk(2,nrc),ureac(nrc)(1:j2)
                write (nttyo,1100) nrk(2,nrc),ureac(nrc)(1:j2)
1100 format(/' * Error - (EQ6/chkinz) Have specified a',/7x,'backward rate law code of ',i2,' for reactant',/7x,a,'. This rate law incorporates an affinity',/7x,'dependence and can only be used for a reactant which',/7x,'is a mineral or solid solution.')

                write (noutpt,1070)
                write (nttyo,1070)
                ier = ier + 1
            end if
        else if (nrk(2,nrc) .eq. 3) then
            ! Linear  rate law.
            if (rkb(1,2,nrc) .lt. 0.) then
                write (noutpt,1090) i,rkb(1,2,nrc),ureac(nrc)(1:j2),nrk(2,nrc)
                write (nttyo,1090) i,rkb(1,2,nrc),ureac(nrc)(1:j2),nrk(2,nrc)
                ier = ier + 1
            end if

            if (jcode(nrc) .gt. 1) then
                if (nrk(1,nrc).eq.3 .or. nrk(1,nrc).eq.1) then
                    write (noutpt,1110) ureac(nrc)(1:j2),nrk(1,nrc)
                    write (nttyo,1110) ureac(nrc)(1:j2),nrk(1,nrc)
1110 format(' * Error - (EQ6/chkinz) The forward rate law',/7x,'code for reactant ',a,' may not be 3 if the',/7x,'backward rate law code is',i2,'. This is permitted',/7x,'only for a reactant which is a pure mineral or',/7x,'a solid solution..')

                    ier = ier + 1
                end if
            end if
        end if
    end do

    ! Matrix variable values.
    !    zvclgi = logarithmic basis variable, starting estimate
    !       or pick-up value
    do kcol = 1,kxt
        if (zvclgi(kcol) .le. -99999.) then
            if (uzveci(kcol)(1:8).ne.uspeca(no2gaq)(1:8) .and.      uzveci(kcol)(1:8).ne.uspeca(nelect)(1:8)) then
                call fmspnx(jlen,uzveci(kcol),uspn56)
                write (noutpt,1200) uspn56(1:jlen)
                write (nttyo,1200) uspn56(1:jlen)
1200 format(/' * Error - (EQ6/chkinz) The basis species ',a,/7x,'has a log number of moles (zvclgi) value less than',' or equal to -99999.',/7x,'on the input file.')

                ier = ier + 1
            end if
        end if
    end do
end subroutine chkinz