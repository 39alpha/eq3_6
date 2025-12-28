subroutine ndiagx(delmax,delvec,eps100,idelmx,iebal,iindx1,irdxc3,jflag,kcarb,kebal,khydr,kmax,ko2gaq,nbasp,nbtmax,nhydr,noutpt,nstmax,nttyo,screwd,uspec)
    !! This subroutine attempts to generate diagnostics if hybrid
    !! Newton-Raphson iteration fails.
    !! This subroutine is called by:
    !!   EQ3NR/eq3nr.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: kmax
    integer :: nbtmax
    integer :: nstmax

    integer :: noutpt
    integer :: nttyo

    integer :: iindx1(kmax)
    integer :: jflag(nstmax)
    integer :: nbasp(nbtmax)
    integer :: idelmx
    integer :: iebal
    integer :: irdxc3
    integer :: kcarb
    integer :: kebal
    integer :: khydr
    integer :: ko2gaq
    integer :: nhydr

    character(len=48) :: uspec(nstmax)

    real(kind=8) :: delvec(kmax)
    real(kind=8) :: delmax
    real(kind=8) :: eps100
    real(kind=8) :: screwd

    ! Local variable declarations.
    integer :: j2
    integer :: nb
    integer :: ns

    integer :: ilnobl

    logical :: qadtst
    logical :: qdltst

    real(kind=8) :: adx
    real(kind=8) :: xx

    ! Check to see if idelmx is zero. If so, diagnostics can not be
    ! generated.
    if (idelmx .eq. 0) then
        write (noutpt,1000)
        write (nttyo,1000)
1000 format(/' * Note - (EQ3NR/ndiagx) No crash diagnostics were',/7x,'generated. Look at the contents of the delvec and beta',/7x,'arrays in the iteration summary for clues.',/)

        go to 999
    end if

    ! Check to see if the ion being adjusted for electrical balance
    ! was crashing to zero. This implies that an ion of opposite
    ! charge must be used to achieve electrical balance.
    qdltst = (delvec(idelmx) + screwd) .le. eps100

    if (idelmx.eq.kebal .and. qdltst) then
        write (noutpt,1010)
        write (nttyo,1010)
1010 format(/' * Note - (EQ3NR/ndiagx) the ion adjusted for',/7x,'electrical balance is crashing to zero. Electrical',/7x,'balancing requires an ion of opposite charge.')

        go to 999
    end if

    ! Check to see if the ion associated with alkalinity balance is
    ! crashing to zero. This implies that non-carbonate alkalinity is
    ! greater than or equal to carbonate alkalinity.
    if (idelmx.eq.kcarb .and. qdltst) then
        nb = iindx1(idelmx)
        ns = nbasp(nb)

        if (jflag(ns) .eq. 7) then
            j2 = ilnobl(uspec(ns)(1:24))
            write (noutpt,1020) uspec(ns)(1:j2)
            write (nttyo,1020) uspec(ns)(1:j2)
1020 format(/' * Note - (EQ3NR/ndiagx) The basis species ',a,' is',/7x,'crashing to zero because the non-carbonate alkalinity',/7x,'exceeds the specified total alkalinity. The total',/7x,"alkalinity can't be used to constrain this species",/7x,'because of excessive interference from OH-, borate,',/7x,'phosphate, organic acids, metal hydroxy complexes, or',/7x,'the like. A total CO2 (IR) or ion chromatographic',/7x,'analysis is required instead.',/)

            go to 999
        end if
    end if

    ! Check to see if fO2 is crashing. This could be due to
    ! a bad combination of constraining the redox state by a
    ! non-fO2 option and divergence to zero of an associated
    ! ion that is constrained by electrical balance.
    qadtst = (delmax - screwd) .le. eps100

    if (idelmx.eq.ko2gaq .and. qadtst .and. iebal.gt.0) then
        if (irdxc3.eq.-1 .and. nbasp(iebal).eq.nhydr) then
            xx = 0.1*screwd
            adx = abs(delvec(khydr))

            if (adx .ge. xx) then
                write (noutpt,1030)
                write (nttyo,1030)
1030 format(/' * Note - (EQ3NR/ndiagx) The fO2 is crashing,',/7x,'probably because a bad electrical balance constraint',/7x,'on H+ is causing the concentration of that ion',/7x,'to crash to zero.',/)

                go to 999
            end if
        end if

        if (irdxc3 .eq. 1) then
            write (noutpt,1040)
            write (nttyo,1040)
1040 format(/' * Note - (EQ3NR/ndiagx) The fO2 is crashing,',/7x,'probably because of a bad constraint on one of the',/7x,'aqueous species appearing in the redox reaction that',/7x,'is being used to constrain the fO2.',/)

            go to 999
        end if
    end if

999 continue
end subroutine ndiagx