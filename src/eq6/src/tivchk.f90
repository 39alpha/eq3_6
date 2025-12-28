subroutine tivchk(deltim,delxi,qtvchk,time1,time0,timemx,tiplol,tiplot,tiprnl,tiprnt,tolxst)
    !! This subroutine checks to make sure that the calculated
    !! time does not exceed any specified limits such as the maximum
    !! time. This routine should be called only if delxi is less than
    !! or equal to the minimum step size, dlxmin.
    !! This subroutine is called by:
    !!   EQ6/path.f
    !!   EQ6/eqshel.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    logical :: qtvchk

    real(kind=8) :: deltim
    real(kind=8) :: delxi
    real(kind=8) :: time1
    real(kind=8) :: time0
    real(kind=8) :: timemx
    real(kind=8) :: tiplol
    real(kind=8) :: tiplot
    real(kind=8) :: tiprnl
    real(kind=8) :: tiprnt
    real(kind=8) :: tolxst

    ! Local variable declarations.
    ! None
    qtvchk = .false.

    ! Check the next print point in time.
    if (((time1 - tiprnt)/tiprnt) .gt. tolxst) then
        time1 = tiprnt
        deltim = tiprnt - time0
        qtvchk = .true.
    end if

    ! Check the next print point in log time.
    if (((time1 - tiprnl)/tiprnl) .gt. tolxst) then
        time1 = tiprnl
        deltim = tiprnl - time0
        qtvchk = .true.
    end if

    ! Check the next plot point in time.
    if (((time1 - tiplot)/tiplot) .gt. tolxst) then
        time1 = tiplot
        deltim = tiplot - time0
        qtvchk = .true.
    end if

    ! Check the next plot point in log time.
    if (((time1 - tiplol)/tiplol) .gt. tolxst) then
        time1 = tiplol
        deltim = tiplol - time0
        qtvchk = .true.
    end if

    ! Check the maximum time.
    if (((time1 - timemx)/timemx) .gt. tolxst) then
        time1 = timemx
        deltim = timemx - time0
        qtvchk = .true.
    end if
end subroutine tivchk