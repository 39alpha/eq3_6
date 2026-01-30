subroutine gpress(iopt,jpress,noptmx,noutpt,nptkmx,nttyo,presg,presh,press,pressb,time1,ptk,xi1)
    !! This subroutine computes the pressure (press) as a function of
    !! reaction progress (xi1) or time (time1).
    !! This subroutine is called by:
    !!   EQ6/eq6.f
    !!   EQ6/tpadv.f
    !! Principal input:
    !!   jpress = pressure tracking flag
    !!   ptk    = pressure tracking coefficients
    !!   time1  = time variable
    !!   xi1    = reaction progress variable
    !! Principal output:
    !!   press  = pressure, bars
    implicit none

    ! Calling sequence variable declarations.
    integer :: noptmx
    integer :: nptkmx

    integer :: noutpt
    integer :: nttyo

    integer :: iopt(noptmx)
    integer :: jpress

    real(kind=8) :: ptk(nptkmx)
    real(kind=8) :: presg
    real(kind=8) :: presh
    real(kind=8) :: press
    real(kind=8) :: pressb
    real(kind=8) :: time1
    real(kind=8) :: xi1

    ! Local variable declarations.
    real(kind=8) :: presmx

    ! Maximum allowed pressure (bars):
    data presmx /10000./

    if (jpress .eq. 0) then
        ! The pressure corresponds to the data file reference presssure
        ! curve.
        press = presg
    else if (jpress .eq. 1) then
        ! The pressure corresponds to the 1.013-bar/steam-saturation
        ! presssure curve.
        press = presh
    else if (jpress .eq. 2) then
        ! Constant pressure.
        press = pressb
    else if (jpress .eq. 3) then
        ! Linear tracking in reaction progress.
        press = pressb + ptk(1)*xi1
    else if (jpress .eq. 4) then
        ! Linear tracking in time.
        press = pressb + ptk(1)*time1
    end if

    ! Stop if the pressure isn't greater than 0.
    if (press .le. 0) then
        write (noutpt,1010) press
        write (nttyo,1010) press
1010 format(/' * Error - (EQ6/gpress) The pressure must be',/7x,'greater than zero. The current pressure is ',1pg12.5,' bars.')

        stop
    end if

    ! Stop if press is greater than presmx.
    if (press .gt. presmx) then
        write (noutpt,1020) press,presmx,xi1
        write (nttyo,1020) press,presmx,xi1
1020 format(/' * Error - (EQ6/gpress) The calculated pressure',/7x,'is ',1pg12.5," C, greater than the code's built-in",/7x,'maximum value of ',g12.5,' bars. Other limits',/7x,'associated with the supporting data file may also apply.',/7x,'The current value of reaction progress is ',g12.5,' moles.')

        if (iopt(2) .gt. 0) then
            write (noutpt,1030) time1
            write (nttyo,1030) time1
1030 format(7x,'The current time value is ',g12.5,' seconds.')
        end if

        stop
    end if
end subroutine gpress
