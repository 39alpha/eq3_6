subroutine gtemp(afcnst,al10,iopt,jtemp,noptmx,noutpt,nttkmx,nttyo,rconst,rtcnst,tempc,tempcb,tempk,time1,ttk,xi1)
    !! This subroutine computes the temperature (tempc) as a function of
    !! reaction progress (xi1) or time (time1).
    !! This subroutine is called by:
    !!   EQ6/eq6.f
    !!   EQ6/tpadv.f
    !! Principal input:
    !!   jtemp  = temperature tracking flag
    !!   ttk    = temperature tracking coefficients
    !!   time1  = time variable
    !!   xi1    = reaction progress variable
    !! Principal output:
    !!   tempc  = temperature, C
    !!   tempk  = temperature, K
    implicit none

    ! Calling sequence variable declarations.
    integer :: noptmx
    integer :: nttkmx

    integer :: iopt(noptmx)
    integer :: jtemp
    integer :: noutpt
    integer :: nttyo

    real(kind=8) :: ttk(nttkmx)
    real(kind=8) :: afcnst
    real(kind=8) :: al10
    real(kind=8) :: rconst
    real(kind=8) :: rtcnst
    real(kind=8) :: tempc
    real(kind=8) :: tempcb
    real(kind=8) :: tempk
    real(kind=8) :: time1
    real(kind=8) :: xi1

    ! Local variable declarations.
    real(kind=8) :: tmpmin
    real(kind=8) :: tmpmax

    ! Allowed range of the temperature (C):
    data tmpmin,tmpmax /-273.15,1400./

    if (jtemp .eq. 0) then
        ! Constant temperature.
        tempc = tempcb
    else if (jtemp .eq. 1) then
        ! Linear tracking in Xi.
        tempc = tempcb + ttk(1)*xi1
    else if (jtemp .eq. 2) then
        ! Linear tracking in time.
        tempc = tempcb + ttk(1)*time1
    else if (jtemp .eq. 3) then
        ! Fluid mixing tracking.
        tempc = ( tempcb*ttk(1) + xi1*ttk(2) )/( xi1 + ttk(1) )
    end if

    tempk = tempc + 273.15
    rtcnst = 0.001*rconst*tempk
    afcnst = al10*rtcnst

    ! Stop if tempc is less than tmpmin or greater than tmpmax.
    if (tempc.lt.tmpmin .or. tempc.gt.tmpmax) then
        write (noutpt,1010) tempc,tmpmin,tmpmax,xi1
        write (nttyo,1010) tempc,tmpmin,tmpmax,xi1
1010 format(/' * Error - (EQ6/gtemp) The calculated temperature',/7x,'is ',g10.3," C, outside the code's built-in allowed",/7x,'range of ',g10.3,'-',g10.3,' C. Other limits associated',/7x,'with the supporting data file may also apply.',/7x,'The current value of reaction progress is ',g12.5,' moles.')

        if (iopt(2) .gt. 0) then
            write (noutpt,1020) time1
            write (nttyo,1020) time1
1020 format(7x,'The current time value is ',g10.3,' seconds.')
        end if

        stop
    end if

999 continue
end subroutine gtemp