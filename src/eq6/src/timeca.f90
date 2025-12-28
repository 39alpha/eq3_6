subroutine timeca(deltim,delxi,drir0,iodb,nodbmx,nord,noutpt,nrd1mx,nttyo,prcinf,qriinf,rirec0,time0,time1)
    !! This subroutine calculates the interval (deltim) of model time
    !! that corresponds to the reaction progress interval (delxi).
    !! This subroutine is called by:
    !!   EQ6/chksti.f
    !!   EQ6/eqshel.f
    !!   EQ6/path.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: nodbmx
    integer :: nrd1mx

    integer :: noutpt
    integer :: nttyo

    integer :: iodb(nodbmx)
    integer :: nord

    logical :: qriinf

    real(kind=8) :: drir0(nrd1mx)
    real(kind=8) :: deltim
    real(kind=8) :: delxi
    real(kind=8) :: prcinf
    real(kind=8) :: rirec0
    real(kind=8) :: time0
    real(kind=8) :: time1

    ! Local variable declarations.
    integer :: i
    integer :: n

    real(kind=8) :: dxp

    real(kind=8) :: fctrl

    if (qriinf) then
        write (noutpt,1000)
        write (nttyo,1000)
1000 format(/' * Note - (EQ6/timeca) Have encountered an infinite',/7x,'time interval.')

        time1 = prcinf
        deltim = prcinf
        go to 999
    end if

    deltim = rirec0*delxi

    if (nord .gt. 0) then
        dxp = delxi

        do n = 1,nord
            i = n + 1
            dxp = dxp*delxi
            deltim = deltim + ( drir0(n)*dxp )/fctrl(i)
        end do
    end if

    time1 = time0 + deltim

    if (deltim.le.0. .and. delxi.gt.0.) then
        if (iodb(2) .gt. 0) then
            write (noutpt,1010) deltim,nord
            write (nttyo,1010) deltim,nord
1010 format(/' * Note - (EQ6/timeca) Have encountered a',/7x,'time increment of ',1pe12.5,' seconds calculated for',/7x,'order ',i2,'.')

            write (noutpt,1020) delxi,rirec0
1020 format(7x,'delxi= ',1pe12.5,/7x,'rirec0= ',e12.5)

            if (nord .gt. 0) then
                write (noutpt,1030) (drir0(n), n = 1,nord)
            end if

1030 format(7x,'drir0: ',/12x,1pe12.5,3x,e12.5,3x,e12.5,/12x,e12.5,3x,e12.5,3x,e12.5)
        end if
    end if

999 continue
end subroutine timeca