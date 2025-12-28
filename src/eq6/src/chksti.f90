subroutine chksti(akmat0,drer0,drir0,deltim,delxi,dlxmin,fdri0,fdrr0,iodb,jreac,kly,kmax,nodbmx,kord,nord,noutpt,npts,nrct,nrctmx,nrd1mx,nttyo,prcinf,qriinf,rirec0,smp100,time0,time1,xi0,xi1)
    !! This subroutine checks the sign of the time increment.
    !! This subroutine is called by:
    !!   EQ6/path.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: kmax
    integer :: nodbmx
    integer :: nrctmx
    integer :: nrd1mx

    integer :: noutpt
    integer :: nttyo

    integer :: iodb(nodbmx)
    integer :: jreac(nrctmx)

    integer :: kly
    integer :: kord
    integer :: nord
    integer :: npts
    integer :: nrct

    logical :: qriinf

    real(kind=8) :: akmat0(nrd1mx,nrd1mx)
    real(kind=8) :: drer0(nrd1mx,nrctmx)
    real(kind=8) :: drir0(nrd1mx)
    real(kind=8) :: fdri0(nrd1mx)
    real(kind=8) :: fdrr0(nrd1mx,nrctmx)

    real(kind=8) :: deltim
    real(kind=8) :: delxi
    real(kind=8) :: dlxmin
    real(kind=8) :: rirec0
    real(kind=8) :: prcinf
    real(kind=8) :: smp100
    real(kind=8) :: time0
    real(kind=8) :: time1
    real(kind=8) :: xi0
    real(kind=8) :: xi1

    ! Local variable declarations.
    real(kind=8) :: dx
    real(kind=8) :: dxsv

100 continue
    xi1 = xi0 + delxi
    call timeca(deltim,delxi,drir0,iodb,nodbmx,nord,noutpt,nrd1mx,nttyo,prcinf,qriinf,rirec0,time0,time1)

    if (deltim .le. smp100) then
        dxsv = delxi
        dx = 0.25*delxi
        delxi = max(dx,dlxmin)

        if (iodb(1).gt.0 .or. iodb(5).gt.0) then
            write (noutpt,1000) deltim
            write (nttyo,1000) deltim
1000 format(/' * Warning - (EQ6/chksti) The calculated time',/7x,"increment is ',1pe11.4. This increment can't be zero",/7x,"or negative, and shouldn't be too close to such a",/7x,'condition. The finite difference representation of',/7x,'the inverse rate may be losing accuracy.')
        end if

        if (dxsv .gt. dlxmin) then
            if (iodb(1).gt.0 .or. iodb(5).gt.0) then
                write (noutpt,1010) dxsv,delxi
                write (nttyo,1010) dxsv,delxi
1010 format(7x,'The step size will be cut from ',1pe11.4,' to ',1pe11.4,'.')
            end if
        else
            nord = nord - 1
            kord = nord
            npts = nord + 1
            kly = 6

            if (iodb(1).gt.0 .or. iodb(5).gt.0) then
                write (noutpt,1020) nord
                write (nttyo,1020) nord
1020 format(7x,'The step size is already at the  minimum value.',/7x,'Cutting the order to ',i2,'.')
            end if

            call rderiv(akmat0,drer0,drir0,fdri0,fdrr0,jreac,nord,nrct,nrctmx,nrd1mx)
        end if

        go to 100
    end if
end subroutine chksti