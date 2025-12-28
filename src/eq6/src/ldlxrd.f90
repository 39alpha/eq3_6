subroutine ldlxrd(delxi,dlxmin,dzvc0,fdzv0,iodb,ipndx1,kdim,km1,kmax,kxt,loph,nodbmx,nord,noutpt,nptmax,nrd1mx,nttyo,uphase,zklogu,zvec0)
    !! This subroutine limits delxi when a phase is rapidly disappearing
    !! from the equilibrium system. The particular mechanism here is
    !! only designed to slow things down enough so that other mechanisms
    !! can efficiently locate the phase disappearance boundary.
    !! This mechanism is not employed to increase information density
    !! around the boundary.
    !! This subroutine is called by:
    !!   EQ6/path.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: kmax
    integer :: nodbmx
    integer :: nptmax
    integer :: nrd1mx

    integer :: noutpt
    integer :: nttyo

    integer :: iodb(nodbmx)
    integer :: ipndx1(kmax)

    integer :: kdim
    integer :: km1
    integer :: kxt
    integer :: nord

    character(len=24) :: uphase(nptmax)

    real(kind=8) :: dzvc0(nrd1mx,kmax)
    real(kind=8) :: fdzv0(nrd1mx,kmax)
    real(kind=8) :: loph(nptmax)
    real(kind=8) :: zvec0(kmax)

    real(kind=8) :: delxi
    real(kind=8) :: dlxmin
    real(kind=8) :: zklogu

    ! Local variable declarations.
    integer :: j
    integer :: j2
    integer :: kcol
    integer :: np

    integer :: ilnobl

    character(len=24) :: unamph

    real(kind=8) :: dlxic
    real(kind=8) :: dxsv

    dxsv = delxi
    unamph = 'Error                   '

    do kcol = km1,kxt
        if (fdzv0(1,kcol) .ge. 0.) then
            go to 100
        end if

        do j = 1,nord
            if (dzvc0(j,kcol) .ge. 0.) then
                go to 100
            end if
        end do

        np = ipndx1(kcol)

        if (loph(np) .le. zklogu) then
            go to 100
        end if

        ! Here dlxic is the value of delxi at which approximately
        ! ninety per cent of the existing mass is destroyed.
        dlxic = -0.9*zvec0(kcol)/dzvc0(1,kcol)

        if (dlxic .lt. delxi) then
            delxi = dlxic
            delxi = max(delxi,dlxmin)
            unamph = uphase(np)
        end if

100 continue
    end do

    if (iodb(5) .gt. 0) then
        if (delxi .lt. dxsv) then
            j2 = ilnobl(unamph)
            write (noutpt,1000) dxsv,delxi,unamph(1:j2)
            write (nttyo,1000) dxsv,delxi,unamph(1:j2)
1000 format(/' * Note - (EQ6/ldlxrd) The step size has been cut',' from ',1pe12.5,/7x,'to ',e12.5,' to limit the dissolution',' of ',a,'.')
        end if
    end if
end subroutine ldlxrd