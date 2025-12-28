subroutine fpbdpp(delxi,demop0,dlxmin,dxval0,emop,emop0,eps100,iemop,iodb,iopt,nodbmx,noptmx,nord,nordmx,noutpt,npet,nrd1mx,npetmx,nptmax,nttyo,uphase,xi0,xi1,xval0)
    !! This subroutine finds the boundary at which a phase or phases
    !! disappear from the equilibrium system (ES). The numbers of moles
    !! of phases in the ES are tracked using finite differences.
    !! This subroutine is called by:
    !!   EQ6/path.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: nodbmx
    integer :: noptmx
    integer :: nordmx
    integer :: npetmx
    integer :: nptmax
    integer :: nrd1mx

    integer :: noutpt
    integer :: nttyo

    integer :: iemop(npetmx)
    integer :: iodb(nodbmx)
    integer :: iopt(noptmx)

    integer :: nord
    integer :: npet

    character(len=24) :: uphase(nptmax)

    real(kind=8) :: demop0(nordmx,npetmx)
    real(kind=8) :: dxval0(nrd1mx)
    real(kind=8) :: emop(npetmx)
    real(kind=8) :: emop0(npetmx)

    real(kind=8) :: delxi
    real(kind=8) :: dlxmin
    real(kind=8) :: eps100
    real(kind=8) :: xi0
    real(kind=8) :: xi1
    real(kind=8) :: xval0

    ! Local variable declarations.
    integer :: icount
    integer :: ier
    integer :: ilsign
    integer :: j2
    integer :: kpdis
    integer :: n
    integer :: np
    integer :: npe
    integer :: npedis

    integer :: ilnobl

    character(len=48) :: usearch
    character(len=24) :: unam24

    real(kind=8) :: dxsv
    real(kind=8) :: mxx
    real(kind=8) :: mxx0
    real(kind=8) :: tolsx
    real(kind=8) :: xtargv
    real(kind=8) :: xval

    data usearch /'a product phase vanishes                        '/

    if (nord .le. 0) then
        go to 999
    end if

    ! The search target is not where the number of moles is
    ! zero, but -eps100. The interval for convergence is (-2*eps100,0.).
    ilsign = 1
    xtargv = -eps100
    tolsx = eps100
    dxsv = delxi

    icount = -1
100 continue
    xi1 = xi0 + delxi
    icount = icount + 1

    if (icount .gt. 50) then
        if (iodb(1).gt.0 .or. iodb(5).gt.0) then
            write (noutpt,1000)
1000 format(/3x,'A scan for product phases vanishing',' has failed.',/3x,'Dropping to the minimum step size.')
        end if

        delxi = dlxmin
        go to 999
    end if

    ! Make a Taylor's series expansion of the number of moles of the
    ! phases in the ES.
    call ptaylr(delxi,demop0,emop0,emop,nord,nordmx,npet,npetmx)

    ! Find any disappearing phases. Note that these count only if the
    ! predicted number of moles is less than or equal to zero.
    xval = 1.0
    npedis = 0
    kpdis = 0
    unam24 = ' '

    do npe = 1,npet
        mxx = emop(npe)

        if (mxx .le. 0) then
            mxx0 = emop0(npe)

            if (mxx.lt.mxx0 .and. mxx0.gt.eps100) then
                kpdis = kpdis + 1

                if (mxx .lt. xval) then
                    xval = mxx
                    npedis = npe
                    np = iemop(npe)
                    unam24 = uphase(np)
                end if
            end if
        end if
    end do

    if (kpdis .gt. 0) then
        if (delxi .gt. dlxmin) then
            if (abs(xval - xtargv) .gt. tolsx) then
                xval0 = emop0(npedis)

                do n = 1,nord
                    dxval0(n) = demop0(n,npedis)
                end do

                call search(delxi,dlxmin,dxval0,eps100,ier,ilsign,iodb,nodbmx,nord,noutpt,nrd1mx,nttyo,tolsx,unam24,usearch,xtargv,xval0)

                if (ier .le. 0) then
                    go to 100
                end if

                if (ier .ge. 2) then
                    ! Note: if ier = 1, the returned "safe" value of delxi
                    ! is used.
                    delxi = dlxmin
                end if

                go to 999
            end if
        end if
    end if

    if (iodb(1).gt.0 .or. iodb(5).gt.0 .or. iopt(3).ge.1) then
        if (kpdis .gt. 0) then
            write (noutpt,1010) xi1,delxi
1010 format(/3x,"Taylor's series predict disappearance of the",' following',/3x,'ES phases at Xi= ',1pe11.4,', delxi= ',e11.4,':',/)

            do npe = 1,npet
                mxx = emop(npe)

                if (abs(mxx - xtargv) .le. tolsx) then
                    mxx0 = emop0(npe)

                    if (mxx.lt.mxx0 .and. mxx0.gt.eps100) then
                        np = iemop(npe)
                        j2 = ilnobl(uphase(np))
                        write (noutpt,1020) uphase(np)(1:j2)
1020 format(5x,a)
                    end if
                end if
            end do
        end if
    end if

    if (iopt(3) .ge. 1) then
        delxi = dxsv
    end if

999 continue
end subroutine fpbdpp