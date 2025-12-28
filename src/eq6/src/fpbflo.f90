subroutine fpbflo(al10,delxi,demop0,dlxmin,dxval0,d1emp1,d2emp1,emop,emop0,eps100,fdpe0,iemop,ier,iodb,nodbmx,nord,nordmx,noutpt,npet,npetmx,nptmax,nrd1mx,nttyo,qdump,toldl,uaqsln,ufixf,uphase,xim1,xi0,xi1,xval0,zklogu)
    !! This subroutine limits delxi by approximate position of
    !! significant maxima in the masses of non-aqeuous species that are
    !! in partial equilibrium with the aqueous solution. The numbers of
    !! moles of such phases are tracked using finite differences.
    !! This subroutine is called by:
    !!   EQ6/path.f.
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: nodbmx
    integer :: nordmx
    integer :: npetmx
    integer :: nptmax
    integer :: nrd1mx

    integer :: noutpt
    integer :: nttyo

    integer :: iemop(npetmx)
    integer :: iodb(nodbmx)

    integer :: ier
    integer :: nord
    integer :: npet

    logical :: qdump

    character(len=24) :: uphase(nptmax)
    character(len=24) :: uaqsln
    character(len=8) :: ufixf

    real(kind=8) :: demop0(nordmx,npetmx)
    real(kind=8) :: dxval0(nrd1mx)
    real(kind=8) :: d1emp1(npetmx)
    real(kind=8) :: d2emp1(npetmx)
    real(kind=8) :: emop(npetmx)
    real(kind=8) :: emop0(npetmx)
    real(kind=8) :: fdpe0(nordmx,npetmx)

    real(kind=8) :: al10
    real(kind=8) :: delxi
    real(kind=8) :: dlxmin
    real(kind=8) :: eps100
    real(kind=8) :: toldl
    real(kind=8) :: xim1
    real(kind=8) :: xi0
    real(kind=8) :: xi1
    real(kind=8) :: xval0
    real(kind=8) :: zklogu

    ! Local variable declarations.
    integer :: icount
    integer :: ilsign
    integer :: j2
    integer :: n
    integer :: np
    integer :: npe
    integer :: npej
    integer :: npj

    integer :: ilnobl

    character(len=48) :: usearch
    character(len=24) :: unam24

    real(kind=8) :: dfmrxx
    real(kind=8) :: dfmxx
    real(kind=8) :: dfmxxj
    real(kind=8) :: dpx
    real(kind=8) :: dpx0
    real(kind=8) :: mxx
    real(kind=8) :: mxxu
    real(kind=8) :: mxx0
    real(kind=8) :: tolsx
    real(kind=8) :: xtargv

    real(kind=8) :: texp

    data usearch /'a product phase maximizes                       '/

    if (nord .le. 1) then
        go to 999
    end if

    ! The search target is not where the first derivative is zero, but
    ! -eps100. The interval for convergence is (-2*eps100,0.).
    ilsign = 1
    xtargv = -eps100
    tolsx = eps100
    mxxu = texp(zklogu)

    icount = -1
100 continue
    xi1 = xi0 + delxi
    icount = icount + 1

    if (icount .gt. 50) then
        if (iodb(1).gt.0 .or. iodb(5).gt.0) then
            write (noutpt,1000)
1000 format(/3x,'A scan for product phase masses maximizing',' has failed.',/3x,'Dropping to the minimum step size.')
        end if

        delxi = dlxmin
        go to 999
    end if

    ! Make a Taylor's series expansion of the first derivatives of the
    ! numbers of moles of the phases in the ES.
    call d1ptay(delxi,demop0,d1emp1,nord,nordmx,npet,npetmx)

    ! Make a Taylor's series expansion of the numbers of moles of the
    ! phases in the ES.
    call ptaylr(delxi,demop0,emop0,emop,nord,nordmx,npet,npetmx)

    ! Find any phases whose mole numbers are decreasing. Pick the
    ! one whose mole number would be decreased the most at the
    ! current step size. Some phase types including the aqueous
    ! solution phase and any fictive fugacity-fixing phases are
    ! exempt. Phases whose mole numbers are already sufficiently
    ! small (.le. 10^zklogu at the base point) point are ignored.
    dfmxxj = 0.
    npej = 0
    npj = 0
    unam24 = ' '

    do npe = 1,npet
        np = iemop(npe)

        if (uphase(np)(1:24) .ne. uaqsln(1:24)) then
            if (uphase(np)(1:5) .ne. ufixf(1:5)) then
                ! Not an excluded phase type.
                mxx0 = emop0(npe)

                if (mxx0 .gt. mxxu) then
                    ! The base point mole number is sufficiently large
                    ! to consider.
                    mxx = emop(npe)
                    dfmxx = mxx - mxx0

                    if (dfmxx .lt. 0.) then
                        ! A decrease in the mole number is predicted.
                        dfmrxx = dfmxx/mxx0

                        if (dfmrxx .lt. -0.0001) then
                            ! The predicted decrease in mole number exceeds
                            ! a small relative amount.
                            dpx0 = demop0(1,npe)
                            dpx = d1emp1(npe)

                            if (dpx0.ge.eps100 .and. dpx.le.-eps100) then
                                ! The derivative of the mole number changes
                                ! from positive at the based point to negative
                                ! at the new point.
                                if (fdpe0(1,npe) .lt. -eps100) then
                                    ! The derivative is confirmed by the first
                                    ! finite difference.
                                    if (dfmxx .lt. dfmxxj) then
                                        ! Choose the phase with the largest predicted
                                        ! decrease in mole number.
                                        dfmxxj = dfmxx
                                        npej = npe
                                        npj = np
                                    end if
                                end if
                            end if
                        end if
                    end if
                end if
            end if
        end if
    end do

    if (npej .gt. 0) then
        qdump = .true.

        if (delxi .gt. dlxmin) then
            dpx = d1emp1(npej)

            if (abs(dpx + eps100) .gt. eps100) then
                unam24 = uphase(npj)
                xval0 = demop0(1,npej)

                do n = 1,nord - 1
                    dxval0(n) = demop0(n + 1,npej)
                end do

                if (nord .gt. 0) then
                    dxval0(nord) = 0.
                end if

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

    if (qdump) then
        if (iodb(1).gt.0 .or.iodb(5).gt.0) then
            write (noutpt,1010) xi1,delxi
1010 format(/3x,"Taylor's series predict maxima for the number of",' moles of',/3x,'the following phases at Xi= ',1pe11.4,', delxi= ',e11.4,':',/)

            do npe = 1,npet
                mxx0 = emop0(npe)
                dpx = d1emp1(npe)
                np = iemop(npe)

                if (mxx0 .gt. mxxu) then
                    if (abs(dpx + eps100) .gt. eps100) then
                        if (uphase(np)(1:24) .ne. uaqsln(1:24)) then
                            if (uphase(np)(1:5) .ne. ufixf(1:5)) then
                                j2 = ilnobl(uphase(np))
                                write (noutpt,1020) uphase(np)(1:j2)
1020 format(5x,a)
                            end if
                        end if
                    end if
                end if
            end do
        end if
    end if

999 continue
end subroutine fpbflo