subroutine fpbnpp(affp,affp0,aftarg,daffp0,delxi,dlxmin,dxval0,eps100,iodb,iopt,jpflag,nodbmx,noptmx,nord,nordmx,noutpt,npchk,npt,nptmax,nrd1mx,nttyo,tolaft,tolsat,uphase,xi0,xi1,xval0)
    !! This subroutine finds the phase boundary at which a phase appears
    !! in the equilibrium system (ES). Note that the phase affinities
    !! are tracked using finite differences.
    !! This subroutine is called by:
    !!   EQ6/path.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: nodbmx
    integer :: nordmx
    integer :: noptmx
    integer :: nptmax
    integer :: nrd1mx

    integer :: noutpt
    integer :: nttyo

    integer :: iodb(nodbmx)
    integer :: iopt(noptmx)
    integer :: jpflag(nptmax)
    integer :: npchk(nptmax)

    integer :: nord
    integer :: npt

    character(len=24) :: uphase(nptmax)

    real(kind=8) :: affp(nptmax)
    real(kind=8) :: affp0(nptmax)
    real(kind=8) :: daffp0(nordmx,nptmax)
    real(kind=8) :: dxval0(nrd1mx)

    real(kind=8) :: aftarg
    real(kind=8) :: delxi
    real(kind=8) :: dlxmin
    real(kind=8) :: eps100
    real(kind=8) :: tolaft
    real(kind=8) :: tolsat
    real(kind=8) :: xi0
    real(kind=8) :: xi1
    real(kind=8) :: xval0

    ! Local variable declarations.
    character(len=48) :: usearch
    character(len=24) :: unam24

    integer :: icount
    integer :: ier
    integer :: ilsign
    integer :: j2
    integer :: kpsst
    integer :: n
    integer :: np
    integer :: npsst

    integer :: ilnobl

    real(kind=8) :: afx
    real(kind=8) :: dxsv
    real(kind=8) :: tolsx
    real(kind=8) :: xtargv
    real(kind=8) :: xval

    data usearch /'a phase just supersaturates                     '/

    if (nord .le. 0) then
        go to 999
    end if

    ! The search target is not where the affinity is zero, but aftarg.
    ! The interval for convergence is (tolsat,tolsst).  The target
    ! (aftarg) is the midpoint of this interval.
    ilsign = -1
    xtargv = aftarg
    tolsx = tolaft
    dxsv = delxi

    icount = -1
100 continue
    xi1 = xi0 + delxi
    icount = icount + 1

    if (icount .gt. 50) then
        if (iodb(1).gt.0 .or. iodb(5).gt.0) then
            write (noutpt,1000)
1000 format(/3x,'A scan for product phases just supersaturating',' has failed.',/3x,'Dropping to the minimum step size.')
        end if

        delxi = dlxmin
        go to 999
    end if

    ! Estimate the affinities of formation of the various phases from
    ! Taylor's series expansions.
    call ataylr(delxi,daffp0,nord,nordmx,npt,nptmax,affp0,affp)

    ! Find any predicted supersaturations. Note that these count only
    ! if they exceed the tolerance tolsat.
    xval = 0.
    npsst = 0
    kpsst = 0
    unam24 = ' '

    do np = 1,npt
        afx = affp(np)

        if (afx .ge. tolsat) then
            if (jpflag(np).le.0 .and. jpflag(np).ne.-1) then
                if (npchk(np) .ne. 1) then
                    kpsst = kpsst + 1

                    if (afx .gt. xval) then
                        xval = affp(np)
                        npsst = np
                        unam24 = uphase(np)
                    end if
                end if
            end if
        end if
    end do

    if (kpsst .gt. 0) then
        if (delxi .gt. dlxmin) then
            if (abs(xval - xtargv) .gt. tolsx) then
                xval0 = affp0(npsst)

                do n = 1,nord
                    dxval0(n) = daffp0(n,npsst)
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
        if (kpsst .gt. 0) then
            write (noutpt,1010) xi1,delxi
1010 format(/3x,"Taylor's series predict the appearance of the",' following',/3x,'phases in the ES at Xi= ',1pe11.4,', delxi= ',e11.4,':',/)

            do np = 1,npt
                afx = affp(np)

                if (abs(afx - xtargv) .le. tolsx) then
                    if (jpflag(np).le.0 .and. jpflag(np).ne.-1) then
                        if (npchk(np) .ne. 1) then
                            j2 = ilnobl(uphase(np))
                            write (noutpt,1020) uphase(np)(1:j2)
1020 format(5x,a)
                        end if
                    end if
                end if
            end do
        end if
    end if

    if (iopt(3) .ge. 1) then
        delxi = dxsv
    end if

999 continue
end subroutine fpbnpp