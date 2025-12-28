subroutine phsdrp(d1zvc1,iindx0,iindx1,iodb,ipndx1,iter,kmax,km1,km1s,kxt,kxts,nodbmx,nord,noutpt,npadd,npdel,nptmax,ntry,uphase,zvclgs,zvclg1)
    !! This subroutine picks a phase to drop from the equilibrium
    !! system. Two independent algorithms are used to find candidates,
    !! and a final choice is made from these.
    !! This subroutine is called by:
    !!   EQ6/eqcalc.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: kmax
    integer :: nodbmx
    integer :: nptmax

    integer :: noutpt

    integer :: iindx0(kmax)
    integer :: iindx1(kmax)
    integer :: ipndx1(kmax)
    integer :: iodb(nodbmx)

    integer :: iter
    integer :: km1
    integer :: km1s
    integer :: kxt
    integer :: kxts
    integer :: nord
    integer :: npadd
    integer :: npdel
    integer :: ntry

    character(len=24) :: uphase(nptmax)

    real(kind=8) :: d1zvc1(kmax)
    real(kind=8) :: zvclgs(kmax)
    real(kind=8) :: zvclg1(kmax)

    ! Local variable declarations.
    integer :: j2
    integer :: j3
    integer :: kcol
    integer :: kj
    integer :: krow
    integer :: npdc1
    integer :: npdc1a
    integer :: npdc2
    integer :: np
    integer :: ns

    integer :: ilnobl

    character(len=24) :: updc1
    character(len=24) :: updc1a
    character(len=24) :: updc2

    real(kind=8) :: crit1
    real(kind=8) :: crit1a
    real(kind=8) :: crit2
    real(kind=8) :: criter
    real(kind=8) :: zx

    npdc1 = 0
    npdc2 = 0
    crit1 = 0.
    crit2 = 0.
    updc1 = 'None'
    updc2 = 'None'

    ! Search for a candidate solid according to the first criterion
    ! (rapid decrease in mass).
    npdc1a = 0
    updc1a = 'None'
    crit1a = 0.

    do kcol = km1,kxt
        zx = zvclg1(kcol) - zvclgs(kcol)

        if (zx .ge. -1.0) then
            zx = 0.
        end if

        if (iter .le. 1) then
            zx = zvclgs(kcol) + 10.
        end if

        if (zvclgs(kcol) .le. -100.) then
            zx = zvclgs(kcol)
        end if

        if (zx .lt. crit1) then
            ns = iindx1(kcol)
            np = ipndx1(kcol)
            crit1a = zx
            npdc1a = np
            updc1a = uphase(np)

            if (np .ne. npadd) then
                crit1 = zx
                npdc1 = np
                updc1 = uphase(np)
            end if
        end if
    end do

    if (npdc1.le.0) then
        crit1 = 0.
    end if

    if (npdc1a.gt.0 .and. ntry.eq.2) then
        if (crit1a .lt. (crit1 -1.)) then
            npdc1 = npdc1a
            updc1 = updc1a
            crit1 = crit1a
        end if
    end if

    ! Search for a candidate solid according to the second criterion
    ! (negative derivative of mass).
    if (nord.gt.0 .and. ntry.le.0) then
        crit2 = -100.

        if (kxt.ge.km1 .and. kxts.ge.km1) then
            do kcol = km1,kxt
                kj = 0
                ns = iindx1(kcol)
                np = ipndx1(kcol)

                do krow = km1s,kxts
                    if (iindx0(krow) .eq. ns) then
                        kj = krow
                        go to 100
                    end if
                end do

100 continue
                if (kj .gt. 0) then
                    zx = d1zvc1(kj)

                    if (zx .lt. crit2) then
                        if (np .ne. npadd) then
                            crit2 = zx
                            npdc2 = np
                            updc2 = uphase(np)
                        end if
                    end if
                end if
            end do
        end if
    end if

    if (npdc2.le.0) then
        crit2 = 0.
    end if

    if (iodb(1) .ge. 1) then
        j2 = ilnobl(updc1)
        j3 = ilnobl(updc2)
        write (noutpt,1000) crit1,updc1(1:j2),crit2,updc2(1:j3)
1000 format(/'   --- Phase Drop Search ---',/7x,'Criterion 1. Negative divergence of log mass variable',/11x,'Value= ',e12.5,4x,/11x,'Name= ',a,/7x,'Criterion 2. Negative derivative of log mass variable',/11x,'Value= ',e12.5,4x,/11x,'Name= ',a,/)
    end if

    criter = 0.
    npdel = 0

    if (crit1 .lt. criter) then
        criter = crit1
        npdel = npdc1
    end if

    if (crit2 .lt. criter) then
        criter = crit2
        npdel = npdc2
    end if
end subroutine phsdrp