subroutine nlkffg(axlks,cess,cdrs,iffg,ifrn1,ifrn2,jffg,jpflag,jsflag,mwtsp,narxmx,ncmpr,ndrs,ndrsmx,ndrsr,ness,nessmx,nessr,nffg,nffgmx,nfrn1,nfrn2,ngrn1,ngrn2,noutpt,nphasx,npt,nptmax,nst,nstmax,ntpr,ntprmx,nttyo,qcntmp,uffg,ufixf,uphase,uspec,vosp0,xlkffg)
    !! This subroutine sets up ficitive minerals, each of which is used
    !! to fixing the fugacity of a specified gas. The fugacity is
    !! actually fixed only if the corresponding ficitive mineral is in
    !! equilibrium with the aqueous system. Otherwise, the fugacity may
    !! be less than the specified fugacity value.
    !! This subroutine is called by:
    !!   EQ6/eq6.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: narxmx
    integer :: ndrsmx
    integer :: nessmx
    integer :: nffgmx
    integer :: nptmax
    integer :: nstmax
    integer :: ntprmx

    integer :: noutpt
    integer :: nttyo

    integer :: iffg(nffgmx)
    integer :: jffg(nffgmx)
    integer :: jpflag(nptmax)
    integer :: jsflag(nstmax)
    integer :: ncmpr(2,nptmax)
    integer :: ness(nessmx)
    integer :: nessr(2,nstmax)
    integer :: ndrs(ndrsmx)
    integer :: ndrsr(2,nstmax)
    integer :: nphasx(nstmax)

    integer :: ifrn1
    integer :: ifrn2
    integer :: nffg
    integer :: nfrn1
    integer :: nfrn2
    integer :: ngrn1
    integer :: ngrn2
    integer :: npt
    integer :: nst
    integer :: ntpr

    logical :: qcntmp

    character(len=48) :: uspec(nstmax)
    character(len=24) :: uffg(nffgmx)
    character(len=24) :: uphase(nptmax)
    character(len=8) :: ufixf

    real(kind=8) :: axlks(narxmx,ntprmx,nstmax)
    real(kind=8) :: cess(nessmx)
    real(kind=8) :: cdrs(ndrsmx)
    real(kind=8) :: mwtsp(nstmax)
    real(kind=8) :: vosp0(nstmax)
    real(kind=8) :: xlkffg(nffgmx)

    ! Local variable declarations.
    integer :: i
    integer :: j
    integer :: j2
    integer :: k
    integer :: n
    integer :: nerr
    integer :: nerr1
    integer :: ng
    integer :: nn
    integer :: np
    integer :: nrf1
    integer :: nrf2
    integer :: nr1
    integer :: nr2
    integer :: ns
    integer :: nsl
    integer :: nss

    integer :: ilnobl

    character(len=24) :: ux24

    nerr = 0
    np = npt
    ns = nst
    nfrn1 = nst + 1
    ifrn1 = npt + 1

    do n = 1,nffg
        do nss = ngrn1,ngrn2
            if (uffg(n)(1:24) .eq. uspec(nss)(1:24)) then
                go to 105
            end if
        end do

        j2 = ilnobl(uffg(n))
        write (noutpt,1000) uffg(n)(1:j2)
        write (nttyo,1000) uffg(n)(1:j2)
1000 format (/' * Error - (EQ6/nlkffg) ',a,' was not on the',/7x,'data file, so its fugacity can not be fixed.')

        nerr = nerr + 1
        go to 200

105 continue
        nerr1 = 0

        if (npt .ge. nptmax) then
            j2 = ilnobl(uffg(n))
            write (noutpt,1010) nptmax,uffg(n)(1:j2)
            write (nttyo,1010) nptmax,uffg(n)(1:j2)
1010 format (/' * Error - (EQ6/nlkffg) The maximum ',i3,' phases',/7x,'would be exceeded by by creating a fugacity fixing',' phase',/7x,'for ',a,'. Increase the dimensioning variable nptpar.')

            nerr1 = nerr1 + 1
            nerr = nerr + 1
        end if

        if (nst .ge. nstmax) then
            j2 = ilnobl(uffg(n))
            write (noutpt,1020) nstmax,uffg(n)(1:j2)
            write (nttyo,1020) nstmax,uffg(n)(1:j2)
1020 format (/' * Error - (EQ6/nlkffg) The maximum ',i3,' species',/7x,'would be exceeded by creating a fugacity fixing phase',/7x,'for ',a,'. Increase the dimensioning parameter nstpar.')

            nerr1 = nerr1 + 1
            nerr = nerr + 1
        end if

        if (nerr1 .gt. 0) then
            go to 200
        end if

        np = np + 1
        npt = np
        ns = ns + 1
        nsl = nst
        nst = ns
        iffg(n) = ns
        nphasx(ns) = np
        ncmpr(1,np) = ns
        ncmpr(2,np) = ns
        ng = nss - ngrn1 + 1
        jffg(n) = ng
        jpflag(np) = 0
        jsflag(ns) = 0

        ux24 = ufixf(1:5)
        ux24(6:24) = uspec(nss)(1:19)
        uspec(ns)(1:24) = ux24
        uspec(ns)(25:48) = ux24
        uphase(np) = ux24

        mwtsp(ns) = mwtsp(nss)
        vosp0(ns) = 0.

        nr1 = nessr(1,nss)
        nr2 = nessr(2,nss)
        nrf1 = nessr(2,nsl) + 1
        nrf2 = nrf1 + nr2 - nr1

        if (nrf2 .gt. nessmx) then
            ux24 = uspec(ns)(1:24)
            j2 = ilnobl(ux24)
            write (noutpt,1030) nessmx,uspec(ns)(1:j2)
            write (nttyo,1030) nessmx,uspec(ns)(1:j2)
1030 format (/' * Error - (EQ6/nlkffg) The maximum ',i5,' entries',/7x,'in the cess/ness arrays would be exceeded by creating',/7x,'the species ',a,'. Increase the dimensioning',/7x,'parameter nesspa.')

            nerr = nerr + 1
            go to 115
        end if

        nessr(1,ns) = nrf1
        nessr(2,ns) = nrf2
        k = nr1 - 1

        do nn = nrf1,nrf2
            k = k + 1
            cess(nn) = cess(k)
            ness(nn) = ness(k)
        end do

115 continue
        nr1 = ndrsr(1,nss)
        nr2 = ndrsr(2,nss)
        nrf1 = ndrsr(2,nsl) + 1
        nrf2 = nrf1 + nr2 - nr1

        if (nrf2 .gt. ndrsmx) then
            ux24 = uspec(ns)(1:24)
            j2 = ilnobl(ux24)
            write (noutpt,1040) ndrsmx,uspec(ns)(1:j2)
            write (nttyo,1040) ndrsmx,uspec(ns)(1:j2)
1040 format (/' * Error - (EQ6/nlkffg) The maximum ',i5,' entries',/7x,'in the cdrs/ndrs arrays would be exceeded by creating',/7x,'the species ',a,'. Increase the dimensioning',/7x,'parameter ndrspa.')

            nerr = nerr + 1
            go to 125
        end if

        ndrsr(1,ns) = nrf1
        ndrsr(2,ns) = nrf2
        k = nr1 - 1

        do nn = nrf1,nrf2
            k = k + 1
            cdrs(nn) = cdrs(k)
            ndrs(nn) = ndrs(k)
        end do

        ndrs(nrf1) = ns

125 continue
        do i = 1,narxmx
            do j = 1,ntprmx
                axlks(i,j,ns) = axlks(i,j,nss)
            end do
        end do

        ! Alter log K values of fixed fugacity solid.
        do j = 1,ntprmx
            if (axlks(1,j,ns) .lt. 9999999.) then
                axlks(1,j,ns) = axlks(1,j,ns) + xlkffg(n)
            else
                if (j.eq.ntpr .or. .not.qcntmp) then
                    j2 = ilnobl(uffg(n))
                    write (noutpt,1050) uffg(n)(1:j2),j
                    write (nttyo,1050) uffg(n)(1:j2),j
1050 format (/" * Error - (EQ6/nlkffg) Can't fix the fugacity",/7x,' of ',a,' in temperature range ',i2,' because',/7x,'actual thermodynamic data are lacking.')

                    nerr = nerr + 1
                end if
            end if
        end do

200 continue
    end do

    nfrn2 = nst
    ifrn2 = npt

    if (nerr .gt. 0) then
        stop
    end if
end subroutine nlkffg