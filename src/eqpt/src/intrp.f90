subroutine intrp(aamatr,apr,avgrid,cof,eps100,gmmatr,ipivot,narxmx,narxt,noutpt,ntprmx,ntprt,nttyo,tempc,tempcs,tmpcmx,xvec,yvec)
    !! This subroutine fits interpolating polynomials to data (avgrid)
    !! on a temperature grid (tempc). The grid is divided into ranges.
    !! A separate polynomial is fitted to the data in each range.
    !! Its coefficients are obtained in the cof array. The coefficients
    !! for all ranges are returned in the apr array.
    !! This subroutine is called by:
    !!   EQPT/pcraq.f
    !!   EQPT/pcrsg.f
    !!   EQPT/wrpar.f
    !! Principal input:
    !!   avgrid = array containing the data on the temperature grid
    !! Principal output:
    !!   apr    = array of polynomial coefficients (for all ranges)
    !!   avgrid = array containing the data on the temperature grid
    !!   narxmx = the maximum number of points or coefficients per
    !!              temperature range
    !!   narxt  = the actual number of points or coefficients per
    !!              temperature range
    !!   nptrmx = the maximum number of temperature ranges
    !!   nptrt  = the actual number of temperature ranges
    !! Workspace:
    !!   aamatr = matrix used to calculate the polynomial coefficients
    !!   cof    = array of fitted polynomial coefficients (for a
    !!              single temperature range)
    !!   xvec   = array of scaled temperatures corresponding to the
    !!              data in the yvec array
    !!   yvec   = array of data to be fitted (for a single temperature
    !!              range)
    implicit none

    ! Calling sequence variable declarations.
    integer :: narxmx
    integer :: ntprmx

    integer :: noutpt
    integer :: nttyo

    integer :: ntprt

    integer :: ipivot(narxmx)
    integer :: narxt(ntprmx)

    real(kind=8) :: apr(narxmx,ntprmx)
    real(kind=8) :: avgrid(narxmx,ntprmx)
    real(kind=8) :: tempc(narxmx,ntprmx)
    real(kind=8) :: tempcs(narxmx,ntprmx)
    real(kind=8) :: tmpcmx(ntprmx)
    real(kind=8) :: aamatr(narxmx,narxmx)
    real(kind=8) :: gmmatr(narxmx,narxmx)
    real(kind=8) :: cof(narxmx)
    real(kind=8) :: xvec(narxmx)
    real(kind=8) :: yvec(narxmx)

    real(kind=8) :: eps100

    ! Local variable declarations.
    integer :: ier
    integer :: n
    integer :: nmax
    integer :: npft
    integer :: ntpr

    real(kind=8) :: tvecmx

    ! Initialize apr to 0.
    nmax = narxmx*ntprmx
    call initaz(apr,nmax)

    ! Loop on temperature ranges.
    do ntpr = 1,ntprt
        ! Put all real values (9999999. indicates no data) in the yvec
        ! array. Put the corresponding scaled temperature values in
        ! the xvec array.
        npft = 0

        do n = 1,narxt(ntpr)
            if (avgrid(n,ntpr) .lt. 9999999.) then
                npft = npft + 1
                xvec(npft) = tempcs(n,ntpr)
                yvec(npft) = avgrid(n,ntpr)
            end if
        end do

        ! Here npft is the number of usable values.
        if (npft .le. 0) then
            ! There are no usable values.
            apr(1,ntpr) = 9999999.
            go to 110
        end if

        if (npft .gt. 1) then
            ! Check for constant yvec.
            do n = 2,npft
                if(abs(yvec(n) -yvec(1)) .gt. eps100) then
                    go to 100
                end if
            end do

            npft = 1
100 continue
        end if

        if (npft .le. 1) then
            apr(1,ntpr) = yvec(1)
            go to 110
        end if

        ! Fit the polynomial.
        call polfit(aamatr,cof,gmmatr,ier,ipivot,narxmx,npft,noutpt,nttyo,xvec,yvec)

        if (ier .gt. 0) then
            write (noutpt,1010) ntpr,tempc(1,ntpr),tempc(narxt(ntpr),ntpr)
            write (nttyo,1010) ntpr,tempc(1,ntpr),tempc(narxt(ntpr),ntpr)
1010 format(/' * Error - (EQPT/intrp) Could not compute the',/7x,'coefficients of an interpolating polynomial in',/7x,'temperature range ',i2,' (',f6.2,' - ',f6.2,' C.')

            stop
        end if

        ! Rescale the coefficients.
        tvecmx = tmpcmx(ntpr)

        ! Calling sequence substitutions:
        !   tvecmx for avxmax
        !   cof for avy
        !   cof for avys
        !   npft for nmax
        call rscaly(tvecmx,cof,cof,eps100,npft)

        ! Store the fitted coefficients for this range in the apr array.
        do n = 1,npft
            apr(n,ntpr) = cof(n)
        end do

110 continue
    end do
end subroutine intrp