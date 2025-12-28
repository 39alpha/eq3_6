subroutine tegrid(itgenf,nacdpr,narxt,nerr,noutpt,ntprmx,ntprt,nttyo,nwarn,ustrgr)
    !! This subroutine tests the contents of a "log K" temperature
    !! grid for sparse contents in any range. The response to an
    !! instance of sparseness is determined by the variable
    !! itgenf (-1 = ignore, 0 = warn, 1 = error). This subroutine
    !! need not be called if itgenf = -1.
    !! This suboutine is called by:
    !!   EQPT/rdpar.f
    !!   EQPT/pcraq.f
    !!   EQPT/pcrsg.f
    !! Principal input:
    !!   itgenf = integer flag controlling the response to a detected
    !!              incidence of sparseness (-1 = ignore, 0 = warn,
    !!              1 = error)
    !!   nacdpr = array giving the number of points containing actual
    !!              data in a temperature range for the "log K" grid
    !!              currently being examined
    !!   narxt  = array giving the number of points (actual data plus
    !!              no data) in that range
    !!   ntprt  = the number of ranges in the current "log K" grid
    !!   ustrgr = string describing what the grid represents
    !! Principal output:
    !!   nerr   = cumulative error counter
    !!   nwarn  = cumulative warning counter
    implicit none

    ! Calling sequence variable declarations.
    integer :: ntprmx

    integer :: noutpt
    integer :: nttyo

    integer :: itgenf
    integer :: nerr
    integer :: ntprt
    integer :: nwarn

    integer :: nacdpr(ntprmx)
    integer :: narxt(ntprmx)

    character(len=56) :: ustrgr

    ! Local variable declarations.
    integer :: j2
    integer :: j3
    integer :: j4
    integer :: ntpr

    logical :: qx

    integer :: ilnobl

    character(len=8) :: ux8a
    character(len=8) :: ux8b

    do ntpr = 1,ntprt
        if (narxt(ntpr) .le. 3) then
            ! For a range with up to three grid points, at least two
            ! must contain actual data. Note that a range can not
            ! be composed of fewer than two grid points.
            qx = nacdpr(ntpr) .lt. 2
        else
            ! Otherwise, at least three points containing actual data
            ! are required.
            qx = nacdpr(ntpr) .lt. 3
        end if

        if (qx .and. itgenf.ge.0) then
            write (ux8a,'(i5)') nacdpr(ntpr)
            call lejust(ux8a)
            j2 = ilnobl(ux8a)
            write (ux8b,'(i5)') ntpr
            call lejust(ux8b)
            j3 = ilnobl(ux8b)
            call lejust(ustrgr)
            j4 = ilnobl(ustrgr)

            if (itgenf .eq. 0) then
                nwarn = nwarn + 1

                if (nacdpr(ntpr) .gt. 0) then
                    write (noutpt,1000) ux8a(1:j2),ux8b(1:j3),ustrgr(1:j4)
                    write (nttyo,1000) ux8a(1:j2),ux8b(1:j3),ustrgr(1:j4)
1000 format(/' * Warning - (EQPT/tegrid) The number of grid',' points containing actual',/7x,'data is only ',a,' in',' range ',a,' of the temperature grid representing',/7x,a,'.')
                else
                    write (noutpt,1010) ux8b(1:j3),ustrgr(1:j4)
                    write (nttyo,1010) ux8b(1:j3),ustrgr(1:j4)
1010 format(/' * Warning - (EQPT/tegrid) There are no grid',' points containing actual',/7x,'data in range ',a,' of the temperature grid representing',/7x,a,'.')
                end if
            else
                nerr = nerr + 1

                if (nacdpr(ntpr) .gt. 0) then
                    write (noutpt,1030) ux8a(1:j2),ux8b(1:j3),ustrgr(1:j4)
                    write (nttyo,1030) ux8a(1:j2),ux8b(1:j3),ustrgr(1:j4)
1030 format(/' * Error - (EQPT/tegrid) The number of grid',' points containing actual',/7x,'data is only ',a,' in',' range ',a,' of the temperature grid representing',/7x,a,'.')
                else
                    write (noutpt,1040) ux8b(1:j3),ustrgr(1:j4)
                    write (nttyo,1040) ux8b(1:j3),ustrgr(1:j4)
1040 format(/' * Error - (EQPT/tegrid) There are no grid',' points containing actual',/7x,'data in range ',a,' of the temperature grid representing',/7x,a,'.')
                end if
            end if
        end if
    end do

999 continue
end subroutine tegrid