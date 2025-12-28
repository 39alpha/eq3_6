subroutine rdgrid(ndat0s,ndbmax,ndbptg,ndbptl,qend,qerr,q500nd,udbval,xdbval)
    !! This subroutine reads data from the standard log K temperature
    !! grid. Other kinds of data are read from such grids.
    !! This suboutine is called by:
    !!   EQPT/rdpar.f
    !!   EQPT/pcraq.f
    !!   EQPT/pcrsg.f
    !! Principal input:
    !!   ndbptg = the number of points on the grid
    !!   ndbptl = the maximum number of points on a single line
    !! Principal output:
    !!   qend   = logical flag, = .true. if end-of-file was encountered
    !!   qerr   = logical flag, = .true. if a read format error was
    !!              encountered
    !!   udbval = string array, the equivalent of xdbval
    !!   xdbval = array containing the data read from the grid
    implicit none

    ! Calling sequence variable declarations.
    integer :: ndbmax

    integer :: ndat0s

    integer :: ndbptg
    integer :: ndbptl

    logical :: qend
    logical :: qerr
    logical :: q500nd

    character(len=16) :: udbval(ndbmax)

    real(kind=8) :: xdbval(ndbmax)

    ! Local variable declarations.
    integer :: i
    integer :: ii
    integer :: j
    integer :: jj
    integer :: k
    integer :: n

    character(len=80) :: uline
    character(len=80) :: ulbufa
    character(len=80) :: ulbufb
    character(len=16) :: ux16

    qend = .false.
    qerr = .false.

    ! Note: the following coding generally assumes that a blank
    ! space trails all but the last number on a line. If a blank
    ! space is not found, it is assumed that a number ends after
    ! four decimal places.
    n = 0

    do i = 1,ndbptg/ndbptl
        read (ndat0s,'(a)',end=990,err=995) uline
        ulbufa = uline

        do k = 1,ndbptl
            call lejust(ulbufa)
            ii = index(ulbufa,' ')
            jj = index(ulbufa,'.')

            if (jj .eq. 0) then
                jj = 80
            else
                jj = jj + 5
            end if

            ii = min(ii,jj)

            if (ii .gt. 1) then
                ux16 = ulbufa(1:ii - 1)
                call locase(ux16)
                udbval(n + k) = ux16

                if (ux16(1:7) .eq. 'no_data') then
                    xdbval(n + k) = 9999999.
                else
                    read (ux16,'(f10.4)',err=995) xdbval(n + k)
                end if
            else
                udbval(n + k) = ' '
                xdbval(n + k) = 0.
            end if

            if (k .lt. ndbptl) then
                ulbufb = ulbufa(ii:80)
                ulbufa = ulbufb
            end if
        end do

        ! read (uline,udbfmt,err=995) (xdbval(n + k), k = 1,ndbptl)
        ! Note: udbfmt = '( (5x,6f10.4) )' with the "6" repeat
        ! count replaced by ndbptl.
        n = n + ndbptl
    end do

    j = mod(ndbptg,ndbptl)

    if (j .gt. 0) then
        read (ndat0s,'(a)',end=990,err=995) uline
        ulbufa = uline

        do k = 1,j
            call lejust(ulbufa)
            ii = index(ulbufa,' ')
            jj = index(ulbufa,'.')

            if (jj .eq. 0) then
                jj = 80
            else
                jj = jj + 5
            end if

            ii = min(ii,jj)

            if (ii .gt. 1) then
                ux16 = ulbufa(1:ii - 1)
                call locase(ux16)
                udbval(n + k) = ux16

                if (ux16(1:7) .eq. 'no_data') then
                    xdbval(n + k) = 9999999.
                else
                    read (ux16,'(f10.4)',err=995) xdbval(n + k)
                end if
            else
                udbval(n + k) = ' '
                xdbval(n + k) = 0.
            end if

            if (k .lt. j) then
                ulbufb = ulbufa(ii:80)
                ulbufa = ulbufb
            end if
        end do

        ! read (uline,udbfmt,err=995) (xdbval(n + k), k = 1,j)
        n = n + j
    end if

    if (q500nd) then
        ! Apply the following filter. Treat values of "500.0000" as
        ! "no data". This accommodates older EQ3/6 data files which
        ! use "500.0000" to mean no data.
        do i = 1,ndbptg
            if (xdbval(i).gt.499.9999 .and. xdbval(i).lt.500.0001) then
                udbval(i) = 'no_data'
                xdbval(i) = 9999999.
            end if
        end do
    end if

    go to 999

990 continue
    qend = .true.
    go to 999

995 continue
    qerr = .true.
    go to 999

999 continue
end subroutine rdgrid