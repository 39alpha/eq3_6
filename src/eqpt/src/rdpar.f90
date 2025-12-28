subroutine rdpar(adh,adhh,adhv,aphi,bdh,bdhh,bdhv,bdot,bdoth,bdotv,cco2,dadhh,dadhv,dbdhh,dbdhv,dbdth,dbdtv,dhfe,dvfe,ipch,ipchmx,ipcv,ipcvmx,itgenf,nacdpr,narxmx,narxt,ndat0s,ndbmax,ndbptg,ndbptl,nerr,noutpt,ntprmx,ntprt,nttyo,nwarn,prehw,presg,q500fl,tdamax,tdamin,tempc,uakey,udbfmt,udbval,xdbval,xhfe,xlke,xvfe)
    !! This subroutine reads the data describing the temperature grid
    !! and data for miscellaneous parameters which are represented
    !! on this grid from the DATA0 file. Such parameters include the
    !! pressure, the Debye-Huckel parameters A(gamma,10), B(gamma),
    !! and A(phi), the B-dot parameter of Helgeson (1969), and the
    !! equilibrium constant for the "Eh" reaction. Depending on the
    !! values of the ipch and ipcv flags, pressure derivatives of
    !! various order of these parameters may also be read. If pressure
    !! corrections are available, a grid for the half-width of the
    !! pressure correction band is also read. This subroutine also
    !! reads the nominal temperature limits for the data file. These
    !! limits are written on the DATA1 and DATA1F files. This suboutine
    !! also reads the coefficients needed to calculate the acivity
    !! coefficient of CO2(aq) using the Drummond (1981) equation.
    !! The counter nerr is incremented by one for each error encountered.
    !! The counter nwarn is similarly incremented for each warning.
    !! This suboutine is called by:
    !!   EQPT/eqpt.f
    !! Principal input:
    !!   ipch   = enthalpy functions data grid flag:
    !!              -1 = no enthalpy grids
    !!               0 = enthalpy grids present
    !!               1 = grids present for the enthalpy and its first
    !!                     pressure derivative
    !!               2 = grids present for the enthalpy and its first
    !!                     and second pressure derivatives
    !!   ipcv   = volume functions data flag:
    !!              -1 = no volume grids
    !!               0 = volume grids present
    !!               1 = grids present for the volume and its first
    !!                     pressure derivative
    !!               2 = grids present for the volume and its first
    !!                     and second pressure derivatives
    !!   ndat0s = unit number of the stripped DATA0 file
    !!   narxt  = array of numbers of coefficients in the temperature
    !!              ranges
    !!   ntprt  = the number of temperature ranges on the standard
    !!              temperature grid
    !!   uakey  = string specifying the type of data file ("SEDH" or
    !!              "Pitzer") being processed
    !! Principal output:
    !!   adh    = array of A(gamma,10) values on the standard
    !!              temperature grid
    !!   adhh   = array of A(H) values on the standard temperature
    !!              grid; A(H) contains but is not equal to
    !!              dA(gamma,10)/dT
    !!   adhv   = array of A(V) values on the standard temperature
    !!              grid; A(V) contains but is not equal to
    !!              dA(gamma,10)/dP
    !!   aphi   = array of A(phi) values on the standard temperature
    !!              grid; A(phi) = A(gamma,e)/3 = 2.303A(gamma,10)/3,
    !!              so the pressure derivatives of A(phi) are obtained
    !!              from A(V) and dnA(V)/dPn.
    !!   bdh    = array of B(gamma) values on the standard temperature
    !!              grid
    !!   bdhh   = array of B(H) values on the standard temperature grid
    !!   bdhv   = array of B(V) values on the standard temperature grid
    !!   bdot   = array of B-dot values on the standard temperature
    !!              grid
    !!   bdoth  = array of B-dot(H)values on the standard temperature
    !!              grid
    !!   bdotv  = array of B-dot(V) values on the standard temperature
    !!              grid
    !!   cco2   = array of coefficients for the Drummond (1981)
    !!              equation
    !!   dadhh  = array of dnA(H)/dPn values on the standard temperature
    !!              grid
    !!   dadhv  = array of dnA(V)/dPn values on the standard temperature
    !!              grid
    !!   dbdhh  = array of dnB(H)/dPn values on the standard temperature
    !!              grid
    !!   dbdhv  = array of dnA(V)/dPn values on the standard temperature
    !!              grid
    !!   dbdth  = array of dnB-dot(H)/dPn values on the standard
    !!              temperature grid
    !!   dbdtv  = array of dnB-dot(V)/dPn values on the standard
    !!              temperature grid
    !!   dhfe   = array of pressure derivatives of the enthalpy of
    !!              reaction for the "Eh" reaction on the standard
    !!              temperature grid
    !!   dvfe   = array of pressure derivatives of the volume of
    !!              reaction for the "Eh" reaction on the standard
    !!              temperature grid
    !!   ndbptg = the number of distinct points on the standard
    !!              temperature grid.
    !!   nerr   = cumulative error counter
    !!   nwarn  = cumulative warning counter
    !!   presg  = array of standard pressures on the standard
    !!              temperature grid
    !!   prehw  = array of half-width values for the standard pressure
    !!              envelope on the standard temperature grid
    !!   tdamax = the nominal maximum temperature (C) for which the
    !!              data file is valid
    !!   tdamin = the nominal minimum temperature (C) for which the
    !!              data file is valid
    !!   tempc  = array of temperatures defining the standard
    !!              temperature grid
    !!   udbfmt = the format used when reading a line of data on the
    !!              standard temperature grid
    !!   xlke   = array of log K values for the "Eh" reaction
    !!              on the standard temperature grid
    !!   xhfe   = array of enthalpy of reaction values for the "Eh"
    !!              reaction hon the standard temperature grid
    !!   xvfe   = array of volume of reaction values for the "Eh"
    !!              reaction on the standard temperature grid
    implicit none

    ! Calling sequence variable declarations.
    integer :: ipchmx
    integer :: ipcvmx
    integer :: narxmx
    integer :: ndbmax
    integer :: ntprmx

    integer :: ndat0s
    integer :: noutpt
    integer :: nttyo

    integer :: nacdpr(ntprmx)
    integer :: narxt(ntprmx)

    integer :: ipch
    integer :: ipcv
    integer :: itgenf
    integer :: ndbptg
    integer :: ndbptl
    integer :: nerr
    integer :: ntprt
    integer :: nwarn

    logical :: q500fl

    character(len=16) :: udbval(ndbmax)

    character(len=16) :: udbfmt
    character(len=8) :: uakey

    real(kind=8) :: tdamax
    real(kind=8) :: tdamin

    real(kind=8) :: cco2(5)

    real(kind=8) :: adh(narxmx,ntprmx)
    real(kind=8) :: adhh(narxmx,ntprmx)
    real(kind=8) :: adhv(narxmx,ntprmx)
    real(kind=8) :: aphi(narxmx,ntprmx)
    real(kind=8) :: bdh(narxmx,ntprmx)
    real(kind=8) :: bdhh(narxmx,ntprmx)
    real(kind=8) :: bdhv(narxmx,ntprmx)
    real(kind=8) :: bdot(narxmx,ntprmx)
    real(kind=8) :: bdoth(narxmx,ntprmx)
    real(kind=8) :: bdotv(narxmx,ntprmx)
    real(kind=8) :: dadhh(narxmx,ntprmx,ipchmx)
    real(kind=8) :: dadhv(narxmx,ntprmx,ipcvmx)
    real(kind=8) :: dbdhh(narxmx,ntprmx,ipchmx)
    real(kind=8) :: dbdhv(narxmx,ntprmx,ipcvmx)
    real(kind=8) :: dbdth(narxmx,ntprmx,ipchmx)
    real(kind=8) :: dbdtv(narxmx,ntprmx,ipcvmx)
    real(kind=8) :: dhfe(narxmx,ntprmx,ipchmx)
    real(kind=8) :: dvfe(narxmx,ntprmx,ipcvmx)
    real(kind=8) :: prehw(narxmx,ntprmx)
    real(kind=8) :: presg(narxmx,ntprmx)
    real(kind=8) :: tempc(narxmx,ntprmx)
    real(kind=8) :: xdbval(ndbmax)
    real(kind=8) :: xhfe(narxmx,ntprmx)
    real(kind=8) :: xlke(narxmx,ntprmx)
    real(kind=8) :: xvfe(narxmx,ntprmx)

    ! Local variable declarations.
    integer :: i
    integer :: ipc
    integer :: j2
    integer :: j3
    integer :: n
    integer :: ndbloc
    integer :: ntpr

    integer :: ilnobl

    logical :: qend
    logical :: qerr
    logical :: q500nd

    character(len=24) :: ux24
    character(len=24) :: ux24lc
    character(len=80) :: ulbufa
    character(len=80) :: ulbufb
    character(len=80) :: uline
    character(len=80) :: ustr80

    ! Skip to the 'Temperature limits' line.
    ux24 = 'Temperature limits (degC'
    ux24lc = ux24
    call locase(ux24lc)
100 continue
    read (ndat0s,1000,end=110,err=995) ustr80
1000 format(a)

    call lejust(ustr80)
    call locase(ustr80)

    if (ustr80(1:24) .ne. ux24lc(1:24)) then
        go to 100
    end if

    go to 120

110 continue
    j2 = ilnobl(ux24)
    write (noutpt,1010) ux24(1:j2)
    write (nttyo,1010) ux24(1:j2)
1010 format(/' * Error - (EQPT/rdpar) Unexpectedly encountered',' end-of-file',/7x,'while searching the DATA0 file for the',' find the line containing',/7x,'"',a,'".')

    stop
120 continue

    ! Read in the parameters. First, read the nominal temperature
    ! limits.
    read (ndat0s,1002,end=990,err=995) tdamin,tdamax
1002 format(5x,2f10.4)

    ! Read the standard temperature grid. The format of this grid will
    ! be analyzed and will serve as a template for that of subsequent
    ! data presented on this grid. The original EQ3/6 format was
    ! four points per line. Up to six points per line are now possible.
    ! The format will be determined by counting the number of data
    ! points on the first line.
    !   ndbptl = maximum number of data points per line (1-6)
    !   ndbptg = number of data points in the entire grid
    !   ndbloc = number of lines containing the grid
    ndbptg = narxt(1)

    do ntpr = 2,ntprt
        ndbptg = ndbptg + narxt(ntpr) - 1
    end do

    ux24 = 'Temperatures            '
    ux24lc = ux24
    call locase(ux24lc)
    read (ndat0s,1000,end=990,err=995) ustr80
    call lejust(ustr80)
    call locase(ustr80)

    if (ustr80(1:24) .ne. ux24lc(1:24)) then
        ux24 = 'Temperature grid (degC) '
        ux24lc = ux24
        call locase(ux24lc)
    end if

    if (ustr80(1:24) .ne. ux24lc(1:24)) then
        ux24 = 'Temperature grid (C)    '
        ux24lc = ux24
        call locase(ux24lc)
    end if

    if (ustr80(1:24) .ne. ux24lc(1:24)) then
        j2 = ilnobl(ux24)
        j3 = ilnobl(ustr80)
        j3 = min(j3,24)
        write (noutpt,1020) ustr80(1:j3),ux24(1:j2)
        write (nttyo,1020) ustr80(1:j3),ux24(1:j2)
1020 format(/' * Error - (EQPT/rdpar) Found the string "',a,'"',/7x,'on the DATA0 file where the string "',a,'"',/7x,'was expected.')

        stop
    end if

    ! Analyze the grid format.
    read (ndat0s,1030,end=990,err=995) uline
1030 format(a80)

    backspace(ndat0s)

    ulbufa = uline
    ulbufb = ' '
    ndbptl = 0

    do n = 1,6
        call lejust(ulbufa)
        i = index(ulbufa,' ')

        if (i .le. 1) then
            go to 130
        end if

        ndbptl = ndbptl + 1
        ulbufb = ulbufa(i:80)
        ulbufa = ulbufb
    end do

130 continue

    ndbloc = ndbptg/ndbptl

    if (mod(ndbptg,ndbptl) .gt. 0) then
        ndbloc = ndbloc + 1
    end if

    udbfmt = '( (5x,6f10.4) )'
    write (udbfmt(7:7),'(i1)') ndbptl

    ! The format analysis is complete.
    ! Read the temperatures on the log K temperature grid.
    ! Return the data in the xdbval holding array.
    q500nd = .false.
    call rdgrid(ndat0s,ndbmax,ndbptg,ndbptl,qend,qerr,q500nd,udbval,xdbval)

    if (qend) then
        go to 990
    end if

    if (qerr) then
        go to 995
    end if

    ! read (ndat0s,udbfmt,end=990,err=995) (xdbval(i), i = 1,ndbptg)
    ! Load the data into the tempc array.
    ! Calling sequence substitutions:
    !   tempc for zdbval
    call ldbar2(nacdpr,narxmx,narxt,ndbmax,ntprmx,ntprt,xdbval,tempc)

    ! Read the standard pressure grid.
    ux24 = 'Pressures               '
    ux24lc = ux24
    call locase(ux24lc)
    read (ndat0s,1000,end=990,err=995) ustr80
    call lejust(ustr80)
    call locase(ustr80)

    if (ustr80(1:24) .ne. ux24lc(1:24)) then
        ux24 = 'Pressure grid (bars)    '
        ux24lc = ux24
        call locase(ux24lc)
    end if

    if (ustr80(1:24) .ne. ux24lc(1:24)) then
        j2 = ilnobl(ux24)
        j3 = ilnobl(ustr80)
        j3 = min(j3,24)
        write (noutpt,1020) ustr80(1:j3),ux24(1:j2)
        write (nttyo,1020) ustr80(1:j3),ux24(1:j2)
        stop
    end if

    ! Read the pressures on the log K temperature grid.
    ! Return the data in the xdbval holding array.
    q500nd = .false.
    call rdgrid(ndat0s,ndbmax,ndbptg,ndbptl,qend,qerr,q500nd,udbval,xdbval)

    if (qend) then
        go to 990
    end if

    if (qerr) then
        go to 995
    end if

    ! Load the data into the presg array.
    ! Calling sequence substitutions:
    !   presg for zdbval
    call ldbar2(nacdpr,narxmx,narxt,ndbmax,ntprmx,ntprt,xdbval,presg)

    if (ipcv .ge. 0) then
        ! Read the standard pressure envelope half-width.
        ux24 = 'Pressure envelope half-w'
        ux24lc = ux24
        call locase(ux24lc)
        read (ndat0s,1000,end=990,err=995) ustr80
        call lejust(ustr80)
        call locase(ustr80)

        if (ustr80(1:24) .ne. ux24lc(1:24)) then
            j2 = ilnobl(ux24)
            j3 = ilnobl(ustr80)
            j3 = min(j3,24)
            write (noutpt,1020) ustr80(1:j3),ux24(1:j2)
            write (nttyo,1020) ustr80(1:j3),ux24(1:j2)
            stop
        end if

        ! Read the pressure half-widths on the log K temperature grid.
        ! Return the data in the xdbval holding array.
        q500nd = .false.
        call rdgrid(ndat0s,ndbmax,ndbptg,ndbptl,qend,qerr,q500nd,udbval,xdbval)

        if (qend) then
            go to 990
        end if

        if (qerr) then
            go to 995
        end if

        ! Load the data into the prehw array.
        ! Calling sequence substitutions:
        !   prehw for zdbval
        call ldbar2(nacdpr,narxmx,narxt,ndbmax,ntprmx,ntprt,xdbval,prehw)
    end if

    if (uakey(1:8) .eq. 'SEDH    ') then
        ! Read the grid for the Debye-Huckel A(gamma,10) parameter and
        ! the grids for related parameters.
        ux24 = 'Debye Huckel A (ADH)    '
        ux24lc = ux24
        call locase(ux24lc)
        read (ndat0s,1000,end=990,err=995) ustr80
        call lejust(ustr80)
        call locase(ustr80)

        if (ustr80(1:24) .ne. ux24lc(1:24)) then
            ux24 = 'Debye-Huckel A_gamma (kg'
            ux24lc = ux24
            call locase(ux24lc)
        end if

        if (ustr80(1:24) .ne. ux24lc(1:24)) then
            ux24 = 'Debye-Huckel A_gamma    '
            ux24lc = ux24
            call locase(ux24lc)
        end if

        if (ustr80(1:24) .ne. ux24lc(1:24)) then
            j2 = ilnobl(ux24)
            j3 = ilnobl(ustr80)
            j3 = min(j3,24)
            write (noutpt,1020) ustr80(1:j3),ux24(1:j2)
            write (nttyo,1020) ustr80(1:j3),ux24(1:j2)
            stop
        end if

        ! Read the Debye-Huckel A(gamma,10) values on the log K
        ! temperature grid. Return the data in the xdbval holding array.
        q500nd = q500fl
        call rdgrid(ndat0s,ndbmax,ndbptg,ndbptl,qend,qerr,q500nd,udbval,xdbval)

        if (qend) then
            go to 990
        end if

        if (qerr) then
            go to 995
        end if

        ! Load the data into the adh array.
        ! Calling sequence substitutions:
        !   adh for zdbval
        call ldbar2(nacdpr,narxmx,narxt,ndbmax,ntprmx,ntprt,xdbval,adh)

        if (ipch .ge. 0) then
            ! Read the grid for the Debye-Huckel A(H) parameter.
            ux24 = 'Debye-Huckel A_H (kcal*k'
            ux24lc = ux24
            call locase(ux24lc)
            read (ndat0s,1000,end=990,err=995) ustr80
            call lejust(ustr80)
            call locase(ustr80)

            if (ustr80(1:24) .ne. ux24lc(1:24)) then
                ux24 = 'Debye-Huckel A_H        '
                ux24lc = ux24
                call locase(ux24lc)
            end if

            if (ustr80(1:24) .ne. ux24lc(1:24)) then
                j2 = ilnobl(ux24)
                j3 = ilnobl(ustr80)
                j3 = min(j3,24)
                write (noutpt,1020) ustr80(1:j3),ux24(1:j2)
                write (nttyo,1020) ustr80(1:j3),ux24(1:j2)
                stop
            end if

            ! Read the Debye-Huckel A(H) values on the log K temperature
            ! grid. Return the data in the xdbval holding array.
            q500nd = q500fl
            call rdgrid(ndat0s,ndbmax,ndbptg,ndbptl,qend,qerr,q500nd,udbval,xdbval)

            if (qend) then
                go to 990
            end if

            if (qerr) then
                go to 995
            end if

            ! Load the data into the adhh array.
            ! Calling sequence substitutions:
            !   adhh for zdbval
            call ldbar2(nacdpr,narxmx,narxt,ndbmax,ntprmx,ntprt,xdbval,adhh)

            do ipc = 1,ipch
                read (ndat0s,1000,end=990,err=995) ustr80

                ! Read the Debye-Huckel A(H) pressure derivatives on the
                ! log K temperature grid. Return the data in the xdbval
                ! holding array.
                q500nd = q500fl
                call rdgrid(ndat0s,ndbmax,ndbptg,ndbptl,qend,qerr,q500nd,udbval,xdbval)

                if (qend) then
                    go to 990
                end if

                if (qerr) then
                    go to 995
                end if

                ! Load the data into the dadhh array.
                ! Calling sequence substitutions:
                !   dadhh for zdbval
                !   ipchmx for ipcmax
                call ldbar3(ipc,ipchmx,nacdpr,narxmx,narxt,ndbmax,ntprmx,ntprt,xdbval,dadhh)
            end do
        end if

        if (ipcv .ge. 0) then
            ! Read the grid for the Debye-Huckel A(V) parameter.
            ux24 = 'Debye-Huckel A_V (cm**3*'
            ux24lc = ux24
            call locase(ux24lc)
            read (ndat0s,1000,end=990,err=995) ustr80
            call lejust(ustr80)
            call locase(ustr80)

            if (ustr80(1:24) .ne. ux24lc(1:24)) then
                ux24 = 'Debye-Huckel A_V        '
                ux24lc = ux24
                call locase(ux24lc)
            end if

            if (ustr80(1:24) .ne. ux24lc(1:24)) then
                j2 = ilnobl(ux24)
                j3 = ilnobl(ustr80)
                j3 = min(j3,24)
                write (noutpt,1020) ustr80(1:j3),ux24(1:j2)
                write (nttyo,1020) ustr80(1:j3),ux24(1:j2)
                stop
            end if

            ! Read the Debye-Huckel A(V) values on the log K temperature
            ! grid. Return the data in the xdbval holding array.
            q500nd = q500fl
            call rdgrid(ndat0s,ndbmax,ndbptg,ndbptl,qend,qerr,q500nd,udbval,xdbval)

            if (qend) then
                go to 990
            end if

            if (qerr) then
                go to 995
            end if

            ! Load the data into the adhv array.
            ! Calling sequence substitutions:
            !   adhv for zdbval
            call ldbar2(nacdpr,narxmx,narxt,ndbmax,ntprmx,ntprt,xdbval,adhv)

            do ipc = 1,ipcv
                read (ndat0s,1000,end=990,err=995) ustr80

                ! Read the Debye-Huckel A(V) pressure derivatives on the
                ! log K temperature grid. Return the data in the xdbval
                ! holding array.
                q500nd = q500fl
                call rdgrid(ndat0s,ndbmax,ndbptg,ndbptl,qend,qerr,q500nd,udbval,xdbval)

                if (qend) then
                    go to 990
                end if

                if (qerr) then
                    go to 995
                end if

                ! Load the data into the dadhv array.
                ! Calling sequence substitutions:
                !   dadhv for zdbval
                !   ipcvmx for ipcmax
                call ldbar3(ipc,ipcvmx,nacdpr,narxmx,narxt,ndbmax,ntprmx,ntprt,xdbval,dadhv)
            end do
        end if

        ! Read the grid for the Debye-Huckel B(gamma) parameter and
        ux24 = 'Debye Huckel B (BDH)    '
        ux24lc = ux24
        call locase(ux24lc)
        read (ndat0s,1000,end=990,err=995) ustr80
        call lejust(ustr80)
        call locase(ustr80)

        if (ustr80(1:24) .ne. ux24lc(1:24)) then
            ux24 = 'Debye-Huckel B_gamma (kg'
            ux24lc = ux24
            call locase(ux24lc)
        end if

        if (ustr80(1:24) .ne. ux24lc(1:24)) then
            ux24 = 'Debye-Huckel B_gamma    '
            ux24lc = ux24
            call locase(ux24lc)
        end if

        if (ustr80(1:24) .ne. ux24lc(1:24)) then
            j2 = ilnobl(ux24)
            j3 = ilnobl(ustr80)
            j3 = min(j3,24)
            write (noutpt,1020) ustr80(1:j3),ux24(1:j2)
            write (nttyo,1020) ustr80(1:j3),ux24(1:j2)
            stop
        end if

        ! Read the Debye-Huckel B(gamma) values on the log K temperature
        ! grid. Return the data in the xdbval holding array.
        q500nd = q500fl
        call rdgrid(ndat0s,ndbmax,ndbptg,ndbptl,qend,qerr,q500nd,udbval,xdbval)

        if (qend) then
            go to 990
        end if

        if (qerr) then
            go to 995
        end if

        ! Load the data into the bdh array.
        ! Calling sequence substitutions:
        !   bdh for zdbval
        call ldbar2(nacdpr,narxmx,narxt,ndbmax,ntprmx,ntprt,xdbval,bdh)

        if (ipch .ge. 0) then
            ! Read the grid for the Debye-Huckel B(H) parameter.
            ux24 = 'Debye-Huckel B_H (cal*kg'
            ux24lc = ux24
            call locase(ux24lc)
            read (ndat0s,1000,end=990,err=995) ustr80
            call lejust(ustr80)
            call locase(ustr80)

            if (ustr80(1:24) .ne. ux24lc(1:24)) then
                ux24 = 'Debye-Huckel B_H        '
                ux24lc = ux24
                call locase(ux24lc)
            end if

            if (ustr80(1:24) .ne. ux24lc(1:24)) then
                j2 = ilnobl(ux24)
                j3 = ilnobl(ustr80)
                j3 = min(j3,24)
                write (noutpt,1020) ustr80(1:j3),ux24(1:j2)
                write (nttyo,1020) ustr80(1:j3),ux24(1:j2)
                stop
            end if

            ! Read the Debye-Huckel B(H) values on the log K temperature
            ! grid. Return the data in the xdbval holding array.
            q500nd = q500fl
            call rdgrid(ndat0s,ndbmax,ndbptg,ndbptl,qend,qerr,q500nd,udbval,xdbval)

            if (qend) then
                go to 990
            end if

            if (qerr) then
                go to 995
            end if

            ! Load the data into the bdhh array.
            ! Calling sequence substitutions:
            !   bdhh for zdbval
            call ldbar2(nacdpr,narxmx,narxt,ndbmax,ntprmx,ntprt,xdbval,bdhh)

            do ipc = 1,ipch
                read (ndat0s,1000,end=990,err=995) ustr80

                ! Read the Debye-Huckel B(H) pressure derivatives on the
                ! log K temperature grid. Return the data in the xdbval
                ! holding array.
                q500nd = q500fl
                call rdgrid(ndat0s,ndbmax,ndbptg,ndbptl,qend,qerr,q500nd,udbval,xdbval)

                if (qend) then
                    go to 990
                end if

                if (qerr) then
                    go to 995
                end if

                ! Load the data into the dbdhh array.
                ! Calling sequence substitutions:
                !   dbdhh for zdbval
                !   ipchmx for ipcmax
                call ldbar3(ipc,ipchmx,nacdpr,narxmx,narxt,ndbmax,ntprmx,ntprt,xdbval,dbdhh)
            end do
        end if

        if (ipcv .ge. 0) then
            ! Read the grid for the Debye-Huckel B(V) parameter.
            ux24 = 'Debye-Huckel B_V (cm**2*'
            ux24lc = ux24
            call locase(ux24lc)
            read (ndat0s,1000,end=990,err=995) ustr80
            call lejust(ustr80)
            call locase(ustr80)

            if (ustr80(1:24) .ne. ux24lc(1:24)) then
                ux24 = 'Debye-Huckel B_V        '
                ux24lc = ux24
                call locase(ux24lc)
            end if

            if (ustr80(1:24) .ne. ux24lc(1:24)) then
                j2 = ilnobl(ux24)
                j3 = ilnobl(ustr80)
                j3 = min(j3,24)
                write (noutpt,1020) ustr80(1:j3),ux24(1:j2)
                write (nttyo,1020) ustr80(1:j3),ux24(1:j2)
                stop
            end if

            ! Read the Debye-Huckel B(V) values on the log K temperature
            ! grid. Return the data in the xdbval holding array.
            q500nd = q500fl
            call rdgrid(ndat0s,ndbmax,ndbptg,ndbptl,qend,qerr,q500nd,udbval,xdbval)

            if (qend) then
                go to 990
            end if

            if (qerr) then
                go to 995
            end if

            ! Load the data into the bdhv array.
            ! Calling sequence substitutions:
            !   bdhv for zdbval
            call ldbar2(nacdpr,narxmx,narxt,ndbmax,ntprmx,ntprt,xdbval,bdhv)

            do ipc = 1,ipcv
                read (ndat0s,1000,end=990,err=995) ustr80

                ! Read the Debye-Huckel B(V) pressure derivatives on the
                ! log K temperature grid. Return the data in the xdbval
                ! holding array.
                q500nd = q500fl
                call rdgrid(ndat0s,ndbmax,ndbptg,ndbptl,qend,qerr,q500nd,udbval,xdbval)

                if (qend) then
                    go to 990
                end if

                if (qerr) then
                    go to 995
                end if

                ! Load the data into the dbdhv array.
                ! Calling sequence substitutions:
                !   dbdhv for zdbval
                !   ipcvmx for ipcmax
                call ldbar3(ipc,ipcvmx,nacdpr,narxmx,narxt,ndbmax,ntprmx,ntprt,xdbval,dbdhv)
            end do
        end if

        ! Read the grid for the B-dot parameter and the grids for
        ! related parameters.
        ux24 = 'Bdot                    '
        ux24lc = ux24
        call locase(ux24lc)
        read (ndat0s,1000,end=990,err=995) ustr80
        call lejust(ustr80)
        call locase(ustr80)

        if (ustr80(1:24) .ne. ux24lc(1:24)) then
            ux24 = 'B-dot                   '
            ux24lc = ux24
            call locase(ux24lc)
        end if

        if (ustr80(1:24) .ne. ux24lc(1:24)) then
            j2 = ilnobl(ux24)
            j3 = ilnobl(ustr80)
            j3 = min(j3,24)
            write (noutpt,1020) ustr80(1:j3),ux24(1:j2)
            write (nttyo,1020) ustr80(1:j3),ux24(1:j2)
            stop
        end if

        ! Read the B-dot values on the log K temperature grid.
        ! Return the data in the xdbval holding array.
        q500nd = q500fl
        call rdgrid(ndat0s,ndbmax,ndbptg,ndbptl,qend,qerr,q500nd,udbval,xdbval)

        if (qend) then
            go to 990
        end if

        if (qerr) then
            go to 995
        end if

        ! Load the data into the bdot array.
        ! Calling sequence substitutions:
        !   bdot for zdbval
        call ldbar2(nacdpr,narxmx,narxt,ndbmax,ntprmx,ntprt,xdbval,bdot)

        if (ipch .ge. 0) then
            ! Read the grid for the B-dot(H) parameter.
            ux24 = 'B-dot_H                 '
            ux24lc = ux24
            call locase(ux24lc)
            read (ndat0s,1000,end=990,err=995) ustr80
            call lejust(ustr80)
            call locase(ustr80)

            if (ustr80(1:24) .ne. ux24lc(1:24)) then
                j2 = ilnobl(ux24)
                j3 = ilnobl(ustr80)
                j3 = min(j3,24)
                write (noutpt,1020) ustr80(1:j3),ux24(1:j2)
                write (nttyo,1020) ustr80(1:j3),ux24(1:j2)
                stop
            end if

            ! Read the B-dot(H) values on the log K temperature grid
            ! Return the data in the xdbval holding array.
            q500nd = q500fl
            call rdgrid(ndat0s,ndbmax,ndbptg,ndbptl,qend,qerr,q500nd,udbval,xdbval)

            if (qend) then
                go to 990
            end if

            if (qerr) then
                go to 995
            end if

            ! Load the data into the bdoth array.
            ! Calling sequence substitutions:
            !   bdoth for zdbval
            call ldbar2(nacdpr,narxmx,narxt,ndbmax,ntprmx,ntprt,xdbval,bdoth)

            do ipc = 1,ipch
                read (ndat0s,1000,end=990,err=995) ustr80

                ! Read the B-dot(H) pressure derivatives on the
                ! log K temperature grid. Return the data in the xdbval
                ! holding array.
                q500nd = q500fl
                call rdgrid(ndat0s,ndbmax,ndbptg,ndbptl,qend,qerr,q500nd,udbval,xdbval)

                if (qend) then
                    go to 990
                end if

                if (qerr) then
                    go to 995
                end if

                ! Load the data into the dbdth array.
                ! Calling sequence substitutions:
                !   dbdth for zdbval
                !   ipchmx for ipcmax
                call ldbar3(ipc,ipchmx,nacdpr,narxmx,narxt,ndbmax,ntprmx,ntprt,xdbval,dbdth)
            end do
        end if

        if (ipcv .ge. 0) then
            ! Read the grid for the B-dot(V) parameter.
            ux24 = 'B-dot_V                 '
            ux24lc = ux24
            call locase(ux24lc)
            read (ndat0s,1000,end=990,err=995) ustr80
            call lejust(ustr80)
            call locase(ustr80)

            if (ustr80(1:24) .ne. ux24lc(1:24)) then
                j2 = ilnobl(ux24)
                j3 = ilnobl(ustr80)
                j3 = min(j3,24)
                write (noutpt,1020) ustr80(1:j3),ux24(1:j2)
                write (nttyo,1020) ustr80(1:j3),ux24(1:j2)
                stop
            end if

            ! Read the B-dot(V) values on the log K temperature
            ! grid. Return the data in the xdbval holding array.
            q500nd = q500fl
            call rdgrid(ndat0s,ndbmax,ndbptg,ndbptl,qend,qerr,q500nd,udbval,xdbval)

            if (qend) then
                go to 990
            end if

            if (qerr) then
                go to 995
            end if

            ! Load the data into the bdotv array.
            ! Calling sequence substitutions:
            !   bdotv for zdbval
            call ldbar2(nacdpr,narxmx,narxt,ndbmax,ntprmx,ntprt,xdbval,bdotv)

            do ipc = 1,ipcv
                read (ndat0s,1000,end=990,err=995) ustr80

                ! Read the B-dot(V) pressure derivatives on the
                ! log K temperature grid. Return the data in the xdbval
                ! holding array.
                q500nd = q500fl
                call rdgrid(ndat0s,ndbmax,ndbptg,ndbptl,qend,qerr,q500nd,udbval,xdbval)

                if (qend) then
                    go to 990
                end if

                if (qerr) then
                    go to 995
                end if

                ! Load the data into the dbdtv array.
                ! Calling sequence substitutions:
                !   dbdtv for zdbval
                !   ipcvmx for ipcmax
                call ldbar3(ipc,ipcvmx,nacdpr,narxmx,narxt,ndbmax,ntprmx,ntprt,xdbval,dbdtv)
            end do
        end if

        ! Read the coefficients for computing the activity coefficient
        ! of CO2(aq) from the Drummond (1981) equation.
        ux24 = 'Cco2   (coefficients for'
        ux24lc = ux24
        call locase(ux24lc)
        read (ndat0s,1000,end=990,err=995) ustr80
        call lejust(ustr80)
        call locase(ustr80)

        if (ustr80(1:24) .ne. ux24lc(1:24)) then
            ux24 = 'Cco2                    '
            ux24lc = ux24
            call locase(ux24lc)
        end if

        if (ustr80(1:24) .ne. ux24lc(1:24)) then
            j2 = ilnobl(ux24)
            j3 = ilnobl(ustr80)
            j3 = min(j3,24)
            write (noutpt,1020) ustr80(1:j3),ux24(1:j2)
            write (nttyo,1020) ustr80(1:j3),ux24(1:j2)
            stop
        end if

        read (ndat0s,1060,end=990,err=995) (cco2(n), n = 1,5)
1060 format(5x,f10.4,11x,f12.7,/10x,f5.1,11x,f12.4,/5x,f10.6)
    end if

    if (uakey(1:8) .eq. 'Pitzer  ') then
        ! Read the grid for the Debye-Huckel A(phi) parameter and
        ! the grids for related parameters.
        ux24 = 'Debye Huckel Aphi       '
        ux24lc = ux24
        call locase(ux24lc)
        read (ndat0s,1000,end=990,err=995) ustr80
        call lejust(ustr80)
        call locase(ustr80)

        if (ustr80(1:24) .ne. ux24lc(1:24)) then
            ux24 = 'Debye_Huckel A_phi      '
            ux24lc = ux24
            call locase(ux24lc)
        end if

        if (ustr80(1:24) .ne. ux24lc(1:24)) then
            j2 = ilnobl(ux24)
            j3 = ilnobl(ustr80)
            j3 = min(j3,24)
            write (noutpt,1020) ustr80(1:j3),ux24(1:j2)
            write (nttyo,1020) ustr80(1:j3),ux24(1:j2)
            stop
        end if

        ! Read the Debye-Huckel A(phi) values on the log K temperature
        ! grid. Return the data in the xdbval holding array.
        q500nd = q500fl
        call rdgrid(ndat0s,ndbmax,ndbptg,ndbptl,qend,qerr,q500nd,udbval,xdbval)

        if (qend) then
            go to 990
        end if

        if (qerr) then
            go to 995
        end if

        ! Load the data into the aphi array.
        ! Calling sequence substitutions:
        !   aphi for zdbval
        call ldbar2(nacdpr,narxmx,narxt,ndbmax,ntprmx,ntprt,xdbval,aphi)

        if (ipch .ge. 0) then
            ! Read the grid for the Debye-Huckel A(H) parameter.
            ux24 = 'Debye-Huckel A_H (kcal*k'
            ux24lc = ux24
            call locase(ux24lc)
            read (ndat0s,1000,end=990,err=995) ustr80
            call lejust(ustr80)
            call locase(ustr80)

            if (ustr80(1:24) .ne. ux24lc(1:24)) then
                ux24 = 'Debye-Huckel A_H        '
                ux24lc = ux24
                call locase(ux24lc)
            end if

            if (ustr80(1:24) .ne. ux24lc(1:24)) then
                j2 = ilnobl(ux24)
                j3 = ilnobl(ustr80)
                j3 = min(j3,24)
                write (noutpt,1020) ustr80(1:j3),ux24(1:j2)
                write (nttyo,1020) ustr80(1:j3),ux24(1:j2)
                stop
            end if

            ! Read the Debye-Huckel A(H) values on the log K temperature
            ! grid. Return the data in the xdbval holding array.
            q500nd = q500fl
            call rdgrid(ndat0s,ndbmax,ndbptg,ndbptl,qend,qerr,q500nd,udbval,xdbval)

            if (qend) then
                go to 990
            end if

            if (qerr) then
                go to 995
            end if

            ! Load the data into the adhh array.
            ! Calling sequence substitutions:
            !   adhh for zdbval
            call ldbar2(nacdpr,narxmx,narxt,ndbmax,ntprmx,ntprt,xdbval,adhh)

            do ipc = 1,ipch
                read (ndat0s,1000,end=990,err=995) ustr80

                ! Read the Debye-Huckel A(H) pressure derivatives on the
                ! log K temperature grid. Return the data in the xdbval
                ! holding array.
                q500nd = q500fl
                call rdgrid(ndat0s,ndbmax,ndbptg,ndbptl,qend,qerr,q500nd,udbval,xdbval)

                if (qend) then
                    go to 990
                end if

                if (qerr) then
                    go to 995
                end if

                ! Load the data into the dadhh array.
                ! Calling sequence substitutions:
                !   dadhh for zdbval
                !   ipchmx for ipcmax
                call ldbar3(ipc,ipchmx,nacdpr,narxmx,narxt,ndbmax,ntprmx,ntprt,xdbval,dadhh)
            end do
        end if

        if (ipcv .ge. 0) then
            ! Read the grid for the Debye-Huckel A(V) parameter.
            ux24 = 'Debye-Huckel A_V (cm**3*'
            ux24lc = ux24
            call locase(ux24lc)
            read (ndat0s,1000,end=990,err=995) ustr80
            call lejust(ustr80)
            call locase(ustr80)

            if (ustr80(1:24) .ne. ux24lc(1:24)) then
                ux24 = 'Debye-Huckel A_V        '
                ux24lc = ux24
                call locase(ux24lc)
            end if

            if (ustr80(1:24) .ne. ux24lc(1:24)) then
                j2 = ilnobl(ux24)
                j3 = ilnobl(ustr80)
                j3 = min(j3,24)
                write (noutpt,1020) ustr80(1:j3),ux24(1:j2)
                write (nttyo,1020) ustr80(1:j3),ux24(1:j2)
                stop
            end if

            ! Read the Debye-Huckel A(V) values on the log K temperature
            ! grid. Return the data in the xdbval holding array.
            q500nd = q500fl
            call rdgrid(ndat0s,ndbmax,ndbptg,ndbptl,qend,qerr,q500nd,udbval,xdbval)

            if (qend) then
                go to 990
            end if

            if (qerr) then
                go to 995
            end if

            ! Load the data into the adhv array.
            ! Calling sequence substitutions:
            !   adhv for zdbval
            call ldbar2(nacdpr,narxmx,narxt,ndbmax,ntprmx,ntprt,xdbval,adhv)

            do ipc = 1,ipcv
                read (ndat0s,1000,end=990,err=995) ustr80

                ! Read the Debye-Huckel A(V) pressure derivatives on the
                ! log K temperature grid. Return the data in the xdbval
                ! holding array.
                q500nd = q500fl
                call rdgrid(ndat0s,ndbmax,ndbptg,ndbptl,qend,qerr,q500nd,udbval,xdbval)

                if (qend) then
                    go to 990
                end if

                if (qerr) then
                    go to 995
                end if

                ! Load the data into the dadhv array.
                ! Calling sequence substitutions:
                !   dadhv for zdbval
                !   ipcvmx for ipcmax
                call ldbar3(ipc,ipcvmx,nacdpr,narxmx,narxt,ndbmax,ntprmx,ntprt,xdbval,dadhv)
            end do
        end if
    end if

    ! Read the grid for the log K of the "Eh" reaction and the grids
    ! for related parameters.
    ux24 = 'log k for eh reaction   '
    ux24lc = ux24
    call locase(ux24lc)
    read (ndat0s,1000,end=990,err=995) ustr80
    call lejust(ustr80)
    call locase(ustr80)

    if (ustr80(1:24) .ne. ux24lc(1:24)) then
        ux24 = 'Eh reaction: logKr (2H2O'
        ux24lc = ux24
        call locase(ux24lc)
    end if

    if (ustr80(1:24) .ne. ux24lc(1:24)) then
        ux24 = 'Eh reaction: logKr      '
        ux24lc = ux24
        call locase(ux24lc)
    end if

    if (ustr80(1:24) .ne. ux24lc(1:24)) then
        j2 = ilnobl(ux24)
        j3 = ilnobl(ustr80)
        j3 = min(j3,24)
        write (noutpt,1020) ustr80(1:j3),ux24(1:j2)
        write (nttyo,1020) ustr80(1:j3),ux24(1:j2)
        stop
    end if

    ! Read the "Eh" reaction log K values on the log K temperature
    ! grid. Return the data in the xdbval holding array.
    q500nd = q500fl
    call rdgrid(ndat0s,ndbmax,ndbptg,ndbptl,qend,qerr,q500nd,udbval,xdbval)

    if (qend) then
        go to 990
    end if

    if (qerr) then
        go to 995
    end if

    ! Load the data into the xlke array.
    ! Calling sequence substitutions:
    !   xlke for zdbval
    call ldbar2(nacdpr,narxmx,narxt,ndbmax,ntprmx,ntprt,xdbval,xlke)

    if (ipch .ge. 0) then
        ! Read the grid for the enthalpy of reaction of the "Eh"
        ! reaction.
        ux24 = 'Eh reaction: delh0r (kca'
        ux24lc = ux24
        call locase(ux24lc)
        read (ndat0s,1000,end=990,err=995) ustr80
        call lejust(ustr80)
        call locase(ustr80)

        if (ustr80(1:24) .ne. ux24lc(1:24)) then
            ux24 = 'Eh reaction: delH0r     '
            ux24lc = ux24
            call locase(ux24lc)
        end if

        if (ustr80(1:24) .ne. ux24lc(1:24)) then
            j2 = ilnobl(ux24)
            j3 = ilnobl(ustr80)
            j3 = min(j3,24)
            write (noutpt,1020) ustr80(1:j3),ux24(1:j2)
            write (nttyo,1020) ustr80(1:j3),ux24(1:j2)
            stop
        end if

        ! Read the "Eh" reaction enthalpy function values on the
        ! log K temperature grid. Return the data in the xdbval
        ! holding array.
        q500nd = q500fl
        call rdgrid(ndat0s,ndbmax,ndbptg,ndbptl,qend,qerr,q500nd,udbval,xdbval)

        if (qend) then
            go to 990
        end if

        if (qerr) then
            go to 995
        end if

        ! Load the data into the xhfe array.
        ! Calling sequence substitutions:
        !   xhfe for zdbval
        call ldbar2(nacdpr,narxmx,narxt,ndbmax,ntprmx,ntprt,xdbval,xhfe)

        do ipc = 1,ipch
            read (ndat0s,1000,end=990,err=995) ustr80

            ! Read the "Eh" reaction enthalpy function derivatives on
            ! the log K temperature grid. Return the data in the xdbval
            ! holding array.
            q500nd = q500fl
            call rdgrid(ndat0s,ndbmax,ndbptg,ndbptl,qend,qerr,q500nd,udbval,xdbval)

            if (qend) then
                go to 990
            end if

            if (qerr) then
                go to 995
            end if

            ! Load the data into the dhfe array.
            ! Calling sequence substitutions:
            !   dhfe for zdbval
            !   ipchmx for ipcmax
            call ldbar3(ipc,ipchmx,nacdpr,narxmx,narxt,ndbmax,ntprmx,ntprt,xdbval,dhfe)
        end do
    end if

    if (ipcv .ge. 0) then
        ! Read the grid for the volume of reaction of the "Eh"
        ! reaction.
        ux24 = 'Eh reaction: delv0r (cm*'
        ux24lc = ux24
        call locase(ux24lc)
        read (ndat0s,1000,end=990,err=995) ustr80
        call lejust(ustr80)
        call locase(ustr80)

        if (ustr80(1:24) .ne. ux24lc(1:24)) then
            ux24 = 'Eh reaction: delV0r     '
            ux24lc = ux24
            call locase(ux24lc)
        end if

        if (ustr80(1:24) .ne. ux24lc(1:24)) then
            j2 = ilnobl(ux24)
            j3 = ilnobl(ustr80)
            j3 = min(j3,24)
            write (noutpt,1020) ustr80(1:j3),ux24(1:j2)
            write (nttyo,1020) ustr80(1:j3),ux24(1:j2)
            stop
        end if

        ! Read the "Eh" reaction volume function values on the
        ! log K temperature grid. Return the data in the xdbval
        ! holding array.
        q500nd = q500fl
        call rdgrid(ndat0s,ndbmax,ndbptg,ndbptl,qend,qerr,q500nd,udbval,xdbval)

        if (qend) then
            go to 990
        end if

        if (qerr) then
            go to 995
        end if

        ! Load the data into the xvfe array.
        ! Calling sequence substitutions:
        !   xvfe for zdbval
        call ldbar2(nacdpr,narxmx,narxt,ndbmax,ntprmx,ntprt,xdbval,xvfe)

        do ipc = 1,ipcv
            read (ndat0s,1000,end=990,err=995) ustr80

            ! Read the "Eh" reaction volume function derivatives on
            ! the log K temperature grid. Return the data in the xdbval
            ! holding array.
            q500nd = q500fl
            call rdgrid(ndat0s,ndbmax,ndbptg,ndbptl,qend,qerr,q500nd,udbval,xdbval)

            if (qend) then
                go to 990
            end if

            if (qerr) then
                go to 995
            end if

            ! Load the data into the dvfe array.
            ! Calling sequence substitutions:
            !   dvfe for zdbval
            !   ipcvmx for ipcmax
            call ldbar3(ipc,ipcvmx,nacdpr,narxmx,narxt,ndbmax,ntprmx,ntprt,xdbval,dvfe)
        end do
    end if

    read (ndat0s,1000,end=990,err=995) ustr80

    go to 999

990 continue
    write (noutpt,2000)
    write (nttyo,2000)
2000 format(/' * Error - (EQPT/rdpar) Unexpectedly encountered',/7x,'end-of-file while reading the DATA0 file.')

    stop

995 continue
    write (noutpt,2010)
    write (nttyo,2010)
2010 format(/' * Error - (EQPT/rdpar) Encountered a read format',/7x,'error while reading the DATA0 file.')

    stop

999 continue
end subroutine rdpar