subroutine wr3w7(cspb,fep,iktmax,iodb,iopg,iopr,iopt,itermx,jflagb,jxmod,kxmod,ncompb,newin,nodbmx,nopgmx,noprmx,noptmx,nsq,nsqmax,ntitl,ntitmx,nxmdmx,nxmod,nxtb,nxtmax,rho,tempc,tdspkg,tdspl,tolbt,toldl,tolsat,ubasis,uebal,umemb,uphas1,uphas2,uredox,usolb,uspecb,utitl,xbarb,uxmd24,xlkmod)
    !! This subroutine writes the EQ3NR input file in compact ("W")
    !! format for versions 7.0-7.2. It thus encompasses two version
    !! levels.
    !! This subroutine is called by:
    !!   XCON3/xcon3.f
    implicit none

    ! Calling sequence variable declarations.
    integer :: iktmax
    integer :: nodbmx
    integer :: nopgmx
    integer :: noprmx
    integer :: noptmx
    integer :: nsqmax
    integer :: ntitmx
    integer :: nxmdmx
    integer :: nxtmax

    integer :: iodb(nodbmx)
    integer :: iopg(nopgmx)
    integer :: iopr(noprmx)
    integer :: iopt(noptmx)
    integer :: jflagb(nsqmax)
    integer :: jxmod(nxmdmx)
    integer :: kxmod(nxmdmx)
    integer :: ncompb(nxtmax)

    integer :: itermx
    integer :: newin
    integer :: nsq
    integer :: ntitl
    integer :: nxmod
    integer :: nxtb

    character(len=80) :: utitl(ntitmx)
    character(len=24) :: ubasis(nsqmax)
    character(len=24) :: umemb(iktmax,nxtmax)
    character(len=24) :: uphas1(nsqmax)
    character(len=24) :: uphas2(nsqmax)
    character(len=24) :: usolb(nxtmax)
    character(len=24) :: uspecb(nsqmax)
    character(len=24) :: uxmd24(nxmdmx)
    character(len=24) :: uebal
    character(len=24) :: uredox

    real(kind=8) :: cspb(nsqmax)
    real(kind=8) :: xbarb(iktmax,nxtmax)
    real(kind=8) :: xlkmod(nxmdmx)

    real(kind=8) :: fep
    real(kind=8) :: rho
    real(kind=8) :: tempc
    real(kind=8) :: tdspkg
    real(kind=8) :: tdspl
    real(kind=8) :: tolbt
    real(kind=8) :: toldl
    real(kind=8) :: tolsat

    ! Local variable declarations.
    integer :: i
    integer :: ikb
    integer :: iktb
    integer :: j2
    integer :: n
    integer :: ns
    integer :: nxb
    integer :: ilnobl

    ! Write the new input file.
    ! Title.
    do 110 n = 1,ntitl
        j2 = ilnobl(utitl(n))
        write (newin,1000) utitl(n)(1:j2)
1000 format(a)

110 continue

        if (ntitl .lt. ntitmx) then
            write (newin,1010)
        end if

1010 format('endit.')

        ! Temperature.
        write (newin,1040) tempc
1040 format(5x,'tempc= ',1pe12.5)

        ! Density, total dissolved salts (per kg), and total dissolved
        ! salts (per liter).
        write (newin,1050) rho,tdspkg,tdspl
1050 format (7x,'rho= ',1pe12.5,4x,'tdspkg= ',1pe12.5,5x,'tdspl= ',1pe12.5)

        ! Redox parameter and name of redox species defining a controlling
        ! couple, if any.
        j2 = ilnobl(uredox)
        write (newin,1060) fep,uredox(1:j2)
1060 format (7x,'fep= ',1pe12.5,4x,'uredox= ',a)

        ! Convergence tolerances (tolbt and toldl) and saturation
        ! tolerance (tolsat).
        write (newin,1070) tolbt,toldl,tolsat
1070 format (5x,'tolbt= ',1pe12.5,5x,'toldl= ',1pe12.5,4x,'tolsat= ',1pe12.5)

        ! Maximum number of iterations.
        write (newin,1080) itermx
1080 format (4x,'itermx= ',i2)

        write (newin,1090)
1090 format('*',15x,'1    2    3    4    5    6    7    8    9   10')

        ! Iopt option switches.
        ! Note: iopt(1) = iopt1, etc.
        write (newin,1100) (iopt(i), i = 1,10)
1100 format(2x,'iopt1-10= ',10i5)

        ! Iopg option switches.
        ! Note: iopg(1) = iopg1, etc.
        write (newin,1110) (iopg(i), i = 1,10)
1110 format(2x,'iopg1-10= ',10i5)

        ! Iopr option switches.
        ! Note: iopr(1) = iopr1, etc.
        write (newin,1120) (iopr(i), i = 1,20)
1120 format(2x,'iopr1-10= ',10i5,/1x,'iopr11-20= ',10i5)

        ! Iodb option switches.
        ! Note: iodb(1) = iodb1, etc.
        write (newin,1130) (iodb(i), i = 1,10)
1130 format(2x,'iodb1-10= ',10i5)

        ! Species for electrical balancing.
        j2 = ilnobl(uebal)
        write (newin,1150) uebal(1:j2)
1150 format(5x,'uebal= ',a)

        ! Nxmod options.
        write (newin,1200) nxmod
1200 format(5x,'nxmod= ',i2)

        if (nxmod .gt. 0) then
            do 420 n = 1,nxmod
                j2 = ilnobl(uxmd24(n))
                write (newin,1220) uxmd24(n)(1:j2),jxmod(n),kxmod(n),xlkmod(n)
1220 format(3x,'species= ',a,/6x,'type= ',i2,14x,'option= ',i2,14x,'xlkmod= ',1pe12.5)

420 continue
            end if

            ! Basis species and associated constraints.
            do 440 ns = 1,nsq
                j2 = ilnobl(uspecb(ns))
                write (newin,1240) uspecb(ns)(1:j2)
1240 format('data file master species= ',a)

                j2 = ilnobl(ubasis(ns))
                write (newin,1250) ubasis(ns)(1:j2)
1250 format(3x,'switch with species= ',a)

                write (newin,1260) jflagb(ns),cspb(ns)
1260 format(3x,'jflag= ',i2,3x,'csp= ',1pg12.5)

                if (jflagb(ns).ge.17 .and. jflagb(ns).le.18) then
                    j2 = ilnobl(uphas1(ns))
                    write (newin,1270) uphas1(ns)(1:j2)
1270 format(2x,'co-ion= ',a)
                else if (jflagb(ns) .eq. 19) then
                    j2 = ilnobl(uphas1(ns))
                    write (newin,1272) uphas1(ns)(1:j2)
1272 format(1x,'mineral= ',a)
                else if (jflagb(ns) .eq. 20) then
                    j2 = ilnobl(uphas2(ns))
                    write (newin,1275) uphas1(ns),uphas2(ns)(1:j2)
1275 format(1x,'sol-sol= ',a24,2x,'end-mem= ',a)
                else if (jflagb(ns) .eq. 21) then
                    j2 = ilnobl(uphas1(ns))
                    write (newin,1277) uphas1(ns)(1:j2)
1277 format(5x,'gas= ',a)
                end if

440 continue

                write (newin,1280)
1280 format('endit.')

                ! Mole fractions of solid solutions.
                if (iopt(4) .ge. 2) then
                    write (newin,1300)
1300 format('*   Solid solution compositions')

                    do 460 nxb = 1,nxtb
                        j2 = ilnobl(usolb(nxb))
                        write (newin,1310) usolb(nxb)(1:j2)
1310 format(3x,a)

                        iktb = ncompb(nxb)

                        do 450 ikb = 1,iktb
                            write (newin,1320) umemb(ikb,nxb),xbarb(ikb,nxb)
1320 format(6x,a24,3x,f10.4)

450 continue

                            write (newin,1330)
1330 format(6x,'endit.')

460 continue

                            write (newin,1340)
1340 format(3x,'endit.')
                        end if

999 continue
                    end subroutine wr3w7