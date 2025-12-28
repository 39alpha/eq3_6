subroutine wrpar(aamatr,adh,adhh,adhv,aphi,apr,avgrid,bdh,bdhh,bdhv,bdot,bdoth,bdotv,cco2,cof,dadhh,dadhv,dbdhh,dbdhv,dbdth,dbdtv,dhfe,dvfe,eps100,gmmatr,ipch,ipchmx,ipcv,ipcvmx,ipivot,narxmx,narxt,ndata1,ndat1f,noutpt,ntprmx,ntprt,nttyo,presg,prehw,tdamax,tdamin,tempc,tempcs,tmpcmx,uakey,xhfe,xlke,xvfe,xvec,yvec)
    !! This subroutine takes the data read from the DATA0 file by the
    !! EQPT/rdpar.f, processes it by converting data on a temperature
    !! grid to the equivalent set of coefficients of interpolating
    !! polynomials, and writes this data on the DATA1 and DATA1F files.
    !! This subroutine is called by:
    !!   EQPT/eqpt.f
    !! Principal input:
    !!   adh    = array of A(gamma,10) values on the standard temperature
    !!              grid
    !!   aphi   = array of A(phi) values on the standard temperature grid
    !!   bdh    = array of B(gamma) values on the standard temperature
    !!              grid
    !!   bdot   = array of B-dot values on the standard temperature
    !!              grid
    !!   cco2   = array of coefficients for the Drummond (1981)
    !!              equation
    !!   ndata1 = unit number of the DATA1 file
    !!   ndat1f = unit number of the DATA1F file
    !!   presg  = array of pressures on the standard temperature grid
    !!   tempc  = array of temperatures defining the standard
    !!              temperature grid
    !!   uakey  = string specifying the type of data file ("SEDH" or
    !!              "Pitzer") being processed
    !!   xlke   = array of log K values for the special "Eh" reaction
    !!              on the standard temperature grid
    !! Principal output:
    !!   apr    = array of polynomial coefficients (combined, for all
    !!              temperature ranges)
    implicit none

    ! Calling sequence variable declarations.
    integer :: ipchmx
    integer :: ipcvmx
    integer :: narxmx
    integer :: ntprmx

    integer :: ipivot(narxmx)
    integer :: narxt(ntprmx)

    integer :: ndata1
    integer :: ndat1f
    integer :: noutpt
    integer :: nttyo

    integer :: ipch
    integer :: ipcv
    integer :: ntprt

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
    real(kind=8) :: xhfe(narxmx,ntprmx)
    real(kind=8) :: xlke(narxmx,ntprmx)
    real(kind=8) :: xvfe(narxmx,ntprmx)

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
    integer :: i
    integer :: ipc
    integer :: j2
    integer :: j3
    integer :: n
    integer :: nt
    integer :: ntpr

    integer :: ilnobl

    character(len=10) :: ux10
    character(len=72) :: uterm
    character(len=72) :: utermc
    character(len=56) :: ux56
    character(len=24) :: ux24
    character(len=24) :: ux24a
    character(len=24) :: ux24b

    real(kind=8) :: adtxmn
    real(kind=8) :: adtx11
    real(kind=8) :: p01_8
    real(kind=8) :: txmn
    real(kind=8) :: tx11

    uterm(1:48) = '+-----------------------------------------------'
    uterm(49:72) = '------------------------'
    utermc = uterm
    utermc(1:1) = '*'

    ! Write the temperature limits on the output files.
    ux56 = 'Data file maximum and minimum temperatures (C)'
    write (ndata1) ux56
    j2 = ilnobl(ux56)
    write (ndat1f,1100) ux56(1:j2)
1100 format(a)

    write (ndata1) tdamin,tdamax
    write (ndat1f,1110) tdamin,tdamax
1110 format(2f10.3)

    write (noutpt,1120) tdamin
    write (nttyo,1120) tdamin
1120 format(//' The minimum temperature is ',f10.3,'C')

    write (noutpt,1130) tdamax
    write (nttyo,1130) tdamax
1130 format(' The maximum temperature is ',f10.3,'C')

    ux56 = 'Maximum temperature (C) by range'
    write (ndata1) ux56
    j2 = ilnobl(ux56)
    write (ndat1f,1140) ux56(1:j2)
1140 format(a)

    do ntpr = 1,ntprt
        n = narxt(ntpr)
        write (ndata1) tempc(n,ntpr)
        write (ndat1f,1150) tempc(n,ntpr)
1150 format(f10.3)
    end do

    write (ndat1f,1010) utermc(1:72)

    write (noutpt,1160)
    write (nttyo,1160)
1160 format(/' The maximum temperatures (C) by range are:')

    do ntpr = 1,ntprt
        n = narxt(ntpr)
        write (noutpt,1170) tempc(n,ntpr)
        write (nttyo,1170) tempc(n,ntpr)
1170 format(3x,f10.3)
    end do

    write (noutpt,1000)
1000 format(/)

    tx11 = tempc(1,1)
    txmn = tdamin

    ! Note: on some data files, 0.01C is used as an approximation
    ! to 0C in order to avoid some equation-of-state difficulties.
    ! In such cases, consider 0.01C to be the same as 0C for the
    ! purpose of the following test.
    ux10 = '0.01      '
    read (ux10,'(f10.3)') p01_8
    adtx11 = abs(tx11 - p01_8)
    adtxmn = abs(txmn - p01_8)

    if (adtx11 .le. 1.e-6) then
        tx11 = 0.
    end if

    if (adtxmn .le. 1.e-6) then
        txmn = 0.
    end if

    if ((tx11 - tdamin) .gt. 0.001) then
        write (ux24a,"(f10.3)") tdamin
        call lejust(ux24a)
        j2 = ilnobl(ux24a)
        write (ux24b,"(f10.3)") tempc(1,1)
        call lejust(ux24b)
        j3 = ilnobl(ux24b)
        write (noutpt,1200) ux24a(1:j2),ux24b(1:j3)
        write (nttyo,1200) ux24a(1:j2),ux24b(1:j3)
1200 format(/' * Warning - (EQPT/rdpar) The minimum temperature for',' this data file',/7x,'is less than the lowest temperature in',' the first temperature',/7x,'range. If EQ3NR and EQ6 are',' required to make any calculations',/7x,'between ',a,'C',' and ',a,'C, they will do so by extrapolating',/7x,'the',' data from the first range.')
    end if

    n = narxt(ntprt)

    if ((tdamax - tempc(n,ntprt)) .gt. 0.001) then
        write (ux24a,"(f10.3)") tempc(n,ntprt)
        call lejust(ux24a)
        j2 = ilnobl(ux24a)
        write (ux24b,"(f10.3)") tdamax
        call lejust(ux24b)
        j3 = ilnobl(ux24b)
        write (noutpt,1210) ux24a(1:j2),ux24b(1:j3)
        write (nttyo,1210) ux24a(1:j2),ux24b(1:j3)
1210 format(/' * Warning - (EQPT/rdpar) The maximum temperature for',' this data file',/7x,'is greater than the highest temperature',' in the last temperature',/7x,'range. If EQ3NR and EQ6 are',' required to make any calculations',/7x,'between ',a,'C',' and ',a,'C, they will do so by extrapolating',/7x,'the',' data from the last range.')
    end if

    ! Process and write the data for the standard pressure grid.
    ! Calling sequence substitutions:
    !   presg for avgrid
    call intrp(aamatr,apr,presg,cof,eps100,gmmatr,ipivot,narxmx,narxt,noutpt,ntprmx,ntprt,nttyo,tempc,tempcs,tmpcmx,xvec,yvec)

    ux24 = 'presg'
    j2 = ilnobl(ux24)
    write (ndata1) ux24
    write (ndat1f,1010) ux24(1:j2)
1010 format(a)

    write (noutpt,1020) ux24(1:j2)
1020 format(7x,a)

    do ntpr = 1,ntprt
        nt = narxt(ntpr)
        write (ndata1) (apr(i,ntpr), i = 1,nt)
        write (ndat1f,1030) (apr(i,ntpr), i = 1,nt)
1030 format( 5(1pe16.9) )
    end do

    write (ndat1f,1010) utermc(1:72)

    if (ipcv .ge. 0) then
        ! Process and write the data for the half-width of the standard
        ! pressure envelope.
        ! Calling sequence substitutions:
        !   prehw for avgrid
        call intrp(aamatr,apr,prehw,cof,eps100,gmmatr,ipivot,narxmx,narxt,noutpt,ntprmx,ntprt,nttyo,tempc,tempcs,tmpcmx,xvec,yvec)

        ux24 = 'prehw'
        j2 = ilnobl(ux24)
        write (ndata1) ux24
        write (ndat1f,1010) ux24(1:j2)
        write (noutpt,1020) ux24(1:j2)

        do ntpr = 1,ntprt
            nt = narxt(ntpr)
            write (ndata1) (apr(i,ntpr), i = 1,nt)
            write (ndat1f,1030) (apr(i,ntpr), i = 1,nt)
        end do

        write (ndat1f,1010) utermc(1:72)
    end if

    ! Process and write the next parameters only if processing a
    ! "SEDH" type data file.
    if (uakey(1:8) .eq. 'SEDH    ') then
        ! Process and write the A(gamma,10) data.
        ! Calling sequence substitutions:
        !   adh for avgrid
        call intrp(aamatr,apr,adh,cof,eps100,gmmatr,ipivot,narxmx,narxt,noutpt,ntprmx,ntprt,nttyo,tempc,tempcs,tmpcmx,xvec,yvec)

        ux24 = 'adh'
        j2 = ilnobl(ux24)
        write (ndata1) ux24
        write (ndat1f,1010) ux24(1:j2)
        write (noutpt,1020) ux24(1:j2)

        do ntpr = 1,ntprt
            nt = narxt(ntpr)
            write (ndata1) (apr(i,ntpr), i = 1,nt)
            write (ndat1f,1030) (apr(i,ntpr), i = 1,nt)
        end do

        write (ndat1f,1010) utermc(1:72)

        if (ipch .ge. 0) then
            ! Process and write the A(H) data.
            ! Calling sequence substitutions:
            !   adhh for avgrid
            call intrp(aamatr,apr,adhh,cof,eps100,gmmatr,ipivot,narxmx,narxt,noutpt,ntprmx,ntprt,nttyo,tempc,tempcs,tmpcmx,xvec,yvec)

            ux24 = 'adhh'
            j2 = ilnobl(ux24)
            write (ndata1) ux24
            write (ndat1f,1010) ux24(1:j2)
            write (noutpt,1020) ux24(1:j2)

            do ntpr = 1,ntprt
                nt = narxt(ntpr)
                write (ndata1) (apr(i,ntpr), i = 1,nt)
                write (ndat1f,1030) (apr(i,ntpr), i = 1,nt)
            end do

            write (ndat1f,1010) utermc(1:72)

            ! Process and write the dnA(H)/dPn data.
            do ipc = 1,ipch
                do ntpr = 1,ntprt
                    do n = 1,narxt(ntpr)
                        avgrid(n,ntpr) = dadhh(n,ntpr,ipc)
                    end do
                end do

                call intrp(aamatr,apr,avgrid,cof,eps100,gmmatr,ipivot,narxmx,narxt,noutpt,ntprmx,ntprt,nttyo,tempc,tempcs,tmpcmx,xvec,yvec)

                ux24 = 'dadhh( )'
                write (ux24(7:7),'(i1)') ipc
                j2 = ilnobl(ux24)
                write (ndata1) ux24
                write (ndat1f,1010) ux24(1:j2)
                write (noutpt,1020) ux24(1:j2)

                do ntpr = 1,ntprt
                    nt = narxt(ntpr)
                    write (ndata1) (apr(i,ntpr), i = 1,nt)
                    write (ndat1f,1030) (apr(i,ntpr), i = 1,nt)
                end do

                write (ndat1f,1010) utermc(1:72)
            end do
        end if

        if (ipcv .ge. 0) then
            ! Process and write the A(V) data.
            ! Calling sequence substitutions:
            !   adhv for avgrid
            call intrp(aamatr,apr,adhv,cof,eps100,gmmatr,ipivot,narxmx,narxt,noutpt,ntprmx,ntprt,nttyo,tempc,tempcs,tmpcmx,xvec,yvec)

            ux24 = 'adhv'
            j2 = ilnobl(ux24)
            write (ndata1) ux24
            write (ndat1f,1010) ux24(1:j2)
            write (noutpt,1020) ux24(1:j2)

            do ntpr = 1,ntprt
                nt = narxt(ntpr)
                write (ndata1) (apr(i,ntpr), i = 1,nt)
                write (ndat1f,1030) (apr(i,ntpr), i = 1,nt)
            end do

            write (ndat1f,1010) utermc(1:72)

            ! Process and write the dnA(V)/dPn data.
            do ipc = 1,ipcv
                do ntpr = 1,ntprt
                    do n = 1,narxt(ntpr)
                        avgrid(n,ntpr) = dadhv(n,ntpr,ipc)
                    end do
                end do

                call intrp(aamatr,apr,avgrid,cof,eps100,gmmatr,ipivot,narxmx,narxt,noutpt,ntprmx,ntprt,nttyo,tempc,tempcs,tmpcmx,xvec,yvec)

                ux24 = 'dadhv( )'
                write (ux24(7:7),'(i1)') ipc
                j2 = ilnobl(ux24)
                write (ndata1) ux24
                write (ndat1f,1010) ux24(1:j2)
                write (noutpt,1020) ux24(1:j2)

                do ntpr = 1,ntprt
                    nt = narxt(ntpr)
                    write (ndata1) (apr(i,ntpr), i = 1,nt)
                    write (ndat1f,1030) (apr(i,ntpr), i = 1,nt)
                end do

                write (ndat1f,1010) utermc(1:72)
            end do
        end if

        ! Process and write the B(gamma) data.
        ! Calling sequence substitutions:
        !   bdh for avgrid
        call intrp(aamatr,apr,bdh,cof,eps100,gmmatr,ipivot,narxmx,narxt,noutpt,ntprmx,ntprt,nttyo,tempc,tempcs,tmpcmx,xvec,yvec)

        ux24 = 'bdh'
        j2 = ilnobl(ux24)
        write (ndata1) ux24
        write (ndat1f,1010) ux24(1:j2)
        write (noutpt,1020) ux24(1:j2)

        do ntpr = 1,ntprt
            nt = narxt(ntpr)
            write (ndata1) (apr(i,ntpr), i = 1,nt)
            write (ndat1f,1030) (apr(i,ntpr), i = 1,nt)
        end do

        write (ndat1f,1010) utermc(1:72)

        if (ipch .ge. 0) then
            ! Process and write the B(H) data.
            ! Calling sequence substitutions:
            !   bdhh for avgrid
            call intrp(aamatr,apr,bdhh,cof,eps100,gmmatr,ipivot,narxmx,narxt,noutpt,ntprmx,ntprt,nttyo,tempc,tempcs,tmpcmx,xvec,yvec)

            ux24 = 'bdhh'
            j2 = ilnobl(ux24)
            write (ndata1) ux24
            write (ndat1f,1010) ux24(1:j2)
            write (noutpt,1020) ux24(1:j2)

            do ntpr = 1,ntprt
                nt = narxt(ntpr)
                write (ndata1) (apr(i,ntpr), i = 1,nt)
                write (ndat1f,1030) (apr(i,ntpr), i = 1,nt)
            end do

            write (ndat1f,1010) utermc(1:72)

            ! Process and write the dnB(H)/dPn data.
            do ipc = 1,ipch
                do ntpr = 1,ntprt
                    do n = 1,narxt(ntpr)
                        avgrid(n,ntpr) = dbdhh(n,ntpr,ipc)
                    end do
                end do

                call intrp(aamatr,apr,avgrid,cof,eps100,gmmatr,ipivot,narxmx,narxt,noutpt,ntprmx,ntprt,nttyo,tempc,tempcs,tmpcmx,xvec,yvec)

                ux24 = 'dbdhh( )'
                write (ux24(7:7),'(i1)') ipc
                j2 = ilnobl(ux24)
                write (ndata1) ux24
                write (ndat1f,1010) ux24(1:j2)
                write (noutpt,1020) ux24(1:j2)

                do ntpr = 1,ntprt
                    nt = narxt(ntpr)
                    write (ndata1) (apr(i,ntpr), i = 1,nt)
                    write (ndat1f,1030) (apr(i,ntpr), i = 1,nt)
                end do

                write (ndat1f,1010) utermc(1:72)
            end do
        end if

        if (ipcv .ge. 0) then
            ! Process and write the B(V) data.
            ! Calling sequence substitutions:
            !   bdhv for avgrid
            call intrp(aamatr,apr,bdhv,cof,eps100,gmmatr,ipivot,narxmx,narxt,noutpt,ntprmx,ntprt,nttyo,tempc,tempcs,tmpcmx,xvec,yvec)

            ux24 = 'bdhv'
            j2 = ilnobl(ux24)
            write (ndata1) ux24
            write (ndat1f,1010) ux24(1:j2)
            write (noutpt,1020) ux24(1:j2)

            do ntpr = 1,ntprt
                nt = narxt(ntpr)
                write (ndata1) (apr(i,ntpr), i = 1,nt)
                write (ndat1f,1030) (apr(i,ntpr), i = 1,nt)
            end do

            write (ndat1f,1010) utermc(1:72)

            ! Process and write the dnB(V)/dPn data.
            do ipc = 1,ipcv
                do ntpr = 1,ntprt
                    do n = 1,narxt(ntpr)
                        avgrid(n,ntpr) = dbdhv(n,ntpr,ipc)
                    end do
                end do

                call intrp(aamatr,apr,avgrid,cof,eps100,gmmatr,ipivot,narxmx,narxt,noutpt,ntprmx,ntprt,nttyo,tempc,tempcs,tmpcmx,xvec,yvec)

                ux24 = 'dbdhv( )'
                write (ux24(7:7),'(i1)') ipc
                j2 = ilnobl(ux24)
                write (ndata1) ux24
                write (ndat1f,1010) ux24(1:j2)
                write (noutpt,1020) ux24(1:j2)

                do ntpr = 1,ntprt
                    nt = narxt(ntpr)
                    write (ndata1) (apr(i,ntpr), i = 1,nt)
                    write (ndat1f,1030) (apr(i,ntpr), i = 1,nt)
                end do

                write (ndat1f,1010) utermc(1:72)
            end do
        end if

        ! Process and write the B-dot data.
        ! interpolate and write the polynomial coefficients for bdot
        ! Calling sequence substitutions:
        !   bdot for avgrid
        call intrp(aamatr,apr,bdot,cof,eps100,gmmatr,ipivot,narxmx,narxt,noutpt,ntprmx,ntprt,nttyo,tempc,tempcs,tmpcmx,xvec,yvec)

        ux24 = 'bdot'
        j2 = ilnobl(ux24)
        write (ndata1) ux24
        write (ndat1f,1010) ux24(1:j2)
        write (noutpt,1020) ux24(1:j2)

        do ntpr = 1,ntprt
            nt = narxt(ntpr)
            write (ndata1) (apr(i,ntpr), i = 1,nt)
            write (ndat1f,1030) (apr(i,ntpr), i = 1,nt)
        end do

        write (ndat1f,1010) utermc(1:72)

        if (ipch .ge. 0) then
            ! Process and write the B-dot(H) data.
            ! Calling sequence substitutions:
            !   bdoth for avgrid
            call intrp(aamatr,apr,bdoth,cof,eps100,gmmatr,ipivot,narxmx,narxt,noutpt,ntprmx,ntprt,nttyo,tempc,tempcs,tmpcmx,xvec,yvec)

            ux24 = 'bdoth'
            j2 = ilnobl(ux24)
            write (ndata1) ux24
            write (ndat1f,1010) ux24(1:j2)
            write (noutpt,1020) ux24(1:j2)

            do ntpr = 1,ntprt
                nt = narxt(ntpr)
                write (ndata1) (apr(i,ntpr), i = 1,nt)
                write (ndat1f,1030) (apr(i,ntpr), i = 1,nt)
            end do

            write (ndat1f,1010) utermc(1:72)

            ! Process and write the dnB-dot(H)/dPn data.
            do ipc = 1,ipch
                do ntpr = 1,ntprt
                    do n = 1,narxt(ntpr)
                        avgrid(n,ntpr) = dbdth(n,ntpr,ipc)
                    end do
                end do

                call intrp(aamatr,apr,avgrid,cof,eps100,gmmatr,ipivot,narxmx,narxt,noutpt,ntprmx,ntprt,nttyo,tempc,tempcs,tmpcmx,xvec,yvec)

                ux24 = 'dbdth( )'
                write (ux24(7:7),'(i1)') ipc
                j2 = ilnobl(ux24)
                write (ndata1) ux24
                write (ndat1f,1010) ux24(1:j2)
                write (noutpt,1020) ux24(1:j2)

                do ntpr = 1,ntprt
                    nt = narxt(ntpr)
                    write (ndata1) (apr(i,ntpr), i = 1,nt)
                    write (ndat1f,1030) (apr(i,ntpr), i = 1,nt)
                end do

                write (ndat1f,1010) utermc(1:72)
            end do
        end if

        if (ipcv .ge. 0) then
            ! Process and write the B-dot(V) data.
            ! Calling sequence substitutions:
            !   bdotv for avgrid
            call intrp(aamatr,apr,bdotv,cof,eps100,gmmatr,ipivot,narxmx,narxt,noutpt,ntprmx,ntprt,nttyo,tempc,tempcs,tmpcmx,xvec,yvec)

            ux24 = 'bdotv'
            j2 = ilnobl(ux24)
            write (ndata1) ux24
            write (ndat1f,1010) ux24(1:j2)
            write (noutpt,1020) ux24(1:j2)

            do ntpr = 1,ntprt
                nt = narxt(ntpr)
                write (ndata1) (apr(i,ntpr), i = 1,nt)
                write (ndat1f,1030) (apr(i,ntpr), i = 1,nt)
            end do

            write (ndat1f,1010) utermc(1:72)

            ! Process and write the dnB-dot(V)/dPn data.
            do ipc = 1,ipcv
                do ntpr = 1,ntprt
                    do n = 1,narxt(ntpr)
                        avgrid(n,ntpr) = dbdtv(n,ntpr,ipc)
                    end do
                end do

                call intrp(aamatr,apr,avgrid,cof,eps100,gmmatr,ipivot,narxmx,narxt,noutpt,ntprmx,ntprt,nttyo,tempc,tempcs,tmpcmx,xvec,yvec)

                ux24 = 'dbdtv( )'
                write (ux24(7:7),'(i1)') ipc
                j2 = ilnobl(ux24)
                write (ndata1) ux24
                write (ndat1f,1010) ux24(1:j2)
                write (noutpt,1020) ux24(1:j2)

                do ntpr = 1,ntprt
                    nt = narxt(ntpr)
                    write (ndata1) (apr(i,ntpr), i = 1,nt)
                    write (ndat1f,1030) (apr(i,ntpr), i = 1,nt)
                end do

                write (ndat1f,1010) utermc(1:72)
            end do
        end if

        ! Write the CCO2 coefficients.
        ux24 = 'cco2'
        j2 = ilnobl(ux24)
        write (ndata1) ux24
        write (ndat1f,1010) ux24(1:j2)
        write (noutpt,1020) ux24(1:j2)

        write (ndata1) (cco2(i), i = 1,5)
        write (ndat1f,1030) (cco2(i), i = 1,5)
        write (ndat1f,1010) utermc(1:72)
    end if

    ! Process and write the next parameter only if processing a
    ! "Pitzer" type data file.
    if (uakey(1:8) .eq. 'Pitzer  ') then
        ! Process and write the A(phi) data.
        ! Calling sequence substitutions:
        !   aphi for avgrid
        call intrp(aamatr,apr,aphi,cof,eps100,gmmatr,ipivot,narxmx,narxt,noutpt,ntprmx,ntprt,nttyo,tempc,tempcs,tmpcmx,xvec,yvec)

        ux24 = 'aphi'
        j2 = ilnobl(ux24)
        write (ndata1) ux24
        write (ndat1f,1010) ux24(1:j2)
        write (noutpt,1020) ux24(1:j2)

        do ntpr = 1,ntprt
            nt = narxt(ntpr)
            write (ndata1) (apr(i,ntpr), i = 1,nt)
            write (ndat1f,1030) (apr(i,ntpr), i = 1,nt)
        end do

        write (ndat1f,1010) utermc(1:72)

        if (ipch .ge. 0) then
            ! Process and write the A(H) data.
            ! Calling sequence substitutions:
            !   adhh for avgrid
            call intrp(aamatr,apr,adhh,cof,eps100,gmmatr,ipivot,narxmx,narxt,noutpt,ntprmx,ntprt,nttyo,tempc,tempcs,tmpcmx,xvec,yvec)

            ux24 = 'adhh'
            j2 = ilnobl(ux24)
            write (ndata1) ux24
            write (ndat1f,1010) ux24(1:j2)
            write (noutpt,1020) ux24(1:j2)

            do ntpr = 1,ntprt
                nt = narxt(ntpr)
                write (ndata1) (apr(i,ntpr), i = 1,nt)
                write (ndat1f,1030) (apr(i,ntpr), i = 1,nt)
            end do

            write (ndat1f,1010) utermc(1:72)

            ! Process and write the dnA(H)/dPn data.
            do ipc = 1,ipch
                do ntpr = 1,ntprt
                    do n = 1,narxt(ntpr)
                        avgrid(n,ntpr) = dadhh(n,ntpr,ipc)
                    end do
                end do

                call intrp(aamatr,apr,avgrid,cof,eps100,gmmatr,ipivot,narxmx,narxt,noutpt,ntprmx,ntprt,nttyo,tempc,tempcs,tmpcmx,xvec,yvec)

                ux24 = 'dadhh( )'
                write (ux24(7:7),'(i1)') ipc
                j2 = ilnobl(ux24)
                write (ndata1) ux24
                write (ndat1f,1010) ux24(1:j2)
                write (noutpt,1020) ux24(1:j2)

                do ntpr = 1,ntprt
                    nt = narxt(ntpr)
                    write (ndata1) (apr(i,ntpr), i = 1,nt)
                    write (ndat1f,1030) (apr(i,ntpr), i = 1,nt)
                end do

                write (ndat1f,1010) utermc(1:72)
            end do
        end if

        if (ipcv .ge. 0) then
            ! Process and write the A(V) data.
            ! Calling sequence substitutions:
            !   adhv for avgrid
            call intrp(aamatr,apr,adhv,cof,eps100,gmmatr,ipivot,narxmx,narxt,noutpt,ntprmx,ntprt,nttyo,tempc,tempcs,tmpcmx,xvec,yvec)

            ux24 = 'adhv'
            j2 = ilnobl(ux24)
            write (ndata1) ux24
            write (ndat1f,1010) ux24(1:j2)
            write (noutpt,1020) ux24(1:j2)

            do ntpr = 1,ntprt
                nt = narxt(ntpr)
                write (ndata1) (apr(i,ntpr), i = 1,nt)
                write (ndat1f,1030) (apr(i,ntpr), i = 1,nt)
            end do

            write (ndat1f,1010) utermc(1:72)

            ! Process and write the dnA(V)/dPn data.
            do ipc = 1,ipcv
                do ntpr = 1,ntprt
                    do n = 1,narxt(ntpr)
                        avgrid(n,ntpr) = dadhv(n,ntpr,ipc)
                    end do
                end do

                call intrp(aamatr,apr,avgrid,cof,eps100,gmmatr,ipivot,narxmx,narxt,noutpt,ntprmx,ntprt,nttyo,tempc,tempcs,tmpcmx,xvec,yvec)

                ux24 = 'dadhv( )'
                write (ux24(7:7),'(i1)') ipc
                j2 = ilnobl(ux24)
                write (ndata1) ux24
                write (ndat1f,1010) ux24(1:j2)
                write (noutpt,1020) ux24(1:j2)

                do ntpr = 1,ntprt
                    nt = narxt(ntpr)
                    write (ndata1) (apr(i,ntpr), i = 1,nt)
                    write (ndat1f,1030) (apr(i,ntpr), i = 1,nt)
                end do

                write (ndat1f,1010) utermc(1:72)
            end do
        end if
    end if

    ! Process and write the thermodynamic data for the "Eh" reaction.
    ! Begin with the log K data.
    ! Calling sequence substitutions:
    !   xlke for avgrid
    call intrp(aamatr,apr,xlke,cof,eps100,gmmatr,ipivot,narxmx,narxt,noutpt,ntprmx,ntprt,nttyo,tempc,tempcs,tmpcmx,xvec,yvec)

    ux24 = 'xlke'
    j2 = ilnobl(ux24)
    write (ndata1) ux24
    write (ndat1f,1010) ux24(1:j2)
    write (noutpt,1020) ux24(1:j2)

    do ntpr = 1,ntprt
        nt = narxt(ntpr)
        write (ndata1) (apr(i,ntpr), i = 1,nt)
        write (ndat1f,1030) (apr(i,ntpr), i = 1,nt)
    end do

    write (ndat1f,1010) utermc(1:72)

    if (ipch .ge. 0) then
        ! Process and write the enthalpy of reaction data.
        ! Calling sequence substitutions:
        !   xhfe for avgrid
        call intrp(aamatr,apr,xhfe,cof,eps100,gmmatr,ipivot,narxmx,narxt,noutpt,ntprmx,ntprt,nttyo,tempc,tempcs,tmpcmx,xvec,yvec)

        ux24 = 'xhfe'
        j2 = ilnobl(ux24)
        write (ndata1) ux24
        write (ndat1f,1010) ux24(1:j2)
        write (noutpt,1020) ux24(1:j2)

        do ntpr = 1,ntprt
            nt = narxt(ntpr)
            write (ndata1) (apr(i,ntpr), i = 1,nt)
            write (ndat1f,1030) (apr(i,ntpr), i = 1,nt)
        end do

        write (ndat1f,1010) utermc(1:72)

        ! Process and write the pressure derivatives of the enthalpy
        ! of reaction.
        do ipc = 1,ipch
            do ntpr = 1,ntprt
                do n = 1,narxt(ntpr)
                    avgrid(n,ntpr) = dhfe(n,ntpr,ipc)
                end do
            end do

            call intrp(aamatr,apr,avgrid,cof,eps100,gmmatr,ipivot,narxmx,narxt,noutpt,ntprmx,ntprt,nttyo,tempc,tempcs,tmpcmx,xvec,yvec)

            ux24 = 'dhfe( )'
            write (ux24(7:7),'(i1)') ipc
            j2 = ilnobl(ux24)
            write (ndata1) ux24
            write (ndat1f,1010) ux24(1:j2)
            write (noutpt,1020) ux24(1:j2)

            do ntpr = 1,ntprt
                nt = narxt(ntpr)
                write (ndata1) (apr(i,ntpr), i = 1,nt)
                write (ndat1f,1030) (apr(i,ntpr), i = 1,nt)
            end do

            write (ndat1f,1010) utermc(1:72)
        end do
    end if

    if (ipcv .ge. 0) then
        ! Process and write the volume of reaction data.
        ! Calling sequence substitutions:
        !   xvfe for avgrid
        call intrp(aamatr,apr,xvfe,cof,eps100,gmmatr,ipivot,narxmx,narxt,noutpt,ntprmx,ntprt,nttyo,tempc,tempcs,tmpcmx,xvec,yvec)

        ux24 = 'xvfe'
        j2 = ilnobl(ux24)
        write (ndata1) ux24
        write (ndat1f,1010) ux24(1:j2)
        write (noutpt,1020) ux24(1:j2)

        do ntpr = 1,ntprt
            nt = narxt(ntpr)
            write (ndata1) (apr(i,ntpr), i = 1,nt)
            write (ndat1f,1030) (apr(i,ntpr), i = 1,nt)
        end do

        write (ndat1f,1010) utermc(1:72)

        ! Process and write the pressure derivatives of the volume
        ! of reaction.
        do ipc = 1,ipcv
            do ntpr = 1,ntprt
                do n = 1,narxt(ntpr)
                    avgrid(n,ntpr) = dvfe(n,ntpr,ipc)
                end do
            end do

            call intrp(aamatr,apr,avgrid,cof,eps100,gmmatr,ipivot,narxmx,narxt,noutpt,ntprmx,ntprt,nttyo,tempc,tempcs,tmpcmx,xvec,yvec)

            ux24 = 'dvfe( )'
            write (ux24(7:7),'(i1)') ipc
            j2 = ilnobl(ux24)
            write (ndata1) ux24
            write (ndat1f,1010) ux24(1:j2)
            write (noutpt,1020) ux24(1:j2)

            do ntpr = 1,ntprt
                nt = narxt(ntpr)
                write (ndata1) (apr(i,ntpr), i = 1,nt)
                write (ndat1f,1030) (apr(i,ntpr), i = 1,nt)
            end do

            write (ndat1f,1010) utermc(1:72)
        end do
    end if
end subroutine wrpar