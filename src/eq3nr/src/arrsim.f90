subroutine arrsim(aamatr,acflg,actlg,bbig,cdrs,cjbasp,cnufac,conc,conclg,coval,delvec,dlogxw,eh,ehfac,eps100,gmmatr,iction,iindx1,iodb,ipivot,irdxc3,ixbasp,jcsort,jflag,jjndex,kbt,ker,khydr,kkndex,kmax,kwater,narn1,narn2,nbasp,nbt,nbti,nbtmax,nbw,ncosp,ndecsp,ndrs,ndrsmx,ndrsr,nelect,nern1,nern2,nhydr,nodbmx,no2gaq,noutpt,npass,nredox,nstmax,nttyo,omega,qawfix,rhsvec,ucospi,uspec,xbar,xbarlg,xbarw,xbrwlg,xlke,xlks,zchar,zvclg1)
    !! This subroutine computes starting estimates of species
    !! concentrations that must be evaluated simultaneously. These
    !! include cases of mean activity constraints, cases of equilibrium
    !! constraints, and cases in which the log fO2 is constrained by Eh,
    !! pe, or a redox couple.
    !! This subroutine is called by:
    !!   EQ3NR/arrset.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: kmax
    integer :: nbtmax
    integer :: ndrsmx
    integer :: nodbmx
    integer :: nstmax

    integer :: noutpt
    integer :: nttyo

    integer :: iindx1(kmax)
    integer :: iodb(nodbmx)
    integer :: iction(nbtmax)
    integer :: ipivot(kmax)
    integer :: ixbasp(nbtmax)
    integer :: jcsort(nstmax)
    integer :: jjndex(nbtmax)
    integer :: jflag(nstmax)
    integer :: kkndex(nbtmax)
    integer :: nbasp(nbtmax)
    integer :: ncosp(nbtmax)
    integer :: ndecsp(nbtmax)
    integer :: ndrs(ndrsmx)
    integer :: ndrsr(2,nstmax)

    integer :: kbt
    integer :: ker
    integer :: khydr
    integer :: kwater
    integer :: narn1
    integer :: narn2
    integer :: nbt
    integer :: nbti
    integer :: nbw
    integer :: nelect
    integer :: nern1
    integer :: nern2
    integer :: nhydr
    integer :: no2gaq
    integer :: npass
    integer :: nredox

    logical :: qawfix

    character(len=48) :: ucospi(nbtmax)
    character(len=48) :: uspec(nstmax)

    real(kind=8) :: aamatr(kmax,kmax)
    real(kind=8) :: acflg(nstmax)
    real(kind=8) :: actlg(nstmax)
    real(kind=8) :: cdrs(ndrsmx)
    real(kind=8) :: cjbasp(nbtmax)
    real(kind=8) :: cnufac(nstmax)
    real(kind=8) :: conc(nstmax)
    real(kind=8) :: conclg(nstmax)
    real(kind=8) :: coval(nbtmax)
    real(kind=8) :: delvec(kmax)
    real(kind=8) :: dlogxw(nbtmax)
    real(kind=8) :: gmmatr(kmax,kmax)
    real(kind=8) :: rhsvec(kmax)
    real(kind=8) :: xbar(nstmax)
    real(kind=8) :: xbarlg(nstmax)
    real(kind=8) :: xlks(nstmax)
    real(kind=8) :: zchar(nstmax)
    real(kind=8) :: zvclg1(kmax)

    real(kind=8) :: bbig
    real(kind=8) :: eh
    real(kind=8) :: ehfac
    real(kind=8) :: eps100
    real(kind=8) :: omega
    real(kind=8) :: xbarw
    real(kind=8) :: xbrwlg
    real(kind=8) :: xlke

    ! Local variable declarations.
    integer :: i
    integer :: ib
    integer :: ibt
    integer :: icol
    integer :: ielect
    integer :: ier
    integer :: ihydr
    integer :: io2gaq
    integer :: irdxc3
    integer :: irow
    integer :: jfl
    integer :: j2
    integer :: j3
    integer :: krow
    integer :: n
    integer :: nb
    integer :: nbh
    integer :: nb1
    integer :: nb2
    integer :: nerr
    integer :: nr1
    integer :: nr2
    integer :: ns
    integer :: nsc
    integer :: ns1
    integer :: ns2

    integer :: ilnobl
    integer :: nbasis

    logical :: qpr
    logical :: qx

    character(len=24) :: uaqeq
    character(len=24) :: ublk24
    character(len=24) :: ueh
    character(len=24) :: umolfr
    character(len=24) :: ured
    character(len=24) :: ust1
    character(len=24) :: ust2
    character(len=24) :: ust3
    character(len=8) :: uptgas

    real(kind=8) :: azns
    real(kind=8) :: azns1
    real(kind=8) :: azsum
    real(kind=8) :: cx
    real(kind=8) :: cxx
    real(kind=8) :: dzx
    real(kind=8) :: rx
    real(kind=8) :: xbrwlo
    real(kind=8) :: zp
    real(kind=8) :: zx
    real(kind=8) :: zxo

    real(kind=8) :: coefdr
    real(kind=8) :: texp

    data ublk24 /'                        '/
    data umolfr /'Mole fraction           '/
    data ueh    /'Eh                      '/
    data ured   /'Aqueous redox reaction  '/
    data uaqeq  /'Aqueous equilibrium     '/
    data uptgas /'Gas     '/

    data qpr    /.false./

    ! Save the current value of the mole fraction of water.
    xbrwlo = xbrwlg

    ! Build a structure for a matrix equation. Note that water
    ! and its mole fraction relation may or may not be included
    ! in this structure. Consequently, the kkndex array is redefined
    ! from that built by EQ3NR/dawfix.f.
    ihydr = 0
    ielect = 0
    io2gaq = 0

    ib = 0

    do krow = 1,kbt
        nb = iindx1(krow)
        kkndex(nb) = 0
        ns = nbasp(nb)
        jfl = jflag(ns)
        qx = .false.

        if (ns .eq. narn1) then
            qx = jfl.eq.0 .and. npass.gt.1 .and. bbig.le.0.1
            qx = qx .or. qawfix
        else if (ns.eq.nelect .or. ns.eq.no2gaq) then
            qx = irdxc3 .ne. 0
        else
            qx = jfl.eq.17 .or. jfl.eq.18 .or. jfl.eq.21 .or.    jfl.eq.25 .or. jfl.eq.27
        end if

        if (qx) then
            ib = ib + 1
            jjndex(ib) = nb
            kkndex(nb) = 1

            if (ns .eq. nhydr) then
                ihydr = ib
            end if

            if (ns .eq. nelect) then
                ielect = ib
            end if

            if (ns .eq. no2gaq) then
                io2gaq = ib
            end if
        end if
    end do

    ibt = ib

    ! Quit if no species concentrations are to be solved simultaneously.
    if (ibt .le. 0) then
        go to 999
    end if

    ! Set up the iction array, which supports the jflag = 17, 18, and 21
    ! options.
    do nb = 1,nbt
        iction(nb) = 0
    end do

    do irow = 1,ibt
        nb = jjndex(irow)
        ns = nbasp(nb)
        jfl = jflag(ns)
        ns1 = 0

        if (jfl.eq.17 .or. jfl.eq.18 .or. jfl.eq.21) then
            ns1 = ncosp(nb)

            ! Calling sequence substitutions:
            !   ns1 for ns
            nb1 = nbasis(nbasp,nbt,nbtmax,ns1)

            if (nb1 .eq. 0) then
                j3 = ilnobl(uspec(ns)(1:24))
                j2 = ilnobl(uspec(ns1)(1:24))
                write (noutpt,1000) uspec(ns1)(1:j2),uspec(ns)(1:j3)
                write (nttyo,1000) uspec(ns1)(1:j2),uspec(ns)(1:j3)
1000 format(/' * Error - (EQ3NR/arrsim) The specified',' counterion',/7x,a,' in the',/7x,'jflag = 17, 18, or',' 21 option  for the basis species',/7x,a," isn't in",' the basis set.')

                stop
            end if

            if (kkndex(nb1) .ge. 1) then
                do icol = 1,ibt
                    nb2 = jjndex(icol)

                    if (nb2 .eq. nb1) then
                        iction(nb) = icol
                        go to 110
                    end if
                end do
            end if
        end if

110 continue
    end do

    ! Optional print of the matrix structure.
    if (iodb(3) .ge. 3) then
        write (noutpt,1020)
1020 format(/7x,'--- Structure for Simultaneous Evaluation of',' Starting Estimates ---',/)

        do irow = 1,ibt
            nb1 = jjndex(irow)
            ns1 = nbasp(nb1)
            jfl = jflag(ns1)
            ust1 = uspec(ns1)
            ust3 = ublk24

            if (ns1 .eq. narn1) then
                ust2 = umolfr
            end if

            if (ns1.eq.nelect .or. ns1.eq.no2gaq) then
                if (ncosp(nb1) .eq. 0) then
                    if (irdxc3 .lt. 0) then
                        ust2 = ueh
                    else
                        ust2 = ured
                    end if

                    go to 120
                end if
            end if

            if (jfl .eq. 17) then
                ns = ncosp(nb1)
                ust2 = 'Log activity combination'
                ust3 = uspec(ns)
            else if (jfl .eq. 18) then
                ns = ncosp(nb1)
                ust2 = 'Mean log activity'
                ust3 = uspec(ns)
            else if (jfl .eq. 21) then
                ns = ncosp(nb1)
                ust2 = 'pHCl'
                ust3 = uspec(ns)
            else if (jfl .eq. 27) then
                ust2 = uaqeq
            else
                do nb = 1,nbti
                    nb2 = ndecsp(nb)

                    if (nb1 .eq. nb2) then
                        ust2 = ucospi(nb)(1:24)
                        ust3 = ucospi(nb)(25:48)
                        go to 120
                    end if
                end do
            end if

120 continue
            j2 = ilnobl(ust2)
            write (noutpt,1030) irow,ust1,ust2(1:j2)
1030 format(2x,i5,2x,a24,2x,a)

            j3 = ilnobl(ust3)

            if (ust3(1:24) .ne. ublk24(1:24)) then
                write (noutpt,1040) ust3(1:j3)
            end if

1040 format(37x,'Constraint species= ',a)
        end do

        write (noutpt,1050)
1050 format(//10x,'--- The zvclg1 and acflg Values ---',/)

        do krow = 1,kbt
            nb = iindx1(krow)
            ns = nbasp(nb)

            if (kkndex(nb) .le. 0) then
                write (noutpt,1060) krow,uspec(ns),zvclg1(krow),acflg(ns)
1060 format(2x,i5,2x,a24,2x,f10.4,2x,f10.4)
            end if
        end do

        write (noutpt,1070)
1070 format(/1x)
    end if

    ! Set up the cnufac array.
    do nsc = narn1,narn2
        cnufac(nsc) = 0.
    end do

    nr1 = ndrsr(1,narn1)

    if (jflag(narn1) .eq. 30) then
        cnufac(narn1) = xbarlg(narn1)/cdrs(nr1)
    end if

    do nsc = narn1 + 1,narn2
        nr1 = ndrsr(1,nsc)

        if (jflag(nsc) .eq. 30) then
            cnufac(nsc) = conc(nsc)/cdrs(nr1)
        end if
    end do

    ! Build the matrix.
    do icol = 1,ibt
        do irow = 1,ibt
            aamatr(irow,icol) = 0.
        end do
    end do

    do irow = 1,ibt
        nb = jjndex(irow)
        ns = nbasp(nb)
        jfl = jflag(ns)

        if (ns .eq. narn1) then
            ! Mole fraction of water, based on a first-order expansion:
            !   log x(w)(new) = log x(w)(old)
            !     + sum over s' in irow (d log x(w)/d log m(s'))
            !       * (log m(s')(new) - log m(s')(old))
            ! Compute the dlogxw array (d log xw/d log ms').
            call gdlgxw(cdrs,cjbasp,cnufac,conc,dlogxw,eps100,ixbasp,jcsort,jflag,narn1,narn2,nbasp,nbt,nbtmax,nbw,ndrs,ndrsmx,ndrsr,nern1,nern2,noutpt,nstmax,nttyo,omega,xbar,xbarw)

            do icol = 1,ibt
                nb1 = jjndex(icol)
                aamatr(irow,icol) = -dlogxw(nb1)
            end do

            aamatr(irow,irow) = 1.0
        else if (jfl .eq. 17) then
            ns1 = ncosp(nb)
            azns = abs(zchar(ns))
            azns1 = abs(zchar(ns1))
            aamatr(irow,irow) = azns1
            icol = iction(nb)

            if (icol .gt. 0) then
                zp = zchar(ns)*zchar(ns1)

                if (zp .lt. 0) then
                    aamatr(irow,icol) = azns
                else
                    aamatr(irow,icol) = -azns
                end if
            end if
        else if (jfl .eq. 18) then
            ns1 = ncosp(nb)
            azns = abs(zchar(ns))
            azns1 = abs(zchar(ns1))
            azsum = azns + azns1
            aamatr(irow,irow) = azns1/azsum
            icol = iction(nb)

            if (icol .gt. 0) then
                aamatr(irow,icol) = azns/azsum
            end if
        else if (jfl .eq. 21) then
            ns1 = ncosp(nb)
            azns = abs(zchar(ns))
            azns1 = abs(zchar(ns1))
            aamatr(irow,irow) = -azns1
            icol = iction(nb)

            if (icol .gt. 0) then
                zp = zchar(ns)*zchar(ns1)

                if (zp .lt. 0) then
                    aamatr(irow,icol) = -azns
                else
                    aamatr(irow,icol) = azns
                end if
            end if
        else if (jfl .eq. 27) then
            do icol = 1,ibt
                nb1 = jjndex(icol)
                ns1 = nbasp(nb1)

                ! Calling sequence substitutions:
                !   ns1 for nse
                aamatr(irow,icol) =      coefdr(cdrs,ndrs,ndrsmx,ndrsr,ns1,ns,nstmax)
            end do

            ! Calling sequence substitutions:
            !   ns for nse
            aamatr(irow,irow)    = coefdr(cdrs,ndrs,ndrsmx,ndrsr,ns,ns,nstmax)
        else if (jfl .eq. 25) then
            ns1 = ncosp(nb)

            do icol = 1,ibt
                nb2 = jjndex(icol)
                ns2 = nbasp(nb2)

                ! Calling sequence substitutions:
                !   ns2 for nse
                !   ns1 for ns
                aamatr(irow,icol)      = coefdr(cdrs,ndrs,ndrsmx,ndrsr,ns2,ns1,nstmax)
            end do
        else if (irdxc3 .lt. 0) then
            aamatr(irow,io2gaq) = 1.

            if (ihydr .gt. 0) then
                aamatr(irow,ihydr) = 4.
            end if
        else
            do icol = 1,ibt
                nb1 = jjndex(icol)
                ns1 = nbasp(nb1)

                if (kkndex(nb1).gt.0) then
                    ! Calling sequence substitutions:
                    !   ns1 for nse
                    !   nredox for ns
                    aamatr(irow,icol)        = coefdr(cdrs,ndrs,ndrsmx,ndrsr,ns1,nredox,nstmax)
                end if
            end do
        end if
    end do

    if (iodb(3) .ge. 3) then
        write (noutpt,1100)
1100 format(/10x,'--- Matrix ---',/)

        do irow = 1,ibt
            write (noutpt,1110) (aamatr(irow,icol), icol = 1,ibt)
1110 format(2x,10(f7.2,2x))
        end do

        write (noutpt,1070)
    end if

    ! Build the right-hand-side vector.
    do irow = 1,ibt
        nb = jjndex(irow)
        ns = nbasp(nb)
        jfl = jflag(ns)

        if (ns .eq. narn1) then
            ! Mole fraction of water, based on a first-order expansion:
            !   log x(w)(new) = log x(w)(old)
            !     + sum over s' in irow (d log x(w)/d log m(s'))
            !       * (log m(s')(new) - log m(s')(old))
            rx = xbrwlo

            do krow = 1,kbt
                if (krow .ne. kwater) then
                    nb1 = iindx1(krow)

                    if (kkndex(nb1) .ge. 1) then
                        rx = rx - dlogxw(nb1)*zvclg1(krow)
                    end if
                end if
            end do

            rhsvec(irow) = rx
        else if (jfl .eq. 17) then
            ns1 = ncosp(nb)
            azns = abs(zchar(ns))
            azns1 = abs(zchar(ns1))
            rx = coval(nb) - azns1*acflg(ns)
            zp = zchar(ns)*zchar(ns1)

            if (zp .lt. 0.) then
                rx = rx - azns*acflg(ns1)
            else
                rx = rx + azns*acflg(ns1)
            end if

            icol = iction(nb)

            if (icol .eq. 0) then
                if (zp .lt. 0.) then
                    rx = rx - azns*conclg(ns1)
                else
                    rx = rx + azns*conclg(ns1)
                end if
            end if

            rhsvec(irow) = rx
        else if (jfl .eq. 18) then
            ns1 = ncosp(nb)
            azns = abs(zchar(ns))
            azns1 = abs(zchar(ns1))
            azsum = azns + azns1
            rx = coval(nb) - ( azns1/azsum )*acflg(ns)    - ( azns/azsum )*acflg(ns1)
            icol = iction(nb)

            if (icol .eq. 0) then
                rx = rx - ( azns/azsum )*conclg(ns1)
            end if

            rhsvec(irow) = rx
        else if (jfl .eq. 21) then
            ns1 = ncosp(nb)
            azns = abs(zchar(ns))
            azns1 = abs(zchar(ns1))
            rx = coval(nb) + azns1*acflg(ns)
            zp = zchar(ns)*zchar(ns1)

            if (zp .lt. 0.) then
                rx = rx + azns*acflg(ns1)
            else
                rx = rx - azns*acflg(ns1)
            end if

            icol = iction(nb)

            if (icol .eq. 0) then
                if (zp .lt. 0.) then
                    rx = rx + azns*conclg(ns1)
                else
                    rx = rx - azns*conclg(ns1)
                end if
            end if

            rhsvec(irow) = rx
        else if (jfl .eq. 27) then
            if (iodb(3) .ge. 3) then
                ! Calling sequence substitutions:
                !   noutpt for nf
                call prreac(cdrs,ndrs,ndrsmx,ndrsr,noutpt,ns,nstmax,uspec)
            end if

            rx = xlks(ns)
            nr1 = ndrsr(1,ns)
            nr2 = ndrsr(2,ns)

            do n = nr1,nr2
                ns2 = ndrs(n)
                cx = cdrs(n)

                if (ns2 .eq. narn1) then
                    rx = rx - cx*actlg(narn1)
                else
                    ! Calling sequence substitutions:
                    !   ns2 for ns
                    nb2 = nbasis(nbasp,nbt,nbtmax,ns2)

                    if (kkndex(nb2) .le. 0) then
                        rx = rx - cx*actlg(ns2)
                    else
                        rx = rx - cx*acflg(ns2)
                    end if
                end if
            end do

            rhsvec(irow) = rx
        else if (jfl .eq. 25) then
            ns1 = ncosp(nb)

            if (iodb(3) .ge. 3) then
                ! Calling sequence substitutions:
                !   noutpt for nf
                !   ns1 for ns
                call prreac(cdrs,ndrs,ndrsmx,ndrsr,noutpt,ns1,nstmax,uspec)
            end if

            rx = xlks(ns1)
            nr1 = ndrsr(1,ns1)
            nr2 = ndrsr(2,ns1)

            do n = nr1,nr2
                ns2 = ndrs(n)
                cx = cdrs(n)

                if (ns2 .eq. ns1) then
                    if (uspec(ns1)(25:32) .eq. uptgas(1:8)) then
                        rx = rx - cx*coval(nb)
                    else
                        rx = rx - cx*actlg(ns1)
                    end if
                else if (ns2 .eq. narn1) then
                    rx = rx - cx*actlg(narn1)
                else
                    ! Calling sequence substitutions:
                    !   ns2 for ns
                    nb2 = nbasis(nbasp,nbt,nbtmax,ns2)

                    if (kkndex(nb2) .le. 0) then
                        rx = rx - cx*actlg(ns2)
                    else
                        rx = rx - cx*acflg(ns2)
                    end if
                end if
            end do

            rhsvec(irow) = rx
        else if (irdxc3 .lt. 0) then
            rx = 4.*eh/ehfac
            rx = rx + xlke + 2.*actlg(narn1)
            rx = rx - 4.*acflg(nhydr)
            nbh = iindx1(khydr)

            if (kkndex(nbh) .le. 0) then
                rx = rx - 4.*conclg(nhydr)
            end if

            rhsvec(irow) = rx
        else if (irdxc3 .gt. 0) then
            rx = xlks(nredox)
            nr1 = ndrsr(1,nredox)
            nr2 = ndrsr(2,nredox)

            do n = nr1,nr2
                ns2 = ndrs(n)
                cx = cdrs(n)

                if (ns2 .eq. narn1) then
                    rx = rx - cx*actlg(narn1)
                else
                    ! Calling sequence substitutions:
                    !   ns2 for ns
                    nb2 = nbasis(nbasp,nbt,nbtmax,ns2)

                    if (kkndex(nb2) .le. 0) then
                        rx = rx - cx*actlg(ns2)
                    else
                        rx = rx - cx*acflg(ns2)
                    end if
                end if
            end do

            rhsvec(irow) = rx
        end if
    end do

    if (iodb(3) .ge. 3) then
        write (noutpt,1200)
1200 format(/10x,'--- Right-Hand-Side Vector ---',/)

        do irow = 1,ibt
            nb = jjndex(irow)
            ns = nbasp(nb)
            write (noutpt,1210) irow,uspec(ns),rhsvec(irow)
1210 format(1x,i5,2x,a24,2x,g12.5)
        end do

        write (noutpt,1070)
    end if

    ! Solve the matrix equation.
    ! Calling sequence substitutions:
    !   ibt for kdim
    call msolvr(aamatr,delvec,gmmatr,ier,ipivot,ibt,kmax,noutpt,nttyo,qpr,rhsvec)

    if (ier .gt. 0) then
        write (noutpt,1220)
        write (nttyo,1220)
1220 format(/' * Error - (EQ3NR/arrsim) The matrix equation',/7x,'required by the speciation model appears to be',/7x,'singular. Check the solubility and other constraints',/7x,'on the input file for a violation of the mineralogic',/7x,'phase rule.')

        stop
    end if

    ! Apply change limits and range limits. Then load the results.
    irow = 0

    do krow = 1,kbt
        nb = iindx1(krow)

        if (kkndex(nb) .ge. 1) then
            irow = irow + 1
            zx = delvec(irow)
            zxo = zvclg1(krow)

            if (zxo .gt. -99999.) then
                dzx = zx - zxo

                if (ns1 .eq. narn1) then
                    if (dzx .gt. 0.2) then
                        dzx = 0.2
                    end if

                    if (dzx .lt. -0.2) then
                        dzx = -0.2
                    end if
                else
                    if (dzx .gt. 2.0) then
                        dzx = 2.0
                    end if

                    if (dzx .lt. -2.0) then
                        dzx = -2.0
                    end if
                end if

                zx = zxo + dzx
            end if

            if (ns1.eq.narn1 .and. zx.gt.0.) then
                zx = 0.
            end if

            if (zx .gt. 2.0) then
                zx = 2.0
            end if

            if (zx .lt. -100.) then
                zx = -100.
            end if

            zvclg1(krow) = zx
            ns = nbasp(nb)
            conclg(ns) = zx
        end if
    end do

    if (iodb(3) .ge. 3) then
        write (noutpt,1050)

        do krow = 1,kbt
            nb = iindx1(krow)
            ns = nbasp(nb)

            if (kkndex(nb) .ge. 1) then
                write (noutpt,1060) krow,uspec(ns),zvclg1(krow),acflg(ns)
            end if
        end do

        write (noutpt,1070)
    end if

    ! Check possible effects of large concentrations resulting from
    ! equilibrium constraints.
    i = ibt

    if (ielect.gt.0 .or. io2gaq.gt.0) then
        i = i - 1
    end if

    if (i .le. 0) then
        go to 999
    end if

    nerr = 0

    do krow = 1,kbt
        nb = iindx1(krow)
        ns = nbasp(nb)

        if (kkndex(nb).ge.1 .and. ns.ne.nelect .and. ns.ne.no2gaq) then
            cx = conclg(ns)
            cxx = texp(cx)

            if (cxx .gt. 20.) then
                ker = 1
                j2 = ilnobl(uspec(ns)(1:24))
                write (noutpt,1230) uspec(ns)(1:j2),cxx
                write (nttyo,1230) uspec(ns)(1:j2),cxx
1230 format(/' * Note - (EQ3NR/arrsim) The species ',a,/7x,'appears to have a required molality near ',g12.5)

                if (npass.gt.1 .and. cxx.gt.100.) then
                    nerr = nerr + 1
                end if
            end if
        end if
    end do

    if (nerr .gt. 0) then
        write (noutpt,1240)
        write (nttyo,1240)
1240 format(/' * Error - (EQ3NR/arrsim) Reconsider your choice ',/7x,'of input constraints.')

        ker = 2
        go to 999
    end if

    if (ker .eq. 1) then
        write (noutpt,1250)
        write (nttyo,1250)
1250 format(/' * Warning - (EQ3NR/arrsim) This run may crash',/7x,'because of poor choice of input constraints.')
    end if

999 continue
end subroutine arrsim