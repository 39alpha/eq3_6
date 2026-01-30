subroutine rc3ocf(amu,jpfcmx,ifcphi1,ifcphi2,ifnnn,ifn2n,ifpsi1,ifpsi2,ifzeta,ilcphi1,ilcphi2,ilnnn,iln2n,ilpsi1,ilpsi2,ilzeta,iodb,nmux,nmut,nmutmx,nodbmx,noutpt,nstmax,nttyo,uspec,zchar)
    !! This subroutine recalculates the original coefficients for
    !! the third-order Cphi, psi, and zeta parameters. Coefficients
    !! for parameters originally in mu form are are unchanged.
    !! After all the third-order data are calculated and otherwise
    !! processed, the Cphi data are transformed into the equivalent
    !! C data.
    !! This routine supports the default option to us C, psi,
    !! and zeta directly in the calculation of activity
    !! coefficients. This can be changed to the former treatment
    !! using the USEOLDPITZERMU option string in the input file
    !! title (this changes qhawep from .true. to .false.).
    !! This subroutine is called by:
    !!   EQ3NR/eq3nr.f
    !!   EQ6/eq6.f
    !! Principal input:
    !!   amu    = array of coefficients for computing mu interaction
    !!            parameters
    !!   nmux   = array containing the indices of the species
    !!              corresponding to data in the amu array
    !!   uspec  = array of species names
    !!   zchar  = array of charge numbers
    !! Principal output:
    !!   amu    = modified, array of coefficients for computing
    !!            Cphi, psi, zeta, and some remaining mu parameters
    implicit none

    ! Calling sequence variable declarations.
    integer :: noutpt
    integer :: nttyo

    integer :: jpfcmx
    integer :: nmutmx
    integer :: nodbmx
    integer :: nstmax

    integer :: ifcphi1
    integer :: ifcphi2
    integer :: ifnnn
    integer :: ifn2n
    integer :: ifpsi1
    integer :: ifpsi2
    integer :: ifzeta
    integer :: ilcphi1
    integer :: ilcphi2
    integer :: ilnnn
    integer :: iln2n
    integer :: ilpsi1
    integer :: ilpsi2
    integer :: ilzeta

    integer :: nmut

    integer :: nmux(3,nmutmx)
    integer :: iodb(nodbmx)

    character(len=48) :: uspec(nstmax)

    real(kind=8) :: amu(jpfcmx,nmutmx)
    real(kind=8) :: zchar(nstmax)

    ! Local variable declarations.
    integer :: ilnobl

    integer :: iz1
    integer :: iz3
    integer :: inote
    integer :: j
    integer :: j1
    integer :: j2
    integer :: j3
    integer :: n
    integer :: nmu
    integer :: nmu1
    integer :: nmu2
    integer :: nnote
    integer :: ns1
    integer :: ns2
    integer :: ns3

    real(kind=8) :: az1
    real(kind=8) :: az2
    real(kind=8) :: az3
    real(kind=8) :: cx
    real(kind=8) :: cx1
    real(kind=8) :: cx2
    real(kind=8) :: z1
    real(kind=8) :: z2
    real(kind=8) :: z3

    ! Recalculate the original coefficients for the Cphi, psi, and
    ! zeta interaction parameters from the derived set of coefficients
    ! for conventional mu interaction parameters.
    ! Take advantage of the expected structure of the data for the
    ! mu parameters. These occur in the following order:
    !   cca data (obtained from Cphi data)
    !   aac data (obtained from the same set of Cphi data)
    !   cc'a data (obtained from the psi(cc'a) data)
    !   aa'c data (obtained from the psi(aa'c) data)
    !   nnn data (mu data, stays the same)
    !   nnn' data (mu data, stays the same)
    !   nca data (obtained from zeta data)
    ! The species within a triplet follow the above pattern.
    ! Thus, so see if one species appears twice, it is only
    ! necessary to see if the first and second species are
    ! identical. Thus, cca data will always be stored as
    ! "cca", never as "cac" or "acc".
    ! The strategy here is as follows: store all converted data
    ! in the original amu array, as originally indexed. The nmxi
    ! "pointer" array can retain its structure as well.
    ! The first step is to convert the mu(cca) back to Cphi(ca)
    ! data. The mu(aac) are likewise converted back to Cphi(ac) data.
    ! This is the same set of Cphi data, just available with different
    ! indexing (Cphi(ac) = Cphi(ca)). This duplication allows the amu
    ! array to be used without changing its structure.
    ! The mu(cc'a) data are converted back to psi(cc'a) data,
    ! and the mu(aa'c) data are likewise converted back to
    ! psi(aa'c) data. The two psi data sets here are separate. Note
    ! that the conversion back to psi data uses the recalculated
    ! Cphi data. Each psi datum has two "cognate" Cphi data.
    ! These cognate Cphi data are for the constituent binary
    ! pairs (e.g., the cognates for psi(cc'a) are Cphi(ca) and
    ! Cphi(c'a). If the data for a cognate Cphi are zero, the
    ! corresponding mu data may not be loaded for use during the
    ! run. In that case, the data will not be found among the
    ! recalculated Cphi data. Zero values will then be assumed.
    ! The mu(nnn) and mu(nnn') data were originally in the same form
    ! and remain that way.
    ! The mu(nca) data are converted back to zeta(nca) data. This is
    ! a simple conversion, as the potential cognate parameters
    ! (mu(ncc) and mu(naa)) are not used in the standard Pitzer model
    ! incorporated into EQ3/6 and other geochemical modeling codes.
    inote = 0
    nnote = 0

    ! Make the mu(cca) to Cphi(ca) conversions.
    ifcphi1 = 1
    ilcphi1 = 0

    do nmu = ifcphi1,nmut
        ns1 = nmux(1,nmu)
        ns2 = nmux(2,nmu)
        ns3 = nmux(3,nmu)
        z1 = zchar(ns1)
        z3 = zchar(ns3)
        az1 = abs(z1)
        az3 = abs(z3)

        if ((ns1.eq.ns2 .and. ns2.ne.ns3) .and.    (z1.gt.0 .and. z3.lt.0.)) then
            ! Have an instance of mu(cca) data. Convert the coefficients
            ! for mu(cca) to those for Cphi(ca).
            cx = 6.*sqrt(az3/az1)

            do j = 1,jpfcmx
                amu(j,nmu) = cx*amu(j,nmu)
            end do

            ilcphi1 = nmu
        else
            ! Have finished finding and processing all mu(cca) data.
            go to 110
        end if
    end do

110 continue

    ! Make the mu(aac) to Cphi(ac) conversions.
    ifcphi2 = ilcphi1 + 1
    ilcphi2 = ilcphi1

    do nmu = ifcphi2,nmut
        ns1 = nmux(1,nmu)
        ns2 = nmux(2,nmu)
        ns3 = nmux(3,nmu)
        z1 = zchar(ns1)
        z3 = zchar(ns3)
        az1 = abs(z1)
        az3 = abs(z3)

        if ((ns1.eq.ns2 .and. ns2.ne.ns3) .and.    (z1.lt.0 .and. z3.gt.0.)) then
            ! Have an instance of mu(aac) data. Convert the coefficients
            ! for mu(aac) to those for Cphi(ac).
            cx = 6.*sqrt(az3/az1)

            do j = 1,jpfcmx
                amu(j,nmu) = cx*amu(j,nmu)
            end do

            ilcphi2 = nmu
        else
            ! Have finished finding and processing all mu(aac) data.
            go to 120
        end if
    end do

120 continue

    ! Make the mu(cc'a) to psi(cc'a) conversions.
    ifpsi1 = ilcphi2 + 1
    ilpsi1 = ilcphi2

    do nmu = ifpsi1,nmut
        ns1 = nmux(1,nmu)
        ns2 = nmux(2,nmu)
        ns3 = nmux(3,nmu)
        z1 = zchar(ns1)
        z2 = zchar(ns2)
        z3 = zchar(ns3)
        az1 = abs(z1)
        az2 = abs(z2)
        az3 = abs(z3)

        if ((ns1.ne.ns2 .and. ns2.ne.ns3) .and.    (z1.gt.0. .and. z2.gt.0. .and. z3.lt.0.)) then
            ! Have an instance of mu(cc'a) data. Convert the coefficients
            ! for mu(cc'a) to those for psi(cc'a).
            ! First find the indices of the cognate mu(cca) and
            ! mu(c'c'a) data.
            inote = 0

            nmu1 = 0

            do n = ifcphi1,ilcphi1
                if (nmux(1,n).eq.ns1 .and. nmux(3,n).eq.ns3) then
                    nmu1 = n
                    go to 130
                end if
            end do

            if (iodb(1) .gt. 0) then
                j1 = ilnobl(uspec(ns1)(1:24))
                j2 = ilnobl(uspec(ns2)(1:24))
                j3 = ilnobl(uspec(ns3)(1:24))
                write (noutpt,1000) uspec(ns1)(1:j1),uspec(ns3)(1:j3),uspec(ns1)(1:j1),uspec(ns2)(1:j2),uspec(ns3)(1:j3)
                write (nttyo,1000) uspec(ns1)(1:j1),uspec(ns3)(1:j3),uspec(ns1)(1:j1),uspec(ns2)(1:j2),uspec(ns3)(1:j3)
1000 format(/' * Note - (EQLIBG/rc3ocf) Could not find the',' Cphi data for',/7x,a,', ',a,'. They are needed to',' recalculate',/7x,'the psi data for ',a,', ',a,', ',a,'.',/7x,'Zero values will be assumed.')

                inote = inote + 1
                nnote = nnote + 1
            end if

130 continue

            nmu2 = 0

            do n = ifcphi1,ilcphi1
                if (nmux(1,n).eq.ns2 .and. nmux(3,n).eq.ns3) then
                    nmu2 = n
                    go to 140
                end if
            end do

            if (iodb(1) .gt. 0) then
                j1 = ilnobl(uspec(ns1)(1:24))
                j2 = ilnobl(uspec(ns2)(1:24))
                j3 = ilnobl(uspec(ns3)(1:24))
                write (noutpt,1000) uspec(ns2)(1:j1),uspec(ns3)(1:j3),uspec(ns1)(1:j1),uspec(ns2)(1:j2),uspec(ns3)(1:j3)
                write (nttyo,1000) uspec(ns2)(1:j1),uspec(ns3)(1:j3),uspec(ns1)(1:j1),uspec(ns2)(1:j2),uspec(ns3)(1:j3)
                inote = inote + 1
                nnote = nnote + 1
            end if

140 continue

            do j = 1,jpfcmx
                amu(j,nmu) = 6.*amu(j,nmu)
            end do

            if (nmu1 .gt. 0) then
                cx1 = 0.5*az2/sqrt(az1*az3)

                do j = 1,jpfcmx
                    amu(j,nmu) = amu(j,nmu) - cx1*amu(j,nmu1)
                end do
            end if

            if (nmu2 .gt. 0) then
                cx2 = 0.5*az1/sqrt(az2*az3)

                do j = 1,jpfcmx
                    amu(j,nmu) = amu(j,nmu) - cx2*amu(j,nmu2)
                end do
            end if

            ilpsi1 = nmu
        else
            ! Have finished finding and processing all mu(aac) data.
            go to 150
        end if
    end do

150 continue

    ! Make the mu(aa'c) to psi(aa'c) conversions.
    ifpsi2 = ilpsi1 + 1
    ilpsi2 = ilpsi1

    do nmu = ifpsi2,nmut
        ns1 = nmux(1,nmu)
        ns2 = nmux(2,nmu)
        ns3 = nmux(3,nmu)
        z1 = zchar(ns1)
        z2 = zchar(ns2)
        z3 = zchar(ns3)
        az1 = abs(z1)
        az2 = abs(z2)
        az3 = abs(z3)

        if ((ns1.ne.ns2 .and. ns2.ne.ns3) .and.    (z1.lt.0. .and. z2.lt.0. .and. z3.gt.0.)) then
            ! Have an instance of mu(aa'c) data. Convert the coefficients
            ! for mu(aa'c) to those for psi(aa'c).
            ! First find the indices of the cognate mu(aac) and
            ! mu(a'a'c) data.
            inote = 0

            nmu1 = 0

            do n = ifcphi2,ilcphi2
                if (nmux(1,n).eq.ns1 .and. nmux(3,n).eq.ns3) then
                    nmu1 = n
                    go to 160
                end if
            end do

            if (iodb(1) .gt. 0) then
                j1 = ilnobl(uspec(ns1)(1:24))
                j2 = ilnobl(uspec(ns2)(1:24))
                j3 = ilnobl(uspec(ns3)(1:24))
                write (noutpt,1000) uspec(ns1)(1:j1),uspec(ns3)(1:j3),uspec(ns1)(1:j1),uspec(ns2)(1:j2),uspec(ns3)(1:j3)
                write (nttyo,1000) uspec(ns1)(1:j1),uspec(ns3)(1:j3),uspec(ns1)(1:j1),uspec(ns2)(1:j2),uspec(ns3)(1:j3)
                inote = inote + 1
                nnote = nnote + 1
            end if

160 continue

            nmu2 = 0

            do n = ifcphi2,ilcphi2
                if (nmux(1,n).eq.ns2 .and. nmux(3,n).eq.ns3) then
                    nmu2 = n
                    go to 170
                end if
            end do

            if (iodb(1) .gt. 0) then
                j1 = ilnobl(uspec(ns1)(1:24))
                j2 = ilnobl(uspec(ns2)(1:24))
                j3 = ilnobl(uspec(ns3)(1:24))
                write (noutpt,1000) uspec(ns2)(1:j1),uspec(ns3)(1:j3),uspec(ns1)(1:j1),uspec(ns2)(1:j2),uspec(ns3)(1:j3)
                write (nttyo,1000) uspec(ns2)(1:j1),uspec(ns3)(1:j3),uspec(ns1)(1:j1),uspec(ns2)(1:j2),uspec(ns3)(1:j3)
                inote = inote + 1
                nnote = nnote + 1
            end if

170 continue

            do j = 1,jpfcmx
                amu(j,nmu) = 6.*amu(j,nmu)
            end do

            if (nmu1 .gt. 0) then
                cx1 = 0.5*az2/sqrt(az3*az1)

                do j = 1,jpfcmx
                    amu(j,nmu) = amu(j,nmu) - cx1*amu(j,nmu1)
                end do
            end if

            if (nmu2 .gt. 0) then
                cx2 = 0.5*az1/sqrt(az3*az2)

                do j = 1,jpfcmx
                    amu(j,nmu) = amu(j,nmu) - cx2*amu(j,nmu2)
                end do
            end if

            ilpsi2 = nmu
        else
            ! Have finished finding and processing all mu(aac) data.
            go to 180
        end if
    end do

180 continue

    ! Skip past the mu(nnn) data.
    ifnnn = ilpsi2 + 1
    ilnnn = ilpsi2

    do nmu = ifnnn,nmut
        ns1 = nmux(1,nmu)
        ns2 = nmux(2,nmu)
        ns3 = nmux(3,nmu)
        z1 = zchar(ns1)
        iz1 = nint(z1)

        if ((ns1.eq.ns2 .and. ns2.eq.ns3) .and. iz1.eq.0) then
            ! Have an instance of mu(nnn) data.
            ilnnn = nmu
        else
            ! Have finished finding all mu(nnn) data.
            go to 190
        end if
    end do

190 continue

    ! Skip past the mu(nnn') data.
    ifn2n = ilnnn + 1
    iln2n = ilnnn

    do nmu = ifn2n,nmut
        ns1 = nmux(1,nmu)
        ns2 = nmux(2,nmu)
        ns3 = nmux(3,nmu)
        z1 = zchar(ns1)
        z3 = zchar(ns3)
        iz1 = nint(z1)
        iz3 = nint(z3)

        if ((ns1.eq.ns2 .and. ns1.ne.ns3) .and.    (iz1.eq.0 .and. iz3.eq.0)) then
            ! Have an instance of mu(nnn') data.
            iln2n = nmu
        else
            ! Have finished finding all mu(nnn') data.
            go to 200
        end if
    end do

200 continue

    ifzeta = iln2n + 1
    ilzeta = iln2n

    ! Make the mu(nca) to zeta(nca) conversions.
    do nmu = ifzeta,nmut
        ns1 = nmux(1,nmu)
        ns2 = nmux(2,nmu)
        ns3 = nmux(3,nmu)
        z1 = zchar(ns1)
        z3 = zchar(ns3)
        az1 = abs(z1)
        az3 = abs(z3)

        if ((ns1.ne.ns2 .and. ns2.ne.ns3) .and.    (iz1.eq.0 .and. z2.gt.0. .and. z3.lt.0.)) then
            ! Have an instance of mu(nca) data. Convert the coefficients
            ! for mu(nca) to those for zeta(nca).
            do j = 1,jpfcmx
                amu(j,nmu) = 6.*amu(j,nmu)
            end do

            ilzeta = nmu
        else
            ! Have finished finding and processing all mu(nca) data.
            go to 210
        end if
    end do

210 continue

    if (iodb(1) .gt. 0) then
        ! Check results.
        write (noutpt,1100)
1100 format(/1x,"Recalculated third-order Pitzer coefficients")

        write (noutpt,1110)
1110 format(/1x,"Cphi(ca) data")

        do n = ifcphi1,ilcphi1
            ns1 = nmux(1,n)
            ns3 = nmux(3,n)
            j3 = ilnobl(uspec(ns3)(1:24))
            write (noutpt,1120) uspec(ns1),uspec(ns3)(1:j3)
1120 format(/1x,a24,2x,a)

            write (noutpt,1130)
1130 format(3x,"Cphi:")

            do j = 1,jpfcmx
                write (noutpt,1150) j,amu(j,n)
1150 format (5x,"a",i1," = ",f10.6)
            end do
        end do

        write (noutpt,1200)
1200 format(/1x,"Cphi(ac) data")

        do n = ifcphi2,ilcphi2
            ns1 = nmux(1,n)
            ns3 = nmux(3,n)
            j3 = ilnobl(uspec(ns3)(1:24))
            write (noutpt,1120) uspec(ns1),uspec(ns3)(1:j3)
            write (noutpt,1130)

            do j = 1,jpfcmx
                write (noutpt,1150) j,amu(j,n)
            end do
        end do

        write (noutpt,1310)
1310 format(/1x,"psi(cc'a) data")

        do n = ifpsi1,ilpsi1
            ns1 = nmux(1,n)
            ns2 = nmux(2,n)
            ns3 = nmux(3,n)
            j3 = ilnobl(uspec(ns3)(1:24))
            write (noutpt,1320) uspec(ns1),uspec(ns2),uspec(ns3)(1:j3)
1320 format(/1x,a24,2x,a24,2x,a)

            write (noutpt,1330)
1330 format(3x,"psi:")

            do j = 1,jpfcmx
                write (noutpt,1150) j,amu(j,n)
            end do
        end do

        write (noutpt,1410)
1410 format(/1x,"psi(aa'c) data")

        do n = ifpsi2,ilpsi2
            ns1 = nmux(1,n)
            ns2 = nmux(2,n)
            ns3 = nmux(3,n)
            j3 = ilnobl(uspec(ns3)(1:24))
            write (noutpt,1320) uspec(ns1),uspec(ns2),uspec(ns3)(1:j3)
            write (noutpt,1330)

            do j = 1,jpfcmx
                write (noutpt,1150) j,amu(j,n)
            end do
        end do

        write (noutpt,1510)
1510 format(/1x,"mu(nnn) data")

        do n = ifnnn,ilnnn
            ns1 = nmux(1,n)
            ns2 = nmux(2,n)
            ns3 = nmux(3,n)
            j3 = ilnobl(uspec(ns3)(1:24))
            write (noutpt,1320) uspec(ns1),uspec(ns2),uspec(ns3)(1:j3)
            write (noutpt,1520)
1520 format(3x,"mu:")

            do j = 1,jpfcmx
                write (noutpt,1150) j,amu(j,n)
            end do
        end do

        write (noutpt,1610)
1610 format(/1x,"mu(nnn') data")

        do n = ifn2n,iln2n
            ns1 = nmux(1,n)
            ns2 = nmux(2,n)
            ns3 = nmux(3,n)
            j3 = ilnobl(uspec(ns3)(1:24))
            write (noutpt,1320) uspec(ns1),uspec(ns1),uspec(ns3)(1:j3)
            write (noutpt,1520)

            do j = 1,jpfcmx
                write (noutpt,1150) j,amu(j,n)
            end do
        end do

        write (noutpt,1710)
1710 format(/1x,"zeta(nca) data")

        do n = ifzeta,ilzeta
            ns1 = nmux(1,n)
            ns2 = nmux(2,n)
            ns3 = nmux(3,n)
            j3 = ilnobl(uspec(ns3)(1:24))
            write (noutpt,1320) uspec(ns1),uspec(ns2),uspec(ns3)(1:j3)
            write (noutpt,1720)
1720 format(3x,"zeta:")

            do j = 1,jpfcmx
                write (noutpt,1150) j,amu(j,n)
            end do
        end do
    end if

    ! Transform the Cphi data into the equivalent C data.
    do n = ifcphi1,ilcphi1
        ns1 = nmux(1,n)
        ns3 = nmux(3,n)
        z1 = zchar(ns1)
        z3 = zchar(ns3)
        az1 = abs(z1)
        az3 = abs(z3)
        cx = 1./(2.*sqrt(az1*az3))

        do j = 1,jpfcmx
            amu(j,n) = cx*amu(j,n)
        end do
    end do

    do n = ifcphi2,ilcphi2
        ns1 = nmux(1,n)
        ns3 = nmux(3,n)
        z1 = zchar(ns1)
        z3 = zchar(ns3)
        az1 = abs(z1)
        az3 = abs(z3)
        cx = 1./(2.*sqrt(az1*az3))

        do j = 1,jpfcmx
            amu(j,n) = cx*amu(j,n)
        end do
    end do

999 continue
end subroutine rc3ocf
