subroutine prtmip(acflg,actlg,conclg,ctb,nbaspd,nbt,nbtmax,nelect,nhydr,nhydx,noutpt,nstmax,uspec,zchar)
    !! This subroutine prints  tables of the mean ionic activities and
    !! activity coefficients.
    !! This subroutine is called by:
    !!   EQ3NR/scripx.f
    !!   EQ6/scripz.f
    !! Principal input:
    !!   acflg  = array of log activity coefficients of species
    !!   actlg  = array of log activities of species
    !!   conclg = array of log concentrations of species
    !!   ctb    = array of total molalities of basis species
    !!   nbt    = the number of species in the basis set
    !!   nelect = index of the fictive species aqueous e-
    !!   nhydr  = index of the species aqueous H+
    !!   nhydx  = index of the species aqueous OH-
    !!   uspec  = array of species names
    !!   zchar  = array of species electrical charge numbers
    !! Principal output:
    !!   None
    implicit none

    ! Calling sequence variable declarations.
    integer :: nbtmax
    integer :: nstmax

    integer :: noutpt

    integer :: nbaspd(nbtmax)
    integer :: nbt
    integer :: nelect
    integer :: nhydr
    integer :: nhydx

    character(len=48) :: uspec(nstmax)

    real(kind=8) :: acflg(nstmax)
    real(kind=8) :: actlg(nstmax)
    real(kind=8) :: conclg(nstmax)
    real(kind=8) :: ctb(nbtmax)
    real(kind=8) :: zchar(nstmax)

    ! Local variable declarations.
    integer :: j2
    integer :: j3
    integer :: nb1
    integer :: nb2
    integer :: ns1
    integer :: ns2

    integer :: ilnobl

    real(kind=8) :: apmlg
    real(kind=8) :: cpmlg
    real(kind=8) :: cpmslg
    real(kind=8) :: ctb1
    real(kind=8) :: ctb2
    real(kind=8) :: gpm
    real(kind=8) :: gpmlg
    real(kind=8) :: gpms
    real(kind=8) :: gpmslg
    real(kind=8) :: zxa
    real(kind=8) :: zx1
    real(kind=8) :: zx2

    real(kind=8) :: texp
    real(kind=8) :: tlg

    ! Compute and print a table of stoichiometric ionic properties.
    write (noutpt,1000)
1000 format(//16x,'--- Stoichiometric Mean Ionic Properties ---')

    write (noutpt,1010)
1010 format(/3x,'Species',18x,'Log a(+/-)',2x,'Log m(+/-)',2x,'Log gamma(+/-)',2x,'gamma(+/-)',/)

    do nb1 = 1,nbt
        ns1 = nbaspd(nb1)
        zx1 = zchar(ns1)

        ctb1 = ctb(nb1)

        if (ctb1 .lt. 0.) then
            if (ns1 .eq. nhydx) then
                ns1 = nhydr
                zx1 = zchar(nhydr)
                ctb1 = -ctb1
            end if
        end if

        if (zx1 .le. 0.) then
            go to 120
        end if

        if (ctb1 .le. 0.) then
            go to 120
        end if

        do nb2 = 1,nbt
            ns2 = nbaspd(nb2)
            zx2 = zchar(ns2)

            ctb2 = ctb(nb2)

            if (ctb2 .lt. 0.) then
                if (ns2 .eq. nhydr) then
                    ns2 = nhydx
                    zx2 = zchar(nhydx)
                    ctb2 = -ctb2
                end if
            end if

            if (zx2 .ge. 0.) then
                go to 110
            end if

            if (ctb2 .le. 0.) then
                go to 110
            end if

            if (ns2 .eq. nelect) then
                go to 110
            end if

            zxa = zx1 - zx2
            apmlg = (zx1*actlg(ns2) - zx2*actlg(ns1))/zxa

            cpmslg = (zx1*tlg(ctb2) - zx2*tlg(ctb1))/zxa

            gpmslg = apmlg - cpmslg
            gpms = texp(gpmslg)

            j2 = ilnobl(uspec(ns1)(1:24))
            j3 = ilnobl(uspec(ns2)(1:24))

100 continue
            if ((j2 + j3) .gt. 24) then
                if (j2 .gt. j3) then
                    j2 = j2 - 1
                else
                    j3 = j3 - 1
                end if

                go to 100
            end if

            write (noutpt,1020) uspec(ns1)(1:j2),uspec(ns2)(1:j3),apmlg,cpmslg,gpmslg,gpms
1020 format(1x,a,'/',a,t29,f9.4,3x,f9.4,4x,f9.4,5x,1pe11.4)

110 continue
        end do

120 continue
    end do

    write (noutpt,1030)
1030 format(/3x,'The stoichiometric mean molalities and activity',' coefficients given',/3x,'above are consistent with the',' sensible composition of the',/3x,'aqueous solution.')

    ! Compute and print a table of ionic properties.
    write (noutpt,1100)
1100 format(//21x,'--- Mean Ionic Properties ---')

    write (noutpt,1010)

    do nb1 = 1,nbt
        ns1 = nbaspd(nb1)
        zx1 = zchar(ns1)

        if (zx1 .le. 0.) then
            go to 220
        end if

        if (conclg(ns1) .le. -99999.) then
            go to 220
        end if

        do nb2 = 1,nbt
            ns2 = nbaspd(nb2)
            zx2 = zchar(ns2)

            if (zx2 .ge. 0.) then
                go to 210
            end if

            if (conclg(ns2) .le. -99999.) then
                go to 210
            end if

            if (ns2 .eq. nelect) then
                go to 210
            end if

            zxa = zx1 - zx2
            apmlg = (zx1*actlg(ns2) - zx2*actlg(ns1))/zxa

            cpmlg = (zx1*conclg(ns2) - zx2*conclg(ns1))/zxa

            gpmlg = (zx1*acflg(ns2) - zx2*acflg(ns1))/zxa
            gpm = texp(gpmlg)

            j2 = ilnobl(uspec(ns1)(1:24))
            j3 = ilnobl(uspec(ns2)(1:24))

200 continue
            if ((j2 + j3) .gt. 32) then
                if (j2 .gt. j3) then
                    j2 = j2 - 1
                else
                    j3 = j3 - 1
                end if

                go to 200
            end if

            write (noutpt,1020) uspec(ns1)(1:j2),uspec(ns2)(1:j3),apmlg,cpmlg,gpmlg,gpm
210 continue
        end do

220 continue
    end do

    write (noutpt,1110)
1110 format(/3x,'The mean molalities and activity coefficients given',' above are',/3x,'consistent with the speciation in the model',' employed.',/)
end subroutine prtmip