subroutine prtsia(affsd,jflagd,jflgi,jsflag,narn1,narn2,nbasp,nbaspd,nbt,nbtmax,ndrsd,ndrsmx,ndrsrd,nhydr,noutpt,nrdxsp,nstmax,sidrsp,uspec)
    !! This subroutine prints tables of saturation indices and affinities
    !! for reactions among aqueous species which are not constrained
    !! to be at equilibrium.
    !! This subroutine is called by:
    !!   EQ3NR/scripx.f
    !!   EQ6/scripz.f
    !! Principal input:
    !!   affsd  = array of affinities for reactions corresponding to
    !!              the various species, as these reactions were read
    !!              from the data file
    !!   jflagd = array of species control flags, used with the 'd'
    !!              set of reactions
    !!   jsflag = array of status flags for species
    !!   nbasp  = array of the indices of the species in the active
    !!              basis set
    !!   nbaspd = array of the indices of the species in the data file
    !!              basis set
    !!   nbt    = the number of species in the basis set
    !!   ndrsd  = array of the indices of the species corresponding to
    !!              the reaction coefficients in the cdrsd array
    !!   ndrsrd = array giving the range in the cdrsd and ndrsd arrays
    !!              corresponding to the reaction for a given species
    !!   nhydr  = index of the species aqueous H+
    !!   nrdxsp = index of the redox basis species
    !!   sidrsp = array of saturation indices for reactions
    !!              corresponding to the various species, as these
    !!              reactions were read from the data file
    !!   uspec  = array of species names
    !! Principal output:
    !!   None
    implicit none

    ! Calling sequence variable declarations.
    integer :: nbtmax
    integer :: ndrsmx
    integer :: nstmax

    integer :: noutpt

    integer :: jflagd(nstmax)
    integer :: jflgi(nbtmax)
    integer :: jsflag(nstmax)
    integer :: nbasp(nbtmax)
    integer :: nbaspd(nbtmax)
    integer :: ndrsd(ndrsmx)
    integer :: ndrsrd(2,nstmax)
    integer :: narn1
    integer :: narn2
    integer :: nbt
    integer :: nhydr
    integer :: nrdxsp

    character(len=48) :: uspec(nstmax)

    real(kind=8) :: affsd(nstmax)
    real(kind=8) :: sidrsp(nstmax)

    ! Local variable declarations.
    integer :: j2
    integer :: j3
    integer :: kount
    integer :: n
    integer :: nb
    integer :: nr1
    integer :: nr2
    integer :: nse
    integer :: ns1
    integer :: ns2
    integer :: nt

    integer :: ilnobl

    character(len=24) :: ux24

    character(len=16) :: ux16a
    character(len=16) :: ux16b

    write (noutpt,1000)
1000 format(//6x,'--- Saturation States of Aqueous Reactions Not',' Fixed at Equilibrium ---')

    write (noutpt,1010)
1010 format(/3x,'Reaction',27x,'Log Q/K',4x,'Affinity, kcal',/)

    kount = 0

    do nb = 1,nbt
        ns1 = nbaspd(nb)

        if (ns1.lt.narn1 .or. ns1.gt.narn2) then
            go to 250
        end if

        if (jflagd(ns1) .eq. 27) then
            go to 250
        end if

        if (jflagd(ns1) .eq. 30) then
            go to 250
        end if

        if (jsflag(ns1) .gt. 0) then
            go to 250
        end if

        nr1 = ndrsrd(1,ns1)
        nr2 = ndrsrd(2,ns1)
        nt = nr2 - nr1 + 1

        if (nt .lt. 2) then
            go to 250
        end if

        ! Find the appropriate matching species to complete the
        ! description of the reaction.
        ns2 = 0

        do n = nr1,nr2
            nse = ndrsd(n)

            if (nse.ne.ns1 .and. nse.ne.nrdxsp      .and. nse.ne.narn1 .and. nse.ne.nhydr) then
                ns2 = nse
                go to 130
            end if
        end do

130 continue

        if (ns2 .eq. 0) then
            do n = nr1,nr2
                nse = ndrsd(n)

                if (nse.ne.ns1 .and. nse.eq.narn1) then
                    ns2 = narn1
                    go to 140
                end if
            end do

140 continue
        end if

        if (ns2 .eq. 0) then
            do n = nr1,nr2
                nse = ndrsd(n)

                if (nse.ne.ns1 .and. nse.eq.nhydr) then
                    ns2 = nhydr
                    go to 150
                end if
            end do

150 continue
        end if

        if (ns2 .eq. 0) then
            do n = nr1,nr2
                nse = ndrsd(n)

                if (nse.ne.ns1 .and. nse.eq.nrdxsp) then
                    ns2 = nrdxsp
                    go to 160
                end if
            end do

160 continue
        end if

        ! At this point, ns2 should not be zero.
        if (ns2 .eq. 0) then
            go to 250
        end if

        if (jsflag(ns2) .gt. 0) then
            go to 250
        end if

        ux24 = uspec(ns2)(1:24)

        kount = kount + 1

        j2 = ilnobl(uspec(ns1)(1:24))
        j3 = ilnobl(ux24)

240 continue
        if ((j2 + j3) .gt. 32) then
            if (j2 .gt. j3) then
                j2 = j2 - 1
            else
                j3 = j3 - 1
            end if

            go to 240
        end if

        if (sidrsp(ns1).gt.-9999999.   .and. sidrsp(ns1).lt.9999999.) then
            write (ux16a,'(f10.5)') sidrsp(ns1)
            write (ux16b,'(f10.5)') affsd(ns1)
        else
            ux16a = '    N/A   '
            ux16b = '    N/A   '
        end if

        write (noutpt,1020) uspec(ns1)(1:j2),ux24(1:j3),ux16a(1:10),ux16b(1:10)
1020 format(1x,a,'/',a,t37,a,3x,a)

250 continue
    end do

    if (kount.le.0) then
        write (noutpt,1030)
    end if

1030 format(1x,'None')

    write (noutpt,1040)
1040 format(1x)
end subroutine prtsia
