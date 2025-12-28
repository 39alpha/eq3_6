subroutine prtrdx(ah,ahrc,cdrsd,eh,ehrc,fo2lg,fo2lrc,jflgi,jsflag,narn1,nbasp,nbaspd,nbt,nbtmax,ndrsd,ndrsmx,ndrsrd,nelect,nhydr,no2gaq,noutpt,nstmax,pe,perc,uspec)
    !! This subroutine prints a table of the Eh, pe-, log fO2, and Ah for
    !! the default redox constraint and each aqueous redox couple which
    !! is not required to satisfy the default redox constraint.
    !! Most of the data printed by this subroutine are calculated by
    !! EQLIB/cdardx.f.
    !! This subroutine is called by:
    !!   EQ3NR/scripx.f
    !!   EQ6/scripz.f
    !! Principal input:
    !!   ah     = Default Ah, kcal
    !!   ahrc   = array of couple-specific Ah values, kcal
    !!   cdrsd  = coefficients of reactions in the 'd' set
    !!   eh     = Default Eh, volts
    !!   ehrc   = array of couple-specific Eh values, volts
    !!   fo2lg  = Default log fO2
    !!   fo2lrc = array of couple-specific log fO2 values
    !!   jsflag = array of species status flags
    !!   nbt    = the number of species in the basis set
    !!   ndrsd  = array of the indices of the species corresponding to
    !!              the reaction coefficients in the cdrsd array
    !!   ndrsrd = array giving the range in the cdrsd and ndrsd arrays
    !!              corresponding to the reaction for a given species
    !!   nelect = index of the fictive species aqueous e-
    !!   no2gaq = index of the fictive species aqueous O2(g)
    !!   pe     = Default pe-
    !!   perc   = array of couple-specific pe- values
    !!   uspec  = array of species names
    !! Principal output:
    !!   None
    implicit none

    ! Calling sequence variable declarations.
    integer :: nbtmax
    integer :: ndrsmx
    integer :: nstmax

    integer :: noutpt

    integer :: jflgi(nbtmax)
    integer :: jsflag(nstmax)
    integer :: nbasp(nbtmax)
    integer :: nbaspd(nbtmax)
    integer :: ndrsd(ndrsmx)
    integer :: ndrsrd(2,nstmax)
    integer :: narn1
    integer :: nbt
    integer :: nelect
    integer :: nhydr
    integer :: no2gaq

    character(len=48) :: uspec(nstmax)

    real(kind=8) :: ahrc(nbtmax)
    real(kind=8) :: cdrsd(ndrsmx)
    real(kind=8) :: ehrc(nbtmax)
    real(kind=8) :: fo2lrc(nbtmax)
    real(kind=8) :: perc(nbtmax)
    real(kind=8) :: ah
    real(kind=8) :: eh
    real(kind=8) :: fo2lg
    real(kind=8) :: pe

    real(kind=8) :: coefdr

    ! Local variable declarations.
    integer :: j2
    integer :: j3
    integer :: n
    integer :: nb
    integer :: nr1
    integer :: nr2
    integer :: ns1
    integer :: ns2
    integer :: nse
    integer :: nt

    integer :: ilnobl

    character(len=24) :: ux24

    real(kind=8) :: cx1

    write (noutpt,1000)
1000 format(//16x,'--- Aqueous Redox Reactions ---',//3x,'Couple',27x,'Eh, volts',6x,'pe-',6x,'log fO2',3x,'Ah, kcal',/)

    write (noutpt,1010) eh,pe,fo2lg,ah
1010 format(1x,'DEFAULT',t37,f7.3,3x,1pe11.4,2x,0pf8.3,2x,f8.3)

    do nb = 1,nbt
        if (jflgi(nb) .eq. 30) then
            go to 110
        end if

        if (jflgi(nb) .eq. 27) then
            go to 110
        end if

        ns1 = nbaspd(nb)
        ns2 = nbasp(nb)

        if (jsflag(ns1) .ne. 0) then
            go to 110
        end if

        if (jsflag(ns2) .ne. 0) then
            go to 110
        end if

        nt = ndrsrd(2,ns1) - ndrsrd(1,ns1) + 1

        if (nt .lt. 2) then
            go to 110
        end if

        ! Test the reaction read from the data file to see if it contains
        ! an entry for O2(g).
        ! Calling sequence substitutions:
        !   cdrsd for cdrs
        !   ndrsd for ndrs
        !   ndrsrd for ndrsr
        !   no2gaq for nse
        !   ns1 for ns
        cx1 = coefdr(cdrsd,ndrsd,ndrsmx,ndrsrd,no2gaq,ns1,nstmax)

        if (cx1 .ne. 0.) then
            ux24 = uspec(narn1)
            nr1 = ndrsrd(1,ns1)
            nr2 = ndrsrd(2,ns1)

            do n = nr1,nr2
                nse = ndrsd(n)

                if (nse .eq. 0) then
                    go to 110
                end if

                if (nse.ne.ns1 .and. nse.ne.no2gaq .and. nse.ne.nelect      .and. nse.ne.narn1 .and. nse.ne.nhydr) then
                    ux24 = uspec(nse)
                end if
            end do

            j2 = ilnobl(uspec(ns1)(1:24))
            j3 = ilnobl(ux24)

100 continue
            if ((j2 + j3) .gt. 32) then
                if (j2 .gt. j3) then
                    j2 = j2 - 1
                else
                    j3 = j3 - 1
                end if

                go to 100
            end if

            write (noutpt,1020) uspec(ns1)(1:j2),ux24(1:j3),ehrc(nb),perc(nb),fo2lrc(nb),ahrc(nb)
1020 format(1x,a,'/',a,t37,f7.3,3x,1pe11.4,2x,0pf8.3,2x,f8.3)
        end if

110 continue
    end do

    write (noutpt,1030)
1030 format(/4x,'Couples required to satisfy the default redox',' constraint are not listed.',/)
end subroutine prtrdx