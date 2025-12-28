subroutine cdardx(actlg,actwlg,ah,ahrc,cdrsd,eh,ehfac,ehrc,farad,fo2lg,fo2lrc,jsflag,mosp,nbasp,nbaspd,nbt,nbtmax,ndrsd,ndrsmx,ndrsrd,no2gaq,nstmax,pe,perc,ph,xlke,xlksd)
    !! This subroutine computes the default Eh, pe-, and Ah from the
    !! default log fO2, and also computes the Eh, pe-, log fO2, and Ah
    !! for each aqueous redox couple.
    !! This subroutine is called by:
    !!   EQ3NR/scripx.f
    !!   EQ6/cdappl.f
    !! Principal input:
    !!   actlg  = array of log activities of species
    !!   cdrsd  = coefficients of reactions in the 'd' set
    !!   jsflag = array of species status flags
    !!   nbt    = the number of species in the basis set
    !!   ndrsd  = array of the indices of the species corresponding to
    !!              the reaction coefficients in the cdrsd array
    !!   ndrsrd = array giving the range in the cdrsd and ndrsd arrays
    !!              corresponding to the reaction for a given species
    !!   no2gaq = index of the fictive species aqueous O2(g)
    !!   xlke   = the log K for the "Eh" reaction
    !! Principal output:
    !!   fo2lrc = array of couple-specific log fO2 values
    !!   ehrc   = array of couple-specific Eh values
    !!   perc   = array of couple-specific pe values
    !!   ahrc   = array of couple-specific Ah values
    implicit none

    ! Calling sequence variable declarations.
    integer :: nbtmax
    integer :: ndrsmx
    integer :: nstmax

    integer :: jsflag(nstmax)
    integer :: nbasp(nbtmax)
    integer :: nbaspd(nbtmax)
    integer :: ndrsd(ndrsmx)
    integer :: ndrsrd(2,nstmax)
    integer :: nbt
    integer :: no2gaq

    real(kind=8) :: actlg(nstmax)
    real(kind=8) :: ahrc(nbtmax)
    real(kind=8) :: cdrsd(ndrsmx)
    real(kind=8) :: ehrc(nbtmax)
    real(kind=8) :: fo2lrc(nbtmax)
    real(kind=8) :: mosp(nstmax)
    real(kind=8) :: perc(nbtmax)
    real(kind=8) :: xlksd(nstmax)
    real(kind=8) :: actwlg
    real(kind=8) :: ah
    real(kind=8) :: eh
    real(kind=8) :: ehfac
    real(kind=8) :: farad
    real(kind=8) :: fo2lg
    real(kind=8) :: pe
    real(kind=8) :: ph
    real(kind=8) :: xlke

    real(kind=8) :: coefdr

    ! Local variable declarations.
    integer :: n
    integer :: nb
    integer :: nr1
    integer :: nr2
    integer :: ns1
    integer :: ns2
    integer :: nse
    integer :: nt

    real(kind=8) :: ahfac
    real(kind=8) :: cx1
    real(kind=8) :: ehx
    real(kind=8) :: fo2lgx

    ahfac = 0.001*farad

    do nb = 1,nbt
        fo2lrc(nb) = fo2lg
        ehrc(nb) = eh
        perc(nb) = pe
        ahrc(nb) = ah
        ns1 = nbaspd(nb)
        ns2 = nbasp(nb)

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

        if (cx1 .eq. 0.) then
            go to 110
        end if

        if (xlksd(ns1) .lt. 9999999.) then
            fo2lgx = xlksd(ns1)
            nr1 = ndrsrd(1,ns1)
            nr2 = ndrsrd(2,ns1)

            do n = nr1,nr2
                nse = ndrsd(n)

                if (nse .le. 0) then
                    go to 110
                end if

                if (nse .ne. no2gaq) then
                    if (mosp(nse) .le. 0.) then
                        go to 110
                    end if

                    fo2lgx = fo2lgx - cdrsd(n)*actlg(nse)
                end if
            end do

            fo2lgx = fo2lgx/cx1
            ehx = (ehfac/4.)*(fo2lgx - 4.*ph - 2.*actwlg - xlke)

            fo2lrc(nb) = fo2lgx
            ehrc(nb) = ehx
            perc(nb) = ehx/ehfac
            ahrc(nb) = ahfac*ehx
        end if

110 continue
    end do
end subroutine cdardx