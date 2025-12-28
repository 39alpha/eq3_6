subroutine gafsir(actlg,afcnst,affsd,cdrsd,jflagd,ndrsd,ndrsmx,ndrsrd,nst,nstmax,sidrsp,uspec,xlksd)
    !! This subroutine computes the saturation indices of the reactions
    !! for the destruction of all species. The saturation index is
    !! defined as SI = log Q/K, where Q is the activity product and K is
    !! the equilibrium constant of the reaction.
    !! This subroutine is very similar in function to EQLIB/afcalc.f.
    !! They both compute affinties and saturation indices. However, the
    !! present subroutine does this for all species; the other subroutine
    !! does it for a single specified reaction. The former could be
    !! written so that it calls the latter inside a loop; however, this
    !! is not done for the sake of avoiding the overhead of repeatedly
    !! making such a call.
    !! This subroutine is called by:
    !!   EQLIB/gaffsd.f
    !! Principal input:
    !!   actlg  = array of logarithms of thermodynamic activities of
    !!              species
    !!   cdrsd  = coefficients of reactions in the 'd' set
    !!   ndrsd  = indices of species appearing in reactions in the
    !!              'd' set
    !!   ndrsrd = range pointer array for the ndrsd array
    !!   xbar   = mole fraction of a species in the phase to which
    !!              it belongs
    !!   xlksd  = logarithms of equilibrium constants for reactions
    !!            in the 'd' set
    !! Principal output:
    !!   affsd  = affinities of reactions in the 'd' set for the
    !!              dissociation or dissolution of species
    !!   sidrsp = saturation indices (log Q/K) of reactions in the 'd'
    !!              set for the dissociation or dissolution of species
    implicit none

    ! Calling sequence variable declarations.
    integer :: ndrsmx
    integer :: nstmax

    integer :: jflagd(nstmax)
    integer :: ndrsd(ndrsmx)
    integer :: ndrsrd(2,nstmax)
    integer :: nst

    character(len=48) :: uspec(nstmax)

    real(kind=8) :: actlg(nstmax)
    real(kind=8) :: affsd(nstmax)
    real(kind=8) :: cdrsd(ndrsmx)
    real(kind=8) :: sidrsp(nstmax)
    real(kind=8) :: xlksd(nstmax)
    real(kind=8) :: afcnst

    ! Local variable declarations.
    integer :: n
    integer :: nr1
    integer :: nr2
    integer :: ns
    integer :: nse
    integer :: nt

    real(kind=8) :: sx

    ! Zero the sidrsp and affsd arrays. This sets the correct values
    ! for species which have jflagd values of 30 and for species whose
    ! reactions are identity reactions.
    do ns = 1,nst
        sidrsp(ns) = 0.
        affsd(ns) = 0.
    end do

    do ns = 1,nst
        if (jflagd(ns).ne.30 .and. jflagd(ns).ne.27) then
            nr1 = ndrsrd(1,ns)
            nr2 = ndrsrd(2,ns)
            nt = nr2 - nr1 + 1

            if (nt .lt. 2) then
                go to 120
            end if

            if (xlksd(ns) .le. -9999999.) then
                sidrsp(ns) = 9999999.
                affsd(ns) = 9999999.
            else if (xlksd(ns) .ge. 9999999.) then
                sidrsp(ns) = -9999999.
                affsd(ns) = -9999999.
            else
                sx = -xlksd(ns)

                do n = nr1,nr2
                    nse = ndrsd(n)

                    if (actlg(nse) .gt. -99999.) then
                        sx = sx + cdrsd(n)*actlg(nse)
                    else
                        sx = -9999999.
                        go to 110
                    end if
                end do

110 continue

                sidrsp(ns) = sx
                affsd(ns) = afcnst*sx
            end if
        end if

120 continue
    end do
end subroutine gafsir