subroutine gaffsd(actlg,afcnst,affpd,affsd,cdrsd,jflagd,jpflag,ncmpr,ndrsd,ndrsmx,ndrsrd,npt,nptmax,nst,nstmax,qxknph,sidrph,sidrsp,uphase,uspec,xbar,xlksd)
    !! This subroutine computes affinities and saturation indices based
    !! on reactions in the 'd' set (cdrsd/ndrsd/ndrsrd arrays).
    !! This subroutine is called by:
    !!   EQ3NR/scripx.f
    !!   EQ6/cdappl.f
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
    !!              in the 'd' set
    !!   qxknph = flag indicating if the composition of the phase is
    !!              known; this is the composition which maximizes the
    !!              affinity (the actual composition if the phase is
    !!              in equilibrium with the aqueous solution)
    !! Principal output:
    !!   affsd  = affinities of reactions in the 'd' set for the
    !!              dissociation or dissolution of species
    !!   affpd  = affinities for the dissolution of phases, based on
    !!              reactions in the 'd' set
    !!   sidrsp = saturation indices (log Q/K) of reactions in the 'd'
    !!              set for the dissociation or dissolution of species
    !!   sidrph = saturation indices of phases, based on reactions in
    !!              the 'd' set
    implicit none

    ! Calling sequence variable declarations.
    integer :: ndrsmx
    integer :: nptmax
    integer :: nstmax

    integer :: jflagd(nstmax)
    integer :: jpflag(nptmax)
    integer :: ncmpr(2,nptmax)
    integer :: ndrsd(ndrsmx)
    integer :: ndrsrd(2,nstmax)
    integer :: npt
    integer :: nst

    logical :: qxknph(nptmax)

    character(len=48) :: uspec(nstmax)
    character(len=24) :: uphase(nptmax)

    real(kind=8) :: actlg(nstmax)
    real(kind=8) :: affpd(nptmax)
    real(kind=8) :: affsd(nstmax)
    real(kind=8) :: cdrsd(ndrsmx)
    real(kind=8) :: sidrph(nptmax)
    real(kind=8) :: sidrsp(nstmax)
    real(kind=8) :: xbar(nstmax)
    real(kind=8) :: xlksd(nstmax)
    real(kind=8) :: afcnst

    ! Local variable declarations.
    integer :: np
    integer :: nr1
    integer :: nr2
    integer :: ns
    integer :: nt

    ! Compute affinities and saturation indices for reactions in the
    ! 'd' set.
    call gafsir(actlg,afcnst,affsd,cdrsd,jflagd,ndrsd,ndrsmx,ndrsrd,nst,nstmax,sidrsp,uspec,xlksd)

    ! Compute affinities and saturation indices for phases.
    do np = 1,npt
        affpd(np) = 0.
        sidrph(np) = 0.

        if (jpflag(np) .lt. 2) then
            nr1 = ncmpr(1,np)
            nr2 = ncmpr(2,np)
            nt = nr2 - nr1 + 1

            if (nt .eq. 1) then
                affpd(np) = affsd(nr1)
                sidrph(np) = sidrsp(nr1)
            else
                if (qxknph(np)) then
                    ! The composition of the phase is known.
                    do ns = nr1,nr2
                        sidrph(np) = sidrph(np) + xbar(ns)*sidrsp(ns)
                    end do

                    affpd(np) = afcnst*sidrph(np)
                else
                    ! The composition of the phase is not known.
                    affpd(np) = -9999999.
                    sidrph(np) = -9999999.
                end if
            end if
        end if
    end do
end subroutine gaffsd