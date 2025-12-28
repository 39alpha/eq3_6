subroutine gpheh(acflg,actlg,actwlg,adh,ah,ahmes,ahnbs,conc,eh,ehfac,ehmes,ehnbs,farad,fo2lg,fxi,iopg,mrmlra,nchlor,nhydr,nopgmx,noutpt,nstmax,nttyo,pch,pe,pemes,penbs,ph,phcl,phmes,phnbs,qphcl,qredox,qrho,xlke)
    !! This subroutine computes the pH, redox potential, and pe
    !! (electron activity function) on the operational pH scale used
    !! in the calculations (see iopg(2)), the NBS pH scale,
    !! and the scale on which pH is numerically equal to -log m(H+).
    !! It also computes the pHCl, which is equivalent to pH + pCl.
    !! This subroutine is called by:
    !!   EQ3NR/scripx.f
    !!   EQ6/path.f
    !! Principal input:
    !!   acflg  = array of logarithms of activity coefficients
    !!   actlg  = array of logarithms of aqueous species activities;
    !!              actlg(nchlor) is the log activity of Cl-
    !!   actwlg = the log activity of water
    !!   adh    = the Debye-Huckel A(gamma) parameter
    !!   conc   = array of species concentrations
    !!   ehfac  = the Eh factor
    !!   farad  = the Faraday constant
    !!   iopg   = array of activity coefficient option switches;
    !!              iopg(2) defines the operational pH scale used in
    !!              the computation:
    !!                -1 = Internal scale corresponding to the activity
    !!                       coefficient model defined by iopg(1)
    !!                 0 = The NBS pH scale (log gamma Cl- is defined
    !!                       by the Bates-Guggenheim expression)
    !!                 1 = The Mesmer pH scale (log gamma H+ is defined
    !!                       as zero)
    !!   fo2lg  = log oxygen fugacity
    !!   nchlor = species index of the chloride ion
    !!   nhydr  = species index of the hydrogen ion
    !!   qredox = logical flag, = .true. if there is oxidation-reduction
    !!              in the modeled system
    !!   qrho   = .true. if a solution density value is available
    !!   fxi    = the ionic strength
    !!   xlke   = log K for special reaction for calculating Eh
    !!             from log oxygen fugacity
    !! Principal output:
    !!   ah     = Ah corresponding to the pH used in the computation
    !!   ahmes  = Ah corresponding to the Mesmer pH scale
    !!   ahnbs  = Ah corresponding to the NBS pH scale
    !!   eh     = Eh corresponding to the pH used in the computation
    !!   ehmes  = Eh corresponding to the Mesmer pH scale
    !!   ehnbs  = Eh corresponding to the NBS pH scale
    !!   ph     = pH on scale used in computation
    !!   phcl   = pHCl
    !!   phmes  = pH on the Mesmer scale
    !!   phnbs  = pH on the NBS scale
    !!   pe     = pe function corresonding to Eh
    !!   pemes  = pe function corresonding to the Mesmer pH scale
    !!   penbs  = pe function corresonding to the NBS pH scale
    !!   qhcl   = logical flag, = .true. if pHCl is defined
    implicit none

    ! Calling sequence variable declarations.
    integer :: nopgmx
    integer :: nstmax

    integer :: noutpt
    integer :: nttyo

    integer :: iopg(nopgmx)
    integer :: nchlor
    integer :: nhydr

    logical :: qphcl
    logical :: qredox
    logical :: qrho

    real(kind=8) :: tlg

    real(kind=8) :: acflg(nstmax)
    real(kind=8) :: actlg(nstmax)
    real(kind=8) :: conc(nstmax)
    real(kind=8) :: actwlg
    real(kind=8) :: adh
    real(kind=8) :: ah
    real(kind=8) :: ahmes
    real(kind=8) :: ahnbs
    real(kind=8) :: eh
    real(kind=8) :: ehfac
    real(kind=8) :: ehmes
    real(kind=8) :: ehnbs
    real(kind=8) :: farad
    real(kind=8) :: fo2lg
    real(kind=8) :: fxi
    real(kind=8) :: mrmlra
    real(kind=8) :: pch
    real(kind=8) :: pe
    real(kind=8) :: pemes
    real(kind=8) :: penbs
    real(kind=8) :: ph
    real(kind=8) :: phcl
    real(kind=8) :: phmes
    real(kind=8) :: phnbs
    real(kind=8) :: xlke

    ! Local variable declarations.
    real(kind=8) :: acfnbs
    real(kind=8) :: ahfx
    real(kind=8) :: delacf
    real(kind=8) :: ehfx
    real(kind=8) :: eterm
    real(kind=8) :: eterm0
    real(kind=8) :: eterm1

    ! Compute pH, Eh, and pe for the operational pH scale.
    ph = -actlg(nhydr)
    eh = -99999.
    pe = -99999.
    ah = -99999.
    ehfx = 0.25*ehfac
    ahfx = 0.001*farad

    if (qredox) then
        eterm = fo2lg - 4.*ph - 2.*actwlg - xlke
        eh = ehfx*eterm
        pe = eh/ehfac
        ah = ahfx*eh
    end if

    ! Compute pH, Eh, and pe consistent with the NBS pH scale.
    ! Note that this requires the activity coefficient of the chloride
    ! ion on the operational scale.
    if (iopg(2) .eq. 0)  then
        phnbs = ph
        ehnbs = eh
        penbs = pe
        ahnbs = ah
    else
        phnbs = -99999.
        ehnbs = -99999.
        penbs = -99999.
        ahnbs = -99999.

        if (nchlor .gt. 0) then
            call nbsgam(acfnbs,adh,fxi,nchlor,noutpt,nttyo)
            delacf = acfnbs - acflg(nchlor)
            phnbs = ph + delacf

            if (qredox) then
                eterm0 = eterm + 4.*(ph - phnbs)
                ehnbs = ehfx*eterm0
                penbs = ehnbs/ehfac
                ahnbs = ahfx*ehnbs
            end if
        end if
    end if

    ! Compute pH, Eh, and pe consistent with the Mesmer pH scale.
    ! Mesmer pH is the same as pmH.
    ehmes = -99999.
    pemes = -99999.

    if (iopg(2) .eq. 1) then
        phmes = ph
        ehmes = eh
        pemes = pe
        ahmes = ah
    else
        delacf = - acflg(nhydr)
        phmes = ph - delacf

        if (qredox) then
            eterm1 = eterm + 4.*(ph - phmes)
            ehmes = ehfx*eterm1
            pemes = ehmes/ehfac
            ahmes = ahfx*ehmes
        end if
    end if

    ! Compute pcH.
    pch = 99999.

    if (qrho) then
        pch = phmes - tlg(mrmlra)
    end if

    ! Compute pHCl = pH + pCl
    qphcl = .false.
    phcl = -99999.

    if (nchlor .gt. 0) then
        if (conc(nchlor) .gt. 0.) then
            qphcl = .true.
            phcl = ph - actlg(nchlor)
        end if
    end if
end subroutine gpheh