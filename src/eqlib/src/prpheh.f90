subroutine prpheh(ah,ahmes,ahnbs,eh,ehmes,ehnbs,iopg,nopgmx,noutpt,pch,pe,pemes,penbs,ph,phcl,phmes,phnbs,qphcl,qredox,qrho)
    !! This subroutine prints the pH, Eh (redox potential), and pe
    !! (log electron activity function) on the following pH scales:
    !!   1. The scale corresponding to the activity coefficient
    !!      model equations, unless these are rescaled to be
    !!      consistent with the NBS pH scale or the Mesmer pH scale
    !!      (see iopg(2), discussed below).
    !!   2. The NBS pH scale.
    !!   3. The Mesmer pH scale (pH is numerically equal to -log m(H+)).
    !! Note the following options for the unscaled model pH scale,
    !! which are chosen by the option switch iopg(2):
    !!   = -1   The internal scale is defined for the activity
    !!          coefficient model defined by iopg(1)
    !!   =  0   The modified NBS pH scale (log gamma Cl-
    !!          is defined by the Bates-Guggenheim expression)
    !!   =  2   The Mesmer pH (pmH) scale
    !! This subroutine also prints the pHCl, which is equivalent to
    !! pH + pCl and independent of the choice of pH scale.
    !! The quantities printed by this subroutine are computed by
    !! EQLIB/gpheh.f.
    !! This subroutine is called by:
    !!   EQ3NR/scripx.f
    !!   EQ6/scripz.f
    !! Principal input:
    !!   ah     = Ah corresponding to the pH used in the computation
    !!   ahmes  = Ah corresponding to the phmes
    !!   ahnbs  = Ah corresponding to the phnbs
    !!   pch    = pcH, pH based on the molarity of the hydrogen ion
    !!   ph     = pH on scale used in computation
    !!   phcl   = pHCl
    !!   phmes  = pH on the Mesmer scale (same as the pmH)
    !!   phnbs  = pH on the modified NBS scale
    !!   eh     = Eh corresponding to the pH used in the computation
    !!   ehmes  = Eh corresponding to the phmes
    !!   ehnbs  = Eh corresponding to the phnbs
    !!   pe     = pe function corresonding to Eh
    !!   pemes  = pe function corresonding to ehmes
    !!   penbs  = pe function corresonding to ehnbs
    !!   iopg   = array of activity coefficient option switches
    !!   noutpt = the unit number of the output file
    !!   qrho   = .true. if a solution density value is available
    !!   qpchl  = logical flag, = .true. if pHCl is defined
    !!   qredox = logical flag, = .true. if there is oxidation-reduction
    !!              in the modeled system
    !! Principal output:
    !!   None
    implicit none

    ! Calling sequence variable declarations.
    integer :: nopgmx

    integer :: noutpt

    integer :: iopg(nopgmx)

    logical :: qphcl
    logical :: qredox
    logical :: qrho

    real(kind=8) :: ah
    real(kind=8) :: ahmes
    real(kind=8) :: ahnbs
    real(kind=8) :: eh
    real(kind=8) :: ehmes
    real(kind=8) :: ehnbs
    real(kind=8) :: pch
    real(kind=8) :: pe
    real(kind=8) :: pemes
    real(kind=8) :: penbs
    real(kind=8) :: ph
    real(kind=8) :: phcl
    real(kind=8) :: phmes
    real(kind=8) :: phnbs

    ! Local variable declarations.
    integer :: j2
    integer :: isc

    integer :: ilnobl

    character(len=24) :: uphscs(6)

    data (uphscs(isc), isc = 1,6)  /'Internal pH scale       ','NBS pH scale            ','Mesmer pH (pmH) scale   ','Pitzer pH scale         ','B-dot pH scale          ','Davies pH scale         '/

    isc = 1

    if (iopg(2) .ge. 0) then
        isc = iopg(2) + 2
    else
        if (iopg(1) .eq. 1) then
            isc = 4
        else if (iopg(1) .eq. 0) then
            isc = 5
        else if (iopg(1) .eq. -1) then
            isc = 6
        end if
    end if

    if (qredox) then
        write (noutpt,1000)
1000 format(/11x,'--- The pH, Eh, pe-, and Ah on various pH',' scales ---')

        write (noutpt,1010)
1010 format(/31x,'pH',5x,'Eh, volts',8x,'pe-',7x,'Ah, kcal',/)

        if (isc.ne.2 .and. isc.ne.3) then
            write (noutpt,1020) uphscs(isc),ph,eh,pe,ah
1020 format(1x,a24,2x,f8.4,3x,f8.4,4x,1pe11.4,2x,0pf10.4)
        end if

        write (noutpt,1020) uphscs(2),phnbs,ehnbs,penbs,ahnbs
        write (noutpt,1020) uphscs(3),phmes,ehmes,pemes,ahmes
    else
        write (noutpt,1030)
1030 format(/11x,'--- The pH on various pH scales ---')

        write (noutpt,1040)
1040 format(/40x,'pH',/)

        if (isc.ne.2 .and. isc.ne.3) then
            write (noutpt,1050) uphscs(isc),ph
1050 format(6x,a24,3x,f11.4)
        end if

        write (noutpt,1050) uphscs(2),phnbs
        write (noutpt,1050) uphscs(3),phmes
    end if

    write (noutpt,1100)
1100 format(/1x)

    if (qrho) then
        write (noutpt,1110) pch
1110 format(6x,'pcH= ',f11.4)
    end if

    if (qphcl) then
        write (noutpt,1140) phcl
1140 format(5x,'pHCl= ',f11.4)
    else
        write (noutpt,1150)
1150 format(5x,'The pHCl is undefined because no Cl- is present.')
    end if

    write (noutpt,1160)
1160 format(1x)

    j2 = ilnobl(uphscs(isc))
    write (noutpt,1190) uphscs(isc)(1:j2)
1190 format(/3x,'The single ion activities and activity coefficients',' listed below',/3x,'are consistent with the ',a,'.',/)
end subroutine prpheh