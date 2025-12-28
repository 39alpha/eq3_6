subroutine prtvpa(abar,acfw,acfwlg,actw,actwlg,a3bar,fje,fjest,fo2,fo2lg,fxi,fxist,iopg,mlmrra,mrmlra,nopgmx,noutpt,osc,oscst,qrho,rhoc,rhowc,sigmam,sigmst,tdsglw,tdspkc,tdsplc,vosol,wfh2o,wftds,woh2o,wosol,wotds,xbarw,xbrwlg)
    !! This subroutine prints a table of various aqueous solution
    !! parameters.
    !! This subroutine is called by:
    !!   EQ3NR/scripx.f
    !!   EQ6/scripz.f
    !! Principal input:
    !!   abar   = average hard core diameter
    !!   acfw   = activity coefficient of water
    !!   acfwlg = log activity coefficient of water
    !!   actw   = activity of water
    !!   actwlg = log activity of water
    !!   a3bar  = average cubed hard core diameter
    !!   fje    = the ionic asymmetry (the 3rd-order electrostatic
    !!              moment function J)
    !!   fjest  = stoichiometric ionic asymmetry
    !!   fo2    = oxygen fugacity
    !!   fo2lg  = log oxygen fugacity
    !!   fxi    = the ionic strength (the 2nd-order electrostatic
    !!              moment function I)
    !!   fxist  = stoichiometric ionic strength
    !!   mlmrra = molality/molarity ratio
    !!   mrmlra = molarity/molality ratio
    !!   osc    = osmotic coefficient
    !!   oscst  = stoichiometric osmotic coefficient
    !!   qrho   = .true. if a solution density value is available
    !!   rhoc   = the solution density (g/mL)
    !!   rhowc  = the solution density (g/L)
    !!   sigmam = sum of solute molalities
    !!   sigmst = sum of stoichiometric solute molalities
    !!   tdsglw = TDS (total dissolved solutes) (g/L)
    !!   tdspkc = TDS (total dissolved solutes) (mg/kg.sol)
    !!   tdsplc = TDS (total dissolved solutes) (mg/L)
    !!   vosol  = Volume (L) of aqueous solution
    !!   wfh2o  = weight (mass) fraction of water in aqueous solution
    !!   wftds  = weight (mass) fraction of solutes in aqueous solution
    !!   woh2o  = weight (mass) of solvent water
    !!   wosol  = weight (mass) of aqueous solution
    !!   wotds  = weight (mass) of total dissolved solutes
    !!   xbarw  = mole fraction of water in aqueous solution
    !!   xbrwlg = log mole fraction of water in aqueous solution
    !! Principal output:
    !!   None
    implicit none

    ! Calling sequence variable declarations.
    integer :: nopgmx

    integer :: noutpt

    integer :: iopg(nopgmx)

    logical :: qrho

    real(kind=8) :: abar
    real(kind=8) :: acfw
    real(kind=8) :: acfwlg
    real(kind=8) :: actw
    real(kind=8) :: actwlg
    real(kind=8) :: a3bar
    real(kind=8) :: fje
    real(kind=8) :: fjest
    real(kind=8) :: fo2
    real(kind=8) :: fo2lg
    real(kind=8) :: fxi
    real(kind=8) :: fxist
    real(kind=8) :: osc
    real(kind=8) :: oscst
    real(kind=8) :: sigmam
    real(kind=8) :: sigmst
    real(kind=8) :: vosol
    real(kind=8) :: wfh2o
    real(kind=8) :: wftds
    real(kind=8) :: woh2o
    real(kind=8) :: wosol
    real(kind=8) :: wotds
    real(kind=8) :: xbarw
    real(kind=8) :: xbrwlg

    real(kind=8) :: mlmrra
    real(kind=8) :: mrmlra
    real(kind=8) :: rhoc
    real(kind=8) :: rhowc
    real(kind=8) :: tdsglw
    real(kind=8) :: tdspkc
    real(kind=8) :: tdsplc

    ! Local variable declarations.
    integer :: j1
    integer :: j2
    integer :: j3

    integer :: ilnobl

    character(len=16) :: ux16a
    character(len=16) :: ux16b
    character(len=16) :: ux16c

    if (fo2lg .gt. -99999.) then
        write (ux16a,'(1pg12.5)') fo2
        call lejust(ux16a)
        j1 = ilnobl(ux16a)
        write (ux16b,'(1pg12.5)') fo2lg
        call lejust(ux16b)
        j2 = ilnobl(ux16b)
        write (noutpt,1000) ux16a(1:j1),ux16b(1:j2)
1000 format(//20x,'Oxygen fugacity= ',a,' bars',/16x,'Log oxygen fugacity= ',a)
    else
        write (noutpt,1010)
1010 format(/8x,'Oxygen fugacity is not defined.')
    end if

    write (ux16a,'(1pg12.5)') actw
    call lejust(ux16a)
    j1 = ilnobl(ux16a)
    write (ux16b,'(1pg12.5)') actwlg
    call lejust(ux16b)
    j2 = ilnobl(ux16b)
    write (noutpt,1020) ux16a(1:j1),ux16b(1:j2)
1020 format(/18x,'Activity of water= ',a,/14x,'Log activity of water= ',a)

    write (ux16a,'(1pg12.5)') xbarw
    call lejust(ux16a)
    j1 = ilnobl(ux16a)
    write (ux16b,'(1pg12.5)') xbrwlg
    call lejust(ux16b)
    j2 = ilnobl(ux16b)
    write (noutpt,1030) ux16a(1:j1),ux16b(1:j2)
1030 format(/13x,'Mole fraction of water= ',a,/9x,'Log mole fraction of water= ',a)

    write (ux16a,'(1pg12.5)') acfw
    call lejust(ux16a)
    j1 = ilnobl(ux16a)
    write (ux16b,'(1pg12.5)') acfwlg
    call lejust(ux16b)
    j2 = ilnobl(ux16b)
    write (noutpt,1040) ux16a(1:j1),ux16b(1:j2)
1040 format(/6x,'Activity coefficient of water= ',a,/2x,'Log activity coefficient of water= ',a)

    write (ux16a,'(1pg12.5)') osc
    call lejust(ux16a)
    j1 = ilnobl(ux16a)
    write (ux16b,'(1pg12.5)') oscst
    call lejust(ux16b)
    j2 = ilnobl(ux16b)
    write (noutpt,1050) ux16a(1:j1),ux16b(1:j2)
1050 format(/16x,'Osmotic coefficient= ',a,/1x,'Stoichiometric osmotic coefficient= ',a)

    write (ux16a,'(1pg12.5)') sigmam
    call lejust(ux16a)
    j1 = ilnobl(ux16a)
    write (ux16b,'(1pg12.5)') sigmst
    call lejust(ux16b)
    j2 = ilnobl(ux16b)
    write (noutpt,1060) ux16a(1:j1),ux16b(1:j2)
1060 format(/18x,'Sum of molalities= ',a,/3x,'Sum of stoichiometric molalities= ',a)

    write (ux16a,'(1pg12.5)') fxi
    call lejust(ux16a)
    j1 = ilnobl(ux16a)
    write (ux16b,'(1pg12.5)') fxist
    call lejust(ux16b)
    j2 = ilnobl(ux16b)
    write (noutpt,1070) ux16a(1:j1),ux16b(1:j2)
1070 format(/17x,'Ionic strength (I)= ',a,' molal',/6x,'Stoichiometric ionic strength= ',a,' molal')

    write (ux16a,'(1pg12.5)') fje
    call lejust(ux16a)
    j1 = ilnobl(ux16a)
    write (ux16b,'(1pg12.5)') fjest
    call lejust(ux16b)
    j2 = ilnobl(ux16b)
    write (noutpt,1080) ux16a(1:j1),ux16b(1:j2)
1080 format(/16x,'Ionic asymmetry (J)= ',a,' molal',/5x,'Stoichiometric ionic asymmetry= ',a,' molal')

    if (iopg(1) .eq. 2) then
        write (ux16a,'(1pg12.5)') abar
        call lejust(ux16a)
        j1 = ilnobl(ux16a)
        write (ux16b,'(1pg12.5)') a3bar
        call lejust(ux16b)
        j2 = ilnobl(ux16b)
        write (noutpt,1090) ux16a(1:j1),ux16b(1:j2)
1090 format(/31x,'Abar = ',a,/30x,'A3bar = ',a)
    end if

    write (ux16a,'(1pg12.5)') woh2o
    call lejust(ux16a)
    j1 = ilnobl(ux16a)
    write (ux16b,'(1pg12.5)') wotds
    call lejust(ux16b)
    j2 = ilnobl(ux16b)
    write (ux16c,'(1pg12.5)') wosol
    call lejust(ux16c)
    j3 = ilnobl(ux16c)
    write (noutpt,1100) ux16a(1:j1),ux16b(1:j2),ux16c(1:j3)
1100 format(/23x,'Solvent mass= ',a,' g',/17x,'Solutes (TDS) mass= ',a,' g',/14x,'Aqueous solution mass= ',a,' g')

    if (qrho) then
        write (ux16a,'(1pg12.5)') vosol
        call lejust(ux16a)
        j1 = ilnobl(ux16a)
        write (noutpt,1110) ux16a(1:j1)
1110 format(/12x,'Aqueous solution volume= ',a,' L')
    end if

    write (ux16a,'(1pg12.5)') wfh2o
    call lejust(ux16a)
    j1 = ilnobl(ux16a)
    write (ux16b,'(1pg12.5)') wftds
    call lejust(ux16b)
    j2 = ilnobl(ux16b)
    write (noutpt,1120) ux16a(1:j1),ux16b(1:j2)
1120 format(/19x,'Solvent fraction= ',a,' kg.H2O/kg.sol',/20x,'Solute fraction= ',a,' kg.tds/kg.sol')

    write (ux16a,'(1pg12.5)') tdspkc
    call lejust(ux16a)
    j1 = ilnobl(ux16a)
    write (noutpt,1140) ux16a(1:j1)
1140 format(/6x,'Total dissolved solutes (TDS)= ',a,' mg/kg.sol')

    if (qrho) then
        write (ux16a,'(1pg12.5)') tdsplc
        call lejust(ux16a)
        j1 = ilnobl(ux16a)
        write (ux16b,'(1pg12.5)') tdsglw
        call lejust(ux16b)
        j2 = ilnobl(ux16b)
        write (noutpt,1210) ux16a(1:j1),ux16b(1:j2)
1210 format(32x,'TDS= ',a,' mg/L',/32x,'TDS= ',a,' g/L')

        write (ux16a,'(1pg12.5)') rhoc
        call lejust(ux16a)
        j1 = ilnobl(ux16a)
        write (ux16b,'(1pg12.5)') rhowc
        call lejust(ux16b)
        j2 = ilnobl(ux16b)
        write (noutpt,1220) ux16a(1:j1),ux16b(1:j2)
1220 format(/19x,'Solution density= ',a,' g/mL',/19x,'Solution density= ',a,' g/L')

        write (ux16a,'(1pg12.5)') mrmlra
        call lejust(ux16a)
        j1 = ilnobl(ux16a)
        write (ux16b,'(1pg12.5)') mlmrra
        call lejust(ux16b)
        j2 = ilnobl(ux16b)
        write (noutpt,1240) ux16a(1:j1),ux16b(1:j2)
1240 format(/18x,'Molarity/molality= ',a,' kg.H2O/L',/18x,'Molality/molarity= ',a,' L/kg.H2O')
    end if

    write (noutpt,1300)
1300 format(1x)
end subroutine prtvpa