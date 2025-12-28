subroutine prtrct(afrc1,aft1,imchmx,imech,iopt,modr,morr,noptmx,noutpt,nrct,nrctmx,nrk,rk,rreacn,rreac1,rrelr1,sfcar,ureac,wodr,wodrt,worr,worrt,xi1,xistsv)
    !! This subroutine prints tables of data for irreversible reactions
    !! and corresponding reaction rates.
    !! This subroutine is called by:
    !!   EQ6/scripz.f
    !! Principal input:
    !! Principal output:
    !!   None
    implicit none

    ! Calling sequence variable declarations.
    integer :: imchmx
    integer :: noptmx
    integer :: nrctmx

    integer :: noutpt

    integer :: imech(2,nrctmx)
    integer :: iopt(noptmx)
    integer :: nrk(2,nrctmx)

    integer :: nrct

    character(len=24) :: ureac(nrctmx)

    real(kind=8) :: afrc1(nrctmx)
    real(kind=8) :: modr(nrctmx)
    real(kind=8) :: morr(nrctmx)
    real(kind=8) :: rk(imchmx,2,nrctmx)
    real(kind=8) :: rreacn(nrctmx)
    real(kind=8) :: rreac1(nrctmx)
    real(kind=8) :: rrelr1(nrctmx)
    real(kind=8) :: sfcar(nrctmx)
    real(kind=8) :: wodr(nrctmx)
    real(kind=8) :: worr(nrctmx)

    real(kind=8) :: aft1
    real(kind=8) :: wodrt
    real(kind=8) :: worrt
    real(kind=8) :: xi1
    real(kind=8) :: xistsv

    ! Local variable declarations.
    integer :: i
    integer :: j2
    integer :: nrc

    integer :: ilnobl

    character(len=16) :: ux16

    write (noutpt,1000)
1000 format(//21x,'--- Reactant Summary ---',//)

    if (xi1 .le. xistsv) then
        write (noutpt,1010)
    end if

1010 format(5x,'Definitions and conventions',//10x,'Delta x = x now - x at start',/10x,'Affinity is + for forward direction (destruction),',/10x,'            - for reverse direction (formation)',/10x,'Rates are + for forward direction (destruction),',/10x,'          - for reverse direction (formation)',/)

    write (noutpt,1020)
1020 format(/5x,'Reactant',18x,'Moles',5x,'Delta moles',4x,'Mass, g',4x,'Delta mass, g',/)

    do nrc = 1,nrct
        write (noutpt,1030) ureac(nrc),morr(nrc),modr(nrc),worr(nrc),wodr(nrc)
1030 format(2x,a24,4(2x,1pe11.4))
    end do

    if (nrct .le. 0) then
        write (noutpt,1040)
    end if

1040 format(2x,'None')

    if (nrct .gt. 0) then
        write (noutpt,1050) worrt,wodrt
1050 format(//20x,'Mass remaining= ',1pe11.4,' grams',/20x,'Mass destroyed= ',e11.4,' grams')
    end if

    if (nrct .gt. 0) then
        if (iopt(2) .le. 0) then
            write (noutpt,1070)
1070 format(//5x,'Reactant',19x,'Affinity',5x,'Rel. Rate',/32x,'kcal/mol',6x,'mol/mol',/)

            do nrc = 1,nrct
                write (noutpt,1090) ureac(nrc),afrc1(nrc),rrelr1(nrc)
1090 format(2x,a24,2x,f13.4,2x,1pe11.4)
            end do
        else
            write (noutpt,1100)
1100 format(//5x,'Reactant',17x,'Rel. Rate',6x,'Rate',9x,'Rate',/31x,'mol/mol',6x,'mol/s',6x,'mol/s/cm2',/)

            do nrc = 1,nrct
                write (noutpt,1110) ureac(nrc),rrelr1(nrc),rreac1(nrc),rreacn(nrc)
1110 format(2x,a24,4(2x,1pe11.4))
            end do

            write (noutpt,1120)
1120 format(//5x,'Reactant',19x,'Affinity',3x,'Surface Area',/32x,'kcal/mol',7x,'cm2',/)

            do nrc = 1,nrct
                write (noutpt,1090) ureac(nrc),afrc1(nrc),sfcar(nrc)
            end do

            write (noutpt,1130)
1130 format(//5x,'Reactant',16x,'Rate Constants, mol/s/cm2',/)

            do nrc = 1,nrct
                write (noutpt,1090) ureac(nrc)

                if (nrk(1,nrc) .ge. 1) then
                    write (noutpt,1140) (rk(i,1,nrc), i = 1,imech(1,nrc))
1140 format(20x,'Forward ',4(2x,1pe11.4),/28x,4(2x,e11.4))
                end if

                if (nrk(2,nrc) .ge. 1) then
                    write (noutpt,1150) (rk(i,2,nrc), i = 1,imech(2,nrc))
1150 format(20x,'Backward',4(2x,1pe11.4),/28x,4(2x,e11.4))
                end if
            end do
        end if

        write (ux16,'(f13.4)') aft1
        call lejust(ux16)
        j2 = ilnobl(ux16)
        write (noutpt,1160) ux16(1:j2)
1160 format(//3x,'Affinity of the overall irreversible reaction= ',a,' kcal.',/3x,'Contributions from irreversible reactions',' with no thermodynamic data',/3x,'are not included.',/)
    end if
end subroutine prtrct