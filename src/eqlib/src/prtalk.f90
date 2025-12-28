subroutine prtalk(alki,alk1,alk2,mrmlra,noutpt,ntf1t,ntf2t,qrho,rho,tempc,wfh2o)
    !! This subroutine writes a table of computed alkalinity parameters.
    !! This subroutine is called by:
    !!   EQ3NR/scripx.f
    !!   EQ6/scripz.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: noutpt

    integer :: ntf1t
    integer :: ntf2t

    real(kind=8) :: alki
    real(kind=8) :: alk1
    real(kind=8) :: alk2
    real(kind=8) :: mrmlra
    real(kind=8) :: rho
    real(kind=8) :: tempc
    real(kind=8) :: wfh2o

    ! Local variable declarations.
    integer :: j2

    integer :: ilnobl

    logical :: qrho

    character(len=16) :: ux16

    real(kind=8) :: alkicv
    real(kind=8) :: alkicw
    real(kind=8) :: alkihv
    real(kind=8) :: alkihw
    real(kind=8) :: alkipl
    real(kind=8) :: alk1cv
    real(kind=8) :: alk1cw
    real(kind=8) :: alk1hv
    real(kind=8) :: alk1hw
    real(kind=8) :: alk1pl
    real(kind=8) :: alk2cv
    real(kind=8) :: alk2cw
    real(kind=8) :: alk2hv
    real(kind=8) :: alk2hw
    real(kind=8) :: alk2pl

    ! Check the temperature. Alkalinity of any kind is only defined
    ! in the temperature range 0-50 C.
    if (tempc.lt.0. .or. tempc.gt.50.) then
        write (ux16,'(f12.3)') tempc
        call lejust(ux16)
        j2 = ilnobl(ux16)
        write (noutpt,1000) ux16(1:j2)
1000 format(//6x,'Alkalinity is not defined at ',a,' C.',/6x,'It is only defined in the temperature range 0-50 C.',/)

        go to 999
    end if

    ! Write a table of the HCO3-CO3-OH total alkalinity.
    write (noutpt,1010)
1010 format(/11x,'--- HCO3-CO3-OH Total Alkalinity ---',/)

    ! Test for the required coefficients.
    if (ntf1t .le. 0) then
        write (noutpt,1020)
1020 format(6x,'Coefficients to compute this alkalinity',' are not available.',/)

        go to 110
    end if

    alk1cw = 50000*wfh2o*alk1
    alk1hw = 1.2192*alk1cw

    if (qrho) then
        ! Have density data, include volumetric measures.
        alk1pl = alk1*mrmlra
        alk1cv = alk1cw*rho
        alk1hv = 1.2192*alk1cv

        write (noutpt,1030) alk1,alk1pl,alk1cw,alk1hw,alk1cv,alk1hv
1030 format(16x,1pg12.5,' eq/kg.H2O',/16x,1pg12.5,' eq/L',/16x,1pg12.5,' mg/kg.sol CaCO3',/16x,1pg12.5,' mg/kg.sol HCO3-',/16x,1pg12.5,' mg/L CaCO3',/16x,1pg12.5,' mg/L HCO3-',/)
    else
        ! Have no density data, do not include volumetric measures.
        write (noutpt,1040) alk1,alk1cw,alk1hw
1040 format(16x,1pg12.5,' eq/kg.H2O',/16x,1pg12.5,' mg/kg.sol CaCO3',/16x,1pg12.5,' mg/kg.sol HCO3-',/)
    end if

110 continue

    ! Write a table of the extended total alkalinity.
    write (noutpt,1050)
1050 format(/11x,'--- Extended Total Alkalinity ---',/)

    ! Test for the required coefficients.
    if (ntf2t .le. 0) then
        write (noutpt,1020)
        go to 120
    end if

    alk2cw = 50000*wfh2o*alk2
    alk2hw = 1.2192*alk2cw

    if (qrho) then
        ! Have density data, include volumetric measures.
        alk2pl = alk2*mrmlra
        alk2cv = alk2cw*rho
        alk2hv = 1.2192*alk2cv

        write (noutpt,1030) alk2,alk2pl,alk2cw,alk2hw,alk2cv,alk2hv
    else
        ! Have no density data, do not include volumetric measures.
        write (noutpt,1040) alk2,alk2cw,alk2hw
    end if

120 continue
    if (alki .gt. 0.) then
        ! Write a table of the input alkalinity. This is relevant in
        ! EQ3NR, but not in EQ6.
        write (noutpt,1060)
1060 format(/11x,'--- Input Total Alkalinity ---',/)

        alkicw = 50000*wfh2o*alki
        alkihw = 1.2192*alkicw

        if (qrho) then
            ! Have density data, include volumetric measures.
            alkipl = alki*mrmlra
            alkicv = alkicw*rho
            alkihv = 1.2192*alkicv
            write (noutpt,1030) alki,alkipl,alkicw,alkihw,alkicv,alkihv
        else
            ! Have no density data, do not include volumetric measures.
            write (noutpt,1040) alki,alkicw,alkihw
        end if
    end if

999 continue
end subroutine prtalk