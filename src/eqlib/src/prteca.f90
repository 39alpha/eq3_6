subroutine prteca(cteaq,mrmlra,nct,nctmax,noutpt,ppmwe,qrho,rho,uelem)
    !! This subroutine prints a table of the elemental composition of the
    !! aqueous solution.
    !! This subroutine is called by:
    !!   EQ3NR/scripx.f
    !!   EQ6/scripz.f
    !! Principal input:
    !!   cteaq  = array of molalities of the chemical elements
    !!   ppmwe  = array of ppms by weight of the chemical elements
    !!   rho    = the density of the aqueous solution
    !!   uelem  = array of names of chemical elements
    !! Principal output:
    !!   None
    implicit none

    ! Calling sequence variable declarations.
    integer :: nctmax

    integer :: noutpt

    integer :: nct

    logical :: qrho

    character(len=8) :: uelem(nctmax)

    real(kind=8) :: cteaq(nctmax)
    real(kind=8) :: ppmwe(nctmax)
    real(kind=8) :: mrmlra
    real(kind=8) :: rho

    ! Local variable declarations.
    integer :: nc

    real(kind=8) :: molrex
    real(kind=8) :: ppmve

    write (noutpt,1000)
1000 format(/11x,'--- Elemental Composition of the Aqueous Solution',' ---')

    if (qrho) then
        ! Have density data, include mg/L and Molarity.
        write (noutpt,1010)
1010 format(/3x,'Element',8x,'mg/L',7x,'mg/kg.sol',4x,'Molarity',5x,'Molality',/)

        do nc = 1,nct
            ppmve = ppmwe(nc)*rho
            molrex = cteaq(nc)*mrmlra
            write (noutpt,1020) uelem(nc),ppmve,ppmwe(nc),molrex,cteaq(nc)
1020 format(5x,a8,1x,1pe12.5,1x,1pe12.5,1x,1pe12.5,1x,1pe12.5)
        end do
    else
        ! Have no density data, so do not include mg/L and Molarity.
        write (noutpt,1030)
1030 format(/3x,'Element',6x,'mg/kg.sol',4x,'Molality',/)

        do nc = 1,nct
            write (noutpt,1040) uelem(nc),ppmwe(nc),cteaq(nc)
1040 format(5x,a8,1x,1pe12.5,1x,1pe12.5)
        end do
    end if

    write (noutpt,1050)
1050 format(1x)
end subroutine prteca