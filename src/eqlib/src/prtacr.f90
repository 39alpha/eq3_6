subroutine prtacr(actlg,iopr,jsflag,nbaspd,nbt,nbtmax,nelect,nhydr,noprmx,no2gaq,noutpt,nstmax,uspec,zchar)
    !! This subroutine prints a table of cation/H+ activity ratios,
    !! anion-H+ activity products, and neutral species activities.
    !! Only aqueous basis species are involved. The level of printing
    !! is controlled by the print control flag iopr(5):
    !!    0 = Don't print
    !!    1 = Print cation/H+ activity ratios only
    !!    2 = Print cation/H+ activity ratios and anion-H+ activity
    !!          products
    !!    3 = Print cation/H+ activity ratios, anion-H+ activity
    !!          products, and neutral species activities
    !! This subroutine is called by:
    !!   EQ3NR/scripx.f
    !!   EQ6/scripz.f
    !! Principal input:
    !!   actlg  = array of log activities of species
    !!   iopr   = array of print control options
    !!   jsflag = array of species status flags
    !!   nbt    = the number of species in the basis set
    !!   nelect = index of the fictive species aqueous e-
    !!   nhydr  = index of the species aqueous H+
    !!   no2gaq = index of the fictive species aqueous O2(g)
    !!   uspec  = array of species names
    !!   zchar  = array of species electrical charge numbers
    !! Principal output:
    !!   None
    implicit none

    ! Calling sequence variable declarations.
    integer :: nbtmax
    integer :: noprmx
    integer :: nstmax

    integer :: noutpt

    integer :: iopr(noprmx)
    integer :: jsflag(nstmax)
    integer :: nbaspd(nbtmax)
    integer :: nbt
    integer :: nelect
    integer :: nhydr
    integer :: no2gaq

    character(len=48) :: uspec(nstmax)

    real(kind=8) :: actlg(nstmax)
    real(kind=8) :: zchar(nstmax)

    ! Local variable declarations.
    integer :: iz
    integer :: j2
    integer :: nb
    integer :: ns

    integer :: ilnobl

    real(kind=8) :: actlh
    real(kind=8) :: actrat
    real(kind=8) :: zx

    if (iopr(5) .le. 0) then
        go to 999
    end if

    if (iopr(5) .eq. 1) then
        write (noutpt,1000)
1000 format(/16x,'--- Cation/H+ Activity Ratios ---',/)
    else if (iopr(5) .eq. 2) then
        write (noutpt,1010)
1010 format(/16x,'--- Ion-H+ Activity Ratios ---',/)
    else if (iopr(5) .eq. 3) then
        write (noutpt,1020)
1020 format(/6x,'--- Ion-H+ Activity Ratios and',' Neutral Species Activities ---',/)
    end if

    actlh = actlg(nhydr)

    do nb = 1,nbt
        ns = nbaspd(nb)

        if (ns .eq. nhydr) then
            go to 100
        end if

        if (ns .eq. no2gaq) then
            go to 100
        end if

        if (ns .eq. nelect) then
            go to 100
        end if

        if (jsflag(ns) .ge. 2) then
            go to 100
        end if

        zx = zchar(ns)
        actrat = actlg(ns)

        if (zx .eq. 0.) then
            if (iopr(5) .ge. 3) then
                j2 = ilnobl(uspec(ns)(1:16))
                write (noutpt,1040) uspec(ns)(1:j2),actrat
1040 format(3x,'Log ( a(',a,') )',t44,'= ',f10.5)
            end if
        else
            actrat = actrat - zx*actlh
            j2 = ilnobl(uspec(ns)(1:16))
            iz = nint(abs(zx))

            if (zx .lt. 0.) then
                if (iopr(5) .ge. 2) then
                    write (noutpt,1050) uspec(ns)(1:j2),iz,actrat
1050 format(3x,'Log ( a(',a,') x a(H+)xx',i2,' )',t44,'= ',f10.5)
                end if
            else
                write (noutpt,1060) uspec(ns)(1:j2),iz,actrat
1060 format(3x,'Log ( a(',a,') / a(H+)xx',i2,' )',t44,'= ',f10.5)
            end if
        end if

100 continue
    end do

    write (noutpt,1070)
1070 format(1x)

999 continue
end subroutine prtacr