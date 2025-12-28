subroutine prtaqs(acflg,actlg,conc,conclg,iopr,jcsort,narn1,narn2,noprmx,noutpt,nstmax,uspec)
    !! This subroutine prints a table of the concentrations, activities,
    !! and activity coefficients of the aqueous solute species. The
    !! species are listed in decreasing order of concentration. The
    !! level of printing is controlled by the print control flag
    !! iopr(4):
    !!   -3  = Omit species with molalities < 1.e-8
    !!   -2 =  Omit species with molalities < 1.e-12
    !!   -1 =  Omit species with molalities < 1.e-20
    !!    0 =  Omit species with molalities < 1.e-100
    !!    1 =  Include all species
    !! This subroutine is called by:
    !!   EQ3NR/scripx.f
    !!   EQ6/scripz.f
    !! Principal input:
    !!   acflg  = array of species activity coefficients
    !!   act    = array of species activities
    !!   actlg  = array of species log activities
    !!   conc   = array of species concentrations
    !!   conclg = array of species log concentrations
    !!   iopr   = array of print control options
    !!   jcsort = array of species indices, in order of increasing
    !!              concentration, but with sorting restricted to within
    !!              phase ranges
    !!   narn1  = start of the range of aqueous species
    !!   narn2  = end of the range of aqueous species
    !!   uspec  = array of names of species
    !! Principal output:
    !!   None
    implicit none

    ! Calling sequence variable declarations.
    integer :: noprmx
    integer :: nstmax

    integer :: noutpt

    integer :: iopr(noprmx)
    integer :: jcsort(nstmax)
    integer :: narn1
    integer :: narn2

    character(len=48) :: uspec(nstmax)

    real(kind=8) :: acflg(nstmax)
    real(kind=8) :: actlg(nstmax)
    real(kind=8) :: conc(nstmax)
    real(kind=8) :: conclg(nstmax)

    real(kind=8) :: texp

    ! Local variable declarations.
    integer :: ipr
    integer :: n
    integer :: ns
    integer :: nss

    real(kind=8) :: ltrunc(-3:0)
    real(kind=8) :: mx

    data (ltrunc(n), n = -3,0)/-8.,-12.,-20.,-100./

    ipr = iopr(4)

    write (noutpt,1000)
1000 format(//16x,'--- Distribution of Aqueous Solute Species ---',/)

    write (noutpt,1010)
1010 format(4x,'Species',18x,'Molality',4x,'Log Molality',3x,'Log Gamma',2x,'Log Activity',/)

    if (iopr(4) .ge. 1) then
        do nss = narn1 + 1,narn2
            n = narn2 - nss + narn1
            ns = jcsort(n)
            write (noutpt,1020) uspec(ns),conc(ns),conclg(ns),acflg(ns),actlg(ns)
1020 format(1x,a24,2x,1pe11.4,2x,0pf11.4,2x,f11.4,2x,f11.4)
        end do
    else
        do nss = narn1 + 1,narn2
            n = narn2 - nss + narn1
            ns = jcsort(n)

            if (conclg(ns) .lt. ltrunc(ipr)) then
                go to 110
            end if

            write (noutpt,1020) uspec(ns),conc(ns),conclg(ns),acflg(ns),actlg(ns)
        end do

        go to 999
110 continue
        mx = texp(ltrunc(ipr))
        write (noutpt,1030) mx
1030 format(/4x,'Species with molalities less than ',1pe10.3,' are',' not listed.')
    end if

    write (noutpt,1040)
1040 format(1x)

999 continue
end subroutine prtaqs