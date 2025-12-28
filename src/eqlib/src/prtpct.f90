subroutine prtpct(conc,csts,ctb,iopr,jcsort,jflag,narn1,narn2,nbasp,nbaspd,nbt,nbtmax,nelect,nhydr,noprmx,no2gaq,noutpt,nstmax,nsts,nstsmx,nstsr,uspec)
    !! This subroutine prints tables giving the percentages of species
    !! making up solute mass totals in the aqueous solution. The
    !! level of printing is controlled by the print control flag
    !! iopr(6):
    !!    0 = Don't print any tables
    !!    1 = Print tables including 99% of all contributing species
    !!    2 = Print tables including all contributing species
    !! This subroutine is called by:
    !!   EQ3NR/scripx.f
    !!   EQ6/scripz.f
    !! Principal input:
    !!   conc   = array of species concentrations
    !!   csts   = array of stoichiometric mass balance factors
    !!   ctb    = array of total molalities of data file basis
    !!              species
    !!   iopr   = array of print control options
    !!   jcsort = array of species indices, in order of increasing
    !!              concentration, but with sorting restricted to within
    !!              phase ranges
    !!   narn1  = start of the range of aqueous species
    !!   narn2  = end of the range of aqueous species
    !!   nbasp  = array of the species indices in the active basis set
    !!   nbaspd = array of the species indices in the data file basis
    !!              set
    !!   nbt    = number of species in the basis set
    !!   nelect = index of the fictive species aqueous e-
    !!   nhydr  = index of the species aqueous H+
    !!   no2gaq = index of the fictive species aqueous O2(g)
    !!   nsts   = array of indices of species corresponding to the
    !!              coefficients in the csts array
    !!   nstsr  = array giving the range in the csts and nsts arrays
    !!              corresponding to a given species
    !!   nwater = index of the species liquid water
    !!   uspec  = array of species names
    !! Principal output:
    !!   None
    implicit none

    integer :: nbtmax
    integer :: noprmx
    integer :: nstmax
    integer :: nstsmx

    integer :: noutpt

    integer :: iopr(noprmx)
    integer :: jcsort(nstmax)
    integer :: jflag(nstmax)
    integer :: nbasp(nbtmax)
    integer :: nbaspd(nbtmax)
    integer :: nsts(nstsmx)
    integer :: nstsr(2,nstmax)
    integer :: narn1
    integer :: narn2
    integer :: nbt
    integer :: nelect
    integer :: nhydr
    integer :: no2gaq

    character(len=48) :: uspec(nstmax)

    real(kind=8) :: conc(nstmax)
    real(kind=8) :: csts(nstsmx)
    real(kind=8) :: ctb(nbtmax)

    ! Local variable declarations.
    integer :: j2
    integer :: nb
    integer :: ns
    integer :: ns1
    integer :: ns2
    integer :: nss
    integer :: nssi

    integer :: ilnobl

    real(kind=8) :: apct
    real(kind=8) :: cx
    real(kind=8) :: cxn
    real(kind=8) :: cxs
    real(kind=8) :: cxt
    real(kind=8) :: fx
    real(kind=8) :: px
    real(kind=8) :: pxs
    real(kind=8) :: pxt
    real(kind=8) :: coefst

    if (iopr(6) .le. -1) then
        go to 999
    end if

    if (iopr(6) .le. 0) then
        write (noutpt,1000)
1000 format(//6x,'--- Major Species by Contribution to Aqueous',' Mass Balances ---',/)
    else
        write (noutpt,1010)
1010 format(//6x,'--- Species by Contribution to Aqueous Mass',' Balances ---',/)
    end if

    do nb = 1,nbt
        ns1 = nbaspd(nb)
        ns2 = nbasp(nb)

        if (ns1.lt.narn1 .or. ns1.gt.narn2) then
            go to 110
        end if

        if (conc(ns2).le.0. .or. jflag(ns2).eq.30) then
            go to 110
        end if

        if (ns1.eq.no2gaq .or. ns1.eq.nelect) then
            go to 110
        end if

        if (ns1.eq.narn1 .or. ns1.eq.nhydr) then
            go to 110
        end if

        j2 = ilnobl(uspec(ns1)(1:24))

        if (iopr(6) .eq. 0) then
            write (noutpt,1020) uspec(ns1)(1:j2)
1020 format(/' Species Accounting for 99% or More of Aqueous ',a,//5x,'Species',19x,'Factor',4x,'Molality',5x,'Per Cent',/)
        else
            write (noutpt,1030) uspec(ns1)(1:j2)
1030 format(/' Species Accounting for Total Aqueous ',a,//5x,'Species',19x,'Factor',4x,'Molality',5x,'Per Cent',/)
        end if

        cxt = ctb(nb)
        pxt = 100.

        if (cxt .gt. 0.) then
            cxn = pxt/cxt
        else
            cxn = 0.
        end if

        cxs = 0.
        pxs = 0.

        do nss = narn1,narn2
            nssi = narn2 - nss + narn1
            ns = jcsort(nssi)

            if (conc(ns) .gt. 0.) then
                fx = coefst(csts,nsts,nstsmx,nstsr,nb,ns,nstmax)

                if (fx .gt. 0.) then
                    cx = fx*conc(ns)
                    cxs = cxs + cx
                    px = cx*cxn
                    pxs = pxs + px
                    apct = abs(px)

                    if (apct .gt. 999.99) then
                        write (noutpt,1040) uspec(ns),fx,conc(ns),px
1040 format(3x,a24,3x,f6.2,3x,1pe11.4,3x,e10.3)
                    else if (apct .lt. 1.00) then
                        write (noutpt,1040) uspec(ns),fx,conc(ns),px
                    else
                        write (noutpt,1050) uspec(ns),fx,conc(ns),px
1050 format(3x,a24,3x,f6.2,3x,1pe11.4,3x,0pf7.2)
                    end if

                    if (iopr(6).le.0 .and. pxs.ge.99.) then
                        go to 105
                    end if
                end if
            end if
        end do

105 continue
        write (noutpt,1060)
1060 format(' - - - - - - - - - - - - - - - - - - - - - - - - - - -',' - - - - -')

        if (iopr(6) .eq. 0) then
            write (noutpt,1070) cxs,pxs
1070 format(3x,'Subtotal',28x,1pe11.4,3x,0pf7.2,/)
        else
            write (noutpt,1080) cxt,pxt
1080 format(3x,'Total',31x,1pe11.4,3x,0pf7.2,/)
        end if

110 continue
    end do

999 continue
end subroutine prtpct