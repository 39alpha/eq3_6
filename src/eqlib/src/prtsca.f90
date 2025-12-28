subroutine prtsca(ctb,jflag,jsflag,mrmlra,mwtsp,narn1,narn2,nbasp,nbaspd,nbt,nbtmax,nelect,nhydr,nhydx,noutpt,no2gaq,nstmax,qrho,rho,uspec,wfh2o)
    !! This subroutine computes and prints a table of the sensible
    !! composition of the aqueous solution in terms of mass balance
    !! totals for component (basis) species that generally correspond to
    !! analyzable solutes. This set of basis species is generally not
    !! identical to either the active basis set or the 'd' basis set,
    !! apart from the fact that it doesn't include the solvent, water.
    !! This subroutine is called by:
    !!   EQ3NR/scripx.f
    !!   EQ6/scripz.f
    !! Principal input:
    !!   ctb    = array of molalities of the component (basis) species
    !!   qrho   = .true. if a solution density value is available
    !!   rho    = the density of the aqueous solution
    !!   uspec  = array of names of aqueous species
    !! Principal output:
    !!   None
    implicit none

    ! Calling sequence variable declarations.
    integer :: nbtmax
    integer :: nstmax

    integer :: noutpt

    integer :: jflag(nstmax)
    integer :: jsflag(nstmax)
    integer :: nbasp(nbtmax)
    integer :: nbaspd(nbtmax)
    integer :: narn1
    integer :: narn2
    integer :: nbt
    integer :: nelect
    integer :: nhydr
    integer :: nhydx
    integer :: no2gaq

    logical :: qrho

    character(len=48) :: uspec(nstmax)

    real(kind=8) :: ctb(nbtmax)
    real(kind=8) :: mwtsp(nstmax)
    real(kind=8) :: mrmlra
    real(kind=8) :: rho
    real(kind=8) :: wfh2o

    ! Local variable declarations.
    integer :: nb
    integer :: ns1
    integer :: ns2

    real(kind=8) :: ctbh
    real(kind=8) :: ctboh
    real(kind=8) :: molrx
    real(kind=8) :: ppmv
    real(kind=8) :: ppmw

    write (noutpt,1000)
1000 format(/11x,'--- Sensible Composition of the',' Aqueous Solution ---')

    if (qrho) then
        ! Have density data, include mg/L and Molarity.
        write (noutpt,1010)
1010 format(/3x,'Species',20x,'mg/L',7x,'mg/kg.sol',4x,'Molarity',5x,'Molality',/)

        do nb = 1,nbt
            ns1 = nbaspd(nb)
            ns2 = nbasp(nb)

            if (ns1.ge.narn1 .and. ns1.le.narn2) then
                if (jsflag(ns2).lt.2 .and. jflag(ns2).ne.30) then
                    if (ns1.ne.narn1 .and. ns1.ne.nelect .and.          ns1.ne.no2gaq) then
                        if (ns1.eq.nhydr .or. ns1.eq.nhydx) then
                            if (ctb(nb) .ge. 0.) then
                                ppmw = 1000.*ctb(nb)*mwtsp(ns1)*wfh2o
                                ppmv = ppmw*rho
                                molrx = ctb(nb)*mrmlra
                                write (noutpt,1020) uspec(ns1),ppmv,ppmw,molrx,ctb(nb)
1020 format(1x,a24,1x,1pe12.5,1x,1pe12.5,1x,1pe12.5,1x,1pe12.5)
                            else if (ns1 .eq. nhydr) then
                                ctboh = -ctb(nb)
                                ppmw = 1000.*ctboh*mwtsp(nhydx)*wfh2o
                                ppmv = ppmw*rho
                                molrx = ctb(nb)*mrmlra
                                write (noutpt,1020) uspec(nhydx),ppmv,ppmw,molrx,ctboh
                            else
                                ctbh = -ctb(nb)
                                ppmw = 1000.*ctbh*mwtsp(nhydr)*wfh2o
                                ppmv = ppmw*rho
                                molrx = ctb(nb)*mrmlra
                                write (noutpt,1020) uspec(nhydr),ppmv,ppmw,molrx,ctbh
                            end if
                        else
                            ppmw = 1000.*ctb(nb)*mwtsp(ns1)*wfh2o
                            ppmv = ppmw*rho
                            molrx = ctb(nb)*mrmlra
                            write (noutpt,1020) uspec(ns1),ppmv,ppmw,molrx,ctb(nb)
                        end if
                    end if
                end if
            end if
        end do
    else
        ! Have no density data, so do not include mg/L and Molarity.
        write (noutpt,1030)
1030 format(/3x,'Species',18x,'mg/kg.sol',4x,'Molality',/)

        do nb = 1,nbt
            ns1 = nbaspd(nb)
            ns2 = nbasp(nb)

            if (ns1.ge.narn1 .and. ns1.le.narn2) then
                if (jsflag(ns2).lt.2 .and. jflag(ns2).ne.30) then
                    if (ns1.ne.narn1 .and. ns1.ne.nelect .and.          ns1.ne.no2gaq) then
                        if (ns1.eq.nhydr .or. ns1.eq.nhydx) then
                            if (ctb(nb) .ge. 0.) then
                                ppmw = 1000.*ctb(nb)*mwtsp(ns1)*wfh2o
                                write (noutpt,1040) uspec(ns1),ppmw,ctb(nb)
1040 format(1x,a24,1x,1pe12.5,1x,1pe12.5)
                            else if (ns1 .eq. nhydr) then
                                ctboh = -ctb(nb)
                                ppmw = 1000.*ctboh*mwtsp(nhydx)*wfh2o
                                write (noutpt,1040) uspec(nhydx),ppmw,ctboh
                            else
                                ctbh = -ctb(nb)
                                ppmw = 1000.*ctbh*mwtsp(nhydr)*wfh2o
                                write (noutpt,1040) uspec(nhydr),ppmw,ctbh
                            end if
                        else
                            ppmw = 1000.*ctb(nb)*mwtsp(ns1)*wfh2o
                            write (noutpt,1040) uspec(ns1),ppmw,ctb(nb)
                        end if
                    end if
                end if
            end if
        end do
    end if

    write (noutpt,1050)
1050 format(/3x,'The above data have physical significance, but',' some may be',/3x,'inconsistent with certain analytical',' methods or reporting schemes.',/)
end subroutine prtsca