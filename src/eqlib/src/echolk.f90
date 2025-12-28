subroutine echolk(axlks,cdrs,ilevel,jsflag,narxmx,ndrs,ndrsmx,ndrsr,nf,nst,ntprmx,nstmax,press,tempc,uspec,xlks)
    !! This subroutine prints the species and reactions that are active
    !! in the current problem, along with the log K values that
    !! correspond to the reactions. Optionally, the coefficients
    !! of the interpolating polynomials may also be printed.
    !! This subroutine is called by:
    !!   EQ3NR/eq3nr.f
    !!   EQ6/echoz.f
    !!   EQ6/eq6.f
    !! Principal input:
    !!   ilevel = print level switch
    !!              =  0  no print
    !!              =  1  print species and reactions only
    !!              =  2  also print equilibrium constants
    !!              =  3  also print the coefficients of the
    !!                    interpolating polynomials
    !!   nf     = unit number of the file to write on
    !! Principal output:
    !!   None
    implicit none

    ! Calling sequence variable declarations.
    integer :: narxmx
    integer :: ndrsmx
    integer :: ntprmx
    integer :: nstmax

    integer :: jsflag(nstmax)
    integer :: ndrs(ndrsmx)
    integer :: ndrsr(2,nstmax)
    integer :: ilevel
    integer :: nf
    integer :: nst

    character(len=48) :: uspec(nstmax)

    real(kind=8) :: axlks(narxmx,ntprmx,nstmax)
    real(kind=8) :: cdrs(ndrsmx)
    real(kind=8) :: xlks(nstmax)
    real(kind=8) :: press
    real(kind=8) :: tempc

    ! Local variable declarations.
    integer :: i
    integer :: j
    integer :: ns

    if (ilevel .le. 0) then
        go to 999
    end if

    write (nf,2005)
2005 format(//' --- Listing of Species and Reactions ---',/)

    if (ilevel .ge. 2) then
        write (nf,2040) tempc,press
    end if

2040 format(7x,'temperature= ',f10.3,'C',/7x,'pressure= ',g13.6,' bars',/)

    do ns = 1,nst
        if (jsflag(ns) .le. 0) then
            write (nf,2050)
2050 format(' --------------------------------------------------')

            call prreac(cdrs,ndrs,ndrsmx,ndrsr,nf,ns,nstmax,uspec)

            if (ilevel .ge. 2) then
                write (nf,2060) xlks(ns)
            end if

2060 format(/10x,'log K= ',f12.4,/)

            if (ilevel .ge. 3) then
                do j = 1,ntprmx
                    write (nf,2070) j,(axlks(i,j,ns), i = 1,narxmx)
2070 format(/3x,'Coefficients for range ',i2,':',/3x,3(2x,g13.6),/3x,3(2x,g13.6),/)
                end do
            end if
        end if
    end do

    write (nf,2050)

999 continue
end subroutine echolk