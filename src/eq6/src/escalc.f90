subroutine escalc(csts,iindx1,jcsort,kbt,kmax,moph,mosp,mtb,mtb0,nbaspd,nbt,nbtmax,ncmpr,noutpt,npt,nptmax,nstmax,nsts,nstsmx,nstsr,qprflg,uspec)
    !! This subroutine recomputes the mass balance totals for the
    !! Equilibrium System (ES). Normally, this is done immediately
    !! after a shift of material from the ES to the Physically Removed
    !! System (PRS). The calculation is performed by summing over what
    !! remains in the ES, so it is not sensitive to the nature of
    !! the material that was transferred.
    !! This subroutine is called by:
    !!   EQ6/dumpdp.f
    !!   EQ6/pshfta.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: kmax
    integer :: nbtmax
    integer :: nptmax
    integer :: nstmax
    integer :: nstsmx

    integer :: noutpt

    integer :: iindx1(kmax)
    integer :: jcsort(nstmax)
    integer :: nbaspd(nbtmax)
    integer :: ncmpr(2,nptmax)
    integer :: nsts(nstsmx)
    integer :: nstsr(2,nstmax)

    integer :: kbt
    integer :: nbt
    integer :: npt

    logical :: qprflg

    character(len=48) :: uspec(nstmax)

    real(kind=8) :: csts(nstsmx)
    real(kind=8) :: moph(nptmax)
    real(kind=8) :: mosp(nstmax)
    real(kind=8) :: mtb(nbtmax)
    real(kind=8) :: mtb0(nbtmax)

    ! Local variable declarations.
    integer :: kcol
    integer :: n
    integer :: nb
    integer :: np
    integer :: nrn1
    integer :: nrn2
    integer :: nr1
    integer :: nr2
    integer :: ns
    integer :: nss

    call copyaa(mtb,mtb0,nbt)

    do nb = 1,nbt
        mtb(nb) = 0.
    end do

    do np = 1,npt
        if (moph(np) .gt. 0.) then
            nrn1 = ncmpr(1,np)
            nrn2 = ncmpr(2,np)

            do nss = nrn1,nrn2
                ns = jcsort(nss)

                if (mosp(ns) .gt. 0.) then
                    nr1 = nstsr(1,ns)
                    nr2 = nstsr(2,ns)

                    do n = nr1,nr2
                        nb = nsts(n)
                        mtb(nb) = mtb(nb) + csts(n)*mosp(ns)
                    end do
                end if
            end do
        end if
    end do

    if (qprflg) then
        write (noutpt,1000)
1000 format(/'   --- Mass Balance Totals ---',//3x,'Species',20x,'Old',10x,'New',/)

        do kcol = 1,kbt
            nb = iindx1(kcol)
            ns = nbaspd(nb)
            write (noutpt,1010) uspec(ns),mtb0(nb),mtb(nb)
1010 format(1x,a24,2(2x,1pe12.5))
        end do

        write (noutpt,1020)
1020 format(1x)
    end if

    call copyaa(mtb,mtb0,nbt)
end subroutine escalc