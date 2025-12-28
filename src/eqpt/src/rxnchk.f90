subroutine rxnchk(cdrs,cess,mtotr,nbt,nbtmx1,nbtmx2,nco,nct,nctmax,nerr,noutpt,ns,nsb,nttyo,uelem,uspec,zchar)
    !! This suboutine checks the reaction associated with the ns-th
    !! species for mass and charge balance. If an imbalance is found,
    !! a message is written to the screen and output files.
    !! This suboutine is called by:
    !!   EQPT/pcraq.f
    !!   EQPT/pcrsg.f
    !! Principal input:
    !!   cdrs   = array of reaction coefficients
    !!   cess   = array of elemental composition coefficients
    !!   nbt    = the number of basis species
    !!   nct    = the number of chemical elements
    !!   nerr   = cumulative error counter
    !!   uelem  = array of chemical element names
    !!   uspec  = array of species names
    !!   zchar  = array of electrical charge numbers
    !! Principal output:
    !!   nerr   = cumulative error counter
    !! Workspace:
    !!   mtotr  = array of mass balance residuals
    implicit none

    ! Calling sequence variable declarations.
    integer :: nbtmx1
    integer :: nbtmx2

    integer :: noutpt
    integer :: nttyo

    integer :: nbt
    integer :: nco
    integer :: nct
    integer :: nctmax
    integer :: nerr
    integer :: ns
    integer :: nsb

    character(len=24) :: uspec(nbtmx1)
    character(len=8) :: uelem(nctmax)
    character(len=8) :: ux8

    real(kind=8) :: cdrs(nbtmx2,nbtmx1)
    real(kind=8) :: cess(nctmax,nbtmx1)
    real(kind=8) :: mtotr(nctmax)
    real(kind=8) :: zchar(nbtmx1)

    ! Local variable declarations.
    logical :: qbadr
    logical :: qo2gx

    integer :: j2
    integer :: j3
    integer :: nbt1
    integer :: nc
    integer :: nfile
    integer :: nlim
    integer :: nse

    integer :: ilnobl

    real(kind=8) :: mx
    real(kind=8) :: rxntol
    real(kind=8) :: ztotr

    data rxntol /1.e-6/

    qbadr = .false.

    qo2gx = uspec(nsb) .eq. 'O2(g)'

    if (qo2gx) then
        cess(nco,nsb) = 2.
    end if

    nlim = ns - 1
    nlim = min(nlim,nbt)
    nbt1 = nbt + 1

    ztotr = cdrs(nbt1,ns)*zchar(ns)

    do nse = 1,nlim
        ztotr = ztotr + cdrs(nse,ns)*zchar(nse)
    end do

    if (abs(ztotr) .gt. rxntol) then
        qbadr = .true.
    end if

    do nc = 1,nct
        mx = cdrs(nbt1,ns)*cess(nc,ns)

        do nse = 1,nlim
            mx = mx + cdrs(nse,ns)*cess(nc,nse)
        end do

        mtotr(nc) = mx

        if (abs(mx) .gt. rxntol) then
            qbadr = .true.
        end if
    end do

    if (qo2gx) then
        cess(nco,nsb) = 0.
    end if

    if (qbadr) then
        ! Have a reaction with one or more imbalances.
        j2 = ilnobl(uspec(ns))
        write (noutpt,1000) uspec(ns)(1:j2)
        write (nttyo,1000) uspec(ns)(1:j2)
1000 format(/' * Error - (EQPT/rxnchk) The reaction for the',' destruction of',/7x,a,' has the following imbalances:',/)

        if (abs(ztotr) .gt. rxntol) then
            ux8 = 'charge'
            j3 = ilnobl(ux8)
            write (noutpt,1010) ztotr,ux8(1:j3)
            write (nttyo,1010) ztotr,ux8(1:j3)
1010 format(9x,1pe12.5,' in ',a)
        end if

        do nc = 1,nct
            if (abs(mtotr(nc)) .gt. rxntol) then
                ux8 = uelem(nc)
                j3 = ilnobl(ux8)
                write (noutpt,1010) mtotr(nc),ux8(1:j3)
                write (nttyo,1010) mtotr(nc),ux8(1:j3)
            end if
        end do

        write (noutpt,1020)
        write (nttyo,1020)
1020 format(/9x,'The reaction is:')

        nfile = noutpt
        call prrecy(cdrs,nbtmx1,nbtmx2,nbt,ns,nfile,uspec)
        nfile = nttyo
        call prrecy(cdrs,nbtmx1,nbtmx2,nbt,ns,nfile,uspec)
        nerr = nerr + 1
    end if
end subroutine rxnchk