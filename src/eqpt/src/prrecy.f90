subroutine prrecy(cdrs,nbtmx1,nbtmx2,nbt,ns,nfile,uspec)
    !! This subroutine writes the reaction associated with the ns-th
    !! species on the file whose unit number is nfile. This subroutine
    !! is virtually identical in function to EQLIB/prreac.f, but
    !! deals with data that is structured somewhat differently.
    !! This subroutine is called by:
    !!   EQPT/rxnchk.f
    !! Principal input:
    !!   cdrs   = array of reaction coefficients
    !!   ns     = species whose reaction is to be printed
    !!   nfile  = unit number of the file to write on
    !!   uspec  = array of species names
    !! Principal input:
    !!   None
    implicit none

    ! Calling sequence variable declarations.
    integer :: nbtmx1
    integer :: nbtmx2

    integer :: nbt
    integer :: nfile
    integer :: ns

    character(len=24) :: uspec(nbtmx1)

    real(kind=8) :: cdrs(nbtmx2,nbtmx1)

    ! Local variable declarations.
    integer :: j2
    integer :: nj
    integer :: nb
    integer :: nbt1

    integer :: ilnobl

    real(kind=8) :: cx

    nbt1 = nbt + 1
    cx = -cdrs(nbt1,ns)
    j2 = ilnobl(uspec(ns))
    write (nfile,1000) cx,uspec(ns)(1:j2)
1000 format(/21x,f10.5,2x,a)

    do nb = 1,nbt
        cx = -cdrs(nb,ns)

        if (cx .gt. 0.) then
            j2 = ilnobl(uspec(nb))
            write (nfile,1010) cx,uspec(nb)(1:j2)
1010 format(18x,' + ',f10.5,2x,a)
        end if
    end do

    write (nfile,1020)
1020 format(/31x,'='/)

    nj = 0

    do nb = 1,nbt
        cx = cdrs(nb,ns)

        if (cx .gt. 0.) then
            j2 = ilnobl(uspec(nb))

            if (nj .le. 0) then
                write (nfile,1030) cx,uspec(nb)(1:j2)
1030 format(21x,f10.5,2x,a)

                nj = 1
            else
                write (nfile,1010) cx,uspec(nb)(1:j2)
            end if
        end if
    end do
end subroutine prrecy