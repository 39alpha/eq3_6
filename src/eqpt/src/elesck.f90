subroutine elesck(cessi,nbtmx1,nctmax,ncts,nentei,nerr,noutpt,ns,nttyo,qblkes,qzeres,uessi,usblkf,uspec)
    !! This subroutine conducts tests on the elemental composition
    !! specified for a species. It detects any blank and duplicate
    !! element names in the composition and any zero-valued
    !! stoichiometric coefficients.
    !! Note: any blank names are replaced by the string '<blank>'.
    !! This subroutine is called by:
    !!   EQPT/pcraq.f
    !!   EQPT/pcrsg.f
    implicit none

    ! Calling sequence variable declarations.
    integer :: nbtmx1
    integer :: nctmax

    integer :: noutpt
    integer :: nttyo

    integer :: nentei(nctmax)

    integer :: ncts
    integer :: nerr
    integer :: ns

    logical :: qblkes
    logical :: qzeres

    character(len=24) :: uspec(nbtmx1)
    character(len=24) :: usblkf
    character(len=8) :: uessi(nctmax)

    real(kind=8) :: cessi(nctmax)

    ! Calling sequence variable declarations.
    integer :: j2
    integer :: j3
    integer :: j4
    integer :: n
    integer :: ncount

    integer :: ilnobl

    logical :: qdupes

    character(len=8) :: ux8

    ! Check for blank element names.
    qblkes = .false.
    ncount = 0

    do n = 1,ncts
        j3 = ilnobl(uessi(n))

        if (j3 .le. 0) then
            ncount = ncount + 1
            uessi(n) = '<blank>'
        end if
    end do

    if (ncount .gt. 0) then
        j2 = ilnobl(uspec(ns))
        j4 = ilnobl(usblkf)
        write (noutpt,1000) uspec(ns)(1:j2),usblkf(1:j4)
        write (nttyo,1000) uspec(ns)(1:j2),usblkf(1:j4)
1000 format(/' * Error - (EQPT/elesck) The species ',a,' appearing',/7x,'on the data file in the ',a,' superblock',' has a specified',/7x,'composition with blank chemical',' element names',/7x,'in the following positions:',/)

        do n = 1,ncts
            if (uessi(n)(1:7) .eq. '<blank>') then
                write (ux8,'(i5)') n
                call lejust(ux8)
                j4 = ilnobl(ux8)
                write (noutpt,1010) ux8(1:j4)
                write (nttyo,1010) ux8(1:j4)
1010 format(9x,a)
            end if
        end do

        qblkes = .true.
        nerr = nerr + 1
    end if

    ! Check for a chemical element appearing more than once in the
    ! composition.
    call nelcck(nctmax,ncts,nentei,nerr,qdupes,uessi)

    if (qdupes)  then
        j2 = ilnobl(uspec(ns))
        j4 = ilnobl(usblkf)
        write (noutpt,1040) uspec(ns)(1:j2),usblkf(1:j4)
        write (nttyo,1040) uspec(ns)(1:j2),usblkf(1:j4)
1040 format(/' * Error - (EQPT/elesck) The species ',a,' appearing',/7x,'on the data file in the ',a,' superblock',' has a specified',/7x,'elemental composition for which',' there is more than one entry',/7x,'for the following',' element(s):',/)

        do n = 1,ncts
            if (nentei(n) .gt. 1) then
                write (ux8,'(i5)') nentei(n)
                call lejust(ux8)
                j3 = ilnobl(ux8)
                j4 = ilnobl(uessi(n))
                write (noutpt,1050) uessi(n)(1:j4),ux8(1:j3)
                write (nttyo,1050) uessi(n)(1:j4),ux8(1:j3)
1050 format(9x,a,' (',a,' entries)')
            end if
        end do

        nerr = nerr + 1
    end if

    ! Check for zero-valued stoichiometric coefficients.
    qzeres = .false.

    do n = 1,ncts
        if (cessi(n) .eq. 0.) then
            qzeres = .true.
        end if
    end do

    if (qzeres) then
        j2 = ilnobl(uspec(ns))
        j4 = ilnobl(usblkf)
        write (noutpt,1080) uspec(ns)(1:j2),usblkf(1:j4)
        write (nttyo,1080) uspec(ns)(1:j2),usblkf(1:j4)
1080 format(/' * Error - (EQPT/elesck) The species ',a,' appearing',/7x,'on the data file in the ',a,' superblock',' has an elemental',/7x,'composition with zero-valued',' coefficients for the following elements:',/)

        do n = 1,ncts
            if (cessi(n) .eq. 0.) then
                ux8 = uessi(n)
                j3 = ilnobl(ux8)
                write (noutpt,1090) ux8(1:j3)
                write (nttyo,1090) ux8(1:j3)
1090 format(9x,a)
            end if
        end do

        nerr = nerr + 1
    end if
end subroutine elesck
