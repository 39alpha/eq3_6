subroutine gakey(ndat0s,noutpt,nttyo,uakey)
    !! This suboutine scans the DATA0 file to determine the aqueous
    !! species activity coefficient model (e.g., Pitzer, Simple
    !! Extended Debye-Huckel) associated with this file.
    !! This suboutine is called by:
    !!   EQPT/eqpt.f
    !! Principal input:
    !!   ndat0s = unit number of the stripped DATA0 file
    !! Principal output:
    !!   uakey  = the type of data file ("SEDH" or "Pitzer")
    implicit none

    ! Calling sequence variable declarations.
    integer :: ndat0s
    integer :: noutpt
    integer :: nttyo

    character(len=8) :: uakey

    ! Local variable declarations.
    integer :: i
    integer :: icount
    integer :: n

    character(len=80) :: uline

    ! Determine whether or not this is a Pitzer-type data file.
    ! Look for patterns indicating the presence of Pitzer interaction
    ! coefficients.
    icount = 0

    do n = 1,500
        read (ndat0s,1000,end=100,err=995) uline
1000 format(a)

        i = index(uline,'beta0 =')
        icount = icount + i
        i = index(uline,'beta1 =')
        icount = icount + i
        i = index(uline,'beta2 =')
        icount = icount + i
        i = index(uline,'cphi =')
        icount = icount + i
        i = index(uline,'Beta0 =')
        icount = icount + i
        i = index(uline,'Beta1 =')
        icount = icount + i
        i = index(uline,'Beta2 =')
        icount = icount + i
        i = index(uline,'Cphi =')
        icount = icount + i
        i = index(uline,'beta(0)')
        icount = icount + i
        i = index(uline,'beta(1)')
        icount = icount + i
        i = index(uline,'beta(2)')
        icount = icount + i
        i = index(uline,'Cphi')
        icount = icount + i
        i = index(uline,'C(phi)')
        icount = icount + i

        if (icount .ge. 24) then
            go to 100
        end if
    end do

100 continue

    if (icount .ge. 8) then
        uakey = 'Pitzer'
    else
        uakey = 'SEDH'
    end if

    rewind(ndat0s)
    go to 999

    ! Write a message for any read error other than end-of-file.
995 continue
    write (noutpt,2010)
    write (nttyo,2010)
2010 format(/' * Error - (EQPT/gakey) Encountered a read format',' error while',/7x,'scanning the DATA0 file to determine',' the associated aqueous',/7x,'activity coefficient model.')

    stop

999 continue
end subroutine gakey