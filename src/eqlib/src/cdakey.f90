subroutine cdakey(iopg,nopgmx,noutpt,nttyo,udakey,udatfi)
    !! This subroutine checks the flag string "udakey" read from the data
    !! file and checks it against a key-list to test whether or not
    !! the data file used is consistent with the specified model for
    !! aqueous species activity coefficients (determined by iopg(1)).
    !! The string "udatfi" identifies the data file. The legal
    !! combinations are:
    !!       Option                 iopg(1)      Legal key
    !!    Davies equation             -1          SEDH
    !!    B-dot equation               0          SEDH
    !!    Pitzer's equations           1          Pitzer
    !!    HC + DHC equations           2          SEDH
    !! This subroutine is called by:
    !!   EQ3NR/eq3nr.f
    !!   EQ6/eq6.f
    !!  Principal input:
    !!  Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: nopgmx

    integer :: noutpt
    integer :: nttyo

    integer :: iopg(nopgmx)

    character(len=8) :: udakey
    character(len=8) :: udatfi

    ! Local variable declarations.
    integer :: j2
    integer :: j3
    integer :: j4
    integer :: j5

    integer :: ilnobl

    character(len=8) :: uiakey
    character(len=8) :: usedh
    character(len=8) :: upitz
    character(len=8) :: ux8

    data usedh  /'SEDH    '/,upitz/'Pitzer  '/

    ! Set the correct key.
    if (iopg(1) .eq. -1) then
        uiakey = usedh
    else if (iopg(1) .eq. 0) then
        uiakey = usedh
    else if (iopg(1) .eq. 1) then
        uiakey = upitz
    else if (iopg(1) .eq. 2) then
        uiakey = usedh
    else
        write (ux8,'(i5)') iopg(1)
        call lejust(ux8)
        j2 = ilnobl(ux8)
        write (noutpt,1000) ux8(1:j2)
        write (nttyo,1000) ux8(1:j2)
1000 format(/' * Error - (EQLIB/cdakey) The iopg(1) value of ',a,' specified on',/7x,"the input file doesn't correspond to an",' available model for the',/7x,'activity coefficients of',' the aqueous species.')

        stop
    end if

    if (udakey(1:8) .ne. uiakey(1:8)) then
        ! Error, bad combination.
        write (ux8,'(i5)') iopg(1)
        call lejust(ux8)
        j2 = ilnobl(ux8)
        j3 = ilnobl(udatfi)
        j4 = ilnobl(uiakey)
        j5 = ilnobl(udakey)
        write (noutpt,1010) ux8(1:j2),udatfi(1:j3),uiakey(1:j4),udakey(1:j5)
        write (nttyo,1010) ux8(1:j2),udatfi(1:j3),uiakey(1:j4),udakey(1:j5)
1010 format(/' * Error - (EQLIB/cdakey) The iopg(1) value of ',a,' specified on',/7x,"the input file doesn't correspond to the",' type of model for the',/7x,'activity coefficients of the',' aqueous species required to use the',/7x,a,' data file.',' The iopg1(1) value calls for using a ',a,' type model,',/7x,'while the data file is consistent with a ',a,' type',' model. Change',/7x,'the iopg(1) value on the input file',' or specify a compatible data file.')

        stop
    end if
end subroutine cdakey