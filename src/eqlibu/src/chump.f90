subroutine chump(eps,irang,noutpt,nttyo)
    !! This subroutine tests the floating point precision and exponent
    !! range for adequacy for EQ3/6 calculations. This is essentially
    !! a trap to insure that most floating point variables are real*8
    !! and, if the computer is a VAX, that the code has been compiled
    !! with the 'G_FLOATING' option (which increases the exponent range).
    !! The machine epsilon and exponent range are respectively tested
    !! against epstst and irgtst. The former should be no greater than
    !! 1.e-14. The latter should be at least 100, and the recommended
    !! value is 300.
    !! This subroutine is called by:
    !!    EQLIBU/flpars.f
    !! Input:
    !!   noutpt = the unit number of the output file
    !!   nttyo  = the unit number of the screen file
    !! Output:
    !!   eps    = the real*8 machine epsilon
    !!   irang  = the real*8 exponent range
    implicit none

    ! Calling sequence variable declarations.
    integer :: irang
    integer :: noutpt
    integer :: nttyo

    real(kind=8) :: eps

    ! Local variable declarations.
    integer :: irgtst

    logical :: qwarn

    real(kind=8) :: epstst

    ! WARNING: Under no circumstances should the value of irgtst be
    ! changed to a value less than 100. If decreased to from the
    ! recommended value of 300, make a corresponding adjustment of the
    ! overflow truncation limit in EQLIBU/texp.f or you will lose the
    ! floating-point overflow protection built into this software.
    data irgtst /300/

    data epstst /1.e-14/

    qwarn = .false.

    if (eps .gt. epstst) then
        write (noutpt,1007) eps,epstst
        write (nttyo,1007) eps,epstst
1007 format(/' * Warning - (EQLIBU/chump) Have insufficient',' floating point',/7x,'epsilon= ',1pe12.5,'. This parameter',' should be no larger',/7x,'than ',1pe12.5,'. You need a',' computer with better floating',/7x,'point characteristics',' to run this software.')

        qwarn=.true.
    end if

    if (irang .lt. irgtst) then
        write (noutpt,1008) irang,irgtst
        write (nttyo,1008) irang,irgtst
1008 format(/' * Warning - (EQLIBU/chump) Have insufficient',' floating point',/7x,'exponent range= +/- ',i4,'. This should',' be at least +/- ',i4,/7x,'to run this software. You need',' a computer with better',/7x,'floating point characteristics.',' Alternatively, you need to recompile',/7x,'with a special',' compiler option.')

        if (irang.ge.37 .and. irang.le.39) then
            write (noutpt,1010)
            write (nttyo,1010)
1010 format(/' * Note - (EQLIBU/chump) This computer appears to',' be a VAX.',/7x,'If you run under VMS, try recompiling with',/7x,'the G_FLOATING option.',/7x,'This will increase the',' exponent range to +/- 308. Some older',/7x,'VAXes do not',' have the necessary hardware to implement this option.',/7x,'Also, there appears to be no equivalent to the',' G_FLOATING option',/7x,'on VAXes running the ULTRIX',' operating system. If you are',/7x,'running a VAX with',' ULTRIX, activate the call to the system subroutine TRPFPE',/7x,'in this subroutine (EQLIBU/chump.f). This call is',' normally commented.',/7x,'out. You might wish to consider',' getting a computer which has better',/7x,'floating point',' characteristics.')
        end if

        qwarn=.true.
    end if

    if (qwarn) then
        write (noutpt,1020)
        write (nttyo,1020)
1020 format(/' * Warning - (EQLIBU/chump) Will continue execution.',' You may encounter',/7x,'floating point overflows, especially',' if you are running EQ3NR',/7x,'or EQ6. You may be able to',' make some EQ3NR and EQ6 problems run',/7x,'by specifying',' judicious basis switches on the input file.',/)
    end if

    ! om   BEGIN_VAX_DEPENDENT_CODING
    ! om
    ! om     The following VAX system call traps and repairs floating point
    ! om     overflows when running under the ULTRIX operating system.
    ! om     This call may have be to be activated (by de-commenting it) if
    ! om     you are running on a VAX with ULTRIX. On some older VAXES, you
    ! om     may have to substitute "TRAPOV" for "TRPFPE".
    ! om
    !        call trpfpe(0,1.d+30)
    ! om
    ! om   END_VAX_DEPENDENT_CODING
end subroutine chump