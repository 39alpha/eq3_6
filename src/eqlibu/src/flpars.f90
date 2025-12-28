subroutine flpars(eps100,irang,noutpt,nttyo,smp100)
    !! This subroutine obtains the following real*8 parameters: the
    !! real*8 machine epsilon (eps), the smallest positive real*8
    !! number (smpos), and the exponent range (irang). It returns
    !! eps100 = 100*eps and smp100 = 100*smpos, not of eps and smpos.
    !! This subroutine calls EQLIBU/chump.f to test the adequacy of the
    !! machine epsilon and the exponent range. If inadequacies are
    !! noted, warnings will be written to the output and screen files.
    !! This subroutine is called by:
    !!   EQPT/eqpt.f
    !!   EQ3NR/eq3nr.f
    !!   EQ6/eq6.f
    !! Input:
    !!   noutpt = unit number of the output file
    !!   nttyo  = unit number of the screen file
    !! Output:
    !!   eps100 = 100*eps
    !!   smp100 = 100*smpos
    !!   irang  = the real*8 exponent range, int( -log10(smpos ))
    !! Local:
    !!   eps    = the real*8 machine epsilon
    !!   smpos  = the smallest positive real*8 number
    implicit none

    ! Calling sequence variable declarations.
    integer :: noutpt
    integer :: nttyo

    integer :: irang

    real(kind=8) :: eps100
    real(kind=8) :: smp100

    ! Local variable declarations.
    real(kind=8) :: eps
    real(kind=8) :: smpos
    real(kind=8) :: x

    data x / 1.0 /

    ! Get eps (real*8 machine epsilon).
    eps = epsilon(x)

    ! Get smpos(real*8 smallest positive number).
    smpos = tiny(x)

    eps100 = 100*eps
    smp100 = 100*smpos

    ! Compute the exponent range.
    irang = int(-log10(smpos))

    ! Clear the IEEE flag for floating-point underflow, if such a
    ! flag is present, to avoid getting an unnecessary system
    ! warning message. Underflow is a normal condition in EQ3/6.
    ! Make any porting changes in the EQLIBU subroutine cliefu.
    ! Do not make the porting changes here.
    ! Note: the calculation of smpos above involves forcing an
    ! underflow.
    call cliefu()

    ! Test the adequacy of the machine epsilon and the exponent range.
    call chump(eps,irang,noutpt,nttyo)
end subroutine flpars