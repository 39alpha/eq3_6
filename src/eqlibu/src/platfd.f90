subroutine platfd(uplatc,uplatm)
    !! This subroutine sets the platform designator strings that are
    !! written on the output and screen files of various codes.
    !! This subroutine is called by:
    !!   EQPT/eqpt.f
    !!   EQ3NR/eq3nr.f
    !!   EQ6/eq6.f
    !! Input:
    !!   None
    !! Output:
    !!   uplatc = platform category (e.g., UNIX, PC, MAC)
    !!   uplatm = platform machine (e.g., SPARC, SGI, Pentium, 486PC)
    implicit none

    ! Calling sequence variable declarations.
    character(len=8) :: uplatc
    character(len=8) :: uplatm

    !      Local variable declarations.
    !        None
    ! om   BEGIN_UNIX_DEPENDENT_CODE
    ! om
    !        uplatc = 'UNIX'
    ! om
    !        uplatm = 'SPARC'
    ! xx     uplatm = 'SGI'
    ! xx     uplatm = 'HP-UX'
    ! xx     uplatm = 'AIX'
    ! xx     uplatm = 'Ultrix'
    ! om
    ! om   END_UNIX_DEPENDENT_CODE
    ! om   BEGIN_PC_DEPENDENT_CODE
    ! om
    uplatc = 'PC'

    ! om
    ! xx     uplatm = 'Pentium II'
    ! xx     uplatm = 'K6'
    ! xx     uplatm = 'K6-2'
    uplatm = 'PC'

    ! xx     uplatm = 'Pentium Pro'
    ! xx     uplatm = 'P5'
    ! xx     uplatm = 'Pentium'
    ! xx     uplatm = '486PC'
    ! xx     uplatm = '386PC'
    ! om
    ! om   END_PC_DEPENDENT_CODE
    ! om   BEGIN_MAC_DEPENDENT_CODE
    ! om
    !        uplatc = 'MAC'
    ! om
    !        uplatm = 'MAC'
    ! om
    ! om   END_MAC_DEPENDENT_CODE
    ! om   BEGIN_VAX_DEPENDENT_CODE
    ! om
    !        uplatc = 'VAX/VMS'
    ! om
    !        uplatm = 'VAX'
    ! om
    ! om   END_VAX_DEPENDENT_CODE
end subroutine platfd