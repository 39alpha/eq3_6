subroutine cliefu()
    !! This subroutine clears the IEEE flag for floating-point underflow,
    !! if such a flag is present, to avoid getting an unnecessary system
    !! warning message. Underflow is a normal condition in EQ3/6.
    !! The presence of an IEEE flag for underflow is platform and
    !! compiler dependent. If the flag is present, the need to clear
    !! it may also depend on the platform and compiler. Currently,
    !! the flag is present and needs to be cleared on Sun SPARCstations
    !! in the case of Sun's Fortran 77 compiler. This is not the case
    !! with any Fortran 90 compiler that has been seen to date,
    !! including Sun's.
    !! This subroutine is called by:
    !!   EQPT/eqpt.f
    !!   EQ3NR/eq3nr.f
    !!   EQ6/eq6.f
    !! Input:
    !!   None
    !! Output:
    !!   None
    implicit none

    !      Calling sequence variable declarations.
    !        None
    ! om   BEGIN_SPARC_DEPENDENT_CODE
    ! om     BEGIN_F77_DEPENDENT_CODE
    ! om
    ! om       Here cieeef() is a C subroutine that is distributed with
    ! om       EQ3/6. The object file is usually attached to the library
    ! om       (.a) file for the EQLIBU library, for convenience in linking.
    ! om
    !          call cieeef()
    !          go to 999
    ! om
    ! om     END_F77_DEPENDENT_CODE
    ! om   END_SPARC_DEPENDENT_CODE
999 continue
end subroutine cliefu