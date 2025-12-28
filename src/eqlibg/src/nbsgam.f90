subroutine nbsgam(acfnbs,adh,fxi,nchlor,noutpt,nttyo)
    !! This subroutine calculates the log activity coefficient of the
    !! chloride ion according to the NBS pH convention (e.g., Covington,
    !! Bates, and Durst, 1985). The convention itself may be extended
    !! outside the specified limit of 0.1 molal on the ionic strength.
    !! Reference:
    !!   Covington, A.K., Bates, R.G., and Durst, R.A., 1985,
    !!   Definition of pH scales, standard reference values, measure-
    !!   ment of pH and related terminology (recommendations, 1984):
    !!   Pure and Applied Chemistry, v. 57, p. 533-542.
    !! This subroutine is called by:
    !!   EQLIB/gpheh.f
    !!   EQLIBG/gcoeff.f
    !! Principal input:
    !!   adh    = Debye-Huckel A(gamma) parameter
    !!   fxi    = the ionic strength (the 2nd-order electrostatic
    !!              moment function I)
    !!   nchlor = the species index of the chloride ion
    !! Principal output:
    !!   acfnbs = log gamma(Cl-), according to the Bates-Guggenheim
    !!            convention
    implicit none

    ! Calling sequence variable declarations.
    integer :: nchlor
    integer :: noutpt
    integer :: nttyo

    real(kind=8) :: acfnbs
    real(kind=8) :: adh
    real(kind=8) :: fxi

    ! Local variable declarations.
    real(kind=8) :: fxisqt

    ! Test the species index of the chloride ion.
    if (nchlor .le. 0) then
        write (noutpt,1000)
        write (nttyo ,1000)
1000 format(/' * Error - (EQLIBG/nbsgam) Have no index for the',/7x,"chloride ion. Can't use the extended NBS pH scale.")

        stop
    end if

    ! Evaluate the NBS expression for the molal activity coefficient
    ! of the chloride ion.
    fxisqt = sqrt(fxi)
    acfnbs = - ( adh * fxisqt ) / ( 1.0 + ( 1.5 * fxisqt ) )
end subroutine nbsgam