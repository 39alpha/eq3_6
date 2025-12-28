subroutine nactop(iopg,nopgmx,noutpt,nttyo,uactop)
    !! This subroutine sets the name of the option ("uactop") for
    !! computing the activity coefficients of aqueous species. It
    !! also sets associated logical flags concerning the generic type
    !! of activity coefficient option. The variable "iopg(1)" determines
    !! the exact activity coefficient model.
    !! This subroutine is called by:
    !!   EQ3NR/eq3nr.f
    !!   EQ6/eq6.f
    !! Principal input:
    !!   iopg   = array of activity coefficient option switches
    !! Principal output:
    !!   uactop = string describing the activity coefficient model
    !!              corresponding to iopg(1)
    implicit none

    ! Calling sequence variable declarations.
    integer :: nopgmx

    integer :: iopg(nopgmx)
    integer :: noutpt
    integer :: nttyo

    character(len=32) :: uactop

    ! Local variable declarations.
    !   None
    if (iopg(1) .eq. -1) then
        ! The Davies equation.
        uactop = 'the Davies equation'
    else if (iopg(1) .eq. 0) then
        ! The B-dot equation.
        uactop = 'the B-dot equation'
    else if (iopg(1) .eq. 1) then
        ! Pitzer's equations.
        uactop = "Pitzer's equations"
    else if (iopg(1) .eq. 2) then
        ! The HC + DHC equations.
        uactop = 'the HC + DH equations'
    else
        ! Have an error.
        write (noutpt,1000) iopg(1)
        write (nttyo,1000) iopg(1)
1000 format(/' * Error - (EQLIBG/nactop) Have iopg(1) = ',i3,'.',/7x,'This does not correspond to a valid option for computing',/7x,'the activity coefficients of aqueous species.')

        stop
    end if
end subroutine nactop