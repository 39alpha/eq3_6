subroutine gspidx(ier,n,nat,natmax,uaqsp,unams)
    !! Get the index (n) of the aqueous species whose name is unams.
    !! This subroutine is called by:
    !!   EQPT/rdpz2.f
    !!   EQPT/rdpz3.f
    !! Principal input:
    !!   nat    = the number of aqueous species
    !!   uaqsp  = array of names of aqueous species
    !!   unams  = the name of an aqueous species whose index n is
    !!              desired
    !! Principal output:
    !!   ier    = error flag
    !!   n      = the index of the species unams in the uaqsp array
    implicit none

    ! Calling sequence variable declarations.
    integer :: natmax

    integer :: ier
    integer :: n
    integer :: nat

    character(len=24) :: uaqsp(natmax)
    character(len=24) :: unams

    ! Local variable declarations.
    ier = 0

    do n = 1,nat
        if (unams(1:24) .eq. uaqsp(n)(1:24)) then
            go to 100
        end if
    end do

    ! The species was not found. Set the error flag.
    ier = 1
    n = 0

100 continue
end subroutine gspidx