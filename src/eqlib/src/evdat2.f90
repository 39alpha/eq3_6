subroutine evdat2(arr,narxmx,narxt,ntpr,ntprmx,prop,tempc)
    !! This subroutine evaluates a thermodynamic property as a function
    !! of temperature, using an interpolating polynomial whose
    !! coefficients are stored in a 2D array arr. The second dimension
    !! of this array corresponds to a temperature range.
    !! Compare with EQLIB/evdat3.f, in which arr is a 3D array, and
    !! EQLIB/evdat4.f, in which this is a 4D array.
    !! This subroutine is called by:
    !!   EQLIB/alters.f
    !!   EQLIB/evdata.f
    !!   EQ6/eq6.f
    !! Principal input:
    !!   arr    = two dimensional array of polynomial coefficients
    !!              describing some thermodynamic function
    !!   tempc  = temperature, C
    !!   ntpr   = temperature range flag
    !!   narxmx = first dimension of the arr array, the maximum number
    !!              of coefficients per temperature range
    !!   narxt  = array of numbers of coefficients in temperature
    !!              ranges
    !!   ntprmx = second dimension of the arr array, the number
    !!              of temperature ranges.
    !! Principal output:
    !!   prop   = the calculated property
    implicit none

    ! Calling sequence variable declarations.
    integer :: narxmx
    integer :: ntprmx

    integer :: narxt(ntprmx)

    integer :: ntpr

    real(kind=8) :: arr(narxmx,ntprmx)
    real(kind=8) :: prop
    real(kind=8) :: tempc

    ! Local variable declarations.
    integer :: n
    integer :: nn
    integer :: nt

    ! Evaluate the polynomial.
    prop = 0.
    nt = narxt(ntpr)

    do nn = 1,nt
        n = nt + 1 - nn
        prop = arr(n,ntpr) + tempc*prop
    end do
end subroutine evdat2