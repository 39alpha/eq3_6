subroutine evdat3(arr,k,nmax,narxmx,narxt,ntpr,ntprmx,prop,tempc)
    !! This subroutine evaluates a thermodynamic property as a function
    !! of temperature, using an interpolating polynomial whose
    !! coefficients are stored in a 3D array arr. The second dimension
    !! of this array corresponds to a temperature range.
    !! Compare with EQLIB/evdat2.f, in which arr is a 2D array, and
    !! EQLIB/evdat4.f, in which this is a 4D array.
    !! This subroutine is called by:
    !!   EQLIB/alters.f
    !!   EQLIB/evdatr.f
    !! Principal input:
    !!   k      = index for the third dimension of array arr
    !!   arr    = two dimensional array of polynomial coefficients
    !!              describing some thermodynamic function
    !!   tempc  = temperature, C
    !!   ntpr   = temperature range flag
    !!   narxmx = first dimension of the arr array, the maximum number
    !!              of coefficients per temperature range
    !!   narxt  = array of numbers of coefficients in temperature
    !!              ranges
    !!   nmax   = third dimension of the arr array (can be ipchmx,
    !!              ipcvmx, ngtmax, nmtmax, or ngtmax)
    !! Prinicpal output:
    !!   prop   = the calculated property
    implicit none

    ! Calling sequence variable declarations.
    integer :: nmax
    integer :: narxmx
    integer :: ntprmx

    integer :: narxt(ntprmx)

    integer :: k
    integer :: ntpr

    real(kind=8) :: arr(narxmx,ntprmx,nmax)
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
        prop = arr(n,ntpr,k) + tempc*prop
    end do
end subroutine evdat3