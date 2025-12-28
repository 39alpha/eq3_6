subroutine evdat4(arr,ipc,ipcxmx,k,nmax,narxmx,narxt,ntpr,ntprmx,prop,tempc)
    !! This subroutine evaluates a thermodynamic property as a function
    !! of temperature, using an interpolating polynomial whose
    !! coefficients are stored in a 4D array arr. The second dimension
    !! of this array corresponds to a temperature range. The third
    !! dimension of this array usually corresponds to an order parameter,
    !! as for pressure correction.
    !! Compare with EQLIB/evdat3.f, in which arr is a 3D array, and
    !! EQLIB/evdat3.f, in which this is a 3D array.
    !! This subroutine is called by:
    !!   EQLIB/evdatr.f
    !! Principal input:
    !!   k      = the index for the fourth dimension of the arr array
    !!   arr    = two dimensional array of polynomial coefficients
    !!              describing some thermodynamic function
    !!   ipcxmx = the third dimension of the arr array
    !!   ipc    = the index for the third dimension of the arr array
    !!   tempc  = temperature, C
    !!   ntpr   = temperature range flag
    !!   narxmx = first dimension of the arr array, the number
    !!              of coefficients per temperature range
    !!   narxt  = array of numbers of coefficients in temperature
    !!              ranges
    !!   ntprmx = second dimension of the arr array, the number
    !!              of temperature ranges.
    !!   nmax   = third dimension of the arr array (can be ngtmax,
    !!              nmtmax, or ngtmax)
    !! Prinicpal output:
    !!   prop   = the calculated property array
    implicit none

    ! Calling sequence variable declarations.
    integer :: ipcxmx
    integer :: nmax
    integer :: narxmx
    integer :: ntprmx

    integer :: narxt(ntprmx)

    integer :: ipc
    integer :: k
    integer :: ntpr

    real(kind=8) :: arr(narxmx,ntprmx,ipcxmx,nmax)
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
        prop = arr(n,ntpr,ipc,k) + tempc*prop
    end do
end subroutine evdat4