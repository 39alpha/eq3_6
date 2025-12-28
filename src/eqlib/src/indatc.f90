subroutine indatc(arr,nad1,narx_asv,narxt,ntpr_asv,ntprt,ux24)
    !! This subroutine reads from the unformatted data file (unit number
    !! nad1) the 2D array arr, which contains coefficients of
    !! interpolating polynomials representing a thermodynamic property
    !! as a function of temperature. There are ntprt temperature
    !! ranges, and a set of narxt(ntpr) coefficients for the ntpr-th
    !! temperature range. When this subroutine is called, another array
    !! name is usually substituted for arr.
    !! Compare with EQLIB/indatd.f, in which arr is a 3D array.
    !! This subroutine is called by:
    !!   EQLIB/indata.f
    !! Principal input:
    !!   nad1   = unit number of the data file
    !! Principal output:
    !!   arr    = 2D array of polynomial coefficients
    implicit none

    ! Calling sequence variable declarations.
    integer :: narx_asv
    integer :: ntpr_asv

    integer :: nad1

    integer :: narxt(ntpr_asv)

    integer :: ntprt

    real(kind=8) :: arr(narx_asv,ntpr_asv)

    character(len=24) :: ux24

    ! Local variable declarations.
    integer :: n
    integer :: ntpr

    ! The content of ux24 may be useful in debugging, but has no
    ! other usage.
    read (nad1) ux24

    do ntpr = 1,ntprt
        read (nad1) (arr(n,ntpr), n = 1,narxt(ntpr))
    end do
end subroutine indatc