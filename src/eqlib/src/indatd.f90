subroutine indatd(arr,ipc,ipcx_asv,nad1,narxt,narx_asv,ntprt,ntpr_asv,ux24)
    !! This subroutine reads from the unformatted data file (unit number
    !! nad1) the 3D array arr, which contains coefficients of
    !! interpolating polynomials representing a thermodynamic property
    !! as a function of temperature. There are ntprt temperature
    !! ranges, and a set of narxt(ntpr) coefficients for the ntpr-th
    !! temperature range for the ipc-th entity. When this subroutine is
    !! called, another array name is usually substituted for arr and
    !! a corresponding dimensioning variable is substituted for ipcx_asv.
    !! Compare with EQLIB/indatc.f, in which arr is a 2D array.
    !! This subroutine is called by:
    !!   EQLIB/indata.f
    !! Principal input:
    !!   nad1   = unit number of the data file
    !! Principal output:
    !!   arr    = 3D array of polynomial coefficients
    implicit none

    ! Calling sequence variable declarations.
    integer :: ipcx_asv
    integer :: narx_asv
    integer :: ntpr_asv

    integer :: nad1

    integer :: narxt(ntpr_asv)

    integer :: ipc
    integer :: ntprt

    real(kind=8) :: arr(narx_asv,ntpr_asv,ipcx_asv)

    character(len=24) :: ux24

    ! Local variable declarations.
    integer :: n
    integer :: ntpr

    ! The content of ux24 may be useful in debugging, but has no
    ! other usage.
    read (nad1) ux24

    do ntpr = 1,ntprt
        read (nad1) (arr(n,ntpr,ipc), n = 1,narxt(ntpr))
    end do
end subroutine indatd