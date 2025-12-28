subroutine ldbar3(ipc,ipcmax,nacdpr,narxmx,narxt,ndbmax,ntprmx,ntprt,xdbval,zdbval)
    !! This subroutine loads data on the temperature grid from the
    !! one-dimensional holding array (xdbval) into the ipc-th part
    !! of the three-dimensional proper data array (represented here by
    !! the dummy name zdbval).
    !! This subroutine is called by:
    !!   EQPT/rdpar.f
    !!   EQPT/pcraq.f
    !!   EQPT/pcrsg.f
    !! Principal input:
    !!   ipc    = the fixed dimension of the proper array
    !!   narxt  = array of numbers of coefficients in the temperature
    !!              ranges
    !!   ntprt  = the number of temperature ranges on the standard
    !!              temperature grid
    !!   xdbval = holding array (1-dimensional)
    !! Principal output:
    !!   nacdpr = array containging the number of actual data points
    !!              by range on the "log K" temperature grid; excludes
    !!   zdbval = dummy name for the proper data array (3-dimensional)
    implicit none

    ! Calling sequence variable declarations.
    integer :: ipcmax
    integer :: narxmx
    integer :: ndbmax
    integer :: ntprmx

    integer :: nacdpr(ntprmx)
    integer :: narxt(ntprmx)

    integer :: ipc
    integer :: ntprt

    real(kind=8) :: xdbval(ndbmax)
    real(kind=8) :: zdbval(narxmx,ntprmx,ipcmax)

    ! Local variable declarations.
    integer :: i
    integer :: k
    integer :: n
    integer :: ntpr

    ! Load the data from the xdbval array into the zdbval array.
    ! Note that a point which comprises the boundary between two
    ! ranges is a member of both ranges.
    ntpr = 1
    k = 0

    do n = 1,narxt(1)
        i = n
        zdbval(n,ntpr,ipc) = xdbval(i)

        if (xdbval(i) .lt. 9999999.) then
            k = k + 1
        end if
    end do

    nacdpr(ntpr) = k

    do ntpr = 2,ntprt
        k = 0
        zdbval(1,ntpr,ipc) = xdbval(i)

        if (xdbval(i) .lt. 9999999.) then
            k = k + 1
        end if

        do n = 2,narxt(ntpr)
            i = i + 1
            zdbval(n,ntpr,ipc) = xdbval(i)

            if (xdbval(i) .lt. 9999999.) then
                k = k + 1
            end if
        end do

        nacdpr(ntpr) = k
    end do
end subroutine ldbar3