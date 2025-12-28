subroutine etoo2(cdrsi,dhfe,dhfs,dvfe,dvfs,ipch,ipchmx,ipcv,ipcvmx,narxmx,narxt,nbtmx1,ndrsts,ns,ntprmx,ntprt,udrsi,xhfe,xhfs,xlke,xlks,xvfe,xvfs)
    !! This subroutine converts the reaction for the ns-th species from
    !! one written in terms of e- to one written in terms of O2(g).
    !! The "Eh" reaction is used to do this.
    !! This subroutine is called by:
    !!   EQPT/pcraq.f
    !!   EQPT/pcrsg.f
    !! Input:
    !!   cdrsi  = array of coefficients for the original reaction
    !!              associated with the ns-th species
    !!   nbtmx1 = the maximum number of basis species plus 1
    !!   ndrsts = number of species in the reaction
    !!   ns     = index of the species whose reaction is to be rewritten
    !!   udrsi  = names of species in the original reaction
    !!   xlke   = the array of log K values for the "Eh reaction"
    !!   xlks   = the array of log K values for the original reaction
    !! Output:
    !!   cdrsi  = array of coefficients for a reaction, rewritten
    !!              in terms of e- instead of O2(g)
    !!   udrsi  = names of species in the rewritten reaction
    !!   xlks   = the array of log K values for the rewritten reaction
    implicit none

    ! Calling sequence variable declarations.
    integer :: ipchmx
    integer :: ipcvmx
    integer :: narxmx
    integer :: nbtmx1
    integer :: ntprmx

    integer :: narxt(ntprmx)

    integer :: ipch
    integer :: ipcv
    integer :: ndrsts
    integer :: ns
    integer :: ntprt

    character(len=24) :: udrsi(nbtmx1)

    real(kind=8) :: cdrsi(nbtmx1)
    real(kind=8) :: dhfe(narxmx,ntprmx,ipchmx)
    real(kind=8) :: dhfs(narxmx,ntprmx,ipchmx,nbtmx1)
    real(kind=8) :: dvfe(narxmx,ntprmx,ipcvmx)
    real(kind=8) :: dvfs(narxmx,ntprmx,ipcvmx,nbtmx1)
    real(kind=8) :: xhfe(narxmx,ntprmx)
    real(kind=8) :: xhfs(narxmx,ntprmx,nbtmx1)
    real(kind=8) :: xlke(narxmx,ntprmx)
    real(kind=8) :: xlks(narxmx,ntprmx,nbtmx1)
    real(kind=8) :: xvfe(narxmx,ntprmx)
    real(kind=8) :: xvfs(narxmx,ntprmx,nbtmx1)

    ! Local variable declarations.
    integer :: i
    integer :: ipc
    integer :: kdrsts
    integer :: n
    integer :: nelect
    integer :: nhydr
    integer :: ntpr
    integer :: nwater

    real(kind=8) :: factor

    ! Search for e-, H+, and H2O in the reaction.
    nelect = 0
    nhydr = 0
    nwater = 0

    do n = 2,ndrsts
        if (udrsi(n)(1:3) .eq. 'e- ') then
            nelect = n
        end if

        if (udrsi(n)(1:3) .eq. 'H+ ') then
            nhydr = n
        end if

        if (udrsi(n)(1:4) .eq. 'H2O ') then
            nwater = n
        end if
    end do

    if (nelect .eq. 0) then
        go to 999
    end if

    ! Re-write the reaction, using the Eh reaction:
    !   H2O = 4H+ + 4e- + O2(g)
    factor = -cdrsi(nelect)/4.

    udrsi(nelect) = 'O2(g)'
    cdrsi(nelect) = factor

    if (nhydr .gt. 0) then
        cdrsi(nhydr) = cdrsi(nhydr) + 4.*factor
    else
        ndrsts = ndrsts + 1
        nhydr = ndrsts
        udrsi(nhydr) = 'H+'
        cdrsi(nhydr) = 4.*factor
    end if

    if (nwater .gt. 0) then
        cdrsi(nwater) = cdrsi(nwater) - 2.*factor
    else
        ndrsts = ndrsts + 1
        nwater = ndrsts
        udrsi(nwater) = 'H2O'
        cdrsi(nwater) = -2.*factor
    end if

    ! Clear any zeroes in the reaction coefficients.
    kdrsts = ndrsts

    do n = 1,ndrsts
        if (abs(cdrsi(n)) .le. 1.e-6) then
            if (n .le. kdrsts - 1) then
                do i = n,kdrsts - 1
                    cdrsi(i) = cdrsi(i + 1)
                    udrsi(i) = udrsi(i + 1)
                end do
            end if

            kdrsts = kdrsts - 1
        end if

        if (n .ge. kdrsts) then
            go to 150
        end if
    end do

150 continue
    ndrsts = kdrsts

    ! Modify the log K values.
    do ntpr = 1,ntprt
        do n = 1,narxt(ntpr)
            if (xlks(n,ntpr,ns).lt.9999999. .and.      xlke(n,ntpr).lt.9999999.) then
                xlks(n,ntpr,ns) = xlks(n,ntpr,ns) + factor*xlke(n,ntpr)
            else
                xlks(n,ntpr,ns) = 9999999.
            end if
        end do
    end do

    if (ipch .ge. 0) then
        do ntpr = 1,ntprt
            do n = 1,narxt(ntpr)
                if (xhfs(n,ntpr,ns).lt.9999999. .and.        xhfe(n,ntpr).lt.9999999.) then
                    xhfs(n,ntpr,ns) = xhfs(n,ntpr,ns) + factor*xhfe(n,ntpr)
                else
                    xhfs(n,ntpr,ns) = 9999999.
                end if
            end do
        end do

        do ipc = 1,ipch
            do ntpr = 1,ntprt
                do n = 1,narxt(ntpr)
                    if (dhfs(n,ntpr,ipc,ns).lt.9999999. .and.          dhfe(n,ntpr,ipc).lt.9999999.) then
                        dhfs(n,ntpr,ipc,ns) = dhfs(n,ntpr,ipc,ns)          + factor*dhfe(n,ntpr,ipc)
                    else
                        dhfs(n,ntpr,ipc,ns) = 9999999.
                    end if
                end do
            end do
        end do
    end if

    if (ipcv .ge. 0) then
        do ntpr = 1,ntprt
            do n = 1,narxt(ntpr)
                if (xvfs(n,ntpr,ns).lt.9999999. .and.        xvfe(n,ntpr).lt.9999999.) then
                    xvfs(n,ntpr,ns) = xvfs(n,ntpr,ns) + factor*xvfe(n,ntpr)
                else
                    xvfs(n,ntpr,ns) = 9999999.
                end if
            end do
        end do

        do ipc = 1,ipcv
            do ntpr = 1,ntprt
                do n = 1,narxt(ntpr)
                    if (dvfs(n,ntpr,ipc,ns).lt.9999999. .and.          dvfe(n,ntpr,ipc).lt.9999999.) then
                        dvfs(n,ntpr,ipc,ns) = dvfs(n,ntpr,ipc,ns)          + factor*dvfe(n,ntpr,ipc)
                    else
                        dvfs(n,ntpr,ipc,ns) = 9999999.
                    end if
                end do
            end do
        end do
    end if

999 continue
end subroutine etoo2