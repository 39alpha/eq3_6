subroutine pcorrm(adh,adhh,adhv,al10,aphi,avcnst,bdh,bdhh,bdhv,bdot,bdoth,bdotv,dadhh,dadhv,dbdhh,dbdhv,dbdth,dbdtv,dhfe,dvfe,iopg,ipch,ipchmx,ipcv,ipcvmx,nopgmx,presg,press,rcnstv,tempk,xhfe,xlke,xvfe)
    !! This subroutine makes pressure corrections for miscellaneous
    !! thermodynamic functions, such as Debye-Huckel parameters.
    !! It normally corrects for pressures off the standard T-P grid
    !! (e.g., from "presg" to "press"). However, it can be used to
    !! correct from any pressure to another by substituting the value
    !! of the former for "presg". This is done is EQ6 to correct for
    !! changing pressure along an isothermal reaction path. In order
    !! to do this, the pressure derivatives of the parameters in
    !! question are also corrected to the pressure of interest. The
    !! derivatives of highest order are necessarily treated as constants,
    !! so there is no correction in this case.
    !! Note the following technical references:
    !!   Ananthaswamy, J., and Atkinson, G., 1984, Thermodynamics of
    !!     concentrated electrolyte mixtures. 4. Pitzer-Debye-Huckel
    !!     limiting slopes for water from 0 to 100 C and from 1 atm
    !!     to 1 kbar, Journal of Chemical and Engineering Data,
    !!     v. 29, p. 81-87.
    !!   Bradley, D. J., and Pitzer, K. S., 1979, Thermodynamics of
    !!     electrolytes. 12. Dielectric properties of waer and
    !!     Debye-Huckel parameters to 350 C and 1 kbar, Journal
    !!     of Physical Chemistry, v. 83, p. 1599-1603.
    !!   Helgeson, H. C., and Kirkham, D. H., 1974b, Theoretical
    !!     prediction of the thermodynamic behavior of aqueous
    !!     electrolytes at high pressures and temperatures: II.
    !!     Debye-Huckel parameters for activity coefficients and
    !!     relative partial molal properties, American Journal
    !!     of Science, v. 274, p. 1199-1261.
    !! EQLIB/pcorrx.f performs the same function for the standard state
    !! thermodynamic functions, such as the equilibrium constants.
    !! This subroutine is called by:
    !!   EQ3NR/eq3nr.f
    !!   EQ6/eq6.f
    !! Principal input:
    !!  avcnst = 2.303 RT, with R in units of bar-cm3/mol-K
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: ipchmx
    integer :: ipcvmx
    integer :: nopgmx

    integer :: iopg(nopgmx)

    integer :: ipch
    integer :: ipcv

    real(kind=8) :: dadhh(ipchmx)
    real(kind=8) :: dadhv(ipcvmx)
    real(kind=8) :: dbdhh(ipchmx)
    real(kind=8) :: dbdhv(ipcvmx)
    real(kind=8) :: dbdth(ipchmx)
    real(kind=8) :: dbdtv(ipcvmx)
    real(kind=8) :: dhfe(ipchmx)
    real(kind=8) :: dvfe(ipcvmx)

    real(kind=8) :: adh
    real(kind=8) :: adhh
    real(kind=8) :: adhv
    real(kind=8) :: al10
    real(kind=8) :: aphi
    real(kind=8) :: avcnst
    real(kind=8) :: bdh
    real(kind=8) :: bdhh
    real(kind=8) :: bdhv
    real(kind=8) :: bdot
    real(kind=8) :: bdoth
    real(kind=8) :: bdotv
    real(kind=8) :: presg
    real(kind=8) :: press
    real(kind=8) :: rcnstv
    real(kind=8) :: tempk
    real(kind=8) :: xhfe
    real(kind=8) :: xlke
    real(kind=8) :: xvfe

    ! Local variable declarations.
    integer :: i
    integer :: ipc
    integer :: n

    real(kind=8) :: a22rt
    real(kind=8) :: avrt
    real(kind=8) :: dp
    real(kind=8) :: xx

    real(kind=8) :: fctrl

    ! a22rt  = 2(2.303)RT (R in P-V units)
    ! avrt   = RT (R in P-V units)
    a22rt = 2*avcnst
    avrt = rcnstv*tempk

    ! Calculate the pressure difference.
    dp = press - presg

    ! This section addresses various parameters at the Gibbs energy
    ! level (e.g., A(gamma,10), B(gamma)). It is divided into two
    ! parts, one for first-order corrections only, the other for
    ! corrections of order 2 or higher. This division is made for
    ! speed, as often the data will only be available for first-order
    ! corrections.
    if (ipcv .eq. 0) then
        ! Make only first-order corrections.
        ! The A(gamma,10) and A(phi) parameters. Note that A(V) (adhv)
        ! is not dA(gamma,10)/dp nor dA(phi)/dp, though it does contain
        ! these derivatives. Also, Helgeson and Kirkham's (1974b)
        ! definition of A(V) is not identical to Bradley and Pitzer's
        ! (1979). See the discussion of Ananthaswamy and Atkinson (1984).
        ! Here if iopg(1) is less than or equal to 0, we follow Helgeson
        ! and Kirkham's definition, otherwise, Bradley and Pitzer's.
        if (iopg(1) .le. 0) then
            ! Helgeson and Kirkham (1974b), eq 48:
            !   A(V) = -2(2.303)RT (dA(gamma,10)/dp)
            ! Compare Ananthaswamy and Atkinson (1984), table 1:
            !   A(V) = -6RT (dA(phi)/dp)
            xx = -adhv*dp/a22rt
            adh = adh + xx
            aphi = adh*al10/3.
        else
            ! Ananthaswamy and Atkinson (1984), table 1:
            !   A(V) = -4RT (dA(phi)/dp)
            xx = -adhv*dp/(4.*avrt)
            aphi = aphi + xx
            adh= 3.*aphi/al10
        end if

        ! The B(gamma) parameter.
        !   Helgeson and Kirkham (1974b), eq 49:
        !     B(V) = 2(2.303)RT (dB(gamma)/dp)
        ! Note that B(gamma) is carried multiplied by 10-8,
        ! while B(V) is carried multiplied by 10-6.
        xx = 0.01*bdhv*dp/a22rt
        bdh = bdh + xx

        ! The B-dot parameter.
        !   by analogy to Helgeson and Kirkham (1974b), eq 49:
        !     B-dot(V) = 2(2.303)RT (dB-dot/dp)
        xx = bdotv*dp/a22rt
        bdot = bdot + xx

        ! The log K for the "Eh" reaction.
        xx = -xvfe*dp/avcnst
        xlke = xlke + xx
    end if

    if (ipcv .gt. 0) then
        ! Make only corrections of order 2 or higher.
        ! The A(gamma,10) and A(phi) parameters.
        xx = adhv*dp

        do ipc = 1,ipcv
            n = ipc + 1
            xx = xx - ( dadhv(ipc)*(dp**n) )/fctrl(n)
        end do

        if (iopg(1) .le. 0) then
            xx = -xx/a22rt
            adh = adh + xx
            aphi = adh*al10/3.
        else
            xx = -xx/(4.*avrt)
            aphi = aphi + xx
            adh= 3.*aphi/al10
        end if

        ! The B(gamma) parameter.
        xx = bdhv*dp

        do ipc = 1,ipcv
            n = ipc + 1
            xx = xx - ( dbdhv(ipc)*(dp**n) )/fctrl(n)
        end do

        xx = 0.01*xx/a22rt
        bdh = bdh + xx

        ! Correct B-dot to the current pressure, using a correction of
        ! order 2 or higher.
        xx = bdotv*dp

        do ipc = 1,ipcv
            n = ipc + 1
            xx = xx - ( dbdtv(ipc)*(dp**n) )/fctrl(n)
        end do

        xx = xx/a22rt
        bdot = bdot + xx

        ! The log K for the "Eh" reaction.
        xx = -xvfe*dp

        do ipc = 1,ipcv
            n = ipc + 1
            xx = xx - ( dvfe(ipc)*(dp**n) )/fctrl(n)
        end do

        xx = xx/avcnst
        xlke = xlke + xx
    end if

    ! This section addresses enthalpy-related parameters.
    if (ipch .gt. 0) then
        ! The A(H) parameter.
        xx = 0.

        do ipc = 1,ipch
            xx = xx  + ( dadhh(ipc)*(dp**ipc) )/fctrl(ipc)
        end do

        adhh = adhh + xx

        do ipc = 1,ipch - 1
            ! The pressure derivatives of the A(H) parameter.
            xx = 0.

            do i = ipc + 1,ipch
                n = i - ipc
                xx = xx  + ( dadhh(i)*(dp**n) )/fctrl(n)
            end do

            dadhh(ipc) = dadhh(ipc) + xx
        end do

        ! The B(H) parameter.
        xx = 0.

        do ipc = 1,ipch
            xx = xx  + ( dbdhh(ipc)*(dp**ipc) )/fctrl(ipc)
        end do

        bdhh = bdhh + xx

        do ipc = 1,ipch - 1
            ! The pressure derivatives of the B(H) parameter.
            xx = 0.

            do i = ipc + 1,ipch
                n = i - ipc
                xx = xx  + ( dbdhh(i)*(dp**n) )/fctrl(n)
            end do

            dbdhh(ipc) = dbdhh(ipc) + xx
        end do

        ! The B-dot(H) parameter.
        xx = 0.

        do ipc = 1,ipch
            xx = xx  + ( dbdth(ipc)*(dp**ipc) )/fctrl(ipc)
        end do

        bdoth = bdoth + xx

        do ipc = 1,ipch - 1
            ! The pressure derivatives of the B-dot(H) parameter.
            xx = 0.

            do i = ipc + 1,ipch
                n = i - ipc
                xx = xx  + ( dbdth(i)*(dp**n) )/fctrl(n)
            end do

            dbdth(ipc) = dbdth(ipc) + xx
        end do

        ! The enthalpy of the "Eh" reaction.
        xx = 0.

        do ipc = 1,ipch
            xx = xx  + ( dhfe(ipc)*(dp**ipc) )/fctrl(ipc)
        end do

        xhfe = xhfe + xx

        do ipc = 1,ipch - 1
            ! The pressure derivatives of the enthalpy of the "Eh" reaction.
            xx = 0.

            do i = ipc + 1,ipch
                n = i - ipc
                xx = xx  + ( dhfe(i)*(dp**n) )/fctrl(n)
            end do

            dhfe(ipc) = dhfe(ipc) + xx
        end do
    end if

    ! This section addresses volume-related parameters.
    if (ipcv .gt. 0) then
        ! The A(V) parameter.
        xx = 0.

        do ipc = 1,ipcv
            xx = xx  + ( dadhv(ipc)*(dp**ipc) )/fctrl(ipc)
        end do

        adhv = adhv + xx

        do ipc = 1,ipcv - 1
            ! The pressure derivatives of the A(V) parameter.
            xx = 0.

            do i = ipc + 1,ipcv
                n = i - ipc
                xx = xx  + ( dadhv(i)*(dp**n) )/fctrl(n)
            end do

            dadhv(ipc) = dadhv(ipc) + xx
        end do

        ! The B(V) parameter.
        xx = 0.

        do ipc = 1,ipcv
            xx = xx  + ( dbdhv(ipc)*(dp**ipc) )/fctrl(ipc)
        end do

        bdhv = bdhv + xx

        do ipc = 1,ipcv - 1
            ! The pressure derivatives of the B(V) parameter.
            xx = 0.

            do i = ipc + 1,ipcv
                n = i - ipc
                xx = xx  + ( dbdhv(i)*(dp**n) )/fctrl(n)
            end do

            dbdhv(ipc) = dbdhv(ipc) + xx
        end do

        ! The B-dot(V) parameter.
        xx = 0.

        do ipc = 1,ipcv
            xx = xx  + ( dbdtv(ipc)*(dp**ipc) )/fctrl(ipc)
        end do

        bdotv = bdotv + xx

        do ipc = 1,ipcv - 1
            ! The pressure derivatives of the B-dot(V) parameter.
            xx = 0.

            do i = ipc + 1,ipcv
                n = i - ipc
                xx = xx  + ( dbdtv(i)*(dp**n) )/fctrl(n)
            end do

            dbdtv(ipc) = dbdtv(ipc) + xx
        end do

        ! The volume of the "Eh" reaction.
        xx = 0.

        do ipc = 1,ipcv
            xx = xx  + ( dvfe(ipc)*(dp**ipc) )/fctrl(ipc)
        end do

        xvfe = xvfe + xx

        do ipc = 1,ipcv - 1
            ! The pressure derivatives of the volume of the "Eh" reaction.
            xx = 0.

            do i = ipc + 1,ipcv
                n = i - ipc
                xx = xx  + ( dvfe(i)*(dp**n) )/fctrl(n)
            end do

            dvfe(ipc) = dvfe(ipc) + xx
        end do
    end if
end subroutine pcorrm