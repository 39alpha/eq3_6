subroutine garmat(act,afrc1,aimatr,al10,armatr,cdac,cdacb,cdrs,csigma,csts,delvec,dlogxw,dvjdte,eact,eps100,fkrc,gmmatr,hact,iact,idirec,iindx1,iktmax,imchmx,imech,ipivot,jcode,jreac,jtemp,kbt,kdim,kmax,mmmatr,morr,mwtrc,nbasp,nbt,nbtmax,ncmpr,ndac,ndacb,ndact,ndctmx,ndrs,ndrsmx,ndrsr,noutpt,nptmax,nrct,nrctmx,nrct1,nrk,nrndex,nsk,nstmax,nsts,nstsmx,nstsr,nttkmx,nttyo,nxridx,nxrtmx,rirec1,rk,rkb,rreac1,rrelr1,rrxfi1,rtcnst,rxbar,sfcar,sgmatr,ssfcar,tempc,tempcb,tempk,ttk,ureac,whcfac,xi1,xlks,xxmatr,xymatr)
    !! This subroutine calculates the Jacobian matrix J[r] (armatr)
    !! used by the higher-order (stiff) ODE integrator.
    !! This subroutine is called by:
    !!   EQ6/path.f
    !! Principal input:
    !! Principal output:
    !!   armatr = the Jacobian matrix J[r]
    implicit none

    ! Calling sequence variable declarations.
    integer :: iktmax
    integer :: imchmx
    integer :: kmax
    integer :: nbtmax
    integer :: ndctmx
    integer :: ndrsmx
    integer :: nrctmx
    integer :: nptmax
    integer :: nstmax
    integer :: nstsmx
    integer :: nttkmx
    integer :: nxrtmx

    integer :: noutpt
    integer :: nttyo

    integer :: jtemp
    integer :: kbt
    integer :: kdim
    integer :: nbt
    integer :: nrct
    integer :: nrct1

    integer :: iact(imchmx,2,nrctmx)
    integer :: idirec(nrctmx)
    integer :: iindx1(kmax)
    integer :: imech(2,nrctmx)
    integer :: ipivot(kmax)
    integer :: jcode(nrctmx)
    integer :: jreac(nrctmx)
    integer :: nbasp(nbtmax)
    integer :: ndac(ndctmx,imchmx,2,nrctmx)
    integer :: ncmpr(2,nptmax)
    integer :: ndacb(nbt,imchmx,2,nrct)
    integer :: ndact(imchmx,2,nrctmx)
    integer :: ndrs(ndrsmx)
    integer :: ndrsr(2,nstmax)
    integer :: nrk(2,nrctmx)
    integer :: nrndex(nrctmx)
    integer :: nsk(nrctmx)
    integer :: nsts(nstsmx)
    integer :: nstsr(2,nstmax)
    integer :: nxridx(nrctmx)

    character(len=24) :: ureac(nrctmx)

    real(kind=8) :: act(nstmax)
    real(kind=8) :: afrc1(nrctmx)
    real(kind=8) :: aimatr(kmax,kmax)
    real(kind=8) :: armatr(nrct1,nrct1)
    real(kind=8) :: cdac(ndctmx,imchmx,2,nrctmx)
    real(kind=8) :: cdacb(nbt,imchmx,2,nrct)
    real(kind=8) :: cdrs(ndrsmx)
    real(kind=8) :: csigma(imchmx,2,nrctmx)
    real(kind=8) :: csts(nstsmx)
    real(kind=8) :: delvec(kmax)
    real(kind=8) :: dlogxw(nbtmax)
    real(kind=8) :: dvjdte(nrct)
    real(kind=8) :: eact(imchmx,2,nrctmx)
    real(kind=8) :: fkrc(nrctmx)
    real(kind=8) :: gmmatr(kmax,kmax)
    real(kind=8) :: hact(imchmx,2,nrctmx)
    real(kind=8) :: mmmatr(nrct,kmax)
    real(kind=8) :: morr(nrctmx)
    real(kind=8) :: mwtrc(nrctmx)
    real(kind=8) :: rk(imchmx,2,nrctmx)
    real(kind=8) :: rkb(imchmx,2,nrctmx)
    real(kind=8) :: rreac1(nrctmx)
    real(kind=8) :: rrelr1(nrctmx)
    real(kind=8) :: rrxfi1(imchmx,nrctmx)
    real(kind=8) :: rxbar(iktmax,nxrtmx)
    real(kind=8) :: sfcar(nrctmx)
    real(kind=8) :: sgmatr(nrct,kmax)
    real(kind=8) :: ssfcar(nrctmx)
    real(kind=8) :: ttk(nttkmx)
    real(kind=8) :: xlks(nstmax)
    real(kind=8) :: xxmatr(kmax,nrct)
    real(kind=8) :: xymatr(nrct,nrct1)

    real(kind=8) :: al10
    real(kind=8) :: eps100
    real(kind=8) :: rirec1
    real(kind=8) :: rtcnst
    real(kind=8) :: tempc
    real(kind=8) :: tempcb
    real(kind=8) :: tempk
    real(kind=8) :: whcfac
    real(kind=8) :: xi1

    ! Local variable declarations.
    integer :: i
    integer :: id
    integer :: idj
    integer :: idn
    integer :: ik
    integer :: j
    integer :: k
    integer :: kcol
    integer :: kk
    integer :: krow
    integer :: n
    integer :: nb
    integer :: np
    integer :: nrc
    integer :: nr1
    integer :: nr2
    integer :: ns
    integer :: nxr

    real(kind=8) :: ax
    real(kind=8) :: cx
    real(kind=8) :: cxs
    real(kind=8) :: dxy
    real(kind=8) :: rt2
    real(kind=8) :: sx
    real(kind=8) :: vx
    real(kind=8) :: xx

    real(kind=8) :: coefst

    ! Zero the matrix J[r] (armatr).
    do j = 1,nrct1
        do k = 1,nrct1
            armatr(j,k) = 0.
        end do
    end do

    ! Calculate the matrix M (mmmatr). This incorporates the
    ! dependency of rate laws on the thermodynamic activities of
    ! species (currently all such species must be of type aqueous).
    call gmmmat(act,afrc1,cdac,cdacb,cdrs,csigma,eps100,fkrc,idirec,iindx1,iktmax,imchmx,imech,jcode,jreac,kbt,kdim,kmax,mmmatr,nbasp,nbt,nbtmax,ncmpr,ndac,ndacb,ndact,ndctmx,ndrs,ndrsmx,ndrsr,noutpt,nptmax,nrct,nrctmx,nrk,nrndex,nstmax,nttyo,nxridx,nxrtmx,rk,rtcnst,rxbar,sfcar,ureac,xlks)

    ! Now calculate the matrix SIGMA (sgmatr) from matrix M and the
    ! array W-squiggle (dlogxw) array.
    do nrc = 1,nrct
        do k = 1,kmax
            sgmatr(nrc,k) = 0.
        end do
    end do

    do nrc = 1,nrct
        if (jreac(nrc) .le. 0) then
            sx = 0.

            do kk = 2,kbt
                sx = sx + mmmatr(nrc,kk)
            end do

            sgmatr(nrc,1) = al10*(mmmatr(nrc,1)*dlogxw(1) - sx)

            do k = 2,kbt
                sgmatr(nrc,k) = al10*(mmmatr(nrc,1)*dlogxw(k)      + mmmatr(nrc,k))
            end do
        end if
    end do

    ! Calculate the matrix PHI (aimatr), which is the inverse of
    ! the Jacobian matrix J[z] used by the algebraic equation solver.
    call ginvrt(aimatr,delvec,gmmatr,ipivot,kdim,kmax)

    ! Emulate the matrix multiplication PHI*THETA. Put the result
    ! in a scratch matrix (xxmatr). Note that xxmatr must have zero
    ! rows for reactants that are not pure minerals or solid solutions.
    ! Those are currently the only types that may have rate dependencies
    ! on activities. That might change in the future. Also, in the row-
    ! column multiplication, contributions are included only over the
    ! index range (1, kbt), not (1, kdim), which would be the general
    ! case. That is because the SIGMA array has no non-zero columns
    ! beyond the kbt-th. That in turn results because all of the
    ! currently allowed rate dependencies on thermodynamic activities
    ! are restricted to those of aqueous species.
    do krow = 1,kmax
        do nrc = 1,nrct
            xxmatr(krow,nrc) = 0.
        end do
    end do

    do krow = 1,kdim
        do nrc = 1,nrct
            if (jreac(nrc) .le. 0) then
                if (jcode(nrc) .eq. 0) then
                    ! Pure mineral.
                    np = nrndex(nrc)
                    ns = ncmpr(1,np)
                    xx = 0.

                    do kcol = 1,kbt
                        nb = iindx1(kcol)
                        cx = coefst(csts,nsts,nstsmx,nstsr,nb,ns,nstmax)
                        xx = xx + aimatr(krow,kcol)*cx
                    end do

                    xxmatr(krow,nrc) = xx
                else if (jcode(nrc) .eq. 1) then
                    ! Solid solution.
                    nxr = nxridx(nrc)
                    np = nrndex(nrc)
                    nr1 = ncmpr(1,np)
                    nr2 = ncmpr(2,np)
                    xx = 0.

                    do kcol = 1,kbt
                        nb = iindx1(kcol)
                        cx = 0.
                        ik = 0

                        do ns = nr1,nr2
                            ik = ik + 1
                            cxs = coefst(csts,nsts,nstsmx,nstsr,nb,ns,nstmax)
                            cx = cx + rxbar(ik,nxr)*cxs
                        end do

                        xx = xx + aimatr(krow,kcol)*cx
                    end do

                    xxmatr(krow,nrc) = xx
                else
                    xxmatr(krow,nrc) = 0.
                end if
            end if
        end do
    end do

    ! Calculate the product SIGMA*PHI*THETA, storing the result in
    ! scratch matrix xymatr. This operation does not put anything
    ! in the last (nrct1-th) column of xymatr. That column is
    ! associated with explicit time dependency of the rate
    ! equations. It is only used if the r vector of "rates
    ! represented by finite differences" is composed of relative
    ! rates plus the inverse rate.
    do nrc = 1,nrct
        do j = 1,nrct1
            xymatr(nrc,j) = 0.
        end do
    end do

    do nrc = 1,nrct
        if (jreac(nrc) .le. 0) then
            do j = 1,nrct
                xx = 0.

                do k = 1,kbt
                    xx = xx + sgmatr(nrc,k)*xxmatr(k,j)
                end do

                xymatr(nrc,j) = xx
            end do
        end if
    end do

    ! Emulate the matrix addition SIGMA*PHI*THETA + ZETA.
    ! The ZETA matrix addresses the dependencies of rates on
    ! surface areas. This coordinates with EQ6\csfar.f, which
    ! evaluates surface area models. Here again, the last (nrct1-th)
    ! column of xymatr is not used.
    do nrc = 1,nrct
        if (jreac(nrc) .le. 0) then
            id = idirec(nrc)

            if (nrk(id,nrc) .eq. 1) then
                ! The rate law is TST-type, therefore there is a surface
                ! area dependency.
                if (jcode(nrc).eq.0 .or. jcode(nrc).eq.1) then
                    np = nrndex(nrc)
                    n = ndrsr(1,ns)
                    cx = cdrs(n)
                else
                    write (noutpt,"(' Illegal jcode value in garmat.f')")
                    write (nttyo,"(' Illegal jcode value in garmat.f')")
                    stop
                end if

                sx = 0.

                if (nsk(nrc) .eq. 0) then
                    ! Constant surface area.
                    continue
                else if (nsk(nrc) .eq. 1) then
                    ! Constant specific surface area.
                    if (sfcar(nrc) .gt. 0.) then
                        sx = rreac1(nrc)*cx*ssfcar(nrc)*mwtrc(nrc)/sfcar(nrc)
                    end if
                else if (nsk(nrc) .eq. 2) then
                    ! Constant particle number surface area growth law.
                    if (morr(nrc) .gt. 0.) then
                        sx = 2.*rreac1(nrc)*cx/(3.*morr(nrc))
                    end if
                else
                    write (noutpt,"(' Illegal nsk value in garmat.f')")
                    write (nttyo,"(' Illegal nsk value in garmat.f')")
                    stop
                end if

                xymatr(nrc,nrc) = xymatr(nrc,nrc) + sx
            end if
        end if
    end do

    ! Overlay for variable temperature. This coordinates with
    ! EQ6\evratc.f, which calculates rate constants as a function
    ! of temperature, and with EQ6\gtemp.f, which calculates the
    ! temperature as a specified function of either Xi or t.
    do nrc = 1,nrct
        dvjdte(nrc) = 0.
    end do

    if (jtemp .gt. 0) then
        ! Calculate temperature derivatives of the rates, following
        ! through the temperature dependency of the rate constants.
        ! This is needed only if the temperature is a function of
        ! Xi or t.
        rt2 = rtcnst*tempk

        do nrc = 1,nrct
            if (jreac(nrc) .le. 0) then
                id = idirec(nrc)

                if (nrk(id,nrc) .ge. 2) then
                    do i = 1,imech(id,nrc)
                        if (iact(i,id,nrc) .eq. 0) then
                            ! No temperature dependence.
                            continue
                        else if (iact(i,id,nrc) .eq. 1) then
                            ! Constant activation energy.
                            dvjdte(nrc) = dvjdte(nrc)            + rrxfi1(i,nrc)*eact(i,id,nrc)/rt2
                        else if (iact(i,id,nrc) .eq. 2) then
                            ! Constant activation enthalpy.
                            dvjdte(nrc) = dvjdte(nrc)            + rrxfi1(i,nrc)*eact(i,id,nrc)/rt2
                        end if
                    end do
                end if
            end if
        end do
    end if

    if (jtemp .eq. 0) then
        ! Constant temperature.
        continue
    else if (jtemp .eq. 1) then
        ! Temperature is a linear function of Xi.
        ! tempc = tempcb + ttk(1)*Xi
        do nrc = 1,nrct
            id = idirec(nrc)

            if (jreac(nrc).le.0 .and. nrk(id,nrc).ge.2) then
                dxy = dvjdte(nrc)*ttk(1)

                do j = 1,nrct
                    idj = idirec(j)

                    if (jreac(j).le.0 .and. nrk(idj,j).ge.2) then
                        xymatr(nrc,j) = xymatr(nrc,j) + dxy
                    end if
                end do
            end if
        end do
    else if (jtemp .eq. 2) then
        ! Temperature is a linear function of time.
        ! tempc = tempcb + ttk(1)*time1
        do nrc = 1,nrct
            id = idirec(nrc)

            if (jreac(nrc).le.0 .and. nrk(id,nrc).ge.2) then
                dxy = dvjdte(nrc)*ttk(1)
                xymatr(nrc,nrct1) = xymatr(nrc,nrct1) + dxy
            end if
        end do
    else if (jtemp .eq. 3) then
        ! Fluid mixing tracking.
        ! tempc = (tempcb*ttk(1) + Xi*ttk(2))/(Xi + ttk(1))
        do nrc = 1,nrct
            id = idirec(nrc)

            if (jreac(nrc).le.0 .and. nrk(id,nrc).ge.2) then
                dxy = dvjdte(nrc)*((ttk(2) - tempc)/(xi1 + ttk(1)))

                do j = 1,nrct
                    idj = idirec(j)

                    if (jreac(j).le.0 .and. nrk(idj,j).ge.2) then
                        xymatr(nrc,j) = xymatr(nrc,j) + dxy
                    end if
                end do
            end if
        end do
    end if

    ! Multiply xymatr by the scalar w (whcfac).
    do nrc = 1,nrct
        do j = 1,nrct1
            xymatr(nrc,j) = whcfac*xymatr(nrc,j)
        end do
    end do

    ! Emulate the matrix multiplication V*xymatr, which yields
    ! the equivalent of J[r] + I. Note that V is nrct1 by nrct,
    ! while xymatr is nrct by nrct1. Thus, J[r] and the I in question
    ! are nrct1 by nrct1.
    do nrc = 1,nrct
        id = idirec(nrc)

        if (jreac(nrc) .le. 0) then
            if (nrk(id,nrc) .eq. 1) then
                ! The current row of V corresponds to a reactant constrained
                ! by a specified relative rate. V(nrc,nrc) = 1.0 and all
                ! other elements of this row are zero.
                armatr(nrc,j) = xymatr(nrc,j)
            else if (nrk(id,nrc).ge.2) then
                ! The current row of V corresponds to a reactant constrained
                ! by an absolute rate.
                do j = 1,nrct
                    ax = 0.

                    do n = 1,nrct
                        idn = idirec(n)

                        if (jreac(n) .le. 0) then
                            if (nrk(idn,n) .ge. 2) then
                                if (n .eq. nrc) then
                                    vx = rirec1*(1.0 - rrelr1(nrc))
                                else
                                    vx = -rirec1*rrelr1(nrc)
                                end if

                                ax = ax + vx*xymatr(n,j)
                            end if
                        end if

                        armatr(nrc,j) = ax
                    end do
                end do
            end if
        end if
    end do

    ! Now get the nrct1-th row.
    vx = -rirec1**2

    do j = 1,nrct
        ax = 0.

        do n = 1,nrct
            idn = idirec(n)

            if (jreac(n) .le. 0) then
                if (nrk(idn,n) .ge. 2) then
                    ax = ax + vx*xymatr(n,j)
                end if
            end if

            armatr(nrct1,j) = ax
        end do
    end do

    ! Finish getting J[r] by substracting the unit matrix.
    do nrc = 1,nrct1
        armatr(nrc,nrc) = armatr(nrc,nrc) - 1.0
    end do
end subroutine garmat