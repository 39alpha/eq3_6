subroutine gmmmat(act,afrc1,cdac,cdacb,cdrs,csigma,eps100,fkrc,idirec,iindx1,iktmax,imchmx,imech,jcode,jreac,kbt,kdim,kmax,mmmatr,nbasp,nbt,nbtmax,ncmpr,ndac,ndacb,ndact,ndctmx,ndrs,ndrsmx,ndrsr,noutpt,nptmax,nrct,nrctmx,nrk,nrndex,nstmax,nttyo,nxridx,nxrtmx,rk,rtcnst,rxbar,sfcar,ureac,xlks)
    !! This subroutine calculates the matrix M (mmmatr) used in turn to
    !! calculate the Jacobian matrix J[r] (armatr) used by the higher-
    !! order (stiff) ODE integrator. The dependence of J[r] on exact
    !! kinetic rate law expressions (coded in EQ6\crrate.f) is handled
    !! through the M matrix. This matrix expresses the dependency of the
    !! rate laws on the thermodynamic activities of species, which
    !! presently must all be of type aqueous.
    !!   The present subroutine must be tightly coordinated with
    !! EQ6\ccrate.f in order to maintain proper consistency and thus
    !! the correctness of the calculated M matrix. If such coordination
    !! and consistency is compromised, the higher-order (stiff) ODE
    !! solver may fail to converge.
    !! This subroutine is called by:
    !!   EQ6/garmat.f
    !! Principal input:
    !! Principal output:
    !!   mmmatr = the matrix M, which is used to calculate the Jacobian
    !!              matrix J[r] used by the higher-order (stiff) ODE
    !!              integrator
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
    integer :: nxrtmx

    integer :: noutpt
    integer :: nttyo

    integer :: kbt
    integer :: kdim
    integer :: nbt
    integer :: nrct

    integer :: idirec(nrctmx)
    integer :: iindx1(kmax)
    integer :: imech(2,nrctmx)
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
    integer :: nxridx(nrctmx)

    character(len=24) :: ureac(nrctmx)

    real(kind=8) :: act(nstmax)
    real(kind=8) :: afrc1(nrctmx)
    real(kind=8) :: cdac(ndctmx,imchmx,2,nrctmx)
    real(kind=8) :: cdacb(nbt,imchmx,2,nrct)
    real(kind=8) :: cdrs(ndrsmx)
    real(kind=8) :: csigma(imchmx,2,nrctmx)
    real(kind=8) :: fkrc(nrctmx)
    real(kind=8) :: mmmatr(nrct,kmax)
    real(kind=8) :: rk(imchmx,2,nrctmx)
    real(kind=8) :: sfcar(nrctmx)
    real(kind=8) :: rxbar(iktmax,nxrtmx)
    real(kind=8) :: xlks(nstmax)

    real(kind=8) :: eps100
    real(kind=8) :: rtcnst

    ! Local variable declarations.
    integer :: i
    integer :: ik
    integer :: j2
    integer :: j3
    integer :: k
    integer :: n
    integer :: nb
    integer :: np
    integer :: nr1
    integer :: nr2
    integer :: nrc
    integer :: ns
    integer :: nse
    integer :: nxr

    integer :: ilnobl

    logical :: qovstp

    character(len=8) :: ux8

    real(kind=8) :: affr
    real(kind=8) :: aprod
    real(kind=8) :: cx
    real(kind=8) :: cxs
    real(kind=8) :: efx
    real(kind=8) :: fs
    real(kind=8) :: mijk
    real(kind=8) :: skfq
    real(kind=8) :: thetij
    real(kind=8) :: vij
    real(kind=8) :: xk
    real(kind=8) :: xlk

    real(kind=8) :: coefdr
    real(kind=8) :: texp

    ! First zero the matrix M.
    do nrc = 1,nrct
        do k = 1,kmax
            mmmatr(nrc,k) = 0.
        end do
    end do

    ! Loop on reactants (rows of M).
    do nrc = 1,nrct
        affr = afrc1(nrc)
        qovstp = affr.lt.0. .and. nrk(2,nrc).eq.0 .and. jreac(nrc).le.0

        if (jreac(nrc).le.0 .and.    (abs(affr).gt.eps100 .or. qovstp))then
            ! Have an active reactant. It is not exhausted or formally
            ! saturated as indicated by its jreac flag. Also, it either
            ! has a finite affinity or the overstep condition pertains.
            ! Compute the corresponding row of M. For a reactant not
            ! satisfying these conditions, the rate is taken to be zero,
            ! and the corresponding row of M should be left filled with
            ! zeroes.
            fs = fkrc(nrc)*sfcar(nrc)

            if (idirec(nrc) .eq. 1) then
                ! The rate law expression used follows the net forward rate
                ! format.
                if (nrk(1,nrc) .eq. 1) then
                    ! Have a specified relative rate (arbitrary kinetics).
                    ! Leave the current row of M filled with zeroes.
                    continue
                else if (nrk(1,nrc) .eq. 2) then
                    ! Have the TST-like rate equation (multi-term). Compute
                    ! the contribution of each successive term to the current
                    ! row of M and accumulate all contributions.
                    ! First get the equilibrium constant K (xk).
                    np = nrndex(nrc)
                    xlk = 0.

                    if (jcode(nrc) .eq. 0) then
                        ns = ncmpr(1,np)
                        xlk = xlks(ns)
                    else if (jcode(nrc) .eq. 1) then
                        nxr = nxridx(nrc)
                        ik = 0
                        nr1 = ncmpr(1,np)
                        nr2 = ncmpr(2,np)

                        do ns = nr1,nr2
                            ik = ik + 1
                            xlk = xlk + rxbar(ik,nxr)*xlks(ns)
                        end do
                    else
                        write (noutpt,"(' Illegal jcode value in gmmmat.f')")
                        write (nttyo,"(' Illegal jcode value in gmmmat.f')")
                        stop
                    end if

                    xk = texp(xlk)

                    ! Loop over terms in the rate law.
                    do i = 1,imech(1,nrc)
                        ! Calculate the kinetic activity product (aprod).
                        aprod = 1.

                        do n = 1,ndact(i,1,nrc)
                            ns = ndac(n,i,1,nrc)
                            aprod = aprod*act(ns)**cdac(n,i,1,nrc)
                        end do

                        skfq = rk(i,1,nrc)*fs*aprod
                        efx = affr/(csigma(i,1,nrc)*rtcnst)
                        vij = skfq*(1.0 - exp(-efx))

                        do k = 1,kbt
                            nb = iindx1(k)
                            nse = nbasp(nb)

                            if (jcode(nrc) .eq. 0) then
                                ns = ncmpr(1,np)
                                cx = coefdr(cdrs,ndrs,ndrsmx,ndrsr,nse,ns,nstmax)
                            else if (jcode(nrc) .eq. 1) then
                                nxr = nxridx(nrc)
                                cx = 0.
                                ik = 0
                                nr1 = ncmpr(1,np)
                                nr2 = ncmpr(2,np)

                                do ns = nr1,nr2
                                    ik = ik + 1
                                    cxs = coefdr(cdrs,ndrs,ndrsmx,ndrsr,nse,ns,nstmax)
                                    cx = cx + rxbar(ik,nxr)*cxs
                                end do
                            else
                                write (noutpt,"(' Illegal jcode value in gmmmat.f')")
                                write (nttyo,"(' Illegal jcode value in gmmmat.f')")
                                stop
                            end if

                            thetij = -skfq*exp(-efx)/csigma(i,1,nrc)
                            mijk = -vij*ndacb(nse,i,1,nrc) + cx*thetij
                            mmmatr(nrc,k) = mmmatr(nrc,k) + mijk
                        end do
                    end do
                else if (nrk(1,nrc) .eq. 3) then
                    ! Have a specified rate (no affinity dependence).
                    ! Leave the current row of M filled with zeroes.
                else
                    ! Unrecognized rate law code.
                    j2 = ilnobl(ureac(nrc))
                    write (ux8,'(i5)') nrk(1,nrc)
                    call lejust(ux8)
                    j3 = ilnobl(ux8)
                    write (noutpt,1000) ureac(nrc)(1:j2),ux8(1:j3)
                    write (nttyo,1000) ureac(nrc)(1:j2),ux8(1:j3)
1000 format(/' * Error - (EQ6/gmmmat) The reactant ',a,' has an',/7x,'unrecognized forward rate law code',' of ',a,'.')

                    stop
                end if
            else
                ! The rate law expression used follows the net backward rate
                ! format.
                if (nrk(2,nrc) .eq. 1) then
                    ! Have specified the case of instantaneous equilibrium.
                    ! For now, leave the current row of M filled with zeroes.
                    continue
                else if (nrk(2,nrc) .eq. 2) then
                    ! Have the TST-like rate equation (multi-term). Compute
                    ! the contribution of each successive term to the current
                    ! row of M and accumulate all contributions.
                    ! First get the equilibrium constant K (xk).
                    np = nrndex(nrc)
                    xlk = 0.

                    if (jcode(nrc) .eq. 0) then
                        ns = ncmpr(1,np)
                        xlk = xlks(ns)
                    else if (jcode(nrc) .eq. 1) then
                        nxr = nxridx(nrc)
                        ik = 0
                        nr1 = ncmpr(1,np)
                        nr2 = ncmpr(2,np)

                        do ns = nr1,nr2
                            ik = ik + 1
                            xlk = xlk + rxbar(ik,nxr)*xlks(ns)
                        end do
                    else
                        write (noutpt,"(' Illegal jcode value in gmmmat.f')")
                        write (nttyo,"(' Illegal jcode value in gmmmat.f')")
                        stop
                    end if

                    xk = texp(xlk)

                    ! Loop over terms in the rate law.
                    do i = 1,imech(2,nrc)
                        ! Calculate the kinetic activity product (aprod).
                        aprod = 1.

                        do n = 1,ndact(i,2,nrc)
                            ns = ndac(n,i,2,nrc)
                            aprod = aprod*act(ns)**cdac(n,i,2,nrc)
                        end do

                        skfq = rk(i,2,nrc)*fs*aprod
                        efx = affr/(csigma(i,2,nrc)*rtcnst)
                        vij = skfq*(1.0 - exp(efx))

                        do k = 1,kbt
                            nb = iindx1(k)
                            nse = nbasp(nb)

                            if (jcode(nrc) .eq. 0) then
                                ns = ncmpr(1,np)
                                cx = coefdr(cdrs,ndrs,ndrsmx,ndrsr,nse,ns,nstmax)
                            else if (jcode(nrc) .eq. 1) then
                                nxr = nxridx(nrc)
                                cx = 0.
                                ik = 0
                                nr1 = ncmpr(1,np)
                                nr2 = ncmpr(2,np)

                                do ns = nr1,nr2
                                    ik = ik + 1
                                    cxs = coefdr(cdrs,ndrs,ndrsmx,ndrsr,nse,ns,nstmax)
                                    cx = cx + rxbar(ik,nxr)*cxs
                                end do
                            else
                                write (noutpt,"(' Illegal jcode value in gmmmat.f')")
                                write (nttyo,"(' Illegal jcode value in gmmmat.f')")
                                stop
                            end if

                            thetij = -skfq*exp(efx)/(csigma(i,2,nrc))
                            mijk = -vij*ndacb(nse,i,2,nrc) + cx*thetij
                            mmmatr(nrc,k) = mmmatr(nrc,k) + mijk
                        end do
                    end do
                else if (nrk(2,nrc) .eq. 3) then
                    ! Have a specified rate (no affinity dependence).
                    ! Leave the current row of M filled with zeroes.
                else
                    ! Unrecognized rate law code.
                    j2 = ilnobl(ureac(nrc))
                    write (ux8,'(i5)') nrk(2,nrc)
                    call lejust(ux8)
                    j3 = ilnobl(ux8)
                    write (noutpt,1020) ureac(nrc)(1:j2),ux8(1:j3)
                    write (nttyo,1020) ureac(nrc)(1:j2),ux8(1:j3)
1020 format(/' * Error - (EQ6/gmmmat) The reactant ',a,' has an',/7x,'unrecognized backward rate law code',' of ',a,'.')

                    stop
                end if
            end if
        end if
    end do
end subroutine gmmmat