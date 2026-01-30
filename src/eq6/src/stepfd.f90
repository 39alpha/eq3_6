subroutine stepfd(acflg,acflg0,affp0,affp,afrc0,afrc1,aw0,aw1,delxi,dxsm00,eh0,eh1,emop,emop0,emos,emos0,fdafm1,fdaf0,fdarm1,fdar0,fdawm1,fdaw0,fdehm1,fdeh0,fdlim,fdo2m1,fdo20,fdpem1,fdpe0,fdphm1,fdph0,fdrem1,fdre0,fdrim1,fdri0,fdrrm1,fdrr0,fdsem1,fdse0,fdzvm1,fdzv0,fje,fje0,fo2lg0,fo2lg1,fxi,fxi0,iemop,iemop0,iemos,iemos0,iindx0,iindx1,iodb,iopt,ipndx0,ipndx1,jcode,jreac,jreac0,jpflag,kdim,kdim0,kmax,km1,km10,kmt,kmt0,kord,kx1,kx10,kxt,kxt0,modr,modr0,moph,moph0,morr,morr0,mosp,mosp0,mtb,mtb0,nbt,nbtmax,ncmpe,ncmpe0,nodbmx,noptmx,nordmx,noutpt,npet,npetmx,npet0,npt,nptmax,npts,nrct,nrctmx,nrd1mx,nrndex,nset,nsetmx,nset0,nstmax,ph0,ph1,qredox,rirec0,rirec1,rreac0,rreac1,rrelr0,rrelr1,sfcar,sfcar0,sigmam,sigmm0,tempc,tempc0,time0,time1,uzvec0,uzvec1,xi0,xi1,xim1,xirct,xirct0,zvclg0,zvclg1,zvec0,zvec1)
    !! This subroutine saves information at the current point and
    !! computes finite differences at the current point of reaction
    !! progress.
    !! This subroutine is called by:
    !!   EQ6/path.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: kmax
    integer :: nbtmax
    integer :: nodbmx
    integer :: noptmx
    integer :: nordmx
    integer :: npetmx
    integer :: nptmax
    integer :: nrctmx
    integer :: nrd1mx
    integer :: nsetmx
    integer :: nstmax

    integer :: noutpt

    integer :: iemop(npetmx)
    integer :: iemop0(npetmx)
    integer :: iemos(nsetmx)
    integer :: iemos0(nsetmx)
    integer :: iindx0(kmax)
    integer :: iindx1(kmax)
    integer :: iodb(nodbmx)
    integer :: iopt(noptmx)
    integer :: ipndx0(kmax)
    integer :: ipndx1(kmax)
    integer :: jcode(nrctmx)
    integer :: jreac(nrctmx)
    integer :: jreac0(nrctmx)
    integer :: jpflag(nptmax)
    integer :: ncmpe(2,npetmx)
    integer :: ncmpe0(2,npetmx)
    integer :: nrndex(nrctmx)

    integer :: kdim
    integer :: kdim0
    integer :: km1
    integer :: km10
    integer :: kmt
    integer :: kmt0
    integer :: kord
    integer :: kx1
    integer :: kx10
    integer :: kxt
    integer :: kxt0
    integer :: nbt
    integer :: npet
    integer :: npet0
    integer :: npt
    integer :: npts
    integer :: nrct
    integer :: nset
    integer :: nset0

    logical :: qredox

    character(len=48) :: uzvec0(kmax)
    character(len=48) :: uzvec1(kmax)

    real(kind=8) :: acflg(nstmax)
    real(kind=8) :: acflg0(nstmax)
    real(kind=8) :: affp(nptmax)
    real(kind=8) :: affp0(nptmax)
    real(kind=8) :: afrc0(nrctmx)
    real(kind=8) :: afrc1(nrctmx)
    real(kind=8) :: dxsm00(nrd1mx)
    real(kind=8) :: emop(npetmx)
    real(kind=8) :: emop0(npetmx)
    real(kind=8) :: emos(nsetmx)
    real(kind=8) :: emos0(nsetmx)
    real(kind=8) :: fdafm1(nordmx,nptmax)
    real(kind=8) :: fdaf0(nordmx,nptmax)
    real(kind=8) :: fdarm1(nordmx,nrctmx)
    real(kind=8) :: fdar0(nordmx,nrctmx)
    real(kind=8) :: fdawm1(nordmx)
    real(kind=8) :: fdaw0(nordmx)
    real(kind=8) :: fdehm1(nordmx)
    real(kind=8) :: fdeh0(nordmx)
    real(kind=8) :: fdo2m1(nordmx)
    real(kind=8) :: fdo20(nordmx)
    real(kind=8) :: fdpem1(nordmx,npetmx)
    real(kind=8) :: fdpe0(nordmx,npetmx)
    real(kind=8) :: fdphm1(nordmx)
    real(kind=8) :: fdph0(nordmx)
    real(kind=8) :: fdrem1(nordmx,nrctmx)
    real(kind=8) :: fdre0(nordmx,nrctmx)
    real(kind=8) :: fdrim1(nrd1mx)
    real(kind=8) :: fdri0(nrd1mx)
    real(kind=8) :: fdrrm1(nrd1mx,nrctmx)
    real(kind=8) :: fdrr0(nrd1mx,nrctmx)
    real(kind=8) :: fdsem1(nordmx,nsetmx)
    real(kind=8) :: fdse0(nordmx,nsetmx)
    real(kind=8) :: fdzvm1(nrd1mx,kmax)
    real(kind=8) :: fdzv0(nrd1mx,kmax)

    real(kind=8) :: modr(nrctmx)
    real(kind=8) :: modr0(nrctmx)
    real(kind=8) :: moph(nptmax)
    real(kind=8) :: moph0(nptmax)
    real(kind=8) :: morr(nrctmx)
    real(kind=8) :: morr0(nrctmx)
    real(kind=8) :: mosp(nstmax)
    real(kind=8) :: mosp0(nstmax)
    real(kind=8) :: mtb(nbtmax)
    real(kind=8) :: mtb0(nbtmax)
    real(kind=8) :: rreac0(nrctmx)
    real(kind=8) :: rreac1(nrctmx)
    real(kind=8) :: rrelr0(nrctmx)
    real(kind=8) :: rrelr1(nrctmx)
    real(kind=8) :: sfcar(nrctmx)
    real(kind=8) :: sfcar0(nrctmx)
    real(kind=8) :: xirct(nrctmx)
    real(kind=8) :: xirct0(nrctmx)
    real(kind=8) :: zvclg0(kmax)
    real(kind=8) :: zvclg1(kmax)
    real(kind=8) :: zvec0(kmax)
    real(kind=8) :: zvec1(kmax)

    real(kind=8) :: aw0
    real(kind=8) :: aw1
    real(kind=8) :: delxi
    real(kind=8) :: eh0
    real(kind=8) :: eh1
    real(kind=8) :: fdlim
    real(kind=8) :: fje
    real(kind=8) :: fje0
    real(kind=8) :: fo2lg0
    real(kind=8) :: fo2lg1
    real(kind=8) :: fxi
    real(kind=8) :: fxi0
    real(kind=8) :: ph0
    real(kind=8) :: ph1
    real(kind=8) :: rirec0
    real(kind=8) :: rirec1
    real(kind=8) :: sigmam
    real(kind=8) :: sigmm0
    real(kind=8) :: tempc
    real(kind=8) :: tempc0
    real(kind=8) :: time0
    real(kind=8) :: time1
    real(kind=8) :: xi0
    real(kind=8) :: xi1
    real(kind=8) :: xim1

    ! Local variable declarations.
    integer :: itrunc
    integer :: j
    integer :: k
    integer :: kcol
    integer :: n
    integer :: nmax
    integer :: np
    integer :: npe
    integer :: nrc
    integer :: nse

    real(kind=8) :: dfx
    real(kind=8) :: dfxl
    real(kind=8) :: dfxu

    if (npts .eq. 1) then
        ! Zero all finite differences.
        nmax = nrd1mx*kmax
        call initaz(fdzvm1,nmax)
        call initaz(fdzv0,nmax)

        nmax = nrd1mx*nrctmx
        call initaz(fdrrm1,nmax)
        call initaz(fdrr0,nmax)

        if (iopt(2) .gt. 0) then
            call initaz(fdrim1,nrd1mx)
            call initaz(fdri0,nrd1mx)
        end if

        call initaz(fdawm1,nordmx)
        call initaz(fdaw0,nordmx)
        call initaz(fdehm1,nordmx)
        call initaz(fdeh0,nordmx)
        call initaz(fdo2m1,nordmx)
        call initaz(fdo20,nordmx)
        call initaz(fdphm1,nordmx)
        call initaz(fdph0,nordmx)

        if (iopt(2) .gt. 0) then
            nmax = nordmx*nrctmx
            call initaz(fdrem1,nmax)
            call initaz(fdre0,nmax)
        end if

        nmax = nordmx*nptmax
        call initaz(fdafm1,nmax)
        call initaz(fdaf0,nmax)

        nmax = nordmx*npetmx
        call initaz(fdpem1,nmax)
        call initaz(fdpe0,nmax)

        nmax = nordmx*nsetmx
        call initaz(fdsem1,nmax)
        call initaz(fdse0,nmax)

        nmax = nordmx*nrctmx
        call initaz(fdarm1,nmax)
        call initaz(fdar0,nmax)

        go to 200
    end if

    ! Calculate the dxsm00 vector.
    do n = 1,kord
        j = kord - n + 1
        dxsm00(j + 1) = dxsm00(j) + delxi
    end do

    dxsm00(1) = delxi

    ! Overflow protection is automatically activated in the code blocks
    ! below in order to run on VAX machines with a small exponent range
    ! (+/- 38) for real*8.
    itrunc = kord

    ! Compute finite differences for the elements of the z vector.
    ! These are mostly the number of moles of basis species.
    do kcol = 1,kdim
        fdzv0(1,kcol) = (zvec1(kcol) - zvec0(kcol))/delxi

        do j = 2,kord + 1
            dfxu = fdlim*dxsm00(j)
            dfxl = -dfxu
            k = j - 1
            dfx = fdzv0(k,kcol) - fdzvm1(k,kcol)
            dfx = min(dfxu,dfx)
            dfx = max(dfxl,dfx)

            if (abs(dfx) .ge. dfxu) then
                itrunc = min(itrunc,k)
            end if

            fdzv0(j,kcol) = dfx/dxsm00(j)
        end do
    end do

    ! Compute finite differences for relative rates.
    do nrc = 1,nrct
        if (jreac(nrc).eq.1 .or. jreac(nrc).eq.2) then
            ! The relative rate must be a constant zero, because no
            ! reactant mass remains:
            !   jreac =  1: exhausted
            !   jreac =  2: saturated; any remaining reactant mass is
            !                 converted to the corresponding product
            !                 phase, so the "reactant" is effectively
            !                 exhausted.
            ! Hence the corresponding finite differences must also be
            ! zero. Avoid calculating finite differences from relative
            ! rate values that may by non-zero owing to convergence
            ! tolerances and the like.
            do j = 1,kord + 1
                fdrr0(j,nrc) = 0.
            end do
        else
            ! The relative rate is calculated for the reactant, which
            ! is actively reacting according to a rate law:
            !   jreac =  0: set to react
            !   jreac = -1: saturated, but the remaining reactant mass
            !                 continues to react irreversibly
            ! Hence calculate the finite differences.
            fdrr0(1,nrc) = (rrelr1(nrc) - rrelr0(nrc))/delxi

            do j = 2,kord + 1
                dfxu = fdlim*dxsm00(j)
                dfxl = -dfxu
                k = j - 1
                dfx = fdrr0(k,nrc) - fdrrm1(k,nrc)
                dfx = min(dfxu,dfx)
                dfx = max(dfxl,dfx)

                if (abs(dfx) .ge. dfxu) then
                    itrunc = min(itrunc,k)
                end if

                fdrr0(j,nrc) = dfx/dxsm00(j)
            end do
        end if
    end do

    if (iopt(2) .gt. 0) then
        ! Compute finite differences for the inverse rate.
        fdri0(1) = (rirec1 - rirec0)/delxi

        do j = 2,kord + 1
            dfxu = fdlim*dxsm00(j)
            dfxl = -dfxu
            k = j - 1
            dfx = fdri0(k) - fdrim1(k)
            dfx = min(dfxu,dfx)
            dfx = max(dfxl,dfx)

            if (abs(dfx) .ge. dfxu) then
                itrunc = min(itrunc,k)
            end if

            fdri0(j) = dfx/dxsm00(j)
        end do
    end if

    ! Compute finite differences for the pH.
    fdph0(1) = (ph1 - ph0)/delxi

    do j = 2,kord
        dfxu = fdlim*dxsm00(j)
        dfxl = -dfxu
        k = j - 1
        dfx = fdph0(k) - fdphm1(k)
        dfx = min(dfxu,dfx)
        dfx = max(dfxl,dfx)

        if (abs(dfx) .ge. dfxu) then
            itrunc = min(itrunc,k)
        end if

        fdph0(j) = dfx/dxsm00(j)
    end do

    ! Compute finite differences for the Eh.
    if (qredox) then
        fdeh0(1) = (eh1 - eh0)/delxi

        do j = 2,kord
            dfxu = fdlim*dxsm00(j)
            dfxl = -dfxu
            k = j - 1
            dfx = fdeh0(k) - fdehm1(k)
            dfx = min(dfxu,dfx)
            dfx = max(dfxl,dfx)

            if (abs(dfx) .ge. dfxu) then
                itrunc = min(itrunc,k)
            end if

            fdeh0(j) = dfx/dxsm00(j)
        end do
    end if

    ! Compute finite differences for the log fO2.
    if (qredox) then
        fdo20(1) = (fo2lg1 - fo2lg0)/delxi

        do j = 2,kord
            dfxu = fdlim*dxsm00(j)
            dfxl = -dfxu
            k = j - 1
            dfx = fdo20(k) - fdo2m1(k)
            dfx = min(dfxu,dfx)
            dfx = max(dfxl,dfx)

            if (abs(dfx) .ge. dfxu) then
                itrunc = min(itrunc,k)
            end if

            fdo20(j) = dfx/dxsm00(j)
        end do
    end if

    ! Compute finite differences for the activity of water.
    fdaw0(1) = (aw1 - aw0)/delxi

    do j = 2,kord
        dfxu = fdlim*dxsm00(j)
        dfxl = -dfxu
        k = j - 1
        dfx = fdaw0(k) - fdawm1(k)
        dfx = min(dfxu,dfx)
        dfx = max(dfxl,dfx)

        if (abs(dfx) .ge. dfxu) then
            itrunc = min(itrunc,k)
        end if

        fdaw0(j) = dfx/dxsm00(j)
    end do

    ! Compute finite differences for the affinities of the
    ! various phases.
    do np = 1,npt
        fdaf0(1,np) = (affp(np) - affp0(np))/delxi

        do j = 2,kord
            dfxu = fdlim*dxsm00(j)
            dfxl = -dfxu
            k = j - 1
            dfx = fdaf0(k,np) - fdafm1(k,np)
            dfx = min(dfxu,dfx)
            dfx = max(dfxl,dfx)

            if (abs(dfx) .ge. dfxu) then
                itrunc = min(itrunc,k)
            end if

            fdaf0(j,np) = dfx/dxsm00(j)
        end do
    end do

    ! Compute finite differences for the number of moles of phases
    ! present in the equilibrium system.
    do npe = 1,npet
        fdpe0(1,npe) = (emop(npe) - emop0(npe))/delxi

        do j = 2,kord
            dfxu = fdlim*dxsm00(j)
            dfxl = -dfxu
            k = j - 1
            dfx = fdpe0(k,npe) - fdpem1(k,npe)
            dfx = min(dfxu,dfx)
            dfx = max(dfxl,dfx)

            if (abs(dfx) .ge. dfxu) then
                itrunc = min(itrunc,k)
            end if

            fdpe0(j,npe) = dfx/dxsm00(j)
        end do
    end do

    ! Compute finite differences for the number of moles of
    ! selected species present in the equilibrium system. These
    ! species include H2O(l) only for the aqueous phase and
    ! all species of the non-aqueous phases.
    do nse = 1,nset
        fdse0(1,nse) = (emos(nse) - emos0(nse))/delxi

        do j = 2,kord
            dfxu = fdlim*dxsm00(j)
            dfxl = -dfxu
            k = j - 1
            dfx = fdse0(k,nse) - fdsem1(k,nse)
            dfx = min(dfxu,dfx)
            dfx = max(dfxl,dfx)

            if (abs(dfx) .ge. dfxu) then
                itrunc = min(itrunc,k)
            end if

            fdse0(j,nse) = dfx/dxsm00(j)
        end do
    end do

    if (iopt(2) .gt. 0) then
        ! Compute finite differences for reaction rates.
        do nrc = 1,nrct
            if (jreac(nrc).eq.1 .or. jreac(nrc).eq.2) then
                ! The reaction rate must be a constant zero, because no
                ! reactant mass remains:
                !   jreac =  1: exhausted
                !   jreac =  2: saturated; any remaining reactant mass is
                !                 converted to the corresponding product
                !                 phase, so the "reactant" is effectively
                !                 exhausted.
                ! Hence the corresponding finite differences must also be
                ! zero. Avoid calculating finite differences from reaction
                ! rate values that may by non-zero owing to convergence
                ! tolerances and the like.
                do j = 1,kord
                    fdre0(j,nrc) = 0.
                end do
            else
                ! The reaction rate is calculated for the reactant, which
                ! is actively reacting according to a rate law:
                !   jreac =  0: set to react
                !   jreac = -1: saturated, but the remaining reactant mass
                !                 continues to react irreversibly
                ! Hence calculate the finite differences.
                fdre0(1,nrc) = (rreac1(nrc) - rreac0(nrc))/delxi

                do j = 2,kord
                    dfxu = fdlim*dxsm00(j)
                    dfxl = -dfxu
                    k = j - 1
                    dfx = fdre0(k,nrc) - fdrem1(k,nrc)
                    dfx = min(dfxu,dfx)
                    dfx = max(dfxl,dfx)

                    if (abs(dfx) .ge. dfxu) then
                        itrunc = min(itrunc,k)
                    end if

                    fdre0(j,nrc) = dfx/dxsm00(j)
                end do
            end if
        end do
    end if

    ! Compute finite differences for the affinities of irreversible
    ! reactions.
    ! Zero the affinities of exhausted reactants which correspond
    ! to product minerals present in the ES. This will prevent
    ! non-zero values due to convergence tolerances, etc., from
    ! adversely affecting the corresponding finite differences for
    ! these reactants.
    do nrc = 1,nrct
        if (jreac(nrc) .eq. 1) then
            if (jcode(nrc).eq.0 .or. jcode(nrc) .eq. 1) then
                ! Note that this is done only for reactants which are
                ! pure minerals or solid solutions.
                np = nrndex(nrc)

                if (jpflag(np) .eq. -1) then
                    afrc0(nrc) = 0.
                    afrc1(nrc) = 0.
                end if
            end if
        end if
    end do

    ! Now compute the finite differences.
    do nrc = 1,nrct
        if (jreac(nrc).eq.-1 .or. jreac(nrc).eq.2) then
            ! The affinity is zero by definition:
            !   jreac = -1: saturated, but the remaining reactant mass
            !                 continues to react irreversibly
            !   jreac =  2: saturated; any remaining reactant mass is
            !                 converted to the corresponding product phase,
            !                 so the "reactant" is effectively exhausted.
            ! Hence the corresponding finite differences must also be zero.
            ! The corresponding computed affinities may be non-zero, owing
            ! to finite convergence tolerances. Avoid calculating finite
            ! differences from such computed affinity values.
            do j = 1,kord
                fdar0(j,nrc) = 0.
            end do
        else
            ! The affinity is generally a non-zero value calculated from
            ! a rate law:
            !   jreac =  0: set to react
            !   jreac =  1: exhausted
            ! Hence calculate the finite differences.
            fdar0(1,nrc) = (afrc1(nrc) - afrc0(nrc))/delxi

            do j = 2,kord
                dfxu = fdlim*dxsm00(j)
                dfxl = -dfxu
                k = j - 1
                dfx = fdar0(k,nrc) - fdarm1(k,nrc)
                dfx = min(dfxu,dfx)
                dfx = max(dfxl,dfx)

                if (abs(dfx) .ge. dfxu) then
                    itrunc = min(itrunc,k)
                end if

                fdar0(j,nrc) = dfx/dxsm00(j)
            end do
        end if
    end do

    if (itrunc .lt. kord) then
        if (iodb(1) .ge. 1) then
            write (noutpt,1000) itrunc
        end if

1000 format(/' * Note - (EQ6/stepfd) Cutting the order to ',i2,/7x,'to stay within the finite difference limit.')

        kord = itrunc
    end if

    ! Make the current point the new base point.
200 continue
    km10 = km1
    kmt0 = kmt
    kx10 = kx1
    kxt0 = kxt
    kdim0 = kdim

    xim1 = xi0
    xi0 = xi1

    time0 = time1
    tempc0 = tempc

    ph0 = ph1
    eh0 = eh1
    fo2lg0 = fo2lg1
    aw0 = aw1

    fje0 = fje
    fxi0 = fxi
    sigmm0 = sigmam

    call copyia(iindx1,iindx0,kdim)
    call copyia(ipndx1,ipndx0,kdim)
    call copyca(uzvec1,uzvec0,kdim)

    call copyaa(zvclg1,zvclg0,kdim)
    call copyaa(zvec1,zvec0,kdim)

    call copyaa(moph,moph0,nptmax)
    call copyaa(mosp,mosp0,nstmax)
    call copyaa(acflg,acflg0,nstmax)

    call copyaa(sfcar,sfcar0,nrct)
    call copyia(jreac,jreac0,nrct)
    call copyaa(xirct,xirct0,nrct)
    call copyaa(morr,morr0,nrct)
    call copyaa(modr,modr0,nrct)

    call copyaa(mtb,mtb0,nbt)

    call copyaa(affp,affp0,npt)

    npet0 = npet
    nset0 = nset

    call copyia(iemop,iemop0,npet)
    call copyaa(emop,emop0,npet)

    nmax = 2*npetmx
    call copyia(ncmpe,ncmpe0,nmax)

    call copyia(iemos,iemos0,nset)
    call copyaa(emos,emos0,nset)

    call copyaa(afrc1,afrc0,nrct)

    call copyaa(rrelr1,rrelr0,nrct)

    if (iopt(2) .gt. 0) then
        rirec0 = rirec1
        call copyaa(rreac1,rreac0,nrct)
    end if

    if (npts .gt. 2) then
        ! Save the old finite differences.
        nmax = nrd1mx*kmax
        call copyaa(fdzv0,fdzvm1,nmax)

        nmax = nrd1mx*nrctmx
        call copyaa(fdrr0,fdrrm1,nmax)

        if (iopt(2) .gt. 0) then
            call copyaa(fdri0,fdrim1,nrd1mx)
        end if

        call copyaa(fdaw0,fdawm1,nordmx)
        call copyaa(fdeh0,fdehm1,nordmx)
        call copyaa(fdo20,fdo2m1,nordmx)
        call copyaa(fdph0,fdphm1,nordmx)

        nmax = nordmx*nptmax
        call copyaa(fdaf0,fdafm1,nmax)

        nmax = nordmx*npetmx
        call copyaa(fdpe0,fdpem1,nmax)

        nmax = nordmx*nsetmx
        call copyaa(fdse0,fdsem1,nmax)

        nmax = nordmx*nrctmx
        call copyaa(fdar0,fdarm1,nmax)

        if (iopt(2) .gt. 0) then
            nmax = nordmx*nrctmx
            call copyaa(fdre0,fdrem1,nmax)
        end if
    end if
end subroutine stepfd
