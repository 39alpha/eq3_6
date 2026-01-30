subroutine ncmpex(acflg,act,actlg,cdrs,cegexs,cgexj,conc,conclg,cpgexs,egexjc,egexjf,egexs,eps100,fo2,fo2lg,fsort,fugac,fugalg,iern1,iern2,ietmax,ifrn1,ifrn2,igas,igstak,iindx1,ilrn1,ilrn2,imrn1,imrn2,istack,ixrn1,ixrn2,jcsort,jern1,jern2,jetmax,jflag,jgext,jgsort,jgstak,jjsort,jpflag,jsflag,jsitex,jssort,jstack,kbt,kdim,kelect,kmax,km1,ko2gaq,kwater,kxt,loph,losp,lsort,mgext,mrgexs,mtb,moph,mosp,narn1,narn2,nbasp,nbt,nbtmax,ncmpr,ndrs,ndrsmx,ndrsr,nelect,nern1,nern2,netmax,ngexsa,ngext,ngrn1,ngrn2,ngt,ngtmax,noutpt,no2gaq,nphasx,npt,nptmax,nst,nstmax,nttyo,omega,omeglg,press,qxbarw,q6mode,ugexj,ugexmo,uphase,uspec,xbar,xbarlg,xbarw,xbarwc,xbrwlc,xbrwlg,xlks,zchar,zgexj,zvclg1,zvec1)
    use iso_fortran_env, only: dp => real64
    !! This subroutine computes all parameters necessary to write the
    !! Jacobian matrix from the zvclg1 array. It thus "expands" the
    !! basis set variable data.
    !! This subroutine is called by:
    !!   EQLIB/newton.f
    !!   EQLIB/nrstep.f
    !!   EQ3NR/arrset.f
    !!   EQ6/eqcalc.f
    !!   EQ6/exivar.f
    !!   EQ6/optmzr.f
    !!   EQ6/path.f
    !! Principal input:
    !!   acflg  = array of logarithms of activity coefficients
    !!   cdrs   = array of reaction coefficients
    !!   ncmpr  = array giving the range in arrays corresponding to
    !!              species of those species which belong to a given
    !!              phase
    !!   ndrs   = array parallel to cdrs giving the index of the
    !!              corresponding species
    !!   ndrsr  = array giving the range in the cdrs/ndrs arrays
    !!              corresonding to the reaction associated with a
    !!              given species
    !!   press  = pressure, bars
    !!   qxbarw = flag controlling iterative improvement of xbarw in
    !!              the present subroutine:
    !!              .false. = no iterative improvement
    !!              .true.  = iterative improvement (relevant to EQ6 only)
    !!   q6mode = flag denoting usage for EQ3NR or EQ6:
    !!              .false. = EQ3NR
    !!              .true.  = EQ6NR
    !!   xbarlg = array of logarithms of mole fractions of species;
    !!            this is primarily an output of this subroutine, but
    !!            xbarlg(narn1), the mole fraction of water, is used
    !!            as an input if q6mode is .false.
    !!   xlks   = array of equilibrium constants
    !!   zvec1  = array of master variables
    !!   zvclg1 = array of logarithms of master variables
    !! Principal output:
    !!   act    = array of species activities
    !!   actlg  = array of logarithms of species activities
    !!   conc   = array of species concentrations
    !!   conclg = array of logarithms of species concentrations
    !!   fo2    = the oxygen fugacity
    !!   fo2lg  = logarithm of the oxygen fugacity
    !!   fugac  = array of gas fugacities
    !!   fugalg = array of logarithms of gas fugacities
    !!   moph   = array of numbers of moles of phases
    !!   loph   = array of logarithms of numbers of moles of phases
    !!   mosp   = array of numbers of moles of species
    !!   losp   = array of logarithms of numbers of moles of species
    !!   xbar   = array of mole fractions of species
    !!   xbarlg = array of logarithms of mole fractions of species
    !!              (but see above under 'principal input")
    implicit none

    ! Calling sequence variable declarations.
    integer :: ietmax
    integer :: jetmax
    integer :: kmax
    integer :: nbtmax
    integer :: ndrsmx
    integer :: netmax
    integer :: ngtmax
    integer :: nptmax
    integer :: nstmax

    integer :: noutpt
    integer :: nttyo

    integer :: igstak(ngtmax)
    integer :: iindx1(kmax)
    integer :: istack(nstmax)
    integer :: jcsort(nstmax)
    integer :: jern1(jetmax,netmax)
    integer :: jern2(jetmax,netmax)
    integer :: jflag(nstmax)
    integer :: jgext(netmax)
    integer :: jgsort(ngtmax)
    integer :: jgstak(ngtmax)
    integer :: jjsort(nstmax)
    integer :: jpflag(nptmax)
    integer :: jsflag(nstmax)
    integer :: jsitex(nstmax)
    integer :: jssort(nstmax)
    integer :: jstack(nstmax)
    integer :: nbasp(nbtmax)
    integer :: ncmpr(2,nptmax)
    integer :: ndrs(ndrsmx)
    integer :: ndrsr(2,nstmax)
    integer :: ngexsa(ietmax,jetmax,netmax)
    integer :: ngext(jetmax,netmax)
    integer :: nphasx(nstmax)

    integer :: iern1
    integer :: iern2
    integer :: ifrn1
    integer :: ifrn2
    integer :: igas
    integer :: ilrn1
    integer :: ilrn2
    integer :: imrn1
    integer :: imrn2
    integer :: ixrn1
    integer :: ixrn2
    integer :: kbt
    integer :: kdim
    integer :: kelect
    integer :: km1
    integer :: ko2gaq
    integer :: kwater
    integer :: kxt
    integer :: narn1
    integer :: narn2
    integer :: nbt
    integer :: nelect
    integer :: nern1
    integer :: nern2
    integer :: ngrn1
    integer :: ngrn2
    integer :: ngt
    integer :: no2gaq
    integer :: npt
    integer :: nst

    logical :: qxbarw
    logical :: q6mode

    character(len=48) :: uspec(nstmax)
    character(len=24) :: ugexmo(netmax)
    character(len=24) :: uphase(nptmax)
    character(len=8) :: ugexj(jetmax,netmax)

    real(kind=8) :: acflg(nstmax)
    real(kind=8) :: act(nstmax)
    real(kind=8) :: actlg(nstmax)
    real(kind=8) :: cdrs(ndrsmx)
    real(kind=8) :: cegexs(ietmax,jetmax,netmax)
    real(kind=8) :: cgexj(jetmax,netmax)
    real(kind=8) :: conc(nstmax)
    real(kind=8) :: conclg(nstmax)
    real(kind=8) :: cpgexs(ietmax,jetmax,netmax)
    real(kind=8) :: egexjc(jetmax,netmax)
    real(kind=8) :: egexjf(jetmax,netmax)
    real(kind=8) :: egexs(ietmax,jetmax,netmax)
    real(kind=8) :: fsort(ngtmax)
    real(kind=8) :: fugac(ngtmax)
    real(kind=8) :: fugalg(ngtmax)
    real(kind=8) :: loph(nptmax)
    real(kind=8) :: losp(nstmax)
    real(kind=8) :: lsort(nstmax)
    real(kind=8) :: mgext(jetmax,netmax)
    real(kind=8) :: mrgexs(ietmax,jetmax,netmax)
    real(kind=8) :: mtb(nbtmax)
    real(kind=8) :: moph(nptmax)
    real(kind=8) :: mosp(nstmax)
    real(kind=8) :: xbar(nstmax)
    real(kind=8) :: xbarlg(nstmax)
    real(kind=8) :: xlks(nstmax)
    real(kind=8) :: zchar(nstmax)
    real(kind=8) :: zgexj(jetmax,netmax)
    real(kind=8) :: zvclg1(kmax)
    real(kind=8) :: zvec1(kmax)

    real(kind=8) :: eps100
    real(kind=8) :: fo2
    real(kind=8) :: fo2lg
    real(kind=8) :: omega
    real(kind=8) :: omeglg
    real(kind=8) :: press
    real(kind=8) :: xbarw
    real(kind=8) :: xbarwc
    real(kind=8) :: xbrwlc
    real(kind=8) :: xbrwlg

    ! Local variable declarations.
    integer :: ix
    integer :: iy
    integer :: je
    integer :: k
    integer :: kcol
    integer :: n
    integer :: nb
    integer :: ne
    integer :: ng
    integer :: np
    integer :: nrr1
    integer :: nrr2
    integer :: nr1
    integer :: nr2
    integer :: ns
    integer :: nss
    integer :: nt
    integer :: nxbarw

    character(len=48) :: ux48
    character(len=48) :: uy48
    character(len=24) :: ux24
    character(len=24) :: uy24

    real(kind=8) :: ajx
    real(kind=8) :: ax
    real(kind=8) :: axw
    real(kind=8) :: axwfnc
    real(kind=8) :: axwmax
    real(kind=8) :: axwmxo
    real(kind=8) :: cx
    real(kind=8) :: cxs
    real(kind=8) :: cxw
    real(kind=8) :: cxx
    real(kind=8) :: dxw
    real(kind=8) :: fx
    real(kind=8) :: lcx
    real(kind=8) :: lx
    real(kind=8) :: mw
    real(kind=8) :: mx
    real(kind=8) :: sigmmc
    real(kind=8) :: sx
    real(kind=8) :: wconst
    real(kind=8) :: wcnstl
    real(kind=8) :: xx
    real(kind=8) :: xy
    real(kind=8) :: zx

    real(kind=8) :: coefdr
    real(kind=8) :: texp
    real(kind=8) :: tlg

    ! Note: the following statements don't really do anything except
    ! cause the compiler not to complain that igas, uspec, uphase, and
    ! press are not used.
    ix = igas
    iy = ix
    igas = iy

    ux48 = uspec(1)
    uy48 = ux48
    uspec(1) = uy48

    ux24 = uphase(1)
    uy24 = ux24
    uphase(1) = uy24

    xx = press
    xy = xx
    press = xy

    ! Initialize the master variables from their logarithms.
    ! Clear the concentrations, activities, and numbers of moles
    ! of the aqueous species.
    do ns = narn1,narn2
        conc(ns) = 0.
        act(ns) = 0.
        mosp(ns) = 0.
        xbar(ns) = 0.
    end do

    do ns = narn1,narn2
        conclg(ns) = -99999.
        actlg(ns) = -99999.
        losp(ns) = -99999.
        xbarlg(ns) = -99999.
    end do

    ! Clear the concentrations, activities, etc., of the generic ion
    ! exchange species. However, do not clear the number of moles
    ! variables (mosp). These will be used to initiate the expansion.
    do ns = nern1,nern2
        conc(ns) = 0.
        act(ns) = 0.
        xbar(ns) = 0.
    end do

    do ns = nern1,nern2
        conclg(ns) = -99999.
        actlg(ns) = -99999.
        losp(ns) = -99999.
        xbarlg(ns) = -99999.
    end do

    ! Clear the concentrations, activities, fugacities, and numbers
    ! of moles of the gas species.
    do ns = ngrn1,ngrn2
        conc(ns) = 0.
        act(ns) = 0.
        mosp(ns) = 0.
        xbar(ns) = 0.
    end do

    do ng = 1,ngt
        fugac(ng) = 0.
    end do

    do ns = ngrn1,ngrn2
        conclg(ns) = -99999.
        actlg(ns) = -99999.
        losp(ns) = -99999.
        xbarlg(ns) = -99999.
    end do

    do ng = 1,ngt
        fugalg(ng) = -99999.
    end do

    do k = 1,kdim
        zx = zvclg1(k)
        zvec1(k) = texp(zx)
    end do

    if (kwater .le. 0) then
        write (noutpt,1000)
        write (nttyo,1000)
1000 format(/' * Error - (EQLIB/ncmpex) Programming error trap:',' Have not yet',/7x,'implemented coding for the case of water',' not in the basis set.')

        stop
    end if

    if (q6mode) then
        ! In EQ6, the mass of solvent water varies.
        wconst = omega/zvec1(kwater)
        wcnstl = tlg(wconst)
    else
        ! In EQ3NR, the mass of solvent water is fixed at 1.0 kg.
        wconst = 1.0
        wcnstl = 0.
    end if

    ! Compute the number of moles, concentrations, and activities of
    ! the basis species. The following loop assumes that all such
    ! species are normal aqueous solute species. Corrections must
    ! be made after this loop for the species water, aqueous O2(g),
    ! and aqueous e-, if they are in the basis set.
    do k = 1,kbt
        nb = iindx1(k)
        ns = nbasp(nb)
        lx = zvclg1(k)
        mx = zvec1(k)
        losp(ns) = lx
        mosp(ns) = mx

        if (ns.ge.narn1 .and. ns.le.narn2) then
            conclg(ns) = wcnstl + lx
            conc(ns) = wconst*mx
            ax = conclg(ns) + acflg(ns)
            actlg(ns) = ax
            act(ns) = texp(ax)
        else if (ns.ge.nern1 .and. ns.le.nern2) then
            conclg(ns) = wcnstl + lx
            conc(ns) = wconst*mx

            ! Mole fractions and activities of generic ion exchange species
            ! will be calculated later. For some exchange models (e.g.,
            ! Gapon), these quantities could be calculated here. However,
            ! for others (e.g., Vanselow), that is not possible, at least
            ! in the general case.
        else
            write (noutpt,1010)
            write (nttyo,1010)
1010 format(/' * Error - (EQLIB/ncmpex) Programming error trap:',' Have not yet',/7x,'implemented coding for the case of a',' basis species which is neither',/7x,'an aqueous species',' nor an ion-exchanger species.')

            stop
        end if
    end do

    ! The following is a return point for correcting the mole fraction
    ! of water in EQ6. Here nxbarw is an iteration counter.
    nxbarw = 0
    xbarw = xbarwc
    xbrwlg= xbrwlc
    axwmxo = 0.
200 continue

    ! Make corrections for water if it is a basis species. The
    ! following coding assumes that water is the first species in the
    ! aqueous solution (i.e., it has species index narn1).
    if (kwater .gt. 0) then
        if (q6mode) then
            xbarlg(narn1) = xbrwlg
            xbar(narn1) = xbarw
        else
            xbrwlg = zvclg1(kwater)
            xbarlg(narn1) = xbrwlg
            xbarw = texp(xbrwlg)
        end if

        ax = xbarlg(narn1) + acflg(narn1)
        actlg(narn1) = ax
        act(narn1) = texp(ax)
    end if

    ! Make corrections for aqueous O2(g), if it is a basis species.
    if (ko2gaq .gt. 0) then
        fo2lg = zvclg1(ko2gaq)
        fo2 = zvec1(ko2gaq)

        ! Note- act(no2gaq) and actlg(no2gaq) will contain the oxygen
        ! fugacity and its logarithm to simplify looping over mass action
        ! expressions. Technically, the aqueous O2(g) species has no
        ! activity.
        actlg(no2gaq) = fo2lg
        act(no2gaq) = fo2
    end if

    ! Make corrections for aqueous e-, if it is a basis species.
    if (kelect .gt. 0) then
        actlg(nelect) = zvclg1(kelect)
        act(nelect) = zvec1(kelect)
    end if

    ! Compute the mole fraction and activity of water, if this
    ! species is not a basis species.
    if (kwater .le. 0) then
        nr1 = ndrsr(1,narn1)
        nr2 = ndrsr(2,narn1)
        cxs = cdrs(nr1)
        cxx = -xlks(narn1) + cxs*acflg(narn1)

        do n = nr1 + 1,nr2
            nss = ndrs(n)
            cxx = cxx + cdrs(n)*actlg(narn1)
        end do

        cxx = -cxx/cxs
        xbarlg(narn1) = cxx
        xbar(narn1) = texp(cxx)
        ax = cxx + acflg(narn1)
        actlg(narn1) = ax
        act(narn1) = texp(ax)
    end if

    ! Compute concentrations and activities of the non-basis
    ! aqueous solute species.
    do ns = narn1 + 1,narn2
        if (jflag(ns) .eq. 30) then
            if (jsflag(ns) .le. 0) then
                nr1 = ndrsr(1,ns)
                nr2 = ndrsr(2,ns)
                cxs = cdrs(nr1)
                cxx = -xlks(ns) + cxs*acflg(ns)

                do n = nr1 + 1,nr2
                    nss = ndrs(n)
                    cxx = cxx + cdrs(n)*actlg(nss)
                end do

                cxx = -cxx/cxs
                conclg(ns) = cxx
                conc(ns) = texp(cxx)
                ax = cxx + acflg(ns)
                actlg(ns) = ax
                act(ns) = texp(ax)
            end if
        end if
    end do

    if (kwater.gt.0 .and. q6mode .and. qxbarw) then
        ! In EQ6 mode, check to see that the mole fraction of water is
        ! sufficiently well determined. If not, iteratively improve it.
        ! This is done using a 1-variable Newton-Raphson method. This
        ! is sometimes necessary because (1) to calculate the
        ! concentrations of the dependent aqueous species, one may have
        ! to have the mole fraction of water (xbarw) and (2) vice
        ! versa. The process is initiated by using an "old" value
        ! of the mole fraction of water (that of xbarwc on entering
        ! the present subroutine).
        ! Recalculate the mole fraction of water using the
        ! updated concentrations of all dependent aqueous species.
        ! The new value is stored as "xbarwc".
        ! Note: here the jcsort array, which provides for a sorted
        ! summation to calculate Sigma m (sigmmc), has not been updated.
        ! This is extremely unlikely to cause a problem here. An unsorted
        ! calculation would almost certainly be adequate.
        call csigm(conc,jcsort,narn1,narn2,nstmax,sigmmc)
        xbarwc = omega/(omega + sigmmc)

        ! Calculate the alpha residual (axw) and the max norm (axwmax).
        axw = xbarwc - xbarw
        axwmax = abs(axw)

        ! Calculate the improvement function (axwfnc).
        axwfnc = 0.

        if (axwmxo .gt. 0.) then
            axwfnc = (axwmxo -axwmax)/axwmxo
        end if

        axwmxo = axwmax

        if (nxbarw .le. 0) then
            write (noutpt,1200)
        end if

1200 format(/1x)

        write (noutpt,1210) nxbarw,xbarw,axwmax,axwfnc
1210 format(2x,'iter= ',i3,2x,'xbarw= ',f12.10,2x,'axwmax= ',1pe10.3,2x,'axwfnc= ',e10.3)

        if (axwmax .gt. eps100) then
            if (nxbarw .le. 20) then
                ! Calculate an improved value of xbarw.
                nxbarw = nxbarw + 1
                sx = 0.

                do ns = narn1 + 1,narn2
                    if (jflag(ns) .eq. 30) then
                        if (jsflag(ns) .le. 0) then
                            nr1 = ndrsr(1,ns)
                            nr2 = ndrsr(2,ns)
                            cxs = cdrs(nr1)

                            ! Calling sequence substitutions:
                            !   narn1 for nse
                            cxw = coefdr(cdrs,ndrs,ndrsmx,ndrsr,narn1,ns,nstmax)
                            sx = sx + cxw*conc(ns)/cxs
                        end if
                    end if
                end do

                ! Here ajx is the 1 x 1 Jacobian, and dxw is (initially) the
                ! raw correction term.
                ajx = (xbarwc/omega)*sx - 1.0
                dxw = -axw/ajx

                ! Apply limits to the correction.
                xx = xbarw + dxw

                if (xx .gt. 1.0) then
                    dxw = 0.5*(1.0 - xbarw)
                end if

                if (xx .le. 0.0) then
                    dxw = -0.5*xbarw
                end if

                dxw = min(dxw,0.05_dp)
                dxw = max(dxw,-0.05_dp)

                ! Make the correction and try again.
                xbarw = xbarw + dxw
                xbrwlg = tlg(xbarw)
                go to 200
            end if
        end if
    end if

    ! Make corrections for water. The following coding assumes that
    ! water is the first species in the aqueous solution (i.e., it
    ! has species index narn1).
    conc(narn1) = omega
    conclg(narn1) = omeglg

    ! Make corrections for aqueous O2(g).
    if (no2gaq .gt. 0) then
        conc(no2gaq) = 0.
        conclg(no2gaq) = -99999.
    end if

    ! Make corrections for aqueous e-.
    if (nelect .gt. 0) then
        conc(nelect) = 0.
        conclg(nelect) = -99999.
    end if

    ! Compute the numbers of moles of the non-basis aqueous solute
    ! species.
    do ns = narn1 + 1,narn2
        if (jflag(ns).eq.30 .and. jsflag(ns).le.0) then
            losp(ns) = conclg(ns) - wcnstl
            mosp(ns) = conc(ns)/wconst
        end if
    end do

    ! Make corrections for water. The following coding assumes that
    ! water is the first species in the aqueous solution (i.e., it
    ! has species index narn1).
    if (.not.q6mode) then
        mosp(narn1) = omega
        losp(narn1) = omeglg
    end if

    ! Make corrections for aqueous O2(g).
    if (no2gaq .gt. 0) then
        mosp(no2gaq) = 0.
        losp(no2gaq) = -99999.
    end if

    ! Make corrections for aqueous e-.
    if (nelect .gt. 0) then
        mosp(nelect) = 0.
        losp(nelect) = -99999.
    end if

    ! Compute the mole fractions and activities of the basis and
    ! non-basis ion-exchanger species, and the numbers of moles
    ! of the non-basis species.
    call ncmpve(acflg,act,actlg,cdrs,cgexj,eps100,iern1,iern2,ietmax,jern1,jern2,jetmax,jflag,jgext,jsflag,losp,mgext,moph,mosp,nbasp,nbt,nbtmax,ndrs,ndrsmx,ndrsr,netmax,noutpt,nptmax,nstmax,nttyo,ugexj,uphase,uspec,xbar,xbarlg,xlks)

    ! Compute the concentrations (mol/kg.H2O) and numbers of moles of
    ! the non-basis ion-exchanger species.
    do np = iern1,iern2
        if (moph(np) .gt. 0.) then
            ne = np - iern1 + 1

            do je = 1,jgext(ne)
                nrr1 = jern1(je,ne)
                nrr2 = jern2(je,ne)

                do ns = nrr1,nrr2
                    if (jflag(ns).eq.30 .and. jsflag(ns).le.0) then
                        lcx = wcnstl + losp(ns)
                        cx = wconst*mosp(ns)
                        conclg(ns) = lcx
                        conc(ns) = cx
                    end if
                end do
            end do
        end if
    end do

    ! Compute the equivalent fractions (egexs) and mole ratios (mrgexs)
    ! of exchanger species of generic ion exchanger phases.
    call gegexs(cegexs,cgexj,egexjc,egexjf,egexs,iern1,iern2,ietmax,jern1,jetmax,jgext,moph,mosp,mrgexs,netmax,ngexsa,ngext,noutpt,nptmax,nstmax,nttyo,zchar,zgexj)

    if (q6mode) then
        ! Compute the numbers of moles and concentrations of the
        ! species belonging to non-aqueous phases present in the
        ! equilibrium sytem (ES).
        if (kxt .ge. km1) then
            do kcol = km1,kxt
                ns = iindx1(kcol)
                losp(ns) = zvclg1(kcol)
                mosp(ns) = zvec1(kcol)
                conclg(ns) = wcnstl + losp(ns)
                conc(ns) = wconst*mosp(ns)
            end do
        end if
    end if

    ! Sort species according to log masses.
    call sortsp(iern1,iern2,istack,jcsort,jern1,jern2,jgext,jsitex,jetmax,jjsort,jssort,jstack,losp,lsort,ncmpr,nern1,nern2,netmax,noutpt,nphasx,npt,nptmax,nst,nstmax,nttyo)

    ! Compute the number of moles of water. The following coding
    ! assumes that the aqueous solution is the first phase (i.e.,
    ! has phase index 1).
    moph(1) = 0.
    loph(1) = -99999.
    mw = 0

    do nss = narn1,narn2
        ns = jcsort(nss)
        mw = mw + mosp(ns)
    end do

    moph(1) = mw
    loph(1) = tlg(mw)

    if (q6mode) then
        ! Compute the mole fractions of the aqueous species. The following
        ! coding assumes that water is the first species in the aqueous
        ! solution (i.e., has species index narn1).
        if (mw .gt. 0.) then
            do ns = narn1,narn2
                xx = mosp(ns)/mw
                xbar(ns) = xx
                xbarlg(ns) = tlg(xx)
            end do

            xbarw = xbar(narn1)
            xbrwlg = xbarlg(narn1)
        end if
    end if

    if (q6mode) then
        ! Compute the number of moles of non-aqueous phases. The following
        ! coding assumes that the aqueous solution is the first phase
        ! (i.e., has phase index 1).
        ! Skip the exchanger phases. Start with the pure liquid phases.
        do np = ilrn1,ilrn2
            moph(np) = 0.
            loph(np) = -99999.

            if (jpflag(np) .le. 0) then
                ns = ncmpr(1,np)
                mx = mosp(ns)
                moph(np) = mx
                loph(np) = tlg(mx)
            end if
        end do

        ! Pure mineral phases.
        do np = imrn1,imrn2
            moph(np) = 0.
            loph(np) = -99999.

            if (jpflag(np) .le. 0) then
                ns = ncmpr(1,np)
                mx = mosp(ns)
                moph(np) = mx
                loph(np) = tlg(mx)
            end if
        end do

        ! Fixed fugacity phases.
        do np = ifrn1,ifrn2
            moph(np) = 0.
            loph(np) = -99999.

            if (jpflag(np) .le. 0) then
                ns = ncmpr(1,np)
                mx = mosp(ns)
                moph(np) = mx
                loph(np) = tlg(mx)
            end if
        end do

        ! Solid solutions.
        do np = ixrn1,ixrn2
            moph(np) = 0.
            loph(np) = -99999.

            if (jpflag(np) .le. 0) then
                mx = 0
                nr1 = ncmpr(1,np)
                nr2 = ncmpr(2,np)

                do nss = nr1,nr2
                    ns = jcsort(nss)
                    mx = mx + mosp(ns)
                end do

                moph(np) = mx
                loph(np) = tlg(mx)
            end if
        end do
    end if

    if (q6mode) then
        ! Clear the activities of species belonging to non-aqueous
        ! solution phases. The following coding assumes that the aqueous
        ! solution is the first phase (i.e., has the phase index 1).
        do np = ixrn1,ixrn2
            nr1 = ncmpr(1,np)
            nr2 = ncmpr(2,np)

            do ns = nr1,nr2
                act(ns) = 0.
                actlg(ns) = -99999.
            end do
        end do
    end if

    if (q6mode) then
        ! Compute the mole fractions and activities of the species which
        ! belong to non-aqueous solution phases in the equilibrium system.
        do np = ixrn1,ixrn2
            mx = moph(np)

            if (mx .gt. 0.) then
                nr1 = ncmpr(1,np)
                nr2 = ncmpr(2,np)

                do ns = nr1,nr2
                    xx = mosp(ns)/mx
                    xbar(ns) = xx
                    xbarlg(ns) = tlg(xx)
                    ax = xbarlg(ns) + acflg(ns)
                    actlg(ns) = ax
                    act(ns) = texp(ax)
                end do
            end if
        end do
    end if

    ! Compute the fugacities of the gas species.
    do ns = ngrn1,ngrn2
        ng = ns - ngrn1 + 1

        if (jsflag(ns) .lt. 2) then
            nr1 = ndrsr(1,ns)
            nr2 = ndrsr(2,ns)
            nt = nr2 - nr1 + 1

            if (nt .ge. 2) then
                cxs = cdrs(nr1)
                fx = -xlks(ns)

                do n = nr1 + 1,nr2
                    nss = ndrs(n)
                    cxx = cdrs(n)
                    fx = fx + cxx*actlg(nss)
                end do

                fx = -fx/cxs
                fugalg(ng) = fx
                fugac(ng) = texp(fx)
            else
                fugalg(ng) = actlg(ns)
                fugac(ng) = act(ns)
            end if
        end if
    end do

    ! Sort all gas species according to equilibrium fugacities.
    ! Put their indices in sorted order in the jgsort array.
    ! Calling sequence substitutions:
    !   fsort for asort
    !   fugac for aval
    !   jgsort for jsort
    !   igstak for istack
    !   jgstak for jstack
    !   ngtmax for nmax
    !   ngt for nval
    ! Caution: the jgsort array from the last call is recycled as a
    ! good starting point. Set jgsort(1) to 0 to make a sort starting
    ! from scratch.
    call qsortw(fsort,fugac,igstak,jgsort,jgstak,ngtmax,noutpt,nttyo,ngt)
end subroutine ncmpex
