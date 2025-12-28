subroutine betas(acflg,actlg,afcnst,alpha,amtb,bbig,beta,betamx,bneg,cdrs,conc,conclg,coval,csts,eh,ehfac,fo2lg,ibetmx,iebal,iindx1,irdxc3,jcsort,jflag,jsflag,jssort,kbt,kdim,kelect,khydr,kmax,km1,ko2gaq,kwater,kxt,mtb,mosp,narn1,narn2,nbasp,nbtmax,ncosp,ndrs,ndrsmx,ndrsr,nelect,nern1,nern2,nhydr,noutpt,no2gaq,nredox,nst,nstmax,nsts,nstsmx,nstsr,ntfx,ntfxmx,ntfxt,nttyo,omega,qredox,q6mode,tfx,ubbig,ubneg,ubetmx,uspec,uzvec1,weight,xbrwlg,xlke,xlks,zchar)
    !! This subroutine computes the Newton-Raphson residual functions.
    !! This subroutine is called by:
    !!   EQLIB/newton.f
    !!   EQLIB/nrstep.f
    !!   EQ3NR/arrset.f
    !!   EQ6/eqcalc.f
    !!   EQ6/optmzr.f
    !! Principal input:
    !!   acflg  = array of log activity coefficients
    !!   q6mode = flag denoting usage for EQ3NR or EQ6:
    !!              .false. = EQ3NR
    !!              .true.  = EQ6NR
    !! Principal output:
    !!   alpha  = array of raw Newton-Raphson residual functions
    !!   beta   = array of normalized Newton-Raphson residual functions
    implicit none

    ! Calling sequence variable declarations.
    integer :: kmax
    integer :: nbtmax
    integer :: ndrsmx
    integer :: nstmax
    integer :: nstsmx
    integer :: ntfxmx

    integer :: noutpt
    integer :: nttyo

    integer :: iindx1(kmax)
    integer :: jcsort(nstmax)
    integer :: jflag(nstmax)
    integer :: jsflag(nstmax)
    integer :: jssort(nstmax)
    integer :: nbasp(nbtmax)
    integer :: ncosp(nbtmax)
    integer :: ndrs(ndrsmx)
    integer :: ndrsr(2,nstmax)
    integer :: nsts(nstsmx)
    integer :: nstsr(2,nstmax)
    integer :: ntfx(ntfxmx)

    integer :: ibetmx
    integer :: iebal
    integer :: irdxc3
    integer :: kbt
    integer :: kdim
    integer :: kelect
    integer :: khydr
    integer :: km1
    integer :: ko2gaq
    integer :: kwater
    integer :: kxt
    integer :: narn1
    integer :: narn2
    integer :: nelect
    integer :: nern1
    integer :: nern2
    integer :: nhydr
    integer :: no2gaq
    integer :: nredox
    integer :: nst
    integer :: ntfxt

    logical :: qredox
    logical :: q6mode

    character(len=48) :: uspec(nstmax)
    character(len=48) :: uzvec1(kmax)
    character(len=48) :: ubbig
    character(len=48) :: ubneg
    character(len=48) :: ubetmx

    real(kind=8) :: acflg(nstmax)
    real(kind=8) :: actlg(nstmax)
    real(kind=8) :: alpha(kmax)
    real(kind=8) :: amtb(nbtmax)
    real(kind=8) :: beta(kmax)
    real(kind=8) :: cdrs(ndrsmx)
    real(kind=8) :: conc(nstmax)
    real(kind=8) :: conclg(nstmax)
    real(kind=8) :: coval(nbtmax)
    real(kind=8) :: csts(nstsmx)
    real(kind=8) :: mtb(nbtmax)
    real(kind=8) :: mosp(nstmax)
    real(kind=8) :: tfx(ntfxmx)
    real(kind=8) :: weight(nstmax)
    real(kind=8) :: xlks(nstmax)
    real(kind=8) :: zchar(nstmax)

    real(kind=8) :: afcnst
    real(kind=8) :: bbig
    real(kind=8) :: betamx
    real(kind=8) :: bneg
    real(kind=8) :: eh
    real(kind=8) :: ehfac
    real(kind=8) :: fo2lg
    real(kind=8) :: omega
    real(kind=8) :: xbrwlg
    real(kind=8) :: xlke

    ! Local variable declarations.
    integer :: jfl
    integer :: jlen
    integer :: kcol
    integer :: krow
    integer :: n
    integer :: nb
    integer :: nr1
    integer :: nr2
    integer :: ns
    integer :: nsc
    integer :: nse
    integer :: nss
    integer :: ns1

    character(len=56) :: uspn56
    character(len=8) :: unone
    character(len=8) :: uptgas

    real(kind=8) :: abeta
    real(kind=8) :: af
    real(kind=8) :: alkc
    real(kind=8) :: atot
    real(kind=8) :: atx
    real(kind=8) :: ax
    real(kind=8) :: aznse
    real(kind=8) :: azns1
    real(kind=8) :: azsum
    real(kind=8) :: belect
    real(kind=8) :: bo2gaq
    real(kind=8) :: bx
    real(kind=8) :: bwater
    real(kind=8) :: cde
    real(kind=8) :: ctot
    real(kind=8) :: cx
    real(kind=8) :: c1
    real(kind=8) :: dx
    real(kind=8) :: fx
    real(kind=8) :: si
    real(kind=8) :: sigmmc
    real(kind=8) :: sigza
    real(kind=8) :: sigzc
    real(kind=8) :: sigzi
    real(kind=8) :: sigzm
    real(kind=8) :: sx
    real(kind=8) :: tx
    real(kind=8) :: xbarwc
    real(kind=8) :: xbrwlc
    real(kind=8) :: zp

    real(kind=8) :: coefst
    real(kind=8) :: tlg

    data unone  /'None    '/
    data uptgas /'Gas     '/

    do krow = 1,kdim
        alpha(krow) = 0.
        beta(krow) = 0.
    end do

    if (q6mode) then
        go to 210
    end if

    ! Compute residuals for EQ3NR.
    do krow = 1,kbt
        nb = iindx1(krow)
        nse = nbasp(nb)
        jfl = jflag(nse)

        if (nse.ge.narn1 .and. nse.le.narn2) then
            ! Aqueous species.
            if (krow.eq.kwater .and. jfl.eq.0) then
                ! The residual for water is based on the equation which
                ! defines the mole fraction of water.
                call csigm(conc,jcsort,narn1,narn2,nstmax,sigmmc)
                xbarwc = omega/(omega + sigmmc)
                xbrwlc = tlg(xbarwc)
                ax = xbrwlc - xbrwlg
                alpha(krow) = ax
                beta(krow) = ax
            else if (nb .eq. iebal) then
                ! Charge balance.
                call gszm(conc,jcsort,narn1,narn2,nstmax,sigza,sigzc,sigzi,sigzm,zchar)
                alpha(krow) = sigzi

                if (sigzm .gt. 0.) then
                    beta(krow) = sigzi/sigzm
                end if
            else if (jfl .eq. 17) then
                ! Log activity combination.
                ns1 = ncosp(nb)
                aznse = abs(zchar(nse))
                azns1 = abs(zchar(ns1))
                ax = coval(nb)/azns1
                zp = zchar(nse)*zchar(ns1)

                if (zp .lt. 0.) then
                    ax = ax - ( aznse/azns1 )*actlg(ns1)
                else
                    ax = ax + ( aznse/azns1 )*actlg(ns1)
                end if

                dx = ax - actlg(nse)
                alpha(krow) = dx
                beta(krow) = dx
            else if (jfl .eq. 18) then
                ! Mean log activity.
                ns1 = ncosp(nb)
                aznse = abs(zchar(nse))
                azns1 = abs(zchar(ns1))
                azsum = aznse + azns1
                ax = ( azsum*coval(nb) - aznse*actlg(ns1) )/azns1
                dx = ax/azns1 - actlg(nse)
                alpha(krow) = dx
                beta(krow) = dx
            else if (jfl .eq. 21) then
                ! pHCl.
                ns1 = ncosp(nb)
                aznse = abs(zchar(nse))
                azns1 = abs(zchar(ns1))
                ax = -coval(nb)/azns1
                zp = zchar(nse)*zchar(ns1)

                if (zp .lt. 0.) then
                    ax = ax - ( aznse/azns1 )*actlg(ns1)
                else
                    ax = ax + ( aznse/azns1 )*actlg(ns1)
                end if

                dx = ax - actlg(nse)
                alpha(krow) = dx
                beta(krow) = dx
            else if (jfl.eq.25 .or. jfl.eq.27) then
                ! Heterogeneous or homogeneous equilibrium with a specified
                ! species (computing fO2 from an aqueous redox couple is not
                ! done here).
                if (jfl .eq. 25) then
                    nsc = ncosp(nb)
                else
                    nsc = nse
                end if

                cx = xlks(nsc)
                nr1 = ndrsr(1,nsc)
                nr2 = ndrsr(2,nsc)

                do n = nr1,nr2
                    nss = ndrs(n)

                    if (nss .eq. nse) then
                        cde = cdrs(n)

                        if (nse .ne. no2gaq) then
                            cx = cx - cde*acflg(nse)
                        end if
                    else if (nss .eq. nsc) then
                        if (uspec(nsc)(25:32) .eq. uptgas(1:8)) then
                            cx = cx - cdrs(n)*coval(nb)
                        end if
                    else if (nss .eq. no2gaq) then
                        cx = cx - cdrs(n)*fo2lg
                    else
                        cx = cx - cdrs(n)*actlg(nss)
                    end if
                end do

                cx = cx/cde
                c1 = conclg(nse)

                if (nse .eq. no2gaq) then
                    c1 = fo2lg
                end if

                dx = cx - c1
                alpha(krow) = dx
                beta(krow) = dx
            else if (nse .eq. no2gaq) then
                ! log fO2.
                if (irdxc3 .eq. -1) then
                    ! Eh residual (Note- if a pe- value was input, it has been
                    ! converted to an Eh value by EQLIB/setup.f). Here it is
                    ! assumed that water is the first species in the aqueous
                    ! phase (has the index narn1).
                    fx = 4.*eh/ehfac  + xlke + 2.*actlg(narn1)        - 4.*actlg(nhydr)
                    dx = fx - fo2lg
                    alpha(krow) = dx
                    beta(krow) = dx
                else if (irdxc3 .eq. 1) then
                    ! Cross-linking (homogeneous aqueous redox) equilibrium.
                    cx = xlks(nredox)
                    nr1 = ndrsr(1,nredox)
                    nr2 = ndrsr(2,nredox)

                    do n = nr1,nr2
                        nss = ndrs(n)

                        if (nss .eq. no2gaq) then
                            cde = cdrs(n)
                        else
                            cx = cx - cdrs(n)*actlg(nss)
                        end if
                    end do

                    cx = cx/cde
                    dx = cx - fo2lg
                    alpha(krow) = dx
                    beta(krow) = dx
                else
                    ! Log fO2 is directly specified (irdxc3 = 0).
                    alpha(krow) = 0.
                    beta(krow) = 0.
                end if
            else if (jfl.ge.7 .and. jfl.le.11) then
                ! Alkalinity balance.
                call calk(alkc,conc,nstmax,ntfx,ntfxmx,ntfxt,tfx)
                atot = coval(nb)
                dx = alkc - atot
                alpha(krow) = dx
                beta(krow) = dx/atot
            else if (jfl .eq. 16) then
                ! Log activity.
                dx = conclg(nse) + acflg(nse) - coval(nb)
                alpha(krow) = dx
                beta(krow) = dx
            else if (jfl.eq.19 .or. jfl.eq.20) then
                ! pX (including pH).
                dx = conclg(nse) + acflg(nse) + coval(nb)
                alpha(krow) = dx
                beta(krow) = dx
            else if (jfl.eq.22 .or. jfl.eq.23) then
                ! pmX (including pmH).
                dx = conclg(nse) + coval(nb)
                alpha(krow) = dx
                beta(krow) = dx
            else if (jfl.ge.0 .and. jfl.le.3) then
                ! Mass balance.
                sx = 0.

                do nss = narn1,narn2
                    ns = jcsort(nss)
                    weight(ns) = coefst(csts,nsts,nstsmx,nstsr,nb,ns,nstmax)
                    sx = sx + weight(ns)*conc(ns)
                end do

                ! do nss = nern1,nern2
                !   ns = jcsort(nss)
                !   weight(ns) = coefst(csts,nsts,nstsmx,nstsr,nb,ns,nstmax)
                !   sx = sx + weight(ns)*conc(ns)
                ! enddo
                ctot = coval(nb)
                dx = sx - ctot
                alpha(krow) = dx
                beta(krow) = dx/ctot
            else
                ! Have found a bad jflag value.
                ! Calling sequence substitutions:
                !   uzvec1(krow) for unam48
                call fmspnm(jlen,uzvec1(krow),uspn56)
                write (noutpt,1000) jfl,uspn56(1:jlen)
                write (nttyo,1000) jfl,uspn56(1:jlen)
1000 format(/' * Error - (EQLIB/betas) Programming error trap:',/7x,'Have encountered a bad jflag value of ',i3,' for',/7x,'the species ',a,'.')

                stop
            end if
        else if (nse.ge.nern1 .and. nse.le.nern2) then
            ! Generic ion exchange species.
            if (jfl .eq. 0) then
                ! Mass balance.
                sx = 0.

                ! do nss = narn1,narn2
                !   ns = jcsort(nss)
                !   weight(ns) = coefst(csts,nsts,nstsmx,nstsr,nb,ns,nstmax)
                !   sx = sx + weight(ns)*conc(ns)
                ! enddo
                do nss = nern1,nern2
                    ns = jcsort(nss)
                    weight(ns) = coefst(csts,nsts,nstsmx,nstsr,nb,ns,nstmax)
                    sx = sx + weight(ns)*conc(ns)
                end do

                ctot = coval(nb)
                dx = sx - ctot
                alpha(krow) = dx
                beta(krow) = dx/ctot
            else
                ! Have found a bad jflag value.
                ! Calling sequence substitutions:
                !   uzvec1(krow) for unam48
                call fmspnm(jlen,uzvec1(krow),uspn56)
                write (noutpt,1000) jfl,uspn56(1:jlen)
                write (nttyo,1000) jfl,uspn56(1:jlen)
            end if
        end if
    end do

    go to 300

210 continue

    ! Compute residuals for EQ6.
    ! Compute mass balance elements.
    do krow = 1,kbt
        nb = iindx1(krow)
        nse = nbasp(nb)

        if ((nse.ne.no2gaq .and. nse.ne.nelect) .or. qredox) then
            ! Note: this calculation is not made if the current basis
            ! species is O2(g) or e-, and the current problem has no redox
            ! aspect. In the absence of a redox aspect, the log fO2 or the
            ! log a(e-) is treated as a "dead" iteration variable.
            sx = 0.

            do nss = 1,nst
                ns = jssort(nss)
                weight(ns) = coefst(csts,nsts,nstsmx,nstsr,nb,ns,nstmax)

                if (weight(ns) .ne. 0.) then
                    if (mosp(ns) .gt. 0.) then
                        tx = weight(ns)*mosp(ns)
                        sx = sx + tx
                    end if
                end if
            end do

            sx = sx - mtb(nb)
            alpha(krow) = sx
            amtb(nb) = mtb(nb)

            if (krow.eq.kwater .or. krow.eq.khydr .or.      krow.eq.ko2gaq .or. krow.eq.kelect) then
                ! For the species H2O, H+, OH-, O2(g,aq), and e-, compute
                ! "absolute" mass balance totals defined as:
                !   sum of positive terms + abs(sum of negative terms)
                ! Use this to calculate the corresponding beta residuals.
                sx = 0.

                do nss = 1,nst
                    ns = jssort(nss)
                    weight(ns) = coefst(csts,nsts,nstsmx,nstsr,nb,ns,nstmax)
                    tx = weight(ns)*mosp(ns)
                    atx = abs(tx)
                    sx = sx + atx
                end do

                amtb(nb) = sx
            end if

            beta(krow) = alpha(krow)/amtb(nb)
        end if
    end do

    ! Compute mass action elements.
    do krow = km1,kxt
        ns = iindx1(krow)
        call afcalc(actlg,af,afcnst,cdrs,jflag,jsflag,ndrs,ndrsmx,ndrsr,ns,nstmax,si,xlks)
        alpha(krow) = si
        beta(krow) = si
    end do

    ! Get characteristic residual parameters.
    !   betamx = max norm of beta
    !   ubetmx = name of the corresponding species
    !   bbig   = value of largest (positive) mass balance residual
    !   ubbig  = name of the corresponding species
    !   bneg   = value of largest (negative) mass balance residual
    !   ubneg  = name of the corresponding species
300 continue
    betamx = 0.
    ibetmx = 0
    ubetmx = unone

    do kcol = 1,kdim
        abeta = abs(beta(kcol))

        if (abeta .gt. betamx) then
            betamx = abeta
            ibetmx = kcol
            ubetmx = uzvec1(kcol)
        end if
    end do

    bbig = 0.
    bneg = 0.
    ubbig = unone
    ubneg = unone

    if (.not.q6mode) then
        if (kwater .gt. 0) then
            bwater = beta(kwater)
            beta(kwater) = 0.
        end if

        if (kelect .gt. 0) then
            belect = beta(kelect)
            beta(kelect) = 0.
        end if

        if (ko2gaq .gt. 0) then
            bo2gaq = beta(ko2gaq)
            beta(ko2gaq) = 0.
        end if
    end if

    do kcol = 1,kbt
        nb = iindx1(kcol)
        nse = nbasp(nb)

        if (jflag(nse) .le. 15) then
            bx = beta(kcol)

            if (bx .lt. bneg) then
                bneg = bx
                ubneg = uzvec1(kcol)
            else if (bx .gt. bbig) then
                bbig = bx
                ubbig = uzvec1(kcol)
            end if
        end if
    end do

    if (.not.q6mode) then
        if (kwater .gt. 0) then
            beta(kwater) = bwater
        end if

        if (kelect .gt. 0) then
            beta(kelect) = belect
        end if

        if (ko2gaq .gt. 0) then
            beta(ko2gaq) = bo2gaq
        end if
    end if
end subroutine betas