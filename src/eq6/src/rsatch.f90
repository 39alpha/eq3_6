subroutine rsatch(csts,egers,egexs,iern1,ietmax,iindx1,iktmax,iopt,ipndx1,jcode,jern1,jern2,jetmax,jgext,jpflag,jreac,kmax,km1,kmt,kx1,kxt,loph,losp,moph,morr,mosp,mrgers,mtb,mtb0,nbaspd,nbtmax,ncmpr,nern1,nern2,nert,nertmx,netmax,ngext,noptmx,noutpt,nptmax,nrct,nrctmx,nrk,nrndex,nstmax,nsts,nstsmx,nstsr,nttyo,nxridx,nxrt,nxrtmx,qreq,rxbar,tolxsf,uphase,ureac,uspec,xbar,xbarlg,zvclg1,zvec1)
    !! This subroutine tests reactants for saturation.
    !! This subroutine is called by:
    !!   EQ6/path.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: ietmax
    integer :: iktmax
    integer :: jetmax
    integer :: kmax
    integer :: nbtmax
    integer :: nertmx
    integer :: netmax
    integer :: noptmx
    integer :: nptmax
    integer :: nrctmx
    integer :: nstmax
    integer :: nstsmx
    integer :: nxrtmx

    integer :: noutpt
    integer :: nttyo

    integer :: iindx1(kmax)
    integer :: iopt(noptmx)
    integer :: ipndx1(kmax)
    integer :: jcode(nrctmx)
    integer :: jern1(jetmax,netmax)
    integer :: jern2(jetmax,netmax)
    integer :: jgext(netmax)
    integer :: jpflag(nptmax)
    integer :: jreac(nrctmx)
    integer :: nbaspd(nbtmax)
    integer :: ncmpr(2,nptmax)
    integer :: ngext(jetmax,netmax)
    integer :: nrk(2,nrctmx)
    integer :: nrndex(nrctmx)
    integer :: nsts(nstsmx)
    integer :: nstsr(2,nstmax)
    integer :: nxridx(nrctmx)

    integer :: iern1
    integer :: km1
    integer :: kmt
    integer :: kx1
    integer :: kxt
    integer :: nern1
    integer :: nern2
    integer :: nert
    integer :: nrct
    integer :: nxrt

    logical :: qreq

    character(len=48) :: uspec(nstmax)
    character(len=24) :: uphase(nptmax)
    character(len=24) :: ureac(nrctmx)

    real(kind=8) :: csts(nstsmx)
    real(kind=8) :: egexs(ietmax,jetmax,netmax)
    real(kind=8) :: loph(nptmax)
    real(kind=8) :: egers(ietmax,jetmax,netmax)
    real(kind=8) :: losp(nstmax)
    real(kind=8) :: moph(nptmax)
    real(kind=8) :: morr(nrctmx)
    real(kind=8) :: mosp(nstmax)
    real(kind=8) :: mrgers(ietmax,jetmax,nertmx)
    real(kind=8) :: mtb(nbtmax)
    real(kind=8) :: mtb0(nbtmax)
    real(kind=8) :: rxbar(iktmax,nxrtmx)
    real(kind=8) :: xbar(nstmax)
    real(kind=8) :: xbarlg(nstmax)
    real(kind=8) :: zvclg1(kmax)
    real(kind=8) :: zvec1(kmax)

    real(kind=8) :: tolxsf

    ! Local variable declarations.
    integer :: ie
    integer :: ik
    integer :: je
    integer :: jlen
    integer :: j2
    integer :: kcol
    integer :: kount
    integer :: n
    integer :: nb
    integer :: ne
    integer :: ner
    integer :: np
    integer :: npp
    integer :: nrc
    integer :: nrn1
    integer :: nrn2
    integer :: nr1
    integer :: nr2
    integer :: ns
    integer :: nss
    integer :: nt
    integer :: nxr

    integer :: ilnobl

    character(len=56) :: uspn56

    real(kind=8) :: dce
    real(kind=8) :: dm
    real(kind=8) :: dmx
    real(kind=8) :: dx
    real(kind=8) :: ecomp
    real(kind=8) :: ed
    real(kind=8) :: lx
    real(kind=8) :: lxx
    real(kind=8) :: mx
    real(kind=8) :: mxx
    real(kind=8) :: xcomp
    real(kind=8) :: xd

    real(kind=8) :: texp
    real(kind=8) :: tlg

    qreq = .false.

    ! Pure minerals.
    do nrc = 1,nrct
        if (jcode(nrc).eq.0 .and. nrk(2,nrc).eq.0) then
            if (jreac(nrc).eq.0 .or. jreac(nrc).eq.-1) then
                np = nrndex(nrc)

                if (jpflag(np) .gt. 0) then
                    go to 230
                end if

                jreac(nrc) = 0

                if (jpflag(np) .ne. -1) then
                    go to 230
                end if

                jreac(nrc) = -1

                if (iopt(1) .eq. 0) then
                    jreac(nrc) = 2
                    qreq = .true.
                    j2 = ilnobl(ureac(nrc))
                    write (noutpt,1000) ureac(nrc)(1:j2)
1000 format(/' Reactant ',a,' has become saturated and any',' remaining',/' unreacted mass has been transferred to',' the ES.',/)

                    dm = morr(nrc)
                    morr(nrc) = 0.
                    mx = moph(np) + dm
                    lx = tlg(mx)
                    moph(np) = mx
                    loph(np) = lx
                    ns = ncmpr(1,np)
                    mosp(ns) = mx
                    losp(ns) = lx

                    do kcol = km1,kmt
                        if (ns .eq. iindx1(kcol)) then
                            zvclg1(kcol) = losp(ns)
                            zvec1(kcol) = mx
                            nr1 = nstsr(1,ns)
                            nr2 = nstsr(2,ns)

                            do n = nr1,nr2
                                nb = nsts(n)
                                dce = csts(n)*dm
                                mtb(nb) = mtb(nb) + dce
                                mtb0(nb) = mtb0(nb) + dce
                            end do

                            go to 230
                        end if
                    end do

                    j2 = ilnobl(ureac(nrc))
                    write (noutpt,1010) ureac(nrc)(1:j2)
                    write (nttyo,1010) ureac(nrc)(1:j2)
1010 format(/' * Error - (EQ6/rsatch) The reactant ',a,'is saturated,',/7x,"but it isn't in the ES.")

                    stop
                end if
            end if
        end if

230 continue
    end do

    ! Solid solutions.
    if (nxrt .le. 0) then
        go to 340
    end if

    do nrc = 1,nrct
        if (jcode(nrc).eq.1 .and. nrk(2,nrc).eq.0) then
            if (jreac(nrc).eq.0 .or. jreac(nrc).eq.-1) then
                np = nrndex(nrc)

                if (jpflag(np) .ge. 2) then
                    go to 330
                end if

                jreac(nrc) = 0

                if (jpflag(np) .ne. -1) then
                    go to 330
                end if

                jreac(nrc) = -1

                if (iopt(1) .ne. 0) then
                    go to 330
                end if

                xcomp = 0.
                nrn1 = ncmpr(1,np)
                nrn2 = ncmpr(2,np)
                nt = nrn2 - nrn1 + 1
                nxr = nxridx(nrc)
                ik = 0

                do ns = nrn1,nrn2
                    ik = ik + 1
                    xd = xbar(ns) - rxbar(ik,nxr)
                    xcomp = xd*xd + xcomp
                end do

                xcomp = sqrt(xcomp/nt)

                if (xcomp .gt. tolxsf) then
                    j2 = ilnobl(ureac(nrc))
                    write (noutpt,1020) uphase(np)(1:j2)
                    write (nttyo,1020) uphase(np)(1:j2)
1020 format(/' * Note - (EQ6/rsatch) The solid solution',' reactant ',a,/7x,'has become saturated, but it',' differs in composition from the',/7x,'corresponding'," product phase. It can't be moved into the ES.")

                    go to 330
                end if

                jreac(nrc) = 2
                qreq = .true.
                write (noutpt,1000) ureac(nrc)
                dm = morr(nrc)
                morr(nrc) = 0.
                mx = moph(np) + dm
                moph(np) = mx
                loph(np) = tlg(mx)

                kount = 0

                do kcol = kx1,kxt
                    ns = iindx1(kcol)
                    npp = ipndx1(kcol)

                    if (npp .eq. np) then
                        kount = kount + 1
                        lxx = loph(np) + xbarlg(ns)
                        losp(ns) = lxx
                        mxx = texp(lxx)
                        mosp(ns) = mxx
                        zvclg1(kcol) = lxx
                        zvec1(kcol) = texp(lxx)
                        dmx = dm*xbar(ns)
                        nr1 = nstsr(1,ns)
                        nr2 = nstsr(2,ns)

                        do n = nr1,nr2
                            nb = nsts(n)
                            dce = csts(n)*dmx
                            mtb(nb) = mtb(nb) + dce
                            mtb0(nb) = mtb0(nb) + dce
                        end do
                    end if
                end do

                if (kount .le. 0) then
                    j2 = ilnobl(ureac(nrc))
                    write (noutpt,1010) ureac(nrc)(1:j2)
                    write (nttyo,1010) ureac(nrc)(1:j2)
                    stop
                end if
            end if
        end if

330 continue
    end do

340 continue

    ! Generic ion exchangers.
    if (nert .le. 0) then
        go to 440
    end if

    do nrc = 1,nrct
        if (jcode(nrc).eq.5 .and. nrk(2,nrc).eq.0) then
            if (jreac(nrc).eq.0 .or. jreac(nrc).eq.-1) then
                np = nrndex(nrc)

                if (jpflag(np) .ge. 2) then
                    go to 430
                end if

                jreac(nrc) = 0

                if (jpflag(np) .ne. -1) then
                    go to 430
                end if

                jreac(nrc) = -1

                ne = np - iern1 + 1
                ner = nxridx(nrc)

                ecomp = 0.

                do je = 1,jgext(ne)
                    do ie = 1,ngext(je,ne)
                        ed = egexs(ie,je,ne) - egers(ie,je,ner)
                        ecomp = ed*ed + ecomp
                    end do
                end do

                ecomp = sqrt(ecomp/nt)

                if (ecomp .gt. tolxsf) then
                    j2 = ilnobl(ureac(nrc))
                    write (noutpt,1040) uphase(np)(1:j2)
                    write (nttyo,1040) uphase(np)(1:j2)
1040 format(/' * Note - (EQ6/rsatch) The generic ion',' exchanger reactant'/7x,a,' has become saturated, but',' it differs in composition',/7x,'from the corresponding'," product phase. It can't be moved into the ES.")

                    go to 440
                end if

                jreac(nrc) = 2
                qreq = .true.
                write (noutpt,1000) ureac(nrc)
                dm = morr(nrc)
                morr(nrc) = 0.
                mx = moph(np) + dm
                moph(np) = mx
                loph(np) = tlg(mx)

                kount = 0

                ! Create matrix index range markers for exchange species.
                do kcol = 2,km1 - 1
                    npp = ipndx1(kcol)

                    if (npp .eq. np) then
                        kount = kount + 1
                        ns = iindx1(kcol)

                        do je = 1,jgext(ne)
                            if (ns.ge.jern1(je,ne) .and. ns.le.jern2(je,ne)) then
                                go to 410
                            end if
                        end do

                        ! Calling sequence substitutions:
                        !   uspec(ns) for unam48
                        call fmspnm(jlen,uspec(ns),uspn56)
                        write (noutpt,1050) uspn56(1:jlen)
                        write (nttyo,1050) uspn56(1:jlen)
1050 format(/' * Error - (EQ6/rsatch) Programming error'," trap: Can't determine",/7x,'the exchange site index',' je for the species',/7x,a,'.')

                        stop

410 continue
                        do ie = 1,ngext(je,ne)
                            nss = ie + jern1(je,ne) -1

                            if (ns .eq. nss) then
                                go to 420
                            end if
                        end do

                        ! Calling sequence substitutions:
                        !   uspec(ns) for unam48
                        call fmspnm(jlen,uspec(ns),uspn56)
                        write (noutpt,1060) uspn56(1:jlen)
                        write (nttyo,1060) uspn56(1:jlen)
1060 format(/' * Error - (EQ6/rsatch) Programming error'," trap: Can't determine",/7x,'the exchange species',' index ie of the species',/7x,a,'.')

                        stop

420 continue
                        dx = dm*mrgers(ie,je,ner)
                        mxx = mosp(ns) + dx
                        mosp(ns) = mxx
                        lxx = tlg(mxx)
                        losp(ns) = lxx
                        zvclg1(kcol) = lxx
                        zvec1(kcol) = mxx
                        nr1 = nstsr(1,ns)
                        nr2 = nstsr(2,ns)

                        if (je .le. 1) then
                            ! Have the first site. Increment the mass balance
                            ! totals in a straightforward manner.
                            do n = nr1,nr2
                                nb = nsts(n)
                                dce = dx*csts(n)
                                mtb(nb) = mtb(nb) + dce
                                mtb0(nb) = mtb0(nb) + dce
                            end do
                        else
                            ! Have a site beyond the first. Increment the mass
                            ! balance totals in the usual manner, except do not
                            ! increment here the total for the exchanger
                            ! substrate. The complete increment for the substrate
                            ! is obtained by considering only one site.
                            do n = nr1,nr2
                                nb = nsts(n)
                                ns = nbaspd(nb)

                                if (ns.lt.nern1 .or. ns.gt.nern2) then
                                    dce = dx*csts(n)
                                    mtb(nb) = mtb(nb) + dce
                                    mtb0(nb) = mtb0(nb) + dce
                                end if
                            end do
                        end if
                    end if
                end do

                if (kount .le. 0) then
                    j2 = ilnobl(ureac(nrc))
                    write (noutpt,1010) ureac(nrc)(1:j2)
                    write (nttyo,1010) ureac(nrc)(1:j2)
                    stop
                end if
            end if
        end if

430 continue
    end do

440 continue
end subroutine rsatch
