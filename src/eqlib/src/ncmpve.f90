subroutine ncmpve(acflg,act,actlg,cdrs,cgexj,eps100,iern1,iern2,ietmax,jern1,jern2,jetmax,jflag,jgext,jsflag,losp,mgext,moph,mosp,nbasp,nbt,nbtmax,ndrs,ndrsmx,ndrsr,netmax,noutpt,nptmax,nstmax,nttyo,ugexj,uphase,uspec,xbar,xbarlg,xlks)
    !! This subroutine computes part of the "expansion" of the
    !! basis set variable data that yields the "total" system
    !! description. This part is the expansion giving the mole
    !! fractions and number of moles of components of generic ion
    !! exchangers, for all exchange models. This part of the expansion
    !! is necessarily iterative in the general case. The overall
    !! expansion is conducted by EQLIB/ncmpex.f.
    !! This subroutine is called by:
    !!   EQLIB/ncmpex.f
    !! Principal input:
    !!   acflg  = array of logarithms of activity coefficients
    !!   actlg  = array of logarithms of species activities
    !!              (input consists of valid results for aqueous
    !!              basis species; output consists of the same
    !!              for generic ion exchanger species)
    !!   cdrs   = array of reaction coefficients
    !!   mosp   = array of numbers of moles of species
    !!              (input contains valid results for basis species
    !!              and carry-over values for non-basis species)
    !!   moph   = array of numbers of moles of phases
    !!   ndrs   = array parallel to cdrs giving the index of the
    !!              corresponding species
    !!   ndrsr  = array giving the range in the cdrs/ndrs arrays
    !!              corresonding to the reaction associated with a
    !!              given species
    !!   xlks   = array of equilibrium constants
    !! Principal output:
    !!   act    = array of species activities
    !!              (output consists of the subset for generic exchanger
    !!              species)
    !!   actlg  = array of logarithms of species activities
    !!   mgext  = array of the sum of the number of moles of exchanger
    !!              species in given sites of generic ion exchange
    !!              phases
    !!   mosp   = array of numbers of moles of species
    !!              (output consists of valid results for non-basis
    !!              species; valid results for basis species were
    !!              input and are retained)
    !!   losp   = array of logarithms of numbers of moles of species
    !!   xbar   = array of mole fractions of species
    !!   xbarlg = array of logarithms of mole fractions of species
    implicit none

    ! Calling sequence variable declarations.
    integer :: ietmax
    integer :: jetmax
    integer :: nbtmax
    integer :: ndrsmx
    integer :: netmax
    integer :: nptmax
    integer :: nstmax

    integer :: noutpt
    integer :: nttyo

    integer :: jern1(jetmax,netmax)
    integer :: jern2(jetmax,netmax)
    integer :: jflag(nstmax)
    integer :: jgext(netmax)
    integer :: jsflag(nstmax)
    integer :: nbasp(nbtmax)
    integer :: ndrs(ndrsmx)
    integer :: ndrsr(2,nstmax)

    integer :: iern1
    integer :: iern2
    integer :: nbt

    character(len=48) :: uspec(nstmax)
    character(len=24) :: uphase(nptmax)
    character(len=8) :: ugexj(jetmax,netmax)

    real(kind=8) :: acflg(nstmax)
    real(kind=8) :: act(nstmax)
    real(kind=8) :: actlg(nstmax)
    real(kind=8) :: cdrs(ndrsmx)
    real(kind=8) :: cgexj(jetmax,netmax)
    real(kind=8) :: losp(nstmax)
    real(kind=8) :: mgext(jetmax,netmax)
    real(kind=8) :: moph(nptmax)
    real(kind=8) :: mosp(nstmax)
    real(kind=8) :: xbar(nstmax)
    real(kind=8) :: xbarlg(nstmax)
    real(kind=8) :: xlks(nstmax)

    real(kind=8) :: eps100

    ! Local variable declarations.
    integer :: icycle
    integer :: icycmx
    integer :: iterve
    integer :: itvemx
    integer :: je
    integer :: j2
    integer :: j3
    integer :: j4
    integer :: n
    integer :: nb
    integer :: ndlxfc
    integer :: ne
    integer :: np
    integer :: nrr1
    integer :: nrr2
    integer :: nrsxfc
    integer :: nr1
    integer :: nr2
    integer :: ns
    integer :: nse
    integer :: nsi
    integer :: nss

    integer :: ilnobl
    integer :: nbasis

    logical :: qprint

    character(len=8) :: ux8

    real(kind=8) :: aajx
    real(kind=8) :: ardlx
    real(kind=8) :: ardlxo
    real(kind=8) :: aresx
    real(kind=8) :: aresxo
    real(kind=8) :: cgxj
    real(kind=8) :: delx
    real(kind=8) :: delxl
    real(kind=8) :: delxu
    real(kind=8) :: dlxfnc
    real(kind=8) :: dm
    real(kind=8) :: cxs
    real(kind=8) :: lax
    real(kind=8) :: lxx
    real(kind=8) :: mtxj
    real(kind=8) :: mtxjc
    real(kind=8) :: resx
    real(kind=8) :: rhsx
    real(kind=8) :: rsxfnc
    real(kind=8) :: resxo
    real(kind=8) :: xbarmx
    real(kind=8) :: xbarsm
    real(kind=8) :: xx

    real(kind=8) :: texp
    real(kind=8) :: tlg

    ! Data statements.
    data icycmx /10/
    data itvemx /100/

    qprint = .false.

    ! Loop over all generic ion exchanger phases.
    do np = iern1,iern2
        ne = np - iern1 + 1

        if (moph(np) .gt. 0.) then
            ! Have found an exchanger phase present in the system.
            ! Loop over all sites.
            do je = 1,jgext(ne)
                nrr1 = jern1(je,ne)
                nrr2 = jern2(je,ne)
                cgxj = cgexj(je,ne)

                ! Determine which species is the basis species for the
                ! current site of the current exchanger phase.
                nse = 0.

                do ns = nrr1,nrr2
                    nb = nbasis(nbasp,nbt,nbtmax,ns)

                    if (nb.gt.0 .and. jsflag(ns).le.0) then
                        nse = ns
                    end if
                end do

                ! Here mtxj (saved as mgext(je,ne)) is the sum of the
                ! numbers of moles of the exchanger species for the current
                ! site. This is equivalent to the sum of the numbers of
                ! moles of the corresponding exchange ions. Depending on
                ! the ion exchange model in question, this sum may not be
                ! simply proportional to the number of moles of the exchanger
                ! phase. More generally, it may be a complex function of
                ! the site composition. In this case, the value of mtxj
                ! must be calculated using an iterative process.
                ! Calculate an initial estimate from mostly inherited
                ! values (use number of moles variables for all exchanger
                ! species for the current site). Only the value for the
                ! species in the basis set is actually current.
                mtxj = 0.

                do ns = nrr1,nrr2
                    mtxj = mtxj + mosp(ns)
                end do

                ! Optional output.
                if (qprint) then
                    j2 = ilnobl(ugexj(je,ne))
                    j3 = ilnobl(uphase(np))
                    write (noutpt,1040) ugexj(je,ne)(1:j2),uphase(np)(1:j3)
1040 format(/' Pre-Newton-Raphson expansion for site ',a,/' of phase ',a,':',/)
                end if

                icycle = -1
                resx = 0.

                ! Loop back here for another cycle of pre-Newton-Raphson
                ! correction.
110 continue
                icycle = icycle + 1

                ! Make a tentative expansion to estimate the mole fractions,
                ! activities, and numbers of moles of all species belonging
                ! to the current site.
                call ncmpvh(acflg,act,actlg,cdrs,cgxj,jflag,jsflag,losp,mosp,mtxj,nbasp,nbt,nbtmax,ndrs,ndrsmx,ndrsr,nrr1,nrr2,nstmax,xbar,xbarlg,xlks)

                ! Check the sum of the mole fractions.
                xbarsm = 0.

                do ns = nrr1,nrr2
                    xbarsm = xbarsm + xbar(ns)
                end do

                ! Calculate the residual function.
                resxo = resx
                resx = xbarsm - 1.0

                ! Calculate the improvement functions.
                aresxo = abs(resxo)
                aresx = abs(resx)

                if (aresxo .gt. 0.) then
                    rsxfnc = (aresxo - aresx)/aresxo
                else
                    rsxfnc = 0.
                end if

                ! Optional output.
                if (qprint) then
                    if (icycle .ge. 0) then
                        write (noutpt,1080) icycle,mtxj,resx,rsxfnc
1080 format(2x,'iter= ',i2,', mtxj= ',1pe10.3,', resx= ',1e10.3,', rsxfnc= ',e10.3)
                    end if
                end if

                if (resx.gt.0.2 .or. resx.lt.-0.2) then
                    if (icycle .lt. icycmx) then
                        ! Make a new estimate of mtxj to use in a new
                        ! expansion calculation.
                        ! Normalize the mole fractions and recompute the
                        ! activities.
                        do ns = nrr1,nrr2
                            xx = xbar(ns)/xbarsm
                            xbar(ns) = xx
                            xbarlg(ns) = tlg(xx)
                            actlg(ns) = cgxj*(xbarlg(ns) + acflg(ns))
                        end do

                        ! Find the species with the largest mole fraction.
                        xbarmx = 0.
                        nsi = 0

                        do ns = nrr1,nrr2
                            if (xbar(ns) .gt. xbarmx) then
                                xbarmx = xbar(ns)
                                nsi = ns
                            end if
                        end do

                        ! Assume that mole fraction is correct. If the species
                        ! with the largest mole fraction is the basis species,
                        ! use its mole fraction to re-estimate mtxj. If it is
                        ! not the basis species, use the mass action relation
                        ! to calculate the mole fraction of the basis species,
                        ! and use that in turn to re-estimate mtxj.
                        if (nsi .eq. nse) then
                            xx = xbar(nse)
                        else
                            nr1 = ndrsr(1,nsi)
                            nr2 = ndrsr(2,nsi)
                            lax = xlks(nsi)

                            do n = nr1,nr2
                                nss = ndrs(n)

                                if (nss .eq. nse) then
                                    cxs = cdrs(n)
                                else
                                    lax = lax - cdrs(n)*actlg(nss)
                                end if
                            end do

                            lax = lax/cxs
                            lxx = (lax/cgxj) - acflg(nse)
                            xx = texp(lxx)
                        end if

                        mtxjc = mosp(nse)/xx

                        if (mtxjc .lt. mosp(nse)) then
                            ! The value of mtxjc is less than the theoretical
                            ! minimum of mosp(nse) (the calculated mole fraction
                            ! xx is greater than 1.0). Move mtxj halfway to its
                            ! theoretical lower limit.
                            dm = mosp(nse) - mtxj
                            mtxj = mtxj + 0.5*dm
                        else
                            ! The value of mtxjc is okay. Use it.
                            mtxj = mtxjc
                        end if

                        go to 110
                    end if
                end if

                ! Optional output.
                if (qprint) then
                    write (noutpt,1090)
1090 format(1x)
                end if

                ! Begin Newton-Raphson iteration. This presumes that a
                ! starting value of mtxj is already provided. The initial
                ! residual is calculated, even if it is already known.
                ! This is intended to make the N-R block as independent
                ! as possible from the above pre-N-R block.
                ! Optional output.
                if (qprint) then
                    j2 = ilnobl(ugexj(je,ne))
                    j3 = ilnobl(uphase(np))
                    write (noutpt,1140) ugexj(je,ne)(1:j2),uphase(np)(1:j3)
1140 format(/' Newton-Raphson expansion for site ',a,/' of phase ',a,':',/)
                end if

                iterve = 0
                resx = 0.
                delx = 0.
                resxo = 0.
                rsxfnc = 0.
                ardlx = 0.
                ardlxo = 0.
                dlxfnc = 0.
                nrsxfc = 0
                ndlxfc = 0

                ! Loop back here to continue more Newton-Raphson
                ! iterations.
120 continue

                ! Make a tentative expansion to estimate the mole fractions,
                ! activities, and numbers of moles of all species belonging
                ! to the current site.
                call ncmpvh(acflg,act,actlg,cdrs,cgxj,jflag,jsflag,losp,mosp,mtxj,nbasp,nbt,nbtmax,ndrs,ndrsmx,ndrsr,nrr1,nrr2,nstmax,xbar,xbarlg,xlks)

                ! Estimate the mole fractions and activities of the
                ! basis species for the current site.
                ! Check the sum of the mole fractions.
                xbarsm = 0.

                do ns = nrr1,nrr2
                    xbarsm = xbarsm + xbar(ns)
                end do

                ! Calculate the residual function.
                resxo = resx
                resx = xbarsm - 1.0

                ! Calculate the improvement functions.
                aresxo = abs(resxo)
                aresx = abs(resx)

                if (aresxo .gt. 0.) then
                    rsxfnc = (aresxo - aresx)/aresxo
                else
                    rsxfnc = 0.
                end if

                ardlxo = ardlx
                ardlx = abs(delx/mtxj)

                if (ardlxo .gt. 0.) then
                    dlxfnc = (ardlxo - ardlx)/ardlxo
                else
                    dlxfnc = 0.
                end if

                if (rsxfnc .ge. 0.) then
                    nrsxfc = 0
                else
                    nrsxfc = nrsxfc + 1
                end if

                if (dlxfnc .ge. 0.) then
                    ndlxfc = 0
                else
                    ndlxfc = ndlxfc + 1
                end if

                ! Optional output.
                if (qprint) then
                    if (iterve .ge. 0) then
                        write (noutpt,1150) iterve,mtxj,resx,rsxfnc
1150 format(2x,'iter= ',i2,', mtxj= ',1pe10.3,', resx= ',1e10.3,', rsxfnc= ',e10.3)
                    end if
                end if

                ! Test for convergence.
                if (aresx .le. eps100) then
                    go to 170
                end if

                ! Test for the maximum number of iterations.
                if (iterve .eq. itvemx) then
                    j2 = ilnobl(ugexj(je,ne))
                    j3 = ilnobl(uphase(np))
                    write (ux8,'(i5)') itvemx
                    call lejust(ux8)
                    j4 = ilnobl(ux8)
                    write (noutpt,1170) ugexj(je,ne)(1:j2),uphase(np)(1:j3),ux8(1:j4)
                    write (nttyo,1170) ugexj(je,ne)(1:j2),uphase(np)(1:j3),ux8(1:j4)
1170 format(/' * Warning - (EQLIB/ncmpve) Newton-Raphson',' iteration failed to'/7x,'converge during an attempt',' to expand the system description',/7x,'for site ',a,' of phase ',a,'. It reached',/7x,'the limit of ',a,' iterations.')

                    write (nttyo,1172) aresx
                    write (noutpt,1172) aresx
1172 format(/9x,'abs. residual= ',1pe12.5)

                    go to 170
                end if

                ! Test for non-convergent behavior.
                if ((nrsxfc.ge.10 .and. ndlxfc.ge.10) .or.        nrsxfc.ge.20 .or. ndlxfc.ge.20) then
                    j2 = ilnobl(ugexj(je,ne))
                    j3 = ilnobl(uphase(np))
                    write (noutpt,1180) ugexj(je,ne)(1:j2),uphase(np)(1:j3)
                    write (nttyo,1180) ugexj(je,ne)(1:j2),uphase(np)(1:j3)
1180 format(/' * Warning - (EQLIB/ncmpve) Newton-Raphson',' iteration failed to'/7x,'converge during an attempt',' to expand the system description',/7x,'for site ',a,' of phase ',a,'. Non-convergent',/7x,'behavior was',' detected.')

                    go to 170
                end if

                ! Do another iteration.
                iterve = iterve + 1

                ! Set the right-hand-side.
                rhsx = -resx

                ! Calculate the derivative (1 x 1 Jacobian) of the
                ! residual function.
                aajx = 0.

                do ns = nrr1,nrr2
                    if (ns .eq. nse) then
                        aajx = aajx + xbar(ns)
                    else
                        nr1 = ndrsr(1,ns)
                        nr2 = ndrsr(2,ns)

                        do n = nr1,nr2
                            nss = ndrs(n)

                            if (nss .eq. nse) then
                                aajx = aajx - cdrs(n)*xbar(ns)/cdrs(nr1)
                                go to 150
                            end if
                        end do

150 continue
                    end if
                end do

                aajx = -aajx/mtxj

                if (aajx .eq. 0.) then
                    j2 = ilnobl(ugexj(je,ne))
                    j3 = ilnobl(uphase(np))
                    write (noutpt,1190) ugexj(je,ne)(1:j2),uphase(np)(1:j3)
                    write (nttyo,1190) ugexj(je,ne)(1:j2),uphase(np)(1:j3)
1190 format(/' * Error - (EQLIB/ncmpve) Have encountered',' a zero 1 x 1 Jacobian'/7x,'during Newton-Raphson',' iteration to expand the description of',/7x,'site ',a,' of phase ',a,'.')

                    stop
                end if

                ! Solve for the Newton-Raphson correction term.
                delx = rhsx/aajx

                ! Apply under-relaxation limits.
                delxu = 9999.*mtxj
                delxl = -0.9999*mtxj
                delx = min(delx,delxu)
                delx = max(delx,delxl)

                ! Apply a special under-relaxation limit.
                mtxjc = mtxj + delx

                if (mtxjc .lt. mosp(nse)) then
                    ! The value of mtxjc is less than the theoretical
                    ! minimum of mosp(nse).  Set up to move mtxj halfway
                    ! to its theoretical lower limit.
                    delx = 0.5*(mosp(nse) - mtxj)
                end if

                ! Apply the correction term.
                mtxj = mtxj + delx
                go to 120

                ! Save the sum of the number of moles of exchanger species
                ! for the current site of the current exchanger phase.
170 continue
                mgext(je,ne) = mtxj

                ! Optional output.
                if (qprint) then
                    write (noutpt,1200)
1200 format(1x)
                end if
            end do
        end if
    end do
end subroutine ncmpve