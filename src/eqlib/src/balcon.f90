subroutine balcon(aamatr,aamgex,al10,cdrs,cjbasp,cnufac,conc,dmlge,eps100,ggmgex,iern1,ietmax,iimgex,iindx1,ixbasp,jcsort,jern1,jern2,jetmax,jjsort,jsitex,kbt,kmax,krow,narn1,narn2,nbasp,nbtmax,ndrs,ndrsmx,ndrsr,nern1,nern2,netmax,noutpt,nphasx,nstmax,nttyo,rhsgex,uspec,weight,xbar)
    !! This subroutine computes a row of the EQ3NR Jacobian matrix
    !! for one of the following:
    !!   Mass balance constraint
    !!   Charge balance (charge adjustment) constraint
    !!   Alkalinity balance constraint
    !! This subroutine is not used to compute any part of the EQ6
    !! Jacobian matrix.
    !! This subroutine is called by:
    !!   EQLIB/matrix.f
    !! Input:
    !!   conc   = array of species concentrations
    !! Output:
    !!   aamatr = the Jacobian matrix
    implicit none

    ! Calling sequence variable declarations.
    integer :: noutpt
    integer :: nttyo

    integer :: ietmax
    integer :: jetmax
    integer :: kmax
    integer :: nbtmax
    integer :: ndrsmx
    integer :: netmax
    integer :: nstmax

    integer :: iimgex(ietmax)
    integer :: iindx1(kmax)
    integer :: ipvgex(ietmax)
    integer :: ixbasp(nbtmax)
    integer :: jcsort(nstmax)
    integer :: jern1(jetmax,netmax)
    integer :: jern2(jetmax,netmax)
    integer :: jjsort(nstmax)
    integer :: jsitex(nstmax)
    integer :: nbasp(nbtmax)
    integer :: ndrs(ndrsmx)
    integer :: ndrsr(2,nstmax)
    integer :: nphasx(nstmax)

    integer :: iern1
    integer :: kbt
    integer :: krow
    integer :: narn1
    integer :: narn2
    integer :: nern1
    integer :: nern2

    character(len=48) :: uspec(nstmax)

    real(kind=8) :: aamatr(kmax,kmax)
    real(kind=8) :: aamgex(ietmax,ietmax)
    real(kind=8) :: cdrs(ndrsmx)
    real(kind=8) :: cjbasp(nbtmax)
    real(kind=8) :: cnufac(nstmax)
    real(kind=8) :: conc(nstmax)
    real(kind=8) :: dmlge(ietmax)
    real(kind=8) :: ggmgex(ietmax,ietmax)
    real(kind=8) :: rhsgex(ietmax)
    real(kind=8) :: weight(nstmax)
    real(kind=8) :: xbar(nstmax)

    real(kind=8) :: al10
    real(kind=8) :: eps100

    ! Local variable declarations.
    integer :: icol
    integer :: idim
    integer :: ier
    integer :: irow
    integer :: je
    integer :: jlen1
    integer :: jlen2
    integer :: kcol
    integer :: nb
    integer :: nbb
    integer :: ne
    integer :: np
    integer :: nrr1
    integer :: nrr2
    integer :: nr1
    integer :: nsc
    integer :: nsi
    integer :: nsj
    integer :: nss

    logical :: qej
    logical :: qmj
    logical :: qpr
    logical :: qrhsgz
    logical :: qwj

    character(len=56) :: usp156
    character(len=56) :: usp256

    real(kind=8) :: axfc
    real(kind=8) :: cgxj
    real(kind=8) :: cxcc
    real(kind=8) :: cxjc
    real(kind=8) :: dmlgej
    real(kind=8) :: sumx

    real(kind=8) :: coefdr

    data qpr    /.false./

    nb = iindx1(krow)
    nsi = nbasp(nb)

    if (nsi.gt.narn1 .and. nsi.le.narn2) then
        ! Have an aqueous species mass balance.
        ! Build columns 1 through kbt.
        do kcol = 1,kbt
            nbb = iindx1(kcol)
            nsj = nbasp(nbb)
            sumx = 0.

            ! First loop over aqueous species.
            do nss = narn1,narn2
                nsc = jcsort(nss)

                if (weight(nsc) .ne. 0.) then
                    if (cnufac(nsc) .ne. 0.) then
                        ! Calling sequence substitutions:
                        !   nsj for nse
                        !   nsc for ns
                        cxjc = coefdr(cdrs,ndrs,ndrsmx,ndrsr,nsj,nsc,nstmax)
                        sumx = sumx + weight(nsc)*cxjc*cnufac(nsc)
                    end if
                end if
            end do

            ! Do not loop over generic ion exchanger species here.
            ! They do not contribute to non-exchanger mass balances
            ! in EQ3NR. This is by design. The situation is different
            ! in EQ6.
            ! Complete the calculation of the current matrix element.
            if (nsj .eq. narn1) then
                aamatr(krow,kcol) = al10*sumx
            else
                aamatr(krow,kcol) = al10*(weight(nsj)*conc(nsj) + sumx)
            end if
        end do
    else if (nsi.ge.nern1 .and. nsi.le.nern2) then
        ! Have a generic ion exchanger species mass balance corresponding
        ! to a given exchanger phase and site.
        cgxj = cjbasp(nb)

        je = jsitex(nsi)
        np = nphasx(nsi)
        ne = np - iern1 + 1

        ! Build columns 1 through kbt.
        do kcol = 1,kbt
            nbb = iindx1(kcol)
            nsj = nbasp(nbb)
            qej = nsj.ge.nern1 .and. nsj.le.nern2
            qwj = nsj .eq. narn1
            qmj = nsj.gt.narn1 .and. nsj.le.narn2

            ! Calculate the dmlge vector for the current site/exchanger
            ! phase and the current basis species. In essence,
            !   dmlge(i,j) is partial m(i)/partial log m(j),
            ! where i denotes a non-basis species for the current site/
            ! exchanger phase and j denotes the basis species for the
            ! current column index kcol. Theoretically, dmlge(j,j) could
            ! also be included, but it can be obtained from a
            ! straightforward calculation (as dmlgej below), in contrast
            ! to the others, which must be obtained simultaneously.
            ! Do not loop over aqueous species as non-basis species.
            ! They can not contribute to this kind of mass balance,
            ! because they do not contain the exchanger ligand.
            ! Loop over all relevant generic ion exchanger species.
            ! Generate an entry only for an active non-basis species.
            ! In the present context of EQ3NR, this includes only
            ! those exchanger species belonging to the site and phase
            ! associated with the current row. The do loop immediately
            ! below is over all species in the current exchanger phase
            ! (e.g., runs over all species on all sites). That is done
            ! to permit the use of the jcsort array. The sorting in that
            ! array is over the entire phase, not by individual site.
            ! The weight array filters out contributions from species
            ! not in the current site in all calculations except those
            ! involving the idim variable and the iimgex array. There,
            ! it is necessary to check the site to which a species
            ! belongs.
            idim = 0
            nrr1 = jern1(je,ne)
            nrr2 = jern2(je,ne)

            do nss = nrr1,nrr2
                nsc = jjsort(nss)

                if (weight(nsc) .ne. 0.) then
                    if (cnufac(nsc) .ne. 0.) then
                        ! Calling sequence substitutions:
                        !   nsj for nse
                        !   nsc for ns
                        cxjc = coefdr(cdrs,ndrs,ndrsmx,ndrsr,nsj,nsc,nstmax)

                        if (cxjc .ne. 0.) then
                            idim = idim + 1
                            iimgex(idim) = nsc
                        end if
                    end if
                end if
            end do

            do irow = 1,idim
                nsc = iimgex(irow)
                nr1 = ndrsr(1,nsc)
                cxcc = cdrs(nr1)

                ! Calling sequence substitutions:
                !   nsj for nse
                !   nsc for ns
                cxjc = coefdr(cdrs,ndrs,ndrsmx,ndrsr,nsj,nsc,nstmax)

                if (qej) then
                    ! The current column is for a generic ion exchanger
                    ! species.
                    axfc = xbar(nsc)*(1.0 + (cxjc/cxcc))

                    do icol = 1,idim
                        aamgex(irow,icol) = -axfc
                    end do

                    aamgex(irow,irow) = 1.0 + aamgex(irow,irow)
                    rhsgex(irow) = al10*(-(cxjc/cxcc)*conc(nsc)        + axfc*conc(nsj))
                else if (qwj) then
                    ! The current column is for solvent water.
                    axfc = xbar(nsc)

                    do icol = 1,idim
                        aamgex(irow,icol) = -axfc
                    end do

                    aamgex(irow,irow) = 1.0 + aamgex(irow,irow)
                    rhsgex(irow) = -al10*(cxjc/cxcc)*conc(nsc)
                else if (qmj) then
                    ! The current column is for an aqueous solute species.
                    axfc = xbar(nsc)*(1.0 + (cxjc/cxcc))

                    do icol = 1,idim
                        aamgex(irow,icol) = -axfc
                    end do

                    aamgex(irow,irow) = 1.0 + aamgex(irow,irow)
                    rhsgex(irow) = -al10*(cxjc/(cxcc*cgxj))*conc(nsc)
                else
                    ! Calling sequence substitutions:
                    !   jlen1 for jlen
                    !   uspec(nsi) for unam48
                    !   usp156 for uspn56
                    call fmspnx(jlen1,uspec(nsi),usp156)

                    ! Calling sequence substitutions:
                    !   jlen2 for jlen
                    !   uspec(nsj) for unam48
                    !   usp256 for uspn56
                    call fmspnx(jlen2,uspec(nsj),usp256)

                    write (noutpt,1110) usp156(1:jlen1),usp256(1:jlen2)
                    write (nttyo,1110) usp156(1:jlen1),usp256(1:jlen2)
1110 format(/' * Error - (EQLIB/balcon) Programming',' error trap: In calculating',/7x,'the Jacobian',' matrix entry for the ',a,' row',/7x,' and the ',a,' column, the type of the column',/7x,' species',' does not correspond to a programmed possibility.')

                    stop
                end if
            end do

            ! Check for a zero rhsgex vector. Note that if idim is zero,
            ! qrhsgz will remain set to .true.
            qrhsgz = .true.

            do irow = 1,idim
                if (rhsgex(irow) .ne. 0.) then
                    qrhsgz = .false.
                    go to 300
                end if
            end do

300 continue

            if (qrhsgz) then
                ! Have a zero right-hand-side vector. Set all dmlge elements
                ! to zero.
                do irow = 1,idim
                    dmlge(irow) = 0.
                end do
            else
                ! Solve the matrix equation. If EQLIBU/msolvr.f can't solve
                ! the matrix, it is because the matrix is either zero
                ! (ier = 1) or non-zero, but computationally singular
                ! (ier = 2).
                ! Calling sequence substitutions:
                !   aamgex for aamatr
                !   dmlge for delvec
                !   ggmgex for ggmatr
                !   ipvgex for ipivot
                !   idim for kdim
                !   ietmax for kmax
                !   rhsgex for rhsvec
                call msolvr(aamgex,dmlge,ggmgex,ier,ipvgex,idim,ietmax,noutpt,nttyo,qpr,rhsgex)

                if (ier .ne. 0) then
                    ! The matrix is zero or it is computationally singular.
                    ! Calling sequence substitutions:
                    !   jlen1 for jlen
                    !   uspec(nsi) for unam48
                    !   usp156 for uspn56
                    call fmspnx(jlen1,uspec(nsi),usp156)

                    ! Calling sequence substitutions:
                    !   jlen2 for jlen
                    !   uspec(nsj) for unam48
                    !   usp256 for uspn56
                    call fmspnx(jlen2,uspec(nsj),usp256)

                    write (noutpt,1120) usp156(1:jlen1),usp256(1:jlen2)
                    write (nttyo,1120) usp156(1:jlen1),usp256(1:jlen2)
1120 format(/' * Error - (EQLIB/balcon) Have encountered a',' zero or computationally',/7x,'singular matrix while',' trying to calculate the Jacobian matrix entry',/7x,'for the ',a,' row and the ',a,' column.')

                    stop
                end if
            end if

            if (qrhsgz) then
                sumx = 0.
            else
                sumx = 0.

                do irow = 1,idim
                    nsc = iimgex(irow)
                    sumx = sumx + weight(nsc)*dmlge(irow)
                end do
            end if

            if (nsj .ne. narn1) then
                dmlgej = al10*conc(nsj)
                sumx = sumx + weight(nsj)*dmlgej
            end if

            aamatr(krow,kcol) = sumx
        end do
    end if
end subroutine balcon