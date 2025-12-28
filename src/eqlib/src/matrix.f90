subroutine matrix(aamatr,al10,bpx,cdrs,cdrtw,cdrw,cjbasp,cnufac,conc,csts,dlogxw,eps100,ibpxmx,iebal,iern1,ietmax,iindx1,ipndx1,irdxc3,ixbasp,ixrn1,jcsort,jern1,jern2,jetmax,jflag,jjsort,jsitex,kbt,kction,kdim,kelect,khydr,kmax,kmt,km1,ko2gaq,kwater,kxt,kx1,mosp,narn1,narn2,nbasp,nbt,nbtmax,nbw,ncosp,ndrs,ndrsmx,ndrsr,nelect,nern1,nern2,netmax,noutpt,no2gaq,nphasx,nredox,nst,nstmax,nsts,nstsmx,nstsr,ntfx,ntfxmx,ntfxt,nttyo,nxtmax,omega,qredox,q6mode,tfx,ugexmo,uspec,weight,xbar,xbarw,xbarwc,zchar)
    !! This subroutine computes the Jacobian matrix J[z] (aamatr).
    !! This matrix is used by the algebraic equation solver.
    !! This subroutine is called by:
    !!   EQLIB/nrstep.f
    !!   EQ6/optmzr.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: ibpxmx
    integer :: ietmax
    integer :: jetmax
    integer :: kmax
    integer :: nbtmax
    integer :: ndrsmx
    integer :: netmax
    integer :: nstmax
    integer :: nstsmx
    integer :: ntfxmx
    integer :: nxtmax

    integer :: noutpt
    integer :: nttyo

    integer :: iindx1(kmax)
    integer :: ipndx1(kmax)
    integer :: ixbasp(nbtmax)
    integer :: jcsort(nstmax)
    integer :: jern1(jetmax,netmax)
    integer :: jern2(jetmax,netmax)
    integer :: jflag(nstmax)
    integer :: jjsort(nstmax)
    integer :: jsitex(nstmax)
    integer :: kction(nbtmax)
    integer :: nbasp(nbtmax)
    integer :: ncosp(nbtmax)
    integer :: ndrs(ndrsmx)
    integer :: ndrsr(2,nstmax)
    integer :: nphasx(nstmax)
    integer :: nsts(nstsmx)
    integer :: nstsr(2,nstmax)
    integer :: ntfx(ntfxmx)

    integer :: iebal
    integer :: iern1
    integer :: irdxc3
    integer :: ixrn1
    integer :: kbt
    integer :: kdim
    integer :: kelect
    integer :: khydr
    integer :: kmt
    integer :: km1
    integer :: ko2gaq
    integer :: kwater
    integer :: kxt
    integer :: kx1
    integer :: narn1
    integer :: narn2
    integer :: nbt
    integer :: nbw
    integer :: nelect
    integer :: nern1
    integer :: nern2
    integer :: no2gaq
    integer :: nredox
    integer :: nst
    integer :: ntfxt

    logical :: qredox
    logical :: qvansc
    logical :: q6mode

    character(len=48) :: uspec(nstmax)
    character(len=24) :: ugexmo(netmax)

    real(kind=8) :: aamatr(kmax,kmax)
    real(kind=8) :: bpx(ibpxmx,nxtmax)
    real(kind=8) :: cdrs(ndrsmx)
    real(kind=8) :: cdrtw(nstmax)
    real(kind=8) :: cdrw(nstmax)
    real(kind=8) :: cjbasp(nbtmax)
    real(kind=8) :: cnufac(nstmax)
    real(kind=8) :: conc(nstmax)
    real(kind=8) :: csts(nstsmx)
    real(kind=8) :: dlogxw(nbtmax)
    real(kind=8) :: mosp(nstmax)
    real(kind=8) :: tfx(ntfxmx)
    real(kind=8) :: weight(nstmax)
    real(kind=8) :: xbar(nstmax)
    real(kind=8) :: zchar(nstmax)

    real(kind=8) :: al10
    real(kind=8) :: eps100
    real(kind=8) :: omega
    real(kind=8) :: xbarw
    real(kind=8) :: xbarwc

    ! Local pseudo-module.
    ! This group contains scratch arrays for computing data (the
    ! quantities 1/2.3026 d m(A)/d m(B), where A is a non-basis
    ! generic ion exchanger species and B is a basis species. Here
    ! x(B) replaces m(B) if B represents solvent water. These
    ! quantities are required to compute the Jacobian matrix. The
    ! quantities for all A sharing the same site of an exchanger
    ! must be determined simultaneously. To do this, it is necessary
    ! to solve a (usually small) matrix equation for each site.
    ! If there are n non-basis exchanger species for a given site,
    ! the matrix is n x n.
    integer :: isv_ietmax

    SAVE isv_ietmax

    integer, dimension(:), allocatable :: iimgex
    integer, dimension(:), allocatable :: ipvgex

    SAVE iimgex,ipvgex

    real(kind=8), dimension(:), allocatable :: dmlge
    real(kind=8), dimension(:), allocatable :: rhsgex

    SAVE dmlge,rhsgex

    real(kind=8), dimension(:,:), allocatable :: aamgex
    real(kind=8), dimension(:,:), allocatable :: ggmgex

    SAVE aamgex,ggmgex

    ! Local variable declarations.
    integer :: i
    integer :: ix
    integer :: iy
    integer :: j
    integer :: jlen
    integer :: jfl
    integer :: j2
    integer :: kcol
    integer :: krow
    integer :: n
    integer :: nb
    integer :: nbb
    integer :: nb2
    integer :: ne
    integer :: nmax
    integer :: np
    integer :: np1
    integer :: nr1
    integer :: ns
    integer :: nsc
    integer :: nsi
    integer :: nsj
    integer :: nss
    integer :: ns1
    integer :: ns2
    integer :: nx

    integer :: ilnobl

    character(len=56) :: uspn56

    real(kind=8) :: aax
    real(kind=8) :: aznsi
    real(kind=8) :: azns1
    real(kind=8) :: cvansf
    real(kind=8) :: cxcc
    real(kind=8) :: cxf
    real(kind=8) :: cxi1
    real(kind=8) :: cxjc
    real(kind=8) :: cxji
    real(kind=8) :: cxor
    real(kind=8) :: cx2r
    real(kind=8) :: cx21
    real(kind=8) :: stx
    real(kind=8) :: sumx
    real(kind=8) :: xx
    real(kind=8) :: zp

    real(kind=8) :: coefdr
    real(kind=8) :: coefst

    ! Note: the following statements don't really do anything except
    ! cause the compiler not to complain that kelect is not used.
    ix = kelect
    iy = ix
    kelect = iy

    if (netmax .gt. 0) then
        ! Have one or more generic ion exchange phases in the current
        ! model. Allocate or reallocate scratch arrays as needed.
        if (.not.ALLOCATED(dmlge)) then
            ! The arrays are not allocated. Zero the saved array size
            ! variable. Note that only one array is tested to see if it
            ! is allocated. It is assumed that all arrays in the set are
            ! either allocated or not.
            isv_ietmax = 0
        else
            ! The arrays are allocated. Check to see if the array size
            ! variable hs changed. If so, deallocate the corresponding
            ! arrays and zero the corresponding saved size variable.
            if (ietmax .ne. isv_ietmax) then
                DEALLOCATE (dmlge)
                DEALLOCATE (rhsgex)
                DEALLOCATE (aamgex,ggmgex)
                DEALLOCATE (iimgex,ipvgex)
                isv_ietmax = 0
            end if
        end if

        ! At this point, the saved array size value is zero if the
        ! corresponding arrays need to be allocated.
        if (isv_ietmax .eq. 0) then
            ALLOCATE (dmlge(ietmax))
            ALLOCATE (rhsgex(ietmax))
            ALLOCATE (aamgex(ietmax,ietmax),ggmgex(ietmax,ietmax))
            ALLOCATE (iimgex(ietmax),ipvgex(ietmax))
            isv_ietmax = ietmax
        end if

        ! Zero these arrays.
        call initaz(dmlge,ietmax)
        call initaz(rhsgex,ietmax)
        nmax = ietmax*ietmax
        call initaz(aamgex,nmax)
        call initaz(ggmgex,nmax)
        call initiz(iimgex,ietmax)
        call initiz(ipvgex,ietmax)
    end if

    ! Zero the matrix aamatr.
    do i = 1,kdim
        do j = 1,kdim
            aamatr(i,j) = 0.
        end do
    end do

    ! Zero the cnufac array.
    do ns = 1,nst
        cnufac(ns) = 0.
    end do

    ! Set up the cnufac array.
    do ns = narn1,narn2
        nr1 = ndrsr(1,ns)

        if (jflag(ns) .eq. 30) then
            cnufac(ns) = -conc(ns)/cdrs(nr1)
        end if
    end do

    do ns = nern1,nern2
        nr1 = ndrsr(1,ns)

        if (jflag(ns) .eq. 30) then
            cnufac(ns) = -conc(ns)/cdrs(nr1)
        end if
    end do

    if (q6mode) then
        go to 210
    end if

    ! Write the matrix for EQ3NR.
    ! Write rows 1 through kbt.
    do krow = 1,kbt
        nb = iindx1(krow)
        nsi = nbasp(nb)
        jfl = jflag(nsi)

        if (nsi.ge.narn1 .and. nsi.le.narn2) then
            if (krow.eq.kwater .and. jfl.eq.0) then
                ! Water, mole fraction equation. Note that if this
                ! block is executed, solvent water is a basis species,
                ! log x(w) is a primary iteration variable, and that
                ! cnufac(narn1) = 0.
                xx = xbarwc/omega

                do kcol = 1,kbt
                    if (kcol .eq. kwater) then
                        ! Have the column for water itself.
                        sumx = 0.

                        do nss = narn1,narn2
                            nsc = jcsort(nss)
                            sumx = sumx + cdrw(nsc)*cnufac(nsc)
                        end do

                        aamatr(krow,kcol) = -xx*sumx - 1.0
                    else
                        nbb = iindx1(kcol)
                        nsj = nbasp(nbb)

                        if (nsj.gt.narn1 .and. nsj.le.narn2) then
                            ! Have a column for an aqueous solute basis species.
                            sumx = 0.

                            do nss = narn1,narn2
                                nsc = jcsort(nss)

                                ! Calling sequence substitutions:
                                !   nsc for ns
                                !   nsj for nse
                                cxjc = coefdr(cdrs,ndrs,ndrsmx,ndrsr,nsj,nsc,nstmax)
                                sumx = sumx + cxjc*cnufac(nsc)
                            end do

                            aamatr(krow,kcol) = -xx*(conc(nsj) + sumx)
                        end if
                    end if
                end do
            else if (nb .eq. iebal) then
                ! Charge balance (charge adjustment).
                do ns = narn1,narn2
                    weight(ns) = zchar(ns)
                end do

                do ns = nern1,nern2
                    weight(ns) = zchar(ns)
                end do

                call balcon(aamatr,aamgex,al10,cdrs,cjbasp,cnufac,conc,dmlge,eps100,ggmgex,iern1,ietmax,iimgex,iindx1,ixbasp,jcsort,jern1,jern2,jetmax,jjsort,jsitex,kbt,kmax,krow,narn1,narn2,nbasp,nbtmax,ndrs,ndrsmx,ndrsr,nern1,nern2,netmax,noutpt,nphasx,nstmax,nttyo,rhsgex,uspec,weight,xbar)
            else if (jfl .eq. 17) then
                ! Log activity combination.
                ns1 = ncosp(nb)
                kcol = kction(nb)
                aamatr(krow,krow) = -1.
                aznsi = abs(zchar(nsi))
                azns1 = abs(zchar(ns1))
                zp = zchar(nsi)*zchar(ns1)

                if (zp .lt. 0.) then
                    aamatr(krow,kcol) = -aznsi/azns1
                else
                    aamatr(krow,kcol) = aznsi/azns1
                end if
            else if (jfl .eq. 18) then
                ! Log mean activity.
                ns1 = ncosp(nb)
                kcol = kction(nb)
                aamatr(krow,krow) = -1.
                aznsi = abs(zchar(nsi))
                azns1 = abs(zchar(ns1))
                aamatr(krow,kcol) = -aznsi/azns1
            else if (jfl .eq. 21) then
                ! pHCl.
                ns1 = ncosp(nb)
                kcol = kction(nb)
                aamatr(krow,krow) = -1.
                aznsi = abs(zchar(nsi))
                azns1 = abs(zchar(ns1))
                zp = zchar(nsi)*zchar(ns1)

                if (zp .lt. 0.) then
                    aamatr(krow,kcol) = -aznsi/azns1
                else
                    aamatr(krow,kcol) = aznsi/azns1
                end if
            else if (jfl .eq. 25) then
                ! Heterogeneous equilibrium.
                ns1 = ncosp(nb)

                ! Calling sequence substitutions:
                !   nsi for nse
                !   ns1 for ns
                cxi1 = coefdr(cdrs,ndrs,ndrsmx,ndrsr,nsi,ns1,nstmax)

                do kcol = 1,kbt
                    if (kcol.ne.krow .and. kcol.ne.kwater) then
                        nb2 = iindx1(kcol)
                        ns2 = nbasp(nb2)

                        ! Calling sequence substitutions:
                        !   ns2 for nse
                        cx21 = coefdr(cdrs,ndrs,ndrsmx,ndrsr,ns2,ns1,nstmax)
                        aamatr(krow,kcol) = -cx21/cxi1
                    end if
                end do

                aamatr(krow,krow) = -1.0
            else if (nsi .eq. no2gaq) then
                ! Log fO2.
                if (irdxc3 .eq. 0) then
                    ! Identity.
                    aamatr(krow,ko2gaq) = 1.
                else if (irdxc3  .eq. -1) then
                    ! Eh residual (Note- if a pe- value was input, it has been
                    ! converted to an Eh value by EQ3NR/setup.f).
                    aamatr(krow,kwater) = 2.0
                    aamatr(krow,khydr) = aamatr(krow,khydr) - 4.0
                    aamatr(krow,krow) = aamatr(krow,krow) - 1.0
                else if (irdxc3  .eq. 1) then
                    ! Cross-linking (homogeneous aqueous redox) equilibrium.
                    aamatr(krow,ko2gaq) = -1.

                    ! Calling sequence substitutions:
                    !   no2gaq for nse
                    !   nredox for ns
                    cxor = coefdr(cdrs,ndrs,ndrsmx,ndrsr,no2gaq,nredox,nstmax)

                    do kcol = 1,kbt
                        if (kcol.ne.ko2gaq .and. kcol.ne.kwater) then
                            nb2 = iindx1(kcol)
                            ns2 = nbasp(nb2)

                            ! Calling sequence substitutions:
                            !   ns2 for nse
                            !   nredox for ns
                            cx2r = coefdr(cdrs,ndrs,ndrsmx,ndrsr,ns2,nredox,nstmax)
                            aamatr(krow,kcol) = -cx2r/cxor
                        end if
                    end do
                else
                    ! Error.
                    write (noutpt,1000) irdxc3
                    write (nttyo,1000) irdxc3
1000 format(/' * Error - (EQLIB/matrix) Programming error',' trap:',/7x,'Have encountered illegal irdxc3 value = ',i5,'.')

                    stop
                end if
            else if (jfl.ge.0 .and. jfl.le.3) then
                ! Mass balance.
                do ns = narn1,narn2
                    weight(ns) = coefst(csts,nsts,nstsmx,nstsr,nb,ns,nstmax)
                end do

                do ns = nern1,nern2
                    weight(ns) = 0.
                end do

                call balcon(aamatr,aamgex,al10,cdrs,cjbasp,cnufac,conc,dmlge,eps100,ggmgex,iern1,ietmax,iimgex,iindx1,ixbasp,jcsort,jern1,jern2,jetmax,jjsort,jsitex,kbt,kmax,krow,narn1,narn2,nbasp,nbtmax,ndrs,ndrsmx,ndrsr,nern1,nern2,netmax,noutpt,nphasx,nstmax,nttyo,rhsgex,uspec,weight,xbar)
            else if (jfl.ge.7 .and. jfl.le.11) then
                ! Alkalinity balance.
                ! Calling sequence substitutions.
                !   weight for array
                !   nstmax for nmax
                do ns = 1,nst
                    weight(ns) = 0.
                end do

                do n = 1,ntfxt
                    ns = ntfx(n)
                    weight(ns) = tfx(n)
                end do

                call balcon(aamatr,aamgex,al10,cdrs,cjbasp,cnufac,conc,dmlge,eps100,ggmgex,iern1,ietmax,iimgex,iindx1,ixbasp,jcsort,jern1,jern2,jetmax,jjsort,jsitex,kbt,kmax,krow,narn1,narn2,nbasp,nbtmax,ndrs,ndrsmx,ndrsr,nern1,nern2,netmax,noutpt,nphasx,nstmax,nttyo,rhsgex,uspec,weight,xbar)
            else if (jfl .eq. 16) then
                ! Log activity.
                aamatr(krow,krow) = 1.0
            else if (jfl.eq.19 .or. jfl.eq.20) then
                ! pX (including pH).
                aamatr(krow,krow) = 1.0
            else if (jfl.eq.22 .or. jfl.eq.23) then
                ! pmX (including pmH).
                aamatr(krow,krow) = 1.0
            else if (jfl .eq. 27) then
                ! Aqueous homogeneous equilibrium.
                ns1 = nsi

                ! Calling sequence substitutions:
                !   nsi for nse
                !   ns1 for ns
                cxi1 = coefdr(cdrs,ndrs,ndrsmx,ndrsr,nsi,ns1,nstmax)

                do kcol = 1,kbt
                    if (kcol.ne.krow .and. kcol.ne.kwater) then
                        nb2 = iindx1(kcol)
                        ns2 = nbasp(nb2)

                        ! Calling sequence substitutions:
                        !   ns2 for nse
                        !   ns1 for ns
                        cx21 = coefdr(cdrs,ndrs,ndrsmx,ndrsr,ns2,ns1,nstmax)
                        aamatr(krow,kcol) = -cx21/cxi1
                    end if
                end do

                aamatr(krow,krow) = -1.0
            else
                ! Have found a bad jflag value.
                ! Calling sequence substitutions:
                !   uspec(nsi) for unam48
                call fmspnx(jlen,uspec(nsi),uspn56)
                write (noutpt,1010) jfl,uspn56(1:jlen)
                write (nttyo,1010) jfl,uspn56(1:jlen)
1010 format(/' * Error - (EQLIB/matrix) Programming error trap:',/7x,'Have encountered a bad jflag value of ',i4,' for',/7x,a,'.')

                stop
            end if
        else if (nsi.ge.nern1 .and. nsi.le.nern2) then
            ! Generic ion exchange species.
            if (jfl .eq. 0) then
                ! Mass balance.
                do ns = narn1,narn2
                    weight(ns) = 0.
                end do

                do ns = nern1,nern2
                    weight(ns) = coefst(csts,nsts,nstsmx,nstsr,nb,ns,nstmax)
                end do

                call balcon(aamatr,aamgex,al10,cdrs,cjbasp,cnufac,conc,dmlge,eps100,ggmgex,iern1,ietmax,iimgex,iindx1,ixbasp,jcsort,jern1,jern2,jetmax,jjsort,jsitex,kbt,kmax,krow,narn1,narn2,nbasp,nbtmax,ndrs,ndrsmx,ndrsr,nern1,nern2,netmax,noutpt,nphasx,nstmax,nttyo,rhsgex,uspec,weight,xbar)
            else
                ! Have found a bad jflag value.
                ! Calling sequence substitutions:
                !   uspec(nsi) for unam48
                call fmspnx(jlen,uspec(nsi),uspn56)
                write (noutpt,1010) jfl,uspn56(1:jlen)
                write (nttyo,1010) jfl,uspn56(1:jlen)
            end if
        end if
    end do

    go to 999

210 continue

    ! Write the matrix for EQ6.
    ! Compute the dlogxw array (d log x(w)/d log m(s')).
    call gdlgxw(cdrs,cjbasp,cnufac,conc,dlogxw,eps100,ixbasp,jcsort,jflag,narn1,narn2,nbasp,nbt,nbtmax,nbw,ndrs,ndrsmx,ndrsr,nern1,nern2,noutpt,nstmax,nttyo,omega,xbar,xbarw)

    ! Change the dlogxw array from derivatives with respect to
    ! molalities to ones with respect to numbers of moles.
    ! Currently, the definition is:
    !   dlogxw = d log x(w)/d log m(s'), s' != w
    ! The new definition will be:
    !   dlogxw = d log x(w)/d log n(s'), s' includes w
    ! There is now a derivative with respect to n(w). There was
    ! none with respect to m(w).
    ! Compute dlogxw(w) (d log x(w)/d log n(w)). For basis species
    ! which are aqueous solute species, it happens that:
    !   d log x(w)/d log n(s') = d log x(w)/d log m(s')
    sumx = 0.

    do nb = 1,nbt
        if (nb .ne. nbw) then
            ns = nbasp(nb)
            sumx = sumx + dlogxw(nb)
        end if
    end do

    dlogxw(nbw) = -sumx

    ! Redefine the cnufac array in terms of numbers of moles in place of
    ! concentrations.
    do ns = narn1,narn2
        nr1 = ndrsr(1,ns)

        if (jflag(ns) .eq. 30) then
            cnufac(ns) = -mosp(ns)/cdrs(nr1)
        end if
    end do

    do ns = nern1,narn2
        nr1 = ndrsr(1,ns)

        if (jflag(ns) .eq. 30) then
            cnufac(ns) = -mosp(ns)/cdrs(nr1)
        end if
    end do

    ! Build the rows corresponding to mass balance expressions for
    ! the active basis species in the 'd' set.
    do krow = 1,kbt
        nb = iindx1(krow)
        nsi = nbasp(nb)

        if (nsi.eq.no2gaq .or. nsi.eq.nelect) then
            if (.not.qredox) then
                aamatr(krow,krow) = 1.0
                go to 250
            end if
        end if

        do ns = 1,nst
            weight(ns) = 0.
        end do

        ! Mass balance.
        do ns = narn1,narn2
            weight(ns) = coefst(csts,nsts,nstsmx,nstsr,nb,ns,nstmax)
        end do

        do ns = nern1,nern2
            weight(ns) = coefst(csts,nsts,nstsmx,nstsr,nb,ns,nstmax)
        end do

        do kcol = km1,kxt
            ns = iindx1(kcol)
            weight(ns) = coefst(csts,nsts,nstsmx,nstsr,nb,ns,nstmax)
        end do

        ! Build column kwater.
        aax = 0.

        do nss = narn1,narn2
            nsc = jcsort(nss)
            aax = aax - weight(nsc)*cnufac(nsc)    *(cdrtw(nsc) + cdrw(nsc)*dlogxw(nbw))
        end do

        aax = aax + weight(narn1)*mosp(narn1)
        aamatr(krow,kwater) = al10*aax

        ! Build columns 1 to kbt, except kwater.
        do kcol = 1,kbt
            if (kcol .ne. kwater) then
                nbb = iindx1(kcol)
                nsj = nbasp(nbb)

                ! Does the operational basis species associated with the
                ! current column belong to a Vanselow ion exchange phase?
                qvansc = .false.

                if (nsj.ge.nern1 .and. nsj.le.nern2) then
                    np = nphasx(nsj)
                    ne = np - iern1 + 1
                    j2 = ilnobl(ugexmo(ne))
                    qvansc = ugexmo(ne)(1:j2).eq.'Vanselow' .or.        ugexmo(ne)(1:9).eq.'Vanselow-'
                end if

                aax = 0.

                do nss = narn1,narn2
                    nsc = jcsort(nss)

                    ! Calling sequence substitutions:
                    !   nsc for ns
                    !   nsj for nse
                    cxjc = coefdr(cdrs,ndrs,ndrsmx,ndrsr,nsj,nsc,nstmax)
                    aax = aax + weight(nsc)*cnufac(nsc)        *(cxjc + cdrw(nsc)*dlogxw(nbb))
                end do

                ! KK
                do nss = nern1,nern2
                    nsc = jcsort(nss)

                    ! Calling sequence substitutions:
                    !   nsc for ns
                    !   nsj for nse
                    cxjc = coefdr(cdrs,ndrs,ndrsmx,ndrsr,nsj,nsc,nstmax)

                    ! KKKKKKKKK
                    cxf = cxjc

                    if (qvansc) then
                        ! The "Vanselow factor" cvansf calculated below is
                        ! a measure of the extent to which the mole fractions
                        ! of two exchanger species (the one, non-basis, the
                        ! other, basis) in a mass action equation for the
                        ! Vanselow ion exchange model can not be simply
                        ! replaced by corresponding molalities. This factor
                        ! is non-zero if the two exchanging ions have different
                        ! charges.
                        nr1 = ndrsr(1,nsc)
                        cxcc = cdrs(nr1)
                        cvansf = -cxcc + (cxjc + cxcc)*xbar(nsj)
                        cxf = cvansf
                    end if

                    ! KKKKKKKKK
                    aax = aax + weight(nsc)*cnufac(nsc)        *(cxf + cdrw(nsc)*dlogxw(nbb))

                    ! Old          aax = aax + weight(nsc)*cnufac(nsc)
                    ! Old $        *(cxjc + cdrw(nsc)*dlogxw(nbb))
                end do

                aax = aax + weight(nsj)*mosp(nsj)
                aamatr(krow,kcol) = al10*aax
            end if
        end do

        ! Build the mineral columns, including those for solid solution
        ! species.
        do kcol = km1,kxt
            nsj = iindx1(kcol)
            aamatr(krow,kcol) = al10*weight(nsj)*mosp(nsj)
        end do

250 continue
    end do

    ! Build the rows corresponding to mass action expressions
    ! for pure minerals.
    do krow = km1,kmt
        nsi = iindx1(krow)

        do kcol = 1,kbt
            nbb = iindx1(kcol)
            nsj = nbasp(nbb)

            if (nsj .eq. narn1) then
                aamatr(krow,kcol) = -cdrtw(nsi) + cdrw(nsi)*dlogxw(nbb)
            else
                ! Calling sequence substitutions:
                !   nsi for ns
                !   nsj for nse
                cxji = coefdr(cdrs,ndrs,ndrsmx,ndrsr,nsj,nsi,nstmax)
                aamatr(krow,kcol) = cxji + cdrw(nsi)*dlogxw(nbb)
            end if
        end do
    end do

    ! Build the rows corresponding to mass action expressions
    ! for species belonging to solid solutions.
    if (kxt .ge. kx1) then
        do krow = kx1,kxt
            nsi = iindx1(krow)
            np = ipndx1(krow)
            nx = np - ixrn1 + 1
            stx = bpx(1,nx)

            ! Build columns 1 through kbt.
            do kcol = 1,kbt
                nbb = iindx1(kcol)
                nsj = nbasp(nbb)

                if (nsj .eq. narn1) then
                    aamatr(krow,kcol) = -cdrtw(nsi) + cdrw(nsi)*dlogxw(nbb)
                else
                    ! Calling sequence substitutions:
                    !   nsi for ns
                    !   nsj for nse
                    cxji = coefdr(cdrs,ndrs,ndrsmx,ndrsr,nsj,nsi,nstmax)
                    aamatr(krow,kcol) = cxji + cdrw(nsi)*dlogxw(nbb)
                end if
            end do

            ! Build columns kx1 through kxt.
            nr1 = ndrsr(1,nsi)
            cxcc = cdrs(nr1)

            do kcol = kx1,kxt
                ns1 = iindx1(kcol)
                np1 = ipndx1(kcol)

                if (np1 .eq. np) then
                    if (ns1 .eq. nsi) then
                        aamatr(krow,kcol) = cxcc*stx*( 1. - xbar(nsi) )
                    else
                        aamatr(krow,kcol) = - cxcc*stx*xbar(ns1)
                    end if
                end if
            end do
        end do
    end if

999 continue
end subroutine matrix