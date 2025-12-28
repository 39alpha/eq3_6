subroutine raff(acflg,actlg,afcnst,affp,afrc1,bpx,cdrs,cgexj,ibpxmx,ibpxt,iern1,ietmax,iktmax,ixrn1,ixrn2,jcode,jern1,jern2,jetmax,jflag,jgext,jpflag,jsflag,jsol,ncmpr,ndrs,ndrsmx,ndrsr,nertmx,net,netmax,ngext,noutpt,nptmax,nrct,nrctmx,nrndex,nstmax,nttyo,nxridx,nxrtmx,nxtmax,rxbar,uphase,uspec,wfac,xbar,xbarlg,xgers,xlks)
    !! This subroutine calculates the affinities of irreversible
    !! reactions (afrc1).
    !! This subroutine is called by:
    !!   EQ6/path.f
    !! Principal input:
    !! Principal output:
    implicit none

    include 'eqlib/eqlpar.h'

    ! Calling sequence variable declarations.
    integer :: ibpxmx
    integer :: ietmax
    integer :: iktmax
    integer :: jetmax
    integer :: ndrsmx
    integer :: nertmx
    integer :: netmax
    integer :: nptmax
    integer :: nrctmx
    integer :: nstmax
    integer :: nxrtmx
    integer :: nxtmax

    integer :: noutpt
    integer :: nttyo

    integer :: ibpxt(nxtmax)
    integer :: jcode(nrctmx)
    integer :: jern1(jetmax,netmax)
    integer :: jern2(jetmax,netmax)
    integer :: jflag(nstmax)
    integer :: jgext(netmax)
    integer :: jpflag(nptmax)
    integer :: jsflag(nstmax)
    integer :: jsol(nxtmax)
    integer :: ncmpr(2,nptmax)
    integer :: ndrs(ndrsmx)
    integer :: ndrsr(2,nstmax)
    integer :: ngext(jetmax,netmax)
    integer :: nrndex(nrctmx)
    integer :: nxridx(nrctmx)

    integer :: iern1
    integer :: ixrn1
    integer :: ixrn2
    integer :: net
    integer :: nrct

    character(len=48) :: uspec(nstmax)
    character(len=24) :: uphase(nptmax)

    real(kind=8) :: acflg(nstmax)
    real(kind=8) :: actlg(nstmax)
    real(kind=8) :: affp(nptmax)
    real(kind=8) :: afrc1(nrctmx)
    real(kind=8) :: bpx(ibpxmx,nxtmax)
    real(kind=8) :: cdrs(ndrsmx)
    real(kind=8) :: cgexj(jetmax,netmax)
    real(kind=8) :: rxbar(iktmax,nxrtmx)
    real(kind=8) :: wfac(iktmax,nxtmax)
    real(kind=8) :: xbar(nstmax)
    real(kind=8) :: xbarlg(nstmax)
    real(kind=8) :: xgers(ietmax,jetmax,nertmx)
    real(kind=8) :: xlks(nstmax)

    real(kind=8) :: afcnst

    ! Local variable declarations with fixed global dimensioning.
    integer :: jflgx(iet_par,net_par)
    integer :: jsflgx(iet_par,net_par)

    real(kind=8) :: acflge(iet_par,jet_par)
    real(kind=8) :: actlge(iet_par,jet_par)
    real(kind=8) :: xbarle(iet_par,jet_par)
    real(kind=8) :: xbare(iet_par,jet_par)

    ! Local variable declarations with variable global dimensioning.
    integer :: isv_iktmax

    SAVE isv_iktmax

    real(kind=8), dimension(:), allocatable :: acflgs
    real(kind=8), dimension(:), allocatable :: actlgs
    real(kind=8), dimension(:), allocatable :: xbarls
    real(kind=8), dimension(:), allocatable :: xbars

    SAVE acflgs,actlgs,xbarls,xbars

    ! Local variable declarations.
    integer :: ie
    integer :: ik
    integer :: je
    integer :: jpflgx
    integer :: k
    integer :: ne
    integer :: ner
    integer :: np
    integer :: nrc
    integer :: nr1
    integer :: nr2
    integer :: ns
    integer :: nxr

    real(kind=8) :: af
    real(kind=8) :: affpr
    real(kind=8) :: si
    real(kind=8) :: xx

    real(kind=8) :: tlg

    ! Allocate or reallocate local work arrays as needed.
    if (.not.ALLOCATED(acflgs)) then
        ! Local work arrays are not allocated. Zero the saved
        ! array size variables. Note that only one array is tested
        ! to see if it is allocated. It is assumed that all local
        ! work arrays are either allocated or not.
        isv_iktmax = 0
    else
        ! Local work arrays are allocated. Check to see if any of the
        ! array size variables have changed. If so, deallocate
        ! the corresponding local work arrays and zero the corresponding
        ! saved size variables.
        if (iktmax .ne. isv_iktmax) then
            DEALLOCATE(acflgs,actlgs,xbarls,xbars)
            isv_iktmax = 0
        end if
    end if

    ! At this point, the saved array size values are zero if the
    ! corresponding arrays need to be allocated.
    if (isv_iktmax .eq. 0) then
        ALLOCATE(acflgs(iktmax),actlgs(iktmax))
        ALLOCATE(xbarls(iktmax),xbars(iktmax))
        isv_iktmax = iktmax
    end if

    ! Zero the contents of the local work arrays.
    do k = 1,iktmax
        acflgs(k) = 0.
        actlgs(k) = -99999.
        xbars(k) = 0.
        xbarls(k) = -99999.
    end do

    ! Recall the following jreac reactant status flag conventions:
    !   jreac =  0: set to react
    !   jreac =  1: exhausted
    !   jreac = -1: saturated, but the remaining reactant mass
    !                 continues to react irreversibly
    !   jreac =  2: saturated; the status of any remaining reactant
    !                 mass is changed to that of a product phase
    ! When jreac = -1 or 2, the affinity is zero by definition.
    ! However, the corresponding computed affinities returned from the
    ! present subroutine may be non-zero, owing to finite convergence
    ! tolerances. Avoid calculating finite differences from such
    ! computed affinity values in EQ6/stepfd.f. However, do not set
    ! the corresponding afrc1 values to zero in the present subroutine.
    do nrc = 1,nrct
        if (jcode(nrc) .eq. 0) then
            ! Pure minerals.
            np = nrndex(nrc)
            afrc1(nrc) = -affp(np)
        else if (jcode(nrc) .eq. 1) then
            ! Solid solutions.
            np = nrndex(nrc)
            nxr = nxridx(nrc)

            nr1 = ncmpr(1,np)
            nr2 = ncmpr(2,np)

            ik = 0

            do ns = nr1,nr2
                ik = ik + 1
                acflgs(ik) = acflg(ns)
                actlgs(ik) = actlg(ns)
                xbarls(ik) = xbarlg(ns)
                xbars(ik) = xbar(ns)
                xx = rxbar(ik,nxr)
                xbar(ns) = xx
                xbarlg(ns) = tlg(xx)
            end do

            ! Calling sequence substitutions:
            !   acflg for acflgc
            call lambda(acflg,afcnst,bpx,ibpxmx,ibpxt,iktmax,ixrn1,ixrn2,jsol,ncmpr,noutpt,np,nptmax,nstmax,nttyo,nxtmax,wfac,xbar,xbarlg,uphase,uspec)

            do ns = nr1,nr2
                actlg(ns) = xbarlg(ns) + acflg(ns)
            end do

            affpr = 0.
            ik = 0

            do ns = nr1,nr2
                ik = ik + 1
                xx = xbar(ns)

                if (xx .gt. 0.) then
                    call afcalc(actlg,af,afcnst,cdrs,jflag,jsflag,ndrs,ndrsmx,ndrsr,ns,nstmax,si,xlks)

                    if (af .gt. -9999999.) then
                        affpr = affpr + xx*af
                    else
                        affpr = -9999999.
                        go to 100
                    end if
                end if
            end do

100 continue

            afrc1(nrc) = -affpr

            ik = 0

            do ns = nr1,nr2
                ik = ik + 1
                acflg(ns) = acflgs(ik)
                actlg(ns) = actlgs(ik)
                xbarlg(ns) = xbarls(ik)
                xbar(ns) = xbars(ik)
            end do
        else if (jcode(nrc) .eq. 2) then
            ! Special reactants.
            afrc1(nrc) = 9999999.
        else if (jcode(nrc) .eq. 3) then
            ! Aqueous species.
            afrc1(nrc) = 9999999.
        else if (jcode(nrc) .eq. 4) then
            ! Gases.
            afrc1(nrc) = 9999999.
        else if (jcode(nrc) .eq. 5) then
            ! Generic ion exchangers.
            np = nrndex(nrc)
            ner = nxridx(nrc)
            ne = np - iern1 + 1

            nr1 = ncmpr(1,np)
            nr2 = ncmpr(2,np)

            jpflgx = jpflag(np)
            jpflag(np) = 0

            do je = 1,jgext(ne)
                ns = jern1(je,ne) - 1

                do ie = 1,ngext(je,ne)
                    ns = ns + 1
                    jsflgx(ie,je) = jsflag(ns)
                    jsflag(ns) = 0
                    jflgx(ie,je) = jflag(ns)
                    jflag(ns) = 0
                    acflge(ie,je) = acflg(ns)
                    actlge(ie,je) = actlg(ns)
                    xbarle(ie,je) = xbarlg(ns)
                    xbare(ie,je) = xbar(ns)
                    xx = xgers(ie,je,ner)
                    xbar(ns) = xx
                    xbarlg(ns) = tlg(xx)
                end do
            end do

            ! Calculate activity coefficients. Recall that these may
            ! non-zero even in the case of ideal exchangers, due to
            ! site-mixing.
            ! Calling sequence substitutions:
            !   acflg for acflgc
            call lamgex(acflg,cgexj,jern1,jern2,jetmax,jgext,net,netmax,nstmax,xbarlg)

            do ns = nr1,nr2
                actlg(ns) = xbarlg(ns) + acflg(ns)
            end do

            affpr = 0.

            do je = 1,jgext(ne)
                ns = jern1(je,ne) - 1

                do ie = 1,ngext(je,ne)
                    ns = ns + 1
                    xx = xbar(ns)

                    if (xx .gt. 0.) then
                        call afcalc(actlg,af,afcnst,cdrs,jflag,jsflag,ndrs,ndrsmx,ndrsr,ns,nstmax,si,xlks)

                        if (af .gt. -9999999.) then
                            affpr = affpr + xx*af
                        else
                            affpr = -9999999.
                            go to 110
                        end if
                    end if
                end do
            end do

110 continue

            afrc1(nrc) = -affpr

            jpflag(np) = jpflgx

            do je = 1,jgext(ne)
                ns = jern1(je,ne) - 1

                do ie = 1,ngext(je,ne)
                    ns = ns + 1
                    jsflag(ns) = jsflgx(ie,je)
                    jflag(ns) = jflgx(ie,je)
                    acflg(ns) = acflge(ie,je)
                    actlg(ns) = actlge(ie,je)
                    xbarlg(ns) = xbarle(ie,je)
                    xbar(ns) = xbare(ie,je)
                end do
            end do
        end if
    end do
end subroutine raff