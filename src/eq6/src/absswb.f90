subroutine absswb(adhfs,adhfsx,advfs,advfsx,avcnst,axhfs,axhfsx,axlks,axlksx,axvfs,axvfsx,cdrs,cdrtw,cdrsx,cdrw,csts,dhfs,dvfs,eps100,ibswx,iindx1,iodb,ipch,ipchmx,ipcv,ipcvmx,jcsort,jflag,jsflag,kbt,kmax,mosp,mtb,narn1,narn2,narxmx,narxt,nbasp,nbaspd,nbaspx,nbt,nbtmax,nbw,ndrs,ndrsmx,ndrsr,ndrsrx,ndrsx,nelect,nhydr,nodbmx,no2gaq,noutpt,nst,nstmax,nsts,nstsmx,nstsr,nswtch,ntpr,ntprmx,nttyo,presg,press,qbassw,qbswx,tempc,uspec,uzvec1,weight,xhfs,xvfs,xlks)
    !! This subroutine carries out automatic basis switching to optimize
    !! optimize the basis set after Newton-Raphson iteration has
    !! converged at the current point of reaction progress. Such
    !! optimization tends to improve the code numerics (finite
    !! difference stability, matrix conditioning, etc.) This subroutine
    !! is similar to EQLIB/absswa.f, which carries out automatic basis
    !! switching as a pre-Newton-Raphson optimization technique.
    !! This subroutine is called by:
    !!   EQ6/path.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: ipchmx
    integer :: ipcvmx
    integer :: kmax
    integer :: narxmx
    integer :: nbtmax
    integer :: ndrsmx
    integer :: nodbmx
    integer :: nstmax
    integer :: nstsmx
    integer :: ntprmx

    integer :: noutpt
    integer :: nttyo

    integer :: ibswx(nbtmax)
    integer :: iindx1(kmax)
    integer :: iodb(nodbmx)
    integer :: jcsort(nstmax)
    integer :: jflag(nstmax)
    integer :: jsflag(nstmax)
    integer :: narxt(ntprmx)
    integer :: nbasp(nbtmax)
    integer :: nbaspd(nbtmax)
    integer :: nbaspx(nbtmax)
    integer :: ndrs(ndrsmx)
    integer :: ndrsr(2,nstmax)
    integer :: ndrsrx(2,nstmax)
    integer :: ndrsx(ndrsmx)
    integer :: nsts(nstsmx)
    integer :: nstsr(2,nstmax)

    integer :: ipch
    integer :: ipcv
    integer :: kbt
    integer :: narn1
    integer :: narn2
    integer :: nbt
    integer :: nbw
    integer :: nelect
    integer :: nhydr
    integer :: no2gaq
    integer :: nst
    integer :: nswtch
    integer :: ntpr

    logical :: qbassw
    logical :: qbswx

    character(len=48) :: uspec(nstmax)
    character(len=48) :: uzvec1(kmax)

    real(kind=8) :: adhfs(narxmx,ntprmx,ipchmx,nstmax)
    real(kind=8) :: adhfsx(narxmx,ntprmx,ipchmx,nstmax)
    real(kind=8) :: advfs(narxmx,ntprmx,ipcvmx,nstmax)
    real(kind=8) :: advfsx(narxmx,ntprmx,ipcvmx,nstmax)
    real(kind=8) :: axhfs(narxmx,ntprmx,nstmax)
    real(kind=8) :: axhfsx(narxmx,ntprmx,nstmax)
    real(kind=8) :: axlks(narxmx,ntprmx,nstmax)
    real(kind=8) :: axlksx(narxmx,ntprmx,nstmax)
    real(kind=8) :: axvfs(narxmx,ntprmx,nstmax)
    real(kind=8) :: axvfsx(narxmx,ntprmx,nstmax)
    real(kind=8) :: cdrs(ndrsmx)
    real(kind=8) :: cdrsx(ndrsmx)
    real(kind=8) :: cdrtw(nstmax)
    real(kind=8) :: cdrw(nstmax)
    real(kind=8) :: csts(nstsmx)
    real(kind=8) :: dhfs(ipchmx,nstmax)
    real(kind=8) :: dvfs(ipcvmx,nstmax)
    real(kind=8) :: mosp(nstmax)
    real(kind=8) :: mtb(nbtmax)
    real(kind=8) :: weight(nstmax)
    real(kind=8) :: xhfs(nstmax)
    real(kind=8) :: xlks(nstmax)
    real(kind=8) :: xvfs(nstmax)

    real(kind=8) :: avcnst
    real(kind=8) :: presg
    real(kind=8) :: press
    real(kind=8) :: tempc

    real(kind=8) :: eps100

    ! Local variable declarations.
    integer :: jlen1
    integer :: jlen2
    integer :: jlen3
    integer :: krow
    integer :: nb
    integer :: nb1
    integer :: nb2
    integer :: ns
    integer :: nse
    integer :: nsi
    integer :: nsj
    integer :: ns1
    integer :: ns2

    character(len=56) :: usp156
    character(len=56) :: usp256
    character(len=56) :: usp356

    real(kind=8) :: api
    real(kind=8) :: apj
    real(kind=8) :: cx
    real(kind=8) :: fx1
    real(kind=8) :: fx2
    real(kind=8) :: rx
    real(kind=8) :: wb1
    real(kind=8) :: wb2
    real(kind=8) :: wsi

    real(kind=8) :: coefdr
    real(kind=8) :: coefst

    ! Initialize the number of basis switches made.
    nswtch = 0

    ! Pick candidates for automatic basis switching. The indices of
    ! candidates are stored in the array ibswx. Note that the method
    ! employed here is similar to that in EQLIB/abswpk.f. However,
    ! there are some differences. Here, basis switching is not being
    ! used as a method of pre-Newton-Raphson optimization. Hence, there
    ! are no residual functions used in defining criteria for making
    ! a switch. the relative contributions to mass balances are used
    ! instead.
    qbswx = .false.

    do nb = 1,nbt
        ibswx(nb) = 0
        nse = nbaspd(nb)
        nsj = nbasp(nb)

        if (nse.ne.narn1 .and. nse.ne.nhydr .and. nse.ne.no2gaq    .and.nse.ne.nelect) then
            ! If the current species is not H2O, H+, O2(g,aq), or e-,
            ! consider a switch.
            ! Get weights for the current mass balance.
            do ns = narn1,narn2
                weight(ns) = coefst(csts,nsts,nstsmx,nstsr,nb,ns,nstmax)
            end do

            ! Screen out species in the mass balance that are not linked
            ! by the current reaction with the current basis species.
            ! Do this by seting the corresponding weights to zero.
            do ns = narn1,narn2
                if (jflag(ns) .eq. 30) then
                    ! Calling sequence substitutions:
                    !   nsj for nse
                    cx = coefdr(cdrs,ndrs,ndrsmx,ndrsr,nsj,ns,nstmax)

                    if (cx .eq. 0.) then
                        weight(ns) = 0.
                    end if
                end if
            end do

            ! Find a candidate.
            call fbassw(jcsort,jflag,mosp,narn1,narn2,nse,nsi,nsj,nstmax,weight,wsi)

            ! If a candidate was found, apply a filter.
            if (nsi .gt. 0) then
                ! Use the current candidate species if its contribution
                ! to the mass balance exceeds that of the current basis
                ! species by a reasonable factor, here 10.
                api = weight(nsi)*mosp(nsi)
                apj = weight(nsj)*mosp(nsj)
                rx = api/apj

                if (rx .gt. 10.) then
                    ibswx(nb) = nsi
                    qbswx = .true.
                end if
            end if
        end if
    end do

    ! If there are no candidates for basis switching, exit the
    ! present subroutine.
    if (.not.qbswx) then
        go to 999
    end if

    if (iodb(1) .ge. 2) then
        write (noutpt,1000)
1000 format(/5x,'--- Candidate Basis Switches ---',/)

        do nb = 1,nbt
            nse = nbaspd(nb)
            nsj = nbasp(nb)
            nsi = ibswx(nb)

            if (nsi .gt. 0) then
                ! Calling sequence substitutions:
                !   jlen1 for jlen
                !   uspec(nsj) for unam48
                !   usp156 for uspn56
                call fmspnx(jlen1,uspec(nsj),usp156)

                ! Calling sequence substitutions:
                !   jlen2 for jlen
                !   uspec(nsi) for unam48
                !   usp256 for uspn56
                call fmspnx(jlen2,uspec(nsi),usp256)

                ! Calling sequence substitutions:
                !   jlen3 for jlen
                !   uspec(nse) for unam48
                !   usp356 for uspn56
                call fmspnx(jlen3,uspec(nse),usp356)

                write (noutpt,1010) usp156(1:jlen1),usp256(1:jlen2),usp356(1:jlen3)
1010 format(/3x,'Could replace ',a,' in the active basis set',' with',/3x,a,' as the species associated with the mass',' balance',/3x,'of ',a,'.')
            end if
        end do

        write (noutpt,1020)
1020 format(1x)
    end if

    ! Resolve any conflicts in candidate basis switches. The algorithm
    ! used here is similar to but not identical to that employed by
    ! EQLIB/gabswx.f, which resolves conflicts in the context of using
    ! automatic basis switching as a form of pre-Newton-Raphson
    ! optimization.
    do nb1 = 1,nbt - 1
        ns1 = ibswx(nb1)

        if (ns1 .gt. 0) then
            do nb2 = nb1 + 1,nbt
                ns2 = ibswx(nb2)

                if (ns2 .eq. ns1) then
                    ! Calling sequence substitutions:
                    !   nb1 for nb
                    wb1 = coefst(csts,nsts,nstsmx,nstsr,nb1,ns,nstmax)
                    fx1 = wb1*mosp(ns1)/mtb(nb1)

                    ! Calling sequence substitutions:
                    !   nb2 for nb
                    wb2 = coefst(csts,nsts,nstsmx,nstsr,nb2,ns,nstmax)
                    fx2 = wb2*mosp(ns2)/mtb(nb2)

                    if (fx1 .gt. fx2) then
                        ibswx(nb2) = 0
                    else
                        ibswx(nb1) = 0
                        go to 100
                    end if
                end if
            end do
        end if

100 continue
    end do

    nswtch = 0

    do nb = 1,nbt
        nsi = ibswx(nb)

        if (nsi .gt. 0) then
            nswtch = nswtch + 1
        end if
    end do

    if (nswtch .le. 0) then
        if (iodb(1) .ge. 2) then
            write (noutpt,1070)
        end if

1070 format(/10x,'No switches will be made.',/)

        go to 999
    end if

    if (iodb(1) .ge. 2) then
        write (noutpt,1100)
1100 format(/5x,'--- Automatic Basis Switches ---',/)

        do nb = 1,nbt
            nse = nbaspd(nb)
            nsj = nbasp(nb)
            nsi = ibswx(nb)

            if (nsi .gt. 0) then
                ! Calling sequence substitutions:
                !   jlen1 for jlen
                !   uspec(nsj) for unam48
                !   usp156 for uspn56
                call fmspnx(jlen1,uspec(nsj),usp156)

                ! Calling sequence substitutions:
                !   jlen2 for jlen
                !   uspec(nsi) for unam48
                !   usp256 for uspn56
                call fmspnx(jlen2,uspec(nsi),usp256)

                ! Calling sequence substitutions:
                !   jlen3 for jlen
                !   uspec(nse) for unam48
                !   usp356 for uspn56
                call fmspnx(jlen3,uspec(nse),usp356)

                write (noutpt,1110) usp156(1:jlen1),usp256(1:jlen2),usp356(1:jlen3)
1110 format(/3x,'Will replace ',a,' in the active basis set',' with',/3x,a,' as the species associated with the mass',' balance',/3x,'of ',a,'.')
            end if
        end do
    end if

    call autosw(adhfs,adhfsx,advfs,advfsx,axhfs,axhfsx,axlks,axlksx,axvfs,axvfsx,cdrs,cdrsx,eps100,ibswx,iindx1,ipch,ipchmx,ipcv,ipcvmx,jflag,jsflag,kbt,kmax,narn1,narxmx,nbasp,nbaspd,nbaspx,nbt,nbtmax,nbw,ndrs,ndrsmx,ndrsx,ndrsr,ndrsrx,noutpt,nst,nstmax,ntprmx,nttyo,qbassw,uspec,uzvec1)

    ! Recompute the cdrw array.
    call gcdrw(cdrs,cdrw,narn1,ndrs,ndrsmx,ndrsr,nst,nstmax)

    ! Recompute the cdrtw array.
    call gcdrtw(cdrs,cdrtw,narn1,narn2,ndrs,ndrsmx,ndrsr,nelect,no2gaq,nst,nstmax)

    ! Update the thermodynamic data to correspond to the new
    ! active basis set. First, recompute the log K, etc., data for
    ! the various reactions.
    call evdatr(adhfs,advfs,axhfs,axlks,axvfs,dhfs,dvfs,ipch,ipchmx,ipcv,ipcvmx,narxmx,narxt,nst,nstmax,ntpr,ntprmx,tempc,xhfs,xlks,xvfs)

    ! Then make pressure corrections to these thermodynamic data.
    if (ipcv .ge. 0) then
        call pcorrx(avcnst,dhfs,dvfs,ipch,ipchmx,ipcv,ipcvmx,nbasp,nbt,nbtmax,ndrsr,nst,nstmax,presg,press,xhfs,xlks,xvfs)
    end if

999 continue
end subroutine absswb