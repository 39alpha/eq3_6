subroutine absswa(adhfs,adhfsx,advfs,advfsx,avcnst,axhfs,axhfsx,axlks,axlksx,axvfs,axvfsx,beta,cdrs,cdrsx,cdrtw,cdrw,csts,dhfs,dvfs,efac,eps100,ibswx,iebal,iindx1,iodb,ipch,ipchmx,ipcv,ipcvmx,jcsort,jflag,jsflag,jssort,kbt,kmax,mosp,narn1,narn2,narxmx,narxt,nbasp,nbaspd,nbaspx,nbt,nbtmax,nbw,ncosp,ndrs,ndrsmx,ndrsr,ndrsrx,ndrsx,nelect,nhydr,nodbmx,no2gaq,noutpt,nst,nstmax,nsts,nstsmx,nstsr,nswtch,ntpr,ntprmx,nttyo,presg,press,qbassw,qbswx,q6mode,tempc,uspec,uzvec1,weight,xvfs,xlks,xhfs)
    !! This subroutine carries out automatic basis switching as a
    !! pre-Newton-Raphson optimization technique (to reduce the
    !! magnitude of very large mass balance residuals). It is similar
    !! to EQ6/absswb.f, which carries out automatic basis switching to
    !! optimize the basis set after the system has been solved at the
    !! latest point of reaction progress.
    !! This subroutine is called by:
    !!   EQ3NR/arrset.f
    !!   EQ6/optmzr.f
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
    integer :: jssort(nstmax)
    integer :: narxt(ntprmx)
    integer :: nbasp(nbtmax)
    integer :: nbaspd(nbtmax)
    integer :: nbaspx(nbtmax)
    integer :: ncosp(nbtmax)
    integer :: ndrs(ndrsmx)
    integer :: ndrsr(2,nstmax)
    integer :: ndrsrx(2,nstmax)
    integer :: ndrsx(ndrsmx)
    integer :: nsts(nstsmx)
    integer :: nstsr(2,nstmax)

    integer :: iebal
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
    logical :: q6mode

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
    real(kind=8) :: beta(kmax)
    real(kind=8) :: dhfs(ipchmx,nstmax)
    real(kind=8) :: dvfs(ipcvmx,nstmax)
    real(kind=8) :: efac(nbtmax)
    real(kind=8) :: mosp(nstmax)
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
    integer :: jfl
    integer :: jlen1
    integer :: jlen2
    integer :: jlen3
    integer :: j2
    integer :: nb
    integer :: nb1
    integer :: nb2
    integer :: ncount
    integer :: nse
    integer :: nsi
    integer :: nsj
    integer :: nss
    integer :: ns2

    integer :: ilnobl

    character(len=56) :: usp156
    character(len=56) :: usp256
    character(len=56) :: usp356
    character(len=8) :: ux8

    ! Initialize the number of basis switches made.
    nswtch = 0

    ! Pick candidates for automatic basis switching. This involves
    ! a search similar to the one above, but with a variety of special
    ! constraints. The indices of candidates are stored in the array
    ! ibswx.
    call abswpk(beta,cdrs,csts,efac,ibswx,iebal,iindx1,jcsort,jflag,jssort,kbt,kmax,mosp,narn1,narn2,nbasp,nbaspd,nbt,nbtmax,ndrs,ndrsmx,ndrsr,nelect,nhydr,no2gaq,nstmax,nsts,nstsmx,nstsr,qbswx,q6mode,weight)

    if (.not.qbswx) then
        nswtch = 0
        go to 999
    end if

    if (iodb(3) .ge. 1) then
        write (noutpt,1000)
1000 format(/5x,'--- Candidate Basis Switches ---')

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

    if (.not.q6mode) then
        ! Running EQ3NR. Wipe out any proposed basis switches in which
        ! the species to be switched out is the other species involved
        ! in a jflag = 17 18, or 21 option.
        ncount = 0

        do nb1 = 1,nbt
            nse = nbaspd(nb1)
            nsj = nbasp(nb1)
            jfl = jflag(nse)

            if (jfl.eq.17 .or. jfl.eq.18 .or. jfl.eq.21) then
                ns2 = ncosp(nb1)

                do nb2 = 1,nbt
                    nss = nbasp(nb2)

                    if (ns2 .eq. nss) then
                        if (ibswx(nb2) .ne. 0) then
                            if (iodb(3) .ge. 2) then
                                ! Calling sequence substitutions:
                                !   jlen2 for jlen
                                !   uspec(ns2) for unam48
                                !   usp256 for uspn56
                                call fmspnx(jlen2,uspec(ns2),usp256)

                                ! Calling sequence substitutions:
                                !   jlen3 for jlen
                                !   uspec(nse) for unam48
                                !   usp356 for uspn56
                                call fmspnx(jlen3,uspec(nse),usp356)

                                write (ux8,'(i5)') jfl
                                call lejust(ux8)
                                j2 = ilnobl(ux8)

                                write (noutpt,1050) usp256(1:jlen2),ux8(1:j2),usp356(1:jlen3)
                                write (nttyo,1050) usp256(1:jlen2),ux8(1:j2),usp356(1:jlen3)
1050 format(/" * Note - (EQLIB/absswa) Can't switch ",' the species ',a,/7x,'out of the active basis set',' because it is tied up in a jflag = ',a,' option',/7x,'for ',a,'.')

                                ncount = ncount + 1
                            end if

                            ibswx(nb2) = 0
                        end if

                        go to 100
                    end if
                end do

100 continue
            end if
        end do

        if (ncount .gt. 0) then
            write (noutpt,1020)
        end if
    end if

    ! Resolve any conflicts in candidate basis switches.
    call gabswx(beta,ibswx,iindx1,kbt,kmax,nbt,nbtmax)

    ! Count the number of switches to make.
    nswtch = 0

    do nb = 1,nbt
        nsi = ibswx(nb)

        if (nsi .gt. 0) then
            nswtch = nswtch + 1
        end if
    end do

    if (nswtch .le. 0) then
        if (iodb(3) .ge. 1) then
            write (noutpt,1070)
        end if

1070 format(/10x,'No switches will be made.',/)

        go to 999
    end if

    if (iodb(3) .ge. 1) then
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

        write (noutpt,1020)
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
end subroutine absswa
