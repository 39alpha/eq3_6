subroutine tpadv(aadh,aadhh,aadhv,aaphi,abdh,abdhh,abdhv,abdoth,abdot,abdotv,adadhh,adadhv,adbdhh,adbdhv,adbdth,adbdtv,adh,adhfe,adhh,adhv,adhfs,adhfsd,advfe,advfs,advfsd,afcnst,al10,amu,aslm,aphi,aprehw,apresg,apresh,apx,avcnst,axhfe,axhfs,axhfsd,axlke,axlks,axlksd,axvfe,axvfs,axvfsd,bdh,bdhh,bdhv,bdot,bdoth,bdotv,dadhh,dadhv,dbdhh,dbdhv,dbdth,dbdtv,dhfe,dhfs,dhfsd,dvfe,dvfs,dvfsd,eact,ehfac,farad,hact,iact,iapxmx,iktmax,imchmx,imech,iopg,iopt,ipbtmx,ipch,ipchmx,ipcv,ipcvmx,ixrn1,ixrn2,jpfcmx,jpress,jptffl,jsol,jtemp,narxmx,narxt,narxth,nbasp,nbaspd,nbt,nbtd,nbtmax,ncmpr,ndrsr,ndrsrd,nmut,nmutmx,nopgmx,noptmx,noutpt,nptkmx,nptmax,nrct,nrctmx,nrk,nslt,nsltmx,nst,nstmax,ntpr,ntprmx,ntprt,nttkmx,nttyo,nweope,nwndpc,nxt,nxtmax,pmu,presg,presh,press,pressb,pressd,pslamn,ptk,rcnstv,rconst,rk,rkb,rtcnst,tempc,tempcb,tempcd,tempcu,tempk,time1,trkb,ttk,uphase,uspec,wfac,xhfe,xhfs,xhfsd,xi1,xlke,xlks,xlksd,xvfe,xvfs,xvfsd)
    !! This subroutine changes the temperature and pressure and
    !! recomputes all temperature-pressure dependent data. Here tempcd
    !! and pressd are the last temperature and pressure, respectively,
    !! at which the data were evaluated.
    !! This subroutine is called by:
    !!   EQ6/path.f
    !!   EQ6/eqshel.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: iapxmx
    integer :: iktmax
    integer :: imchmx
    integer :: ipbtmx
    integer :: ipchmx
    integer :: ipcvmx
    integer :: jpfcmx
    integer :: narxmx
    integer :: nbtmax
    integer :: nmutmx
    integer :: nopgmx
    integer :: noptmx
    integer :: nptkmx
    integer :: nptmax
    integer :: nrctmx
    integer :: nsltmx
    integer :: nstmax
    integer :: ntprmx
    integer :: nttkmx
    integer :: nxtmax

    integer :: noutpt
    integer :: nttyo

    integer :: iact(imchmx,2,nrctmx)
    integer :: imech(2,nrctmx)
    integer :: iopg(nopgmx)
    integer :: iopt(noptmx)
    integer :: jsol(nxtmax)
    integer :: narxt(ntprmx)
    integer :: narxth(2)
    integer :: nbasp(nbtmax)
    integer :: nbaspd(nbtmax)
    integer :: ncmpr(2,nptmax)
    integer :: ndrsr(2,nstmax)
    integer :: ndrsrd(2,nstmax)
    integer :: nrk(2,nrctmx)

    integer :: ipch
    integer :: ipcv
    integer :: ixrn1
    integer :: ixrn2
    integer :: jpress
    integer :: jptffl
    integer :: jtemp
    integer :: nbt
    integer :: nbtd
    integer :: nmut
    integer :: nrct
    integer :: nslt
    integer :: nst
    integer :: ntpr
    integer :: ntprt
    integer :: nweope
    integer :: nwndpc
    integer :: nxt

    character(len=48) :: uspec(nstmax)
    character(len=24) :: uphase(nptmax)

    real(kind=8) :: tempcu(ntprmx)

    real(kind=8) :: aadh(narxmx,ntprmx)
    real(kind=8) :: aadhh(narxmx,ntprmx)
    real(kind=8) :: aadhv(narxmx,ntprmx)
    real(kind=8) :: aaphi(narxmx,ntprmx)
    real(kind=8) :: abdh(narxmx,ntprmx)
    real(kind=8) :: abdhh(narxmx,ntprmx)
    real(kind=8) :: abdhv(narxmx,ntprmx)
    real(kind=8) :: abdot(narxmx,ntprmx)
    real(kind=8) :: abdoth(narxmx,ntprmx)
    real(kind=8) :: abdotv(narxmx,ntprmx)
    real(kind=8) :: adadhh(narxmx,ntprmx,ipchmx)
    real(kind=8) :: adadhv(narxmx,ntprmx,ipcvmx)
    real(kind=8) :: adbdhh(narxmx,ntprmx,ipchmx)
    real(kind=8) :: adbdhv(narxmx,ntprmx,ipcvmx)
    real(kind=8) :: adbdth(narxmx,ntprmx,ipchmx)
    real(kind=8) :: adbdtv(narxmx,ntprmx,ipcvmx)
    real(kind=8) :: adhfe(narxmx,ntprmx,ipchmx)
    real(kind=8) :: adhfs(narxmx,ntprmx,ipchmx,nstmax)
    real(kind=8) :: adhfsd(narxmx,ntprmx,ipchmx,nstmax)
    real(kind=8) :: advfe(narxmx,ntprmx,ipcvmx)
    real(kind=8) :: advfs(narxmx,ntprmx,ipcvmx,nstmax)
    real(kind=8) :: advfsd(narxmx,ntprmx,ipcvmx,nstmax)
    real(kind=8) :: aprehw(narxmx,ntprmx)
    real(kind=8) :: apresg(narxmx,ntprmx)
    real(kind=8) :: apresh(5,2)
    real(kind=8) :: axhfe(narxmx,ntprmx)
    real(kind=8) :: axhfs(narxmx,ntprmx,nstmax)
    real(kind=8) :: axhfsd(narxmx,ntprmx,nstmax)
    real(kind=8) :: axlke(narxmx,ntprmx)
    real(kind=8) :: axlks(narxmx,ntprmx,nstmax)
    real(kind=8) :: axlksd(narxmx,ntprmx,nstmax)
    real(kind=8) :: axvfe(narxmx,ntprmx)
    real(kind=8) :: axvfs(narxmx,ntprmx,nstmax)
    real(kind=8) :: axvfsd(narxmx,ntprmx,nstmax)

    real(kind=8) :: dadhh(ipchmx)
    real(kind=8) :: dadhv(ipcvmx)
    real(kind=8) :: dbdhh(ipchmx)
    real(kind=8) :: dbdhv(ipcvmx)
    real(kind=8) :: dbdth(ipchmx)
    real(kind=8) :: dbdtv(ipcvmx)
    real(kind=8) :: dhfe(ipchmx)
    real(kind=8) :: dhfs(ipchmx,nstmax)
    real(kind=8) :: dhfsd(ipchmx,nstmax)
    real(kind=8) :: dvfe(ipcvmx)
    real(kind=8) :: dvfs(ipcvmx,nstmax)
    real(kind=8) :: dvfsd(ipcvmx,nstmax)
    real(kind=8) :: xhfs(nstmax)
    real(kind=8) :: xhfsd(nstmax)
    real(kind=8) :: xlks(nstmax)
    real(kind=8) :: xlksd(nstmax)
    real(kind=8) :: xvfs(nstmax)
    real(kind=8) :: xvfsd(nstmax)

    real(kind=8) :: adh
    real(kind=8) :: adhh
    real(kind=8) :: adhv
    real(kind=8) :: aphi
    real(kind=8) :: bdh
    real(kind=8) :: bdhh
    real(kind=8) :: bdhv
    real(kind=8) :: bdot
    real(kind=8) :: bdoth
    real(kind=8) :: bdotv
    real(kind=8) :: prehw
    real(kind=8) :: presg
    real(kind=8) :: presh
    real(kind=8) :: press
    real(kind=8) :: xhfe
    real(kind=8) :: xlke
    real(kind=8) :: xvfe

    real(kind=8) :: amu(jpfcmx,nmutmx)
    real(kind=8) :: aslm(jpfcmx,0:ipbtmx,nsltmx)
    real(kind=8) :: pmu(nmutmx)
    real(kind=8) :: pslamn(0:ipbtmx,nsltmx)

    real(kind=8) :: apx(iapxmx,nxtmax)
    real(kind=8) :: wfac(iktmax,nxtmax)

    real(kind=8) :: eact(imchmx,2,nrctmx)
    real(kind=8) :: hact(imchmx,2,nrctmx)
    real(kind=8) :: ptk(nptkmx)
    real(kind=8) :: rk(imchmx,2,nrctmx)
    real(kind=8) :: rkb(imchmx,2,nrctmx)
    real(kind=8) :: trkb(imchmx,2,nrctmx)
    real(kind=8) :: ttk(nttkmx)

    real(kind=8) :: afcnst
    real(kind=8) :: al10
    real(kind=8) :: avcnst
    real(kind=8) :: ehfac
    real(kind=8) :: farad
    real(kind=8) :: pressb
    real(kind=8) :: pressd
    real(kind=8) :: rcnstv
    real(kind=8) :: rconst
    real(kind=8) :: rtcnst
    real(kind=8) :: tempc
    real(kind=8) :: tempcb
    real(kind=8) :: tempcd
    real(kind=8) :: tempk
    real(kind=8) :: time1
    real(kind=8) :: xi1

    ! Local variable declarations.
    integer :: nwpclm

    logical :: qnewp
    logical :: qnewt

    real(kind=8) :: dp
    real(kind=8) :: dt
    real(kind=8) :: pxl
    real(kind=8) :: pxu
    real(kind=8) :: toldt
    real(kind=8) :: toldp

    ! nwpclm  = maximum number of repeated warnings regarding problems
    !            associated with pressure corrections
    ! nwndpc = number of warnings of no data to support pressure
    !            corrections
    ! nweope = number of warnings of excursions outside the recommended
    !            pressure envelope
    data nwpclm /5/

    ! toldt = tolerance for recalculating temperature-dependent data (C)
    ! toldp = tolerance for recalculating pressure-dependent data (bars)
    data toldt / 1.e-4 /, toldp /1.e-4 /

    ! Compute the new temperature.
    call gtemp(afcnst,al10,iopt,jtemp,noptmx,noutpt,nttkmx,nttyo,rconst,rtcnst,tempc,tempcb,tempk,time1,ttk,xi1)

    dt = tempc - tempcd
    qnewt = .false.

    if (abs(dt) .gt. toldt) then
        ! Determine the corresponding temperature range flag.
        call gntpr(ntpr,ntprmx,ntprt,tempc,tempcu)

        qnewt = .true.

        ! Compute thermodynamic data for the current temperature.
        call evdata(aadh,aadhh,aadhv,aaphi,abdh,abdhh,abdhv,abdot,abdoth,abdotv,adadhh,adadhv,adbdhh,adbdhv,adbdth,adbdtv,adh,adhfe,adhh,adhv,adhfs,adhfsd,advfe,advfs,advfsd,afcnst,al10,amu,aslm,aphi,aprehw,apresg,apresh,apx,avcnst,axhfe,axhfs,axhfsd,axlke,axlks,axlksd,axvfe,axvfs,axvfsd,bdh,bdhh,bdhv,bdot,bdoth,bdotv,dadhh,dadhv,dbdhh,dbdhv,dbdth,dbdtv,dhfe,dhfs,dhfsd,dvfe,dvfs,dvfsd,ehfac,farad,iapxmx,iktmax,iopg,iopt,ipbtmx,ipch,ipchmx,ipcv,ipcvmx,ixrn1,ixrn2,jpfcmx,jptffl,jsol,narxmx,narxt,narxth,ncmpr,nmut,nmutmx,nopgmx,noptmx,noutpt,nptmax,nslt,nsltmx,nst,nstmax,ntpr,ntprmx,nttyo,nxt,nxtmax,pmu,prehw,presg,presh,press,pslamn,rconst,rcnstv,rtcnst,tempc,tempk,uphase,uspec,wfac,xhfe,xhfs,xhfsd,xlke,xlks,xlksd,xvfe,xvfs,xvfsd)

        ! Compute kinetic data for the current temperature.
        call evratc(eact,hact,iact,imchmx,imech,nrct,nrctmx,nrk,rk,rkb,rtcnst,tempk,trkb)
    end if

    ! Compute the new pressure.
    call gpress(iopt,jpress,noptmx,noutpt,nptkmx,nttyo,presg,presh,press,pressb,time1,ptk,xi1)

    ! Skip pressure corrections if the pressure is tracking on the
    ! data file reference pressure curve.
    if (jpress .le. 0) then
        go to 990
    end if

    ! Skip pressure corrections if the temperature and pressure
    ! haven't changed significantly.
    dp = press - pressd
    qnewp = abs(dp) .gt. toldp

    if (.not.qnewt .and. .not.qnewp) then
        go to 990
    end if

    ! Skip pressure corrections if the pressure happens to match the
    ! data file reference pressure curve value.
    dp = press - presg

    if (abs(dp) .le. toldp) then
        go to 990
    end if

    ! Make pressure corrections, if the requisite data are available.
    if (ipcv .lt. 0) then
        ! There are no data to support needed pressure corrections.
        if (nwndpc .le. nwpclm) then
            write (noutpt,1000) press,presg,dp
            write (nttyo,1000) press,presg,dp
1000 format(/' * Warning - (EQ6/tpadv) The supporting data file',/7x,'contains no data to support making thermodynamic',/7x,'pressure corrections. No such corrections will be made.',/7x,'The current pressure is ',1pg12.5,' bars, the standard',/7x,'grid pressure is ',g12.5,' bars, and the pressure',/7x,'difference is ',g12.5,' bars.')

            if (nwndpc .eq. nwpclm) then
                write (noutpt,1010)
                write (nttyo,1010)
1010 format(/' This warning will not be repeated during the',/7x,'rest of this run.')
            end if

            nwndpc = nwndpc + 1
        end if
    else
        if (abs(dp) .gt. prehw) then
            ! The pressure is outside the recommended envelope.
            pxu = presg + prehw
            pxl = presg - prehw

            if (pxl .le. 0.) then
                pxl = 0.
            end if

            write (noutpt,1020) press,pxl,pxu,tempc
            write (nttyo,1020) press,pxl,pxu,tempc
1020 format(/' * Warning - (EQ6/tpadv) The current pressure',/7x,'of ',1pg12.5,' bars is outside the recommended',/7x,'pressure envelope of ',g12.5,' to ',g12.5,' bars',/7x,'at ',0pf6.2,' C.')

            if (nweope .eq. nwpclm) then
                write (noutpt,1010)
                write (nttyo,1010)
            end if

            nweope = nweope + 1
        end if

        ! Make pressure corrections to the thermodynamic data.
        call pcorrx(avcnst,dhfs,dvfs,ipch,ipchmx,ipcv,ipcvmx,nbasp,nbt,nbtmax,ndrsr,nst,nstmax,presg,press,xhfs,xlks,xvfs)

        ! Calling sequence substitutions:
        !   dhfsd for dhfs
        !   dvfsd for dvfs
        !   nbaspd for nbasp
        !   nbtd for nbt
        !   ndrsrd for ndrsr
        !   xhfsd for xhfs
        !   xlksd for xlks
        !   xvfsd for xvfs
        call pcorrx(avcnst,dhfsd,dvfsd,ipch,ipchmx,ipcv,ipcvmx,nbaspd,nbtd,nbtmax,ndrsrd,nst,nstmax,presg,press,xhfsd,xlksd,xvfsd)

        ! Make pressure corrections to the kinetic data.
        !   Presently there are no such corrections.
    end if

990 continue
    tempcd = tempc
    pressd = press
end subroutine tpadv
