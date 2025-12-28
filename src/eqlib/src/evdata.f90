subroutine evdata(aadh,aadhh,aadhv,aaphi,abdh,abdhh,abdhv,abdot,abdoth,abdotv,adadhh,adadhv,adbdhh,adbdhv,adbdth,adbdtv,adh,adhfe,adhh,adhv,adhfs,adhfsd,advfe,advfs,advfsd,afcnst,al10,amu,aslm,aphi,aprehw,apresg,apresh,apx,avcnst,axhfe,axhfs,axhfsd,axlke,axlks,axlksd,axvfe,axvfs,axvfsd,bdh,bdhh,bdhv,bdot,bdoth,bdotv,dadhh,dadhv,dbdhh,dbdhv,dbdth,dbdtv,dhfe,dhfs,dhfsd,dvfe,dvfs,dvfsd,ehfac,farad,iapxmx,iktmax,iopg,iopt,ipbtmx,ipch,ipchmx,ipcv,ipcvmx,ixrn1,ixrn2,jpfcmx,jptffl,jsol,narxmx,narxt,narxth,ncmpr,nmut,nmutmx,nopgmx,noptmx,noutpt,nptmax,nslt,nsltmx,nst,nstmax,ntpr,ntprmx,nttyo,nxt,nxtmax,pmu,prehw,presg,presh,press,pslamn,rconst,rcnstv,rtcnst,tempc,tempk,uphase,uspec,wfac,xhfe,xhfs,xhfsd,xlke,xlks,xlksd,xvfe,xvfs,xvfsd)
    !! This subroutine evaluates thermodynamic properties as a function
    !! of temperature and the standard grid pressure. Pressure
    !! corrections are made by EQLIB/pcorrx.f.
    !! This subroutine is called by:
    !!   EQ3NR/eq3nr.f
    !!   EQ6/eq6.f
    !!   EQ6/tpadv.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: iapxmx
    integer :: iktmax
    integer :: ipbtmx
    integer :: ipchmx
    integer :: ipcvmx
    integer :: jpfcmx
    integer :: narxmx
    integer :: nmutmx
    integer :: nopgmx
    integer :: noptmx
    integer :: nptmax
    integer :: nsltmx
    integer :: nstmax
    integer :: ntprmx
    integer :: nxtmax

    integer :: noutpt
    integer :: nttyo

    integer :: iopg(nopgmx)
    integer :: iopt(noptmx)
    integer :: jsol(nxtmax)
    integer :: narxt(ntprmx)
    integer :: narxth(2)
    integer :: ncmpr(2,nptmax)

    integer :: ipch
    integer :: ipcv
    integer :: ixrn1
    integer :: ixrn2
    integer :: jptffl
    integer :: nmut
    integer :: nslt
    integer :: nst
    integer :: ntpr
    integer :: nxt

    character(len=48) :: uspec(nstmax)
    character(len=24) :: uphase(nptmax)

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

    real(kind=8) :: afcnst
    real(kind=8) :: al10
    real(kind=8) :: avcnst
    real(kind=8) :: ehfac
    real(kind=8) :: farad
    real(kind=8) :: rconst
    real(kind=8) :: rcnstv
    real(kind=8) :: rtcnst
    real(kind=8) :: tempc
    real(kind=8) :: tempk

    ! Local variable declarations.
    integer :: ipc
    integer :: ntprh

    real(kind=8) :: prop

    ! Compute temperature dependent constants.
    rtcnst = 0.001*rconst*tempk
    afcnst = al10*rtcnst
    avcnst = al10*rcnstv*tempk
    ehfac = (al10*rconst*tempk)/farad

    ! Compute the standard grid pressure.
    ! Calling sequence substitutions:
    !   apresg for arr
    !   presg for prop
    call evdat2(apresg,narxmx,narxt,ntpr,ntprmx,presg,tempc)

    ! Compute the 1.013-bar/steam-saturation curve pressure at the
    ! initial temperature.
    if (tempc .le. 100.) then
        ntprh = 1
    else
        ntprh = 2
    end if

    ! Calling sequence substitutions:
    !   apresh for arr
    !   5 for narxmx
    !   narxth for narxt
    !   ntprh for ntpr
    !   2 for ntprmx
    !   presh for prop
    call evdat2(apresh,5,narxth,ntprh,2,presh,tempc)

    if (ipcv .ge. 0) then
        ! Compute the half-width of the recommended pressure envelope.
        ! Calling sequence substitutions:
        !   aprehw for arr
        !   prehw for prop
        call evdat2(aprehw,narxmx,narxt,ntpr,ntprmx,prehw,tempc)
    end if

    if (iopg(1) .le. 0) then
        ! Davies equation or B-dot equation.
        ! Compute the Debye-Huckel A(gamma,10) parameter.
        ! Calling sequence substitutions:
        !   aadh for arr
        !   adh for prop
        call evdat2(aadh,narxmx,narxt,ntpr,ntprmx,adh,tempc)

        aphi = adh*al10/3.
    end if

    if (iopg(1) .eq. 1) then
        ! Pitzer's equations.
        ! Compute the Debye-Huckel A(phi) parameter.
        ! Calling sequence substitutions:
        !   aaphi for arr
        !   aphi for prop
        call evdat2(aaphi,narxmx,narxt,ntpr,ntprmx,aphi,tempc)
        adh= 3.*aphi/al10
    end if

    if (iopg(1) .eq. 2) then
        ! HC + DH equations.
        ! Compute the Debye-Huckel A(gamma,10) parameter.
        ! Calling sequence substitutions:
        !   aadh for arr
        !   adh for prop
        call evdat2(aadh,narxmx,narxt,ntpr,ntprmx,adh,tempc)

        aphi = adh*al10/3.
    end if

    if (ipch .ge. 0) then
        ! Compute the Debye-Huckel A(H) parameter.
        ! Calling sequence substitutions:
        !   aadhh for arr
        !   adhh for prop
        call evdat2(aadhh,narxmx,narxt,ntpr,ntprmx,adhh,tempc)

        do ipc = 1,ipch
            ! Compute the pressure derivatives of the Debye-Huckel A(H)
            ! parameter.
            ! Calling sequence substitutions:
            !   adadhh for arr
            !   ipc for k
            !   ipchmx for nmax
            call evdat3(adadhh,ipc,ipchmx,narxmx,narxt,ntpr,ntprmx,prop,tempc)
            dadhh(ipc) = prop
        end do
    end if

    if (ipcv .ge. 0) then
        ! Compute the Debye-Huckel A(V) parameter.
        ! Calling sequence substitutions:
        !   aadhv for arr
        !   adhv for prop
        call evdat2(aadhv,narxmx,narxt,ntpr,ntprmx,adhv,tempc)

        do ipc = 1,ipcv
            ! Compute the pressure derivatives of the Debye-Huckel A(V)
            ! parameter.
            ! Calling sequence substitutions:
            !   adadhv for arr
            !   ipc for k
            !   ipcvmx for nmax
            call evdat3(adadhv,ipc,ipcvmx,narxmx,narxt,ntpr,ntprmx,prop,tempc)
            dadhv(ipc) = prop
        end do
    end if

    if (iopg(1).eq.0 .or. iopg(1).eq.2) then
        ! B-dot equation or HC + DH equations.
        ! Compute the Debye-Huckel B(gamma) parameter.
        ! Calling sequence substitutions:
        !   abdh for arr
        !   bdh for prop
        call evdat2(abdh,narxmx,narxt,ntpr,ntprmx,bdh,tempc)

        if (ipch .ge. 0) then
            ! Compute the Debye-Huckel B(H) parameter.
            ! Calling sequence substitutions:
            !   abdhh for arr
            !   bdhh for prop
            call evdat2(abdhh,narxmx,narxt,ntpr,ntprmx,bdhh,tempc)

            do ipc = 1,ipch
                ! Compute the pressure derivatives of the Debye-Huckel B(H)
                ! parameter.
                ! Calling sequence substitutions:
                !   adbdhh for arr
                !   ipc for k
                !   ipchmx for nmax
                call evdat3(adbdhh,ipc,ipchmx,narxmx,narxt,ntpr,ntprmx,prop,tempc)
                dbdhh(ipc) = prop
            end do
        end if

        if (ipcv .ge. 0) then
            ! Compute the Debye-Huckel B(V) parameter.
            ! Calling sequence substitutions:
            !   abdhv for arr
            !   bdhv for prop
            call evdat2(abdhv,narxmx,narxt,ntpr,ntprmx,bdhv,tempc)

            do ipc = 1,ipcv
                ! Compute the pressure derivatives of the Debye-Huckel B(V)
                ! parameter.
                ! Calling sequence substitutions:
                !   adbdhv for arr
                !   ipc for k
                !   ipcvmx for nmax
                call evdat3(adbdhv,ipc,ipcvmx,narxmx,narxt,ntpr,ntprmx,prop,tempc)
                dbdhv(ipc) = prop
            end do
        end if
    end if

    if (iopg(1) .eq. 0) then
        ! B-dot equation.
        ! Compute the Helgeson (1969) B-dot parameter.
        ! Calling sequence substitutions:
        !   abdot for arr
        !   bdot for prop
        call evdat2(abdot,narxmx,narxt,ntpr,ntprmx,bdot,tempc)

        if (ipch .ge. 0) then
            ! Compute the Debye-Huckel B-dot(H) parameter.
            ! Calling sequence substitutions:
            !   abdoth for arr
            !   bdoth for prop
            call evdat2(abdoth,narxmx,narxt,ntpr,ntprmx,bdoth,tempc)

            do ipc = 1,ipch
                ! Compute the pressure derivatives of the Debye-Huckel
                ! B-dot(H) parameter.
                ! Calling sequence substitutions:
                !   adbdth for arr
                !   ipc for k
                !   ipchmx for nmax
                call evdat3(adbdth,ipc,ipchmx,narxmx,narxt,ntpr,ntprmx,prop,tempc)
                dbdth(ipc) = prop
            end do
        end if

        if (ipcv .ge. 0) then
            ! Compute the Debye-Huckel B-dot(V) parameter.
            ! Calling sequence substitutions:
            !   abdotv for arr
            !   bdotv for prop
            call evdat2(abdotv,narxmx,narxt,ntpr,ntprmx,bdotv,tempc)

            do ipc = 1,ipcv
                ! Compute the pressure derivatives of the Debye-Huckel
                ! B-dot(V) parameter.
                ! Calling sequence substitutions:
                !   adbdtv for arr
                !   ipc for k
                !   ipcvmx for nmax
                call evdat3(adbdtv,ipc,ipcvmx,narxmx,narxt,ntpr,ntprmx,prop,tempc)
                dbdtv(ipc) = prop
            end do
        end if
    end if

    ! Compute the log K for the 'Eh' reaction.
    if (axlke(1,ntpr) .lt. 9999999.) then
        ! Calling sequence substitutions:
        !   axlke for arr
        !   xlke for prop
        call evdat2(axlke,narxmx,narxt,ntpr,ntprmx,xlke,tempc)
    else
        xlke = 9999999.
    end if

    if (ipch .ge. 0) then
        ! Compute the enthalpy of reaction for the 'Eh' reaction.
        if (axhfe(1,ntpr) .lt. 9999999.) then
            ! Calling sequence substitutions:
            !   axhfe for arr
            !   xhfe for prop
            call evdat2(axhfe,narxmx,narxt,ntpr,ntprmx,xhfe,tempc)
        else
            xhfe = 9999999.
        end if

        do ipc = 1,ipch
            ! Compute the pressure derivatives of the enthalpy of
            ! reaction for the 'Eh' reaction.
            ! Calling sequence substitutions:
            !   adhfe for arr
            !   ipc for k
            !   ipchmx for nmax
            call evdat3(adhfe,ipc,ipchmx,narxmx,narxt,ntpr,ntprmx,prop,tempc)
            dhfe(ipc) = prop
        end do
    end if

    if (ipcv .ge. 0) then
        ! Compute the volume of reaction for the 'Eh' reaction.
        if (axvfe(1,ntpr) .lt. 9999999.) then
            ! Calling sequence substitutions:
            !   axvfe for arr
            !   xvfe for prop
            call evdat2(axvfe,narxmx,narxt,ntpr,ntprmx,xvfe,tempc)
        else
            xvfe = 9999999.
        end if

        do ipc = 1,ipcv
            ! Compute the pressure derivatives of the volume of
            ! reaction for the 'Eh' reaction.
            ! Calling sequence substitutions:
            !   advfe for arr
            !   ipc for k
            !   ipcvmx for nmax
            call evdat3(advfe,ipc,ipcvmx,narxmx,narxt,ntpr,ntprmx,prop,tempc)
            dvfe(ipc) = prop
        end do
    end if

    ! Compute the log K values for all reactions as they are
    ! currently written.
    call evdatr(adhfs,advfs,axhfs,axlks,axvfs,dhfs,dvfs,ipch,ipchmx,ipcv,ipcvmx,narxmx,narxt,nst,nstmax,ntpr,ntprmx,tempc,xhfs,xlks,xvfs)

    ! Compute the log K values for all reactions as they were
    ! written on the data file.
    ! Calling sequence substitutions:
    !   adhfsd for adhfs
    !   advfsd for advfs
    !   axhfsd for axhfs
    !   axlksd for axlks
    !   axvfsd for avhfs
    !   dhfsd for dhfs
    !   dvfsd for dvfs
    !   xhfsd for xhfs
    !   xlksd for xlks
    !   xvfsd for xvfs
    call evdatr(adhfsd,advfsd,axhfsd,axlksd,axvfsd,dhfsd,dvfsd,ipch,ipchmx,ipcv,ipcvmx,narxmx,narxt,nst,nstmax,ntpr,ntprmx,tempc,xhfsd,xlksd,xvfsd)

    ! Compute remaining temperature-dependent aqueous species activity
    ! coefficient parameters.
    if (iopg(1) .eq. 1) then
        ! Computer the S-lambda(n) and mu coefficients for Pitzer's
        ! equations.
        call evptzc(amu,aslm,ipbtmx,jpfcmx,jptffl,nmut,nmutmx,noutpt,nttyo,nslt,nsltmx,pmu,pslamn,tempc)
    end if

    ! Compute temperature dependent solid solution activity
    ! coefficient parameters.
    if (iopt(4) .ge. 1) then
        call wterm(apx,iapxmx,iktmax,ixrn1,ixrn2,jsol,ncmpr,noutpt,nptmax,nstmax,nttyo,nxt,nxtmax,press,tempk,uphase,uspec,wfac)
    end if
end subroutine evdata