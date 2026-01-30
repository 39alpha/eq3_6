subroutine satchk(acflg,act,actlg,afcnst,affp,affs,apx,bpx,cdrs,eps100,iindx1,iodb,iopt,iapxmx,ibpxmx,iktmax,ixrn1,jflag,jpflag,jsflag,jsol,kmax,km1,kpsat,kpsst,kxt,nbasp,nbt,nbtmax,ncmpr,ndrs,ndrsmx,ndrsr,nodbmx,noptmx,noutpt,npchk,npt,nptmax,nstmax,nttyo,nxrn1,nxrn2,nxtmax,qxknph,sidrph,sidrsp,tolsat,uphase,uspec,wfac,xbar,xbarlg,xlks)
    !! This subroutine checks for newly saturated phases and
    !! supersaturated phases.
    !! This subroutine is called by:
    !!   EQ6/eqphas.f
    !! Principal input:
    !! Principal output:
    !!   affp   = array of affinities associated with phases
    !!   affs   = array of affinities associated with species
    !!   jpflag = status flag array for phases
    !!   kpsat  = number of newly saturated phases
    !!   kpsst  = number of supersaturated phases
    !!   sidrph = array of saturation indices (log Q/K) associated
    !!              with phases
    !!   sidrsp = array of saturation indices (log Q/K) associated
    !!              with species
    implicit none

    ! Calling sequence variable declarations.
    integer :: iapxmx
    integer :: ibpxmx
    integer :: iktmax
    integer :: kmax
    integer :: nbtmax
    integer :: ndrsmx
    integer :: nodbmx
    integer :: noptmx
    integer :: nptmax
    integer :: nstmax
    integer :: nxtmax

    integer :: noutpt
    integer :: nttyo

    integer :: iindx1(kmax)
    integer :: iodb(nodbmx)
    integer :: iopt(noptmx)
    integer :: jflag(nstmax)
    integer :: jpflag(nptmax)
    integer :: jsflag(nstmax)
    integer :: jsol(nxtmax)
    integer :: nbasp(nbtmax)
    integer :: ncmpr(2,nptmax)
    integer :: ndrs(ndrsmx)
    integer :: ndrsr(2,nstmax)
    integer :: npchk(nptmax)

    integer :: ixrn1
    integer :: km1
    integer :: kpsat
    integer :: kpsst
    integer :: kxt
    integer :: nbt
    integer :: npt
    integer :: nxrn1
    integer :: nxrn2

    logical :: qxknph(nptmax)

    character(len=48) :: uspec(nstmax)
    character(len=24) :: uphase(nptmax)

    real(kind=8) :: acflg(nstmax)
    real(kind=8) :: act(nstmax)
    real(kind=8) :: actlg(nstmax)
    real(kind=8) :: affp(nptmax)
    real(kind=8) :: affs(nstmax)
    real(kind=8) :: apx(iapxmx,nxtmax)
    real(kind=8) :: bpx(ibpxmx,nxtmax)
    real(kind=8) :: cdrs(ndrsmx)
    real(kind=8) :: sidrph(nptmax)
    real(kind=8) :: sidrsp(nstmax)
    real(kind=8) :: wfac(iktmax,nxtmax)
    real(kind=8) :: xbar(nstmax)
    real(kind=8) :: xbarlg(nstmax)
    real(kind=8) :: xlks(nstmax)

    real(kind=8) :: afcnst
    real(kind=8) :: eps100
    real(kind=8) :: tolsat

    ! Local variable declarations.
    integer :: ier
    integer :: j2
    integer :: kcol
    integer :: nb
    integer :: np
    integer :: nr1
    integer :: nr2
    integer :: ns
    integer :: nt

    integer :: ilnobl

    logical :: qtest

    character(len=24) :: uaqsln
    character(len=24) :: ugas

    real(kind=8) :: af
    real(kind=8) :: si

    data uaqsln /'Aqueous solution        '/     ugas   /'Gas                     '/

    ! XX   This double-checking of the jsflag and jpflag arrays shouldn't be
    ! XX   necessary.
    !      Double-check flags on status of non-aqueous phases which are
    !      in the matrix.
    do np = 1,npt
        if (uphase(np)(1:24) .ne. uaqsln(1:24)) then
            if (jpflag(np) .le. 0) then
                jpflag(np) = 0
                nr1 = ncmpr(1,np)
                nr2 = ncmpr(2,np)

                do ns = nr1,nr2
                    if (jsflag(ns) .le. 0) then
                        jsflag(ns) = 0

                        do kcol = km1,kxt
                            if (ns .eq. iindx1(kcol)) then
                                jsflag(ns) = -1
                                jpflag(np) = -1
                                go to 30
                            end if
                        end do
                    end if

30 continue
                end do
            end if
        end if
    end do

    kpsat = 0
    kpsst = 0

    ! Compute affinities.
    do np = 1,npt
        nr1 = ncmpr(1,np)
        nr2 = ncmpr(2,np)
        nt = nr2 - nr1 + 1

        if (uphase(np)(1:24) .eq. ugas(1:24)) then
            ! Case of a gas phase.
            affp(np) = -9999999.
            sidrph(np) = -9999999.

            do ns = nr1,nr2
                affs(ns) = -9999999.
                sidrsp(ns) = -9999999.
            end do
        else if (nt .eq. 1) then
            ! Case of a pure phase.
            ns = nr1

            if (jpflag(np) .le. 1) then
                call afcalc(actlg,af,afcnst,cdrs,jflag,jsflag,ndrs,ndrsmx,ndrsr,ns,nstmax,si,xlks)
                affs(ns) = af
                affp(np) = af
            else
                af = -9999999.
                si = -9999999.
            end if

            affs(ns) = af
            affp(np) = af
            sidrsp(ns) = si
            sidrph(np) = si
        else if (uphase(np)(1:24) .eq. uaqsln(1:24)) then
            ! Case of an aqueous solution.
            if (jpflag(np) .le. 1) then
                affp(np) = 0.
                sidrph(np) = 0.

                do nb = 1,nbt
                    ns = nbasp(nb)

                    if (jsflag(ns) .le. 1) then
                        call afcalc(actlg,af,afcnst,cdrs,jflag,jsflag,ndrs,ndrsmx,ndrsr,ns,nstmax,si,xlks)
                        affs(ns) = af
                        sidrsp(ns) = si
                        affp(np) = affp(np) + xbar(ns)*affs(ns)
                        sidrph(np) = sidrph(np) + xbar(ns)*sidrsp(ns)
                    else
                        affs(ns) = -9999999.
                        sidrsp(ns) = -9999999.
                    end if
                end do
            else
                affp(np) = -9999999.
                sidrph(np) = -9999999.
            end if
        else if (iopt(4) .ge. 1) then
            ! If not doing non-aqueous solutions, skip.
            if (jpflag(np) .eq. -1) then
                ! Case of non-aqueous solutions in the ES.
                qxknph(np) = .true.
                affp(np) = 0.
                sidrph(np) = 0.
                qtest = .false.

                do ns = nr1,nr2
                    if (jsflag(ns) .le. 1) then
                        call afcalc(actlg,af,afcnst,cdrs,jflag,jsflag,ndrs,ndrsmx,ndrsr,ns,nstmax,si,xlks)
                        affs(ns) = af
                        sidrsp(ns) = si
                        affp(np) = affp(np) + xbar(ns)*affs(ns)
                        sidrph(np) = sidrph(np) + xbar(ns)*sidrsp(ns)
                        qtest = .true.
                    else
                        affs(ns) = -9999999.
                        sidrsp(ns) = -9999999.
                    end if
                end do

                if (.not.qtest) then
                    affp(np) = -9999999.
                    sidrph(np) = -9999999.
                end if
            else if (jpflag(np) .le. 1) then
                ! Case of non-aqueous solutions not in the ES.
                call hpsat(acflg,act,actlg,afcnst,affp,affs,apx,bpx,cdrs,eps100,iapxmx,ibpxmx,ier,iktmax,ixrn1,jflag,jpflag,jsflag,jsol,ncmpr,ndrs,ndrsmx,ndrsr,noutpt,np,nptmax,nstmax,nttyo,nxrn1,nxrn2,nxtmax,sidrsp,sidrph,uphase,uspec,wfac,xbar,xbarlg,xlks)

                if (ier .le. 0) then
                    qxknph(np) = .true.
                else
                    qxknph(np) = .false.
                    j2 = ilnobl(uphase(np))
                    write (noutpt,2005) uphase(np)(1:j2)
                    write (nttyo,2005) uphase(np)(1:j2)
2005 format(/' * Warning - (EQ6/satchk) A hypothetical',' affinity calculation',/7x,'failed for solid',' solution ',a,'.')

                    affp(np) = -9999999.
                    sidrph(np) = -9999999.

                    do ns = nr1,nr2
                        affs(ns) = -9999999.
                        sidrsp(ns) = -9999999.
                    end do
                end if
            end if
        end if

        ! Set the jpflag array to mark supersaturated and newly saturated
        ! phases.
        if (jpflag(np).le.0 .and. jpflag(np).ne.-1 .and.    npchk(np) .ne. 1) then
            if (abs(affp(np)) .le. tolsat) then
                jpflag(np) = -10
                kpsat = kpsat + 1
            else if (affp(np) .gt. tolsat) then
                jpflag(np) = -2
                kpsst = kpsst + 1
            end if
        end if

        ! Set the jsflag array to mark the corresponding species.
        if (uphase(np)(1:24) .eq. uaqsln(1:24)) then
            do nb = 1,nbt
                ns = nbasp(nb)

                if (jsflag(ns).le.0 .and. jsflag(ns).ne.-1) then
                    if (abs(affs(ns)) .le. tolsat) then
                        jsflag(ns) = -10
                    else if (affs(ns) .gt. tolsat) then
                        jsflag(ns) = -2
                    end if
                end if
            end do
        else
            do ns = nr1,nr2
                if (jsflag(ns).le.0 .and. jsflag(ns).ne.-1) then
                    if (abs(affs(ns)) .le. tolsat) then
                        jsflag(ns) = -10
                    else if (affs(ns) .gt. tolsat) then
                        jsflag(ns) = -2
                    end if
                end if
            end do
        end if
    end do
end subroutine satchk
