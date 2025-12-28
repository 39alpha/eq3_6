subroutine shftph(emop,emop0,emos,emos0,fdpe0,fdpem1,fdse0,fdsem1,iemop,iemos,iern1,iern2,ietmax,iindx1,imrn1,imrn2,ipndx1,ixrn1,ixrn2,jern1,jetmax,jgext,jpflag,jsflag,kbt,kmax,km1,kmt,kx1,kxt,loph,losp,moph,mosp,mprph,mprsp,mrgexs,nbtmax,ncmpe,ncmpr,netmax,ngext,nordmx,noutpt,np,npet,npetmx,nptmax,nsetmx,nstmax,nttyo,qshftd,qtotsh,uphase,xbar,xbarlg,zklgmn,zklogl,zvclg0,zvclg1,zvec0,zvec1)
    !! This subroutine shifts mass of a single phase from the equilibrium
    !! system (ES) to the physically removed subsystem (PRS), which is
    !! a system conceptually out of contact with the aqueous solution.
    !! Here np is the index of the desired phase.
    !! Presently, the following kinds of phases may be shifted:
    !!   Generic ion exchange phases
    !!   Pure minerals
    !!   Solid solutions
    !! A shift may be total or partial. If qtotsh = .true., the shift
    !! is total and the entire mineral phase is transferred. Otherwise,
    !! the shift is partial and the shift transfers only some of the
    !! mass. At least 10**zklgmn moles of mass will remain in the ES.
    !! If the mass is less than that to start with, no shift takes
    !! place. Otherwise, the ES mass is reduced by multiplying the
    !! existing mass by a factor of 10**-zklogl.
    !! This subroutine decrements the ES mass of the phase and
    !! correspondingly increments its PRS mass. However, it does
    !! not recalculate the component total masses for the ES. After
    !! this subroutine has been called, once or in a sequence, a call
    !! must be made to EQ6/escalc.f, which recalculates the component
    !! mass balance totals.
    !! If a total shift is made, it will subsequently be necessary to
    !! update the iemop indexing array and associated arrays used in
    !! tracking the masses of phases in the ES.
    !! If a shift is made, qshftd is returned as .true.; otherwise,
    !! it is returned as .false.
    !! The way the shift is made presumes that the calculation is
    !! currently at the "zero" point, the value of Xi from which the
    !! calculation will then step out.
    !! This subroutine is called by:
    !!   EQ6/dumpdp.f
    !!   EQ6/pshfta.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: ietmax
    integer :: jetmax
    integer :: kmax
    integer :: nbtmax
    integer :: netmax
    integer :: nordmx
    integer :: npetmx
    integer :: nptmax
    integer :: nsetmx
    integer :: nstmax

    integer :: noutpt
    integer :: nttyo

    integer :: iemop(npetmx)
    integer :: iemos(nsetmx)
    integer :: iindx1(kmax)
    integer :: ipndx1(kmax)
    integer :: jern1(jetmax,netmax)
    integer :: jgext(netmax)
    integer :: jpflag(nptmax)
    integer :: jsflag(nstmax)
    integer :: ncmpe(2,npetmx)
    integer :: ncmpr(2,nptmax)
    integer :: ngext(jetmax,netmax)

    integer :: iern1
    integer :: iern2
    integer :: imrn1
    integer :: imrn2
    integer :: ixrn1
    integer :: ixrn2
    integer :: kbt
    integer :: km1
    integer :: kmt
    integer :: kx1
    integer :: kxt
    integer :: np
    integer :: npet

    logical :: qshftd
    logical :: qtotsh

    character(len=24) :: uphase(nptmax)

    real(kind=8) :: emop(npetmx)
    real(kind=8) :: emop0(npetmx)
    real(kind=8) :: emos(nsetmx)
    real(kind=8) :: emos0(nsetmx)
    real(kind=8) :: fdpe0(nordmx,npetmx)
    real(kind=8) :: fdpem1(nordmx,npetmx)
    real(kind=8) :: fdse0(nordmx,nsetmx)
    real(kind=8) :: fdsem1(nordmx,nsetmx)
    real(kind=8) :: loph(nptmax)
    real(kind=8) :: losp(nstmax)
    real(kind=8) :: mrgexs(ietmax,jetmax,netmax)
    real(kind=8) :: moph(nptmax)
    real(kind=8) :: mosp(nstmax)
    real(kind=8) :: mprph(nptmax)
    real(kind=8) :: mprsp(nstmax)
    real(kind=8) :: xbar(nstmax)
    real(kind=8) :: xbarlg(nstmax)
    real(kind=8) :: zvclg0(kmax)
    real(kind=8) :: zvclg1(kmax)
    real(kind=8) :: zvec0(kmax)
    real(kind=8) :: zvec1(kmax)

    real(kind=8) :: zklgmn
    real(kind=8) :: zklogl

    ! Local variable declarations.
    integer :: ie
    integer :: iptype
    integer :: je
    integer :: j2
    integer :: kcol
    integer :: kstart
    integer :: kend
    integer :: n
    integer :: ne
    integer :: npe
    integer :: nr1
    integer :: nr2
    integer :: ns
    integer :: nse

    integer :: ilnobl

    real(kind=8) :: delmp
    real(kind=8) :: delms
    real(kind=8) :: lgexpl
    real(kind=8) :: lnew
    real(kind=8) :: lnew1
    real(kind=8) :: lold
    real(kind=8) :: lrx
    real(kind=8) :: mgexpl
    real(kind=8) :: mnew
    real(kind=8) :: mnew1
    real(kind=8) :: mold
    real(kind=8) :: mrx
    real(kind=8) :: zx

    real(kind=8) :: texp
    real(kind=8) :: tlg

    qshftd = .false.

    ! Set the minimum mass of an exchanger phase that must remain in
    ! the ES.
    lgexpl = -14.0
    mgexpl = texp(lgexpl)

    ! Get the type of phase.
    iptype = 0

    if (np.ge.imrn1 .and. np.le.imrn2) then
        ! Have a pure mineral.
        iptype = 1
    else if (np.ge.ixrn1 .and. np.le.ixrn2) then
        ! Have a solid solution.
        iptype = 2
    else if (np.ge.iern1 .and. np.le.iern2) then
        ! Have a generic ion exchange phase.
        iptype = 3
    end if

    if (iptype .le. 0) then
        j2 = ilnobl(uphase(np))
        write (noutpt,1000) np,uphase(np)(1:j2)
        write (nttyo,1000) np,uphase(np)(1:j2)
1000 format(/' * Error - (EQ6/shftph) Programming error trap: The',' phase',/7x,a," is not of a recognized type and can't have",' mass shifted',/7x,'from the ES to the PRS.')

        stop
    end if

    ! Save the current numbers of moles values.
    lold = loph(np)
    mold = moph(np)

    ! Calculate the shift increment.
    if (qtotsh) then
        ! Do a total shift.
        if (iptype .eq. 3) then
            !          Have a generic ion exchanger phase.
            ! XX       Presently a generic ion exchange phase can not be
            ! XX       totally shifted.
            if (mold .le. mgexpl) then
                write (noutpt,1010) uphase(np)(1:j2),mold
1010 format(3x,'Not shifted: ',a,/15x,'Mass= ',1pe12.5,' moles')

                go to 999
            end if

            lnew = lgexpl
            mnew = mgexpl
            delmp = mold - mgexpl
            jpflag(np) = 0
            j2 = ilnobl(uphase(np))
            write (noutpt,1020) uphase(np)(1:j2),mold,mnew
1020 format(7x,'Shifted: ',a,/    15x,'Old mass= ',1pe12.5,' moles, new mass= ',1pe12.5,' moles')
        else
            ! Have a pure mineral or a solid solution.
            lnew = -99999.
            mnew = 0.
            delmp = mold
            jpflag(np) = 0
            j2 = ilnobl(uphase(np))
            write (noutpt,1020) uphase(np)(1:j2),mold,mnew
        end if
    else
        ! Do a partial shift.
        if (lold .le. zklgmn) then
            ! Do not shift if the mass is already small.
            j2 = ilnobl(uphase(np))
            write (noutpt,1010) uphase(np)(1:j2),mold
            go to 999
        end if

        zx = lold - zklogl
        lnew = max(zx,zklgmn)
        mnew = texp(lnew)
        delmp = mold - mnew
        j2 = ilnobl(uphase(np))
        write (noutpt,1020) uphase(np)(1:j2),mold,mnew
    end if

    ! Complete the shift calculations.
    qshftd = .true.

    ! Reset the primary mass variables.
    loph(np) = lnew
    moph(np) = mnew

    ! Increment the PRS mass.
    mprph(np) = mprph(np) + delmp

    if (iptype.eq.1 .or. iptype.eq.2) then
        ! Have a pure mineral or a solid solution.
        nr1 = ncmpr(1,np)
        nr2 = ncmpr(2,np)

        do ns = nr1,nr2
            delms = delmp*xbar(ns)
            mprsp(ns) = mprsp(ns) + delms

            if (qtotsh) then
                ! Do a total shift.
                lnew1 = -99999.
                jsflag(ns) = 0
            else
                ! Do a partial shift.
                lnew1 = lnew + xbarlg(ns)
            end if

            mnew1 = texp(lnew1)
            losp(ns) = lnew1
            mosp(ns) = mnew1
        end do
    else if (iptype .eq. 3) then
        ! Have a generic ion exchange phase.
        ne = np - iern1 + 1

        do je = 1,jgext(ne)
            ns = jern1(je,ne) - 1

            do ie = 1,ngext(je,ne)
                ns = ns + 1
                mrx = mrgexs(ie,je,ne)
                lrx = tlg(mrx)
                delms = delmp*mrx
                mprsp(ns) = mprsp(ns) + delms

                if (qtotsh) then
                    !              Do a total shift of the component species.
                    ! XX           Presently a generic ion exchange phase can not be
                    ! XX           totally shifted; hence, neither can the component
                    ! XX           species.
                    lnew1 = lgexpl
                else
                    ! Do a partial shift.
                    lnew1 = lnew + lrx
                end if

                mnew1 = texp(lnew1)
                losp(ns) = lnew1
                mosp(ns) = mnew1
            end do
        end do
    end if

    ! Reset the variables used for tracking mole numbers. If a phase
    ! is totally purged from the ES, it will subsequently be necessary
    ! to update the associated index arrays and reload the
    ! corresponding data.
    do npe = 1,npet
        if (np .eq. iemop(npe)) then
            go to 110
        end if
    end do

    npe = 0
110 continue

    if (npe .gt. 0) then
        mnew = moph(np)
        emop(npe) = mnew
        emop0(npe) = mnew

        ! Zero the associated finite differences as they are
        ! no longer valid. The derivatives (demop0) previously
        ! derived from them remain valid and are not zeroed here.
        ! Zeroing the finite differences here will effectively
        ! drop the order for the mole number tracking of the
        ! affected phase.
        do n = 1,nordmx
            fdpe0(n,npe) = 0.
            fdpem1(n,npe) = 0.
        end do

        nr1 = ncmpe(1,npe)
        nr2 = ncmpe(2,npe)

        do nse = nr1,nr2
            ns = iemos(nse)
            mnew1 = mosp(ns)
            emos(nse) = mnew1
            emos0(nse) = mnew1

            ! Zero the associated finite differences as they are
            ! no longer valid. The derivatives (demos0) previously
            ! derived from them remain valid and are not zeroed here.
            ! Zeroing the finite differences here will effectively
            ! drop the order for the mole number tracking of the
            ! affected species.
            do n = 1,nordmx
                fdse0(n,nse) = 0.
                fdsem1(n,nse) = 0.
            end do
        end do
    end if

    ! Find the corresponding matrix positions. As coded here, there
    ! need not be any such.
    kstart = 0
    kend  = -1

    if (iptype .eq. 1) then
        ! Have a pure mineral.
        do kcol = km1,kmt
            if (ipndx1(kcol) .eq. np) then
                kstart = kcol
                go to 130
            end if
        end do

130 continue
        if (kstart .gt. 0) then
            kend = kstart
        end if
    else if (iptype .eq. 2) then
        ! Have a solid solution.
        do kcol = kx1,kxt
            if (ipndx1(kcol) .eq. np) then
                kstart = kcol
                go to 150
            end if
        end do

150 continue
        if (kstart .gt. 0) then
            do kcol = kxt,kstart,-1
                if (ipndx1(kcol) .eq. np) then
                    kend = kcol
                    go to 160
                end if
            end do

160 continue
        end if
    else if (iptype .eq. 3) then
        ! Have a generic ion exchange phase.
        do kcol = 1,kbt
            if (ipndx1(kcol) .eq. np) then
                kstart = kcol
                go to 170
            end if
        end do

170 continue
        if (kstart .gt. 0) then
            do kcol = kbt,kstart,-1
                if (ipndx1(kcol) .eq. np) then
                    kend = kcol
                    go to 180
                end if
            end do

180 continue
        end if
    end if

    do kcol = kstart,kend
        ns = iindx1(kcol)
        lnew1 = losp(ns)
        mnew1 = mosp(ns)
        zvclg0(kcol) = lnew1
        zvclg1(kcol) = lnew1
        zvec0(kcol) = mnew1
        zvec1(kcol) = mnew1
    end do

999 continue
end subroutine shftph