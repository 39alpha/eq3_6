subroutine dumpdp(csts,demop0,demos0,emop,emop0,emos,emos0,fdpe0,fdpem1,fdse0,fdsem1,iemop,iemos,iern1,iern2,ietmax,iindx0,iindx1,imrn1,imrn2,iodb,ipndx0,ipndx1,ixrn1,ixrn2,jcsort,jern1,jetmax,jgext,jpflag,jsflag,kbt,kdim,kdim0,kmax,km1,km10,kmt,kmt0,kord,kstep,kx1,kx10,kxt,kxt0,loph,losp,moph,mosp,mprph,mprsp,mrgexs,mtb,mtb0,nbasp,nbaspd,nbt,nbtmax,ncmpe,ncmpr,ndelay,netmax,ngext,nodbmx,nordmx,noutpt,npet,npetmx,npet0,npt,nptmax,npts,nset,nsetmx,nset0,nstmax,nsts,nstsmx,nstsr,nttyo,qbye,uaqsln,ufixf,uspec,uphase,uzvec0,uzvec1,xbar,xbarlg,zklgmn,zklogl,zvclg0,zvclg1,zvec0,zvec1)
    !! This subroutine transfers the entire mass of eligible phases
    !! in the Equilibrium System (ES) to the Physically Removed System
    !! (PRS). The eligible phases exclude the aqueous solution phase
    !! and any fictive fugacity-fixing phases. The ES phase membership
    !! and the corresponding finite-difference representations are
    !! changed accordingly.
    !! This subroutine is called by:
    !!   EQ6/path.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: ietmax
    integer :: jetmax
    integer :: kmax
    integer :: nbtmax
    integer :: netmax
    integer :: nodbmx
    integer :: nordmx
    integer :: npetmx
    integer :: nptmax
    integer :: nsetmx
    integer :: nstmax
    integer :: nstsmx

    integer :: noutpt
    integer :: nttyo

    integer :: iemop(npetmx)
    integer :: iemos(nsetmx)
    integer :: iindx0(kmax)
    integer :: iindx1(kmax)
    integer :: iodb(nodbmx)
    integer :: ipndx0(kmax)
    integer :: ipndx1(kmax)
    integer :: jcsort(nstmax)
    integer :: jern1(jetmax,netmax)
    integer :: jgext(netmax)
    integer :: jpflag(nptmax)
    integer :: jsflag(nstmax)
    integer :: nbasp(nbtmax)
    integer :: nbaspd(nbtmax)
    integer :: ncmpe(2,npetmx)
    integer :: ncmpr(2,nptmax)
    integer :: ngext(jetmax,netmax)
    integer :: nsts(nstsmx)
    integer :: nstsr(2,nstmax)

    integer :: iern1
    integer :: iern2
    integer :: imrn1
    integer :: imrn2
    integer :: ixrn1
    integer :: ixrn2
    integer :: kbt
    integer :: kdim
    integer :: kdim0
    integer :: km1
    integer :: km10
    integer :: kmt
    integer :: kmt0
    integer :: kord
    integer :: kstep
    integer :: kx1
    integer :: kx10
    integer :: kxt
    integer :: kxt0
    integer :: nbt
    integer :: ndelay
    integer :: npet
    integer :: npet0
    integer :: npt
    integer :: npts
    integer :: nset
    integer :: nset0

    logical :: qbye

    character(len=48) :: uspec(nstmax)
    character(len=48) :: uzvec0(kmax)
    character(len=48) :: uzvec1(kmax)
    character(len=24) :: uphase(nptmax)
    character(len=24) :: uaqsln
    character(len=8) :: ufixf

    real(kind=8) :: csts(nstsmx)
    real(kind=8) :: demop0(nordmx,npetmx)
    real(kind=8) :: demos0(nordmx,nsetmx)
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
    real(kind=8) :: moph(nptmax)
    real(kind=8) :: mosp(nstmax)
    real(kind=8) :: mprph(nptmax)
    real(kind=8) :: mprsp(nstmax)
    real(kind=8) :: mrgexs(ietmax,jetmax,netmax)
    real(kind=8) :: mtb(nbtmax)
    real(kind=8) :: mtb0(nbtmax)
    real(kind=8) :: xbar(nstmax)
    real(kind=8) :: xbarlg(nstmax)
    real(kind=8) :: zvclg0(kmax)
    real(kind=8) :: zvclg1(kmax)
    real(kind=8) :: zvec0(kmax)
    real(kind=8) :: zvec1(kmax)

    real(kind=8) :: zklgmn
    real(kind=8) :: zklogl

    ! Local variable declarations.
    integer :: ier
    integer :: j2
    integer :: kcol
    integer :: nmax
    integer :: np
    integer :: npe
    integer :: nr1
    integer :: nr2
    integer :: ns
    integer :: nse

    integer :: ilnobl

    logical :: qprflg
    logical :: qshftd
    logical :: qtotsh

    ! Do a total shift of any eligible ES phase whose mass has begun
    ! to decrease.
    qtotsh = .true.
    qshftd = .false.

    do npe = 1,npet
        if (fdpe0(1,npe) .lt. 0.) then
            np = iemop(npe)

            if (uphase(np)(1:24) .ne. uaqsln(1:24)) then
                if (uphase(np)(1:5) .ne. ufixf(1:5)) then
                    j2 = ilnobl(uphase(np))
                    write (noutpt,1000) uphase(np)(1:j2)
1000 format(/' --- ',a,' is no longer precipitating ---')

                    if (iodb(1) .ge. 1) then
                        write (noutpt,1010) emos(npe),fdpe0(1,npe)
1010 format(/'   Mass= ',1pe12.5,' moles',/'   Two-point relative rate= ',1pe12.5,' mol/mol',/)
                    end if

                    call shftph(emop,emop0,emos,emos0,fdpe0,fdpem1,fdse0,fdsem1,iemop,iemos,iern1,iern2,ietmax,iindx1,imrn1,imrn2,ipndx1,ixrn1,ixrn2,jern1,jetmax,jgext,jpflag,jsflag,kbt,kmax,km1,kmt,kx1,kxt,loph,losp,moph,mosp,mprph,mprsp,mrgexs,nbtmax,ncmpe,ncmpr,netmax,ngext,nordmx,noutpt,np,npet,npetmx,nptmax,nsetmx,nstmax,nttyo,qshftd,qtotsh,uphase,xbar,xbarlg,zklgmn,zklogl,zvclg0,zvclg1,zvec0,zvec1)

                    jpflag(np) = 0
                    nr1 = ncmpr(1,np)
                    nr2 = ncmpr(2,np)

                    do ns = nr1,nr2
                        if (jsflag(ns) .eq. -10) then
                            jsflag(ns) = 0
                        end if
                    end do
                end if
            end if
        end if
    end do

    if (qshftd) then
        qbye = .true.
        kord = 0
        ndelay = 1
        npts = 1
        qprflg = iodb(10) .ge. 1

        call escalc(csts,iindx1,jcsort,kbt,kmax,moph,mosp,mtb,mtb0,nbaspd,nbt,nbtmax,ncmpr,noutpt,npt,nptmax,nstmax,nsts,nstsmx,nstsr,qprflg,uspec)

        call miidxz(ier,iindx1,ipndx1,jpflag,jsflag,kbt,kdim,kmax,km1,kmt,kx1,kxt,losp,ncmpr,noutpt,npt,nptmax,nstmax,nttyo,uspec,uzvec1,zvclg1,zvec1)

        if (ier .gt. 0) then
            write (noutpt,1020)
            write (nttyo,1020)
1020 format(/' * Error - (EQ6/dumpdp) Programming error trap:',"Can't recover from",/7x,'having exceeded the maximum ',i4,' elements of the iindx1 array. This should',/7x,'not be',' possible in the present subroutine.')

            stop
        end if

        kdim0 = kdim
        kxt0 = kxt
        kmt0 = kmt
        km10 = km1
        kx10 = kx1

        do kcol = 1,kdim
            iindx0(kcol) = iindx1(kcol)
            ipndx0(kcol) = ipndx1(kcol)
            uzvec0(kcol) = uzvec1(kcol)
            zvclg0(kcol) = zvclg1(kcol)
            zvec0(kcol) = zvec1(kcol)
        end do

        ! The ES phase assemblage has changed. Re-set the index arrays
        ! associated with finite-difference description of the numbers
        ! of mole of phases and species in the Equilibrium System (ES).
        call iiemop(iemop,iemos,iindx1,ipndx1,jsflag,kdim,kmax,ncmpe,ncmpr,noutpt,npet,npetmx,npt,nptmax,nset,nsetmx,nstmax,nttyo,uaqsln,uspec,uphase)

        npet0 = npet
        nset0 = nset

        nmax = nordmx*npetmx
        call initaz(fdpe0,nmax)
        call initaz(fdpem1,nmax)
        call initaz(demop0,nmax)

        nmax = nordmx*nsetmx
        call initaz(fdse0,nmax)
        call initaz(fdsem1,nmax)
        call initaz(demos0,nmax)

        do npe = 1,npet
            np = iemop(npe)
            emop0(npe) = moph(np)
            emop(npe) = moph(np)
            nr1 = ncmpe(1,npe)
            nr2 = ncmpe(2,npe)

            do nse = nr1,nr2
                ns = iemos(nse)
                emos0(nse) = mosp(ns)
                emos(nse) = mosp(ns)
            end do
        end do
    end if
end subroutine dumpdp
