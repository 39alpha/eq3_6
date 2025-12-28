subroutine goback(acflg,acflg0,emop,emop0,emos,emos0,fje,fje0,fxi,fxi0,iemop,iemop0,iemos,iemos0,iindx0,iindx1,ipndx0,ipndx1,jpflag,jsflag,jreac,jreac0,kdim,kdim0,kmax,km1,km10,kmt,kmt0,kx1,kx10,kxt,kxt0,loph,losp,moph,moph0,mosp,mosp0,ncmpe,ncmpe0,npet,npetmx,npet0,npt,nptmax,nrct,nrctmx,nset,nsetmx,nset0,nst,nstmax,qreq,qriinf,sigmam,sigmm0,uzvec0,uzvec1,xi0,xi1)
    !! This subroutine sets up to go back to the previous point of
    !! reaction progress.
    !! This subroutine is called by:
    !!   EQ6/path.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: kmax
    integer :: npetmx
    integer :: nptmax
    integer :: nrctmx
    integer :: nsetmx
    integer :: nstmax

    integer :: iemop(npetmx)
    integer :: iemop0(npetmx)
    integer :: iemos(nsetmx)
    integer :: iemos0(nsetmx)
    integer :: iindx0(kmax)
    integer :: iindx1(kmax)
    integer :: ipndx0(kmax)
    integer :: ipndx1(kmax)
    integer :: jreac(nrctmx)
    integer :: jreac0(nrctmx)
    integer :: jpflag(nptmax)
    integer :: jsflag(nstmax)
    integer :: ncmpe(2,npetmx)
    integer :: ncmpe0(2,npetmx)

    integer :: kdim
    integer :: kdim0
    integer :: km1
    integer :: km10
    integer :: kmt
    integer :: kmt0
    integer :: kx1
    integer :: kx10
    integer :: kxt
    integer :: kxt0
    integer :: npet
    integer :: npet0
    integer :: npt
    integer :: nrct
    integer :: nset
    integer :: nset0
    integer :: nst

    logical :: qreq
    logical :: qriinf

    character(len=48) :: uzvec0(kmax)
    character(len=48) :: uzvec1(kmax)

    real(kind=8) :: acflg(nstmax)
    real(kind=8) :: acflg0(nstmax)
    real(kind=8) :: emop(npetmx)
    real(kind=8) :: emop0(npetmx)
    real(kind=8) :: emos(nsetmx)
    real(kind=8) :: emos0(nsetmx)
    real(kind=8) :: loph(nptmax)
    real(kind=8) :: losp(nstmax)
    real(kind=8) :: moph(nptmax)
    real(kind=8) :: moph0(nptmax)
    real(kind=8) :: mosp(nstmax)
    real(kind=8) :: mosp0(nstmax)

    real(kind=8) :: fje
    real(kind=8) :: fje0
    real(kind=8) :: fxi
    real(kind=8) :: fxi0
    real(kind=8) :: sigmam
    real(kind=8) :: sigmm0
    real(kind=8) :: xi0
    real(kind=8) :: xi1

    ! Local variable declarations.
    integer :: ileft
    integer :: kcol
    integer :: np
    integer :: npe
    integer :: nrc
    integer :: ns
    integer :: nse

    xi1 = xi0

    qreq = .false.
    qriinf = .false.

    kdim = kdim0
    kxt = kxt0
    kmt = kmt0
    km1 = km10
    kx1 = kx10

    do kcol = 1,kdim
        uzvec1(kcol) = uzvec0(kcol)
        iindx1(kcol) = iindx0(kcol)
        ipndx1(kcol) = ipndx0(kcol)
    end do

    do np = 1,npt
        moph(np) = moph0(np)
    end do

    ! Note that the loop is unrolled.
    ileft = (npt/8)*8

    do np = 1,ileft,8
        loph(np) = -99999.
        loph(np + 1) = -99999.
        loph(np + 2) = -99999.
        loph(np + 3) = -99999.
        loph(np + 4) = -99999.
        loph(np + 5) = -99999.
        loph(np + 6) = -99999.
        loph(np + 7) = -99999.
        jpflag(np) = max(jpflag(np),0)
        jpflag(np + 1) = max(jpflag(np + 1),0)
        jpflag(np + 2) = max(jpflag(np + 2),0)
        jpflag(np + 3) = max(jpflag(np + 3),0)
        jpflag(np + 4) = max(jpflag(np + 4),0)
        jpflag(np + 5) = max(jpflag(np + 5),0)
        jpflag(np + 6) = max(jpflag(np + 6),0)
        jpflag(np + 7) = max(jpflag(np + 7),0)
    end do

    do np = ileft + 1,npt
        loph(np) = -99999.
        jpflag(np) = max(jpflag(np),0)
    end do

    fje = fje0
    fxi = fxi0
    sigmam = sigmm0

    do ns = 1,nst
        mosp(ns) = mosp0(ns)
        acflg(ns) = acflg0(ns)
    end do

    ! Note that the loop is unrolled.
    ileft = (nst/8)*8

    do ns = 1,ileft,8
        losp(ns) = -99999.
        losp(ns + 1) = -99999.
        losp(ns + 2) = -99999.
        losp(ns + 3) = -99999.
        losp(ns + 4) = -99999.
        losp(ns + 5) = -99999.
        losp(ns + 6) = -99999.
        losp(ns + 7) = -99999.
        jsflag(ns) = max(jsflag(ns),0)
        jsflag(ns + 1) = max(jsflag(ns + 1),0)
        jsflag(ns + 2) = max(jsflag(ns + 2),0)
        jsflag(ns + 3) = max(jsflag(ns + 3),0)
        jsflag(ns + 4) = max(jsflag(ns + 4),0)
        jsflag(ns + 5) = max(jsflag(ns + 5),0)
        jsflag(ns + 6) = max(jsflag(ns + 6),0)
        jsflag(ns + 7) = max(jsflag(ns + 7),0)
    end do

    do ns = ileft + 1,nst
        losp(ns) = -99999.
        jsflag(ns) = max(jsflag(ns),0)
    end do

    if (kxt .ge. km1) then
        do kcol = km1,kxt
            ns = iindx1(kcol)
            np = ipndx1(kcol)
            jpflag(np) = -1
            jsflag(ns) = -1
        end do
    end if

    do npe = 1,npet0
        iemop(npe) = iemop0(npe)
        emop(npe) = emop0(npe)
        ncmpe(1,npe) = ncmpe0(1,npe)
        ncmpe(2,npe) = ncmpe0(2,npe)
    end do

    npet = npet0

    do nse = 1,nset0
        iemos(nse) = iemos0(nse)
        emos(nse) = emos0(nse)
    end do

    nset = nset0

    do nrc = 1,nrct
        jreac(nrc) = jreac0(nrc)
    end do
end subroutine goback