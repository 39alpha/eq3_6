subroutine pshfta(csts,emop,emop0,emos,emos0,fdpe0,fdpem1,fdse0,fdsem1,iemop,iemos,iern1,iern2,ietmax,iindx1,imrn1,imrn2,iodb,ipndx1,ixrn1,ixrn2,jcsort,jern1,jetmax,jgext,jpflag,jsflag,kbt,km1,kmax,kmt,kx1,kxt,loph,losp,moph,mosp,mprph,mprsp,mrgexs,mtb,mtb0,nbasp,nbaspd,nbt,nbtmax,ncmpe,ncmpr,netmax,ngext,nodbmx,nordmx,noutpt,npet,npetmx,npt,nptmax,nsetmx,nstmax,nsts,nstsmx,nstsr,nttyo,uaqsln,ufixf,uphase,uspec,xbar,xbarlg,zklgmn,zklogl,zvclg0,zvclg1,zvec0,zvec1)
    !! This subroutine oversees the partial transfer of eligible phases
    !! (as a group) from the Equilibrium System (ES) to the Physically
    !! Removed System (PRS). The eligible phases exclude the aqueous
    !! solution phase and any fictive fugacity-fixing phases. Partial
    !! transfer means that some mass of each phase transferred remains
    !! in the ES.
    !! This subroutine is called by:
    !!   EQ6/dumpdp.f
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
    integer :: iindx1(kmax)
    integer :: iodb(nodbmx)
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
    integer :: km1
    integer :: kmt
    integer :: kx1
    integer :: kxt
    integer :: nbt
    integer :: npet
    integer :: npt

    character(len=48) :: uspec(nstmax)
    character(len=24) :: uphase(nptmax)
    character(len=24) :: uaqsln
    character(len=8) :: ufixf

    real(kind=8) :: csts(nstsmx)
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
    real(kind=8) :: zvclg0(kmax)
    real(kind=8) :: zvclg1(kmax)
    real(kind=8) :: zvec0(kmax)
    real(kind=8) :: zvec1(kmax)
    real(kind=8) :: xbar(nstmax)
    real(kind=8) :: xbarlg(nstmax)

    real(kind=8) :: zklgmn
    real(kind=8) :: zklogl

    ! Local variable declarations.
    logical :: qprflg
    logical :: qshftd
    logical :: qtotsh

    integer :: np
    integer :: npe
    integer :: nshftd

    nshftd = 0
    qtotsh = .false.

    write (noutpt,1000)
1000 format(' --- Shifting (partial) ES solid(s) to the PRS ---',/)

    ! Shift all ES phases except the aqueous solution phase and
    ! any fictive fugacity-fixing phases.
    do npe = 1,npet
        np = iemop(npe)

        if (uphase(np)(1:24) .ne. uaqsln(1:24)) then
            if (uphase(np)(1:5) .ne. ufixf(1:5)) then
                call shftph(emop,emop0,emos,emos0,fdpe0,fdpem1,fdse0,fdsem1,iemop,iemos,iern1,iern2,ietmax,iindx1,imrn1,imrn2,ipndx1,ixrn1,ixrn2,jern1,jetmax,jgext,jpflag,jsflag,kbt,kmax,km1,kmt,kx1,kxt,loph,losp,moph,mosp,mprph,mprsp,mrgexs,nbtmax,ncmpe,ncmpr,netmax,ngext,nordmx,noutpt,np,npet,npetmx,nptmax,nsetmx,nstmax,nttyo,qshftd,qtotsh,uphase,xbar,xbarlg,zklgmn,zklogl,zvclg0,zvclg1,zvec0,zvec1)

                if (qshftd) then
                    nshftd = nshftd + 1
                end if
            end if
        end if
    end do

    if (nshftd .gt. 0) then
        ! Recompute the composition of the ES.
        qprflg = iodb(10) .ge. 1
        call escalc(csts,iindx1,jcsort,kbt,kmax,moph,mosp,mtb,mtb0,nbaspd,nbt,nbtmax,ncmpr,noutpt,npt,nptmax,nstmax,nsts,nstsmx,nstsr,qprflg,uspec)
    end if
end subroutine pshfta
