subroutine hpsat(acflg,act,actlg,afcnst,affp,affs,apx,bpx,cdrs,eps100,iapxmx,ibpxmx,ier,iktmax,ixrn1,jflag,jpflag,jsflag,jsol,ncmpr,ndrs,ndrsmx,ndrsr,noutpt,np,nptmax,nstmax,nttyo,nxrn1,nxrn2,nxtmax,sidrsp,sidrph,uphase,uspec,wfac,xbar,xbarlg,xlks)
    !! This subroutine calculates the most stable (least soluble)
    !! composition of a given solid solution, given the composition
    !! of the aqueous phase it is in equilibrium with.
    !! This subroutine is called by:
    !!   EQ3NR/scripx.f
    !!   EQ6/satchk.f
    !! Principal input:
    !!   apx    = array of temperature-independent coefficients for
    !!              solid solution activity coefficients
    !!   bpx    = array of site-mixing parameters for solid solution
    !!              activity coefficients
    !!   ixrn1  = start of solid solution phases in the list of phases
    !!   ixrn2  = end of solid solution phases in the list of phases
    !!   jsol   = array of activity coefficient model indices for solid
    !!              solutions
    !!   ncmpr  = species range pointer array for phases
    !!   uphase = array of phase names
    !!   uspec  = array of species names
    !!   wterm  = array of temperature-dependent coefficients for
    !!              solid solution activity coefficients
    !! Principal output:
    !!   acflg  = array of activity coefficients
    !!   act    = array of activities
    !!   actlg  = array of log activity values
    !!   affp   = array of affinities for phases
    !!   affs   = array of affinities for species
    !!   sidrph = array of saturation indices for phases
    !!   sidrsp = array of saturation indices for species
    !!   xbar   = array of mole fractions
    !!   xbarlg = array of log mole fraction values
    !! Local:
    !!   sisppu = saturation index of a pure end-member, normalized
    !!              so that one mole of end-member is destroyed in the
    !!              corresponding reaction
    implicit none

    ! Calling sequence variable declarations.
    integer :: iapxmx
    integer :: ibpxmx
    integer :: iktmax
    integer :: ndrsmx
    integer :: nptmax
    integer :: nstmax
    integer :: nxtmax

    integer :: noutpt
    integer :: nttyo

    integer :: jflag(nstmax)
    integer :: jpflag(nptmax)
    integer :: jsflag(nstmax)
    integer :: jsol(nxtmax)
    integer :: ncmpr(2,nptmax)
    integer :: ndrs(ndrsmx)
    integer :: ndrsr(2,nstmax)

    integer :: ier
    integer :: ixrn1
    integer :: np
    integer :: nxrn1
    integer :: nxrn2

    character(len=48) :: uspec(nstmax)
    character(len=24) :: uphase(nptmax)

    real(kind=8) :: acflg(nstmax)
    real(kind=8) :: actlg(nstmax)
    real(kind=8) :: act(nstmax)
    real(kind=8) :: affp(nptmax)
    real(kind=8) :: affs(nstmax)
    real(kind=8) :: apx(iapxmx,nxtmax)
    real(kind=8) :: bpx(ibpxmx,nxtmax)
    real(kind=8) :: cdrs(ndrsmx)
    real(kind=8) :: sidrsp(nstmax)
    real(kind=8) :: sidrph(nptmax)
    real(kind=8) :: wfac(iktmax,nxtmax)
    real(kind=8) :: xbar(nstmax)
    real(kind=8) :: xbarlg(nstmax)
    real(kind=8) :: xlks(nstmax)

    real(kind=8) :: afcnst
    real(kind=8) :: eps100

    ! Local variable declarations.
    integer :: nactc
    integer :: nrn1
    integer :: nr1
    integer :: nr2
    integer :: ns
    integer :: nx
    integer :: ns2

    real(kind=8) :: af
    real(kind=8) :: cx
    real(kind=8) :: si
    real(kind=8) :: siph
    real(kind=8) :: sisp
    real(kind=8) :: sisppu
    real(kind=8) :: stx
    real(kind=8) :: stxi
    real(kind=8) :: stxm1
    real(kind=8) :: xl
    real(kind=8) :: xx
    real(kind=8) :: cx2

    real(kind=8) :: texp
    real(kind=8) :: tlg

    ier = 0
    nx = np - ixrn1 + 1

    nr1 = ncmpr(1,np)
    nr2 = ncmpr(2,np)

    affp(np) = -9999999.
    sidrph(np) = -9999999.

    ! Get the affinities and saturation indices of the pure
    ! end-members.
    do ns = nr1,nr2
        xbar(ns) = 1.0
        xbarlg(ns) = 0.
        acflg(ns) = 0.
        act(ns) = 1.0
        actlg(ns) = 0.

        call afcalc(actlg,af,afcnst,cdrs,jflag,jsflag,ndrs,ndrsmx,ndrsr,ns,nstmax,si,xlks)
        affs(ns) = af
        sidrsp(ns) = si
    end do

    ! Count the number of active components.
    nactc = 0

    do ns = nr1,nr2
        if (jsflag(ns) .le. 0) then
            nactc = nactc + 1
        end if
    end do

    ! Check the value of the site-mixing parameter.
    if (bpx(1,nx) .le. 0.) then
        bpx(1,nx) = 1.0
    end if

    stx = bpx(1,nx)

    ! Set coefficients/exponents based on the site-mixing parameter
    if (abs(stx - 1.0) .gt. eps100) then
        stxi = 1./bpx(1,nx)
        stxm1 = stx - 1.0
    else
        stx = 1.0
        stxi = 1.0
        stxm1 = 0.0
    end if

    ! Calculate the mole fractions for the case of an ideal molecular-
    ! or site-mixing solution.
    do ns = nr1,nr2
        if (jsflag(ns) .le. 0) then
            nrn1 = ndrsr(1,ns)
            cx = (-1.0/cdrs(nrn1))
            xbar(ns) = 0.0

            do ns2 = nr1,nr2
                if (ns2 .ne. ns) then
                    nrn1 = ndrsr(1,ns2)
                    cx2 = (-1.0/cdrs(nrn1))
                    xx = cx2*stxi*sidrsp(ns2)
                    xx = xx - cx*stxi*sidrsp(ns)
                    xbar(ns) = xbar(ns) + texp(xx)
                end if
            end do

            xbar(ns) = 1.0 / (1.0 + xbar(ns))
            xbarlg(ns) = tlg(xbar(ns))
        end if
    end do

    ! Calculate the activities, activity coefficients, affinities,
    ! and saturation indices for the case of an ideal molecular- or
    ! site-mixing solution.
    siph = 0.

    do ns = nr1,nr2
        if (jsflag(ns) .le. 0) then
            nrn1 = ndrsr(1,ns)
            cx = (-1.0/cdrs(nrn1))
            sisppu = cx*sidrsp(ns)
            xx = xbar(ns)
            xl = xbarlg(ns)
            sisp = sisppu - stx*xl
            acflg(ns) = stxm1*xl
            act(ns) = xx**stx
            actlg(ns) = stx*xl
            sisp = sisp/cx
            sidrsp(ns) = sisp
            affs(ns) = afcnst*sisp
            siph = siph + xx*sisp
        end if
    end do

    sidrph(np) = siph
    affp(np) = afcnst*siph
end subroutine hpsat
