subroutine stpkmd(cbsri,cbsr1,cdac,cesri,cesr1,csigma,eact,fkrc,hact,iact,ibsrti,ibsrt1,iesrti,iesrt1,iktmax,imchmx,imech,iopt,ixrti,jcode,jreac,modr,morr,morrw1,nbt1mx,nctmax,ndact,ndctmx,noptmx,noutpt,nprob,nrct,nrctmx,nrk,nsk,nsrt,nsrtmx,ntitl1,ntitmx,nttyo,nxrt,nxrtmx,rkb,rxbari,sfcar,ssfcar,trkb,ubsri,ubsr1,ucxri,udac,uesri,uesr1,ureac,ureac1,utitl1,vreac)
    !! This subroutine sets up certain variables and arrays for writing
    !! the top half of an EQ6 input file when one of the advanced pickup
    !! file options is selected. Presently there is only one such option,
    !! which is to replace the set of reactants by a special reactant
    !! defined as the aqueous solution at the last point of reaction
    !! progress. This subroutine is somewhat analogous to EQ3NR/stpk36.f,
    !! which sets up certain variables and arrays for the advanced
    !! pickup file options for EQ3NR.
    !! This subroutine should only be called if iopt(1) is greater than
    !! zero.
    !! This subroutine is called by:
    !!   EQ6/path.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: iktmax
    integer :: imchmx
    integer :: nbt1mx
    integer :: nctmax
    integer :: ndctmx
    integer :: noptmx
    integer :: nsrtmx
    integer :: ntitmx
    integer :: nrctmx
    integer :: nxrtmx

    integer :: noutpt
    integer :: nttyo

    integer :: iact(imchmx,2,nrctmx)
    integer :: ibsrti(nsrtmx)
    integer :: iesrti(nsrtmx)
    integer :: imech(2,nrctmx)
    integer :: iopt(noptmx)
    integer :: ixrti(nxrtmx)
    integer :: jcode(nrctmx)
    integer :: jreac(nrctmx)
    integer :: ndact(imchmx,2,nrctmx)
    integer :: nrk(2,nrctmx)
    integer :: nsk(nrctmx)

    integer :: ibsrt1
    integer :: iesrt1
    integer :: nprob
    integer :: nrct
    integer :: nsrt
    integer :: ntitl1
    integer :: nxrt

    character(len=80) :: utitl1(ntitmx)
    character(len=24) :: ubsri(nbt1mx,nsrtmx)
    character(len=24) :: ubsr1(nbt1mx)
    character(len=24) :: ucxri(iktmax,nxrtmx)
    character(len=24) :: udac(ndctmx,imchmx,2,nrctmx)
    character(len=24) :: ureac(nrctmx)
    character(len=8) :: uesri(nctmax,nsrtmx)
    character(len=8) :: uesr1(nctmax)

    character(len=24) :: ureac1

    real(kind=8) :: cbsri(nbt1mx,nsrtmx)
    real(kind=8) :: cbsr1(nbt1mx)
    real(kind=8) :: cdac(ndctmx,imchmx,2,nrctmx)
    real(kind=8) :: cesri(nctmax,nsrtmx)
    real(kind=8) :: cesr1(nctmax)
    real(kind=8) :: csigma(imchmx,2,nrctmx)
    real(kind=8) :: eact(imchmx,2,nrctmx)
    real(kind=8) :: fkrc(nrctmx)
    real(kind=8) :: hact(imchmx,2,nrctmx)
    real(kind=8) :: modr(nrctmx)
    real(kind=8) :: morr(nrctmx)
    real(kind=8) :: rkb(imchmx,2,nrctmx)
    real(kind=8) :: rxbari(iktmax,nxrtmx)
    real(kind=8) :: sfcar(nrctmx)
    real(kind=8) :: ssfcar(nrctmx)
    real(kind=8) :: trkb(imchmx,2,nrctmx)
    real(kind=8) :: vreac(nrctmx)

    real(kind=8) :: morrw1

    ! Local variable declarations.
    integer :: n
    integer :: nmax

    if (iopt(20) .le. 0) then
        go to 999
    end if

    ! Set up the initial part of the new main title.
    call initcb(utitl1,ntitmx)
    ntitl1 = 25
    utitl1(1) = 'EQ6 input file name= sample.6i'
    utitl1(2) = 'Description= "Sample"'
    utitl1(3) = 'Version level= 8.0'
    utitl1(4) = 'Revised mm/dd/yy    Revisor= Username'
    utitl1(6) = '  This is a sample EQ6 input file written as an EQ6 PICKUP file.'
    utitl1(7) = 'The EQ6 input file used to generate this output is identified in'
    utitl1(8) = 'the second title given below.'
    utitl1(10) = '  You are expected to modify this EQ6 input file tomeet your own needs.'

    ! Delete all current reactants.
    nrct = 0
    call initiz(jcode,nrctmx)
    call initiz(jreac,nrctmx)
    call initiz(nsk,nrctmx)

    nsrt = 0
    call initiz(ibsrti,nsrtmx)
    call initiz(iesrti,nsrtmx)

    nmax = nbt1mx*nsrtmx
    call initaz(cbsri,nmax)
    call initcb(ubsri,nmax)

    nmax = nctmax*nsrtmx
    call initaz(cesri,nmax)
    call initcb(uesri,nmax)

    call initiz(ixrti,nxrtmx)

    nxrt = 0
    nmax = iktmax*nxrtmx
    call initaz(rxbari,nmax)
    call initcb(ucxri,nmax)

    call initcb(ureac,nrctmx)
    call initaz(fkrc,nrctmx)
    call initaz(sfcar,nrctmx)
    call initaz(ssfcar,nrctmx)

    nmax = 2*nrctmx
    call initiz(imech,nmax)
    call initiz(nrk,nmax)

    nmax = imchmx*2*nrctmx
    call initaz(csigma,nmax)
    call initaz(rkb,nmax)
    call initaz(trkb,nmax)
    call initaz(eact,nmax)
    call initaz(hact,nmax)
    call initiz(iact,nmax)
    call initiz(ndact,nmax)

    nmax = ndctmx*imchmx*2*nrctmx
    call initaz(cdac,nmax)
    call initcb(udac,nmax)

    ! Write an EQ6 input file in which "Fluid 2" is the sole reactant
    ! added to the equilibrium system (ES) containing "Fluid 1" and
    ! possible other phases. "Fluid 2" is defined as a special
    ! reactant and appears in the top half of the file. It is the
    ! aqueous solution from the end of the first reaction-path problem
    ! of one or more such problems included on the current EQ6 input
    ! file. "Fluid 1" is the aqueous solution from the end of the
    ! current such problem, and the whole ES described is that at
    ! the end of this problem.
    ! If this option is chosen for the first (perhaps the only)
    ! problem on the EQ6 input file, the pickup file produced
    ! is set up for mixing of "Fluid 2" with an ES containing an
    ! identical aqueous solution. This would normally be not very
    ! useful as a whole, thought the top and bottom halves might be
    ! separately useful.
    utitl1(12) = '  This sample file has "Fluid 2" mixing into an equilibrium system (ES)'
    utitl1(13) = 'which contains "Fluid 1". "Fluid 2" is the aqueous solution at the end'
    utitl1(14) = 'of the first EQ6 problem. "Fluid 1" is the aqueous solution at the'
    utitl1(15) = 'end of the current EQ6 problem, and the ES is the corresponding ES.'

    if (nprob .eq. 1) then
        utitl1(17) = '  Warning: This EQ6 input file is set up to mix Fl"uid 2" with an ES'
        utitl1(18) = 'containing the identical aqueous solution. You may wish to change'
        utitl1(19) = 'the bottom half of this file so that this fluid mixes with something else.'
    end if

    ! Note: "Fluid 2" as a special reactant is defined in EQ6/path.f.
    nrct = 1
    nsrt = 1
    ureac(1) = ureac1
    jcode(1) = 2
    jreac(1) = 0
    morr(1) = morrw1
    modr(1) = 0.
    nrk(1,1) = 1
    rkb(1,1,1) = 1.0

    vreac(1) = 1.0

    iesrti(1) = iesrt1

    do n = 1,iesrt1
        uesri(n,1) = uesr1(n)
        cesri(n,1) = cesr1(n)
    end do

    ibsrti(1) = ibsrt1

    do n = 1,ibsrt1
        ubsri(n,1) = ubsr1(n)
        cbsri(n,1) = cbsr1(n)
    end do

999 continue
end subroutine stpkmd