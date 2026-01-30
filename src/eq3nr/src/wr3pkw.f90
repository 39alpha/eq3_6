subroutine wr3pkw(electr,cgexj,ietmax,iopg,jetmax,jflgi,jgext,kbt,kct,kdim,kmax,kmt,kxmod,kxt,mtbi,mtbaqi,mwtges,nbti,nbtmax,net,netmax,newin,ngexrt,nobswt,nopgmx,nsbswt,ntitl2,ntitmx,nxmdmx,nxmod,qgexsh,pressi,tempci,tgexp,ubmtbi,ugexj,ugexmo,ugexp,ugexr,uhfgex,uobsw,usbsw,utitl2,uvfgex,uxkgex,uxmod,uzveci,xhfgex,xlkgex,xlkmod,xvfgex,zgexj,zvclgi)
    !! This subroutine writes the pickup file in compact ("W") format.
    !! This file is used to communicate data from EQ3NR to EQ6 (it
    !! comprises the bottom half of an EQ6 input file).
    !! This subroutine is a near clone of parts of XCON6/wr6w8.f and
    !! EQ6/wr6pkw.f.
    !! The calling sequence of this subroutine is identical to that of
    !! EQ3NR/wr3pkd.f.
    !! This subroutine is called by:
    !!   EQ3NR/eq3nr.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: ietmax
    integer :: jetmax
    integer :: kmax
    integer :: nbtmax
    integer :: netmax
    integer :: nopgmx
    integer :: ntitmx
    integer :: nxmdmx

    integer :: newin

    integer :: iopg(nopgmx)
    integer :: jgext(netmax)
    integer :: jflgi(nbtmax)
    integer :: kxmod(nxmdmx)
    integer :: ngexrt(jetmax,netmax)

    integer :: kbt
    integer :: kct
    integer :: kdim
    integer :: kmt
    integer :: kxt
    integer :: nbti
    integer :: net
    integer :: nobswt
    integer :: nsbswt
    integer :: ntitl2
    integer :: nxmod

    logical :: qgexsh

    character(len=80) :: utitl2(ntitmx)
    character(len=56) :: ugexr(ietmax,jetmax,netmax)
    character(len=48) :: ubmtbi(nbtmax)
    character(len=48) :: uobsw(2,nbtmax)
    character(len=48) :: usbsw(2,nbtmax)
    character(len=48) :: uxmod(nxmdmx)
    character(len=48) :: uzveci(kmax)
    character(len=24) :: ugexmo(netmax)
    character(len=24) :: ugexp(netmax)
    character(len=8) :: ugexj(jetmax,netmax)
    character(len=8) :: uhfgex(ietmax,jetmax,netmax)
    character(len=8) :: uvfgex(ietmax,jetmax,netmax)
    character(len=8) :: uxkgex(ietmax,jetmax,netmax)

    real(kind=8) :: cgexj(jetmax,netmax)
    real(kind=8) :: mtbaqi(nbtmax)
    real(kind=8) :: mtbi(nbtmax)
    real(kind=8) :: mwtges(netmax)
    real(kind=8) :: tgexp(netmax)
    real(kind=8) :: xhfgex(ietmax,jetmax,netmax)
    real(kind=8) :: xlkgex(ietmax,jetmax,netmax)
    real(kind=8) :: xlkmod(nxmdmx)
    real(kind=8) :: xvfgex(ietmax,jetmax,netmax)
    real(kind=8) :: zgexj(jetmax,netmax)
    real(kind=8) :: zvclgi(kmax)

    real(kind=8) :: electr
    real(kind=8) :: pressi
    real(kind=8) :: tempci

    ! Local variable declarations.
    integer :: i
    integer :: je
    integer :: j2
    integer :: j3
    integer :: j4
    integer :: kcol
    integer :: kprs
    integer :: krow
    integer :: n
    integer :: nbi
    integer :: ne

    integer :: ilnobl

    character(len=8) :: uendit

    data uendit /'endit.  '/

    ! The items in this block correspond to things defined in the
    ! first part of XCON6/wr6w8.f or EQ6/wr6pkw.f.
    j3 = ilnobl(uendit)

    ! Comment header for iopt, iopr, iodb, and iopg option switches.
1200 format('*',15x,'1    2    3    4    5    6    7    8    9   10')

1020 format(a)

    ! Old title.
    write (newin,2590)
2590 format('* Title of previous run')

    do n = 1,ntitl2
        j2 = ilnobl(utitl2(n))
        write (newin,1020) utitl2(n)(1:j2)
    end do

    if (ntitl2 .lt. ntitmx) then
        write (newin,1020) uendit(1:j3)
    end if

    ! Special basis switches.
    write (newin,2610)
2610 format('* Special basis switches')

    write (newin,2630) nsbswt
2630 format(4x,'nsbswt= ',i3)

    do n = 1,nsbswt
        j2 = ilnobl(usbsw(1,n))
        write (newin,2650) usbsw(1,n)(1:j2)
2650 format('species= ',a)

        j2 = ilnobl(usbsw(2,n))
        write (newin,2670) usbsw(2,n)(1:j2)
2670 format(2x,'switch with= ',a)
    end do

    ! Original temperature.
    write (newin,2680)
2680 format('* Original temperature and pressure')

    write (newin,2690) tempci
2690 format(4x,'tempci= ',1pe12.5)

    ! Original pressure.
    write (newin,2710) pressi
2710 format(4x,'pressi= ',1pe12.5)

    ! Ion exchanger creation.
    write (newin,2790)
2790 format('* Ion exchangers')

    write (newin,2810) net
2810 format(7x,'net= ',i3)

    do ne = 1,net
        j2 = ilnobl(ugexp(ne))
        write (newin,2840) ugexp(ne)(1:j2)
2840 format(5x,'ugexp= ',a)

        write (newin,2860) mwtges(ne)
2860 format(4x,'mwtges= ',1pe12.5)

        j3 = ilnobl(ugexmo(ne))
        write (newin,2870) ugexmo(ne)(1:j3)
2870 format(4x,'ugexmo= ',a)

        write (newin,2880) tgexp(ne)
2880 format(5x,'tgexp= ',1pe12.5)

        write (newin,2900) jgext(ne)
2900 format(5x,'jgext= ',i3)

        do je = 1,jgext(ne)
            j3 = ilnobl(ugexj(je,ne))
            write (newin,2940) ugexj(je,ne)(1:j3)
2940 format(5x,'ugexj= ',a)

            write (newin,2960) cgexj(je,ne),zgexj(je,ne)
2960 format(5x,'cgexj= ',1pe12.5,5x,'zgexj= ',e12.5)

            write (newin,3010) ngexrt(je,ne)
3010 format(4x,'ngexrt= ',i3)

            do n = 1,ngexrt(je,ne)
                j2 = ilnobl(ugexr(n,je,ne))
                write (newin,3040) ugexr(n,je,ne)(1:j2)
3040 format(5x,'ugexr= ',a)

                j4 = ilnobl(uxkgex(n,je,ne))
                write (newin,3060) xlkgex(n,je,ne),uxkgex(n,je,ne)(1:j4)
3060 format(4x,'xlkgex= ',1pe12.5,5x,'units= ',a)

                j4 = ilnobl(uhfgex(n,je,ne))
                write (newin,3070) xhfgex(n,je,ne),uhfgex(n,je,ne)(1:j4)
3070 format(4x,'xhfgex= ',1pe12.5,5x,'units= ',a)

                j4 = ilnobl(uvfgex(n,je,ne))
                write (newin,3080) xvfgex(n,je,ne),uvfgex(n,je,ne)(1:j4)
3080 format(4x,'xvfgex= ',1pe12.5,5x,'units= ',a)
            end do
        end do
    end do

    ! Number of nxmod options.
    write (newin,3110)
3110 format('* Alter/suppress options')

    write (newin,3130) nxmod
3130 format(5x,'nxmod= ',i2)

    ! Nxmod options.
    do n = 1,nxmod
        j2 = ilnobl(uxmod(n))
        write (newin,3160) uxmod(n)(1:j2),kxmod(n),xlkmod(n)
3160 format(3x,'species= ',a,/4x,'option= ',i2,14x,'xlkmod= ',1pe12.5)
    end do

    ! Iopg options.
    ! Note: iopg(1) = iopg1, etc.
    write (newin,3172)
3172 format('* Iopg options')

    write (newin,1200)

    write (newin,3180) (iopg(i), i = 1,20)
3180 format(2x,'iopg1-10= ',10i5,/1x,'iopg11-20= ',10i5)

    ! Index limits.
    write (newin,3210)
3210 format('* Matrix index limits')

    kprs = 0
    write (newin,3230) kct,kbt,kmt,kxt,kdim,kprs
3230 format(7x,'kct= ',i2,17x,'kbt= ',i2,17x,'kmt= ',i2,/ 7x,'kxt= ',i2,16x,'kdim= ',i2,16x,'kprs= ',i2)

    ! Species for which mass balances are defined.
    !   ubmtbi = names of the data file basis species for which mass
    !              balances are defined
    !   jflgi  = jflag input for basis species
    !      0 = Retain as an active basis species
    !     30 = Convert to a dependent species; fold the mass balance
    !            total for this species into the mass balance totals
    !            of basis species which remain active
    write (newin,3330)
3330 format('* Data file basis species and jflgi values')

    do nbi = 1,nbti
        write (newin,3350) ubmtbi(nbi),jflgi(nbi)
3350 format(3x,a48,3x,i2)
    end do

    ! Mass balance totals.
    write (newin,3400)
3400 format('* Mass balance totals (ES and aqueous solution)')

    do nbi = 1,nbti
        write (newin,3420) mtbi(nbi),mtbaqi(nbi)
3420 format(3x,1pe25.15,3x,e25.15)
    end do

    write (newin,3440) electr
3440 format(9x,'Electrical imbalance= ',1pe25.15)

    ! Ordinary basis switches.
    write (newin,3450)
3450 format('* Ordinary basis switches')

    write (newin,3460) nobswt
3460 format(4x,'nobswt= ',i3)

    do n = 1,nobswt
        j2 = ilnobl(uobsw(1,n))
        write (newin,2650) uobsw(1,n)(1:j2)
        j2 = ilnobl(uobsw(2,n))
        write (newin,2670) uobsw(2,n)(1:j2)
    end do

    ! Matrix column variables.
    write (newin,3470)
3470 format('* Matrix species or entities')

    do krow = 1,kdim
        j2 = ilnobl(uzveci(krow))
        write (newin,3490) uzveci(krow)(1:j2)
3490 format(3x,a)
    end do

    ! Matrix variable values.
    write (newin,3510)
3510 format('* Values of matrix variables')

    do kcol = 1,kdim
        write (newin,3530) zvclgi(kcol)
3530 format(3x,1pe22.15)
    end do

    ! There is no PRS in EQ3NR, so the block of data for phases and
    ! species in the PRS is not written here.
end subroutine wr3pkw
