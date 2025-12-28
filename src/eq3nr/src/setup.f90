subroutine setup(coval,eh,ehfac,ier,irdxc3,itdsf3,jflag,mwtsp,narn1,nbaspd,nbt,nbtmax,noutpt,nstmax,nttyo,pe,rho,tdspkg,tdspl,uspec)
    !! This subroutine converts input coval data which are not on the
    !! molal concentration scale to that scale.
    !! This subroutine is called by:
    !!   EQ3NR/eq3nr.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: nbtmax
    integer :: nstmax

    integer :: noutpt
    integer :: nttyo

    integer :: jflag(nstmax)
    integer :: nbaspd(nbtmax)

    integer :: ier
    integer :: irdxc3
    integer :: itdsf3
    integer :: narn1
    integer :: nbt

    character(len=48) :: uspec(nstmax)

    real(kind=8) :: coval(nbtmax)
    real(kind=8) :: mwtsp(nstmax)

    real(kind=8) :: eh
    real(kind=8) :: ehfac
    real(kind=8) :: pe
    real(kind=8) :: rho
    real(kind=8) :: tdspkg
    real(kind=8) :: tdspl

    ! Local variable declarations.
    integer :: jfl
    integer :: jlen
    integer :: nb
    integer :: ns

    character(len=56) :: uspn56

    real(kind=8) :: wfh2o

    ier = 0

    if (itdsf3 .ge. 1) then
        tdspkg = tdspl/rho
    end if

    wfh2o = (1000. - 0.001*tdspkg)/1000.

    if (irdxc3 .eq. -2) then
        ! Convert input pe to Eh.
        eh = pe*ehfac
        irdxc3 = -1
    end if

    do nb = 1,nbt
        ns = nbaspd(nb)
        jfl = jflag(ns)

        if (ns .eq. narn1) then
            go to 100
        end if

        if (jfl .eq. 0) then
            go to 100
        end if

        if (jfl .eq. 30) then
            go to 100
        end if

        if (jfl .eq. 7) then
            go to 100
        end if

        if (jfl.eq.16 .or. jfl.eq.17 .or. jfl.eq.18) then
            go to 100
        end if

        if (jfl.eq.19 .or. jfl.eq.20 .or. jfl.eq.21) then
            go to 100
        end if

        if (jfl.eq.22 .or. jfl.eq.23) then
            go to 100
        end if

        if (jfl.eq.25 .or. jfl.eq.27) then
            go to 100
        end if

        if (jfl .eq.- 1) then
            go to 100
        end if

        if (jfl .eq. 1) then
            ! Convert molarity to molality.
            coval(nb) = coval(nb)/(rho*wfh2o)
            jflag(ns) = 0
        else if (jfl .eq. 2) then
            ! Convert mg/L to molality.
            coval(nb) = coval(nb)*1.e-3/(mwtsp(ns)*rho*wfh2o)
            jflag(ns) = 0
        else if (jfl .eq. 3) then
            ! Convert mg/kg.sol to molality.
            coval(nb) = coval(nb)*1.e-3/(mwtsp(ns)*wfh2o)
            jflag(ns) = 0
        else if (jfl .eq. 8) then
            ! Convert alkalinity in eq/L to eq/kg.H2O.
            coval(nb) = coval(nb)/(rho*wfh2o)
            jflag(ns) = 7
        else if (jfl .eq. 9) then
            ! Convert alkalinity in eq/kg.sol to eq/kg.H2O.
            coval(nb) = coval(nb)/wfh2o
            jflag(ns) = 7
        else if (jfl .eq. 10) then
            ! Convert alkalinity in mg/L CaCO3 to eq/kg.H2O.
            coval(nb) = coval(nb)/(50000.*rho*wfh2o)
            jflag(ns) = 7
        else if (jfl .eq. 11) then
            ! Convert alkalinity in mg/L HCO3- to eq/kg.H2O.
            coval(nb) = coval(nb)/(60960*rho*wfh2o)
            jflag(ns) = 7
        else
            ! No match.
            ! Calling sequence substitutions:
            !   uspec(ns) for unam48
            call fmspnx(jlen,uspec(ns),uspn56)
            write (noutpt,1000) jfl,uspn56(1:jlen)
            write (nttyo,1000) jfl,uspn56(1:jlen)
1000 format(/' * Error - (EQ3NR/setup) An undefined jflag value',' of ',i3,/7x,'was specified on the input file for the',' species',/7x,a,'.')

            ier = 1
        end if

100 continue
    end do
end subroutine setup