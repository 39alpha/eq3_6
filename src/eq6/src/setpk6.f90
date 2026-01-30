subroutine setpk6(actwlg,awmax,awmaxi,awmin,awmini,eh,ehmax,ehmaxi,ehmin,ehmini,fo2lg,iindx1,jflag,jflgi,kbt,kdim,kmax,kprs,mprph,mprphi,mprsp,mprspi,mtb,mtbi,mtbaq,mtbaqi,nbasp,nbaspd,nbaspi,nbti,nbtmax,ncmpr,nobswt,noutpt,nprpmx,nprpti,nprsmx,nprsti,npt,nptmax,nttyo,nstmax,o2max,o2maxi,o2min,o2mini,ph,phmax,phmaxi,phmin,phmini,prcinf,press,pressi,tempc,tempci,time1,timemx,timmxi,tistti,ubmtbi,uobsw,uphase,uprphi,uprspi,uspec,uzveci,uzvec1,xi1,ximax,ximaxi,xistti,zvclgi,zvclg1)
    !! This subroutine sets up certain variables and arrays for writing
    !! on the pickup file. This subroutine must be called prior to
    !! writing a pickup file.
    !! This subroutine is called by:
    !!   EQ6/path.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: kmax
    integer :: nbtmax
    integer :: nprpmx
    integer :: nprsmx
    integer :: nptmax
    integer :: nstmax

    integer :: noutpt
    integer :: nttyo

    integer :: iindx1(kmax)
    integer :: jflag(nstmax)
    integer :: jflgi(nbtmax)
    integer :: nbasp(nbtmax)
    integer :: nbaspd(nbtmax)
    integer :: nbaspi(nbtmax)
    integer :: ncmpr(2,nptmax)

    integer :: kbt
    integer :: kdim
    integer :: kprs
    integer :: nbti
    integer :: nobswt
    integer :: nprpti
    integer :: nprsti
    integer :: npt

    character(len=48) :: ubmtbi(nbtmax)
    character(len=48) :: uobsw(2,nbtmax)
    character(len=48) :: uprspi(nprsmx)
    character(len=48) :: uspec(nstmax)
    character(len=48) :: uzveci(kmax)
    character(len=48) :: uzvec1(kmax)
    character(len=24) :: uphase(nptmax)
    character(len=24) :: uprphi(nprpmx)

    real(kind=8) :: mprph(nptmax)
    real(kind=8) :: mprphi(nprpmx)
    real(kind=8) :: mprsp(nstmax)
    real(kind=8) :: mprspi(nprsmx)
    real(kind=8) :: mtb(nbtmax)
    real(kind=8) :: mtbi(nbtmax)
    real(kind=8) :: mtbaq(nbtmax)
    real(kind=8) :: mtbaqi(nbtmax)
    real(kind=8) :: zvclgi(kmax)
    real(kind=8) :: zvclg1(kmax)

    real(kind=8) :: actwlg
    real(kind=8) :: awmax
    real(kind=8) :: awmaxi
    real(kind=8) :: awmin
    real(kind=8) :: awmini
    real(kind=8) :: eh
    real(kind=8) :: ehmax
    real(kind=8) :: ehmaxi
    real(kind=8) :: ehmin
    real(kind=8) :: ehmini
    real(kind=8) :: fo2lg
    real(kind=8) :: o2max
    real(kind=8) :: o2maxi
    real(kind=8) :: o2min
    real(kind=8) :: o2mini
    real(kind=8) :: ph
    real(kind=8) :: phmax
    real(kind=8) :: phmaxi
    real(kind=8) :: phmin
    real(kind=8) :: phmini
    real(kind=8) :: prcinf
    real(kind=8) :: press
    real(kind=8) :: pressi
    real(kind=8) :: tempc
    real(kind=8) :: tempci
    real(kind=8) :: time1
    real(kind=8) :: timemx
    real(kind=8) :: timmxi
    real(kind=8) :: tistti
    real(kind=8) :: xi1
    real(kind=8) :: ximax
    real(kind=8) :: ximaxi
    real(kind=8) :: xistti

    ! Local variable declarations.
    integer :: j2
    integer :: kcol
    integer :: krow
    integer :: n
    integer :: nb
    integer :: nbi
    integer :: nerr
    integer :: nmax
    integer :: np
    integer :: nr1
    integer :: nr2
    integer :: ns
    integer :: nsi
    integer :: ns1
    integer :: ns2

    integer :: ilnobl

    character(len=8) :: ux8
    character(len=24) :: uaqsln

    data uaqsln /'Aqueous solution        '/

    nerr = 0

    ! Starting values of reaction progress and time.
    xistti = xi1
    tistti = time1

    ! Maximum values of reaction progress and time.
    if (xi1 .ge. ximax) then
        ximaxi = prcinf
    end if

    if (time1 .ge. timemx) then
        timmxi = prcinf
    end if

    ! Minimum and maximum values of other quantities.
    if (pH .le. phmin) then
        phmini = -prcinf
    end if

    if (pH .ge. phmax) then
        phmaxi = prcinf
    end if

    if (eh .le. ehmin) then
        ehmini = -prcinf
    end if

    if (eh .ge. ehmax) then
        ehmaxi = prcinf
    end if

    if (fo2lg .le. o2min) then
        o2mini = -prcinf
    end if

    if (fo2lg .ge. o2max) then
        o2maxi = prcinf
    end if

    if (actwlg .le. awmin) then
        awmini = -prcinf
    end if

    if (actwlg .ge. awmax) then
        awmaxi = prcinf
    end if

    ! Original temperature.
    tempci = tempc

    ! Original pressure.
    pressi = press

    ! Mass balance totals and jflag values.
    call initaz(mtbi,nbtmax)
    call initaz(mtbaqi,nbtmax)
    call initiz(jflgi,nbtmax)

    do krow = 1,kbt
        nb = iindx1(krow)
        ns = nbaspd(nb)
        nsi = nbaspi(nb)
        nbi = krow
        ubmtbi(nbi) = uspec(nsi)
        mtbi(nbi) = mtb(nb)
        mtbaqi(nbi) = mtbaq(nb)
        jflgi(nbi) = jflag(ns)
    end do

    nbti = kbt

    ! Ordinary basis switching directives.
    nmax = 2*nbtmax
    call initcb(uobsw,nmax)

    n = 0

    do kcol = 1,kbt
        nb = iindx1(kcol)
        ns1 = nbaspd(nb)
        ns2 = nbasp(nb)

        if (ns1 .ne. ns2) then
            n = n + 1
            uobsw(1,n) = uspec(ns1)
            uobsw(2,n) = uspec(ns2)
        end if
    end do

    nobswt = n

    ! Matrix column variables and corresponding values.
    call initcb(uzveci,kmax)
    call initaz(zvclgi,kmax)

    do kcol = 1,kdim
        uzveci(kcol) = uzvec1(kcol)
        zvclgi(kcol) = zvclg1(kcol)
    end do

    ! Phases and species in the PRS.
    kprs = 0
    nprpti = 0
    nprsti = 0

    call initcb(uprphi,nprpmx)
    call initcb(uprspi,nprsmx)
    call initaz(mprphi,nprpmx)
    call initaz(mprspi,nprsmx)

    do np = 1,npt
        if (uphase(np)(1:24) .ne. uaqsln(1:24)) then
            if (mprph(np) .gt. 0.) then
                nprpti = nprpti + 1

                if (nprpti .gt. nprpmx) then
                    write (ux8,'(i5)') nprpmx
                    call lejust(ux8)
                    j2 = ilnobl(ux8)
                    write (noutpt,1000) ux8(1:j2)
                    write (nttyo,1000) ux8(1:j2)
1000 format(/' * Error - (EQ6/setpk6) Have too many phases',' in the physically removed system',/7x,'(PRS) to write',' them all on the pickup file. The code is only',' dimensioned',/7x,'to allow ',a,' such phases on the',' input and pickup files. Increase',/7x,'the',' dimensioning parameter nprppa.')

                    nerr = nerr + 1
                end if

                uprphi(nprpti) = uphase(np)
                mprphi(nprpti) = mprph(np)
                nr1 = ncmpr(1,np)
                nr2 = ncmpr(2,np)

                do ns = nr1,nr2
                    if (mprsp(ns) .gt. 0.) then
                        nprsti = nprsti + 1

                        if (nprsti .gt. nprsmx) then
                            write (ux8,'(i5)') nprsmx
                            call lejust(ux8)
                            j2 = ilnobl(ux8)
                            write (noutpt,1010) ux8
                            write (nttyo,1010) ux8
1010 format(/' * Error - (EQ6/setpk6) Have too many',' species in the physically removed system',/7x,'(PRS) to write them all on the pickup file.',' The code is only dimensioned',/7x,'to allow ',a,' such species on the input and pickup files.',' Increase',/7x,'the dimensioning parameter nprspa.')

                            nerr = nerr + 1
                        end if

                        uprspi(nprsti) = uspec(ns)
                        mprspi(nprsti) = mprsp(ns)
                    end if
                end do
            end if
        end if
    end do

    if (nprpti .gt. 0) then
        kprs = 1
    end if

    if (nerr .gt. 0) then
        stop
    end if
end subroutine setpk6
