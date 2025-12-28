subroutine miidxz(ier,iindx1,ipndx1,jpflag,jsflag,kbt,kdim,kmax,km1,kmt,kx1,kxt,losp,ncmpr,noutpt,npt,nptmax,nstmax,nttyo,uspec,uzvec1,zvclg1,zvec1)
    !! This subroutine modifies the matrix indexing originally read in
    !! from the input file. Modification occurs whenever a phase
    !! boundary is crossed.
    !! This subroutine is called by:
    !!   EQ6/dumpdp.f
    !!   EQ6/eqcalc.f
    !!   EQ6/setffg.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: kmax
    integer :: nptmax
    integer :: nstmax

    integer :: iindx1(kmax)
    integer :: ipndx1(kmax)
    integer :: jpflag(nptmax)
    integer :: jsflag(nstmax)
    integer :: ncmpr(2,nptmax)

    integer :: ier
    integer :: kbt
    integer :: kdim
    integer :: km1
    integer :: kmt
    integer :: kx1
    integer :: kxt
    integer :: noutpt
    integer :: npt
    integer :: nttyo

    character(len=48) :: uspec(nstmax)
    character(len=48) :: uzvec1(kmax)

    real(kind=8) :: losp(nstmax)
    real(kind=8) :: zvclg1(kmax)
    real(kind=8) :: zvec1(kmax)

    ! Local sequence variable declarations.
    integer :: jlen
    integer :: kcol
    integer :: nr1
    integer :: nr2
    integer :: np
    integer :: ns
    integer :: nt

    character(len=56) :: uspn56

    real(kind=8) :: lxx

    real(kind=8) :: texp

    ier = 0

    kcol = kbt

    ! Do non-aqueous phases with only one species.
    do np = 1,npt
        if (jpflag(np) .eq. -1) then
            nr1 = ncmpr(1,np)
            nr2 = ncmpr(2,np)
            nt = nr2 - nr1 + 1

            if (nt .eq. 1) then
                ns = nr1
                kcol = kcol + 1

                if (kcol .le. kmax) then
                    iindx1(kcol) = ns
                    ipndx1(kcol) = np
                    uzvec1(kcol) = uspec(ns)
                    jsflag(ns) = -1
                    lxx = losp(ns)
                    zvclg1(kcol) = lxx
                    zvec1(kcol) = texp(lxx)
                else
                    ! Calling sequence substitutions:
                    !   uspec(ns) for unam48
                    call fmspnm(jlen,uspec(ns),uspn56)

                    write (noutpt,1000) kmax,uspn56(1:jlen)
                    write (nttyo,1000) kmax,uspn56(1:jlen)
1000 format(/' * Warning - (EQ6/miidxz) Have exceeded the',' maximum ',i4,' elements of the',/7x,'iindx1 array',' while trying to add ',a,/7x,'as a member of the',' equilibrium system. Increase the dimensioning',/7x,'parameter kpar.')

                    ier = 1
                end if
            end if
        end if
    end do

    km1 = kbt + 1
    kmt = kcol

    ! Do non-aqueous phases with more than one species.
    do np = 2,npt
        if (jpflag(np) .eq. -1) then
            nr1 = ncmpr(1,np)
            nr2 = ncmpr(2,np)
            nt = nr2 - nr1 + 1

            if (nt .ge. 2) then
                do ns = nr1,nr2
                    if (jsflag(ns) .eq. -1) then
                        kcol = kcol + 1

                        if (kcol .le. kmax) then
                            iindx1(kcol) = ns
                            ipndx1(kcol) = np
                            uzvec1(kcol) = uspec(ns)
                            jsflag(ns) = -1
                            lxx = losp(ns)
                            zvclg1(kcol) = lxx
                            zvec1(kcol) = texp(lxx)
                        else
                            ! Calling sequence substitutions:
                            !   uspec(ns) for unam48
                            call fmspnm(jlen,uspec(ns),uspn56)

                            write (noutpt,1000) uspn56(1:jlen)
                            write (nttyo,1000) uspn56(1:jlen)
                            stop
                        end if
                    end if
                end do
            end if
        end if
    end do

    kx1 = kmt + 1
    kxt = kcol
    kdim = kcol
end subroutine miidxz