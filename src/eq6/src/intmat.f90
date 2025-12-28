subroutine intmat(iaqsln,iindx1,ipndx1,kbt,kdim,kelect,khydr,khydx,kmax,km1,kmt,ko2gaq,kwater,kx1,kxt,narn1,narn2,nbasp,nbt,nbti,nbtmax,ncmpr,nelect,nern1,nern2,nhydr,nhydx,nobswt,noutpt,no2gaq,nphasx,npt,nptmax,nstmax,nttyo,qloffg,ubmtbi,ufixf,uobsw,uphase,uspec,uzveci,uzvec1,zvclgi,zvclg1,zvec1)
    !! This routine interprets matrix variables read from the input
    !! file. It constructs the master variable array zvclgi and
    !! builds iindx1, the master variable index array, deleting any
    !! left over fixed fugacity phases. Note that iindx1(kcol) gives
    !! the basis index for kcol = 1,kbt. For kcol = km1,kxt, it gives
    !! the species index. this subroutine also constructs the ipndx1
    !! array. Note that ipndx1(kcol) gives the phase index for
    !! kcol = 1,kxt.
    !! This routine also resolves any issues between the matrix variable
    !! input (which may imply instances of ordinary basis switching) and
    !! the ordinary basis switching input.
    !! This subroutine is called by:
    !!   EQ6/eq6.f
    !! Principal input:
    !! Principal output:
    !!   qloffg = logical flag, = .true. if left over fixed fugacity
    !!              phases are in the system
    implicit none

    ! Calling sequence variable declarations.
    integer :: kmax
    integer :: nbtmax
    integer :: nptmax
    integer :: nstmax

    integer :: noutpt
    integer :: nttyo

    integer :: iindx1(kmax)
    integer :: ipndx1(kmax)
    integer :: nbasp(nbtmax)
    integer :: ncmpr(2,nptmax)
    integer :: nphasx(nstmax)

    integer :: iaqsln
    integer :: kbt
    integer :: kdim
    integer :: kelect
    integer :: khydr
    integer :: khydx
    integer :: km1
    integer :: kmt
    integer :: ko2gaq
    integer :: kwater
    integer :: kx1
    integer :: kxt
    integer :: narn1
    integer :: narn2
    integer :: nbt
    integer :: nbti
    integer :: nelect
    integer :: nern1
    integer :: nern2
    integer :: nhydr
    integer :: nhydx
    integer :: nobswt
    integer :: no2gaq
    integer :: npt

    logical :: qloffg

    character(len=48) :: ubmtbi(nbtmax)
    character(len=48) :: uobsw(2,nbtmax)
    character(len=48) :: uspec(nstmax)
    character(len=48) :: uzveci(kmax)
    character(len=48) :: uzvec1(kmax)
    character(len=24) :: uphase(nptmax)
    character(len=8) :: ufixf

    real(kind=8) :: zvclgi(kmax)
    real(kind=8) :: zvclg1(kmax)
    real(kind=8) :: zvec1(kmax)

    ! Local variable declarations.
    integer :: jlena
    integer :: jlenb
    integer :: jlenc
    integer :: jlend
    integer :: jlen2
    integer :: j2
    integer :: kbts
    integer :: kcol
    integer :: km1s
    integer :: kmts
    integer :: krow
    integer :: kx1s
    integer :: kxts
    integer :: n
    integer :: nb
    integer :: nbi
    integer :: nerr
    integer :: nn
    integer :: np
    integer :: nr1
    integer :: nr2
    integer :: ns
    integer :: ns1
    integer :: ns2
    integer :: nt

    integer :: ilnobl

    character(len=56) :: uspa56
    character(len=56) :: uspb56
    character(len=56) :: uspc56
    character(len=56) :: uspd56
    character(len=56) :: uspn56
    character(len=24) :: ux24
    character(len=8) :: ux8

    real(kind=8) :: lxx

    real(kind=8) :: texp

    nerr = 0
    qloffg = .false.
    km1s = km1
    kmts = kmt
    kx1s = kx1
    kxts = kxt
    kbts = kbt

    ! Aqueous species and generic ion exchanger species.
    kcol = 0

    do krow = 1,kbts
        do nb = 1,nbt
            ns = nbasp(nb)

            if (ns.ge.nern1 .and. ns.le.nern2) then
                ns = ns + 1
            end if

            if (uzveci(krow)(1:48) .eq. uspec(ns)(1:48)) then
                kcol = kcol + 1
                iindx1(kcol) = nb

                if (ns.ge.narn1 .and. ns.le.narn2) then
                    ! Have an aqueous species.
                    ipndx1(kcol) = iaqsln
                else if (ns.ge.nern1 .and. ns.le.nern2) then
                    ! Have a generic ion exchanger species.
                    ipndx1(kcol) = nphasx(ns)
                end if

                uzvec1(kcol) = uspec(ns)
                lxx = zvclgi(krow)
                zvclg1(kcol) = lxx
                zvec1(kcol) = texp(lxx)
                go to 100
            end if
        end do

        ! Did not find the listed active basis species in the
        ! "data file" set. It may be a species that is to be
        ! switched in by an ordinary basis switch.
        if (krow .le. nbti) then
            nbi = krow

            do ns2 = narn1,narn2
                if (uzveci(krow)(1:48) .eq. uspec(ns2)(1:48)) then
                    ! Calling sequence substitutions:
                    !   jlenb for jlen2
                    !   uspb56 for uspn56
                    call fmspnx(jlenb,uzveci(krow),uspb56)

                    do nb = 1,nbt
                        ns1 = nbasp(nb)

                        if (ubmtbi(nbi)(1:48) .eq. uspec(ns1)(1:48)) then
                            ! A basis switch is implied, with the species whose
                            ! index is ns2 being switched in for the one whose
                            ! index is ns1. The actual switch will be made later.
                            ! For the moment, put the name of the species to be
                            ! switched out (the "data file" basis species) in the
                            ! current position in the uzvec1 array.
                            ! Calling sequence substitutions:
                            !   jlena for jlen2
                            !   uspa56 for uspn56
                            call fmspnx(jlena,ubmtbi(nbi),uspa56)

                            kcol = kcol + 1
                            iindx1(kcol) = nb
                            ipndx1(kcol) = iaqsln
                            uzvec1(kcol) = uspec(ns1)
                            lxx = zvclgi(krow)
                            zvclg1(kcol) = lxx
                            zvec1(kcol) = texp(lxx)

                            ! Check for explicit basis switch directive.
                            do n = 1,nobswt
                                if (uobsw(1,n)(1:48) .eq. uspec(ns1)(1:48)) then
                                    if (uobsw(2,n)(1:48) .eq. uspec(ns2)(1:48)) then
                                        ! The n-th ordinary basis switching directive
                                        ! covers this pair. This is the normally expected
                                        ! case.
                                        go to 100
                                    else
                                        ! The first species is set to be switched with
                                        ! a different second species.
                                        ! Calling sequence substitutions:
                                        !   jlend for jlen2
                                        !   uspd56 for uspn56
                                        call fmspnx(jlend,uobsw(2,n),uspd56)

                                        write (noutpt,1010) uspb56(1:jlenb),uspa56(1:jlena),uspd56(1:jlend)
                                        write (nttyo,1010) uspb56(1:jlenb),uspa56(1:jlena),uspd56(1:jlend)
1010 format(/' * Error - (EQ6/intmat) The species ',a,' must be in',/7x,'the active basis set for',' the system described on the input',/7x,'file. This implies it should replace the',' "data file" basis',/7x,'species ',a,'. However, an explicit replacement by',' the species',/7x,a,' is called for by an',' ordinary basis switch specified',/7x,'on the input file.')

                                        nerr = nerr + 1
                                        go to 100
                                    end if
                                else if (uobsw(2,n)(1:48) .eq. uspec(ns2)(1:48)) then
                                    ! The second species is set to be switched with a
                                    ! different first species.
                                    !   Calling sequence substitutions:
                                    !     jlenc for jlen2
                                    !     uspc56 for uspn56
                                    call fmspnx(jlenc,uobsw(1,n),uspc56)

                                    write (noutpt,1020) uspb56(1:jlenb),uspa56(1:jlena),uspc56(1:jlenc)
                                    write (nttyo,1020) uspb56(1:jlenb),uspa56(1:jlena),uspc56(1:jlenc)
1020 format(/' * Error - (EQ6/intmat) The species ',a,' must be in',/7x,'the active basis set for',' the system described on the input',/7x,'file. This implies it should replace the',' "data file" basis',/7x,'species ',a,'. However, an explicit replacing of',' the species',/7x,a,' is called for by an',' ordinary basis switch specified',/7x,'on the input file.')

                                    nerr = nerr + 1
                                    go to 100
                                end if
                            end do

                            ! The implied ordinary basis switch does not match
                            ! any specified switch. Create one.
                            write (noutpt,1030) uspa56(1:jlena),uspb56(1:jlenb)
                            write (nttyo,1030) uspa56(1:jlena),uspb56(1:jlenb)
1030 format(/' * Note - (EQ6/intmat) An implied',' ordinary basis switch',/7x,'replacing ',a,' with ',a,' was found.',/7x,'An explicit switch',' will be constructed.')

                            nn = nobswt + 1

                            if (nn .gt. nbtmax) then
                                write (ux8,'(i8)') nbtmax
                                call lejust(ux8)
                                j2 = ilnobl(ux8)
                                write (noutpt,1040) ux8(1:j2)
                                write (nttyo,1040) ux8(1:j2)
1040 format(/' * Error - (EQ6/intmat) Cannot construct',' an explicit ordinary',/7x,'basis switch because',' this would exceed the limit of ',a,' switches',/7x,'(which is equal to the dimensioned limit',' for basis species).')

                                nerr = nerr + 1
                            end if

                            nobswt = nn
                            uobsw(1,nn) = uspec(ns1)
                            uobsw(2,nn) = uspec(ns2)
                            go to 100
                        end if
                    end do
                end if
            end do
        end if

        ! Cannot identify the active basis species specified in the
        ! uzveci array on the input file.
        write (noutpt,1060) uspb56(1:jlenb)
        write (nttyo,1060) uspb56(1:jlenb)
1060 format(/' * Error - (EQ6/intmat) The species ',a,' is in the',/7x,'active basis set for the system described on the input',' file,',/7x,"but it isn't on the data file and wasn't created",' by',/7x,'an input file directive.')

        nerr = nerr + 1

100 continue
    end do

    kbt = kcol

    ! Pure phases.
    km1 = kcol + 1

    if (kmts .lt. km1s) then
        go to 120
    end if

    qloffg = .false.

    do krow = km1s,kmts
        do np = 1,npt
            nr1 = ncmpr(1,np)
            nr2 = ncmpr(2,np)
            nt = nr2 - nr1 + 1

            if (nt .eq. 1) then
                ns = nr1

                if (uzveci(krow)(1:48) .eq. uspec(ns)(1:48)) then
                    kcol = kcol + 1
                    iindx1(kcol) = ns
                    ipndx1(kcol) = np
                    uzvec1(kcol) = uspec(ns)
                    lxx = zvclgi(krow)
                    zvclg1(kcol) = lxx
                    zvec1(kcol) = texp(lxx)
                    go to 110
                end if
            end if
        end do

        if (uzveci(krow)(1:5) .eq. ufixf(1:5)) then
            ux24 = uzveci(krow)(6:24)
            j2 = ilnobl(ux24)
            write (noutpt,1070) ux24(1:j2)
            write (nttyo,1070) ux24(1:j2)
1070 format(/' * Note - (EQ6/intmat) A left over fixed fugacity',' phase for ',a,/7x,'is in the system described on the',' input file. It will be purged.')

            qloffg = .true.
        else
            call fmspnm(jlen2,uzveci(krow),uspn56)
            write (noutpt,1080) uspn56(1:jlen2)
            write (nttyo,1080) uspn56(1:jlen2)
1080 format(/' * Error - (EQ6/intmat) The species ',a,/7x,'is in the system described on the input file, but it'," isn't",/7x,"on the data file and wasn't created",' following an input file directive.')

            nerr = nerr + 1
        end if

110 continue
    end do

120 continue
    kmt = kcol

    ! Non-aqueous solutions.
    kx1 = kcol + 1

    do krow = kx1s,kxts
        do np = 1,npt
            nr1 = ncmpr(1,np)
            nr2 = ncmpr(2,np)
            nt = nr2 - nr1 + 1

            if (nt.gt.1 .and.     uzveci(krow)(25:48).eq.uphase(np)(1:24)) then
                go to 130
            end if
        end do

        call fmspnm(jlen2,uzveci(krow),uspn56)
        write (noutpt,1080) uspn56(1:jlen2)
        write (nttyo,1080) uspn56(1:jlen2)
        nerr = nerr + 1
        go to 140

130 continue
        do ns = nr1,nr2
            if (uzveci(krow)(1:24) .eq. uspec(ns)(1:24)) then
                kcol = kcol + 1
                iindx1(kcol) = ns
                ipndx1(kcol) = np
                uzvec1(kcol) = uspec(ns)
                lxx = zvclgi(krow)
                zvclg1(kcol) = lxx
                zvec1(kcol) = texp(lxx)
                go to 140
            end if
        end do

        call fmspnm(jlen2,uzveci(krow),uspn56)
        write (noutpt,1080) uspn56(1:jlen2)
        write (nttyo,1080) uspn56(1:jlen2)
        nerr = nerr + 1

140 continue
    end do

    kxt = kcol
    kdim = kxt

    kwater = 0

    if (narn1 .ne. 0) then
        do kcol = 1,kbt
            nb = iindx1(kcol)
            ns = nbasp(nb)

            if (ns .eq. narn1) then
                kwater = kcol
                go to 160
            end if
        end do

        call fmspnx(jlen2,uspec(narn1),uspn56)
        write (noutpt,1200) uspn56(1:jlen2)
        write (nttyo,1200) uspn56(1:jlen2)
1200 format(/' * Error - (EQ6/intmat) The species ',a," isn't",/7x,'on the input file, as is required.')

        nerr = nerr + 1
    end if

160 continue

    khydr = 0

    if (nhydr .ne. 0) then
        do kcol = 1,kbt
            nb = iindx1(kcol)
            ns = nbasp(nb)

            if (ns .eq. nhydr) then
                khydr = kcol
                go to 170
            end if
        end do

        if (nhydx .eq. 0) then
            call fmspnx(jlen2,uspec(nhydr),uspn56)
            write (noutpt,1200) uspn56(1:jlen2)
            write (nttyo,1200) uspn56(1:jlen2)
            nerr = nerr + 1
        end if
    end if

170 continue

    khydx = 0

    if (nhydx .ne. 0) then
        do kcol = 1,kbt
            nb = iindx1(kcol)
            ns = nbasp(nb)

            if (ns .eq. nhydx) then
                khydx = kcol
                go to 180
            end if
        end do

        if (nhydr .eq. 0) then
            call fmspnx(jlen2,uspec(nhydx),uspn56)
            write (noutpt,1200) uspn56(1:jlen2)
            write (nttyo,1200) uspn56(1:jlen2)
            nerr = nerr + 1
        end if
    end if

180 continue

    ko2gaq = 0

    if (no2gaq .ne. 0) then
        do kcol = 1,kbt
            nb = iindx1(kcol)
            ns = nbasp(nb)

            if (ns .eq. no2gaq) then
                ko2gaq = kcol
                go to 190
            end if
        end do

        if (nelect .eq. 0) then
            call fmspnm(jlen2,uspec(no2gaq),uspn56)
            write (noutpt,1200) uspn56(1:jlen2)
            write (nttyo,1200) uspn56(1:jlen2)
            nerr = nerr + 1
        end if
    end if

190 continue

    kelect = 0

    if (nelect .ne. 0) then
        do kcol = 1,kbt
            nb = iindx1(kcol)
            ns = nbasp(nb)

            if (ns .eq. nelect) then
                kelect = kcol
                go to 200
            end if
        end do

        if (no2gaq .eq. 0) then
            call fmspnm(jlen2,uspec(nelect),uspn56)
            write (noutpt,1200) uspn56(1:jlen2)
            write (nttyo,1200) uspn56(1:jlen2)
            nerr = nerr + 1
        end if
    end if

200 continue

    if (nerr .gt. 0) then
        stop
    end if
end subroutine intmat