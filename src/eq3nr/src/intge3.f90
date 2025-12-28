subroutine intge3(cgexp,cgexpi,ier,iern1,iern2,ietmax,jern1,jern2,jetmax,jgext,jgexti,net,neti,netmax,ngexpi,ngexti,noutpt,nptmax,nstmax,nttyo,ugexj,ugexji,ugexp,ugexpi,ugexsi,uphase,uspec,xbar,xbarlg,xgexsi)
    !! This subroutine interprets any concentrations and compositions
    !! of generic exchange phases read from the EQ3NR input file. Here
    !! concentration means moles of exchanger phase per kg H2O, and
    !! composition means a set of mole fractions of the component
    !! species on each site of the exchanger phase.
    !! This subroutine is called by:
    !!   EQ3NR/eq3nr.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: ietmax
    integer :: jetmax
    integer :: netmax
    integer :: nptmax
    integer :: nstmax

    integer :: noutpt
    integer :: nttyo

    integer :: jern1(jetmax,netmax)
    integer :: jern2(jetmax,netmax)
    integer :: jgext(netmax)
    integer :: jgexti(netmax)
    integer :: ngexpi(netmax)
    integer :: ngexti(jetmax,netmax)
    integer :: ier
    integer :: iern1
    integer :: iern2
    integer :: net
    integer :: neti

    character(len=48) :: uspec(nstmax)
    character(len=24) :: ugexp(netmax)
    character(len=24) :: ugexpi(netmax)
    character(len=24) :: ugexsi(ietmax,jetmax,netmax)
    character(len=24) :: uphase(nptmax)
    character(len=8) :: ugexj(jetmax,netmax)
    character(len=8) :: ugexji(jetmax,netmax)

    real(kind=8) :: cgexp(netmax)
    real(kind=8) :: cgexpi(netmax)
    real(kind=8) :: xbar(nstmax)
    real(kind=8) :: xbarlg(nstmax)
    real(kind=8) :: xgexsi(ietmax,jetmax,netmax)

    ! Local variable declarations.
    integer :: iei
    integer :: je
    integer :: jeei
    integer :: jei
    integer :: jj
    integer :: jj2
    integer :: jj3
    integer :: j2
    integer :: j3
    integer :: j4
    integer :: j5
    integer :: j6
    integer :: nd2
    integer :: ne
    integer :: neei
    integer :: nei
    integer :: np
    integer :: nr1
    integer :: nr2
    integer :: ns

    integer :: ilnobl

    character(len=24) :: ugexpd
    character(len=24) :: ustr
    character(len=24) :: uxp
    character(len=24) :: uxs
    character(len=8) :: ugexjd
    character(len=8) :: uxj
    character(len=8) :: ux8a
    character(len=8) :: ux8b

    real(kind=8) :: xsum
    real(kind=8) :: xx

    real(kind=8) :: tlg

    ier = 0

    if (neti .le. 0) then
        go to 999
    end if

    if (net .le. 0) then
        call lejust(ugexpi(1))
        j4 = ilnobl(ugexpi(1))

        if (neti .eq. 1) then
            nei = 1
            write (noutpt,1000) ugexpi(nei)(1:j4)
            write (nttyo,1000) ugexpi(nei)(1:j4)
1000 format(/' * Error - (EQ3NR/intge3) No generic exchanger',' phases were declared',/7x,'on the input file, but a ',' concentration and composition are',/7x,'specified for',' the generic exchanger phase ',a,'.',/7x,'Specifying a',' concentration and  composition does not',/7x,'automatically create such a phase.',/7x,'Other input',' file directives are required.')

            ier = ier + 1
        else
            nei = 1
            write (noutpt,1010) neti,ugexpi(nei)(1:j4)
            write (nttyo,1010) neti,ugexpi(nei)(1:j4)
1010 format(/' * Error - (EQ3NR/intge3) No generic exchanger',' phases were declared',/7x,'on the input file, but',' concentrations and compositions are specified',/7x,'for ',i3,' such phases, including ',a,'.',/7x,'Specifying a',' concentration and  composition does not',/7x,'automatically create such a phase.',/7x,'Other input',' file directives are required.')

            ier = ier + 1
        end if
    end if

    if (ier .gt. 0) then
        go to 999
    end if

    ! Left-justify all input exchanger phase names and site names. Make
    ! sure that no such names are blank. Assign defaults if necessary.
    do nei = 1,neti
        call lejust(ugexpi(nei))
        j4 = ilnobl(ugexpi(nei))

        if (j4 .le. 0) then
            ! Calling sequence substitutions:
            !   nei for ne
            call adgexp(nei,noutpt,nttyo,ugexpd)
            ugexpi(nei) = ugexpd
        end if

        do jei = 1,jgexti(nei)
            call lejust(ugexji(jei,nei))
            j3 = ilnobl(ugexji(jei,nei))

            if (j3 .le. 0) then
                ! Calling sequence substitutions:
                !   jei for je
                call adgexj(jei,noutpt,nttyo,ugexjd)
                ugexji(jei,nei) = ugexjd
            end if
        end do
    end do

    ! Make sure that the names of exchanger phases for which
    ! concentrations and compositions are specified are unique.
    ! Make sure that the names of the sites of these phases are
    ! unique for each such phase.
    do nei = 1,neti
        j4 = ilnobl(ugexpi(nei))

        if (nei .gt. 1) then
            do neei = 1,nei - 1
                if (ugexpi(nei)(1:24) .eq. ugexpi(neei)(1:24)) then
                    write (noutpt,1030) ugexpi(nei)(1:j4),nei,neei
                    write (nttyo,1030) ugexp(nei)(1:j4),nei,neei
1030 format (/' * Error - (EQ3NR/intge3) "',a,'" is the name',' specified',/7x,'for exchanger phases ',i3,' and ',i3,' on the list of such phases',/7x,'whose concentrations',' and compositions are given on the input file.',/7x,'Only one concentration and composition may be',' entered per',/7x,'exchanger phase.')

                    ier = ier + 1
                end if
            end do
        end if

        do jei = 1,jgexti(nei)
            j3 = ilnobl(ugexji(jei,nei))

            if (jei .gt. 1) then
                do jeei = 1,jei - 1
                    if (ugexji(jei,nei)(1:8) .eq. ugexji(jeei,nei)(1:8)) then
                        write (noutpt,1040) ugexji(jei,nei)(1:j3),jei,jeei,ugexpi(nei)(1:j4)
                        write (nttyo,1040) ugexji(jei,nei)(1:j3),jei,jeei,ugexpi(nei)(1:j4)
1040 format (/' * Error - (EQ3NR/intge3) "',a,'" is the',' name specified',/7x,'for sites ',i3,' and ',i3,' of',' exchanger phase ',a,/7x,', whose composition is',' given on the input file.',/7x,'A site name must be',' unique for a given exchanger phase.')

                        ier = ier + 1
                    end if
                end do
            end if
        end do
    end do

    if (ier .gt. 0) then
        go to 999
    end if

    ! Interpret the input data.
    do nei = 1,neti
        uxp = ugexpi(nei)
        j4 = ilnobl(uxp)
        ne = 0

        do np = iern1,iern2
            if (uxp(1:24) .eq. uphase(np)(1:24)) then
                go to 100
            end if
        end do

        write (noutpt,1070) uxp(1:j4)
        write (nttyo,1070) uxp(1:j4)
1070 format(/' * Error - (EQ3NR/intge3) The generic exchanger',' phase',/7x,a,', whose composition is specified on the input',/7x,"file, doesn't match the name of any such phase declared",/7x,'on the input file.')

        ier = ier + 1
        go to 200

100 continue
        ne = np - iern1 + 1
        ngexpi(nei) = np
        cgexp(ne) = cgexpi(nei)

        if (jgexti(nei) .ne. jgext(ne)) then
            write (ux8a,'(i5)') jgexti(nei)
            call lejust(ux8a)
            j5 = ilnobl(ux8a)
            write (ux8b,'(i5)') jgext(ne)
            call lejust(ux8b)
            j6 = ilnobl(ux8b)
            write (noutpt,1080) uxp(1:j4),ux8a(1:j5),ux8b(1:j6)
            write (nttyo,1080) uxp(1:j4),ux8a(1:j5),ux8b(1:j6)
1080 format(/' * Error - (EQ3NR/intge3) The generic exchanger',' phase ',a,/7x,'has a specified composition including ',a,' exchange site(s), but',/7x,'this phase actually has ',a,' exchange site(s). The number of such',/7x,'sites much match.')

            ier = ier + 1
            go to 200
        end if

        do jei = 1,jgexti(nei)
            uxj = ugexji(jei,nei)
            j3 = ilnobl(uxj)

            do je = 1,jgext(ne)
                if (uxj(1:8) .eq. ugexj(je,ne)(1:8)) then
                    go to 110
                end if
            end do

            write (noutpt,1110) uxj(1:j3),uxp(1:j4)
            write (nttyo,1110) uxj(1:j3),uxp(1:j4)
1110 format(/' * Error - (EQ3NR/intge3) Site ',a,' of generic',' exchanger phase',/7x,a,', whose composition is specified',' on',/7x,"the input file, doesn't match the name of any",' site',/7x,'of this phase as declared on the input file.')

            ier = ier + 1
            go to 190

110 continue
            xsum = 0.

            do iei = 1,ngexti(jei,nei)
                uxs = ugexsi(iei,jei,nei)
                j2 = ilnobl(uxs)

                ! Compose the exchanger species name, by combining the
                ! name of the exchangeable species (the aqueous species)
                ! with the site name. The algorithm for this follows that
                ! used in EQLIB/intexi.f, except that the phase part of the
                ! name is not employed here.
                jj2 = j2
                jj3 = j3
                nd2 = 0
250 continue
                jj = jj2 + jj3 + 1

                if (jj .gt. 24) then
                    if (nd2 .lt. 2) then
                        jj2 = jj2 - 1
                        nd2 = nd2 + 1
                        go to 250
                    else
                        jj3 = jj3 - 1
                        nd2 = 0
                        go to 250
                    end if
                end if

                ustr = ' '
                ustr(1:jj2) = ugexsi(iei,jei,nei)(1:jj2)
                ustr(jj2 + 1:jj2 + 1) = ' '
                ustr(jj2 + 2:jj) = ugexji(jei,nei)(1:jj3)

                nr1 = jern1(je,ne)
                nr2 = jern2(je,ne)

                do ns = nr1,nr2
                    if (ustr(1:24) .eq. uspec(ns)(1:24)) then
                        go to 130
                    end if
                end do

                write (noutpt,1120) uxs(1:j2),uxj(1:j3),uxp(1:j4)
                write (nttyo,1120) uxs(1:j2),uxj(1:j3),uxp(1:j4)
1120 format(/' * Error - (EQ3NR/intge3) The species ',a,' for',/7x,'site ',a,' of generic exchanger phase ',a,',',/7x,'whose composition is specified on the input file,'," doesn't match",/7x,'the name of any species for this',' site of this phase',/7x,'as declared on the input file.')

                ier = ier + 1
                go to 180

130 continue
                xx = xgexsi(iei,jei,nei)

                if (xx .ge. 0.) then
                    xsum = xsum + xx
                    xbar(ns) = xx
                    xbarlg(ns) = tlg(xx)
                else
                    write (noutpt,1130) uxs(1:j2),uxj(1:j3),uxp(1:j4),xx
                    write (nttyo,1130) uxs(1:j2),uxj(1:j3),uxp(1:j4),xx
1130 format(/' * Error - (EQ3NR/intge3) The species ',a,' for ',/7x,'site ',a,' of generic exchanger phase ',a,/7x,'has an illegal negative mole fraction of ',1pe12.5,' specified',/7x,'on the input file.')

                    ier = ier + 1
                end if

180 continue
            end do

            ! Normalize the mole fractions on the current site.
            nr1 = jern1(je,ne)
            nr2 = jern2(je,ne)

            if (xsum .gt. 0.) then
                do ns = nr1,nr2
                    xx = xbar(ns)/xsum
                    xbar(ns) = xx
                    xbarlg(ns) = tlg(xx)
                end do
            end if

190 continue
        end do

200 continue
    end do

999 continue
end subroutine intge3