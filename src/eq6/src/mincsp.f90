subroutine mincsp(cdrsd,jpflag,jsflag,nbaspd,nbtd,nbtmax,ncmpra,ndrsd,ndrsmx,ndrsrd,nmrn1a,nmrn2a,noutpt,npta,nptmax,nstmax,nttyo,nxopex,nxopmx,nxopt,nxpemx,uspeca,uxcat,uxopex,uxopt)
    !! This subroutine executes the mineral subset-selection suppression
    !! (nxopt) options.
    !! This subroutine is called by:
    !!   EQ6/eq6.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: nbtmax
    integer :: ndrsmx
    integer :: nptmax
    integer :: nstmax
    integer :: nxopmx
    integer :: nxpemx

    integer :: jpflag(nptmax)
    integer :: jsflag(nstmax)
    integer :: ncmpra(2,nptmax)
    integer :: nbaspd(nbtmax)
    integer :: ndrsd(ndrsmx)
    integer :: ndrsrd(2,nstmax)

    integer :: nbtd
    integer :: nmrn1a
    integer :: nmrn2a
    integer :: noutpt
    integer :: npta
    integer :: nttyo
    integer :: nxopex
    integer :: nxopt

    character(len=48) :: uspeca(nstmax)
    character(len=24) :: uxcat(nxopmx)
    character(len=24) :: uxopex(nxpemx)
    character(len=8) :: uxopt(nxopmx)

    real(kind=8) :: cdrsd(ndrsmx)

    ! Local variable declarations.
    integer :: i
    integer :: j
    integer :: jlen
    integer :: j2
    integer :: j3
    integer :: n
    integer :: nb
    integer :: nerr
    integer :: np
    integer :: np1m1
    integer :: ns
    integer :: nse

    integer :: ilnobl
    integer :: nphase

    character(len=56) :: uspn56
    character(len=24) :: ux24
    character(len=8) :: ux8

    real(kind=8) :: cx

    real(kind=8) :: coefdr

    if (nxopt .le. 0) then
        go to 999
    end if

    nerr = 0

    i = 0

    do n = 1,nxopt
        ux8 = uxopt(n)
        call locase(ux8)

        if (ux8(1:4) .eq. 'all ') then
            i = i + 1
        end if
    end do

    if (i .gt. 1) then
        write (noutpt,1010)
        write (nttyo,1010)
1010 format(/' * Error - (EQ6/mincsp) The "all" mineral subset-','selection suppression',/7x,'option has been specified more',' than once.')

        nerr = nerr + 1
    end if

    do n = 1,nxopt
        ux8 = uxopt(n)
        call locase(ux8)

        if (ux8(1:4) .eq. 'all ') then
            j2 = ilnobl(uxcat(n))

            if (j2 .gt. 0) then
                write (noutpt,1020) uxcat(n)(1:j2)
                write (nttyo,1020) uxcat(n)(1:j2)
1020 format(/' * Error - (EQ6/mincsp) A category may not be',' specified for an',/7x,'"all" mineral subset-selection',' suppression option. Here "',a,'"',/7x,'is specified',' as defining a category.')

                nerr = nerr + 1
            end if
        end if
    end do

    do n = 1,nxopt
        ux8 = uxopt(n)
        call locase(ux8)

        if (ux8(1:8).eq.'alwith ' .or. ux8(1:8).eq.'allwith') then
            ux24 = uxcat(n)
            j2 = ilnobl(ux24)

            if (j2 .le. 0) then
                write (noutpt,1030)
                write (nttyo,1030)
1030 format(/' * Error - (EQ6/mincsp) A category must be',' specified for each',/7x,'"allwith" mineral subset-','selection suppression option. Have an',/7x,'instance',' where this is lacking.')

                nerr = nerr + 1
            else
                do j = n + 1,nxopt
                    ux8 = uxopt(j)
                    call locase(ux8)

                    if (ux8(1:8).eq.'alwith ' .or. ux8(1:8).eq.'allwith') then
                        if (uxcat(j)(1:24) .eq. ux24(1:24)) then
                            write (noutpt,1040) ux24(1:j2)
                            write (nttyo,1040) ux24(1:j2)
1040 format(/' * Error - (EQ6/mincsp) The string "',a,'"',' is repeated in defining',/7x,'a category for an',' "allwith" mineral subset-selection suppression',/7x,'option. It may only be so used once.')

                            nerr = nerr + 1
                        end if
                    end if
                end do
            end if
        end if
    end do

    if (nerr .gt. 0) then
        stop
    end if

    ns = nmrn1a

    ! Calling sequence substitutions:
    !   ncmpra for ncmpr
    !   npta for npt
    np1m1 = nphase(ncmpra,npta,nptmax,ns) - 1

    do n = 1,nxopt
        ux8 = uxopt(n)
        call locase(ux8)

        if (ux8(1:4) .eq. 'all ') then
            write (noutpt,1110)
            write (nttyo,1110)
1110 format(/' * Note - (EQ6/mincsp) Executing the "all" mineral',' subset-selection',/7x,'suppression option.')

            if (nxopex .gt. 0) then
                write (noutpt,1120)
                write (nttyo,1120)
1120 format(/9x,'Minerals that are exceptions are:',/)
            end if

            np = np1m1

            do ns = nmrn1a,nmrn2a
                np = np + 1

                if (jsflag(ns) .le. 0) then
                    do i = 1,nxopex
                        if (uspeca(ns)(1:24) .eq. uxopex(i)(1:24)) then
                            ! Calling sequence substitutions:
                            !   uspeca(ns) for unam48
                            call fmspnm(jlen,uspeca(ns),uspn56)
                            write (noutpt,1130) uspn56(1:jlen)
                            write (nttyo,1130) uspn56(1:jlen)
1130 format(11x,a)

                            go to 200
                        end if
                    end do

                    jsflag(ns) = 1
                    jpflag(np) = 1
                end if

200 continue
            end do
        else if (ux8(1:8).eq.'alwith ' .or. ux8(1:8).eq.'allwith') then
            do nb = 1,nbtd
                nse = nbaspd(nb)

                if (uspeca(nse)(1:24) .eq. uxcat(n)(1:24)) then
                    j3 = ilnobl(uxcat(n))
                    write (noutpt,1140) uxcat(n)(1:j3)
                    write (nttyo,1140) uxcat(n)(1:j3)
1140 format(/' * Note - (EQ6/mincsp) Executing the "alwith ',a,'" mineral',/7x,'subset-selection suppression option.')

                    if (nxopex .gt. 0) then
                        write (noutpt,1120)
                        write (nttyo,1120)
                    end if

                    np = np1m1

                    do ns = nmrn1a,nmrn2a
                        np = np + 1

                        if (jsflag(ns) .le. 0) then
                            ! Calling sequence substitutions:
                            !   cdrsd for cdrs
                            !   ndrsd for ndrs
                            !   ndrsrd for ndrsr
                            cx =  coefdr(cdrsd,ndrsd,ndrsmx,ndrsrd,nse,ns,nstmax)

                            if (cx .ne. 0.) then
                                if (nxopex .gt. 0) then
                                    do i = 1,nxopex
                                        if (uspeca(ns)(1:24) .eq. uxopex(i)(1:24)) then
                                            ! Calling sequence substitutions:
                                            !   uspeca(ns) for unam48
                                            call fmspnm(jlen,uspeca(ns),uspn56)
                                            write (noutpt,1130) uspn56(1:jlen)
                                            write (nttyo,1130) uspn56(1:jlen)
                                            go to 210
                                        end if
                                    end do
                                end if

                                jsflag(ns) = 1
                                jpflag(np) = 1
                            end if
                        end if

210 continue
                    end do

                    go to 220
                end if
            end do

            j3 = ilnobl(uxcat(n))
            write (noutpt,1150) uxcat(n)(1:j3)
            write (nttyo,1150) uxcat(n)(1:j3)
1150 format(" * Error - (EQ6/mincsp) Don't recognize",' "',a,'" as',/7x,'an argument to the "alwith" mineral',' subset-selection',/7x,'suppression option.')

            nerr = nerr + 1

220 continue
        else
            j2 = ilnobl(uxopt(n))
            j3 = ilnobl(uxcat(n))
            write (noutpt,1160) uxopt(n)(1:j2),uxcat(n)(1:j3)
            write (nttyo,1160) uxopt(n)(1:j2),uxcat(n)(1:j3)
1160 format(/" * Error - (EQ6/mincsp) Don't recognize",' "',a,1x,a,'"',/7x,'as a mineral subset-selection ',' suppression option.')

            nerr = nerr + 1
        end if
    end do

    write (noutpt,1170)
    write (nttyo,1170)
1170 format(1x)

    if (nerr .gt. 0) then
        stop
    end if

    ! r * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
999 continue
end subroutine mincsp