subroutine intbs6(jflag,jflgi,kmax,narn1a,narn2a,nbaspd,nbtd,nbti,nbtmax,ndrsrd,ndecsp,noutpt,nsta,nstmax,nttyo,uspeca,ubmtbi)
    !! This subroutine interprets the data file basis species read from
    !! the input file. "Data file" species to be created according to
    !! directives read from the input file are ignored. This subroutine
    !! also sets up the jflag array.
    !! This subroutine is called by:
    !!   EQ6/eq6.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: kmax
    integer :: nbtmax
    integer :: nstmax

    integer :: jflag(nstmax)
    integer :: jflgi(nbtmax)
    integer :: nbaspd(nbtmax)
    integer :: ndecsp(nbtmax)
    integer :: ndrsrd(2,nstmax)

    integer :: narn1a
    integer :: narn2a
    integer :: nbtd
    integer :: nbti
    integer :: noutpt
    integer :: nsta
    integer :: nttyo

    character(len=56) :: uspn56
    character(len=48) :: uspeca(nstmax)
    character(len=48) :: ubmtbi(nbtmax)

    ! Local variable declarations.
    integer :: jlen
    integer :: j2
    integer :: nb
    integer :: nbi
    integer :: nerr
    integer :: ns
    integer :: nt

    integer :: ilnobl
    integer :: nbasis

    character(len=8) :: ux8

    nerr = 0

    do ns = 1,nsta
        jflag(ns) = 0
    end do

    do ns = narn1a,narn2a
        jflag(ns) = 30
    end do

    do nb = 1,nbtd
        ns = nbaspd(nb)
        nt = ndrsrd(2,ns) - ndrsrd(1,ns) + 1

        if (nt .lt. 2) then
            jflag(ns) = -1
        end if
    end do

    do nbi = 1,nbti
        do ns = 1,nsta
            if (uspeca(ns)(1:48) .eq. ubmtbi(nbi)(1:48)) then
                if (jflgi(nbi).ne.0 .and. jflgi(nbi).ne.30) then
                    ! Calling sequence substitutions:
                    !   ubmtbi(nbi) for unam48
                    call fmspnx(jlen,ubmtbi(nbi),uspn56)
                    write (ux8,'(i5)') jflgi(nbi)
                    call lejust(ux8)
                    j2 = ilnobl(ux8)
                    write (noutpt,1000) uspn56(1:jlen),ux8(1:j2)
                    write (nttyo,1000) uspn56(1:jlen),ux8(1:j2)
1000 format(/' * Error - (EQ6/intbs6) The species ',a,/7x,'has an has illegal jflag value of ',a,' on the',' input file',/7x,'(the only legal values are 0 and 30).')

                    nerr = nerr + 1
                end if

                jflag(ns) = jflgi(nbi)

                ! Calling sequence substitutions:
                !   nbaspd for nbasp
                !   nbtd for nbt
                nb = nbasis(nbaspd,nbtd,nbtmax,ns)

                if (nb .eq. 0) then
                    ! Add a species to the working basis set.
                    nbtd = nbtd + 1

                    if (nbtd .gt. nbtmax) then
                        ! Calling sequence substitutions:
                        !   ubmtbi(nbi) for unam48
                        call fmspnx(jlen,ubmtbi(nbi),uspn56)
                        write (ux8,'(i5)') nbtmax
                        call lejust(ux8)
                        j2 = ilnobl(ux8)
                        write (noutpt,1005) ux8(1:j2),uspn56(1:jlen)
                        write (nttyo,1005) ux8(1:j2),uspn56(1:jlen)
1005 format(/' * Error - (EQ6/intbs6) The maximum ',a,/7x,'basis species have been exceeded in',' interpreting the',/7x,'input file while processing',' data file basis species',/7x,a,'. Increase the',' dimensioning parameter nbtpar.')

                        nerr = nerr + 1
                        go to 115
                    end if

                    nb = nbtd
                    nbaspd(nb) = ns
                end if

                ndecsp(nbi) = nb
                go to 115
            end if
        end do

        ! Calling sequence substitutions:
        !   ubmtbi(nbi) for unam48
        call fmspnx(jlen,ubmtbi(nbi),uspn56)
        write (noutpt,1010) uspn56(1:jlen)
        write (nttyo,1010) uspn56(1:jlen)
1010 format(/' * Note - (EQ6/intbs6) The species "',a,'"',/7x,'appears as a basis species on the input file, but it'," wasn't",/7x,'read from the data file. If it is a species',' to be created by the code,',/7x,'such as a generic ion',' exchanger species, there is no problem.')

115 continue
    end do

    if (nerr .gt. 0) then
        stop
    end if
end subroutine intbs6