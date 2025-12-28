subroutine intbs3(covali,ier,jflag,jflgi,narn1a,narn2a,nbaspd,nbtd,nbti,nbtmax,ndrsrd,ndecsp,noutpt,nrdxsp,nsta,nstmax,nttyo,uspeca,uspeci)
    !! This subroutine interprets the basis species listed on the input
    !! file. It sets up the jflag arrays.
    !! This subroutine is called by:
    !!   EQ3NR/eq3nr.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: nbtmax
    integer :: nstmax

    integer :: jflag(nstmax)
    integer :: jflgi(nbtmax)
    integer :: nbaspd(nbtmax)
    integer :: ndecsp(nbtmax)
    integer :: ndrsrd(2,nstmax)

    integer :: ier
    integer :: narn1a
    integer :: narn2a
    integer :: nbtd
    integer :: nbti
    integer :: noutpt
    integer :: nrdxsp
    integer :: nsta
    integer :: nt
    integer :: nttyo

    character(len=48) :: uspeca(nstmax)
    character(len=48) :: uspeci(nbtmax)

    real(kind=8) :: covali(nbtmax)

    ! Local variable declarations.
    integer :: j2
    integer :: nb
    integer :: nbi
    integer :: ns

    integer :: ilnobl
    integer :: nbasis

    ier = 0

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

    jflag(narn1a) = 0

    if (nrdxsp .gt. 0) then
        jflag(nrdxsp) = 0
    end if

    do nbi = 1,nbti
        j2 = ilnobl(uspeci(nbi)(1:24))

        do ns = 1,nsta
            if (uspeca(ns)(1:24) .eq. uspeci(nbi)(1:24)) then
                if (ndecsp(nbi) .gt. 0) then
                    write (noutpt,1000) uspeci(nbi)(1:j2)
                    write (nttyo,1000) uspeci(nbi)(1:j2)
1000 format(/' * Error - (EQ3NR/intbs3) More than one',' constraint was specified',/7x,'on the input file for',/7x,'the species ',a,'.')

                    ier = ier + 1
                else
                    if (jflgi(nbi).le.15 .and. covali(nbi).le.0.) then
                        jflgi(nbi) = -1
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
                            write (noutpt,1005) nbtmax,uspeci(nbi)(1:j2)
                            write (nttyo,1005) nbtmax,uspeci(nbi)(1:j2)
1005 format(/' * Error - (EQ3NR/intbs3) The maximum ',i7,' basis species have been exceeded',/7x,'while',' trying to load ',a,'. Increase the dimensioning',/7x,'parameter nbtpar.')

                            ier = ier + 1
                            go to 115
                        end if

                        nb = nbtd
                        nbaspd(nb) = ns
                    end if

                    ndecsp(nbi) = nb
                end if

                go to 115
            end if
        end do

        write (noutpt,1010) uspeci(nbi)(1:j2)
        write (nttyo,1010) uspeci(nbi)(1:j2)
1010 format(/' * Error - (EQ3NR/intbs3) The species ',a,' is',/7x,"on the input file, but it isn't on the supporting",' data file.')

        ier = ier + 1

115 continue
    end do
end subroutine intbs3