subroutine flgset(axlksd,iopt,jflag,jpflag,jsflag,kxmod,narn1a,narn2a,narxmx,nbaspd,nbtd,nbtmax,ncmpra,ncta,ndrsd,ndrsmx,ndrsrd,noptmx,noutpt,npta,nptmax,nrdxsp,nsta,nstmax,ntpr,ntprmx,nttyo,nxmdmx,nxmod,uphasa,uptypa,uspeca,uxmod)
    !! This subroutine sets up the status arrays jpflag and jsflag.
    !! These flags denote the statuses, respectively, of phases
    !! and species. The relevant values and their meanings are
    !! as follows:
    !!    jpflag:
    !!      = 0   The species is present or potentially present
    !!      = 1   The physical presence of the phase in the model
    !!              is suppressed; associated variables, such as a
    !!              reaction affinity, may be calculated
    !!      = 2   The phase is completely ignored by the code
    !!    jsflag:
    !!      = 0   The species is present or potentially present
    !!      = 1   The physical presence of the species in the model
    !!              is suppressed; associated variables, such as a
    !!              reaction affinity, may be calculated
    !!      = 2   The species is completely ignored by the code
    !! This subroutine is called by:
    !!   EQ3NR/eq3nr.f
    !!   EQ6/eq6.f
    !! Principal input:
    !! Principal output:
    !!   jpflag = phase status flag array
    !!   jsflag = species status flag array
    implicit none

    ! Calling sequence variable declarations.
    integer :: narxmx
    integer :: nbtmax
    integer :: ndrsmx
    integer :: noptmx
    integer :: nptmax
    integer :: nstmax
    integer :: ntprmx
    integer :: nxmdmx

    integer :: noutpt
    integer :: nttyo

    integer :: iopt(noptmx)
    integer :: jflag(nstmax)
    integer :: jpflag(nptmax)
    integer :: jsflag(nstmax)
    integer :: kxmod(nxmdmx)
    integer :: nbaspd(nbtmax)
    integer :: ncmpra(2,nptmax)
    integer :: ndrsd(ndrsmx)
    integer :: ndrsrd(2,nstmax)

    integer :: narn1a
    integer :: narn2a
    integer :: nbtd
    integer :: ncta
    integer :: npta
    integer :: nrdxsp
    integer :: nsta
    integer :: ntpr
    integer :: nxmod

    character(len=48) :: uspeca(nstmax)
    character(len=48) :: uxmod(nxmdmx)
    character(len=24) :: uphasa(nptmax)
    character(len=24) :: uptypa(nptmax)

    real(kind=8) :: axlksd(narxmx,ntprmx,nstmax)

    ! Local variable declarations.
    integer :: jlen
    integer :: n
    integer :: nb
    integer :: ncount
    integer :: nn
    integer :: np
    integer :: nr1
    integer :: nr2
    integer :: ns
    integer :: nse
    integer :: nss
    integer :: nt

    integer :: nbasis

    logical :: qheadr

    character(len=56) :: uspn56
    character(len=24) :: uptsld
    character(len=24) :: uptliq
    character(len=24) :: ux24
    character(len=24) :: uy24

    data uptsld /'Solid                   '/
    data uptliq /'Liquid                  '/

    ! Note: the following statements don't really do anything except
    ! cause the compiler not to complain that uptliq is not used.
    ! are not used.
    ux24 = uptliq
    uy24 = ux24
    uptliq = uy24

    ! Zero the jsflag and jpflag arrays.
    do ns = 1,nsta
        jsflag(ns) = 0
    end do

    do np = 1,npta
        jpflag(np) = 0
    end do

    ! Set jsflag to 2 if jflag = -1 for a species. If not, if the
    ! species is an active dependent species, set jsflag to 2 if
    ! jsflag = -1 for any other species appearing in its reaction.
    do ns = 1,nsta
        nr1 = ndrsrd(1,ns)
        nr2 = ndrsrd(2,ns)
        nt = nr2 - nr1 + 1

        if (jflag(ns) .eq. -1) then
            jsflag(ns) = 2
        else if (ns.lt.narn1a .or. ns.gt.narn2a    .or. jflag(ns) .eq. 30) then
            if (nt .ge. 2) then
                do n = nr1,nr2
                    nss = ndrsd(n)

                    if (jflag(nss) .eq. -1) then
                        jsflag(ns) = 2
                        go to 130
                    end if
                end do
            end if
        end if

130 continue
    end do

    ! Now provide an exception to the above. Scan auxiliary basis
    ! species with jflag .le. 1 and ensure that any basis species that
    ! they are linked to do not have jsflag = 2. If so, reset jsflag
    ! to 1. Do this recursively, to ensure that no species is skipped.
    ! The idea here is to keep these strict basis species in the
    ! compressed species set produced by EQLIB\cmpdat.f, so that the
    ! one-to-one association with a chemical element is preserved.
    nn = max(ncta + 1,nrdxsp + 1)
    qheadr = .true.
140 continue
    ncount = 0

    do nb = nn,nbtd
        ns = nbaspd(nb)

        if (jsflag(ns) .le. 1) then
            nr1 = ndrsrd(1,ns)
            nr2 = ndrsrd(2,ns)

            do n = nr1,nr2
                nse = ndrsd(n)

                if (jsflag(nse) .ge. 2) then
                    jsflag(nse) = 1
                    ncount = ncount + 1

                    if (qheadr) then
                        write (noutpt,1000)
                        write (nttyo,1000)
1000 format(/' Preserving the following basis species',' in the model',/' to maintain linkage to strict',' basis species associated',/' one-to-one with',' active chemical elements:',/)
                    end if

                    qheadr = .false.

                    ! Calling sequence substitutions:
                    !   uspeca(nse) for unam48
                    call fmspnx(jlen,uspeca(nse),uspn56)

                    write (noutpt,1010) uspn56(1:jlen)
                    write (nttyo,1010) uspn56(1:jlen)
1010 format('   ',a)
                end if
            end do
        end if
    end do

    if (ncount .gt. 0) then
        go to 140
    end if

    ! Set jsflag to 2 if there exist no log K data in the current
    ! temperature range, unless the species happens to be in the
    ! basis set.
    do ns = 1,nsta
        if (axlksd(1,ntpr,ns) .ge. 9999999.) then
            ! Calling sequence substitutions:
            !   nbaspd for nbasp
            !   nbtd for nbt
            nb = nbasis(nbaspd,nbtd,nbtmax,ns)

            if (nb .eq. 0) then
                jsflag(ns) = 2
            end if
        end if
    end do

    ! Execute the nxmod suppression options.
    call supprs(kxmod,jpflag,jsflag,ncmpra,noutpt,npta,nptmax,nsta,nstmax,nttyo,nxmdmx,nxmod,uphasa,uspeca,uxmod)

    ! Suppress any phase with zero species.
    do np = 1,npta
        nt = ncmpra(2,np) - ncmpra(1,np) + 1

        if (nt .le. 0) then
            jpflag(np) = 2
        end if
    end do

    ! Exercise the option to hard suppress all solid solutions.
    if (iopt(4) .le. 0) then
        do np = 1,npta
            if (uptypa(np)(1:24) .eq. uptsld(1:24)) then
                nr1 = ncmpra(1,np)
                nr2 = ncmpra(2,np)
                nt = nr2 - nr1 + 1

                if (nt .ge. 2) then
                    jpflag(np) = 2

                    do ns = nr1,nr2
                        jsflag(ns) = 2
                    end do
                end if
            end if
        end do
    end if
end subroutine flgset