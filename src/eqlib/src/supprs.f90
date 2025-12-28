subroutine supprs(kxmod,jpflag,jsflag,ncmpra,noutpt,npta,nptmax,nsta,nstmax,nttyo,nxmdmx,nxmod,uphasa,uspeca,uxmod)
    !! This subroutine suppresses phases/species as directed by what is
    !! on the input file. Here uxmod is the name of the associated
    !! species and kxmod = -1. EQLIB/alters.f handles the log K alter
    !! function (kxmod = 0, 1, or 2). Suppression of a phase results in
    !! suppression of all of its component species.
    !! The string uxmod is 48 characters in length. This can contain
    !! a phase-specific species name in which the first 24 characters
    !! contain the species name proper and the second 24 characters
    !! contain the phase name. If the second 24 characters are blank,
    !! then every species or phase whose name matches what is in the
    !! first 24 characters will be suppressed. Note that it is not
    !! possible to specify a name in the second 24 character field
    !! with nothing in the first such field, as the contents of the
    !! uxmod field read from the input file is automatically left-
    !! adjusted when its contents are read. To specify a phase name
    !! only, put the string 'Phase_name_only:' in the first 24
    !! characters, where a species name proper would ordinarily
    !! appear, and put the phase name in the second 24 characters.
    !! In practice, this should not be necessary, as phase names
    !! should be unique from species names, save for the case of
    !! a pure phase, for which the species name and the phase
    !! name are identical.
    !! This subroutine is called by:
    !!   EQLIB/flgset.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: nptmax
    integer :: nstmax
    integer :: nxmdmx

    integer :: noutpt
    integer :: nttyo

    integer :: kxmod(nxmdmx)
    integer :: jsflag(nstmax)
    integer :: jpflag(nptmax)
    integer :: ncmpra(2,nptmax)
    integer :: npta
    integer :: nsta
    integer :: nxmod

    character(len=24) :: uphasa(nptmax)
    character(len=48) :: uspeca(nstmax)
    character(len=48) :: uxmod(nxmdmx)

    ! Local variable declarations.
    integer, parameter :: nlpmax = 25,nlsmax = 25

    integer :: jlen
    integer :: jlenx
    integer :: j2
    integer :: n
    integer :: nchar
    integer :: nerr
    integer :: nhitp
    integer :: nhitpl
    integer :: nhits
    integer :: nhitsl
    integer :: nlist
    integer :: nn
    integer :: np
    integer :: nr1
    integer :: nr2
    integer :: ns

    integer :: ilnobl

    logical :: qponly

    character(len=24), dimension(:), allocatable :: uusupp
    character(len=56), dimension(:), allocatable :: uusups

    character(len=56) :: uspn56
    character(len=56) :: ux56
    character(len=48) :: unam48
    character(len=48) :: ux48
    character(len=24) :: ublk24
    character(len=8) :: ufix
    character(len=8) :: ux8

    data ublk24 /'                        '/
    data ufix   /'fix     '/

    nerr = 0

    if (nxmod .le. 0) then
        go to 999
    end if

    ! Allocate arrays for lists of user-suppressed phases and species.
    ALLOCATE (uusupp(nlpmax))
    ALLOCATE (uusups(nlsmax))

    nhitpl = 0
    nhitsl = 0

    do n = 1,nxmod
        if (kxmod(n) .ge. 0) then
            go to 150
        end if

        unam48 = uxmod(n)
        nchar = 48

        qponly = .false.

        if (unam48(1:16) .eq. 'Phase_name_only:') then
            qponly = .true.
            ux48 = unam48(17:48)
            call lejust(ux48)
            j2 = ilnobl(ux48)

            if (j2 .gt. 24) then
                write (noutpt,1010) ux48(1:j2)
                write (nttyo,1010) ux48(1:j2)
1010 format(/' * Error - (EQLIB/supprs) Have encountered ',' an input file directive',/7x,'to suppress the phase',' "',a,'".',/7x,'This name exceeds the allowed',' 24 characters.')

                nerr = nerr + 1
                go to 150
            end if

            unam48(1:24) = ux48(1:24)
            unam48(25:48) = ublk24(1:24)
        end if

        j2 = ilnobl(unam48(1:24))
        call fmspnm(jlen,unam48,uspn56)

        if (unam48(25:48) .eq. ublk24(1:24)) then
            nchar = nchar - 24
        end if

        if (unam48(1:24) .eq. ublk24(1:24)) then
            nchar = nchar - 24
        end if

        if (nchar .eq. 48) then
            if (unam48(1:24) .eq. unam48(25:48)) then
                ! Have specified a pure phase species like Albite (Albite).
                ! Set up to suppress this as the pure phase.
                unam48(25:48) = ublk24(1:24)
                qponly = .true.
                nchar = 24
            end if
        end if

        if (nchar .eq. 0) then
            write (noutpt,1020)
            write (nttyo,1020)
1020 format(/' * Error - (EQLIB/supprs) Have encountered a blank',' uxmod input for an',/7x,'nxmod suppress option.')

            nerr = nerr + 1
            go to 150
        end if

        nhitp = 0

        if (nchar .eq. 24) then
            ! Look in the list of phases.
            do np = 1,npta
                if (unam48(1:24) .eq. uphasa(np)) then
                    nhitp = nhitp + 1
                    nhitpl = nhitpl + 1

                    if (nhitp .le. nlpmax) then
                        uusupp(nhitpl) = uphasa(np)
                    end if

                    if (jpflag(np) .le. 0) then
                        jpflag(np) = 1
                    end if

                    nr1 = ncmpra(1,np)
                    nr2 = ncmpra(2,np)

                    do ns = nr1,nr2
                        if (jsflag(ns) .le. 0) then
                            jsflag(ns) = 1
                        end if
                    end do

                    if (qponly) then
                        go to 150
                    end if

                    go to 110
                end if
            end do
        end if

        if (qponly) then
            ! Was only looking to suppress a phase, and did not find it.
            write (noutpt,1050) unam48(1:j2)
            write (nttyo,1050) unam48(1:j2)
1050 format(/" * Warning - (EQLIB/supprs) Can't find the phase ",'"',a,'",',/7x,'which is specified in an nxmod suppress',' option. Check to make sure',/7x,'that a phase of this',' name appears on the supporting data file.')

            go to 150
        end if

110 continue

        ! Look in the list of species.
        nhits = 0

        do ns = 1,nsta
            if (unam48(1:nchar) .eq. uspeca(ns)(1:nchar)) then
                nhits = nhits + 1

                if (jsflag(ns) .lt. 1) then
                    jsflag(ns) = 1
                end if

                if (uspeca(ns)(1:24).ne.uspeca(ns)(25:48)) then
                    call fmspnm(jlenx,uspeca(ns),ux56)
                    nhitsl = nhitsl + 1

                    if (nhitsl .le. nlsmax) then
                        uusups(nhitsl) = ux56
                    end if
                end if
            end if
        end do

        if (nhits .gt. 0) then
            go to 150
        end if

        ! The species to be suppressed was not found.
        if (nchar.eq.48) then
            write (noutpt,1120) uspn56(1:jlen)
            write (nttyo,1120) uspn56(1:jlen)
1120 format(/" * Warning - (EQLIB/supprs) Can't find the",' species "',a,'"',/7x,'which is specified in an nxmod',' suppress option. Check to make sure',/7x,'that a',' species of this name appears on the supporting',' data file.')
        else if (nhitp .le. 0) then
            ! The name could also have referred to a phase, but no
            ! such phase was found.
            write (noutpt,1130) uspn56(1:jlen)
            write (nttyo,1130) uspn56(1:jlen)
1130 format(/" * Warning - (EQLIB/supprs) Can't find an",' entity "',a,'"',/7x,'which is specified in an nxmod',' suppress option. Check to make sure',/7x,'that a',' phase or species of this name appears on the supporting',' data file.')
        end if

        ! See if an attempt was made to suppress a fictive fugacity
        ! fixing phase.
        if (unam48(1:8) .eq. ufix(1:8)) then
            write (noutpt,1150) unam48(1:j2)
            write (nttyo,1150) unam48(1:j2)
1150 format(/' * Note - (EQLIB/supprs) The phase "',a,'" specified',/7x,'in an nxmod suppress option is a',' fictive fugacity-fixing phase.',/7x,"Such a phase can't",' be suppressed.')
        end if

150 continue
    end do

    ! Write lists of user-suppressed phases and species.
    nlist = min(nhitpl,nlpmax)

    if (nlist .gt. 0) then
        if (nlist .eq. 1) then
            write (noutpt,1200)
            write (nttyo,1200)
1200 format(/' The following phase has been user-suppressed:',/)
        else
            write (noutpt,1210)
            write (nttyo,1210)
1210 format(/' The following phases have been user-suppressed:',/)
        end if

        do n = 1,nlist
            j2 = ilnobl(uusupp(n))
            write (noutpt,1220) uusupp(n)(1:j2)
            write (nttyo,1220) uusupp(n)(1:j2)
1220 format(4x,a)
        end do

        if (nhitpl .gt. nlist) then
            nn = nhitpl - nlist
            write (ux8,'(i8)') nn
            call lejust(ux8)
            j2 = ilnobl(ux8)
            write (noutpt,1230) ux8(1:j2)
            write (nttyo,1230) ux8(1:j2)
1230 format(6x,'plus ',a,' others')
        end if
    end if

    nlist = min(nhitsl,nlsmax)

    if (nlist .gt. 0) then
        if (nlist .eq. 1) then
            write (noutpt,1250)
            write (nttyo,1250)
1250 format(/' The following species has been user-suppressed:',/)
        else
            write (noutpt,1260)
            write (nttyo,1260)
1260 format(/' The following species have been user-suppressed:',/)
        end if

        do n = 1,nlist
            j2 = ilnobl(uusups(n))
            write (noutpt,1220) uusups(n)(1:j2)
            write (nttyo,1220) uusups(n)(1:j2)
        end do

        if (nhitsl .gt. nlist) then
            nn = nhitsl - nlist
            write (ux8,'(i8)') nn
            call lejust(ux8)
            j2 = ilnobl(ux8)
            write (noutpt,1230) ux8(1:j2)
            write (nttyo,1230) ux8(1:j2)
        end if
    end if

    if (nerr .gt. 0) then
        stop
    end if

    ! Deallocate arrays for lists of user-suppressed phases and species.
    DEALLOCATE (uusupp)
    DEALLOCATE (uusups)

999 continue
end subroutine supprs