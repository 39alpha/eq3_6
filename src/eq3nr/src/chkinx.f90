subroutine chkinx(cdrs,coval,ier,irdxc3,jflag,jsflag,narn1,narn2,nbasp,nbt,nbtmax,ncosp,ndrs,ndrsmx,ndrsr,nelect,nhydr,nhydx,noutpt,no2gaq,nredox,nstmax,nttyo,tempc,ucospi,uredox,uspec,zchar)
    !! This subroutine checks the code input for various kinds of
    !! errors and inconsistencies. Here ier accumulates the
    !! number of errors caught.
    !! This subroutine is called by:
    !!   EQ3NR/eq3nr.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: nbtmax
    integer :: ndrsmx
    integer :: nstmax

    integer :: noutpt
    integer :: nttyo

    integer :: jflag(nstmax)
    integer :: jsflag(nstmax)

    integer :: nbasp(nbtmax)
    integer :: ncosp(nbtmax)
    integer :: ndrs(ndrsmx)
    integer :: ndrsr(2,nstmax)

    integer :: ier
    integer :: irdxc3
    integer :: narn1
    integer :: narn2
    integer :: nbt
    integer :: nelect
    integer :: nhydr
    integer :: nhydx
    integer :: no2gaq
    integer :: nredox

    character(len=48) :: ucospi(nbtmax)
    character(len=48) :: uspec(nstmax)
    character(len=24) :: uredox

    real(kind=8) :: cdrs(ndrsmx)
    real(kind=8) :: coval(nbtmax)
    real(kind=8) :: zchar(nstmax)

    real(kind=8) :: tempc

    ! Local variable declarations.
    integer :: iredox
    integer :: jfl
    integer :: jfl2
    integer :: jlen
    integer :: j2
    integer :: j3
    integer :: j4
    integer :: n
    integer :: nalk
    integer :: nb
    integer :: nbb
    integer :: nb1
    integer :: nb2
    integer :: nr1
    integer :: nr2
    integer :: ns
    integer :: nss
    integer :: ns1
    integer :: ns2
    integer :: nt

    integer :: ilnobl
    integer :: nbasis

    character(len=56) :: uspn56
    character(len=24) :: ux24
    character(len=16) :: ux16

    real(kind=8) :: cx
    real(kind=8) :: zp

    real(kind=8) :: coefdr

    ier = 0

    ! Check usage of the log activity combination (jflag= 17) and
    ! mean log activity (jflag= 18) options.
    do nb = 1,nbt
        ns = nbasp(nb)
        jfl = jflag(ns)
        j2 = ilnobl(uspec(ns)(1:24))

        if (jfl.eq.17 .or. jfl.eq.18) then
            ns1 = ncosp(nb)

            if (jsflag(ns) .ne. 0) then
                ier = ier + 1
                write (noutpt,1000) jfl,uspec(ns)(1:j2)
                write (nttyo,1000) jfl,uspec(ns)(1:j2)
1000 format(/" * Error - (EQ3NR/chkinx) Can't apply the jflag= ",i3,' option to',/7x,a,', because this species is',' suppressed.')
            end if

            if (zchar(ns) .eq. 0.) then
                ier = ier + 1
                write (noutpt,1010) jfl,uspec(ns)(1:j2)
                write (nttyo,1010) jfl,uspec(ns)(1:j2)
1010 format(/" * Error - (EQ3NR/chkinx) Can't apply the jflag= ",i3,' option to',/7x,a,', because this species has',' zero charge.')
            end if

            ! Calling sequence substitutions:
            !   ns1 for ns
            nb1 = nbasis(nbasp,nbt,nbtmax,ns1)

            if (nb1 .le. 0) then
                ier = ier + 1
                j3 = ilnobl(ucospi(nb)(1:24))
                write (noutpt,1020) jfl,uspec(ns)(1:j2),ucospi(nb)(1:j3)
                write (nttyo,1020) jfl,uspec(ns)(1:j2),ucospi(nb)(1:j3)
1020 format(/" * Error - (EQ3NR/chkinx) Can't apply the jflag= ",i3,' option to',/7x,a,", because this species isn't",' in the basis set.')
            else
                if (jflag(ns1).le.15 .and. coval(nb1).le.0.) then
                    ier = ier + 1
                    j3 = ilnobl(uspec(ns1)(1:24))
                    write (noutpt,1030) jfl,uspec(ns)(1:j2),uspec(ns1)(1:j3)
                    write (nttyo,1030) jfl,uspec(ns)(1:j2),uspec(ns1)(1:j3)
1030 format(/" * Error - (EQ3NR/chkinx) Can't apply the",' jflag= ',i3,' option to',/7x,a,', because the',' specified counter-ion ',a,/7x,"isn't present.")
                else if (jsflag(ns1) .gt. 0) then
                    ier = ier + 1
                    j3 = ilnobl(uspec(ns1)(1:24))
                    write (noutpt,1040) jfl,uspec(ns)(1:j2),uspec(ns1)(1:j3)
                    write (nttyo,1040) jfl,uspec(ns)(1:j2),uspec(ns1)(1:j3)
1040 format(/" * Error - (EQ3NR/chkinx) Can't apply the",' jflag= ',i3,' option to',/7x,a,', because the',' specified counter-ion ',a,/7x,'is suppressed.')
                end if

                if (zchar(ns1) .eq. 0.) then
                    ier = ier + 1
                    j3 = ilnobl(uspec(ns1)(1:24))
                    write (noutpt,1050) jfl,uspec(ns)(1:j2),uspec(ns1)(1:j3)
                    write (nttyo,1050) jfl,uspec(ns)(1:j2),uspec(ns1)(1:j3)
1050 format(/" * Error - (EQ3NR/chkinx) Can't apply the",' jflag= ',i3,' option to',/7x,a,', because the',' specified counter-ion ',a,/7x,'has zero charge.')
                end if

                if (jfl .eq. 18) then
                    zp = zchar(ns)*zchar(ns1)

                    if (zp .gt. 0.) then
                        ier = ier + 1
                        j3 = ilnobl(uspec(ns1)(1:24))
                        write (noutpt,1060) uspec(ns1)(1:j3),uspec(ns)(1:j2)
                        write (nttyo,1060) jfl,uspec(ns1)(1:j3),uspec(ns)(1:j2)
1060 format(/" * Error - (EQ3NR/chkinx) Can't apply the",' jflag= ',i3,' option',/7x,'to the species ',a,' because the specified counter-ion ',a,/7x,'has',' the same charge sign.')
                    end if
                end if
            end if
        end if
    end do

    ! Check usage of the pHCl (jflag= 21) option. Note that other checks
    ! have been made previously in EQ3NR/intnsp.f, where the counter-ion
    ! has been identified (Cl- for H+, and vice versa).
    do nb = 1,nbt
        ns = nbasp(nb)
        jfl = jflag(ns)
        j2 = ilnobl(uspec(ns)(1:24))

        if (jfl .eq. 21) then
            ns1 = ncosp(nb)

            if (jsflag(ns) .ne. 0) then
                ! The species is suppressed.
                ier = ier + 1
                write (noutpt,1000) jfl,uspec(ns)(1:j2)
                write (nttyo,1000) jfl,uspec(ns)(1:j2)
            end if

            ! Calling sequence substitutions:
            !   ns1 for ns
            if (jsflag(ns1) .gt. 0) then
                ! The counter-ion is suppressed.
                ier = ier + 1
                j3 = ilnobl(uspec(ns1)(1:24))
                write (noutpt,1040) jfl,uspec(ns)(1:j2),uspec(ns1)(1:j3)
                write (nttyo,1040) jfl,uspec(ns)(1:j2),uspec(ns1)(1:j3)
            end if
        end if
    end do

    ! Check usage of the log activity (jflag= 16) and pX (jflag= 19)
    ! options.
    do nb = 1,nbt
        ns = nbasp(nb)
        jfl = jflag(ns)
        j2 = ilnobl(uspec(ns)(1:24))

        if (jfl.eq.16 .or. jfl.eq.19) then
            if (jsflag(ns) .ne. 0) then
                ! The species is suppressed.
                ier = ier + 1
                write (noutpt,1000) jfl,uspec(ns)(1:j2)
                write (nttyo,1000) jfl,uspec(ns)(1:j2)
            end if

            if (jfl .eq. 16) then
                if (coval(nb) .gt. 3.0) then
                    ier = ier + 1
                    write (ux16,'(1pg12.5)') coval(nb)
                    call lejust(ux16)
                    j4 = ilnobl(ux16)
                    write (noutpt,1120) uspec(ns)(1:j2),ux16(1:j4)
                    write (nttyo,1120) uspec(ns)(1:j2),ux16(1:j4)
1120 format(/' * Error - (EQ3NR/chkinx) The log activity of ',a,' is specified as ',/7x,a,' on the input file. This',' is unrealistically high.',/7x,'Is the sign correct?')
                end if
            else if (jfl .eq. 19) then
                if (coval(nb) .lt. -3.0) then
                    ier = ier + 1
                    write (ux16,'(1pg12.5)') coval(nb)
                    call lejust(ux16)
                    j4 = ilnobl(ux16)
                    write (noutpt,1130) uspec(ns)(1:j2),ux16(1:j4)
                    write (nttyo,1130) uspec(ns)(1:j2),ux16(1:j4)
1130 format(/' * Error - (EQ3NR/chkinx) The pX of ',a,' is specified as ',/7x,a,' on the input file. This',' is unrealistically low.',/7x,'Is the sign correct?')
                end if
            end if
        end if
    end do

    ! Check usage of the pH (jflag= 20) option.
    do nb = 1,nbt
        ns = nbasp(nb)
        jfl = jflag(ns)
        j2 = ilnobl(uspec(ns)(1:24))

        if (jfl .eq. 20) then
            if (jsflag(ns) .ne. 0) then
                ! The species is suppressed.
                ier = ier + 1
                write (noutpt,1000) jfl,uspec(ns)(1:j2)
                write (nttyo,1000) jfl,uspec(ns)(1:j2)
            end if

            if (ns .ne. nhydr) then
                ier = ier + 1
                write (noutpt,1150) uspec(ns)(1:j2)
                write (nttyo,1150) uspec(ns)(1:j2)
1150 format(/' * Error - (EQ3NR/chkinx) A pH value has been',' entered on the input file',/7x,'for ',a,'. A pH value',' can be used as a constaint only for H+.',/7x,'Use the',' pX option (jflag= 19) to enter analogous values for',' other species.')
            end if

            if (coval(nb) .gt. 17.0) then
                ier = ier + 1
                write (ux16,'(1pg12.5)') coval(nb)
                call lejust(ux16)
                j4 = ilnobl(ux16)
                write (noutpt,1160) ux16(1:j4)
                write (nttyo,1160) ux16(1:j4)
1160 format(/' * Error - (EQ3NR/chkinx) The pH value of ',a,' entered on the',/7x,'input file is unrealistically high.')
            else if (coval(nb) .lt. -3.0) then
                ier = ier + 1
                write (ux16,'(1pg12.5)') coval(nb)
                call lejust(ux16)
                j4 = ilnobl(ux16)
                write (noutpt,1170) ux16(1:j4)
                write (nttyo,1170) ux16(1:j4)
1170 format(/' * Error - (EQ3NR/chkinx) The pH value of ',a,' entered on the',/7x,'input file is unrealistically low.')
            end if
        end if
    end do

    ! Check usage of the pmH (jflag= 22) option.
    do nb = 1,nbt
        ns = nbasp(nb)
        jfl = jflag(ns)
        j2 = ilnobl(uspec(ns)(1:24))

        if (jfl .eq. 22) then
            if (jsflag(ns) .ne. 0) then
                ! The species is suppressed.
                ier = ier + 1
                write (noutpt,1000) jfl,uspec(ns)(1:j2)
                write (nttyo,1000) jfl,uspec(ns)(1:j2)
            end if

            if (ns .ne. nhydr) then
                ier = ier + 1
                write (noutpt,1250) uspec(ns)(1:j2)
                write (nttyo,1250) uspec(ns)(1:j2)
1250 format(/' * Error - (EQ3NR/chkinx) A pmH value has been',' entered on the input file',/7x,'for ',a,'. A pmH value',' can be used as a constaint only for H+.',/7x,'Use the',' pmX option (jflag= 23) to enter analogous values for',' other species.')
            end if

            if (coval(nb) .gt. 17.0) then
                ier = ier + 1
                write (ux16,'(1pg12.5)') coval(nb)
                call lejust(ux16)
                j4 = ilnobl(ux16)
                write (noutpt,1260) ux16(1:j4)
                write (nttyo,1260) ux16(1:j4)
1260 format(/' * Error - (EQ3NR/chkinx) The pmH value of ',a,' entered on the',/7x,'input file is unrealistically high.')
            else if (coval(nb) .lt. -3.0) then
                ier = ier + 1
                write (ux16,'(1pg12.5)') coval(nb)
                call lejust(ux16)
                j4 = ilnobl(ux16)
                write (noutpt,1270) ux16(1:j4)
                write (nttyo,1270) ux16(1:j4)
1270 format(/' * Error - (EQ3NR/chkinx) The pmH value of ',a,' entered on the',/7x,'input file is unrealistically low.')
            end if
        end if
    end do

    ! Check usage of the pmX (jflag= 23) option.
    do nb = 1,nbt
        ns = nbasp(nb)
        jfl = jflag(ns)
        j2 = ilnobl(uspec(ns)(1:24))

        if (jfl .eq. 23) then
            if (jsflag(ns) .ne. 0) then
                ! The species is suppressed.
                ier = ier + 1
                write (noutpt,1000) jfl,uspec(ns)(1:j2)
                write (nttyo,1000) jfl,uspec(ns)(1:j2)
            end if

            if (coval(nb) .lt. -3.0) then
                ier = ier + 1
                write (ux16,'(1pg12.5)') coval(nb)
                call lejust(ux16)
                j4 = ilnobl(ux16)
                write (noutpt,1290) uspec(ns)(1:j2),ux16(1:j4)
                write (nttyo,1290) uspec(ns)(1:j2),ux16(1:j4)
1290 format(/' * Error - (EQ3NR/chkinx) The pmX of ',a,' is specified as ',/7x,a,' on the input file. This',' is unrealistically low.',/7x,'Is the sign correct?')
            end if
        end if
    end do

    ! Check the irdxc3 .eq. -3 option for the needed input constraints.
    if (irdxc3 .eq. -3) then
        if (jflag(no2gaq) .ne. 25) then
            ier = ier + 1
            write (noutpt,1540)
            write (nttyo,1540)
1540 format(/' * Error - (EQ3NR/chkinx) The value of irdxc3 on',' the input file',/7x,'is -3, but the jflag value for O2(g)',' is not 25, as is required.')
        end if
    end if

    ! Check the irdxc3 .eq. 1 option to make sure it implies a redox
    ! reaction.
    if (irdxc3 .ne. 1) then
        go to 110
    end if

    ! Calling sequence substitutions:
    !   no2gaq for nse
    !   nredox for ns
    cx = coefdr(cdrs,ndrs,ndrsmx,ndrsr,no2gaq,nredox,nstmax)

    if (cx .eq. 0.) then
        ier = ier + 1
        j2 = ilnobl(uredox)
        write (noutpt,1570) uredox(1:j2)
        write (nttyo,1570) uredox(1:j2)
1570 format(/' * Error - (EQ3NR/chkinx) The reaction for the',' species ',a,/7x,'is specified on the input file to',' define the default redox condition.',/7x,'However, this',' is not a redox reaction.')

        ! Calling sequence substitutions:
        !   noutpt for nf
        !   nredox for ns
        call prreac(cdrs,ndrs,ndrsmx,ndrsr,noutpt,nredox,nstmax,uspec)

        ! Calling sequence substitutions:
        !   nttyo for nf
        !   nredox for ns
        call prreac(cdrs,ndrs,ndrsmx,ndrsr,nttyo,nredox,nstmax,uspec)
    end if

    ! Check for the needed input constraints on the irdxc3 .eq. 1
    ! option.
    jfl = jflag(nredox)

    if (jfl .le. 15) then
        ! Calling sequence substitutions:
        !   nredox for ns
        iredox = nbasis(nbasp,nbt,nbtmax,nredox)

        if (coval(iredox) .le. 0.) then
            ier = ier + 1
            j2 = ilnobl(uredox)
            write (noutpt,1580) uredox(1:j2)
            write (nttyo,1580) uredox(1:j2)
1580 format(/' * Error - (EQ3NR/chkinx) The species ',a," can't",/7x,'be used to define a redox couple because this species',"isn't present.")
        end if
    else if (jfl.eq.27 .or. jfl.eq.30) then
        ier = ier + 1
        write (noutpt,1590) uredox,jfl
        write (nttyo,1590) uredox,jfl
1590 format(/' * Error - (EQ3NR/chkinx) The reaction for the',' species ',a,/7x,'is specified on the input file to',' define the default redox condition.',/7x,'However, this',' species has an illegal jflag value of ',i3,'.')
    end if

    nr1 = ndrsr(1,nredox)
    nr2 = ndrsr(2,nredox)

    do n = nr1,nr2
        nss = ndrs(n)

        if (nss.ne.nredox .and. nss.ne.narn1 .and. nss.ne.nhydr  .and. nss.ne.nhydx .and. nss.ne.no2gaq .and. nss.ne.nelect) then
            go to 100
        end if
    end do

    go to 110

    ! Calling sequence substitutions:
    !   nss for ns
100 continue
    nbb = nbasis(nbasp,nbt,nbtmax,nss)

    if (jflag(nss).le.15 .and. coval(nbb).le.0.) then
        ier = ier + 1
        j2 = ilnobl(uspec(nss))
        write (noutpt,1580) uspec(nss)(1:j2)
        write (nttyo,1580) uspec(nss)(1:j2)
    end if

110 continue

    ! Check for negative values for non-negative coval quantities.
    do nb = 1,nbt
        ns = nbasp(nb)
        jfl = jflag(ns)

        if (jfl .le. 15) then
            if (coval(nb).lt.0. .and. ns.ne.no2gaq) then
                ier = ier + 1
                j2 = ilnobl(uspec(ns)(1:24))
                write (ux16,'(1pg12.5)') coval(nb)
                call lejust(ux16)
                j4 = ilnobl(ux16)
                write (noutpt,1610) ux16(1:j4),uspec(ns)(1:j2)
                write (nttyo,1610) ux16(1:j4),uspec(ns)(1:j2)
1610 format(/' * Error - (EQ3NR/chkinx) Have read an illegal',' negative concentration',/7x,'(coval) value of ',a,' from the input file for ',a,'.')
            end if
        end if
    end do

    ! Check for the alkalinity constraint on a species other than HCO3-
    ! or CO3--.
    do nb = 1,nbt
        ns = nbasp(nb)
        jfl = jflag(ns)

        if (jfl.ge.7 .and. jfl.le.11) then
            if (uspec(ns)(1:6).ne.'HCO3- ' .and.      uspec(ns)(1:6).ne.'CO3-- ') then
                ier = ier + 1
                j2 = ilnobl(uspec(ns)(1:24))
                write (noutpt,1620) uspec(ns)(1:j2)
                write (nttyo,1620) uspec(ns)(1:j2)
1620 format(/' * Error - (EQ3NR/chkinx) An illegal alkalinity',/7x,'constraint on the species ',a,' was read from the',/7x,'input file. This constraint can only be placed on',/7x,'HCO3- or CO3--.')
            end if
        end if
    end do

    ! Check for an alkalinity constraint on more than one species.
    nalk = 0

    do nb = 1,nbt
        ns = nbasp(nb)
        jfl = jflag(ns)

        if (jfl.ge.7 .and. jfl.le.11) then
            nalk = nalk + 1

            if (nalk .eq. 1) then
                ux24 = uspec(ns)(1:24)
            else
                ier = ier + 1
                j2 = ilnobl(ux24)
                j3 = ilnobl(uspec(ns)(1:24))
                write (noutpt,1630) ux24(1:j2),uspec(ns)(1:j3)
                write (nttyo,1630) ux24(1:j2),uspec(ns)(1:j3)
1630 format(/' * Error - (EQ3NR/chkinx) An alkalinity',' constraint may only',/7x,'be applied to one species in',' in a given problem. One was applied to',/7x,a,', another to ',a,'. This code only accepts as input',/7x,'the total alkalinity, not partitioned forms, such',' as the HCO3,',/7x,'CO3, and OH alkalinities.')
            end if
        end if
    end do

    if (tempc.lt.0. .or. tempc.gt.50.) then
        do nb = 1,nbt
            ns = nbasp(nb)
            jfl = jflag(ns)

            if (jfl.ge.7 .and. jfl.le.11) then
                ier = ier + 1
                j2 = ilnobl(uspec(ns)(1:24))
                write (ux16,'(1pg8.3)') tempc
                call lejust(ux16)
                j4 = ilnobl(ux16)
                write (noutpt,1650) uspec(ns)(1:j2),ux16(1:j4)
                write (nttyo,1650) uspec(ns)(1:j2),ux16(1:j4)
1650 format(/' * Error - (EQ3NR/chkinx) An alkalinity',' constraint may not',/7x,'be applied to the species ',a,' at ',a,' C. Alkalinity is',/7x,'defined in EQ3/6 only',' in the range 0-50 C.',/)
            end if
        end do
    end if

    ! Check the ncosp values for heterogeneous reaction constraints.
    do nb = 1,nbt
        ns = nbasp(nb)
        jfl = jflag(ns)

        if (jfl .eq. 25) then
            ns1 = ncosp(nb)

            ! Calling sequence substitutions:
            !   ns for nse
            !   ns1 for ns
            cx = coefdr(cdrs,ndrs,ndrsmx,ndrsr,ns,ns1,nstmax)

            if (cx.eq.0. .or. jsflag(ns1).gt.0) then
                ier = ier + 1
                j2 = ilnobl(uspec(ns)(1:24))
                write (noutpt,1670) uspec(ns)(1:j2)
                write (nttyo,1670) uspec(ns)(1:j2)
1670 format(/' * Error - (EQ3NR/chkinx) The species ',a," can't be constrained by",/7x,'the following reaction:',/)

                ! Calling sequence substitutions:
                !   noutpt for nf
                !   ns1 for ns
                call prreac(cdrs,ndrs,ndrsmx,ndrsr,noutpt,ns1,nstmax,uspec)

                ! Calling sequence substitutions:
                !   nttyo for nf
                !   ns1 for ns
                call prreac(cdrs,ndrs,ndrsmx,ndrsr,nttyo,ns1,nstmax,uspec)
            end if
        end if
    end do

    ! Check for illegal jflag value of 27 or 30 for any species in the
    ! strict basis set.
    do nb = 1,nbt
        ns = nbasp(nb)
        jfl = jflag(ns)

        if (jfl.eq.27 .or. jfl.eq.30) then
            nt = ndrsr(2,ns) - ndrsr(1,ns) + 1

            if (nt .lt. 2) then
                ier = ier + 1
                j2 = ilnobl(uspec(ns)(1:24))
                write (noutpt,1700) jfl,uspec(ns)(1:j2)
                write (nttyo,1700) jfl,uspec(ns)(1:j2)
1700 format(/' * Error - (EQ3NR/chkinx) Have an illegal jflag',' value of ',i3,/7x,'for ',a,'.')
            end if
        end if
    end do

    ! Check to see that no heterogeneous reaction constraint has
    ! been applied more than once.
    do nb1 = 1,nbt
        ns1 = nbasp(nb1)
        jfl = jflag(ns1)

        if (jfl .eq. 25) then
            do nb2 = 1,nbt
                if (nb2 .ne. nb1) then
                    ns2 = nbasp(nb2)
                    jfl2 = jflag(ns2)

                    if (jfl2 .eq. 25) then
                        if (ncosp(nb1) .eq. ncosp(nb2)) then
                            ier = ier + 1
                            write (noutpt,1720)
                            write (nttyo,1720)
1720 format(/' * Error - (EQ3NR/chkinx) The following',' heterogeneous reaction constraint',/7x,'has been',' used more than once:')
                        end if
                    end if
                end if
            end do
        end if
    end do

    ! Check to see that jflag(no2gaq) .ne. 25 unless irdxc3 .eq. -3.
    if (irdxc3 .ne. -3) then
        jfl = jflag(no2gaq)

        if (jfl .eq. 25) then
            ier = ier + 1
            write (noutpt,1750) jfl
            write (nttyo,1750) jfl
1750 format(/' * Error - (EQ3NR/chkinx) Have a jflag value of ',i3,' on the input file',/7x,'for "aqueous" O2(g), but irdxc3',' is not -3 as is required.')
        end if
    end if

    ! Check to see that only legal jflag values are applied to species
    ! other than aqueous species.
    do nb = 1,nbt
        ns = nbasp(nb)

        if (ns.lt.narn1 .or. ns.gt.narn2) then
            jfl = jflag(ns)

            if (jfl.ne.0 .and. jfl.ne.30) then
                ier = ier + 1

                ! Calling sequence substitutions:
                !   uspec(ns) for unam48
                call fmspnm(jlen,uspec(ns),uspn56)
                write (noutpt,1800) jfl,uspn56(1:jlen)
                write (nttyo,1800) jfl,uspn56(1:jlen)
1800 format(/' * Error - (EQ3NR/chkinx) An jflag value of ',i3,' has been specified',/7x,'for ',a,'. This value may not',' be applied to non-aqueous species.')
            end if
        end if
    end do

999 continue
end subroutine chkinx