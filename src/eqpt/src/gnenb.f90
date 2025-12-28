subroutine gnenb(ipbt_asv,ikt_asv,jpdblo,jpfc_asv,nap_asv,nat_asv,nazt_asv,nbt_asv,nct_asv,ndat0s,ngt_asv,nlt_asv,nmt_asv,noutpt,npt_asv,npx2_asv,npx3_asv,nsb,nst_asv,nttyo,nxt_asv,uakey)
    !! This subroutine makes a first pass through the DATA0 file to
    !! determine the necessary dimensioning of arrays. The arrays are
    !! then allocated back in the main program, and the data file is
    !! then processed.
    !! The returned array allocation size variables generally end
    !! in "_asv." For example, the number of chemical elements is
    !! normally nct. The associated allocation size variable is
    !! nct_asv.
    !! Note: although for example nxt_asv is described as the number
    !! of solid solution phases on the data file, it is actually an
    !! array allocation size variable which might have to be utilized
    !! in allocating an array, even if the number of items in the set
    !! is zero. Thus, the minimum value of one is generally imposed
    !! on any allocation size variable.
    !! Principal input:
    !!   ipbt_asv = the number of parameters in a Pitzer alpha set
    !!   ndat0s   = unit number of the stripped DATA0 file
    !! Principal output:
    !!   nap_asv  = maximum number of distinct sets of Pitzer alpha
    !!                parameters
    !!   nazt_asv = the number of aqueous species on the data file
    !!                for which hard core diameters are specified
    !!   nbt_asv  = the number of basis species on the data file
    !!   nct_asv  = the number of chemical elements on the data file
    !!   ikt_asv  = the maximum number of end-member component species
    !!                in any solid solution on the data file
    !!   jpdblo   = integer flag denoting the Pitzer data block
    !!                organization; -1 = classical, 0 = newlockile
    !!   jpfc_asv = the number of terms in the temperature function
    !!                used to represent Pitzer interaction parameters
    !!   nat_asv  = the number of aqueous species on the data file
    !!   ngt_asv  = the number of gas species on the data file
    !!   nlt_asv  = the number of pure liquid species on the data file
    !!   nmt_asv  = the number of gas species on the data file
    !!   npt_asv  = the number of phases of all types on the data file
    !!   npx2_asv = the number of pairs of species not of the same
    !!                charge sign for which Pitzer parameters are
    !!                defined; typically, one of the pair is a cation
    !!                and the other is an anion, but one or both
    !!                species may also be electrically neutral
    !!   npx3_asv = the number of triplets of species corresponding to
    !!                aqueous electrolyte mixtures for which Pitzer
    !!                parameters are defined; generally, no more than
    !!                two of these may have an electrical charge number
    !!                that is postive, negative, or zero
    !!   nst_asv  = the number of species of all types on the data file
    !!   nxt_asv  = the number of solid-solution phases on the data file
    implicit none

    ! Calling sequence variable declarations.
    integer :: ipbt_asv
    integer :: ikt_asv
    integer :: jpfc_asv
    integer :: nap_asv
    integer :: nat_asv
    integer :: nazt_asv
    integer :: nbt_asv
    integer :: nct_asv
    integer :: ngt_asv
    integer :: nlt_asv
    integer :: nmt_asv
    integer :: npt_asv
    integer :: npx2_asv
    integer :: npx3_asv
    integer :: nst_asv
    integer :: nxt_asv

    integer :: ndat0s
    integer :: noutpt
    integer :: nttyo

    integer :: jpdblo
    integer :: nsb

    character(len=8) :: uakey

    ! Local variable declarations.
    integer, parameter :: nap_par = 500

    integer :: i
    integer :: ier
    integer :: j
    integer :: ja
    integer :: j2
    integer :: j3
    integer :: j4
    integer :: j5
    integer :: ikt_sum
    integer :: ikt_x
    integer :: n
    integer :: nmax

    integer :: ilnobl

    logical :: qalpha

    character(len=80) :: udastr
    character(len=80) :: uline
    character(len=80) :: ux80
    character(len=16) :: uaqu
    character(len=16) :: uaux
    character(len=16) :: uele
    character(len=16) :: ubas
    character(len=16) :: ubdot
    character(len=16) :: ugas
    character(len=16) :: uliq
    character(len=16) :: umin
    character(len=16) :: uref
    character(len=16) :: usso
    character(len=16) :: uterm
    character(len=16) :: ustr16
    character(len=8) :: ux8

    real(kind=8) :: palpha(ipbt_asv,nap_par)
    real(kind=8) :: palphi(ipbt_asv)

    real(kind=8) :: aax
    real(kind=8) :: var

    data uaqu   /'aqueous species '/
    data uaux   /'auxiliary basis '/
    data ubas   /'basis species   '/
    data ubdot  /'bdot parameters '/
    data uele   /'elements        '/
    data ugas   /'gases           '/
    data uliq   /'liquids         '/
    data umin   /'solids          '/
    data uref   /'references      '/
    data usso   /'solid solutions '/
    data uterm  /'+---------------'/

    ! Initialize some variables to zero.
    ikt_asv = 0
    nazt_asv = 0
    nap_asv = 0
    npx2_asv = 0
    npx3_asv = 0

    if (uakey(1:8) .eq. 'SEDH    ') then
        ! Have a Simple Extended Debye-Huckel model for the activity
        ! coefficients of aqueous species. Determine the number of aqueous
        ! species for which hard core diameters are specified.
120 continue
        read (ndat0s,1000,end=990,err=995) uline
1000 format(a)

        if (uline(1:16) .ne. ubdot(1:16)) then
            go to 120
        end if

        ! Skip the terminator line.
        read (ndat0s,1000,end=990,err=995) uline

        ! Tally lines until next terminator.
        n = 0
130 continue
        read (ndat0s,1000,end=990,err=995) uline

        if (uline(1:16) .eq. uterm(1:16)) then
            go to 140
        end if

        n = n + 1
        go to 130

140 continue
        nazt_asv = n
    end if

    if (uakey(1:8) .eq. 'Pitzer  ') then
        ! Have a Pitzer model for the activity coefficients of aqueous
        ! species.
        if (jpdblo .eq. -1) then
            ! Have the "classical" Pitzer data block organization.
            ! Determine the number of pairs of ions corresponding to pure
            ! aqueous electrolytes for which Pitzer parameters are defined
            ! ('single-salt parameters').
150 continue
            read (ndat0s,1000,end=990,err=995) uline

            if (uline(1:16) .ne. 'single-salt para') then
                go to 150
            end if

            ! Skip the terminator line.
            read (ndat0s,1000,end=990,err=995) uline

            ! Determine the number of such pairs of ions (npx2_asv) by
            ! counting the number of terminator lines prior to encountering
            ! the 'mixture term parameters' line.
            ! Also determine the number of distinct sets of Pitzer alpha
            ! parameters (nap_asv).
            n = 0
            ja = 0
160 continue
            read (ndat0s,1000,end=990,err=995) uline

            if (uline(1:16) .eq. 'mixture term par') then
                backspace(ndat0s)
                go to 180
            end if

            i = index(uline,'alpha1')

            if (i .gt. 0) then
                ! Found the line with the alpha parameters.
                read (uline,1020,end=990,err=995) palphi(1),palphi(2)
1020 format(18x,2(16x,f5.1))

                if (ja .gt. 1) then
                    do j = 1,ja
                        if (abs(palpha(1,j) - palphi(1)).le.1.e-12            .and. abs(palpha(2,j) - palphi(2)).le.1.e-12) then
                            go to 170
                        end if
                    end do
                end if

                ! Not in the set; add it.
                ja = ja + 1

                if (ja .gt. nap_par) then
                    write (ux8,'(i5)') nap_par
                    call lejust(ux8)
                    j2 = ilnobl(ux8)
                    write (noutpt,1030) ux8(1:j2)
                    write (nttyo,1030) ux8(1:j2)
1030 format(/' * Error - (EQPT/gnenb) Have overflowed the',' palpha array while',/7x,'attempting to find the number',' of distinct sets of Pitzer alpha',/7x,'parameters.',' Increase the dimensioning parameter nap_par in this',/7x,'subroutine from its present value of ',a,'.')

                    stop
                end if

                palpha(1,ja) = palphi(1)
                palpha(2,ja) = palphi(2)
            end if

170 continue

            if (uline(1:16) .eq. uterm(1:16)) then
                n = n + 1
            end if

            go to 160

180 continue
            npx2_asv = n
            nap_asv = ja

            ! Determine the number of triplets of ions corresponding to
            ! aqueous electrolyte mixtures for which Pitzer parameters are
            ! defined ('mixture term parameters').
190 continue
            read (ndat0s,1000,end=990,err=995) uline

            if (uline(1:16) .ne. 'mixture term par') then
                go to 190
            end if

            ! Skip the terminator line.
            read (ndat0s,1000,end=990,err=995) uline

            ! Skip the E-theta flag line.
            read (ndat0s,1000,end=990,err=995) uline

            ! Determine the number such triplets of by counting the number
            ! of terminator lines prior to encountering the 'elements' line.
            n = 0
200 continue
            read (ndat0s,1000,end=990,err=995) uline

            if (uline(1:16) .eq. uele(1:16)) then
                backspace(ndat0s)
                go to 210
            end if

            if (uline(1:16) .eq. uterm(1:16)) then
                n = n + 1
            end if

            go to 200

210 continue
            npx3_asv = n
        else if (jpdblo .eq. 0) then
            ! Have the "new" Pitzer data block organization.
            ! First determine the number of distinct sets of Pitzer alpha
            ! parameters (nap_asv). It is possible that no values are
            ! input here, in which case only the default set holds.
            ! Load the default set into the palpha array.
            nmax = ipbt_asv*nap_par
            call initaz(palpha,nmax)
            palpha(1,1) = 2.0
            palpha(2,1) = 12.0
            palpha(1,2) = 1.4
            palpha(2,2) = 12.0
            ja = 2
            nap_asv = ja

            ! Skip down to the "ca combinations" superblock. This is
220 continue
            read (ndat0s,1000,end=990,err=995) uline

            if (uline(1:15) .ne. 'ca combinations') then
                go to 220
            end if

            ! Skip the terminator line.
            read (ndat0s,1000,end=990,err=995) uline

230 continue
            read (ndat0s,1000,end=990,err=995) uline
            j3 = ilnobl(uline)

            if (uline(1:26) .eq. "cc'a and aa'c combinations") then
                backspace(ndat0s)
                go to 250
            end if

            qalpha = .false.
            ux80 = uline
            call lejust(ux80)

            ! Find the next set of alpha lines. Look for "alpha(1)".
            ! That marks the start of a set.
            ustr16 = 'alpha(1) ='
            j4 = 10

            if (ux80(1:j4) .eq. ustr16(1:j4)) then
                ! Found a line with the alpha(1) parameter.
                j5 = index(ux80,'=')
                udastr = ux80(j5 + 1:80)
                call g1dat(ier,noutpt,nttyo,udastr,var)

                if (ier .gt. 0) then
                    write (noutpt,1110) uline(1:j3)
                    write (nttyo,1110) uline(1:j3)
1110 format(/' * Error - (EQPT/gnenb) Have found a line',' starting with:',/7x,'"',a,'"',/7x,'containing an',' expected numerical field that could not be read.',/7x,'This occurred while scanning the data file to',' determine the number',/7x,'of distinct sets of',' Pitzer alpha parameters.')

                    stop
                end if

                palphi(1) = var
                qalpha = .true.

                do i = 2,ipbt_asv
                    read (ndat0s,1000,end=990,err=995) uline
                    j3 = ilnobl(uline)
                    ux80 = uline
                    call lejust(ux80)
                    ustr16 = 'alpha( ) ='
                    j4 = 10
                    write (ustr16(7:7),'(i1)') i

                    if (ux80(1:j4) .eq. ustr16(1:j4)) then
                        ! Found a line with the expected alpha(i) parameter.
                        j5 = index(ux80,'=')
                        udastr = ux80(j5 + 1:80)
                        call g1dat(ier,noutpt,nttyo,udastr,var)

                        if (ier .le. 0) then
                            palphi(i) = var
                        else
                            write (noutpt,1110) uline(1:j3)
                            write (nttyo,1110) uline(1:j3)
                            stop
                        end if
                    else
                        ! Did not find a line with the expected alpha(i)
                        ! parameter.
                        ustr16 = 'beta(0) ='
                        j4 = 9

                        if (ux80(1:j4) .ne. ustr16(1:j4)) then
                            ! Did not find a beta(0) line either. This
                            ! condition is invalid.
                            write (ux8,'(i1)') i
                            write (noutpt,1120) uline(1:j3),ux8(1:1)
                            write (nttyo,1120) uline(1:j3),ux8(1:1)
1120 format(/' * Error - (EQPT/gnenb) Have found a line',' starting with:',/7x,'"',a,'"',/7x,'that contains',' neither the expected alpha(',a,') data nor the',/7x,'possible beta(0) data that could mark the end',' of an abbreviated',/7x,'alpha set. This occurred',' while scanning the data file',/7x,'to determine',' the number of distinct sets of Pitzer',/7x,'alpha parameters.')

                            stop
                        end if
                    end if
                end do
            end if

            if (qalpha) then
                do j = 1,ja
                    do i = 1,ipbt_asv
                        aax = abs(palpha(i,j) - palphi(i))

                        if (aax .gt. 1.e-12) then
                            go to 240
                        end if
                    end do

                    ! The current alpha combination is already in
                    ! the known set.
                    go to 230
240 continue
                end do

                ! The current alpha set is not in the set; add it.
                ja = ja + 1

                if (ja .gt. nap_par) then
                    ! Have overflowed the palpha array, which was statically
                    ! dimensioned.
                    write (ux8,'(i5)') nap_par
                    call lejust(ux8)
                    j2 = ilnobl(ux8)
                    write (noutpt,1030) ux8(1:j2)
                    write (nttyo,1030) ux8(1:j2)
                    stop
                end if

                do i = 1,ipbt_asv
                    palpha(i,ja) = palphi(i)
                end do
            end if

            go to 230

250 continue
            nap_asv = ja

            rewind(ndat0s)

            ! Skip back down to the "ca combinations" superblock.
260 continue
            read (ndat0s,1000,end=990,err=995) uline

            if (uline(1:15) .ne. 'ca combinations') then
                go to 260
            end if

            ! Skip the terminator line.
            read (ndat0s,1000,end=990,err=995) uline

            ! Determine the number of pair combinations (npx2_asv) by
            ! counting the number of terminator lines prior to
            ! encountering the "cc'a and aa'c combinations" line.
            ! Correct this estimate for other pair superblock
            ! headers.
            n = 0
270 continue
            read (ndat0s,1000,end=990,err=995) uline

            if (uline(1:26) .eq. "cc'a and aa'c combinations") then
                backspace(ndat0s)
                go to 280
            end if

            if (uline(1:24) .eq. "cc' and aa' combinations") then
                n = n - 1
            end if

            if (uline(1:22) .eq. "nc and na combinations") then
                n = n - 1
            end if

            if (uline(1:15) .eq. "nn combinations") then
                n = n - 1
            end if

            if (uline(1:16) .eq. "nn' combinations") then
                n = n - 1
            end if

            if (uline(1:16) .eq. uterm(1:16)) then
                n = n + 1
            end if

            go to 270

280 continue
            npx2_asv = n

            ! Skip down to the "cc'a and aa'c combinations" superblock.
290 continue
            read (ndat0s,1000,end=990,err=995) uline

            if (uline(1:26) .ne. "cc'a and aa'c combinations") then
                go to 290
            end if

            ! Skip the terminator line.
            read (ndat0s,1000,end=990,err=995) uline

            ! Determine the number of triplet combinations (npx3_asv) by
            ! counting the number of terminator lines prior to
            ! encountering the "cc'a and aa'c combinations" line.
            ! Correct this estimate for other triplet superblock
            ! headers.
            n = 0
300 continue
            read (ndat0s,1000,end=990,err=995) uline

            if (uline(1:16) .eq. uele(1:16)) then
                backspace(ndat0s)
                go to 310
            end if

            if (uline(1:16) .eq. "nca combinations") then
                n = n - 1
            end if

            if (uline(1:16) .eq. "nnn' combination") then
                n = n - 1
            end if

            if (uline(1:16) .eq. uterm(1:16)) then
                n = n + 1
            end if

            go to 300

310 continue
            npx3_asv = n
        end if
    end if

    ! Find the 'elements' line.
350 continue
    read (ndat0s,1000,end=990,err=995) uline

    if (uline(1:16) .ne. uele(1:16)) then
        go to 350
    end if

    ! Skip the terminator line.
    read (ndat0s,1000,end=990,err=995) uline

    ! Tally lines until next terminator.
    n = 0
360 continue
    read (ndat0s,1000,end=990,err=995) uline

    if (uline(1:16) .eq. uterm(1:16)) then
        go to 370
    end if

    n = n + 1
    go to 360

370 continue
    nct_asv = n

    ! Skip to the 'basis species' line.
400 continue
    read (ndat0s,1000,end=990,err=995) uline

    if (uline(1:16) .ne. ubas(1:16)) then
        go to 400
    end if

    ! Skip the terminator line.
    read (ndat0s,1000,end=990,err=995) uline

    ! Determine the number of basis species by counting the number
    ! of terminator lines prior to encountering the 'aqueous species'
    ! line. In this count, note that an "extra" terminator line
    ! is encountered after the "auxiliary basis species" line.
    ! The number of strict basis species is also determined, by
    ! noting the number of basis species prior to encountering the
    ! "auxiliary basis species" line.
    n = 0
410 continue
    read (ndat0s,1000,end=990,err=995) uline

    if (uline(1:16) .eq. uaux(1:16)) then
        nsb = n
    end if

    if (uline(1:16) .eq. uaqu(1:16)) then
        backspace(ndat0s)
        go to 420
    end if

    if (uline(1:16) .eq. uterm(1:16)) then
        n = n + 1
    end if

    go to 410

420 continue
    nbt_asv = n - 1

    ! Skip to the 'aqueous species' line.
450 continue
    read (ndat0s,1000,end=990,err=995) uline

    if (uline(1:16) .ne. uaqu(1:16)) then
        go to 450
    end if

    ! Skip the terminator line.
    read (ndat0s,1000,end=990,err=995) uline

    ! Determine the number of non-basis aqueous species by counting the
    ! number of terminator lines prior to encountering the 'solids'
    ! line. The total number of aqueous species is the sum of the
    ! basis and non-basis aqueous species.
    n = 0
460 continue
    read (ndat0s,1000,end=990,err=995) uline

    if (uline(1:16) .eq. umin(1:16)) then
        backspace(ndat0s)
        go to 470
    end if

    if (uline(1:16) .eq. uterm(1:16)) then
        n = n + 1
    end if

    go to 460

470 continue
    nat_asv = nbt_asv + n

    ! Skip to the 'solids' line.
500 continue
    read (ndat0s,1000,end=990,err=995) uline

    if (uline(1:16) .ne. umin(1:16)) then
        go to 500
    end if

    ! Skip the terminator line.
    read (ndat0s,1000,end=990,err=995) uline

    ! Determine the number of pure mineral species by counting the
    ! number of terminator lines prior to encountering the 'liquids'
    ! line.
    n = 0
510 continue
    read (ndat0s,1000,end=990,err=995) uline

    if (uline(1:16) .eq. uliq(1:16)) then
        backspace(ndat0s)
        go to 520
    end if

    if (uline(1:16) .eq. uterm(1:16)) then
        n = n + 1
    end if

    go to 510

520 continue
    nmt_asv = n

    ! Skip to the 'liquids' line.
550 continue
    read (ndat0s,1000,end=990,err=995) uline

    if (uline(1:16) .ne. uliq(1:16)) then
        go to 550
    end if

    ! Skip the terminator line.
    read (ndat0s,1000,end=990,err=995) uline

    ! Determine the number of pure liquid species by counting the
    ! number of terminator lines prior to encountering the 'gases'
    ! line.
    n = 0
560 continue
    read (ndat0s,1000,end=990,err=995) uline

    if (uline(1:16) .eq. ugas(1:16)) then
        backspace(ndat0s)
        go to 570
    end if

    if (uline(1:16) .eq. uterm(1:16)) then
        n = n + 1
    end if

    go to 560

    ! Note: water as a pure liquid phase is not explicitly on the
    ! DATA0 file. It is added by EQ3NR or EQ6 by cloning the liquid
    ! water species belonging to the aqueous solution phase. The pure
    ! liquid water phase must be accounted for in nlt_asv, so the
    ! actual count of pure liquid phases explicitly on the data file
    ! is incremented by one.
570 continue
    nlt_asv = n + 1

    ! Skip to the 'gases' line.
600 continue
    read (ndat0s,1000,end=990,err=995) uline

    if (uline(1:16) .ne. ugas(1:16)) then
        go to 600
    end if

    ! Skip the terminator line.
    read (ndat0s,1000,end=990,err=995) uline

    ! Determine the number of gas species by counting the number of
    ! terminator lines prior to encountering the 'solid solutions'
    ! line.
    n = 0
610 continue
    read (ndat0s,1000,end=990,err=995) uline

    if (uline(1:16) .eq. usso(1:16)) then
        backspace(ndat0s)
        go to 620
    end if

    if (uline(1:16) .eq. uterm(1:16)) then
        n = n + 1
    end if

    go to 610

620 continue
    ngt_asv = n

    ! Skip to the 'solid solutions' line.
650 continue
    read (ndat0s,1000,end=990,err=995) uline

    if (uline(1:16) .ne. usso(1:16)) then
        go to 650
    end if

    ! Skip the terminator line.
    read (ndat0s,1000,end=990,err=995) uline

    ! Determine the number of solid solution phases by counting the
    ! number of terminator lines prior to encountering the 'references'
    ! line.
    !   ikt_x   = the number of components in the current solid
    !               solution
    !   ikt_sum = the total number of components in all solid solutions
    n = 0
    ikt_x = 0
    ikt_sum = 0

660 continue
    read (ndat0s,1000,end=990,err=995) uline

    if (uline(1:16) .eq. uref(1:16)) then
        backspace(ndat0s)
        go to 670
    end if

    if (ikt_x .le. 0) then
        ! Have not yet found a 'component' line.
        i = index(uline,'component')

        if (i .eq. 5) then
            ! Have found the line with the number of components.
            ! Note that there is a two-component minimum here.
            read (uline,'(1x,i2)',err=995) ikt_x
            ikt_x = max(2,ikt_x)
            ikt_asv = max(ikt_asv,ikt_x)
            ikt_sum = ikt_sum + ikt_x
        end if
    end if

    if (uline(1:16) .eq. uterm(1:16)) then
        ! Have found the end of the current solid solution block.
        n = n + 1

        if (ikt_x .le. 0) then
            ! No line defining the number of components was found.
            ikt_x = 2
            ikt_asv = max(ikt_asv,ikt_x)
            ikt_sum = ikt_sum + ikt_x
        end if

        ikt_x = 0
    end if

    go to 660

670 continue
    nxt_asv = n

    ! Get the total number of phases. Be sure to account for the
    ! aqueous solution and the gas phase. Water as a pure liquid
    ! phase has been implicitly added to nlt_asv.
    npt_asv = 2 + nmt_asv + nlt_asv + nxt_asv

    ! Get the total number of species.
    nst_asv = nat_asv + nmt_asv + nlt_asv + ngt_asv + ikt_sum

    ! Each dimensioning parameter should have a minimum value of 1,
    ! even if the corresponding number of items is 0.
    ikt_asv = max(1,ikt_asv)
    nap_asv = max(1,nap_asv)
    nat_asv = max(1,nat_asv)
    nazt_asv = max(1,nazt_asv)
    nbt_asv = max(1,nbt_asv)
    nct_asv = max(1,nct_asv)
    ngt_asv = max(1,ngt_asv)
    nlt_asv = max(1,nlt_asv)
    nmt_asv = max(1,nmt_asv)
    npt_asv = max(1,npt_asv)
    npx2_asv = max(1,npx2_asv)
    npx3_asv = max(1,npx3_asv)
    nst_asv = max(1,nst_asv)
    nxt_asv = max(1,nxt_asv)

    ! Rewind the DATA0 file and exit.
    rewind(ndat0s)
    go to 999

    ! Write a message for any read error.
990 continue
    write (noutpt,2000)
    write (nttyo,2000)
2000 format(/' * Error - (EQPT/gnenb) Unexpectedly encountered',' end-of-file',/7x,'while scanning the DATA0 file to determine',' the necessary array',/7x,'dimensions.')

    stop

995 continue
    write (noutpt,2010)
    write (nttyo,2010)
2010 format(/' * Error - (EQPT/gnenb) Encountered a read format',' error while',/7x,'scanning the DATA0 file to determine',' the necessary array dimensions.')

    stop

999 continue
end subroutine gnenb