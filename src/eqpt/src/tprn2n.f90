subroutine tprn2n(amun2n,apsi,innpr,in2pr,in2ntr,ipbtmx,jpfcmx,natmax,ncvn2n,nerr,nnnpr,nn2pr,nn2ntr,noutpt,npx3mx,npx3t,nttyo,nwarn,pcvn2n,qpdnn,qpdn2,qpdn2n,uaqsp,utripl)
    !! Test and process the Pitzer data for nnn' (neutral, same neutral,
    !! different neutral) triplets read from the DATA0 file. Find and
    !! flag errors, such as duplication of data (e.g., two data blocks
    !! for the same repeated neutral-distinct neutral triplet). The
    !! conventional primitive Pitzer parameters are identical to the
    !! observable parameters read from the data file:
    !!   mu(nnn') -> mu(nnn')
    !! Note that there are two mu parameters mu(nnn') and mu(n'n'n)
    !! for each corresponding lambda parameter lambda(nn'). Both are
    !! handled in one set.
    !! Check the coverage of entered Pitzer data against all possible
    !! nnn' triplets that can be composed of the neutral species present
    !! on the data file.
    !! This subroutine is called by:
    !!   EQPT/eqpt.f
    !! Principal input:
    !!   nat    = the number of aqueous species
    !!   uaqsp  = array of names of aqueous species
    !! Principal output:
    !!   nerr   = cumulative error counter
    !!   nwarn  = cumulative warning counter
    implicit none

    ! Calling sequence variable declarations.
    integer :: ipbtmx
    integer :: jpfcmx
    integer :: natmax
    integer :: nnnpr
    integer :: nn2ntr
    integer :: nn2pr
    integer :: npx3mx

    integer :: noutpt
    integer :: nttyo

    integer :: innpr(2,nnnpr)
    integer :: in2pr(nn2pr)
    integer :: in2ntr(2,nn2ntr)

    integer :: ncvn2n
    integer :: nerr
    integer :: npx3t
    integer :: nwarn

    logical :: qpdnn(nnnpr)
    logical :: qpdn2(nn2pr)
    logical :: qpdn2n(nn2ntr)

    character(len=24) :: uaqsp(natmax)
    character(len=24) :: utripl(3,npx3mx)

    real(kind=8) :: amun2n(jpfcmx,nn2ntr)
    real(kind=8) :: apsi(jpfcmx,npx3mx)

    real(kind=8) :: pcvn2n

    ! Local variable declarations.
    integer :: i
    integer :: ii
    integer :: j
    integer :: jj
    integer :: jt
    integer :: jtripl
    integer :: j2
    integer :: j3
    integer :: j4
    integer :: j5
    integer :: j6
    integer :: k
    integer :: n
    integer :: ncount
    integer :: ndupl
    integer :: nlistl
    integer :: nn
    integer :: nodatc

    integer :: ilnobl

    logical :: qfound1
    logical :: qfound2
    logical :: qfound3

    character(len=80) :: ustr80
    character(len=24) :: unam1
    character(len=24) :: unam2
    character(len=24) :: unam3
    character(len=8) :: ux8

    ! Limit on the list of triplets for which no data were found.
    data nlistl / 20 /

    ! Initialize the data arrays.
    do n = 1,nn2ntr
        qpdn2n(n) = .false.

        do j = 1,jpfcmx
            amun2n(j,n) = 0.
        end do
    end do

    ! Check the entered data for nnn' triplets.
    nodatc = 0

    do n = 1,nn2ntr
        i = in2ntr(1,n)
        j = in2ntr(1,n)
        k = in2ntr(2,n)
        unam1 = uaqsp(i)
        unam2 = uaqsp(j)
        unam3 = uaqsp(k)
        j2 = ilnobl(unam1)
        j3 = ilnobl(unam2)
        j4 = ilnobl(unam3)

        ! Search for unam1, unam2, unam3 in the utripl array.
        ! That array corresponds to the species triplets blocks.
        call srch33(jtripl,unam1,unam2,unam3,utripl,npx3mx,npx3t)

        if (jtripl .le. 0) then
            ! No data block was found on the DATA0 file.
            ! Note that qpdn2n(n) is left with a value of .false.
            nodatc = nodatc + 1
        else
            ! Have found an entry.
            qpdn2n(n) = .true.

            ! Store the data.
            do j = 1,jpfcmx
                amun2n(j,n) = apsi(j,jtripl)
            end do

            ! Check for duplicate data sets.
            ndupl = 0

            do jt = jtripl + 1,npx3t
                if (unam1(1:24) .eq. utripl(1,jt)(1:24)) then
                    if (unam2(1:24) .eq. utripl(2,jt)(1:24)) then
                        if (unam3(1:24) .eq. utripl(3,jt)(1:24)) then
                            ndupl = ndupl + 1
                        end if
                    end if
                end if
            end do

            if (ndupl .gt. 0) then
                if (ndupl .eq. 1) then
                    write (noutpt,1010) unam1(1:j2),unam2(1:j3),unam3(1:j4)
                    write (nttyo,1010) unam1(1:j2),unam2(1:j3),unam3(1:j4)
1010 format(/' * Error - (EQPT/tprn2n) Have found a',' duplicate data block on the DATA0 file',/7x,"for the nnn' triplet ",a,', ',a,', ',a,'.')
                else
                    ux8 = ' '
                    write (ux8,'(i5)') ndupl
                    call lejust(ux8)
                    j6 = ilnobl(ux8)
                    write (noutpt,1020) ux8(1:j6),unam1(1:j2),unam2(1:j3),unam3(1:j4)
                    write (nttyo,1020) ux8(1:j6),unam1(1:j2),unam2(1:j3),unam3(1:j4)
1020 format(/' * Error - (EQPT/tprn2n) Have found ',a,' duplicate data blocks on the DATA0 file',/7x,"for the nnn' triplet ",a,', ',a,', ',a,'.')
                end if

                nerr = nerr + ndupl
            end if

            ! Check for the presence of the corresponding lamda(nn),
            ! lambda(n'n'), and lambda(nn')  data. The presence of these
            ! data on the data file is required.
            qfound1 = .false.

            do nn = 1,nn2pr
                ii = in2pr(nn)

                if (uaqsp(ii)(1:24) .eq. unam1(1:24)) then
                    if (qpdn2(nn)) then
                        qfound1 = .true.
                        go to 130
                    end if
                end if
            end do

130 continue

            qfound2 = .false.

            do nn = 1,nn2pr
                ii = in2pr(nn)

                if (uaqsp(ii)(1:24) .eq. unam3(1:24)) then
                    if (qpdn2(nn)) then
                        qfound2 = .true.
                        go to 140
                    end if
                end if
            end do

140 continue

            qfound3 = .false.

            if (unam1 .lt. unam3) then
                do nn = 1,nnnpr
                    ii = innpr(1,nn)

                    if (uaqsp(ii)(1:24) .eq. unam1(1:24)) then
                        jj = innpr(2,nn)

                        if (uaqsp(jj)(1:24) .eq. unam3(1:24)) then
                            if (qpdnn(nn)) then
                                qfound3 = .true.
                                go to 150
                            end if
                        end if
                    end if
                end do

150 continue
            else
                do nn = 1,nnnpr
                    ii = innpr(1,nn)

                    if (uaqsp(ii)(1:24) .eq. unam3(1:24)) then
                        jj = innpr(2,nn)

                        if (uaqsp(jj)(1:24) .eq. unam1(1:24)) then
                            if (qpdnn(nn)) then
                                qfound3 = .true.
                                go to 160
                            end if
                        end if
                    end if
                end do

160 continue
            end if

            if (.not.(qfound1 .and. qfound2 .and. qfound3)) then
                ! Certain required data were not found on the DATA0 file.
                write (noutpt,1030) unam1(1:j2),unam2(1:j3),unam3(1:j4)
                write (nttyo,1030) unam1(1:j2),unam2(1:j3),unam3(1:j4)
1030 format(/' * Error - (EQPT/tprn2n) Have data on the DATA0'," file for the nnn'",/7x,'triplet ',a,', ',a,', ',a,', but'," don't have the",/7x,'required data for the following',' neutral-neutral pair(s):',/)

                if (.not.qfound1) then
                    write (noutpt,1040) unam1(1:j2),unam1(1:j2)
                    write (nttyo,1040) unam1(1:j2),unam1(1:j2)
1040 format(9x,a,', ',a)
                end if

                if (.not.qfound2) then
                    write (noutpt,1040) unam3(1:j4),unam3(1:j4)
                    write (nttyo,1040) unam3(1:j4),unam3(1:j4)
                end if

                if (.not.qfound3) then
                    write (noutpt,1040) unam1(1:j2),unam3(1:j4)
                    write (nttyo,1040) unam1(1:j2),unam3(1:j4)
                end if

                nerr = nerr + 1
            end if
        end if
    end do

    if (nodatc .gt. 0) then
        write (noutpt,1050)
        write (nttyo,1050)
1050 format(/' * Warning - (EQPT/tprn2n) Did not find a data',' block on the DATA0 file',/7x,'for any of the following'," nnn' triplets:",/)

        ncount = 0

        do n = 1,nn2ntr
            if (.not.qpdn2n(n)) then
                ncount = ncount + 1
                i = in2ntr(1,n)
                j = in2ntr(1,n)
                k = in2ntr(2,n)
                unam1 = uaqsp(i)
                unam2 = uaqsp(j)
                unam3 = uaqsp(k)
                j2 = ilnobl(unam1)
                j3 = ilnobl(unam2)
                j4 = ilnobl(unam3)
                ustr80 = unam1(1:j2) // ', ' // unam2(1:j3) //      ', ' // unam3(1:j4)
                j5 = ilnobl(ustr80)
                write (noutpt,1060) ustr80(1:j5)
                write (nttyo,1060) ustr80(1:j5)
1060 format(9x,a)

                if (ncount .eq. nlistl) then
                    go to 200
                end if
            end if
        end do

200 continue

        nn = nodatc - ncount

        if (nn .gt. 0) then
            write (ux8,'(i5)') nn
            call lejust(ux8)
            j3 = ilnobl(ux8)
            write (noutpt,1070) ux8(1:j3)
            write (nttyo,1070) ux8(1:j3)
1070 format(/9x,'plus ',a,' others')
        end if

        write (noutpt,1080)
        write (nttyo,1080)
1080 format(1x)

        nwarn = nwarn + 1
    end if

    ncvn2n = nn2ntr - nodatc
    pcvn2n = (100.*ncvn2n)/float(nn2ntr)
end subroutine tprn2n