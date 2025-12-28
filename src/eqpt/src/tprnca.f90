subroutine tprnca(amunca,apsi,inapr,incatr,incpr,ipbtmx,jpfcmx,natmax,ncvnca,nerr,nnapr,nncatr,nncpr,noutpt,npx3mx,npx3t,nttyo,nwarn,pcvnca,qpdna,qpdnca,qpdnc,uaqsp,utripl)
    !! Test and process the Pitzer data for nca (neutral, cation,
    !! anion) triplets read from the DATA0 file. Find and flag errors,
    !! such as duplication of data (e.g., two data blocks for the
    !! same nca triplet). Calculate the conventional primitive Pitzer
    !! parameters from the observable compound parameters read from
    !! the data file:
    !!   zeta(nca) -> mu(nca)
    !! Check the coverage of entered Pitzer data against all possible
    !! nca triplets that can be composed of the neutral species, cations,
    !! and anions present on the data file.
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
    integer :: nnapr
    integer :: nncatr
    integer :: nncpr
    integer :: npx3mx

    integer :: noutpt
    integer :: nttyo

    integer :: inapr(2,nnapr)
    integer :: incatr(3,nncatr)
    integer :: incpr(2,nncpr)

    integer :: ncvnca
    integer :: nerr
    integer :: npx3t
    integer :: nwarn

    logical :: qpdna(nnapr)
    logical :: qpdnca(nncatr)
    logical :: qpdnc(nncpr)

    character(len=24) :: uaqsp(natmax)
    character(len=24) :: utripl(3,npx3mx)

    real(kind=8) :: amunca(jpfcmx,nncatr)
    real(kind=8) :: apsi(jpfcmx,npx3mx)

    real(kind=8) :: pcvnca

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

    character(len=80) :: ustr80
    character(len=24) :: unam1
    character(len=24) :: unam2
    character(len=24) :: unam3
    character(len=8) :: ux8

    ! Limit on the list of triplets for which no data were found.
    data nlistl / 20 /

    ! Initialize the data arrays.
    do n = 1,nncatr
        qpdnca(n) = .false.

        do j = 1,jpfcmx
            amunca(j,n) = 0.
        end do
    end do

    ! Check the entered data for nca triplets.
    nodatc = 0

    do n = 1,nncatr
        i = incatr(1,n)
        j = incatr(2,n)
        k = incatr(3,n)
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
            ! Note that qpdnca(n) is left with a value of .false.
            nodatc = nodatc + 1
        else
            ! Have found an entry.
            qpdnca(n) = .true.

            ! Store the data.
            ! Note: here "cphi" is really zeta.
            do j = 1,jpfcmx
                amunca(j,n) = apsi(j,jtripl)/6.
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
1010 format(/' * Error - (EQPT/tprnca) Have found a',' duplicate data block on the DATA0 file',/7x,'for the nca triplet ',a,', ',a,', ',a,'.')
                else
                    ux8 = ' '
                    write (ux8,'(i5)') ndupl
                    call lejust(ux8)
                    j6 = ilnobl(ux8)
                    write (noutpt,1020) ux8(1:j6),unam1(1:j2),unam2(1:j3),unam3(1:j4)
                    write (nttyo,1020) ux8(1:j6),unam1(1:j2),unam2(1:j3),unam3(1:j4)
1020 format(/' * Error - (EQPT/tprnca) Have found ',a,' duplicate data blocks on the DATA0 file',/7x,'for the nca triplet ',a,', ',a,', ',a,'.')
                end if

                nerr = nerr + ndupl
            end if

            ! Check for the presence of the corresponding lamda(nc) and
            ! lambda(na) data. The presence of these data on the data
            ! file is required.
            qfound1 = .false.

            do nn = 1,nncpr
                ii = incpr(1,nn)

                if (uaqsp(ii)(1:24) .eq. unam1(1:24)) then
                    jj = incpr(2,nn)

                    if (uaqsp(jj)(1:24) .eq. unam2(1:24)) then
                        if (qpdnc(nn)) then
                            qfound1 = .true.
                            go to 130
                        end if
                    end if
                end if
            end do

130 continue

            qfound2 = .false.

            do nn = 1,nnapr
                ii = inapr(1,nn)

                if (uaqsp(ii)(1:24) .eq. unam1(1:24)) then
                    jj = inapr(2,nn)

                    if (uaqsp(jj)(1:24) .eq. unam3(1:24)) then
                        if (qpdna(nn)) then
                            qfound2 = .true.
                            go to 140
                        end if
                    end if
                end if
            end do

140 continue

            if (.not.(qfound1 .and. qfound2)) then
                ! Certain required data were not found on the DATA0 file.
                write (noutpt,1030) unam1(1:j2),unam2(1:j3),unam3(1:j4)
                write (nttyo,1030) unam1(1:j2),unam2(1:j3),unam3(1:j4)
1030 format(/' * Error - (EQPT/tprnca) Have data on the DATA0',' file for the nca',/7x,'triplet ',a,', ',a,', ',a,', but'," don't have the",/7x,'required data for the following',' neutral-ion pair(s):',/)

                if (.not.qfound1) then
                    write (noutpt,1040) unam1(1:j2),unam2(1:j3)
                    write (nttyo,1040) unam1(1:j2),unam2(1:j3)
1040 format(9x,a,', ',a)
                end if

                if (.not.qfound2) then
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
1050 format(/' * Warning - (EQPT/tprnca) Did not find a data',' block on the DATA0 file',/7x,'for any of the following',' nca triplets:',/)

        ncount = 0

        do n = 1,nncatr
            if (.not.qpdnca(n)) then
                ncount = ncount + 1
                i = incatr(1,n)
                j = incatr(2,n)
                k = incatr(3,n)
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

    ncvnca = nncatr - nodatc
    pcvnca = (100.*ncvnca)/float(nncatr)
end subroutine tprnca