subroutine tprnn(abeta,acphi,alamnn,innpr,in2pr,ipbtmx,jpfcmx,natmax,ncvnn,nerr,nnnpr,nn2pr,noutpt,npx2mx,npx2t,nttyo,nwarn,pcvnn,qpdnn,qpdn2,uaqsp,upair)
    !! Test and process the Pitzer data for nn' (neutral, different
    !! neutral) pairs read from the DATA0 file. Find and flag errors,
    !! such as duplication of data (e.g., two data blocks for the same
    !! nn' pair). The conventional primitive Pitzer parameters are
    !! identical to the observable parameters read from the data file:
    !!   lambda(nn') -> lambda(nn')
    !! Check the coverage of entered Pitzer data against all possible
    !! nn' pairs that can be composed of the neutral species present
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
    integer :: nn2pr
    integer :: npx2mx

    integer :: noutpt
    integer :: nttyo

    integer :: innpr(2,nnnpr)
    integer :: in2pr(nn2pr)

    integer :: ncvnn
    integer :: nerr
    integer :: npx2t
    integer :: nwarn

    logical :: qfound1
    logical :: qfound2

    logical :: qpdnn(nnnpr)
    logical :: qpdn2(nn2pr)

    character(len=24) :: uaqsp(natmax)
    character(len=24) :: upair(2,npx2mx)

    real(kind=8) :: abeta(jpfcmx,0:ipbtmx,npx2mx)
    real(kind=8) :: acphi(jpfcmx,npx2mx)
    real(kind=8) :: alamnn(jpfcmx,0:ipbtmx,nnnpr)

    real(kind=8) :: pcvnn

    ! Local variable declarations.
    integer :: i
    integer :: ii
    integer :: j
    integer :: jp
    integer :: jpair
    integer :: j2
    integer :: j3
    integer :: j4
    integer :: j5
    integer :: n
    integer :: ncount
    integer :: ndupl
    integer :: nlistl
    integer :: nn
    integer :: nodatc

    integer :: ilnobl

    character(len=56) :: ustr56
    character(len=24) :: unam1
    character(len=24) :: unam2
    character(len=8) :: ux8

    ! Limit on the list of pairs for which no data were found.
    data nlistl / 20 /

    ! Initialize the data arrays.
    do n = 1,nnnpr
        qpdnn(n) = .false.

        do j = 1,jpfcmx
            do i = 0,ipbtmx
                alamnn(j,i,n) = 0.
            end do
        end do
    end do

    ! Check the entered data for nn' pairs.
    nodatc = 0

    do n = 1,nnnpr
        i = innpr(1,n)
        j = innpr(2,n)
        unam1 = uaqsp(i)
        unam2 = uaqsp(j)
        j2 = ilnobl(unam1)
        j3 = ilnobl(unam2)

        ! Search for unam1, unam2 in the upair array.
        ! That array corresponds to the species pairs blocks.
        call srch22(jpair,unam1,unam2,upair,npx2mx,npx2t)

        if (jpair .gt. 0) then
            ! Have found an entry.
            qpdnn(n) = .true.

            ! Store the data.
            do j = 1,jpfcmx
                do i = 0,ipbtmx
                    alamnn(j,i,n) = abeta(j,i,jpair)
                end do
            end do

            ! Check for duplicate data sets.
            ndupl = 0

            do jp = jpair + 1,npx2t
                if (unam1(1:24) .eq. upair(1,jp)(1:24)) then
                    if (unam2(1:24) .eq. upair(2,jp)(1:24)) then
                        ndupl = ndupl + 1
                    end if
                end if
            end do

            if (ndupl .gt. 0) then
                if (ndupl .eq. 1) then
                    write (noutpt,1010) unam1(1:j2),unam2(1:j3)
                    write (nttyo,1010) unam1(1:j2),unam2(1:j3)
1010 format(/' * Error - (EQPT/tprnn) Have found a',' duplicate data block on the DATA0 file',/7x,"for the nn' pair ",a,', ',a,'.')
                else
                    ux8 = ' '
                    write (ux8,'(i5)') ndupl
                    call lejust(ux8)
                    j5 = ilnobl(ux8)
                    write (noutpt,1020) ux8(1:j5),unam1(1:j2),unam2(1:j3)
                    write (nttyo,1020) ux8(1:j5),unam1(1:j2),unam2(1:j3)
1020 format(/' * Error - (EQPT/tprnn) Have found ',a,' duplicate data blocks on the DATA0 file',/7x,"for the nn' pair ",a,', ',a,'.')
                end if

                nerr = nerr + ndupl
            end if
        else
            ! No data block was found on the DATA0 file.
            ! Note that qpdnn(n) is left with a value of .false.
            nodatc = nodatc + 1
        end if

        ! Check for the presence of the corresponding lamda(nn),
        ! and lambda(n'n') data. The presence of these data on
        ! the  data file is required.
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

            if (uaqsp(ii)(1:24) .eq. unam2(1:24)) then
                if (qpdn2(nn)) then
                    qfound2 = .true.
                    go to 140
                end if
            end if
        end do

140 continue

        if (qpdnn(n) .and. .not.(qfound1 .and. qfound2)) then
            ! Certain required data were not found on the DATA0 file.
            write (noutpt,1030) unam1(1:j2),unam2(1:j3)
            write (nttyo,1030) unam1(1:j2),unam2(1:j3)
1030 format(/' * Error - (EQPT/tprnn) Have data on the DATA0'," file for the nn'",/7x,'pair ',a,', ',a,', but'," don't have the",/7x,'required data for the following',' neutral-neutral pair(s):',/)

            if (.not.qfound1) then
                write (noutpt,1040) unam1(1:j2),unam1(1:j2)
                write (nttyo,1040) unam1(1:j2),unam1(1:j2)
1040 format(9x,a,', ',a)
            end if

            if (.not.qfound2) then
                write (noutpt,1040) unam2(1:j3),unam2(1:j3)
                write (nttyo,1040) unam2(1:j3),unam2(1:j3)
            end if

            nerr = nerr + 1
        end if
    end do

    if (nodatc .gt. 0) then
        write (noutpt,1050)
        write (nttyo,1050)
1050 format(/' * Warning - (EQPT/tprnn) Did not find a data',' block on the DATA0 file',/7x,'for any of the following'," nn' pairs:",/)

        ncount = 0

        do n = 1,nnnpr
            if (.not.qpdnn(n)) then
                ncount = ncount + 1
                i = innpr(1,n)
                j = innpr(2,n)
                unam1 = uaqsp(i)
                unam2 = uaqsp(j)
                j2 = ilnobl(unam1)
                j3 = ilnobl(unam2)
                ustr56 = unam1(1:j2) // ', ' // unam2(1:j3)
                j4 = ilnobl(ustr56)
                write (noutpt,1060) ustr56(1:j4)
                write (nttyo,1060) ustr56(1:j4)
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

    ncvnn = nnnpr - nodatc
    pcvnn = (100.*ncvnn)/float(nnnpr)
end subroutine tprnn