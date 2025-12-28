subroutine tprna(abeta,alamna,inapr,ipbtmx,jpfcmx,natmax,ncvna,nerr,nnapr,noutpt,npx2mx,npx2t,nttyo,nwarn,pcvna,qpdna,uaqsp,upair)
    !! Test and process the na (neutral, anion) pair Pitzer data read
    !! from the DATA0 file. Find and flag errors, such as duplication of
    !! data (e.g., two data blocks for the same na pair). The
    !! conventional primitive Pitzer parameters are identical to the
    !! observable parameters read from the data file:
    !!   lambda(na) -> lambda(na)
    !! Check the coverage of entered Pitzer data against all possible
    !! na pairs that can be composed of the aqueous neutral species and
    !! anions present on the data file.
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
    integer :: npx2mx

    integer :: noutpt
    integer :: nttyo

    integer :: inapr(2,nnapr)

    integer :: ncvna
    integer :: nerr
    integer :: npx2t
    integer :: nwarn

    logical :: qpdna(nnapr)

    character(len=24) :: uaqsp(natmax)
    character(len=24) :: upair(2,npx2mx)

    real(kind=8) :: abeta(jpfcmx,0:ipbtmx,npx2mx)
    real(kind=8) :: alamna(jpfcmx,0:ipbtmx,nnapr)

    real(kind=8) :: pcvna

    ! Local variable declarations.
    integer :: i
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
    do n = 1,nnapr
        qpdna(n) = .false.

        do j = 1,jpfcmx
            do i = 0,ipbtmx
                alamna(j,i,n) = 0.
            end do
        end do
    end do

    ! Check the entered data for na pairs.
    nodatc = 0

    do n = 1,nnapr
        i = inapr(1,n)
        j = inapr(2,n)
        unam1 = uaqsp(i)
        unam2 = uaqsp(j)

        ! Search for unam1, unam2 in the upair array.
        ! That array corresponds to the species pairs blocks.
        call srch22(jpair,unam1,unam2,upair,npx2mx,npx2t)

        if (jpair .gt. 0) then
            ! Have found an entry.
            qpdna(n) = .true.

            ! Store the data.
            do j = 1,jpfcmx
                do i = 0,ipbtmx
                    alamna(j,i,n) = abeta(j,i,jpair)
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
                j2 = ilnobl(unam1)
                j3 = ilnobl(unam2)

                if (ndupl .eq. 1) then
                    write (noutpt,1010) unam1(1:j2),unam2(1:j3)
                    write (nttyo,1010) unam1(1:j2),unam2(1:j3)
1010 format(/' * Error - (EQPT/tprna) Have found a',' duplicate data block on the DATA0 file',/7x,'for the na pair ',a,', ',a,'.')
                else
                    ux8 = ' '
                    write (ux8,'(i5)') ndupl
                    call lejust(ux8)
                    j5 = ilnobl(ux8)
                    write (noutpt,1020) ux8(1:j5),unam1(1:j2),unam2(1:j3)
                    write (nttyo,1020) ux8(1:j5),unam1(1:j2),unam2(1:j3)
1020 format(/' * Error - (EQPT/tprna) Have found ',a,' duplicate data blocks on the DATA0 file',/7x,'for the na pair ',a,', ',a,'.')
                end if

                nerr = nerr + ndupl
            end if
        else
            ! No data block was found on the DATA0 file.
            ! Note that qpdna(n) is left with a value of .false.
            nodatc = nodatc + 1
        end if
    end do

    if (nodatc .gt. 0) then
        write (noutpt,1040)
        write (nttyo,1040)
1040 format(/' * Warning - (EQPT/tprna) Did not find a data',' block on the DATA0 file',/7x,'for any of the following',' na pairs:',/)

        ncount = 0

        do n = 1,nnapr
            if (.not.qpdna(n)) then
                ncount = ncount + 1
                i = inapr(1,n)
                j = inapr(2,n)
                unam1 = uaqsp(i)
                unam2 = uaqsp(j)
                j2 = ilnobl(unam1)
                j3 = ilnobl(unam2)
                ustr56 = unam1(1:j2) // ', ' // unam2(1:j3)
                j4 = ilnobl(ustr56)
                write (noutpt,1050) ustr56(1:j4)
                write (nttyo,1050) ustr56(1:j4)
1050 format(9x,a)

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
            write (noutpt,1060) ux8(1:j3)
            write (nttyo,1060) ux8(1:j3)
1060 format(/9x,'plus ',a,' others')
        end if

        write (noutpt,1070)
        write (nttyo,1070)
1070 format(1x)

        nwarn = nwarn + 1
    end if

    ncvna = nnapr - nodatc
    pcvna = (100.*ncvna)/float(nnapr)
end subroutine tprna