subroutine tprn2(abeta,acphi,alamn2,amun3,in2pr,ipbtmx,jpfcmx,natmax,ncvn2,nerr,nn2pr,nn3tr,noutpt,npx2mx,npx2t,nttyo,nwarn,pcvn2,qpdn2,uaqsp,upair)
    !! Test and process Pitzer data read from the data file that
    !! pertain to a single neutral species (nn and nnn combinations,
    !! such as CO2(aq)-CO2(aq) and CO2(aq)-CO2(aq)-CO2(aq). Find and
    !! flag errors, such as duplication of data (e.g., two data
    !! blocks for the same pair composed of a repeated neutral).
    !! The conventional primitive Pitzer parameters are identical
    !! to the observable parameters read from the data file:
    !!   lambda(nn) -> lambda(nn)
    !!   mu(nnn)    -> mu(nnn)
    !! Check the coverage of entered Pitzer data against all possible
    !! pairs that can be composed by repeating a neutral species present
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
    integer :: nn2pr
    integer :: nn3tr
    integer :: npx2mx

    integer :: noutpt
    integer :: nttyo

    integer :: in2pr(nn2pr)

    integer :: ncvn2
    integer :: nerr
    integer :: npx2t
    integer :: nwarn

    logical :: qpdn2(nn2pr)

    character(len=24) :: uaqsp(natmax)
    character(len=24) :: upair(2,npx2mx)

    real(kind=8) :: abeta(jpfcmx,0:ipbtmx,npx2mx)
    real(kind=8) :: acphi(jpfcmx,npx2mx)
    real(kind=8) :: alamn2(jpfcmx,0:ipbtmx,nn2pr)
    real(kind=8) :: amun3(jpfcmx,nn3tr)

    real(kind=8) :: pcvn2

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
    character(len=8) :: ux8

    ! Limit on the list of pairs for which no data were found.
    data nlistl / 20 /

    ! Initialize the data arrays.
    do n = 1,nn2pr
        qpdn2(n) = .false.

        do j = 1,jpfcmx
            do i = 0,ipbtmx
                alamn2(j,i,n) = 0.
            end do
        end do
    end do

    ! Note: nn3tr = nn2pr.
    do n = 1,nn3tr
        do j = 1,jpfcmx
            amun3(j,n) = 0.
        end do
    end do

    ! Check the entered data for nn pairs.
    nodatc = 0

    do n = 1,nn2pr
        i = in2pr(n)
        unam1 = uaqsp(i)

        ! Search for unam1 in the upair array.
        call srch22(jpair,unam1,unam1,upair,npx2mx,npx2t)

        if (jpair .gt. 0) then
            ! Have found an entry.
            qpdn2(n) = .true.

            ! Store the data.
            do j = 1,jpfcmx
                do i = 0,ipbtmx
                    alamn2(j,i,n) = abeta(j,i,jpair)
                end do
            end do

            do j = 1,jpfcmx
                amun3(j,n)   = acphi(j,jpair)
            end do

            ! Check for duplicate data sets.
            ndupl = 0

            do jp = jpair + 1,npx2t
                if (unam1(1:24) .eq. upair(1,jp)(1:24)) then
                    if (unam1(1:24) .eq. upair(2,jp)(1:24)) then
                        ndupl = ndupl + 1
                    end if
                end if
            end do

            if (ndupl .gt. 0) then
                j2 = ilnobl(unam1)

                if (ndupl .eq. 1) then
                    write (noutpt,1010) unam1(1:j2),unam1(1:j2)
                    write (nttyo,1010) unam1(1:j2),unam1(1:j2)
1010 format(/' * Error - (EQPT/tprn2) Have found a',' duplicate data block on the DATA0 file',/7x,'for the nn pair ',a,', ',a,'.')
                else
                    ux8 = ' '
                    write (ux8,'(i5)') ndupl
                    call lejust(ux8)
                    j5 = ilnobl(ux8)
                    write (noutpt,1020) ux8(1:j5),unam1(1:j2),unam1(1:j2)
                    write (nttyo,1020) ux8(1:j5),unam1(1:j2),unam1(1:j2)
1020 format(/' * Error - (EQPT/tprn2) Have found ',a,' duplicate data blocks on the DATA0 file',/7x,'for the nn pair ',a,', ',a,'.')
                end if

                nerr = nerr + ndupl
            end if
        else
            ! No data block was found on the DATA0 file.
            ! Note that qpdn2(n) is left with a value of .false.
            nodatc = nodatc + 1
        end if
    end do

    if (nodatc .gt. 0) then
        write (noutpt,1050)
        write (nttyo,1050)
1050 format(/' * Warning - (EQPT/tprn2) Did not find a data',' block on the DATA0 file',/7x,'for any of the following',' nn pairs:',/)

        ncount = 0

        do n = 1,nn2pr
            if (.not.qpdn2(n)) then
                ncount = ncount + 1
                i = in2pr(n)
                unam1 = uaqsp(i)
                j2 = ilnobl(unam1)
                ustr56 = unam1(1:j2) // ', ' // unam1(1:j2)
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

    ncvn2 = nn2pr - nodatc
    pcvn2 = (100.*ncvn2)/float(nn2pr)
end subroutine tprn2