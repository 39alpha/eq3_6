subroutine tprca(abeta,acphi,alamca,alpha,alphca,amua2c,amuc2a,icapr,ipbtmx,jpfcmx,natmax,na2ctr,ncapr,ncvca,nc2atr,nerr,noutpt,npx2mx,npx2t,nttyo,nwarn,pcvca,qpdca,uaqsp,upair,zaqsp)
    !! Test and process the Pitzer data for ca (cation, anion) pairs
    !! read from the DATA0 file. Find and flag errors, such as duplication
    !! of data (e.g., two data blocks for the same ca pair). Calculate
    !! the conventional primitive Pitzer parameters from the observable
    !! compound parameters read from the data file:
    !!   beta(n)(ca) -> lambda(n)(ca)   (n = 0,2)
    !!   Cphi(ca)    -> mu(cca) and mu(aac)
    !! Check the coverage of entered Pitzer data against all possible
    !! ca pairs that can be composed of the aqueous cations and anions
    !! present on the data file.
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
    integer :: na2ctr
    integer :: ncapr
    integer :: nc2atr
    integer :: npx2mx

    integer :: noutpt
    integer :: nttyo

    integer :: icapr(2,ncapr)

    integer :: ncvca
    integer :: nerr
    integer :: npx2t
    integer :: nwarn

    logical :: qpdca(ncapr)

    character(len=24) :: uaqsp(natmax)
    character(len=24) :: upair(2,npx2mx)

    real(kind=8) :: abeta(jpfcmx,0:ipbtmx,npx2mx)
    real(kind=8) :: acphi(jpfcmx,npx2mx)
    real(kind=8) :: alamca(jpfcmx,0:ipbtmx,ncapr)
    real(kind=8) :: alphca(ipbtmx,ncapr)
    real(kind=8) :: alpha(2,npx2mx)
    real(kind=8) :: amua2c(jpfcmx,na2ctr)
    real(kind=8) :: amuc2a(jpfcmx,nc2atr)
    real(kind=8) :: zaqsp(natmax)

    real(kind=8) :: pcvca

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

    real(kind=8) :: fmuc2a
    real(kind=8) :: fmua2c
    real(kind=8) :: z1
    real(kind=8) :: z2

    ! Limit on the list of pairs for which no data were found.
    data nlistl / 20 /

    ! Initialize the data arrays.
    do n = 1,ncapr
        qpdca(n) = .false.

        do i = 1,ipbtmx
            alphca(i,n) = 0.
        end do
    end do

    do n = 1,ncapr
        do i = 0,ipbtmx
            do j = 1,jpfcmx
                alamca(j,i,n) = 0.
            end do
        end do
    end do

    do n = 1,ncapr
        do j = 1,jpfcmx
            amuc2a(j,n) = 0.
            amua2c(j,n) = 0.
        end do
    end do

    ! Check the entered data for ca pairs.
    nodatc = 0

    do n = 1,ncapr
        i = icapr(1,n)
        j = icapr(2,n)
        unam1 = uaqsp(i)
        unam2 = uaqsp(j)
        z1 = zaqsp(i)
        z2 = zaqsp(j)

        ! Search for unam1, unam2 in the upair array.
        ! That array corresponds to the species pairs blocks.
        call srch22(jpair,unam1,unam2,upair,npx2mx,npx2t)

        if (jpair .gt. 0) then
            ! Have found an entry.
            qpdca(n) = .true.

            ! Store the data.
            do i = 1,ipbtmx
                alphca(i,n) = alpha(i,jpair)
            end do

            do i = 0,ipbtmx
                do j = 1,jpfcmx
                    alamca(j,i,n) = abeta(j,i,jpair)
                end do
            end do

            fmuc2a = sqrt(-z1/z2)/6.
            fmua2c = sqrt(-z2/z1)/6.

            do j = 1,jpfcmx
                amuc2a(j,n)   = acphi(j,jpair)*fmuc2a
                amua2c(j,n)   = acphi(j,jpair)*fmua2c
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
1010 format(/' * Error - (EQPT/tprca) Have found a duplicate',' data block on the DATA0 file',/7x,'for the ca pair ',a,', ',a,'.')
                else
                    ux8 = ' '
                    write (ux8,'(i5)') ndupl
                    call lejust(ux8)
                    j5 = ilnobl(ux8)
                    write (noutpt,1020) ux8(1:j5),unam1(1:j2),unam2(1:j3)
                    write (nttyo,1020) ux8(1:j5),unam1(1:j2),unam2(1:j3)
1020 format(/' * Error - (EQPT/tprca) Have found ',a,' duplicate data blocks on the DATA0 file',/7x,'for the ca pair ',a,', ',a,'.')
                end if

                nerr = nerr + ndupl
            end if
        else
            ! No data block was found on the DATA0 file.
            ! Note that qpdca(n) is left with a value of .false.
            nodatc = nodatc + 1
        end if
    end do

    if (nodatc .gt. 0) then
        write (noutpt,1040)
        write (nttyo,1040)
1040 format(/' * Warning - (EQPT/tprca) Did not find a data',' block on the DATA0 file',/7x,'for any of the following',' ca pairs:',/)

        ncount = 0

        do n = 1,ncapr
            if (.not.qpdca(n)) then
                ncount = ncount + 1
                i = icapr(1,n)
                j = icapr(2,n)
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

    ncvca = ncapr - nodatc
    pcvca = (100.*ncvca)/float(ncapr)
end subroutine tprca