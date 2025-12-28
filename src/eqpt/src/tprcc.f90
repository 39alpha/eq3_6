subroutine tprcc(alamcc,atheta,iccpr,ipbtmx,jpfcmx,natmax,nccpr,ncvcc,nerr,noutpt,npx3mx,nthdt,nttyo,nwarn,pcvcc,qpdcc,uaqsp,uthdtr)
    !! Test and process the Pitzer data for cc' (cation, different
    !! cation) pairs read from the DATA0 file. Find and flag errors,
    !! such as duplication of data (e.g., two data blocks for the same
    !! cc' pair). Calculate the conventional primitive Pitzer parameters
    !! from the observable compound parameters read from the data file:
    !!   theta(cc') -> lambda(cc')
    !! Check the coverage of entered Pitzer data against all possible
    !! cation-distinct cation pairs that can be composed of the aqueous
    !! cations present on the data file.
    !! This subroutine is complementary to tpraa.f, which does the
    !! same thing for aa' (anion, different anion) pairs.
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
    integer :: nccpr
    integer :: npx3mx

    integer :: noutpt
    integer :: nttyo

    integer :: iccpr(2,nccpr)

    integer :: ncvcc
    integer :: nerr
    integer :: nthdt
    integer :: nwarn

    logical :: qpdcc(nccpr)

    character(len=24) :: uaqsp(natmax)
    character(len=24) :: uthdtr(3,npx3mx)

    real(kind=8) :: alamcc(jpfcmx,0:ipbtmx,nccpr)
    real(kind=8) :: atheta(jpfcmx,npx3mx)

    real(kind=8) :: pcvcc

    ! Local variable declarations.
    integer :: i
    integer :: j
    integer :: jth
    integer :: jthpr
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

    logical :: qmatch
    logical :: qzerov

    character(len=56) :: ustr56
    character(len=24) :: unam1
    character(len=24) :: unam2
    character(len=8) :: ux8

    ! Limit on the list of pairs for which no data were found.
    data nlistl / 20 /

    ! Initialize the data arrays.
    do n = 1,nccpr
        qpdcc(n) = .false.

        do j = 1,jpfcmx
            do i = 0,ipbtmx
                alamcc(j,i,n) = 0.
            end do
        end do
    end do

    ! Check the entered data for cc' pairs.
    nodatc = 0

    do n = 1,nccpr
        i = iccpr(1,n)
        j = iccpr(2,n)
        unam1 = uaqsp(i)
        unam2 = uaqsp(j)

        jthpr = 0
        ndupl = 0

        ! Search for unam1, unam2 in the species triplet blocks.
        ! This is actually done by searching the uthdtr array,
        ! not the utripl array.
        do jth = 1,nthdt
            if (unam1(1:24) .eq. uthdtr(1,jth)(1:24)) then
                if (unam2(1:24) .eq. uthdtr(2,jth)(1:24)) then
                    jthpr = jth
                    go to 110
                end if
            end if
        end do

110 continue

        if (jthpr .gt. 0) then
            ! Have found an entry in the triplet blocks.
            qpdcc(n) = .true.

            ! Store the data.
            do j = 1,jpfcmx
                alamcc(j,0,n) = atheta(j,jthpr)

                do i = 1,ipbtmx
                    alamcc(j,i,n) = 0.
                end do
            end do

            ! Search for duplicates in the species triplets blocks.
            ! Ignore duplications in the triplet blocks if the
            ! values are all zeros (take this to mean "no data input")
            ! or if the values all match the first-encountered
            ! values (exact duplication).
            do jth = jthpr + 1,nthdt
                if (unam1(1:24) .eq. uthdtr(1,jth)(1:24)) then
                    if (unam2(1:24) .eq. uthdtr(2,jth)(1:24)) then
                        qzerov = .true.

                        do j = 1,jpfcmx
                            if (atheta(j,jth) .ne. 0.) then
                                qzerov = .false.
                            end if
                        end do

                        qmatch = .true.

                        do j = 1,jpfcmx
                            if ((atheta(j,jth) - atheta(j,jthpr)) .gt. 1.e-6) then
                                qmatch = .false.
                            end if
                        end do

                        if (.not.qzerov .and. .not.qmatch) then
                            ndupl = ndupl + 1
                        end if
                    end if
                end if
            end do

            if (ndupl .gt. 0) then
                j2 = ilnobl(unam1)
                j3 = ilnobl(unam2)

                if (ndupl .eq. 1) then
                    write (noutpt,1010) unam1(1:j2),unam2(1:j3)
                    write (nttyo,1010) unam1(1:j2),unam2(1:j3)
1010 format(/' * Error - (EQPT/tprcc) Have found a',' duplicate data block on the DATA0 file',/7x,"for the cc' pair ",a,', ',a,'.')
                else
                    ux8 = ' '
                    write (ux8,'(i5)') ndupl
                    call lejust(ux8)
                    j5 = ilnobl(ux8)
                    write (noutpt,1020) ux8(1:j5),unam1(1:j2),unam2(1:j3)
                    write (nttyo,1020) ux8(1:j5),unam1(1:j2),unam2(1:j3)
1020 format(/' * Error - (EQPT/tprcc) Have found ',a,' duplicate data blocks on the DATA0 file',/7x,"for the cc' pair ",a,', ',a,'.')
                end if

                nerr = nerr + ndupl
            end if
        end if

        if (jthpr.le.0) then
            ! No data block was found on the DATA0 file.
            ! Note that qpdcc(n) is left with a value of .false.
            nodatc = nodatc + 1
        end if
    end do

    if (nodatc .gt. 0) then
        write (noutpt,1040)
        write (nttyo,1040)
1040 format(/' * Warning - (EQPT/tprcc) Did not find a data',' block on the DATA0 file',/7x,'for any of the following'," cc' pairs:",/)

        ncount = 0

        do n = 1,nccpr
            if (.not.qpdcc(n)) then
                ncount = ncount + 1
                i = iccpr(1,n)
                j = iccpr(2,n)
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

    ncvcc = nccpr - nodatc
    pcvcc = (100.*ncvcc)/float(nccpr)

999 continue
end subroutine tprcc