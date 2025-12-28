subroutine tpraac(amuaac,amua2c,apsi,iaactr,icapr,ipbtmx,jpfcmx,naactr,natmax,na2ctr,ncapr,ncvaac,nerr,noutpt,npx3mx,npx3t,nttyo,nwarn,pcvaac,qpdaac,qpdca,uaqsp,utripl,zaqsp)
    !! Test and process the Pitzer data for aa'c (anion, different
    !! anion, cation) triplets read from the DATA0 file. Find and flag
    !! errors, such as duplication of data (e.g., two data blocks for
    !! the same aa'c triplet). Calculate the conventional primitive
    !! Pitzer parameters from the observable compound parameters read
    !! from the data file:
    !!   psi(aa'c) -> mu(aa'c)
    !! Note: mu(aa'c) depends on cphi(ca) and cphi(c'a) as well as
    !! psi(aa'c).
    !! Check the coverage of entered Pitzer data against all possible
    !! aa'c triplets that can be composed of the cations and anions
    !! present on the data file.
    !! This subroutine is complementary to tprcca.f, which does the
    !! same thing for cc'a (cation, different cation, anion) triplets.
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
    integer :: naactr
    integer :: natmax
    integer :: na2ctr
    integer :: ncapr
    integer :: npx3mx

    integer :: noutpt
    integer :: nttyo

    integer :: iaactr(3,naactr)
    integer :: icapr(2,ncapr)

    integer :: ncvaac
    integer :: nerr
    integer :: npx3t
    integer :: nwarn

    logical :: qpdaac(naactr)
    logical :: qpdca(ncapr)

    character(len=24) :: uaqsp(natmax)
    character(len=24) :: utripl(3,npx3mx)

    real(kind=8) :: apsi(jpfcmx,npx3mx)
    real(kind=8) :: amuaac(jpfcmx,naactr)
    real(kind=8) :: amua2c(jpfcmx,na2ctr)
    real(kind=8) :: zaqsp(natmax)

    real(kind=8) :: pcvaac

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

    real(kind=8), dimension(:), allocatable :: amu1
    real(kind=8), dimension(:), allocatable :: amu2

    real(kind=8) :: z1
    real(kind=8) :: z2

    ! Limit on the list of triplets for which no data were found.
    data nlistl / 20 /

    ! Allocate work space arrays.
    ALLOCATE(amu1(jpfcmx))
    ALLOCATE(amu2(jpfcmx))

    ! Initialize the data arrays.
    do n = 1,naactr
        qpdaac(n) = .false.

        do j = 1,jpfcmx
            amuaac(j,n) = 0.
        end do
    end do

    ! Check the entered data for aa'c triplets.
    nodatc = 0

    do n = 1,naactr
        i = iaactr(1,n)
        j = iaactr(2,n)
        k = iaactr(3,n)
        unam1 = uaqsp(i)
        unam2 = uaqsp(j)
        unam3 = uaqsp(k)
        j2 = ilnobl(unam1)
        j3 = ilnobl(unam2)
        j4 = ilnobl(unam3)
        z1 = zaqsp(i)
        z2 = zaqsp(j)

        ! Search for unam1, unam2, unam3 in the utripl array.
        ! That array corresponds to the species triplets blocks.
        call srch33(jtripl,unam1,unam2,unam3,utripl,npx3mx,npx3t)

        if (jtripl .le. 0) then
            ! No data block was found on the DATA0 file.
            ! Note that qpdaac(n) is left with a value of .false.
            nodatc = nodatc + 1
        else
            ! Have found an entry.
            qpdaac(n) = .true.

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
1010 format(/' * Error - (EQPT/tpraac) Have found a',' duplicate data block',/7x,'on the DATA0 file'," for the aa'c triplet",/7x,a,', ',a,', ',a,'.')
                else
                    ux8 = ' '
                    write (ux8,'(i5)') ndupl
                    call lejust(ux8)
                    j6 = ilnobl(ux8)
                    write (noutpt,1020) ux8(1:j6),unam1(1:j2),unam2(1:j3),unam3(1:j4)
                    write (nttyo,1020) ux8(1:j6),unam1(1:j2),unam2(1:j3),unam3(1:j4)
1020 format(/' * Error - (EQPT/tpraac) Have found ',a,' duplicate data blocks',/7x,'on the DATA0 file'," for the aa'c triplet",/7x,a,', ',a,', ',a,'.')
                end if

                nerr = nerr + ndupl
            end if

            ! Get the data for the two constituent binary systems.
            qfound1 = .false.

            do j = 1,jpfcmx
                amu1(j) = 0.
            end do

            do nn = 1,ncapr
                ii = icapr(1,nn)

                if (uaqsp(ii)(1:24) .eq. unam3(1:24)) then
                    jj = icapr(2,nn)

                    if (uaqsp(jj)(1:24) .eq. unam1(1:24)) then
                        if (qpdca(nn)) then
                            qfound1 = .true.

                            do j = 1,jpfcmx
                                amu1(j) = amua2c(j,nn)
                            end do

                            go to 130
                        end if
                    end if
                end if
            end do

130 continue

            qfound2 = .false.

            do j = 1,jpfcmx
                amu2(j) = 0.
            end do

            do nn = 1,ncapr
                ii = icapr(1,nn)

                if (uaqsp(ii)(1:24) .eq. unam3(1:24)) then
                    jj = icapr(2,nn)

                    if (uaqsp(jj)(1:24) .eq. unam2(1:24)) then
                        if (qpdca(nn)) then
                            qfound2 = .true.

                            do j = 1,jpfcmx
                                amu2(j) = amua2c(j,nn)
                            end do

                            go to 140
                        end if
                    end if
                end if
            end do

140 continue

            if (qfound1 .and. qfound2) then
                do j = 1,jpfcmx
                    amuaac(j,n) = (apsi(j,jtripl) +        3.0*((z2/z1)*amu1(j) + (z1/z2)*amu2(j)))/6.
                end do
            else
                ! Certain required data were not found on the DATA0 file.
                write (noutpt,1030) unam1(1:j2),unam2(1:j3),unam3(1:j4)
                write (nttyo,1030) unam1(1:j2),unam2(1:j3),unam3(1:j4)
1030 format(/' * Error - (EQPT/tpraac) Have data on the DATA0'," file for the aa'c",/7x,'triplet ',a,', ',a,', ',a,', but'," don't have the",/7x,'required data for the following',' cation-anion pair(s):',/)

                if (.not.qfound1) then
                    write (noutpt,1040) unam3(1:j2),unam1(1:j4)
                    write (nttyo,1040) unam3(1:j2),unam1(1:j4)
1040 format(9x,a,', ',a)
                end if

                if (.not.qfound2) then
                    write (noutpt,1040) unam3(1:j3),unam1(1:j4)
                    write (nttyo,1040) unam3(1:j3),unam1(1:j4)
                end if

                nerr = nerr + 1
            end if
        end if
    end do

    if (nodatc .gt. 0) then
        write (noutpt,1060)
        write (nttyo,1060)
1060 format(/' * Warning - (EQPT/tpraac) Did not find a data',' block on the DATA0 file',/7x,'for any of the following'," aa'c triplets:",/)

        ncount = 0

        do n = 1,naactr
            if (.not.qpdaac(n)) then
                ncount = ncount + 1
                i = iaactr(1,n)
                j = iaactr(2,n)
                k = iaactr(3,n)
                unam1 = uaqsp(i)
                unam2 = uaqsp(j)
                unam3 = uaqsp(k)
                j2 = ilnobl(unam1)
                j3 = ilnobl(unam2)
                j4 = ilnobl(unam3)
                ustr80 = unam1(1:j2) // ', ' // unam2(1:j3) // ', ' //      unam3(1:j4)
                j5 = ilnobl(ustr80)
                write (noutpt,1070) ustr80(1:j5)
                write (nttyo,1070) ustr80(1:j5)
1070 format(9x,a)

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
            write (noutpt,1080) ux8(1:j3)
            write (nttyo,1080) ux8(1:j3)
1080 format(/9x,'plus ',a,' others')
        end if

        write (noutpt,1090)
        write (nttyo,1090)
1090 format(1x)

        nwarn = nwarn + 1
    end if

    ncvaac = naactr - nodatc
    pcvaac = (100.*ncvaac)/float(naactr)

    ! Deallocate work space arrays.
    DEALLOCATE(amu1)
    DEALLOCATE(amu2)
end subroutine tpraac