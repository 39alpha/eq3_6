subroutine thetck(athetx,jpfcmx,na,nc,nerr,nn,noutpt,nttyo,n1,n2,n3,unam1,unam2,unam3)
    !! This suboutine checks data that were read from the S-theta
    !! fields of the DATA0 file.
    !! This suboutine is called by:
    !!   EQPT/rdpz3.f
    !! Principal input:
    !!   na     = number of anions in a species triplet
    !!   nc     = number of cations in a species triplet
    !!   nn     = number of neutral species in a species triplet
    !!   n1     = index of the first species in a triplet
    !!   n2     = index of the second species in a triplet
    !!   n3     = index of the third species in a triplet
    !!   unam1  = first name in a species triplet
    !!   unam2  = second name in a species triplet
    !!   unam3  = third name in a species triplet
    !! Principal output:
    !!   nerr   = error counter
    implicit none

    ! Calling sequence variable declarations.
    integer :: jpfcmx

    integer :: noutpt
    integer :: nttyo

    integer :: na
    integer :: nc
    integer :: nn
    integer :: nerr
    integer :: n1
    integer :: n2
    integer :: n3

    character(len=24) :: unam1
    character(len=24) :: unam2
    character(len=24) :: unam3

    real(kind=8) :: athetx(jpfcmx)

    ! Local variable declarations.
    integer :: j
    integer :: j2
    integer :: j3
    integer :: j4

    integer :: ilnobl

    logical :: qzethx

    ! Check for erroneous data in the thetax data.
    qzethx = .true.

    do j = 1,jpfcmx
        if (athetx(j) .ne. 0.) then
            qzethx = .false.
        end if
    end do

    if (.not.qzethx) then
        if (nn .eq. 3) then
            ! Note: the unique neutral is assumed to be the third
            ! species in the triplet.
            if (n1.eq.n2 .and. n3.ne.n1) then
                ! Have an nnn' combination.
                j2 = ilnobl(unam1)
                j3 = ilnobl(unam2)
                j4 = ilnobl(unam3)
                write (noutpt,1400) unam1(1:j2),unam2(1:j3),unam3(1:j4)
                write (nttyo,1400) unam1(1:j2),unam2(1:j3),unam3(1:j4)
1400 format(/' * Error - (EQPT/thetck) Have found an illegal',' data block for the',/7x,'species triplet ',a,', ',a,', ',a,' among the blocks for',/7x,'"mixture parameters".',' Have non-zero input for lambda data in the',/7x,'"theta" field. Only mu data may be specified (in',' the "psi" field)',/7x,'in a block of this type. This'," is an nnn' combination. Enter the",/7x,'corresponding'," lambda(nn') data",' in a "single-salt parameters" block,',/7x,'using the "beta0" field. Leave the'," mu(nnn') data",' in the present',/7x,'"mixture parameters" block (in the',' "psi" field). Enter the',/7x,"corresponding mu(n'n'n)",' data in a second "mixture parameters" block',/7x,'(in the "psi" field).')

                nerr = nerr + 1
            end if
        end if

        if (nn.eq.1 .and. nc.eq.1 .and. na.eq.1) then
            ! Have an nca combination.
            j2 = ilnobl(unam1)
            j3 = ilnobl(unam2)
            j4 = ilnobl(unam3)
            write (noutpt,1410) unam1(1:j2),unam2(1:j3),unam3(1:j4)
            write (nttyo,1410) unam1(1:j2),unam2(1:j3),unam3(1:j4)
1410 format(/' * Error - (EQPT/thetck) Have found an illegal',' data block for the',/7x,'species triplet ',a,', ',a,', ',a,' among the blocks for',/7x,'"mixture parameters".',' Have non-zero input for lambda data in the',/7x,'"theta"',' field. Only zeta data may be specified (in the "psi"',' field)',/7x,'in a block of this type. This is an nca',' combination. Enter the',/7x,'corresponding lambda(nc)',' and lambda(na) data in separate',/7x,'"single-salt',' parameters" blocks, using the "beta0" fields.',/7x,'Leave the zeta(nca) data in the present "mixture',' parameters"',/7x,'block (in the "psi" field).')

            nerr = nerr + 1
        end if
    end if
end subroutine thetck