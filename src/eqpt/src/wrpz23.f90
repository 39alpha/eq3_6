subroutine wrpz23(alphca,alamaa,alamca,alamcc,alamna,alamnc,alamnn,alamn2,amuaac,amua2c,amucca,amuc2a,amunca,amun2n,amun3,ipbtmx,iaapr,icapr,iccpr,inapr,incpr,innpr,in2pr,iaactr,ia2ctr,iccatr,ic2atr,incatr,in2ntr,in3tr,jpdblo,jpfcmx,natmax,naapr,ncapr,nccpr,nnapr,nncpr,nnnpr,nn2pr,naactr,na2ctr,nccatr,nc2atr,nncatr,nn2ntr,nn3tr,ndata1,ndat1f,noutpt,nttyo,uaqsp,uethfl)
    !! This subroutine writes the processed Pitzer parameters (in
    !! conventional lambda and mu form) on the DATA1 and DATA1F files.
    !! This subroutine is called by:
    !!   EQPT/eqpt.f
    !! Principal input:
    !!   alamaa  = array of coefficients for computing lamda parameters
    !!               for aa' pairs
    !!   alamca  = array of coefficients for computing lamda parameters
    !!               for ca pairs
    !!   alamcc  = array of coefficients for computing lamda parameters
    !!               for cc' pairs
    !!   alamna  = array of coefficients for computing lamda parameters
    !!               for na pairs
    !!   alamnc  = array of coefficients for computing lamda parameters
    !!               for nc pairs
    !!   alamnn  = array of coefficients for computing lamda parameters
    !!               for nn' pairs
    !!   alamn2  = array of coefficients for computing lamda parameters
    !!               for nn pairs
    !!   amuaac  = array of coefficients for computing mu parameters
    !!               for aa'c triplets
    !!   amua2c  = array of coefficients for computing mu parameters
    !!               for aac triplets
    !!   amucca  = array of coefficients for computing mu parameters
    !!               for cc'a triplets
    !!   amuc2a  = array of coefficients for computing mu parameters
    !!               for cca triplets
    !!   amunca  = array of coefficients for computing mu parameters
    !!               for nca triplets
    !!   amun2n  = array of coefficients for computing mu parameters
    !!               for nnn' triplets
    !!   amun3   = array of coefficients for computing mu parameters
    !!               for nnn triplets
    !!   iaapr   = species index array for aa' pairs
    !!   icapr   = species index array for ca pairs
    !!   iccpr   = species index array for cc' pairs
    !!   inapr   = species index array for na pairs
    !!   incpr   = species index array for nc pairs
    !!   innpr   = species index array for nn' pairs
    !!   in2pr   = species index array for nn pairs
    !!   iaactr  = species index array for aa'c triplets
    !!   ia2ctr  = species index array for aac triplets
    !!   iccatr  = species index array for cc'a triplets
    !!   ic2atr  = species index array for cca triplets
    !!   incatr  = species index array for nca triplets
    !!   in2ntr  = species index array for nn'n triplets
    !!   naapr   = number of aa' pairs
    !!   ncapr   = number of ca pairs
    !!   nccpr   = number of cc' pairs
    !!   nnapr   = number of na pairs
    !!   nncpr   = number of nc pairs
    !!   nnnpr   = number of nn' pairs
    !!   nn2pr   = number of nn pairs
    !!   naactr  = number of aa'c triplets
    !!   na2ctr  = number of aac triplets
    !!   nccatr  = number of cc'a triplets
    !!   nc2atr  = number of cca triplets
    !!   nncatr  = number of nca triplets
    !!   nn2ntr  = number of nn'n triplets
    !!   ndata1 = unit number of the DATA1 file
    !!   ndat1f = unit number of the DATA1F file
    !!   uethfl = E-theta flag string
    !!   upair  = array of names in species pairs
    !!   utripl = array of names in species triplets
    !! Principal output:
    !!   none
    implicit none

    ! Calling sequence variable declarations.
    integer :: ipbtmx
    integer :: jpfcmx
    integer :: natmax

    integer :: ndata1
    integer :: ndat1f
    integer :: noutpt
    integer :: nttyo

    integer :: jpdblo
    integer :: naapr
    integer :: ncapr
    integer :: nccpr
    integer :: nnapr
    integer :: nncpr
    integer :: nnnpr
    integer :: nn2pr
    integer :: naactr
    integer :: na2ctr
    integer :: nccatr
    integer :: nc2atr
    integer :: nncatr
    integer :: nn2ntr
    integer :: nn3tr

    integer :: iaapr(2,naapr)
    integer :: icapr(2,ncapr)
    integer :: iccpr(2,nccpr)
    integer :: inapr(2,nnapr)
    integer :: incpr(2,nncpr)
    integer :: innpr(2,nnnpr)
    integer :: in2pr(nn2pr)
    integer :: iaactr(3,naactr)
    integer :: ia2ctr(2,na2ctr)
    integer :: iccatr(3,nccatr)
    integer :: ic2atr(2,nc2atr)
    integer :: incatr(3,nncatr)
    integer :: in2ntr(2,nn2ntr)
    integer :: in3tr(nn3tr)

    character(len=24) :: uaqsp(natmax)
    character(len=8) :: uethfl

    real(kind=8) :: alphca(ipbtmx,ncapr)

    real(kind=8) :: alamaa(jpfcmx,0:ipbtmx,naapr)
    real(kind=8) :: alamca(jpfcmx,0:ipbtmx,ncapr)
    real(kind=8) :: alamcc(jpfcmx,0:ipbtmx,nccpr)
    real(kind=8) :: alamna(jpfcmx,0:ipbtmx,nnapr)
    real(kind=8) :: alamnc(jpfcmx,0:ipbtmx,nncpr)
    real(kind=8) :: alamnn(jpfcmx,0:ipbtmx,nnnpr)
    real(kind=8) :: alamn2(jpfcmx,0:ipbtmx,nn2pr)

    real(kind=8) :: amuaac(jpfcmx,naactr)
    real(kind=8) :: amua2c(jpfcmx,na2ctr)
    real(kind=8) :: amucca(jpfcmx,nccatr)
    real(kind=8) :: amuc2a(jpfcmx,nc2atr)
    real(kind=8) :: amunca(jpfcmx,nncatr)
    real(kind=8) :: amun2n(jpfcmx,nn2ntr)
    real(kind=8) :: amun3(jpfcmx,nn3tr)

    ! Local variable declarations.
    integer :: i
    integer :: j
    integer :: jpfc
    integer :: j2
    integer :: j3
    integer :: j4
    integer :: k
    integer :: n

    integer :: ilnobl

    character(len=16) :: ux16

    logical :: qnzdat

    character(len=80) :: ustr
    character(len=24) :: ublk24
    character(len=24) :: ustp24
    character(len=24) :: uend24
    character(len=24) :: unam1
    character(len=24) :: unam2
    character(len=24) :: unam3
    character(len=72) :: uterm
    character(len=72) :: utermc

    real(kind=8) :: zero

    data uend24 / 'endit.                  '  /
    data ublk24 / '                        '  /
    data ustp24 / 'stop.                   '  /
    data zero   /    0.0    /

    uterm(1:48) = '+-----------------------------------------------'
    uterm(49:72) = '------------------------'
    utermc = uterm
    utermc(1:1) = '*'

    j3 = ilnobl(uend24)
    j4 = ilnobl(ustp24)

    if (jpdblo .eq. -1) then
        ! Old Pitzer data block organization.
        ! Write the E-lambda flag.
        ustr = 'E-lambda flag = ' // uethfl
        write (ndata1) ustr
        j2 = ilnobl(ustr)
        write (ndat1f,1010) ustr(1:j2)
1010 format('*',/a,/'*')
    else
        ! New Pitzer data block organization.
        ! Do not write the E-lambda flag. It is implicitly "on" in the
        ! case of the new Pitzer data block organization.
        continue
    end if

    ! Write the lambda set for ca pairs.
    write (ndat1f,1100)
1100 format("* Coefficients for lamda(ca) parameters")

    write (ndat1f,1110) utermc
1110 format(a72)

    do n = 1,ncapr
        ! Are there any non-zero data?
        qnzdat = .false.

        do i = 0,ipbtmx
            do jpfc = 1,jpfcmx
                if (alamca(jpfc,i,n) .ne. 0.) then
                    qnzdat = .true.
                    go to 120
                end if
            end do
        end do

120 continue

        if (qnzdat) then
            ! Write a block.
            i = icapr(1,n)
            j = icapr(2,n)
            unam1 = uaqsp(i)
            unam2 = uaqsp(j)

            write (ndata1) unam1,unam2
            j2 = ilnobl(unam2)
            write (ndat1f,1120) unam1,unam2(1:j2)
1120 format(a24,2x,a)

            if (jpdblo .eq. -1) then
                ! Old Pitzer data block organization.
                write (ndata1) alamca(1,0,n),alamca(1,1,n),alamca(1,2,n)
                write (ndat1f,1130) alamca(1,0,n),alamca(1,1,n),alamca(1,2,n)
1130 format(3x,'lambda0 = ',f9.5,3x,'lambda1 = ',f9.5,3x,'lambda2 = ',f9.5)

                write (ndata1) alphca(1,n),alphca(2,n)
                write (ndat1f,1140) alphca(1,n),alphca(2,n)
1140 format(25x,'alpha1 = ',f5.1,7x,'alpha2 = ',f5.1)

                write (ndata1) alamca(2,0,n),alamca(3,0,n)
                write (ndat1f,1150) alamca(2,0,n),alamca(3,0,n)
1150 format(4x,'dl0/dt = ',1pe10.3,2x,'d2l0/dt2 = ',1pe10.3)

                write (ndata1) alamca(2,1,n),alamca(3,1,n)
                write (ndat1f,1160) alamca(2,1,n),alamca(3,1,n)
1160 format(4x,'dl1/dt = ',1pe10.3,2x,'d2l1/dt2 = ',1pe10.3)

                write (ndata1) alamca(2,2,n),alamca(3,2,n)
                write (ndat1f,1170) alamca(2,2,n),alamca(3,2,n)
1170 format(4x,'dl2/dt = ',1pe10.3,2x,'d2l2/dt2 = ',1pe10.3)
            else
                ! New Pitzer data block organization.
                ux16 = 'alpha( )'
                j2 = 8

                do i = 1,ipbtmx
                    write (ux16(7:7),'(i1)') i
                    write (ndata1) ux16,alphca(i,n)
                    write (ndat1f,'(2x,a," = ",1pe12.5)')        ux16(1:j2),alphca(i,n)
                end do

                ux16 = 'lambda( )'
                j2 = 9

                do i = 0,ipbtmx
                    write (ux16(8:8),'(i1)') i
                    write (ndata1) ux16
                    write (ndat1f,'(2x,a,":")') ux16(1:j2)

                    do j = 1,jpfcmx
                        write (ndata1) alamca(j,i,n)
                        write (ndat1f,'(4x,"a",i1," = ",1pe15.8)')          j,alamca(j,i,n)
                    end do
                end do
            end if

            ! Write the block terminator.
            write (ndata1) uterm
            write (ndat1f,1180) uterm
1180 format(a72)
        end if
    end do

    ! Write the lambda set for nn pairs.
    write (ndat1f,1200)
1200 format("* Coefficients for lamda(nn) parameters")

    write (ndat1f,1110) utermc

    do n = 1,nn2pr
        ! Are there any non-zero data?
        qnzdat = .false.

        do i = 0,ipbtmx
            do jpfc = 1,jpfcmx
                if (alamn2(jpfc,i,n) .ne. 0.) then
                    qnzdat = .true.
                    go to 130
                end if
            end do
        end do

130 continue

        if (qnzdat) then
            ! Write a block.
            i = in2pr(n)
            unam1 = uaqsp(i)
            unam2 = unam1

            write (ndata1) unam1,unam2
            j2 = ilnobl(unam2)
            write (ndat1f,1120) unam1,unam2(1:j2)

            if (jpdblo .eq. -1) then
                ! Old Pitzer data block organization.
                write (ndata1) alamn2(1,0,n),alamn2(1,1,n),alamn2(1,2,n)
                write (ndat1f,1130) alamn2(1,0,n),alamn2(1,1,n),alamn2(1,2,n)

                ! No alpha parameters for this kind of pair.
                write (ndata1) zero,zero
                write (ndat1f,1140) zero, zero

                write (ndata1) alamn2(2,0,n),alamn2(3,0,n)
                write (ndat1f,1150) alamn2(2,0,n),alamn2(3,0,n)

                write (ndata1) alamn2(2,1,n),alamn2(3,1,n)
                write (ndat1f,1160) alamn2(2,1,n),alamn2(3,1,n)

                write (ndata1) alamn2(2,2,n),alamn2(3,2,n)
                write (ndat1f,1170) alamn2(2,2,n),alamn2(3,2,n)
            else
                ! New Pitzer data block organization.
                ux16 = 'lambda(0)'
                j2 = 9
                write (ndata1) ux16
                write (ndat1f,'(2x,a,":")') ux16(1:j2)

                do j = 1,jpfcmx
                    write (ndata1) alamn2(j,0,n)
                    write (ndat1f,'(4x,"a",i1," = ",1pe15.8)') j,alamn2(j,0,n)
                end do
            end if

            ! Write the block terminator.
            write (ndata1) uterm
            write (ndat1f,1180) uterm
        end if
    end do

    ! Write the lambda set for nc pairs.
    write (ndat1f,1210)
1210 format("* Coefficients for lamda(nc) parameters")

    write (ndat1f,1110) utermc

    do n = 1,nncpr
        ! Are there any non-zero data?
        qnzdat = .false.

        do i = 0,ipbtmx
            do jpfc = 1,jpfcmx
                if (alamnc(jpfc,i,n) .ne. 0.) then
                    qnzdat = .true.
                    go to 140
                end if
            end do
        end do

140 continue

        if (qnzdat) then
            ! Write a block.
            i = incpr(1,n)
            j = incpr(2,n)
            unam1 = uaqsp(i)
            unam2 = uaqsp(j)

            write (ndata1) unam1,unam2
            j2 = ilnobl(unam2)
            write (ndat1f,1120) unam1,unam2(1:j2)

            if (jpdblo .eq. -1) then
                ! Old Pitzer data block organization.
                write (ndata1) alamnc(1,0,n),alamnc(1,1,n),alamnc(1,2,n)
                write (ndat1f,1130) alamnc(1,0,n),alamnc(1,1,n),alamnc(1,2,n)

                ! No alpha parameters for this kind of pair.
                write (ndata1) zero,zero
                write (ndat1f,1140) zero, zero

                write (ndata1) alamnc(2,0,n),alamnc(3,0,n)
                write (ndat1f,1150) alamnc(2,0,n),alamnc(3,0,n)

                write (ndata1) alamnc(2,1,n),alamnc(3,1,n)
                write (ndat1f,1160) alamnc(2,1,n),alamnc(3,1,n)

                write (ndata1) alamnc(2,2,n),alamnc(3,2,n)
                write (ndat1f,1170) alamnc(2,2,n),alamnc(3,2,n)
            else
                ! New Pitzer data block organization.
                ux16 = 'lambda(0)'
                j2 = 9
                write (ndata1) ux16
                write (ndat1f,'(2x,a,":")') ux16(1:j2)

                do j = 1,jpfcmx
                    write (ndata1) alamnc(j,0,n)
                    write (ndat1f,'(4x,"a",i1," = ",1pe15.8)') j,alamnc(j,0,n)
                end do
            end if

            ! Write the block terminator.
            write (ndata1) uterm
            write (ndat1f,1180) uterm
        end if
    end do

    ! Write the lambda set for na pairs.
    write (ndat1f,1220)
1220 format("* Coefficients for lamda(na) parameters")

    write (ndat1f,1110) utermc

    do n = 1,nnapr
        ! Are there any non-zero data?
        qnzdat = .false.

        do i = 0,ipbtmx
            do jpfc = 1,jpfcmx
                if (alamna(jpfc,i,n) .ne. 0.) then
                    qnzdat = .true.
                    go to 150
                end if
            end do
        end do

150 continue

        if (qnzdat) then
            ! Write a block.
            i = inapr(1,n)
            j = inapr(2,n)
            unam1 = uaqsp(i)
            unam2 = uaqsp(j)

            write (ndata1) unam1,unam2
            j2 = ilnobl(unam2)
            write (ndat1f,1120) unam1,unam2(1:j2)

            if (jpdblo .eq. -1) then
                ! Old Pitzer data block organization.
                write (ndata1) alamna(1,0,n),alamna(1,1,n),alamna(1,2,n)
                write (ndat1f,1130) alamna(1,0,n),alamna(1,1,n),alamna(1,2,n)

                ! No alpha parameters for this kind of pair.
                write (ndata1) zero,zero
                write (ndat1f,1140) zero, zero

                write (ndata1) alamna(2,0,n),alamna(3,0,n)
                write (ndat1f,1150) alamna(2,0,n),alamna(3,0,n)

                write (ndata1) alamna(2,1,n),alamna(3,1,n)
                write (ndat1f,1160) alamna(2,1,n),alamna(3,1,n)

                write (ndata1) alamna(2,2,n),alamna(3,2,n)
                write (ndat1f,1170) alamna(2,2,n),alamna(3,2,n)
            else
                ! New Pitzer data block organization.
                ux16 = 'lambda(0)'
                j2 = 9
                write (ndata1) ux16
                write (ndat1f,'(2x,a,":")') ux16(1:j2)

                do j = 1,jpfcmx
                    write (ndata1) alamna(j,0,n)
                    write (ndat1f,'(4x,"a",i1," = ",1pe15.8)') j,alamna(j,0,n)
                end do
            end if

            ! Write the block terminator.
            write (ndata1) uterm
            write (ndat1f,1180) uterm
        end if
    end do

    ! Write the lambda set for nn' pairs.
    write (ndat1f,1230)
1230 format("* Coefficients for lamda(nn') parameters")

    write (ndat1f,1110) utermc

    do n = 1,nnnpr
        ! Are there any non-zero data?
        qnzdat = .false.

        do i = 0,ipbtmx
            do jpfc = 1,jpfcmx
                if (alamnn(jpfc,i,n) .ne. 0.) then
                    qnzdat = .true.
                    go to 160
                end if
            end do
        end do

160 continue

        if (qnzdat) then
            ! Write a block.
            i = innpr(1,n)
            j = innpr(2,n)
            unam1 = uaqsp(i)
            unam2 = uaqsp(j)

            write (ndata1) unam1,unam2
            j2 = ilnobl(unam2)
            write (ndat1f,1120) unam1,unam2(1:j2)

            if (jpdblo .eq. -1) then
                ! Old Pitzer data block organization.
                write (ndata1) alamnn(1,0,n),alamnn(1,1,n),alamnn(1,2,n)
                write (ndat1f,1130) alamnn(1,0,n),alamnn(1,1,n),alamnn(1,2,n)

                ! No alpha parameters for this kind of pair.
                write (ndata1) zero,zero
                write (ndat1f,1140) zero, zero

                write (ndata1) alamnn(2,0,n),alamnn(3,0,n)
                write (ndat1f,1150) alamnn(2,0,n),alamnn(3,0,n)

                write (ndata1) alamnn(2,1,n),alamnn(3,1,n)
                write (ndat1f,1160) alamnn(2,1,n),alamnn(3,1,n)

                write (ndata1) alamnn(2,2,n),alamnn(3,2,n)
                write (ndat1f,1170) alamnn(2,2,n),alamnn(3,2,n)
            else
                ! New Pitzer data block organization.
                ux16 = 'lambda(0)'
                j2 = 9
                write (ndata1) ux16
                write (ndat1f,'(2x,a,":")') ux16(1:j2)

                do j = 1,jpfcmx
                    write (ndata1) alamnn(j,0,n)
                    write (ndat1f,'(4x,"a",i1," = ",1pe15.8)') j,alamnn(j,0,n)
                end do
            end if

            ! Write the block terminator.
            write (ndata1) uterm
            write (ndat1f,1180) uterm
        end if
    end do

    write (ndata1) uend24,ublk24,ublk24
    write (ndat1f,1270) uend24(1:j3)
1270 format(a)

    ! Write the lambda set for cc' pairs.
    write (ndat1f,1320)
1320 format("* Coefficients for lamda(cc') parameters")

    write (ndat1f,1110) utermc

    do n = 1,nccpr
        ! Are there any non-zero data?
        qnzdat = .false.

        do i = 0,ipbtmx
            do jpfc = 1,jpfcmx
                if (alamcc(jpfc,i,n) .ne. 0.) then
                    qnzdat = .true.
                    go to 170
                end if
            end do
        end do

170 continue

        if (qnzdat) then
            ! Write a block.
            i = iccpr(1,n)
            j = iccpr(2,n)
            unam1 = uaqsp(i)
            unam2 = uaqsp(j)

            write (ndata1) unam1,unam2
            j2 = ilnobl(unam2)
            write (ndat1f,1120) unam1,unam2(1:j2)

            if (jpdblo .eq. -1) then
                ! Old Pitzer data block organization.
                write (ndata1) alamcc(1,0,n),alamcc(1,1,n),alamcc(1,2,n)
                write (ndat1f,1130) alamcc(1,0,n),alamcc(1,1,n),alamcc(1,2,n)

                ! No alpha parameters for this kind of pair.
                write (ndata1) zero,zero
                write (ndat1f,1140) zero, zero

                write (ndata1) alamcc(2,0,n),alamcc(3,0,n)
                write (ndat1f,1150) alamcc(2,0,n),alamcc(3,0,n)

                write (ndata1) alamcc(2,1,n),alamcc(3,1,n)
                write (ndat1f,1160) alamcc(2,1,n),alamcc(3,1,n)

                write (ndata1) alamcc(2,2,n),alamcc(3,2,n)
                write (ndat1f,1170) alamcc(2,2,n),alamcc(3,2,n)
            else
                ! New Pitzer data block organization.
                ux16 = 'lambda(0)'
                j2 = 9
                write (ndata1) ux16
                write (ndat1f,'(2x,a,":")') ux16(1:j2)

                do j = 1,jpfcmx
                    write (ndata1) alamcc(j,0,n)
                    write (ndat1f,'(4x,"a",i1," = ",1pe15.8)') j,alamcc(j,0,n)
                end do
            end if

            ! Write the block terminator.
            write (ndata1) uterm
            write (ndat1f,1180) uterm
        end if
    end do

    ! Write the lambda set for aa' pairs.
    write (ndat1f,1330)
1330 format("* Coefficients for lamda(aa') parameters")

    write (ndat1f,1110) utermc

    do n = 1,naapr
        ! Are there any non-zero data?
        qnzdat = .false.

        do i = 0,ipbtmx
            do jpfc = 1,jpfcmx
                if (alamaa(jpfc,i,n) .ne. 0.) then
                    qnzdat = .true.
                    go to 180
                end if
            end do
        end do

180 continue

        if (qnzdat) then
            ! Write a block.
            i = iaapr(1,n)
            j = iaapr(2,n)
            unam1 = uaqsp(i)
            unam2 = uaqsp(j)

            write (ndata1) unam1,unam2
            j2 = ilnobl(unam2)
            write (ndat1f,1120) unam1,unam2(1:j2)

            if (jpdblo .eq. -1) then
                ! Old Pitzer data block organization.
                write (ndata1) alamaa(1,0,n),alamaa(1,1,n),alamaa(1,2,n)
                write (ndat1f,1130) alamaa(1,0,n),alamaa(1,1,n),alamaa(1,2,n)

                ! No alpha parameters for this kind of pair.
                write (ndata1) zero,zero
                write (ndat1f,1140) zero, zero

                write (ndata1) alamaa(2,0,n),alamaa(3,0,n)
                write (ndat1f,1150) alamaa(2,0,n),alamaa(3,0,n)

                write (ndata1) alamaa(2,1,n),alamaa(3,1,n)
                write (ndat1f,1160) alamaa(2,1,n),alamaa(3,1,n)

                write (ndata1) alamaa(2,2,n),alamaa(3,2,n)
                write (ndat1f,1170) alamaa(2,2,n),alamaa(3,2,n)
            else
                ! New Pitzer data block organization.
                ux16 = 'lambda(0)'
                j2 = 9
                write (ndata1) ux16
                write (ndat1f,'(2x,a,":")') ux16(1:j2)

                do j = 1,jpfcmx
                    write (ndata1) alamaa(j,0,n)
                    write (ndat1f,'(4x,"a",i1," = ",1pe15.8)') j,alamaa(j,0,n)
                end do
            end if

            ! Write the block terminator.
            write (ndata1) uterm
            write (ndat1f,1180) uterm
        end if
    end do

    write (ndata1) uend24,ublk24,ublk24
    write (ndat1f,1270) uend24(1:j3)

    ! Write the mu set for cca triplets.
    write (ndat1f,1340)
1340 format("* Coefficients for mu(cca) parameters")

    write (ndat1f,1110) utermc

    do n = 1,nc2atr
        ! Are there any non-zero data?
        qnzdat = .false.

        do jpfc = 1,jpfcmx
            if (amuc2a(jpfc,n) .ne. 0.) then
                qnzdat = .true.
                go to 220
            end if
        end do

220 continue

        if (qnzdat) then
            ! Write a block.
            i = ic2atr(1,n)
            k = ic2atr(2,n)
            unam1 = uaqsp(i)
            unam2 = unam1
            unam3 = uaqsp(k)

            write (ndata1) unam1,unam2,unam3
            j2 = ilnobl(unam3)
            write (ndat1f,1370) unam1,unam2,unam3(1:j2)
1370 format(a24,2x,a24,2x,a)

            if (jpdblo .eq. -1) then
                ! Old Pitzer data block organization.
                write (ndata1) amuc2a(1,n),amuc2a(2,n),amuc2a(3,n)
                write (ndat1f,1380) amuc2a(1,n),amuc2a(2,n),amuc2a(3,n)
1380 format(5x,'mummx = ',f9.5,2x,'dmummx/dt = ',1pe10.3,2x,'d2mummx/dt2 = ',1pe10.3)
            else
                ! New Pitzer data block organization.
                ux16 = 'mu'
                j2 = 2
                write (ndata1) ux16
                write (ndat1f,'(2x,a,":")') ux16(1:j2)

                do j = 1,jpfcmx
                    write (ndata1) amuc2a(j,n)
                    write (ndat1f,'(4x,"a",i1," = ",1pe15.8)') j,amuc2a(j,n)
                end do
            end if

            ! Write the block terminator.
            write (ndata1) uterm
            write (ndat1f,1180) uterm
        end if
    end do

    ! Write the mu set for aac triplets.
    write (ndat1f,1430)
1430 format("* Coefficients for mu(aac) parameters")

    write (ndat1f,1110) utermc

    do n = 1,na2ctr
        ! Are there any non-zero data?
        qnzdat = .false.

        do jpfc = 1,jpfcmx
            if (amua2c(jpfc,n) .ne. 0.) then
                qnzdat = .true.
                go to 230
            end if
        end do

230 continue

        if (qnzdat) then
            ! Write a block.
            i = ia2ctr(1,n)
            k = ia2ctr(2,n)
            unam1 = uaqsp(i)
            unam2 = unam1
            unam3 = uaqsp(k)

            write (ndata1) unam1,unam2,unam3
            j2 = ilnobl(unam3)
            write (ndat1f,1370) unam1,unam2,unam3(1:j2)

            if (jpdblo .eq. -1) then
                ! Old Pitzer data block organization.
                write (ndata1) amua2c(1,n),amua2c(2,n),amua2c(3,n)
                write (ndat1f,1380) amua2c(1,n),amua2c(2,n),amua2c(3,n)
            else
                ! New Pitzer data block organization.
                ux16 = 'mu'
                j2 = 2
                write (ndata1) ux16
                write (ndat1f,'(2x,a,":")') ux16(1:j2)

                do j = 1,jpfcmx
                    write (ndata1) amua2c(j,n)
                    write (ndat1f,'(4x,"a",i1," = ",1pe15.8)') j,amua2c(j,n)
                end do
            end if

            ! Write the block terminator.
            write (ndata1) uterm
            write (ndat1f,1180) uterm
        end if
    end do

    write (ndata1) uend24,ublk24,ublk24
    write (ndat1f,1270) uend24(1:j3)

    ! Write the mu set for cc'a triplets.
    write (ndat1f,1470)
1470 format("* Coefficients for mu(cc'a) parameters")

    write (ndat1f,1110) utermc

    do n = 1,nccatr
        ! Are there any non-zero data?
        qnzdat = .false.

        do jpfc = 1,jpfcmx
            if (amucca(jpfc,n) .ne. 0.) then
                qnzdat = .true.
                go to 240
            end if
        end do

240 continue

        if (qnzdat) then
            ! Write a block.
            i = iccatr(1,n)
            j = iccatr(2,n)
            k = iccatr(3,n)
            unam1 = uaqsp(i)
            unam2 = uaqsp(j)
            unam3 = uaqsp(k)

            write (ndata1) unam1,unam2,unam3
            j2 = ilnobl(unam3)
            write (ndat1f,1370) unam1,unam2,unam3(1:j2)

            if (jpdblo .eq. -1) then
                ! Old Pitzer data block organization.
                write (ndata1) amucca(1,n),amucca(2,n),amucca(3,n)
                write (ndat1f,1380) amucca(1,n),amucca(2,n),amucca(3,n)
            else
                ! New Pitzer data block organization.
                ux16 = 'mu'
                j2 = 2
                write (ndata1) ux16
                write (ndat1f,'(2x,a,":")') ux16(1:j2)

                do j = 1,jpfcmx
                    write (ndata1) amucca(j,n)
                    write (ndat1f,'(4x,"a",i1," = ",1pe15.8)') j,amucca(j,n)
                end do
            end if

            ! Write the block terminator.
            write (ndata1) uterm
            write (ndat1f,1180) uterm
        end if
    end do

    ! Write the mu set for aa'c triplets.
    write (ndat1f,1490)
1490 format("* Coefficients for mu(aa'c) parameters")

    write (ndat1f,1110) utermc

    do n = 1,naactr
        ! Are there any non-zero data?
        qnzdat = .false.

        do jpfc = 1,jpfcmx
            if (amuaac(jpfc,n) .ne. 0.) then
                qnzdat = .true.
                go to 260
            end if
        end do

260 continue

        if (qnzdat) then
            ! Write a block.
            i = iaactr(1,n)
            j = iaactr(2,n)
            k = iaactr(3,n)
            unam1 = uaqsp(i)
            unam2 = uaqsp(j)
            unam3 = uaqsp(k)

            write (ndata1) unam1,unam2,unam3
            j2 = ilnobl(unam3)
            write (ndat1f,1370) unam1,unam2,unam3(1:j2)

            if (jpdblo .eq. -1) then
                ! Old Pitzer data block organization.
                write (ndata1) amuaac(1,n),amuaac(2,n),amuaac(3,n)
                write (ndat1f,1380) amuaac(1,n),amuaac(2,n),amuaac(3,n)
            else
                ! New Pitzer data block organization.
                ux16 = 'mu'
                j2 = 2
                write (ndata1) ux16
                write (ndat1f,'(2x,a,":")') ux16(1:j2)

                do j = 1,jpfcmx
                    write (ndata1) amuaac(j,n)
                    write (ndat1f,'(4x,"a",i1," = ",1pe15.8)') j,amuaac(j,n)
                end do
            end if

            ! Write the block terminator.
            write (ndata1) uterm
            write (ndat1f,1180) uterm
        end if
    end do

    ! Write the mu set for nnn triplets.
    write (ndat1f,1510)
1510 format("* Coefficients for mu(nnn) parameters")

    write (ndat1f,1110) utermc

    do n = 1,nn3tr
        ! Are there any non-zero data?
        qnzdat = .false.

        do jpfc = 1,jpfcmx
            if (amun3(jpfc,n) .ne. 0.) then
                qnzdat = .true.
                go to 270
            end if
        end do

270 continue

        if (qnzdat) then
            ! Write a block.
            i = in3tr(n)
            unam1 = uaqsp(i)
            unam2 = unam1
            unam3 = unam1

            write (ndata1) unam1,unam2,unam3
            j2 = ilnobl(unam3)
            write (ndat1f,1370) unam1,unam2,unam3(1:j2)

            if (jpdblo .eq. -1) then
                ! Old Pitzer data block organization.
                write (ndata1) amun3(1,n),amun3(2,n),amun3(3,n)
                write (ndat1f,1380) amun3(1,n),amun3(2,n),amun3(3,n)
            else
                ! New Pitzer data block organization.
                ux16 = 'mu'
                j2 = 2
                write (ndata1) ux16
                write (ndat1f,'(2x,a,":")') ux16(1:j2)

                do j = 1,jpfcmx
                    write (ndata1) amun3(j,n)
                    write (ndat1f,'(4x,"a",i1," = ",1pe15.8)') j,amun3(j,n)
                end do
            end if

            ! Write the block terminator.
            write (ndata1) uterm
            write (ndat1f,1180) uterm
        end if
    end do

    ! Write the mu set for nnn' triplets.
    write (ndat1f,1520)
1520 format("* Coefficients for mu(nnn') parameters")

    write (ndat1f,1110) utermc

    do n = 1,nn2ntr
        ! Are there any non-zero data?
        qnzdat = .false.

        do jpfc = 1,jpfcmx
            if (amun2n(jpfc,n) .ne. 0.) then
                qnzdat = .true.
                go to 280
            end if
        end do

280 continue

        if (qnzdat) then
            ! Write a block.
            i = in2ntr(1,n)
            k = in2ntr(2,n)
            unam1 = uaqsp(i)
            unam2 = unam1
            unam3 = uaqsp(k)

            write (ndata1) unam1,unam2,unam3
            j2 = ilnobl(unam3)
            write (ndat1f,1370) unam1,unam2,unam3(1:j2)

            if (jpdblo .eq. -1) then
                ! Old Pitzer data block organization.
                write (ndata1) amun2n(1,n),amun2n(2,n),amun2n(3,n)
                write (ndat1f,1380) amun2n(1,n),amun2n(2,n),amun2n(3,n)
            else
                ! New Pitzer data block organization.
                ux16 = 'mu'
                j2 = 2
                write (ndata1) ux16
                write (ndat1f,'(2x,a,":")') ux16(1:j2)

                do j = 1,jpfcmx
                    write (ndata1) amun2n(j,n)
                    write (ndat1f,'(4x,"a",i1," = ",1pe15.8)') j,amun2n(j,n)
                end do
            end if

            ! Write the block terminator.
            write (ndata1) uterm
            write (ndat1f,1180) uterm
        end if
    end do

    ! Write the mu set for nca triplets.
    write (ndat1f,1530)
1530 format("* Coefficients for mu(nca) parameters")

    write (ndat1f,1110) utermc

    do n = 1,nncatr
        ! Are there any non-zero data?
        qnzdat = .false.

        do jpfc = 1,jpfcmx
            if (amunca(jpfc,n) .ne. 0.) then
                qnzdat = .true.
                go to 290
            end if
        end do

290 continue

        if (qnzdat) then
            ! Write a block.
            i = incatr(1,n)
            j = incatr(2,n)
            k = incatr(3,n)
            unam1 = uaqsp(i)
            unam2 = uaqsp(j)
            unam3 = uaqsp(k)

            write (ndata1) unam1,unam2,unam3
            j2 = ilnobl(unam3)
            write (ndat1f,1370) unam1,unam2,unam3(1:j2)

            if (jpdblo .eq. -1) then
                ! Old Pitzer data block organization.
                write (ndata1) amunca(1,n),amunca(2,n),amunca(3,n)
                write (ndat1f,1380) amunca(1,n),amunca(2,n),amunca(3,n)
            else
                ! New Pitzer data block organization.
                ux16 = 'mu'
                j2 = 2
                write (ndata1) ux16
                write (ndat1f,'(2x,a,":")') ux16(1:j2)

                do j = 1,jpfcmx
                    write (ndata1) amunca(j,n)
                    write (ndat1f,'(4x,"a",i1," = ",1pe15.8)') j,amunca(j,n)
                end do
            end if

            ! Write the block terminator.
            write (ndata1) uterm
            write (ndat1f,1180) uterm
        end if
    end do

    write (ndata1) uend24,ublk24,ublk24
    write (ndat1f,1270) uend24(1:j3)

    write (ndata1) ustp24,ublk24,ublk24
    write (ndat1f,1270) ustp24(1:j4)
end subroutine wrpz23