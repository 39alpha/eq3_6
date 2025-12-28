subroutine bldspc(iaapr,icapr,iccpr,inapr,incpr,innpr,in2pr,iaactr,ia2ctr,iccatr,ic2atr,incatr,in2ntr,in3tr,jassan,jassca,jassne,nat,natmax,naapr,ncapr,nccpr,nnapr,nncpr,nnnpr,nn2pr,naactr,na2ctr,nncatr,nccatr,nc2atr,nn2ntr,nn3tr,uaqsp,zaqsp)
    !! Build index arrays for the pair and triplet combinations of
    !! aqueous solute species for use with Pitzer parameters.
    !! All pairs and triplets are distinct in the sense that the
    !! ordering of the members is not relevant; e.g., ca = ac,
    !! cc' = c'c, nca = anc = can = cna = acn = nac.
    !! This subroutine is called by:
    !!   EQPT/eqpt.f
    !! Principal input:
    !!   naapr  = number of aa' pairs
    !!   ncapr  = number of ca pairs
    !!   nccpr  = number of cc' pairs
    !!   nnapr  = number of na pairs
    !!   nncpr  = number of nc pairs
    !!   nnnpr  = number of nn' pairs
    !!   nn2pr  = number of nn pairs
    !!   naactr = number of aa'c triplets
    !!   na2ctr = number of aac triplets
    !!   nccatr = number of cc'a triplets
    !!   nc2atr = number of cca triplets
    !!   nncatr = number of nca triplets
    !!   nn3tr = number of nnn triplets
    !!   nn2ntr = number of nnn' triplets
    !!   nat    = the number of aqueous species
    !!   uaqsp  = array of names of aqueous species
    !!   zaqsp  = array of electrical charages of aqueous species
    !! Principal output:
    !!   iaapr  = array of indices of species in aa' pairs
    !!   icapr  = array of indices of species in ca pairs
    !!   iccpr  = array of indices of species in cc' pairs
    !!   inapr  = array of indices of species in na pairs
    !!   incpr  = array of indices of species in nc pairs
    !!   innpr  = array of indices of species in nn' pairs
    !!   in2pr  = array of indices of species in nn pairs
    !!   iaactr = array of indices of species in aa'c triplets
    !!   ia2ctr = array of indices of species in aac triplets
    !!   iccatr = array of indices of species in cc'a triplets
    !!   ic2atr = array of indices of species in cca triplets
    !!   incatr = array of indices of species in nca triplets
    !!   in3tr = array of indices of species in nnn triplets
    !!   in2ntr = array of indices of species in nnn' triplets
    implicit none

    ! Calling sequence variable declarations.
    integer :: natmax

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

    integer :: jassan
    integer :: jassca
    integer :: jassne
    integer :: nat

    character(len=24) :: uaqsp(natmax)

    real(kind=8) :: zaqsp(natmax)

    ! Local variable declarations.
    integer, dimension(:), allocatable :: iworka
    integer, dimension(:), allocatable :: iworkc
    integer, dimension(:), allocatable :: iworkn

    character(len=24), dimension(:), allocatable :: uworka
    character(len=24), dimension(:), allocatable :: uworkc
    character(len=24), dimension(:), allocatable :: uworkn

    integer :: i
    integer :: j
    integer :: k
    integer :: jn
    integer :: jc
    integer :: ja
    integer :: n
    integer :: npr

    ! Allocate work arrays.
    ALLOCATE(iworkn(jassne))
    ALLOCATE(iworkc(jassca))
    ALLOCATE(iworka(jassan))

    ALLOCATE(uworkn(jassne))
    ALLOCATE(uworkc(jassca))
    ALLOCATE(uworka(jassan))

    ! Build scratch lists for cations, anions, and neutrals.
    ! Do not include solvent water or O2(g) in the neutrals list.
    ! Do not include e- on the anions list.
    jn = 0
    jc = 0
    ja = 0

    do n = 2,nat
        if (zaqsp(n) .eq. 0.) then
            if (uaqsp(n)(1:6) .ne. 'O2(g) ') then
                jn = jn + 1
                uworkn(jn) = uaqsp(n)
                iworkn(jn) = n
            end if
        else if (zaqsp(n) .gt. 0.) then
            jc = jc + 1
            uworkc(jc) = uaqsp(n)
            iworkc(jc) = n
        else
            if (uaqsp(n)(1:3) .ne. 'e- ') then
                ja = ja + 1
                uworka(ja) = uaqsp(n)
                iworka(ja) = n
            end if
        end if
    end do

    ! Construct the index array for ca pairs.
    n = 0

    do i = 1,jc
        do j = 1,ja
            n = n + 1
            icapr(1,n) = iworkc(i)
            icapr(2,n) = iworka(j)
        end do
    end do

    ! Construct the index array for cc' pairs (c' != c).
    n = 0

    do i = 1,jc
        do j = i + 1,jc
            n = n + 1

            ! Store the pair in alphabetical order.
            if (uworkc(i) .lt. uworkc(j)) then
                iccpr(1,n) = iworkc(i)
                iccpr(2,n) = iworkc(j)
            else
                iccpr(1,n) = iworkc(j)
                iccpr(2,n) = iworkc(i)
            end if
        end do
    end do

    ! Construct the index array for anion-distinct anion pairs.
    ! Ignore pairs where an anion is paired with itself, such as
    ! Cl-, Cl-, because the associated Pitzer coefficients are
    ! either implicitly undefined or explicitly set to zero.
    ! Ignore pairs where an anion is paired with itself.
    n = 0

    do i = 1,ja
        do j = i + 1,ja
            n = n + 1

            ! Store the pair in alphabetical order.
            if (uworka(i) .lt. uworka(j)) then
                iaapr(1,n) = iworka(i)
                iaapr(2,n) = iworka(j)
            else
                iaapr(1,n) = iworka(j)
                iaapr(2,n) = iworka(i)
            end if
        end do
    end do

    ! Construct the index array for nn pairs.
    n = 0

    do i = 1,jn
        n = n + 1

        ! Store the "pair."
        in2pr(n) = iworkn(i)
    end do

    ! Construct the index array for nn' pairs.
    n = 0

    do i = 1,jn
        do j = i + 1,jn
            n = n + 1

            ! Store the pair in alphabetical order.
            if (uworkn(i) .le. uworkn(j)) then
                innpr(1,n) = iworkn(i)
                innpr(2,n) = iworkn(j)
            else
                innpr(1,n) = iworkn(j)
                innpr(2,n) = iworkn(i)
            end if
        end do
    end do

    ! Construct the index array for nc pairs.
    n = 0

    do i = 1,jn
        do j = 1,jc
            n = n + 1
            incpr(1,n) = iworkn(i)
            incpr(2,n) = iworkc(j)
        end do
    end do

    ! Construct the index array for na pairs.
    n = 0

    do i = 1,jn
        do j = 1,ja
            n = n + 1
            inapr(1,n) = iworkn(i)
            inapr(2,n) = iworka(j)
        end do
    end do

    ! Construct the index array for cca triplets. Note that we will
    ! generally just use the icapr array whenever the ic2atr array
    ! is required.
    do n = 1,nc2atr
        ic2atr(1,n) = icapr(1,n)
        ic2atr(2,n) = icapr(2,n)
    end do

    ! Construct the index array for aac triplets. Note that we will
    ! generally just use the icapr array whenever the ia2ctr array
    ! is required.
    do n = 1,na2ctr
        ia2ctr(1,n) = icapr(2,n)
        ia2ctr(2,n) = icapr(1,n)
    end do

    ! Construct the index array for nnn triplets. Note that we will
    ! generally just use the in2pr array whenever the ic2atr array
    ! is required.
    do n = 1,nn3tr
        in3tr(n) = in2pr(n)
    end do

    ! Construct the index array for nnn' triplets. Note that we could
    ! just use the innpr array (with some manipulation to account for
    ! the fact that there are two triplets nnn' and n'n'n for each pair
    ! nn') whenever the in2ntr array is required. See comments above
    ! for the ic2atr array. Also note that nn2ntr = 2*nn2pr.
    n = 0

    do npr = 1,nn2pr
        ! Store the two triplet combinations corresponding to a given
        ! pair together.
        n = n + 1
        in2ntr(1,n) = innpr(1,npr)
        in2ntr(2,n) = innpr(2,npr)
        n = n + 1
        in2ntr(1,n) = innpr(2,npr)
        in2ntr(2,n) = innpr(1,npr)
    end do

    ! Construct the index array for cc'a triplets.
    n = 0

    do i = 1,jc
        do j = i + 1,jc
            do k = 1,ja
                n = n + 1

                ! Store the two cations in alphabetical order.
                if (uworkc(i) .le. uworkc(j)) then
                    iccatr(1,n) = iworkc(i)
                    iccatr(2,n) = iworkc(j)
                else
                    iccatr(1,n) = iworkc(j)
                    iccatr(2,n) = iworkc(i)
                end if

                iccatr(3,n) = iworka(k)
            end do
        end do
    end do

    ! Construct the index array for aac triplets.
    n = 0

    do i = 1,ja
        do j = i + 1,ja
            do k = 1,jc
                n = n + 1

                ! Store the two anions in alphabetical order.
                if (uworka(i) .le. uworka(j)) then
                    iaactr(1,n) = iworka(i)
                    iaactr(2,n) = iworka(j)
                else
                    iaactr(1,n) = iworka(j)
                    iaactr(2,n) = iworka(i)
                end if

                iaactr(3,n) = iworkc(k)
            end do
        end do
    end do

    ! Construct the index array for nca triplets.
    n = 0

    do i = 1,jn
        do j = 1,jc
            do k = 1,ja
                n = n + 1
                incatr(1,n) = iworkn(i)
                incatr(2,n) = iworkc(j)
                incatr(3,n) = iworka(k)
            end do
        end do
    end do

    ! Deallocate work arrays.
    DEALLOCATE(iworka)
    DEALLOCATE(iworkc)
    DEALLOCATE(iworkn)

    DEALLOCATE(uworka)
    DEALLOCATE(uworkc)
    DEALLOCATE(uworkn)
end subroutine bldspc