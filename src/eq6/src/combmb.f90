subroutine combmb(cdrs,iindx1,ipndx1,jflag,kbt,kdim,km1,kmax,kmt,kx1,kxt,mtb,mtbaq,ndrsmx,nbasp,nbt,nbtmax,ndrs,ndrsr,noutpt,nstmax,nttyo,uspec,uzvec1,zvclg1,zvec1)
    !! This subroutine combines mass balances totals so that the total
    !! mass of an active auxiliary basis species whose jflag value
    !! is 30 is combined into the total masses of the other basis
    !! species which appear in its reaction. This species is then
    !! removed from the active basis set.
    !! This subroutine is called by:
    !!   EQ6/eq6.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: kmax
    integer :: ndrsmx
    integer :: nbtmax
    integer :: nstmax

    integer :: noutpt
    integer :: nttyo

    integer :: iindx1(kmax)
    integer :: ipndx1(kmax)
    integer :: jflag(nstmax)
    integer :: nbasp(nbtmax)
    integer :: ndrs(ndrsmx)
    integer :: ndrsr(2,nstmax)

    integer :: kbt
    integer :: kdim
    integer :: km1
    integer :: kmt
    integer :: kx1
    integer :: kxt
    integer :: nbt

    character(len=48) :: uspec(nstmax)
    character(len=48) :: uzvec1(kmax)

    real(kind=8) :: cdrs(ndrsmx)
    real(kind=8) :: mtb(nbtmax)
    real(kind=8) :: mtbaq(nbtmax)
    real(kind=8) :: zvclg1(kmax)
    real(kind=8) :: zvec1(kmax)

    ! Local variable declarations with global dimensions.
    integer :: isv_kmax

    SAVE isv_kmax

    integer, dimension(:), allocatable :: iindxs
    integer, dimension(:), allocatable :: ipndxs

    SAVE iindxs,ipndxs

    real(kind=8), dimension(:), allocatable :: zvclgs

    SAVE zvclgs

    ! Local variable declarations.
    integer :: jlen
    integer :: k
    integer :: kbts
    integer :: km1s
    integer :: kmts
    integer :: krow
    integer :: kx1s
    integer :: kxts
    integer :: n
    integer :: nb
    integer :: nbb
    integer :: nerr
    integer :: np
    integer :: nr1
    integer :: nr2
    integer :: ns
    integer :: nss
    integer :: nt

    integer :: nbasis

    character(len=56) :: uspn56

    real(kind=8) :: cx
    real(kind=8) :: mx
    real(kind=8) :: lx

    real(kind=8) :: texp

    ! Allocate or reallocate local work arrays as needed.
    if (.not.ALLOCATED(iindxs)) then
        ! Local work arrays are not allocated. Zero the saved
        ! array size variables. Note that only one array is tested
        ! to see if it is allocated. It is assumed that all local
        ! work arrays are either allocated or not.
        isv_kmax = 0
    else
        ! Local work arrays are allocated. Check to see if any of the
        ! array size variables have changed. If so, deallocate
        ! the corresponding local work arrays and zero the corresponding
        ! saved size variables.
        if (kmax .ne. isv_kmax) then
            DEALLOCATE(iindxs,ipndxs,zvclgs)
            isv_kmax = 0
        end if
    end if

    ! At this point, the saved array size values are zero if the
    ! corresponding arrays need to be allocated.
    if (isv_kmax .eq. 0) then
        ALLOCATE(iindxs(kmax),ipndxs(kmax))
        ALLOCATE(zvclgs(kmax))
        isv_kmax = kmax
    end if

    ! Zero the contents of the local work arrays.
    do k = 1,kmax
        iindxs(k) = 0
        ipndxs(k) = 0
    end do

    do k = 1,kmax
        zvclgs(k) = 0.
    end do

    nerr = 0
    kbts = kbt
    km1s = km1
    kmts = kmt
    kx1s = kx1
    kxts = kxt

    do krow = 1,kdim
        iindxs(krow) = iindx1(krow)
        ipndxs(krow) = ipndx1(krow)
        zvclgs(krow) = zvclg1(krow)
    end do

    kbt = 0

    do krow = 1,kbts
        nb = iindxs(krow)
        np = ipndxs(krow)
        ns = nbasp(nb)

        if (jflag(ns) .eq. 0) then
            kbt = kbt + 1
            uzvec1(kbt) = uspec(ns)
            iindx1(kbt) = nb
            ipndx1(kbt) = np
            zvclg1(kbt) = zvclgs(krow)
        else
            nr1 = ndrsr(1,ns)
            nr2 = ndrsr(2,ns)
            nt = nr2 - nr1 + 1

            if (nt .lt. 2) then
                ! Calling sequence substitutions:
                !   uspec(ns) for unam48
                call fmspnm(jlen,uspec(ns),uspn56)
                write (noutpt,1000) uspn56(1:jlen),jflag(ns)
                write (nttyo,1000) uspn56(1:jlen),jflag(ns)
1000 format(/' * Error- (combmb) The species ',a,/7x,'has jflag value of ',i2,", but it's a strict basis",' species.',/7x,"Its mass balance can't be combined into",' that of another basis species.')

                nerr = nerr + 1
                go to 100
            end if

            do n = nr1 + 1,nr2
                cx = -cdrs(n)/cdrs(nr1)
                nss = ndrs(n)

                ! Calling sequence substitutions:
                !   nss for ns
                nbb = nbasis(nbasp,nbt,nbtmax,nss)
                mtb(nbb) = mtb(nbb) + cx*mtb(nb)
                mtbaq(nbb) = mtbaq(nbb) + cx*mtbaq(nb)
            end do

            mtb(nb) = 0.
            mtbaq(nb) = 0.
        end if

100 continue
    end do

    kdim = kbt
    km1 = kdim + 1
    kmt = kdim

    if (kmts .ge. km1s) then
        do krow = km1s,kmts
            kmt = kmt + 1
            ns = iindxs(krow)
            uzvec1(kmt) = uspec(ns)
            iindx1(kmt) = ns
            ipndx1(kmt) = ipndxs(krow)
            zvclg1(kmt) = zvclgs(krow)
        end do

        kdim = kmt
    end if

    kx1 = kdim + 1
    kxt = kdim

    if (kxts .ge. kx1s) then
        do krow = kx1s,kxts
            kxt = kxt + 1
            ns = iindxs(krow)
            uzvec1(kxt) = uspec(ns)
            iindx1(kxt) = ns
            ipndx1(kxt) = ipndxs(krow)
            zvclg1(kxt) = zvclgs(krow)
        end do

        kdim = kxt
    end if

    do krow = 1,kdim
        lx = zvclg1(krow)
        mx = texp(lx)
        zvec1(krow) = mx
    end do

    if (nerr .gt. 0) then
        stop
    end if
end subroutine combmb