subroutine chdxtz(delxi,dlxmx0,dzvc0,iodb,iopt,jordlm,kdim,kmax,km1,kord,kxt,nodbmx,noptmx,nordmx,nordzs,noutpt,nrd1mx,nsscmx,scalim,scfczs,sscrew,qmin,smp100,uzvec1,zklogu,zvec0,zvclg0)
    !! This subroutine chooses a step size and order according to the
    !! Taylor's series accuracy criterion, examining the z vector and
    !! its associated smoothed or average derivatives. Subroutine
    !! chdxtz.f performs the same function for the r vector and its
    !! associated average derivatives.
    !! This subroutine is called by:
    !!   EQ6/path.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: kmax
    integer :: nodbmx
    integer :: noptmx
    integer :: nordmx
    integer :: nrd1mx
    integer :: nsscmx

    integer :: noutpt

    integer :: iodb(nodbmx)
    integer :: iopt(noptmx)

    integer :: jordlm
    integer :: kdim
    integer :: km1
    integer :: kord
    integer :: kxt
    integer :: nordzs

    logical :: qmin

    character(len=48) :: uzvec1(kmax)

    real(kind=8) :: dzvc0(nrd1mx,kmax)
    real(kind=8) :: sscrew(nsscmx)
    real(kind=8) :: zvclg0(kmax)
    real(kind=8) :: zvec0(kmax)

    real(kind=8) :: delxi
    real(kind=8) :: dlxmx0
    real(kind=8) :: scalim
    real(kind=8) :: scfczs
    real(kind=8) :: smp100
    real(kind=8) :: zklogu

    ! Local variable declarations with global dimensions.
    integer :: isv_nordmx

    SAVE isv_nordmx

    integer, dimension(:), allocatable :: kkk

    SAVE kkk

    character(len=24), dimension(:), allocatable :: uscal

    SAVE uscal

    real(kind=8), dimension(:), allocatable :: scmax
    real(kind=8), dimension(:), allocatable :: scraw

    SAVE scmax,scraw

    ! Local variable declarations.
    integer :: i
    integer :: ip1
    integer :: j2
    integer :: j3
    integer :: k
    integer :: kcol
    integer :: kcoli

    integer :: ilnobl

    character(len=24) :: ubas
    character(len=24) :: uzt
    character(len=24) :: uzn

    real(kind=8) :: adx
    real(kind=8) :: delchg
    real(kind=8) :: dlxii
    real(kind=8) :: fctr
    real(kind=8) :: fctrmn
    real(kind=8) :: lxx0
    real(kind=8) :: mxx0
    real(kind=8) :: relchg
    real(kind=8) :: xip1

    real(kind=8) :: fctrl

    data ubas /'Basis variable          '/

    ! Allocate or reallocate local work arrays as needed.
    if (.not.ALLOCATED(kkk)) then
        ! Local work arrays are not allocated. Zero the saved
        ! array size variables. Note that only one array is tested
        ! to see if it is allocated. It is assumed that all local
        ! work arrays are either allocated or not.
        isv_nordmx = 0
    else
        ! Local work arrays are allocated. Check to see if any of the
        ! array size variables have changed. If so, deallocate
        ! the corresponding local work arrays and zero the corresponding
        ! saved size variables.
        if (nordmx .ne. isv_nordmx) then
            DEALLOCATE(kkk,uscal,scmax,scraw)
            isv_nordmx = 0
        end if
    end if

    ! At this point, the saved array size values are zero if the
    ! corresponding arrays need to be allocated.
    if (isv_nordmx .eq. 0) then
        ALLOCATE(kkk(0:nordmx))
        ALLOCATE(uscal(0:nordmx))
        ALLOCATE(scmax(0:nordmx),scraw(0:nordmx))
        isv_nordmx = nordmx
    end if

    ! Zero the contents of the local work arrays.
    do k = 0,nordmx
        kkk(k) = 0
    end do

    do k = 0,nordmx
        uscal(k) = ' '
    end do

    do k = 0,nordmx
        scmax(k) = 0.
        scraw(k) = 0.
    end do

    if (kord .le. 0) then
        nordzs = 0
        scfczs = dlxmx0/delxi
        go to 999
    end if

    ! Set up the regular choice of orders. Here i is the order, which
    ! technically runs from 0 to jordlm. Here jordlm + 1 is the hidden
    ! extra order, which is used to constrain the error in using order
    ! jordlm.
    do i = 0,jordlm
        scmax(i) = 0.
        scraw(i) = 0.
    end do

    ! For i = 0, take delxi = dlxmx0. This keeps delxi greater
    ! than or equal to dlxmx0, except within searches, as for phase
    ! boundaries.
    scmax(0) = dlxmx0/delxi
    scraw(0) = scmax(0)
    kkk(0) = 0
    uscal(0) = 'None'

    ! Process the higher orders.
    do i = 1,jordlm
        ip1 = i + 1
        dlxii = (delxi**(ip1))/fctrl(ip1)
        fctrmn = 1.e+38
        kkk(i) = 0
        uscal(i) = 'None'

        ! Constrain the step size to satisfy the tolerance on the
        ! fractional error of matrix variables.
        do kcol = 1,kdim
            qmin = kcol.ge.km1 .and. kcol.le.kxt
            lxx0 = zvclg0(kcol)

            ! Skip the constraint for a mineral with a very small mass.
            if (.not.qmin .or. lxx0.gt.zklogu) then
                adx = abs(dzvc0(ip1,kcol))
                delchg = adx*dlxii

                if (delchg .ne. 0.) then
                    ! Taylor's series in a linear quantity (e.g., mass).
                    mxx0 = zvec0(kcol)

                    if (mxx0 .ge. smp100) then
                        relchg = delchg/mxx0
                        fctr = sscrew(1)/relchg
                    else
                        relchg = 1.e+38
                        fctr = 1.e+38
                    end if

                    if (fctr .lt. fctrmn) then
                        fctrmn = fctr
                        uscal(i)(1:24) = ubas(1:24)
                        kkk(i) = kcol
                    end if
                end if
            end if
        end do

        xip1 = ip1
        scmax(i) = fctrmn**(1./xip1)
        scraw(i) = scmax(i)

        if (scmax(i) .gt. scalim) then
            scmax(i) = scalim
        end if
    end do

    if (iodb(5) .ge. 2) then
        write (noutpt,1000)
1000 format(/' --- Order and Scale Factor Data (chdxtz) ---',/)

        do i = 0,jordlm
            uzt = uscal(i)
            kcoli = kkk(i)

            if (uzt(1:24) .eq. ubas(1:24)) then
                uzn = uzvec1(kcoli)(1:24)
            else
                uzn = ' '
            end if

            write (noutpt,1010) i,scmax(i),scraw(i)
1010 format('   Order= ',i2,', scale= ',f10.5,', raw scale= ',f10.5)

            j2 = ilnobl(uzt)
            j3 = ilnobl(uzn)

            if (j3 .gt. 0) then
                write (noutpt,1020) uzt(1:j2),uzn(1:j3)
1020 format(5x,'Constrained by ',a,1x,a)
            else
                write (noutpt,1030) uzt(1:j2)
1030 format(5x,'Constrained by ',a)
            end if
        end do
    end if

    ! Choose the order that gives the highest scale factor.
    scfczs = scmax(0)
    nordzs = 0

    do i = 1,jordlm
        if (scmax(i) .gt. scfczs) then
            scfczs = scmax(i)
            nordzs = i
        end if
    end do

999 continue
end subroutine chdxtz