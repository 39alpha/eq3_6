subroutine chdxgr(delxi,dlxmx0,fdri0,fdrr0,iodb,jordlm,jreac,kord,nodbmx,nordmx,nordr,noutpt,nrct,nrctmx,nrd1mx,nsscmx,scalim,scfcr,sscrew,qriinf,rirec0,rrelr0,ureac)
    !! This subroutine chooses a step size and order according to the
    !! Gear accuracy criterion, examining the r vector and its
    !! associated finite differences. Subroutine chdxgz.f performs
    !! the same function for the z vector and its associated finite
    !! differences.
    !! This subroutine is called by:
    !!   EQ6/path.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: nodbmx
    integer :: nrctmx
    integer :: nordmx
    integer :: nrd1mx
    integer :: nsscmx

    integer :: noutpt

    integer :: iodb(nodbmx)
    integer :: jreac(nrctmx)

    integer :: jordlm
    integer :: kord
    integer :: nordr
    integer :: nrct

    logical :: qriinf

    character(len=24) :: ureac(nrctmx)

    real(kind=8) :: fdri0(nrd1mx)
    real(kind=8) :: fdrr0(nrd1mx,nrctmx)
    real(kind=8) :: rrelr0(nrctmx)
    real(kind=8) :: sscrew(nsscmx)

    real(kind=8) :: delxi
    real(kind=8) :: dlxmx0
    real(kind=8) :: rirec0
    real(kind=8) :: scalim
    real(kind=8) :: scfcr

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
    integer :: nrc

    integer :: ilnobl

    character(len=24) :: uinv
    character(len=24) :: urel
    character(len=24) :: uzt
    character(len=24) :: uzn

    real(kind=8) :: adx
    real(kind=8) :: delchg
    real(kind=8) :: dlxii
    real(kind=8) :: fctr
    real(kind=8) :: fctrmn
    real(kind=8) :: relchg
    real(kind=8) :: rxx0
    real(kind=8) :: xip1

    data uinv /'Inverse rate            '/,urel /'Relative rate           '/

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
    do k = 1,nordmx
        kkk(k) = 0
    end do

    do k = 1,nordmx
        uscal(k) = ' '
    end do

    do k = 1,nordmx
        scmax(k) = 0.
        scraw(k) = 0.
    end do

    if (kord .le. 0) then
        nordr = 0
        scfcr= dlxmx0/delxi
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
        dlxii = delxi**(ip1)
        fctrmn = 1.e+38
        kkk(i) = 0
        uscal(i) = 'None'

        ! Get results based on the inverse rate.
        if (.not.qriinf) then
            rxx0 = rirec0

            if (rxx0 .ne. 0.) then
                adx = abs(fdri0(ip1))
                delchg = adx*dlxii
                relchg = delchg/rxx0

                if (relchg .ne. 0.) then
                    fctr = sscrew(3)/relchg

                    if (fctr .lt. fctrmn) then
                        fctrmn = fctr
                        uscal(i)(1:24) = uinv(1:24)
                        kkk(i) = 0
                    end if
                end if
            end if
        end if

        ! Get results based on relative rates.
        do nrc = 1,nrct
            if (jreac(nrc).eq.0 .or. jreac(nrc).eq.-1)  then
                rxx0 = abs(rrelr0(nrc))

                if (rxx0 .ne. 0.) then
                    adx = abs(fdrr0(ip1,nrc))
                    delchg = adx*dlxii

                    if (delchg .ne. 0.) then
                        fctr = sscrew(3)/delchg

                        if (fctr .lt. fctrmn) then
                            fctrmn = fctr
                            uscal(i)(1:24) = urel(1:24)
                            kkk(i) = nrc
                        end if
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
1000 format(/' --- Order and Scale Factor Data (chdxgr) ---',/)

        do i = 0,jordlm
            uzt = uscal(i)
            nrc = kkk(i)

            if (uzt(1:24) .eq. urel(1:24)) then
                uzn = ureac(nrc)
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
    scfcr = scmax(0)
    nordr = 0

    do i = 1,jordlm
        if (scmax(i) .gt. scfcr) then
            scfcr = scmax(i)
            nordr = i
        end if
    end do

999 continue
end subroutine chdxgr