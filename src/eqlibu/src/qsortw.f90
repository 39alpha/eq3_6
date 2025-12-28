subroutine qsortw(asort,aval,istack,jsort,jstack,nmax,noutpt,nttyo,nval)
    !! This subroutine sorts the first nval elements of the array
    !! aval, using the "quicksort" method of C.A.R. Hoare. The
    !! sorted elements of aval are placed in the array asort,
    !! and the corresponding indices are placed in the array jsort.
    !! If jsort(1) is 0, the sort starts from scratch. Otherwise,
    !! the jsort array is presumed to contain a prior result which
    !! should be a good starting point for the present sorting. Arrays
    !! istack and jstack are scratch arrays. Here they are used
    !! as stacks containing the left and right indices, respectively,
    !! of lists to be sorted. All the arrays are dimensioned 1:nmax.
    !! EQLIBU/ssortw.f is a possible alternative to this subroutine.
    !! It uses the shell sort method.
    !! Caution: to be safe, set the first element of the actual array
    !! referenced here as jsort to 0 before calling this routine. If
    !! this is not done, be sure that nothing has changed the actual
    !! array since the previous call for this array. In particular,
    !! the length of the array (here represented as nmax should be the
    !! same. If it has changed, there is a very good possibility that
    !! one or more of the elements of this array are invalid. The present
    !! subroutine does only minimal validity checking. Note that with
    !! regard to previous calls, include calls to either qsortw.f or
    !! ssortw.f, as they are interchangeable in function. Note that
    !! even if only one of these subroutines is employed, it is generally
    !! used to sort various arrays. Thus, it is not a useful check to
    !! save the previous value of nmax for comparison.
    !! This subroutine is called by:
    !!   EQLIB/ncmpex.f
    !!   EQLIB/sortsp.f
    !! Input:
    !!   aval   = array of values to be sorted
    !!   istack = scratch integer array
    !!   jstack = scratch integer array
    !!   nmax   = dimension of the asort, aval, istack, jsort, and
    !!            jstack arrays
    !!   nval   = the number of elements in the aval array to be sorted,
    !!            beginning with the first element
    !! Output:
    !!   asort  = array of sorted values, corresponding to the range
    !!            composed of the first nval elements of the aval array
    !!   jsort  = array of sorted indices, corresponding to the values
    !!            in the asort array
    implicit none

    ! Calling sequence variable declarations:
    integer :: nmax

    integer :: noutpt
    integer :: nttyo

    integer :: istack(nmax)
    integer :: jsort(nmax)
    integer :: jstack(nmax)
    integer :: nval

    real(kind=8) :: asort(nmax)
    real(kind=8) :: aval(nmax)

    ! Local variable declarations.
    integer :: i
    integer :: iend
    integer :: ileft
    integer :: ipt
    integer :: jend
    integer :: jpt
    integer :: jx
    integer :: j2
    integer :: j3
    integer :: nstack

    integer :: ilnobl
    integer :: jiarmn
    integer :: jiarmx

    logical :: qx

    character(len=8) :: ux8
    character(len=8) :: ux8a

    real(kind=8) :: aax

    if (jsort(1) .ne. 0) then
        ! Perform minimal checking of the incoming jsort array.
        qx = .false.
        i = jiarmn(jsort,nval)

        if (i .ne. 1) then
            write (ux8,'(i5)') i
            call lejust(ux8)
            j2 = ilnobl(ux8)
            write (noutpt,1000) ux8(1:j2)
            write (nttyo,1000) ux8(1:j2)
1000 format(/' * Warning - (EQLIBU/qsortw) Programming error',' trap: The minimum',/7x,'element of the incoming jsort',' array is ',a,', not 1 as is required.',/7x,'Will ignore',' this jsort array and create a new one from scratch.')

            qx = .true.
        end if

        i = jiarmx(jsort,nval)

        if (i .ne. nval) then
            write (ux8,'(i5)') i
            call lejust(ux8)
            j2 = ilnobl(ux8)
            write (ux8a,'(i5)') nval
            call lejust(ux8a)
            j3 = ilnobl(ux8a)
            write (noutpt,1010) ux8(1:j2),ux8a(1:j3)
            write (nttyo,1010) ux8(1:j2),ux8a(1:j3)
1010 format(/' * Warning - (EQLIBU/qsortw) Programming error',' trap: The maximum',/7x,'element of the incoming jsort',' array is ',a,', not nval (',a,'),',/7x,'as is required.',' Will ignore this jsort array and create a new one',/7x,'from scratch.')

            qx = .true.
        end if

        if (qx) then
            jsort(1) = 0
        end if
    end if

    ! Set up the entire array as one list.
    nstack = 1
    istack(1) = 1
    jstack(1) = nval
    ileft = (nval/8)*8

    if (jsort(1) .eq. 0) then
        ! Start from scratch. Note that the loop is unrolled.
        do jpt = 1,ileft,8
            jsort(jpt) = jpt
            asort(jpt) = aval(jpt)
            jsort(jpt + 1) = jpt + 1
            asort(jpt + 1) = aval(jpt + 1)
            jsort(jpt + 2) = jpt + 2
            asort(jpt + 2) = aval(jpt + 2)
            jsort(jpt + 3) = jpt + 3
            asort(jpt + 3) = aval(jpt + 3)
            jsort(jpt + 4) = jpt + 4
            asort(jpt + 4) = aval(jpt + 4)
            jsort(jpt + 5) = jpt + 5
            asort(jpt + 5) = aval(jpt + 5)
            jsort(jpt + 6) = jpt + 6
            asort(jpt + 6) = aval(jpt + 6)
            jsort(jpt + 7) = jpt + 7
            asort(jpt + 7) = aval(jpt + 7)
        end do

        do jpt = ileft + 1,nval
            jsort(jpt) = jpt
            asort(jpt) = aval(jpt)
        end do
    else
        ! Use the pre-existing jsort array as a starting point.
        ! Note that the loop is unrolled.
        do jpt = 1,ileft,8
            asort(jpt) = aval(jsort(jpt))
            asort(jpt + 1) = aval(jsort(jpt + 1))
            asort(jpt + 2) = aval(jsort(jpt + 2))
            asort(jpt + 3) = aval(jsort(jpt + 3))
            asort(jpt + 4) = aval(jsort(jpt + 4))
            asort(jpt + 5) = aval(jsort(jpt + 5))
            asort(jpt + 6) = aval(jsort(jpt + 6))
            asort(jpt + 7) = aval(jsort(jpt + 7))
        end do

        do jpt = ileft + 1,nval
            asort(jpt) = aval(jsort(jpt))
        end do
    end if

    ! The label below is a return point for sorting more lists.
    ! If the stack of lists is empty the sort is finished.
120 continue
    if (nstack .le. 0) then
        go to 999
    end if

    ! Remove the last list from the stack.
    iend = istack(nstack)
    jend = jstack(nstack)
    nstack = nstack - 1
    ipt = iend
    jpt = jend

    ! Sort that list. The label below is a return point to continue
    ! sorting it.
130 continue
    if (ipt .ge. jpt) then
        go to 180
    end if

    ! Move the left pointer as far as possible.
140 continue
    if (asort(ipt) .gt. asort(jpt)) then
        go to 150
    end if

    ipt = ipt + 1

    if (ipt .ge. jpt) then
        go to 180
    end if

    go to 140

    ! Have found an out of order element. Exchange elements.
150 continue
    aax = asort(ipt)
    asort(ipt) = asort(jpt)
    asort(jpt) = aax
    jx = jsort(ipt)
    jsort(ipt) = jsort(jpt)
    jsort(jpt) = jx

    ! Move the right pointer as far as possible.
160 continue
    if (asort(ipt) .gt. asort(jpt)) then
        go to 170
    end if

    jpt = jpt - 1

    if (ipt .ge. jpt) then
        go to 180
    end if

    go to 160

    ! Have found an out of order element. Exchange elements.
170 continue
    aax = asort(ipt)
    asort(ipt) = asort(jpt)
    asort(jpt) = aax
    jx = jsort(ipt)
    jsort(ipt) = jsort(jpt)
    jsort(jpt) = jx
    go to 130

    ! Here jpt is the dividing point. Split the list in two, and
    ! sort each independently. Do not sort empty lists.
180 continue
    jx = jpt - 1

    if (jx .gt. iend) then
        nstack = nstack + 1
        istack(nstack) = iend
        jstack(nstack) = jx
    end if

    jx = jpt + 1

    if (jx .ge. jend) then
        go to 120
    end if

    nstack = nstack + 1
    istack(nstack) = jx
    jstack(nstack) = jend
    go to 120

999 continue
end subroutine qsortw