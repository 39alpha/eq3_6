subroutine lindep(aamatr,eps100,irow1,irow2,jcol1,jcol2,kmax,qldep)
    !! This subroutine determines whether the rows of a submatrix of
    !! matrix aamatr are linearly dependent or not. The submatrix runs
    !! from rows irow1 to irow2 and columns jcol1 to jcol2. If there is
    !! linear dependence, qldep is returned with a value of ".true.".
    !! This subroutine factors all rows if possible. It is assumed that
    !! the calling subroutine may use these rows if qldep is returned
    !! with a value of ".true.".
    !! This subroutine is called by:
    !!   EQ3NR/dawfix.f
    !!   EQ6/eqcalc.f
    !!   EQ6/jgibbs.f
    !! Input:
    !!   aamatr = the matrix to be tested
    !!   eps100 = 100 times the real*8 machine epsilon
    !!   irow1  = start of range of rows to examine
    !!   irow2  = end of range of rows to examine
    !!   jcol1  = start of range of columns to examine
    !!   jcol2  = end of range of columns to examine
    !!   kmax   = dimension of aamatr
    !! Output:
    !!   qldep  = flag variable
    !!              .true.  = linear dependence
    !!              .false. = no linear dependence
    implicit none

    ! Calling sequence variable declarations.
    integer :: kmax

    integer :: irow1
    integer :: irow2
    integer :: jcol1
    integer :: jcol2

    logical :: qldep

    real(kind=8) :: aamatr(kmax,kmax)
    real(kind=8) :: eps100

    ! Local variable declarations.
    integer :: i
    integer :: idim
    integer :: ipiv
    integer :: irow
    integer :: j
    integer :: jcol
    integer :: jdim
    integer :: nrowi
    integer :: nrowp

    real(kind=8) :: aax
    real(kind=8) :: avmax
    real(kind=8) :: av
    real(kind=8) :: factor

    idim = irow2 - irow1 + 1
    jdim = jcol2 - jcol1 + 1

    do i = irow1,irow2
        avmax = 0.

        do j = jcol1,jcol2
            av = abs(aamatr(i,j))

            if (av .gt. avmax) then
                avmax = av
            end if
        end do

        if (avmax .eq. 0.) then
            ! Have a row containing all zeros.
            qldep = .true.
            go to 999
        end if
    end do

    if (idim .gt. jdim) then
        ! Have linear dependence because the number of rows is greater
        ! than the number of columns.
        qldep = .true.
        go to 999
    end if

    ! Note: nrowp is the number of rows processed; irow is the
    ! index of the row currently being tested.
    nrowp = 0

    do jcol = jcol1,jcol2
        irow = irow1 + nrowp
        avmax = 0.
        ipiv = irow

        do i = irow,irow2
            av = abs(aamatr(i,jcol))

            if (av .gt. avmax) then
                avmax = av
                ipiv = i
            end if
        end do

        if (avmax .gt. 0.) then
            ! The present column is non-empty for the present and following
            ! rows. This suffices to establish the existence of another
            ! independent row.
            nrowp = nrowp + 1

            if (ipiv .ne. irow) then
                ! Pivot- exchange rows irow and ipiv.
                do j = jcol,jcol2
                    aax = aamatr(irow,j)
                    aamatr(irow,j) = aamatr(ipiv,j)
                    aamatr(ipiv,j) = aax
                end do
            end if

            ! Move down the present column from the current row to the last
            ! row. If a non-zero entry is found, multiply the row by the
            ! factor necessary to set this entry to unity. Note that
            ! entries in any preceding columns of the submatrix are zero
            ! at this point.
            do i = irow,irow2
                if (aamatr(i,jcol) .ne. 0.) then
                    factor = 1./aamatr(i,jcol)

                    do j = jcol,jcol2
                        aamatr(i,j) = factor*aamatr(i,j)
                    end do
                end if
            end do

            do i = irow + 1,irow2
                if (aamatr(i,jcol) .ne. 0.) then
                    ! Have a following row whose entry in the current column is
                    ! not zero (hence it has a value of unity). Make it zero by
                    ! subtracting the irow-th row (whose entry in the current
                    ! column is also unity) from the following row now under
                    ! consideration.
                    do j = jcol,jcol2
                        aamatr(i,j) = aamatr(i,j) - aamatr(irow,j)

                        if (abs(aamatr(i,j)) .le. eps100) then
                            aamatr(i,j) = 0.
                        end if
                    end do
                end if
            end do
        end if
    end do

    ! Count the number of independent rows.
    nrowi = 0

    do i = irow1,irow2
        do j = jcol1,jcol2
            if (aamatr(i,j) .ne. 0.) then
                nrowi = nrowi + 1
                go to 100
            end if
        end do

100 continue
    end do

    qldep = nrowi .lt. idim

999 continue
end subroutine lindep