subroutine initlv(qarray,nmax,qvalue)
    !! This subroutine initializes the logical array qarray to qvalue
    !! over the first nmax positions. Normally, nmax would be the
    !! dimension of the 1D qarray. However, nmax could be less than the
    !! true dimension. Also, qarray could actually have more than one
    !! dimension. For example, if the array is really qa, having
    !! dimensions of i1 and i2, this subroutine could be used by calling
    !! it in the following manner: call initlv(qa,(i1*i2),qvalue).
    !! Use this subroutine with caution if nmax is not the product of the
    !! true (declared) dimensions.
    !! This subroutine is called by:
    !!   Any
    !! Input:
    !!   qarray = the array
    !!   nmax   = dimension or pseudo-dimension of array
    !!   qvalue = the initialization value (.true. or .false.)
    !! Output:
    !!   qarray = qarray, with the first nmax positions set to the
    !!            initialization value
    implicit none

    ! Calling sequence variable declarations.
    integer :: nmax

    logical :: qarray(nmax)
    logical :: qvalue

    ! Local variable declarations.
    integer :: i
    integer :: ileft

    logical :: qv

    ! Note the assignment of the value to a local variable.
    ! Note that the loop is unrolled.
    qv = qvalue
    ileft = (nmax/8)*8

    do i = 1,ileft,8
        qarray(i) = qv
        qarray(i + 1) = qv
        qarray(i + 2) = qv
        qarray(i + 3) = qv
        qarray(i + 4) = qv
        qarray(i + 5) = qv
        qarray(i + 6) = qv
        qarray(i + 7) = qv
    end do

    do i = ileft + 1,nmax
        qarray(i) = qv
    end do
end subroutine initlv