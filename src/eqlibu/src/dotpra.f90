real(kind=8) function dotpra(array1,array2,nmax)
    !! This subroutine computes the dot product of two real*8 arrays
    !! array1 and array2.
    !! This subroutine is nearly identical to EQLIBU/ddot.f, which has
    !! the same function, but a different calling sequence.
    !! This subroutine is called by:
    !!   EQLIBU/sgeco.f
    !! Input:
    !!   array1 =  real*8 vector with nmax elements
    !!   array2 =  real*8 vector with nmax elements
    !! Output:
    !!   dotpra = the dot product
    implicit none

    ! Calling sequence variable declarations.
    integer :: nmax

    real(kind=8) :: array1(nmax)
    real(kind=8) :: array2(nmax)

    ! Local variable declarations.
    integer :: i
    integer :: ileft

    real(kind=8) :: dx

    ! Note the use of a local variable (ax) within the loop.
    ! Note also that the loop is unrolled.
    dx = 0.
    ileft = (nmax/8)*8

    do i = 1,ileft,8
        dx = dx + array1(i)*array2(i) + array1(i + 1)*array2(i + 1)  + array1(i + 2)*array2(i + 2) + array1(i + 3)*array2(i + 3)  + array1(i + 4)*array2(i + 4) + array1(i + 5)*array2(i + 5)  + array1(i + 6)*array2(i + 6) + array1(i + 7)*array2(i + 7)
    end do

    do i = ileft + 1,nmax
        dx = dx + array1(i)*array2(i)
    end do

    dotpra = dx
end function dotpra