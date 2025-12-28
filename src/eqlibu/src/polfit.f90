subroutine polfit(aamatr,cof,gmmatr,ier,ipivot,npfmax,npft,noutpt,nttyo,xvec,yvec)
    !! This subroutine fits an exact polynomial through npft points.
    !! Each point is an x,y pair, where the values of x are in the
    !! xvec array, those of y in the yvec array. The method is to solve
    !! the matrix equation (A)(c) = (y) where:
    !!           | 1  x(1)  x(1)**2  ... |
    !!           | 1  x(2)  x(2)**2  ... |
    !!   A =     | 1  x(3)  x(3)**2  ... |
    !!           | ...                   |
    !!           | 1  x(n)  x(n)**2  ... |
    !!  and solve for the coefficients (c), which give the
    !!  representation:
    !!    y = c(1) + c(2)*x + c(3)*x**2 + ... + c(n)*x**n-1
    !! It is assumed that the xvec array contains distinct values
    !! (i.e., there are no duplicates). The matrix A is the
    !! array aamatr, the matrix dimension n is npft, and the
    !! vector c is the array cof.
    !! If the xvec array was scaled using EQLIBU/scalx1.f, the resulting
    !! yvec array should be rescaled using EQLIBU/rscaly.f.
    !! This subroutine is called by:
    !!   EQPT/intrp.f
    !!   Any
    !! Input:
    !!   xvec   = the x vector
    !!   yvec   = the y vector
    !!   npft   = number of x,y points
    !! Output:
    !!   cof    = the array of coefficients fo the fitted polynomial
    !!   ier    = error flag (returned by EQLIBU/msolvr.f):
    !!              =  0  Okay
    !!              =  1  Encountered a zero matrix
    !!              =  2  Encountered a non-zero, computationally
    !!                      singular matrix
    implicit none

    ! Calling sequence variable declarations.
    integer :: npfmax

    integer :: ipivot(npfmax)
    integer :: ier
    integer :: noutpt
    integer :: npft
    integer :: nttyo

    real(kind=8) :: aamatr(npfmax,npfmax)
    real(kind=8) :: gmmatr(npfmax,npfmax)
    real(kind=8) :: cof(npfmax)
    real(kind=8) :: xvec(npfmax)
    real(kind=8) :: yvec(npfmax)

    ! Local variable declarations.
    integer :: i
    integer :: j

    logical :: qpr

    real(kind=8) :: aax

    data qpr    /.false./

    ! Initialize error flag.
    ier = 0

    ! Set up a matrix equation.
    do j = 1,npft
        cof(j) = 0.
        aax = 1.
        aamatr(j,1) = 1.

        do i = 2,npft
            aax = aax*xvec(j)
            aamatr(j,i) = aax
        end do
    end do

    ! Solve the matrix equation.
    ! Calling sequence substitutions:
    !   cof for delvec
    !   npft for kdim
    !   npfmax for kmax
    !   yvec for rhsvec
    call msolvr(aamatr,cof,gmmatr,ier,ipivot,npft,npfmax,noutpt,nttyo,qpr,yvec)
end subroutine polfit