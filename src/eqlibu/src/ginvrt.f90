subroutine ginvrt(aimatr,delvec,gmmatr,ipivot,kdim,kmax)
    !! This subroutine inverts the matrix aamatr. The inverted matrix
    !! is returned in aimatr. Thus
    !!   [aamatr][aimatr] = I
    !! where I is the unit matrix (1's on the diagonal, 0's elsewhere).
    !! This subroutine assumes that aamatr has already been factored
    !! into its L-U decomposition, gmmatr. The inverse matrix aimatr
    !! is obtained by solving a series of matrix equations, each using
    !! the same decomposed matrix but with a different unit vector as
    !! the right hand side. The original aamatr matrix itself is not
    !! used here.
    !! A companion subroutine minvrt.f has the same function as the
    !! present routine. However, it does its own factoring of aamatr
    !! to obtain gmmatr.
    !! This subroutine uses the following subroutines, which are
    !! adaptations of 1979 Linpack subroutines of the same names:
    !!   EQLIBU/dgesl.f
    !! The calling sequences have been preserved to facilitate
    !! substituting equivalents optimized for specific platforms,
    !! for anyone desiring to do so.
    !! This subroutine is called by:
    !!   EQ6/garmat.f
    !! Input:
    !!   gmmatr = the L-U decomposition of the matrix to be inverted
    !!   kdim   = the actual dimention of the matrix
    !!   kmax   = the maximum dimension of the matrix
    !! Workspace:
    !!   delvec = at one point, the right-hand-side vector; at another,
    !!              the solution vector
    !!   ipivot = the pivot vector
    !! Output:
    !!   aimatr = the inverted matrix
    implicit none

    ! Calling sequence variable declarations.
    integer :: kmax

    integer :: ipivot(kmax)
    integer :: kdim

    real(kind=8) :: aimatr(kmax,kmax)
    real(kind=8) :: gmmatr(kmax,kmax)
    real(kind=8) :: delvec(kmax)

    ! Local variable declarations.
    integer :: j
    integer :: k

    ! Invert the matrix aamatr by looping over columns for the inverted
    ! matrix (aimatr).
    do k = 1,kdim
        ! Put the corresponding unit vector in delvec. This serves to
        ! set the right-hand-side vector.
        do j = 1,kdim
            delvec(j) = 0.
        end do

        delvec(k) = 1.

        ! Solve for the unknown vector (delvec). Note that delvec
        ! initially contains the right-hand-side vector. Normally
        ! separate arrays would be used for the two vectors. This
        ! example of bad programming is made here to retain consistency
        ! with the original Linpack routine.
        call dgesl(gmmatr,kmax,kdim,ipivot,delvec)

        ! Load the result in the inverted matrix.
        do j = 1,kdim
            aimatr(j,k) = delvec(j)
        end do
    end do

999 continue
end subroutine ginvrt