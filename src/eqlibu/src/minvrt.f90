subroutine minvrt(aamatr,aimatr,delvec,gmmatr,ier,ipivot,kdim,kmax,noutpt,nttyo,qpr)
    !! This subroutine inverts the matrix aamatr. The inverted matrix
    !! is returned in aimatr. Thus
    !!   [aamatr][aimatr] = I
    !! where I is the unit matrix (1's on the diagonal, 0's elsewhere).
    !! The method is L-U decomposition of the matrix, in combination
    !! with the subsequent solution of a series of matrix equations,
    !! each using the same decomposed matrix but with a different unit
    !! vector as the right hand side. The original aamatr matrix is
    !! preserved. The workspace vector gmmatr holds the L-U
    !! decomposition.
    !! A companion subroutine ginvrt.f has the same function as the
    !! present routine. However, it assumes that aamatr has already
    !! been factored into gmmatr.
    !! This subroutine uses the following subroutines, which are
    !! adaptations of 1979 Linpack subroutines of the same names:
    !!   EQLIBU/dgefa.f
    !!   EQLIBU/dgesl.f
    !! The calling sequences have been preserved to facilitate
    !! substituting equivalents optimized for specific platforms,
    !! for anyone desiring to do so.
    !! This subroutine is called by:
    !!   EQ6/path.f
    !! Input:
    !!   aamatr = the matrix
    !!   kdim   = the actual dimention of the matrix
    !!   kmax   = the maximum dimension of the matrix
    !!   noutpt = the unit number of the output file
    !!   nttyo  = the unit number of the screen file
    !!   qpr    = logical switch to print messages from this subroutine
    !! Workspace:
    !!   delvec = at one point, the right-hand-side vector; at another,
    !!              the solution vector
    !!   gmmatr = work space for a copy of the matrix; holds the
    !!              L-U decomposition
    !!   ipivot = the pivot vector
    !! Output:
    !!   aimatr = the inverted matrix
    !!   ier    = error flag:
    !!              =  1  Encountered a zero matrix
    !!              =  2  Encountered a non-zero, computationally
    !!                      singular matrix
    implicit none

    ! Calling sequence variable declarations.
    integer :: kmax

    integer :: noutpt
    integer :: nttyo

    integer :: ipivot(kmax)
    integer :: ier
    integer :: kdim

    logical :: qpr

    real(kind=8) :: aamatr(kmax,kmax)
    real(kind=8) :: aimatr(kmax,kmax)
    real(kind=8) :: gmmatr(kmax,kmax)
    real(kind=8) :: delvec(kmax)

    ! Local variable declarations.
    integer :: info
    integer :: j
    integer :: k

    ! Copy aamatr to gmmatr.
    do k = 1,kdim
        do j = 1,kdim
            gmmatr(j,k) = aamatr(j,k)
        end do
    end do

    ! OR
    ! do k = 1,kdim
    !   call copyaa(aamatr(1,k),gmmatr(1,k),kdim)
    ! enddo
    ! OR
    ! integer nmax
    ! nmax = kmax*kdim
    ! call copyaa(aamatr,gmmatr,nmax)
    ! Avoid trying to decompose a zero matrix.
    do k = 1,kdim
        do j = 1,kdim
            if (gmmatr(j,k) .ne. 0.) then
                go to 100
            end if
        end do
    end do

    ier = 1

    if (qpr) then
        write (noutpt,1000)
        write (nttyo,1000)
1000 format(/' * Note - (EQLIBU/minvrt) Have a zero matrix to',' invert.')
    end if

    go to 999

    ! Factor the matrix (make the L-U decomposition).
100 continue
    call dgefa(gmmatr,kmax,kdim,ipivot,info)

    ! Test info.
    if (info .ne. 0) then
        ier = 2

        if (qpr) then
            write (noutpt,1010)
            write (nttyo,1010)
1010 format(/' * Note - (EQLIBU/minvrt) Have a computationally',' singular matrix.')
        end if

        go to 999
    end if

    ! Complete the inversion. Loop over columns for the inverted
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
end subroutine minvrt