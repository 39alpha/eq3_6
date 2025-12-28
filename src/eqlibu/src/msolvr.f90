subroutine msolvr(aamatr,delvec,gmmatr,ier,ipivot,kdim,kmax,noutpt,nttyo,qpr,rhsvec)
    !! This subroutine solves the matrix equation:
    !!   (aamatr)(delvec) = (rhsvec)
    !! The method employed is L-U decomposition. The original aamatr
    !! matrix is preserved. The workspace vector gmmatr holds the
    !! L-U decomposition.
    !! This subroutine uses
    !! the following subroutines, which are adaptations of 1979 Linpack
    !! subroutines of the same names:
    !!   EQLIBU/dgefa.f
    !!   EQLIBU/dgesl.f
    !! The calling sequences have been preserved to facilitate
    !! substituting equivalents optimized for specific platforms,
    !! for anyone desiring to do so.
    !! This subroutine is called by:
    !!   EQLIBU/polfit.f
    !!   EQLIBU/minvrt.f
    !!   EQLIB/nrstep.f
    !!   EQ3NR/arrsim.f
    !!   EQ6/jgibbs.f
    !! Input:
    !!   aamatr = the matrix
    !!   kdim   = the actual dimention of the matrix
    !!   kmax   = the maximum dimension of the matrix
    !!   noutpt = the unit number of the output file
    !!   nttyo  = the unit number of the screen file
    !!   rhsvec = the right hand side vector
    !!   qpr    = logical switch to print messages from this subroutine
    !! Workspace:
    !!   gmmatr = work space for a copy of the matrix; holds the
    !!              L-U decomposition
    !!   ipivot = the pivot vector
    !! Output:
    !!   delvec = the solution vector
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
    real(kind=8) :: gmmatr(kmax,kmax)
    real(kind=8) :: rhsvec(kmax)
    real(kind=8) :: delvec(kmax)

    ! Local variable declarations.
    integer :: info
    integer :: j
    integer :: k

    ! Initialize error flag.
    ier = 0

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
1000 format(/' * Note - (EQLIBU/msolvr) Have a zero matrix to',' factor.')
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
1010 format(/' * Note - (EQLIBU/msolvr) Have a computationally',' singular matrix.')
        end if

        go to 999
    end if

    ! Copy rhsvec to delvec.
    do j = 1,kdim
        delvec(j) = rhsvec(j)
    end do

    ! OR
    ! call copyaa(rhsvec,delvec,kdim)
    ! Solve for the unknown vector (delvec). Note that delvec initially
    ! contains a copy of the right-hand-side vector.
    call dgesl(gmmatr,kmax,kdim,ipivot,delvec)

999 continue
end subroutine msolvr