subroutine dgesl(gmmatr,kmax,kdim,ipivot,delvec)
    !! This subroutine solves the real*8 system (gmmatr)*(x) = (b),
    !! where the matrix gmmatr has previously been factored by
    !! EQLIBU/dgefa.f. It is expected that singularity of the matrix
    !! has been trapped by that subroutine. Initially, the array delvec
    !! contains the right-hand-side vector (b); this is replaced by
    !! the solution vector (x). This subroutine is an adaptation of the
    !! 1979 Linpack subroutine of the same name. The size and order in
    !! the original calling sequence has been preserved:
    !! (a,lda,n,ipvt,b) = (gmmatr,kmax,kdim,ipivot,delvec)
    !! This subroutine uses the following pseudo-Linpack BLAS (Basic
    !! Linear Algebra Subsystem) subroutines:
    !!   EQLIBU/daxpy.f
    !! This subroutine is called by:
    !!   EQLIBU/msolvr.f
    !! Input:
    !!   gmmatr = the matrix
    !!   kmax   = the declared (maximum) dimension of the matrix
    !!   kdim   = the used dimension of the matrix
    !!   ipivot = the pivot vector (from EQLIBU/dgefa.f)
    !!   delvec = initially, the right hand side vector (b)
    !! Output:
    !!   delvec = the solution vector (x).
    implicit none

    ! Calling sequence variable declarations.
    integer :: kmax
    integer :: kdim

    integer :: ipivot(kmax)

    real(kind=8) :: gmmatr(kmax,kmax)
    real(kind=8) :: delvec(kmax)

    ! Local variable declarations.
    integer :: il
    integer :: k
    integer :: kb

    real(kind=8) :: tx

    ! First solve (L)*(y) = (b).
    if (kdim .ge. 2) then
        do k = 1, kdim - 1
            il = ipivot(k)
            tx = delvec(il)

            if (il .ne. k) then
                delvec(il) = delvec(k)
                delvec(k) = tx
            end if

            call daxpy(kdim - k,tx,gmmatr(k + 1,k),1,delvec(k + 1),1)
        end do
    end if

    ! Now solve  (U)*(x) = (y).
    do kb = 1, kdim
        k = kdim + 1 - kb
        delvec(k) = delvec(k)/gmmatr(k,k)
        tx = -delvec(k)
        call daxpy(k - 1,tx,gmmatr(1,k),1,delvec(1),1)
    end do
end subroutine dgesl