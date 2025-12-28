subroutine dgefa(gmmatr,kmax,kdim,ipivot,info)
    !! This subroutine factors the real*8 array gmmatr, using the method
    !! of L-U decomposition. It is an adaptation of the 1979 Linpack
    !! subroutine of the same name. The size and order of the original
    !! calling sequence has been preserved:
    !! (a,lda,n,ipvt,info) = (gmmatr,kmax,kdim,ipivot,info)
    !! Note: detection of singularity using the test corresponding to
    !! the condition info = -1 was not implemented in the original
    !! Linpack version of the subroutine and may be missing from other
    !! modern equivalents.
    !! This subroutine uses the following pseudo-Linpack BLAS (Basic
    !! Linear Algebra Subsystem) subroutines:
    !!   EQLIBU/idamax.f
    !!   EQLIBU/daxpy.f
    !!   EQLIBU/dscal.f
    !! This subroutine is called by:
    !!   EQLIBU/msolvr.f
    !! Input:
    !!   gmmatr = the matrix (unfactored)
    !!   kmax   = the declared (maximum) dimension of the matrix
    !!   kdim   = the used dimension of the matrix
    !!   ipivot = the pivot vector (used by EQLIBU/dgesl.f)
    !! Output:
    !!   gmmatr = the matrix, factored
    !!   info   = error flag:
    !!              = -1      singularity
    !!              =  0      no singularity
    !!              =  k > 0  singularity
    implicit none

    ! Calling sequence variable declarations.
    integer :: kmax
    integer :: info
    integer :: kdim

    integer :: ipivot(kmax)

    real(kind=8) :: gmmatr(kmax,kmax)

    ! Local variable declarations.
    integer :: j
    integer :: k
    integer :: kp1
    integer :: il
    integer :: idamax

    real(kind=8) :: adiv
    real(kind=8) :: div
    real(kind=8) :: tx

    ! Gaussian elimination with partial pivoting.
    info = 0

    if (kdim .ge. 2) then
        do k = 1,kdim - 1
            kp1 = k + 1

            ! Find il = pivot index.
            il = idamax(kdim - k + 1,gmmatr(k,k),1) + k - 1
            ipivot(k) = il

            ! Zero pivot implies this column already triangularized.
            if (gmmatr(il,k) .ne. 0.) then
                ! Interchange if necessary.
                if (il .ne. k) then
                    tx = gmmatr(il,k)
                    gmmatr(il,k) = gmmatr(k,k)
                    gmmatr(k,k) = tx
                end if

                ! Compute multipliers.
                div = gmmatr(k,k)

                ! Test for potential zero divide.
                adiv=abs(div)

                if (adiv.le.0.) then
                    info = -1
                    go to 999
                end if

                tx = -1./div
                call dscal(kdim - k,tx,gmmatr(k + 1,k),1)

                ! Row elimination with column indexing.
                do j = kp1,kdim
                    tx = gmmatr(il,j)

                    if (il .ne. k) then
                        gmmatr(il,j) = gmmatr(k,j)
                        gmmatr(k,j) = tx
                    end if

                    call daxpy(kdim - k,tx,gmmatr(k + 1,k),1,gmmatr(k + 1,j),1)
                end do

                go to 100
            end if

            info = k
100 continue
        end do
    end if

    ipivot(kdim) = kdim

    if (gmmatr(kdim,kdim) .eq. 0.) then
        info = kdim
    end if

999 continue
end subroutine dgefa