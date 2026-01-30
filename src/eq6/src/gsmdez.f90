subroutine gsmdez(delxia,dzvc0,dzvc0s,kdim,kmax,nord,nrd1mx)
    !! This subroutine computes smoothed or averaged derivatives of the
    !! elements of the z vector. The actual derivatives (of various
    !! orders) are averaged over the interval (-delxia,+delxia). The
    !! base point (point 0) is at the center of this interval.
    !! Subroutine gsmder.f performs the same function for elements of
    !! the r vector.
    !! This subroutine is called by:
    !!   EQ6/path.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: kmax
    integer :: nrd1mx

    integer :: nord
    integer :: kdim

    real(kind=8) :: dzvc0(nrd1mx,kmax)
    real(kind=8) :: dzvc0s(nrd1mx,kmax)

    real(kind=8) :: delxia

    ! Local variable declarations.
    integer :: k
    integer :: kcol
    integer :: n
    integer :: nmax

    real(kind=8) :: dx
    real(kind=8) :: dzx

    ! Zero the averaged derivative arrays.
    nmax = nrd1mx*kmax
    call initaz(dzvc0s,nmax)

    ! The averaged derivatives are calculated over the interval
    ! (-delxia,+delxia), about the base point (point 0).
    dx = 2.*delxia

    ! Calculate average derivatives for elements of the z vector.
    do kcol = 1,kdim
        do n = 1,nord
            dzx = dzvc0(nord,kcol)

            do k = nord - 1,n,-1
                dzx = dzvc0(k,kcol) + (dzx*dx/(nord - n + 1))
            end do

            dzvc0s(n,kcol) = dzx
        end do
    end do
end subroutine gsmdez
