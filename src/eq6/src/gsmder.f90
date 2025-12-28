subroutine gsmder(delxia,drer0,drer0s,drir0,drir0s,jreac,nord,nrct,nrctmx,nrd1mx)
    !! This subroutine computes smoothed or averaged derivatives of the
    !! elements of the r vector. The actual derivatives (of various
    !! orders) are averaged over the interval (-delxia,+delxia). The
    !! base point (point 0) is at the center of this interval.
    !! Subroutine gsmdez.f performs the same function for elements of
    !! the z vector.
    !! This subroutine is called by:
    !!   EQ6/path.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: nrctmx
    integer :: nrd1mx

    integer :: jreac(nrctmx)

    integer :: nord
    integer :: nrct

    real(kind=8) :: drer0(nrd1mx,nrctmx)
    real(kind=8) :: drer0s(nrd1mx,nrctmx)
    real(kind=8) :: drir0(nrd1mx)
    real(kind=8) :: drir0s(nrd1mx)

    real(kind=8) :: delxia

    ! Local variable declarations.
    integer :: k
    integer :: n
    integer :: nmax
    integer :: nrc

    real(kind=8) :: drx
    real(kind=8) :: dx

    ! Zero the averaged derivative arrays.
    do n = 1,nrd1mx
        drir0s(n) = 0.
    end do

    nmax = nrd1mx*nrctmx
    call initaz(drer0s,nmax)

    ! The averaged derivatives are calculated over the interval
    ! (-delxia,+delxia), about the base point (point 0).
    dx = 2.*delxia

    ! Calculate average derivatives for the inverse rate.
    do n = 1,nord
        drx = drir0(nord)

        do k = nord - 1,n,-1
            drx = drir0(k) + (drx*dx/(nord - n + 1))
        end do

        drir0s(n) = drx
    end do

    ! Calculate average derivatives for relative rates.
    do nrc = 1,nrct
        if (jreac(nrc).eq.0 .or. jreac(nrc).eq.-1)  then
            do n = 1,nord
                drx = drer0(nord,nrc)

                do k = nord - 1,n,-1
                    drx = drer0(k,nrc) + (drx*dx/(nord - n + 1))
                end do

                drer0s(n,nrc) = drx
            end do
        end if
    end do

999 continue
end subroutine gsmder