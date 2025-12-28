subroutine caft1(afrc1,aft1,nrct,nrctmx,rrelr1)
    !! This subroutine computes the total affinity (aft1).
    !! This subroutine is called by:
    !!   EQ6/path.f
    !! Principal input:
    !! Principal output:
    implicit none

    integer :: nrctmx

    integer :: nrct

    real(kind=8) :: aft1
    real(kind=8) :: afrc1(nrctmx)
    real(kind=8) :: rrelr1(nrctmx)

    ! Local variable declarations.
    integer :: nrc

    aft1 = 0.

    do nrc = 1,nrct
        if (rrelr1(nrc) .ne. 0.) then
            if (afrc1(nrc) .lt. 9999999.) then
                aft1 = aft1 + abs(afrc1(nrc)*rrelr1(nrc))
            else
                ! Have a reactant with an "infinite" affinity. Set the
                ! total affinity to "infinity" also.
                aft1 = 9999999.
                go to 100
            end if
        end if
    end do

100 continue
end subroutine caft1