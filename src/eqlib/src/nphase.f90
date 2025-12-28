integer function nphase(ncmpr,npt,nptmax,ns)
    !! This subroutine returns the index of the phase which contains
    !! the ns-th species. If a match is not found, a zero value is
    !! returned.
    !! This subroutine is called by:
    !!   EQ6/mincsp.f
    !! Input:
    !!   ncmpr  = array giving the start and the end of the range of the
    !!            species belonging to a given phase
    !!   npt    = the number of phases
    !!   ns     = the index of the desired species
    !! Output:
    !!   nphase = the index of the phase containing the ns-th species
    implicit none

    ! Calling sequence variable declarations.
    integer :: nptmax

    integer :: ncmpr(2,nptmax)
    integer :: npt
    integer :: ns

    ! Local variable declarations.
    integer :: np
    integer :: nr1
    integer :: nr2

    nphase = 0

    do np = 1,npt
        nr1 = ncmpr(1,np)
        nr2 = ncmpr(2,np)

        if (ns .le. nr2) then
            if (ns .ge. nr1) then
                nphase = np
                go to 999
            end if
        end if
    end do

999 continue
end function nphase