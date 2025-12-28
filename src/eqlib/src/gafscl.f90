subroutine gafscl(cdrsd,cscale,ndrsmx,ndrsrd,nst,nstmax)
    !! This subroutine calculates the cscale array of affinity scaling
    !! factors. Affinity scaling is used to help make decisions on
    !! which phases are the best choices to precipitate when there
    !! are multiple supersaturations.
    !! This subroutine is called by:
    !!   EQ6/eq6.f
    !! Principal input:
    !!   cdrsd  = array of reaction coefficients ('d' set)
    !! Principal output:
    !!   cscale = array of affinity scaling factors
    implicit none

    ! Calling sequence variable declarations.
    integer :: ndrsmx
    integer :: nstmax

    integer :: ndrsrd(2,nstmax)
    integer :: nst

    real(kind=8) :: cdrsd(ndrsmx)
    real(kind=8) :: cscale(nstmax)

    ! Local variable declarations.
    integer :: n
    integer :: ns
    integer :: nr1
    integer :: nr2

    real(kind=8) :: cx

    do ns = 1,nst
        cx = 0.
        nr1 = ndrsrd(1,ns)
        nr2 = ndrsrd(2,ns)

        do n = nr1 + 1,nr2
            cx = cx + abs(cdrsd(n))
        end do

        if (cx .le. 0.) then
            cx = 1.0
        end if

        cscale(ns) = cx
    end do
end subroutine gafscl