real(kind=8) function coefel(cess,ness,nessmx,nessr,nc,ns,nstmax)
    !! This subroutine finds the coefficient of the nc-th element
    !! in the composition of the ns-th species.
    !! This subroutine is called by:
    !!   None
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: nessmx
    integer :: nstmax

    integer :: ness(nessmx)
    integer :: nessr(2,nstmax)
    integer :: nc
    integer :: ns

    real(kind=8) :: cess(nessmx)

    ! Local variable declarations.
    integer :: n
    integer :: n1
    integer :: n2

    coefel = 0.
    n1 = nessr(1,ns)
    n2 = nessr(2,ns)

    do n = n1,n2
        if (nc .eq. ness(n)) then
            coefel = cess(n)
            go to 999
        end if
    end do

999 continue
end function coefel