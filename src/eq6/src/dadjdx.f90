subroutine dadjdx(delxi,dlxmin,iodb,nodbmx,noutpt,qadjdx)
    !! This subroutine determines whether or not delxi can be reduced
    !! to satisfy some criterion, such as the pH not exceeding the
    !! requested maximum value.
    !! This subroutine is called by:
    !!   EQ6/path.f
    !! Principal input:
    !! Principal output:
    implicit none

    integer :: nodbmx

    integer :: noutpt

    integer :: iodb(nodbmx)

    logical :: qadjdx

    real(kind=8) :: delxi
    real(kind=8) :: dlxmin

    ! Local variable declarations.
    ! None
    qadjdx = .false.

    if (delxi .le. dlxmin) then
        ! The step size is already at the minimum value.
        ! Do not cut it.
        if (iodb(1) .gt. 0) then
            write (noutpt,1100)
        end if

1100 format(3x,'The step size will not be cut because it is',' already at the',/5x,'minimum value.',/)
    else
        ! Set up to go back and cut the step size to satisfy the
        ! accuracy criterion for calculating the point at which the
        ! event of concern occurs.
        qadjdx = .true.
    end if
end subroutine dadjdx