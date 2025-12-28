subroutine gakmat(akmat0,dxsm00,nord,nrd1mx)
    !! This subroutine computes the akmat0 matrix, which allows
    !! calculation of estimates of derivatives from corresponding finite
    !! differences. Note: this subroutine presumes that the diagonal
    !! elements have been previously set to factorials and that the
    !! lower triangle elements have previously been set to zero. These
    !! parts of the matrix are invariant.
    !! This routine can also be used to calculate the akmat1 matrix from
    !! the dlxsm1 array, provided that in the calling sequence akmat1
    !! is substituted for akmat0 and dlxsm1 for dxsm00. The akmat0 matrix
    !! and the dxsm00 array are used with finite differences employed
    !! to generate predictor functions, while akmat1 and dlxsm1 are used
    !! are used with finite differences employed to generate corrector
    !! functions.
    !! This subroutine is called by:
    !!   EQ6/path.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: nrd1mx

    integer :: nord

    real(kind=8) :: akmat0(nrd1mx,nrd1mx)
    real(kind=8) :: dxsm00(nrd1mx)

    ! Local variable declarations.
    integer :: i
    integer :: im1
    integer :: ip1
    integer :: j
    integer :: jm1

    if (nord .ge. 2) then
        i = 1
        ip1 = i + 1

        do j = ip1,nord
            jm1 = j - 1
            akmat0(i,j) = akmat0(i,jm1)*dxsm00(jm1)
        end do

        do i = 2,nord - 1
            ip1 = i + 1
            im1 = i - 1

            do j = ip1,nord
                jm1 = j - 1
                akmat0(i,j) = akmat0(i,jm1)*dxsm00(jm1) + akmat0(im1,jm1)*i
            end do
        end do
    end if
end subroutine gakmat