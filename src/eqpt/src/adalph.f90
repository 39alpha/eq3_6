subroutine adalph(alpha,ipbtmx,iz1,iz2,npx2,npx2mx)
    !! This subroutine assigns standard values of the Pitzer
    !! alpha parameters for the npx2-th species pair. The standard
    !! values depend on the charge combination.
    !! This subroutine is called by:
    !!   EQPT/rdpca.f
    !! Principal input:
    !!   ndat0s = unit number of the stripped DATA0 file
    !! Principal output:
    !!   nazt   = the number of specified hard core diameters
    !!   uazp   = array of aqueous species names used to specify
    !!              hard core diamters on the data file
    !!   azero  = array of corresponding hard core diameters
    !!   insgf  = array of corresponding neutral species
    !!              activity coefficient flags
    implicit none

    ! Calling sequence variable declarations.
    integer :: ipbtmx
    integer :: npx2mx

    integer :: iz1
    integer :: iz2
    integer :: npx2

    real(kind=8) :: alpha(ipbtmx,npx2mx)

    ! Local variable declarations.
    integer :: i
    integer :: iaz1
    integer :: iaz2

    ! Assign standard alpha values for the npx2-th species pair,
    ! according the charge combination.
    iaz1 = abs(iz1)
    iaz2 = abs(iz2)

    do i = 1,ipbtmx
        alpha(i,npx2) = 0.
    end do

    if (iaz1.gt.0 .or. iaz2.gt.0) then
        if (iaz1.eq.1 .or. iaz2.eq.1) then
            alpha(1,npx2) = 2.
        else
            alpha(1,npx2) = 1.4
            alpha(2,npx2) = 12.0
        end if
    end if
end subroutine adalph