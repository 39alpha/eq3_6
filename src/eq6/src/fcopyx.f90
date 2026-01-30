subroutine fcopyx(nf1,nf2,nllnmx,ulinex)
    !! This subroutine appends the contents of the file whose unit
    !! number is nf1 to the file whose unit number is nf2. The line
    !! length is assumed to be nllnmx characters. The second file
    !! must already be open and correctly positioned.
    !! This subroutine is called by:
    !!   EQ6/eq6.f
    !! Input:
    !!   nf1    = unit number of the first file
    !!   nf2    = unit number of the second file
    !!   nllnmx = the maximum character length of a line on the
    !!            scrambled and descrambled files
    !!   ulinex = variable holding a line
    !! Output:
    !!   None
    implicit none

    ! Calling sequence variable declarations.
    integer :: nf1
    integer :: nf2

    integer :: nllnmx

    integer :: j1

    character(len=nllnmx) :: ulinex

    rewind nf1
100 continue
    read (nf1,'(a)',end=110) ulinex
    j1 = len_trim(ulinex)
    write (nf2,'(a)') ulinex(1:j1)
    go to 100

110 continue
end subroutine fcopyx
