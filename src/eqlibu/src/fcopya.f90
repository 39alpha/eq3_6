subroutine fcopya(nf1,nf2)
    !! This subroutine appends the contents of the file whose unit number
    !! is nf1 to the file whose unit number is nf2. The line length
    !! is assumed to be 128 characters. The second file  must already
    !! be open and correctly positioned.
    !! This subroutine is called by:
    !!   EQ6/eq6.f
    !! Input:
    !!   nf1    = unit number of the first file
    !!   nf2    = unit number of the second file
    !! Output:
    !!   None
    implicit none

    ! Calling sequence variable declarations.
    integer :: nf1
    integer :: nf2

    ! Local variable declarations.
    integer :: j

    character(len=128) :: uline

    rewind nf1

    do j = 1,10000
        read (nf1,1000,end=999) uline
1000 format(a128)

        write (nf2,1000) uline
    end do

999 continue
end subroutine fcopya