subroutine dscram(nf1,nf2)
    !! This subroutine unscrambles a file of tables whose lines are
    !! interspersed, but which are marked 'a', 'b', 'c', etc., in
    !! column one. The contents of the scrambled file are copied to
    !! the unscrambled file as the unscrambling takes place. The
    !! unscrambled file must already be open. The record length of
    !! this file will is set at 128 characters.
    !! This subroutine is called by:
    !!   EQ6/eq6.f
    !! Input:
    !!   nf1    = unit number of the scrambled file
    !!   nf2   = unit number of the unscrambled file
    !! Output:
    !!   None
    implicit none

    ! Calling sequence variable declarations.
    integer :: nf1
    integer :: nf2

    ! Local variable declarations.
    integer :: i
    integer :: i1
    integer :: i2
    integer :: j

    character(len=128) :: uline
    character(len=1) :: uc
    character(len=1) :: ux

    i1 = ichar('a')
    i2 = ichar('Z')

    if (i2 .lt. i1) then
        i1 = ichar('A')
        i2 = ichar('z')
    end if

    rewind nf2

    do i = i1,i2
        uc = char(i)
        rewind nf1

        do j = 1,10000
            read (nf1,1000,end=100) ux,uline
1000 format(a1,a128)

            if (ux .eq. uc) then
                write (nf2,1005) uline
            end if

1005 format(a128)
        end do

100 continue
    end do
end subroutine dscram