subroutine dscrax(nf1,nf2,nllnmx,ulinex)
    !! This subroutine descrambles a file of tables whose lines are
    !! interspersed, but which are marked 'a', 'b', 'c', etc., in
    !! column one. The contents of the scrambled file are copied to
    !! the descrambled file as the descrambling takes place. The
    !! descrambled file must already be open. The maximum line length
    !! of this file is nllnmx characters.
    !! This subroutine is called by:
    !!   EQ6/eq6.f
    !! Input:
    !!   nf1    = unit number of the scrambled file
    !!   nf2    = unit number of the descrambled file
    !!   nllnmx = the maximum character length of a line on the
    !!            descrambled files; the maximum character length
    !!            for the scrambled file is nllnmx + 1
    !!   ulinex = variable holding a line
    !! Output:
    !!   None
    implicit none

    ! Calling sequence variable declarations.
    integer :: nf1
    integer :: nf2

    integer :: nllnmx

    character(len=nllnmx) :: ulinex

    ! Local variable declarations.
    integer :: i
    integer :: i1
    integer :: i2
    integer :: j1

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
100 continue
        read (nf1,'(a1,a)',end=110) ux,ulinex

        if (ux .eq. uc) then
            j1 = len_trim(ulinex)
            write (nf2,'(a)') ulinex(1:j1)
        end if

        go to 100

110 continue
    end do
end subroutine dscrax