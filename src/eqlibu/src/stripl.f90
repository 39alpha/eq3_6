subroutine stripl(nin,nout)
    !! This subroutine copies the file whose unit number is "nin" to
    !! that whose unit number is "nout". Lines beginning with an asterix
    !! are not copied. Lines exceeding a length of 80 characters are
    !! truncated.
    !! This subroutine is called by:
    !!   EQPT/ofiles.f
    !!   EQ3NR/eq3nr.f
    !!   EQ6/eq6.f
    !! Input:
    !!   nin    = unit number of the file to be stripped
    !!   nout   = unit number of the stripped file
    !! Output:
    !!   None
    implicit none

    ! Calling sequence variable declarations.
    integer :: nin
    integer :: nout

    ! Local variable declarations.
    integer :: j2

    integer :: ilnobl

    character(len=80) :: uibuf
    character(len=1) :: ustar
    character(len=1) :: ux

    data ustar  /'*'/

10 continue
    read (nin,1000,end=999) uibuf
1000 format(a80)

    ! Skip lines with an asterisk in column 1.
    ux = uibuf(1:1)

    if (ux .eq. ustar) then
        go to 10
    end if

    j2 = ilnobl(uibuf)

    if (j2 .le. 0) then
        j2 = 1
    end if

    write (nout,1010) uibuf(1:j2)
1010 format(a)

    go to 10

999 continue
end subroutine stripl