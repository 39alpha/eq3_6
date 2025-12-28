subroutine g1dat(ier,noutpt,nttyo,udastr,var)
    !! This subroutine reads a number (var) from a string (udastr)
    !! This subroutine is called by:
    !!   EQPT/gnenb.f
    !!   EQPT/rdpca.f
    !!   EQPT/rdpth.f
    !!   EQPT/rdpni.f
    !!   EQPT/rdpn2.f
    !!   EQPT/rdpnn.f
    !! Principal input:
    !!   udastr = a string containing one numerical data field
    !! Principal output:
    !!   ier    = error flag
    !!   var    = the number contained in that field
    implicit none

    ! Calling sequence variable declarations.
    integer :: noutpt
    integer :: nttyo

    integer :: ier

    character(len=80) :: udastr

    real(kind=8) :: var

    ! Local variable declarations.
    integer :: jj
    integer :: jlen

    integer :: ilnobl

    var = 0.
    ier = 0

    call lejust(udastr)
    jj = index(udastr,' ')

    if (jj .eq. 0) then
        jj = 81
    end if

    if (jj .lt. 80) then
        udastr(jj:80) = ' '
    end if

    jlen = ilnobl(udastr)

    if (jlen .gt. 25) then
        write (noutpt,1000) udastr(1:jlen)
        write (nttyo,1000) udastr(1:jlen)
1000 format(/' * Error - (EQPT/g1dat) Have found a data',' field in the string:',/7x,'"',a,'"',/7x,'that appears',' to exceed the allowed 25 characters.')

        ier = 1
    else
        read (udastr,1010,err=995) var
1010 format(e25.18)
    end if

    go to 999

995 continue
    ier = 1

999 continue
end subroutine g1dat