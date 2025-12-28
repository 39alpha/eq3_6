subroutine wrazp(azero,insgf,nazt,naztmx,ndata1,ndat1f,noutpt,nttyo,uazp)
    !! This subroutine writes on the DATA1 and DATA1F files the "bdot"
    !! data read from the DATA0 file by EQPT/rdazp.f.
    !!   EQPT/eqpt.f
    !! Principal input:
    !!   ndata1 = unit number of the DATA1 file
    !!   ndat1f = unit number of the DATA1F file
    !!   nazt   = the number of elements in the uazp array
    !!   uazp   = array containing lines of data
    !!   azero  = array of corresponding hard core diameters
    !!   insgf  = array of corresponding neutral species
    !!              activity coefficient flags
    !! Principal output:
    !!   None
    implicit none

    ! Calling sequence variable declarations.
    integer :: naztmx

    integer :: ndata1
    integer :: ndat1f
    integer :: noutpt
    integer :: nttyo

    integer :: nazt

    integer :: insgf(naztmx)

    character(len=24) :: uazp(naztmx)

    real(kind=8) :: azero(naztmx)

    ! Local variable declarations.
    integer :: ineu
    integer :: naz

    character(len=72) :: uterm
    character(len=72) :: utermc
    character(len=24) :: uendit

    real(kind=8) :: zero

    data uendit / 'endit.' /

    uterm(1:48) = '+-----------------------------------------------'
    uterm(49:72) = '------------------------'
    utermc = uterm
    utermc(1:1) = '*'

    ! Write the azero ('bdot') data on the DATA1 and DATA1F files.
    do naz = 1,nazt
        write (ndata1) uazp(naz),azero(naz),insgf(naz)
        write (ndat1f,1010) uazp(naz),azero(naz),insgf(naz)
1010 format(a24,2x,f7.1,2x,i2)
    end do

    ! Write the block terminator.
    zero = 0.
    ineu = 0
    write (ndata1) uendit,zero,ineu
    write (ndat1f,1010) uendit,zero,ineu
    write (ndat1f,1020) utermc(1:72)
1020 format(a)
end subroutine wrazp