subroutine intinx(ier,ixrn1,ixrn2,ncmpr,ncmpri,noutpt,npnxp,nptmax,nstmax,nttyo,nxicmx,nxti,nxtimx,umemi,uphase,usoli,uspec,xbar,xbari,xbarlg)
    !! This subroutine interpets the solid solution compositions
    !! specified on the input file.
    !! This subroutine is called by:
    !!   EQ3NR/eq3nr.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: nptmax
    integer :: nstmax
    integer :: nxicmx
    integer :: nxtimx

    integer :: noutpt
    integer :: nttyo

    integer :: ncmpr(2,nptmax)
    integer :: ncmpri(2,nxtimx)
    integer :: npnxp(nxtimx)
    integer :: ier
    integer :: ixrn1
    integer :: ixrn2
    integer :: nxti

    character(len=48) :: uspec(nstmax)
    character(len=24) :: umemi(nxicmx)
    character(len=24) :: uphase(nptmax)
    character(len=24) :: usoli(nxtimx)

    real(kind=8) :: xbar(nstmax)
    real(kind=8) :: xbari(nxicmx)
    real(kind=8) :: xbarlg(nstmax)

    ! Local variable declarations.
    integer :: j2
    integer :: j3
    integer :: np
    integer :: nr1
    integer :: nr1i
    integer :: nr2
    integer :: nr2i
    integer :: ns
    integer :: nxi
    integer :: nxic

    integer :: ilnobl

    character(len=24) :: uxp
    character(len=24) :: uxs

    real(kind=8) :: xx

    real(kind=8) :: tlg

    ier = 0

    if (nxti .le. 0) then
        go to 999
    end if

    do nxi = 1,nxti
        npnxp(nxi) = 0
        uxp = usoli(nxi)

        do np = ixrn1,ixrn2
            if (uxp(1:24) .eq. uphase(np)(1:24)) then
                go to 100
            end if
        end do

        ier = ier + 1
        j2 = ilnobl(uxp)
        write (noutpt,1000) uxp(1:j2)
        write (nttyo,1000) uxp(1:j2)
1000 format(/' * Error - (EQ3NR/intinx) The solid solution ',a,/7x,'is on the input file, but not on the data file.')

        go to 130

100 continue
        npnxp(nxi) = np
        nr1i = ncmpri(1,nxi)
        nr2i = ncmpri(2,nxi)
        nr1 = ncmpr(1,np)
        nr2 = ncmpr(2,np)

        do nxic = nr1i,nr2i
            uxs = umemi(nxic)

            do ns = nr1,nr2
                if (uxs(1:24) .eq. uspec(ns)(1:24)) then
                    go to 110
                end if
            end do

            ier = ier + 1
            j2 = ilnobl(uxs)
            j3 = ilnobl(uxp)
            write (noutpt,1010) uxs(1:j2),uxp(1:j3)
            write (nttyo,1010) uxs(1:j2),uxp(1:j3)
1010 format(/' * Error - (EQ3NR/intinx) The solid solution',/7x,'species ',a,' (',a,') is on the'    /7x,'input file, but not on the data file.')

            go to 120

110 continue

            xx = xbari(nxic)
            xbar(ns) = xx
            xbarlg(ns) = tlg(xx)
120 continue
        end do
    end do

130 continue

999 continue
end subroutine intinx