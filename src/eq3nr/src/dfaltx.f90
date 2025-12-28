subroutine dfaltx(itermx,rho,scamas,tdspkg,tdspl,tolbt,toldl,tolspf)
    !! This subroutine sets the defaults for various run parameters. In
    !! some cases, it forces the parameters to take on certain values,
    !! or to fall in certain ranges.
    !! This subroutine is called by:
    !!   EQ3NR/eq3nr.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: itermx

    real(kind=8) :: rho
    real(kind=8) :: scamas
    real(kind=8) :: tdspkg
    real(kind=8) :: tdspl
    real(kind=8) :: tolbt
    real(kind=8) :: toldl
    real(kind=8) :: tolspf

    ! Local variable declarations.
    !   None
    if (itermx .le. 0) then
        itermx = 200
    end if

    if (rho .le. 0.) then
        rho = 1.0
    end if

    if (tdspkg .lt. 0.) then
        tdspkg = 0.
    end if

    if (tdspl .lt. 0.) then
        tdspl = 0.
    end if

    if (scamas .le. 0.) then
        scamas = 1.0
    end if

    if (tolbt .le. 0.) then
        tolbt = 1.e-6
    end if

    if (tolbt .lt. 1.e-10) then
        tolbt = 1.e-10
    end if

    if (tolbt .gt. 1.e-2) then
        tolbt = 1.e-2
    end if

    if (toldl .le. 0.) then
        toldl = 1.e-6
    end if

    if (toldl .lt. 1.e-10) then
        toldl = 1.e-10
    end if

    if (toldl .gt. 1.e-2) then
        toldl = 1.e-2
    end if

    if (tolspf .le. 0.) then
        tolspf = 0.0005
    end if

    if (tolspf .lt. 0.01) then
        tolspf = 0.01
    end if

    if (tolspf .gt. 0.00005) then
        tolspf = 0.00005
    end if
end subroutine dfaltx