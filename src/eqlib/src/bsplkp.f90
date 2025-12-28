subroutine bsplkp(axlks,narxmx,nbasp,nbt,nbtmax,ndrs,ndrsmx,ndrsr,nstmax,ntprmx)
    !! This subroutine looks at each active auxiliary basis species.
    !! It changes the corresponding log K polynomial coefficients so
    !! that log K is fixed at a value of -9999999. if any other species
    !! in the corresponding dissociation reaction is not in the model.
    !! This subroutine is called by:
    !!   EQ3NR/eq3nr.f
    !!   EQ6/eq6.f
    !! Principal input:
    !!   nbasp  = array of indices of species in the active basis set
    !! Principal output:
    !!   axlks  = array of coefficients for computing equilbrium
    !!              constants as a function of temperature
    implicit none

    ! Calling sequence variable declarations.
    integer :: narxmx
    integer :: nbtmax
    integer :: ndrsmx
    integer :: nstmax
    integer :: ntprmx

    integer :: nbasp(nbtmax)
    integer :: ndrs(ndrsmx)
    integer :: ndrsr(2,nstmax)
    integer :: nbt

    real(kind=8) :: axlks(narxmx,ntprmx,nstmax)

    ! Local variable declarations.
    integer :: i
    integer :: j
    integer :: n
    integer :: nb
    integer :: ns
    integer :: nr1
    integer :: nr2
    integer :: nt1
    integer :: nse

    do nb = 1,nbt
        ns = nbasp(nb)
        nr1 = ndrsr(1,ns)
        nr2 = ndrsr(2,ns)
        nt1 = nr2 - nr1 + 1

        if (nt1 .ge. 2) then
            do n = nr1,nr2
                nse = ndrs(n)

                if (nse .eq. 0) then
                    do j = 1,ntprmx
                        do i = 2,narxmx
                            axlks(i,j,ns) = 0.
                        end do

                        axlks(1,j,ns) = -9999999.
                    end do

                    go to 130
                end if
            end do
        end if

130 continue
    end do
end subroutine bsplkp