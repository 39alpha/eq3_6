subroutine jflaux(jflag,nbaspd,nbtd,nbtmax,ndrsd,ndrsmx,ndrsrd,nstmax)
    !! This subroutine sets jflag to -1 for auxiliary basis species that
    !! can not appear in the model.
    !! This subroutine is called by:
    !!   EQ3NR/eq3nr.f
    !!   EQ6/eq6.f
    !! Principal input:
    !!   jflag  = array of species control flags
    !!   nbaspd = array of indices of species in the data file basis
    !!              set
    !!   nbtd   = number of species in the data file basis set
    !! Principal output:
    !!   jflag  = array of species control flags
    implicit none

    ! Calling sequence variable declarations.
    integer :: nbtmax
    integer :: ndrsmx
    integer :: nstmax

    integer :: jflag(nstmax)
    integer :: nbaspd(nbtmax)
    integer :: ndrsd(ndrsmx)
    integer :: ndrsrd(2,nstmax)
    integer :: nbtd

    ! Local variable declarations.
    integer :: n
    integer :: nb
    integer :: ncount
    integer :: nr1
    integer :: nr2
    integer :: ns
    integer :: nse

    ! Note that looping is required to accommodate multi-level
    ! dependencies in the auxiliary basis set. If jflag is set
    ! to -1 for any auxiliary basis species, another pass is required
    ! to insure that jflag is also set to -1 for any higher-level
    ! auxiliary basis species.
100 continue
    ncount = 0

    do nb = 1,nbtd
        ns = nbaspd(nb)
        nr1 = ndrsrd(1,ns)
        nr2 = ndrsrd(2,ns)

        if (jflag(ns) .eq. 30) then
            do n = nr1 + 1,nr2
                nse = ndrsd(n)

                if (jflag(nse) .eq. -1) then
                    jflag(ns) = -1
                    ncount = ncount + 1
                    go to 110
                end if
            end do
        end if

110 continue
    end do

    if (ncount .gt. 0) then
        go to 100
    end if
end subroutine jflaux