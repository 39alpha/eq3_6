subroutine gcdrw(cdrs,cdrw,narn1,ndrs,ndrsmx,ndrsr,nst,nstmax)
    !! This subroutine computes the cdrw array. Each element of this
    !! array contains the reaction coefficient for liquid water in the
    !! corresponding reaction.
    !! Do not confuse this array with the cdrtw array computed by
    !! EQLIB/gcdrtw.f.
    !! This subroutine is called by:
    !!   EQ3NR/eq3nr.f
    !!   EQ6/path.f
    !!   EQ6/absswb.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: ndrsmx
    integer :: nstmax

    integer :: ndrs(ndrsmx)
    integer :: ndrsr(2,nstmax)

    integer :: narn1
    integer :: nst

    real(kind=8) :: cdrs(ndrsmx)
    real(kind=8) :: cdrw(nstmax)

    ! Local variable declarations.
    integer :: n
    integer :: nn
    integer :: nr1
    integer :: nr2
    integer :: ns

    do ns = 1,nst
        cdrw(ns) = 0.
        nr1 = ndrsr(1,ns)
        nr2 = ndrsr(2,ns)

        do n = nr1,nr2
            nn = ndrs(n)

            if (nn .eq. narn1) then
                cdrw(ns) = cdrs(n)
                go to 100
            end if
        end do

100 continue
    end do
end subroutine gcdrw