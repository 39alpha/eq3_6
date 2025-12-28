subroutine gdlgxw(cdrs,cjbasp,cnufac,conc,dlogxw,eps100,ixbasp,jcsort,jflag,narn1,narn2,nbasp,nbt,nbtmax,nbw,ndrs,ndrsmx,ndrsr,nern1,nern2,noutpt,nstmax,nttyo,omega,xbar,xbarw)
    !! This subroutine computes the array dlogxw, which contains the
    !! partial derivatives:
    !!   d log x(w)/d log m(s')
    !! where s' denotes a basis species other than water. Allowed
    !! basis species may be aqueous or generic ion exchanger species.
    !! This subroutine is called by:
    !!   EQLIB/matrix.f
    !!   EQ3NR/arrsim.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: nbtmax
    integer :: ndrsmx
    integer :: nstmax

    integer :: noutpt
    integer :: nttyo

    integer :: ixbasp(nbtmax)
    integer :: jcsort(nstmax)
    integer :: jflag(nstmax)
    integer :: nbasp(nbtmax)
    integer :: ndrs(ndrsmx)
    integer :: ndrsr(2,nstmax)
    integer :: narn1
    integer :: narn2
    integer :: nbt
    integer :: nbw
    integer :: nern1
    integer :: nern2

    real(kind=8) :: cdrs(ndrsmx)
    real(kind=8) :: cjbasp(nbtmax)
    real(kind=8) :: cnufac(nstmax)
    real(kind=8) :: conc(nstmax)
    real(kind=8) :: dlogxw(nbtmax)
    real(kind=8) :: xbar(nstmax)
    real(kind=8) :: eps100
    real(kind=8) :: omega
    real(kind=8) :: xbarw

    ! Local variable declarations.
    integer :: nb
    integer :: nsc
    integer :: nse
    integer :: nss

    real(kind=8) :: cw
    real(kind=8) :: cxec
    real(kind=8) :: dx
    real(kind=8) :: sumw
    real(kind=8) :: xx

    real(kind=8) :: coefdr

    if (nbw .eq. 0) then
        write (noutpt,1000)
        write (nttyo,1000)
1000 format(/' * Error - (EQLIB/dxlgxw) Programming error trap:',/7x,'Water is not in the basis set, so the dlogxw array',/7x,'should not be calculated.')

        stop
    end if

    ! Make sure that the "conc" of water is zero. This will be restored
    ! at the end of this subroutine. This practice allows avoiding "if"
    ! testing in the following loops. Note that there is no actual
    ! element of the dlogxw array for water (basis index nbw). However,
    ! there is a sum used in the denominator of the other array elements
    ! that is based on water, and this resembles the sums accumulated
    ! for the other elements. Therefore, dlogxw(nbw) will be used here
    ! to accumulate that sum. It will be reset to zero after the sum
    ! is used.
    cw = conc(narn1)
    conc(narn1) = 0.

    ! Compute the sums needed to calculate the dlogxw array.
    ! Store them in this array.
    do nb = 1,nbt
        nse = nbasp(nb)
        dlogxw(nb) = conc(nse)
    end do

    do nss = narn1,narn2
        nsc = jcsort(nss)

        if (jflag(nsc) .eq. 30) then
            do nb = 1,nbt
                nse = nbasp(nb)

                ! Calling sequence substitutions:
                !   nsc for ns
                cxec = coefdr(cdrs,ndrs,ndrsmx,ndrsr,nse,nsc,nstmax)

                ! Here cnufac(nsc) = -conc(nsc)/cdrs(ndrsr(1,nsc)), or
                ! the negative of the molality of a dependent species
                ! divided by its (intrinsically negative) reaction
                ! coefficient in its associated reaction.
                dlogxw(nb) = dlogxw(nb) + cxec*cnufac(nsc)
            end do
        end if
    end do

    ! Compute the dlogxw array from the sums. The common denominator
    ! (dx) for the elements corresponding to basis species other than
    ! water is computed from the sum for water.
    ! Do not protect here against a value of xbarw that is ridiculously
    ! small, say on the order of 1.e-30. Such a value seems to be
    ! necessary to obtaining good results in pre-Newton-Raphson
    ! optimization. Only protect against a zero divide below, where
    ! dx (derived from xbarw) appears in a denominator.
    xx = xbarw/omega
    sumw = dlogxw(nbw)
    dlogxw(nbw) = 0.
    dx = 1. + xx*sumw

    ! Protect against the case in which dx is zero. This is designed
    ! to prevent a fatal floating point exception in pre-Newton-Raphson
    ! optimization.
    if (dx.ge.0. .and.  dx.lt.eps100) then
        dx = eps100
    end if

    if (dx.lt.0. .and.  dx.gt.-eps100) then
        dx = -eps100
    end if

    do nb = 1,nbt
        if (nb .ne. nbw) then
            dlogxw(nb) = -xx*dlogxw(nb)/dx

            if (ixbasp(nb) .ge. 1) then
                ! Here nb refers to a basis species whose thermodynamic
                ! activity is defined by its mole fraction.
                nse = nbasp(nb)
                dlogxw(nb) = dlogxw(nb)*cjbasp(nb)*(1.0 - xbar(nse))
            end if
        end if
    end do

    ! Restore the "conc" of water.
    conc(narn1) = cw
end subroutine gdlgxw