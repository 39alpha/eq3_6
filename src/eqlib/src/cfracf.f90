subroutine cfracf(cdrs,csts,efac,jcsort,jflag,jssort,kmax,mosp,narn1,narn2,nbasp,nbaspd,nbt,nbtmax,ndrs,ndrsmx,ndrsr,nern1,nern2,nfac,nst,nstmax,nsts,nstsmx,nstsr,q6mode,weight)
    !! This subroutine determines the dominant species and the associated
    !! exponents required for continued fraction corrections. The data
    !! are returned without any filtering as to their applicability.
    !! particularly when the calling prgram is EQ3NR. For example, data
    !! may be returned for mass balances for H+ or O2(g,aq), when
    !! no such balances are used.
    !! This subroutine is called by:
    !!   EQ3NR/arrset.f
    !!   EQ6/optmzr.f
    !! Principal input:
    !!   mosp   = array containing the number of moles variables for the
    !!              species; the conc array may be substituted for this
    !!              in some calls
    !!   jcsort = array of species indices in sorted order within phase
    !!              ranges
    !!   jssort = array of species indices in sorted order, ignoring
    !!              phase ranges
    !!   q6mode = flag denoting usage for EQ3NR or EQ6:
    !!              .false. = EQ3NR
    !!              .true.  = EQ6NR
    !! Principal output:
    !!   efac   = array of exponents for continued fraction corrections
    !!   nfac   = array of indices of dominant species
    implicit none

    ! Calling sequence variable declarations.
    integer :: kmax
    integer :: nbtmax
    integer :: ndrsmx
    integer :: nstmax
    integer :: nstsmx

    integer :: jcsort(nstmax)
    integer :: jflag(nstmax)
    integer :: jssort(nstmax)
    integer :: nbasp(nbtmax)
    integer :: nbaspd(nbtmax)
    integer :: ndrs(ndrsmx)
    integer :: ndrsr(2,nstmax)
    integer :: nfac(nbtmax)
    integer :: nsts(nstsmx)
    integer :: nstsr(2,nstmax)

    integer :: narn1
    integer :: narn2
    integer :: nbt
    integer :: nern1
    integer :: nern2
    integer :: nst

    logical :: q6mode

    real(kind=8) :: cdrs(ndrsmx)
    real(kind=8) :: csts(nstsmx)
    real(kind=8) :: efac(nbtmax)
    real(kind=8) :: mosp(nstmax)
    real(kind=8) :: weight(nstmax)

    ! Local variable declarations.
    integer :: nb
    integer :: nr1
    integer :: ns
    integer :: nse
    integer :: nsi
    integer :: nss

    real(kind=8) :: cx
    real(kind=8) :: stx
    real(kind=8) :: wsi
    real(kind=8) :: coefdr
    real(kind=8) :: coefst

    ! Initialize the nfac and efac arrays.
    do nb = 1,nbt
        nfac(nb) = 0
        efac(nb) = 1.0
    end do

    ! Loop over all mass balances, used or not.
    do nb = 1,nbt
        nse = nbasp(nb)

        if (jflag(nse) .eq. 0) then
            do ns = 1,nst
                weight(ns) = 0.
            end do

            do nss = narn1,narn2
                ns = jcsort(nss)
                weight(ns) = coefst(csts,nsts,nstsmx,nstsr,nb,ns,nstmax)
            end do

            if (q6mode) then
                do nss = nern1,nern2
                    ns = jcsort(nss)
                    weight(ns) = coefst(csts,nsts,nstsmx,nstsr,nb,ns,nstmax)
                end do
            end if

            ! Find the species (nfac) that makes the largest contribution
            ! to the mass balance. Get the exponent (efac) required for the
            ! continued fraction correction.
            call fdomsp(jssort,mosp,nsi,nst,nstmax,weight,wsi)

            ! Get the exponent (efac) required for the continued fraction
            ! correction.
            nfac(nb) = nsi
            ns = nbaspd(nb)

            if (nse .eq. ns) then
                efac(nb) = 1./wsi
            else
                cx = coefdr(cdrs,ndrs,ndrsmx,ndrsr,nse,ns,nstmax)
                nr1 = ndrsr(1,ns)
                stx = -cx/cdrs(nr1)
                efac(nb) = 1./(stx*wsi)
            end if
        end if
    end do
end subroutine cfracf