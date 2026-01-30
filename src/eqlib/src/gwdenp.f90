subroutine gwdenp(adwipp,bdwipp,jcsort,mlmrra,mosp,mrmlra,mwtsp,narn1,narn2,nstmax,qdwipp,rhoc,rhowc,tdsgks,tdsglw,tdspkc,tdsplc,tempc,vosol,wfh2o,wftds,wkgwi,woh2o,wosol,wotds)
    !! This subroutine gets the weights of solute, total dissolved
    !! solutes, and aqueous solution, and obtains the density of
    !! the aqueous solution by evaluating the WIPP brine density
    !! model. It also computes some related parameters, such
    !! as the solution volume and the molarity/molality ratio.
    !! This subroutine is called by:
    !!   EQ3NR/eq3nr.f
    !!   EQ6/path.f
    !! Principal input:
    !!   adwipp = 'a' constant in the WIPP brine density model
    !!   bdwipp = 'b' constant in the WIPP brine density model
    !!   jcsort = array of species indices in order of increasing
    !!              concentration, organized into blocks, one per phase
    !!   mosp   = numbers of moles of species
    !!   mwtsp  = molecular weights of species
    !!   mwtsp  = molecular weights of species
    !! Principal output:
    !!   mlmrra = molality/molarity ratio
    !!   mrmlra = molarity/molality ratio
    !!   qdwipp = .true. if results from the WIPP brine density model
    !!              are available
    !!   rhoc   = the solution density (g/mL), calculated from the WIPP
    !!              brine density model
    !!   rhowc  = the solution density (g/L), calculated from the WIPP
    !!              brine density model
    !!   tdsglw = TDS (total dissolved solutes) (g/L), calculated
    !!              from the WIPP brine density model
    !!   tdsgks = TDS (total dissolved solutes) (g/kg.sol)
    !!   tdspkc = TDS (total dissolved solutes) (mg/kg.sol)
    !!   tdsplc = TDS (total dissolved solutes) (mg/L), calculated
    !!              from the WIPP brine density model
    !!   vosol  = Volume (L) of aqueous solution
    !!   wfh2o  = weight (mass) fraction of water in aqueous solution
    !!   wftds  = weight (mass) fraction of solutes in aqueous solution
    !!   woh2o  = weight (mass) of solvent water
    !!   wosol  = weight (mass) of aqueous solution
    !!   wotds  = weight (mass) of total dissolved solutes
    implicit none

    ! Calling sequence variable declarations.
    integer :: nstmax

    integer :: jcsort(nstmax)

    integer :: narn1
    integer :: narn2

    logical :: qdwipp

    real(kind=8) :: mosp(nstmax)
    real(kind=8) :: mwtsp(nstmax)

    real(kind=8) :: adwipp
    real(kind=8) :: bdwipp
    real(kind=8) :: mlmrra
    real(kind=8) :: mrmlra
    real(kind=8) :: rhoc
    real(kind=8) :: rhowc
    real(kind=8) :: tdsgks
    real(kind=8) :: tdsglw
    real(kind=8) :: tdspkc
    real(kind=8) :: tdsplc
    real(kind=8) :: tempc
    real(kind=8) :: vosol
    real(kind=8) :: wfh2o
    real(kind=8) :: wftds
    real(kind=8) :: wkgwi
    real(kind=8) :: woh2o
    real(kind=8) :: wosol
    real(kind=8) :: wotds

    ! Local variable declarations.
    ! None
    ! Calculate the weights (actually, the masses) of total dissolved
    ! solutes, the solvent, and the solution, and also some related
    ! quantities.
    ! Get the total dissolved solute mass (wotds, g).
    call ctds(jcsort,mosp,mwtsp,narn1,narn2,nstmax,wotds)

    ! Get the solvent mass (fixed at 1000 grams in EQ3NR, variable
    ! in EQ6). The following works for either code.
    woh2o = mosp(narn1)*mwtsp(narn1)

    ! Get the mass of aqueous solution.
    wosol =  woh2o + wotds

    ! Get the weight fraction of solvent and the weight fraction
    ! of total dissolved solutes.
    wfh2o = woh2o/wosol
    wftds = wotds/wosol

    ! Get the inverse of the weight of solvent, in kg.
    wkgwi = 1000./woh2o

    ! Get the TDS (g/kg.solution)
    tdsgks = 1000.*wftds

    ! Get the TDS (mg/kg.solution)
    tdspkc = 1000.*tdsgks

    ! Calculate the solution density (rhowc) from the WIPP brine
    ! density model if the temperature lies in the range 20-30 C.
    ! This model is a fit to NaCl solution data having the form:
    !   density (g/L) = a + b x TDS (g/L)
    ! where a = 1000.96 and b = 0.639963. Note that this formulation
    ! requires knowing TDS in g/L (tdsglw), which prevents a direct
    ! evaluation using the above form because the density is not
    ! available in the appropriate units. Also calculate the molarity/
    ! molality conversion ratio (mrmlra, molarity = mrmlra x molality)
    ! and its inverse.
    qdwipp = .false.
    rhowc = 0.
    vosol = 0.
    tdsglw = 0.
    mrmlra = 0.
    mlmrra = 0.
    rhoc = 0.
    tdsplc = 0.

    if (tempc.ge.20. .and. tempc.le.30.) then
        ! Note that wfh2o corresponds to Nw, the "solvent fraction"
        ! in the density calculations writeup. Here wftds corresponds
        ! to "(1 - Nw)".
        qdwipp = .true.
        rhowc = adwipp/(1.0 - bdwipp*wftds)
        vosol = wosol/rhowc
        tdsglw = rhowc*wftds
        mrmlra = 0.001*wfh2o*rhowc

        if (mrmlra .gt. 0.) then
            mlmrra = 1./mrmlra
        else
            mlmrra = 1.e+38
        end if

        rhoc = 0.001*rhowc
        tdsplc = 1000.*tdsglw
    end if
end subroutine gwdenp
