subroutine sortsp(iern1,iern2,istack,jcsort,jern1,jern2,jgext,jsitex,jetmax,jjsort,jssort,jstack,losp,lsort,ncmpr,nern1,nern2,netmax,noutpt,nphasx,npt,nptmax,nst,nstmax,nttyo)
    !! This subroutine sorts the species in order of increasing mass.
    !! This subroutine is called by:
    !!   EQLIB/ncmpex.f
    !! Principal input:
    !!   jern1  = array marking the start of a site range for a
    !!              generic ion exchange phase
    !!   jern2  = array marking the end of a site range for a
    !!              generic ion exchange phase
    !!   iern1  = start of the phase range for generic ion exchange
    !!              phases
    !!   iern2  = end of the phase range for generic ion exchange
    !!              phases
    !!   jgext  = array giving the number of sites in generic
    !!              generic ion exchange phases
    !!   jsitex = site index array, gives the site to which a
    !!              species belongs
    !!   losp   = array of logarithms of numbers of moles of species
    !!   ncmpr  = species index array, gives the first and last
    !!              species in a given phase
    !!   nern1  = start of the species range for generic ion exchange
    !!              phases
    !!   nern2  = end of the species range for generic ion exchange
    !!              phases
    !!   nphasx = phase index array, gives the phase to which a
    !!              species belongs
    !!   nst    = number of species
    !! Principal output:
    !!   jcsort = species index array, arranged in order of increasing
    !!              mass, but arranged according to phase
    !!   jjsort = species index array, arranged in order of increasing
    !!              mass, but arranged according to site and phase
    !!   jssort = species index array, arranged in order of increasing
    !!              mass
    !! Scratch arrays:
    !!   lsort, istack, and jstack
    implicit none

    ! Calling sequence variable declarations:
    integer :: jetmax
    integer :: netmax
    integer :: nptmax
    integer :: nstmax

    integer :: noutpt
    integer :: nttyo

    integer :: iern1
    integer :: iern2
    integer :: nern1
    integer :: nern2

    integer :: istack(nstmax)
    integer :: jcsort(nstmax)
    integer :: jern1(jetmax,netmax)
    integer :: jern2(jetmax,netmax)
    integer :: jgext(netmax)
    integer :: jjsort(nstmax)
    integer :: jsitex(nstmax)
    integer :: jssort(nstmax)
    integer :: jstack(nstmax)
    integer :: ncmpr(2,nptmax)
    integer :: nphasx(nstmax)

    integer :: npt
    integer :: nst

    real(kind=8) :: losp(nstmax)
    real(kind=8) :: lsort(nstmax)

    ! Local variable declarations.
    integer :: je
    integer :: ne
    integer :: np
    integer :: nrr1
    integer :: nrr2
    integer :: ns
    integer :: nss
    integer :: nsi

    ! Sort all species according to numbers of moles. There indices
    ! in sorted order in the jssort array. Species belonging to
    ! different phases are mixed together in this sort.
    ! Calling sequence substitutions:
    !   lsort for asort
    !   losp for aval
    !   jssort for jsort
    !   nstmax for nmax
    !   nst for nval
    ! Caution: the jssort array from the last call is recycled as a
    ! good starting point. Set jssort(1) to 0 to make a sort starting
    ! from scratch.
    call qsortw(lsort,losp,istack,jssort,jstack,nstmax,noutpt,nttyo,nst)

    ! Compute the jcsort array using the jssort array already computed.
    ! In this array, the species belonging to a given phase are grouped
    ! together within the specified phase limits. The scratch array
    ! istack is used here to represent the next available species
    ! position within the range for any given phase. Note that istack
    ! is dimensioned 1:nstmax. Here it must be of length at least
    ! nptmax, which is necessarily no greater than nstmax.
    do np = 1,npt
        istack(np) = ncmpr(1,np)
    end do

    do nss = 1,nst
        ns = jssort(nss)
        np = nphasx(ns)
        nsi = istack(np)
        jcsort(nsi) = ns
        istack(np) = istack(np) + 1
    end do

    ! Compute the jjsort array using the jcsort array already computed.
    ! In this array, the species belonging to a given phase are grouped
    ! together within the specified site limits. Presently only generic
    ! ion exchange phases have sites in the present context. In the
    ! future, other types of phases with distinct sites may be added
    ! to the list. Phases without distinct sites are treated as having
    ! "one" site.
    ! The scratch array istack is used here to represent the next
    ! available species position within the range for any given site
    ! of the phase currently being processed. any given phase. Note
    ! that istack is dimensioned 1:nstmax. Here it must have a length
    ! of at least the maximum number of species in any possible phase,
    ! which is necessarily no greater than nstmax.
    do ns = 1,nst
        jjsort(ns) = jcsort(ns)
    end do

    do ns = nern1,nern2
        jjsort(ns) = 0
    end do

    do np = iern1,iern2
        ne = np - iern1 + 1

        do je = 1,jgext(ne)
            istack(je) = jern1(je,ne)
        end do

        nrr1 = ncmpr(1,np)
        nrr2 = ncmpr(2,np)

        do nss = nrr1,nrr2
            ns = jcsort(nss)
            je = jsitex(ns)
            nsi = istack(je)
            jjsort(nsi) = ns
            istack(je) = istack(je) + 1
        end do
    end do
end subroutine sortsp