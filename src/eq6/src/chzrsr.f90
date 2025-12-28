subroutine chzrsr(cbsr,elecsr,eps100,jcode,nbasp,nbt,nbtmax,nbt1mx,noutpt,nrct,nrctmx,nrndex,nsrtmx,nstmax,nttyo,ureac,uspec,zchar)
    !! This subroutine checks the reactions for special reactants to
    !! ensure that they satisfy charge balance. Reactions made by
    !! EQ6/makrsr.f should satisfy this condition, but reactions
    !! composed by the user may not. All reactions are checked by this
    !! subroutine. All special reactants are presumed to be uncharged.
    !! This subroutine is called by:
    !!   EQ6/eq6.f
    !! Principal input:
    !!   jcode  = flag array denoting the types of reactants
    !!   nbasp  = array of indices of species in the basis set
    !!   nbt    = the number of basis species
    !!   nbtmax = the maximum number of basis species
    !!   nbt1mx = the maximum number of basis species plus 1
    !!   nrctmx = the maximum number of reactants
    !!   nsrtmx = the maximum number of special reactants
    !!   ureac  = array of names of reactants
    !!   uspec  = array of names of species
    !!   zchar  = array of electrical charges of species
    !! Principal output:
    !!   cbsr   = array of coefficients for reactions for special
    !!              reactants
    !!   elecsr = array of electrical imbalances of reactions for
    !!              special reactants
    implicit none

    ! Calling sequence variable declarations.
    integer :: nbtmax
    integer :: nbt1mx
    integer :: nrctmx
    integer :: nsrtmx
    integer :: nstmax

    integer :: noutpt
    integer :: nttyo

    integer :: jcode(nrctmx)
    integer :: nbasp(nbtmax)
    integer :: nrndex(nrctmx)
    integer :: nbt
    integer :: nrct

    character(len=48) :: uspec(nstmax)
    character(len=24) :: ureac(nrctmx)

    real(kind=8) :: cbsr(nbt1mx,nsrtmx)
    real(kind=8) :: elecsr(nsrtmx)
    real(kind=8) :: zchar(nstmax)
    real(kind=8) :: eps100

    ! Local variable declarations.
    integer :: nf

    integer :: j2
    integer :: nb
    integer :: nrc
    integer :: ns
    integer :: nsr

    integer :: ilnobl

    logical :: qprrsr

    character(len=24) :: ux24

    real(kind=8) :: acx
    real(kind=8) :: aztx
    real(kind=8) :: tlersr
    real(kind=8) :: ztx

    tlersr = nbt*eps100

    do nrc = 1,nrct
        if (jcode(nrc) .eq. 2) then
            nsr = nrndex(nrc)
            ztx = 0.
            aztx = 0.

            do nb = 1,nbt
                ns = nbasp(nb)
                ztx = ztx + cbsr(nb,nsr)*zchar(ns)
                aztx = aztx + cbsr(nb,nsr)*abs(zchar(ns))
            end do

            elecsr(nsr) = ztx

            qprrsr = .false.

            if (abs(ztx) .gt. tlersr) then
                acx = abs(ztx/aztx)

                if (acx .gt. 1.e-6) then
                    j2 = ilnobl(ureac(nrc))
                    write (noutpt,1000) ureac(nrc)(1:j2)
                    write (nttyo,1000) ureac(nrc)(1:j2)
1000 format(/' * Warning - (EQ6/chzrsr) The reaction for the',' special reactant',/7x,a," isn't well charge balanced.",/7x,'The reaction is:',/)

                    qprrsr = .true.
                else if (acx .gt. 1.e-10) then
                    j2 = ilnobl(ureac(nrc))
                    write (noutpt,1010) ureac(nrc)(1:j2)
                    write (nttyo,1010) ureac(nrc)(1:j2)
1010 format(/' * Warning - (EQ6/chzrsr) The reaction for the',' special reactant',/7x,a," isn't charge balanced to",' high precision.',/7x,'The reaction is:',/)

                    qprrsr = .true.
                end if
            end if

            if (qprrsr) then
                nf = noutpt

                ! Calling sequence substitutions:
                !   nbasp for nbaspd
                !   nbt for nbtd
                call prrsr(cbsr,jcode,nbasp,nbt,nbtmax,nbt1mx,nf,noutpt,nrc,nrctmx,nrndex,nsrtmx,nstmax,nttyo,ureac,uspec)

                ! Calling sequence substitutions:
                !   nbasp for nbaspd
                !   nbt for nbtd
                nf = nttyo
                call prrsr(cbsr,jcode,nbasp,nbt,nbtmax,nbt1mx,nf,noutpt,nrc,nrctmx,nrndex,nsrtmx,nstmax,nttyo,ureac,uspec)

                ux24 = ' '
                write (ux24,'(1pe12.5)') ztx
                j2 = ilnobl(ux24)
                write (noutpt,1050) ux24(1:j2)
                write (nttyo,1050) ux24(1:j2)
1050 format(/7x,'The imbalance is ',a,'. This will be factored',' into the computed',/7x,'charge balance a a compensating',' offset. The offset at any point',/7x,'will be',' proportional to the extent of reaction. The pH and',/7x,'redox will be properly calculated, as they are',' computed entirely',/7x,'from mass balance relations for',' the data file basis species.',/7x,'This version of EQ6',' does not compute any of these properties using',/7x,'the charge balance equation.')
            end if
        end if
    end do
end subroutine chzrsr