subroutine prrsr(cbsr,jcode,nbaspd,nbtd,nbtmax,nbt1mx,nf,noutpt,nrc,nrctmx,nrndex,nsrtmx,nstmax,nttyo,ureac,uspec)
    !! This subroutine writes the reaction for the nsr-th special
    !! reactant on the file whose unit number is nf.
    !! This subroutine is similar in function to EQLIB/prreac.f, which
    !! writes an ordinary reaction.
    !! This subroutine is called by:
    !!   EQ6/chzrsr.f
    !!   EQ6/ckfrsr.f
    !!   EQ6/echoz.f
    !! Principal input:
    !!   cbsr   = array of coefficients for reactions for special
    !!              reactants
    !!   nbtmax = the maximum number of basis species
    !!   nbt1mx = the maximum number of basis species plus 1
    !!   nf     = the unit number of the file on which to write the
    !!              reaction
    !!   nrc    = the reactant index of the special reactant whose
    !!              reaction is to be printed
    !!   nsrtmx = the maximum number of special reactants
    !!   ureac  = array of names of reactants
    !!   uspec  = array of names of species
    !! Principal output:
    !!   None
    implicit none

    ! Calling sequence variable declarations.
    integer :: nf
    integer :: noutpt
    integer :: nttyo

    integer :: nbtmax
    integer :: nbt1mx
    integer :: nrctmx
    integer :: nsrtmx
    integer :: nstmax

    integer :: jcode(nrctmx)
    integer :: nbaspd(nbtmax)
    integer :: nrndex(nrctmx)

    integer :: nbtd
    integer :: nrc

    character(len=48) :: uspec(nstmax)
    character(len=24) :: ureac(nrctmx)

    real(kind=8) :: cbsr(nbt1mx,nsrtmx)

    ! Local variable declarations.
    integer :: jlen
    integer :: j2
    integer :: nb
    integer :: ns
    integer :: nsr

    integer :: ilnobl

    logical :: qfirst

    character(len=56) :: uspn56

    real(kind=8) :: cx

    ! Make sure that the reactant is a special reactant.
    if (jcode(nrc) .ne. 2) then
        j2 = ilnobl(ureac(nrc))
        write (noutpt,1000) ureac(nrc)(1:j2)
        write (nttyo,1000) ureac(nrc)(1:j2)
1000 format(/' * Error - (EQ6/prrsr.f) Programming error trap: this',' subroutine was called',/7x,'to print the reaction for the',' special reactant ',a,'. However,',/7x,'this reactant is not',' a special reactant.')

        stop
    end if

    ! Get the special reactant index.
    nsr = nrndex(nrc)

    write (nf,1020)
1020 format(1x)

    ! The reactant.
    cx = -cbsr(nbt1mx,nsr)
    j2 = ilnobl(ureac(nrc))
    write (nf,1030) cx,ureac(nrc)(1:j2)
1030 format(6x,1pe22.15,2x,a)

    ! Reactants exclusive of the special reactant.
    do nb = 1,nbtd
        cx = cbsr(nb,nsr)

        if (cx .lt. 0.) then
            cx = -cx
            ns = nbaspd(nb)

            ! Calling sequence substitutions:
            !   uspec(ns) for unam48
            call fmspnx(jlen,uspec(ns),uspn56)
            write (nf,1040) cx,uspn56(1:jlen)
1040 format(4x,'+ ',1pe22.15,2x,a)
        end if
    end do

    write (nf,1050)
1050 format(18x,'==')

    ! Products.
    qfirst = .true.

    do nb = 1,nbtd
        cx = cbsr(nb,nsr)

        if (cx .gt. 0.) then
            ns = nbaspd(nb)

            ! Calling sequence substitutions:
            !   uspec(ns) for unam48
            call fmspnx(jlen,uspec(ns),uspn56)

            if (qfirst) then
                write (nf,1030) cx,uspn56(1:jlen)
                qfirst = .false.
            else
                write (nf,1040) cx,uspn56(1:jlen)
            end if
        end if
    end do
end subroutine prrsr
