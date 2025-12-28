subroutine ckfrsr(cbsr,csts,jcode,jflag,nbaspd,nbtd,nbtmax,nbt1mx,noutpt,nrct,nrctmx,nrndex,nsrtmx,nstmax,nsts,nstsmx,nstsr,nttyo,ureac,uspec)
    !! This subroutine checks the reactions for special reactants. As
    !! needed, it rewrites these reactions to eliminate any basis
    !! species for which jflag = 30. Such species are not in the active
    !! basis set, and there are no corresonding mass balances.
    !! This subroutine is called by:
    !!   EQ6/eq6.f
    !! Principal input:
    !!   cbsr   = array of coefficients for reactions for special
    !!              reactants
    !!   nbtmax = the maximum number of basis species
    !!   nbt1mx = the maximum number of basis species plus 1
    !!   nsrtmx = the maximum number of special reactants
    !!   ureac  = array of names of reactants
    !!   uspec  = array of names of species
    !! Principal output:
    !!   cbsr   = array of coefficients for reactions for special
    !!              reactants (modified if necessary)
    implicit none

    ! Calling sequence variable declarations.
    integer :: nbtmax
    integer :: nbt1mx
    integer :: nrctmx
    integer :: nsrtmx
    integer :: nstmax
    integer :: nstsmx

    integer :: noutpt
    integer :: nttyo

    integer :: jcode(nrctmx)
    integer :: jflag(nstmax)
    integer :: nbaspd(nbtmax)
    integer :: nrndex(nrctmx)
    integer :: nsts(nstsmx)
    integer :: nstsr(2,nstmax)

    integer :: nbtd
    integer :: nrct

    character(len=48) :: uspec(nstmax)
    character(len=24) :: ureac(nrctmx)

    real(kind=8) :: cbsr(nbt1mx,nsrtmx)
    real(kind=8) :: csts(nstsmx)

    ! Local variable declarations.
    integer :: nf

    integer :: jlen
    integer :: j2
    integer :: n
    integer :: nb
    integer :: nbb
    integer :: nrc
    integer :: nr1
    integer :: nr2
    integer :: ns
    integer :: nsr

    integer :: ilnobl

    logical :: qcaught

    character(len=56) :: uspn56

    real(kind=8) :: cx

    do nrc = 1,nrct
        if (jcode(nrc) .eq. 2) then
            nsr = nrndex(nrc)
            qcaught = .false.

            do nb = 1,nbtd
                ns = nbaspd(nb)
                cx = cbsr(nb,nsr)

                if (cx.ne.0. .and. jflag(ns).eq.30) then
                    if (.not.qcaught) then
                        ! Write the original reaction.
                        j2 = ilnobl(ureac(nrc))
                        write (noutpt,1000) ureac(nrc)(1:j2)
                        write (nttyo,1000) ureac(nrc)(1:j2)
1000 format(/1x,'The reaction for special reactant ',a,' is written in terms of',/1x,'one or more species',' which are not in the active basis set. The',/1x,'currently written reaction is:')

                        nf = noutpt
                        call prrsr(cbsr,jcode,nbaspd,nbtd,nbtmax,nbt1mx,nf,noutpt,nrc,nrctmx,nrndex,nsrtmx,nstmax,nttyo,ureac,uspec)
                        write (noutpt,1010)
1010 format(/1x,'The following species are not in the',' active basis set and will be eliminated',/1x,'from the reaction:',/)

                        nf = nttyo
                        call prrsr(cbsr,jcode,nbaspd,nbtd,nbtmax,nbt1mx,nf,noutpt,nrc,nrctmx,nrndex,nsrtmx,nstmax,nttyo,ureac,uspec)
                        write (nttyo,1010)

                        qcaught = .true.
                    end if

                    ! Calling sequence substitutions:
                    !   uspec(ns) for unam48
                    call fmspnx(jlen,uspec(ns),uspn56)

                    write (noutpt,1020) uspn56(1:jlen)
                    write (nttyo,1020) uspn56(1:jlen)
1020 format(3x,a)

                    ! The reaction contains a non-zero reaction coefficient for
                    ! a species which is not in the active basis set. Rewrite
                    ! the reaction so that this species does not appear.
                    nr1 = nstsr(1,ns)
                    nr2 = nstsr(2,ns)

                    do n = nr1,nr2
                        nbb = nsts(n)
                        cbsr(nbb,nsr) = cbsr(nbb,nsr) + csts(n)*cx
                    end do

                    cbsr(nb,nsr) = 0.
                end if
            end do

            if (qcaught) then
                ! Write the modified reaction.
                write (noutpt,1070)
                write (nttyo,1070)
1070 format(//1x,'The modified reaction is:')

                nf = noutpt
                call prrsr(cbsr,jcode,nbaspd,nbtd,nbtmax,nbt1mx,nf,noutpt,nrc,nrctmx,nrndex,nsrtmx,nstmax,nttyo,ureac,uspec)
                write (noutpt,1030)
1030 format(/1x)

                nf = nttyo
                call prrsr(cbsr,jcode,nbaspd,nbtd,nbtmax,nbt1mx,nf,noutpt,nrc,nrctmx,nrndex,nsrtmx,nstmax,nttyo,ureac,uspec)
                write (nttyo,1030)
            end if
        end if
    end do
end subroutine ckfrsr