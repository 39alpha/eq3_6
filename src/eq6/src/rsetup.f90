subroutine rsetup(atwt,cbsr,cesr,iern1,ietmax,iindx1,iktmax,jcode,jern1,jern2,jetmax,jgext,kbt,kmax,mwtges,mwtrc,mwtsp,nbaspd,nbt,nbtmax,nbt1mx,ncmpr,nct,nctmax,nertmx,netmax,ngext,noutpt,nptmax,nrct,nrctmx,nrndex,nsrtmx,nstmax,nsts,nstsmx,nstsr,nttyo,nxridx,nxrtmx,rxbar,ureac,uspec,vosp0,vreac,xgers)
    !! This subroutine assigns the molecular weights and molar volumes
    !! of the reactants. Volumes are treated only for solid reactants
    !! for later comparison with the volume of solid products.
    !!    nrct = total number of reactants
    !!    ureac = name of reactant
    !!    mwtrc = molecular weight of reactant
    !!    vreac = molar volume of (solid) reactant
    !!    jcode = reactant type code
    !!            0   Mineral
    !!            1   Solid solution
    !!            2   Special reactant
    !!            3   Aqueous species
    !!            4   Gas
    !!            5   Generic ion exchanger
    !! This subroutine is called by:
    !!   EQ6/eq6.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: ietmax
    integer :: iktmax
    integer :: jetmax
    integer :: kmax
    integer :: nbtmax
    integer :: nbt1mx
    integer :: nctmax
    integer :: nertmx
    integer :: netmax
    integer :: nptmax
    integer :: nrctmx
    integer :: nsrtmx
    integer :: nstmax
    integer :: nstsmx
    integer :: nxrtmx

    integer :: noutpt
    integer :: nttyo

    integer :: iindx1(kmax)
    integer :: jcode(nrctmx)
    integer :: jern1(jetmax,netmax)
    integer :: jern2(jetmax,netmax)
    integer :: jgext(netmax)
    integer :: nbaspd(nbtmax)
    integer :: ncmpr(2,nptmax)
    integer :: ngext(jetmax,netmax)
    integer :: nrndex(nrctmx)
    integer :: nsts(nstsmx)
    integer :: nstsr(2,nstmax)
    integer :: nxridx(nrctmx)

    integer :: iern1
    integer :: kbt
    integer :: nbt
    integer :: nct
    integer :: nrct

    character(len=48) :: uspec(nstmax)
    character(len=24) :: ureac(nrctmx)

    real(kind=8) :: atwt(nctmax)
    real(kind=8) :: cbsr(nbt1mx,nsrtmx)
    real(kind=8) :: cesr(nctmax,nsrtmx)
    real(kind=8) :: mwtges(netmax)
    real(kind=8) :: mwtrc(nrctmx)
    real(kind=8) :: mwtsp(nstmax)
    real(kind=8) :: rxbar(iktmax,nxrtmx)
    real(kind=8) :: vosp0(nstmax)
    real(kind=8) :: vreac(nrctmx)
    real(kind=8) :: xgers(ietmax,jetmax,nertmx)

    ! Local variable declarations.
    integer :: ie
    integer :: ik
    integer :: je
    integer :: j2
    integer :: jlen
    integer :: kcol
    integer :: n
    integer :: nb
    integer :: nbb
    integer :: nc
    integer :: ne
    integer :: ner
    integer :: np
    integer :: nrc
    integer :: nrn1
    integer :: nrn2
    integer :: nr1
    integer :: nr2
    integer :: ns
    integer :: nsr
    integer :: nss
    integer :: nxr

    integer :: ilnobl

    logical :: qstop

    character(len=56) :: uspn56

    qstop = .false.

    do nrc = 1,nrct
        if (jcode(nrc) .eq. 0) then
            ! Pure mineral reactants.
            np = nrndex(nrc)
            ns = ncmpr(1,np)
            mwtrc(nrc) = mwtsp(ns)
            vreac(nrc) = vosp0(ns)

            nr1 = nstsr(1,ns)
            nr2 = nstsr(1,ns)

            do n = nr1,nr2
                nb = nsts(n)

                do kcol = 1,kbt
                    nbb = iindx1(kcol)

                    if (nbb .eq. nb) then
                        go to 100
                    end if
                end do

                nss = nbaspd(nb)
                j2 = ilnobl(ureac(nrc))

                ! Calling sequence substitutions:
                !   uspec(nss) for unam48
                call fmspnx(jlen,uspec(nss),uspn56)

                write (noutpt,2100) ureac(nrc)(1:j2),uspn56(1:jlen)
                write (nttyo,2100) ureac(nrc)(1:j2),uspn56(1:jlen)
2100 format(/' * Error - (EQ6/rsetup) The reactant ',a,' is',' composed of',/7x,'the basis species ',a,', which',/7x,'is not in the active basis set.')

                qstop = .true.
100 continue
            end do
        else if (jcode(nrc) .eq. 1) then
            ! Solid solution reactants.
            np = nrndex(nrc)
            nrn1 = ncmpr(1,np)
            nrn2 = ncmpr(2,np)
            nxr = nxridx(nrc)
            mwtrc(nrc) = 0.
            vreac(nrc) = 0.

            ik = 0

            do ns = nrn1,nrn2
                ik = ik + 1
                mwtrc(nrc) = mwtrc(nrc) + rxbar(ik,nxr)*mwtsp(ns)
                vreac(nrc) = vreac(nrc) + rxbar(ik,nxr)*vosp0(ns)

                nr1 = nstsr(1,ns)
                nr2 = nstsr(1,ns)

                do n = nr1,nr2
                    nb = nsts(n)

                    do kcol = 1,kbt
                        nbb = iindx1(kcol)

                        if (nbb .eq. nb) then
                            go to 110
                        end if
                    end do

                    nss = nbaspd(nb)
                    j2 = ilnobl(ureac(nrc))

                    ! Calling sequence substitutions:
                    !   uspec(nss) for unam48
                    call fmspnx(jlen,uspec(nss),uspn56)

                    write (noutpt,2100) ureac(nrc)(1:j2),uspn56(1:jlen)
                    write (nttyo,2100) ureac(nrc)(1:j2),uspn56(1:jlen)
                    qstop = .true.
110 continue
                end do
            end do
        else if (jcode(nrc) .eq. 2) then
            ! Special reactants.
            mwtrc(nrc) = 0.
            nsr = nrndex(nrc)

            do nc = 1,nct
                if (cesr(nc,nsr) .ne. 0.) then
                    mwtrc(nrc) = mwtrc(nrc) + atwt(nc)*cesr(nc,nsr)
                end if
            end do

            vreac(nrc) = 0.

            do nb = 1,nbt
                if (cbsr(nb,nsr) .ne. 0.) then
                    do kcol = 1,kbt
                        nbb = iindx1(kcol)

                        if (nbb .eq. nb) then
                            go to 120
                        end if
                    end do

                    nss = nbaspd(nb)
                    j2 = ilnobl(ureac(nrc))

                    ! Calling sequence substitutions:
                    !   uspec(nss) for unam48
                    call fmspnx(jlen,uspec(nss),uspn56)

                    write (noutpt,2100) ureac(nrc)(1:j2),uspn56(1:jlen)
                    write (nttyo,2100) ureac(nrc)(1:j2),uspn56(1:jlen)
                    qstop = .true.
                end if

120 continue
            end do
        else if (jcode(nrc).eq.3 .or. jcode(nrc).eq.4) then
            ! Aqueous species reactants and gas reactants.
            ns = nrndex(nrc)
            mwtrc(nrc) = mwtsp(ns)
            vreac(nrc) = 0.

            nr1 = nstsr(1,ns)
            nr2 = nstsr(1,ns)

            do n = nr1,nr2
                nb = nsts(n)

                do kcol = 1,kbt
                    nbb = iindx1(kcol)

                    if (nbb .eq. nb) then
                        go to 130
                    end if
                end do

                nss = nbaspd(nb)
                j2 = ilnobl(ureac(nrc))

                ! Calling sequence substitutions:
                !   uspec(nss) for unam48
                call fmspnx(jlen,uspec(nss),uspn56)

                write (noutpt,2100) ureac(nrc)(1:j2),uspn56(1:jlen)
                write (nttyo,2100) ureac(nrc)(1:j2),uspn56(1:jlen)
                qstop = .true.
130 continue
            end do
        else if (jcode(nrc) .eq. 5) then
            ! Generic ion exchangers.
            np = nrndex(nrc)
            ner = nxridx(nrc)
            ne = np - iern1 + 1

            mwtrc(nrc) = mwtges(ne)
            vreac(nrc) = 0.

            do je = 1,jgext(ne)
                ns = jern1(je,ne) - 1

                do ie = 1,ngext(je,ne)
                    ns = ns + 1
                    mwtrc(nrc) = mwtrc(nrc) + xgers(ie,je,ner)*mwtsp(ns)
                end do
            end do

            do je = 1,jgext(ne)
                do ns = jern1(je,ne),jern2(je,ne)
                    nr1 = nstsr(1,ns)
                    nr2 = nstsr(1,ns)

                    do n = nr1,nr2
                        nb = nsts(n)

                        do kcol = 1,kbt
                            nbb = iindx1(kcol)

                            if (nbb .eq. nb) then
                                go to 140
                            end if
                        end do

                        nss = nbaspd(nb)
                        j2 = ilnobl(ureac(nrc))

                        ! Calling sequence substitutions:
                        !   uspec(nss) for unam48
                        call fmspnx(jlen,uspec(nss),uspn56)

                        write (noutpt,2100) ureac(nrc)(1:j2),uspn56(1:jlen)
                        write (nttyo,2100) ureac(nrc)(1:j2),uspn56(1:jlen)
                        qstop = .true.

140 continue
                    end do
                end do
            end do
        end if
    end do

    if (qstop) then
        stop
    end if
end subroutine rsetup