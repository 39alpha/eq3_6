subroutine reacts(cbsr,csts,delxi,drer0,iern1,ietmax,iktmax,iodb,jcode,jetmax,jgext,jreac,modr,modr0,morr,morr0,mrgers,mtb,mtb0,nbaspd,nbt,nbtmax,nbt1mx,ncmpr,nern1,nern2,nertmx,netmax,ngext,nodbmx,nord,noutpt,nptmax,nrct,nrctmx,nrd1mx,nrndex,nsrtmx,nstmax,nsts,nstsmx,nstsr,nttyo,nxridx,nxrtmx,rrelr0,rxbar,ureac,xirct,xirct0)
    !! This subroutine computes the destroyed and current masses of the
    !! reactants and the current mass balance totals for the equilibrium
    !! system.
    !!    morr   = moles of reactant remaining
    !!    modr   = moles of irreversibly destroyed reactant
    !!    jreac  = reactant control switch
    !!       = -1  The reactant has saturated but continues to
    !!             be available for irreversible reaction
    !!       =  0  The reactant is currently reacting
    !!       =  1  The reactant has been exhausted
    !!       =  2  The reactant has saturated and any remaining mass
    !!             has been transferred to the equilibrium system
    !!             (This occurs only if iopt(1) = 0)
    !! This subroutine is called by:
    !!   EQ6/eqshel.f
    !!   EQ6/path.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: ietmax
    integer :: iktmax
    integer :: jetmax
    integer :: nbtmax
    integer :: nbt1mx
    integer :: nertmx
    integer :: netmax
    integer :: nodbmx
    integer :: nptmax
    integer :: nrctmx
    integer :: nrd1mx
    integer :: nsrtmx
    integer :: nstmax
    integer :: nstsmx
    integer :: nxrtmx

    integer :: iodb(nodbmx)
    integer :: jcode(nrctmx)
    integer :: jgext(netmax)
    integer :: jreac(nrctmx)
    integer :: nbaspd(nbtmax)
    integer :: ncmpr(2,nptmax)
    integer :: ngext(jetmax,netmax)
    integer :: nrndex(nrctmx)
    integer :: nsts(nstsmx)
    integer :: nstsr(2,nstmax)
    integer :: nxridx(nrctmx)

    integer :: iern1
    integer :: nbt
    integer :: nern1
    integer :: nern2
    integer :: nord
    integer :: noutpt
    integer :: nrct
    integer :: nttyo

    character(len=24) :: ureac(nrctmx)

    real(kind=8) :: cbsr(nbt1mx,nsrtmx)
    real(kind=8) :: csts(nstsmx)
    real(kind=8) :: drer0(nrd1mx,nrctmx)
    real(kind=8) :: modr(nrctmx)
    real(kind=8) :: modr0(nrctmx)
    real(kind=8) :: morr(nrctmx)
    real(kind=8) :: morr0(nrctmx)
    real(kind=8) :: mrgers(ietmax,jetmax,nertmx)
    real(kind=8) :: mtb(nbtmax)
    real(kind=8) :: mtb0(nbtmax)
    real(kind=8) :: rrelr0(nrctmx)
    real(kind=8) :: rxbar(iktmax,nxrtmx)
    real(kind=8) :: xirct(nrctmx)
    real(kind=8) :: xirct0(nrctmx)

    real(kind=8) :: delxi

    ! Local variable declarations.
    integer :: ie
    integer :: ik
    integer :: je
    integer :: j2
    integer :: n
    integer :: nb
    integer :: ne
    integer :: ner
    integer :: np
    integer :: nsr
    integer :: nrc
    integer :: nrn1
    integer :: nrn2
    integer :: nr1
    integer :: nr2
    integer :: ns
    integer :: nxr

    integer :: ilnobl

    real(kind=8) :: dlxrct
    real(kind=8) :: dmodr
    real(kind=8) :: dx

    ! Reset morr, modr, xirct, and mtb to values at the previous step.
    do nb = 1,nbt
        mtb(nb) = mtb0(nb)
    end do

    do nrc = 1,nrct
        morr(nrc) = morr0(nrc)
        modr(nrc) = modr0(nrc)
        xirct(nrc) = xirct0(nrc)
    end do

    if (delxi .le. 0.) then
        go to 999
    end if

    ! Increment morr, modr, xirct, and mtb to correspond to the
    ! advance in reaction progress.
    if (iodb(1) .ge. 3) then
        write (noutpt,1000)
    end if

1000 format(' --- Advancement of Irreversible Reactions ---',/7x,'Reaction',14x,'dlxrct',7x,'dmodr')

    do nrc = 1,nrct
        if (jreac(nrc).eq.0 .or. jreac(nrc).eq.-1) then
            call integr(delxi,dlxrct,drer0,nord,nrc,nrctmx,nrd1mx,rrelr0)

            ! XX       It is assumed that each reactant has a reaction coefficient
            ! XX       of -1. Below it is implied that dlxrct is multiplied by -(-1).
            dmodr = dlxrct

            if (dmodr .gt. morr(nrc)) then
                ! The amount destroyed in the current step slightly exceeds
                ! the amount remaining.
                dmodr = morr(nrc)

                ! XX         It is assumed that each reactant has a reaction coefficient
                ! XX         of -1. Below it is implied that dlxrct is multiplied
                ! XX         by -(-1).
                dlxrct = dmodr
            end if

            if (iodb(1) .ge. 3) then
                write (noutpt,1010) ureac(nrc),dlxrct,dmodr
            end if

1010 format(3x,a24,2(3x,1pe12.5))

            if (dlxrct .ne. 0.) then
                xirct(nrc) = xirct(nrc) + dlxrct
                modr(nrc) = modr(nrc) + dmodr
                morr(nrc) = morr(nrc) - dmodr

                if (jcode(nrc) .eq. 0) then
                    ! Pure mineral.
                    np = nrndex(nrc)
                    ns = ncmpr(1,np)
                    nr1 = nstsr(1,ns)
                    nr2 = nstsr(2,ns)

                    do n = nr1,nr2
                        nb = nsts(n)
                        mtb(nb) = mtb(nb) + dmodr*csts(n)
                    end do
                else if (jcode(nrc) .eq. 1) then
                    ! Solid solution.
                    nxr = nxridx(nrc)
                    np = nrndex(nrc)
                    nrn1 = ncmpr(1,np)
                    nrn2 = ncmpr(2,np)
                    ik = 0

                    do ns = nrn1,nrn2
                        ik = ik + 1
                        nr1 = nstsr(1,ns)
                        nr2 = nstsr(2,ns)
                        dx = dmodr*rxbar(ik,nxr)

                        do n = nr1,nr2
                            nb = nsts(n)
                            mtb(nb) = mtb(nb) + dx*csts(n)
                        end do
                    end do
                else if (jcode(nrc) .eq. 2) then
                    ! Special reactant.
                    nsr = nrndex(nrc)

                    do nb = 1,nbt
                        mtb(nb) = mtb(nb) + dmodr*cbsr(nb,nsr)
                    end do
                else if (jcode(nrc).eq.3 .or. jcode(nrc).eq.4) then
                    ! Aqueous species or gas.
                    ns = nrndex(nrc)
                    nr1 = nstsr(1,ns)
                    nr2 = nstsr(2,ns)

                    do n = nr1,nr2
                        nb = nsts(n)
                        mtb(nb) = mtb(nb) + dmodr*csts(n)
                    end do
                else if (jcode(nrc) .eq. 5) then
                    ! Generic ion exchanger.
                    ner = nxridx(nrc)
                    np = nrndex(nrc)
                    ne = np - iern1 + 1
                    ns = ncmpr(1,np) - 1

                    do je = 1,jgext(ne)
                        do ie = 1,ngext(je,ne)
                            dx = dmodr*mrgers(ie,je,ner)
                            ns = ns + 1
                            nr1 = nstsr(1,ns)
                            nr2 = nstsr(2,ns)

                            if (je .le. 1) then
                                ! Have the first site. Increment the mass balance
                                ! totals in a straightforward manner.
                                do n = nr1,nr2
                                    nb = nsts(n)
                                    mtb(nb) = mtb(nb) + dx*csts(n)
                                end do
                            else
                                ! Have a site beyond the first. Increment the mass
                                ! balance totals in the usual manner, except do not
                                ! increment here the total for the exchanger
                                ! substrate. The complete increment for the substrate
                                ! is obtained by considering only one site.
                                do n = nr1,nr2
                                    nb = nsts(n)
                                    ns = nbaspd(nb)

                                    if (ns.lt.nern1 .or. ns.gt.nern2) then
                                        mtb(nb) = mtb(nb) + dx*csts(n)
                                    end if
                                end do
                            end if
                        end do
                    end do
                else
                    j2 = ilnobl(ureac(nrc))
                    write (noutpt,64) jcode(nrc),ureac(nrc)(1:j2)
                    write (nttyo,64) jcode(nrc),ureac(nrc)(1:j2)
64 format(/' * Error - (EQ6/reacts) Programming error',' trap: Have unknown',/7x,'reactant type code (jcode)',' value of ',i5,' for reactant ',/7x,a,'.')

                    stop
                end if
            end if
        end if
    end do

999 continue
end subroutine reacts