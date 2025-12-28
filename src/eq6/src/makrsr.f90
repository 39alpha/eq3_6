subroutine makrsr(cbsri,cesri,cess,eps100,ibsrti,iesrti,jcode,nbt1mx,nct,nctmax,ness,nessmx,nessr,noutpt,nrct,nrctmx,nsrtmx,nstmax,nttyo,ubsri,uelem,uesri,ureac,uspec,zchar)
    !! This subroutine makes a reaction for a special reactant. The
    !! reaction is written in terms of the special reactant and strict
    !! basis species only, and is charge balanced.
    !! This subroutine is called by:
    !!   EQ6/eq6.f
    !! Principal input:
    !!   cesri  = array of coefficients for chemical elements composing
    !!              special reactants
    !!   ibsrti = array of numbers of species in reactions for special
    !!              reactants
    !!   iesrti = array of numbers of chemical elements composing special
    !!              reactants
    !!   jcode  = flag array denoting the types of reactants
    !!   nbt1mx = the maximum number of basis species plus 1
    !!   nct    = the number of chemical elements
    !!   nctmax = the maximum number of chemical elements
    !!   ness   = array of indices of chemical elements composing
    !!              species
    !!   nessmx = the maximum number of entries in the cess or ness array
    !!   nessr  = pointer array giving the ranges in the cess and ness
    !!              arrays giving the compositions of species
    !!   nrctmx = the maximum number of reactants
    !!   nsrtmx = the maximum number of special reactants
    !!   uelem  = array of names of chemical elements
    !!   uesri  = array of names of chemical elements composing special
    !!              reactants
    !!   ureac  = array of names of reactants
    !!   uspec  = array of names of species
    !!   zchar  = array of electrical charges of species
    !! Principal output:
    !!   cbsri  = array of coefficients for reactions for special
    !!              reactants
    !!   ibsrti = array of numbers of species in reactions for special
    !!              reactants
    !!   ubsri  = array of names of species in reactions for special
    !!              reactants
    implicit none

    ! Calling sequence variable declarations.
    integer :: nbt1mx
    integer :: nctmax
    integer :: nessmx
    integer :: nrctmx
    integer :: nsrtmx
    integer :: nstmax

    integer :: noutpt
    integer :: nttyo

    integer :: ibsrti(nsrtmx)
    integer :: iesrti(nsrtmx)
    integer :: jcode(nrctmx)
    integer :: ness(nessmx)
    integer :: nessr(2,nstmax)

    integer :: nct
    integer :: nrct

    character(len=48) :: uspec(nstmax)
    character(len=24) :: ubsri(nbt1mx,nsrtmx)
    character(len=24) :: ureac(nrctmx)
    character(len=8) :: uelem(nctmax)
    character(len=8) :: uesri(nctmax,nsrtmx)

    real(kind=8) :: cbsri(nbt1mx,nsrtmx)
    real(kind=8) :: cess(nessmx)
    real(kind=8) :: cesri(nctmax,nsrtmx)
    real(kind=8) :: zchar(nstmax)

    real(kind=8) :: eps100

    ! Local variable declarations.
    integer :: j2
    integer :: j3
    integer :: n
    integer :: nc
    integer :: nch
    integer :: nci
    integer :: ncih
    integer :: ncio
    integer :: nco
    integer :: nerr
    integer :: nn
    integer :: nrc
    integer :: nr1
    integer :: nr2
    integer :: ns
    integer :: nsr

    integer :: ilnobl

    real(kind=8) :: ctxo
    real(kind=8) :: ctxh
    real(kind=8) :: cx
    real(kind=8) :: cxe
    real(kind=8) :: cxh
    real(kind=8) :: cxo
    real(kind=8) :: ztx

    nerr = 0

    nco = 0

    do nc = 1,nct
        if (uelem(nc)(1:2) .eq. 'O ') then
            nco = nc
            go to 100
        end if
    end do

100 continue

    nch = 0

    do nc = 1,nct
        if (uelem(nc)(1:2) .eq. 'H ') then
            nch = nc
            go to 110
        end if
    end do

110 continue

    nsr = 0

    do nrc = 1,nrct
        if (jcode(nrc) .eq. 2) then
            nsr = nsr + 1

            if (ibsrti(nsr) .eq. 0) then
                cbsri(1,nsr) = -1.0
                ubsri(1,nsr) = ureac(nrc)
                n = 1

                ncio = 0

                do nci = 1,iesrti(nsr)
                    if (uesri(nci,nsr)(1:2) .eq. 'O ') then
                        ncio = nci
                        go to 120
                    end if
                end do

120 continue

                ncih = 0

                do nci = 1,iesrti(nsr)
                    if (uesri(nci,nsr)(1:2) .eq. 'H ') then
                        ncih = nci
                        go to 130
                    end if
                end do

130 continue

                ctxo = 0.
                ctxh = 0.

                if (ncio .gt. 0) then
                    ctxo = -cesri(ncio,nsr)
                end if

                if (ncih .gt. 0) then
                    ctxh = -cesri(ncih,nsr)
                end if

                ztx = 0.

                ! Map all chemical elements except H and O.
                do nci = 1,iesrti(nsr)
                    if (nci.ne.ncio .and. nci.ne.ncih) then
                        do nc = 1,nct
                            if (uesri(nci,nsr)(1:8) .eq. uelem(nc)(1:8)) then
                                n = n + 1
                                ns = nc
                                nr1 = nessr(1,ns)
                                nr2 = nessr(2,ns)

                                cxe = 0.
                                cxo = 0.
                                cxh = 0.

                                do nn = nr1,nr2
                                    if (ness(nn) .eq. nco) then
                                        cxo = cess(nn)
                                    else if (ness(nn) .eq. nch) then
                                        cxh = cess(nn)
                                    else
                                        cxe = cess(nn)
                                    end if
                                end do

                                if (cxe .eq. 0.) then
                                    if (cxh.ne.0. .and. cxo.eq.0.) then
                                        cxe = cxh
                                        cxh = 0.
                                    else if (cxo.ne.0. .and. cxh.eq.0.) then
                                        cxe = cxo
                                        cxo = 0.
                                    else if (cxo.ne.0. .and. cxh.ne.0.) then
                                        cxe = cxo
                                        cxo = 0.
                                    end if
                                end if

                                cx = cesri(nci,nsr)/cxe
                                cbsri(n,nsr) = cx
                                ubsri(n,nsr) = uspec(ns)(1:24)
                                ctxo = ctxo + cxo*cx
                                ctxh = ctxh + cxh*cx
                                ztx = ztx + zchar(ns)*cx
                                go to 140
                            end if
                        end do

                        j2 = ilnobl(uesri(nci,nsr))
                        j3 = ilnobl(ureac(nrc))
                        write (noutpt,1000) uesri(nci,nsr)(1:j2),ureac(nrc)(1:j3)
                        write (nttyo,1000) uesri(nci,nsr)(1:j2),ureac(nrc)(1:j3)
1000 format(/" * Error - (EQ6/makrsr) Can't map the",/7x,'chemical element ',a,', which appears in the'          /7x,'composition for the special reactant ',a,',',/7x,'into a corresponding basis species to use in',/7x,'composing reaction for that reactant.')

                        nerr = nerr + 1
140 continue
                    end if
                end do

                ! Get H+ from charge balance.
                if (abs(ztx) .gt. eps100) then
                    n = n + 1
                    cx = -ztx
                    cbsri(n,nsr) = cx
                    ubsri(n,nsr) = 'H+ '
                    ztx = ztx + cx
                    ctxh = ctxh + cx
                end if

                ! Get H2O from H balance.
                if (abs(ctxh) .gt. eps100) then
                    n = n + 1
                    cx = -0.5*ctxh
                    cbsri(n,nsr) = cx
                    ubsri(n,nsr) = 'H2O '
                    ctxo = ctxo + cx
                    ctxh = ctxh + 2.*cx
                end if

                ! Get O2(g) from O balance.
                if (abs(ctxo) .gt. eps100) then
                    n = n + 1
                    cx = -0.5*ctxo
                    cbsri(n,nsr) = cx
                    ubsri(n,nsr) = 'O2(g)'
                    ctxo = ctxo + 2.*cx
                end if

                ! Check the balances on H, O, and charge.
                if (abs(ctxh).gt.eps100 .or. abs(ctxo).gt.eps100 .or.        abs(ztx).gt.eps100) then
                    j2 = ilnobl(ureac(nrc))
                    write (noutpt,1110) ureac(nrc)(1:j2)
                    write (nttyo,1110) ureac(nrc)(1:j2)
1110 format(/' * Error - (EQ6/makrsr) The reaction which',/7x,'was composed for the special reactant "',a,'"',/7x,'fails to satisfy all balance conditions.')

                    nerr = nerr + 1
                end if

                ibsrti(nsr) = n
            end if
        end if
    end do

    if (nerr .gt. 0) then
        stop
    end if
end subroutine makrsr