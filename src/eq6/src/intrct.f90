subroutine intrct(cbsr,cbsri,cesr,cesri,egers,egersi,ibsrti,iern1,iern2,iesrti,ietmax,igerti,iktmax,imrn1,imrn2,ixrn1,ixrn2,ixrti,jcode,jetmax,jgerti,jgext,narn1,narn2,nbaspd,nbt,nbtmax,nbt1mx,ncmpr,nct,nctmax,nertmx,netmax,ngexsa,ngext,ngrn1,ngrn2,noutpt,nptmax,nrct,nrctmx,nrndex,nsrtmx,nstmax,nttyo,nxridx,nxrtmx,rxbar,rxbari,ubsri,ucxri,uelem,uesri,ugerji,ugermo,ugersi,ugexj,ugexmo,uphase,ureac,uspec,xgers,xgersi)
    !! This subroutine assigns phase or species indices corresponding to
    !! reactants.
    !! This subroutine is called by:
    !!   EQ6/eq6.f
    !! Principal input:
    !! Principal output:
    !!   jcode  = flag which determines the type of reactant:
    !!              0 =  Mineral
    !!              1 =  Solid solution
    !!              2 =  Special reactant
    !!              3 =  Aqueous species
    !!              4 =  Gas
    !!              5 =  Generic ion exchanger
    implicit none

    ! Calling sequence variable declarations.
    integer :: ietmax
    integer :: iktmax
    integer :: jetmax
    integer :: nbtmax
    integer :: nbt1mx
    integer :: nctmax
    integer :: nertmx
    integer :: netmax
    integer :: nptmax
    integer :: nrctmx
    integer :: nsrtmx
    integer :: nstmax
    integer :: nxrtmx

    integer :: ibsrti(nsrtmx)
    integer :: iesrti(nsrtmx)
    integer :: igerti(jetmax,nertmx)
    integer :: ixrti(nxrtmx)
    integer :: jcode(nrctmx)
    integer :: jgerti(nertmx)
    integer :: jgext(netmax)
    integer :: ncmpr(2,nptmax)
    integer :: nbaspd(nbtmax)
    integer :: ngexsa(ietmax,jetmax,netmax)
    integer :: ngext(jetmax,netmax)
    integer :: nrndex(nrctmx)
    integer :: nxridx(nrctmx)

    integer :: iern1
    integer :: iern2
    integer :: imrn1
    integer :: imrn2
    integer :: ixrn1
    integer :: ixrn2
    integer :: narn1
    integer :: narn2
    integer :: nbt
    integer :: nct
    integer :: ngrn1
    integer :: ngrn2
    integer :: noutpt
    integer :: nrct
    integer :: nttyo

    character(len=48) :: uspec(nstmax)
    character(len=24) :: ubsri(nbt1mx,nsrtmx)
    character(len=24) :: ucxri(iktmax,nxrtmx)
    character(len=24) :: ugermo(nertmx)
    character(len=24) :: ugersi(ietmax,jetmax,nertmx)
    character(len=24) :: ugexmo(netmax)
    character(len=24) :: uphase(nptmax)
    character(len=24) :: ureac(nrctmx)
    character(len=8) :: uelem(nctmax)
    character(len=8) :: uesri(nctmax,nsrtmx)
    character(len=8) :: ugerji(jetmax,nertmx)
    character(len=8) :: ugexj(jetmax,netmax)

    real(kind=8) :: cesr(nctmax,nsrtmx)
    real(kind=8) :: cesri(nctmax,nsrtmx)
    real(kind=8) :: cbsr(nbt1mx,nsrtmx)
    real(kind=8) :: cbsri(nbt1mx,nsrtmx)
    real(kind=8) :: egers(ietmax,jetmax,netmax)
    real(kind=8) :: egersi(ietmax,jetmax,nertmx)
    real(kind=8) :: rxbar(iktmax,nxrtmx)
    real(kind=8) :: rxbari(iktmax,nxrtmx)
    real(kind=8) :: xgers(ietmax,jetmax,netmax)
    real(kind=8) :: xgersi(ietmax,jetmax,nertmx)

    ! Local variable declarations.
    integer :: ie
    integer :: iei
    integer :: ik
    integer :: iki
    integer :: je
    integer :: jei
    integer :: j2
    integer :: j3
    integer :: j4
    integer :: j5
    integer :: j6
    integer :: nb
    integer :: nbi
    integer :: nc
    integer :: nci
    integer :: ne
    integer :: ner
    integer :: nerr
    integer :: np
    integer :: nrc
    integer :: nr1
    integer :: nr2
    integer :: ns
    integer :: nsr
    integer :: nss
    integer :: nt
    integer :: nxr

    integer :: ilnobl

    character(len=24) :: urn
    character(len=24) :: ux
    character(len=24) :: uaq
    character(len=24) :: umn
    character(len=24) :: ugs
    character(len=24) :: uss
    character(len=24) :: uge
    character(len=16) :: ux16
    character(len=8) :: uxd
    character(len=8) :: ux8a
    character(len=8) :: ux8b

    real(kind=8) :: ex
    real(kind=8) :: exc
    real(kind=8) :: exo
    real(kind=8) :: rx
    real(kind=8) :: rxo
    real(kind=8) :: rxsum

    data uaq    /'aqueous species         '/,umn    /'minerals                '/,ugs    /'gases                   '/,uss    /'solid solutions         '/,uge    /'generic ion exchangers  '/

    if (nrct .eq. 0) then
        go to 999
    end if

    nerr = 0
    nsr = 0
    nxr = 0
    ner = 0

    do nrc = 1,nrct
        urn = ureac(nrc)

        if (jcode(nrc) .eq. 0) then
            ! Pure minerals.
            do np = imrn1,imrn2
                if (urn(1:24) .eq. uphase(np)(1:24)) then
                    nrndex(nrc) = np
                    go to 200
                end if
            end do

            j2 = ilnobl(urn)
            j3 = ilnobl(umn)
            write (noutpt,1000) urn(1:j2),umn(1:j3)
            write (nttyo,1000) urn(1:j2),umn(1:j3)
1000 format(/' * Error - (EQ6/intrct) The reactant ',a," isn't among",/7x,'the loaded ',a,'.')

            nerr = nerr + 1
        else if (jcode(nrc) .eq. 1) then
            ! Solid solutions.
            nxr = nxr + 1

            do np = ixrn1,ixrn2
                if (urn(1:24) .eq. uphase(np)(1:24)) then
                    nrndex(nrc) = np
                    nxridx(nrc) = nxr
                    nr1 = ncmpr(1,np)
                    nr2 = ncmpr(2,np)
                    nt = nr2 - nr1 + 1

                    ! Make sure that the mole fractions sum to unity.
                    rxsum = 0.

                    do iki = 1,ixrti(nxr)
                        rxsum = rxsum + rxbari(iki,nxr)
                    end do

                    if (abs(rxsum - 1.0) .gt. 1.e-3) then
                        j3 = ilnobl(urn)
                        write (ux16,'(g11.4)') rxsum
                        call lejust(ux16)
                        j4 = ilnobl(ux16)
                        write (noutpt,1002) urn(1:j3),ux16(1:j4)
                        write (nttyo,1002) urn(1:j3),ux16(1:j4)
1002 format(/' * Warning - (EQ6/intrct) The mole fractions',' given for the composition',/7x,'of the solid',' solution reactant ',a,' sum to',/7x,a,', not unity.',' They will be renormalized.')

                        ! Renormalize the mole fractions.
                        write (noutpt,1004)
                        write (nttyo,1004)
1004 format(/5x,'Renormalizing the Input Mole Fractions',//9x,'Component',19x,'Old x',9x,'New x',/)

                        do iki = 1,ixrti(nxr)
                            rxo = rxbari(iki,nxr)
                            rx = rxo/rxsum
                            write (noutpt,1006) ucxri(iki,nxr),rxo,rx
                            write (nttyo,1006) ucxri(iki,nxr),rxo,rx
1006 format(7x,a24,3x,1pe11.4,3x,e11.4)

                            rxbari(iki,nxr) = rx
                        end do
                    end if

                    ! Decode the solid solution composition.
                    do iki = 1,ixrti(nxr)
                        ux = ucxri(iki,nxr)
                        ns = nr1 - 1

                        do ik = 1,nt
                            ns = ns + 1

                            if (uspec(ns)(1:24) .eq. ux(1:24)) then
                                rxbar(ik,nxr) = rxbari(iki,nxr)
                                go to 100
                            end if
                        end do

                        j2 = ilnobl(ux)
                        j3 = ilnobl(urn)
                        write (noutpt,1010) ux(1:j2),urn(1:j3)
                        write (nttyo,1010) ux(1:j2),urn(1:j3)
1010 format(/' * Error - (EQ6/intrct) The species ',a," doesn't appear",/7x,'on the data file as a',' component of ',a,'.')

                        nerr = nerr + 1

100 continue
                    end do

                    go to 200
                end if
            end do

            j2 = ilnobl(urn)
            j3 = ilnobl(uss)
            write (noutpt,1000) urn(1:j2),uss(1:j3)
            write (nttyo,1000) urn(1:j2),uss(1:j3)
            nerr = nerr + 1
        else if (jcode(nrc) .eq. 2) then
            ! Special reactants.
            nsr = nsr + 1
            nrndex(nrc) = nsr
            nxridx(nrc) = nsr

            ! Decode the chemical composition.
            do nci = 1,iesrti(nsr)
                uxd = uesri(nci,nsr)

                do nc = 1,nct
                    if (uelem(nc)(1:8) .eq. uxd(1:8)) then
                        cesr(nc,nsr) = cesri(nci,nsr)
                        go to 110
                    end if
                end do

                j2 = ilnobl(uxd)
                j3 = ilnobl(urn)
                write (noutpt,1020) uxd(1:j2),urn(1:j3)
                write (nttyo,1020) uxd(1:j2),urn(1:j3)
1020 format(/' * Error - (EQ6/intrct) The element ',a,' in special reactant',/7x,a," isn't on the data file.")

                nerr = nerr + 1

110 continue
            end do

            ! Decode the chemical reaction.
            ux = ubsri(1,nsr)

            if (ux(1:24) .eq. ureac(nrc)(1:24)) then
                cbsr(nbt1mx,nsr) = cbsri(1,nsr)
            else
                j2 = ilnobl(ux)
                j3 = ilnobl(urn)
                write (noutpt,1030) ux(1:j2),urn(1:j3)
                write (nttyo,1030) ux(1:j2),urn(1:j3)
1030 format(' * Error - (EQ6/intrct) The species ',a,' appears',' first',/7x,'in the reaction associated with the special',' reactant',/7x,a,'. The special reactant itself must',' appear first.')

                nerr = nerr + 1
            end if

            do nbi = 2,ibsrti(nsr)
                ux = ubsri(nbi,nsr)

                do nb = 1,nbt
                    ns = nbaspd(nb)

                    if (uspec(ns)(1:24) .eq. ux(1:24)) then
                        cbsr(nb,nsr) = cbsri(nbi,nsr)
                        go to 120
                    end if
                end do

                j2 = ilnobl(ux)
                j3 = ilnobl(urn)
                write (noutpt,1040) ux(1:j2),urn(1:j3)
                write (nttyo,1040) ux(1:j2),urn(1:j3)
1040 format(' * Error - (EQ6/intrct) The species ',a,' which',' appears',/7x,'in the reaction associated with the',' special reactant',/7x,a," isn't in the data file',' basis set.")

                nerr = nerr + 1
120 continue
            end do
        else if (jcode(nrc) .eq. 3) then
            ! Aqueous species.
            do ns = narn1,narn2
                if (urn(1:24) .eq. uspec(ns)(1:24)) then
                    nrndex(nrc) = ns
                    go to 200
                end if
            end do

            j2 = ilnobl(urn)
            j3 = ilnobl(uaq)
            write (noutpt,1000) urn(1:j2),uaq(1:j3)
            write (nttyo,1000) urn(1:j2),uaq(1:j3)
            nerr = nerr + 1
        else if (jcode(nrc) .eq. 4) then
            ! Gases.
            do ns = ngrn1,ngrn2
                if (urn(1:24) .eq. uspec(ns)(1:24)) then
                    nrndex(nrc) = ns
                    go to 200
                end if
            end do

            j2 = ilnobl(urn)
            j3 = ilnobl(ugs)
            write (noutpt,1000) urn(1:j2),ugs(1:j3)
            write (nttyo,1000) urn(1:j2),ugs(1:j3)
            nerr = nerr + 1
        else if (jcode(nrc) .eq. 5) then
            ! Generic ion exchangers.
            ner = ner + 1

            do np = iern1,iern2
                if (urn(1:24) .eq. uphase(np)(1:24)) then
                    nxridx(nrc) = ner
                    nrndex(nrc) = np
                    ne = np - iern1 + 1

                    if (ugermo(ner)(1:24) .ne. ugexmo(ne)(1:24)) then
                        j2 = ilnobl(urn)
                        j5 = ilnobl(ugermo(ner))
                        j6 = ilnobl(ugexmo(ne))
                        write (noutpt,1042) urn(1:j2),ugermo(ner)(1:j5),ugexmo(ne)(1:j6)
                        write (nttyo,1042) urn(1:j2),ugermo(ner)(1:j5),ugexmo(ne)(1:j6)
1042 format(/' * Error - (EQ6/intrct) The generic exchanger',' phase reactant',/7x,a,' has the specified exchange',' model',/7x,'"',a,'", but the defined exchange model',' for this',/7x,'exchanger phase is "',a,'".')

                        nerr = nerr + 1
                    end if

                    if (jgerti(ner) .ne. jgext(ne)) then
                        j2 = ilnobl(urn)
                        write (ux8a,'(i5)') jgerti(ner)
                        call lejust(ux8a)
                        j5 = ilnobl(ux8a)
                        write (ux8b,'(i5)') jgext(ne)
                        call lejust(ux8b)
                        j6 = ilnobl(ux8b)
                        write (noutpt,1044) urn(1:j2),ux8a(1:j5),ux8b(1:j6)
                        write (nttyo,1044) urn(1:j2),ux8a(1:j5),ux8b(1:j6)
1044 format(/' * Error - (EQ6/intrct) The generic exchanger',' phase reactant',/7x,a,' has a specified composition',' including',/7x,a,' exchange site(s), but this phase',' actually has ',a,' exchange',/7x,'site(s). The',' number of such sites must match.')

                        nerr = nerr + 1
                    end if

                    do jei = 1,jgerti(ner)
                        ! Make sure that the exchange fractions on the current
                        ! site sum to unity.
                        exc = 0.

                        do iei = 1,igerti(jei,ner)
                            exc = exc + egersi(iei,jei,ner)
                        end do

                        if (abs(exc - 1.0) .gt. 1.e-3) then
                            j3 = ilnobl(urn)
                            j4 = ilnobl(ugerji(jei,ner))
                            write (ux16,'(g11.4)') exc
                            call lejust(ux16)
                            j5 = ilnobl(ux16)
                            write (noutpt,1050) ugerji(1,j5),urn(1:j3),ux16(1:j5)
                            write (nttyo,1050) ugerji(1,j5),urn(1:j3),ux16(1:j5)
1050 format(/' * Warning - (EQ6/intrct) The exchange',' fractions given for the composition',/7x,'of',' exchange site ',a,' of the generic ion exchanger',/7x,'reactant ',a,' sum to ',a,', not unity.',/7x,'They will be renormalized.')

                            ! Renormalize the exchange fractions.
                            write (noutpt,1060)
                            write (nttyo,1060)
1060 format(/5x,'Renormalizing the Input Exchange',' Fractions',//9x,'Component',19x,'Old e',9x,'New e',/)

                            do iei = 1,igerti(jei,ner)
                                exo = egersi(iei,jei,ner)
                                ex = exo/exc
                                write (noutpt,1006) ugersi(iei,jei,ner),exo,ex
                                write (nttyo,1006) ugersi(iei,jei,ner),exo,ex
                                egersi(iei,jei,ner) = ex
                            end do
                        end if

                        ! Decode the composition of the current exchange site.
                        do je = 1,jgext(ne)
                            if (ugexj(je,ne)(1:8) .eq. ugerji(jei,ner)(1:8)) then
                                go to 150
                            end if
                        end do

                        j4 = ilnobl(ugerji(jei,ner))
                        write (noutpt,1070) ugerji(jei,ner)(1:j4),urn(1:j3)
                        write (nttyo,1070) ugerji(jei,ner)(1:j4),urn(1:j3)
1070 format(/' * Error - (EQ6/intrct) The exchange site ',a,' specified for',/7x,'the generic ion exchanger',' reactant ',a,'does not',/7x,'match any site',' defined for that exchanger.')

                        nerr = nerr + 1
                        go to 170

150 continue
                        do iei = 1,igerti(jei,ner)
                            ux = ugersi(iei,jei,ner)

                            ! Match with the exchange species name. This contains
                            ! the site name.
                            ns = ncmpr(1,np)

                            do je = 1,jgext(ne)
                                do ie = 1,ngext(je,ne)
                                    ns = ns + 1

                                    if (uspec(ns)(1:24) .eq. ux(1:24)) then
                                        egers(ie,je,ne) = egersi(iei,jei,ner)
                                        go to 160
                                    end if
                                end do
                            end do

                            ! No match was found above. Match with the name of
                            ! the aqueous species which exchanges. This does
                            ! not contain the site name.
                            do je = 1,jgext(ne)
                                do ie = 1,ngext(je,ne)
                                    nss = ngexsa(ie,je,ne)

                                    if (nss .gt. 0) then
                                        if (uspec(nss)(1:24) .eq. ux(1:24)) then
                                            egers(ie,je,ne) = egersi(iei,jei,ner)
                                            go to 160
                                        end if
                                    end if
                                end do
                            end do

                            j2 = ilnobl(ux)
                            j3 = ilnobl(urn)
                            write (noutpt,1080) ux(1:j2),urn(1:j3)
                            write (nttyo,1080) ux(1:j2),urn(1:j3)
1080 format(/' * Error - (EQ6/intrct) The species ',a,' specified for',/7x,'the generic ion exchanger',' reactant ',a,' does not',/7x,'match any species',' defined for that exchanger.')

                            nerr = nerr + 1

160 continue
                        end do

170 continue
                    end do

                    go to 200
                end if
            end do

            j2 = ilnobl(urn)
            j3 = ilnobl(uge)
            write (noutpt,1000) urn(1:j2),uge(1:j3)
            write (nttyo,1000) urn(1:j2),uge(1:j3)
            nerr = nerr + 1
        end if

200 continue
    end do

    if (nerr .gt. 0) then
        stop
    end if

999 continue
end subroutine intrct