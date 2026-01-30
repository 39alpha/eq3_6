subroutine indats(adhfsa,advfsa,axhfsa,axlksa,axvfsa,cdrsa,cdrsv,cessa,cessv,ipch,ipch_asv,ipcv,ipcv_asv,mwtspa,nad1,narxt,narx_asv,nata,nata_asv,nbta,nbta_asv,nbta1_asv,nbtafd,ncmpra,ncta,ncta_asv,ndrsa,ndrsa_asv,ndrsn,ndrsra,nerr,nessa,nessa_asv,nessn,nessra,ngta,ngta_asv,nlta,nlta_asv,nmta,nmta_asv,noutpt,np,npta_asv,ns,nsta_asv,ntprt,ntpr_asv,nttyo,uaqsln,ubasp,udrsv,uelema,uendit,uessv,uphasa,uphasv,uptgas,uptliq,uptsld,uptypa,usblkf,uspeca,vosp0a,zchara)
    !! This subroutine reads a superblock of species blocks from the
    !! supporting data file "data1".
    !! This subroutine is called by:
    !!   EQLIB/indata.f
    !! Principal input:
    !!   nad1   = unit number of the data file
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: ipch_asv
    integer :: ipcv_asv
    integer :: narx_asv
    integer :: nata_asv
    integer :: nbta_asv
    integer :: nbta1_asv
    integer :: ncta_asv
    integer :: ndrsa_asv
    integer :: nessa_asv
    integer :: ngta_asv
    integer :: nlta_asv
    integer :: nmta_asv
    integer :: npta_asv
    integer :: nsta_asv
    integer :: ntpr_asv

    integer :: noutpt
    integer :: nttyo

    integer :: narxt(ntpr_asv)
    integer :: ncmpra(2,npta_asv)
    integer :: ndrsa(ndrsa_asv)
    integer :: ndrsra(2,nsta_asv)
    integer :: nessa(nessa_asv)
    integer :: nessra(2,nsta_asv)

    integer :: ipch
    integer :: ipcv
    integer :: nad1
    integer :: nata
    integer :: nbta
    integer :: nbtafd
    integer :: ncta
    integer :: ndrsn
    integer :: nerr
    integer :: nessn
    integer :: ngta
    integer :: nlta
    integer :: nmta
    integer :: np
    integer :: ns
    integer :: ntprt

    character(len=48) :: ubasp(nbta_asv)
    character(len=48) :: uspeca(nsta_asv)
    character(len=24) :: udrsv(nbta1_asv)
    character(len=24) :: uphasa(npta_asv)
    character(len=24) :: uptypa(npta_asv)
    character(len=24) :: uphasv
    character(len=24) :: uaqsln
    character(len=24) :: uptsld
    character(len=24) :: uptliq
    character(len=24) :: uptgas
    character(len=24) :: usblkf
    character(len=8) :: uelema(ncta_asv)
    character(len=8) :: uessv(ncta_asv)
    character(len=8) :: uendit

    real(kind=8) :: adhfsa(narx_asv,ntpr_asv,ipch_asv,nsta_asv)
    real(kind=8) :: advfsa(narx_asv,ntpr_asv,ipcv_asv,nsta_asv)
    real(kind=8) :: axhfsa(narx_asv,ntpr_asv,nsta_asv)
    real(kind=8) :: axlksa(narx_asv,ntpr_asv,nsta_asv)
    real(kind=8) :: axvfsa(narx_asv,ntpr_asv,nsta_asv)
    real(kind=8) :: cdrsa(ndrsa_asv)
    real(kind=8) :: cdrsv(nbta1_asv)
    real(kind=8) :: cessa(nessa_asv)
    real(kind=8) :: cessv(ncta_asv)
    real(kind=8) :: mwtspa(nsta_asv)
    real(kind=8) :: vosp0a(nsta_asv)
    real(kind=8) :: zchara(nsta_asv)

    ! Local variable declarations.
    integer :: ipc
    integer :: j2
    integer :: j3
    integer :: j4
    integer :: j5
    integer :: j6
    integer :: j7
    integer :: n
    integer :: nb
    integer :: nc
    integer :: nctsv
    integer :: ndrstv
    integer :: nf
    integer :: ntpr

    integer :: ilnobl

    character(len=24) :: uelect
    character(len=24) :: uhx
    character(len=24) :: uspecv
    character(len=24) :: ux24

    data uelect /'e-                      '/

    ! The label below defines a loop for reading the blocks in the
    ! superblock.
100 continue

    ! Read a species block.
    ! uspecv = species name
    ! nctsv  = number of chemical elements comprising the species
    ! ndrstv = number of reaction coefficients in the reaction
    !          associated with the species
    read (nad1) uspecv,nctsv,ndrstv
    j2 = ilnobl(uspecv)
    j5 = ilnobl(usblkf)

    ! Check for the end of the current superblock.
    if (uspecv(1:8) .eq. uendit(1:8)) then
        go to 999
    end if

    ! If the species is a pure mineral or pure liquid, set the phase
    ! name to match the species name.
    if (usblkf(1:24).eq.uptsld(1:24) .or. usblkf(1:24).eq.uptliq(1:24)) then
        uphasv(1:24) = uspecv(1:24)
    end if

    j3 = ilnobl(uphasv)

    ! Have another species.
    ns = ns + 1

    ! Load the species name array. Note that the last 24 characters
    ! contain the name of the associated phase.
    uspeca(ns)(1:24) = uspecv(1:24)
    uspeca(ns)(25:48) = uphasv(1:24)

    ! mwtspa = molecular weight
    ! zchara = electrical charge
    ! vosp0a = molar volume (not read).
    read (nad1) mwtspa(ns),zchara(ns)

    vosp0a(ns) = 0.

    ! cessv  = chemical element stoichiometric coefficient
    ! uessv  = name of corresponding chemical element.
    if (nctsv .gt. 0) then
        read (nad1) (cessv(n),uessv(n), n = 1,nctsv)
    end if

    ! Decode the elemental composition.
    nessra(1,ns) = nessn + 1

    if (nctsv .gt. 0) then
        do n = 1,nctsv
            do nc = 1,ncta
                if (uessv(n)(1:8) .eq. uelema(nc)(1:8)) then
                    nessn = nessn + 1
                    cessa(nessn) = cessv(n)
                    nessa(nessn) = nc
                    go to 120
                end if
            end do

            j4 = ilnobl(uessv(n))
            write (noutpt,1010) uspecv(1:j2),uphasv(1:j3),usblkf(1:j5),uessv(n)(1:j4)
            write (nttyo,1010) uspecv(1:j2),uphasv(1:j3),usblkf(1:j5),uessv(n)(1:j4)
1010 format(/' * Error - (EQLIB/indats) The species',/7x,a,' (',a,')',/7x,'appearing on the data file in the ',a,' superblock',/7x,'is listed as being composed of an unrecognized',/7x,'chemical element ',a,'.')

            nerr = nerr + 1
120 continue
        end do

        nessra(2,ns) = nessn
    else
        ! Have a species comprised of no chemical elements. This is
        ! permitted only for the following:
        !   O2(g) in a phase other than "gas"
        !   e-    in any phase
        if (uspecv(1:5) .eq. 'O2(g)') then
            if (usblkf(1:24) .eq. uptgas(1:24)) then
                write (noutpt,1020) uspecv(1:j2),uphasv(1:j3),usblkf(1:j5)
                write (nttyo,1020) uspecv(1:j2),uphasv(1:j3),usblkf(1:j5)
1020 format(/' * Error - (EQLIB/indats) The species',/7x,a,' (',a,')',/7x,'appearing on the data file in the ',a,' superblock',/7x,'is comprised of no chemical elements.')

                nerr = nerr + 1
            end if
        else if (uspecv(1:2) .eq. uelect(1:2)) then
            continue
        else
            write (noutpt,1020) uspecv(1:j2),uphasv(1:j3),usblkf(1:j5)
            write (nttyo,1020) uspecv(1:j2),uphasv(1:j3),usblkf(1:j5)
            nerr = nerr + 1
        end if

        ! A species with no elements is treated as being comprised of
        ! one "null" element.
        nessn = nessn + 1
        cessa(nessn) = 0.
        nessa(nessn) = 0
        nessra(2,ns) = nessn
    end if

    ! Read the associated reaction data, if any.
    ! Read the reaction, if any.
    if (ndrstv .gt. 0) then
        read (nad1) (cdrsv(n),udrsv(n), n = 1,ndrstv)
    end if

    ! Decode the associated reaction into a temporary format. The
    ! species associated with the reaction is moved if necessary so
    ! that it is the first species to appear in the reaction. Its
    ! species index is stored in the ndrsa array. The other species
    ! are presumed to be basis species. If one of these species has
    ! not been already loaded into the ubasp array, it is loaded at
    ! this time. The basis number of the species is stored in the
    ! ndrsa array. It will be converted to the corresponding species
    ! number after all species are loaded. This scheme permits the
    ! loading of a reaction which is written in terms of a basis
    ! species which itself has not yet been loaded.
    ndrsra(1,ns) = ndrsn + 1

    if (ndrstv .le. 0) then
        ! A species with no reaction listed on the data file is treated
        ! as a basis species. An associated reaction with one "null"
        ! species is created. The actual reaction is formally treated
        ! as one in which one mole of the associated species is destroyed
        ! and one mole of the same species is created. This is termed
        ! an "identity reaction." However, the reaction stored is not
        ! the identity reaction but the one with the "null" species as
        ! the only species.
        nbta = nbta + 1
        ubasp(nbta) = uspeca(ns)
        ndrsn = ndrsn + 1
        cdrsa(ndrsn) = 0.
        ndrsa(ndrsn) = 0
        ndrsra(2,ns) = ndrsn

        do ntpr = 1,ntpr_asv
            do n = 1,narx_asv
                axlksa(n,ntpr,ns) = 0.
            end do
        end do
    else if (ndrstv .eq. 1) then
        write (noutpt,1030) uspecv(1:j2),uphasv(1:j3),usblkf(1:j5)
        write (nttyo,1030) uspecv(1:j2),uphasv(1:j3),usblkf(1:j5)
1030 format(/' * Error - (EQLIB/indats) The species'  /7x,a,' (',a,')',/7x,'appearing on the data file in the ',a,' superblock',/7x,'has a reaction with only one species in it.')

        nerr = nerr + 1
    else
        ! Make sure that the species associated with the reaction is
        ! the first species in the reaction.
        do n = 1,ndrstv
            if (udrsv(n)(1:24) .eq. uspecv(1:24)) then
                ndrsn = ndrsn + 1
                cdrsa(ndrsn) = cdrsv(n)
                ndrsa(ndrsn) = ns
                go to 160
            end if
        end do

        write (noutpt,1040) uspecv(1:j2),uphasv(1:j3),usblkf(1:j5)
        write (nttyo,1040) uspecv(1:j2),uphasv(1:j3),usblkf(1:j5)
1040 format(/' * Error - (EQLIB/indats) The species',/7x,a,' (',a,')',/7x,'appearing on the data file in the ',a,' superblock',/7x,'does not appear in its associated reaction.')

        j4 = ilnobl(udrsv(1))
        write (noutpt,1050) udrsv(1)(1:j4),uphasv(1:j3),uspecv(1:j2),uphasv(1:j3)
        write (nttyo,1050) udrsv(1)(1:j4),uphasv(1:j3),uspecv(1:j2),uphasv(1:j3)
1050 format(/9x,'Is ',a,' (',a,') a doppelganger of',/12x,a,' (',a,')?')

        udrsv(1) = uspecv
        nerr = nerr + 1

160 continue
        nf = 0

        do n = 1,ndrstv
            if (udrsv(n)(1:24).eq.uspecv(1:24) .and. nf.eq.0) then
                nf = 1

                if (ns .le. nbtafd) then
                    ! The species is in the range of aqueous species that are
                    ! intended to be basis species. Add it to the basis set.
                    nbta = nbta + 1
                    ubasp(nbta) = udrsv(n)
                end if
            else
                ndrsn = ndrsn + 1
                cdrsa(ndrsn) = cdrsv(n)

                do nb = 1,nbta
                    if (ubasp(nb)(1:24) .eq. udrsv(n)(1:24)) then
                        ndrsa(ndrsn) = nb
                        go to 180
                    end if
                end do

                ! The reaction is written in terms of another species that
                ! is not yet in the basis set. Add it to that set.
                nbta = nbta + 1
                ubasp(nbta) = udrsv(n)
                ndrsa(ndrsn) = nbta
            end if

180 continue
        end do

        ndrsra(2,ns) = ndrsn

        ! Read the interpolating polynomial coefficients for the log K
        ! of the reaction.
        read (nad1) ux24
        j6 = ilnobl(ux24)
        uhx = 'Log K'
        j7 = ilnobl(uhx)

        if (uhx(1:j7) .ne. ux24(1:j6)) then
            write (noutpt,1052) ux24(1:j6),uhx(1:j7)
            write (nttyo,1052) ux24(1:j6),uhx(1:j7)
1052 format(/' * Error - (EQLIB/indats) Have read the tag',' string "',a,'"',/7x,'where "',a,'" was expected.')
        end if

        do ntpr = 1,ntprt
            read (nad1) (axlksa(n,ntpr,ns), n = 1,narxt(ntpr))
        end do
    end if

    if (ipch .ge. 0) then
        ! Read the interpolating polynomial coefficients for the enthalpy
        ! function and its pressure derivatives. For a basis species, the
        ! enthalpy function is the standard partial molar enthalpy of
        ! formation; for all other species, the enthalpy function is
        ! the standard partial molar enthalpy of reaction.
        read (nad1) ux24
        j6 = ilnobl(ux24)
        uhx = 'Enthalpy'
        j7 = ilnobl(uhx)

        if (uhx(1:j7) .ne. ux24(1:j6)) then
            write (noutpt,1052) ux24(1:j6),uhx(1:j7)
            write (nttyo,1052) ux24(1:j6),uhx(1:j7)
        end if

        do ntpr = 1,ntprt
            read (nad1) (axhfsa(n,ntpr,ns), n = 1,narxt(ntpr))
        end do

        do ipc = 1,ipch
            read (nad1) ux24
            j6 = ilnobl(ux24)

            if (ipc .eq. 1) then
                uhx = 'dEnthalpy/dP'
            else
                uhx = 'd( )Enthalpy/dP( )'
                write (uhx(3:3),'(i1)') ipc,ipc
                write (uhx(17:17),'(i1)') ipc,ipc
            end if

            j7 = ilnobl(uhx)

            if (uhx(1:j7) .ne. ux24(1:j6)) then
                write (noutpt,1052) ux24(1:j6),uhx(1:j7)
                write (nttyo,1052) ux24(1:j6),uhx(1:j7)
            end if

            do ntpr = 1,ntprt
                read (nad1) (adhfsa(n,ntpr,ipc,ns), n = 1,narxt(ntpr))
            end do
        end do
    end if

    if (ipcv .ge. 0) then
        ! Read the interpolating polynomial coefficients for the volume
        ! function and its pressure derivatives. For a basis species, the
        ! volume function is the standard partial molar enthalpy of
        ! formation; for all other species, the volume function is
        ! the standard partial molar volume of reaction.
        read (nad1) ux24
        j6 = ilnobl(ux24)
        uhx = 'Volume'
        j7 = ilnobl(uhx)

        if (uhx(1:j7) .ne. ux24(1:j6)) then
            write (noutpt,1052) ux24(1:j6),uhx(1:j7)
            write (nttyo,1052) ux24(1:j6),uhx(1:j7)
        end if

        do ntpr = 1,ntprt
            read (nad1) (axvfsa(n,ntpr,ns), n = 1,narxt(ntpr))
        end do

        do ipc = 1,ipcv
            read (nad1) ux24
            j6 = ilnobl(ux24)

            if (ipc .eq. 1) then
                uhx = 'dVolume/dP'
            else
                uhx = 'd( )Volume/dP( )'
                write (uhx(3:3),'(i1)') ipc,ipc
                write (uhx(15:15),'(i1)') ipc,ipc
            end if

            j7 = ilnobl(uhx)

            if (uhx(1:j7) .ne. ux24(1:j6)) then
                write (noutpt,1052) ux24(1:j6),uhx(1:j7)
                write (nttyo,1052) ux24(1:j6),uhx(1:j7)
            end if

            do ntpr = 1,ntprt
                read (nad1) (advfsa(n,ntpr,ipc,ns), n = 1,narxt(ntpr))
            end do
        end do
    end if

    if (usblkf(1:24) .eq. uaqsln(1:24)) then
        nata = nata + 1
    end if

    if (usblkf(1:24) .eq. uptsld(1:24)) then
        nmta = nmta + 1
        np = np + 1
        uphasa(np) = uspecv
        uptypa(np) = uptsld
        ncmpra(1,np) = ns
        ncmpra(2,np) = ns
    end if

    if (usblkf(1:24) .eq. uptliq(1:24)) then
        nlta = nlta + 1
        np = np + 1
        uphasa(np) = uspecv
        uptypa(np) = uptliq
        ncmpra(1,np) = ns
        ncmpra(2,np) = ns
    end if

    if (usblkf(1:24) .eq. uptgas(1:24)) then
        ngta = ngta + 1
    end if

    ! Loop back to read another species block.
    go to 100

999 continue
end subroutine indats
