subroutine cmpdat(adhfs,adhfsd,advfs,advfsd,amu,amua,apx,apxa,aslm,aslma,atwt,atwta,axhfs,axhfsd,axlks,axlksd,axvfs,axvfsd,azero,azeroa,bpx,bpxa,cdrs,cdrsd,cess,cessa,iapxmx,iapxt,iapxta,iaqsla,iaqsln,ibpxmx,ibpxt,ibpxta,ilrn1,ilrn2,imrn1,imrn2,ipbtmx,ipch,ipchmx,ipcv,ipcvmx,insgf,insgfa,iopg,ixrn1,ixrn1a,ixrn2,ixrn2a,jflag,jpfcmx,jpflag,jsflag,jsitex,jsol,jsola,mwtsp,mwtspa,nalpaa,nalpha,napmax,napt,napta,narn1,narn1a,narn2,narn2a,narxmx,nat,natmax,nbasp,nbaspd,nbmap,nbt,nbtd,nbti,nbtmax,nchlor,ncmap,ncmpr,ncmpra,nct,ncta,nctmax,ndecsp,ndrs,ndrsd,ndrsmx,ndrsr,ndrsrd,ness,nessa,nessmx,nessr,nessra,ngrn1,ngrn1a,ngrn2,ngrn2a,ngt,nlrn1,nlrn1a,nlrn2,nlrn2a,nlt,nmrn1,nmrn1a,nmrn2,nmrn2a,nmt,nmut,nmuta,nmutmx,nmux,nmuxa,nopgmx,nslt,nsltmx,nslta,nslx,nslxa,nphasx,npt,npta,nptmax,nsmap,nsta,nst,nstmax,ntf1,ntf1a,ntf1mx,ntf1t,ntf1ta,ntf2,ntf2a,ntf2mx,ntf2t,ntf2ta,ntprmx,nxrn1,nxrn1a,nxrn2,nxrn2a,nxt,nxtmax,palpaa,palpha,qchlor,tf1,tf1a,tf2,tf2a,uelem,uelema,uphasa,uphase,uptypa,uptype,uspec,uspeca,vosp0,vosp0a,zchar,zchara)
    !! This sburoutine compresses the data read from the data file. The
    !! objective is to write data arrays (and associated variables)
    !! which exclude data belonging to phases and species that are not
    !! required by the current problem. Such phases are marked by
    !! the condition jpflag = 2. If a phase is not needed, none of its
    !! constituent species are needed. If a phase is needed, any
    !! constituent species that are not needed are marked by the
    !! condition jsflag = 2.
    !! The data read from the data file are stored in arrays (and a
    !! few important variables) which are characterized by an 'a'
    !! suffix, but which otherwise have the same names as the working
    !! arrays (and associated variables) which are created by this
    !! compression. Thus, the number of species read from the data
    !! file is nsta. The number remaining after compression is nst.
    !! The array of reaction coefficients read from the data file is
    !! cdrsd. The array created by the compression is cdrs.
    !! During compression, the basis set is compressed from nbtd elements
    !! in the nbaspd array to nbt elements in the nbasp array. Here nbasp
    !! and nbt are basically used as workspace. At the end of this
    !! subroutine, the nbaspd is set equal to the compressed nbasp array
    !! and nbtd is set equal to nbt.
    !! This subroutine is called by:
    !!   EQ3NR/eq3nr.f
    !!   EQ6/eq6.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: iapxmx
    integer :: ibpxmx
    integer :: ipbtmx
    integer :: ipchmx
    integer :: ipcvmx
    integer :: jpfcmx
    integer :: napmax
    integer :: narxmx
    integer :: natmax
    integer :: nbtmax
    integer :: nctmax
    integer :: ndrsmx
    integer :: nessmx
    integer :: nmutmx
    integer :: nopgmx
    integer :: nptmax
    integer :: nsltmx
    integer :: nstmax
    integer :: ntf1mx
    integer :: ntf2mx
    integer :: ntprmx
    integer :: nxtmax

    integer :: iapxt(nxtmax)
    integer :: iapxta(nxtmax)
    integer :: ibpxt(nxtmax)
    integer :: ibpxta(nxtmax)
    integer :: insgf(natmax)
    integer :: insgfa(natmax)
    integer :: iopg(nopgmx)
    integer :: jflag(nstmax)
    integer :: jpflag(nptmax)
    integer :: jsflag(nstmax)
    integer :: jsitex(nstmax)
    integer :: jsol(nxtmax)
    integer :: jsola(nxtmax)
    integer :: nalpaa(nsltmx)
    integer :: nalpha(nsltmx)
    integer :: nbasp(nbtmax)
    integer :: nbaspd(nbtmax)
    integer :: nbmap(nbtmax)
    integer :: ncmap(nctmax)
    integer :: ncmpr(2,nptmax)
    integer :: ncmpra(2,nptmax)
    integer :: ndecsp(nbtmax)
    integer :: ndrs(ndrsmx)
    integer :: ndrsd(ndrsmx)
    integer :: ndrsr(2,nstmax)
    integer :: ndrsrd(2,nstmax)
    integer :: ness(nessmx)
    integer :: nessa(nessmx)
    integer :: nessr(2,nstmax)
    integer :: nessra(2,nstmax)
    integer :: nmux(3,nmutmx)
    integer :: nmuxa(3,nmutmx)
    integer :: nphasx(nstmax)
    integer :: nslx(2,nsltmx)
    integer :: nslxa(2,nsltmx)
    integer :: nsmap(nstmax)
    integer :: ntf1(ntf1mx)
    integer :: ntf1a(ntf1mx)
    integer :: ntf2(ntf2mx)
    integer :: ntf2a(ntf2mx)

    integer :: iaqsla
    integer :: iaqsln
    integer :: ilrn1
    integer :: ilrn2
    integer :: imrn1
    integer :: imrn2
    integer :: ipch
    integer :: ipcv
    integer :: ixrn1
    integer :: ixrn1a
    integer :: ixrn2
    integer :: ixrn2a
    integer :: napt
    integer :: napta
    integer :: narn1
    integer :: narn1a
    integer :: narn2
    integer :: narn2a
    integer :: nat
    integer :: nbt
    integer :: nbtd
    integer :: nbti
    integer :: nchlor
    integer :: nct
    integer :: ncta
    integer :: ngrn1
    integer :: ngrn1a
    integer :: ngrn2
    integer :: ngrn2a
    integer :: ngt
    integer :: nlrn1
    integer :: nlrn1a
    integer :: nlrn2
    integer :: nlrn2a
    integer :: nlt
    integer :: nmrn1
    integer :: nmrn1a
    integer :: nmrn2
    integer :: nmrn2a
    integer :: nmt
    integer :: nmut
    integer :: nmuta
    integer :: nslt
    integer :: nslta
    integer :: npt
    integer :: npta
    integer :: nsta
    integer :: nst
    integer :: ntf1t
    integer :: ntf1ta
    integer :: ntf2t
    integer :: ntf2ta
    integer :: nxrn1
    integer :: nxrn1a
    integer :: nxrn2
    integer :: nxrn2a
    integer :: nxt

    logical :: qchlor

    character(len=8) :: uelem(nctmax)
    character(len=8) :: uelema(nctmax)
    character(len=48) :: uspec(nstmax)
    character(len=48) :: uspeca(nstmax)
    character(len=24) :: uphasa(nptmax)
    character(len=24) :: uphase(nptmax)
    character(len=24) :: uptypa(nptmax)
    character(len=24) :: uptype(nptmax)

    real(kind=8) :: adhfs(narxmx,ntprmx,ipchmx,nstmax)
    real(kind=8) :: adhfsd(narxmx,ntprmx,ipchmx,nstmax)
    real(kind=8) :: advfs(narxmx,ntprmx,ipcvmx,nstmax)
    real(kind=8) :: advfsd(narxmx,ntprmx,ipcvmx,nstmax)
    real(kind=8) :: apx(iapxmx,nxtmax)
    real(kind=8) :: apxa(iapxmx,nxtmax)
    real(kind=8) :: atwt(nctmax)
    real(kind=8) :: atwta(nctmax)
    real(kind=8) :: axhfs(narxmx,ntprmx,nstmax)
    real(kind=8) :: axhfsd(narxmx,ntprmx,nstmax)
    real(kind=8) :: axlks(narxmx,ntprmx,nstmax)
    real(kind=8) :: axlksd(narxmx,ntprmx,nstmax)
    real(kind=8) :: axvfs(narxmx,ntprmx,nstmax)
    real(kind=8) :: axvfsd(narxmx,ntprmx,nstmax)
    real(kind=8) :: azero(natmax)
    real(kind=8) :: azeroa(natmax)

    real(kind=8) :: amu(jpfcmx,nmutmx)
    real(kind=8) :: amua(jpfcmx,nmutmx)
    real(kind=8) :: aslm(jpfcmx,0:ipbtmx,nsltmx)
    real(kind=8) :: aslma(jpfcmx,0:ipbtmx,nsltmx)
    real(kind=8) :: bpx(ibpxmx,nxtmax)
    real(kind=8) :: bpxa(ibpxmx,nxtmax)
    real(kind=8) :: cdrs(ndrsmx)
    real(kind=8) :: cdrsd(ndrsmx)
    real(kind=8) :: cess(nessmx)
    real(kind=8) :: cessa(nessmx)
    real(kind=8) :: mwtsp(nstmax)
    real(kind=8) :: mwtspa(nstmax)
    real(kind=8) :: palpaa(ipbtmx,napmax)
    real(kind=8) :: palpha(ipbtmx,napmax)
    real(kind=8) :: tf1(ntf1mx)
    real(kind=8) :: tf1a(ntf1mx)
    real(kind=8) :: tf2(ntf2mx)
    real(kind=8) :: tf2a(ntf2mx)
    real(kind=8) :: vosp0(nstmax)
    real(kind=8) :: vosp0a(nstmax)
    real(kind=8) :: zchar(nstmax)
    real(kind=8) :: zchara(nstmax)

    ! Local variable declarations.
    integer :: i
    integer :: ipc
    integer :: j
    integer :: n
    integer :: na
    integer :: naa
    integer :: nap
    integer :: napa
    integer :: nb
    integer :: nba
    integer :: nbi
    integer :: nc
    integer :: nca
    integer :: nd
    integer :: ne
    integer :: nmu
    integer :: nmua
    integer :: nqt
    integer :: nsa
    integer :: np
    integer :: npa
    integer :: nrn1a
    integer :: nrn2a
    integer :: nr1
    integer :: nr1a
    integer :: nr2
    integer :: nr2a
    integer :: ns
    integer :: nsl
    integer :: nsla
    integer :: ns1
    integer :: ns1a
    integer :: ns2
    integer :: ns2a
    integer :: ns3
    integer :: ns3a
    integer :: nx
    integer :: nxa

    ! As long as the aqueous phase itself is not eliminated by
    ! compression, keep the chloride ion from being eliminated
    ! even if it is not present in the system to be modeled.
    ! This allows the interpretation of NBS pH, which requires
    ! the activity coefficient of the chloride ion.
    qchlor = .false.

    if (jpflag(iaqsla) .le. 0) then
        if (jsflag(nchlor) .eq. 2) then
            jsflag(nchlor) = 0
            qchlor = .true.
        end if
    end if

    ! The nsmap array is a pointer array. Here nsmap(nsa) = ns,
    ! where nsa is the index of a species before compression and ns
    ! is that of the same species after compression. If the species
    ! is eliminated by the compression, nsmap(nsa) = 0.
    do nsa = 1,nsta
        nsmap(nsa) = 0
    end do

    np = 0
    ns = 0
    ne = 0
    nd = 0

    do npa = 1,npta
        if (jpflag(npa) .lt. 2) then
            np = np + 1
            jpflag(np) = jpflag(npa)
            uphase(np) = uphasa(npa)
            uptype(np) = uptypa(npa)
            nrn1a = ncmpra(1,npa)
            nrn2a = ncmpra(2,npa)

            ! Loop over all the species belonging to this phase.
            ! At least one species should have jsflag .lt. 2.
            ncmpr(1,np) = ns + 1

            do nsa = nrn1a,nrn2a
                if (jsflag(nsa) .lt. 2) then
                    ns = ns + 1
                    nsmap(nsa) = ns
                    jsflag(ns) = jsflag(nsa)
                    uspec(ns) = uspeca(nsa)
                    mwtsp(ns) = mwtspa(nsa)
                    zchar(ns) = zchara(nsa)
                    vosp0(ns) = vosp0a(nsa)
                    jsitex(ns) = 1
                    nphasx(ns) = np

                    nr1a = nessra(1,nsa)
                    nr2a = nessra(2,nsa)
                    nessr(1,ns) = ne + 1

                    do n = nr1a,nr2a
                        ne = ne + 1
                        ness(ne) = nessa(n)
                        cess(ne) = cessa(n)
                    end do

                    nessr(2,ns) = ne

                    nr1a = ndrsrd(1,nsa)
                    nr2a = ndrsrd(2,nsa)
                    ndrsr(1,ns) = nd + 1

                    do n = nr1a,nr2a
                        nd = nd + 1
                        ndrs(nd) = ndrsd(n)
                        cdrs(nd) = cdrsd(n)
                    end do

                    ndrsr(2,ns) = nd

                    do j = 1,ntprmx
                        do i = 1,narxmx
                            axlks(i,j,ns) = axlksd(i,j,nsa)
                        end do
                    end do

                    if (ipch .ge. 0) then
                        do j = 1,ntprmx
                            do i = 1,narxmx
                                axhfs(i,j,ns) = axhfsd(i,j,nsa)
                            end do
                        end do

                        do ipc = 1,ipch
                            do j = 1,ntprmx
                                do i = 1,narxmx
                                    adhfs(i,j,ipc,ns) = adhfsd(i,j,ipc,nsa)
                                end do
                            end do
                        end do
                    end if

                    if (ipcv .ge. 0) then
                        do j = 1,ntprmx
                            do i = 1,narxmx
                                axvfs(i,j,ns) = axvfsd(i,j,nsa)
                            end do
                        end do

                        do ipc = 1,ipcv
                            do j = 1,ntprmx
                                do i = 1,narxmx
                                    advfs(i,j,ipc,ns) = advfsd(i,j,ipc,nsa)
                                end do
                            end do
                        end do
                    end if
                end if
            end do

            ncmpr(2,np) = ns
        end if
    end do

    npt = np
    nst = ns

    ! Suppress the aqueous chloride ion if it is not present.
    if (qchlor) then
        ns = nsmap(nchlor)
        jsflag(ns) = 2
    end if

    ! Update the jflag array.
    do nsa = 1,nsta
        ns = nsmap(nsa)

        if (ns .gt. 0) then
            jflag(ns) = jflag(nsa)
        end if
    end do

    do ns = nst + 1,nsta
        jflag(ns) = 0
    end do

    ! Update the species indices stored in the ndrs array.
    do n = 1,nd
        nsa = ndrs(n)

        if (nsa .gt. 0.) then
            ndrs(n) = nsmap(nsa)
        end if
    end do

    ! Update species range pointers.
    ! Calling sequence substitutions:
    !   narn1 for nlim1
    !   narn2 for nlim2
    !   narn1a for nlim1a
    !   narn2a for nlim2a
    !   nat for ntot
    call cnvndx(narn1,narn1a,narn2,narn2a,nsmap,nstmax,nat)

    ! Calling sequence substitutions:
    !   nmrn1 for nlim1
    !   nmrn2 for nlim2
    !   nmrn1a for nlim1a
    !   nmrn2a for nlim2a
    !   nmt for ntot
    call cnvndx(nmrn1,nmrn1a,nmrn2,nmrn2a,nsmap,nstmax,nmt)

    ! Calling sequence substitutions:
    !   nlrn1 for nlim1
    !   nlrn2 for nlim2
    !   nlrn1a for nlim1a
    !   nlrn2a for nlim2a
    !   nlt for ntot
    call cnvndx(nlrn1,nlrn1a,nlrn2,nlrn2a,nsmap,nstmax,nlt)

    ! Calling sequence substitutions:
    !   ngrn1 for nlim1
    !   ngrn2 for nlim2
    !   ngrn1a for nlim1a
    !   ngrn2a for nlim2a
    !   ngt for ntot
    call cnvndx(ngrn1,ngrn1a,ngrn2,ngrn2a,nsmap,nstmax,ngt)

    ! Note- nxt is the number of solid solution phases, not the
    ! number of solid solution species. Hence the use of nqt below
    ! for the number of species in solid solution phases. This
    ! variable (nqt) is not used elsewhere in this code.
    ! Calling sequence substitutions:
    !   nxrn1 for nlim1
    !   nxrn2 for nlim2
    !   nxrn1a for nlim1a
    !   nxrn2a for nlim2a
    !   nqt for ntot
    call cnvndx(nxrn1,nxrn1a,nxrn2,nxrn2a,nsmap,nstmax,nqt)

    ! Get the phase index of the aqueous solution.
    iaqsln = nphasx(narn1)

    ! Update phase range pointers for pure liquids, pure solids,
    ! and solid solutions.
    if (nlrn1 .le. 0) then
        ilrn1 = 0
        ilrn2 = -1
    else
        ilrn1 = nphasx(nlrn1)
        ilrn2 = nphasx(nlrn2)
    end if

    if (nmrn1 .le. 0) then
        imrn1 = 0
        imrn2 = -1
    else
        imrn1 = nphasx(nmrn1)
        imrn2 = nphasx(nmrn2)
    end if

    if (nxrn1 .le. 0) then
        nxt = 0
        ixrn1 = 0
        ixrn2 = -1
    else
        ixrn1 = nphasx(nxrn1)
        ixrn2 = nphasx(nxrn2)
        nxt = ixrn2 - ixrn1 + 1
    end if

    ! Compress the basis set. Use the nbmap array as a pointer array
    ! analogous to the nsmap array.
    nb = 0

    do nba = 1,nbtd
        nbmap(nba) = 0
        nsa = nbaspd(nba)
        ns = nsmap(nsa)

        if (ns .gt. 0) then
            nb = nb + 1
            nbmap(nba) = nb
            nbasp(nb) = ns
        end if
    end do

    nbt = nb

    ! Reset the ndecsp basis pointer array.
    do nbi = 1,nbti
        nba = ndecsp(nbi)

        if (nba .gt. 0) then
            ndecsp(nbi) = nbmap(nba)
        end if
    end do

    ! Compress the set of chemical elements. Use the ncmap
    ! array as a pointer array analogous to the nsmap array.
    nc = 0

    do nca = 1,ncta
        ncmap(nca) = 0

        do nba = 1,nbtd
            nsa = nbaspd(nba)
            ns = nsmap(nsa)

            if (ns .gt. 0) then
                nr1 = nessra(1,nsa)
                nr2 = nessra(2,nsa)

                do n = nr1,nr2
                    if (nessa(n) .eq. nca) then
                        nc = nc + 1
                        ncmap(nca) = nc
                        uelem(nc) = uelema(nca)
                        atwt(nc) = atwta(nca)
                        go to 460
                    end if
                end do
            end if
        end do

460 continue
    end do

    nct = nc

    ! Update the element indices stored in the ness array.
    do n = 1,ne
        nca = ness(n)

        if (nca .gt. 0) then
            ness(n) = ncmap(nca)
        end if
    end do

    ! Redefine the nbaspd array so that it matches the compressed
    ! nbasp array.
    do nb = 1,nbtd
        nbaspd(nb) = 0
    end do

    do nb = 1,nbt
        nbaspd(nb) = nbasp(nb)
    end do

    nbtd = nbt

    !      Create compressed activity coefficient data for solution
    !      phases other than the aqueous solution.
    ! XX   Note- the scheme below is not entirely satisfactory, as it
    ! XX   does not account for the effects of compression of component
    ! XX   species. The problem is that the apx-apxa arrays do not have
    ! XX   the data explcitly tied to the corresponding pairs or larger
    ! XX   groups of species. Revisit this issue when solid solutions
    ! XX   are redone.
    nxa = 0
    nx = 0

    do npa = ixrn1a,ixrn2a
        nxa = npa - ixrn1a + 1

        if (jpflag(npa) .lt. 2) then
            nx = nx + 1
            jsol(nx) = jsola(nxa)

            iapxt(nx) = iapxta(nxa)

            do i = 1,iapxta(nx)
                apx(i,nx) = apxa(i,nxa)
            end do

            ibpxt(nx) = ibpxta(nxa)

            do i = 1,ibpxta(nx)
                bpx(i,nx) = bpxa(i,nxa)
            end do
        end if
    end do

    ! Create compressed arrays for titration factors.
    ntf1t = 0

    do n = 1,ntf1ta
        nsa = ntf1a(n)
        ns = nsmap(nsa)

        if (ns .gt. 0) then
            ntf1t = ntf1t + 1
            ntf1(ntf1t) = ns
            tf1(ntf1t) = tf1a(n)
        end if
    end do

    ntf2t = 0

    do n = 1,ntf2ta
        nsa = ntf2a(n)
        ns = nsmap(nsa)

        if (ns .gt. 0) then
            ntf2t = ntf2t + 1
            ntf2(ntf2t) = ns
            tf2(ntf2t) = tf2a(n)
        end if
    end do

    ! Create compressed arrays for hard core diameters, neutral
    ! solute species flags, and hydration numbers.
    na = 0

    do nsa = narn1a,narn2a
        ns = nsmap(nsa)

        if (ns .gt. 0) then
            naa = nsa - narn1a + 1
            na = na + 1
            azero(na) = azeroa(naa)
            insgf(na) = insgfa(naa)
        end if
    end do

    if (iopg(1) .eq. 1) then
        ! Create compressed arrays for the parameters belonging to
        ! Pitzer's equations.
        nsl = 0

        do nsla = 1,nslta
            ns1a = nslxa(1,nsla)
            ns2a = nslxa(2,nsla)
            ns1 = nsmap(ns1a)
            ns2 = nsmap(ns2a)
            n = ns1*ns2

            if (n .gt. 0) then
                nsl = nsl + 1

                nalpha(nsl) = nalpaa(nsla)

                nslx(1,nsl) = ns1
                nslx(2,nsl) = ns2

                do i = 0,ipbtmx
                    do j = 1,jpfcmx
                        aslm(j,i,nsl) = aslma(j,i,nsla)
                    end do
                end do
            end if
        end do

        nslt = nsl

        ! Copy, but don't compress, the Pitzer alpha parameters.
        do napa = 1,napta
            nap = napa

            do i = 1,ipbtmx
                palpha(i,nap) = palpaa(i,napa)
            end do
        end do

        napt = napta

        nmu = 0

        do nmua = 1,nmuta
            ns1a = nmuxa(1,nmua)
            ns2a = nmuxa(2,nmua)
            ns3a = nmuxa(3,nmua)
            ns1 = nsmap(ns1a)
            ns2 = nsmap(ns2a)
            ns3 = nsmap(ns3a)
            n = ns1*ns2*ns3

            if (n .gt. 0) then
                nmu = nmu + 1
                nmux(1,nmu) = ns1
                nmux(2,nmu) = ns2
                nmux(3,nmu) = ns3

                do j = 1,jpfcmx
                    amu(j,nmu) = amua(j,nmua)
                end do
            end if
        end do

        nmut = nmu
    end if
end subroutine cmpdat