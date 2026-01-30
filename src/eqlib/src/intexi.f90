subroutine intexi(al10,axhfs,axlks,axvfs,cegexs,cess,cdrs,cgexj,cpgexs,egexjf,iern1,iern2,ietmax,jern1,jern2,jetmax,jflag,jgext,jpflag,jsflag,jsitex,kern1,kern2,ketmax,kgexsa,mwtges,mwtsp,narn1,narn2,narxmx,narxt,nbasp,nbt,nbtmax,ncmpr,ndrs,ndrsmx,ndrsr,ness,nessmx,nessr,netmax,net,nern1,nern2,ngexro,ngexrt,ngexsa,ngexso,ngext,noutpt,nphasx,npt,nptmax,nst,nstmax,ntprt,ntprmx,nttyo,nvetmx,rconst,tgexp,ugexj,ugexmo,ugexmv,ugexp,ugexr,ugexs,ugexsr,uhfgex,uphase,uspec,uvfgex,uxkgex,xhfgex,xlkgex,xvfgex,zchar,zgexj)
    !! This subroutine inteprets data read from the input file for the
    !! generic ion exchange model. It composes the necessary generic
    !! exchange phases and species and sets up the associated reactions
    !! and thermodynamic properties.
    !! "Ion exchange" may be strictly or loosely interpreted here,
    !! depending on the models chosen. In the case of the Gapon and
    !! Vanselow models, the interpretation is strict. For example, a
    !! site must be negatively or positively charged, only ions of
    !! the opposite sign may occupy it, and the formal exchange
    !! capacity is exactly 100% utilized, save for a negligible
    !! trace of unfilled positions. The end-member exchanger species
    !! are electrically neutral, save for the trace unfilled site
    !! species. In contrast, the "Site-mixing" model allows more
    !! flexibility, including charged exchanger end-members, and
    !! overloading or underloading of the formal exchange capacities.
    !! Each exchanger species is specific to a given site (type of
    !! site), and these species are organized in the species list such
    !! that those belonging to a given site are listed contiguously.
    !! An exchanger with only one site is treated in the same manner
    !! as a phase composed of constituent species. In the case of an
    !! exchanger with more than one site, each site is treated as a
    !! kind of subphase, composed of its own constituent species.
    !! To illustrate, an exchange phase with substrate Z and two exchange
    !! sites is represented by:
    !!   (Na+,K+)3(Fe+++,Al+++)2(Z)
    !! Here a pair of parentheses represents a site. Z represents the
    !! substrate. The site it occupies is not an exchange site. There is
    !! always one mole of substrate occupying this site. The first
    !! exchange site, here called site A, contains a mixture of Na+
    !! and K+ ions. There are 3 moles of this site per mole of substrate.
    !! This site has an intrinsic charge (charge when the site is not
    !! occupied) of -1. Thus, 3 moles of monovalent cations per mole of
    !! substrate are required to attain electrical neutrality on this
    !! site. The second exchange site, here called site B, contains a
    !! mixture of Fe+++ and Al+++ ions. There are 2 moles of this site
    !! per mole of substrate, and it has an intrinsic charge of -3.
    !! The species for this example may be written as:
    !!   site A:
    !!     (Na+)3()2()
    !!     (K+)3()2()
    !!   site B:
    !!     ()3(Fe+++)2()
    !!     ()3(Al+++)2()
    !!   substrate site:
    !!     ()3()2(Z)
    !! Here an empty pair of parentheses denotes an empty site.
    !! However, an empty site is a definite part of the species.
    !! An empty site makes no contribution to the molecular weight of
    !! of such a species. Thus, the molecular weight of (Na+)3()2()
    !! is 3 times the atomic weight of Na. The mass in grams of the
    !! complete exchanger phase is the sum of the masses of all of
    !! the constitutent species. The number of moles of the exchanger
    !! phase is equal to the number of moles of the substrate.
    !! One could conceivably define an exchanger with a non-unit
    !! stoichiometry for the substrate component; e.g., something
    !! like:
    !!   (Na+,K+)3(Fe+++,Al+++)2(Z')2
    !! where Z' is a redefined substrate for the previous example.
    !! However, there do not seem to be any advantages to introducing
    !! or even allowing such a complexity. One can always define a
    !! stoichiometric with unit substrate, in this example by taking
    !! Z = 2*Z'.
    !! The activity of such a species in the ideal case (implicitly in
    !! the site mixing sense) is the mole fraction on the site for which
    !! the species is defined, raised to the power equal to the
    !! stoichiometric factor for the site. In the case of (Na+)3()2(),
    !! the relation is:
    !!   a{(Na+)3()2()} = x{(Na+)3()2()}**3.
    !! The origin of this relationship is as follows. Write a pseudo-
    !! reaction in which the (Na+)3()2() species goes to a form of
    !! itself, but with the stoichiometry redefined so that there is
    !! one mole of the site of interest per mole of the exchanger
    !! species:
    !!   (Na+)3()2() = 3(Na+)()[2/3]()[1/3]
    !! Because this is a pseudo-reaction, reflecting only a change
    !! in components, there is no change in the Gibbs energy, enthalpy,
    !! or entropy. Alternatively, the equilibrium constant is unity.
    !! Thus,
    !!   a{(Na+)3()2()} = a{(Na+)()[2/3]()[1/3]}**3
    !! In the ideal case, it is clear that:
    !!   a{(Na+)()[2/3]()[1/3]} = x{(Na+)()[2/3]()[1/3]}
    !! Thus,
    !!   a{(Na+)3()2()} = x{(Na+)()[2/3]()[1/3]}**3
    !! The mole fractions of the Na+ species are not affected by the
    !! stoichiometric arbitrariness in how the species are defined.
    !! Hence:
    !!   x{(Na+)3()2()} = x{(Na+)()[2/3]()[1/3]}
    !! Subsitution of this into the equaiton immediately above then
    !! yields:
    !!   a{(Na+)3()2()} = x{(Na+)3()2()}**3.
    !! This explains the origin of the site number appearing as an
    !! exponent in the calculation of the activity of a component.
    !! The above kinds of species are not end-members. End-members for
    !! the above exchanger phase would be:
    !!   (Na+)3(Fe+++)2(Z)
    !!   (Na+)3(Al+++)2(Z)
    !!   (K+)3(Fe+++)2(Z)
    !!   (K+)3(Al+++)2(Z)
    !! Note that for example:
    !!   a((Na+)3(Fe+++)2(Z)) = a((Na+)3()2())**3 a(()3(Fe+++)2())**2
    !! Note that a factor of a(Z)**1 is missing from the right hand
    !! side. This is because this factor has a fixed value of 1 by
    !! definition. This would be true even if we allowed for the
    !! formal possibility of an arbitrary number of moles of substrate
    !! in the exchanger phase, because there is no mixing on the
    !! substrate site. Thus, a(Z) = 1.
    !! This subroutine is called by:
    !!   EQ3NR/eq3nr.f
    !!   EQ6/eq6.f
    !! Principal input:
    !!   ugexj  = array of exchanger site names
    !!   ugexmo = array of strings identifying exchange models, which
    !!              are used for composing ion exchange species and
    !!              corresponding reactions; examples include:
    !!                'Gapon'    = Gapon model
    !!                'Vanselow' = Vanselow model
    !!                'Site-mixing' = the general site-mixing model
    !!   ugexmv = array of valid strings identifying exchange models
    !!   ugexp  = array of exchanger phase names
    !!   ugexr  = array of strings containing compact representations
    !!              of the exchange reactions; e.g., 'Na+ = Ca++' for a
    !!              reaction in which Na+ on the exchanger is replaced
    !!              by Ca++. One may make specifications such as
    !!              'Na+ = ' in which case the ion goes into solution
    !!              leaving a bare substrate. All reactions are
    !!              normalized to the exchange (or loss) of one
    !!              equivalent. The exact form of the reaction is
    !!              otherwise dependent on the mixing law specifed in
    !!              the element of the ugexmo array for the current
    !!              exchanger.
    !!   ugexs  = array of exchange ions names
    !!   ugexsr = array of exchange ions names extracted from the
    !!              ugexr array
    !! Principal input/output:
    !!   uphase = array of phase names (extended to include the names
    !!              of generic exchange phases)
    !!   uspec  = array of species names (extended to include the names
    !!              of species belonging to generic exchange phases)
    !! Principal output:
    !!   cegexs = array of coefficients giving the number of equivalents
    !!              per mole of each exchanger species.
    !!   cpgexs = array of coefficients giving the number of moles of
    !!              exchanger substrate per mole of each exchanger
    !!              species.
    !!   egexjf = array of formal exchange capacities of the sites of
    !!              the exchangers, defined in terms of equivalents per
    !!              mole of exchanger substrate. The formal exchange
    !!              capacity of an exchanger phase is the sum of these
    !!              for all its exchange sites.
    !!   jern1  = array giving the start of the range in the species
    !!              list corresponding to species in the je-th site
    !!              of the ne-th exchanger.
    !!   jern2  = array giving the end of the range in the species
    !!              list corresponding to species in the je-th site
    !!              of the ne-th exchanger.
    !!   jgext  = array giving the number of exchange sites in each of
    !!              the ion exchange phases
    implicit none

    ! Calling sequence variable declarations.
    integer :: ietmax
    integer :: jetmax
    integer :: ketmax
    integer :: narxmx
    integer :: nbtmax
    integer :: ndrsmx
    integer :: nessmx
    integer :: netmax
    integer :: nptmax
    integer :: nstmax
    integer :: ntprmx
    integer :: nvetmx

    integer :: noutpt
    integer :: nttyo

    integer :: jern1(jetmax,netmax)
    integer :: jern2(jetmax,netmax)
    integer :: jflag(nstmax)
    integer :: jgext(netmax)
    integer :: jpflag(nptmax)
    integer :: jsflag(nstmax)
    integer :: jsitex(nstmax)
    integer :: kern1(netmax)
    integer :: kern2(netmax)
    integer :: kgexsa(ketmax,netmax)
    integer :: narxt(ntprmx)
    integer :: nbasp(nbtmax)
    integer :: ncmpr(2,nptmax)
    integer :: ndrs(ndrsmx)
    integer :: ndrsr(2,nstmax)
    integer :: ness(nessmx)
    integer :: nessr(2,nstmax)
    integer :: ngexro(ietmax,jetmax,netmax)
    integer :: ngexrt(jetmax,netmax)
    integer :: ngext(jetmax,netmax)
    integer :: ngexsa(ietmax,jetmax,netmax)
    integer :: ngexso(ietmax,jetmax,netmax)
    integer :: nphasx(nstmax)

    integer :: iern1
    integer :: iern2
    integer :: narn1
    integer :: narn2
    integer :: nbt
    integer :: net
    integer :: nern1
    integer :: nern2
    integer :: npt
    integer :: nst
    integer :: ntprt

    character(len=56) :: ugexr(ietmax,jetmax,netmax)
    character(len=48) :: uspec(nstmax)
    character(len=24) :: ugexmo(netmax)
    character(len=24) :: ugexmv(nvetmx)
    character(len=24) :: ugexp(netmax)
    character(len=24) :: ugexs(ietmax,jetmax,netmax)
    character(len=24) :: ugexsr(2,ietmax,jetmax,netmax)
    character(len=24) :: uphase(nptmax)
    character(len=8) :: ugexj(jetmax,netmax)
    character(len=8) :: uhfgex(ietmax,jetmax,netmax)
    character(len=8) :: uvfgex(ietmax,jetmax,netmax)
    character(len=8) :: uxkgex(ietmax,jetmax,netmax)

    real(kind=8) :: axhfs(narxmx,ntprmx,nstmax)
    real(kind=8) :: axlks(narxmx,ntprmx,nstmax)
    real(kind=8) :: axvfs(narxmx,ntprmx,nstmax)
    real(kind=8) :: cegexs(ietmax,jetmax,netmax)
    real(kind=8) :: cess(nessmx)
    real(kind=8) :: cdrs(ndrsmx)
    real(kind=8) :: cgexj(jetmax,netmax)
    real(kind=8) :: cpgexs(ietmax,jetmax,netmax)
    real(kind=8) :: egexjf(jetmax,netmax)
    real(kind=8) :: mwtges(netmax)
    real(kind=8) :: mwtsp(nstmax)
    real(kind=8) :: tgexp(netmax)
    real(kind=8) :: xhfgex(ietmax,jetmax,netmax)
    real(kind=8) :: xlkgex(ietmax,jetmax,netmax)
    real(kind=8) :: xvfgex(ietmax,jetmax,netmax)
    real(kind=8) :: zchar(nstmax)
    real(kind=8) :: zgexj(jetmax,netmax)

    real(kind=8) :: al10
    real(kind=8) :: rconst

    ! Local variable declarations.
    integer :: i
    integer :: ie
    integer :: iee
    integer :: iej
    integer :: ieo
    integer :: itot
    integer :: j
    integer :: je
    integer :: jee
    integer :: jj
    integer :: jj2
    integer :: jj3
    integer :: j2
    integer :: j3
    integer :: j4
    integer :: j5
    integer :: j6
    integer :: k
    integer :: ke
    integer :: kee
    integer :: ke1
    integer :: n
    integer :: nb
    integer :: nd2
    integer :: ne
    integer :: nee
    integer :: nerr
    integer :: nn
    integer :: np
    integer :: nrf1
    integer :: nrf2
    integer :: nr1
    integer :: nr2
    integer :: ns
    integer :: nsba
    integer :: nse
    integer :: nsl
    integer :: nspect
    integer :: nss
    integer :: nsse
    integer :: nve
    integer :: nvet

    integer :: ilnobl
    integer :: nbasis

    logical :: qmoerr
    logical :: qsperr

    character(len=56) :: ustr56
    character(len=24) :: ugexpd
    character(len=24) :: ustr1
    character(len=24) :: ustr1e
    character(len=24) :: ustr2
    character(len=8) :: ugexjd

    real(kind=8) :: afhx
    real(kind=8) :: afvx
    real(kind=8) :: afxb
    real(kind=8) :: afx0
    real(kind=8) :: arcnst
    real(kind=8) :: bfx
    real(kind=8) :: cfactr
    real(kind=8) :: cgx
    real(kind=8) :: cfxi
    real(kind=8) :: cfxz
    real(kind=8) :: cx
    real(kind=8) :: efx
    real(kind=8) :: tfxb
    real(kind=8) :: tfx0
    real(kind=8) :: xx
    real(kind=8) :: zprod
    real(kind=8) :: zpx
    real(kind=8) :: zsi
    real(kind=8) :: zxj
    real(kind=8) :: zxjt

    ! Initialize some constants.
    nerr = 0
    arcnst = -al10*0.001*rconst
    tfx0 = 273.15

    ! Count the number of valid exchange model options.
    nvet = 0

    do nve = 1,nvetmx
        if (ugexmv(nve)(1:6) .ne. 'ERROR ') then
            nvet = nvet + 1
        end if
    end do

    if (net.gt.0 .and. nvet.le.0) then
        write (noutpt,1000)
        write (nttyo,1000)
1000 format (/' * Error - (EQLIB/intexi) Programming error trap:',/7x,'There are no programmed valid strings representing',/7x,'exchange models for the generic ion phases. Check the',/7x,'data statement which initializes the ugexmv array. This',/7x,'is located in the EQLIB INCLUDE file eqlo8d.h. If the',/7x,'INCLUDE files have been "pre-stuffed" in your source',/7x,'code, this data statement may be found in the main',/7x,'programs for EQ3NR and EQ6.')

        stop
    end if

    ! Pass 1. Loop on exchanger phases. Check the input. Provide
    ! defaults as necessary. build a list of exchanger species, and
    ! map exchange ions to their aqueous species counterparts.
    ! The exchange phases and species will be created (in the
    ! sense of being added to the general lists of phases and species)
    ! below in pass 2.
    do ne = 1,net
        call lejust(ugexp(ne))
        j4 = ilnobl(ugexp(ne))

        if (j4 .le. 0) then
            ! Have a blank exchanger name. Assign a default name.
            call adgexp(ne,noutpt,nttyo,ugexpd)
            ugexp(ne) = ugexpd
            j4 = ilnobl(ugexp(ne))
        end if

        ! Make sure that the phase name is unique among all exchangers.
        do nee = 1,ne - 1
            if (ugexp(ne)(1:24) .eq. ugexp(nee)(1:24)) then
                write (noutpt,1010) ugexp(ne)(1:j4),ne,nee
                write (nttyo,1010) ugexp(ne)(1:j4),ne,nee
1010 format (/' * Error - (EQLIB/intexi) "',a,'" is the name',' specified',/7x,'for exchanger phases ',i3,' and ',i3,'. Exchanger phase names,',/7x,'including any assigned',' defaults, must be unique.')

                nerr = nerr + 1
            end if
        end do

        ! Check the prescribed exchange model.
        call lejust(ugexmo(ne))

        ! Provide a default.
        if (ugexmo(ne)(1:3) .eq. '   ') then
            ugexmo(ne) = ugexmv(1)
            write (noutpt,1012) ugexp(ne)(1:j4),ugexmo(ne)(1:j2)
            write (nttyo,1012) ugexp(ne)(1:j4),ugexmo(ne)(1:j2)
1012 format (/' * Warning - (EQLIB/intexi) No exchange model',' was species for',/7x,'the exchanger phase ',a,'. The ',a,' model has been',/7x,'assigned as a default.')
        end if

        j5 = ilnobl(ugexmo(ne))

        do nve = 1,nvetmx
            j6 = ilnobl(ugexmv(nve))

            if (ugexmo(ne)(1:j5) .eq. ugexmv(nve)(1:j6)) then
                go to 130
            end if
        end do

        write (noutpt,1020) ugexmo(ne)(1:j2),ugexp(ne)(1:j4)
        write (nttyo,1020) ugexmo(ne)(1:j2),ugexp(ne)(1:j4)
1020 format (/" * Error - (EQLIB/intexi) Don't recognize the",' exchange model type',/7x,'"',a,'" specified for the',' generic exchanger phase ',a,'.',/7x,'Valid choices',' include the following:',/)

        do nve = 1,nvetmx
            if (ugexmv(nve)(1:6) .ne. 'ERROR ') then
                j6 = ilnobl(ugexmv(nve))
                write (noutpt,1022) ugexmv(nve)(1:j6)
                write (nttyo,1022) ugexmv(nve)(1:j6)
1022 format(9x,a)
            end if
        end do

        j6 = ilnobl(ugexmv(1))
        write (noutpt,1024) ugexmv(1)(1:j6)
        write (nttyo,1024) ugexmv(1)(1:j6)
1024 format(/7x,'A blank input defaults to "',a,'".')

        nerr = nerr + 1

130 continue
        if (mwtges(ne) .le. 0) then
            mwtges(ne) = 100.
            write (noutpt,1030) ugexp(ne)(1:j4)
            write (nttyo,1030) ugexp(ne)(1:j4)
1030 format (/" * Warning - (EQLIB/intexi) Don't have a valid",' input for the',/7x,'molecular weight of the substrate for',' exchange phase'/7x,a,'. Assigning a default value of',' 100 grams/mol.')
        end if

        if (tgexp(ne) .eq. 0.) then
            tgexp(ne) = 25.0
            write (noutpt,1040) ugexp(ne)(1:j4)
            write (nttyo,1040) ugexp(ne)(1:j4)
1040 format (/" * Warning - (EQLIB/intexi) Don't have a valid",' input for the',/7x,'temperature reference required for',' the thermodynamic data for',/7x,'exchange phase ',a,'. Assigning a default value of 25C.')
        end if

        cfactr = 1./(arcnst*(tgexp(ne) + 273.15))

        ! Loop on sites.
        do je = 1,jgext(ne)
            call lejust(ugexj(je,ne))
            j3 = ilnobl(ugexj(je,ne))

            ! Calculate the formally declared exchange capacity (in
            ! equivalents per mole of exchanger substrate) of the
            ! current site. The sign of this is opposite to that of the
            ! charge on the site itself.
            egexjf(je,ne) = -cgexj(je,ne)*zgexj(je,ne)

            ! Check the name of the current site.
            if (j3 .eq. 0) then
                ! No name was given. Assign one (e.g., "S(1)" to site 1,
                ! "S(2)" to site 2).
                call adgexj(je,noutpt,nttyo,ugexjd)
                ugexj(je,ne) = ugexjd
                j3 = ilnobl(ugexj(je,ne))
            end if

            ! Make sure that the site name is unique among all sites for
            ! the current exchanger phase.
            do jee = 1,je - 1
                if (ugexj(je,ne)(1:8) .eq. ugexj(jee,nee)(1:8)) then
                    write (noutpt,1060) ugexj(je,ne)(1:j3),je,jee,ugexp(ne)(1:j4)
                    write (nttyo,1060) ugexj(je,ne)(1:j3),je,jee,ugexp(ne)(1:j4)
1060 format (/' * Error - (EQLIB/intexi) "',a,'" is the name',' of',/7x,'sites ',i3,' and ',i3,' of exchanger phase ',a,'.',/7x,'Site names, including any assigned',' defaults, must be unique',/7x,'for a given',' exchanger phase.')

                    nerr = nerr + 1
                end if
            end do

            ! Check the electrical charge of the site against the specified
            ! exchange model.
            if (ugexmo(ne)(1:5).eq.'Gapon' .or.      ugexmo(ne)(1:6).eq.'Gapon-' .or.      ugexmo(ne)(1:8).eq.'Vanselow' .or.      ugexmo(ne)(1:9).eq.'Vanselow-') then
                if (zgexj(je,ne) .eq. 0.) then
                    write (noutpt,1070) ugexj(je,ne)(1:j3),ugexp(ne)(1:j4),ugexmo(ne)(1:j5)
                    write (nttyo,1070) ugexj(je,ne)(1:j3),ugexp(ne)(1:j4),ugexmo(ne)(1:j5)
1070 format (/' * Error - (EQLIB/intexi) The electrical',' charge of site ',a,''/7x,'of exchanger phase ',a,' is zero. This is not valid',/7x,'for the ',a,' exchange model.')

                    nerr = nerr + 1
                end if
            end if

            ! Loop on condensed reactions read from the input file. Map
            ! these to species corresponding to ions on the sites.
            qsperr = .false.

            do n = 1,ngexrt(je,ne)
                k = index(ugexr(n,je,ne),' == ')

                if (k .gt. 0) then
                    ustr56 = ugexr(n,je,ne)(1:k - 1)
                    call lejust(ustr56)
                    ustr1 = ustr56(1:24)
                    ustr56 = ugexr(n,je,ne)(k + 4:56)
                    call lejust(ustr56)
                    ustr2 = ustr56(1:24)
                else
                    k = index(ugexr(n,je,ne),' = ')

                    if (k .gt. 0) then
                        ustr56 = ugexr(n,je,ne)(1:k - 1)
                        call lejust(ustr56)
                        ustr1 = ustr56(1:24)
                        ustr56 = ugexr(n,je,ne)(k + 3:56)
                        call lejust(ustr56)
                        ustr2 = ustr56(1:24)
                    else
                        k = index(ugexr(n,je,ne),'==')

                        if (k .gt. 0) then
                            ustr56 = ugexr(n,je,ne)(1:k - 1)
                            call lejust(ustr56)
                            ustr1 = ustr56(1:24)
                            ustr56 = ugexr(n,je,ne)(k + 2:56)
                            call lejust(ustr56)
                            ustr2 = ustr56(1:24)
                        else
                            k = index(ugexr(n,je,ne),'=')

                            if (k .gt. 0) then
                                ustr56 = ugexr(n,je,ne)(1:k - 1)
                                call lejust(ustr56)
                                ustr1 = ustr56(1:24)
                                ustr56 = ugexr(n,je,ne)(k + 1:56)
                                call lejust(ustr56)
                                ustr2 = ustr56(1:24)
                            else
                                j2 = ilnobl(ugexr(n,je,ne))
                                write (noutpt,1080) ugexr(n,je,ne)(1:j2),ugexj(je,ne)(1:j3),ugexp(ne)(1:j4)
                                write (nttyo,1080) ugexr(n,je,ne)(1:j2),ugexj(je,ne)(1:j3),ugexp(ne)(1:j4)
1080 format (/" * Error - (EQLIB/intexi) Can't decode",' the condensed reaction string',/7x,'"',a,'",',' which is given for site ',a,' of exchange',/7x,'phase ',a,'. The two species in the string',' must be',/7x,'separated by " == ", " = ", "==",',' or "=".')

                                ustr1 = 'Error'
                                ustr2 = 'Error'
                                nerr = nerr + 1
                                qsperr = .true.
                            end if
                        end if
                    end if
                end if

                ! Prettify the strings in the ugexr array. Use "__" to
                ! indicate a bare exchange site.
                if (ustr1(1:3) .eq. '   ') then
                    ustr1 = '__'
                end if

                if (ustr1(1:3) .eq. '_  ') then
                    ustr1 = '__'
                end if

                if (ustr2(1:3) .eq. '   ') then
                    ustr2 = '__'
                end if

                if (ustr2(1:3) .eq. '_  ') then
                    ustr2 = '__'
                end if

                ugexsr(1,n,je,ne) = ustr1
                ugexsr(2,n,je,ne) = ustr2
                j2 = ilnobl(ustr1)
                ugexr(n,je,ne)(j2 + 1:j2 + 3) = ' = '
                ugexr(n,je,ne)(j2 + 4:56) = ustr2

                ! Make sure the two species in the reaction are not identical.
                if (ustr1(1:6) .ne. 'Error ') then
                    if (ustr1(1:24) .eq. ustr2(1:24)) then
                        j2 = ilnobl(ugexr(n,je,ne))

                        if (ustr1(1:3) .eq. '__ ') then
                            write (noutpt,1090) ugexr(n,je,ne)(1:j2),ugexj(je,ne)(1:j3),ugexp(ne)(1:j4)
                            write (nttyo,1090) ugexr(n,je,ne)(1:j2),ugexj(je,ne)(1:j3),ugexp(ne)(1:j4)
1090 format (/' * Error - (EQLIB/intexi) The condensed',' exchange reaction "',a,'"',/7x,'for site ',a,' of exchange phase ',a,' is an',/7x,'identity'," reaction. Such a reaction shouldn't appear on the",/7x,'input file.')
                        else
                            write (noutpt,1100) ugexr(n,je,ne)(1:j2),ugexj(je,ne)(1:j3),ugexp(ne)(1:j4)
                            write (nttyo,1100) ugexr(n,je,ne)(1:j2),ugexj(je,ne)(1:j3),ugexp(ne)(1:j4)
1100 format (/' * Error - (EQLIB/intexi) The condensed',' exchange reaction "',a,'"',/7x,'for site ',a,' of exchange phase ',a,' is invalid.',/7x,'The two'," species in the reaction can't be identical.")
                        end if

                        nerr = nerr + 1
                        qsperr = .true.
                    end if
                end if
            end do

            if (qsperr) then
                go to 200
            end if

            ! Add a condensed reaction which involves one of the
            ! exchangeable species and the bare site species, if none
            ! such was included on the input.
            do n = 1,ngexrt(je,ne)
                ustr1 = ugexsr(1,n,je,ne)
                ustr2 = ugexsr(2,n,je,ne)

                if (ustr1(1:3) .eq. '__ ') then
                    go to 140
                end if

                if (ustr2(1:3) .eq. '__ ') then
                    go to 140
                end if
            end do

            n = ngexrt(je,ne) + 1
            ngexrt(je,ne) = n
            ustr1 = ugexsr(1,1,je,ne)
            ustr2 = '__'
            ugexsr(1,n,je,ne) = ustr1
            ugexsr(2,n,je,ne) = ustr2
            j2 = ilnobl(ustr1)
            ugexr(n,je,ne) = ustr1
            ugexr(n,je,ne)(j2 + 1:j2 + 3) = ' = '
            ugexr(n,je,ne)(j2 + 4:56) = ustr2
            xlkgex(n,je,ne) = -12.0
            xhfgex(n,je,ne) = 0.
            xvfgex(n,je,ne) = 0.
            uxkgex(n,je,ne) = 'LogK/eq'
            uhfgex(n,je,ne) = 'kcal/eq'
            uvfgex(n,je,ne) = 'cm3/eq'
140 continue

            ! Add an identity reaction for the bare site species.
            ! This is done to simplify the indexing relations between
            ! exchanger species and corresponding reactions by making
            ! them one to one. This simplification is important, because
            ! the species must be created in a special order so that the
            ! corresponding reaction properties for dissocation to the
            ! bare site species can be calculated from the specified
            ! exchange data without having to resort to matrix equations.
            n = ngexrt(je,ne) + 1
            ngexrt(je,ne) = n
            ustr1 = '__'
            ustr2 = '__'
            ugexsr(1,n,je,ne) = ustr1
            ugexsr(2,n,je,ne) = ustr2
            j2 = ilnobl(ustr1)
            ugexr(n,je,ne) = ustr1
            ugexr(n,je,ne)(j2 + 1:j2 + 3) = ' = '
            ugexr(n,je,ne)(j2 + 4:56) = ustr2
            xlkgex(n,je,ne) = 0.
            xhfgex(n,je,ne) = 0.
            xvfgex(n,je,ne) = 0.
            uxkgex(n,je,ne) = 'LogK/eq'
            uhfgex(n,je,ne) = 'kcal/eq'
            uvfgex(n,je,ne) = 'cm3/eq'

            ! Put the names of the ions on the current site
            ! in the ugexs array.
            nspect = 2
            ustr1 = ugexsr(1,1,je,ne)
            ustr2 = ugexsr(2,1,je,ne)
            ugexs(1,je,ne) = ustr1
            ugexs(2,je,ne) = ustr2

            ! The last reaction is the identity reaction for the bare site
            ! species, and it is guaranteed that a previous reaction
            ! involves this species, so looping from 2,ngexrt(je,ne) - 1
            ! is sufficient.
            do n = 2,ngexrt(je,ne) - 1
                do i = 1,2
                    ustr1 = ugexsr(i,n,je,ne)

                    do nn = 1,nspect
                        if (ustr1(1:24) .eq. ugexs(nn,je,ne)) then
                            go to 150
                        end if
                    end do

                    nspect = nspect + 1
                    ugexs(nspect,je,ne) = ustr1
150 continue
                end do
            end do

            if (nspect .ne. ngexrt(je,ne)) then
                write (noutpt,1110) nspect,ngexrt(je,ne),ugexj(je,ne)(1:j3),ugexp(ne)(1:j4)
                write (nttyo,1110) nspect,ngexrt(je,ne),ugexj(je,ne)(1:j3),ugexp(ne)(1:j4)
1110 format (/' * Error - (EQLIB/intexi) Have ',i3,' species',' and ',i3,' reactions for',/7x,'site ',a,' of exchange',' phase ',a,' after automatic',/7x,'additions to deal with',' the bare site species. The number',/7x,'of species',' (including the bare site species) must be equal to',/7x,'the number of reactions. Check the input. If you have',' listed for',/7x,'a site for example two exchange',' reactions such as "Na+ = K+" and',/7x,'"Rb+ = Cs+,"',' you must include a third one to complete the linkage,',/7x,'such as "Na+ = Rb+".')

                nerr = nerr + 1
            end if

            ! Map each exchange species to its aqueous species counterpart
            ! (the exchangeable species).
            do ie = 1,nspect
                if (ugexs(ie,je,ne)(1:3) .eq. '__ ') then
                    ngexsa(ie,je,ne) = 0
                    go to 170
                end if

                j2 = ilnobl(ugexs(ie,je,ne))

                do nss = narn1,narn2
                    if (ugexs(ie,je,ne)(1:24) .eq. uspec(nss)(1:24)) then
                        ngexsa(ie,je,ne) = nss
                        go to 170
                    end if
                end do

                write (noutpt,1130) ugexs(ie,je,ne)(1:j2),ugexj(je,ne)(1:j3),ugexp(ne)(1:j4)
                write (nttyo,1130) ugexs(ie,je,ne)(1:j2),ugexj(je,ne)(1:j3),ugexp(ne)(1:j4)
1130 format (/' * Error - (EQLIB/intexi) The exchange species ',a,/7x,'specified for site ',a,' of exchange phase ',a,/7x,'does not correspond to any aqueous species read from',' the data file.')

                nerr = nerr + 1
                qsperr = .true.

170 continue
            end do

            if (qsperr) then
                go to 200
            end if

            ! Test the electrical charges of the exchangeable species
            ! and the charge of the site against the exchanger model.
            if (ugexmo(ne)(1:5).eq.'Gapon' .or.      ugexmo(ne)(1:6).eq.'Gapon-' .or.      ugexmo(ne)(1:8).eq.'Vanselow' .or.      ugexmo(ne)(1:9).eq.'Vanselow-') then
                do ie = 1,nspect
                    nss = ngexsa(ie,je,ne)

                    if (nss .gt. 0) then
                        zpx = zchar(nss)*zgexj(je,ne)

                        if (zpx .ge. 0.) then
                            write (noutpt,1150) ugexs(ie,je,ne)(1:j2),zchar(nss),ugexj(je,ne)(1:j3),ugexp(ne)(1:j4),ugexmo(ne)(1:j5)
                            write (nttyo,1150) ugexs(ie,je,ne)(1:j2),zchar(nss),ugexj(je,ne)(1:j3),ugexp(ne)(1:j4),ugexmo(ne)(1:j5)
1150 format (/' * Error - (EQLIB/intexi) The electrical',' charge of ',a,/7x,'is ',f5.1,'. This species'," can't be exchanged onto site ",a,/7x,'of exchanger',' phase ',a,' because the site',/7x,'has the same',' charge sign or zero charge. Opposite charge',/7x,'signs are required for the ',a,' exchange',' model.')

                            nerr = nerr + 1
                        end if
                    end if
                end do
            end if

            ! Convert the input thermodynamic parameters for exchange or
            ! dissociation reactions to standard units.
            do ie = 1,nspect
                call lejust(uxkgex(ie,je,ne))

                if (uxkgex(ie,je,ne)(1:3) .eq. '   ') then
                    uxkgex(ie,je,ne) = 'LogK/eq'
                else if (uxkgex(ie,je,ne)(1:8) .eq. 'kcal/eq ') then
                    xlkgex(ie,je,ne) = cfactr*xlkgex(ie,je,ne)
                    uxkgex(ie,je,ne) = 'LogK/eq'
                else if (uxkgex(ie,je,ne)(1:8) .eq. 'kJ/eq   ') then
                    xlkgex(ie,je,ne) = xlkgex(ie,je,ne)/4.184
                    xlkgex(ie,je,ne) = cfactr*xlkgex(ie,je,ne)
                    uxkgex(ie,je,ne) = 'LogK/eq'
                else if (uxkgex(ie,je,ne)(1:8) .ne. 'LogK/eq ') then
                    j2 = ilnobl(ugexr(ie,je,ne))
                    j5 = ilnobl(uxkgex(ie,je,ne))
                    write (noutpt,1170) uxkgex(ie,je,ne)(1:j5),ugexr(ie,je,ne)(1:j2),ugexj(je,ne)(1:j3),ugexp(ne)(1:j4)
                    write (nttyo,1170) uxkgex(ie,je,ne)(1:j5),ugexr(ie,je,ne)(1:j2),ugexj(je,ne)(1:j3),ugexp(ne)(1:j4)
1170 format (/" * Error - (EQLIB/intexi) Don't recognize the",' string "',a,'", which is',/7x,'given to define the',' units for the xlkgex (log K) input for the',/7x,'reaction "',a,'" on site ',a,' of exchanger phase',/7x,a,'. The string must be blank (defaults to',' "LogK/eq"),',/7x,'"LogK/eq", "kcal/eq", or "kJ/eq".')

                    nerr = nerr + 1
                end if

                call lejust(uhfgex(ie,je,ne))

                if (uhfgex(ie,je,ne)(1:3) .eq. '   ') then
                    uhfgex(ie,je,ne) = 'kcal/eq'
                else if (uhfgex(ie,je,ne)(1:8) .eq. 'kJ/eq   ') then
                    xhfgex(ie,je,ne) = xhfgex(ie,je,ne)/4.184
                    uhfgex(ie,je,ne) = 'kcal/eq'
                else if (uhfgex(ie,je,ne)(1:8) .ne. 'kcal/eq ') then
                    j2 = ilnobl(ugexr(ie,je,ne))
                    j5 = ilnobl(uhfgex(ie,je,ne))
                    write (noutpt,1180) uhfgex(ie,je,ne)(1:j5),ugexr(ie,je,ne)(1:j2),ugexj(je,ne)(1:j3),ugexp(ne)(1:j4)
                    write (nttyo,1180) uhfgex(ie,je,ne)(1:j5),ugexr(ie,je,ne)(1:j2),ugexj(je,ne)(1:j3),ugexp(ne)(1:j4)
1180 format (/" * Error - (EQLIB/intexi) Don't recognize the",' string "',a,'", which is',/7x,'given to define the',' units for the xhfgex (enthalpy) input for the',/7x,'reaction "',a,'" on site ',a,' of exchanger phase',/7x,a,'. The string must be blank (defaults to',' "kcal/eq"),',/7x,'"kcal/eq", or "kJ/eq".')

                    nerr = nerr + 1
                end if

                call lejust(uvfgex(ie,je,ne))

                if (uvfgex(ie,je,ne)(1:3) .eq. '   ') then
                    uvfgex(ie,je,ne) = 'cm3/eq'
                else if (uvfgex(ie,je,ne)(1:8) .ne. 'cm3/eq  ') then
                    j2 = ilnobl(ugexr(ie,je,ne))
                    j5 = ilnobl(uvfgex(ie,je,ne))
                    write (noutpt,1190) uvfgex(ie,je,ne)(1:j5),ugexr(ie,je,ne)(1:j2),ugexj(je,ne)(1:j3),ugexp(ne)(1:j4)
                    write (nttyo,1190) uvfgex(ie,je,ne)(1:j5),ugexr(ie,je,ne)(1:j2),ugexj(je,ne)(1:j3),ugexp(ne)(1:j4)
1190 format (/" * Error - (EQLIB/intexi) Don't recognize the",' string "',a,'", which is',/7x,'given to define the',' units for the xvfgex (volume) input for the',/7x,'reaction "',a,' on site ',a,' of exchanger phase',/7x,a,'. The string must be blank (defaults to',' "cm3/eq")',/7x,'or "cm3/eq".')

                    nerr = nerr + 1
                end if
            end do

200 continue
        end do
    end do

    if (nerr .gt. 0) then
        stop
    end if

    np = npt
    ns = nst
    nern1 = nst + 1
    iern1 = npt + 1

    ! Pass 2. Loop on exchanger phases. Formally create exchanger phases
    ! and species, and comute the corresponding reactions and
    ! thermodynamic data for such species.
    do ne = 1,net
        j4 = ilnobl(ugexp(ne))

        if (npt .ge. nptmax) then
            write (noutpt,1240) nptmax,ugexp(ne)(1:j2)
            write (nttyo,1240) nptmax,ugexp(ne)(1:j2)
1240 format (/' * Error - (EQLIB/intexi) The maximum ',i3,' phases would be',/7x,'exceeded by creating the exchange',' phase ',a,'. Increase',/7x,'the dimensioning parameter',' nptpar.')

            stop
        end if

        ! Add the current exchange phase to the list of phases.
        np = np + 1
        npt = np
        jpflag(np) = 0
        uphase(np) = ugexp(ne)
        ncmpr(1,np) = ns + 1
        tfxb = tgexp(ne)

        ! Loop on sites.
        nsba = 0

        do je = 1,jgext(ne)
            j3 = ilnobl(ugexj(je,ne))

            nspect = ngexrt(je,ne)

            ! Change the condensed reactions if necessary so that any
            ! association reactions are converted to dissociation
            ! reactions (e.g., for Na+, from "__ = Na+" to "Na+ = __").
            ! Do not modify the last reaction, which is the identity
            ! reaction for the bare site species.
            do ie = 1,nspect - 1
                if (ugexsr(1,ie,je,ne)(1:3) .eq. '__ ') then
                    xlkgex(ie,je,ne) = -xlkgex(ie,je,ne)
                    xhfgex(ie,je,ne) = -xhfgex(ie,je,ne)
                    xvfgex(ie,je,ne) = -xvfgex(ie,je,ne)
                    ugexsr(1,ie,je,ne)(1:24) = ugexsr(2,ie,je,ne)(1:24)
                    ugexsr(2,ie,je,ne)(1:24) = '__ '
                end if
            end do

            ! Now convert all exchange reactions in this set to
            ! dissociation reactions. To give an idea of what is happening,
            ! a single exchange reaction "Ca++ = Na+" read from the input
            ! file would now be expanded to the following condensed
            ! reaction set:
            !   1. Ca++ = Na+
            !   2. Ca++ = __
            !   3. __ = __
            ! The corresponding species (corresponding in the set of having
            ! matching indices) would be:
            !   1. Ca++
            !   2. Na+
            !   3. __
            ! Note that there is no formal correspondence in that species 2
            ! (Na+) does not appear in reaction 2 (Ca++ = __). Also, the
            ! base site species ("__") need not appear last, depending
            ! on what condensed reactions were read from the input file.
            ! The identity reaction for this species is guaranteed to
            ! appear last, as it can't be read from the input file, and
            ! is created in this position by the current subroutine in the
            ! "Pass 1" section above.
            ! The first step is to convert the set of existing condensed
            ! reactions into a set of linearly independent dissociation
            ! reactions. In the example, the set of reactions becomes:
            !   1. Na+ = __
            !   2. Ca++ = __
            !   3. __ = __
            ! Note that it is guaranteed that the last reaction is the
            ! identity reaction for the bare site species and that at
            ! least one of the remaining reactions is a dissociation
            ! reaction. Tests performed above also guarantee that the
            ! desired conversion is possible. Basically, the process
            ! requires going through the reactions in an order which
            ! insures that for each exchange reaction there currently
            ! exists a dissociation reaction for one of the two species.
            ! First, loop through the reactions, and mark all those which
            ! are in the desired form (the identity reaction for the bare
            ! site species and all reactions in dissociation format) by
            ! setting a value of 1 in the corresponding element of the
            ! ngexro array. It is guaranteed that there is at least one
            ! reaction in dissociation format. Then repeatedly loop
            ! through the remaining reactions, each time finding one to
            ! transform from exchange to dissociation format by combining
            ! it with a reaction in the processed set. Mark it as
            ! transformed using the ngexro array. Repeated looping
            ! is required to insure that all such reactions can be
            ! transformed. Stop when no reactions remain in exchange
            ! format.
            ! The ngexro array will later be used to create a mapping
            ! between the species and the reactions.
            do iee = 1,nspect
                ngexro(iee,je,ne) = 0
            end do

            itot = 1
            ngexro(nspect,je,ne) = 1

            do ie = 1,nspect - 1
                ustr1 = ugexsr(1,ie,je,ne)
                ustr2 = ugexsr(2,ie,je,ne)

                if (ustr2(1:3) .eq. '__ ') then
                    itot = itot + 1
                    ngexro(ie,je,ne) = 1
                end if
            end do

210 continue
            if (itot .eq. nspect) then
                go to 220
            end if

            do ie = 1,nspect - 1
                if (ngexro(ie,je,ne) .eq. 0) then
                    ustr1 = ugexsr(1,ie,je,ne)
                    ustr2 = ugexsr(2,ie,je,ne)

                    do iee = 1,nspect - 1
                        if (ngexro(iee,je,ne) .eq. 1) then
                            ustr1e = ugexsr(1,iee,je,ne)

                            if (ustr1e(1:24) .eq. ustr2(1:24)) then
                                itot = itot + 1
                                ngexro(ie,je,ne) = 1
                                ugexsr(1,ie,je,ne) = ustr1
                                ugexsr(2,ie,je,ne) = '__'
                                j2 = ilnobl(ustr1)
                                ugexr(ie,je,ne) = ustr1
                                ugexr(ie,je,ne)(j2 + 1:j2 + 3) = ' = '
                                ugexr(ie,je,ne)(j2 + 4:56) = '__'
                                xx = xlkgex(ie,je,ne) + xlkgex(iee,je,ne)
                                xlkgex(ie,je,ne) = xx
                                xx = xhfgex(ie,je,ne) + xhfgex(iee,je,ne)
                                xhfgex(ie,je,ne) = xx
                                xx = xvfgex(ie,je,ne) + xvfgex(iee,je,ne)
                                xvfgex(ie,je,ne) = xx
                                go to 210
                            else if (ustr1e(1:24) .eq. ustr1(1:24)) then
                                itot = itot + 1
                                ngexro(ie,je,ne) = 1
                                ugexsr(1,ie,je,ne) = ustr2
                                ugexsr(2,ie,je,ne) = '__'
                                j2 = ilnobl(ustr2)
                                ugexr(ie,je,ne) = ustr2
                                ugexr(ie,je,ne)(j2 + 1:j2 + 3) = ' = '
                                ugexr(ie,je,ne)(j2 + 4:56) = '__'
                                xx = -xlkgex(ie,je,ne) + xlkgex(iee,je,ne)
                                xlkgex(ie,je,ne) = xx
                                xx = -xhfgex(ie,je,ne) + xhfgex(iee,je,ne)
                                xhfgex(ie,je,ne) = xx
                                xx = -xvfgex(ie,je,ne) + xvfgex(iee,je,ne)
                                xvfgex(ie,je,ne) = xx
                                go to 210
                            end if
                        end if
                    end do
                end if
            end do

            write (noutpt,1250) ugexj(je,ne)(1:j3),ugexp(ne)(1:j4)
            write (nttyo,1250) ugexj(je,ne)(1:j3),ugexp(ne)(1:j4)
1250 format(/' * Error - (EQLIB/intexi) Programming error trap:'," Couldn't transform",/7x,'one or more  condensed reactions',' for site ',a,' of exchange',/7x,'phase ',a,' from',' exchange format to dissociation',/7x,'format. Check the',' responsible coding.')

            stop

220 continue

            ! Now set up the ngexro array as a pointer array giving the
            ! index of the reaction corresponding to a given species.
            itot = 0

            do ie = 1,nspect
                ngexro(ie,je,ne) = 0
                ustr1 = ugexs(ie,je,ne)

                do iee = 1,nspect
                    if (ustr1(1:24) .eq. ugexsr(1,iee,je,ne)(1:24)) then
                        itot = itot + 1
                        ngexro(ie,je,ne) = iee
                        go to 230
                    end if
                end do

                write (noutpt,1260) ugexs(ie,je,ne)(1:j2),ugexj(je,ne)(1:j3),ugexp(ne)(1:j4)
                write (nttyo,1260) ugexs(ie,je,ne)(1:j2),ugexj(je,ne)(1:j3),ugexp(ne)(1:j4)
1260 format(/' * Error - (EQLIB/intexi) Programming error trap:'," Couldn't find",/7x,'a condensed identity or dissociation',' reaction for ',a,/7x,'on site ',a,' of exchange phase ',a,'. Check the',/7x,'responsible coding.')

                stop

230 continue
            end do

            ! Order the species for creation, using the ngexso array.
            ! Use the existing order, except that the bare site species
            ! is created first.
            itot = 0

            do ie = 1,nspect
                if (ugexs(ie,je,ne)(1:3) .eq. '__ ') then
                    itot = itot + 1
                    ngexso(itot,je,ne) = ie
                    go to 240
                end if
            end do

240 continue

            do ie = 1,nspect
                if (ugexs(ie,je,ne)(1:3) .ne. '__ ') then
                    itot = itot + 1
                    ngexso(itot,je,ne) = ie
                end if
            end do

            ! Create the species for the current site.
            jern1(je,ne) = ns + 1

            ! Loop on ions on sites. Note that the bare exchanger species
            ! will be created first.
            do ieo = 1,nspect
                ie = ngexso(ieo,je,ne)
                j2 = ilnobl(ugexs(ie,je,ne))

                if (nst .ge. nstmax) then
                    write (noutpt,1290) nstmax,ugexs(ie,je,ne)(1:j2),ugexj(je,ne)(1:j3),ugexp(ne)(1:j4)
                    write (nttyo,1290) nstmax,ugexs(ie,je,ne)(1:j2),ugexj(je,ne)(1:j3),ugexp(ne)(1:j4)
1290 format (/' * Error - (EQLIB/intexi) The maximum ',i4,' species would be',/7x,'exceeded by creating an',' exchange species for ',a,/7x,'on site ',a,' of exchange phase ',a,'. Increase the',/7x,'dimensioning parameter nstpar.')

                    stop
                end if

                ns = ns + 1
                nsl = nst
                nst = ns
                jsflag(ns) = 0

                jsitex(ns) = je
                nphasx(ns) = np

                ! If the current species is the bare one for the current
                ! site, add it to the set of strict basis species.
                nss = ngexsa(ie,je,ne)

                if (nss .eq. 0) then
                    nbt = nbt + 1

                    if (nbt .gt. nbtmax) then
                        write (noutpt,1300) nbtmax,ugexs(ie,je,ne)(1:j2),ugexj(je,ne)(1:j3),ugexp(ne)(1:j4)
                        write (nttyo,1300) nbtmax,ugexs(ie,je,ne)(1:j2),ugexj(je,ne)(1:j3),ugexp(ne)(1:j4)
1300 format (/' * Error - (EQLIB/intexi) The maximum,',i3,' basis species would be',/7x,'exceeded by',' adding a species for ',a,' on site ',a,/7x,'of',' exchange phase ',a,'. Increase the dimensioning',/7x,'parameter nbtpar.')

                        stop
                    else
                        nbasp(nbt) = ns
                        jflag(ns) = 0
                        nsba = ns
                    end if
                end if

                ! Get the index of the corresponding aqueous species
                ! (the exchangeable species). Compute the stoichiometric
                ! relationship between the exchanger species and the
                ! corresponding aqueous species. This relationship depends
                ! on the specified exchange model.
                nss = ngexsa(ie,je,ne)

                if (nss .gt. 0) then
                    zsi = zchar(nss)
                else
                    zsi = 0.
                end if

                zxj = zgexj(je,ne)
                zprod = zsi*zxj

                ! Note: cfxi below is the reaction coefficient for the
                ! exchangeable ion in the dissociation reaction for the
                ! corresponding exchanger species. The dissociation reaction
                ! is always written such that one mole of exchanger species
                ! is dissociated.
                j5 = ilnobl(ugexmo(ne))

                if (nss .gt. 0) then
                    if (ugexmo(ne)(1:5).eq.'Gapon' .or.          ugexmo(ne)(1:6).eq.'Gapon-' .or.          ugexmo(ne)(1:8).eq.'Vanselow' .or.          ugexmo(ne)(1:9).eq.'Vanselow-') then
                        ! Check electrical charges of the exchangeable ion and
                        ! the exchange site.
                        qmoerr = .false.

                        if (zsi .eq. 0.) then
                            write (noutpt,1320) ugexs(ie,je,ne)(1:j2),ugexj(je,ne)(1:j3),ugexp(ne)(1:j4),ugexmo(ne)(1:j5)
                            write (nttyo,1320) ugexs(ie,je,ne)(1:j2),ugexj(je,ne)(1:j3),ugexp(ne)(1:j4),ugexmo(ne)(1:j5)
1320 format (/" * Error - (EQLIB/intexi) Can't create",' a generic ion exchange species',/7x,'for ',a,' on site ',a,' of exchange phase',/7x,a,' for the ',a,' model because',/7x,'the',' exchangeable ion has no electrical charge.')

                            nerr = nerr + 1
                            qmoerr = .true.
                        end if

                        if (zxj .eq. 0.) then
                            write (noutpt,1330) ugexs(ie,je,ne)(1:j2),ugexj(je,ne)(1:j3),ugexp(ne)(1:j4),ugexmo(ne)(1:j5)
                            write (nttyo,1330) ugexs(ie,je,ne)(1:j2),ugexj(je,ne)(1:j3),ugexp(ne)(1:j4),ugexmo(ne)(1:j5)
1330 format (/" * Error - (EQLIB/intexi) Can't create",' a generic ion exchange species',/7x,'for ',a,' on site ',a,' of exchange phase',/7x,a,' for the ',a,' model because',/7x,'the',' the exchange site has no electrical charge.')

                            nerr = nerr + 1
                            qmoerr = .true.
                        end if

                        if (qmoerr) then
                            go to 270
                        end if

                        if (zprod .gt. 0.) then
                            write (noutpt,1340) ugexs(ie,je,ne)(1:j2),ugexj(je,ne)(1:j3),ugexp(ne)(1:j4),ugexmo(ne)(1:j5)
                            write (nttyo,1340) ugexs(ie,je,ne)(1:j2),ugexj(je,ne)(1:j3),ugexp(ne)(1:j4),ugexmo(ne)(1:j5)
1340 format (/" * Error - (EQLIB/intexi) Can't create",' a generic ion exchange species',/7x,'for ',a,' on site ',a,' of exchange phase',/7x,a,' for the ',a,' model because the exchangeable',/7x,'ion and',' the exchange site have the same electrical',' charge sign.')

                            nerr = nerr + 1
                            go to 270
                        end if
                    end if
                end if

                cgx = cgexj(je,ne)
                zxjt = cgx*zxj

                ! cfxi = the number of moles of exchangeable ion per
                !          mole of the exchanger species
                ! cfxz = the number of moles of exchanger substrate per
                !          mole of exchanger species
                if (nss .gt. 0) then
                    ! Have other than the bare site species.
                    if (ugexmo(ne)(1:j5).eq.'Gapon' .or.          ugexmo(ne)(1:6).eq.'Gapon-') then
                        ! Gapon (Gapon-?) model.
                        !   For:
                        !     cgexj(je,ne) = 1
                        !     zgexj(je,ne) = -1
                        !     "Na+ = Ca++"
                        !   the species are:
                        !     Na-Z and Ca[1/2]-Z
                        cfxi = -cgx*zxj/zsi
                        cfxz = -cfxi*zsi/zxjt
                    else if (ugexmo(ne)(1:j5).eq.'Vanselow' .or.          ugexmo(ne)(1:9).eq.'Vanselow-') then
                        ! Vanselow (Vanselow-?) model.
                        !   For:
                        !     cgexj(je,ne) = 1
                        !     zgexj(je,ne) = -1
                        !     "Na+ = Ca++"
                        !   the species are:
                        !     Na-Z and Ca-Z2
                        cfxi = -cgx*zxj*abs(zsi)/zsi
                        cfxz = -cfxi*zsi/zxjt
                    else if (ugexmo(ne)(1:j5) .eq. 'Site-mixing') then
                        ! Site-mixing model.
                        cfxi = cgx
                        cfxz = 1.0
                    else
                        write (noutpt,1350) ugexmo(ne)(1:j5),ugexp(ne)(1:j4)
                        write (nttyo,1350) ugexmo(ne)(1:j5),ugexp(ne)(1:j4)
1350 format (/' * Error - (EQLIB/intexi.f) Programming'," error trap: Don't recognize the",/7x,'exchange',' model type"',a,'" specified for the exchanger phase',/7x,a,'. Valid choices include "Gapon" and "Vanselow".',/7x,'Invalid choices should have been trapped above in',' this subroutine.',/7x,'Check the coding at the',' present point. Now trying to compute',/7x,'factors',' for the compositions and reactions of species',' belonging',/7x,'to this exchanger phase.')

                        stop
                    end if
                else
                    ! Have the bare site species.
                    cfxi = cgx
                    cfxz = 1.0
                end if

                ! Here efx is the number of equivalents of exchangeable
                ! ion appearing in the dissociation reaction.
                if (nss .gt. 0) then
                    ! Have other than the bare site species.
                    zchar(ns) = cfxi*zsi - zxjt
                    mwtsp(ns) = cfxi*mwtsp(nss)
                    efx = abs(cfxi*zsi)
                else
                    ! Have the bare site species.
                    zchar(ns) = zxjt
                    mwtsp(ns) = 0.
                    efx = 0.
                end if

                ! Compose the full name of the exchanger species. The phase
                ! part of the name is the exchange phase name. The species
                ! part is a concatenation of the species part of the name of
                ! the ion released in an exchange reaction (the aqueous
                ! ion) and the site name. A blank space is included in the
                ! concatenation between the two parts. The concatenation is
                ! trimmed as necessary to fit into 24 characters.
                ! Approximately two characters of the exchange ion name are
                ! trimmed for each character of the site name. The full name
                ! of the species should should then look like
                ! "Na+ S(1)                Exchanger(A)". This algorithm
                ! matches that in EQLIB/intgex.f. As formatted by
                ! EQLIBU/fmspnm.f or EQLIBU/fmspnx.f, the name in this
                ! example would appear as "Na+ S(1) Exchanger(A)".
                ! Warning: the stoichiometry of the species can't be deduced
                ! from the name alone. The the number of formula units of
                ! exchangeable Na+ in site S(1) in the above example could be
                ! 1.0, 2.0, or some other number. The same is true for the
                ! number of moles of the site itself per mole of exchanger.
                ! Thus, "Na+ S(1) Exchanger(A)" has one stoichiometry for
                ! the Vanselow model, but a different one for the Gapon model.
                jj2 = j2
                jj3 = j3
                nd2 = 0
250 continue
                jj = jj2 + jj3 + 1

                if (jj .gt. 24) then
                    if (nd2 .lt. 2) then
                        jj2 = jj2 - 1
                        nd2 = nd2 + 1
                        go to 250
                    else
                        jj3 = jj3 - 1
                        nd2 = 0
                        go to 250
                    end if
                end if

                uspec(ns) = ' '
                uspec(ns)(1:jj2) = ugexs(ie,je,ne)(1:jj2)
                uspec(ns)(jj2 + 1:jj2 + 1) = ' '
                uspec(ns)(jj2 + 2:jj) = ugexj(je,ne)(1:jj3)
                uspec(ns)(25:48) = ugexp(ne)(1:24)

                ! Set up the elemental composition.
                nrf1 = nessr(2,nsl) + 1

                if (nss .gt. 0) then
                    nr1 = nessr(1,nss)
                    nr2 = nessr(2,nss)
                    nrf2 = nrf1 + nr2 - nr1
                else
                    nrf2 = nrf1
                end if

                if (nrf2 .gt. nessmx) then
                    write (noutpt,1380) nessmx,ugexs(ie,je,ne)(1:j2),ugexj(je,ne)(1:j3),ugexp(ne)(1:j4)
                    write (nttyo,1380) nessmx,ugexs(ie,je,ne)(1:j2),ugexj(je,ne)(1:j3),ugexp(ne)(1:j4)
1380 format (/' * Error - (EQLIB/intexi) The maximum ',i5,' entries in the',/7x,'cess/ness arrays would be',' exceeded by creating the',/7x,'exchange species',' for ',a,' on site',/7x,a,' of exchange phase ',a,'.',' Increase the',/7x,'dimensioning parameter nesspa.')

                    stop
                end if

                nessr(1,ns) = nrf1
                nessr(2,ns) = nrf2

                if (nss .gt. 0) then
                    k = nr1 - 1

                    do n = nrf1,nrf2
                        k = k + 1
                        cess(n) = cfxi*cess(k)
                        ness(n) = ness(k)
                    end do
                else
                    n = nrf1
                    cess(n) = 0.
                    ness(n) = 0
                end if

                ! Set up the corresponding reaction.
                nrf1 = ndrsr(2,nsl) + 1

                if (nss .gt. 0) then
                    ! Case of a non-bare site species.
                    nrf2 = nrf1 + 2
                else
                    ! Case of the bare site species.
                    nrf2 = nrf1
                end if

                if (nrf2 .gt. ndrsmx) then
                    write (noutpt,1390) ndrsmx,ugexs(ie,je,ne)(1:j2),ugexj(je,ne)(1:j3),ugexp(ne)(1:j4)
                    write (nttyo,1390) ndrsmx,ugexs(ie,je,ne)(1:j2),ugexj(je,ne)(1:j3),ugexp(ne)(1:j4)
1390 format (/' * Error - (EQLIB/intexi) The maximum ',i5,' entries in the',/7x,'cdrs/ndrs arrays would be',' exceeded by creating the',/7x,'exchange species',' for ',a,'on site',/7x,a,' of exchange phase ',a,'.',' Increase the',/7x,'dimensioning parameter ndrspa.')

                    stop
                end if

                ndrsr(1,ns) = nrf1
                ndrsr(2,ns) = nrf2

                if (nss .gt. 0) then
                    ! Case of a non-bare site species. Create a reaction in
                    ! which the species dissociates to the bare site species.
                    n = nrf1
                    cdrs(n) = -1.0
                    ndrs(n) = ns
                    n = n + 1
                    cdrs(n) = cfxi
                    ndrs(n) = nss
                    n = n + 1
                    cdrs(n) = cfxz
                    ndrs(n) = nsba
                else
                    ! Case of the bare site species. Create a null reaction.
                    ! The bare site species then becomes a strict basis
                    ! species, though it doesn't actually formally correspond
                    ! to a chemical element.
                    n = nrf1
                    cdrs(n) = 0.
                    ndrs(n) = 0
                end if

                if (nss .eq. 0) then
                    go to 270
                end if

                ! Set up to the thermodynamic properties of the reaction
                ! just created. This reaction should match one of the
                ! existing condensed reactions.
                iee = ngexro(ie,je,ne)

                ! Calculate the thermodynamic properties of the reaction.
                do j = 1,ntprt
                    if (efx .ne. 0.) then
                        afxb = efx*xlkgex(iee,je,ne)
                        afhx = efx*xhfgex(iee,je,ne)
                        afvx = efx*xvfgex(iee,je,ne)
                    else
                        afxb = xlkgex(iee,je,ne)
                        afhx = xhfgex(iee,je,ne)
                        afvx = xvfgex(iee,je,ne)
                    end if

                    ! Here afxb is the log K value and tfxb the temperature (C)
                    ! at the base temperature. Use the van't Hoff relation to
                    ! find afx0, the log K at 0C. Approximate the van't Hoff
                    ! relation across the total temperature range using the
                    ! standard power series.
                    bfx = -afhx/(arcnst*tfx0)
                    afx0 = afxb + bfx*(tfxb/(tfxb + tfx0))
                    axlks(1,j,ns) = afx0
                    xx = -bfx

                    do i = 2,narxt(j)
                        xx = -xx/tfx0
                        axlks(i,j,ns) = xx
                    end do

                    axhfs(1,j,ns) = afhx
                    axvfs(1,j,ns) = afvx

                    do i = 2,narxt(j)
                        axhfs(i,j,ns) = 0.
                        axvfs(i,j,ns) = 0.
                    end do
                end do

                ! Is the aqueous species which goes onto the substrate a
                ! basis species? If not, the reaction must be rewritten
                ! in terms of equivalent basis species.
                nrf1 = ndrsr(1,ns)
                nrf2 = ndrsr(2,ns)

                do n = nrf1 + 1,nrf2
                    nse = ndrs(n)

                    if (nse.ge.narn1 .and. nse.le.narn2) then
                        ! Calling sequence substitutions:
                        !   nse for ns
                        nb = nbasis(nbasp,nbt,nbtmax,nse)

                        if (nb .eq. 0) then
                            nbt = nbt + 1

                            if (nbt .gt. nbtmax) then
                                j6 = ilnobl(uspec(nse))
                                write (noutpt,1400) nbtmax,uspec(nse)(1:j6),ugexs(ie,je,ne)(1:j2),ugexj(je,ne)(1:j3),ugexp(ne)(1:j4)
                                write (nttyo,1400) nbtmax,uspec(nse)(1:j6),ugexs(ie,je,ne)(1:j2),ugexj(je,ne)(1:j3),ugexp(ne)(1:j4)
1400 format (/' * Error - (EQLIB/intexi) The maximum,',i3,' basis species would be',/7x,'exceeded by',' adding ',a,' in order to',/7x,'accomodate the',' required setup for species ',a,/7x,'on site ',a,' of exchange phase ',a,'.',/7x,'Increase the',' dimensioning parameter nbtpar.')

                                stop
                            else
                                nbasp(nbt) = nse
                                jflag(nse) = 30
                            end if
                        end if
                    end if
                end do

270 continue
            end do

            ! Modify the ngexsa array if necessary so that the species
            ! index of the bare site species appears in the first position
            ! for the current site of the current exchange phase.
            do ie = 1,nspect
                nss = ngexsa(ie,je,ne)

                if (nss .eq. 0) then
                    do i = 1,ie - 1
                        iee = ie - i
                        ngexsa(iee + 1,je,ne) = ngexsa(iee,je,ne)
                    end do

                    ngexsa(1,je,ne) = 0
                    go to 280
                end if
            end do

280 continue

            jern2(je,ne) = nst
            ngext(je,ne) = jern2(je,ne) - jern1(je,ne) + 1
        end do

        ncmpr(2,np) = nst
    end do

    nern2 = nst
    iern2 = npt

    ! Pass 3. Loop on exchanger phases. Change the ordering of the
    ! condensed reactions and associated data to match that of the
    ! corresponding created species. Then eliminate the identity
    ! reaction for the bare site species. This puts these data
    ! in the desired form for writing on a pickup file.
    ! In the example discussed above (input of a single exchange
    ! reaction, Ca++ = Na+), the expanded set of condensed reactions
    ! was:
    !       1. Na+ = __
    !       2. Ca++ = __
    !       3. __ = __
    ! The actual species created were processed in the order
    ! __, Ca++, Na+. The bare site species is always processed first.
    ! The other species are processed in order of appearance in the
    ! condensed reactions read from the input file. Thus, the order of
    ! creation of the corresponding reactions from the condensed
    ! counterparts was:
    !       1. (3. __ = __)
    !       2. (2. Ca++ = __)
    !       3. (1. Na+ = __)
    ! Here the set of condensed reactions is changed to:
    !       1. Ca++ = __
    !       2. Na+ = __
    ! As input, this is expanded to:
    !       1. Ca++ = __
    !       2. Na+ = __
    !       3. __ = __
    ! The actual species are then created in the same order as before,
    ! (__, Ca++, Na+).
    ! Note that less setup work is required for the case of the input
    ! of two dissociation reactions than for the case of the input of
    ! the single equivalent exchange reaction.
    do ne = 1,net
        ! Loop on sites.
        do je = 1,jgext(ne)
            nspect = ngexrt(je,ne)

            ! Re-arrange the condensed reactions so that their order
            ! matches that of the local species list. This eliminates
            ! the need for the ngexro pointer array to find the
            ! condensed reaction for a given species in the local
            ! species list.
            do ie = 1,nspect
                iee = ngexro(ie,je,ne)

                if (iee .ne. ie) then
                    ! Exchange positions.
                    ustr56 = ugexr(iee,je,ne)
                    ugexr(iee,je,ne) = ugexr(ie,je,ne)
                    ugexr(ie,je,ne) = ustr56
                    xx = xlkgex(iee,je,ne)
                    xlkgex(iee,je,ne) = xlkgex(ie,je,ne)
                    xlkgex(ie,je,ne) = xx
                    xx = xhfgex(iee,je,ne)
                    xhfgex(iee,je,ne) = xhfgex(ie,je,ne)
                    xhfgex(ie,je,ne) = xx
                    xx = xvfgex(iee,je,ne)
                    xvfgex(iee,je,ne) = xvfgex(ie,je,ne)
                    xvfgex(ie,je,ne) = xx

                    do iej = 1,nspect
                        if (ngexro(iej,je,ne) .eq. ie) then
                            go to 310
                        end if
                    end do

310 continue
                    ngexro(iej,je,ne) = iee
                    ngexro(ie,je,ne) = ie
                end if
            end do

            do ie = 1,nspect
                iee = ngexso(ie,je,ne)

                if (iee .ne. ie) then
                    ! Exchange positions.
                    ustr56 = ugexr(iee,je,ne)
                    ugexr(iee,je,ne) = ugexr(ie,je,ne)
                    ugexr(ie,je,ne) = ustr56
                    xx = xlkgex(iee,je,ne)
                    xlkgex(iee,je,ne) = xlkgex(ie,je,ne)
                    xlkgex(ie,je,ne) = xx
                    xx = xhfgex(iee,je,ne)
                    xhfgex(iee,je,ne) = xhfgex(ie,je,ne)
                    xhfgex(ie,je,ne) = xx
                    xx = xvfgex(iee,je,ne)
                    xvfgex(iee,je,ne) = xvfgex(ie,je,ne)
                    xvfgex(ie,je,ne) = xx

                    do iej = 1,nspect
                        if (ngexso(iej,je,ne) .eq. ie) then
                            go to 320
                        end if
                    end do

320 continue
                    ngexso(iej,je,ne) = iee
                    ngexso(ie,je,ne) = ie
                end if
            end do

            ! Remove the identity reaction from the set of condensed
            ! reactions.
            do ie = 1,nspect - 1
                iee = ie + 1
                ugexr(ie,je,ne) = ugexr(iee,je,ne)
                xlkgex(ie,je,ne) = xlkgex(iee,je,ne)
                xhfgex(ie,je,ne) = xhfgex(iee,je,ne)
                xvfgex(ie,je,ne) = xvfgex(iee,je,ne)
            end do

            ie = nspect
            ugexr(ie,je,ne) = '__ = __ '
            xlkgex(ie,je,ne) = 0.
            xhfgex(ie,je,ne) = 0.
            xvfgex(ie,je,ne) = 0.
            nspect = nspect - 1
            ngexrt(je,ne) = nspect
        end do
    end do

    ! Calculate the cegexs array. This is an array of coefficients
    ! giving the number of equivalents per mole of each exchanger
    ! species. This supports the calculation of equivalent fractions.
    do ne = 1,net
        do je = 1,jgext(ne)
            ns = jern1(je,ne) - 1

            do ie = 1,ngext(je,ne)
                ns = ns + 1

                if (uspec(ns)(1:3) .eq. '__ ') then
                else
                    ! Have the bare site species.
                    cegexs(ie,je,ne) = 0.

                    ! Have other than the bare site species.
                    nr1 = ndrsr(1,ns)
                    nr2 = ndrsr(2,ns)
                    cx = 0.

                    do n = nr1 + 1,nr2
                        nse = ndrs(n)

                        ! There are contributions only from the aqueous
                        ! species. There is none from the bare site species,
                        ! which also appears in the reaction.
                        if (nse.ge.narn1 .and. nse.le.narn2) then
                            cx = cx + cdrs(n)*zchar(nse)
                        end if
                    end do

                    cegexs(ie,je,ne) = cx
                end if
            end do
        end do
    end do

    ! Calculate the cpgexs array. This is an array of coefficients
    ! giving the number of moles of exchanger substrate (Z) per mole of
    ! exchanger species (e.g., Na-Z, Ca-Z2). This supports the
    ! calculation of the number of moles of the exchanger phase.
    do ne = 1,net
        j5 = ilnobl(ugexmo(ne))

        if (ugexmo(ne)(1:j5).eq.'Gapon' .or.    ugexmo(ne)(1:6).eq.'Gapon-' .or.    ugexmo(ne)(1:j5) .eq. 'Site-mixing') then
            ! Gapon (Gapon-?) or Site-mixing model.
            do je = 1,jgext(ne)
                do ie = 1,ngext(je,ne)
                    cpgexs(ie,je,ne) = 1.0
                end do
            end do
        else if (ugexmo(ne)(1:j5).eq.'Vanselow' .or.    ugexmo(ne)(1:9).eq.'Vanselow-') then
            ! Vanselow (Vanselow-?) model.
            do je = 1,jgext(ne)
                do ie = 1,ngext(je,ne)
                    cpgexs(ie,je,ne) = cegexs(ie,je,ne)/egexjf(je,ne)
                end do
            end do
        end if
    end do

    ! Set up the kgexsa, kern1, and kern2 pointer arrays. The kgexsa
    ! array contains the indices of the exchange species (e.g., Na+).
    ! The latter two arrays point to the start and end of the range
    ! in kgexsa corresponding to a given generic ion exchanger phase.
    ! Here kgexsa(ke,ne) is the index of the ke-th species exchanging
    ! on the ne-th generic ion exchanger phase.
    ke = 0

    do ne = 1,net
        ke1 = ke + 1
        kern1(ne) = ke1

        do je = 1,jgext(ne)
            do ie = 1,ngext(je,ne)
                nss = ngexsa(ie,je,ne)

                ! Has this exchange species already been loaded because
                ! because it appears in another site?
                do kee = ke1,ke
                    nsse = kgexsa(kee,ne)

                    if (nsse .eq. nss) then
                        go to 400
                    end if
                end do

                ! Have an exchange species which has not already been
                ! loaded. Load it.
                ke = ke + 1
                kgexsa(ke,ne) = nss
400 continue
            end do
        end do

        kern2(ne) = ke
    end do

990 continue
    if (nerr .gt. 0) then
        stop
    end if
end subroutine intexi
