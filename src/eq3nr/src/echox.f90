subroutine echox(azero,cdrs,covali,eh,fo2lg,iebal3,iodb,iopg,iopr,iopt,irdxc3,itdsf3,itermx,ixrn1,ixrn2,jflag,jflgi,jpres3,narn1,narn2,nat,nata,natmax,nbt,nbta,nbti,nbtmax,nct,ncta,nctmax,ncmpr,ncosp,ndecsp,ndrs,ndrsmx,ndrsr,nelect,ngt,ngta,ngtmax,nhydr,nhydx,njfmax,nlt,nlta,nltmax,nmt,nmta,nmtmax,nodbmx,nopgmx,noprmx,noptmx,noutpt,no2gaq,npt,npta,nptmax,nredox,nst,nsta,nstmax,nxt,nxta,nxti,nxtmax,pe,presg,press,rho,scamas,tdspkg,tdspl,tempc,tolbt,toldl,tolspf,uactop,ucospi,uebal,ujflls,uphase,uredox,uspec,uspeci,xbar)
    !! This subroutine echoes the problem input defining the
    !! compositional constraints on the aqueous solution being modeled
    !! and the options and tolerances selected.
    !! This subroutine is called by:
    !!   EQ3NR/eq3nr.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: natmax
    integer :: nbtmax
    integer :: nctmax
    integer :: ndrsmx
    integer :: ngtmax
    integer :: njfmax
    integer :: nltmax
    integer :: nmtmax
    integer :: nodbmx
    integer :: nopgmx
    integer :: noprmx
    integer :: noptmx
    integer :: nptmax
    integer :: nstmax
    integer :: nxtmax

    integer :: noutpt

    integer :: iodb(nodbmx)
    integer :: iopg(nopgmx)
    integer :: iopr(noprmx)
    integer :: iopt(noptmx)
    integer :: jflag(nstmax)
    integer :: jflgi(nbtmax)
    integer :: ncmpr(2,nptmax)
    integer :: ncosp(nbtmax)
    integer :: ndecsp(nbtmax)
    integer :: ndrs(ndrsmx)
    integer :: ndrsr(2,nstmax)

    integer :: iebal3
    integer :: irdxc3
    integer :: itdsf3
    integer :: itermx
    integer :: ixrn1
    integer :: ixrn2
    integer :: jpres3
    integer :: narn1
    integer :: narn2
    integer :: nat
    integer :: nata
    integer :: nbt
    integer :: nbta
    integer :: nbti
    integer :: nct
    integer :: ncta
    integer :: nelect
    integer :: ngt
    integer :: ngta
    integer :: nhydr
    integer :: nhydx
    integer :: nlt
    integer :: nlta
    integer :: nmt
    integer :: nmta
    integer :: no2gaq
    integer :: npt
    integer :: npta
    integer :: nredox
    integer :: nst
    integer :: nsta
    integer :: nxt
    integer :: nxta
    integer :: nxti

    character(len=48) :: ucospi(nbtmax)
    character(len=48) :: uspec(nstmax)
    character(len=48) :: uspeci(nbtmax)
    character(len=32) :: ujflls(0:njfmax)
    character(len=32) :: uactop
    character(len=24) :: uphase(nptmax)
    character(len=24) :: uebal
    character(len=24) :: uredox

    real(kind=8) :: azero(natmax)
    real(kind=8) :: cdrs(ndrsmx)
    real(kind=8) :: covali(nbtmax)
    real(kind=8) :: xbar(nstmax)

    real(kind=8) :: dp
    real(kind=8) :: eh
    real(kind=8) :: fo2lg
    real(kind=8) :: pe
    real(kind=8) :: presg
    real(kind=8) :: press
    real(kind=8) :: rho
    real(kind=8) :: scamas
    real(kind=8) :: tdspkg
    real(kind=8) :: tdspl
    real(kind=8) :: tempc
    real(kind=8) :: tolbt
    real(kind=8) :: toldl
    real(kind=8) :: tolspf

    ! Local variable declarations.
    integer :: jfl
    integer :: j2
    integer :: j3
    integer :: j4
    integer :: n
    integer :: na
    integer :: nbi
    integer :: nb1
    integer :: np
    integer :: nr1
    integer :: nr2
    integer :: ns
    integer :: nss

    integer :: ilnobl

    character(len=24) :: ux24

    j2 = ilnobl(uactop)
    write (noutpt,1000) uactop(1:j2)
1000 format(/' The activity coefficients of aqueous species will be',/' calculated using ',a,'.'/)

    if (iopr(3) .gt. 1) then
        write (noutpt,1010)
1010 format(/7x,'Species',17x,'HC Diameter',/)

        do ns = narn1 + 1,narn2
            na = ns - narn1 + 1
            write (noutpt,1020) uspec(ns),azero(na)
1020 format(5x,a24,3x,f7.3)
        end do
    end if

    write (noutpt,1030) tempc
1030 format(/' Temperature= ',f6.2,' C',/)

    write (noutpt,1032) jpres3
1032 format(/' jpres3=  ',i3,' (Pressure option switch)')

    if (jpres3 .eq. 0) then
        write (noutpt,1034) press
1034 format(/'   Pressure= ',1pg12.5,' bars (data file reference',' curve value)',/)
    else if (jpres3 .eq. 1) then
        dp = press - presg
        write (noutpt,1036) press,presg,dp
1036 format(/'   Pressure= ',1pg12.5,' bars (1.013-bar/steam-','saturation curve value)',/'   Data file reference curve pressure= ',g12.5,' bars',/'   Pressure difference= ',g12.5,' bars',/)
    else if (jpres3 .eq. 2) then
        dp = press - presg
        write (noutpt,1038) press,presg,dp
1038 format(/'   Pressure= ',1pg12.5,' bars (specified value)',/'   Data file reference curve pressure= ',g12.5,' bars',/'   Pressure difference= ',g12.5,' bars',/)
    end if

    ! Print a table showing for phases, species, and groups thereof,
    ! the number of each of entity on the data base, the number the
    ! software is dimensioned for, and the number appearing in the
    ! current problem.
    call prtntt(nat,nata,natmax,nbt,nbta,nbtmax,nct,ncta,nctmax,ngt,ngta,ngtmax,nlt,nlta,nltmax,nmt,nmta,nmtmax,noutpt,npt,npta,nptmax,nst,nsta,nstmax,nxt,nxta,nxtmax)

    write (noutpt,1040) (iopt(n), n = 1,10)
1040 format(/' iopt(1)=  ',i2,' (Used only by EQ6)',/' iopt(2)=  ',i2,' (Used only by EQ6)',/' iopt(3)=  ',i2,' (Used only by EQ6)',/' iopt(4)=  ',i2,' (Solid solutions)',/' iopt(5)=  ',i2,' (Used only by EQ6)',/' iopt(6)=  ',i2,' (Used only by EQ6)',/' iopt(7)=  ',i2,' (Not used)',/' iopt(8)=  ',i2,' (Not used)',/' iopt(9)=  ',i2,' (Not used)',/' iopt(10)= ',i2,' (Not used)')

    write (noutpt,1045) (iopt(n), n = 11,18)
1045 format(' iopt(11)= ',i2,' (Auto basis switching, in',' pre-Newton-Raphson optimization)',/' iopt(12)= ',i2,' (Used only by EQ6)',/' iopt(13)= ',i2,' (Not used)',/' iopt(14)= ',i2,' (Not used)',/' iopt(15)= ',i2,' (Used only by EQ6)',/' iopt(16)= ',i2,' (Not used)',/' iopt(17)= ',i2,' (pickup file options)',/' iopt(18)= ',i2,' (Used only by EQ6)',/' iopt(19)= ',i2,' (Advanced EQ3NR pickup file options)',/)

    write (noutpt,1050) (iopg(n), n = 1,2)
1050 format(/'   iopg(1)=  ',i2,' (Aqueous species activity',' coefficient model)',/'   iopg(2)=  ',i2,' (pH scale)',/)

    write (noutpt,1060) (iopr(n), n = 1,17)
1060 format(/' iopr(1)=  ',i2,' (List all species)',/' iopr(2)=  ',i2,' (List all reactions)',/' iopr(3)=  ',i2,' (List HC diameters)',/' iopr(4)=  ',i2,' (Aqueous species concentration print',' cut-off)',/' iopr(5)=  ',i2,' (Ion/H+ activity ratios)',/' iopr(6)=  ',i2,' (Mass balance percentages)',/' iopr(7)=  ',i2,' (Affinity print cut-off)',/' iopr(8)=  ',i2,' (Fugacities)',/' iopr(9)=  ',i2,' (Mean molal activity coefficients)',/' iopr(10)= ',i2,' (Pitzer coefficients tabulation)',/' iopr(11)= ',i2,' (Not used)',/' iopr(12)= ',i2,' (Not used)',/' iopr(13)= ',i2,' (Not used)',/' iopr(14)= ',i2,' (Not used)',/' iopr(15)= ',i2,' (Not used)',/' iopr(16)= ',i2,' (Not used)',/' iopr(17)= ',i2,' (pickup file format)',/)

    write (noutpt,1070) (iodb(n), n = 1,7)
1070 format(/'   iodb(1)=  ',i2,' (General diagnostics)',/'   iodb(2)=  ',i2,' (Used only by EQ6)',/'   iodb(3)=  ',i2,' (pre-Newton-Raphson optimization',' iterations)',/'   iodb(4)=  ',i2,' (Newton-Raphson iterations)',/'   iodb(5)=  ',i2,' (Used only by EQ6)',/'   iodb(6)=  ',i2,' (Hypothetical affinity iterations)',/'   iodb(7)=  ',i2,' (Used only by EQ6)',/)

    jfl = jflag(no2gaq)
    write (noutpt,1080) irdxc3
1080 format(/' irdxc3=  ',i3,' (Default redox constraint switch)')

    if (irdxc3 .eq. 1) then
        if (uredox(1:7).ne.'O2(aq) ' .and.    uredox(1:7).ne.'H2(aq) ') then
            nr1 = ndrsr(1,nredox)
            nr2 = ndrsr(2,nredox)

            do n = nr1,nr2
                nss = ndrs(n)

                if (nss.ne.nredox .and. nss.ne.narn1 .and. nss.ne.nhydr      .and. nss.ne.nhydx .and. nss.ne.no2gaq .and. nss.ne.nelect) then
                    go to 100
                end if
            end do

100 continue
        else
            nss = narn1
        end if

        j2 = ilnobl(uredox)
        j3 = ilnobl(uspec(nss)(1:24))
        write (noutpt,1100) uredox(1:j2),uspec(nss)(1:j3)
1100 format(/'   The default redox state is constrained by the',/'   ',a,'/',a,' couple.',/)
    else if (irdxc3 .eq. 0) then
        write (noutpt,1110) fo2lg
1110 format(/'   The default redox state is constrained by',' Log fO2 = ',g12.5,' (log bars).',/)
    else if (irdxc3 .eq. -1) then
        write (noutpt,1120) eh
1120 format(/'   The default redox state is constrained by',' Eh = ',f8.5,' volts.',/)
    else if (irdxc3 .eq. -2) then
        write (noutpt,1130) pe
1130 format(/'   The default redox state is constrained by',' pe- = ',1pe12.5,'.',/)
    else if (irdxc3 .eq. -3) then
        write (noutpt,1140)
1140 format(/'   The default redox state is controlled by a',/' heterogeneous reaction (see below).',/)
    end if

    write (noutpt,1142) iebal3
1142 format(/' iebal3=  ',i3,' (Electrical balancing option switch)')

    if (iebal3 .le. 0) then
        write (noutpt,1144)
1144 format(/'   No electrical balancing adjustment will be made.',/'   The imbalance will be calculated.',/)
    else if (iebal3 .eq. 1) then
        j2 = ilnobl(uebal)
        write (noutpt,1146) uebal(1:j2)
1146 format(/'   The species ',a,' will be adjusted to',/'   achieve electrical balance.',/)
    end if

    write (noutpt,1200) rho
1200 format(/' Solution density = ',f8.5,' g/ml',/)

    write (noutpt,1202) itdsf3
1202 format(/' itdsf3=  ',i3,' (Total dissolved solutes option',' switch)')

    if (itdsf3 .le. 0) then
        write (noutpt,1204) tdspkg
1204 format( /'   Total dissolved salts = ',f10.2,' mg/kg.sol',/)
    else
        write (noutpt,1206) tdspl
1206 format(/'   Total dissolved salts = ',f10.2,' mg/L',/)
    end if

    write (noutpt,1210) tolbt,toldl,tolspf
1210 format(/' tolbt  = ',1pe12.5,' (convergence tolerance on',' residual functions)',/' toldl  = ',e12.5,' (convergence',' tolerance on correction terms)',/' tolspf = ',e12.5,' (saturation print flag tolerance, does not affect',/25x,'convergence)',/)

    write (noutpt,1220) itermx
1220 format(/' itermx = ',i3,' (maximum number of iterations)',/)

    write (noutpt,1230) scamas
1230 format(/' scamas = ',1pe12.5,' (scale factor for aqueous',' solution',/25x,'mass written on the pickup file)',/)

    write (noutpt,1300)
1300 format(/21x,'--- Original Input Constraints ---',//5x,'Species',20x,'coval   jflag   Type of Input',/)

    do nbi = 1,nbti
        jfl = jflgi(nbi)
        nb1 = ndecsp(nbi)

        if (nb1 .gt. 0) then
            j2 = ilnobl(ujflls(jfl))

            if (jfl.eq.17 .or. jfl.eq.18) then
                write (noutpt,1310) uspeci(nbi),covali(nbi),jfl,ujflls(jfl)(1:j2)
1310 format(2x,a24,2x,1pe12.5,2x,i2,2x,a)

                j3 = ilnobl(ucospi(nbi)(1:24))
                write (noutpt,1315) ucospi(nbi)(1:j3)
1315 format(39x,'Counterion= ',a)
            else if (jfl .eq. 21) then
                write (noutpt,1310) uspeci(nbi),covali(nbi),jfl,ujflls(jfl)(1:j2)
            else if (jfl .eq. 25) then
                write (noutpt,1310) uspeci(nbi),covali(nbi),jfl,ujflls(jfl)(1:j2)
                ns = ncosp(nb1)
                j3 = ilnobl(ucospi(nbi)(1:24))
                ux24 = ucospi(nbi)(25:48)
                j4 = ilnobl(ux24)

                if (j4 .le. 0) then
                    ux24 = ucospi(nbi)(1:24)
                    j4 = j3
                end if

                write (noutpt,1320) ucospi(nbi)(1:j3),ux24(1:j4)
1320 format(46x,'Species= ',a,/48x,'Phase= ',a)

                ! Calling sequence substitutions:
                !   noutpt for nf
                call prreac(cdrs,ndrs,ndrsmx,ndrsr,noutpt,ns,nstmax,uspec)
                write (noutpt,1322)
1322 format(1x)
            else if (jfl.eq.27 .or. jfl.eq.30) then
                write (noutpt,1323) uspeci(nbi),jfl,ujflls(jfl)(1:j2)
1323 format(2x,a24,16x,i2,2x,a)
            else if (uspeci(nbi)(1:6) .eq. 'O2(g) ') then
                write (noutpt,1325) uspeci(nbi),covali(nbi),jfl
1325 format(2x,a24,2x,1pe12.5,2x,i2,2x,'Log fO2')
            else
                if (jfl .ge. 0) then
                    write (noutpt,1310) uspeci(nbi),covali(nbi),jfl,ujflls(jfl)(1:j2)
                else
                    write (noutpt,1327) uspeci(nbi)
1327 format(2x,a24,16x,'Not present in the model')
                end if
            end if
        end if
    end do

    write (noutpt,1340)
1340 format(1x)

    if (nxti .gt. 0) then
        write (noutpt,1600)
1600 format(/21x,'--- Input Solid Solution Compositions ---',/)

        do np = ixrn1,ixrn2
            nr1 = ncmpr(1,np)
            nr2 = ncmpr(2,np)

            do ns = nr1,nr2
                if (xbar(ns) .ne. 0.) then
                    go to 200
                end if
            end do

            go to 190

200 continue
            j2 = ilnobl(uphase(np))
            write (noutpt,1610) uphase(np)(1:j2)
1610 format(/6x,a,/)

            write (noutpt,1620)
1620 format(30x,'Mole Fraction')

            do ns = nr1,nr2
                write (noutpt,1630) uspec(ns),xbar(ns)
1630 format(12x,a24,3x,f6.4)
            end do

190 continue
        end do

        write (noutpt,1340)
    end if

    write (noutpt,1900)
1900 format(/' - - - - - - - - - - - - - - - - - - - - - - - - ','- - - - - - ',/)
end subroutine echox