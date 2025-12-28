subroutine prtgex(acflg,actlg,affpd,affsd,cegexs,conc,egexjc,egexjf,egexpa,egexpc,egexs,egexw,iern1,iern2,ietmax,jern1,jern2,jetmax,jgext,kern1,kern2,ketmax,kgexsa,moph,mosp,netmax,ngexsa,ngext,noutpt,np,nptmax,nstmax,sidrph,sidrsp,tolspf,ugexj,ugexmo,uspec,uphase,xbar,xbarlg,xgexw,wkgwi)
    !! This subroutine prints tables describing the state and properties
    !! of the np-th phase (a generic ion exchanger). It is analogous to
    !! EQLIB/prtsso.f, which prints the same tables for solid solution
    !! phases.
    !! This subroutine is called by:
    !!   EQ3NR/scripx.f
    !!   EQ6/scripz.f
    !! Principal input:
    !!   acflg  = array of log activity coefficients
    !!   actlg  = array of log activities of species
    !!   affpd  = array of phase affinities
    !!   affsd  = array of species affinities
    !!   sidrph = array of phase saturation indices
    !!   sidrsp = array of species saturation indices
    !!   tolspf = saturation print flag tolerance, used to flag those
    !!              phases which are close to saturation
    !!   uphase = array of phase names
    !!   xbar   = array of mole fractions of species
    !!   xbarlg = array of log mole fractions of species
    !! Principal output:
    !!   None
    implicit none

    ! Calling sequence variable declarations.
    integer :: ietmax
    integer :: jetmax
    integer :: ketmax
    integer :: netmax
    integer :: nptmax
    integer :: nstmax

    integer :: noutpt

    integer :: jern1(jetmax,netmax)
    integer :: jern2(jetmax,netmax)
    integer :: jgext(netmax)
    integer :: kern1(netmax)
    integer :: kern2(netmax)
    integer :: kgexsa(ketmax,netmax)
    integer :: ngexsa(ietmax,jetmax,netmax)
    integer :: ngext(jetmax,netmax)

    integer :: iern1
    integer :: iern2
    integer :: np

    character(len=48) :: uspec(nstmax)
    character(len=24) :: ugexmo(netmax)
    character(len=24) :: uphase(nptmax)
    character(len=8) :: ugexj(jetmax,netmax)

    real(kind=8) :: acflg(nstmax)
    real(kind=8) :: actlg(nstmax)
    real(kind=8) :: affpd(nptmax)
    real(kind=8) :: affsd(nstmax)
    real(kind=8) :: cegexs(ietmax,jetmax,netmax)
    real(kind=8) :: conc(nstmax)
    real(kind=8) :: egexjc(jetmax,netmax)
    real(kind=8) :: egexjf(jetmax,netmax)
    real(kind=8) :: egexpa(netmax)
    real(kind=8) :: egexpc(netmax)
    real(kind=8) :: egexs(ietmax,jetmax,netmax)
    real(kind=8) :: egexw(ketmax,netmax)
    real(kind=8) :: moph(nptmax)
    real(kind=8) :: mosp(nstmax)
    real(kind=8) :: sidrph(nptmax)
    real(kind=8) :: sidrsp(nstmax)
    real(kind=8) :: xbar(nstmax)
    real(kind=8) :: xbarlg(nstmax)
    real(kind=8) :: xgexw(ketmax,netmax)

    real(kind=8) :: tolspf
    real(kind=8) :: wkgwi

    ! Local variable declarations.
    integer :: ie
    integer :: je
    integer :: j2
    integer :: j3
    integer :: j4
    integer :: ke
    integer :: nr1
    integer :: nr2
    integer :: ne
    integer :: ns
    integer :: nss

    integer :: ilnobl

    character(len=24) :: ux24

    real(kind=8) :: aafx
    real(kind=8) :: afx
    real(kind=8) :: cxp
    real(kind=8) :: elx
    real(kind=8) :: ex
    real(kind=8) :: expkg
    real(kind=8) :: mxp
    real(kind=8) :: xlx
    real(kind=8) :: xx

    real(kind=8) :: tlg

    mxp = moph(np)
    j2 = ilnobl(uphase(np))

    if (np.lt.iern1 .or. np.gt. iern2) then
        write (noutpt,1000) uphase(np)(1:j2)
        write (noutpt,1000) uphase(np)(1:j2)
1000 format(/' * Error - (EQLIB/prtgex) Programming error trap:',/7x,"Can't write an output table for ",a," because",/7x,"it isn't a generic ion exchanger.")

        stop
    end if

    ne = np - iern1 + 1

    write (noutpt,1010) uphase(np)(1:j2)
1010 format(/16x,'--- ',a,' ---',/)

    cxp = mxp*wkgwi
    j3 = ilnobl(ugexmo(ne))
    write (noutpt,1020) ugexmo(ne)(1:j3),mxp,cxp
1020 format(/5x,'Exchanger model= ',a,/5x,'Amount of substrate (Z)= ',1pe11.4,' mol',/5x,'Concentration of substrate= ',e11.4)

    if (egexpc(ne) .ne. 0.) then
        ex = egexpc(ne)*mxp
        expkg = ex*wkgwi
        write (noutpt,1024) egexpc(ne),ex,expkg
1024 format(/5x,'Apparent cation exchange capacity= ',1pe11.4,' eq/mol',/5x,'Total exchangeable cations= ',e11.4,' eq',/5x,'Total exchangeable cations= ',e11.4,' eq/kg.H2O')
    end if

    if (egexpa(ne) .ne. 0.) then
        ex = egexpa(ne)*mxp
        expkg = ex*wkgwi
        write (noutpt,1030) egexpa(ne),ex,expkg
1030 format(/5x,'Apparent anion exchange capacity= ',1pe11.4,' eq/mol',/5x,'Total exchangeable anions= ',e11.4,' eq',/5x,'Total exchangeable anions= ',e11.4,' eq/kg.H2O')
    end if

    do je = 1,jgext(ne)
        j4 = ilnobl(ugexj(je,ne))
        ex = egexjc(je,ne)*mxp
        expkg = ex*wkgwi
        write (noutpt,1040) ugexj(je,ne)(1:j4),egexjf(je,ne),egexjc(je,ne),ex,expkg
1040 format(/7x,'Site: ',a,//9x,'Formal exchange capacity= ',1pe11.4,' eq/mol',/9x,'Apparent exchange capacity= ',e11.4,' eq/mol',/9x,'Total exchangeable ions= ',e11.4,' eq',/9x,'Total exchangeable ions= ',e11.4,' eq/kg.H2O',/)

        ! Print equivalents and equivalents/kg.H2O.
        write (noutpt,1070)
1070 format(/13x,'Component',19x,'eq',10x,'eq/kg.H2O',/)

        ns = jern1(je,ne) - 1

        do ie = 1,ngext(je,ne)
            ns = ns + 1
            ex = cegexs(ie,je,ne)*mosp(ns)
            expkg = ex*wkgwi
            nss = ngexsa(ie,je,ne)

            if (nss .le. 0) then
                ux24 = '__ '
            else
                ux24 = uspec(nss)
            end if

            if (nss .gt. 0) then
                write (noutpt,1080) ux24,ex,expkg
            end if

1080 format(10x,a24,3x,1pe11.4,3x,e11.4)
        end do

        ! Print equivalent fractions.
        write (noutpt,1110)
1110 format(//13x,'Component',20x,'e',11x,'Log e',/)

        ns = jern1(je,ne) - 1

        do ie = 1,ngext(je,ne)
            ns = ns + 1
            ex = egexs(ie,je,ne)
            nss = ngexsa(ie,je,ne)

            if (nss .gt. 0) then
                ux24 = uspec(nss)
                elx = tlg(ex)
                write (noutpt,1130) ux24,ex,elx
1130 format(10x,a24,3x,1pe11.4,3x,0pf9.4)
            end if
        end do

        write (noutpt,1090)
1090 format(1x)
    end do

    write (noutpt,1140)
1140 format(/6x,'Whole-Phase Equivalent Fractions',//4x,'Component',20x,'e',11x,'Log e',/)

    do ke = kern1(ne),kern2(ne)
        nss = kgexsa(ke,ne)

        if (nss .gt. 0) then
            ex = egexw(ke,ne)
            elx = tlg(ex)
            write (noutpt,1160) uspec(nss),ex,elx
1160 format(1x,a24,3x,1pe11.4,3x,0pf9.4)
        end if
    end do

    write (noutpt,1170)
1170 format(//6x,'Whole-Phase Mole Fractions',//4x,'Component',20x,'x',11x,'Log x',/)

    do ke = kern1(ne),kern2(ne)
        nss = kgexsa(ke,ne)

        if (nss .gt. 0) then
            xx = xgexw(ke,ne)
            xlx = tlg(xx)
            write (noutpt,1160) uspec(nss),xx,xlx
        end if
    end do

    ! Print mole fractions, activity coefficients, and activities.
    write (noutpt,1180)
1180 format(//6x,'Thermodynamic Summary',//4x,'Component',20x,'x',11x,'Log x',3x,'Log lambda',2x,'Log activity',/)

    do je = 1,jgext(ne)
        nr1 = jern1(je,ne)
        nr2 = jern2(je,ne)

        do ns = nr1,nr2
            if (xbar(ns) .gt. 0.) then
                write (noutpt,1190) uspec(ns),xbar(ns),xbarlg(ns),acflg(ns),actlg(ns)
1190 format(1x,a24,3x,1pe11.4,3(3x,0pf9.4))
            end if
        end do

        write (noutpt,1090)
    end do

    ! Print saturation states and affinities.
    write (noutpt,1200)
1200 format(/4x,'Exchanger',21x,'Log Q/K',9x,'Aff, kcal',4x,'State',/)

    afx = affpd(np)
    aafx = abs(afx)

    if (aafx .le. tolspf) then
        write (noutpt,1210) uphase(np),sidrph(np),affpd(np)
1210 format(1x,a24,2x,2(3x,f13.4),3x,'SATD')
    else if (afx .gt. tolspf) then
        write (noutpt,1220) uphase(np),sidrph(np),affpd(np)
1220 format(1x,a24,2x,2(3x,f13.4),3x,'SSATD')
    else
        write (noutpt,1230) uphase(np),sidrph(np),affpd(np)
1230 format(1x,a24,2x,2(3x,f13.4))
    end if

    write (noutpt,1090)

    do je = 1,jgext(ne)
        nr1 = jern1(je,ne)
        nr2 = jern2(je,ne)

        do ns = nr1,nr2
            if (xbar(ns) .gt. 0.) then
                afx = affsd(ns)
                aafx = abs(afx)

                if (aafx .le. tolspf) then
                    write (noutpt,1240) uspec(ns),sidrsp(ns),affsd(ns)
1240 format(3x,a24,2(3x,f13.4),3x,'SATD')
                else if (afx .gt. tolspf) then
                    write (noutpt,1250) uspec(ns),sidrsp(ns),affsd(ns)
1250 format(3x,a24,2(3x,f13.4),3x,'SSATD')
                else
                    write (noutpt,1260) uspec(ns),sidrsp(ns),affsd(ns)
1260 format(3x,a24,2(3x,f13.4))
                end if
            end if
        end do

        write (noutpt,1090)
    end do
end subroutine prtgex