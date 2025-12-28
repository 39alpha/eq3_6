subroutine chksar(afrc0,afrcp,dafrc0,delxi,dlxmin,dxval0,eps100,iodb,jreac,nodbmx,noutpt,nord,nordmx,nrct,nrctmx,nrd1mx,nrk,nttyo,tolsar,ureac,xi0,xi1,xval0)
    !! This subroutine checks the signs of the affinities of the
    !! reactants. It finds the point of reaction progress at which the
    !! affinity of any irreversible reaction changes sign. The reactant
    !! affinities are tracked using finite differences.
    !! This mechanism is used for two purposes. One is to insure that
    !! a rate law for say the forward direction is not extrapolated
    !! to the region of unfavorable affinity when a distinct rate
    !! law is specified for the opposite direction. The second purpose
    !! of this mechanism is to find the point at which a formerly
    !! exhausted reactant must be reactivated so that it can be
    !! precipitated according to a specified rate law (not a partial
    !! equilibrium condition).
    !! Note the following:
    !!   jreac = reactant status flag:
    !!              0 = set to react
    !!             -1 = saturated, but the remaining reactant mass
    !!                    continues to react irreversibly
    !!              1 = exhausted
    !!              2 = saturated; the status of any remaining reactant
    !!                    mass is changed to that of a product phase
    !! When jreac = -1 or 2, the affinity is zero by definition. In this
    !! case, finite-difference expressions of the affinity should be
    !! avoided because the affinity values used to build them may
    !! reflect somewhat random non-zero values consistent with the
    !! applied convergence tolerances. The corresponding finite
    !! differences may therefore behave wildly.
    !! This subroutine is called by:
    !!   EQ6/path.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: nodbmx
    integer :: nordmx
    integer :: nrctmx
    integer :: nrd1mx

    integer :: noutpt
    integer :: nttyo

    integer :: iodb(nodbmx)
    integer :: jreac(nrctmx)
    integer :: nrk(2,nrctmx)

    integer :: nord
    integer :: nrct

    character(len=24) :: ureac(nrctmx)

    real(kind=8) :: afrc0(nrctmx)
    real(kind=8) :: afrcp(nrctmx)
    real(kind=8) :: dafrc0(nordmx,nrctmx)
    real(kind=8) :: dxval0(nrd1mx)

    real(kind=8) :: delxi
    real(kind=8) :: dlxmin
    real(kind=8) :: eps100
    real(kind=8) :: tolsar
    real(kind=8) :: xi0
    real(kind=8) :: xi1
    real(kind=8) :: xval0

    ! Local variable declarations.
    integer :: icount
    integer :: ier
    integer :: ilsign
    integer :: j2
    integer :: krzero
    integer :: n
    integer :: nrc
    integer :: nrzero

    integer :: ilnobl

    character(len=48) :: usearch
    character(len=24) :: unam24

    real(kind=8) :: aaxxp
    real(kind=8) :: atgsar
    real(kind=8) :: axx
    real(kind=8) :: axx0
    real(kind=8) :: axxp
    real(kind=8) :: dxp
    real(kind=8) :: tolsx
    real(kind=8) :: xtargv
    real(kind=8) :: xval

    real(kind=8) :: fctrl

    data usearch /'a reactant affinity changes sign                '/

    if (nord .le. 0) then
        go to 999
    end if

    ! The search target is not where a reactant affinity is zero. The
    ! magnitude of the search target is atgsar = 0.5*tolsar. The search
    ! tolerance has the same value. Thus, in the case of a reactant
    ! affinity going from positive to negative values, the search target
    ! is -atgsar and the interval for convergence is (-tolsar,0.). In
    ! the case of a reactant affinity going from negative to positive
    ! values, the search target is atgsar, and the interval for
    ! convergence is (0.,tolsar).
    atgsar = 0.5*tolsar
    tolsx = atgsar

    icount = -1
100 continue
    xi1 = xi0 + delxi
    icount = icount + 1

    if (icount .gt. 50) then
        if (iodb(1).gt.0 .or. iodb(5).gt.0) then
            write (noutpt,1000)
1000 format(/3x,'A scan for reactant affinities changing sign',' has failed.',/3x,'Dropping to the minimum step size.')
        end if

        delxi = dlxmin
        go to 999
    end if

    ! Estimate the reactant affinities from Taylor's series expansions.
    ! Avoid using Taylor's series for reactants whose affinities must
    ! be zero.
    do nrc = 1,nrct
        afrcp(nrc) = 0.

        if (jreac(nrc).eq.0 .or. jreac(nrc).eq. 1) then
            axx = afrc0(nrc)
            dxp = 1.

            do n = 1,nord
                dxp = dxp*delxi
                axx = axx + ( dafrc0(n,nrc)/fctrl(n) )*dxp
            end do

            afrcp(nrc) = axx
        end if
    end do

    ! Find any cross-overs. Note that there are two kinds, positive to
    ! negative, and negative to positive.
    xval = 0.
    nrzero = 0
    krzero = 0
    unam24 = ' '

    do nrc = 1,nrct
        if (jreac(nrc).eq.0 .or. jreac(nrc).eq. 1) then
            axx0 = afrc0(nrc)
            axxp = afrcp(nrc)

            if (axx0 .gt. 0.) then
                if (axxp .lt. -tolsar) then
                    ! A case has been found in which the cross over tolerance
                    ! is exceeded for the positive to negative case.
                    if (nrk(2,nrc) .ne. 0) then
                        ! The case is not one in which precipitation is to be
                        ! governed by instantaneous equilibrium. In that case,
                        ! let other controls apply.
                        krzero = krzero + 1
                        aaxxp = -axxp

                        if (aaxxp .gt. xval) then
                            xval = aaxxp
                            nrzero = nrc
                            unam24 = ureac(nrc)
                            ilsign = 1
                            xtargv = -atgsar
                        end if
                    end if
                end if
            else if (axx0 .lt. 0.) then
                if (axxp .gt. tolsar) then
                    ! A case has been found in which the cross over tolerance
                    ! is exceeded for the negative to positive case.
                    if (jreac(nrc) .ne. 1) then
                        ! The case is not one in which any new forward reaction
                        ! is impossible because the reactant is exhausted.
                        ! That case may be ignored.
                        krzero = krzero + 1
                        aaxxp = axxp

                        if (aaxxp .gt. xval) then
                            xval = aaxxp
                            nrzero = nrc
                            unam24 = ureac(nrc)
                            ilsign = -1
                            xtargv = atgsar
                        end if
                    end if
                end if
            end if
        end if
    end do

    if (krzero .gt. 0) then
        if (delxi .gt. dlxmin) then
            if (abs(xval - xtargv) .gt. tolsx) then
                xval0 = afrc0(nrzero)

                do n = 1,nord
                    dxval0(n) = dafrc0(n,nrzero)
                end do

                call search(delxi,dlxmin,dxval0,eps100,ier,ilsign,iodb,nodbmx,nord,noutpt,nrd1mx,nttyo,tolsx,unam24,usearch,xtargv,xval0)

                if (ier .le. 0) then
                    go to 100
                end if

                if (ier .ge. 2) then
                    ! Note: if ier = 1, the returned "safe" value of delxi
                    ! is used.
                    delxi = dlxmin
                end if

                go to 999
            end if
        end if
    end if

    if (iodb(1).gt.0 .or. iodb(5).gt.0) then
        ! Print any occurrences of reactant affinity sign changes.
        if (krzero .gt. 0) then
            do nrc = 1,nrct
                if (jreac(nrc).eq.0 .or. jreac(nrc).eq. 1) then
                    axx0 = afrc0(nrc)
                    axxp = afrcp(nrc)

                    if (abs(axxp) .le. tolsx) then
                        if (axx0.gt.tolsx .and. dafrc0(1,nrc).lt.0.) then
                            if (nrk(2,nrc) .ne. 0) then
                                j2 = ilnobl(ureac(nrc))
                                write (noutpt,1010) ureac(nrc)(1:j2),xi1,delxi
1010 format(/3x,"Taylor's series predict that the",' affinity for reactant',/3x,a,' has decreased to',' nearly zero at Xi= ',1pe12.5,',',/3x,'delxi= ',1pe12.5,'.')
                            end if
                        else if (axx0.lt.tolsx .and. dafrc0(1,nrc).gt.0.) then
                            j2 = ilnobl(ureac(nrc))
                            write (noutpt,1020) ureac(nrc)(1:j2),xi1,delxi
1020 format(/3x,"Taylor's series predict that the",' affinity for reactant',/3x,a,' has increased to',' nearly zero at Xi= ',1pe12.5,',',/3x,'delxi= ',1pe12.5,'.')
                        end if
                    end if
                end if
            end do
        end if
    end if

999 continue
end subroutine chksar