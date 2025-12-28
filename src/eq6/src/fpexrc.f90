subroutine fpexrc(delxi,dlxmin,drer0,dxval0,eps100,iodb,jreac,morr,morr0,nodbmx,nord,noutpt,nrct,nrctmx,nrd1mx,nttyo,qdump,rrelr0,ureac,xi0,xi1,xval0)
    !! This subroutine finds the point of reaction progress at which a
    !! reactant becomes exhausted.
    !! This subroutine is called by:
    !!   EQ6/path.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: nodbmx
    integer :: nrctmx
    integer :: nrd1mx

    integer :: noutpt
    integer :: nttyo

    integer :: iodb(nodbmx)
    integer :: jreac(nrctmx)

    integer :: nord
    integer :: nrct

    logical :: qdump

    character(len=24) :: ureac(nrctmx)

    real(kind=8) :: drer0(nrd1mx,nrctmx)
    real(kind=8) :: dxval0(nrd1mx)
    real(kind=8) :: morr(nrctmx)
    real(kind=8) :: morr0(nrctmx)
    real(kind=8) :: rrelr0(nrctmx)

    real(kind=8) :: delxi
    real(kind=8) :: dlxmin
    real(kind=8) :: dlxrct
    real(kind=8) :: dxp
    real(kind=8) :: eps100
    real(kind=8) :: mxx
    real(kind=8) :: xi0
    real(kind=8) :: xi1
    real(kind=8) :: xval0

    ! Local variable declarations.
    integer :: icount
    integer :: ier
    integer :: ilsign
    integer :: j2
    integer :: kexrc
    integer :: n
    integer :: nexrc
    integer :: nordp1
    integer :: nrc

    integer :: ilnobl

    character(len=48) :: usearch
    character(len=24) :: unam24

    real(kind=8) :: dxsv
    real(kind=8) :: tolsx
    real(kind=8) :: xtargv
    real(kind=8) :: xval

    real(kind=8) :: fctrl

    data usearch /'a reactant is exhausted                         '/

    dxsv = delxi

    ! The search target is not where the number of moles remaining is
    ! zero, but -eps100. The interval for convergence is (-2*eps100,0).
    nordp1 = nord + 1
    ilsign = 1
    xtargv = -eps100
    tolsx = eps100

    icount = -1
100 continue
    xi1 = xi0 + delxi
    icount = icount + 1

    if (icount .gt. 50) then
        if (iodb(1).gt.0 .or. iodb(5).gt.0) then
            write (noutpt,1000)
1000 format(/3x,'A scan for reactant masses becoming exhausted',' has failed.',/3x,'Dropping to the minimum step size.')
        end if

        delxi = dlxmin
        go to 990
    end if

    ! Estimate the number of moles remaining of each reactant from
    ! Taylor's series expansions.
    do nrc = 1,nrct
        dlxrct = 0.

        if (jreac(nrc).eq.0 .or. jreac(nrc).eq.-1) then
            dlxrct = rrelr0(nrc)*delxi
            dxp = delxi

            do n = 1,nord
                dxp = dxp*delxi
                dlxrct = dlxrct + ( drer0(n,nrc)/fctrl(n + 1) )*dxp
            end do
        end if

        ! XX     The following assumes that each reactant has a reaction
        ! XX     coefficient of 1.
        morr(nrc) = morr0(nrc) - dlxrct
    end do

    ! Find any cases of exhausted reactants.
    xval = 1.0
    nexrc = 0
    kexrc = 0
    unam24 = ' '

    do nrc = 1,nrct
        mxx = morr(nrc)

        if (mxx .le. 0.) then
            if (jreac(nrc).eq.0 .or. jreac(nrc).eq.-1) then
                kexrc = kexrc + 1

                if (mxx .lt. xval) then
                    xval = morr(nrc)
                    nexrc = nrc
                    unam24 = ureac(nrc)
                end if
            end if
        end if
    end do

    if (kexrc .gt. 0) then
        if (delxi .gt. dlxmin) then
            if (abs(xval - xtargv) .gt. tolsx) then
                xval0 = morr0(nexrc)

                ! XX         The following assumes that each reactant has a reaction
                ! XX         coefficient of 1.
                dxval0(1) = -rrelr0(nexrc)

                do n = 2,nordp1
                    dxval0(n) = -drer0(n - 1,nexrc)
                end do

                ! Calling sequence substitutions:
                !   nordp1 for nord
                call search(delxi,dlxmin,dxval0,eps100,ier,ilsign,iodb,nodbmx,nordp1,noutpt,nrd1mx,nttyo,tolsx,unam24,usearch,xtargv,xval0)

                if (ier .le. 0) then
                    go to 100
                end if

                if (ier .ge. 2) then
                    ! Note: if ier = 1, the returned "safe" value of delxi
                    ! is used.
                    delxi = dlxmin
                end if

                go to 990
            end if
        end if
    end if

    if (iodb(1).gt.0 .or. iodb(5).gt.0) then
        if (kexrc .gt. 0) then
            write (noutpt,1010) xi1,delxi
1010 format(/3x,"Taylor's series predict exhaustion of the",' following',/3x,'reactants at Xi= ',1pe11.4,', delxi= ',e11.4,':',/)

            do nrc = 1,nrct
                mxx = morr(nrc)

                if (abs(mxx - xtargv) .le. tolsx) then
                    if (jreac(nrc).eq.0 .or. jreac(nrc).eq.-1) then
                        j2 = ilnobl(ureac(nexrc))
                        write (noutpt,1020) ureac(nexrc)(1:j2)
1020 format(5x,a)
                    end if
                end if
            end do
        end if
    end if

990 continue
    if (dxsv .gt. delxi) then
        qdump = .false.
    end if

999 continue
end subroutine fpexrc