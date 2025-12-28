subroutine prtsso(acflg,actlg,affpd,affsd,ixrn1,ixrn2,jsol,jsomax,ncmpr,noutpt,np,nptmax,nstmax,nxtmax,sidrph,sidrsp,tolspf,uspec,uphase,uxtype,xbar,xbarlg)
    !! This subroutine prints tables describing the state and properties
    !! of the np-th phase (a solid solution).
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
    integer :: jsomax
    integer :: nptmax
    integer :: nstmax
    integer :: nxtmax

    integer :: noutpt

    integer :: jsol(nxtmax)
    integer :: ncmpr(2,nptmax)
    integer :: ixrn1
    integer :: ixrn2
    integer :: np

    character(len=48) :: uspec(nstmax)
    character(len=32) :: uxtype(jsomax)
    character(len=24) :: uphase(nptmax)

    real(kind=8) :: acflg(nstmax)
    real(kind=8) :: actlg(nstmax)
    real(kind=8) :: affpd(nptmax)
    real(kind=8) :: affsd(nstmax)
    real(kind=8) :: sidrph(nptmax)
    real(kind=8) :: sidrsp(nstmax)
    real(kind=8) :: xbar(nstmax)
    real(kind=8) :: xbarlg(nstmax)

    real(kind=8) :: tolspf

    ! Local variable declarations.
    integer :: j2
    integer :: j3
    integer :: k
    integer :: nr1
    integer :: nr2
    integer :: ns
    integer :: nx

    integer :: ilnobl

    real(kind=8) :: aafx
    real(kind=8) :: afx

    j2 = ilnobl(uphase(np))

    if (np.lt.ixrn1 .or. np.gt. ixrn2) then
        write (noutpt,1000) uphase(np)(1:j2)
        write (noutpt,1000) uphase(np)(1:j2)
1000 format(/' * Error - (EQLIB/prtsso) Programming error trap:',/7x,"Can't write an output table for ",a," because",/7x,"it isn't a solid solution.")

        stop
    end if

    nr1 = ncmpr(1,np)
    nr2 = ncmpr(2,np)

    nx = np - ixrn1 + 1
    k = jsol(nx)

    j3 = ilnobl(uxtype(k))
    write (noutpt,1010) uphase(np)(1:j2),uxtype(k)(1:j3)
1010 format(/16x,'--- ',a,' ---'//3x,a,/)

    ! Print mole fractions, activity coefficients, and activities.
    write (noutpt,1020)
1020 format(4x,'Component',20x,'x',11x,'Log x',3x,'Log lambda',2x,'Log activity',/)

    do ns = nr1,nr2
        if (xbar(ns) .le. 0.) then
            write (noutpt,1030) uspec(ns),xbar(ns)
        else
            write (noutpt,1030) uspec(ns),xbar(ns),xbarlg(ns),acflg(ns),actlg(ns)
1030 format(1x,a24,3x,1pe11.4,3(3x,0pf9.4))
        end if
    end do

    write (noutpt,1040)
1040 format(1x)

    ! Print saturation states and affinities.
    write (noutpt,1100)
1100 format(/4x,'Mineral',23x,'Log Q/K',9x,'Aff, kcal',4x,'State',/)

    afx = affpd(np)
    aafx = abs(afx)

    if (aafx .le. tolspf) then
        write (noutpt,1110) uphase(np),sidrph(np),affpd(np)
1110 format(1x,a24,2x,2(3x,f13.4),3x,'SATD')
    else if (afx .gt. tolspf) then
        write (noutpt,1120) uphase(np),sidrph(np),affpd(np)
1120 format(1x,a24,2x,2(3x,f13.4),3x,'SSATD')
    else
        write (noutpt,1130) uphase(np),sidrph(np),affpd(np)
1130 format(1x,a24,2x,2(3x,f13.4))
    end if

    do ns = nr1,nr2
        if (xbar(ns) .gt. 0.) then
            afx = affsd(ns)
            aafx = abs(afx)

            if (aafx .le. tolspf) then
                write (noutpt,1140) uspec(ns),sidrsp(ns),affsd(ns)
1140 format(3x,a24,2(3x,f13.4),3x,'SATD')
            else if (afx .gt. tolspf) then
                write (noutpt,1150) uspec(ns),sidrsp(ns),affsd(ns)
1150 format(3x,a24,2(3x,f13.4),3x,'SSATD')
            else
                write (noutpt,1160) uspec(ns),sidrsp(ns),affsd(ns)
1160 format(3x,a24,2(3x,f13.4))
            end if
        end if
    end do

    write (noutpt,1070)
1070 format(1x)
end subroutine prtsso