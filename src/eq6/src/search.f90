subroutine search(delxi,dlxmin,dxval0,eps100,ier,ilsign,iodb,nodbmx,nord,noutpt,nrd1mx,nttyo,tolsx,unam24,usearch,xtargv,xval0)
    !! This subroutine finds the value of delxi at which occurs an event
    !! whose type is described by the string in the usearch variable.
    !! This subroutine is called by:
    !!   EQ6/chksar.f
    !!   EQ6/chksir.f
    !!   EQ6/chksrr.f
    !!   EQ6/chktmx.f
    !!   EQ6/chktpl.f
    !!   EQ6/chktpr.f
    !!   EQ6/ckawmn.f
    !!   EQ6/ckawmx.f
    !!   EQ6/ckehmn.f
    !!   EQ6/ckehmx.f
    !!   EQ6/cko2mn.f
    !!   EQ6/cko2mx.f
    !!   EQ6/ckphmn.f
    !!   EQ6/ckphmx.f
    !!   EQ6/fpbdpp.f
    !!   EQ6/fpbflo.f
    !!   EQ6/fpbnpp.f
    !!   EQ6/fpexrc.f
    !! Principal input:
    !!   unam24 = name of entity, one of whose properties is involved in
    !!              the search
    !!   ilsign = expected sign of the residual function at the
    !!              left boundary (delxi = 0)
    !!   tolsx  = the convergence tolerance
    !! Principal output:
    !!   ier     = error flag
    !!               = 1  the search did not converge
    !!               = 2  event being searched for appears to have
    !!                      occurred prior to the interval being
    !!                      searched
    !!               = 3  event being searched for doesn't appear to
    !!                      occur in the interval being searched
    implicit none

    ! Calling sequence variable declarations.
    integer :: nodbmx
    integer :: nrd1mx

    integer :: noutpt
    integer :: nttyo

    integer :: iodb(nodbmx)

    integer :: ier
    integer :: ilsign
    integer :: nord

    character(len=48) :: usearch
    character(len=24) :: unam24

    real(kind=8) :: dxval0(nrd1mx)

    real(kind=8) :: delxi
    real(kind=8) :: dlxmin
    real(kind=8) :: eps100
    real(kind=8) :: tolsx
    real(kind=8) :: xtargv
    real(kind=8) :: xval0

    ! Local variable declarations.
    integer :: iter
    integer :: j2
    integer :: j3

    integer :: ilnobl

    real(kind=8) :: x(2)
    real(kind=8) :: y(2)

    real(kind=8) :: adx
    real(kind=8) :: ares
    real(kind=8) :: ares0
    real(kind=8) :: resfnc
    real(kind=8) :: resx
    real(kind=8) :: slope
    real(kind=8) :: xleft
    real(kind=8) :: xnew
    real(kind=8) :: xright

    ier = 0

    ! The following local variables are signficant:
    !   iter   = iteration counter
    !   resfnc = convergence function
    iter = 0
    resfnc = 0.

    j2 = ilnobl(usearch)
    j3 = ilnobl(unam24)

    if (iodb(7) .ge. 1) then
        write (noutpt,1000) usearch(1:j2)
1000 format(/3x,'--- Search for where ',a,' ---',/)

        if (j3 .gt. 0) then
            write (noutpt,1010) unam24(1:j3)
        end if

1010 format(9x,'The entity involved is ',a,'.',/)
    end if

    ! Save the entering value of delxi in xright.
    !   xleft  = the left boundary of the interval being examined
    !   xright = the right boundary
    xleft = 0.
    xright = delxi

    ! Check the left boundary point.
    !   The current left point is ( x(1),y(1) ).
    !   The current right point is ( x(2),y(2) ).
    x(1) = xleft
    delxi = xleft

    call sfncge(delxi,xval0,xtargv,dxval0,nord,nrd1mx,resx)

    y(1) = resx
    ares0 = abs(resx)

    if (iodb(7) .ge. 1) then
        write (noutpt,1020) iter,x(1),resx,resfnc
    end if

1020 format(3x,'iter= ',i2,' delxi= ',1pe11.4,' residual= ',e11.4,' resfnc=',e11.4)

    ! Check for a direct hit on the left boundary.
    if (ares0 .le. 0.) then
        delxi = dlxmin
        go to 999
    end if

    if ( (resx.lt.0 .and. ilsign.gt.0) .or.  (resx.gt.0 .and. ilsign.lt.0) ) then
        write (noutpt,1030) usearch(1:j2)
        write (nttyo,1030) usearch(1:j2)
1030 format(/' * Note - (EQ6/search) A search for where ',a,/7x,' indicates that the event has been stepped over.')

        if (j3 .gt. 0) then
            write (noutpt,1040) unam24(1:j3)
            write (nttyo,1040) unam24(1:j3)
1040 format(7x,'The entity involved is ',a,'.')
        end if

        ier = 2
        delxi = dlxmin
        write (noutpt,1050) delxi
        write (nttyo,1050) delxi
1050 format(7x,'The step size will be set to ',1pe11.4,'.')

        go to 999
    end if

    ! Check the right boundary point.
100 continue
    x(2) = xright
    delxi = xright

    call sfncge(delxi,xval0,xtargv,dxval0,nord,nrd1mx,resx)

    y(2) = resx
    ares = abs(resx)

    if (iodb(7) .ge. 1) then
        write (noutpt,1020) iter,x(2),resx,resfnc
    end if

    ! Check for a direct hit on the right boundary.
    if (ares .le. 0.) then
        delxi = xright
        go to 999
    end if

    ares = min(ares0,ares)

    if ( (resx.lt.0 .and. ilsign.lt.0) .or.  (resx.gt.0 .and. ilsign.gt.0) ) then
        write (noutpt,1120) usearch(1:j2)
        write (nttyo,1120) usearch(1:j2)
1120 format(/' * Note - (EQ6/search) A search for where ',a,/7x,'indicates that the event being searched for does',/7x,'not take place in the interval being examined.')

        if (j3 .gt. 0) then
            write (noutpt,1040) unam24(1:j3)
            write (nttyo,1040) unam24(1:j3)
        end if

        write (noutpt,1130)
        write (nttyo,1130)
1130 format(7x,'The step size will not be decreased.')

        ier = 3
        delxi = xright
        go to 999
    end if

    ! Find a new point. Each new point will replace either the
    ! left point or the right point, depending on the value of
    ! the residual. At least one iteration is required. The
    ! secant method is the principal algorithm for generating a
    ! new point. However, if convergence is slow, interval halving
    ! may occasionally be used instead.
110 continue
    iter = iter + 1
    ares0 = ares

    if (iodb(7) .ge. 2) then
        write (noutpt,1140) x(1),x(2)
    end if

1140 format(7x,'Left point= ',1pe11.4,', Right point= ',e11.4)

    if (iter.le.4 .or. resfnc.ge.0.5 .or. mod(iter,4).ne.0) then
        ! Secant method.
        if (iodb(7) .ge. 1) then
            write (noutpt,1150)
        end if

1150 format(3x,'Secant method')

        slope =  (y(2) - y(1) )/( x(2) - x(1) )
        xnew = x(1) - (y(1)/slope)
    else
        ! Interval halving method.
        if (iodb(7) .ge. 1) then
            write (noutpt,1160)
        end if

1160 format(3x,'Interval halving')

        xnew = 0.5*( x(1) + x(2) )
    end if

    delxi = xnew

    call sfncge(delxi,xval0,xtargv,dxval0,nord,nrd1mx,resx)

    ares = abs(resx)
    resfnc = (ares0 - ares)/ares0

    if (iodb(7) .ge. 1) then
        write (noutpt,1020) iter,xnew,resx,resfnc
    end if

    ! Note there are multiple convergence tests.
    if (ares .le. tolsx) then
        go to 999
    end if

    if (xnew .gt. 0.) then
        ! Check for underflow of the most recent correction.
        adx = abs((xnew - x(1))/xnew)

        if (adx .le. eps100) then
            go to 999
        end if
    end if

    if (abs(xnew - x(1)) .le. dlxmin) then
        ! Check for correction within the limit of the minimum
        ! step size.
        delxi = max(delxi,dlxmin)
        go to 999
    end if

    if (iter .ge. 50) then
        delxi = max(x(1),dlxmin)
        write (noutpt,1170) usearch(1:j2)
        write (nttyo,1170) usearch(1:j2)
1170 format(/' * Note - (EQ6/search) A search for where ',a,/7x,'has failed to converge within the maximum number of',' iterations.')

        if (j3 .gt. 0) then
            write (noutpt,1040) unam24(1:j3)
            write (nttyo,1040) unam24(1:j3)
        end if

        write (noutpt,1180) delxi
        write (nttyo,1180) delxi
1180 format(7x,'The step size will be set to a "safe" value of ',1pe11.4)

        ier = 1
        go to 999
    end if

    if ( (resx.lt.0 .and. ilsign.lt.0 ) .or.  (resx.gt.0 .and. ilsign.gt.0) ) then
        x(1) = xnew
        y(1) = resx
    else
        x(2) = xnew
        y(2) = resx
    end if

    go to 110

999 continue
end subroutine search