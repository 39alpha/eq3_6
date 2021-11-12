      subroutine chktmx(delxi,dlxmin,dlxtmx,drir0,dxval0,eps100,
     $ iodb,nodbmx,nord,noutpt,nrd1mx,nttyo,prcinf,qdump,qriinf,
     $ rirec0,timemx,time0,time1,tolxst,xi0,xi1,xval0)
c
c     This subroutine checks to see that the requested maximum value of
c     time is not exceeded.
c
c     This subroutine is called by:
c
c       EQ6/path.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c
c     Principal output:
c
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer nodbmx,nrd1mx
c
      integer noutpt,nttyo
c
      integer iodb(nodbmx)
c
      integer nord
c
      logical qdump,qriinf
c
      real*8 drir0(nrd1mx),dxval0(nrd1mx)
c
      real*8 delxi,dlxmin,dlxtmx,eps100,prcinf,rirec0,timemx,time0,
     $ time1,tolxst,xi0,xi1,xval0,xtargv
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer icount,ier,ilsign,n,nordp1
c
      character*48 usearch
      character*24 unam24
c
      real*8 dxp,dxsv,tolsx
c
      real*8 fctrl
c
c-----------------------------------------------------------------------
c
      data usearch /'time reaches the requested maximum value        '/
      data unam24  /'                       '/
c
c-----------------------------------------------------------------------
c
      dxsv = delxi
c
      nordp1 = nord + 1
      dlxtmx = prcinf
      if (timemx .le. 0.) go to 999
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     The search target is the requested limit on the time. The interval
c     for convergence is ((1. - tolxst)*timemx,(1. + tolxst)*timemx)).
c
      ilsign = -1
      xtargv = timemx
      tolsx = tolxst*timemx
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      icount = -1
  100 xi1 = xi0 + delxi
      icount = icount + 1
      if (icount .gt. 10) then
        if (iodb(1).gt.0 .or. iodb(5).gt.0) then
          write (noutpt,1000)
 1000     format(/3x,'A scan for where time reaches the requested',
     $    ' maximum value',/3x,'has failed. Dropping to the minimum',
     $    ' step size.')
        endif
        delxi = dlxmin
        go to 990
      endif
c
c     Estimate the time from a Taylor's series expansion.
c
      if (qriinf) then
        time1 = prcinf
      else
        time1 = time0 + rirec0*delxi
        dxp = delxi
        do n = 1,nord
          dxp = dxp*delxi
          time1 = time1 + ( drir0(n)*dxp )/fctrl(n + 1)
        enddo
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (delxi .gt. dlxmin) then
        if (time0 .gt. (1. - tolxst)*timemx) then
c
c         The time at the base point is too close to the requested
c         maximum value.
c
          delxi = dlxmin
          dlxtmx = delxi
          if (iodb(1).gt.0 .or. iodb(5).gt.0) then
            write (noutpt,1010) time0,timemx,xi0,delxi
 1010       format(/3x,'The base point time of ',1pe11.4,' seconds',
     $      ' is already very close',/3x,'to the requested value of ',
     $      e11.4,' seconds',/3x,'at Xi= ',e11.4,'. Have set delxi',
     $      ' equal to the minimum',/7x,'value of ',e11.4,'.')
            go to 100
          endif
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (delxi .gt. dlxmin) then
        if ((time1 - timemx) .gt. tolsx) then
c
          xval0 = time0
          dxval0(1) = rirec0
          do n = 2,nordp1
            dxval0(n) = drir0(n - 1)
          enddo
c
c         Calling sequence substitutions:
c           nordp1 for nord
c
          call search(delxi,dlxmin,dxval0,eps100,ier,ilsign,iodb,
     $    nodbmx,nordp1,noutpt,nrd1mx,nttyo,tolsx,unam24,usearch,
     $    xtargv,xval0)
c
          dlxtmx = delxi
          if (ier .le. 0) go to 100
          if (ier .ge. 2) then
c
c           Note: if ier = 1, the returned "safe" value of delxi
c           is used.
c
            delxi = dlxmin
          endif
          dlxtmx = delxi
          go to 990
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iodb(1).gt.0 .or. iodb(5).gt.0) then
        if (abs(time1 - timemx) .le. tolsx) then
          write (noutpt,1020) timemx,xi1,delxi
 1020     format(/3x,"Taylor's series predict that the model time is",
     $    ' at the requested',/3x,'maximum value of ',1pe11.4,
     $    ' seconds at Xi= ',e11.4,',',/3x,'delxi= ',e11.4,'.')
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  990 if (dxsv .gt. delxi) qdump = .false.
c
  999 continue
      end
