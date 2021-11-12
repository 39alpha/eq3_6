      subroutine cko2pl(delxi,dlxmin,do20,dxo0pl,dxo1pl,dxval0,
     $ eps100,iodb,nodbmx,nord,noutpt,nrd1mx,nttyo,fo2lg0,
     $ fo2lg1,o20plo,o21plo,prcinf,qdump,tolxsu,xi0,xi1,xval0)
c
c     This subroutine checks to see that the next log fO2-based plot
c     point is not exceeded. Because the log fO2 might be decreasing
c     or increasing, two potential target points (o20plo and o21plo)
c     must be addressed.
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
      logical qdump
c
      real*8 do20(nrd1mx),dxval0(nrd1mx)
c
      real*8 delxi,dlxmin,dxo0pl,dxo1pl,eps100,prcinf,fo2lg0,o20plo,
     $ fo2lg1,o21plo,tolxsu,xi0,xi1,xval0,xtargv
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer icount,ier,ilsign,n
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
      data usearch /'log fO2 reaches a plot point value              '/
      data unam24  /'                       '/
c
c-----------------------------------------------------------------------
c
      dxsv = delxi
c
      dxo0pl = prcinf
      dxo1pl = prcinf
      if (nord .le. 0) go to 999
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     The first search target is the currently defined log fO2-based
c     plot point associated with the lesser log fO2 value.
c
      ilsign = +1
      xtargv = o20plo
      tolsx = tolxsu
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      icount = -1
  100 xi1 = xi0 + delxi
      icount = icount + 1
      if (icount .gt. 10) then
        if (iodb(1).gt.0 .or. iodb(5).gt.0) then
          write (noutpt,1000)
 1000     format(/3x,'A scan for where the log fO2 reaches the lesser',
     $    ' of two plot point values',/3x,'has failed. Dropping to',
     $    ' the minimum step size.')
        endif
        delxi = dlxmin
        go to 110
      endif
c
c     Estimate the log fO2 from a Taylor's series expansion.
c
      fo2lg1 = fo2lg0
      dxp = 1.
      do n = 1,nord
        dxp = dxp*delxi
        fo2lg1 = fo2lg1 + ( do20(n)/fctrl(n) )*dxp
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (delxi .gt. dlxmin) then
        if (fo2lg0 .le. (1. + tolxsu)*o20plo) then
c
c         The log fO2 at the base point is too close to the requested
c         lesser plot point value.
c
          delxi = dlxmin
          dxo0pl = delxi
          if (iodb(1).gt.0 .or. iodb(5).gt.0) then
            write (noutpt,1010) fo2lg0,o20plo,xi0,delxi
 1010       format(/3x,'The base point log fO2 of ',1pe11.4,' is',
     $      ' already very close',/3x,'to the requested lesser plot',
     $      ' point value of',e11.4,/3x,'at Xi= ',e11.4,'. Have set',
     $      ' delxi equal to the minimum',/7x,'value of ',e11.4,'.')
            go to 100
          endif
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (delxi .gt. dlxmin) then
        if ((fo2lg1 - o20plo) .lt. -tolsx) then
c
          xval0 = fo2lg0
          do n = 1,nord
            dxval0(n) = do20(n)
          enddo
c
          call search(delxi,dlxmin,dxval0,eps100,ier,ilsign,iodb,
     $    nodbmx,nord,noutpt,nrd1mx,nttyo,tolsx,unam24,usearch,
     $    xtargv,xval0)
c
          dxo0pl = delxi
          if (ier .le. 0) go to 100
          if (ier .ge. 2) then
c
c           Note: if ier = 1, the returned "safe" value of delxi
c           is used.
c
            delxi = dlxmin
          endif
          dxo0pl = delxi
          go to 110
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iodb(1).gt.0 .or. iodb(5).gt.0) then
        if (abs(fo2lg1 - o20plo) .le. tolsx) then
          write (noutpt,1020) o20plo,xi1,delxi
 1020     format(/3x,"Taylor's series predict that the log fO2 is",
     $    ' at the requested',/3x,'lesser plot point value of ',
     $    1pe11.4,' at Xi= ',e11.4,',',/3x,'delxi= ',e11.4,'.')
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  110 continue
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     The second search target is the currently defined log fO2-based
c     plot point associated with the greater log fO2 value.
c
      ilsign = -1
      xtargv = o21plo
      tolsx = tolxsu
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      icount = -1
  200 xi1 = xi0 + delxi
      icount = icount + 1
      if (icount .gt. 10) then
        if (iodb(1).gt.0 .or. iodb(5).gt.0) then
          write (noutpt,1100)
 1100     format(/3x,'A scan for where the log fO2 reaches the greater',
     $    ' of two plot point values',/3x,'has failed. Dropping to',
     $    ' the minimum step size.')
        endif
        delxi = dlxmin
        go to 210
      endif
c
c     Estimate the log fO2 from a Taylor's series expansion.
c
      fo2lg1 = fo2lg0
      dxp = 1.
      do n = 1,nord
        dxp = dxp*delxi
        fo2lg1 = fo2lg1 + ( do20(n)/fctrl(n) )*dxp
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (delxi .gt. dlxmin) then
        if (fo2lg0 .ge. (1. - tolxsu)*o21plo) then
c
c         The log fO2 at the base point is too close to the requested
c         greater plot point value.
c
          delxi = dlxmin
          dxo1pl = delxi
          if (iodb(1).gt.0 .or. iodb(5).gt.0) then
            write (noutpt,1110) fo2lg0,o21plo,xi0,delxi
 1110       format(/3x,'The base point log fO2 of ',1pe11.4,' is',
     $      ' already very close',/3x,'to the requested greater plot',
     $      ' point value of',e11.4,/3x,'at Xi= ',e11.4,'. Have set',
     $      ' delxi equal to the minimum',/7x,'value of ',e11.4,'.')
            go to 200
          endif
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (delxi .gt. dlxmin) then
        if ((fo2lg1 - o21plo) .gt. tolsx) then
c
          xval0 = fo2lg0
          do n = 1,nord
            dxval0(n) = do20(n)
          enddo
c
          call search(delxi,dlxmin,dxval0,eps100,ier,ilsign,iodb,
     $    nodbmx,nord,noutpt,nrd1mx,nttyo,tolsx,unam24,usearch,
     $    xtargv,xval0)
c
          dxo1pl = delxi
          if (ier .le. 0) go to 200
          if (ier .ge. 2) then
c
c           Note: if ier = 1, the returned "safe" value of delxi
c           is used.
c
            delxi = dlxmin
          endif
          dxo1pl = delxi
          go to 210
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iodb(1).gt.0 .or. iodb(5).gt.0) then
        if (abs(fo2lg1 - o21plo) .le. tolsx) then
          write (noutpt,1120) o21plo,xi1,delxi
 1120     format(/3x,"Taylor's series predict that the log fO2 is",
     $    ' at the requested',/3x,'greater plot point value of ',
     $    1pe11.4,' at Xi= ',e11.4,',',/3x,'delxi= ',e11.4,'.')
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  210 continue
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  990 if (dxsv .gt. delxi) qdump = .false.
c
  999 continue
      end
