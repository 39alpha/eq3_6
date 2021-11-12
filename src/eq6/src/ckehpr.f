      subroutine ckehpr(delxi,dlxmin,deh0,dxe0pr,dxe1pr,dxval0,
     $ eps100,iodb,nodbmx,nord,noutpt,nrd1mx,nttyo,eh0,eh1,
     $ eh0prn,eh1prn,prcinf,qdump,tolxsu,xi0,xi1,xval0)
c
c     This subroutine checks to see that the next Eh-based print
c     point is not exceeded. Because the Eh might be decreasing or
c     increasing, two potential target points (eh0prn and eh1prn)
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
      real*8 deh0(nrd1mx),dxval0(nrd1mx)
c
      real*8 delxi,dlxmin,dxe0pr,dxe1pr,eps100,prcinf,eh0,eh0prn,eh1,
     $ eh1prn,tolxsu,xi0,xi1,xval0,xtargv
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
      data usearch /'Eh reaches a print point value                  '/
      data unam24  /'                       '/
c
c-----------------------------------------------------------------------
c
      dxsv = delxi
c
      dxe0pr = prcinf
      dxe1pr = prcinf
      if (nord .le. 0) go to 999
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     The first search target is the currently defined Eh-based print
c     point associated with the lesser Eh value.
c
      ilsign = +1
      xtargv = eh0prn
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
 1000     format(/3x,'A scan for where the Eh reaches the lesser of',
     $    ' two print point values',/3x,'has failed. Dropping to the',
     $    ' minimum step size.')
        endif
        delxi = dlxmin
        go to 110
      endif
c
c     Estimate the Eh from a Taylor's series expansion.
c
      eh1 = eh0
      dxp = 1.
      do n = 1,nord
        dxp = dxp*delxi
        eh1 = eh1 + ( deh0(n)/fctrl(n) )*dxp
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (delxi .gt. dlxmin) then
        if (eh0 .le. (1. + tolxsu)*eh0prn) then
c
c         The Eh at the base point is too close to the requested
c         lesser print point value.
c
          delxi = dlxmin
          dxe0pr = delxi
          if (iodb(1).gt.0 .or. iodb(5).gt.0) then
            write (noutpt,1010) eh0,eh0prn,xi0,delxi
 1010       format(/3x,'The base point Eh of ',1pe11.4,' v is already',
     $      ' very close',/3x,'to the requested lesser print point',
     $      ' value of',e11.4,/3x,'at Xi= ',e11.4,'. Have set delxi',
     $      ' equal to the minimum',/7x,'value of ',e11.4,'.')
            go to 100
          endif
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (delxi .gt. dlxmin) then
        if ((eh1 - eh0prn) .lt. -tolsx) then
c
          xval0 = eh0
          do n = 1,nord
            dxval0(n) = deh0(n)
          enddo
c
          call search(delxi,dlxmin,dxval0,eps100,ier,ilsign,iodb,
     $    nodbmx,nord,noutpt,nrd1mx,nttyo,tolsx,unam24,usearch,
     $    xtargv,xval0)
c
          dxe0pr = delxi
          if (ier .le. 0) go to 100
          if (ier .ge. 2) then
c
c           Note: if ier = 1, the returned "safe" value of delxi
c           is used.
c
            delxi = dlxmin
          endif
          dxe0pr = delxi
          go to 110
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iodb(1).gt.0 .or. iodb(5).gt.0) then
        if (abs(eh1 - eh0prn) .le. tolsx) then
          write (noutpt,1020) eh0prn,xi1,delxi
 1020     format(/3x,"Taylor's series predict that the Eh is",
     $    ' at the requested',/3x,'lesser print point value of ',
     $    1pe11.4,' v at Xi= ',e11.4,',',/3x,'delxi= ',e11.4,'.')
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  110 continue
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     The second search target is the currently defined Eh-based print
c     point associated with the greater Eh value.
c
      ilsign = -1
      xtargv = eh1prn
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
 1100     format(/3x,'A scan for where the Eh reaches the greater of',
     $    ' two print point values',/3x,'has failed. Dropping to the',
     $    ' minimum step size.')
        endif
        delxi = dlxmin
        go to 210
      endif
c
c     Estimate the Eh from a Taylor's series expansion.
c
      eh1 = eh0
      dxp = 1.
      do n = 1,nord
        dxp = dxp*delxi
        eh1 = eh1 + ( deh0(n)/fctrl(n) )*dxp
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (delxi .gt. dlxmin) then
        if (eh0 .ge. (1. - tolxsu)*eh1prn) then
c
c         The Eh at the base point is too close to the requested
c         greater print point value.
c
          delxi = dlxmin
          dxe1pr = delxi
          if (iodb(1).gt.0 .or. iodb(5).gt.0) then
            write (noutpt,1110) eh0,eh1prn,xi0,delxi
 1110       format(/3x,'The base point Eh of ',1pe11.4,' v is already',
     $      ' very close',/3x,'to the requested greater print point',
     $      ' value of',e11.4,/3x,'at Xi= ',e11.4,'. Have set delxi',
     $      ' equal to the minimum',/7x,'value of ',e11.4,'.')
            go to 200
          endif
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (delxi .gt. dlxmin) then
        if ((eh1 - eh1prn) .gt. tolsx) then
c
          xval0 = eh0
          do n = 1,nord
            dxval0(n) = deh0(n)
          enddo
c
          call search(delxi,dlxmin,dxval0,eps100,ier,ilsign,iodb,
     $    nodbmx,nord,noutpt,nrd1mx,nttyo,tolsx,unam24,usearch,
     $    xtargv,xval0)
c
          dxe1pr = delxi
          if (ier .le. 0) go to 200
          if (ier .ge. 2) then
c
c           Note: if ier = 1, the returned "safe" value of delxi
c           is used.
c
            delxi = dlxmin
          endif
          dxe1pr = delxi
          go to 210
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iodb(1).gt.0 .or. iodb(5).gt.0) then
        if (abs(eh1 - eh1prn) .le. tolsx) then
          write (noutpt,1120) eh1prn,xi1,delxi
 1120     format(/3x,"Taylor's series predict that the Eh is",
     $    ' at the requested',/3x,'greater print point value of ',
     $    1pe11.4,' v at Xi= ',e11.4,',',/3x,'delxi= ',e11.4,'.')
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
