      subroutine aaaepl()
c
c     This is a dummy subroutine designed to lead-off the source code
c     of the EPCLIB library.
c
c     EPCLIB: EQ3/6 PC interface software library
c     EQ3/6 version 8.0a
c
c     Last revised 09/07/11 by TJW
c
c-----------------------------------------------------------------------
c
c     See the readme.txt file that came with this software for further
c     information, including contact information and references.
c
c-----------------------------------------------------------------------
c
c     The routines in this library are called by each other and by the
c     EQ3/6 interface software consisting of the codes RUNEQPT, RUNEQ3,
c     RUNEQ6, XCIF3, and XCIF6.
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       None
c
c     Principal output:
c
c       None
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
c       None
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
c       None
c
c-----------------------------------------------------------------------
c
      end
      subroutine ckferr(noutpt,qerr)
c
c     This subroutine checks the output file for EQ3/6 error
c     messages. The variable qerr is returned with a value of .true.
c     if any are detected.
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer noutpt
c
      logical qerr
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,j
c
      character(len=80) uline
c
c-----------------------------------------------------------------------
c
      qerr = .false.
      rewind(noutpt)
c
  100 read (noutpt,1000,end=999,err=110) uline
 1000 format(a80)
      i = index(uline,'* error')
      j = index(uline,'* Error')
      if (i .gt. 0) go to 110
      if (j .gt. 0) go to 110
      go to 100
c
  110 qerr = .true.
c
  999 continue
      end
      subroutine ckfwar(noutpt,qwarn)
c
c     This subroutine checks the output file for EQ3/6 warning
c     messages. The variable qwarn is returned with a value of .true.
c     if any are detected.
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer noutpt
c
      logical qwarn
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,j
c
      character(len=80) uline
c
c-----------------------------------------------------------------------
c
      qwarn = .false.
      rewind(noutpt)
c
  100 read (noutpt,1000,end=999,err=999) uline
 1000 format(a80)
      i = index(uline,'* warning')
      j = index(uline,'* Warning')
      if (i .gt. 0) go to 110
      if (j .gt. 0) go to 110
      go to 100
c
  110 qwarn = .true.
c
  999 continue
      end
      subroutine ggammo(iopg1,nfile,qgcoef)
c
c     This subroutine determines the activity coefficient model
c     type specified on an EQ3NR or EQ6 input file.
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer nfile
c
      integer iopg1
c
      logical qgcoef
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i
c
      character(len=80) uline
      character(len=8) uendit,uinfor
c
c-----------------------------------------------------------------------
c
      data uendit/'endit.  '/
c
c-----------------------------------------------------------------------
c
      qgcoef = .false.
      rewind(nfile)
c
c     Determine the problem input format.
c
      read (nfile,1000,end=999,err=999) uline
 1000 format(a80)
      backspace nfile
      uinfor = 'W'
      if (uline(1:8) .eq. '|-------') uinfor = 'D'
c
c     Scan the problem input.
c
      if (uinfor(1:1) .eq. 'W') then
c
c       Compact (W) format.
c
c       Skip past the title.
c
  100   read (nfile,1000,end=999,err=999) uline
        if (uline(1:8). eq. uendit(1:8)) go to 110
        go to 100
c
c       Find the iopg1 option.
c
  110   read (nfile,1000,end=999,err=999) uline
c
        i = index(uline,'  iopg1-10=     0')
        if (i .eq. 1) then
          qgcoef = .true.
          iopg1 = 0
          go to 999
        endif
c
        i = index(uline,'  iopg1-10=    -1')
        if (i .eq. 1) then
          qgcoef = .true.
          iopg1 = -1
          go to 999
        endif
c
        i = index(uline,'  iopg1-10=     1')
        if (i .eq. 1) then
          qgcoef = .true.
          iopg1 = 1
          go to 999
        endif
c
        go to 110
c
      else
c
c       Menu-style (D) format.
c
c       Skip the title.
c
        read (nfile,1000,end=999,err=999) uline
        read (nfile,1000,end=999,err=999) uline
        read (nfile,1000,end=999,err=999) uline
c
  120   read (nfile,1000,end=999,err=999) uline
        if (uline(1:8) .eq. '|-------') go to 130
        go to 120
c
c       Find the iopg1 option.
c
  130   read (nfile,1000,end=999,err=999) uline
c
        i = index(uline,'|  [x] ( 0) The B-dot equation')
        if (i .eq. 1) then
          qgcoef = .true.
          iopg1 = 0
          go to 999
        endif
c
        i = index(uline,'|  [X] ( 0) The B-dot equation')
        if (i .eq. 1) then
          qgcoef = .true.
          iopg1 = 0
          go to 999
        endif
c
        i = index(uline,'|  [*] ( 0) The B-dot equation')
        if (i .eq. 1) then
          qgcoef = .true.
          iopg1 = 0
          go to 999
        endif
c
        i = index(uline,'|  [x] (-1) The Davies equation')
        if (i .eq. 1) then
          qgcoef = .true.
          iopg1 = -1
          go to 999
        endif
c
        i = index(uline,'|  [X] (-1) The Davies equation')
        if (i .eq. 1) then
          qgcoef = .true.
          iopg1 = -1
          go to 999
        endif
c
        i = index(uline,'|  [*] (-1) The Davies equation')
        if (i .eq. 1) then
          qgcoef = .true.
          iopg1 = -1
          go to 999
        endif
c
        i = index(uline,"|  [x] ( 1) Pitzer's equation")
        if (i .eq. 1) then
          qgcoef = .true.
          iopg1 = 1
          go to 999
        endif
c
        i = index(uline,"|  [X] ( 1) Pitzer's equation")
        if (i .eq. 1) then
          qgcoef = .true.
          iopg1 = 1
          go to 999
        endif
c
        i = index(uline,"|  [*] ( 1) Pitzer's equation")
        if (i .eq. 1) then
          qgcoef = .true.
          iopg1 = 1
          go to 999
        endif
c
        go to 130
c
      endif
c
  999 continue
c
      end
      integer function ifnobl(ustr)
c
c     This subroutine finds the position of the first non-blank
c     character in the string ustr.
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       ustr   = the input string variable
c
c     Principal output:
c
c       ifnobl = the position of the first non-blank character
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      character(len=*) ustr
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer j,nchars
c
c-----------------------------------------------------------------------
c
c     Get the length of the string variable.
c
      nchars = len(ustr)
c
c     Find the first non-blank character.
c
      ifnobl = 0
      do j = 1,nchars
        if (ustr(j:j) .ne. ' ') then
          ifnobl = j
          go to 999
        endif
      enddo
c
  999 continue
      end
      integer function ilnobl(ustr)
c
c     This subroutine finds the position of the first non-blank
c     character in the string ustr.
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       ustr   = the input string variable
c
c     Principal output:
c
c       ilnobl = the position of the first non-blank character
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      character(len=*) ustr
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer j,nchars
c
c-----------------------------------------------------------------------
c
c     Get the length of the string variable.
c
      nchars = len(ustr)
c
c     Find the first non-blank character.
c
      ilnobl = 0
      do j = nchars,1,-1
        if (ustr(j:j) .ne. ' ') then
          ilnobl = j
          go to 999
        endif
      enddo
c
  999 continue
      end
      subroutine kilfil(ufile)
c
c     This subroutine kills the file ufile.
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       ufile  = the name of the file to be killed
c
c     Principal output:
c
c       none
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      character(len=*) ufile
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer j2,nlu
c
      integer ilnobl
c
      logical qex,qopen
c
c-----------------------------------------------------------------------
c
      inquire(file=ufile,exist=qex)
      if (qex) then
        do nlu=19,0,-1
          inquire(unit=nlu,opened=qopen)
          if (.not.qopen) go to 100
        enddo
c
        j2 =ilnobl(ufile)
        write (6,1000) ufile(1:j2)
 1000   format(" * Error - (EPCLIB\kilfil) Couldn't find an available",
     $  ' logical',/7x,'unit number to use in killing ',a,'.',/)
        stop
c
  100   continue
        open(nlu,file=ufile,err=110)
        close(nlu,status='delete',err=110)
      endif
  110 continue
c
      inquire(file=ufile,exist=qex)
      if (qex) then
        j2 =ilnobl(ufile)
        j2 =ilnobl(ufile)
        write (6,1010) ufile(1:j2)
 1010   format(' * Error - (EPCLIB\kilfil) The file ',a,' still',
     $  /7x,'exists after attempting to kill it.',/)
        stop
      endif
c
      end
      subroutine lejust(ustr)
c
c     This subroutine left-justifies the non-blank portion of the
c     string ustr.
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       ustr   = the input string variable
c
c     Principal output:
c
c       ustr   = the output string variable
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      character(len=*) ustr
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer j,jj,jbl,j1,nchars
      integer ifnobl
c
c-----------------------------------------------------------------------
c
c     Get the length of the string variable.
c
      nchars = len(ustr)
c
c     Get the position of the first non-blank character and the number
c     of blanks on the left-hand-side.
c
      j1 = ifnobl(ustr)
      jbl = j1 - 1
c
      if (jbl .gt. 0) then
        do jj = j1,nchars
          j = jj - jbl
          ustr(j:j) = ustr(jj:jj)
        enddo
        do j = nchars - jbl,nchars
          ustr(j:j) = ' '
        enddo
      endif
c
      end
      subroutine locase(ustr)
c
c     This subroutine converts a string from upper case to lower case.
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      character(len=*) ustr
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer idel,j,nchars
c
      character(len=1) u1
c
c-----------------------------------------------------------------------
c
      idel = ichar('A') - ichar('a')
      if (idel .ne. 0) then
        nchars = len(ustr)
        do j = 1,nchars
          u1 = ustr(j:j)
          if (u1.ge.'A' .and. u1.le.'Z')
     $    ustr(j:j) = char(ichar(u1) - idel)
        enddo
      endif
c
      end
      subroutine texdir(udir,qex)
c
c     This subroutine tests directory udir to see if it exists. The
c     variable qex is returned with a value of .true. if this is so.
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      logical qex
c
      character(len=*) udir
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
c     None
c
c-----------------------------------------------------------------------
c
c     This test is presently disabled. The Fortran inquire statement
c     doesn't seem to work correctly with directories. The DOS
c     "if exist" construct seems to work in some versions of Windows,
c     but not all.
c
      qex = .true.
c
      end
      subroutine xchstr(nfile,qchstr,uchstr)
c
c     This subroutine determines the presence of a cheat string
c     (uchstr) in the main title of an EQ3NR or EQ6 input file.
c     The presence of a cheat string in a comment line should
c     be ignored.
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer nfile
c
      logical qchstr
c
      character(len=*) uchstr
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,j2
c
      integer ilnobl
c
      character(len=80) uline,ux80
      character(len=8) uendit,uinfor
c
c-----------------------------------------------------------------------
c
      data uendit/'endit.  '/
c
c-----------------------------------------------------------------------
c
      qchstr = .false.
      j2 = ilnobl(uchstr)
      rewind(nfile)
c
c     Determine the problem input format.
c
      read (nfile,1000,end=999,err=999) uline
 1000 format(a80)
      backspace nfile
      uinfor = 'W'
      if (uline(1:8) .eq. '|-------') uinfor = 'D'
c
c     Scan the problem input.
c
      if (uinfor(1:1) .eq. 'W') then
c
c       Compact (W) format.
c
  100   read (nfile,1000,end=999,err=999) uline
c
c       Skip the current line if it is a comment line.
c
        if (uline(1:1) .eq. "*") go to 100
c
        ux80 = uline
        call locase(ux80)
        i = index(ux80,uchstr(1:j2))
        if (i .gt. 0) then
          qchstr = .true.
          go to 999
        endif
        if (uline(1:8). eq. uendit(1:8)) go to 999
        go to 100
c
      else
c
c       Menu-style (D) format.
c
c       Skip the title header block.
c
        read (nfile,1000,end=999,err=999) uline
        read (nfile,1000,end=999,err=999) uline
        read (nfile,1000,end=999,err=999) uline
c
  110   read (nfile,1000,end=999,err=999) uline
        if (uline(1:1) .eq. "*") go to 110
        ux80 = uline
        call locase(ux80)
        i = index(ux80,uchstr(1:j2))
        if (i .gt. 0) then
          qchstr = .true.
          go to 999
        endif
        if (uline(1:8) .eq. '|-------') go to 999
        go to 110
c
      endif
c
  999 continue
c
      end
