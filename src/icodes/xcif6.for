      program xcif6
c
c     XCIF6: XCON6 interface software
c     EQ3/6 version 8.0a
c
c     Last revised 09/07/11 by TJW
c
c-----------------------------------------------------------------------
c
c     See the readme.txt file that came with this software for further
c     information, including contact information and references.
c
c
c-----------------------------------------------------------------------
c
c     This is a PC utility program to run XCON6, the EQ6 input file
c     reformatter. It is written in Fortran with Lahey compiler
c     extensions, and utilizes DOS system calls. It is a functional
c     near-equivalent of the Linux/UNIX shell script xcif36, when
c     that script is called as "xcif6".
c
c     This program can be adapted to run under Linux/UNIX. Look for
c     system commands run using "call system".
c
c     This code makes the necessary alterations in the IXCON control
c     file (to which both XCON3 and XCON6 respond) to pass on the
c     desired version level (e.g., '6.0', '7.0', 7.2', or 8.0') and
c     format (compact "W" or menu-style "D").
c
c     The input file or files to be run are specified on the command
c     line following the old version level, the new version level,
c     and the new format ("W" or "D"). Leading pathnames may be
c     specified. The original input files are not preserved.
c
c     To use this code, the XCON6 executable must be in the code
c     directory defined in the environment variable "EQ36CO". The code
c     directory should normally be:
c
c       c:\eq3_6v8.0a\bin
c
c     Substitute "d:" for "c:" if the d: drive is used instead. The
c     EQ36CO variable is normally set by the eq36cfg.bat file.
c
c     This code requires modules from the EPCLIB library.
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer nargpa,ninppa,nlchpa
      parameter (nargpa = 24,ninppa = 150,nlchpa = 80)
c
      integer nscrch,nttyo,nxcon
      integer i,i2,itot,j,jj,j2,k,kk,n,nargmx,nargv,nferr,nfnoex,
     $ nfrun,nfwrex,nlchmx
c
      integer iargc,ilnobl
c
      logical qerr,qex
c
      character(len=80) ucodir,ucodrn,ucomln,uline,ulscr,upathn,uxc6c
      character(len=80) uferr(ninppa),ufnoex(ninppa),
     $ ufrun(ninppa),ufwrex(ninppa)
      character(len=80) uargv(1:nargpa),uinput(1:ninppa)
      character(len=8) ux8
      character(len=3) unewvl,uoldvl,uver
      character(len=1) u1,u2,u3,u4,ufor,unewf
c
c-----------------------------------------------------------------------
c
      data nscrch /14/,nttyo  /6/,nxcon  /8/
c
      data ucodrn /"EQ36CO"/
c
c-----------------------------------------------------------------------
c
c     Set dimensioning variables.
c
      nargmx = nargpa
      nlchmx = nlchpa
c
c-----------------------------------------------------------------------
c
c     Get the number of command line arguments.
c
      nargv = iargc()
      nargv = min(nargv,nargmx)
c
      do n = 1,nargmx
        uargv(n) = " "
      enddo
c
c     Get the command line arguments.
c
      do n = 1,nargv
        call getarg(n,uargv(n))
      enddo
c
c     Check the number of arguments on the command line.
c
      if (nargv .lt. 4) then
        write (nttyo,1000)
 1000   format(' Usage: xcif6 oldversion newversion newformat',
     $  ' inputfile(s)',
     $  /'   oldversion = old version level, "6", "7", "7.2", or "8"',
     $  /'   newversion = new version level, "6", "7", "7.2", or "8"',
     $  /'   newformat = format, "W" (compact) or "D" (menu-style)',
     $  /'   inputfiles must end in .6i')
        go to 999
      endif
c
c-----------------------------------------------------------------------
c
      write (nttyo,1010)
 1010 format(' Running XCIF6',/)
c
c-----------------------------------------------------------------------
c
c     Get and check the code directory.
c
      call getenv(ucodrn,ucodir)
      j = ilnobl(ucodir)
      k = ilnobl(ucodrn)
      if (j .le. 0) then
        write (nttyo,1020) ucodrn(1:k)
 1020   format (/' * Error - The ',a," environment variable isn't set.",
     $  /7x,'This is normally set by the eq36cfg.bat file.',
     $  /7x,'Please refer to the installation instructions.',/)
        go to 999
      endif
c
c     Calling sequence substitutions:
c       ucodir for ufile
c
      call texdir(ucodir,qex)
      if (.not.qex) then
        write (nttyo,1040) ucodir(1:j)
 1040   format (/' * Error - The EQ3/6 code directory ',a,
     $  /7x,"doesn't exist.",/)
        go to 999
      endif
      write (nttyo,1050) ucodir(1:j)
 1050 format('   The code directory is ',a)
c
c     Get the XCON6 executable.
c
      uxc6c = ucodir(1:j) // '\xcon6.exe'
c
c     Check the XCON6 executable.
c
      inquire(file=uxc6c,exist=qex)
      if (.not.qex) then
        write (nttyo,1060) ucodir(1:j),ucodrn(1:k)
 1060   format (/" * Error - The executable code xcon6.exe doesn't",
     $  /7x,'exist in the EQ3/6 code directory ',a,
     $  /7x,'defined in the environment variable ',a,'.',/)
        go to 999
      endif
c
c-----------------------------------------------------------------------
c
c     Get and test the default for the old version level.
c
      uver = uargv(1)
      if (uver(1:3) .eq. '6.0') then
        uoldvl = '6.0'
      elseif (uver(1:3) .eq. '6  ') then
        uoldvl = '6.0'
      elseif (uver(1:3) .eq. '7.0') then
        uoldvl = '7.0'
      elseif (uver(1:3) .eq. '7  ') then
        uoldvl = '7.0'
      elseif (uver(1:3) .eq. '7.2') then
        uoldvl = '7.2'
      elseif (uver(1:3) .eq. '8.0') then
        uoldvl = '8.0'
      elseif (uver(1:3) .eq. '8  ') then
        uoldvl = '8.0'
      else
        j = ilnobl(uver)
        write (nttyo,1061) uver(1:j)
 1061   format(/" * Error - Don't recognize default old version level",
     $  ' "',a,'".',/)
        write (nttyo,1000)
      endif
c
c-----------------------------------------------------------------------
c
c     Get and test the new version level.
c
      uver = uargv(2)
      if (uver(1:3) .eq. '6.0') then
        unewvl = '6.0'
      elseif (uver(1:3) .eq. '6  ') then
        unewvl = '6.0'
      elseif (uver(1:3) .eq. '7.0') then
        unewvl = '7.0'
      elseif (uver(1:3) .eq. '7  ') then
        unewvl = '7.0'
      elseif (uver(1:3) .eq. '7.2') then
        unewvl = '7.2'
      elseif (uver(1:3) .eq. '8.0') then
        unewvl = '8.0'
      elseif (uver(1:3) .eq. '8  ') then
        unewvl = '8.0'
      else
        j = ilnobl(uver)
        write (nttyo,1062) uver(1:j)
 1062   format(/" * Error - Don't recognize new version level",
     $  ' "',a,'".',/)
        write (nttyo,1000)
      endif
c
c-----------------------------------------------------------------------
c
c     Get and test the new format.
c
      ufor = uargv(3)
      if (ufor(1:1) .eq. 'W') then
        unewf = 'W'
      elseif (ufor(1:1) .eq. 'D') then
        unewf = 'D'
      else
        j = ilnobl(ufor)
        write (nttyo,1063) ufor(1:j)
 1063   format(/" * Error - Don't recognize new format",' "',a,'".',/)
        write (nttyo,1000)
      endif
c
c-----------------------------------------------------------------------
c
c     Set up the IXCON control file.
c
      call kilfil('ixcon')
      open(nxcon,file='ixcon',status='new')
      write (nxcon,1065)
 1065 format("  This file contains the user-controlled options for",
     $ /"  the EQ3/6 input file reformatters XCON3 and XCON6.",
     $ /"  The code will scan the old input file to determine",
     $ /"  its format. It will also try to determine its",
     $ /"  version level. If it can't do so, the version level",
     $ /"  must be set using the default in this control file.",/)
c
      u1 = ' '
      u2 = ' '
      u3 = ' '
      u4 = ' '
      if (uoldvl(1:3) .eq. '6.0') then
        u1 = 'X'
      elseif (uoldvl(1:3) .eq. '7.0') then
        u2 = 'X'
      elseif (uoldvl(1:3) .eq. '7.2') then
        u3 = 'X'
      elseif (uoldvl(1:3) .eq. '8.0') then
        u4 = 'X'
      endif
      write (nxcon,1067) u1,u2,u3,u4
 1067 format("OLD INPUT FILE VERSION LEVEL (DEFAULT ONLY)",
     $ /" |",a1,"| 6.0 (versions 6.0-6.1)",
     $ /" |",a1,"| 7.0 (versions 7.0-7.1)",
     $ /" |",a1,"| 7.2 (version 7.2)",
     $ /" |",a1,"| 8.0 (version 8.0)")
c
      u1 = ' '
      u2 = ' '
      if (unewf(1:1) .eq. 'W') then
        u1 = 'X'
      elseif (unewf(1:1) .eq. 'D') then
        u2 = 'X'
      endif
      write (nxcon,1068) u1,u2
 1068 format("NEW INPUT FILE FORMAT",
     $ /" |",a1,"| W   (compact)",
     $ /" |",a1,"| D   (menu-style)")
c
      u1 = ' '
      u2 = ' '
      u3 = ' '
      u4 = ' '
      if (unewvl(1:3) .eq. '6.0') then
        u1 = 'X'
      elseif (unewvl(1:3) .eq. '7.0') then
        u2 = 'X'
      elseif (unewvl(1:3) .eq. '7.2') then
        u3 = 'X'
      elseif (unewvl(1:3) .eq. '8.0') then
        u4 = 'X'
      endif
      write (nxcon,1069) u1,u2,u3,u4
 1069 format("NEW INPUT FILE VERSION LEVEL",
     $ /" |",a1,"| 6.0 (versions 6.0-6.1)",
     $ /" |",a1,"| 7.0 (versions 7.0-7.1)",
     $ /" |",a1,"| 7.2 (version 7.2)",
     $ /" |",a1,"| 8.0 (version 8.0)",/"END")
c
      close(nxcon)
c
c-----------------------------------------------------------------------
c
c     Initialize counters:
c
c       nfrun  = number of input files converted
c       nfwrex = number of input files not converted because they have
c                  the wrong file name extension
c       nferr  = number of input files which were not converted
c                  due to XCON6 errors
c       nfnoex = number of input files which were not converted because
c                  they don't exist
c
      nfrun = 0
      nfwrex = 0
      nferr = 0
      nfnoex = 0
c
c-----------------------------------------------------------------------
c
c     Process the desired input files.
c
c     Loop on the remaining command line arguments. These may be
c     individual files or sets defined using metacharacters.
c
      do n = 4,nargv
c
c       Get the leading pathname in the current argument, if any.
c
        upathn = ' '
        uline = uargv(n)
        k = 0
  122   jj = index(uline,'\')
        if (jj .eq. 0) go to 124
        k = k + jj
        upathn = uargv(n)(1:k)
        ulscr = uline(jj + 1:80)
        uline = ulscr
        go to 122
  124   continue
c
c       Are there any metacharacters? If so, expand the argument
c       by building a list of the implied files.
c
        j = scan(uargv(n),'*?[]')
        if (j .eq. 0) then
          itot = 1
          uinput(1) = uargv(n)
          call locase(uinput(1))
        else
          j = ilnobl(uargv(n))
          call kilfil('e36scr')
c PC
          uline = 'dir /l /b ' // uargv(n)(1:j) // ' > e36scr'
c endPC
c
c Linux
c         uline = 'ls -1 ' // uargv(n)(1:j) // ' > e36scr'
c endLinux
c
          call system(uline)
c
          inquire(file='e36scr',exist=qex)
          if (qex) then
            open(nscrch,file='e36scr',status='old')
            itot = 0
            do i = 1,ninppa + 1
              read (nscrch,1074,end=125) uline
 1074         format(a80)
              if (i .gt. ninppa) then
                write (ux8,'(i5)') ninppa
                call lejust(ux8)
                j2 = ilnobl(ux8)
                write (nttyo,1075) uargv(n)(1:j),ux8(1:j2)
 1075           format(/' * Error - The command line argument ',a,
     $          /7x,'expands to more than the allowed ',a,' files.',/)
                stop
              endif
              itot = i
              call lejust(uline)
              k = ilnobl(upathn)
              jj = ilnobl(uline)
              uinput(i) = upathn(1:k) // uline(1:jj)
            enddo
  125       close(nscrch,status='delete')
          endif
          if (itot .eq. 0) then
            write (nttyo,1076) uargv(n)(1:j)
 1076       format(/' * Warning - The command line argument ',a,
     $      /7x,"doesn't expand to include any files.",/)
            go to 195
          endif
        endif
c
c       Loop over all files implied by the current argument.
c
        do i = 1,itot
          write (nttyo,1078)
 1078     format(/' --------------------------------------------------')
          write (nttyo,1072)
 1072     format(1x)
c
c         Get the name of the current input file.
c
          j = ilnobl(uinput(i))
          write (nttyo,1079) uinput(i)(1:j)
 1079     format('   Processing ',a)
c
c         Does the current input file exist?
c
          if (itot .eq. 1) then
            inquire(file=uinput(i),exist=qex)
            if (.not.qex) then
              write (nttyo,1090) uinput(1)(1:j)
 1090         format (/' * Error - The input file "',a,'"',
     $        /7x,"doesn't exist. It will not be converted.",/)
              nfnoex = nfnoex + 1
              ufnoex(nfnoex) = uinput(i)(1:j)
              go to 190
            endif
          endif
c
c         Does the current input file have the proper file name
c         extension?
c
          k = index(uinput(i),'.6i')
          if (k .eq. 0) then
            write (nttyo,1092) uinput(i)(1:j)
 1092       format (/" * Error - Can't run XCIF6 on input file ",'"',a,
     $      '"',
     $      /7x,"because it doesn't have a .6i file name extension.",/)
            kk = index(uinput(i),'.3i')
            if (kk .gt. 0) then
              write (nttyo,1095)
 1095         format (7x,"It's an EQ3NR input file.")
            endif
            nfwrex = nfwrex + 1
            ufwrex(nfwrex) = uinput(i)(1:j)
            go to 190
          endif
          i2 = ilnobl(uinput(i))
          if ((k + 2) .ne. i2) then
            write (nttyo,1092) uinput(i)(1:j)
            go to 190
          endif
c
          call kilfil('input')
          call kilfil('inputs')
          call kilfil('newin')
c
c         Copy the specifed input file into the current directory as
c         "input".
c
c PC
          uline = 'copy ' // uinput(i)(1:j) // ' input > c:\nul'
c endPC
c
c Linux
c         uline = 'ln -s ' // uinput(i)(1:j) // ' input'
c endLinux
c
          call system(uline)
c
c         Run XCON6.
c
          call system(uxc6c)
c
c         Delete the input file.
c
          call kilfil('input')
c
c         Note: there is no attempt here to delete empty files, as in
c         the UNIX shell script counterpart of this program. This is
c         because DOS treats files a bit differently than UNIX. An empty
c         file does not have a zero byte size.
c
c         Was a newin file written?
c
          inquire(file='newin',exist=qex)
          qerr = .not.qex
c
c         If there were errors, write a message.
c
          if (qerr) then
            write (nttyo,1130)
 1130       format(/" * Error - XCON6 couldn't write a reformatted",
     $      /7x,'version of this input file.',/)
          endif
c
c         Set remaining counters.
c
          j = ilnobl(uinput(i))
          if (.not.qerr) then
            nfrun = nfrun + 1
            ufrun(nfrun) = uinput(i)(1:j)
          else
            nferr = nferr + 1
            uferr(nferr) = uinput(i)(1:j)
          endif
c
          if (.not.qerr) then
c
c           Replace the original input file.
c
            call kilfil(uinput(i))
c
c PC
            uline = 'ren newin ' // uinput(i)(1:j)
c endPC
c
c Linux
c           uline = 'mv newin ' // uinput(i)(1:j)
c endLinux
c
            call system(uline)
c
          else
c
c           Delete the newin file.
c
            call kilfil('newin')
c
          endif
c
  190     continue
        enddo
  195   continue
      enddo
c
c-----------------------------------------------------------------------
c
c     Delete the IXCON control file.
c
      call kilfil('ixcon')
c
c-----------------------------------------------------------------------
c
c     Write some statistics.
c
      if (nfrun .gt. 0) then
        write (nttyo,1078)
        write (nttyo,1140)
 1140   format(/'  The following input files were converted:',/)
        do n = 1,nfrun
          j2 = ilnobl(ufrun(n))
          j2 = min(j2,74)
          write (nttyo,1150) ufrun(n)(1:j2)
 1150     format(4x,a)
        enddo
      endif
c
      if (nfwrex .gt. 0) then
        write (nttyo,1078)
        write (nttyo,1155)
 1155   format(/'  The following input files were not converted',
     $  /'  because they have the wrong file name extension:',/)
        do n = 1,nfwrex
          j2 = ilnobl(ufwrex(n))
          j2 = min(j2,74)
          write (nttyo,1150) ufwrex(n)(1:j2)
        enddo
      endif
c
      if (nferr .gt. 0) then
        write (nttyo,1078)
        write (nttyo,1160)
 1160   format(/'  The following input files were not converted due to',
     $  /'  errors encountered while running XCON6:',/)
        do n = 1,nferr
          j2 = ilnobl(uferr(n))
          j2 = min(j2,74)
          write (nttyo,1150) uferr(n)(1:j2)
        enddo
      endif
c
      if (nfnoex .gt. 0) then
        write (nttyo,1078)
        write (nttyo,1170)
 1170   format(/"  The following input files don't exist:",/)
        do n = 1,nfnoex
          j2 = ilnobl(ufnoex(n))
          j2 = min(j2,74)
          write (nttyo,1150) ufnoex(n)(1:j2)
        enddo
      endif
c
c-----------------------------------------------------------------------
c
  999 write (nttyo,1078)
      write (nttyo,1210)
 1210 format(' All done')
      end
