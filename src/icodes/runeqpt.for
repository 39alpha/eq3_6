      program runeqpt
c
c     RUNEQPT: EQPT interface software
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
c     This is a PC utility program to run EQPT, the EQ3/6 data file
c     preprocessor. It is written in Fortran with Lahey compiler
c     extensions, and utilizes DOS system calls. It is a functional
c     near-equivalent of the Linux/UNIX shell script runeqpt.
c
c     This program can be adapted to run under Linux/UNIX. Look for
c     system commands run using "call system".
c
c     To use this code, the data0 files must be in the current
c     directory. The EQPT executable must be in the code directory
c     defined in the environment variable "EQ36CO". The code directory
c     should normally be:
c
c       c:\eq3_6v8.0a\bin
c
c     Substitute "d:" for "c:" if the d: drive is used instead. The
c     EQ36CO variable is normally set by the eq36cfg.bat file.
c
c     The files produced by EQPT are renamed after the parent data0
c     file and appear in the current directory.
c
c     This code requires modules from the EPCLIB library.
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer nargpa,nlchpa,ndafpa
      parameter (nargpa = 24,nlchpa = 80,ndafpa = 100)
c
      integer noutpt,nscrch,nttyo
      integer j2,j3,n,nargmx,nargv,nerr,nf,nferr,nflist,nfnoex,nfrun,
     $ nfwarn,nlchmx,ntot
c
      integer iargc,ilnobl
c
      logical qerr,qex,qwarn
c
      character(len=80) ucodir,ucomln,udata0,ueqptc,ufile,uline,
     $ ulscr,ucodrn
      character(len=80) ufrun(ndafpa),ufwarn(ndafpa),uferr(ndafpa),
     $ ufnoex(ndafpa)
      character(len=80) uargv(1:nargpa)
      character(len=24) uflist(1:8)
      character(len=8) ukey(ndafpa)
      character(len=8) ux8
c
c-----------------------------------------------------------------------
c
      data noutpt /10/,nscrch /14/,nttyo  /6/
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
      if (nargv .lt. 1) then
        write (nttyo,1000)
 1000   format(' Usage: runeqpt datafilekey(s)',
     $  /'   datafilekeys = cmp, sup, 1kb, hmw, ypf, and so forth'
     $  /'   "all" matches all represented in the current directory')
        go to 999
      endif
c
      j2 = ilnobl(uargv(1))
      if (uargv(1)(1:j2).eq.'help' .or. uargv(1)(1:j2).eq.'HELP') then
        write (nttyo,1000)
        go to 999
      endif
c
c-----------------------------------------------------------------------
c
      write (nttyo,1010)
 1010 format(' Running RUNEQPT',/)
c
c     Get and check the code directory.
c
      call getenv(ucodrn,ucodir)
      j2 = ilnobl(ucodir)
      j3 = ilnobl(ucodrn)
c
      if (j2 .le. 0) then
        write (nttyo,1020) ucodrn(1:j3)
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
        write (nttyo,1040) ucodir(1:j2)
 1040   format (/' * Error - The EQ3/6 code directory ',a,
     $  /7x,"doesn't exist.",/)
        go to 999
      endif
      write (nttyo,1050) ucodir(1:j2)
 1050 format('   The code directory is ',a)
c
c     Get the EQPT executable.
c
      ueqptc = ucodir(1:j2) // '\eqpt.exe'
c
c     Check the EQPT executable.
c
      inquire(file=ueqptc,exist=qex)
      if (.not.qex) then
        write (nttyo,1060) ucodir(1:j2),ucodrn(1:j3)
 1060   format (/" * Error - The executable code eqpt.exe doesn't",
     $  /7x,'exist in the EQ3/6 code directory ',a,
     $  /7x,'defined in the environment variable ',a,'.',/)
        go to 999
      endif
c
c-----------------------------------------------------------------------
c
c     Set up the list of keys representing data files to be processed.
c
      if (uargv(1)(1:4).ne.'all ' .and. uargv(1)(1:4).ne.'ALL ') then
c
c       Individual keys are listed.
c
        ntot = nargv
        if (ntot .gt. ndafpa) then
          write (ux8,'(i5)') ndafpa
          call lejust(ux8)
          j2 = ilnobl(ux8)
          write (nttyo,1150) ux8(1:j2)
 1150     format(/' * Error - The command line contains data file',
     $    /7x,'keys for more than the allowed ',a,' files.',/)
          stop
        endif
        do n = 1,nargv
          ukey(n) = uargv(n)
        enddo
c
c       Test the individual keys specified.
c
        nerr = 0
        do n = 1,nargv
          j2 = ilnobl(ukey(n))
          if (j2 .gt. 3) then
            nerr = nerr + 1
            write (nttyo,1152) ukey(n)(1:j2)
 1152       format(/' * Error - The specified data file key "',a,'"',
     $      /7x,'is not valid. Use a three-letter string.',/)
          endif
        enddo
c
        if (nerr .gt. 0) stop
c
      else
c
c       The "all" option was specified.
c
        if (nargv .gt. 1) then
          write (nttyo,1160)
 1160     format (/' * Error - Extra arguments are present following',
     $    /7x,'the "all" option. No processing will be done.',/)
          write (nttyo,1000)
          go to 999
        endif
c
        call kilfil('e36scr')
c
c PC
        uline = 'dir /l /b ' // 'data0.* ' // '> e36scr'
c endPC
c
c Linux
c       uline = 'ls -1 ' // 'data0.* ' // '> e36scr'
c endLinux
c
        call system(uline)
c
        inquire(file='e36scr',exist=qex)
        if (qex) then
          open(nscrch,file='e36scr',status='old')
          ntot = 0
          do n = 1,ndafpa + 1
            read (nscrch,1170,end=110) uline
 1170       format(a80)
            if (n .gt. ndafpa) then
              write (ux8,'(i5)') ndafpa
              call lejust(ux8)
              j2 = ilnobl(ux8)
              write (nttyo,1180) ux8(1:j2)
 1180         format(/' * Error - The "all" command line option',
     $        /7x,'expands to more than the allowed ',a,' files.',/)
              stop
            endif
            ntot = n
            call lejust(uline)
            ukey(n) = uline(7:9)
          enddo
  110     close(nscrch,status='delete')
        else
          write (nttyo,1190)
 1190     format(/" * Error - Can't find any data0 files in the",
     $    /7x,'current directory for the "all" option.',/)
          stop
        endif
      endif
c
c-----------------------------------------------------------------------
c
c     Delete the generic file that EQPT reads, if it exists.
c
      call kilfil('data0')
c
c     Delete the generic files that EQPT writes, if any copies already
c     exist.
c
      call kilfil('data1')
      call kilfil('data1f')
      call kilfil('output')
      call kilfil('slist')
c
c-----------------------------------------------------------------------
c
c     Initialize counters:
c
c       nfrun  = number of data files successfully processed with
c                  no errors or warnings
c       nfwarn = number of data files successfully processed with
c                  no errors, but with warnings
c       nferr  = number of data files that could not be successfully
c                  processed due to errors
c       nfnoex = number of data files which could not be successfully
c                  processed because they don't exist
c
      nfrun = 0
      nfwarn = 0
      nferr = 0
      nfnoex = 0
c
c-----------------------------------------------------------------------
c
c     Process the desired data files.
c
      do n = 1,ntot
        write (nttyo,1200)
 1200   format(/' --------------------------------------------------')
        write (nttyo,1210)
 1210   format(1x)
c
c       Get the name of the current data0 file.
c
  120   udata0 = 'data0.' // ukey(n)(1:3)
        j2 = ilnobl(udata0)
        write (nttyo,1220) udata0(1:j2)
 1220   format('   Processing ',a)
c
c       Does the current data0 file exist?
c
        j2 = ilnobl(udata0)
c
        inquire(file=udata0,exist=qex)
        if (.not.qex) then
          write (nttyo,1230) udata0(1:j2)
 1230     format (/' * Error - The data file "',a,'"',
     $    /7x,"doesn't exist. It will not be processed.",/)
          nfnoex = nfnoex + 1
          ufnoex(nfnoex) = udata0(1:j2)
          go to 390
        endif
c
c       Copy the specifed data0 file into the current directory as
c       "data0".
c
c PC
        uline = 'copy ' // udata0(1:j2) // ' data0 > c:\nul'
c endPC
c
c Linux
c       uline = 'ln -s ' // udata0(1:j2)
c endLinux
c
        call system(uline)
c
c       Delete any old copies of files to be created.
c
        ufile = 'data1.' // ukey(n)(1:3)
        call kilfil(ufile)
        ufile = 'data1f.' // ukey(n)(1:3)
        call kilfil(ufile)
        ufile = 'output.' // ukey(n)(1:3)
        call kilfil(ufile)
        ufile = 'slist.' // ukey(n)(1:3)
        call kilfil(ufile)
c
c       Run EQPT.
c
        call system(ueqptc)
c
c       Delete the data0 file (the copy of the data0 file now being
c       processed).
c
        call kilfil('data0')
c
c       Note: there is no attempt here to delete empty files, as in
c       the UNIX shell script counterpart of this program. This is
c       because DOS treats files a bit differently than UNIX. An empty
c       file does not have a zero byte size.
c
c-----------------------------------------------------------------------
c
c       Check for errors and warnings in processing the current
c       data file. Warnings will be checked for only if no errors
c       are found.
c
        qerr = .false.
        qwarn = .false.
c
        inquire(file='output',exist=qex)
        if (.not.qex) qerr = .true.
        inquire(file='data1',exist=qex)
        if (.not.qex) qerr = .true.
c
        if (.not.qerr) then
          open(noutpt,file='output',status='old')
          call ckferr(noutpt,qerr)
          if (.not.qerr) call ckfwar(noutpt,qwarn)
          close(noutpt)
        endif
c
c-----------------------------------------------------------------------
c
c       If there were errors, destroy the data1 file that was produced.
c
        if (qerr) then
          call kilfil('data1')
          write (nttyo,1270)
 1270     format(/' * Note - EQPT errors were encountered while',
     $    ' processing the',/7x,'current DATA0 file. The DATA1 file',
     $    ' has been deleted',/7x,'to prevent possible misuse.',/)
        endif
c
c-----------------------------------------------------------------------
c
c       Set remaining counters.
c
        j2 = ilnobl(udata0)
        if (.not.qerr) then
          if (.not.qwarn) then
            nfrun = nfrun + 1
            ufrun(nfrun) = udata0(1:j2)
          else
            nfwarn = nfwarn + 1
            ufwarn(nfwarn) = udata0(1:j2)
          endif
        else
          nferr = nferr + 1
          uferr(nferr) = udata0(1:j2)
        endif
c
c       Rename the data1 file.
c
        nflist = 0
        inquire(file='data1',exist=qex)
        if (qex) then
c
c PC
          uline = 'ren data1 data1.' // ukey(n)(1:3)
c endPC
c
c Linux
c         uline = 'mv data1 data1.' // ukey(n)(1:3)
c endLinux
c
          call system(uline)
c
          nflist = nflist + 1
          uflist(nflist) = 'data1.' // ukey(n)(1:3)
        endif
c
c       Rename the other output files.
c
        inquire(file='data1f',exist=qex)
        if (qex) then
c
c PC
          uline = 'ren data1f data1f.' // ukey(n)(1:3)
c endPC
c
c Linux
c         uline = 'mv data1f data1f.' // ukey(n)(1:3)
c endLinux
c
          call system(uline)
c
          nflist = nflist + 1
          uflist(nflist) = 'data1f.' // ukey(n)(1:3)
        endif
        inquire(file='output',exist=qex)
        if (qex) then
c
c PC
          uline = 'ren output output.' // ukey(n)(1:3)
c endPC
c
c Linux
c         uline = 'mv output output.' // ukey(n)(1:3)
c endLinux
c
          call system(uline)
c
          nflist = nflist + 1
          uflist(nflist) = 'output.' // ukey(n)(1:3)
        endif
        inquire(file='slist',exist=qex)
        if (qex) then
c
c PC
          uline = 'ren slist slist.' // ukey(n)(1:3)
c endPC
c
c Linux
c         uline = 'mv slist slist.' // ukey(n)(1:3)
c endLinux
c
          call system(uline)
c
          nflist = nflist + 1
          uflist(nflist) = 'slist.' // ukey(n)(1:3)
        endif
c
        write (nttyo,1310)
 1310   format(/'  The following output files were written:')
        do nf = 1,nflist
          j2 = ilnobl(uflist(nf))
          write (nttyo,1320) uflist(nf)(1:j2)
 1320     format('    ',a)
        enddo
c
  390   continue
c
      enddo
c
c-----------------------------------------------------------------------
c
c     Write some statistics.
c
      if (nfrun .gt. 0) then
        write (nttyo,1200)
        write (nttyo,1350)
 1350   format(/'  The following data files were successfully',
     $  ' processed',/'  with no errors or warnings:',/)
        do n = 1,nfrun
          j2 = ilnobl(ufrun(n))
          j2 = min(j2,74)
          write (nttyo,1360) ufrun(n)(1:j2)
 1360     format(4x,a)
        enddo
      endif
c
      if (nfwarn .gt. 0) then
        write (nttyo,1200)
        write (nttyo,1370)
 1370   format(/'  The following data files were successfully',
     $  ' processed',/'  with no errors, but with warnings:',/)
        do n = 1,nfwarn
          j2 = ilnobl(ufwarn(n))
          j2 = min(j2,74)
          write (nttyo,1360) ufwarn(n)(1:j2)
        enddo
      endif
c
      if (nferr .gt. 0) then
        write (nttyo,1200)
        write (nttyo,1380)
 1380   format(/'  The following data files were processed, but EQPT',
     $  /'  error messages were generated:',/)
        do n = 1,nferr
          j2 = ilnobl(uferr(n))
          j2 = min(j2,74)
          write (nttyo,1360) uferr(n)(1:j2)
        enddo
      endif
c
      if (nfnoex .gt. 0) then
        write (nttyo,1200)
        write (nttyo,1390)
 1390   format(/"  The following data files don't exist:",/)
        do n = 1,nfnoex
          j2 = ilnobl(ufnoex(n))
          j2 = min(j2,74)
          write (nttyo,1360) ufnoex(n)(1:j2)
        enddo
      endif
c
c-----------------------------------------------------------------------
c
  999 write (nttyo,1200)
      write (nttyo,1400)
 1400 format(' All done')
      end
