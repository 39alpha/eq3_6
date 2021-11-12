      program runeq3
c
c     RUNEQ3: EQ3NR interface software
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
c     This is a PC utility program to run EQ3NR, the EQ3/6
c     speciation-solubility code. It is written in Fortran with Lahey
c     compiler extensions, and utilizes DOS system calls. It is
c     a functional near-equivalent of the Linux/UNIX shell script
c     runeq36, when that script is called as "runeq3".
c
c     This program can be adapted to run under Linux/UNIX. Look for
c     system commands run using "call system".
c
c     To use this software, the EQ3NR executable must be in the code
c     directory defined in the environment variable "EQ36CO". The code
c     directory should normally be:
c
c       c:\eq3_6v8.0b\bin
c
c     Substitute "d:" for "c:" if the d: drive is used instead. The
c     EQ36CO variable is normally set by the eq36cfg.bat file.
c
c     The supporting data1 data files must be in the data file
c     directory defined in the environment variable "EQ36DA".
c     The data file directory should normally be:
c
c       c:\eq3_6v8.0b\bin
c
c     Substitute "d:" for "c:" if the d: drive is used instead. The
c     EQ36DA variable is normally set by the eq36cfg.bat file.
c
c     The desired data file is specified by a data file key given as
c     the first argument on the command line. A list of programmed
c     data file keys is kept in the uvkey array, which is also defined
c     below in a data statement. The parameter nkeypa defines the number
c     of programmed keys. A note will be generated if a key is specified
c     which is not on this list. The purpose of programmed keys is to
c     allow this code to screen out EQ3NR input files which are
c     incompatible with the selected data file. If the selected data
c     file's key is not on the programmed list, no such screening
c     is carried out.
c
c     The input file or files to be run are specified on the command
c     line following the data file key. If no leading pathname is given,
c     an input file is assumed to be in the current directory. The files
c     produced by EQ3NR corresponding to each input file are renamed
c     after that file and appear in the current directory.
c
c     This code requires modules from the EPCLIB library.
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer nargpa,ninppa,nlchpa,nkeypa
      parameter (nargpa = 24,ninppa = 150,nlchpa = 80,nkeypa = 45)
c
      integer nfile,noutpt,nscrch,nttyo
c
      integer i,i2,ikey,iopg1,itot,j,jj,j2,k,kk,n,nargmx,nargv,nf,
     $ nferr,nfindf,nflist,nfnoex,nfnwr,nfrun,nfwrex,nlchmx
      integer igkey(1:nkeypa)
c
      integer iargc,ilnobl
c
      logical qerr,qex,qgascr,qgcoef
c
      character(len=80) ucodir,ucodrn,ucomln,udadir,udadrn,udata1,
     $ ueq3c,ufile,uline,ulscr,upathn
      character(len=80) uferr(ninppa),ufindf(ninppa),ufnoex(ninppa),
     $ ufrun(ninppa),ufwrex(ninppa)
      character(len=80) uargv(1:nargpa),uinput(1:ninppa)
      character(len=24) uflist(1:8),ufnwr(1:8)
      character(len=8) uname,ux8
      character(len=8) uvkey(1:nkeypa)
c
c-----------------------------------------------------------------------
c
      data noutpt /10/,nfile  /12/,nscrch /14/,nttyo  /6/
c
      data ucodrn /"EQ36CO"/
      data udadrn /"EQ36DA"/
c
      data (uvkey(n), n = 1,nkeypa) /'dav','com','cmp','ymp','cm1',
     $ 'cm2','cm3','sup','nea','chv','jph','cv1','cv2','cv3','phr',
     $ 'skb','wat','shv','sat','500','750','1kb','2kb','3kb','4kb',
     $ '5kb','bdt','hmw','ypf','ypa','yp1','yp2','yp3','fwe','gmo',
     $ 'smw','fch','frz','pit','pt1','pt2','pt3','pt4','pt5','fmt'/
      data (igkey(n), n = 1,nkeypa) /-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     $ 0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1/
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
      if (nargv .lt. 2) then
        write (nttyo,1000)
 1000   format(' Usage: runeq3 datafilekey inputfile(s)',
     $  /'   datafilekeys = cmp, sup, 1kb, hmw, ypf, and so forth',
     $  /'   inputfiles must end in .3i')
        go to 999
      endif
c
c-----------------------------------------------------------------------
c
      write (nttyo,1010)
 1010 format(' Running RUNEQ3',/)
c
c-----------------------------------------------------------------------
c
c     Get and check the code directory.
c
      call getenv(ucodrn,ucodir)
      j = ilnobl(ucodir)
      k = ilnobl(ucodrn)
      if (ucodir(1:5) .eq. 'ERROR') then
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
c     Get the EQ3NR executable.
c
      ueq3c = ucodir(1:j) // '\eq3nr.exe'
c
c     Check the EQ3NR executable.
c
      inquire(file=ueq3c,exist=qex)
      if (.not.qex) then
        write (nttyo,1060) ucodir(1:j),ucodrn(1:k)
 1060   format (/" * Error - The executable code eq3nr.exe doesn't",
     $  /7x,'exist in the EQ3/6 code directory ',a,
     $  /7x,'defined in the environment variable ',a,'.',/)
        go to 999
      endif
c
c-----------------------------------------------------------------------
c
c     Get and check the data file directory.
c
      call getenv(udadrn,udadir)
      j = ilnobl(udadir)
      k = ilnobl(udadrn)
      if (udadir(1:5) .eq. 'ERROR') then
        write (nttyo,1020) udadrn(1:k)
        go to 999
      endif
c
c     Calling sequence substitutions:
c       udadir for ufile
c
      call texdir(udadir,qex)
      if (.not.qex) then
        write (nttyo,1063) udadir(1:j)
 1063   format (/' * Error - The EQ3/6 data file directory ',a,
     $  /7x,"doesn't exist.",/)
        go to 999
      endif
      write (nttyo,1067) udadir(1:j)
 1067 format('   The data file directory is ',a)
c
c     Check the specified data file.
c
      udata1 = udadir(1:j) // '\data1.' // uargv(1)(1:3)
c
      inquire(file=udata1,exist=qex)
      if (.not.qex) then
        write (nttyo,1069) uargv(1),udadir(1:j),udadrn(1:k)
 1069   format (/' * Error - The data file "data1.',a3,'" ',"doesn't",
     $  /7x,'exist in the EQ3/6 data file directory ',a,
     $  /7x,'defined in the environment variable ',a,'.',/)
        go to 999
      endif
c
c     Copy the specifed data1 file into the current directory as
c     "data1".
c
      j = ilnobl(udata1)
      write (nttyo,1070) udata1(1:j)
 1070 format('   The supporting data file is ',a)
      call kilfil('data1')
c
c PC
      uline = 'copy ' // udata1(1:j) // ' data1' // ' > c:\nul'
c endPC
c
c Linux
c     uline = 'ln -s ' // udata1(1:j) // ' data1'
c endLinux
c
      call system(uline)
c
c-----------------------------------------------------------------------
c
c     Is the data file key in the list of programmed keys?
c
      ikey = 0
      do i = 1,nkeypa
        if (uargv(1)(1:4) .eq. uvkey(i)(1:4)) then
          ikey = i
          go to 120
        endif
      enddo
      j = ilnobl(uargv(1))
      write (nttyo,1071) uargv(1)(1:j)
 1071 format (/" * Note - Don't recognize the data file key ",
     $ '"',a,'".',
     $ /7x,"Won't be able to screen out possible mismatches",
     $ /7x,'between the data file and the input file or files.',/)
  120 continue
c
c-----------------------------------------------------------------------
c
c     Delete any existing copies of the generic files that EQ3NR reads
c     or writes, if any copies already exist.
c
      call kilfil('input')
      call kilfil('output')
      call kilfil('pickup')
c
c-----------------------------------------------------------------------
c
c     Initialize counters:
c
c       nfrun  = number of input files run without generating any
c                  EQ3NR error messages
c       nfwrex = number of input files not run because they have the
c                  wrong file name extension
c       nfindf = number of input files requiring a different kind
c                  of data file
c       nferr  = number of input files for which EQ3NR error messages
c                  are generated
c       nfnoex = number of input files which are not run because they
c                  don't exist
c
      nfrun = 0
      nfwrex = 0
      nfindf = 0
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
      do n = 2,nargv
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
c
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
     $        /7x,"doesn't exist. It will not be run.",/)
              nfnoex = nfnoex + 1
              ufnoex(nfnoex) = uinput(i)(1:j)
              go to 190
            endif
          endif
c
c         Does the current input file have the proper file name
c         extension?
c
          k = index(uinput(i),'.3i')
          if (k .eq. 0) then
            write (nttyo,1092) uinput(i)(1:j)
 1092       format (/" * Error - Can't run EQ3NR on input file ",'"',a,
     $      '"',/7x,"because it doesn't have a .3i file name",
     $      " extension.",/)
            kk = index(uinput(i),'.6i')
            if (kk .gt. 0) then
              write (nttyo,1095)
 1095         format (7x,"It's an EQ6 input file.",/)
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
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c         Check options on the input file.
c
          open(nfile,file='input',status='old')
c
          qgascr = .false.
          if (ikey .gt. 0) then
c
c           Determine the activity coefficient model type specified
c           on the input file.
c
            call ggammo(iopg1,nfile,qgcoef)
c
c           Test the activity coefficient formalism specified on the
c           on the input file to see if it is compatible with that
c           on the data file.
c
            if (qgcoef) then
              if (iopg1 .eq. -1) iopg1 = 0
              if (iopg1 .ne. igkey(ikey)) then
                write (nttyo,1093)
 1093           format(/" * Note - This input file requires an",
     $          " activity coefficient",/7x,"model that isn't",
     $          " compatible with the specified",/7x,"data file.",
     $          " It will not be run.",/)
                j = ilnobl(uinput(i))
                nfindf = nfindf + 1
                ufindf(nfindf) = uinput(i)(1:j)
                qgascr = .true.
              endif
            else
              write (nttyo,1094)
 1094         format(/" * Note - Can't determine the activity",
     $        " coefficient",/7x,"option on this input file.",/)
            endif
          else
            write (nttyo,1096)
 1096       format(/" * Note - Can't check the activity",
     $      " coefficient option",/7x,"for this input file",
     $      " because the specified data file",/7x,"keystring",
     $      " is unknown to this interface software.",/)
          endif
c
          if (qgascr) then
            close(nfile,status='delete')
            go to 190
          else
            close(nfile)
          endif
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c         Get the name of the input file, less the pathname and file
c         name extension.
c
          k = ilnobl(upathn)
          uline = uinput(i)(k + 1:80)
          k = index(uline,'.3i')
          uname = uline(1:k - 1)
c
c         Delete any old copies of files to be created.
c
          j = ilnobl(uname)
          ufile = uname(1:j) // '.3o'
          call kilfil(ufile)
          ufile = uname(1:j) // '.3p'
          call kilfil(ufile)
c
c         Run EQ3NR.
c
          call system(ueq3c)
c
c         Delete the input file (the copy of the current input file).
c
          call kilfil('input')
c
c         Check for errors in running the current input file.
c         Warnings will be ignored here.
c
          qerr = .false.
c
          inquire(file='output',exist=qex)
          if (.not.qex) qerr = .true.
c
          if (.not.qerr) then
            open(noutpt,file='output',status='old')
            call ckferr(noutpt,qerr)
            close(noutpt)
          endif
c
c         If there were errors, write a warning message.
c
          if (qerr) then
            write (nttyo,1130)
 1130       format(/' * Warning - EQ3NR has encountered errors in',
     $      /7x,'running this input file. Check the output file',
     $      /7x,'for details.',/)
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
c         Rename the output files.
c
          nflist = 0
          nfnwr = 0
          j = ilnobl(uname)
          inquire(file='output',exist=qex)
          if (qex) then
c
c PC
            uline = 'ren output ' // uname(1:j) // '.3o'
c endPC
c
c Linux
c           uline = 'mv output ' // uname(1:j) // '.3o'
c endLinux
c
            call system(uline)
c
            nflist = nflist + 1
            uflist(nflist) = uname(1:j) // '.3o'
          else
            nfnwr = nfnwr + 1
            ufnwr(nfnwr) = uname(1:j) // '.3o'
          endif
          inquire(file='pickup',exist=qex)
          if (qex) then
c
c PC
            uline = 'ren pickup ' // uname(1:j) // '.3p'
c endPC
c
c Linux
            uline = 'mv pickup ' // uname(1:j) // '.3p'
c endLinux
c
            call system(uline)
c
            nflist = nflist + 1
            uflist(nflist) = uname(1:j) // '.3p'
          else
            nfnwr = nfnwr + 1
            ufnwr(nfnwr) = uname(1:j) // '.3p'
          endif
c
          write (nttyo,1135)
 1135     format(/'  The following output files were written:')
          do nf = 1,nflist
            j2 = ilnobl(uflist(nf))
            write (nttyo,1137) uflist(nf)(1:j2)
 1137       format('    ',a)
          enddo
c
          if (nfnwr .gt. 0) then
            write (nttyo,1138)
 1138       format(/'  The following output files were not written:')
            do nf = 1,nfnwr
              j2 = ilnobl(ufnwr(nf))
              write (nttyo,1137) ufnwr(nf)(1:j2)
            enddo
          endif
c
  190     continue
        enddo
  195   continue
      enddo
c
c-----------------------------------------------------------------------
c
c     Delete the data1 file.
c
      call kilfil('data1')
c
c-----------------------------------------------------------------------
c
c     Write some statistics.
c
      if (nfrun .gt. 0) then
        write (nttyo,1078)
        write (nttyo,1140)
 1140   format(/'  The following input files were run without',
     $  /'  generating any EQ3NR error messages:',/)
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
 1155   format(/'  The following input files were not run because they',
     $  /"  don't have the correct file name extension:",/)
        do n = 1,nfwrex
          j2 = ilnobl(ufwrex(n))
          j2 = min(j2,74)
          write (nttyo,1150) ufwrex(n)(1:j2)
        enddo
      endif
c
      if (nfindf .gt. 0) then
        write (nttyo,1078)
        write (nttyo,1157)
 1157   format(/'  The following input files were not run because',
     $  /"  they require an activity coefficient model that isn't",
     $  /'  compatible with the specified data file:',/)
        do n = 1,nfindf
          j2 = ilnobl(ufindf(n))
          j2 = min(j2,74)
          write (nttyo,1150) ufindf(n)(1:j2)
        enddo
      endif
c
      if (nferr .gt. 0) then
        write (nttyo,1078)
        write (nttyo,1160)
 1160   format(/'  The following input files were run, but EQ3NR error',
     $  /'  messages were generated:',/)
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
 1170   format(/"  The following input files were not run because",
     $  /"  they don't exist:",/)
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
