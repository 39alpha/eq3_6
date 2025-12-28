subroutine aaaxc6(ustxc6,uvexc6)
    !! XCON6: EQ6 Input File Converter
    !! EQ3/6 version 8.0 R43a (Patched 10/01/2009)
    !! Last revised 04/24/02 by TJW
    !! Copyright (c) 1993, 1995, 1997, 2002 The Regents of the University
    !! of California, Lawrence Livermore National Laboratory. All rights
    !! reserved.
    !! This work was produced at the University of California,
    !! Lawrence Livermore National Laboratory (UC LLNL) under
    !! contract no. W-7405-ENG-48 between the U.S. Department of
    !! Energy (DOE) and The Regents of the University of California
    !! (University) for the operation of UC LLNL. Copyright is
    !! reserved to the University for purposes of controlled
    !! dissemination, commercialization through formal licensing,
    !! or other disposition under terms of Contract 48; DOE
    !! policies, regulations, and orders; and U.S. statutes.
    !!                        DISCLAIMER
    !! This computer code was prepared as an account of work
    !! sponsored by an agency of the United States Government.
    !! Neither the United States Government nor the University of
    !! California nor any of their employees, makes any warranty,
    !! express or implied, or assumes any liability or responsibility
    !! for the accuracy, completeness, or usefulness of any
    !! information, apparatus, product, or process disclosed, or
    !! represents that its use would not infringe privately-owned
    !! rights. Reference herein to any specific commercial products,
    !! process, or service by trade name, trademark, manufacturer,
    !! or otherwise, does not necessarily constitute or imply its
    !! endorsement, recommendation, or favoring by the United States
    !! Government or the University of California. The views and
    !! opinions of authors expressed herein do not necessarily state
    !! or reflect those of the United States government or the
    !! University of California, and shall not be used for
    !! advertising or product endorsement purposes.
    !! See the readme.txt file that came with this software for further
    !! information, including contact information and references.
    !! This is a utility program to convert EQ6 input files from
    !! one version level to another and from one format to another.
    !! This program is most conveniently run using the XCIF6 utility
    !! discussed below. It may also be run directly.
    !! The old input file is known as 'input'. It is preserved in the
    !! conversion. The new input file is known as 'newin'. The version
    !! levels are:
    !!   '6.0' = versions 6.0-6.1 (3245.0288, 3245.0888)
    !!   '7.0' = versions 7.0-7.1
    !!   '7.2' = version 7.2
    !!   '8.0' = version 8.0
    !! The formats are:
    !!   'W' = compact
    !!   'D' = menu-style
    !! Version level '8.0' is presently not accessible. 'D' format does
    !! not exist for version level '6.0'.
    !! To specify the choice of version level and format of the new
    !! input file, use the control file IXCON. The same control file
    !! also governs the EQ3NR input file converter XCON3. A copy of
    !! the contents of this file can be found in the comments in
    !! EQLIBU/rddixc.f.
    !! This code automatically determines the format of the old input
    !! file. It may have trouble determining the version level, however.
    !! It will first look for a string of the form 'Version level= X.X'
    !! or 'Version number= X.X' in the first title of the first problem
    !! on this file. Here 'X.X' denotes the version level. If none is
    !! found, the code will use a default value specified on the IXCON
    !! file. If no IXCON file is present, the code will use a default
    !! specified in the main program.
    !! If a version level marker as described above exists on the old
    !! input file, it will be modified on the new input file to give
    !! its correct version level. If no version level marker is on the
    !! old input file, one will be added to the new input file.
    !! The version level and format of the new input file may match
    !! those of the old input file. The new input file may be identical
    !! to the old one, or made somewhat more pretty.
    !! As noted above, it is more convenient to run XCON6 using the
    !! interface software XCIF6 than to run XCON6 directly. On a UNIX
    !! system, XCIF6 is a symbolic link to the shell script xcif36RXX,
    !! where RXX is a revision number (R01, R02, R03, etc.). On a
    !! PC system, XCIF6 is a Fortran program (XCIF6.EXE; the source
    !! code file is XCIF6.FOR) containing DOS system calls. XCIF6
    !! replaces the specified input files. XCIF6 also creates the
    !! IXCON control file used by XCON6, and deletes it prior to
    !! completing execution. To convert a set of EQ6 input files from
    !! version level 7.2 ('W' or 'D' format) to version level 8.0,
    !! 'W' format, enter:
    !!   xcif6 7.2 8.0 W *.6i
    !! To convert the same files to 'D' format instead, enter:
    !!   xcif6 7.2 8.0 D *.6i
    !! XCON3 and XCIF3 are analogues of XCON6 and XCIF6, respectively.
    !! They perform the analogous conversions of EQ3NR input files.
    !! This code requires subroutines from the following EQ3/6 libraries:
    !!   EQLIBU
    !! This subroutine is designed to ensure that the copyright statement
    !! and legal and other disclaimers appear at the beginning of the
    !! concatenated source code for XCON6. Concatenation is normally
    !! alphabetical, hence the form of the name. This subroutine returns
    !! the stage and version numbers of this computer code.
    !! This subroutine is called by:
    !!   XCON6/xcon6.f
    !! Principal input:
    !!   None
    !! Principal output:
    !!   ustxc6 = XCON6 stage number
    !!   uvexc6 = XCON6 version number
    implicit none

    ! Calling sequence variable declarations.
    character(len=8) :: ustxc6
    character(len=8) :: uvexc6

    ! Local variable declarations.
    character(len=8) :: ust
    character(len=8) :: uve

    data ust /'R43a    '/
    data uve /'8.0     '/

    ustxc6 = ust
    uvexc6 = uve
end subroutine aaaxc6