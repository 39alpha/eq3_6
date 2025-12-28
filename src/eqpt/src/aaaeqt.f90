subroutine aaaeqt(usteqt,uveeqt)
    !! EQPT: EQ3/6 Data File Preprocessor Code
    !! EQ3/6 version 8.0a R43a (Patched 10/01/2009)
    !! Last revised 04/24/02 by TJW
    !! Copyright (c) 1987, 1990-1993, 1995, 1997, 2002 The Regents of
    !! the University of California, Lawrence Livermore National
    !! Laboratory. All rights reserved.
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
    !! This code is a data file preprocessor for EQ3/6. It reads a
    !! primary data file (a DATA0 file) and writes a secondary,
    !! unformatted data file (a DATA1 file), which can be read by
    !! EQ3NR and EQ6. EQPT conducts a number of checks on the data as
    !! it processes it. For example, it checks reactions for mass
    !! balance and electrical balance. The code also fits interpolating
    !! polynomials to temperature-dependent data that are represented
    !! on a DATA0 file by values on a temperature grid. The coefficients
    !! of these polynomials are written on the DATA1 file in place of
    !! the original gridded data.
    !! This code requires subroutines from the following EQ3/6 libraries:
    !!   EQLIBU
    !! This subroutine is designed to ensure that the copyright statement
    !! and legal and other disclaimers appear at the beginning of the
    !! concatenated source code for EQPT. Concatenation is normally
    !! alphabetical, hence the form of the name. This subroutine returns
    !! the stage and version numbers of this computer code.
    !! This subroutine is called by:
    !!   EQPT/eqpt.f
    !! Input:
    !!   None
    !! Output:
    !!   usteqt = EQPT stage number
    !!   uveeqt = EQPT version number
    implicit none

    ! Calling sequence variable declarations.
    character(len=8) :: usteqt
    character(len=8) :: uveeqt

    ! Local variable declarations.
    character(len=8) :: ust
    character(len=8) :: uve

    data ust /'R43a    '/
    data uve /'8.0a    '/

    usteqt = ust
    uveeqt = uve
end subroutine aaaeqt