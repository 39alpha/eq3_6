subroutine aaaelg(ustelg,uveelg)
    !! EQLIBG: EQ3/6 Library - Activity Coefficients
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
    !! This library supports the EQ3NR and EQ6 codes. It contains
    !! subroutines that deal with calculations specific to activity
    !! coefficient models. The EQLIB library contains all other
    !! subroutines that address elements of geochemical modeling
    !! calculations that are common to EQ3NR and EQ6. The EQLIBU
    !! library contains subroutines with generic utility and
    !! mathematical functions.
    !! This subroutine is designed to ensure that the copyright statement
    !! and legal and other disclaimers appear at the beginning of the
    !! concatenated source code for EQLIBG. Concatenation is normally
    !! alphabetical, hence the form of the name. This subroutine is
    !! generally called by the main program of any code that utilizes
    !! subroutines from this library. It returns the stage and version
    !! numbers of this library.
    !! This subroutine is called by:
    !!   EQ3NR/eq3nr.f
    !!   EQ6/eq6.f
    !! Input:
    !!   None
    !! Output:
    !!   ustelg = EQLIBG stage number
    !!   uveelg = EQLIBG version number
    implicit none

    ! Calling sequence variable declarations.
    character(len=8) :: ustelg
    character(len=8) :: uveelg

    ! Local variable declarations.
    character(len=8) :: ust
    character(len=8) :: uve

    data ust /'R43a    '/
    data uve /'8.0a    '/

    ustelg = ust
    uveelg = uve
end subroutine aaaelg