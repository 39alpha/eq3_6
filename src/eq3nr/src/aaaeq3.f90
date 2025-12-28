subroutine aaaeq3(usteq3,uveeq3)
    !! EQ3NR: EQ3NR Speciation-Solubility Code
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
    !! This code computes chemical models of aqueous solutions. Problems
    !! are composed using mostly analytical data. Thermodynamic data
    !! are used in the calculations. This code is similar in function
    !! to the codes in the WATEQ series. It can be used by itself to
    !! determine aqueous speciation and solution-mineral saturation
    !! relations. It can also be used to initialize EQ6 calculations.
    !! The user defines the run parameters on the file called input.
    !! The file called data1 provides the basic supporting thermodynamic,
    !! etc., data. EQ3NR writes a pickup file and an output file called
    !! output. The pickup file is used to pass data to EQ6. This file
    !! is the second half of an EQ6 input file.
    !! This code requires subroutines from the following EQ3/6 libraries:
    !!   EQLIB
    !!   EQLIBG
    !!   EQLIBU
    !! This subroutine is designed to ensure that the copyright statement
    !! and legal and other disclaimers appear at the beginning of the
    !! concatenated source code for EQ3NR. Concatenation is normally
    !! alphabetical, hence the form of the name. This subroutine returns
    !! the stage and version numbers of this computer code.
    !! This subroutine is called by:
    !!   EQ3NR/eq3nr.f
    !! Principal input:
    !!   None
    !! Principal output:
    !!   usteq3 = EQ3NR stage number
    !!   uveeq3 = EQ3NR version number
    implicit none

    ! Calling sequence variable declarations.
    character(len=8) :: usteq3
    character(len=8) :: uveeq3

    ! Local variable declarations.
    character(len=8) :: ust
    character(len=8) :: uve

    data ust /'R43a    '/
    data uve /'8.0a    '/

    usteq3 = ust
    uveeq3 = uve
end subroutine aaaeq3