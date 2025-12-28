subroutine aaaeq6(usteq6,uveeq6)
    !! EQ6: EQ6 Reaction-Path Code
    !! EQ3/6 version 8.0a R43a (Patched 10/01/2009)
    !! Last revised 04/24/02 by TJW
    !! Copyright (c) 1987, 1990-1993, 1995, 1997, 2002 The Regents
    !! of the University of California, Lawrence Livermore National
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
    !! This code computes reaction-path models of water/rock and water/
    !! rock/waste systems. A model of an aqueous solution must first
    !! be created using the EQ3NR code. Problems are composed by choosing
    !! substances to react with this water. Such a substance can include
    !! another aqueous solution. Thermodynamic data are used in the
    !! calculations. Kinetic data of some sort is also required. The
    !! code predicts rates of reaction, changes in fluid composition
    !! (including pH and Eh), and the precipitation and possible
    !! re-dissolution of secondary phases. The temperature and pressure
    !! may vary along a reaction path. This code primarily models the
    !! chemistry in a small box. It can also model the reaction of a
    !! packet of water (the first packet) traversing a reactive medium.
    !! This system differs in that secondary minerals are left behind
    !! and can't re-dissolve. The fugacities of selected gases may be
    !! fixed to simulate reaction in the presence of a large gas volume,
    !! such as the atmosphere. This code is similar in function to such
    !! codes as CHILLER, REACT, and SOLMINEQ.88.
    !! The user defines the run parameters on the file called input.
    !! The file called DATA1 provides the basic supporting thermodynamic,
    !! etc., data. EQ6 writes a pickup file and an output file called
    !! output. The pickup file is an EQ6 input file which can be used
    !! to restart a calculation. This code also writes data to the
    !! tab file. This is an ASCII file containing tables which can be
    !! used in graphics post-processing. This code also creates a tabx
    !! (scrambled tab) file. Normally this is just used at the end of
    !! a run to write the tab file. However, an option exists to use it
    !! to create a tab file that spans more than one run.
    !! This code requires subroutines from the following EQ3/6 libraries:
    !!   EQLIB
    !!   EQLIBG
    !!   EQLIBU
    !! This subroutine is designed to ensure that the copyright statement
    !! and legal and other disclaimers appear at the beginning of the
    !! concatenated source code for EQ6. Concatenation is normally
    !! alphabetical, hence the form of the name. This subroutine returns
    !! the stage and version numbers of this computer code.
    !! This subroutine is called by:
    !!   EQ6/eq6.f
    !! Input:
    !!   None
    !! Output:
    !!   usteq6 = EQ6 stage number
    !!   uveeq6 = EQ6 version number
    implicit none

    ! Calling sequence variable declarations.
    character(len=8) :: usteq6
    character(len=8) :: uveeq6

    ! Local variable declarations.
    character(len=8) :: ust
    character(len=8) :: uve

    data ust /'R43a    '/
    data uve /'8.0a    '/

    usteq6 = ust
    uveeq6 = uve
end subroutine aaaeq6