! eqlo8.h
!     Strings, arrays, and such for the menu-style ("D") input format
!     for INPUT file options common to EQ3NR and EQ6 in Version 8.
!     INCLUDE file is referenced by EQ6, EQ3NR, XCON3, and XCON6.
!       Nxmod alter/supress options:
!         ukxm   = array of valid option strings for nxmod options
!       Exchange model options for generic ion exchangers:
!         ugexmv = array of valid exchange model strings for generic
!                    ion exchangers
!       Units options for thermodynamic parameters of generic ion
!         exchangers:
!         uxfuni = array of strings for units options for
!                    thermodynamic functions of generic ion
!                    exchangers
!       Major option switch (iopt, iopr, iodb, and iopg) arrays:
!         jptxpa = maximum number of choices per iopt option switch
!         jprxpa = maximum number of choices per iopr print  option
!                    switch
!         jdbxpa = maximum number of choices per iodb debugging print
!                    option switch
!         jpgxpa = maximum number of choices per iopg activity
!           coefficient option switch
!         nptxpa = maximum number of iopt option switches for which
!                  strings are defined
!         nprxpa = maximum number of print option switches for which
!                    strings are defined
!         ndbxpa = maximum number of debugging print option switches
!                    for which strings are defined
!         npgxpa = maximum number of activity coefficient option
!                    switches for which strings are defined
!         uopttx = array of title strings for iopt option switches
!         uoprtx = array of title strings for iopr option switches
!         uodbtx = array of title strings for iodb option switches
!         uopgtx = array of title strings for iopg option switches
!         uoptcx = array of code relevancy strings for iopt option
!                    switches
!         uoprcx = array of code relevancy strings for iopr option
!                    switches
!         uodbcx = array of code relevancy strings for iodb option
!                    switches
!         uopgcx = array of code relevancy strings for iopg option
!                    switches
!         uoptox = array of option strings for iopt option switches
!         uoprox = array of option strings for iopr option switches
!         uodbox = array of option strings for iodb option switches
!         uopgox = array of option strings for iopg option switches
!         ioptox = array of option values for iopt option switches
!         ioprox = array of option values for iopr option switches
!         iodbox = array of option values for iodb option switches
!         iopgox = array of option values for iopg option switches
!     The strings are defined in the EQLIB include file eqlo8d.h.
integer :: jptxpa
integer :: jprxpa
integer :: jdbxpa
integer :: jpgxpa
integer :: nptxpa
integer :: nprxpa
integer :: ndbxpa
integer :: npgxpa
integer :: nvet_par
parameter (jptxpa = 5,jprxpa = 5,jdbxpa = 5,jpgxpa = 5,nptxpa = 20,nprxpa = 20,ndbxpa = 20,npgxpa = 20,nvet_par = 10)

integer :: ioptox
integer :: ioprox
integer :: iodbox
integer :: iopgox

common /eqlo8i/ ioptox(jptxpa,nptxpa),ioprox(jprxpa,nprxpa),iodbox(jdbxpa,ndbxpa),iopgox(jpgxpa,npgxpa)

character(len=72) :: uopttx
character(len=72) :: uoprtx
character(len=72) :: uodbtx
character(len=72) :: uopgtx
character(len=72) :: uoptox
character(len=72) :: uoprox
character(len=72) :: uodbox
character(len=72) :: uopgox
character(len=24) :: uoptcx
character(len=24) :: uoprcx
character(len=24) :: uodbcx
character(len=24) :: uopgcx
character(len=24) :: ugexmv
character(len=16) :: ukxm
character(len=8) :: uxfuni

common /eqlo8c/ uopttx(nptxpa),uoprtx(nprxpa),uodbtx(ndbxpa),uopgtx(npgxpa),uoptox(jptxpa,nptxpa),uoprox(jprxpa,nprxpa),uodbox(jdbxpa,ndbxpa),uopgox(jpgxpa,npgxpa),uoptcx(nptxpa),uoprcx(nprxpa),uodbcx(ndbxpa),uopgcx(npgxpa),ugexmv(nvet_par),ukxm(-1:2),uxfuni(4)

