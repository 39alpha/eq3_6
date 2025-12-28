! eqlj8.h
!     Strings, arrays, and such for the menu-style ("D") input format
!     for INPUT file options exclusive to EQ3NR in Version 8. This
!     INCLUDE file is referenced by EQ3NR and XCON3.
!       Jflgi options:
!         njfxpa = the largest jflgi value
!         ujf3   = array of strings corresponding to jflgi values
!       Note: a distinct, smaller set of jflgi option strings (ujf6)
!       is used by EQ6.
!     The strings are defined in the EQLIB include file eqlj8d.h.
integer :: njfxpa
parameter (njfxpa = 30)

character(len=16) :: ujf3

common /eqlj8c/ ujf3(-1:njfxpa)

