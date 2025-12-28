! eqlk8.h
!     Strings, arrays, and such for the menu-style ("D") input format
!     for INPUT file options exclusive to EQ6 in Version 8. This INCLUDE
!     file is referenced by EQ6, EQ3NR, and XCON6.
!       Reactant type strings:
!         urcjco = array of strings corresponding to jcode reactant
!                    type flags
!       Reactant status strings:
!         urcjre = array of strings corresponding to jreac reactant
!                    status flags
!       Rate law strings:
!         urcnrk = array of strings corresponding to nrk rate law codes
!       Rate direction header strings:
!         urcrld = array of strings corresponding to rate law directions
!                    (forward and backward)
!       Rate direction sign strings:
!         urcrls = array of strings corresponding to signs of rate law
!                    directions (+ and -)
!       Relative rate header strings:
!         urcrel = array of strings corresponding to relative rates of
!                    various orders
!       Nxopt mineral subset-selection suppression option strings:
!         uxopti = array of strings corresponding to valid options
!       Jflgi strings for EQ6:
!         ujf6   = array of jflgi option strings for EQ6
!       Note: a distinct, larger set of jflgi option strings (ujf3)
!       is used by EQ3NR.
!     The strings are defined in the EQLIB include file eqlk8d.h.
character(len=24) :: urcrld
character(len=24) :: urcjco
character(len=24) :: urcjre
character(len=24) :: urcnrk
character(len=24) :: urcrel
character(len=16) :: ujf6
character(len=16) :: uxopti
character(len=2) :: urcrls

common /eqlk8c/ urcrld(2),urcjco(0:5),urcjre(-1:2),urcnrk(-1:4,2),urcrel(3,2),ujf6(2),uxopti(4),urcrls(2)

