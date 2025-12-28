! eqltf1.h
!     This is a pseudo-data file for the titration factors which are
!     required to compute the HCO3-CO3-OH total alkalinity in normal
!     (e.g., dilute) groundwaters. The species considered here are
!     the three anions and their 1:1 ion pairs with Na+, K+, Mg++,
!     and Ca++. Addition of higher order ion pairs, such as
!     Mg4(OH)4++++ or Ca(OH)2(aq), is not recommended, as these
!     are not signficant in dilute groundwaters. Do not add other
!     kinds of species which may contribute to total alkalinity,
!     such as acetate, other organics, borate, silicate, phosphate,
!     or transition metal-hydroxy complexes. Instead, add any such
!     needed species to the pseudo-data file in the companion EQLIB
!     INCLUDE file eqltf2.h.
!     Note that the list of species terminates with an "endit."
!     This pseudo-data file is referenced by:
!       EQ3NR/eq3nr.f
!       EQ6/eq6.f
integer, parameter :: ntf1a_par = 15

character(len=32), dimension(ntf1a_par) :: utf1xd = (/ '"CaCO3(aq)           ", 2.  ','"CaHCO3+             ", 1.  ','"CaOH+               ", 1.  ','"CO3--               ", 2.  ','"HCO3-               ", 1.  ','"KCO3-               ", 2.  ','"KHCO3(aq)           ", 1.  ','"KOH(aq)             ", 1.  ','"MgCO3(aq)           ", 2.  ','"MgHCO3+             ", 1.  ','"MgOH+               ", 1.  ','"NaCO3-              ", 2.  ','"NaHCO3(aq)          ", 1.  ','"NaOH(aq)            ", 1.  ','"OH-                 ", 1.  '/)

