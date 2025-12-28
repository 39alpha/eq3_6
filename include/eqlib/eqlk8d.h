! eqlk8d.h
!     Strings, arrays, and such for the menu-style ("D") input format
!     for INPUT file options exclusive to EQ6 in Version 8. This INCLUDE
!     INCLUDE file is referenced by EQ6, EQ3NR, and XCON6.
!     The relevant arrays are declared in the EQLIB INCLUDE file
!     eqlk8.h.
!     Reactant type strings.
data urcjco(0) /'Pure mineral            '/
data urcjco(1) /'Solid solution          '/
data urcjco(2) /'Special reactant        '/
data urcjco(3) /'Aqueous species         '/
data urcjco(4) /'Gas species             '/
data urcjco(5) /'Generic ion exchanger   '/

! Reactant status strings.
data urcjre(-1) /'Saturated, reacting     '/
data urcjre(0)  /'Reacting                '/
data urcjre(1)  /'Exhausted               '/
data urcjre(2)  /'Saturated, not reacting '/

! Rate law strings.
data urcnrk(-1,1) /'Use backward rate law   '/
data urcnrk(0,1)  /'Illegal value           '/
data urcnrk(1,1)  /'Relative rate equation  '/
data urcnrk(2,1)  /'TST rate equation       '/
data urcnrk(3,1)  /'Linear rate equation    '/
data urcnrk(-1,2) /'Use forward rate law    '/
data urcnrk(0,2)  /'Partial equilibrium     '/
data urcnrk(1,2)  /'Relative rate equation  '/
data urcnrk(2,2)  /'TST rate equation       '/
data urcnrk(3,2)  /'Linear rate equation    '/

! Rate direction header strings.
data urcrld(1) /'Forward rate law        '/
data urcrld(2) /'Backward rate law       '/

! Rate direction sign strings.
data urcrls(1) /'+'/
data urcrls(2) /'-'/

! Relative rate header strings.
data urcrel(1,1) /'dXi(n)/dXi (mol/mol)    '/
data urcrel(2,1) /'d2Xi(n)/dXi2 (mol/mol2) '/
data urcrel(3,1) /'d3Xi(n)/dXi3 (mol/mol3) '/
data urcrel(1,2) /'-dXi(n)/dXi (mol/mol)   '/
data urcrel(2,2) /'-d2Xi(n)/dXi2 (mol/mol2)'/
data urcrel(3,2) /'-d3Xi(n)/dXi3 (mol/mol3)'/

! Nxopt mineral subset-selection suppression option strings:
data uxopti(1) /'None    '/
data uxopti(2) /'All     '/
data uxopti(3) /'Alwith  '/
data uxopti(4) /'Allwith '/

! Jflgi strings for EQ6.
data ujf6(1) /'Moles           '/
data ujf6(2) /'Make non-basis  '/

