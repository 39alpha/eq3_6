
 EQ3/6, Version 8.0a (EQ3/6-V8-REL-V8.0a-PC)
 EQPT Data File Preprocessor Code (EQ3/6-V8-EQPT-EXE-R43a-PC)
 Supported by the following EQ3/6 libraries:
   EQLIBU (EQ3/6-V8-EQLIBU-LIB-R43a-PC)

 Copyright (c) 1987, 1990-1993, 1995, 1997, 2002 The Regents of the
 University of California, Lawrence Livermore National Laboratory.
 All rights reserved.

 This work is subject to additional statements and
 disclaimers which may be found in the README.txt file
 included in the EQ3/6 software transmittal package.


 Run  10:35:00  04Dec2021


           Number of elements =    86
           Number of basis species =   189

 data0.geo.R2
 
 THERMODYNAMIC DATABASE
 After data0.ymp.R5
 
 This data file is consistent with the Rimstidt paradigm for quartz
 solubility, which makes SiO2(aq) more stable than in the case of
 the Fournier paradigm.
 
 Fixed second range log K grids for Sepiolite and SiO2(am).
 Replaced YMP Laumontite with SUPCRT92 version (Rimstidt-consistent)
 Added Beidellite-ss, Nontronite-ss, and Montmorillonite_ss
 Renamed Saponite-tri to Saponite-ss
 
 -----------------------
 This EQ3/6 database file is an update to the previous version of the
 thermodynamic database data0.ymp.R2 DTN: MO0210SPATHDYN.000. The changes
 result from both correcting errors and adding new data. The changes in
 this version of the database (data0.ymp.R5) relative to the immediately
 preceding version (data0.ymp.R4) are:
 
  1. The log K values were corrected for the following neptunium solids:
       Na3NpF8, Np, Np2O5, NpBr3, NpBr4, NpC0.91, NpCl3, NpCl4, NpF3, NpF4,
       NpF5, NpF6, NpI3,NpN, NpO2, NpOBr2, and NpOCl2.
  2. The log K values were corrected for the following plutonium solids:
       Cs3PuCl6, CsPu2Cl7, Pu, Pu2C3, Pu2O3,Pu3C2, PuBr3, PuC0.84, PuCl3, PuF3,
       PuF6, PuI3, PuO1.61, PuO2, PuOBr, PuOCl, PuOF, and PuOI.
  3. The log K values were corrected for the following minerals: Analcime,
       Chabazite, Clinoptilolite, Erionite, Laumontite, Mesolite, Phillipsite,
       and Stellerite.
  4. The log K values were corrected for HF2- [[DTN SN0410T0510404.001]
       HF2-_Dissoc1_CFJC_fix.xls
  5. The stoichiometry of the of the Pu solid species PuO2(OH)2:H2O which
       was incorrectly labeled as PuO2(OH)2:2H2O in data0.ymp.R2
       (DTN:MO0210SPATHDYN.000). The datablock with the updated Pu
       stoichiometry was obtained from data0.ymp.R3 (DTN MO0312SPATDMIF.000).
       The corrected stoichiometry is consistent with that of the source
       (01lem/nea).
  6. The log K value at 25 deg.C for Syngenite was updated to -7.4484.
       Source was data0.ypf.R0 (DTN: SN0302T0510102.002)
  7. Updated the H2O stoichiometry for REE phosphate solids. Fixed the
        misinterpretation of "xH2O" from 95spa/bru as "10H2O".
  8. Replacement of NiCO3(s) data block by the updated one in data0.ypf.R2.
  9. Updated PO2, PO3, PO4, P2O7, and PO3F species affected by the thermo-
       dynamic key data discrepancy. Two PO2F2 species were added.
 10. Eskolaite (Cr2O3) and CrO2 data were updated; reactions were written in
       terms of Cr+++, not CrO4-- or other Cr basis species.
 11. The reaction for the species CrO3Cl- was rewritten using CrO4-- instead
       of Cr+++.
 12. Addition of the phases Palygorskite, Amorphous Sepiolite, and
       Antigorite(am)
 13. The log K values of metal solids update tabulated in the R4 report but
       not included in that version of the database were added in this version.
 14. Addition of the phases NiMoO4 and Kogarkoite (Na3SO4F).
 15. Addition of SnO2(am), Sn(OH)5- and Sn(OH)6--.
 16. Reactions for P2O7 solids rewritten in terms of H2P2O7-- instead of
       HPO4--.
 17. PuCl++ data block removed, data de-selected by 03GUI/NEA errata and
       corrigenda.
 18. Update of 'azer0' values in accord with suggested values in
       ANL-WIS-GS-000003 REV 01 report.
 19. Update of aqueous species in 'azer0' list.
 20. Updates that include additions and removal of aqueous Th species
       consistent with the model of 05alt/nec.
 21. Updates of uranyl phases that include the addition of Becquerelite,
       Compreignacite, Weeksite, Boltwoodite, Dehydrated Schoepite,
       Sklodowskite and Studtite.
 22. Addition of the species CaSeO4(aq) and CaSeO4:2H2O.
 23. Removed polynuclear species Ni2OH+++ and Ni4(OH)++++; updated log K
       values for NiOH+, Ni(OH)2(aq), Ni(OH)3- using data from data.ypf.R2
 24. Removed polynuclear species Cr3(OH)4(5+) and Cr2(OH)2++++; replaced
       CrOH++ with log K values from 98Bal/Nor; added Cr(OH)3(aq) and
       Cr(OH)3(am) using log K from 98Bal/Nor.
 25. Adjusted formatting for compatibility with EQ3/6 v. 7 series code.
       Cleaned up comment material in the data blocks.
 26. Removed "500.000" values assigned to mineral volumes; assigned zeros.
 27. Replacement of "UO3:.9H2O(alpha)" by "Schoepite(dehyd,0.9)".
 28. Update of the phosphate solids whitlockite, hydroxyapatite, strengite,
     and fluorapatite.
 +--------------------------------------------------------------------

 Number of logK temperature grid ranges= 2

 Enthalpy functions flag= -1 (Not present)
 Volume functions flag  = -1 (Not present)


 element = O       , atwt =   15.99940
 element = Ag      , atwt =  107.86820
 element = Al      , atwt =   26.98154
 element = Am      , atwt =  243.00000
 element = Ar      , atwt =   39.94800
 element = Au      , atwt =  196.96655
 element = B       , atwt =   10.81100
 element = Ba      , atwt =  137.32700
 element = Be      , atwt =    9.01218
 element = Bi      , atwt =  208.98038
 element = Br      , atwt =   79.90400
 element = Ca      , atwt =   40.07800
 element = Cd      , atwt =  112.41100
 element = Ce      , atwt =  140.11600
 element = Cl      , atwt =   35.45270
 element = Co      , atwt =   58.93320
 element = Cr      , atwt =   51.99610
 element = Cs      , atwt =  132.90545
 element = Cu      , atwt =   63.54600
 element = Dy      , atwt =  162.50000
 element = Er      , atwt =  167.26000
 element = Eu      , atwt =  151.96400
 element = F       , atwt =   18.99840
 element = Fe      , atwt =   55.84500
 element = Fr      , atwt =  223.00000
 element = Ga      , atwt =   69.72300
 element = Gd      , atwt =  157.25000
 element = H       , atwt =    1.00794
 element = As      , atwt =   74.92160
 element = C       , atwt =   12.01070
 element = P       , atwt =   30.97376
 element = He      , atwt =    4.00260
 element = Hf      , atwt =  178.49000
 element = Hg      , atwt =  200.59000
 element = Ho      , atwt =  164.93032
 element = I       , atwt =  126.90447
 element = In      , atwt =  114.81800
 element = K       , atwt =   39.09830
 element = Kr      , atwt =   83.80000
 element = La      , atwt =  138.90550
 element = Li      , atwt =    6.94100
 element = Lu      , atwt =  174.96700
 element = Mg      , atwt =   24.30500
 element = Mn      , atwt =   54.93805
 element = Mo      , atwt =   95.94000
 element = N       , atwt =   14.00674
 element = Na      , atwt =   22.98977
 element = Nb      , atwt =   92.90638
 element = Nd      , atwt =  144.24000
 element = Ne      , atwt =   20.17970
 element = Ni      , atwt =   58.69340
 element = Np      , atwt =  237.00000
 element = Pb      , atwt =  207.20000
 element = Pd      , atwt =  106.42000
 element = Pm      , atwt =  145.00000
 element = Pr      , atwt =  140.90765
 element = Pt      , atwt =  195.07800
 element = Pu      , atwt =  244.00000
 element = Ra      , atwt =  226.00000
 element = Rb      , atwt =   85.46780
 element = Re      , atwt =  186.20700
 element = Rh      , atwt =  102.90550
 element = Rn      , atwt =  222.00000
 element = Ru      , atwt =  101.07000
 element = S       , atwt =   32.06600
 element = Sb      , atwt =  121.76000
 element = Sc      , atwt =   44.95591
 element = Se      , atwt =   78.96000
 element = Si      , atwt =   28.08550
 element = Sm      , atwt =  150.36000
 element = Sn      , atwt =  118.71000
 element = Sr      , atwt =   87.62000
 element = Tb      , atwt =  158.92534
 element = Tc      , atwt =   98.00000
 element = Th      , atwt =  232.03810
 element = Ti      , atwt =   47.86700
 element = Tl      , atwt =  204.38330
 element = Tm      , atwt =  168.93421
 element = U       , atwt =  238.02890
 element = V       , atwt =   50.94150
 element = W       , atwt =  183.84000
 element = Xe      , atwt =  131.29000
 element = Y       , atwt =   88.90585
 element = Yb      , atwt =  173.04000
 element = Zn      , atwt =   65.39000
 element = Zr      , atwt =   91.22400


 The minimum temperature is      0.000C
 The maximum temperature is    100.000C

 The maximum temperatures (C) by range are:
      100.000
      300.000


       presg
       adh
       bdh
       bdot
       cco2
       xlke


 aqueous

     1 H2O
     2 Ag+
     3 Al+++
     4 Am+++
     5 Ar(aq)
     6 Au+
     7 B(OH)3(aq)
     8 Ba++
     9 Be++
    10 Bi+++
    11 Br-
    12 Ca++
    13 Cd++
    14 Ce+++
    15 Cl-
    16 Co++
    17 CrO4--
    18 Cs+
    19 Cu++
    20 Dy+++
    21 Er+++
    22 Eu+++
    23 F-
    24 Fe++
    25 Fr+
    26 Ga+++
    27 Gd+++
    28 H+
    29 H2AsO4-
    30 HCO3-
    31 HPO4--
    32 He(aq)
    33 Hf++++
    34 Hg++
    35 Ho+++
    36 I-
    37 In+++
    38 K+
    39 Kr(aq)
    40 La+++
    41 Li+
    42 Lu+++
    43 Mg++
    44 Mn++
    45 MoO4--
    46 NH3(aq)
    47 Na+
    48 NbO3-
    49 Nd+++
    50 Ne(aq)
    51 Ni++
    52 NpO2+
    53 Pb++
    54 Pd++
    55 Pm+++
    56 Pr+++
    57 Pt++
    58 PuO2+
    59 Ra++
    60 Rb+
    61 ReO4-
    62 Rh+++
    63 Rn(aq)
    64 RuO4--
    65 SO4--
    66 SbO2-
    67 Sc+++
    68 SeO3--
    69 SiO2(aq)
    70 Sm+++
    71 Sn++
    72 Sr++
    73 Tb+++
    74 TcO4-
    75 Th++++
    76 Ti(OH)4(aq)
    77 Tl+
    78 Tm+++
    79 UO2++
    80 VO2+
    81 WO4--
    82 Xe(aq)
    83 Y+++
    84 Yb+++
    85 Zn++
    86 ZrO++
    87 O2(g)
    88 HS-
    89 Acetic_acid(aq)
    90 Formic_acid(aq)
    91 S2--
    92 Ag++
    93 Am++
    94 Am++++
    95 AmO2+
    96 AmO2++
    97 Au+++
    98 BH4-
    99 Br3-
   100 BrO-
   101 BrO3-
   102 BrO4-
   103 CN-
   104 Ce++
   105 Ce++++
   106 ClO-
   107 ClO2-
   108 ClO3-
   109 ClO4-
   110 Co+++
   111 Cr++
   112 Cr+++
   113 CrO4---
   114 Cu+
   115 Dy++
   116 Er++
   117 Ethane(aq)
   118 Eu++
   119 Fe+++
   120 Formaldehyde(aq)
   121 Gd++
   122 H2(aq)
   123 H2O2(aq)
   124 H2AsO3-
   125 H2P2O7--
   126 H2PO2-
   127 HPO3--
   128 HSO5-
   129 HSe-
   130 Hg2++
   131 Ho++
   132 I3-
   133 IO-
   134 IO3-
   135 IO4-
   136 La++
   137 Methane(aq)
   138 Methanol(aq)
   139 Mn+++
   140 MnO4--
   141 N2(aq)
   142 N2O2--
   143 N3-
   144 NO2-
   145 NO3-
   146 Nd++
   147 Np+++
   148 Np++++
   149 NpO2++
   150 O2(aq)
   151 Pm++
   152 Pr++
   153 Pu+++
   154 Pu++++
   155 PuO2++
   156 Rh++
   157 Ru(OH)2++
   158 Ru++
   159 Ru+++
   160 RuO4-
   161 RuO4(aq)
   162 S2O3--
   163 S2O4--
   164 S2O6--
   165 S2O8--
   166 S3--
   167 S3O6--
   168 S4--
   169 S4O6--
   170 S5--
   171 S5O6--
   172 SCN-
   173 SO3--
   174 SeO4--
   175 Sm++
   176 Sn++++
   177 Tb++
   178 TcO++
   179 TcO4--
   180 Tl+++
   181 Tm++
   182 U+++
   183 U++++
   184 UO2+
   185 V++
   186 V+++
   187 VO++
   188 Yb++
   189 Zr++++
   190 (NH4)2Sb2S4(aq)
   191 (NpO2)2(OH)2++
   192 (NpO2)2CO3(OH)3-
   193 (NpO2)3(CO3)6(-6)
   194 (NpO2)3(OH)5+
   195 (PuO2)2(OH)2++
   196 Pu(CO3)5------
   197 (PuO2)3(CO3)6(6-)
   198 (UO2)11(CO3)6(OH)12--
   199 (UO2)2(OH)2++
   200 (UO2)2(PuO2)(CO3)6(6-)
   201 (UO2)2CO3(OH)3-
   202 (UO2)2NpO2(CO3)6(-6)
   203 (UO2)2OH+++
   204 (UO2)3(CO3)6(6-)
   205 (UO2)3(OH)4++
   206 (UO2)3(OH)5+
   207 (UO2)3(OH)7-
   208 (UO2)3O(OH)2(HCO3)+
   209 (UO2)4(OH)7+
   210 (VO)2(OH)2++
   211 Ag(CO3)2---
   212 Ag(HS)2-
   213 AgCO3-
   214 AgCl(aq)
   215 AgCl2-
   216 AgCl3--
   217 AgCl4---
   218 AgF(aq)
   219 AgNO3(aq)
   220 AgO-
   221 AgOH(aq)
   222 Al(SO4)2-
   223 Al13O4(OH)24(7+)
   224 Al2(OH)2++++
   225 Al3(OH)4(5+)
   226 AlF++
   227 AlF2+
   228 AlF3(aq)
   229 AlF4-
   230 AlH2PO4++
   231 AlHPO4+
   232 AlO+
   233 AlO2-
   234 AlOH++
   235 AlSO4+
   236 Am(CO3)2-
   237 Am(CO3)3---
   238 Am(CO3)5(6-)
   239 AmO2(CO3)3----
   240 Am(OH)2+
   241 Am(OH)3(aq)
   242 Am(SO4)2-
   243 AmCO3+
   244 AmCl++
   245 AmF++
   246 AmF2+
   247 AmH2PO4++
   248 AmN3++
   249 AmNO2++
   250 AmNO3++
   251 AmOH++
   252 AmSO4+
   253 AsO3F--
   254 AsO4---
   255 Au(HS)2-
   256 AuCl(aq)
   257 AuCl2-
   258 AuCl3--
   259 AuCl4-
   260 B2O(OH)5-
   261 BF2(OH)2-
   262 BF3OH-
   263 BF4-
   264 BO2-
   265 BaB(OH)4+
   266 BaCO3(aq)
   267 BaCl+
   268 BaF+
   269 BaHCO3+
   270 BaNO3+
   271 BaOH+
   272 BeCl+
   273 BeCl2(aq)
   274 BeF+
   275 BeF2(aq)
   276 BeF3-
   277 BeF4--
   278 BeO(aq)
   279 BeO2--
   280 BeOH+
   281 BiO+
   282 BiO2-
   283 BiOH++
   284 Br2(aq)
   285 CO2(aq)
   286 CO3--
   287 CaB(OH)4+
   288 CaCO3(aq)
   289 CaCl+
   290 CaCl2(aq)
   291 CaCrO4(aq)
   292 CaF+
   293 CaHCO3+
   294 CaHPO4(aq)
   295 CaHSiO3+
   296 CaNO3+
   297 CaOH+
   298 CaP2O7--
   299 CaPO4-
   300 CaSO4(aq)
   301 Cd(CN)2(aq)
   302 Cd(CN)3-
   303 Cd(CN)4--
   304 Cd(CO3)2--
   305 Cd(N3)2(aq)
   306 Cd(N3)3-
   307 Cd(N3)4--
   308 Cd(NH3)++
   309 Cd(NH3)2++
   310 Cd(NH3)4++
   311 Cd(OH)Cl(aq)
   312 Cd(SCN)2(aq)
   313 Cd(SCN)3-
   314 Cd2OH+++
   315 Cd4(OH)4++++
   316 CdBr+
   317 CdBr2(aq)
   318 CdBr3-
   319 CdCN+
   320 CdCO3(aq)
   321 CdCl+
   322 CdCl2(aq)
   323 CdCl3-
   324 CdCl4--
   325 CdF+
   326 CdF2(aq)
   327 CdHCO3+
   328 CdI+
   329 CdI2(aq)
   330 CdI3-
   331 CdI4--
   332 CdN3+
   333 CdNO2+
   334 CdO(aq)
   335 CdO2--
   336 CdOH+
   337 CdP2O7--
   338 CdSCN+
   339 CdSO4(aq)
   340 CdSeO4(aq)
   341 Ce(CO3)2-
   342 Ce(HPO4)2-
   343 Ce(OH)2+
   344 Ce(OH)2++
   345 Ce(OH)3(aq)
   346 Ce(PO4)2---
   347 Ce2(OH)2(6+)
   348 Ce3(OH)5++++
   349 CeCO3+
   350 CeCl++
   351 CeF++
   352 CeF2+
   353 CeH2PO4++
   354 CeHCO3++
   355 CeHPO4+
   356 CeNO3++
   357 CeOH++
   358 CeOH+++
   359 CePO4(aq)
   360 CeSO4+
   361 Co2(OH)3+
   362 Co4(OH)4++++
   363 CoBr2(aq)
   364 CoCl+
   365 CoF+
   366 CoI2(aq)
   367 CoNO3+
   368 CoO(aq)
   369 CoO2--
   370 CoOH+
   371 CoOH++
   372 CoSO4(aq)
   373 CoSeO4(aq)
   374 Cr2O7--
   375 CrOH++
   376 Cr(OH)3(aq)
   377 CrBr++
   378 CrCl++
   379 CrCl2+
   380 CrO+
   381 CrO2-
   382 CrO3Cl-
   383 CsBr(aq)
   384 CsCl(aq)
   385 CsI(aq)
   386 CsOH(aq)
   387 Cu(NH3)2++
   388 Cu(NH3)3++
   389 CuCl(aq)
   390 CuCl+
   391 CuCl2(aq)
   392 CuCl2-
   393 CuCl3-
   394 CuCl3--
   395 CuCl4--
   396 CuF+
   397 CuNH3++
   398 CuO(aq)
   399 CuO2--
   400 CuOH+
   401 CuSO4(aq)
   402 Dy(CO3)2-
   403 Dy(HPO4)2-
   404 Dy(OH)2+
   405 Dy(OH)3(aq)
   406 Dy(OH)4-
   407 Dy(PO4)2---
   408 Dy(SO4)2-
   409 DyCO3+
   410 DyCl++
   411 DyF++
   412 DyH2PO4++
   413 DyHCO3++
   414 DyHPO4+
   415 DyNO3++
   416 DyOH++
   417 DyPO4(aq)
   418 DySO4+
   419 Er(CO3)2-
   420 Er(HPO4)2-
   421 Er(OH)2+
   422 Er(OH)3(aq)
   423 Er(OH)4-
   424 Er(PO4)2---
   425 Er(SO4)2-
   426 ErCO3+
   427 ErCl++
   428 ErF++
   429 ErH2PO4++
   430 ErHCO3++
   431 ErHPO4+
   432 ErNO3++
   433 ErOH++
   434 ErPO4(aq)
   435 ErSO4+
   436 Eu(CO3)2-
   437 Eu(HPO4)2-
   438 Eu(OH)2+
   439 Eu(OH)3(aq)
   440 Eu(PO4)2---
   441 Eu(SO4)2-
   442 EuBr++
   443 EuBr2+
   444 EuBrO3++
   445 EuCO3+
   446 EuCl++
   447 EuCl2+
   448 EuF++
   449 EuF2+
   450 EuF3(aq)
   451 EuH2PO4++
   452 EuHCO3++
   453 EuHPO4+
   454 EuIO3++
   455 EuNO3++
   456 EuOH++
   457 EuPO4(aq)
   458 EuSO4+
   459 Fe(CO3)2--
   460 Fe(OH)3-
   461 Fe(OH)4--
   462 Fe(SO4)2-
   463 Fe2(OH)2++++
   464 Fe3(OH)4(5+)
   465 FeCO3(aq)
   466 FeCl+
   467 FeCl++
   468 FeCl2(aq)
   469 FeF+
   470 FeF++
   471 FeF2+
   472 FeH2PO4+
   473 FeH2PO4++
   474 FeHCO3+
   475 FeNO2++
   476 FeNO3++
   477 FeO(aq)
   478 FeO+
   479 FeO2-
   480 FeOH+
   481 FeOH++
   482 FeSO4(aq)
   483 FeSO4+
   484 Formate
   485 GaO+
   486 GaO2-
   487 GaOH++
   488 Gd(CO3)2-
   489 Gd(HPO4)2-
   490 Gd(OH)2+
   491 Gd(OH)3(aq)
   492 Gd(OH)4-
   493 Gd(PO4)2---
   494 Gd(SO4)2-
   495 GdCO3+
   496 GdCl++
   497 GdF++
   498 GdF2+
   499 GdH2PO4++
   500 GdHCO3++
   501 GdHPO4+
   502 GdNO3++
   503 GdOH++
   504 GdPO4(aq)
   505 GdSO4+
   506 H2CrO4(aq)
   507 H2MoO4(aq)
   508 H2N2O2(aq)
   509 H2PO3-
   510 H2PO3F(aq)
   511 H2PO4-
   512 H2S(aq)
   513 H2S2O3(aq)
   514 H2S2O4(aq)
   515 H2SO3(aq)
   516 H2SO4(aq)
   517 H2Se(aq)
   518 H2SeO3(aq)
   519 H2VO4-
   520 H3AsO4(aq)
   521 H3P2O7-
   522 H3PO2(aq)
   523 H3PO3(aq)
   524 H3PO4(aq)
   525 H3VO4(aq)
   526 H4P2O7(aq)
   527 HAlO2(aq)
   528 HAsO2(aq)
   529 HAsO3F-
   530 HAsO4--
   531 HAsS2(aq)
   532 HBeO2-
   533 HBiO2(aq)
   534 HBrO(aq)
   535 HCdO2-
   536 HClO(aq)
   537 HClO2(aq)
   538 HCoO2-
   539 HCrO2(aq)
   540 HCrO4-
   541 HCuO2-
   542 HF(aq)
   543 HF2-
   544 HFeO2(aq)
   545 HFeO2-
   546 HGaO2(aq)
   547 HHfO2+
   548 HHfO3-
   549 HHgO2-
   550 HIO(aq)
   551 HIO3(aq)
   552 HInO2(aq)
   553 HMnO2-
   554 HMoO4-
   555 HN2O2-
   556 HN3(aq)
   557 HNO2(aq)
   558 HNO3(aq)
   559 HNbO3(aq)
   560 HNiO2-
   561 HO2-
   562 HP2O7---
   563 HPO2F2(aq)
   564 PO2F2-
   565 HPO3F-
   566 HPbO2-
   567 HRuO5-
   568 HS2O3-
   569 HS2O4-
   570 HSO3-
   571 HSO4-
   572 HSbO2(aq)
   573 HScO2(aq)
   574 HSeO3-
   575 HSeO4-
   576 CaSeO4(aq)
   577 HSiO3-
   578 HSnO2-
   579 HTlO2(aq)
   580 HUO2(aq)
   581 HUO2+
   582 HUO3-
   583 HUO4-
   584 HVO4--
   585 HWO4-
   586 HZnO2-
   587 HZrO2+
   588 HZrO3-
   589 HfO++
   590 HfO2(aq)
   591 HfOH+++
   592 HgCl+
   593 HgCl2(aq)
   594 HgCl3-
   595 HgCl4--
   596 HgF+
   597 HgO(aq)
   598 HgOH+
   599 Ho(CO3)2-
   600 Ho(HPO4)2-
   601 Ho(OH)2+
   602 Ho(OH)3(aq)
   603 Ho(PO4)2---
   604 Ho(SO4)2-
   605 HoCO3+
   606 HoF++
   607 HoH2PO4++
   608 HoHCO3++
   609 HoHPO4+
   610 HoNO3++
   611 HoOH++
   612 HoPO4(aq)
   613 HoSO4+
   614 InCl++
   615 InF++
   616 InO+
   617 InO2-
   618 InOH++
   619 KBr(aq)
   620 KCl(aq)
   621 KHSO4(aq)
   622 KI(aq)
   623 KOH(aq)
   624 KP2O7---
   625 KSO4-
   626 La(CO3)2-
   627 La(HPO4)2-
   628 La(OH)2+
   629 La(OH)3(aq)
   630 La(PO4)2---
   631 La(SO4)2-
   632 La2(OH)2++++
   633 La5(OH)9(6+)
   634 LaCO3+
   635 LaCl++
   636 LaF++
   637 LaF2+
   638 LaF3(aq)
   639 LaH2PO4++
   640 LaHCO3++
   641 LaHPO4+
   642 LaNO3++
   643 LaOH++
   644 LaPO4(aq)
   645 LaSO4+
   646 LiCl(aq)
   647 LiOH(aq)
   648 LiSO4-
   649 Lu(CO3)2-
   650 Lu(HPO4)2-
   651 Lu(OH)2+
   652 Lu(OH)3(aq)
   653 Lu(PO4)2---
   654 Lu(SO4)2-
   655 LuCO3+
   656 LuCl++
   657 LuF++
   658 LuH2PO4++
   659 LuHCO3++
   660 LuHPO4+
   661 LuNO3++
   662 LuOH++
   663 LuPO4(aq)
   664 LuSO4+
   665 Mg4(OH)4++++
   666 MgB(OH)4+
   667 MgCO3(aq)
   668 MgCl+
   669 MgF+
   670 MgHCO3+
   671 MgHPO4(aq)
   672 MgHSiO3+
   673 MgOH+
   674 MgP2O7--
   675 MgSO4(aq)
   676 Mn(NO3)2(aq)
   677 Mn2(OH)3+
   678 Mn2OH+++
   679 MnCl+
   680 MnCl3-
   681 MnF+
   682 MnHCO3+
   683 MnNO3+
   684 MnO(aq)
   685 MnO2--
   686 MnO4-
   687 MnOH+
   688 MnSO4(aq)
   689 MnSeO4(aq)
   690 N2H5+
   691 N2H6++
   692 NH4+
   693 NH4SO4-
   694 NH4SbO2(aq)
   695 Na2P2O7--
   696 NaB(OH)4(aq)
   697 NaBr(aq)
   698 NaCO3-
   699 NaCl(aq)
   700 NaF(aq)
   701 NaHCO3(aq)
   702 NaHP2O7--
   703 NaHSiO3(aq)
   704 NaI(aq)
   705 NaOH(aq)
   706 NaP2O7---
   707 NaSO4-
   708 Nd(CO3)2-
   709 Nd(HPO4)2-
   710 Nd(OH)2+
   711 Nd(OH)3(aq)
   712 Nd(OH)4-
   713 Nd(PO4)2---
   714 Nd(SO4)2-
   715 Nd2(OH)2++++
   716 NdCO3+
   717 NdCl++
   718 NdF++
   719 NdH2PO4++
   720 NdHCO3++
   721 NdHPO4+
   722 NdNO3++
   723 NdOH++
   724 NdPO4(aq)
   725 NdSO4+
   726 NiCO3(aq)
   727 NiCrO4(aq)
   728 Ni(NH3)2++
   729 Ni(NH3)6++
   730 Ni(NO3)2(aq)
   731 Ni(OH)2(aq)
   732 Ni(OH)3-
   733 NiBr+
   734 NiCl+
   735 NiF+
   736 NiHP2O7-
   737 NiNO3+
   738 NiO(aq)
   739 NiO2--
   740 NiOH+
   741 NiP2O7--
   742 NiSeO4(aq)
   743 Np(CO3)3---
   744 Np(CO3)4----
   745 Np(CO3)5(6-)
   746 Np(OH)4(aq)
   747 Np(SO4)2(aq)
   748 Np(SCN)+++
   749 Np(SCN)2++
   750 Np(SCN)3+
   751 NpCl+++
   752 NpF+++
   753 NpF2++
   754 NpI+++
   755 NpNO3+++
   756 NpO2(CO3)2--
   757 NpO2(CO3)2---
   758 NpO2(CO3)2OH----
   759 NpO2(CO3)3(5-)
   760 NpO2(CO3)3----
   761 NpO2(HPO4)2--
   762 NpO2(OH)2-
   763 NpO2(SO4)2--
   764 NpO2Cl+
   765 NpO2CO3-
   766 NpO2CO3(aq)
   767 NpO2F(aq)
   768 NpO2F+
   769 NpO2F2(aq)
   770 NpO2H2PO4+
   771 NpO2HPO4-
   772 NpO2HPO4(aq)
   773 NpO2IO3(aq)
   774 NpO2IO3+
   775 NpO2OH(aq)
   776 NpO2OH+
   777 NpO2SO4(aq)
   778 NpO2SO4-
   779 NpOH++
   780 NpOH+++
   781 NpSO4++
   782 OH-
   783 P2O7----
   784 PH4+
   785 PO3F--
   786 PO4---
   787 Pb(BrO3)2(aq)
   788 Pb(ClO3)2(aq)
   789 Pb(HS)2(aq)
   790 Pb(HS)3-
   791 Pb(OH)2(aq)
   792 Pb(OH)3-
   793 Pb(SCN)2(aq)
   794 Pb2OH+++
   795 Pb3(OH)4++
   796 Pb4(OH)4++++
   797 Pb6(OH)8++++
   798 PbBr+
   799 PbBr2(aq)
   800 PbBr3-
   801 PbBrO3+
   802 PbCl+
   803 PbCl2(aq)
   804 PbCl3-
   805 PbCl4--
   806 PbClO3+
   807 PbF+
   808 PbF2(aq)
   809 PbH2PO4+
   810 PbHPO4(aq)
   811 PbI+
   812 PbI2(aq)
   813 PbI3-
   814 PbI4--
   815 PbNO3+
   816 PbO(aq)
   817 PbOH+
   818 PbP2O7--
   819 PbSCN+
   820 Pd(SO4)2--
   821 Pd(SO4)3----
   822 PdCl+
   823 PdCl2(aq)
   824 PdCl3-
   825 PdCl4--
   826 PdO(aq)
   827 PdOH+
   828 PdSO4(aq)
   829 Pm(CO3)2-
   830 Pm(HPO4)2-
   831 Pm(OH)2+
   832 Pm(OH)3(aq)
   833 Pm(PO4)2---
   834 Pm(SO4)2-
   835 PmCO3+
   836 PmCl++
   837 PmF++
   838 PmH2PO4++
   839 PmHCO3++
   840 PmHPO4+
   841 PmNO3++
   842 PmOH++
   843 PmPO4(aq)
   844 PmSO4+
   845 Pr(CO3)2-
   846 Pr(HPO4)2-
   847 Pr(OH)2+
   848 Pr(OH)3(aq)
   849 Pr(PO4)2---
   850 Pr(SO4)2-
   851 PrCO3+
   852 PrCl++
   853 PrF++
   854 PrH2PO4++
   855 PrHCO3++
   856 PrHPO4+
   857 PrNO3++
   858 PrOH++
   859 PrPO4(aq)
   860 PrSO4+
   861 Pt(SO4)2--
   862 Pt(SO4)3----
   863 PtCl+
   864 PtCl2(aq)
   865 PtCl3-
   866 PtCl4--
   867 PtO(aq)
   868 PtOH+
   869 PtSO4(aq)
   870 Pu(SO4)2(aq)
   871 Pu(SO4)2-
   872 PuBr+++
   873 PuCl+++
   874 PuF+++
   875 PuF2++
   876 PuH3PO4++++
   877 PuI++
   878 PuNO3+++
   879 PuO2(CO3)2--
   880 PuO2(CO3)3----
   881 PuO2(CO3)3(5-)
   882 PuO2(OH)2(aq)
   883 PuO2(SO4)2--
   884 PuO2Cl+
   885 PuO2Cl2(aq)
   886 PuO2CO3(aq)
   887 PuO2CO3-
   888 PuO2F+
   889 PuO2F2(aq)
   890 PuO2OH(aq)
   891 PuO2OH+
   892 PuO2SO4(aq)
   893 PuOH++
   894 PuOH+++
   895 PuSCN++
   896 PuSO4+
   897 PuSO4++
   898 RbBr(aq)
   899 RbCl(aq)
   900 RbF(aq)
   901 RbI(aq)
   902 RbOH(aq)
   903 Rh(SO4)2-
   904 Rh(SO4)2--
   905 Rh(SO4)3---
   906 Rh(SO4)3----
   907 RhCl+
   908 RhCl++
   909 RhCl2(aq)
   910 RhCl2+
   911 RhCl3(aq)
   912 RhCl3-
   913 RhCl4-
   914 RhCl4--
   915 RhO(aq)
   916 RhO+
   917 RhOH+
   918 RhOH++
   919 RhSO4(aq)
   920 RhSO4+
   921 Ru(Cl)2+
   922 Ru(Cl)3(aq)
   923 Ru(OH)2+
   924 Ru(OH)2Cl+
   925 Ru(OH)2Cl2(aq)
   926 Ru(OH)2Cl3-
   927 Ru(OH)2Cl4--
   928 Ru(OH)2SO4(aq)
   929 Ru(OH)4(aq)
   930 Ru(SO4)2-
   931 Ru(SO4)2--
   932 Ru(SO4)3---
   933 Ru(SO4)3----
   934 Ru4(OH)12++++
   935 RuCl+
   936 RuCl++
   937 RuCl2(aq)
   938 RuCl2+
   939 RuCl3(aq)
   940 RuCl3-
   941 RuCl4-
   942 RuCl4--
   943 RuCl5--
   944 RuCl6---
   945 RuO(aq)
   946 RuO+
   947 RuOH+
   948 RuOH++
   949 RuSO4(aq)
   950 RuSO4+
   951 S--
   952 S2O5--
   953 SO2(aq)
   954 ScO+
   955 ScO2-
   956 ScOH++
   957 SiF6--
   958 Sm(CO3)2-
   959 Sm(HPO4)2-
   960 Sm(OH)2+
   961 Sm(OH)3(aq)
   962 Sm(OH)4-
   963 Sm(PO4)2---
   964 Sm(SO4)2-
   965 SmCO3+
   966 SmCl++
   967 SmF++
   968 SmH2PO4++
   969 SmHCO3++
   970 SmHPO4+
   971 SmNO3++
   972 SmOH++
   973 SmPO4(aq)
   974 SmSO4+
   975 Sn(OH)2++
   976 Sn(OH)3+
   977 Sn(OH)4(aq)
   978 Sn(OH)5-
   979 Sn(OH)6--
   980 Sn(SO4)2(aq)
   981 SnCl+
   982 SnCl2(aq)
   983 SnCl3-
   984 SnF+
   985 SnF2(aq)
   986 SnF3-
   987 SnO(aq)
   988 SnOH+
   989 SnOH+++
   990 SnSO4++
   991 SrCO3(aq)
   992 SrCl+
   993 SrF+
   994 SrHCO3+
   995 SrHPO4(aq)
   996 SrNO3+
   997 SrOH+
   998 SrP2O7--
   999 SrSO4(aq)
  1000 Tb(CO3)2-
  1001 Tb(HPO4)2-
  1002 Tb(OH)2+
  1003 Tb(OH)3(aq)
  1004 Tb(PO4)2---
  1005 Tb(SO4)2-
  1006 TbCO3+
  1007 TbCl++
  1008 TbF++
  1009 TbH2PO4++
  1010 TbHCO3++
  1011 TbHPO4+
  1012 TbNO3++
  1013 TbOH++
  1014 TbPO4(aq)
  1015 TbSO4+
  1016 TcCO3(OH)2(aq)
  1017 TcCO3(OH)3-
  1018 TcO(OH)2(aq)
  1019 TcO(OH)3-
  1020 TcOOH+
  1021 Th(OH)2++
  1022 Th(OH)3+
  1023 Th(OH)(CO3)4(5-)
  1024 Th(OH)2CO3(aq)
  1025 Th(OH)2(CO3)2--
  1026 Th(OH)3CO3-
  1027 Th(OH)4CO3--
  1028 Th(OH)4(aq)
  1029 Th(OH)4PO4---
  1030 Th(SO4)2(aq)
  1031 Th(SO4)3--
  1032 Th(SO4)4----
  1033 Th2(OH)2(6+)
  1034 Th2(OH)3(5+)
  1035 ThCl+++
  1036 ThCl2++
  1037 ThCl3+
  1038 ThCl4(aq)
  1039 ThF+++
  1040 ThF2++
  1041 ThF3+
  1042 ThF4(aq)
  1043 Th(OH)+++
  1044 ThSO4++
  1045 Ti(OH)3+
  1046 Ti(OH)5-
  1047 TlCl(aq)
  1048 TlCl++
  1049 TlF(aq)
  1050 TlO+
  1051 TlO2-
  1052 TlOH(aq)
  1053 TlOH++
  1054 Tm(CO3)2-
  1055 Tm(HPO4)2-
  1056 Tm(OH)2+
  1057 Tm(OH)3(aq)
  1058 Tm(PO4)2---
  1059 Tm(SO4)2-
  1060 TmCO3+
  1061 TmCl++
  1062 TmF++
  1063 TmH2PO4++
  1064 TmHCO3++
  1065 TmHPO4+
  1066 TmNO3++
  1067 TmOH++
  1068 TmPO4(aq)
  1069 TmSO4+
  1070 U(CO3)4----
  1071 U(CO3)5(6-)
  1072 U(NO3)2++
  1073 U(SCN)2++
  1074 U(SO4)2(aq)
  1075 UBr+++
  1076 UCl+++
  1077 UF+++
  1078 UF2++
  1079 UF3+
  1080 UF4(aq)
  1081 UF5-
  1082 UF6--
  1083 UI+++
  1084 UNO3+++
  1085 UO+
  1086 UO++
  1087 UO2(CO3)2--
  1088 UO2(CO3)3(5-)
  1089 UO2(CO3)3----
  1090 UO2(H2PO4)(H3PO4)+
  1091 UO2(H2PO4)2(aq)
  1092 UO2(IO3)2(aq)
  1093 UO2(N3)2(aq)
  1094 UO2(N3)3-
  1095 UO2(N3)4--
  1096 UO2(SCN)2(aq)
  1097 UO2(SCN)3-
  1098 UO2(SO4)2--
  1099 UO2(aq)
  1100 UO2-
  1101 UO2Br+
  1102 UO2BrO3+
  1103 UO2CO3(aq)
  1104 UO2Cl+
  1105 UO2Cl2(aq)
  1106 UO2ClO3+
  1107 UO2F+
  1108 UO2F2(aq)
  1109 UO2F3-
  1110 UO2F4--
  1111 UO2H2PO4+
  1112 UO2H3PO4++
  1113 UO2HPO4(aq)
  1114 UO2IO3+
  1115 UO2N3+
  1116 UO2NO3+
  1117 UO2OH(aq)
  1118 UO2OH+
  1119 UO2OSi(OH)3+
  1120 UO2PO4-
  1121 UO2S2O3(aq)
  1122 UO2SCN+
  1123 UO2SO3(aq)
  1124 UO2SO4(aq)
  1125 UO3(aq)
  1126 UO3-
  1127 UO4--
  1128 UOH++
  1129 UOH+++
  1130 USCN+++
  1131 USO4++
  1132 V(OH)2+
  1133 V2(OH)2++++
  1134 VO(OH)3(aq)
  1135 VO+
  1136 VO2(HPO4)2---
  1137 VO2F(aq)
  1138 VO2F2-
  1139 VO2H2PO4(aq)
  1140 VO2HPO4-
  1141 VO2SO4-
  1142 VO4---
  1143 VOF+
  1144 VOF2(aq)
  1145 VOH+
  1146 VOH++
  1147 VOOH+
  1148 VOSO4(aq)
  1149 VSO4+
  1150 Y(CO3)2-
  1151 Y(HPO4)2-
  1152 Y(OH)2+
  1153 Y(OH)3(aq)
  1154 Y(OH)4-
  1155 Y(PO4)2---
  1156 Y(SO4)2-
  1157 Y2(OH)2++++
  1158 YCO3+
  1159 YCl++
  1160 YF++
  1161 YF2+
  1162 YF3(aq)
  1163 YH2PO4++
  1164 YHCO3++
  1165 YHPO4+
  1166 YNO3++
  1167 YOH++
  1168 YPO4(aq)
  1169 YSO4+
  1170 Yb(CO3)2-
  1171 Yb(HPO4)2-
  1172 Yb(OH)2+
  1173 Yb(OH)3(aq)
  1174 Yb(OH)4-
  1175 Yb(PO4)2---
  1176 Yb(SO4)2-
  1177 YbCO3+
  1178 YbCl++
  1179 YbF++
  1180 YbH2PO4++
  1181 YbHCO3++
  1182 YbHPO4+
  1183 YbNO3++
  1184 YbOH++
  1185 YbPO4(aq)
  1186 YbSO4+
  1187 Zn(CN)4--
  1188 Zn(N3)2(aq)
  1189 Zn(NH3)++
  1190 Zn(NH3)2++
  1191 Zn(NH3)3++
  1192 Zn(NH3)4++
  1193 Zn(OH)Cl(aq)
  1194 Zn(SCN)2(aq)
  1195 Zn(SCN)4--
  1196 ZnBr+
  1197 ZnBr2(aq)
  1198 ZnBr3-
  1199 ZnCO3(aq)
  1200 ZnCl+
  1201 ZnCl2(aq)
  1202 ZnCl3-
  1203 ZnClO4+
  1204 ZnF+
  1205 ZnH2PO4+
  1206 ZnHCO3+
  1207 ZnHPO4(aq)
  1208 ZnI+
  1209 ZnI2(aq)
  1210 ZnI3-
  1211 ZnI4--
  1212 ZnN3+
  1213 ZnO(aq)
  1214 ZnO2--
  1215 ZnOH+
  1216 ZnSO4(aq)
  1217 ZnSeO4(aq)
  1218 ZrO2(aq)
  1219 ZrOH+++


 minerals

     1 (C4AF)
     2 (C4AH13)
     3 (CAH10)
     4 (C2AH8)
     5 (C4AH19)
     6 (CA2)
     7 (CA)
     8 (C3A)
     9 (C12A7)
    10 (NH4)4NpO2(CO3)3
    11 (UO2)2As2O7
    12 (UO2)2Cl3
    13 (UO2)2P2O7
    14 (UO2)3(AsO4)2
    15 (UO2)3(PO4)2
    16 (UO2)3(PO4)2:4H2O
    17 (UO2)3(PO4)2:6H2O
    18 (VO)3(PO4)2
    19 Acanthite
    20 Afwillite
    21 Silver
    22 Ag3PO4
    23 AgTcO4
    24 Akermanite
    25 Aluminum
    26 Al2(SO4)3
    27 Al2(SO4)3:6H2O
    28 AlF3
    29 Alabandite
    30 Alamosite
    31 Albite
    32 Albite_high
    33 Albite_low
    34 Allite (C3S)
    35 Alstonite
    36 Alum-K
    37 Alunite
    38 Americium
    39 Am(OH)3
    40 Am(OH)3(am)
    41 Am2(CO3)3
    42 Am2C3
    43 Am2O3
    44 AmBr3
    45 AmCl3
    46 AmF3
    47 AmF4
    48 AmH2
    49 AmI3
    50 AmO2
    51 AmOBr
    52 AmOCl
    53 AmOHCO3
    54 AmPO4(am)
    55 Amesite-14A
    56 Amesite-7A
    57 Analcime
    58 Analcime-dehy
    59 Anatase
    60 Andalusite
    61 Andradite
    62 Anglesite
    63 Anhydrite
    64 Annite
    65 Anorthite
    66 Antarcticite
    67 Anthophyllite
    68 Antigorite
    69 Antigorite(am)
    70 Antlerite
    71 Aphthitalite
    72 Aragonite
    73 Arcanite
    74 Arsenolite
    75 Arsenopyrite
    76 Artinite
    77 Arsenic
    78 As2O5
    79 As4O6(cubi)
    80 As4O6(mono)
    81 Gold
    82 Azurite
    83 Boron
    84 B2O3
    85 Barium
    86 Ba(OH)2:8H2O
    87 Ba2Si3O8
    88 Ba2SiO4
    89 Ba2U2O7
    90 Ba3UO6
    91 BaBr2
    92 BaBr2:2H2O
    93 BaCl2
    94 BaCl2:2H2O
    95 BaCl2:H2O
    96 BaCrO4
    97 BaHPO4
    98 BaI2
    99 BaMnO4
   100 BaO
   101 BaS
   102 BaSeO3
   103 BaSeO4
   104 BaSiF6
   105 BaU2O7
   106 BaUO4
   107 Baddeleyite
   108 Barite
   109 Barytocalcite
   110 Bassanite
   111 Beryllium
   112 Be13U
   113 Becquerelite
   114 Beidellite-Ca
   115 Beidellite-H
   116 Beidellite-K
   117 Beidellite-Mg
   118 Beidellite-Na
   119 Bellite (C2S)
   120 Berlinite
   121 Berndtite
   122 Bieberite
   123 Bischofite
   124 Bixbyite
   125 Bloedite
   126 Boehmite
   127 Boltwoodite
   128 Boltwoodite-Na
   129 Borax
   130 Boric_acid
   131 Bornite
   132 Brezinaite
   133 Bromellite
   134 Brucite
   135 Brushite
   136 Bunsenite
   137 Burkeite
   138 Graphite
   139 Calcium
   140 Ca-Al_Pyroxene
   141 Ca2Al2O5:8H2O
   142 Ca2Cl2(OH)2:H2O
   143 Ca2V2O7
   144 Ca3(AsO4)2
   145 Ca3Al2O6
   146 Ca3V2O8
   147 Ca4Al2Fe2O10
   148 Ca4Al2O7:13H2O
   149 Ca4Al2O7:19H2O
   150 Ca4Cl2(OH)6:13H2O
   151 CaAl2O4
   152 CaAl2O4:10H2O
   153 CaAl4O7
   154 CaHfO3
   155 CaSO4:0.5H2O(beta)
   156 CaSeO4:2H2O
   157 CaUO4
   158 CaV2O6
   159 CaZrO3
   160 Cadmoselite
   161 Calcite
   162 Calomel
   163 Carnallite
   164 Cassiterite
   165 SnO2(am)
   166 Cattierite
   167 Cadmium
   168 Cd(BO2)2
   169 Cd(IO3)2
   170 Cd(OH)2
   171 Cd(OH)Cl
   172 Cd3(AsO4)2
   173 Cd3(PO4)2
   174 Cd3(SO4)(OH)4
   175 Cd3(SO4)2(OH)2
   176 CdBr2
   177 CdBr2:4H2O
   178 CdCl2
   179 CdCl2(NH3)2
   180 CdCl2(NH3)4
   181 CdCl2(NH3)6
   182 CdCl2:H2O
   183 CdCr2O4
   184 CdF2
   185 CdI2
   186 CdS
   187 CdSO4
   188 CdSO4:2.667H2O
   189 CdSO4:H2O
   190 CdSeO3
   191 CdSeO4
   192 CdSiO3
   193 Cerium
   194 Ce(OH)3
   195 Ce(OH)3(am)
   196 Ce2(CO3)3:8H2O
   197 Ce2O3
   198 Ce3(PO4)4
   199 CeF3:.5H2O
   200 CeO2
   201 CePO4:H2O
   202 Celadonite
   203 Celestite
   204 Cerussite
   205 Chabazite
   206 Chalcanthite
   207 Chalcedony
   208 Chalcocite
   209 Chalcocyanite
   210 Chalcopyrite
   211 Chamosite-7A
   212 Chlorargyrite
   213 Chloromagnesite
   214 Chromite
   215 Chrysotile
   216 Cinnabar
   217 Claudetite
   218 Clausthalite
   219 Clinochlore-14A
   220 Clinochlore-7A
   221 Clinoptilolite
   222 Clinoptilolite-Ca
   223 Clinoptilolite-Cs
   224 Clinoptilolite-K
   225 Clinoptilolite-NH4
   226 Clinoptilolite-Na
   227 Clinoptilolite-Sr
   228 Clinoptilolite-dehy
   229 Clinozoisite
   230 Cobalt
   231 Co(OH)2
   232 Co2SiO4
   233 Co3(AsO4)2
   234 Co3(PO4)2
   235 CoCl2
   236 CoCl2:2H2O
   237 CoCl2:6H2O
   238 CoCr2O4
   239 CoF2
   240 CoF3
   241 CoFe2O4
   242 CoHPO4
   243 CoO
   244 CoSO4
   245 CoSO4:3Co(OH)2
   246 CoSO4:6H2O
   247 CoSeO3
   248 CoTiO3
   249 CoWO4
   250 Coesite
   251 Coffinite
   252 Colemanite
   253 Compreignacite
   254 Cordierite_anhyd
   255 Cordierite_hydr
   256 Corundum
   257 Cotunnite
   258 Covellite
   259 Chromium
   260 CrCl3
   261 CrF3
   262 CrF4
   263 CrI3
   264 CrO2
   265 CrO3
   266 CrS
   267 Cr-ettringite
   268 Cr-ferrihydrite
   269 Cristobalite(alpha)
   270 Cristobalite(beta)
   271 Crocoite
   272 Cronstedtite-7A
   273 Cesium
   274 Cs2NaAmCl6
   275 Cs2NaPuCl6
   276 Cs2NpBr6
   277 Cs2NpCl6
   278 Cs2PuBr6
   279 Cs2PuCl6
   280 Cs2U2O7
   281 Cs2U4O12
   282 Cs2UO4
   283 Cs3PuCl6
   284 CsPu2Cl7
   285 CsTcO4
   286 CSH:1.7
   287 Copper
   288 Cu3(PO4)2
   289 CuCl2
   290 CuCr2O4
   291 CuSeO3
   292 Cuprite
   293 Daphnite-14A
   294 Daphnite-7A
   295 Dawsonite
   296 Diaspore
   297 Dicalcium_silicate
   298 Diopside
   299 Dolomite
   300 Dolomite-dis
   301 Dolomite-ord
   302 Downeyite
   303 Dysprosium
   304 Dy(OH)3
   305 Dy(OH)3(am)
   306 Dy2(CO3)3
   307 Dy2O3
   308 DyF3:.5H2O
   309 DyPO4:2H2O
   310 Enstatite
   311 Epidote
   312 Epidote-ord
   313 Epsomite
   314 Erbium
   315 Er(OH)3
   316 Er(OH)3(am)
   317 Er2(CO3)3
   318 Er2O3
   319 ErF3:.5H2O
   320 ErPO4:2H2O
   321 Erionite
   322 Eskolaite
   323 Cr(OH)3(am)
   324 Ettringite
   325 Europium
   326 Eu(IO3)3:2H2O
   327 Eu(NO3)3:6H2O
   328 Eu(OH)2.5Cl.5
   329 Eu(OH)2Cl
   330 Eu(OH)3
   331 Eu2(CO3)3:3H2O
   332 Eu2(SO4)3:8H2O
   333 Eu2O3(cubic)
   334 Eu2O3(monoclinic)
   335 Eu3O4
   336 EuBr3
   337 EuCl2
   338 EuCl3:6H2O
   339 EuF3:0.5H2O
   340 EuO
   341 EuOCl
   342 EuOHCO3
   343 EuPO4:H2O
   344 EuS
   345 EuSO4
   346 Eucryptite
   347 Fayalite
   348 Iron
   349 Fe(OH)2
   350 Fe(OH)3
   351 Fe2(MoO4)3
   352 Fe2(SO4)3
   353 FeF2
   354 FeF3
   355 FeO
   356 FeSO4
   357 FeV2O4
   358 Ferrite-Ca
   359 Ferrite-Cu
   360 Ferrite-Dicalcium
   361 Ferrite-Mg
   362 Ferrite-Ni
   363 Ferrite-Zn
   364 Ferroaluminoceladonite
   365 Ferroceladonite
   366 Ferroselite
   367 Ferrosilite
   368 Fluorapatite
   369 Fluorite
   370 Forsterite
   371 Foshagite
   372 Frankdicksonite
   373 Freboldite
   374 Friedl_salt
   375 Gallium
   376 Galena
   377 Gaylussite
   378 Gadolinium
   379 Gd(OH)3
   380 Gd(OH)3(am)
   381 Gd2(CO3)3
   382 Gd2O3
   383 GdF3:.5H2O
   384 GdPO4:2H2O
   385 Gehlenate_Hydrate
   386 Gehlenite
   387 Gibbsite
   388 Gismondine-Na
   389 Gismondine-Ca
   390 Glauberite
   391 Goethite
   392 Greenalite
   393 Grossular
   394 Gypsum
   395 Gyrolite
   396 Halite
   397 Hatrurite
   398 Hausmannite
   399 Heazlewoodite
   400 Hedenbergite
   401 Hematite
   402 Hemicarboaluminate
   403 Hercynite
   404 Herzenbergite
   405 Heulandite
   406 Hexahydrite
   407 Hafnium
   408 HfB2
   409 HfBr2
   410 HfBr4
   411 HfC
   412 HfCl2
   413 HfCl4
   414 HfF2
   415 HfF4
   416 HfI2
   417 HfI4
   418 HfN
   419 HfO2
   420 HfS2
   421 HfS3
   422 Hg2SO4
   423 Hg2SeO3
   424 HgSeO3
   425 Hillebrandite
   426 Holmium
   427 Ho(OH)3
   428 Ho(OH)3(am)
   429 Ho2(CO3)3
   430 Ho2O3
   431 HoF3:.5H2O
   432 HoPO4:2H2O
   433 Hopeite
   434 Huntite
   435 Hydroboracite
   436 Hydrogarnet
   437 Hydromagnesite
   438 Hydrophilite
   439 Hydrotalcite
   440 Hydroxylapatite
   441 Hydrozincite
   442 Iodine
   443 Ice
   444 Illite
   445 Ilmenite
   446 Indium
   447 Jadeite
   448 Jarosite
   449 Jarosite-Na
   450 Potassium
   451 K-Feldspar
   452 K2CO3:1.5H2O
   453 K2O
   454 K2Se
   455 K2UO4
   456 K3H(SO4)2
   457 K4NpO2(CO3)3
   458 K8H4(CO3)6:3H2O
   459 KAl(SO4)2
   460 KBr
   461 KMgCl3:2H2O
   462 KNaCO3:6H2O
   463 KTcO4
   464 Kainite
   465 Kalicinite
   466 Kalsilite
   467 Kaolinite
   468 Karelianite
   469 Katoite
   470 Kieserite
   471 Klockmannite
   472 Krutaite
   473 Kyanite
   474 Lanthanum
   475 La(OH)3
   476 La(OH)3(am)
   477 La2(CO3)3:8H2O
   478 La2O3
   479 LaCl3
   480 LaCl3:7H2O
   481 LaF3:.5H2O
   482 LaPO4:H2O
   483 Lammerite
   484 Lanarkite
   485 Lansfordite
   486 Larnite
   487 Laumontite
   488 Laurite
   489 Lawrencite
   490 Lawsonite
   491 Leonite
   492 Lithium
   493 Li2Se
   494 Li2UO4
   495 Lime
   496 Linnaeite
   497 Litharge
   498 Lopezite
   499 Lutetium
   500 Lu(OH)3
   501 Lu(OH)3(am)
   502 Lu2(CO3)3
   503 Lu2O3
   504 LuF3:.5H2O
   505 LuPO4:0.5H2O
   506 Magnesiochromite
   507 Magnesite
   508 Magnetite
   509 Malachite
   510 Manganite
   511 Manganosite
   512 Margarite
   513 Massicot
   514 Maximum_Microcline
   515 Mayenite
   516 Melanterite
   517 Mercallite
   518 Merwinite
   519 Mesolite
   520 Metacinnabar
   521 Magnesium
   522 Mg2V2O7
   523 MgBr2
   524 MgBr2:6H2O
   525 MgCl2:2H2O
   526 MgCl2:4H2O
   527 MgCl2:H2O
   528 MgOHCl
   529 MgSO4
   530 MgSeO3
   531 MgUO4
   532 MgV2O6
   533 Millerite
   534 Minium
   535 Minnesotaite
   536 Mirabilite
   537 Misenite
   538 Manganese
   539 Mn(OH)2(am)
   540 MnCl2:2H2O
   541 MnCl2:4H2O
   542 MnCl2:H2O
   543 MnO2(gamma)
   544 MnSO4
   545 MnSe
   546 MnSeO3
   547 MnV2O6
   548 Molybdenum
   549 MoO2Cl2
   550 MoSe2
   551 Molysite
   552 Monocarboaluminate
   553 Monohydrocalcite
   554 Monosulphate
   555 Monteponite
   556 Monticellite
   557 Montmorillonite-Ca
   558 Montmorillonite-H
   559 Montmorillonite-K
   560 Montmorillonite-Mg
   561 Montmorillonite-Na
   562 Montroydite
   563 Mordenite
   564 Mordenite-dehy
   565 Morenosite
   566 Muscovite
   567 NH4HSe
   568 NH4TcO4
   569 Sodium
   570 Na2CO3
   571 Na2CO3:7H2O
   572 Na2Cr2O7
   573 Na2CrO4
   574 Na2O
   575 Na2Se
   576 Na2Se2
   577 Na2SiO3
   578 Na2U2O7
   579 Na2UO4(alpha)
   580 Na3H(SO4)2
   581 Kogarkoite
   582 Na3NpF8
   583 Na3NpO2(CO3)2
   584 Na3UO4
   585 Na4Ca(SO4)3:2H2O
   586 Na4SiO4
   587 Na4UO2(CO3)3
   588 Na6Si2O7
   589 NaBr
   590 NaBr:2H2O
   591 NaFeO2
   592 NaNpO2CO3
   593 NaNpO2CO3:3.5H2O
   594 NaTcO4:4H2O
   595 NaUO3
   596 Nahcolite
   597 Nantokite
   598 Natrolite
   599 Natron
   600 Natrosilite
   601 Naumannite
   602 Neodymium
   603 Nd(OH)3
   604 Nd(OH)3(am)
   605 Nd(OH)3(c)
   606 Nd2(CO3)3
   607 Nd2O3
   608 NdF3:.5H2O
   609 NdOHCO3
   610 NdPO4:H2O
   611 Nepheline
   612 Nesquehonite
   613 Nickel
   614 NiCO3
   615 NiCO3:5.5H2O
   616 Ni(OH)2
   617 Ni2P2O7
   618 Ni3(PO4)2
   619 NiCl2
   620 NiCl2:2H2O
   621 NiCl2:4H2O
   622 NiCr2O4
   623 NiF2
   624 NiF2:4H2O
   625 NiSO4
   626 NiSO4:6H2O(alpha)
   627 NiTiO3
   628 NiWO4
   629 Nickelbischofite
   630 NiMoO4
   631 Niter
   632 Nitrobarite
   633 Nontronite-Ca
   634 Nontronite-H
   635 Nontronite-K
   636 Nontronite-Mg
   637 Nontronite-Na
   638 Neptunium
   639 Np(OH)4(am)
   640 Np2C3
   641 Np2O5
   642 NpBr3
   643 NpBr4
   644 NpC0.91
   645 NpCl3
   646 NpCl4
   647 NpF3
   648 NpF4
   649 NpF5
   650 NpF6
   651 NpI3
   652 NpN
   653 NpO2
   654 NpO2(am,hyd)
   655 NpO2(NO3)2:6H2O
   656 NpO2CO3
   657 NpO2OH(am)
   658 NpO2OH(am,aged)
   659 NpO3:H2O
   660 NpOBr2
   661 NpOCl2
   662 Okenite
   663 Orpiment
   664 Otavite
   665 Ottemannite
   666 Oxychloride-Mg
   667 Phosphorus
   668 Paragonite
   669 Paralaurionite
   670 Pargasite
   671 Lead
   672 Pb(H2PO4)2
   673 Pb(IO3)2
   674 Pb(N3)2(mono)
   675 Pb(N3)2(orth)
   676 Pb2Cl2CO3
   677 Pb2Cl5NH4
   678 Pb2O(N3)2
   679 Pb2SiO4
   680 Pb3(PO4)2
   681 Pb3SO6
   682 Pb4Cl2(OH)6
   683 Pb4O(PO4)2
   684 Pb4SO7
   685 PbBr2
   686 PbBrF
   687 PbCO3.PbO
   688 PbF2
   689 PbFCl
   690 PbHPO4
   691 PbI2
   692 PbSO4(NH3)2
   693 PbSO4(NH3)4
   694 PbSeO4
   695 Palladium
   696 Pd(OH)2
   697 Pd4S
   698 PdO
   699 PdS
   700 PdS2
   701 Penroseite
   702 Pentahydrite
   703 Periclase
   704 Perovskite
   705 Petalite
   706 Phillipsite
   707 Phlogopite
   708 Phosgenite
   709 Picromerite
   710 Pirssonite
   711 Plattnerite
   712 Plombierite
   713 Promethium
   714 Pm(OH)3
   715 Pm(OH)3(am)
   716 Pm2(CO3)3
   717 Pm2O3
   718 PmF3:.5H2O
   719 PmPO4:1.5H2O
   720 Polydymite
   721 Polyhalite
   722 Portlandite
   723 Powellite
   724 Praseodymium
   725 Pr(OH)3
   726 Pr(OH)3(am)
   727 Pr2(CO3)3
   728 Pr2O3
   729 PrF3:.5H2O
   730 PrPO4:H2O
   731 Prehnite
   732 Pseudowollastonite
   733 Platinum
   734 PtS
   735 PtS2
   736 Plutonium
   737 Pu(HPO4)2(am,hyd)
   738 Pu(OH)3
   739 Pu2C3
   740 Pu2O3
   741 Pu3C2
   742 PuAs
   743 PuBi
   744 PuBi2
   745 PuBr3
   746 PuC0.84
   747 PuCl3
   748 PuCl3:6H2O
   749 PuCl4
   750 PuF3
   751 PuF6
   752 PuI3
   753 PuN
   754 PuO1.61
   755 PuO2
   756 PuO2(hyd,aged)
   757 PuO2(NO3)2:6H2O
   758 PuO2(OH)2:H2O
   759 PuO2CO3
   760 PuO2OH(am)
   761 PuOBr
   762 PuOCl
   763 PuOF
   764 PuOI
   765 PuP
   766 PuPO4(s,hyd)
   767 PuSb
   768 Pyrite
   769 Pyrolusite
   770 Pyromorphite
   771 Pyromorphite-OH
   772 Pyrophyllite
   773 Pyrrhotite
   774 Quartz
   775 Radium
   776 Ra(NO3)2
   777 RaCl2:2H2O
   778 RaSO4
   779 Rankinite
   780 Rubidium
   781 Rb2UO4
   782 Rhenium
   783 Realgar
   784 Rhodium
   785 Rh2O3
   786 Rhodochrosite
   787 Rhodonite
   788 Ripidolite-14A
   789 Ripidolite-7A
   790 Riversideite
   791 Romarchite
   792 Ruthenium
   793 Ru(OH)3:H2O(am)
   794 RuBr3
   795 RuCl3
   796 RuI3
   797 RuO2
   798 RuO2:2H2O(am)
   799 RuO4
   800 RuSe2
   801 Rutherfordine
   802 Rutile
   803 Sulfur
   804 Sanbornite
   805 Sanidine_high
   806 Saponite-Ca
   807 Saponite-H
   808 Saponite-K
   809 Saponite-Mg
   810 Saponite-Na
   811 Antimony
   812 Sb(OH)3
   813 Sb2O4
   814 Sb2O5
   815 Sb4O6(cubic)
   816 Sb4O6(orthorhombic)
   817 SbBr3
   818 SbCl3
   819 Scandium
   820 Scacchite
   821 Schoepite
   822 Scolecite
   823 Selenium
   824 Se2O5
   825 SeCl4
   826 SeO3
   827 Sellaite
   828 Sepiolite
   829 Sepiolite(am)
   830 Palygorskite
   831 Shcherbinaite
   832 Silicon
   833 SiO2(am)
   834 Siderite
   835 Sillimanite
   836 Sklodowskite
   837 Samarium
   838 Sm(OH)3
   839 Sm(OH)3(am)
   840 Sm2(CO3)3
   841 Sm2(SO4)3
   842 Sm2O3
   843 SmF3:.5H2O
   844 SmPO4:H2O
   845 Smectite-high-Fe-Mg
   846 Smectite-low-Fe-Mg
   847 Smectite_Reykjanes
   848 Smithsonite
   849 Studtite
   850 Tin
   851 Sn(OH)2
   852 Sn(SO4)2
   853 Sn3S4
   854 SnBr2
   855 SnBr4
   856 SnCl2
   857 SnSO4
   858 SnSe
   859 SnSe2
   860 Soddyite
   861 Sphaerocobaltite
   862 Sphalerite
   863 Spinel
   864 Spinel-Co
   865 Spodumene
   866 Strontium
   867 Sr(NO3)2
   868 Sr(NO3)2:4H2O
   869 Sr(OH)2
   870 Sr2SiO4
   871 Sr3(AsO4)2
   872 SrBr2
   873 SrBr2:6H2O
   874 SrBr2:H2O
   875 SrCl2
   876 SrCl2:2H2O
   877 SrCl2:6H2O
   878 SrCl2:H2O
   879 SrCrO4
   880 SrF2
   881 SrHPO4
   882 SrI2
   883 SrO
   884 SrS
   885 SrSeO4
   886 SrSiO3
   887 SrUO4(alpha)
   888 Starkeyite
   889 Stellerite
   890 Stilbite
   891 Stilleite
   892 Strengite
   893 Strontianite
   894 Sylvite
   895 Syngenite
   896 Tachyhydrite
   897 Talc
   898 Tarapacaite
   899 Terbium
   900 Tb(OH)3
   901 Tb(OH)3(am)
   902 Tb2(CO3)3
   903 Tb2O3
   904 TbF3:.5H2O
   905 TbPO4:2H2O
   906 Technetium
   907 Tc2O7
   908 Tc2O7:H2O
   909 TcO2
   910 TcO2:1.6H2O
   911 Tenorite
   912 Tephroite
   913 Thorium
   914 Th(NO3)4:5H2O
   915 Th(SO4)2
   916 Th.75PO4
   917 Th2S3
   918 Th2Se3
   919 Th7S12
   920 ThBr4
   921 ThCl4
   922 ThF4
   923 ThF4:2.5H2O
   924 ThI4
   925 ThO2:2H2O(am)
   926 ThO2(am)
   927 ThS
   928 ThS2
   929 Thenardite
   930 Thermonatrite
   931 Thorianite
   932 Titanium
   933 Ti2O3
   934 Ti3O5
   935 TiB2
   936 TiBr3
   937 TiBr4
   938 TiC
   939 TiCl2
   940 TiCl3
   941 TiF3
   942 TiF4(am)
   943 TiI4
   944 TiN
   945 TiO
   946 TiO(alpha)
   947 Tiemannite
   948 Titanite
   949 Thallium
   950 TlTcO4
   951 Thulium
   952 Tm(OH)3
   953 Tm(OH)3(am)
   954 Tm2(CO3)3
   955 Tm2O3
   956 TmF3:.5H2O
   957 TmPO4:2H2O
   958 Tobermorite
   959 Tremolite
   960 Trevorite
   961 Tridymite
   962 Troilite
   963 Trona-K
   964 Uranium
   965 U(HPO4)2:4H2O
   966 U(OH)2SO4
   967 U(SO3)2
   968 U(SO4)2
   969 U(SO4)2:4H2O
   970 U(SO4)2:8H2O
   971 U2C3
   972 U2F9
   973 U2O2Cl5
   974 U2O3F6
   975 U2S3
   976 U2Se3
   977 U3As4
   978 U3O5F8
   979 U3P4
   980 U3S5
   981 U3Sb4
   982 U3Se4
   983 U3Se5
   984 U4F17
   985 U5O12Cl
   986 UAs
   987 UAs2
   988 UBr2Cl
   989 UBr2Cl2
   990 UBr3
   991 UBr3Cl
   992 UBr4
   993 UBr5
   994 UBrCl2
   995 UBrCl3
   996 UC
   997 UC1.94(alpha)
   998 UCl2F2
   999 UCl2I2
  1000 UCl3
  1001 UCl3F
  1002 UCl3I
  1003 UCl4
  1004 UCl5
  1005 UCl6
  1006 UClF3
  1007 UClI3
  1008 UF3
  1009 UF4
  1010 UF4:2.5H2O
  1011 UF5(alpha)
  1012 UF5(beta)
  1013 UF6
  1014 UH3(beta)
  1015 UI3
  1016 UI4
  1017 UN
  1018 UN1.59(alpha)
  1019 UN1.73(alpha)
  1020 UO2(AsO3)2
  1021 UO2(IO3)2
  1022 UO2(NO3)2
  1023 UO2(NO3)2:2H2O
  1024 UO2(NO3)2:3H2O
  1025 UO2(NO3)2:6H2O
  1026 UO2(NO3)2:H2O
  1027 UO2(OH)2(beta)
  1028 UO2.25
  1029 UO2.25(beta)
  1030 UO2.3333(beta)
  1031 UO2.6667
  1032 UO2Br2
  1033 UO2Br2:3H2O
  1034 UO2Br2:H2O
  1035 UO2BrOH:2H2O
  1036 UO2Cl
  1037 UO2Cl2
  1038 UO2Cl2:3H2O
  1039 UO2Cl2:H2O
  1040 UO2ClOH:2H2O
  1041 UO2F2
  1042 UO2F2:3H2O
  1043 UO2FOH:2H2O
  1044 UO2FOH:H2O
  1045 UO2HPO4:4H2O
  1046 UO2SO3
  1047 UO2SO4
  1048 UO2SO4:2.5H2O
  1049 UO2SO4:3.5H2O
  1050 UO2SO4:3H2O
  1051 UO3(alpha)
  1052 UO3(beta)
  1053 UO3(gamma)
  1054 Schoepite(dehyd,0.9)
  1055 Schoepite(dehyd,1.0)
  1056 UOBr2
  1057 UOBr3
  1058 UOCl
  1059 UOCl2
  1060 UOCl3
  1061 UOF2
  1062 UOF2:H2O
  1063 UOF4
  1064 UOFOH
  1065 UOFOH:.5H2O
  1066 UP
  1067 UP2
  1068 UP2O7
  1069 UPO5
  1070 US
  1071 US1.9
  1072 US2
  1073 US3
  1074 USb
  1075 USb2
  1076 USe
  1077 USe2(alpha)
  1078 USe2(beta)
  1079 USe3
  1080 Umangite
  1081 Uraninite
  1082 Uranophane(alpha)
  1083 Vanadium
  1084 V2O4
  1085 V3O5
  1086 V4O7
  1087 Vaesite
  1088 Tungsten
  1089 WCl2(s)
  1090 WCl4(s)
  1091 WCl5(s)
  1092 WCl6(s)
  1093 WO2Cl2(s)
  1094 WOCl4(s)
  1095 WOF4(s)
  1096 Wairakite
  1097 Weeksite-K
  1098 Weeksite-Na
  1099 Whitlockite
  1100 Wilkmanite
  1101 Witherite
  1102 Wollastonite
  1103 Wurtzite
  1104 Wustite
  1105 Xonotlite
  1106 Yttrium
  1107 Ytterbium
  1108 Yb(OH)3
  1109 Yb(OH)3(am)
  1110 Yb2(CO3)3
  1111 Yb2O3
  1112 YbF3:.5H2O
  1113 YbPO4:2H2O
  1114 Zincite
  1115 Zircon
  1116 Zinc
  1117 Zn(BO2)2
  1118 Zn(ClO4)2:6H2O
  1119 Zn(IO3)2
  1120 Zn(NO3)2:6H2O
  1121 Zn(OH)2(beta)
  1122 Zn(OH)2(epsilon)
  1123 Zn(OH)2(gamma)
  1124 Zn2(OH)3Cl
  1125 Zn2SO4(OH)2
  1126 Zn2SiO4
  1127 Zn2TiO4
  1128 Zn3(AsO4)2
  1129 Zn3O(SO4)2
  1130 Zn5(NO3)2(OH)8
  1131 ZnBr2
  1132 ZnBr2:2H2O
  1133 ZnCO3:H2O
  1134 ZnCl2
  1135 ZnCl2(NH3)2
  1136 ZnCl2(NH3)4
  1137 ZnCl2(NH3)6
  1138 ZnCr2O4
  1139 ZnF2
  1140 ZnI2
  1141 ZnSO4
  1142 ZnSO4:6H2O
  1143 ZnSO4:7H2O
  1144 ZnSO4:H2O
  1145 ZnSeO3:H2O
  1146 Zoisite
  1147 Zirconium
  1148 ZrB2
  1149 ZrC
  1150 ZrCl
  1151 ZrCl2
  1152 ZrCl3
  1153 ZrCl4
  1154 ZrF4(beta)
  1155 ZrH2
  1156 ZrN


 liquids

     1 Bromine
     2 Quicksilver

 * Note - (EQPT/pcrsg) The pure liquids block has
       not been written on the DATA1 and DATA1F files,
       because the EQ3NR and EQ6 codes presently do not
       treat non-aqeuous liquids.


 gases

     1 Ag(g)
     2 Al(g)
     3 Am(g)
     4 AmF3(g)
     5 Argon
     6 B(g)
     7 BF3(g)
     8 Be(g)
     9 Br2(g)
    10 C(g)
    11 CH4(g)
    12 CO(g)
    13 CO2(g)
    14 Ca(g)
    15 Cd(g)
    16 Chlorine
    17 CoCl2(g)
    18 CoCl3(g)
    19 CoF2(g)
    20 CrCl4(g)
    21 Cs(g)
    22 Cu(g)
    23 Fluorine
    24 FeCl2(g)
    25 FeCl3(g)
    26 FeF2(g)
    27 FeF3(g)
    28 H2(g)
    29 H2O(g)
    30 H2O2(g)
    31 H2S(g)
    32 HBr(g)
    33 HCl(g)
    34 HF(g)
    35 HI(g)
    36 HNO3(g)
    37 Helium
    38 Hf(g)
    39 Hg(g)
    40 I2(g)
    41 K(g)
    42 Krypton
    43 Li(g)
    44 Mg(g)
    45 Nitrogen
    46 NH3(g)
    47 NO(g)
    48 NO2(g)
    49 NO3(g)
    50 N2O(g)
    51 N2O3(g)
    52 N2O4(g)
    53 N2O5(g)
    54 Na(g)
    55 Neon
    56 NiCl2(g)
    57 NiF2(g)
    58 O2(g)
    59 Pb(g)
    60 Rb(g)
    61 Radon
    62 RuCl3(g)
    63 RuO3(g)
    64 S2(g)
    65 SO2(g)
    66 Si(g)
    67 SiF4(g)
    68 Sn(g)
    69 Tc(g)
    70 Tc2O7(g)
    71 TcC(g)
    72 TcO(g)
    73 TcS(g)
    74 Th(g)
    75 Ti(g)
    76 TiBr4(g)
    77 TiCl(g)
    78 TiCl2(g)
    79 TiCl3(g)
    80 TiCl4(g)
    81 TiF(g)
    82 TiF2(g)
    83 TiF3(g)
    84 TiF4(g)
    85 TiO(g)
    86 U(g)
    87 U2Cl10(g)
    88 U2Cl8(g)
    89 U2F10(g)
    90 UBr(g)
    91 UBr2(g)
    92 UBr3(g)
    93 UBr4(g)
    94 UBr5(g)
    95 UCl(g)
    96 UCl2(g)
    97 UCl3(g)
    98 UCl4(g)
    99 UCl5(g)
   100 UCl6(g)
   101 UF(g)
   102 UF2(g)
   103 UF3(g)
   104 UF4(g)
   105 UF5(g)
   106 UF6(g)
   107 UI(g)
   108 UI2(g)
   109 UI3(g)
   110 UI4(g)
   111 UO(g)
   112 UO2(g)
   113 UO2Cl2(g)
   114 UO2F2(g)
   115 UO3(g)
   116 UOF4(g)
   117 WCl2(g)
   118 WCl4(g)
   119 WCl6(g)
   120 WF(g)
   121 WF6(g)
   122 WO2Cl2(g)
   123 WOCl4(g)
   124 WOF4(g)
   125 Xenon
   126 Zn(g)
   127 Zr(g)
   128 ZrF4(g)


 solid solutions

 Biotite
     Number of components=   2  Activity coefficient code=   1
            Annite                    Phlogopite
 Carbonate-Calcite
     Number of components=   6  Activity coefficient code=   1
            Calcite                   Magnesite
            Rhodochrosite             Siderite
            Smithsonite               Strontianite
 Chlorite-ss
     Number of components=   2  Activity coefficient code=   1
            Clinochlore-14A           Daphnite-14A
 Clinoptilolite-ss
     Number of components=   6  Activity coefficient code=   1
            Clinoptilolite-Ca         Clinoptilolite-Cs
            Clinoptilolite-K          Clinoptilolite-NH4
            Clinoptilolite-Na         Clinoptilolite-Sr
 Epidote-ss
     Number of components=   2  Activity coefficient code=   1
            Clinozoisite              Epidote
 Garnet-ss
     Number of components=   2  Activity coefficient code=   1
            Andradite                 Grossular
 Orthopyroxene
     Number of components=   2  Activity coefficient code=   1
            Enstatite                 Ferrosilite
 Plagioclase
     Number of components=   2  Activity coefficient code=   1
            Albite_high               Anorthite
 Sanidine-ss
     Number of components=   2  Activity coefficient code=   1
            Albite_high               Sanidine_high
 Beidellite-ss
     Number of components=   5  Activity coefficient code=   1
            Beidellite-Ca             Beidellite-H
            Beidellite-K              Beidellite-Mg
            Beidellite-Na
 Montmorillonite-ss
     Number of components=   5  Activity coefficient code=   1
            Montmorillonite-Ca        Montmorillonite-H
            Montmorillonite-K         Montmorillonite-Mg
            Montmorillonite-Na
 Nontronite-ss
     Number of components=   5  Activity coefficient code=   1
            Nontronite-Ca             Nontronite-H
            Nontronite-K              Nontronite-Mg
            Nontronite-Na
 Saponite-ss
     Number of components=   5  Activity coefficient code=   1
            Saponite-Ca               Saponite-H
            Saponite-K                Saponite-Mg
            Saponite-Na


  Aqueous species hard core diameter coverage:

     1218 aqueous species have hard core diameters specified on the data file
     1218 aqueous solute species are present on this file
    Coverage is 100.00 per cent


 Completed processing the SEDH data file data0.geo.R2.


 No errors were encountered.

 No warnings were encountered.


          Start time = 10:35:00  04Dec2021
            End time = 10:35:00  04Dec2021

           Run time =  0.141     seconds

 Normal exit