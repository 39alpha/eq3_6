 EQPT Species List (SLIST) File:


           Number of elements =    40
           Number of basis species =    69

 data0.ypf.R2
 PITZER THERMODYNAMIC DATABASE INCLUDING ACTINIDES AND TRANSITION METALS
 (09/28/2006)
 
 This database file represents a combination of the Pitzer parameters for
 major cations and anions found in data0.ypf (renamed data0.ypf.R0)(DTN:
 SN0302T0510102.002) and added Pitzer parameters for actinide and transition
 metal species. However, the listed solid and aqueous species in this
 data0.ypf.R2 database are somewhat different from those in data0.ypf.R0.
 Some log K values for solids and aqueous species were actually taken from
 the data0.ymp.R4 database source (DTN: SN0410T0510404.002 ) or its updated
 version data0.ymp.R5 (DTN: MO0608SPAYMPR5.000), but log K data for some solids
 and aqueous species were obtained from other sources in so far as to maintain
 consistency with the newly-acquired Pitzer parameters. Also, additional log K
 data for some solids was estimated based on the agreement between model
 predictions and reported solubility data. Documentation of this database can
 be found in ANL-WIS-GS-000001 and ANL-EBS-MD-000045. This version supersedes
 the preliminary version data0.ypf.R2 (DTN: MO0608SPAYPFR2.000). Since this
 database incorporates data from data0.ypf.R0(DTN: SN0302T0510102.002),
 calculations using data valid to elevated temperatures are possible. However,
 the user must be aware of the database limitations regarding the valid range
 of temperatures and pressures as well as ionic strength for the Pitzer
 parameters assigned to the species in this database. Therefore, the user
 must be familiar with these limitations and should consult the reports
 associated with the database development and/or relevant sources for these
 parameters. Use of these parameters outside their respective ranges of
 validity is inadvisable and is not permitted for applications on the Yucca
 Mountain Project unless specific justification is provided.
 
 This datafile uses the standard form of Pitzer's equations relevant to
 the temperature range. Details on the compilation of these Pitzer
 parameters are given in the document 'In-Drift Precipitates/Salts Model'
 (ANL-EBS-MD-000045) and 'Pitzer Thermodynamic Database for some Actinides and
 Transition Metal Species (data0.ypf.R1)' (ANL-WIS-GS-000001).
 
 Pitzer parameters are represented by the 25C-centric four-term temperature
 function given by:
 
 x(T) = a1 + a2*(1/T - 1/298.15) + a3*ln(T/298.15) + a4*(T - 298.15)
 
 where T is temperature in Kelvin and a1 through a4 denote the temperature
 function fitting coefficients for the temperature-dependent Pitzer
 parameters. The conversion of non-standard or expanded forms of Pitzer
 interaction parameters recently adopted by several workers for highly
 soluble salts to the standard form currently embedded in EQ3/6 Version 8.0
 was conducted using the approach described in 02Rard/Wij.  This conversion
 imposes usage limits on these parameters within a valid range of temperature
 and ionic strength.
 
 Some data blocks provide comments on the gathered Pitzer parameters and
 solid phase solubility data to make the user aware of any convention
 adopted in the data extraction/compilation process. DTN numbers for some
 sources are given at the bottom of this file within the references section.
 Also, spreadsheet and/or other types of source files are includes within
 brackets '[]'.
 
 Below is a list of current corrections to this database in addition to those
 made for the version data0.ypf.R1:
 
  (1) Adding a zero to the beta(2) a1 parameter for K+ - F- ca pair.
  (2) Adding a space after the '=' symbol and the corresponding value
        of the a1 constant in some of the ca parameters for Mn-Br, Ba-Cl,
        Ba-Br, Ba-I, Ba-ClO4, Ba-NO3, Ca-ClO4, Co-Br, Co-I, Co-NO3, Co-SO4,
        Cr-Cl, Cr-NO3, Na-Cr(OH)4, Na-Cr2O2(OH)4, K-Cr2O7, Cu-Cl, Cu-SO4,
        Cu-Br, Fe-Cl, K-CrO4, K-PO4, Mn-Cl, Mn-SO4, Na-ClO4, Na-PO4, Na-CrO4,
        Ni-Cl, Ni-Br, Ni-SO4, Zn-Cl, Zn-Br, Zn-I, Zn-NO3, Gd-NO3 and Zn-SO4.
  (3) Similar fixes to those described in (2) for some cation-cation and
        anion-anion interaction (theta) parameters: Al-Na, Cr-Na, Cr-K, Cu-K,
        and Br-Cr2O7.
  (4) Correction to the Beta(0) and Beta(1) parameters Cs-ClO4
  (5) Correction to the C(phi) parameter for Cr-NO3
  (6) Correction to the C(phi) parameter for Na-CrO4
  (7) Correction to the C(phi) parameter for Fe-SO4
  (8) Correction to the Beta(0)parameter for Zn-ClO4
  (9) Correction to the theta parameter for Na-Co
 (10) Correction to the psi parameter for Na+/B4O5(OH)4--/Cl-
 (11) Correction to the C(phi) parameter for Na-PO4
 (12) Updates/corrections to several references
 (13) Correction to the log K for Th6(OH)15(9+) aqueous complexation reaction
 (14) Correction to the log K for Th4(OH)12++++ aqueous complexation reaction
 (15) Correction to the sign of the log K value for NpO2(CO3)3(5-),
        NpO2(CO3)2---, and NpO2CO3-
 (16) Corrections to the log K values for Fe(OH)2(aq), NiOH+, NiCO3(aq),
        NiHCO3+, and NiPO4-
 (17) Correction to the C(phi) parameter for Co-NO3
 (18) Log K values for Np, U, and Pu solids from data0.ymp.R4 (Np and Pu)
        and DTN SN0410T0510404.001 (U only): Neptunium, Np2O5, NpOCl2, NpCl4,
        NpF4, NpF5, Na3NpF8, Plutonium, Pu2O3, Pu2C3, Pu3C2, PuC0.84, PuF3,
        PuF6, PuO1.61, PuOCl, PuOF, uranium, U(SO4)2, UCl4, UCl6, UF4, UF6,
        UI4, UO2Cl2, UO2F2, UO2SO4, UO3(gamma), and UOCl2.
 (19) Log K values for some zeolite phases: Analcime, Chabazite,
        Na-Clinoptilolite, Erionite, Laumontite, Mesolite, Phillipsite, and
        Stellerite from data0.ymp.R4
 (20) Log K values obtained from DTN SN0410T0510404.001 for the solid phases
        NiCr2O4, CrF3, CrO3, Iron, Fe(OH)2, FeF2, Ferrite-Ca, Ferrite_Cu,
        Ferrite-Dicalcium, Ferrite-Mg, Ferrite-Ni, Co(OH)2, CoCl2, CoCr2O4,
        CoF2, CoFe2O4, CoO, MnSO4, MoO2Cl2, and Witherite
 (21) Correction of the Pu phase PuO2(OH)2:H2O stoichiometry using the
        formula of Lemire (2001 [DIRS 159027]). Corrected data block obtained
        from data0.ymp.R4
 (22) Correction to the log K value at 25C for the aqueous species for
        Cr2O2(OH)4-- to the value of 51.9302 from the previous value of 4.05.
 (23) Removal of ternary parameters (theta) for Cr3+ - Na+ and Al3+ - Na+
        (source 02Chr(2))
 (24) Removal of the phases Cr2(SO4)3.Na2SO4:24H2O and Cr2(SO4)3.K2SO4:24H2O
        (source 02Chr(1))
 (25) Addition of the psi parameter for Ni++ - Na+ - SO4-- from 02Chr(2),
        Table 4.
 (26) Replace Pitzer parameters for 2-2 electrolytes (NiSO4, CuSO4, ZnSO4,
        and MnSO4) by 03Gue/Mou by those in Table I of 74Pit/May.
 (27) Add the recommended beta(2) value of -40.0 by 74Pit/May to the above
        datablocks for which beta(2) is not given in the source Table I
        (see Section 6.2.6.1 of ANL-WIS-GS-000001 model report. Also, the
        source pointer was updated to 74Pit/May instead of 91Pit.
 (28) The binary parameters from 98Alb/Riz for Cu-NO3 were replaced by those
        listed in 88Kim/Fre (Table III)
 (29) Addition of the gas species CH4(g).
 (30) Addition of binary Ni-NO3 Pitzer parameters.
 (31) Removal of the aqueous species Ni(NO3)2(aq), NaF(aq), UO2Cl2(aq),
        and FeCl2(aq).
 (32) Change of the UO3(aq) log K value to 12.15 using data given by 03Gui.
 (33) Correction: addition of a negative sign  to the C(phi) parameter for
        Cu-SO4.
 (34) Replaced all the Pitzer lambda parameters for CO2(aq) with those based
      on the work by 90Cor/dePab, 93Rum/Mau, and 94rum/nic. These new lambda
      parameters are valid to temperatures up to 250 deg.C.
 (35) The lambda (CO2(aq)-K) and lambda (CO2(aq)-HSO4-)Pitzer parameters
      sourced to data0.ypf.R0 were deleted from the database.
 (36) Modified log K value for Calcite solubility.
 (37) Replaced zetas for CO2(aq)-Na-Cl and CO2(aq)-Na-SO4.
 (38) Deleted remaining zeta parameters for CO2(aq) sourced to data0.ypf.R0.
 (39) Removal of Ca-HCO3 binary Pitzer parameters that were based on the work
        by 93He/Mor.
 (40) Added/updated high order parameters from 82pei/pit and 99kon/kon for
        Cl-CO3, Cl-HCO3, Na-Cl-HCO3, and Na-CO3-Cl.
 (41) Update/addition of CaCl2 Pitzer parameter and ion-pair log K values
        after Model 3 of 98ste/fel. Adoption of this model also includes the
        addition of the ion pair CaCl+ and CaCl2(aq).
 (42) Update of log K values for Soda Niter, Thermonatrite, Antarcticite,
        CaCl2:4H2O, CaCl2:2H2O.
 (43) Editorial changes the involve formatting and fixing two entries in the
      'references' section.
 
 BEGIN CONFIGURATION DATA BLOCK
 Do not change the data in this block unless you know what you
 are doing.
 INTERPRET 500 AS NO DATA= YES
 SPARSE GRID RANGE CONDITION=IGNORE
 Pitzer data parameters:
    PITZER DATA BLOCK ORG.= NEW
    PITZER TEMP FUNCTION= LIVERMORE
    NO. OF PITZER TEMP FUNC TERMS= 4
 END CONFIGURATION DATA
 +--------------------------------------------------------------------


 element = O       , atwt =   15.99940
 element = Al      , atwt =   26.98154
 element = Ba      , atwt =  137.32700
 element = Br      , atwt =   79.90400
 element = C       , atwt =   12.01100
 element = Ca      , atwt =   40.07800
 element = Cl      , atwt =   35.45270
 element = Cr      , atwt =   51.99610
 element = Cs      , atwt =  132.90543
 element = F       , atwt =   18.99840
 element = Fe      , atwt =   55.84700
 element = H       , atwt =    1.00794
 element = I       , atwt =  126.90447
 element = K       , atwt =   39.09830
 element = Li      , atwt =    6.94100
 element = Mg      , atwt =   24.30500
 element = Mn      , atwt =   54.93805
 element = Mo      , atwt =   95.94000
 element = N       , atwt =   14.00674
 element = Na      , atwt =   22.98977
 element = P       , atwt =   30.97376
 element = Rb      , atwt =   85.46780
 element = S       , atwt =   32.06600
 element = Si      , atwt =   28.08550
 element = Sr      , atwt =   87.62000
 element = B       , atwt =   10.81100
 element = Am      , atwt =  243.00000
 element = Cm      , atwt =  247.00000
 element = Gd      , atwt =  157.25000
 element = Nd      , atwt =  144.24000
 element = Np      , atwt =  237.00000
 element = Pu      , atwt =  244.00000
 element = Tc      , atwt =   99.00000
 element = Th      , atwt =  232.03810
 element = U       , atwt =  238.02890
 element = Co      , atwt =   58.93320
 element = Cu      , atwt =   63.54600
 element = Ni      , atwt =   58.69340
 element = Zn      , atwt =   65.39000
 element = Zr      , atwt =   91.22400


 aqueous

     1  H2O                       Al+++
     3  Ba++                      Br-
     5  HCO3-                     Ca++
     7  Cl-                       Cr+++
     9  Cs+                       F-
    11  Fe++                      H+
    13  I-                        K+
    15  Li+                       Mg++
    17  Mn++                      MoO4--
    19  NO3-                      Na+
    21  HPO4--                    Rb+
    23  SO4--                     SiO2(aq)
    25  Sr++                      B(OH)3(aq)
    27  Am+++                     Cm+++
    29  Gd+++                     Nd+++
    31  NpO2+                     PuO2+
    33  TcO4-                     Th++++
    35  UO2++                     Co++
    37  Cu++                      Ni++
    39  Zn++                      Zr++++
    41  O2(g)                     ClO4-
    43  CrO4--                    Fe+++
    45  H2(aq)                    HPO3--
    47  HS-                       Mn+++
    49  MnO4--                    NH4+
    51  NO2-                      SO3--
    53  Am++                      Am++++
    55  AmO2+                     AmO2++
    57  BH4-                      Np+++
    59  Np++++                    NpO2++
    61  PO4---                    Pu+++
    63  Pu++++                    PuO2++
    65  TcO++                     U++++
    67  O2(aq)                    OH-
    69  CO3--                     AlO2-
    71  AlOH++                    AlO+
    73  HAlO2(aq)                 AlF++
    75  AlF2+                     AlF3(aq)
    77  AlF4-                     B(OH)4-
    79  B3O3(OH)4-                B4O5(OH)4--
    81  CaB(OH)4+                 CaCO3(aq)
    83  CaOH+                     CaSO4(aq)
    85  Cm(OH)3(aq)               CaCl+
    87  CaCl2(aq)                 CO2(aq)
    89  Cr2O7--                   FeCl+
    91  FeHCO3+                   FeCO3(aq)
    93  Fe(CO3)2--                FeOH+
    95  Fe(OH)2(aq)               FeOH++
    97  Fe(OH)2+                  Fe(OH)3(aq)
    99  Fe(OH)4-                  FeSO4+
   101  Fe(SO4)2-                 FeCl++
   103  FeCl2+                    FeF+
   105  FeF++                     FeF2+
   107  HSO4-                     HSiO3-
   109  H2PO4-                    H3PO4(aq)
   111  MgCO3(aq)                 MgHCO3+
   113  MgOH+                     NH3(aq)
   115  MgB(OH)4+                 Am(CO3)+
   117  Am(CO3)2-                 Am(CO3)3---
   119  Am(CO3)4(5-)              Am(OH)2+
   121  Am(OH)++                  Am(SO4)2-
   123  Am(SO4)+                  AmCl++
   125  AmCl2+                    Am(OH)3(aq)
   127  AmF++                     AmF2+
   129  AmH2PO4++                 AmNO3++
   131  Cm(CO3)+                  Cm(CO3)2-
   133  Cm(CO3)3---               Cm(CO3)4(5-)
   135  Cm(OH)2+                  Cm(OH)++
   137  Cm(SO4)2-                 Cm(SO4)+
   139  CmCl++                    CmCl2+
   141  NpOH+++                   Np(OH)2++
   143  Nd(CO3)+                  Nd(CO3)2-
   145  Nd(CO3)3---               Nd(CO3)4(5-)
   147  Nd(OH)3(aq)               Nd(OH)2+
   149  Nd(OH)++                  Nd(SO4)2-
   151  Nd(SO4)+                  NdCl++
   153  NdCl2+                    Np(OH)3+
   155  Np(OH)4(aq)               Pu(CO3)+
   157  Pu(CO3)2-                 Pu(CO3)3---
   159  Pu(CO3)4(5-)              Pu(OH)3(aq)
   161  Pu(OH)2+                  Pu(OH)++
   163  Pu(SO4)2-                 Pu(SO4)+
   165  PuCl++                    PuCl2+
   167  PuOH+++                   Pu(OH)2++
   169  Pu(OH)3+                  Pu(OH)4(aq)
   171  ThOH+++                   Th(OH)2++
   173  Th(OH)3+                  Th(OH)4(aq)
   175  Th2(OH)2(6+)              Th4(OH)8(8+)
   177  Th4(OH)12++++             Th6(OH)15(9+)
   179  Th(CO3)5(6-)              Th(SO4)2(aq)
   181  Th(SO4)3--                Th(SO4)4----
   183  UOH+++                    U(OH)2++
   185  U(OH)3+                   U(OH)4(aq)
   187  NiCrO4(aq)                H2CrO4(aq)
   189  HCrO4-                    CrBr++
   191  CrOH++                    Cr(OH)2+
   193  Cr(OH)3(aq)               Cr(OH)4-
   195  Cr2O2(OH)4--              NiOH+
   197  Ni(OH)2(aq)               Ni(OH)3-
   199  Ni(OH)4--                 NiBr+
   201  NiCO3(aq)                 NiF+
   203  NiHCO3+                   NiHPO4(aq)
   205  NiH2PO4+                  NiPO4-
   207  NiHP2O7-                  NiP2O7--
   209  NiNO3+                    NpO2OH(aq)
   211  NpO2(OH)2-                NpO2(CO3)3(5-)
   213  NpO2(CO3)2---             NpO2CO3-
   215  PuO2(CO3)3----            PuO2(CO3)2--
   217  PuO2CO3(aq)               UO2CO3(aq)
   219  UO2(CO3)2--               UO2(CO3)3----
   221  (UO2)2CO3(OH)3-           (UO2)11(CO3)6(OH)12--
   223  (UO2)3(CO3)6(6-)          (UO2)3O(OH)2(HCO3)+
   225  (UO2)2NpO2(CO3)6(6-)      (UO2)2(PuO2)(CO3)6(6-)
   227  UO2OH+                    (UO2)2(OH)2++
   229  (UO2)3(OH)4++             (UO2)2OH+++
   231  (UO2)3(OH)5+              (UO2)3(OH)7-
   233  (UO2)4(OH)7+              UO3(aq)
   235  HUO4-                     UO4--
   237  UO2Cl+


 minerals

     1  Albite                    Alunite
     3  Amesite-7A                Amesite-14A
     5  Analcime                  Analcime-dehy
     7  Anhydrite                 Antarcticite
     9  Antigorite(am)            Aragonite
    11  Arcanite                  Artinite
    13  Barite                    Beidellite-Mg
    15  Beidellite-Ca             Beidellite-K
    17  Beidellite-Na             Beidellite-H
    19  Bischofite                Bloedite
    21  Boehmite                  Borax
    23  Boric_acid                KB5O8:4H2O
    25  K2B4O7:4H2O               NaBO2:4H2O
    27  NaB5O8:5H2O               NaBO2:NaCl:2H2O
    29  Brushite                  Burkeite
    31  CaBr2                     Ca2Cl2(OH)2:H2O
    33  Ca4Cl2(OH)6:13H2O         CaCl2
    35  CaCl2:2H2O                CaCl2:4H2O
    37  CaI2                      Calcite
    39  Ca(NO3)2                  Ca(NO3)2:2H2O
    41  Ca(NO3)2:3H2O             Ca(NO3)2:4H2O
    43  Carnallite                Carobbite
    45  Celadonite                Celestite
    47  Chabazite                 Chamosite-7A
    49  Chloromagnesite           Clinoptilolite
    51  Clinoptilolite-dehy       Clinoptilolite-Ca
    53  Clinoptilolite-Cs         Clinoptilolite-K
    55  Clinoptilolite-NH4        Clinoptilolite-Na
    57  Clinoptilolite-Sr         Cristobalite(alpha)
    59  Cronstedtite-7A           Cryolite
    61  Daphnite-14A              Daphnite-7A
    63  Dawsonite                 Dolomite
    65  Epsomite                  Erionite
    67  Ferroaluminoceladonite    Ferroceladonite
    69  Fe2(MoO4)3                FeF3
    71  Fe2(SO4)3                 Fluorapatite
    73  Fluorite                  Gaylussite
    75  Gibbsite                  Glaserite
    77  Glauberite                Goethite
    79  Greenalite                Gypsum
    81  Halite                    Hematite
    83  Hemihydrate               Heulandite
    85  Hexahydrite               Huntite
    87  Hydroxylapatite           Hydromagnesite
    89  Illite                    Jarosite
    91  Jarosite-Na               K-Feldspar
    93  K2CO3                     K2CO3:1.5H2O
    95  K2O                       K2Si4O9
    97  K3H(SO4)2                 K8H4(CO3)6:3H2O
    99  Kainite                   KAlCl4
   101  K2HPO4                    K3AlCl6
   103  K3AlF6                    K3PO4
   105  Kalicinite                KAl(SO4)2
   107  KAl(SO4)2:3H2O            KAl(SO4)2:12H2O
   109  Kaolinite                 KBr
   111  KClO4                     KH2PO4
   113  KI                        Kieserite
   115  KMgCl3:2H2O               KNaCO3:6H2O
   117  KOH                       Labile_Salt
   119  Lansfordite               Laumontite
   121  Leonhardtite              Leonite
   123  Lime                      Magnesite
   125  Magnetite                 Maximum_Microcline
   127  Mercallite                Mesolite
   129  MgBr2                     MgCl2:H2O
   131  MgCl2:2H2O                MgCl2:4H2O
   133  MgI2                      MgMoO4
   135  Mg(NO3)2                  MgOHCl
   137  MgSO4                     Minnesotaite
   139  Mirabilite                Misenite
   141  Molysite                  Montmorillonite-H
   143  Montmorillonite-Na        Montmorillonite-K
   145  Montmorillonite-Ca        Montmorillonite-Mg
   147  Mordenite                 NaBr
   149  NaClO4                    NaI
   151  NaNO2                     NaOH
   153  Na2CO3:7H2O               Na2CrO4
   155  Na2MoO4                   Na2O
   157  Na2SO4(Sol-3)             Na3H(SO4)2
   159  Kogarkoite                Na4Ca(SO4)3:2H2O
   161  Nahcolite                 Natrite
   163  Natron                    Natrolite
   165  Nesquehonite              NH4Cl
   167  NH4ClO4                   NH4I
   169  (NH4)2SO4                 Niter
   171  Nontronite-Mg             Nontronite-Ca
   173  Nontronite-K              Nontronite-Na
   175  Nontronite-H              Oxychloride-Mg
   177  Pentahydrite              Pentasalt
   179  Periclase                 Phillipsite
   181  Picromerite               Pirssonite
   183  Polyhalite                Portlandite
   185  Pyrolusite                Pyrophyllite
   187  Quartz                    Ripidolite-7A
   189  Ripidolite-14A            Saponite-H
   191  Saponite-Na               Saponite-K
   193  Saponite-Ca               Saponite-Mg
   195  Scolecite                 Sellaite
   197  Sepiolite                 Sepiolite(am)
   199  Palygorskite              SiO2(am)
   201  Smectite-high-Fe-Mg       Smectite-low-Fe-Mg
   203  Soda Niter                Strontianite
   205  SrBr2                     SrCl2
   207  SrF2                      SrI2
   209  SrMoO4                    SrO
   211  Sr(OH)2                   Stellerite
   213  Stilbite                  Sylvite
   215  Syngenite                 Tachyhydrite
   217  Talc                      Thenardite
   219  Thermonatrite             Trona
   221  Trona-K                   Villiaumite
   223  Whitlockite               Americium
   225  Am(OH)3                   Am(OH)3(am)
   227  Am2(CO3)3                 Am2O3
   229  AmCl3                     AmF3
   231  AmF4                      AmO2
   233  AmOCl                     AmOHCO3
   235  AmPO4(am)                 Brucite
   237  Chromium                  Eskolaite
   239  Cr(OH)3(am)               CrCl3:6H2O
   241  CrF3                      KFe3(CrO4)2(OH)6
   243  Na2CrO4:4H2O              Tarapacaite
   245  Lopezite                  CrO3
   247  CuCr2O4                   Copper
   249  CuCl2:2H2O                CuSO4:5H2O
   251  CuSO4.Na2SO4:2H2O         CuSO4.K2SO4:6H2O
   253  Cu3(PO4)2                 Iron
   255  Fe(OH)2                   Fe(OH)3
   257  Siderite                  FeF2
   259  FeO                       FeSO4
   261  Melanterite               Ferrite-Ca
   263  Ferrite-Cu                Ferrite-Dicalcium
   265  Ferrite-Mg                Ferrite-Ni
   267  Ferrite-Zn                Gadolinium
   269  Gd(OH)3                   Gd(OH)3(am)
   271  Gd2(CO3)3                 Gd2O3
   273  GdF3:.5H2O                GdPO4:10H2O
   275  Hausmannite               Hydrozincite
   277  Lawrencite                Cobalt
   279  Co(OH)2                   Co2SiO4
   281  Co3(PO4)2                 CoCl2
   283  CoCl2:2H2O                CoCl2:6H2O
   285  CoCr2O4                   CoF2
   287  CoFe2O4                   CoHPO4
   289  CoO                       CoSO4:3Co(OH)2
   291  CoSO4:6H2O                Bieberite
   293  Malachite                 Mg2(OH)3Cl:4H2O
   295  Manganite                 Manganosite
   297  MgUO4                     Manganese
   299  Mn(OH)2(am)               Rhodochrosite
   301  MnCl2:4H2O                MnSO4
   303  Molybdite                 Powellite
   305  MoO2Cl2                   Sodium
   307  Neodymium                 Nd(OH)3
   309  Nd(OH)3(am)               Nd(OH)3(c)
   311  Nd2(CO3)3                 Nd2O3
   313  NdF3:.5H2O                NdOHCO3
   315  NdPO4:10H2O               Nickel
   317  NiCO3                     NiCl2:4H2O
   319  Nickelbischofite          NiMoO4
   321  Ni(OH)2                   Ni2P2O7
   323  Ni3(PO4)2                 Ni4CrO4(OH)6
   325  NiSO4:7H2O                NiCr2O4
   327  NiF2:4H2O                 Neptunium
   329  NpO2(cr)                  Np2O5
   331  Np(OH)4(am)               NpOCl2
   333  NpCl4                     NpF4
   335  NpF5                      NpO2(NO3)2:6H2O
   337  NpO2CO3                   NpO2OH(am)
   339  NpO2OH(am,aged)           NpO3:H2O
   341  Na3NpF8                   Na3NpO2(CO3)2
   343  NaNpO2CO3                 NaNpO2CO3:3.5H2O
   345  KNpO2CO3                  K3NpO2(CO3)2
   347  Plutonium                 PuO2(cr)
   349  Pu(OH)4(am)               Pu(HPO4)2(am,hyd)
   351  Pu(OH)3                   Pu2O3
   353  Pu2C3                     Pu3C2
   355  PuC0.84                   PuCl3:6H2O
   357  PuCl4                     PuF3
   359  PuF6                      PuO1.61
   361  PuO2(NO3)2:6H2O           PuO2(OH)2:H2O
   363  PuO2CO3                   PuOCl
   365  PuOF                      PuPO4(s,hyd)
   367  Technetium                Tc2O7
   369  Tc2O7:H2O                 TcO2
   371  TcO2:1.6H2O               NaTcO4:4H2O
   373  KTcO4                     CsTcO4
   375  CsClO4                    Thorium
   377  Th(OH)4(am)               ThO2(cr)
   379  ThO2(am)                  Th(SO4)2
   381  Th.75PO4                  ThF4
   383  ThF4:2.5H2O               Uranium
   385  UO2(cr)                   U(OH)4(am)
   387  Schoepite                 Ca(UO2)6O4(OH)6:8H2O
   389  Boltwoodite-Na            UO2CO3
   391  U(HPO4)2:4H2O             U(OH)2SO4
   393  U(SO3)2                   U(SO4)2
   395  U(SO4)2:4H2O              U(SO4)2:8H2O
   397  U2O3F6                    U3O5F8
   399  UCl2F2                    UCl3F
   401  UCl4                      UCl6
   403  UClF3                     UF4
   405  UF4:2.5H2O                UF6
   407  UI4                       UO2(NO3)2
   409  UO2(NO3)2:2H2O            UO2(NO3)2:3H2O
   411  UO2(NO3)2:6H2O            UO2(NO3)2:H2O
   413  UO2(OH)2(beta)            UO2.3333(beta)
   415  UO2.6667                  UO2Cl2
   417  UO2Cl2:3H2O               UO2Cl2:H2O
   419  UO2ClOH:2H2O              UO2F2
   421  UO2F2:3H2O                UO2FOH:2H2O
   423  UO2FOH:H2O                UO2HPO4:4H2O
   425  UO2SO3                    UO2SO4
   427  UO2SO4:2.5H2O             UO2SO4:3.5H2O
   429  UO2SO4:3H2O               UO3(alpha)
   431  UO3(beta)                 UO3(gamma)
   433  UO3:.9H2O(alpha)          UOCl2
   435  UOF2                      UOF2:H2O
   437  UOFOH                     UOFOH:.5H2O
   439  UP2O7                     US2
   441  US3                       Uranophane(alpha)
   443  Na2U2O7                   Na2UO4(alpha)
   445  Na4UO2(CO3)3              Witherite
   447  Zinc                      Zincite
   449  Zn(BO2)2                  Zn(ClO4)2:6H2O
   451  Zn(NO3)2:6H2O             Zn(OH)2(beta)
   453  Zn(OH)2(epsilon)          Zn(OH)2(gamma)
   455  Zn2(OH)3Cl                Zn2SO4(OH)2
   457  Zn3O(SO4)2                Zn5(NO3)2(OH)8
   459  ZnBr2:2H2O                ZnCO3:H2O
   461  ZnCl2                     ZnCl2(NH3)2
   463  ZnCl2(NH3)4               ZnCl2(NH3)6
   465  ZnCr2O4                   ZnF2
   467  ZnSO4:6H2O                ZnSO4:7H2O
   469  ZnSO4:H2O                 ZrO2


 liquids

     0  none

 * Note - (EQPT/pcrsg) The pure liquids block has
       not been written on the DATA1 and DATA1F files,
       because the EQ3NR and EQ6 codes presently do not
       treat non-aqeuous liquids.


 gases

     1  CH4(g)                    CO2(g)
     3  H2(g)                     H2O(g)
     5  HBr(g)                    HCl(g)
     7  HF(g)                     HNO3(g)
     9  N2O5(g)                   NO3(g)
    11  O2(g)


 solid solutions

     0  none
