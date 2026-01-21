 EQPT Species List (SLIST) File:


           Number of elements =    18
           Number of basis species =    28

 data0.fmt.R1
 PITZER THERMODYNAMIC DATABASE (12/04/2008)
 data0.fmt.R1
 
  This file contains two revisions to data0.fmt.R0:
 
  ***** Aphi = 0.392, not 0.39 *****
  ***** beta(1) for NaCl = 0.2664, not 0.2644 *****
 
 Revision 1: Harvie et al. (1984) gave Aphi = 0.39. As documented by L.N.
 Plummer, D.L. Parkhurst, G.W. Fleming, and S.A. Dunkle, A Computer
 Program Incorporating Pitzer's Equations for Calculation of Geochemical
 Reactions in Brines, U.S. Geological Survey, Water-Resources Investigations
 Report 88-4153 (see p. 3, the results of Harvie et al. (1984) are actually
 consistent with an Aphi value of 0.392 (see Plummer et al.,1988, p. 3).
 
 Revision 2: the beta(1) for NaCl of 0.2644 given by Harvie et al. (1984)
 appears to be a typographical error. The value of 0.2664 is given in the
 earlier paper by C.E. Harvie and J.H. Weare (1980, The prediction of mineral
 solubilities in natural waters: the Na-K-Mg-Ca-Cl-SO4-H2O system from zero
 to high concentration at 25C, Geochimica et Cosmochimica Acta, v. 44,
 p. 981-997; see Table 1, p. 987). There appears to be nothing of substance
 in the subsequent Weare and company literature describing a revision, just a
 different value appearing in a table in the 1984 paper. The 0.2664 value
 is also cited by Pitzer (1991, Chapter 3, Ion Interaction Approach: Theory
 and Correlation, p. 75-153 in Pitzer, K.S., ed., Activity Coefficients in
 Electrolyte Solutions, 2nd edition, CRC Press, Boca Raton; see Table 2,
 p. 100). This is also the value given by Plummer et al. (1988, p. 150).
 There may be no practical difference between the 0.2644 and 0.2664
 values.
 
  ************************************************************
 
 This datafile is a corrected translation of the "fmt_050405.chemdat"
 datafile used with the WIPP code FMT. The purpose of this datafile
 is to transfer the WIPP geochemistry model from FMT to EQ3/6.
 
 This datafile uses the standard form of Pitzer's equations relevant to
 the temperature range.
 
 Pitzer parameters are represented by the 25C-centric four-term temperature
 function given by:
 
 x(T) = a1 + a2*(1/T - 1/298.15) + a3*ln(T/298.15) + a4*(T - 298.15)
 
 where T is temperature in Kelvin and a1 through a4 denote the temperature
 function fitting coefficients for the temperature-dependent Pitzer
 parameters. This is the form used by the Yucca Mountain Project. The form
 itself is sufficient to accurately describe data up to a temperature of at
 least 250C, although the YMP database is only qualified up to 140C. At the
 present time, this WIPP database is for 25C only. The form of the equation
 is that a1 is the 25C value. Consequently, the a2, a3,and a4 constants
 can simply be assigned a value of zero.
 
 Some data blocks provide comments on the gathered Pitzer parameters and
 solid phase solubility data to make the user aware of any convention
 adopted in the data extraction/compilation process.
 
 The following comments apply to data0.fmt.R0:
 
   Data for O2(aq) and H2(aq) were added from YMP DTN: SN0302T0510102.002
   (data0.ypf.R1). Pitzer data were not included for H2(aq) because they
   were not available from the YMP DTN source. The O2(aq) data were added
   to allow calculations involving the oxygen fugacity, and this can be
   accomplished for essentially non-redox problems by adding a trace amount
   of O2(aq) or H2(aq). Until the requisite set of Pitzer data for H2(aq)
   have been added, this database should not be used to model systems with
   non-trace amounts of that species.
 
   Data for uranium and plutonium were removed, as these elements
   are not included in the present WIPP geochemistry model. These data
   are not needed for the EQ3/6-FMT code comparison study.
 
   The thermodynamic data for CO2(g) are based on data hard-wired into
   FMT v. 2.4.
 
 
 BEGIN CONFIGURATION DATA BLOCK
 Do not change the data in this block unless you know what you
 are doing.
 INTERPRET 500 AS NO DATA= NO
 SPARSE GRID RANGE CONDITION=IGNORE
 Pitzer data parameters:
    PITZER DATA BLOCK ORG.= NEW
    PITZER TEMP FUNCTION= LIVERMORE
    NO. OF PITZER TEMP FUNC TERMS= 4
 END CONFIGURATION DATA
 +--------------------------------------------------------------------


 element = O       , atwt =   15.99940
 element = Am      , atwt =  243.00000
 element = B       , atwt =   10.81100
 element = Br      , atwt =   79.90400
 element = C       , atwt =   12.01100
 element = Ca      , atwt =   40.07800
 element = Cl      , atwt =   35.45270
 element = H       , atwt =    1.00794
 element = K       , atwt =   39.09830
 element = Mg      , atwt =   24.30500
 element = N       , atwt =   14.00674
 element = Na      , atwt =   22.98977
 element = Np      , atwt =  237.00000
 element = P       , atwt =   30.97376
 element = S       , atwt =   32.06600
 element = Th      , atwt =  232.03810
 element = Null-   , atwt =    0.00000
 element = Null+   , atwt =    0.00000


 aqueous

     1  H2O                       Am+++
     3  B(OH)4-                   Br-
     5  HCO3-                     Ca++
     7  Cl-                       H+
     9  K+                        Mg++
    11  NO3-                      Na+
    13  NpO2+                     HPO4--
    15  SO4--                     Th++++
    17  NegIon                    PosIon
    19  O2(g)                     O2(aq)
    21  H2(aq)                    ClO4-
    23  Acetate-                  Citrate---
    25  EDTA----                  Oxalate--
    27  Lactate-                  NH3(aq)
    29  MgOH+                     HSO4-
    31  OH-                       CO3--
    33  CO2(aq)                   CaCO3(aq)
    35  MgCO3(aq)                 B(OH)3(aq)
    37  B3O3(OH)4-                B4O5(OH)4--
    39  CaB(OH)4+                 MgB(OH)4+
    41  H3PO4(aq)                 H2PO4-
    43  PO4---                    AmCO3+
    45  Am(CO3)2-                 Am(CO3)3---
    47  AmOH++                    Am(OH)2+
    49  Am(OH)3(aq)               AmCl++
    51  AmCl2+                    Am(CO3)4(5-)
    53  Am(SO4)2-                 AmSO4+
    55  Th(CO3)5(6-)              Th(OH)3(CO3)-
    57  Th(OH)4(aq)               Th(SO4)2(aq)
    59  Th(SO4)3--                NpO2CO3-
    61  NpO2(CO3)2---             NpO2(CO3)3(5-)
    63  NpO2OH(aq)                NpO2(OH)2-
    65  HAcetate(aq)              H3Citrate(aq)
    67  H2Citrate-                HCitrate--
    69  H4EDTA(aq)                H3EDTA-
    71  H2EDTA--                  HEDTA---
    73  H2Oxalate(aq)             HOxalate-
    75  HLactate(aq)              AmAcetate++
    77  AmCitrate(aq)             AmEDTA-
    79  AmOxalate+                AmLactate++
    81  Th(Acetate)2++            NpO2H2EDTA-
    83  Th(Lactate)2++            ThAcetate+++
    85  ThCitrate+                ThEDTA(aq)
    87  ThOxalate++               ThLactate+++
    89  NpO2Acetate(aq)           NpO2Citrate--
    91  NpO2EDTA---               NpO2Oxalate-
    93  NpO2Lactate(aq)           MgAcetate+
    95  MgCitrate-                MgEDTA--
    97  MgOxalate(aq)             CaAcetate+
    99  CaCitrate-                CaEDTA--
   101  CaOxalate(aq)             NpO2HEDTA--


 minerals

     1  AmOHCO3(c)                Am(OH)3(s)
     3  NaAm(CO3)2.6H2O(c)        AmPO4(c)
     5  ThO2(am)                  Th(SO4)2.9H2O(s)
     7  Th(SO4)2.8H2O(s)          Th(SO4)2.Na2SO4.6H2O
     9  Th(SO4)2.K2SO4.4H2O       Th(SO4)2.2K2SO4.2H2O
    11  2[Th(SO4)2.7/2K2SO4]      NpO2OH(aged)
    13  NpO2OH(am)                2[NaNpO2CO3.7/2H2O]
    15  Na3NpO2(CO3)2             KNpO2CO3
    17  K3NpO2(CO3)2              H2Oxalate.2H2O
    19  NaHOxalate.H2O            Na2Oxalate
    21  Anhydrite                 Aphthitalite/Glaserite
    23  Whewellite                Aragonite
    25  Arcanite                  Bischofite
    27  Bloedite                  Brucite
    29  Burkeite                  Calcite
    31  CaCl2.4H2O                CaOxychloride_A
    33  CaOxychloride_B           Carnallite
    35  Epsomite                  Gaylussite
    37  Glauberite                Gypsum
    39  Halite                    Hexahydrite
    41  Kainite                   Kalicinite
    43  Kieserite                 Leonite
    45  Labile_Salt               Magnesite
    47  Mg2Cl(OH)3.4H2O           Mercallite
    49  Mirabilite                Misenite
    51  Nahcolite                 Natron
    53  Nesquehonite              Picromerite/Schoenite
    55  Pirssonite                Polyhalite
    57  Portlandite               K2CO3.3/2H2O
    59  K8H4(CO3)6.3H2O           KNaCO3.6H2O
    61  K_Trona                   K3H(SO4)2
    63  Na3H(SO4)2                Na2CO3.7H2O
    65  Sylvite                   Syngenite
    67  Tachyhydrite              Thenardite
    69  Thermonatrite             Trona
    71  Borax                     B(OH)3
    73  K-Pentaborate(30C)        K-Tetraborate(30C)
    75  Na_Metaborate             Na_Pentaborate
    77  Teepleite(20C)            Dolomite
    79  Hydromagnesite5424        Hydromagnesite4323


 liquids

     0  none

 * Note - (EQPT/pcrsg) The pure liquids block has
       not been written on the DATA1 and DATA1F files,
       because the EQ3NR and EQ6 codes presently do not
       treat non-aqeuous liquids.


 gases

     1  CH4(g)                    CO2(g)
     3  H2(g)                     H2O(g)
     5  O2(g)


 solid solutions

     0  none
