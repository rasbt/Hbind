/* NOTE: It is strongly urged that you do not make changes to any values in this header file.  Proceed at your own risk.  To access the user adjustable parameters, please look at the params.h header file.  */


#ifndef _DEFS_H
#define _DEFS_H

#include <params.h> /* params.h is the header file which contains user adjustable variables */

/* Enable output of bumps -- listing and mol2 files */
#undef OUTPUT_BUMPS
#undef BUMP_FILES
#undef OUTPUT_ALL_BUMPS /* ignore max cutoff values for bump output */
/*************************************/
/*************************************/

/* maximal number of atoms that can be read in from a pdb-file */
#define  MAX_PDB_ATOMS              8000

/* maximal number of water atoms from a pdb-file, added by Litian 07/24/03 */
#define  MAX_PDB_WATER_ATOMS        800

/* maximal number of water atoms from a seperate file, added by Litian 08/28/03 */
#define  MAX_WATER_ATOMS            5000

/* maximal number of residues which can be read from a pdb-file */
#define  MAX_PDB_RESIDUES           1000

/* maximal length of a line in a pdb-file */
#define  MAX_PDB_LINELENGTH          200

/* maximal length of a line in a template file */
#define  MAX_TEMPLATE_LINELENGTH     200

/* maximal length of a line in a parameter file */
#define  MAX_PARAMETER_LINELENGTH    128

/* maximal number of template-points */
#define  MAX_TEMPLATE_POINTS         512

/* define the number of buckets in the type hash table */
#define  NUMBER_OF_TYPE_HASH_CLASSES  10
#define  BUCKETS_PER_A_PERIMETER      2.0
#define  BUCKETS_PER_A_LENGTH         4.0


/* maximal number of hbonds per template point */
#define  MAX_HBONDS_PER_TEMPLATE_POINT   5

/* maximal number of ligand atoms, also taken as an upper bound for the
   number of interaction points */
#define  MAX_LIGAND_ATOMS                550

/* maximal number of residues of a ligand */
#define  MAX_LIGAND_RESIDUES             50

/* maximal number of atoms in a single residue */
#define  MAX_RESIDUE_ATOMS               14 /* changed from 12 - Sameer 13/Mar/05*/

/* maximal size of distance matrices, i.e. maximal number of interaction
   points to match */
#define  MAX_DIST_SIZE                   10

/* when matching ligand interaction centers onto template points, the
   matching for the hydrophobic points has to be less exact, h-bond
   point matchings are weight 1.0 during the DME, the least-squares-fit,
   and RMSD computations, this value is the relative weight for the
   hydrophobic points */
#define  HPHOB_MATCH_WEIGHT              0.3

/* maximal number of interaction-pairs between ligand and protein */
#define  MAX_INTER_PAIRS                 16

/* maximal number of H-bonds between protein and ligand */
#define  MAX_TOTAL_HBONDS                500
#define  MAX_TOTAL_METAL_HBONDS          100

/* maximal number of salt bridges between protein and ligand */
#define  MAX_TOTAL_SALT_BRIDGES          50

/* maximale length of a compound name */
#define  MAX_COMPOUNDNAME_LEN            128

/* maximal number of compounds in a pts file */
#define  MAX_NUMBER_OF_PTS_COMPOUNDS     7000

/* maximal number of points for a set of compounds in a pts file */
#define  MAX_NUMBER_OF_TOTAL_COMPOUND_INTERACTION_POINTS 200000

/* maximal number of conserved waters in the ligand binding site */
#define  MAX_BINDING_SITE_WATERS         1000

/* maximal number of anchor fragment bumps that are corrected for a ligand */
#define  MAX_ANCHOR_BUMPS                10

/* maximal number of bumps accepted for all side chains */
#define  MAX_SIDE_CHAIN_BUMPS            20

/* maximal number of iterations when trying to get rid of bumps */
#define  MAX_ANCHOR_CORRECTION_STEPS     100

/* maximal number of iterations to get rid of side-chain bumps */
#define  MAX_NUMBER_OF_SIDE_CHAIN_CORRECTIONS 8

/* maximal number of tries to resolve ligand bumps */
#ifndef ROTAMER_SEARCH
#define  MAX_LIGAND_BUMP_RESOLVE_ITERATIONS  10
#endif
#ifdef ROTAMER_SEARCH
#define  MAX_LIGAND_BUMP_RESOLVE_ITERATIONS  30
#define  MAX_UNSAT_SIDE_CHAINS  			 50 /* changed from 10 - Sam 20/Jan/05*/
#endif

/* maximal number of translations for a bumping water */
#define  MAX_NUMBER_WATER_UNBUMP_ITERATIONS  10

/* maximal number of rotamers for a single side-chain */
#define  MAX_NUMBER_OF_ROTAMERS          81

/* the relative weight of matched side-chain atoms during distance geometry
   computations */
#define  SIDE_CHAIN_DG_WEIGHT            0.5

/* the relative weight of matched side-chain atoms during ligand
   transformation */
#define  SIDE_CHAIN_TRANSFORM_WEIGHT     0.5

/* the maximal allowed distance of the ligand centroid from the center
   of the template after transforming the ligand into the binding site */
#define  MAX_CENTER_DISTANCE             8.0

/* maximal number of neighbors to a single atom */
#define  MAX_NEIGHBOR_ATOMS              300

/* distance up to which atoms are considered as being neighbors, i. e. this
   defines the set of atoms checked for intramolecular bumps */
#define  NEIGHBOR_DISTANCE               8.0

/* maximal distance between atom centers when checking for hydro-compl */
#define  HYDRO_DIST                      4.5
#define  HYDRO_DIST_2			20.2545  /* (HYDRO_DIST + 0.0005)**2 */

/* maximal distance to look for polar ligand atoms for interface target
   atoms */
#define  HPHIL_MATCH_DIST                5.0

/* the range of donor acceptor distances in which a H-bond can be formed */
#define  MIN_HBOND_LENGTH_HYDROGEN       1.5

/*Changed to 2.4 for donor-acceptor project
#define  MIN_HBOND_LENGTH */
#define  MIN_HBOND_LENGTH                2.4


#define  MAX_HBOND_LENGTH                3.5


#define  MIN_HBOND_LENGTH_HYDROGEN_2     2.25
/* (MIN_HBOND_LENGTH - 0.0005)**2 */

/*Changed to 2.4^2 for donor-acceptor project
#define  MIN_HBOND_LENGTH_MINUS_TOL_2    6.2475 */
#define  MIN_HBOND_LENGTH_MINUS_TOL_2     5.7595

/*Changed to 2.4^2 for donor-acceptor project
#define  MIN_HBOND_LENGTH_2              6.25 */
#define  MIN_HBOND_LENGTH_2              5.7595

#define  MAX_HBOND_LENGTH_2             12.25
/* (MAX_HBOND_LENGTH + 0.0005)**2 */
#define  MAX_HBOND_LENGTH_PLUS_TOL_2    12.2535


/* added to check metal hbond  12/03/2003 */
#define  MIN_METAL_1_HBOND_LENGTH        2.0
#define  MIN_METAL_2_HBOND_LENGTH        1.7
#define  MAX_METAL_1_HBOND_LENGTH        2.9
#define  MAX_METAL_2_HBOND_LENGTH        2.6

#define  MIN_METAL_1_HBOND_LENGTH_2      4.0
#define  MIN_METAL_2_HBOND_LENGTH_2      2.89
/* (MAX_METAL_1_HBOND_LENGTH + 0.0005)**2 */
#define  MAX_METAL_1_HBOND_LENGTH_PLUS_TOL_2 8.4129
/* (MAX_METAL_2_HBOND_LENGTH + 0.0005)**2 */
#define  MAX_METAL_2_HBOND_LENGTH_PLUS_TOL_2 6.7626
/*#define  INFINITE_DIST                   100.0 */

#define  MAX_TOTAL_METAL_HBOND           100

/* the range for the donor-hydrogen-acceptor angle in a H-bond */
#define  MIN_HYDROGEN_ANGLE              120
#define  MAX_HYDROGEN_ANGLE              180
#define  FAILURE_ANGLE                   300
#define  MIN_PREACCEPTOR_ANGLE            90
#define  MAX_PREACCEPTOR_ANGLE           180

/* the maximal distance between two complementary charged atoms in a
   salt bridge */
#define  MIN_SALT_BRIDGE_LENGTH          2.5
#define  MAX_SALT_BRIDGE_LENGTH          4.5

/* (MIN_SALT_BRIDGE_LENGTH - 0.0005)**2 */
#define  MIN_SALT_BRIDGE_LENGTH_MINUS_TOL_2 6.2475
/* (MAX_SALT_BRIDGE_LENGTH + 0.0005)**2 */
#define  MAX_SALT_BRIDGE_LENGTH_PLUS_TOL_2 20.2545

/* Sameer Jan/09/05 - looser hbond for rotamer search
* #ifdef HBOND_LOOSENING
*	#undef   MAX_HBOND_LENGTH
*	#define  MAX_HBOND_LENGTH           ( 3.5 + HBOND_LOOSENING )
*
*	#undef  MAX_SALT_BRIDGE_LENGTH
*	#define  MAX_SALT_BRIDGE_LENGTH     ( 4.5  + HBOND_LOOSENING )
* #endif
*/
/* the number of checks which can rule out infeasible ligands */
#define  NUMBER_OF_CHECKS                9

/* the number of template point types: DONOR, ACCEPTOR, DONEPTOR, HYDROPHOB */
#define  NUMBER_OF_INTERACTION_TYPES     4

/* the maximal number of side-chain rotations checked in the recursive
   side-chain bump-resolving routine */
#define  MAX_ROTATION_DEPTH              3

/* we do not care during the bump checks for any atoms with a distance
   larger than this */
#define  DONT_CARE_BUMP_DISTANCE         4.5
#define  DONT_CARE_BUMP_DISTANCE_SQ     20.25

/* van der Waals radius of a water molecule */
#define  WATER_RAD                       1.36

/* the hydrophilicity value for water */
#define  WATER_HYDRO                     580

/* the number of bonds that have to be considered during the side-chain
   unbumping, this is an average of two bonds per bump */
#define  MAX_UNBUMP_BONDS        3*MAX_SIDE_CHAIN_BUMPS

/* the maximal number of dependent bonds for a single */
#define  MAX_UNBUMP_DEPENDENCIES 10

#define  NUMBER_OF_FILTERS       7
#define  TRIANGLE_MATCH          0
#define  DME                     1
#define  RMS_DEVI                2
#define  BUMP_ANCHOR             3
#define  BUMP_SIDE_CHAIN         4
#define  SCORING                 5

#define  PASSED                  6

/* indices for coordinate vector */
#define  X            0
#define  Y            1
#define  Z            2

/* interaction types */
#define  NOTHING      16384 /* Changed from 0 to 1024 since sum-conflict of NOTHING + DONOR_HYDROGEN = DONEPTOR_DONEPTOR - Sameer 18/Feb/05*/
#define  ACCEPTOR          1
#define  DONOR             4
#define  HYDROPHOB        16
#define  DONEPTOR         64
#define  DONOR_HYDROGEN  256
#define  METAL_1        1024
#define  METAL_2        4096

/* summed interaction types */
#define  ACCEPTOR_ACCEPTOR         2
#define  DONOR_DONOR               8
#define  HYDROPHOB_HYDROPHOB      32
#define  ACCEPTOR_DONOR            5
#define  ACCEPTOR_DONEPTOR        65
#define  DONOR_DONEPTOR           68
#define  DONEPTOR_DONEPTOR       128
#define  ACCEPTOR_DONOR_HYDROGEN 257
#define  DONEPTOR_DONOR_HYDROGEN 320
#define  ACCEPTOR_METAL_1       1025
#define  ACCEPTOR_METAL_2       4095
#define  DONEPTOR_METAL_1       1028
#define  DONEPTOR_METAL_2       4104

#define  NO_MATCH     0
#define  MATCH        1

#define  BUMP        -1
#define  NO_BUMP      1

#define  FATAL_FAILURE    -1
#define  FAILURE           0
#define  SUCCESS           1
#define  CHECK_AGAIN       2
#define  CHECK_NEXT_BOND   3
#define  RESTART_SKIP_MOL  4
#define  FALSE             0
#define  TRUE              1

#define  UNKNOWN              -1
#define  NO_MORE              -1
#define  UNKNOWN_INDEX         0
#define  UNKNOWN_TARGET_INDEX -1
#define  NO_ANGLE         -999.0

#define  TARGET       1
#define  LIGAND       3

#define  SHORT        0
#define  LONG         1

#define  ARRAY_DUMMY -1

#define  ATOMS_ONLY   1
#define  ALSO_HETERO  2

/* atom levels for side-chain flexibility */
#define  ALPHA    0
#define  BETA     1
#define  GAMMA    2
#define  DELTA    3
#define  EPSILON  4
#define  ZETA     5
#define  ETA      6

#define  MAX_NUMBER_OF_MOL2_ATOMS      500
#define  MAX_NUMBER_OF_MOL2_BONDS      500
#define  MAX_NUMBER_OF_FLEXIBLE_BONDS  500
#define  MAX_NUMBER_OF_CARBON_RINGS     50
#define  MAX_MOL2_LINELENGTH           128
#define  MAX_LEN_MOL2_COMPOUND_NAME    128
#define  MAX_RULES_PER_BONDED_ATOM       5
#define  MAX_NUMBER_OF_FLEX_BOND_RULES  40
#define  MAX_NUMBER_OF_HYD_RULES        10
#define  MAX_NEIGHBORS                  20
#define  MAX_CARBON_RING_SIZE            6


#define  FALSE        0
#define  TRUE         1
#define  NO           0
#define  YES          1

#define  ANY         99

#define  ATOM1        0
#define  ATOM2        1

#define  OFF_ANCHOR   0
#define  ANCHOR       1
#define  COVALENT    99

#define  NO_FLEX_BOND 0
#define  FLEX        42
#define  FLEX_REVERSE 2
#define  RIGID        3

#define  REVERSE      0
#define  STRAIGHT     1
#define  END          2

#define  FIXED        0
#define  ROTATABLE    1

#define  FAILURE      0
#define  SUCCESS      1

#define  DISPLACED    0
#define  CONSERVED    1
#define  POLAR_DISPLACED 2
#define  WATER_OVERLAP_NEEDS_TO_BE_RESOLVED 3

#define  UNKNOWN     -1
#define  SINGLE       1
#define  DOUBLE       2
#define  TRIPLE       3
#define  QUADRUPLE    4
#define  AROMATIC     5
#define  AMIDE        6
#define  DELOCALIZED  7
#define  CYCLE_BOND  10

#define  NO_CYCLE     0
#define  CYCLE        1
#define  UNVISITED    6
#define  VISITED     20
#define  CYCLE_START 30

#define  CLOSE        0
#define  DISTANT      1

#define  JUST_SCORE  0
#define  DOCK_AND_SCORE 1
#define  RECORD_UNSAT 998 /* Sameer 29/Nov/04 */
#define  RECORD_PRE_ROT_SRCH_TERMS 997 /* Sameer 03/Apr/05 */

/* Keep this for now just because the check could cause slightly different results if we used a different value, however this is an arbitrary value */
#define  SMALL_FLOAT  1.0e-5
/*
#define  SMALL_DOUBLE 1.0e-10
*/

/* parameter-file entries */
#define  INVALID_ENTRY_OR_ERROR            0
#define  FLEX_BOND_DEFN_FILE               1
#define  DATA_ROOT                         2
#define  DME_THRESHOLD                     3
#define  RMS_THRESHOLD                     4
#define  ANCHOR_TRANSLATION                5
#define  ANCHOR_OVERLAP                    6
#define  SIDE_CHAIN_OVERLAP                7
#define  INTRA_OVERLAP                     8
#define  INTERMEDIATELY_TOLERATED_OVERLAP  9
#define  FINALLY_TOLERATED_OVERLAP        10
#define  FINALLY_TOLERATED_MAX_BUMP       11
#define  SCORE_CUTOFF                     12
#define  MAX_TEMPLATE_TRIANGLES           13
#define  MATCH_2_KEY_POINTS               14
#define  GROUP_CONFORMERS                 15
#define  RESTART_MOLECULE                 16
#define  IGNORE                           99

/** 18/Jan/05 - Sameer.
 * Use this macro instead of magic
 * number 1.0 in rotate_unbump.c
 */
#define  ALLOWED_INTRA_OVERLAP_DURING_ROTATIONS 1.0

/* batch-file entries */
#define  LOCAL_PARAMETER_FILE 40
#define  PROTEIN              41
#define  TEMPLATE             42
#define  DATABASE             43

/* all atom types */
/* To add a new PDB atom type, add a line below (with the next number),
   add a radius definition in inc/vdwrad.h and add the parsing code to read_pdb.c
   (instructions at 'ADD ATOM') */
/* To add a new MOL2 atom type, add a line below (with the next number),
   add a radius definition in inc/vdwrad.h, add the type parsing code to
   assign_type.c (instructions at 'ADD MOL2 ATOM'), and add the radius
   assignment code to read_mol2.c (also at 'ADD MOL2 ATOM') */
#define  CA        1
#define  CB        2
#define  C         3
#define  O         4
#define  N         5
#define  CD        6
#define  CD1       7
#define  CD2       8
#define  CE        9
#define  CE1      10
#define  CE2      11
#define  CE3      12
#define  CG       13
#define  CG1      14
#define  CG2      15
#define  CH2      16
#define  CZ       17
#define  CZ2      18
#define  CZ3      19
#define  ND1      20
#define  ND2      21
#define  NE       22
#define  NE1      23
#define  NE2      24
#define  NH1      25
#define  NH2      26
#define  NZ       27
#define  OD1      28
#define  OD2      29
#define  OE1      30
#define  OE2      31
#define  OG       32
#define  OG1      33
#define  OH       34
#define  OXT      35
#define  SD       36
#define  SG       37
#define  S        40
#define  P        41
#define  H        42
#define  SI       43
#define  CO       44    /* metal type 2, hbond distance 2.6 */
#define  AG       45
#define  NI       46    /* metal type 1, hbond distance 2.0 */
#define  BR       47
#define  TL       48
#define  OS       49
#define  CL       50
#define  LA       51
#define  K        52    /* metal type 3, hbond distance 2.9 */
#define  FE       53    /* metal type 2, hbond distance 2.6 */
#define  DUMMY    54    /* dummy atom */
#define  MO       55
#define  I        56
#define  AS       57
#define  U        58
#define  F        59
#define  RE       60
#define  RU       61
#define  ER       62
#define  LI       63
#define  MN       64    /* metal type 2, hbond distance 2.6 */
#define  PD       65
#define  CU       66    /* metal type 2, hbond distance 2.6 */
#define  AL       67
#define  MG       68    /* metal type 2, hbond distance 2.6 */
#define  NA       69    /* metal type 3, hbond distance 2.9 */
#define  CR       70
#define  CS       71
#define  TI       72
#define  ZN       73    /* metal type 2, hbond distance 2.6 */
#define  RB       74
#define  V        75
#define  RH       76
#define  GE       77
#define  ZR       78
#define  GA       79
#define  Y_       80
#define  CA_M     81    /* metal type 2, this is metal "CA", not alpha carbon */
#define  AU       82
#define  BE       83
#define  SR       84
#define  BA       85
#define  HG       86
#define  B        87
#define  IN       88
#define  SN       89
#define  SB       90
#define  BI       91
#define  SE       92
#define  TE       93

/***** Added by PCS -- 22-Mar-00 *****/
#define  OC       94  /* Too Close Oxygen - Very small radius */
/* added by Litian He 03/17/2004 */
#define  AD1      95   /* we have these 2 because the crystallographic can't decide O/N in ASN */
#define  AD2      96
#define  AE1      97   /* we have these 2 because the crystallographic can't decide O/N in GLN */
#define  AE2      98

/*************************************/
#define  HETATM  999
#define  WATER  1000

#define  ADDED     0
#define  SP1       1    /*  .1   */
#define  SP2       2    /*  .2   */
#define  SP3       3    /*  .3   */
/*       O         4        .o   */
#define  SP4       5    /*  .4   */
#define  AR        6    /*  .ar  */
#define  CO2       7    /*  .co2 */
#define  AM        8    /*  .am  */
#define  PL3       9    /*  .pl3 */
#define  CAT      10    /*  .cat */
#define  O2       11    /*  .o2  */
#define  TH       12    /*  .th  */
#define  AMBIG    13    /*  .ambiguous  */
#define  HISN     14    /*  .HIS-nitrogen  */
/*       OH       34        .oh  */

/* all residues */
#define  ALA      0
#define  ARG      1
#define  ASN      2
#define  ASP      3
#define  CYS      4
#define  GLN      5
#define  GLU      6
#define  GLY      7
#define  HIS      8
#define  ILE      9
#define  LEU     10
#define  LYS     11
#define  MET     12
#define  PHE     13
#define  PRO     14
#define  SER     15
#define  THR     16
#define  TRP     17
#define  TYR     18
#define  VAL     19
#define  ACE     20
#define  PCA     21

#define  HOH     22

/* this should be defined in stdio.h */
#ifndef FILENAME_MAX
#define FILENAME_MAX    1024
#endif

/* linux and sgi compatibility - usually is defined in netdb.h or sys/param.h */
#ifndef MAXHOSTNAMELEN
#define MAXHOSTNAMELEN 256
#endif

/** RSK - Aug 13 2002 **/
#define MAX_SUBST_NAME_LENGTH 5
#define MAX_SUBSTS          100
#define MAX_SUBST_LEN       256

#endif

/** added by Litian 07/23/2003 **/
/* maximal distance from neighbors to buried atom */
#define  MAX_BURY_DIST               4.5 /* Sameer 02/17/05 - as recommended by MIZ */
/* minimal distance from neighbors to buried atom */
#define  MIN_BURY_DIST               1.9
/* minimal number of neighbors of the buried atom */
#define  MIN_NUMBER_BURY_NEIGHBOR    8

/* In the future rounding should be performed correctly
 * In the event of that happening there is no need to add and
 * subtract 0.0005 from distances.  Doing so makes the code more
 * confusing and why not add 0.001, ....
 */
#define MAX_BURY_DIST_2 	20.25
#define MIN_BURY_DIST_2		 3.61


/* to calculate the score after we got weights from each parameter in
   the scoring funcions defined in src/hbind/scoring/score_define.txt */
#define INITIAL_VALUE                   -1000000
#define MAX_NUM_OF_PARAMETERS           20
#define MAX_NUM_OF_SCORE_FUNCTIONS      67
#define MAX_ENTRY_PER_FUN               6
#define MAX_PARAM_PER_ENTRY             3
#define MAX_LEN_PER_ENTRY               10
#define MAX_DEFINITION_FILE_PATH        100


/* added to rename the different types of atoms
    -- Litian, 09/03/2003 */
#define INITIAL                          0
#define INTERFACIAL_ATOM                -1
#define METAL_DIRECT_HBOND               1
#define SALT_BRIDGE                      2
#define DIRECT_HBOND                     3
#define ONE_WATER_MEDIATED_HBOND         4
#define TWO_WATER_MEDIATED_HBOND         5
#define METAL_ONE_WATER_MEDIATED_HBOND   6
#define METAL_TWO_WATER_MEDIATED_HBOND   7
#define UNSAT_POLAR                      8
#define UNSAT_CHARGE                     9
#define REPULSIVE_POLAR                  10
#define REPULSIVE_CHARGE                 11
#define SATISFIED_SOMEHOW     			 12

#define WATER_HBOND                      21

#define LIGAND_WATER_HBOND               22
#define INTRA_LIGAND_HBOND               23
#define INTRA_LIGAND_SALT_BRIDGE         24

#define TARGET_WATER_HBOND               25
#define INTRA_TARGET_HBOND               26
#define INTRA_TARGET_SALT_BRIDGE         27

/* maximal distance for unsatisfied polar/charge and repulsive interfactions */
#define MAX_UNSAT_REPULSIVE_DIST         3.5

/***** Added by PCS -- 11-Jan-00 *****/
/* Score Weightings */
#define HPHOB_SCORE_WEIGHT               0.586
#define HBOND_SCORE_WEIGHT               2.755
#define SALT_BR_SCORE_WEIGHT             2.755
/*************************************/

/* added to give weight to different interactions
    -- Litian, 08/12/2004
*/
#define WEIGHT_OF_HBONDS                    -0.6727    /* term I*/
#define WEIGHT_OF_METAL_HBONDS              -0.2321    /* term G*/
#define WEIGHT_OF_SALT_BRIDGE               -0.2911    /* term H*/
#define WEIGHT_OF_HPHOB_HPHOB_CONTACT       -0.091     /* term S*/
#define WEIGHT_OF_POLAR_CHARGE              -0.282     /* term I+G+H*/
#define WEIGHT_OF_UNSAT_REPULSIVE           -0.049     /* term J+K+L+M*/
#define WEIGHT_OF_EXPOSED_LIG_HPHOB         -0.398     /* term P*/

/*  -- Litian, 10/13/2004 */
#define POSITIVE_CHARGE                     1.0
#define NEGATIVE_CHARGE                    -1.0
#define MINIMAL_CHARGE                      0.5
#define HYDROPHOB_VALUE_SHIFT             235.0

/* Out of bounds values used for initializations */
#define PRE_ROT_CONFIG 0

#ifdef INTERACTION_OPPORTUNITIES
#define MAX_INTERACTION_OPPORTUNITIES 20
#endif


#if 0
#define SML_JUNK_FLOAT -899.77

#define BIG_JUNK_FLOAT 999.77
#define BIG_JUNK_DOUBLE 999.77
#define SML_JUNK_DOUBLE -899.77
#define BIG_JUNK_INT 777
#define SML_JUNK_INT -666
#define MIN_RMSD_SCORE -999.0 /* Sameer 05/Mar/05*/
#endif


/* For correct calculation of metal-bonds and polar
 * repulsive contacts in presence of a metal, following
 * macros coniditonally compile in some code. If change work,
 * the code should be made part of distribution, unconditionally
 * Sameer - 24/May/2005
 */
#define CORRECT_SF_L_TERM
#define NON_METALBONDED_REPULSIVE
#ifdef CORRECT_SF_L_TERM
#define MAX_TARGET_METAL_ATOMS 20
#endif


#define NO_OF_SCORE_COMPONENT_TERMS 53
