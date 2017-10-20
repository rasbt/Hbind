#ifndef _VDWRAD_
#define _VDWRAD_

/***** PDB file van der Waals Radii Definitions *****/
/***** Small adjustments made by L Kuhn, 7/31/07, to match *****/
/***** data in Table XI: The Minimum Radii of Atoms, in Li & Nussinov (1998)
"A Set of van der Waals and Coulombic Radii of Protein Atoms for Molecular
and Solvent-Accessible Surface Calculation, Packing Evaluation, and Docking",
Proteins 21, 111-127. *****/

/***** Instructions to add a new atom type *****/
/* Add a radius definition below, then add the atom type
   assignment in read_pdb.c (see that file at 'ADD ATOM'),
   and add the atom name definition in inc/defs.h (also at
   'ADD ATOM' in that file */

#define RAD_ACEO    1.26
#define RAD_ACEC    1.62
#define RAD_ACECA   1.73
#define RAD_ALAN    1.43
#define RAD_ALACA   1.73
#define RAD_ALAC    1.62
#define RAD_ALAO    1.26
#define RAD_ALAOXT  1.24
#define RAD_ALACB   1.66
#define RAD_ARGN    1.43
#define RAD_ARGCA   1.73
#define RAD_ARGC    1.62
#define RAD_ARGO    1.26
#define RAD_ARGOXT  1.24
#define RAD_ARGCB   1.69
#define RAD_ARGCG   1.68
#define RAD_ARGCD   1.62
#define RAD_ARGNE   1.34
#define RAD_ARGCZ   1.62
#define RAD_ARGNH1  1.34
#define RAD_ARGNH2  1.34
#define RAD_ASNN    1.43
#define RAD_ASNCA   1.73
#define RAD_ASNC    1.62
#define RAD_ASNO    1.26
#define RAD_ASNOXT  1.24
#define RAD_ASNCB   1.69
#define RAD_ASNCG   1.67
#define RAD_ASNOD1  1.20
#define RAD_ASNND2  1.31
#define RAD_ASNAD1  1.25   /* average of OD1 and ND2 */
#define RAD_ASNAD2  1.25   /* average of OD1 and ND2 */
#define RAD_ASPN    1.43
#define RAD_ASPCA   1.73
#define RAD_ASPC    1.62
#define RAD_ASPO    1.26
#define RAD_ASPOXT  1.24
#define RAD_ASPCB   1.69
#define RAD_ASPCG   1.65
#define RAD_ASPOD1  1.24
#define RAD_ASPOD2  1.24
#define RAD_CYSN    1.43
#define RAD_CYSCA   1.73
#define RAD_CYSC    1.62
#define RAD_CYSO    1.26
#define RAD_CYSOXT  1.24
#define RAD_CYSCB   1.69
#define RAD_CYSSG   1.54
#define RAD_GLNN    1.43
#define RAD_GLNCA   1.73
#define RAD_GLNC    1.62
#define RAD_GLNO    1.26
#define RAD_GLNOXT  1.24
#define RAD_GLNCB   1.69
#define RAD_GLNCG   1.62
#define RAD_GLNCD   1.67
#define RAD_GLNOE1  1.20
#define RAD_GLNNE2  1.31
#define RAD_GLNAE1  1.25   /* average of OE1 and NE2 */
#define RAD_GLNAE2  1.25   /* average of OE1 and NE2 */
#define RAD_GLUN    1.43
#define RAD_GLUCA   1.73
#define RAD_GLUC    1.62
#define RAD_GLUO    1.26
#define RAD_GLUOXT  1.24
#define RAD_GLUCB   1.69
#define RAD_GLUCG   1.62
#define RAD_GLUCD   1.65
#define RAD_GLUOE1  1.24
#define RAD_GLUOE2  1.24
#define RAD_GLYN    1.43
#define RAD_GLYCA   1.68   /* CH2 atom type in Li & Nussinov */
#define RAD_GLYC    1.62
#define RAD_GLYO    1.26
#define RAD_GLYOXT  1.24
#define RAD_HISN    1.43
#define RAD_HISCA   1.73
#define RAD_HISC    1.62
#define RAD_HISO    1.26
#define RAD_HISOXT  1.24
#define RAD_HISCB   1.69
#define RAD_HISCG   1.62
#define RAD_HISND1  1.38
#define RAD_HISCD2  1.50
#define RAD_HISCE1  1.50
#define RAD_HISNE2  1.38
#define RAD_ILEN    1.43
#define RAD_ILECA   1.73
#define RAD_ILEC    1.62
#define RAD_ILEO    1.26
#define RAD_ILEOXT  1.24
#define RAD_ILECB   1.82
#define RAD_ILECG1  1.68
#define RAD_ILECG2  1.66
#define RAD_ILECD1  1.66
#define RAD_LEUN    1.43
#define RAD_LEUCA   1.73
#define RAD_LEUC    1.62
#define RAD_LEUO    1.26
#define RAD_LEUOXT  1.24
#define RAD_LEUCB   1.69
#define RAD_LEUCG   1.82
#define RAD_LEUCD1  1.66
#define RAD_LEUCD2  1.66
#define RAD_LYSN    1.43
#define RAD_LYSCA   1.73
#define RAD_LYSC    1.62
#define RAD_LYSO    1.26
#define RAD_LYSOXT  1.24
#define RAD_LYSCB   1.69
#define RAD_LYSCG   1.68
#define RAD_LYSCD   1.68
#define RAD_LYSCE   1.62
#define RAD_LYSNZ   1.22
#define RAD_METN    1.43
#define RAD_METCA   1.73
#define RAD_METC    1.62
#define RAD_METO    1.26
#define RAD_METOXT  1.24
#define RAD_METCB   1.69
#define RAD_METCG   1.68
#define RAD_METSD   1.67
#define RAD_METCE   1.66
#define RAD_PHEN    1.43
#define RAD_PHECA   1.73
#define RAD_PHEC    1.62
#define RAD_PHEO    1.26
#define RAD_PHEOXT  1.24
#define RAD_PHECB   1.69
#define RAD_PHECG   1.62
#define RAD_PHECD1  1.62
#define RAD_PHECD2  1.62
#define RAD_PHECE1  1.62
#define RAD_PHECE2  1.62
#define RAD_PHECZ   1.62
#define RAD_PRON    1.43
#define RAD_PROCA   1.73
#define RAD_PROC    1.62
#define RAD_PROO    1.26
#define RAD_PROOXT  1.24
#define RAD_PROCB   1.69
#define RAD_PROCG   1.68
#define RAD_PROCD   1.68
#define RAD_SERN    1.43
#define RAD_SERCA   1.73
#define RAD_SERC    1.62
#define RAD_SERO    1.26
#define RAD_SEROXT  1.24
#define RAD_SERCB   1.69
#define RAD_SEROG   1.30
#define RAD_THRN    1.43
#define RAD_THRCA   1.73
#define RAD_THRC    1.62
#define RAD_THRO    1.26
#define RAD_THROXT  1.24
#define RAD_THRCB   1.82
#define RAD_THROG1  1.30
#define RAD_THRCG2  1.66
#define RAD_TRPN    1.43
#define RAD_TRPCA   1.73
#define RAD_TRPC    1.62
#define RAD_TRPO    1.26
#define RAD_TRPOXT  1.24
#define RAD_TRPCB   1.69
#define RAD_TRPCG   1.62
#define RAD_TRPCD1  1.62
#define RAD_TRPCD2  1.62
#define RAD_TRPNE1  1.35
#define RAD_TRPCE2  1.62
#define RAD_TRPCE3  1.62
#define RAD_TRPCZ2  1.62
#define RAD_TRPCZ3  1.62
#define RAD_TRPCH2  1.62
#define RAD_TYRN    1.43
#define RAD_TYRCA   1.73
#define RAD_TYRC    1.62
#define RAD_TYRO    1.26
#define RAD_TYROXT  1.24
#define RAD_TYRCB   1.69
#define RAD_TYRCG   1.62
#define RAD_TYRCD1  1.62
#define RAD_TYRCD2  1.62
#define RAD_TYRCE1  1.62
#define RAD_TYRCE2  1.62
#define RAD_TYRCZ   1.62
#define RAD_TYROH   1.30
#define RAD_VALN    1.43
#define RAD_VALCA   1.73
#define RAD_VALC    1.62
#define RAD_VALO    1.26
#define RAD_VALOXT  1.24
#define RAD_VALCB   1.82
#define RAD_VALCG1  1.66
#define RAD_VALCG2  1.66
#define RAD_H       1.00 /* Not from Li & Nussinov */

/***** MOL2 file van der Waals Radii Definitions *****/

/***** Instructions to add a new atom type *****/
/* Add a radius definition below, then add the atom type
   assignment in assign_type.c (see that file at 'ADD MOL2 ATOM'),
   add the atom radius assignment code in read_mol2.c (at
   'ADD MOL2 ATOM RAD'), and add the atom name definition in
   inc/defs.h (at 'ADD MOL2 ATOM') */
/* NOTE: toneroma 2009_02_02 - New metal handling doesn't use vdwrad for bump check, but instead uses the minimum interaction distance */


#define MRAD_AG       1.550
#define MRAD_AL       1.500  /* Batsanov et al has this value as 1.8 20MAR08 toneroma - previous value 1.500, Amber/Dock5 has this value as 1.17 */
#define MRAD_AS       0.830  /* Batsanov et al 20MAR08 toneroma - previous value 0.830*/
#define MRAD_AU       1.900  /* Batsanov et al 20MAR08 toneroma - previous value none */
#define MRAD_B        1.700  /* Batsanov et al 20MAR08 toneroma - previous value none */
#define MRAD_BA       3.100  /* Batsanov et al 20MAR08 toneroma - previous value none */
#define MRAD_BE       1.900  /* Batsanov et al 20MAR08 toneroma - previous value none */
#define MRAD_BI       1.800  /* Batsanov et al 20MAR08 toneroma - previous value none */
#define MRAD_BR       1.950  /* Batsanov et al has this value as 1.9  20MAR08 toneroma */
#define MRAD_CSP1     1.850
#define MRAD_CSP2     1.850
#define MRAD_CSP3     1.800
#define MRAD_CCAT     1.800
#define MRAD_CAR      1.800
#define MRAD_CDEF     1.850
#define MRAD_CA       1.600  /* Batsanov et alhas this value as 2.5 20MAR08 toneroma */
#define MRAD_CD       1.750  /* Batsanov et al 20MAR08 toneroma - previous value 1.750 */
#define MRAD_CE       1.860
#define MRAD_CL       2.030  /* Batsanov et al has this value as 1.8  20MAR08 toneroma */
#define MRAD_CO       1.130
#define MRAD_CR       1.130
#define MRAD_CS       3.010  /* Batsanov et al 20MAR08 toneroma - previous value 3.010 */
#define MRAD_CU       1.150  /* Batsanov et al 20MAR08 toneroma - previous value 1.150 */
#define MRAD_ER       1.590
#define MRAD_F        1.550  /* Batsanov et al has this value as 1.5 20MAR08 toneroma */
#define MRAD_FE       1.950
#define MRAD_GA       1.550
#define MRAD_GE       2.720  /* Batsanov et al 20MAR08 toneroma - previous value 2.720 */
#define MRAD_H        1.000
#define MRAD_HG       2.000  /* Batsanov et al 20MAR08 toneroma - previous value none */
#define MRAD_I        2.350  /* Batsanov et al has this value as 2.1 20MAR08 toneroma */
#define MRAD_IN       1.900  /* Batsanov et al 20MAR08 toneroma - previous value none */
#define MRAD_K        2.390  /* Batsanov et al has this value as 3.0 20MAR08 toneromav.*/
#define MRAD_LA       1.830
#define MRAD_LI       1.220  /* Batsanov et al has this value as 2.6 20MAR08 toneroma - previous value 1.220, Amber/Dock5 has this value as 1.17 */
#define MRAD_MG       1.500  /* Batsanov et al has this value as 2.2 20MAR08 toneroma - previous value 1.500, Amber/Dock5 has this value as 1.17 */
#define MRAD_MN       1.190
#define MRAD_MO       1.750
#define MRAD_NSP3     1.850
#define MRAD_NSP4     1.850
#define MRAD_NDEF     1.750
#define MRAD_NA       1.600  /* Batsanov et al has this value as 2.8 20MAR08 toneroma */
#define MRAD_NI       1.240
#define MRAD_OSP3     1.650
#define MRAD_ODEF     1.600
#define MRAD_OS       1.580  /* Batsanov et al 20MAR08 toneroma - previous value 1.580 */
#define MRAD_P        2.100  /* Batsanov et al has this value as 1.7 20MAR08 toneroma */
#define MRAD_PD       1.440
#define MRAD_RB       2.650  /* Batsanov et al 20MAR08 toneroma - previous value 2.650 */
#define MRAD_RE       1.300
#define MRAD_RH       1.220
#define MRAD_RU       1.200
#define MRAD_S        2.000  /* Batsanov et al has this value as 1.8 20MAR08 toneroma */
#define MRAD_SB       1.750  /* Batsanov et al 20MAR08 toneroma - previous value none */
#define MRAD_SE       1.850  /* Batsanov et al 20MAR08 toneroma - previous value none */
#define MRAD_SI       2.200  /* Batsanov et al 20MAR08 toneroma - previous value 2.200 */
#define MRAD_SN       2.300  /* Batsanov et al 20MAR08 toneroma - previous value none */
#define MRAD_SR       2.700  /* Batsanov et al 20MAR08 toneroma - previous value none */
#define MRAD_TE       2.000  /* Batsanov et al 20MAR08 toneroma - previous value none */
#define MRAD_TI       1.950
#define MRAD_TL       1.710  /* Batsanov et al 20MAR08 toneroma - previous value 1.710 */
#define MRAD_U        1.750
#define MRAD_V        1.060
#define MRAD_Y_       1.610
#define MRAD_ZN       1.150  /* Batsanov et al 20MAR08 toneroma - previous value 1.150 */
#define MRAD_ZR       1.420
#define MRAD_UNKNOWN  2.000

#endif
