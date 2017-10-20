#ifndef _CALC_SCORE_FROM_TERMS_
#define _CALC_SCORE_FROM_TERMS_

void score_from_terms (float raw_term_counts[], float affinity_scores[], int max_num_of_score_functions );

enum raw_term_indices								    /* alphabetic symbol*/
  { 	NUMBER_OF_LIGAND_NONH_ATOMS = 0,					/* A*/
	TOTAL_SUM_HPHOB_HPHOB,								/* B*/
	TOTAL_DIFF_HPHOB_HPHOB,								/* C*/
	TOTAL_DIFF_HPHOB_HPHIL,								/* D*/
	TOTAL_SUM_HPHIL_HPHIL,								/* E*/
	TOTAL_DIFF_HPHIL_HPHIL,								/* F*/
	NUMBER_OF_METAL_HBONDS,								/* G*/
	NUMBER_OF_SALT_BRIDGES,								/* H*/
	NUMBER_OF_HBONDS,									/* I*/
	NUMBER_OF_UNSAT_POLAR,								/* J*/
	NUMBER_OF_UNSAT_CHARGE,								/* K*/
	NUMBER_OF_REPULSIVE_POLAR,							/* L*/
	NUMBER_OF_REPULSIVE_CHARGE,							/* M*/
	RATIO_OF_INTERFACIAL_LIGAND_HEAVY_ATOMS,			/* N*/
	NUMBER_OF_INTERFACIAL_LIGAND_ATOMS,					/* O*/
	NUMBER_OF_EXPOSED_HYDRO_LIGAND_ATOMS,				/* P*/
	NUMBER_OF_LIGAND_FLEXIBLE_BONDS,					/* Q*/
	NUMBER_OF_ALL_INTERFACIAL_FLEXIBLE_BONDS,			/* R*/
	NUMBER_OF_HPHOB_HPHOB_CONTACT,						/* S*/
	TOTAL_AVG_HPHOB_HPHOB,								/* T*/
	OLD_HBIND_TOTAL_SUM_HPHOB_HPHOB,					/* U*/
	NUMBER_OF_FLEXIBLE_INTERFACIAL_LIGAND_BONDS,		/* V*/
	NUMBER_OF_HPHOB_HPHIL_CONTACT,						/* W*/
	TOTAL_HPHOB_COMP,                                                      /*X*/
	TOTAL_HPHOB_HPHOB_COMP,                                                      /*Y*/
	CONTACT_HPHOB_HPHOB,                                                      /*Z*/
	CONTACT_HPHIL_HPHIL,                                                      /*AA*/
	CONTACT_HPHOB_HPHIL,                                                      /*AB*/
	TOTAL_TARGET_HPHOB_HPHOB,                                                      /*AC*/
	INCREASE_HPHOB_ENVIRON,                                                      /*AD*/
	TOTAL_TARGET_HYDRO,                                                      /*AE*/
	NORM_TARGET_HPHOB_HPHOB,                                                      /*AF*/
	TOTAL_TARGET_HPHOB_CONTACT,                                                      /*AG*/
	NORM_TARGET_HPHOB_CONTACT,                                                      /*AH*/
	TOTAL_LIGAND_HYDRO,                                                      /*AI*/
	HBOND_RATIO,
	NORM_POLAR,
	NORM_UNSAT,
	NORM_HPHOB_COMP,						/* N_HPHOB_COMP*/
	NUMBER_OF_INTRA_TARGET_HBONDS,						/* X*/
	NUMBER_OF_INTRA_TARGET_SALT_BRIDGES,				/* Y*/
	INTERMEDIATE_OVERLAP,								/* Z*/
	TOTAL_OVERLAP,										/* [*/
	NUMBER_OF_UNSAT_LIGAND_POLAR,						/* \ */
	NUMBER_OF_UNSAT_TARGET_POLAR,						/* ]*/
	NUMBER_OF_UNSAT_LIGAND_CHARGE,						/* ^*/
	NUMBER_OF_UNSAT_TARGET_CHARGE,						/* _*/
	NORM_HPHOB_COMP_2,						/* ELECTRO_ENVIRON*/
	TOTAL_AVG_HPHOB_HPHOB_NEW,								/* T_NEW*/
	TOTAL_BURIED_HPHOB,								/* B_HYDRO	*/
	TOTAL_EXPOSED_HPHOB,								/* E_HYDRO	*/
	NORM_HPHOB_COMP_3,								/* TE_HYDRO	*/

  /* Repulsive interaction is an atom pair of P-L. 
   * So following numbers are halves of NUMBER_OF_REPULSIVE_POLAR,
	NUMBER_OF_REPULSIVE_CHARGE resp.
	NUMBER_OF_REPULSIVE_LIGAND_POLAR,					 LL
	NUMBER_OF_REPULSIVE_TARGET_POLAR,					 LT
	NUMBER_OF_REPULSIVE_LIGAND_CHARGE,					 ML
	NUMBER_OF_REPULSIVE_TARGET_CHARGE					 MT
	*/
		
};

#endif
