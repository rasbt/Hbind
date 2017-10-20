#include "sf_weights.h"
#include "calc_score_from_terms.h"

/* OrientScore: Scoring function 26 from Affinity training-testing results */
float aff_26_score (float raw_term_counts[] )
{
#define SF 25
  return ( affinity_wts[SF][REGRESSION_CONST]
	   /*C*/    + affinity_wts[SF][TERM_1] * (raw_term_counts[ TOTAL_DIFF_HPHOB_HPHOB] )
	   /*D*/    + affinity_wts[SF][TERM_2] * (raw_term_counts[ TOTAL_DIFF_HPHOB_HPHIL] )
	   );
  
#undef SF
}


/* AffiScore: Scoring function 62 from Affinity training-testing results */
float aff_62_score (float raw_term_counts[] )
{

#define SF 61
     return (    affinity_wts[SF][REGRESSION_CONST]
    /*T*/       + affinity_wts[SF][TERM_1] * raw_term_counts[ TOTAL_AVG_HPHOB_HPHOB ]
    /*G+H+I*/       + affinity_wts[SF][TERM_2] * ( raw_term_counts[ NUMBER_OF_METAL_HBONDS ]
                                                        + raw_term_counts[ NUMBER_OF_SALT_BRIDGES ]
                                                        + raw_term_counts[ NUMBER_OF_HBONDS ])
    /*J+K*/       + affinity_wts[SF][TERM_3] * (  raw_term_counts[ NUMBER_OF_UNSAT_POLAR ] 
                                           + raw_term_counts[ NUMBER_OF_UNSAT_CHARGE ] )
           );
#undef SF
}


/* AffiScore: Scoring function 62 from Affinity training-testing results */
float aff_64_score (float raw_term_counts[] )
{

#define SF 62
     return (    affinity_wts[SF][REGRESSION_CONST]
    /*T_new*/       + affinity_wts[SF][TERM_1] * raw_term_counts[ TOTAL_HPHOB_HPHOB_COMP ]
    /*G+H+I*/       + affinity_wts[SF][TERM_2] * ( raw_term_counts[ NUMBER_OF_METAL_HBONDS ]
                                                        + raw_term_counts[ NUMBER_OF_SALT_BRIDGES ]
                                                        + raw_term_counts[ NUMBER_OF_HBONDS ])
    /*J+K*/       + affinity_wts[SF][TERM_3] * (  raw_term_counts[ NUMBER_OF_UNSAT_POLAR ] 
                                           + raw_term_counts[ NUMBER_OF_UNSAT_CHARGE ] )
           );
#undef SF
}

/* AffiScore: Scoring function 62 from Affinity training-testing results */
float orient_99_score (float raw_term_counts[] )
{

#define SF 63
     return (    affinity_wts[SF][REGRESSION_CONST]
    /*AE*/       + affinity_wts[SF][TERM_1] * raw_term_counts[ TOTAL_TARGET_HYDRO ]
    /*Z*/       + affinity_wts[SF][TERM_2] * raw_term_counts[ CONTACT_HPHOB_HPHOB ]		 
    /*G+H+I*/       + affinity_wts[SF][TERM_3] * ( raw_term_counts[ NUMBER_OF_METAL_HBONDS ]
                                                        + raw_term_counts[ NUMBER_OF_SALT_BRIDGES ]
                                                        + raw_term_counts[ NUMBER_OF_HBONDS ])
           );
#undef SF
}

/* AffiScore: Scoring function 62 from Affinity training-testing results */
float aff_100_score (float raw_term_counts[] )
{

#define SF 64
     return (    affinity_wts[SF][REGRESSION_CONST]
    /*AE*/       + affinity_wts[SF][TERM_1] * raw_term_counts[ TOTAL_TARGET_HYDRO ]
    /*Z*/       + affinity_wts[SF][TERM_2] * raw_term_counts[ CONTACT_HPHOB_HPHOB ]		 
    /*G+H+I*/       + affinity_wts[SF][TERM_3] * ( raw_term_counts[ NUMBER_OF_METAL_HBONDS ]
                                                        + raw_term_counts[ NUMBER_OF_SALT_BRIDGES ]
                                                        + raw_term_counts[ NUMBER_OF_HBONDS ])
    /*J+K*/       + affinity_wts[SF][TERM_4] * (  raw_term_counts[ NUMBER_OF_UNSAT_POLAR ] 
                                           + raw_term_counts[ NUMBER_OF_UNSAT_CHARGE ] )
           );
#undef SF
}

/* AffiScore: Scoring function 62 from Affinity training-testing results */
float combo_99_score (float raw_term_counts[] )
{

#define SF 65
     return (    affinity_wts[SF][REGRESSION_CONST]
    /*AE*/       + affinity_wts[SF][TERM_1] * raw_term_counts[ TOTAL_TARGET_HYDRO ]
    /*Z*/       + affinity_wts[SF][TERM_2] * raw_term_counts[ CONTACT_HPHOB_HPHOB ]
    /*G+H+I*/       + affinity_wts[SF][TERM_3] * ( raw_term_counts[ NUMBER_OF_METAL_HBONDS ]
                                                        + raw_term_counts[ NUMBER_OF_SALT_BRIDGES ]
                                                        + raw_term_counts[ NUMBER_OF_HBONDS ])
           );
#undef SF
}

/* AffiScore: Scoring function 62 from Affinity training-testing results */
float aff_99_score (float raw_term_counts[] )
{

#define SF 66
     return (    affinity_wts[SF][REGRESSION_CONST]
    /*AE*/       + affinity_wts[SF][TERM_1] * raw_term_counts[ TOTAL_TARGET_HYDRO ]
    /*Z*/       + affinity_wts[SF][TERM_2] * raw_term_counts[ CONTACT_HPHOB_HPHOB ]
    /*G+H+I*/       + affinity_wts[SF][TERM_3] * ( raw_term_counts[ NUMBER_OF_METAL_HBONDS ]
                                                        + raw_term_counts[ NUMBER_OF_SALT_BRIDGES ]
                                                        + raw_term_counts[ NUMBER_OF_HBONDS ])
           );
#undef SF
}

/* AffiScore: Scoring function 116 from Conformation/Orientation training-testing results */
float conf_116_score (float raw_term_counts[] )
{

#define SF 67
  return (    affinity_wts[SF][REGRESSION_CONST]
	      /*AE+AI*/       + affinity_wts[SF][TERM_1] * (raw_term_counts[ TOTAL_TARGET_HYDRO ] 
							    + raw_term_counts[TOTAL_LIGAND_HYDRO])
	      /*Z*/       + affinity_wts[SF][TERM_2] * raw_term_counts[ CONTACT_HPHOB_HPHOB ]
	      /*G+H+I*/       + affinity_wts[SF][TERM_3] * ( raw_term_counts[ NUMBER_OF_METAL_HBONDS ]
							     + raw_term_counts[ NUMBER_OF_SALT_BRIDGES ]
							     + raw_term_counts[ NUMBER_OF_HBONDS ])
	      /*J+K*/       + affinity_wts[SF][TERM_4] * ( raw_term_counts[ NUMBER_OF_UNSAT_POLAR ]
							   + raw_term_counts[ NUMBER_OF_UNSAT_CHARGE ])
	      );
#undef SF
}

void score_from_terms (float raw_term_counts[], float affinity_scores[], int max_num_of_score_functions )
{
  affinity_scores[62] = aff_64_score(raw_term_counts);
  affinity_scores[63] = orient_99_score(raw_term_counts);
}
