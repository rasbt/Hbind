#ifndef  _UNBUMP_TRANSLATE_H
#define  _UNBUMP_TRANSLATE_H

extern void  compute_global_vector ( vector_pt vectors,
				     int       number_of_points,
				     vector_t  sum_vector );

extern int  quadratic_formula ( double  a, 
				double  b, 
				double  c,
				double  *sol );

extern int  compute_translation_coefficient ( vector_t v,
					      vector_t w,
					      double   r,
					      double   *factor );

#endif
