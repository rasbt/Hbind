#ifndef _LEAST_SQUARE_FIT_H
#define _LEAST_SQUARE_FIT_H

extern int least_square_fit ( double  a[3][3], 
			      double  b[3][3], 
			      double  c[3][3],
			      double  weights[3],
			      double  matrix[3][3], 
			      double  vector[3] );

#endif
