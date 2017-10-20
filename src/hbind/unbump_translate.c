#include <math.h>
#include <stdio.h>
#include "types.h"
#include "defs.h"

/*
 *  This function computes a global transformation vector that is 
 *  simply the sum of the different translation vectors needed to get
 *  rid of the single bumps. By taking the sum of the vectors, some kind
 *  of a weighted global vector is computed, since the influence of
 *  more intense bumps is larger than that of small bumps.
 */
void  compute_global_vector ( vector_pt vectors,
			      int       number_of_points,
			      vector_t  sum_vector )
{
  int       i, j;
  
  for ( i = 0; i < 3; i++ )
    sum_vector[i] = 0;
  for ( i = 0; i < number_of_points; i++ )
    for ( j = 0; j < 3; j++ )
      sum_vector[j] += vectors[i][j];
}

/*
 *  Simple function that returns the positive solution of the
 *  quadratic formula  a * x^2 + b * x + c = 0
 */ 
int  quadratic_formula ( double  a, 
			 double  b, 
			 double  c,
			 double  *sol )
{
  double  discr;
  double  sol1, sol2;

  discr = b * b - 4.0 * a * c;
  if ( discr < 0 )
    return FAILURE;
  sol1 = ( (-1.0) * b - sqrt ( discr ) ) / ( 2.0 * a );
  sol2 = ( (-1.0) * b + sqrt ( discr ) ) / ( 2.0 * a );
  if ( sol1 >= 0 )
    *sol = sol1;
  else 
    *sol = sol2;
  return SUCCESS;
}

/*
 *  This function computes the coefficient `alpha' for the
 *  translation vector w to move an atom out of the
 *  bumping zone of another atom. Tails of vectors v and
 *  u are at center of the fixed atom, head of v and tail of
 *  w are at the center of the moving atom, and heads of vectors
 *  w and u are at the target position for the moving atom
 *  at the nearest bump-free point when translating it
 *  along vector w. Length of vector u is r, which is the 
 *  sum of the van der Waals radii of both atoms minus the allowed 
 *  overlap. This function solves the following set of four equations
 *  by the use of the quadratic formula:
 *  (I)-(III)  (v) + alpha * (w) = (u)
 *  (IV)       r = sqrt ( u_x * u_x + u_y * u_y + u_z * u_z )
 *  if there is no solution to the quadratic formula, FAILURE is returned
 */
int  compute_translation_coefficient ( vector_t v,
				       vector_t w,
				       double   r,
				       double   *factor )
{
  double  a, b, c;
  int     result;
  
  a = w[X] * w[X] + w[Y] * w[Y] + w[Z] * w[Z];
  b = 2 * w[X] * v[X] + 2 * w[Y] * v[Y] + 2 * w[Z] * v[Z];
  c = v[X] * v[X] + v[Y] * v[Y] + v[Z] * v[Z] - r * r;

  result = quadratic_formula ( a, b, c, factor );
  return result;
}
