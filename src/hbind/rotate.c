#include <math.h>
#include "defs.h"
#include "eigen.h"

int  rotate ( double  a[3][3],  
	      double  b[3][3],  
	      double  weights[3],
	      double  r[3][3] ) 
{
  double u[3][3],
         eigenvector[6][6],
         omega[6][6],	
         h[6][3],
         k[6][3];
  double eigenvalue[6],	
         d[6],
         x[3];
  double sum, 
         det,
         sqrt2;
  int    rank;
  int    i, j, l;
  
  for ( i = 0; i < 3; i++ )
    for ( j = 0; j < 3; j++ ) 
      {
	u[i][j] = 0.0;
	for ( l = 0; l < 3; l++ )
	  u[i][j] += weights[l] * a[l][i] * b[l][j];
      }
  for ( i = 0; i < 3; i++ )
    for ( j = 0; j < 3; j++ )
      {
	omega[i][j] = 0.0;
	omega[i][j+3] = u[i][j];
	omega[i+3][j] = u[j][i];
      }
  for ( i = 3; i < 6; i++ )
    for ( j = 3; j < 6; j++ )
      omega[i][j] = 0.0;
  if ( eigen ( omega, eigenvector, eigenvalue ) == FAILURE )
    return FAILURE;
  sqrt2 = sqrt ( 2.0 );
  for ( i = 0; i < 6; i++ )
    for ( j = 0; j < 3; j++ )
      {
	h[i][j] = eigenvector[j][5-i] * sqrt2;
	k[i][j] = eigenvector[j+3][5-i] * sqrt2;
      }
  for ( i = 0; i < 6; i++ )
    d[i] = eigenvalue[5-i];
  x[0] = h[0][1] * h[1][2] - h[0][2] * h[1][1];
  x[1] = h[0][2] * h[1][0] - h[0][0] * h[1][2];
  x[2] = h[0][0] * h[1][1] - h[0][1] * h[1][0];
  sum = 0.0;
  for ( i = 0; i < 3; i++ ) 
    sum += ( x[i] - h[2][i] ) * ( x[i] - h[2][i] );
  if ( sum > 1.0E-16 ) 
    for ( i = 0; i < 3; i++ ) 
      {
	h[2][i] = - h[2][i];
	k[2][i] = - k[2][i];
      }
  det = u[0][0] * ( u[1][1] * u[2][2] - u[1][2] * u[2][1] )
    + u[0][1] * ( u[2][0] * u[1][2] - u[1][0] * u[2][2] )
    + u[0][2] * ( u[1][0] * u[2][1] - u[2][0] * u[1][1] );
  if ( fabs ( det ) < 1.0E-8 )
    det = 0.0;
  rank = 0;
  for ( i = 0; i < 3; i++ )
    if ( fabs ( d[i] ) > 1.0E-8 ) 
      rank++;
  if ( det < -1.0E-8 ) 
    {
      if ( d[1] != d[2] ) 
	for ( i = 0; i < 3; i++ )
	  k[2][i] = - k[2][i];
      else
	return FAILURE;
    }
  else 
    if ( rank == 2 ) 
      {
	for ( i = 0; i < 3; i++ )
	  h[2][i] = x[i];
	k[2][0] = k[0][1] * k[1][2] - k[0][2] * k[1][1];
	k[2][1] = k[0][2] * k[1][0] - k[0][0] * k[1][2];
	k[2][2] = k[0][0] * k[1][1] - k[0][1] * k[1][0];
      }
    else 
      if ( rank < 2 ) 
	return FAILURE; 
  for ( i = 0; i < 3; i++ )
    for ( j = 0; j < 3; j++ ) 
      {
	r[i][j] = 0.0;
	for ( l = 0; l < 3; l++ )
	  r[i][j] += k[l][i] * h[l][j];
	if ( fabs ( r[i][j] ) < 1.0E-8 ) 
	  r[i][j] = 0.0;
      }
  return SUCCESS;
}

