#include "defs.h"
#include "rotate.h"

int least_square_fit ( double  a[3][3], 
		       double  b[3][3], 
		       double  c[3][3],
		       double  weights[3],
		       double  matrix[3][3], 
		       double  vector[3] )
{
  double  newa[3][3],	
          newb[3][3];	
  double  ta[3],
          tb[3],	
          help[3];
  double  weightsum;
  int     i, j;

  for( i = 0; i < 3; i++ )
    ta[i] = tb[i] = 0.0;
  weightsum = 0.0;
  for ( i = 0; i < 3; i++ ) 
    {
      weightsum += weights[i];
      for ( j = 0; j < 3; j++ )
	{
	  ta[j] += weights[i] * a[i][j];
	  tb[j] += weights[i] * b[i][j];
	}
    }
  for ( i = 0; i < 3; i++ )
    {
      ta[i] /= weightsum;
      tb[i] /= weightsum;
    }
  for ( i = 0; i < 3; i++ )
    for ( j = 0; j < 3; j++ )
      {
	newa[i][j] = a[i][j] - ta[j];
	newb[i][j] = b[i][j] - tb[j];
      }
  if ( rotate ( newa, 
		newb, 
		weights,
		matrix ) == FAILURE )   
    return FAILURE;
  for ( i = 0; i < 3; i++ ) 
    {
      for ( j = 0; j < 3; j++ )
	c[i][j] = a[i][j] - ta[j];
      for ( j = 0; j < 3; j++ )
	{
	  help[j] = matrix[j][0] * c[i][0];
	  help[j] += matrix[j][1] * c[i][1];
	  help[j] += matrix[j][2] * c[i][2];
	}
      for ( j = 0; j < 3; j++ )
	c[i][j] = help[j] + tb[j];
      vector[i] = tb[i] - ( matrix[i][0] * ta[0] 
		          + matrix[i][1] * ta[1] 
			  + matrix[i][2] * ta[2] );
    }
  return SUCCESS;
}
