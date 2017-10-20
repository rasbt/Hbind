#include <math.h>
#include "defs.h"
#include "basics.h"

void  tred2 ( double  a[6][6],
	      double  z[6][6],
	      double  d[6],
	      double  e[6],
	      double  tolerance )
{
  double  f, g, h, m;
  int     l;
  int     i, j, k;
  
  for ( i = 0; i < 6; i++ )
    for ( j = 0; j < 6; j++ )
      z[i][j] = a[i][j];
  for ( i = 5; i >= 1; i-- ) 
    {
      l = i - 2;
      f = z[i][i-1];
      g = 0.0;
      for ( k = 0; k <= l; k++ ) 
	g += z[i][k] * z[i][k];
      h = g + f * f;	    
      if ( g <= tolerance ) 
	{
	  e[i] = f;
	  h = 0.0;
	}
      else 
	{
	  l++;
	  if ( f >= 0.0 )	
	    g = (-1.0) * sqrt ( h );
	  else		
	    g = sqrt ( h );
	  e[i] = g;
	  h -= f * g;
	  z[i][i-1] = f - g;
	  f = 0.0;
	  for ( j = 0; j <= l; j++ ) 
	    {
	      z[j][i] = z[i][j] / h;
	      g = 0.0;
	      for ( k = 0; k <= j; k++ )
		g += z[j][k] * z[i][k];
	      for ( k = j + 1; k <= l; k++ )
		g += z[k][j] * z[i][k];
	      e[j] = g / h;
	      f += g * z[j][i];
	    }
	  m  = f / ( h + h );
	  for ( j = 0; j <= l; j++ ) 
	    {
	      f = z[i][j];
	      g = e[j] - m * f;
	      e[j] = g;
	      for ( k = 0; k <= j; k++ )
		z[j][k] = z[j][k] 
		  - ( f * e[k] + g * z[i][k] );
	    }
	}
      d[i] = h;
    }
  d[0] = 0.0;
  e[0] = 0.0;
  for ( i = 0; i < 6; i++ ) 
    {
      l = i - 1;
      if ( d[i] != 0.0 )
	for ( j = 0; j <= l; j++ ) 
	  {
	    g = 0.0;
	    for ( k = 0; k <= l; k++ ) 
	      g += ( z[i][k] * z[k][j] );
	    for ( k = 0; k <= l; k++ ) 
	      z[k][j] -= g * z[k][i];
	  }
      d[i] = z[i][i];
      z[i][i] = 1.0;
      for ( j = 0; j <= l; j++ ) 
	{
	  z[i][j] = 0.0;
	  z[j][i] = 0.0;
	}
    }
}

int  tql2 ( double   z[6][6],
	    double   d[6],
	    double   e[6],
	    double   eps )
{
  double h, p, r, b, c, f, g, s;
  int    L, k, m;
  int    i, j;

  for ( i = 1; i < 6; i++ )
    e[i-1] = e[i];
  e[5] = 0.0;
  b = 0.0;
  f = 0.0;
  for ( L = 0; L < 6; L++ ) 
    {
      j = 0;
      h = eps * ( fabs ( d[L] ) + fabs ( e[L] ) );
      if ( b < h ) 
	b = h;
      for ( m = L; compare_double( fabs ( e[m] ), b ) == 1 && m < 6; m++ ) ;
      if ( m != L )
	do 
	  {
	    if ( j >= 30 ) 
	      return FAILURE; 
	    j++;
	    p = ( d[L+1] - d[L] ) / ( 2.0 * e[L] );
	    r = sqrt ( p * p + 1.0 );
	    if ( p < 0.0 )
	      h = p - r;
	    else
	      h = p + r;
	    h = d[L] - e[L] / h;
	    for ( i = L; i < 6; i++ )
	      d[i] -= h;
	    f += h;
	    p = d[m];
	    c = 1.0;
	    s = 0.0;
	    for ( i = m - 1; i >= L; i-- ) 
	      {
		g = c * e[i];
		h = c * p;
		if ( fabs ( p ) >= fabs ( e[i] ) ) 
		  {
		    c = e[i] / p;
		    r = sqrt ( c * c + 1.0 );
		    e[i+1] = s * p * r;
		    s = c / r;
		    c = 1.0 / r;
		  }
		else 
		  {
		    c = p / e[i];
		    r = sqrt ( c * c + 1.0 );
		    e[i+1] = s * e[i] * r;
		    s = 1.0 / r;
		    c *= s;
		  }
		p = c * d[i] - s * g;
		d[i+1] = h + s * ( c * g + s * d[i] );
		for ( k = 0; k < 6; k++ ) 
		  {
		    h = z[k][i+1];
		    z[k][i+1] = s * z[k][i] + c * h;
		    z[k][i] = c * z[k][i] - s * h;
		  }
	      }
	    e[L] = s * p;
	    d[L] = c * p;
	  } 
	while ( compare_double( fabs ( e[L] ), b ) == 1 );
      d[L] += f;
    }
  for ( i = 0; i < 6; i++ ) 
    {
      k = i;
      p = d[i];
      for ( j = i + 1; j < 6; j++ )
	if ( d[j] < p ) 
	  {
	    k = j;
	    p = d[j];
	  }
      if ( k != i ) 
	{
	  d[k] = d[i];
	  d[i] = p;
	  for ( j = 0; j < 6; j++ ) 
	    {
	      p = z[j][i];
	      z[j][i] = z[j][k];
	      z[j][k] = p;
	    }
	}
    }
  return SUCCESS;
}
 
int  eigen ( double a[6][6],
	     double b[6][6],
	     double c[6] )	

{
  double e[6];

  tred2 ( a, b, c, e, 1.0 / HUGE );
  if ( tql2 ( b, c, e, 1.0E-20 ) == FAILURE )
    return FAILURE;
  return SUCCESS;
}

