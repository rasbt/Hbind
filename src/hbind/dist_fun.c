#include <math.h>

float  dist_fun ( float  pos1[3], 
		  float  pos2[3] )
{
  double dx, dy, dz, dist;

  dx = (double) pos1[0] - pos2[0];
  dy = (double) pos1[1] - pos2[1];
  dz = (double) pos1[2] - pos2[2];

  dist = (double) sqrt ( (double) dx * dx + (double) dy * dy + (double) dz * dz );
  return (float) dist;
  
}


double dist_fun_cut ( float  pos1[3], 
	              float  pos2[3] )
{
  double dx, dy, dz;
  int    temp;
  double dist;

  dx = (double) ( pos1[0] - pos2[0] );
  dy = (double) ( pos1[1] - pos2[1] );
  dz = (double) ( pos1[2] - pos2[2] );

  dist = sqrt ( dx * dx + dy * dy + dz * dz );

  
  temp = (int) ( (dist + 0.00005) * 10000 );
  return ( (float) temp / 10000 );
  
}


double dist_fun_double ( double  pos1[3],
			 double  pos2[3] )
{
   double dx, dy, dz;

   dx = pos1[0] - pos2[0];
   dy = pos1[1] - pos2[1];
   dz = pos1[2] - pos2[2];

   return sqrt ( dx * dx + dy * dy + dz * dz );
}


