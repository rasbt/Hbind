#ifndef _DIST_FUN_H
#define _DIST_FUN_H

static inline float
squared_dist(const float *a, const float *b)
{
  return (b[0] - a[0])*(b[0] - a[0]) + (b[1] - a[1])*(b[1] - a[1]) + 
         (b[2] - a[2])*(b[2] - a[2]); 
}

static inline double
squared_dist_double(const double *a, const double *b)
{
  return (b[0] - a[0])*(b[0] - a[0]) + (b[1] - a[1])*(b[1] - a[1]) + 
         (b[2] - a[2])*(b[2] - a[2]); 
}

extern float  dist_fun ( float  pos1[3], float  pos2[3] );

extern double dist_fun_cut ( float  pos1[3], float pos2[3] );

extern double dist_fun_double ( double  pos1[3], double  pos2[3] );

#endif
