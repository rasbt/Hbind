#ifndef _BASICS_H
#define _BASICS_H

#include <stdio.h>

#define MIN_DOUBLE 0.00000001

static inline int compare_double(double d1, double d2)
{
  double my_diff = d1 - d2;
  if(my_diff > MIN_DOUBLE) return 1;
  else if(my_diff < -MIN_DOUBLE) return -1;
  return 0;
}

int hbind_strtod(const char *nptr, double *rv);

FILE *hbind_fopen(const char *path, const char *mode);


#endif

