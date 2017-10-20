
#ifndef DISTANCE_ARRAY_HEADER_INCLUDED
#define DISTANCE_ARRAY_HEADER_INCLUDED

#include <atom.h>

typedef struct{
  atom_pt *atoms;
  int num_atoms;
  int size;
}*dist_bin_pt, dist_bin_t;
  
typedef struct{
  dist_bin_pt bins;
  float box_width;
  float min_corner[3];
  float max_corner[3];
  int num_bins[3];
}*dist_array_pt, dist_array_t;

void
distance_array(dist_array_pt darray, const atom_pt atoms, int num_atoms,
               const float box_width, const float max_dist);

void
free_distance_array(dist_array_pt d);

dist_bin_pt 
get_bin(dist_array_pt darray, float *pt);

#endif
