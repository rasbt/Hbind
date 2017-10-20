#ifndef DISTANCE_MATRICES_HEADER_FILE_INCLUDED
#define DISTANCE_MATRICES_HEADER_FILE_INCLUDED

#include <octree_array.h>

static inline int diff_int(int  a, int  b)
{
  return (a > b ? a - b : b - a);
}

void
init_target_nbr_arrays(atom_pt target_atoms, const int num_atoms);

void
initialize_inter_dist_matrix(const float* targ_positions, 
                             const size_t num_targ_atoms, 
                             const float* lig_positions, 
                             const size_t num_lig_atoms,
                             float *M, const float max_dist);

#endif
