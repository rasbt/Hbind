#ifndef _FIND_ALL_BUMPS_H
#define _FIND_ALL_BUMPS_H

#include <types.h>

int find_all_bumps(dock_feats_pt features, global_data_pt global);

int atoms_overlap(const atom_pt a, const atom_pt b, const float distance,
		  const float overlap_tolerance, float *overlap);

int atoms_overlap2(const atom_pt a, const atom_pt b, const float sq_dist, 
                   const float overlap_tolerance, float *overlap, float *dist);

/*! At the present this method is not used -- if the handling of water is
 * turned back on, then this will need to be verified for correctness
 */
int water_atom_overlap(const atom_pt a, const float distance, 
                       const float overlap_tolerance, atom_pt w);
#endif
