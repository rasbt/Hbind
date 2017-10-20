#include <math.h>
#include <stdio.h>
#include "types.h"
#include "defs.h"
#include "dist_fun.h"
#include "unbump_water.h"
#include <check_complementarity.h>
#include <find_all_bumps.h>


/* Seems to tally total_overlap */
int find_all_bumps(dock_feats_pt features, global_data_pt global)
{
  atom_pt atom;
  atom_pt water;
  atom_pt *neighbors;
  float distance;
  float sq_dist;
  float overlap_tolerance = global->side_chain_overlap;
  float overlap;
  float tmp;
  int number_of_bumps = 0,
           sum,
           dist_int;
  int      i, j, k;

  atom_pt target = global->target_atoms;
  atom_pt ligand = global->ligand->atoms;
  int *target_bumps = global->target_bumps;
  int *ligand_bumps = global->ligand_bumps;
  const int num_targ_atoms = global->number_of_target_atoms;
  features->total_overlap = 0.0;

#ifdef TRACE
  printf("\n");
#endif

  for ( i = 0; i < num_targ_atoms; i++ ){
    for ( j = 0; j < global->ligand->number_of_atoms; j++ ){
      /* Ligand hydrogen atoms are considered if there were not added by HBIND 
       * */
      if(ligand[j].orbit == ADDED) continue;

      sq_dist = global->target_ligand_sq_dists[j*num_targ_atoms + i];
      if(sq_dist > DONT_CARE_BUMP_DISTANCE_SQ) continue;

      if(atoms_overlap2(&target[i], &ligand[j], sq_dist, overlap_tolerance, 
                        &overlap, &distance)){
#ifdef TRACE
        printf("sum=%d, distance=%f, ligand[j].rad=%f, target[i].rad=%f, "
               "tolerated overlap=%f, remaining overlap=%f\n", 
               target[i].act + ligand[j].act, distance, ligand[j].rad, 
               target[i].rad, overlap_tolerance, overlap);
        printf("cumulative overlap = %f\n", features->total_overlap);
        printf("BUMP no. %2d: %3s %s %3s (%7.3f %7.3f %7.3f) <--> %2d %4s %2d "
               "(%7.3f %7.3f %7.3f) : %7.3f (%7.3f)\n", number_of_bumps,
               target[i].name, target[i].residue, target[i].residue_num,
               target[i].pos[X], target[i].pos[Y], target[i].pos[Z],
               j + 1, ligand[j].type_str, ligand[j].fragment,
               ligand[j].pos[X], ligand[j].pos[Y], ligand[j].pos[Z],
               distance, overlap);
#endif
        /* the vector `bump_point' is already full, thus there is no hope for 
         * this ligand ==> bump-check failed */
        if ( number_of_bumps >= MAX_SIDE_CHAIN_BUMPS ) {
          global->number_of_bumps = number_of_bumps;
          return FAILURE;
        }
        target_bumps[number_of_bumps] = i;
        ligand_bumps[number_of_bumps] = j;
        features->total_overlap += overlap;
        number_of_bumps++;
      }
    }
  }
#ifdef TRACE
  printf("cumulative overlap = %f\n", features->total_overlap);
#endif


  overlap_tolerance = global->intra_overlap;
  for(i = 0; i < global->number_of_target_residues; i++){
    if(global->target_intra_overlap[i] != YES) continue;
     
    for(j = 0; j < global->target_residues[i].number_of_atoms; j++){
      atom = &target[global->target_residues[i].start_atom+j];
      neighbors = atom->neighbors;

      for(k = 0; k < atom->num_nbrs + atom->num_added_nbrs; ++k){

        sq_dist = squared_dist(atom->pos, neighbors[k]->pos);
        if(sq_dist > DONT_CARE_BUMP_DISTANCE_SQ) continue;

        tmp = atom->neighbor_dist[k] - overlap_tolerance;
        /* Added the first inequality to keep this consistent with previous
         * hbind code */
        if(sq_dist < tmp * tmp &&
           atoms_overlap2(atom, neighbors[k], sq_dist, overlap_tolerance, 
                          &overlap, &distance)){
          /* the vector `bump_point' is already full, thus there is no hope for 
           * this ligand ==> bump-check failed */
          if ( number_of_bumps >= MAX_SIDE_CHAIN_BUMPS ) {
            global->number_of_bumps = number_of_bumps;
            return FAILURE;
          }
          features->total_overlap += overlap;
          target_bumps[number_of_bumps] = 
            global->target_residues[i].start_atom+j;
          ligand_bumps[number_of_bumps] = (-1) * ((neighbors[k] - target) + 1);
          number_of_bumps++;
#ifdef TRACE
          printf("intra BUMP %s %s %s (%7.3f,%d) with %s %s %s (%7.3f,%d): %f "
                 "(%f), orig = %f \n", atom->name, atom->residue,
                 atom->residue_num, atom->rad, atom->act, neighbors[k]->name, 
                 neighbors[k]->residue, neighbors[k]->residue_num, 
                 neighbors[k]->rad, neighbors[k]->act, distance,
                 atom->rad + neighbors[k]->rad - overlap_tolerance - distance,
                 atom->neighbor_dist[k]);
          printf("cumulative overlap = %f\n", features->total_overlap);
#endif
        }
      }
    }
  }

  global->number_of_bumps = number_of_bumps;
  if(number_of_bumps == 0) return NO_BUMP;
  return BUMP;
}

int 
atoms_overlap(const atom_pt a, const atom_pt b, const float distance, 
              const float overlap_tolerance, float *overlap)
{
  *overlap = 0.0;
  if(is_hbond_interaction(a->act, b->act)){ 
    if(distance < MIN_HBOND_LENGTH) *overlap = MIN_HBOND_LENGTH - distance;
  }else if((a->act + b->act == ACCEPTOR_DONOR_HYDROGEN ||
          a->act + b->act == DONEPTOR_DONOR_HYDROGEN)){
    if(distance < MIN_HBOND_LENGTH_HYDROGEN) 
      *overlap = MIN_HBOND_LENGTH_HYDROGEN - distance;
  }else if(a->act == METAL_1 || b->act == METAL_1){
    if(distance < MIN_METAL_1_HBOND_LENGTH)
      *overlap = MIN_METAL_1_HBOND_LENGTH - distance;
  }else if(a->act == METAL_2 || b->act == METAL_2){
    if(distance < MIN_METAL_2_HBOND_LENGTH)
      *overlap = MIN_METAL_2_HBOND_LENGTH - distance;
  }else if(distance < a->rad + b->rad - overlap_tolerance) 
    *overlap = a->rad + b->rad - overlap_tolerance - distance;

  if(*overlap > 0.0) return 1;
  return 0;
} 

int
atoms_overlap2(const atom_pt a, const atom_pt b, const float sq_dist, 
               const float overlap_tolerance, float *overlap, float *dist)
{
  float tmp;
  *overlap = 0.0;

  if(is_hbond_interaction(a->act, b->act)){ 
    if(sq_dist < MIN_HBOND_LENGTH_2){
      *dist = sqrt(sq_dist);
      *overlap = MIN_HBOND_LENGTH - *dist;
    }
  }else if((a->act + b->act == ACCEPTOR_DONOR_HYDROGEN ||
          a->act + b->act == DONEPTOR_DONOR_HYDROGEN)){
    if(sq_dist < MIN_HBOND_LENGTH_HYDROGEN_2){
      *dist = sqrt(sq_dist);
      *overlap = MIN_HBOND_LENGTH_HYDROGEN - *dist;
    }
  }else if(a->act == METAL_1 || b->act == METAL_1){
    if(sq_dist < MIN_METAL_1_HBOND_LENGTH_2){
      *dist = sqrt(sq_dist);
      *overlap = MIN_METAL_1_HBOND_LENGTH - *dist;
    }
  }else if(a->act == METAL_2 || b->act == METAL_2){
    if(sq_dist < MIN_METAL_2_HBOND_LENGTH_2){
      *dist = sqrt(sq_dist);
      *overlap = MIN_METAL_2_HBOND_LENGTH - *dist;
    }
  }else{
    tmp = a->rad + b->rad - overlap_tolerance;
    if(sq_dist < tmp*tmp){
      *dist = sqrt(sq_dist);
      *overlap = tmp - *dist;
    }
  }

  if(*overlap > 0.0) return 1;
  return 0;
}

int 
water_atom_overlap(const atom_pt a, const float distance, 
                   const float overlap_tolerance, atom_pt w)
{
  /* there is a bump of this water with a polar ligand atom, so
 can be displaced at not cost */
  if((a->act == ACCEPTOR || a->act == DONEPTOR || a->act == DONOR) && 
     distance < MIN_HBOND_LENGTH)
    w->state = POLAR_DISPLACED;
  /* the ligand atom is pretty close to the water, so we just
     assume that the water will be displaced */
  else if(a->act == NOTHING && 
          distance < WATER_RAD / 2.0 + a->rad - overlap_tolerance)
    w->state = DISPLACED;
  /* there is only a relatively small overlap, so rather than assuming we
   * can displace the water, we attempt to "unbump" it */
  else if(a->act == NOTHING && 
          distance < WATER_RAD + a->rad - overlap_tolerance)
    w->state = WATER_OVERLAP_NEEDS_TO_BE_RESOLVED;

  return SUCCESS;
}

