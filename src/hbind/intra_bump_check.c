#include <stdio.h>
#include <string.h>
#include <math.h>
#include <types.h>
#include <mymalloc.h>
#include <dist_fun.h>
#include "bitstrings.h"
#include "trans_rotate.h"
#include <find_all_bumps.h>
#include <debug_funs.h>

int intra_ligand_bump(molecule_pt ligand, float max_overlap, int *rotated_atoms,
                      int number_of_rotated_atoms )
{
  int      check[MAX_NUMBER_OF_MOL2_ATOMS];
  float    distance,
           orig_distance;
  float sq_dist, tmp;
  int relation;
  int      i, j;
  float overlap;
  int sum;

  atom_pt atoms = ligand->atoms;
  int *relations = ligand->relations;
  int number_of_atoms = ligand->number_of_atoms;
  float *orig_atom_pos = ligand->orig_atom_positions;

  for(i = 0; i < ligand->number_of_atoms; i++) check[i] = RIGID;
  for(i = 0; i < number_of_rotated_atoms; i++) check[rotated_atoms[i]] = FLEX;

  for(i = 0; i < number_of_rotated_atoms; i++)
    for(j = 0; j < number_of_atoms; j++){
      if(rotated_atoms[i] == j || check[j] != RIGID || 
         atoms[j].orbit == ADDED || atoms[rotated_atoms[i]].orbit == ADDED) 
        continue;
	
      if(rotated_atoms[i] < j)
        relation = relations[rotated_atoms[i]*number_of_atoms+j];
      else relation = relations[j*number_of_atoms+rotated_atoms[i]];
      if(relation != DISTANT) continue;

      sq_dist = squared_dist(atoms[j].pos, atoms[rotated_atoms[i]].pos);
      if(sq_dist > DONT_CARE_BUMP_DISTANCE_SQ) continue;

      orig_distance = dist_fun(&orig_atom_pos[3*rotated_atoms[i]], 
                               &orig_atom_pos[3*j]);
 
      sum = atoms[rotated_atoms[i]].act + atoms[j].act;
      tmp = orig_distance - max_overlap;
      if((sq_dist < tmp * tmp || sum == ACCEPTOR_DONOR_HYDROGEN || 
          sum == DONEPTOR_DONOR_HYDROGEN) &&
         atoms_overlap2(&atoms[rotated_atoms[i]], &atoms[j], sq_dist, 
                        max_overlap, &overlap, &distance)){
#ifdef TRACE
        printf("Atom %d (%s %s, %5.3f) bumps with %d (%s %s, %5.3f): %5.3f "
               "(orig = %5.3f, illegal overlap: %5.3f)\n", rotated_atoms[i] + 1,
               atoms[rotated_atoms[i]].name, atoms[rotated_atoms[i]].type_str, 
               atoms[rotated_atoms[i]].rad, j + 1, atoms[j].name, 
               atoms[j].type_str, atoms[j].rad, distance, orig_distance, 
               atoms[rotated_atoms[i]].rad + atoms[j].rad - max_overlap - 
               distance);
#endif
        return BUMP;
      }
    }
  return NO_BUMP;
}

void  identify_second_grade_neighbors ( global_data_pt  global )
{
  int  **neighbors;
  int  *relations,
       *number_of_neighbors;
  int  number_of_atoms,
       neighbor;
  int  i, j, k;

  relations = global->ligand->relations;
  number_of_atoms = global->ligand->number_of_atoms;
  neighbors = global->ligand->neighbors;
  number_of_neighbors = global->ligand->number_of_neighbors;
  for ( i = 0; i < number_of_atoms; i++ )
    for ( j = i + 1; j < number_of_atoms; j++ )
      relations[i*number_of_atoms+j] = DISTANT;
  for ( i = 0; i < number_of_atoms; i++ )
    for ( j = 0; j < number_of_neighbors[i]; j++ )
      {
	neighbor = neighbors[i][j];
	if ( i < neighbor )
	  relations[i*number_of_atoms+neighbor] = CLOSE;
	for ( k = 0; k < number_of_neighbors[neighbor]; k++ )
	  if ( i != neighbors[neighbor][k] 
	       && neighbor < neighbors[neighbor][k] )
	    relations[i*number_of_atoms+neighbors[neighbor][k]] = CLOSE;
      }
}

double  bump_check_neighbors(atom_pt  atoms,
			       int      *atoms_to_check,
			       int      number_of_atoms_to_check,
			       float    allowed_overlap )
{
  atom_pt atom;
  float   radius,
          distance;
  float sq_dist;
  float tmp[2];
  double  total_overlap;
  atom_pt *neighbors;
  int     i, index;
  int     sum;

  total_overlap = 0.0;
  for( index = 0; index < number_of_atoms_to_check; index++ ){
    atom = &atoms[atoms_to_check[index]];
    neighbors = atom->neighbors;
    radius = atom->rad;
    for(i = 0; i < atom->num_nbrs + atom->num_added_nbrs; ++i){
      sq_dist = squared_dist(atom->pos, neighbors[i]->pos);
      if(sq_dist > DONT_CARE_BUMP_DISTANCE_SQ) continue;

      sum = atom->act + neighbors[i]->act;
      tmp[0] = atom->neighbor_dist[i] - allowed_overlap;
      tmp[1] = radius + neighbors[i]->rad - allowed_overlap;
      if(sq_dist < tmp[0]*tmp[0] &&
         ((sum != ACCEPTOR_DONOR && sum != ACCEPTOR_DONEPTOR && 
           sum != DONOR_DONEPTOR && sum != DONEPTOR_DONEPTOR && 
           sq_dist < tmp[1]*tmp[1]) ||
          ((sum == ACCEPTOR_DONOR || sum == DONOR_DONEPTOR || 
            sum == ACCEPTOR_DONEPTOR || sum == DONEPTOR_DONEPTOR) && 
           sq_dist < MIN_HBOND_LENGTH_2))){
        distance = sqrtf(sq_dist);
#ifdef TRACE
        printf("intra BUMP %4s %s %s (%5.3f,%s) with %4s %s %s (%5.3f,%s): "
               "%f (%f), orig = %f \n", atom->name, atom->residue, 
               atom->residue_num, atom->rad, act2str(atom->act), 
               neighbors[i]->name, neighbors[i]->residue, 
               neighbors[i]->residue_num, neighbors[i]->rad, 
               act2str(neighbors[i]->act), distance, tmp[1] - distance, 
               atom->neighbor_dist[i]);
#endif
        total_overlap += tmp[1] - distance;
      }
    }
  }
  return total_overlap;
}
 

int check_for_intra_ligand_bump(global_data_pt global,
                                transform_matrix_pt matrix, int flex_bond, 
                                double angle)
{
  int      atom_indices[MAX_NUMBER_OF_MOL2_ATOMS];
  int      result;
  int      i;
  float *atom_positions = 0;
  atom_pt atoms = global->ligand->atoms;
  int count = 0;

  if(global->ligand->bond_directions[flex_bond] == REVERSE ){
    for(i = 0; i < global->ligand->number_of_atoms; i++ )
        /* this atom is located 'right' from the rotation bond, i.e., at the 
         * side with the larger atom index, which is the side pointing away 
         * from the base */
	if(bitstring_get_bit(atoms[i].fragments, flex_bond))
	  atom_indices[count++] = i;
  }else{
    for(i = 0; i < global->ligand->number_of_atoms; i++ )
      /* this atom is located 'left' from the rotation bond */
      if(bitstring_get_bit ( atoms[i].fragments, flex_bond) == 0)
        atom_indices[count++] = i;
  }      

  atom_positions = (float*) mymalloc(3*count * sizeof(float));
  for(i = 0; i < count; i++){
    memcpy(&atom_positions[3*i], atoms[atom_indices[i]].pos, 3*sizeof(float));
    transform_point(atoms[atom_indices[i]].pos, *matrix );
    rotate_single_atom_around_y_axis(atoms[atom_indices[i]].pos, angle);
    transform_point_back(atoms[atom_indices[i]].pos, *matrix);
  }
  result = intra_ligand_bump(global->ligand, global->intra_overlap, 
                             atom_indices, count);
  /* Restore the ligand atom positions as they where when this function was
   * called */
  for(i = 0; i < count; i++)
    memcpy(atoms[atom_indices[i]].pos, &atom_positions[3*i], 3*sizeof(float));

  if(atom_positions) free(atom_positions);
  atom_positions = 0;
  return result;
}
