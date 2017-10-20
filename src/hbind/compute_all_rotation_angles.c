#include <math.h>
#include <float.h>
#include <string.h>
#include "types.h"
#include "rotatable_bonds.h"
#include "trans_rotate.h"
#include "unbump_side_chains.h"
#include "compute_ligand_angles.h"
#include "compute_target_angles.h"

double compute_rotation_angle(atom_pt fixed_atom, atom_pt rotated_atom, 
                              float max_overlap, 
                              transform_matrix_t unbump_matrix)
{
  double angle;
  double orig_angle;
  float rot_atom_orig_pos[3];
  float fix_atom_orig_pos[3];
  double lxz;
  double l;
  double d =  fixed_atom->rad + rotated_atom->rad - max_overlap + 0.0001;

  memcpy(rot_atom_orig_pos, rotated_atom->pos, 3*sizeof(*rot_atom_orig_pos));
  memcpy(fix_atom_orig_pos, fixed_atom->pos, 3*sizeof(*fix_atom_orig_pos));
  
  /* first transform the two points so that a computation of a rotation
     angle is possible */
  transform_point(rotated_atom->pos, unbump_matrix);  
  transform_point(fixed_atom->pos, unbump_matrix);
  lxz = sqrt(((double) rotated_atom->pos[X] *
              (double) rotated_atom->pos[X]) + 
             ((double) rotated_atom->pos[Z] * 
              (double) rotated_atom->pos[Z]));
  l = sqrt(((double) rotated_atom->pos[X] * 
            (double) rotated_atom->pos[X]) + 
           ((double) rotated_atom->pos[Y] * 
            (double) rotated_atom->pos[Y]) + 
           ((double) rotated_atom->pos[Z] * 
            (double) rotated_atom->pos[Z]));

  angle = compute_unbump_angle(fixed_atom->pos[X], d, rotated_atom->pos[Y], 
                               fixed_atom->pos[Y], lxz, l);
  /*printf("compute_all_angle:  angle = %f ",angle);*/
  if(angle == FLT_MAX){
    /* Copy the original positions rather than "move" back */ 
    memcpy(rotated_atom->pos, rot_atom_orig_pos, 3*sizeof(*rot_atom_orig_pos));
    memcpy(fixed_atom->pos, fix_atom_orig_pos, 3*sizeof(*fix_atom_orig_pos));
    return angle;
  }
  /* a negative angle is used to indicate if the final position of
     the side-chain atom is in the negative x-range, in this case
     the angle relative to the (positive) x-axis is 180 - angle */

  if(angle < 0) angle = M_PI + angle;

  if ((rotated_atom->pos[X] * rotated_atom->pos[X]) < SMALL_FLOAT){
    if (rotated_atom->pos[Z] < 0 ) orig_angle = 0.0 - M_PI/2.0;
    else orig_angle = M_PI/2.0;
  }else
    orig_angle = asin ( (double) rotated_atom->pos[Z] / 
		      sqrt ( (double ) 
			     ( ( rotated_atom->pos[X]
				 * rotated_atom->pos[X] )
			       + ( rotated_atom->pos[Z]
				   * rotated_atom->pos[Z] ) ) ) );
  /* if the two angles are too close to each other, we consider them to be
     the same */
  if(fabs(fabs(orig_angle) - fabs(angle)) < 0.001) angle = 0.0;
  else{
	  /* depending on the location of the rotate atom relative to the X axis
	     (pos. vs. neg. Z value) and on the choice of the smallest or the
	     larger rotation, one of four possible ways to compute the rotation
	     angle is chosen */
    if(rotated_atom->pos[Z] < 0) angle = fabs(angle) - fabs(orig_angle);
    else angle = (-1) * (fabs(angle) - fabs(orig_angle));
  }

  /* It turned out when, because in the PDB the position of the
     atoms is given with up to 0.001 Angstrom accuracy, only floats
     are used for the coordinates, a rotation by less than 0.02 rad
     might not reflect a change in the original position. To avoid this
     (which leads to an infinite loop of minimal - but non effective -
     rotations), the bumping side chain is rotated at least by 0.02 rad
     which are 1.146 degrees */
  if(angle > 0 && angle < 0.02) angle = 0.02;
  if(angle < 0 && angle > -0.02) angle = -0.02;

  /* Copy the original positions rather than "move" back */ 
  memcpy(rotated_atom->pos, rot_atom_orig_pos, 3*sizeof(*rot_atom_orig_pos));
  memcpy(fixed_atom->pos, fix_atom_orig_pos, 3*sizeof(*fix_atom_orig_pos));

#if 0
  /* now transform the points back to their original position */
  transform_point_back(rotated_atom->pos, unbump_matrix);
  transform_point_back(fixed_atom->pos, unbump_matrix);
#endif
  return angle;
}

int compute_all_rotation_angles(global_data_pt global)
{
  int   fragment;
  int   i, j;

  for(i = 0; i < global->number_of_target_residues; i++ )
    for(j = 0; j < 5; j++ )
      global->unbump_target_indices[i][j] = UNKNOWN_TARGET_INDEX;
  for(i = 0; i < global->ligand->number_of_flexible_bonds; i++ )
    global->unbump_indices[i] = UNKNOWN_INDEX;

  global->number_of_unbump_entries = 0;
  for(i = 0; i < global->number_of_bumps; i++){
    /* only compute angle for non main-chain atoms, if target_bumps[i] is below
     * zero, the bumping target atom is a water, i.e., the collision cannot be 
     * resolved by a target rotation */
    if(global->target_bumps[i] >= 0 && 
       global->target_atoms[global->target_bumps[i]].type > 5 && 
       global->target_atoms[global->target_bumps[i]].level != HETATM )
      compute_target_angles(global, i);

    /* otherwise this is an intramolecular bump within the target, since in 
     * that case the ligand index is negative */
    if(global->ligand_bumps[i] >= 0){
      fragment = global->ligand->atoms[global->ligand_bumps[i]].fragment;
      /* if the bumping ligand atom is located on the base fragment, then there
       * is no way to resolve this collision by rotating a ligand_side_chain */
      if(global->ligand->fragment_locations[fragment] == OFF_ANCHOR )
        compute_ligand_angles(global, i);
    }
  }
  return SUCCESS;
}
