#include <stdio.h>
#include <float.h>
#include <math.h>
#include "types.h"
#include "err_handle.h"
#include "rotatable_bonds.h"
#include "unbump_side_chains.h"
#include "compute_all_rotation_angles.h"
#include "dist_fun.h"
#include "trans_rotate.h"
#include "rotate_unbump_bonds.h"
#include "intra_bump_check.h"


void  compute_target_angles ( global_data_pt  global,
			      int             bump )
{
  atom_pt             target_atom,
                      ligand_atom,
                      atoms;
  residue_pt          residue;
  transform_matrix_pt unbump_matrix;
  unbump_data_pt      unbump_data;
  int                 **unbump_target_indices;
  int                 atoms_to_rotate[MAX_RESIDUE_ATOMS];
  int                 *unbump_bonds;
  double              angle;
  float               intra_bump_penalty;
  int                 res_type,
                      res_index,
                      bond,
                      index,
                      entry,
                      number_of_atoms_to_rotate,
                      result;
  int                 j;
       
  atoms = global->target_atoms;
  target_atom = &global->target_atoms[global->target_bumps[bump]];
  /* this is an intermolecular bump */
  if(global->ligand_bumps[bump] >= 0)
    ligand_atom = &global->ligand->atoms[global->ligand_bumps[bump]];
  /* this is an intra-target bump, so, the other bumping atom is a target 
     atom, too, however, we are using the ligand_atom variable here, it
     shouldn't make a difference (just cause some confusion, if somebody
     else has to work with this code) */
  else ligand_atom = &global->target_atoms[(-1)*(global->ligand_bumps[bump]+1)];

  unbump_target_indices = global->unbump_target_indices;
  unbump_bonds = global->unbump_bonds;
  res_type = target_atom->res;
  res_index = target_atom->residue_index;
  residue = &global->target_residues[res_index];
  if ( target_atom->level - 1 >  number_of_rotatable_bonds[res_type] )
    bond = number_of_rotatable_bonds[res_type];
  else
    bond = target_atom->level - 1; 

  /* the bumping atom is adjacent to the rotatable bond, so skip this bond */
  if ( rotatable_bonds[res_type][bond] == target_atom->type )
    fprintf(stderr, "******** THIS SHOULDN'T HAPPEN ********\n");

  while(bond > 0){ 
    /* we haven't had this bond yet */
    if(unbump_target_indices[res_index][bond] == UNKNOWN_TARGET_INDEX )
      index = global->number_of_unbump_bonds;
    else index = unbump_target_indices[res_index][bond];

    entry = global->number_of_unbump_entries;
    unbump_matrix = &global->unbump_matrices[entry];
    compute_matrix(&global->target_atoms[residue->start_atom], ligand_atom, 
                   res_type, bond, *unbump_matrix);
      
    /* check, if a rotation around this bond can resolve the bump */
    result = check_resolvability(ligand_atom->pos, target_atom->pos, 
                                 ligand_atom->rad + target_atom->rad 
                                 - global->side_chain_overlap + 0.0001, 
                                 *unbump_matrix );
    if(result != BUMP){ 
      angle = compute_rotation_angle(ligand_atom, target_atom, 
                                     global->side_chain_overlap, 
                                     *unbump_matrix); 
      /* otherwise there was no solution to the equation that
         computes this angle */
      if(angle != FLT_MAX){
        unbump_bonds[index] = res_index * 10 + bond;
        number_of_atoms_to_rotate = 
          compute_atoms_to_rotate(global, index, atoms_to_rotate);

#ifdef TRACE
        printf("target bond %d (%d,%d): angle = %5.3f, atoms = %d\n", index, 
               res_index, bond, angle, number_of_atoms_to_rotate);
#endif
        for ( j = 0; j < number_of_atoms_to_rotate; j++){ 
          transform_point(atoms[atoms_to_rotate[j]].pos, *unbump_matrix); 
          rotate_single_atom_around_y_axis(atoms[atoms_to_rotate[j]].pos,
                                           angle); 
          transform_point_back(atoms[atoms_to_rotate[j]].pos, *unbump_matrix); 
        } 
        intra_bump_penalty = bump_check_neighbors(atoms, atoms_to_rotate, 
                                                  number_of_atoms_to_rotate, 
                                                  global->intra_overlap);
#ifdef TRACE
        printf ( "intra_bump_penalty = %5.3lf\n", intra_bump_penalty );
#endif
	      for ( j = 0; j < number_of_atoms_to_rotate; j++ )
		{
		  transform_point ( atoms[atoms_to_rotate[j]].pos,
				    *unbump_matrix );
		  rotate_single_atom_around_y_axis 
		    ( atoms[atoms_to_rotate[j]].pos,
		      (-1) * angle );
		  transform_point_back ( atoms[atoms_to_rotate[j]].pos,
					 *unbump_matrix );
		}
	      if ( result == NO_BUMP 
		   && intra_bump_penalty < global->intermediate_overlap )
		{
		  /* this bond is ok, so keep the matrix */
		  /* store the index for this bond */
		  unbump_target_indices[res_index][bond] = index;
		  unbump_bonds[index] = res_index * 10 + bond;
		  /* now store all the data for this particular entry */
		  global->unbump_data[bump][index] = 
		    &global->unbump_entries[entry];
		  global->number_of_unbump_entries++;
		  if ( index == global->number_of_unbump_bonds )
		    /* this is a new bond */
		    {
		      global->number_of_unbump_bonds++;
		      if ( global->number_of_unbump_bonds == MAX_UNBUMP_BONDS )
			err_panic2 ( "find_target_rotatable_bonds",
				    "more than MAX_UNBUMP_BONDS");
		    }
		  unbump_data = global->unbump_data[bump][index];
		  unbump_data->matrix = unbump_matrix;
		  unbump_data->angle = angle;
		  unbump_data->penalty = 1.0 + intra_bump_penalty;
		  unbump_data->old_distance = dist_fun ( ligand_atom->pos,
							 target_atom->pos );
		}
	    }
	}
      bond--;
    }
}

