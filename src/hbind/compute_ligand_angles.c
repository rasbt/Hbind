#include <float.h>
#include <math.h>
#include <string.h>
#include <types.h>
#include <mymalloc.h>
#include "err_handle.h"
#include "compute_all_rotation_angles.h"
#include "trans_rotate.h"
#include "unbump_side_chains.h"
#include "dist_fun.h"
#include "intra_bump_check.h"

int reduce_ligand_bond(global_data_pt global, int flex_bond, int *fragment)
{
  molecule_pt  ligand;
  int          bond;
  
  ligand = global->ligand;
  bond = ligand->flexible_bonds[flex_bond];
  if ( *fragment == ligand->bonds[bond].fragment1 )
    {
      *fragment = ligand->bonds[bond].fragment2;
      if ( ligand->fragment_locations[*fragment] == ANCHOR )
	/* this was already the last bond, since the new fragment is part of
	   the anchor */
	return -1;
      flex_bond = ligand->way_to_anchor[*fragment];
    }
  else
    {
      *fragment = ligand->bonds[bond].fragment1;
      if ( ligand->fragment_locations[*fragment] == ANCHOR )
	/* this was already the last bond, since the new fragment is part of
	   the anchor */
	return -1;
      flex_bond = ligand->way_to_anchor[*fragment];
    }
  return flex_bond;
}

void  compute_ligand_matrix ( molecule_pt        ligand,
			      atom_pt            target_atom,
			      int                flex_bond,
			      transform_matrix_t matrix )
{
  int      bond_index,
           base,
           head;
  
  bond_index = ligand->flexible_bonds[flex_bond];
  if ( ligand->bond_directions[flex_bond] == REVERSE )
    {
      base = ligand->bonds[bond_index].atom1;
      head = ligand->bonds[bond_index].atom2;
    }
  else
    {
      base = ligand->bonds[bond_index].atom2;
      head = ligand->bonds[bond_index].atom1;
    }
  /* the rotatable bond is between atoms[base] and atoms[head] */
  compute_transformation_matrix ( target_atom->pos,
				  ligand->atoms[base].pos,
				  ligand->atoms[head].pos,
				  matrix );
}

void  compute_ligand_angles ( global_data_pt  global,
			      int             bump )
{
  molecule_pt         ligand;
  atom_pt             ligand_atom;
  atom_t              target_atom;
  transform_matrix_pt unbump_matrix;
  unbump_data_pt      unbump_data;
  double              angle;
  int                 *unbump_indices,
                      *unbump_bonds;
  int                 bond,
                      flex_bond,
                      fragment,
                      index,
                      entry,
                      result;
  int                 j;

  target_atom.pos = (float*) mymalloc(3*sizeof(float));    

  ligand = global->ligand;
  unbump_indices = global->unbump_indices;
  unbump_bonds = global->unbump_bonds;
  /* the bumping target atom is a water */
  if ( global->target_bumps[bump] < 0 ){
    j = ( (-1) * global->target_bumps[bump] ) - 1;
    memcpy(target_atom.pos, global->waters[j].pos, 3*sizeof(float));
    target_atom.rad = WATER_RAD;
  }else{
    j = global->target_bumps[bump];
    memcpy(target_atom.pos, global->target_atoms[j].pos, 3*sizeof(float));
    target_atom.rad = global->target_atoms[j].rad;
  }
  ligand_atom = &ligand->atoms[global->ligand_bumps[bump]];
  fragment = ligand_atom->fragment;
  flex_bond = ligand->way_to_anchor[ligand_atom->fragment];
  bond = ligand->flexible_bonds[flex_bond];
  if ( global->ligand_bumps[bump] == ligand->bonds[bond].atom1
       || global->ligand_bumps[bump] == ligand->bonds[bond].atom2 )
    /* otherwise we cannot reduce this bond, since it is already the
       last bond before the anchor fragment */
    flex_bond = reduce_ligand_bond ( global,
				     flex_bond,
				     &fragment );

  while ( flex_bond >= 0 )
    {
      if ( unbump_indices[flex_bond] == UNKNOWN_INDEX )
        /* we haven't had this bond yet */
	index = global->number_of_unbump_bonds;
      else
        index = unbump_indices[flex_bond];
      entry = global->number_of_unbump_entries;
      unbump_matrix = &global->unbump_matrices[entry];
      compute_ligand_matrix ( ligand,
			      &target_atom,
			      flex_bond,
			      *unbump_matrix );
      /* check, if a rotation around this bond can resolve the bump */
      result = check_resolvability ( target_atom.pos,
				     ligand_atom->pos,
				     target_atom.rad 
				     + ligand_atom->rad 
				     - global->side_chain_overlap 
				     + 0.0001,
				     *unbump_matrix );
      if ( result == NO_BUMP )
	{
	  angle = 
	    compute_rotation_angle ( &target_atom,
				     ligand_atom,
				     global->side_chain_overlap,
				     *unbump_matrix );
#ifdef TRACE
	  printf ( "ligand bond %d (flex %d): angle = %5.3f\n",
		   index,
		   flex_bond,
		   angle );
#endif
	  if ( angle != FLT_MAX
	       && check_for_intra_ligand_bump ( global,
						unbump_matrix,
						flex_bond,
						angle )
	       == NO_BUMP )
	    /* angle is equal to FLT_MAX, if there was no solution for the
	       set of equations used to compute the rotation angle */
	    {
#ifdef TRACE
	      printf ( "Intra bump passed\n" );
#endif
	      /* this bond is ok, so keep the matrix */
	      unbump_indices[flex_bond] = index;
	      unbump_bonds[index] = (-1) * ( flex_bond + 1);
	      /* now store all the data for this particular entry */
	      global->unbump_data[bump][index] 
		= &global->unbump_entries[entry];
	      global->number_of_unbump_entries++;
	      if ( index == global->number_of_unbump_bonds)
		/* this is a new bond */
		{
		  global->number_of_unbump_bonds++;
		  if ( global->number_of_unbump_bonds == MAX_UNBUMP_BONDS )
		    err_panic2 ( "compute_ligand_angles",
				"more than MAX_UNBUMP_BONDS");
		}
	      unbump_data = global->unbump_data[bump][index];
	      unbump_data->matrix = unbump_matrix;
	      unbump_data->angle = angle;
	      unbump_data->penalty = 1.0;
	      unbump_data->old_distance = dist_fun ( ligand_atom->pos,
						     target_atom.pos );
	    }
#ifdef TRACE
      else printf ( "Intra bump failed\n" );
#endif
    }
    flex_bond = reduce_ligand_bond ( global, flex_bond, &fragment );
  }
  free(target_atom.pos);
}
