#include <stdio.h>
#include "defs.h"
#include "types.h"
#include "err_handle.h"
#include "rotatable_bonds.h"

#define TRAC

int  find_dependent_target_angles ( global_data_pt  global )
{
  int             **unbump_target_indices,
                  **dependents;
  int             *dependencies,
                  *unbump_bonds,
                  *number_of_dependents;
  int             dep[MAX_UNBUMP_BONDS];
  int             res_index,
                  res_type,
                  bond_index,
                  counter,
                  global_counter;
  int             i, j;

  unbump_target_indices = global->unbump_target_indices;
  unbump_bonds = global->unbump_bonds;
  dependencies = global->unbump_dependencies;
  dependents = global->unbump_dependents;
  number_of_dependents = global->number_of_unbump_dependents;
  global_counter = 0;
  for ( i = 0; i < global->number_of_unbump_bonds; i++ )
    {
      counter = 0;
      if ( unbump_bonds[i] >= 0 )
	/* this is a target bond */
	{
	  res_index = unbump_bonds[i] / 10;
	  res_type = global->target_residues[res_index].type;
	  bond_index = unbump_bonds[i] % 10;
	  for ( j = bond_index + 1; 
		j <= number_of_rotatable_bonds[res_type]; 
		j++ )
	    if ( unbump_target_indices[res_index][j] != UNKNOWN_TARGET_INDEX )
	      /* this bond can be used to resolve a bump, so it is one
		 of the bonds that depend on he current one */
	      dep[counter++] = unbump_target_indices[res_index][j];
	}
      if ( counter > 0 )
	{
	  dependents[i] = &dependencies[global_counter];
	  for ( j = 0; j < counter; j++ ) 
	    dependents[i][j] = dep[j];
	  global_counter += counter;
	  number_of_dependents[i] = counter;
	}
      else
	number_of_dependents[i] = 0;
      if ( global_counter >= MAX_UNBUMP_BONDS * MAX_UNBUMP_DEPENDENCIES )
	err_panic2 ( "find_dependent_angles", "unbump_dependencies overflow");
    }
  return global_counter;
}

void  find_dependent_ligand_angles ( global_data_pt  global,
				     int             global_counter )
{
  molecule_pt ligand;
  bond_pt     bonds;
  int         **dependents;
  int         *dependencies,
              *unbump_indices,
              *unbump_bonds,
              *way_to_anchor,
              *anchor_dist,
              *number_of_dependents;
  int         dep[MAX_NUMBER_OF_FLEXIBLE_BONDS][MAX_NUMBER_OF_FLEXIBLE_BONDS];
  int         counter[MAX_NUMBER_OF_FLEXIBLE_BONDS];
  int         fragment,
              bond_index,
              index,
              flex_bond;
  int         i, j;

  unbump_indices = global->unbump_indices;
  unbump_bonds = global->unbump_bonds;
  ligand = global->ligand;
  bonds = ligand->bonds;
  dependencies = global->unbump_dependencies;
  dependents = global->unbump_dependents;
  number_of_dependents = global->number_of_unbump_dependents;
  way_to_anchor = ligand->way_to_anchor;
  anchor_dist = ligand->anchor_dist;
  for ( i = 0; i < global->number_of_unbump_bonds; i++ )
    counter[i] = 0;
  for ( i = 0; i < ligand->number_of_flexible_bonds; i++ )
    if ( unbump_indices[i] != UNKNOWN_INDEX )
      {
	bond_index = ligand->flexible_bonds[i];
	/* first check, which fragment is directed away from the anchor 
	   fragment */
	if ( way_to_anchor[bonds[bond_index].fragment1] == i )
	  /* fragment2 is located towards the anchor, so this is the right
	     way to go */
	  fragment = bonds[bond_index].fragment2;
	else
	  fragment = bonds[bond_index].fragment1;
#ifdef TRACE
	printf ( "starting with flex bond %d -> ", i );
#endif
	/*
	while ( way_to_anchor[fragment] != UNKNOWN )
	*/
	while ( anchor_dist[fragment] > 0 )
	  {
	    flex_bond = way_to_anchor[fragment];
#ifdef TRACE
	    printf ( "%d ", flex_bond );
#endif
	    index = unbump_indices[flex_bond];
	    if ( index != UNKNOWN_INDEX )	      
	      dep[index][(counter[index])++] = unbump_indices[i];
	    bond_index = ligand->flexible_bonds[flex_bond];
	    /* first check, which fragment is directed away from the anchor 
	       fragment */
	    if ( way_to_anchor[bonds[bond_index].fragment1] == flex_bond )
	      /* fragment2 is located away from the anchor, so that's the right
		 way to go */
	      fragment = bonds[bond_index].fragment2;
	    else
	      fragment = bonds[bond_index].fragment1;
	  }	    
#ifdef TRACE
	printf ( "\n" );
#endif
      }
  for ( i = 0; i < global->number_of_unbump_bonds; i++ )
    if ( unbump_bonds[i] < 0 )
      /* the current bond is a ligand bond */
      {
	if ( counter[i] > 0 )
	  {
	    dependents[i] = &dependencies[global_counter];
	    for ( j = 0; j < counter[i]; j++ )
	      dependents[i][j] = dep[i][j];
	    global_counter += counter[i];
	    number_of_dependents[i] = counter[i];
	  }
	else
	  number_of_dependents[i] = 0;
      }
  if ( global_counter >= MAX_UNBUMP_BONDS * MAX_UNBUMP_DEPENDENCIES )
    err_panic2 ( "find_dependent_angles", "unbump_dependencies overflow");
}

void  compute_unbump_dependencies ( global_data_pt  global )
{
  int  counter;
  counter = find_dependent_target_angles ( global );
  find_dependent_ligand_angles ( global, counter );

#ifdef TRACE
  int i, j;
  for(i = 0; i < global->number_of_unbump_bonds; i++){
    printf ( "dependents of bond %d: ", i );
    for ( j = 0; j < global->number_of_unbump_dependents[i]; j++ )
      printf ( "%d ", global->unbump_dependents[i][j] );
    printf ( "\n" );
  }
#endif
}
	      
