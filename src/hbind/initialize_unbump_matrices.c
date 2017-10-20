#include <stdio.h>
#include "defs.h"
#include "types.h"
#include "bitstrings.h"

#define TRAC

void  initialize_probability_matrix ( global_data_pt  global )
{
  double         probability;
  int            count;
  int            b, r;
  unbump_data_pt **unbump_data = global->unbump_data;
  for(b = 0; b < global->number_of_bumps; b++){
    count = 0;
    for(r = 0; r < global->number_of_unbump_bonds; r++)
      if(unbump_data[b][r]) count++;
    if(!count) continue;

    probability = 1.0 / count;
    for(r = 0; r < global->number_of_unbump_bonds; r++)
      if(unbump_data[b][r]) unbump_data[b][r]->probability = probability;
  }
  
#ifdef TRACE
  for ( r = 0; r < global->number_of_unbump_bonds; r++){
    for ( b = 0; b < global->number_of_bumps; b++ ){
      if ( unbump_data[b][r] != NULL )
        printf ( " %6.3lf", unbump_data[b][r]->probability );
      else printf ( "  *.***" );
    }
    printf ( "\n" );
  }
#endif
}
      
int  number_of_rotating_atoms ( global_data_pt  global,
				int             index )
{
  atom_pt  target,
           atoms;
  int      *unbump_bonds;
  int      res_index,
           bond,
           flex_bond,
           offset, 
           count;
  int      i;

  target = global->target_atoms;
  atoms = global->ligand->atoms;
  unbump_bonds = global->unbump_bonds;       
  count = 0;
  if ( unbump_bonds[index] >= 0 )
    /* this is a target bond */
    {
      res_index = unbump_bonds[index] / 10;
      bond = unbump_bonds[index] % 10;
      offset = global->target_residues[res_index].start_atom;
      for ( i = 0; 
	    i < global->target_residues[res_index].number_of_atoms;
	    i++ )
	if ( target[offset+i].level > bond )
	  count++;
    }
  else
    /* this is a ligand bond */
    {      
      flex_bond = ( (-1) * unbump_bonds[index] ) - 1;
      if ( global->ligand->bond_directions[flex_bond] == REVERSE )
	{
	  for ( i = 0; i < global->ligand->number_of_atoms; i++ )
	    if ( bitstring_get_bit ( atoms[i].fragments,
				     flex_bond ) 
	     	 && atoms[i].type != H )
	      /* this atom is located 'right' from the rotation bond, 
		 i.e., at the side with the larger atom index, which is 
		 the side pointing away from the base */
	      count++;
	}
      else
	{
	  for ( i = 0; i < global->ligand->number_of_atoms; i++ )
	    if ( bitstring_get_bit ( atoms[i].fragments,
				     flex_bond ) == 0 
		 && atoms[i].type != H ) 
	      /* this atom is located 'left' from the rotation bond */
	      count++;
	}      
    }
  return count;
}

void initialize_force_matrix(global_data_pt  global)
{
  int rotated_atoms;
  int b, r;
  unbump_data_pt  **unbump_data = global->unbump_data;
  
  for(r = 0; r < global->number_of_unbump_bonds; r++){
    rotated_atoms = number_of_rotating_atoms(global, r);
    for ( b = 0; b < global->number_of_bumps; b++ ){
      if(!unbump_data[b][r]) continue;

      if(unbump_data[b][r]->angle > 0)
        unbump_data[b][r]->force = (double) rotated_atoms * 
          unbump_data[b][r]->angle * unbump_data[b][r]->penalty;
      else
        unbump_data[b][r]->force = (double) rotated_atoms * 
          (-1.0) * unbump_data[b][r]->angle * unbump_data[b][r]->penalty;
    }
  }
#ifdef TRACE
  for(r = 0; r < global->number_of_unbump_bonds; r++){
    for(b = 0; b < global->number_of_bumps; b++){
      if(unbump_data[b][r] != 0) printf(" %6.3lf", unbump_data[b][r]->force);
      else printf ( "     **" );
    }
    printf ( "\n" );
  }
#endif
}

void initialize_unbump_matrices(global_data_pt global)
{
  initialize_probability_matrix(global);
  initialize_force_matrix(global);
}

